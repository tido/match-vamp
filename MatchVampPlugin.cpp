/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    This file copyright 2007 Simon Dixon, Chris Cannam and QMUL.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "MatchVampPlugin.h"

#include "Matcher.h"
#include "MatchFeeder.h"
#include "Path.h"

#include <vamp/vamp.h>
#include <vamp-sdk/PluginAdapter.h>
#include <vamp-sdk/RealTime.h>

#include <vector>
#include <algorithm>

//static int extant = 0;

#ifdef _WIN32
HANDLE
MatchVampPlugin::m_serialisingMutex;
#else
pthread_mutex_t 
MatchVampPlugin::m_serialisingMutex;
#endif

bool
MatchVampPlugin::m_serialisingMutexInitialised = false;

// We want to ensure our freq map / crossover bin in Matcher.cpp are
// always valid with a fixed FFT length in seconds, so must reject low
// sample rates
static float sampleRateMin = 5000.f;

static float defaultStepTime = 0.020;

MatchVampPlugin::MatchVampPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(0),
    m_stepTime(defaultStepTime),
    m_blockSize(0),
    m_serialise(false),
    m_begin(true),
    m_locked(false)
{
    if (inputSampleRate < sampleRateMin) {
        std::cerr << "MatchVampPlugin::MatchVampPlugin: input sample rate "
                  << inputSampleRate << " < min supported rate "
                  << sampleRateMin << ", plugin will refuse to initialise"
                  << std::endl;
    }

    if (!m_serialisingMutexInitialised) {
        m_serialisingMutexInitialised = true;
#ifdef _WIN32
        m_serialisingMutex = CreateMutex(NULL, FALSE, NULL);
#else
        pthread_mutex_init(&m_serialisingMutex, 0);
#endif
    }

    pm1 = 0;
    pm2 = 0;
    feeder = 0;
//    std::cerr << "MatchVampPlugin::MatchVampPlugin(" << this << "): extant = " << ++extant << std::endl;
}

MatchVampPlugin::~MatchVampPlugin()
{
//    std::cerr << "MatchVampPlugin::~MatchVampPlugin(" << this << "): extant = " << --extant << std::endl;

    delete feeder;
    delete pm1;
    delete pm2;

    if (m_locked) {
#ifdef _WIN32
        ReleaseMutex(m_serialisingMutex);
#else
        pthread_mutex_unlock(&m_serialisingMutex);
#endif
        m_locked = false;
    }
}

string
MatchVampPlugin::getIdentifier() const
{
    return "match";
}

string
MatchVampPlugin::getName() const
{
    return "Match Performance Aligner";
}

string
MatchVampPlugin::getDescription() const
{
    return "Calculate alignment between two performances in separate channel inputs";
}

string
MatchVampPlugin::getMaker() const
{
    return "Simon Dixon (plugin by Chris Cannam)";
}

int
MatchVampPlugin::getPluginVersion() const
{
    return 1;
}

string
MatchVampPlugin::getCopyright() const
{
    return "GPL";
}

MatchVampPlugin::ParameterList
MatchVampPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor desc;
    desc.identifier = "serialise";
    desc.name = "Serialise Plugin Invocations";
    desc.description = "Reduce potential memory load at the expense of multiprocessor performance by serialising multi-threaded plugin runs";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    return list;
}

float
MatchVampPlugin::getParameter(std::string name) const
{
    if (name == "serialise") {
        return m_serialise ? 1.0 : 0.0; 
    }
    return 0.0;
}

void
MatchVampPlugin::setParameter(std::string name, float value)
{
    if (name == "serialise") {
        m_serialise = (value > 0.5);
//        std::cerr << "MatchVampPlugin::setParameter: set serialise to " << m_serialise << std::endl;
    }
}

size_t
MatchVampPlugin::getPreferredStepSize() const
{
    return m_inputSampleRate * defaultStepTime;
}

size_t
MatchVampPlugin::getPreferredBlockSize() const
{
    return 2048;
}

void
MatchVampPlugin::createMatchers() const
{
    Matcher::Parameters params(m_inputSampleRate, m_stepTime, m_blockSize);
    pm1 = new Matcher(params, 0);
    pm2 = new Matcher(params, pm1);
    pm1->setOtherMatcher(pm2);
    feeder = new MatchFeeder(pm1, pm2);
}

bool
MatchVampPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (m_inputSampleRate < sampleRateMin) {
        std::cerr << "MatchVampPlugin::MatchVampPlugin: input sample rate "
                  << m_inputSampleRate << " < min supported rate "
                  << sampleRateMin << std::endl;
        return false;
    }
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;
    if (stepSize > blockSize/2 ||
        blockSize != getPreferredBlockSize()) return false;

    m_stepSize = stepSize;
    m_stepTime = float(stepSize) / m_inputSampleRate;
    m_blockSize = blockSize;

    cerr << "step size = " << m_stepSize << ", time = " << m_stepTime << endl;

    createMatchers();
    m_begin = true;
    m_locked = false;

    return true;
}

void
MatchVampPlugin::reset()
{
    delete feeder;
    delete pm1;
    delete pm2;
    feeder = 0;
    pm1 = 0;
    pm2 = 0;

    createMatchers();
    m_begin = true;
    m_locked = false;
}

MatchVampPlugin::OutputList
MatchVampPlugin::getOutputDescriptors() const
{
    OutputList list;

    float outRate = 1.0 / m_stepTime;

    OutputDescriptor desc;
    desc.identifier = "path";
    desc.name = "Path";
    desc.description = "Alignment path";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    list.push_back(desc);

    desc.identifier = "a_b";
    desc.name = "A-B Timeline";
    desc.description = "Timing in performance B corresponding to moments in performance A";
    desc.unit = "sec";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    list.push_back(desc);

    desc.identifier = "b_a";
    desc.name = "B-A Timeline";
    desc.description = "Timing in performance A corresponding to moments in performance B";
    desc.unit = "sec";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    list.push_back(desc);

    desc.identifier = "a_b_divergence";
    desc.name = "A-B Divergence";
    desc.description = "Difference between timings in performances A and B";
    desc.unit = "sec";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    list.push_back(desc);

    desc.identifier = "a_b_temporatio";
    desc.name = "A-B Tempo Ratio";
    desc.description = "Ratio of tempi between performances A and B";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    list.push_back(desc);

    desc.identifier = "a_features";
    desc.name = "A Features";
    desc.description = "Spectral features extracted from performance A";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    list.push_back(desc);

    return list;
}

MatchVampPlugin::FeatureSet
MatchVampPlugin::process(const float *const *inputBuffers,
                         Vamp::RealTime timestamp)
{
    if (m_begin) {
        if (!m_locked && m_serialise) {
            m_locked = true;
#ifdef _WIN32
            WaitForSingleObject(m_serialisingMutex, INFINITE);
#else
            pthread_mutex_lock(&m_serialisingMutex);
#endif
        }
        m_startTime = timestamp;
        m_begin = false;
    }
    
//    std::cerr << timestamp.toString();

    feeder->feed(inputBuffers);

//    std::cerr << ".";
//    std::cerr << std::endl;

    return FeatureSet();
}

MatchVampPlugin::FeatureSet
MatchVampPlugin::getRemainingFeatures()
{
    int x = pm2->getFrameCount() - 1;
    int y = pm1->getFrameCount() - 1;

    Finder *finder = feeder->getFinder();

    std::vector<int> pathx;
    std::vector<int> pathy;

//    std::cerr << "initial x,y = " << x << std::endl;

    while (finder->find(y, x) && ((x > 0) || (y > 0))) {

        pathx.push_back(x);
        pathy.push_back(y);

//        std::cerr << pathx.size() << ": (" << x << "," << y << ")" << std::endl;

        switch (finder->getDistance() & ADVANCE_BOTH){
        case ADVANCE_THIS:  y--; break;
        case ADVANCE_OTHER: x--; break;
        case ADVANCE_BOTH:  x--; y--; break;
        default: // this would indicate a bug, but we wouldn't want to hang
            std::cerr << "WARNING: MatchVampPlugin::getRemainingFeatures: Neither matcher advanced in path backtrack at (" << x << "," << y << ")" << std::endl;
            if (x > y) x--; else y--; break;
        }
    }

    std::reverse(pathx.begin(), pathx.end());
    std::reverse(pathy.begin(), pathy.end());

    int smoothedLen = Path().smooth(pathx, pathy, pathx.size());

    FeatureSet returnFeatures;

    int prevx = 0;
    int prevy = 0;

    for (int i = 0; i < smoothedLen; ++i) {

        int x = pathx[i];
        int y = pathy[i];

        Vamp::RealTime xt = Vamp::RealTime::frame2RealTime
            (x * m_stepSize, lrintf(m_inputSampleRate));
        Vamp::RealTime yt = Vamp::RealTime::frame2RealTime
            (y * m_stepSize, lrintf(m_inputSampleRate));

        Feature feature;
        feature.hasTimestamp = true;
        feature.timestamp = m_startTime + xt;
        feature.values.clear();
        feature.values.push_back(yt.sec + double(yt.nsec)/1.0e9);
        returnFeatures[0].push_back(feature);
        
        if (x != prevx) {

            feature.hasTimestamp = true;
            feature.timestamp = m_startTime + xt;
            feature.values.clear();
            feature.values.push_back(yt.sec + yt.msec()/1000.0);
            returnFeatures[1].push_back(feature);

            Vamp::RealTime diff = yt - xt;
            feature.values.clear();
            feature.values.push_back(diff.sec + diff.msec()/1000.0);
            returnFeatures[3].push_back(feature);

            if (i > 0) {
                int lookback = 100; //!!! arbitrary
                if (lookback > i) lookback = i;
                int xdiff = x - pathx[i-lookback];
                int ydiff = y - pathy[i-lookback];
                if (xdiff != 0 && ydiff != 0) {
                    float ratio = float(ydiff)/float(xdiff);
                    if (ratio < 8 && ratio > (1.0/8)) { //!!! just for now, since we aren't dealing properly with silence yet
                        feature.values.clear();
                        feature.values.push_back(ratio);
                        returnFeatures[4].push_back(feature);
                    }
                }
            }
        }

        if (y != prevy) {
            feature.hasTimestamp = true;
            feature.timestamp = m_startTime + yt;
            feature.values.clear();
            feature.values.push_back(xt.sec + xt.msec()/1000.0);
            returnFeatures[2].push_back(feature);
        }

        prevx = x;
        prevy = y;
    }

    delete feeder;
    delete pm1;
    delete pm2;
    feeder = 0;
    pm1 = 0;
    pm2 = 0;

    if (m_locked) {
#ifdef _WIN32
        ReleaseMutex(m_serialisingMutex);
#else
        pthread_mutex_unlock(&m_serialisingMutex);
#endif
        m_locked = false;
    }

    return returnFeatures;
    

/*
    for (int i = 0; i < smoothedLen; ++i) {
        std::cerr << i << ": [" << pathx[i] << "," << pathy[i] << "]" << std::endl;
    }

    std::cerr << std::endl;
    std::cerr << "File: A" << std::endl;
    std::cerr << "Marks: -1" << std::endl;
    std::cerr << "FixedPoints: true 0" << std::endl;
    std::cerr << "0" << std::endl;
    std::cerr << "0" << std::endl;
    std::cerr << "0" << std::endl;
    std::cerr << "0" << std::endl;
    std::cerr << "File: B" << std::endl;
    std::cerr << "Marks: 0" << std::endl;
    std::cerr << "FixedPoints: true 0" << std::endl;
    std::cerr << "0.02" << std::endl;
    std::cerr << "0.02" << std::endl;

    std::cerr << smoothedLen << std::endl;
    for (int i = 0; i < smoothedLen; ++i) {
        std::cerr << pathx[i] << std::endl;
    }

    std::cerr << smoothedLen << std::endl;
    for (int i = 0; i < smoothedLen; ++i) {
        std::cerr << pathy[i] << std::endl;
    }
*/
}

static Vamp::PluginAdapter<MatchVampPlugin> mvpAdapter;

const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int version,
                                                    unsigned int index)
{
    if (version < 1) return 0;

    switch (index) {
    case  0: return mvpAdapter.getDescriptor();
    default: return 0;
    }
}
