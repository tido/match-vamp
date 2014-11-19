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

#include "Finder.h"

#include "Path.h"

#include <algorithm>

using namespace std;

Finder::Finder(Matcher *pm)
{
    m_m = pm;
    m_duration1 = -1;
    m_duration2 = -1;
} // constructor

Finder::~Finder()
{
}

void
Finder::setDurations(int d1, int d2)
{
    m_duration1 = d1;
    m_duration2 = d2;
}

Matcher::Advance
Finder::getExpandDirection(int row, int col)
{
    return getExpandDirection(row, col, false);
} // getExpandDirection()

Matcher::Advance
Finder::getExpandDirection(int row, int col, bool check)
{
    double min = m_m->getPathCost(row, col);
    
    int bestRow = row;
    int bestCol = col;

    pair<int, int> rowRange = m_m->getRowRange(col);
    if (rowRange.second > row+1) {
        rowRange.second = row+1;	// don't cheat by looking at future :)
    }
    for (int index = rowRange.first; index < rowRange.second; index++) {
        double tmp = m_m->getPathCost(index, col);
        if (tmp < min) {
            min = tmp;
            bestRow = index;
        }
    }

    pair<int, int> colRange = m_m->getColRange(row);
    if (colRange.second > col+1) {
        colRange.second = col+1;	// don't cheat by looking at future :)
    }
    for (int index = colRange.first; index < colRange.second; index++) {
        double tmp = m_m->getPathCost(row, index);
        if (tmp < min) {
            min = tmp;
            bestCol = index;
            bestRow = row;
        }
    }

    if (bestRow == row) {
        if (bestCol == col) {
            return Matcher::AdvanceBoth;
        } else {
            return Matcher::AdvanceThis;
        }
    } else if (bestCol == col) {
        return Matcher::AdvanceOther;
    } else {
        return Matcher::AdvanceNone;
    }

}

void
Finder::recalculatePathCostMatrix(int r1, int c1, int r2, int c2) 
{
    float diagonalWeight = sqrtf(2.f);
    
    int prevRowStart = 0, prevRowStop = 0;

    for (int r = r1; r <= r2; r++) {

        pair<int, int> colRange = m_m->getColRange(r);

        int rowStart = max(c1, colRange.first);
        int rowStop = min(c2 + 1, colRange.second);
        
        for (int c = rowStart; c < rowStop; c++) {

            float newCost = m_m->getDistance(r, c);
            Matcher::Advance dir = Matcher::AdvanceNone;

            if (r > r1) {	// not first row
                double min = -1;
                if ((c > prevRowStart) && (c <= prevRowStop)) {
                    // diagonal from (r-1,c-1)
                    min = m_m->getPathCost(r-1, c-1) + newCost * diagonalWeight;
                    dir = Matcher::AdvanceBoth;
                }
                if ((c >= prevRowStart) && (c < prevRowStop)) {
                    // vertical from (r-1,c)
                    double cost = m_m->getPathCost(r-1, c) + newCost;
                    if ((min < 0) || (cost < min)) {
                        min = cost;
                        dir = Matcher::AdvanceThis;
                    }
                }
                if (c > rowStart) {
                    // horizontal from (r,c-1)
                    double cost = m_m->getPathCost(r, c-1) + newCost;
                    if ((min < 0) || (cost < min)) {
                        min = cost;
                        dir = Matcher::AdvanceOther;
                    }
                }
                
                m_m->setPathCost(r, c, dir, min);

            } else if (c > rowStart) {	// first row
                // horizontal from (r,c-1)
                m_m->setPathCost(r, c, Matcher::AdvanceOther,
                                   m_m->getPathCost(r, c-1) + newCost);
            }
        }

        prevRowStart = rowStart;
        prevRowStop = rowStop;
    }
} 

int
Finder::retrievePath(bool smooth, vector<int> &pathx, vector<int> &pathy)
{
    pathx.clear();
    pathy.clear();

    int ex = m_m->getOtherFrameCount() - 1;
    int ey = m_m->getFrameCount() - 1;

    if (ex < 0 || ey < 0) {
        return 0;
    }
    
    int x = ex;
    int y = ey;

    if (m_duration2 > 0 && m_duration2 < m_m->getOtherFrameCount()) {
        x = m_duration2 - 1;
    }
    if (m_duration1 > 0 && m_duration1 < m_m->getFrameCount()) {
        y = m_duration1 - 1;
    }
    
//    cerr << "before: x = " << x << ", y = " << y << endl;

    if (!m_m->isAvailable(y, x)) {
        // Path did not pass through the expected end point --
        // probably means the pieces are substantially different in
        // the later bits. Reset the expected end point to the end of
        // both files including any trailing silence.
        cerr << "NOTE: Path did not pass through expected end point, inputs are probably significantly different" << endl;
        x = ex;
        y = ey;
    }

    recalculatePathCostMatrix(0, 0, y, x);

//    cerr << "start: x = " << x << ", y = " << y << endl;
    
    while (m_m->isAvailable(y, x) && (x > 0 || y > 0)) {

//        cerr << "x = " << x << ", y = " << y;
        
        pathx.push_back(x);
        pathy.push_back(y);

        switch (m_m->getAdvance(y, x)) {
        case Matcher::AdvanceThis:
//            cerr << ", going down (dist = " << getDistance() << ")" << endl;
            y--;
            break;
        case Matcher::AdvanceOther:
//            cerr << ", going left (dist = " << getDistance() << ")" << endl;
            x--;
            break;
        case Matcher::AdvanceBoth:
//            cerr << ", going diag (dist = " << getDistance() << ")" << endl;
            x--;
            y--;
            break;
        case Matcher::AdvanceNone: // this would indicate a bug, but we wouldn't want to hang
            cerr << "WARNING: Neither matcher advanced in path backtrack at (" << x << "," << y << ")" << endl;
            if (x > y) {
                x--;
            } else {
                y--;
            }
            break;
        }
    }

    if (x > 0 || y > 0) {
        cerr << "WARNING: Ran out of available path at (" << y << "," << x
             << ")!" << endl;
    }
    
    reverse(pathx.begin(), pathx.end());
    reverse(pathy.begin(), pathy.end());

    if (smooth) {
        int smoothedLen = Path().smooth(pathx, pathy, pathx.size());
        return smoothedLen;
    } else {
        return pathx.size();
    }
}


