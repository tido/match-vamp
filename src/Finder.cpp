/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright (c) 2007-2015 Simon Dixon, Chris Cannam, and Queen Mary
    University of London, Copyright (c) 2014-2015 Tido GmbH.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

//#define DEBUG_FINDER 1
//#define PERFORM_ERROR_CHECKS 1

#include "Finder.h"

#include "Path.h"

#include <algorithm>
#include <iomanip>

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
Finder::setMatcher(Matcher *pm)
{
    cerr << "Finder::setMatcher: finder " << this << ", matcher " << pm << endl;
    m_m = pm;
}

void
Finder::setDurations(int d1, int d2)
{
#ifdef DEBUG_FINDER
    cerr << "*** setDurations: " << d1 << ", " << d2 << endl;
#endif
    m_duration1 = d1;
    m_duration2 = d2;
}

bool
Finder::getBestRowCost(int row, int &bestCol, normpathcost_t &min)
{
    if (!m_m->isRowAvailable(row)) return false;
    pair<int, int> colRange = m_m->getColRangeForRow(row);
    if (colRange.first >= colRange.second) return false;
    for (int index = colRange.first; index < colRange.second; index++) {
        normpathcost_t tmp = m_m->getNormalisedPathCost(row, index);
        if (index == colRange.first || tmp < min) {
            min = tmp;
            bestCol = index;
        }
    }
    return true;
}    

bool
Finder::getBestColCost(int col, int &bestRow, normpathcost_t &min)
{
    if (!m_m->isColAvailable(col)) return false;
    pair<int, int> rowRange = m_m->getRowRangeForCol(col);
    if (rowRange.first >= rowRange.second) return false;
    bestRow = rowRange.first;
    for (int index = rowRange.first; index < rowRange.second; index++) {
        normpathcost_t tmp = m_m->getNormalisedPathCost(index, col);
        if (index == rowRange.first) {
            min = tmp;
        } else if (tmp < min) {
            min = tmp;
            bestRow = index;
        }
    }
    return true;
}    

void
Finder::getBestEdgeCost(int row, int col,
                        int &bestRow, int &bestCol,
                        normpathcost_t &min)
{
    min = m_m->getNormalisedPathCost(row, col);
    
    bestRow = row;
    bestCol = col;

    pair<int, int> rowRange = m_m->getRowRangeForCol(col);
    if (rowRange.second > row+1) {
        rowRange.second = row+1;	// don't cheat by looking at future :)
    }
    for (int index = rowRange.first; index < rowRange.second; index++) {
        normpathcost_t tmp = m_m->getNormalisedPathCost(index, col);
        if (tmp < min) {
            min = tmp;
            bestRow = index;
        }
    }

    pair<int, int> colRange = m_m->getColRangeForRow(row);
    if (colRange.second > col+1) {
        colRange.second = col+1;	// don't cheat by looking at future :)
    }
    for (int index = colRange.first; index < colRange.second; index++) {
        normpathcost_t tmp = m_m->getNormalisedPathCost(row, index);
        if (tmp < min) {
            min = tmp;
            bestCol = index;
            bestRow = row;
        }
    }
}

advance_t
Finder::getExpandDirection()
{
    return getExpandDirection(m_m->getFrameCount() - 1,
                              m_m->getOtherFrameCount() - 1);
}

advance_t
Finder::getExpandDirection(int row, int col)
{
    // To determine which direction to expand the search area in, we
    // look at the path costs along the leading edges of the search
    // area (the final row and column within the area). We find the
    // lowest path cost within the final row, and the lowest within
    // the final column, and we compare them. If the row is cheaper
    // then we expand by adding another row next to it; if the column
    // is cheaper then we expand by adding another column next to
    // it. (The overall lowest path cost across the row and column
    // represents the best alignment we have within the entire search
    // area given the data available and the assumption that the piece
    // is not ending yet.)

    int bestRow = row;
    int bestCol = col;
    normpathcost_t bestCost = INVALID_PATHCOST;

//    cerr << "Finder " << this << "::getExpandDirection: ";
    
    getBestEdgeCost(row, col, bestRow, bestCol, bestCost);

//    cerr << "at [" << row << "," << col << "] (cost " << m_m->getPathCost(row, col) << ") blocksize = " << m_m->getBlockSize() << " best is [" << bestRow << "," << bestCol << "] (cost " << bestCost << ")" << endl;
    
    if (bestRow == row) {
        if (bestCol == col) {
            return AdvanceBoth;
        } else {
            return AdvanceThis;
        }
    } else if (bestCol == col) {
        return AdvanceOther;
    } else {
        return AdvanceNone;
    }
}

void
Finder::recalculatePathCostMatrix(int r1, int c1, int r2, int c2) 
{
    int prevRowStart = 0, prevRowStop = 0;

    for (int r = r1; r <= r2; r++) {

        pair<int, int> colRange = m_m->getColRangeForRow(r);

        int rowStart = max(c1, colRange.first);
        int rowStop = min(c2 + 1, colRange.second);
        
        for (int c = rowStart; c < rowStop; c++) {

            advance_t dir = AdvanceNone;
            pathcost_t straightIncrement = m_m->getDistance(r, c);
            pathcost_t diagIncrement = pathcost_t(straightIncrement *
                                                  m_m->getDiagonalWeight());

            if (r > r1) {	// not first row
                pathcost_t min = INVALID_PATHCOST;
                if ((c > prevRowStart) && (c <= prevRowStop)) {
                    // diagonal from (r-1,c-1)
                    min = m_m->getPathCost(r-1, c-1) + diagIncrement;
                    dir = AdvanceBoth;
                }
                if ((c >= prevRowStart) && (c < prevRowStop)) {
                    // vertical from (r-1,c)
                    pathcost_t cost = m_m->getPathCost(r-1, c) + straightIncrement;
                    if ((min == INVALID_PATHCOST) || (cost < min)) {
                        min = cost;
                        dir = AdvanceThis;
                    }
                }
                if (c > rowStart) {
                    // horizontal from (r,c-1)
                    pathcost_t cost = m_m->getPathCost(r, c-1) + straightIncrement;
                    if ((min == INVALID_PATHCOST) || (cost < min)) {
                        min = cost;
                        dir = AdvanceOther;
                    }
                }
                
                m_m->setPathCost(r, c, dir, min);

            } else if (c > rowStart) {	// first row
                // horizontal from (r,c-1)
                m_m->setPathCost(r, c, AdvanceOther,
                                 m_m->getPathCost(r, c-1) + straightIncrement);
            }
        }

        prevRowStart = rowStart;
        prevRowStop = rowStop;
    }
} 

#ifdef PERFORM_ERROR_CHECKS
Finder::ErrorPosition
Finder::checkPathCostMatrix() 
{
    ErrorPosition err;

    int r1 = 0;
    int c1 = 0;
    int r2 = m_m->getFrameCount() - 1;
    int c2 = m_m->getOtherFrameCount() - 1;

    if (r2 < r1 || c2 < c1) {
        return err;
    }

    int prevRowStart = 0, prevRowStop = 0;

    for (int r = r1; r <= r2; r++) {

        pair<int, int> colRange = m_m->getColRangeForRow(r);

        int rowStart = max(c1, colRange.first);
        int rowStop = min(c2 + 1, colRange.second);
        
        for (int c = rowStart; c < rowStop; c++) {

            advance_t dir = AdvanceNone;
            pathcost_t updateTo = INVALID_PATHCOST;
            distance_t distance = m_m->getDistance(r, c);
            pathcost_t straightIncrement = distance;
            pathcost_t diagIncrement = pathcost_t(distance * m_m->getDiagonalWeight());
            err.distance = distance;

            if (r > r1) { // not first row
                pathcost_t min = INVALID_PATHCOST;
                if ((c > prevRowStart) && (c <= prevRowStop)) {
                    // diagonal from (r-1,c-1)
                    min = m_m->getPathCost(r-1, c-1) + diagIncrement;
                    err.prevCost = m_m->getPathCost(r-1, c-1);
                    dir = AdvanceBoth;
                }
                if ((c >= prevRowStart) && (c < prevRowStop)) {
                    // vertical from (r-1,c)
                    pathcost_t cost = m_m->getPathCost(r-1, c) + straightIncrement;
                    if ((min == INVALID_PATHCOST) || (cost < min)) {
                        min = cost;
                        err.prevCost = m_m->getPathCost(r-1, c);
                        dir = AdvanceThis;
                    }
                }
                if (c > rowStart) {
                    // horizontal from (r,c-1)
                    pathcost_t cost = m_m->getPathCost(r, c-1) + straightIncrement;
                    if ((min == INVALID_PATHCOST) || (cost < min)) {
                        min = cost;
                        err.prevCost = m_m->getPathCost(r, c-1);
                        dir = AdvanceOther;
                    }
                }

                updateTo = min;

            } else { // first row

                if (c > rowStart) {
                    // horizontal from (r,c-1)
                    updateTo = m_m->getPathCost(r, c-1) + straightIncrement;
                    err.prevCost = m_m->getPathCost(r, c-1);
                    dir = AdvanceOther;
                }
            }

            if (dir != AdvanceNone) {
                if (m_m->getAdvance(r, c) != dir) {
                    err.type = ErrorPosition::WrongAdvance;
                    err.r = r;
                    err.c = c;
                    err.costWas = m_m->getPathCost(r, c);
                    err.costShouldBe = updateTo;
                    err.advanceWas = m_m->getAdvance(r, c);
                    err.advanceShouldBe = dir;
                    return err;
                }
                if (m_m->getPathCost(r, c) != updateTo) {
                    err.type = ErrorPosition::WrongCost;
                    err.r = r;
                    err.c = c;
                    err.costWas = m_m->getPathCost(r, c);
                    err.costShouldBe = updateTo;
                    err.advanceWas = m_m->getAdvance(r, c);
                    err.advanceShouldBe = dir;
                    return err;
                }
            } else {
                // AdvanceNone should occur only at r = r1, c = c1
                if (r != r1 || c != c1) {
                    err.type = ErrorPosition::NoAdvance;
                    err.r = r;
                    err.c = c;
                    err.costWas = m_m->getPathCost(r, c);
                    err.costShouldBe = updateTo;
                    err.advanceWas = m_m->getAdvance(r, c);
                    err.advanceShouldBe = dir;
                    return err;
                }
            }
        }

        prevRowStart = rowStart;
        prevRowStop = rowStop;
    }

    return err;
}

void
Finder::checkAndReport()
{
    cerr << "Finder: Checking path-cost matrix..." << endl;
    ErrorPosition err = checkPathCostMatrix();
    if (err.type == ErrorPosition::NoError) {
        cerr << "No errors found" << endl;
    } else {
        cerr << "\nWARNING: Checking path-cost matrix returned mismatch:" << endl;
        cerr << "Type: " << err.type << ": ";
        switch (err.type) {
        case ErrorPosition::NoError: break;
        case ErrorPosition::WrongCost: cerr << "WrongCost"; break;
        case ErrorPosition::WrongAdvance: cerr << "WrongAdvance"; break;
        case ErrorPosition::NoAdvance: cerr << "NoAdvance"; break;
        }
        cerr << endl;
        cerr << "At row " << err.r << ", column " << err.c
             << "\nShould be advancing "
             << Matcher::advanceToString(err.advanceShouldBe)
             << ", advance in matrix is "
             << Matcher::advanceToString(err.advanceWas)
             << "\nPrev cost " << err.prevCost
             << " plus distance " << distance_print_t(err.distance)
             << " [perhaps diagonalised] gives "
             << err.costShouldBe << ", matrix contains " << err.costWas
             << endl;
        cerr << "Note: diagonal weight = " << m_m->getDiagonalWeight() << endl;
        cerr << endl;

        int w(4);
        int ww(15);

        cerr << "Distance matrix leading up to this point:" << endl;
        cerr << setprecision(12) << setw(w) << "";
        for (int i = -4; i <= 0; ++i) {
            cerr << setw(ww) << i;
        }
        cerr << endl;
        for (int j = -4; j <= 0; ++j) {
            cerr << setw(w) << j;
            for (int i = -4; i <= 0; ++i) {
                cerr << setw(ww)
                     << distance_print_t(m_m->getDistance(err.r + j, err.c + i));
            }
            cerr << endl;
        }
        cerr << endl;

        cerr << "Cost matrix leading up to this point:" << endl;
        cerr << setw(w) << "";
        for (int i = -4; i <= 0; ++i) {
            cerr << setw(ww) << i;
        }
        cerr << endl;
        for (int j = -4; j <= 0; ++j) {
            cerr << setw(w) << j;
            for (int i = -4; i <= 0; ++i) {
                cerr << setw(ww) << m_m->getPathCost(err.r + j, err.c + i);
            }
            cerr << endl;
        }
        cerr << endl;
    }
}
#endif

pathcost_t
Finder::getOverallCost()
{
    int ex = m_m->getOtherFrameCount() - 1;
    int ey = m_m->getFrameCount() - 1;

    if (ex < 0 || ey < 0) {
        return 0;
    }

    return m_m->getPathCost(ey, ex);
}

int
Finder::retrievePath(bool smooth, vector<int> &pathx, vector<int> &pathy)
{
    pathx.clear();
    pathy.clear();

#ifdef PERFORM_ERROR_CHECKS
    checkAndReport();
#endif

    int ex = m_m->getOtherFrameCount() - 1;
    int ey = m_m->getFrameCount() - 1;

    if (ex < 0 || ey < 0) {
        return 0;
    }
    
    int x = ex;
    int y = ey;

#ifdef DEBUG_FINDER
    cerr << "*** retrievePath: smooth = " << smooth << endl;
    cerr << "*** retrievePath: before: x = " << x << ", y = " << y << endl;
#endif

    if (m_duration2 > 0 && m_duration2 < m_m->getOtherFrameCount()) {
        x = m_duration2 - 1;
    }
    if (m_duration1 > 0 && m_duration1 < m_m->getFrameCount()) {
        y = m_duration1 - 1;
    }

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

#ifdef DEBUG_FINDER
    cerr << "*** retrievePath: start: x = " << x << ", y = " << y << endl;
#endif
    
    while (m_m->isAvailable(y, x) && (x > 0 || y > 0)) {

//        cerr << "x = " << x << ", y = " << y;
        
        pathx.push_back(x);
        pathy.push_back(y);

        switch (m_m->getAdvance(y, x)) {
        case AdvanceThis:
//            cerr << ", going down (dist = " << getDistance() << ")" << endl;
            y--;
            break;
        case AdvanceOther:
//            cerr << ", going left (dist = " << getDistance() << ")" << endl;
            x--;
            break;
        case AdvanceBoth:
//            cerr << ", going diag (dist = " << getDistance() << ")" << endl;
            x--;
            y--;
            break;
        case AdvanceNone: // this would indicate a bug, but we wouldn't want to hang
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
        int smoothedLen = Path().smooth
            (pathx, pathy, static_cast<int>(pathx.size()));
        return smoothedLen;
    } else {
        return static_cast<int>(pathx.size());
    }
}

