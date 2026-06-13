
#ifndef _DR_STATE_H_
#define _DR_STATE_H_

#include "../define.h"
#include <set>
#include <math.h>
#include <Kokkos_Core.hpp>

class Node   {

public:
    realArr xcoor,
        ycoor, head, rptFlag, invertElev,
        initDepth, fullDepth, surDepth, pondedArea, degree, inlet,
        updated, crownElev, inflow, outflow, losses, oldVolume, newVolume,
        fullVolume, overflow, oldDepth, newDepth, oldLatFlow, newLatFlow,
        oldFlowInflow, oldNetInflow, apiExtInflow, SDinflow, SDoutflow,
        SDinflowObs, SDoutflowObs,
        converged, newSurfArea, oldSurfArea, sumdqdh, dYdT, draingeVolume,
        backflow, backflowRatio, outfallFixedStage;
    intArr typee, subIndex, toNode, connectedInlet, numInlet, inletIndex, sealed,
        outfallType, outfallStageSeries, outfallHasFlapGate;
};


class Outfall   {
    public:
    realArr xcoor, ycoor, hasFlapGate, fixedStage, tideCurve, stageSeries, routeTo, vRouted;
    intArr typee;
};

/** Storage-unit shape coefficients (FUNCTIONAL / CYLINDRICAL only). */
class Storage {
public:
    intArr shape;
    realArr a0;
    realArr a1;
    realArr a2;
};

class Link   {
    public:
    realArr rptFlag, qFull,
    offset1, offset2, yFull, wMax, ywMax, aFull,
    rFull, sFull, sMax, s, qc, q0, qLimit, cLossInlet, 
    cLossOutlet, cLossAvg, seepRate, hasFlapGate, oldFlow, 
    newFlow, oldDepth, newDepth, oldVolume, newVolume, surfArea1, 
    surfArea2,dqdh, setting, targetSetting, timeLastSet,
    froude, direction, bypassed, normalFlow, inletControl;
    intArr flowClass;
    intArr typee;
    intArr node1, node2;
    intArr subIndex;
};

class Conduit   {
    public:
    realArr length, modLength, roughness, slope, beta, qMax, 
    a1, a2, q1, q2, q1Old, q2Old, evapLossRate, seepLossRate, 
    capacityLimited, superCritical, hasLosses, fullState, roughFactor;
};

/** TYPE3 pump curve table stored in flat arrays for Kokkos CPU/GPU use. */
class PumpCurves {
public:
    int nCurves;
    int nPts;
    intArr startIdx;
    intArr nPtsPerCurve;
    realArr x;
    realArr y;
};

class Pump {
public:
    intArr pumpCurve;
    intArr type;
    realArr yOn;
    realArr yOff;
    realArr initSetting;
    realArr xMin;
    realArr xMax;
};

/** River stage time series for TIMESERIES outfalls (flat storage). */
class RiverStages {
public:
    int nSeries;
    int nPts;
    intArr startIdx;
    intArr nPtsPerSeries;
    realArr time;
    realArr stage;
};

/** Node-Link adjacency table for node-centric parallel computation.
 *  Allows efficient lookup of all links connected to a given node.
 *  
 *  Usage:
 *    for (int j = adjStart(nodeId); j < adjStart(nodeId+1); j++) {
 *        int linkId = adjLinks(j);
 *        int position = adjPosition(j);  // 0=upstream(node1), 1=downstream(node2)
 *    }
 */
class NodeLinkAdjacency {
public:
    int nNodes;
    int nEntries;  // = nLinks * 2 (each link connects 2 nodes)
    
    // adjStart[i] = starting index in adjLinks for node i
    // adjStart[nNodes] = total number of entries (for boundary)
    intArr adjStart;
    
    // adjLinks[j] = link index for the j-th adjacency entry
    intArr adjLinks;
    
    // adjPosition[j] = 0 if node is link's node1 (upstream), 1 if node2 (downstream)
    intArr adjPosition;
    
    bool isBuilt = false;
};

#endif
