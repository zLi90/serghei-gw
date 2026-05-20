
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
        converged, newSurfArea, oldSurfArea, sumdqdh, dYdT, draingeVolume,
        backflow, backflowRatio;
    intArr typee, subIndex, toNode, connectedInlet, numInlet, inletIndex, sealed;
};


class Outfall   {
    public:
    realArr xcoor, ycoor, hasFlapGate, fixedStage, tideCurve, stageSeries, routeTo, vRouted;
    intArr typee;
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

#endif
