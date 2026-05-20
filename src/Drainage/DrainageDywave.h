#ifndef DRAINAGEDYWAVE_H
#define DRAINAGEDYWAVE_H

#include <cmath>
#include <string.h>
#include <stdio.h>
#include <array>
#include "../define.h"
#include "enums.h"
#include "globals.h"
#include "funcs.h"
#include "xsect.h"


class DrainageDywave   {

    friend struct FindNodeDepthsFunctor;

public:
    double Vrouted = 0.0;
    double extotal = 0.0;
    double outfallInflow = 0.0;
    double outfallDischarge = 0.0;
    double routingStep = 0.0;
    double drainage_two2one = 0.0;
    double drainage_one2two = 0.0;

    /*--------------------------------------
    drainage module exchange flow computation
    ---------------------------------------*/
    // Inlet dimensions (0.5m x 0.2m)
    double inletLength = 0.5;
    double inletWidth = 0.2;
    // Junction (manhole) dimensions - circular with diameter 0.6m, treated as square for grid calculation
    double junctionSize = 0.6;
    
    // Discharge coefficients
    double Cw = 0.6;  // Weir discharge coefficient
    double Co = 0.6;  // Orifice discharge coefficient
    double w = 1.0;
    
    // Option to choose between multi-grid and single-grid calculation method
    // true: use multi-grid method (251017 version)
    // false: use single-grid method (250523 version)
    bool useMultiGrid = false;

    inline void DrainageExchange_solve(State &state, Domain &dom, Node &Tnode)
    {
        // Initialize drainage exchange accumulators
        double drainageTwo2One = 0.0;
        double drainageOne2Two = 0.0;
        
        // Capture class member variables as local variables to avoid conflicts with global variables
        // This ensures correct behavior on both CPU and GPU
        const double localInletLength = inletLength;
        const double localInletWidth = inletWidth;
        const double localJunctionSize = junctionSize;
        const double localCw = Cw;  // Use class member, not global variable
        const double localCo = Co;  // Use class member, not global variable
        const double localw = w;
        const bool localUseMultiGrid = useMultiGrid;
        
    //    printf("inletLength:%f,inletWidth:%f,junctionSize:%f,Cw:%f,Co:%f,w:%f\n",localInletLength,localInletWidth,localJunctionSize,localCw,localCo,localw);
        
            // Multi-grid method (251017 version)
            Kokkos::parallel_reduce(
                Nobjects[NODE], 
                KOKKOS_LAMBDA(int ii, real &localTwo2One, real &localOne2Two) {
                    // Early exit for non-relevant node types
                    // printf("Processing node %d\n", ii);
                    // printf("Tnode.typee(ii)=%d\n", Tnode.typee(ii));
                    if (Tnode.typee(ii) != INLET && 
                        !(Tnode.typee(ii) == JUNCTION && Tnode.sealed(ii) == FALSE && Tnode.connectedInlet(ii) == FALSE) &&
                        !(Tnode.typee(ii) == JUNCTION && Tnode.connectedInlet(ii) == TRUE && Tnode.sealed(ii) == FALSE) &&
                        !(Tnode.typee(ii) == JUNCTION && Tnode.sealed(ii) == TRUE && Tnode.connectedInlet(ii) == FALSE) &&
                        !(Tnode.typee(ii) == JUNCTION && Tnode.sealed(ii) == TRUE && Tnode.connectedInlet(ii) == TRUE)) {
                        return;
                    }
                    const int j = floor(dom.ny - ((Tnode.ycoor(ii) / dom.dxConst) + 1));
                    const int i = floor((Tnode.xcoor(ii) / dom.dxConst));
                    const int iGlobNode = packIndicesUniformGrid((dom.ny + 2 * hc), (dom.nx + 2 * hc), j, i);
                    
                    // Update node head
                    Tnode.head(ii) = Tnode.invertElev(ii) + Tnode.newDepth(ii);
                    const double nodeHead = Tnode.head(ii);
                    
                    // Process different node types directly in lambda
                    const int nodeType = Tnode.typee(ii);
                    const bool isInlet = (nodeType == INLET);
                    const bool isJunctionUnsealedNoInlet = (nodeType == JUNCTION && Tnode.sealed(ii) == FALSE && Tnode.connectedInlet(ii) == FALSE);
                    const bool isJunctionUnsealedWithInlet = (nodeType == JUNCTION && Tnode.connectedInlet(ii) == TRUE && Tnode.sealed(ii) == FALSE);
                    const bool isJunctionSealedNoInlet = (nodeType == JUNCTION && Tnode.sealed(ii) == TRUE && Tnode.connectedInlet(ii) == FALSE);
                    const bool isJunctionSealedWithInlet = (nodeType == JUNCTION && Tnode.sealed(ii) == TRUE && Tnode.connectedInlet(ii) == TRUE);
                    
                    // Determine node coverage dimensions
                    double nodeLength, nodeWidth;
                    
                    if (isInlet) {
                        nodeLength = localInletLength;
                        nodeWidth = localInletWidth;
                    } else if (isJunctionUnsealedNoInlet || isJunctionUnsealedWithInlet || 
                              isJunctionSealedNoInlet || isJunctionSealedWithInlet) {
                        // For junction, use square coverage (circular manhole approximated as square)
                        nodeLength = localJunctionSize;
                        nodeWidth = localJunctionSize;
                    }
                    
                    // For nodes that need multi-grid treatment
                if (localUseMultiGrid) {
                    // Calculate node coverage area
                    const double halfLength = nodeLength / 2.0;
                    const double halfWidth = nodeWidth / 2.0;
                    
                    // Calculate grid indices for node boundaries
                    const double xMin = Tnode.xcoor(ii) - halfLength;
                    const double xMax = Tnode.xcoor(ii) + halfLength;
                    const double yMin = Tnode.ycoor(ii) - halfWidth;
                    const double yMax = Tnode.ycoor(ii) + halfWidth;
                    
                    const int iMin = floor(xMin / dom.dxConst);
                    const int iMax = floor(xMax / dom.dxConst);
                    const int jMin = floor(dom.ny - ((yMax / dom.dxConst) + 1));
                    const int jMax = floor(dom.ny - ((yMin / dom.dxConst) + 1));
                    
                    // Ensure indices are within bounds
                    const int iStart = max(0, iMin);
                    const int iEnd = min((int)(dom.nx + 2 * hc - 1), iMax);
                    const int jStart = max(0, jMin);
                    const int jEnd = min((int)(dom.ny + 2 * hc - 1), jMax);
                    
                    // Calculate average state.h and total area from covered cells
                    double avgH = 0.0;
                    double avgZ = 0.0;
                    int numCells = 0;
                    
                    for (int jj = jStart; jj <= jEnd; jj++) {
                        for (int iii = iStart; iii <= iEnd; iii++) {
                            const int iGlob = packIndicesUniformGrid((dom.ny + 2 * hc), (dom.nx + 2 * hc), jj, iii);
                            
                            // Check if this cell is within node bounds
                            const double cellX = (iii + 0.5) * dom.dxConst;
                            const double cellY = (dom.ny - (jj + 0.5)) * dom.dxConst;
                            
                            if (cellX >= xMin && cellX <= xMax && cellY >= yMin && cellY <= yMax) {
                                avgH += state.h(iGlob);
                                avgZ += state.z(iGlob);
                                numCells++;
                            }
                        }
                    }
                    
                    if (numCells > 0) {
                        avgH /= numCells;
                        avgZ /= numCells;
                        const double totalArea = numCells * dom.cellArea();
                        const double surfaceLevel = avgH + avgZ;
                        
                        // Case 1: Surface water flows to drainage system (drainage)
                        // For inlet: always check; for junction: only if unsealed or connected to inlet
                        bool shouldDrain = false;
                        if (isInlet) {
                            shouldDrain = true;
                        } else if (isJunctionUnsealedNoInlet || isJunctionUnsealedWithInlet) {
                            shouldDrain = true;
                        }
                        
                        if (shouldDrain && nodeHead < surfaceLevel && avgH > TOL1) {
                            double dischargeRate;
                            if (nodeHead < avgZ) {
                            // if (nodeHead < state.z(iGlobNode)) {
                                // Free surface flow
                                if (isInlet) {
                                    dischargeRate = localCw * 2*(localInletLength+localInletWidth) * avgH * sqrt(2 * GRAV * avgH);
                                    // dischargeRate = localCw * w * state.h(iGlobNode) * sqrt(2 * GRAV * state.h(iGlobNode));
                                } else {
                                    // For junction, use appropriate width
                                    dischargeRate = localCw * nodeWidth*4 * avgH * sqrt(2 * GRAV * avgH);
                                    }
                            } else {
                                // Orifice flow
                                dischargeRate = localCo * Amh * sqrt(2 * GRAV * (surfaceLevel - nodeHead));
                            }
                            
                            const double volume = min(dischargeRate * dom.dt, avgH * totalArea);
                            
                            // Distribute volume extraction equally across covered cells
                            const double volumePerCell = volume / numCells;
                            const double hChangePerCell = volumePerCell / dom.cellArea();
                            
                            for (int jj = jStart; jj <= jEnd; jj++) {
                                for (int iii = iStart; iii <= iEnd; iii++) {
                                    const int iGlob = packIndicesUniformGrid((dom.ny + 2 * hc), (dom.nx + 2 * hc), jj, iii);
                                    // const double cellX = (iii + 0.5) * dom.dxConst;
                                    // const double cellY = (dom.ny - (jj + 0.5)) * dom.dxConst;
                                    
                                    // if (cellX >= xMin && cellX <= xMax && cellY >= yMin && cellY <= yMax) {
                                        state.h(iGlob) -= hChangePerCell;
                                        state.Dqss(iGlob) = dischargeRate / numCells;
                                        state.Vqss(iGlob) = volumePerCell;
                                        
                                        // Ensure minimum water depth
                                        if (state.h(iGlob) < TOL12) {
                                            state.h(iGlob) = 0.0;
                                        }
                                    // }
                                }
                            }
                            
                            // Update node inflow
                            Tnode.SDinflow(ii) = volume / dom.dt;
                            localTwo2One += volume;
                        }
                        
                        // Case 2: Drainage system overflows to surface (overflow)
                        // For inlet: always check; for junction: only if unsealed or sealed without inlet
                        bool shouldOverflow = false;
                        if (isInlet) {
                            shouldOverflow = true;
                        } else if (isJunctionUnsealedNoInlet || isJunctionSealedNoInlet) {
                            shouldOverflow = true;
                        }
                        
                        if (shouldOverflow && nodeHead > surfaceLevel) {
                            double dischargeRate;
                            if (Tnode.invertElev(ii) > surfaceLevel) {
                                // Overflow from manhole
                                dischargeRate = localCo * Amh * sqrt(2 * GRAV * max((Tnode.newDepth(ii) - Tnode.fullDepth(ii)), 0.0));
                            } else {
                                // Orifice flow
                                dischargeRate = localCo * Amh * sqrt(2 * GRAV * (nodeHead - surfaceLevel));
                            }
                            
                            // Calculate volume and distribute equally across covered cells
                            const double volume = dischargeRate * dom.dt;
                            const double volumePerCell = volume / numCells;
                            const double hChangePerCell = volumePerCell / dom.cellArea();
                            
                            for (int jj = jStart; jj <= jEnd; jj++) {
                                for (int iii = iStart; iii <= iEnd; iii++) {
                                    const int iGlob = packIndicesUniformGrid((dom.ny + 2 * hc), (dom.nx + 2 * hc), jj, iii);
                                    const double cellX = (iii + 0.5) * dom.dxConst;
                                    const double cellY = (dom.ny - (jj + 0.5)) * dom.dxConst;
                                    
                                    if (cellX >= xMin && cellX <= xMax && cellY >= yMin && cellY <= yMax) {
                                        state.h(iGlob) += hChangePerCell;
                                        state.Dqss(iGlob) = dischargeRate / numCells;
                                        state.Vqss(iGlob) = volumePerCell;
                                    }
                                }
                            }
                            
                            // Update node outflow
                            Tnode.SDoutflow(ii) = dischargeRate;
                            localOne2Two += volume;
                        }
                    }
                }

                else{

                    bool shouldDrain = false;
                    if (isInlet) {
                        shouldDrain = true;
                    } else if (isJunctionUnsealedNoInlet || isJunctionUnsealedWithInlet) {
                        shouldDrain = true;
                    }

                    if (shouldDrain && Tnode.head(ii) < state.h(iGlobNode) + state.z(iGlobNode) && state.h(iGlobNode) > TOL1)
                            {
                                if (Tnode.head(ii) < state.z(iGlobNode))
                                {
                                    state.Dqss(iGlobNode) = localCw * localw * state.h(iGlobNode) * sqrt(2 * GRAV * state.h(iGlobNode)); 
                                }
                                else
                                {
                                    state.Dqss(iGlobNode) = localCo * Amh * sqrt(2 * GRAV * (state.h(iGlobNode) + state.z(iGlobNode) - Tnode.head(ii)));
                                    
                                }
                                state.Vqss(iGlobNode) = min(state.Dqss(iGlobNode) * dom.dt, state.h(iGlobNode) * dom.cellArea());
                                state.h(iGlobNode) -= state.Vqss(iGlobNode) / dom.cellArea();
                                Tnode.SDinflow(ii) = state.Vqss(iGlobNode) / dom.dt;
                                localTwo2One += state.Vqss(iGlobNode);
                            }
                    bool shouldOverflow = false;
                    if (isInlet) {
                        shouldOverflow = true;
                    } else if (isJunctionUnsealedNoInlet || isJunctionSealedNoInlet) {
                        shouldOverflow = true;
                    }
                    
                    if (shouldOverflow && Tnode.head(ii) > state.h(iGlobNode) + state.z(iGlobNode))
                            {
                                // printf("Tnode.invertElev(ii) %f Tnode.newDepth %f state.h(iGlob) %f state.z(iGlob) %f \n", Tnode.invertElev(ii), Tnode.newDepth(ii), state.h(iGlob), state.z(iGlob));
                                if (Tnode.invertElev(ii) > state.h(iGlobNode) + state.z(iGlobNode))
                                {
                                    state.Dqss(iGlobNode) = localCo * Amh * sqrt(2 * GRAV * max((Tnode.newDepth(ii) - Tnode.fullDepth(ii)),0));    
                                }
                                else
                                {
                                    state.Dqss(iGlobNode) = localCo * Amh * sqrt(2 * GRAV * (Tnode.head(ii) - (state.h(iGlobNode) + state.z(iGlobNode))));
                                }
                                state.Vqss(iGlobNode) = state.Dqss(iGlobNode) * dom.dt;
                                state.h(iGlobNode) += state.Vqss(iGlobNode) / dom.cellArea();
                                Tnode.SDoutflow(ii) = state.Dqss(iGlobNode);
                                localOne2Two += state.Vqss(iGlobNode);
                            }
                        
                            if (state.h(iGlobNode) < TOL12) {
                                state.h(iGlobNode) = 0.0;
                            }
                    }
            },
                Kokkos::Sum<real>(drainageTwo2One), Kokkos::Sum<real>(drainageOne2Two));
            
                // Update global drainage totals
                drainage_two2one += drainageTwo2One;
                drainage_one2two += drainageOne2Two;
                
                Kokkos::fence();
    }


public:
/*--------------------------------------
drainage module pipe flow comuputation
---------------------------------------*/
    inline void routing_execute(double tStep, Node& Tnode, Link& Tlink, Conduit& Tconduit, Outfall& Toutfall,
                                SergheiTimers& timers)
    {   

        initSystemInflows(Tnode);
        inletBackflow(Tnode);
        addSystemInflows(Tnode);
        routeFlow(tStep, Tnode, Tlink, Tconduit, Toutfall, timers);
        routingStep = getVariableStep(tStep, Tnode, Tlink, Tconduit);
    }

    // inline void initSystemInflows(Node& Tnode)
    // {
    //     int j;
    //     for (j = 0; j < Nobjects[NODE]; j++)
    //     {
    //         Tnode.oldLatFlow(j)  = Tnode.newLatFlow(j);
    //         Tnode.newLatFlow(j)  = 0.0;
    //     }
    // }
    inline void initSystemInflows(Node& Tnode)
    {
        const int nNode = Nobjects[NODE];
        Kokkos::parallel_for(
            nNode,
            KOKKOS_LAMBDA(int j) {
                Tnode.oldLatFlow(j) = Tnode.newLatFlow(j);
                Tnode.newLatFlow(j) = 0.0;
            });
        Kokkos::fence();  // 若后续立刻依赖更新后的 newLatFlow，建议保留；与项目其它 drainage 处一致
    }
    
    inline void addSystemInflows(Node& Tnode)
    {
        int j;
        for (j = 0; j < Nobjects[NODE]; j++)
        {
            if (Tnode.typee(j) == INLET)//zn251020
            {
                Tnode.newLatFlow(j) += Tnode.backflow(j);//检查井回流量
                // Tnode.newLatFlow(j) += Tnode.SDinflow(j);//地表来水
                Tnode.newLatFlow(j) -= Tnode.SDoutflow(j);//溢流量
            }
            else if (Tnode.typee(j) == JUNCTION)
            {
                //对于连接了inlet且有盖子的manhole，要考虑inlet的入流和回流。
                //当该manhole的overflow<=0时，考虑inlet的入流量；
                if (Tnode.connectedInlet(j) == TRUE && Tnode.sealed(j) == TRUE)
                {
                    if (Tnode.overflow(j) <= 0)
                    {
                        Tnode.newLatFlow(j) += Tnode.SDinflow(Tnode.inletIndex(j));//检查井回流量
                    }
                    //当该manhole的overflow>0时，此时不考虑inlet的入流量；
                    else if (Tnode.overflow(j) > 0)
                    {
                        Tnode.newLatFlow(j) -= Tnode.backflow(Tnode.inletIndex(j));
                    }
                }
                //对于连接了inlet且没有盖子的manhole，要考虑inlet的入流和回流，以及地表来水。
                else if (Tnode.connectedInlet(j) == TRUE && Tnode.sealed(j) == FALSE)
                {
                    if (Tnode.overflow(j) <= 0)
                    {
                        Tnode.newLatFlow(j) += Tnode.SDinflow(Tnode.inletIndex(j));//检查井回流量
                    }
                    //当该manhole的overflow>0时，此时不考虑inlet的入流量；
                    else if (Tnode.overflow(j) > 0)
                    {
                        Tnode.newLatFlow(j) -= Tnode.backflow(Tnode.inletIndex(j));
                    }
                    Tnode.newLatFlow(j) += Tnode.SDinflow(j);
                }
                //对于没有连接inlet且没有盖子的manhole，要考虑地表来水和溢流。
                else if (Tnode.connectedInlet(j) == FALSE && Tnode.sealed(j) == FALSE)
                {
                    Tnode.newLatFlow(j) += Tnode.SDinflow(j);
                    Tnode.newLatFlow(j) -= Tnode.SDoutflow(j);
                }
                //对于没有连接inlet且有盖子的manhole，只考虑溢流。
                else{
                    Tnode.newLatFlow(j) -= Tnode.SDoutflow(j);
                    }
             
            }
            Tnode.SDinflow(j) = 0.0;
            Tnode.SDoutflow(j) = 0.0;
        }
    }

    inline void inletBackflow(Node& Tnode)
    {
        int j;
        for (j = 0; j < Nobjects[NODE]; j++)
        {
            if (Tnode.typee(j) == INLET)//zn251020
            {
                Tnode.backflowRatio(j) = 1/Tnode.numInlet(Tnode.toNode(j));
                Tnode.backflow(j) = Tnode.overflow(Tnode.toNode(j)) * Tnode.backflowRatio(j);//zn251020当前不考虑具体的backflowRatio
            }
        }
    }

inline double getVariableStep(double maxStep, Node& Tnode, Link& Tlink, Conduit& Tconduit)

{
    int    minLink = -1;                
    int    minNode = -1;                
    double tMin, tMinLink, tMinNode;                                   
    tMin = maxStep;

    tMinLink = getLinkStep(tMin, &minLink, Tlink, Tconduit);
    tMinNode = getNodeStep(tMinLink, &minNode,Tnode);

    tMin = tMinLink;
    if ( tMinNode < tMin )
    {
        tMin = tMinNode ;
        minLink = -1;
    }
    if ( tMin < MinRouteStep ) tMin = MinRouteStep;
    return tMin;
}

inline double getLinkStep(double tMin, int *minLink, Link& Tlink, Conduit& Tconduit)
{
    int    i;                          
    int    k;                          
    double q;                         
    double t;                           
    double tLink = tMin;              

    for ( i = 0; i < Nobjects[LINK]; i++ )
    {
        if ( Tlink.typee(i) == CONDUIT )
        {
            k = Tlink.subIndex(i);
            q = fabs(Tlink.newFlow(i));
            if (q <= FUDGE || Tconduit.a1(k) <= FUDGE || Tlink.froude(i) <= 0.01)
                continue;
            t = Tlink.newVolume(i)/ q;
            t = t * Tconduit.modLength(k) / link_getLength(i, Tlink, Tconduit);
            t = t * Tlink.froude(i) / (1.0 + Tlink.froude(i)) * CourantFactor;
            if ( t < tLink )
            {
                tLink = t;
                *minLink = i;
            }
        }
    }
    return tLink;
}

inline double getNodeStep(double tMin, int *minNode, Node& Tnode)
{
    int    i;                           
    double maxDepth;                    
    double dYdT;                       
    double t1;                         
    double tNode = tMin;               

    for ( i = 0; i < Nobjects[NODE]; i++ )
    {
        if ( Tnode.typee(i) == OUTFALL ) continue;
        if ( Tnode.newDepth(i) <= FUDGE) continue;
        if ( Tnode.newDepth(i)  + FUDGE >=
             Tnode.crownElev(i) - Tnode.invertElev(i) ) continue;
        maxDepth = (Tnode.crownElev(i) - Tnode.invertElev(i)) * 0.25;
        if ( maxDepth < FUDGE ) continue;
        dYdT = Tnode.dYdT(i);
        if (dYdT < FUDGE ) continue;

        t1 = maxDepth / dYdT;
        if ( t1 < tNode )
        {
            tNode = t1;
            *minNode = i;
        }
    }
    return tNode;
}


    inline void routeFlow(double routingStep, Node& Tnode, Link& Tlink, Conduit& Tconduit, Outfall& Toutfall,
                          SergheiTimers& timers)
    {
        int j;
        for (j = 0; j < Nobjects[LINK]; j++) link_setOldHydState(j, Tlink, Tconduit, Tnode);
        for (j = 0; j < Nobjects[NODE]; j++) node_setOldHydState(j, Tnode);
        for (j = 0; j < Nobjects[NODE]; j++) node_initFlows(j, routingStep, Tnode);
        if ( Nobjects[LINK] > 0 )
        {
            dynwave_execute(routingStep, Tnode, Tlink, Tconduit, Toutfall, timers);
        }
    }

    inline void link_setOldHydState(int j, Link& Tlink, Conduit& Tconduit, Node& Tnode)

    {
        int k, n1, n2;
        double z;

        Tlink.oldDepth(j)  = Tlink.newDepth(j);
        Tlink.oldFlow(j)   = Tlink.newFlow(j);
        Tlink.oldVolume(j) = Tlink.newVolume(j);
     
        n1 = Tlink.node1(j);
        z = Tnode.invertElev(n1) + Tlink.offset1(j) + Tlink.yFull(j);
        Tnode.crownElev(n1) = max(Tnode.crownElev(n1), z);
        n2 = Tlink.node2(j);
        z = Tnode.invertElev(n2) + Tlink.offset2(j) + Tlink.yFull(j);
        Tnode.crownElev(n2) = max(Tnode.crownElev(n2), z);
        Tlink.flowClass(j) = DRY;
        Tlink.dqdh(j) = 0.0;
    }

    inline void node_setOldHydState(int j, Node& Tnode)
 
    {
        Tnode.oldDepth(j)    = Tnode.newDepth(j);
        Tnode.oldVolume(j)   = Tnode.newVolume(j);
        Tnode.oldFlowInflow(j) = Tnode.inflow(j);
        Tnode.oldNetInflow(j) = Tnode.inflow(j) - Tnode.outflow(j);
    }

    void node_initFlows(int j, double tStep, Node& Tnode)

    {
        if (Tnode.newLatFlow(j) > 0.0)
        {
            Tnode.inflow(j) = Tnode.newLatFlow(j);
        }else{
            Tnode.outflow(j) = Tnode.newLatFlow(j);
        }

        if ( Tnode.newVolume(j) > Tnode.fullVolume(j) )
            Tnode.overflow(j) = (Tnode.newVolume(j) - Tnode.fullVolume(j)) / tStep;
        else Tnode.overflow(j) = 0.0;
    }


    inline void dynwave_execute(double tStep, Node& Tnode, Link& Tlink, Conduit& Tconduit, Outfall& Toutfall,
                                SergheiTimers& timers)
    {
        int converged;
        Steps = 0;
        converged = FALSE;
        initRoutingStep(Tnode, Tlink, Tconduit);
        while ( Steps < MaxTrials )
        {
            initNodeStates(Tnode);
            Kokkos::Timer timerLinkFlows;
            findLinkFlows(Steps,tStep, Tnode, Tlink, Tconduit, Toutfall);
            timers.drainageLinkFlows += timerLinkFlows.seconds();
            timerLinkFlows.reset();
            converged = findNodeDepths(tStep, Tnode, Tlink, Tconduit);
            timers.drainageNodeDepths += timerLinkFlows.seconds();
            Steps++;
            if ( Steps > 1 )
            {
                if ( converged ) break;
                findBypassedLinks(Tnode, Tlink);
            }
        }
        for (int ii = 0; ii<Nobjects[NODE];ii++)
            {
                double v1,v2;
                if (Tnode.typee(ii) == OUTFALL)
                {
                    // printf("the outfall index is :%d\n",ii);
                    // printf("the inflow is: %f\n",Tnode.inflow(ii) );
                    // printf("t is : %f\n",tStep);
                    
                    v1 = Tnode.inflow(ii) * tStep;
                    v2 = Tnode.outflow(ii) * tStep;
                    outfallInflow += v1;
                    outfallDischarge = Tnode.inflow(ii);
                } 
                }
        findLimitedLinks(Tnode,Tlink,Tconduit);
    }

    inline void   findBypassedLinks(Node& Tnode, Link& Tlink)
    {
        int i;
        for (i = 0; i < Nobjects[LINK]; i++)
        {
            if ( Tnode.converged(Tlink.node1(i)) &&
                Tnode.converged(Tlink.node2(i)) )
                Tlink.bypassed(i) = TRUE;
            else Tlink.bypassed(i) = FALSE;
        }
    }


    inline void  findLimitedLinks(Node& Tnode, Link& Tlink, Conduit& Tconduit)

    {
        int    j, n1, n2, k;
        double h1, h2;

        for (j = 0; j < Nobjects[LINK]; j++)
        {
            k = Tlink.subIndex(j);
            Tconduit.capacityLimited(k) = FALSE;
            if ( Tconduit.a1(k) >= Tlink.aFull(j) )
            {

                n1 = Tlink.node1(j);
                n2 = Tlink.node2(j);
                h1 = Tnode.newDepth(n1) + Tnode.invertElev(n1);
                h2 = Tnode.newDepth(n2) + Tnode.invertElev(n2);
                if ( (h1 - h2) > fabs(Tconduit.slope(k)) * Tconduit.length(k) )
                    Tconduit.capacityLimited(k) = TRUE;
            }
        }
    }

     inline void updateNodeDepths(double dt, Node& Tnode)

    {
        for (int ii = 0; ii < Nobjects[NODE]; ii++){
            if ( Tnode.SDinflow(ii) > 0.0 )
            {    
                Tnode.inflow(ii) = Tnode.SDinflow(ii);
            }
            if (Tnode.SDoutflow(ii) > 0.0)
            {    
                Tnode.outflow(ii) = Tnode.SDoutflow(ii);
            }
        if (Tnode.typee(ii) !=OUTFALL ) setNodeDepth(ii, dt, Tnode, Steps);
        };      
    }

    inline void link_setOutfallDepth(int j, Node& Tnode, Link& Tlink, Conduit& Tconduit)
    {
        int     n;                    
        double  z;                
        double  q;                      
        double  yCrit = 0.0;               
        double  yNorm = 0.0;   

        if ( Tnode.typee(Tlink.node2(j)) == OUTFALL )
        {
            n = Tlink.node2(j);
            z = Tlink.offset2(j);
        }
        else if ( Tnode.typee(Tlink.node1(j)) == OUTFALL )
        {
            n = Tlink.node1(j);
            z = Tlink.offset1(j);
        }
        else return;

        if ( Tlink.typee(j) == CONDUIT )
        {

            q = fabs(Tlink.newFlow(j));
            yNorm = link_getYnorm(j, q, Tlink, Tconduit);
            yCrit = xsect_getYcrit(j, q, Tlink);
        }

        outfall_setOutletDepth(n, yNorm, yCrit, z, Tnode);
    }
    void outfall_setOutletDepth(int j, double yNorm, double yCrit, double z, Node& Tnode)

{
    if ( z > 0.0 ) Tnode.newDepth(j) = 0.0;
    else Tnode.newDepth(j) = min(yNorm, yCrit);
}


    static KOKKOS_INLINE_FUNCTION void setNodeDepth(int i, double dt, Node& Tnode, int steps)
    {
        int     canPond;
        int     isPonded;
        int     isSurcharged = FALSE;
        double  dQ;
        double  dV;
        double  dy;
        double  yMax;
        double  yOld;
        double  yLast;
        double  yNew;
        double  yCrown;
        double  surfArea;
        double  denom;
        double  corr;
        double  f;

        canPond = TRUE;
        isPonded = (canPond && Tnode.newDepth(i) > Tnode.fullDepth(i));
        yCrown = Tnode.crownElev(i) - Tnode.invertElev(i);
        yOld = Tnode.oldDepth(i);
        yLast = Tnode.newDepth(i);
        Tnode.overflow(i) = 0.0;
        surfArea = Tnode.newSurfArea(i);
        surfArea = max(surfArea, MINSURFAREA);
        dQ = Tnode.inflow(i) - Tnode.outflow(i);
        dV = 0.5 * (Tnode.oldNetInflow(i) + dQ) * dt;
        if (isPonded) isSurcharged = FALSE;
        else isSurcharged = (yCrown > 0.0 && yLast > yCrown);
        if (!isSurcharged)
        {
            dy = dV / surfArea;
            yNew = yOld + dy;
            if ( !isPonded ) Tnode.oldSurfArea(i) = surfArea;
            if ( steps > 0 )
            {
                yNew = (1.0 - OMEGA) * yLast + OMEGA * yNew;
            }
            if ( isPonded && yNew < Tnode.fullDepth(i) )
                yNew = Tnode.fullDepth(i) - FUDGE;

        }
        else
        {
            corr = 1.0;
            denom = Tnode.sumdqdh(i);
            if ( yLast < 1.25 * yCrown )
            {
                f = (yLast - yCrown) / yCrown;
                denom += (Tnode.oldSurfArea(i)/dt -
                        Tnode.sumdqdh(i)) * exp(-15.0 * f);
            }

            if ( denom == 0.0 ) dy = 0.0;
            else dy = corr * dQ / denom;
            yNew = yLast + dy;
            if ( yNew < yCrown ) yNew = yCrown - FUDGE;
            if ( canPond && yNew > Tnode.fullDepth(i))
                yNew = Tnode.fullDepth(i) + FUDGE;
        }
        if ( yNew < 0 ) yNew = 0.0;
        yMax = Tnode.fullDepth(i);
        if ( yNew > yMax )
        {
            yNew = DrainageDywave::getFloodedDepth(i, canPond, dV, yNew, yMax, dt, Tnode);
        }
        else Tnode.newVolume(i) = DrainageDywave::node_getVolume(i, yNew, Tnode);
        Tnode.dYdT(i) = fabs(yNew - yOld) / dt;
        Tnode.newDepth(i) = yNew;
    }

    static KOKKOS_INLINE_FUNCTION double getFloodedDepth(int i, int canPond, double dV, double yNew,
                        double yMax, double dt, Node& Tnode)
    {
            Tnode.newVolume(i) = max((Tnode.oldVolume(i)+dV), Tnode.fullVolume(i));
            Tnode.overflow(i) = (Tnode.newVolume(i) -
                max(Tnode.oldVolume(i), Tnode.fullVolume(i))) / dt;

        if ( Tnode.overflow(i) < FUDGE ) Tnode.overflow(i) = 0.0;
        return yNew;
    }

    static KOKKOS_INLINE_FUNCTION double node_getVolume(int j, double d, Node const& Tnode)
    {
        if ( Tnode.fullDepth(j) > 0.0 )
            return Tnode.fullVolume(j) * (d / Tnode.fullDepth(j));
        else return 0.0;
    }

    /** View handles copied by value; functor is trivially copyable to device like Kokkos patterns. */
    struct FindNodeDepthsFunctor {
        mutable Node node;
        double dt;
        int stepsSnap;
        KOKKOS_INLINE_FUNCTION void operator()(int ii) const
        {
            if (node.typee(ii) == OUTFALL) {
                return;
            }
            const double yOld = node.newDepth(ii);
            DrainageDywave::setNodeDepth(ii, dt, node, stepsSnap);
            node.converged(ii) = TRUE;
            if (fabs(yOld - node.newDepth(ii)) > 0) {
                node.converged(ii) = FALSE;
            }
        }
    };

    struct FindNodeDepthsConvergedFunctor {
        Node node;
        KOKKOS_INLINE_FUNCTION void operator()(int i, int& lsum) const
        {
            if (node.typee(i) == OUTFALL) {
                return;
            }
            if (node.converged(i) == FALSE) {
                lsum += 1;
            }
        }
    };

    inline int findNodeDepths(double dt, Node& Tnode, Link& Tlink, Conduit& Tconduit)
    {
        for ( int i = 0; i < Nobjects[LINK]; i++ ) link_setOutfallDepth(i, Tnode, Tlink, Tconduit);
        {
            const int nNode = Nobjects[NODE];
            const int stepsSnap = Steps;
            Kokkos::parallel_for(
                nNode,
                FindNodeDepthsFunctor{Tnode, dt, stepsSnap});
            Kokkos::fence();
        }
        
    {
            const int nNode = Nobjects[NODE];
            int nNotConverged = 0;
            Kokkos::parallel_reduce(
                "findNodeDepths_converged",
                nNode,
                FindNodeDepthsConvergedFunctor{Tnode},
                Kokkos::Sum<int>(nNotConverged));
            Kokkos::fence();
            return (nNotConverged == 0) ? TRUE : FALSE;
        }

        // for (int i = 0; i < Nobjects[NODE]; i++)
        // {
        //     if ( Tnode.typee(i) == OUTFALL ) continue;
        //     if (Tnode.converged(i) == FALSE) return FALSE;
        // }
        // return TRUE;
    }

    inline void  initRoutingStep(Node &Tnode, Link &Tlink, Conduit &Tconduit)
    {
        
        for (int i = 0; i < Nobjects[NODE]; i++)
        {
            Tnode.converged(i) = FALSE;
            Tnode.dYdT(i) = 0.0;
        }
        for (int i = 0; i < Nobjects[LINK]; i++)
        {
                Tlink.bypassed(i) = FALSE;
                Tlink.surfArea1(i) = 0.0;
                Tlink.surfArea2(i) = 0.0;
            }
        for (int i = 0; i < Nlinks[CONDUIT]; i++)
            {
                Tconduit.a2(i) = Tconduit.a1(i);
            }

    }

    inline void initNodeStates(Node& Tnode)
    {
        const int nNode = Nobjects[NODE];
        Kokkos::parallel_for(
            nNode,
            KOKKOS_LAMBDA(int i) {
                Tnode.newSurfArea(i) =
                    DrainageDywave::node_getPondedArea(i, Tnode.newDepth(i), Tnode);
                Tnode.inflow(i) = 0.0;
                Tnode.outflow(i) = 0.0;
                const double nlf = Tnode.newLatFlow(i);
                if (nlf >= 0.0) {
                    Tnode.inflow(i) += nlf;
                } else {
                    Tnode.outflow(i) -= nlf;
                }
                Tnode.sumdqdh(i) = 0.0;
            });
        Kokkos::fence();
    }

    /** Host + device: ponded area for node j (d unused; kept for call-site parity). */
    static KOKKOS_INLINE_FUNCTION double node_getPondedArea(int j, double /*d*/, Node const& Tnode)
    {
        double a = Tnode.pondedArea(j);
        return (a <= 0.0) ? 0.0 : a;
    }

KOKKOS_INLINE_FUNCTION double getThetaOfPsiC(double psi)
{
    int    k;
    double theta, theta1, ap, tt, tt23, t3, d;

    if      (psi > 0.90)  theta = 4.17 + 1.12 * (psi - 0.90) / 0.176;
    else if (psi > 0.5)   theta = 3.14 + 1.03 * (psi - 0.5) / 0.4;
    else if (psi > 0.015) theta = 1.2 + 1.94 * (psi - 0.015) / 0.485;
    else                  theta = 0.12103 - 55.5075 * psi +
                                  15.62254 * sqrt(psi);
    theta1 = theta;
    ap     = (2.0*PI) * psi;

    for (k = 1; k <= 40; k++)
    {
        theta    = fabs(theta);
        tt       = theta - sin(theta);
        tt23     = pow(tt, 2./3.);
        t3       = pow(theta, 1./3.);
        d        = ap * theta / t3 - tt * tt23;
        d        = d / ( ap*(2./3.)/t3 - (5./3.)*tt23*(1.0-cos(theta)) );
        theta    = theta - d;
        if ( fabs(d) <= 0.0001 ) return theta;
    }
    return theta1;
} 
KOKKOS_INLINE_FUNCTION double getAcircularC(double psi)
{
    double theta;
    if ( psi >= 1.0 ) return 1.0;
    if ( psi <= 0.0 ) return 0.0;
    if ( psi <= 1.0e-6 )
    {
        theta = pow(124.4797*psi, 3./13.);
        return theta*theta*theta / 37.6911;
    }
    theta = getThetaOfPsiC(psi);
    return (theta - sin(theta)) / (2.0 * PI);
}

inline void findLinkFlows(int Steps,double dt, Node& Tnode, Link& Tlink, Conduit& Tconduit, Outfall& Toutfall)
    {

        Kokkos::parallel_for(
            Nobjects[LINK], KOKKOS_LAMBDA(int ii) {
                int k = Tlink.subIndex(ii);
                int n1 = Tlink.node1(ii);
                int n2 = Tlink.node2(ii);
                int TNtypee1 = Tnode.typee(n1);
                int TNtypee2 = Tnode.typee(n2);
                double TLyFull = Tlink.yFull(ii);
                int TLtypee = Tlink.typee(ii);
                int TLflowClass = Tlink.flowClass(ii);
                double TNnewDepth1 = Tnode.newDepth(n1);
                double TNnewDepth2 = Tnode.newDepth(n2);
                double TCqMax = Tconduit.qMax(k);
                double TCbeta = Tconduit.beta(k);
                double TLsFull = Tlink.sFull(ii);
                double TLsMax = Tlink.sMax(ii);
                double TLaFull = Tlink.aFull(ii);
                double TLqc = Tlink.qc(ii);
                double TLwMax = Tlink.wMax(ii);
                double TNinvertElev1 = Tnode.invertElev(n1);
                double TNinvertElev2 = Tnode.invertElev(n2);
                double TLsurfArea1 = Tlink.surfArea1(ii);
                double TLsurfArea2 = Tlink.surfArea2(ii);
                double length = Tconduit.modLength(k);
                double TLdqdh = Tlink.dqdh(ii);
                double TLfroude = Tlink.froude(ii);
                double TLnewDepth = Tlink.newDepth(ii);
                double TLnewVolume = Tlink.newVolume(ii);
                double TLnewFlow = Tlink.newFlow(ii);
                double TCroughFactor = Tconduit.roughFactor(k);
                double TLinletControl = Tlink.inletControl(ii);
                double TLnormalFlow = Tlink.normalFlow(ii);
                double TLhasFlapGate = Tlink.hasFlapGate(ii);
                double TLdirection = Tlink.direction(ii);
                double TClength = Tconduit.length(k);
                double TLqLimit = Tlink.qLimit(ii);
                double TCfullState = Tconduit.fullState(k);

                double z1, z2;                    
                double h1, h2;                  
                double y1, y2;                  
                double a1, a2;                 
                double r1;                     
                double yMid, rMid, aMid;        
                double aWtd, rWtd;            
                double qLast;                     
                double qOld;                  
                double aOld;                  
                double v;                     
                double rho;                 
                double sigma;               
                double wSlot;                  
                double dq1, dq2, dq3, dq4, dq5, dq6;    
                double denom;           
                double q;                        
                char   isFull = FALSE;            
                char   isClosed = FALSE;
                double omega = 0.5;             
                double flowDepth1;   
                double flowDepth2;   
                double flowDepthMid;  
                double width1;         
                double width2;        
                double widthMid;     
                double surfArea1 = 0.0; 
                double surfArea2 = 0.0; 
                double criticalDepth = 0.0;   
                double normalDepth = 0.0;    
                double fullDepth;      
                double fasnh = 1.0;   
                
                qOld = Tlink.oldFlow(ii);
                qLast = Tconduit.q1(k);
                z1 = Tnode.invertElev(n1);
                z2 = Tnode.invertElev(n2);
                h1 = Tnode.newDepth(n1) + Tnode.invertElev(n1);
                h2 = Tnode.newDepth(n2) + Tnode.invertElev(n2);
                h1 = max(h1, z1);
                h2 = max(h2, z2);
                y1 = h1 - z1;
                y2 = h2 - z2;
                y1 = max(y1, FUDGE);
                y2 = max(y2, FUDGE);
                y1 = min(y1, Tlink.yFull(ii));
                y2 = min(y2, Tlink.yFull(ii));
                aOld = Tconduit.a2(k);
                aOld = max(aOld, FUDGE);

                flowDepth1 = y1;
                flowDepth2 = y2;
                normalDepth = (flowDepth1 + flowDepth2) / 2.0;
                criticalDepth = normalDepth;
                fullDepth = Tlink.yFull(ii);
                if (flowDepth1 >= fullDepth && flowDepth2 >= fullDepth)
                {
                    Tlink.flowClass(ii) = SUBCRITICAL;
                }
                else Tlink.flowClass(ii) = getFlowClassC(&fasnh,qLast,h1,h2,y1,y2, &criticalDepth, &normalDepth, n1, n2,TNtypee1,TNtypee2,TLflowClass,
            TLyFull, TLaFull, TLqc, TLwMax, TCqMax,TCbeta, TLsFull, TLsMax, TLtypee, TNinvertElev1,TNinvertElev2,TNnewDepth1,TNnewDepth2);
            // printf("the flow calss is : %d\n",Tlink.flowClass(ii));   
            switch (Tlink.flowClass(ii))
                {
                    case SUBCRITICAL:
                        flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
                        if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
                        width1 = getWidthC(flowDepth1, TLyFull,TLwMax);
                        width2 =   getWidthC(flowDepth2, TLyFull,TLwMax);
                        widthMid = getWidthC(flowDepthMid, TLyFull,TLwMax);
                        surfArea1 = (width1 + widthMid) * length / 4.;
                        surfArea2 = (widthMid + width2) * length / 4. * fasnh;
                        break;

                    case UP_CRITICAL:
                        flowDepth1 = criticalDepth;
                        if ( normalDepth < criticalDepth ) flowDepth1 = normalDepth;
                        flowDepth1 = max(flowDepth1, FUDGE);
                        h1 = TNinvertElev1 + flowDepth1;
                        flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
                        if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
                        width2   = getWidthC(flowDepth2, TLyFull,TLwMax);
                        widthMid = getWidthC(flowDepthMid, TLyFull,TLwMax);
                        surfArea2 = (widthMid + width2) * length * 0.5;
                        break;

                    case DN_CRITICAL:
                        flowDepth2 = criticalDepth;
                        if ( normalDepth < criticalDepth ) flowDepth2 = normalDepth;
                        flowDepth2 = max(flowDepth2, FUDGE);
                        h2 = TNinvertElev2  + flowDepth2;
                        width1 = getWidthC(flowDepth1,TLyFull,TLwMax);
                        flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
                        if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
                        widthMid = getWidthC(flowDepthMid, TLyFull,TLwMax);
                        surfArea1 = (width1 + widthMid) * length * 0.5;
                        break;

                    case UP_DRY:
                        flowDepth1 = FUDGE;
                        flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
                        if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
                        width1 = getWidthC(flowDepth1, TLyFull,TLwMax);
                        width2 = getWidthC(flowDepth2, TLyFull,TLwMax);
                        widthMid = getWidthC(flowDepthMid, TLyFull,TLwMax);
                        surfArea2 = (widthMid + width2) * length / 4.;
                        surfArea1 = (width1 + widthMid) * length / 4.;
                        break;

                    case DN_DRY:
                        flowDepth2 = FUDGE;
                        flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
                        if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
                        width1 = getWidthC(flowDepth1, TLyFull,TLwMax);
                        width2 = getWidthC(flowDepth2, TLyFull,TLwMax);
                        widthMid = getWidthC(flowDepthMid, TLyFull,TLwMax);
                        surfArea1 = (widthMid + width1) * length / 4.;
                        surfArea2 = (width2 + widthMid) * length / 4.;
                        break;

                    case DRY:
                        surfArea1 = FUDGE * length / 2.0;
                        surfArea2 = surfArea1;
                        break;
                }

                Tlink.surfArea1(ii) = surfArea1;
                Tlink.surfArea2(ii) = surfArea2;
                y1 = flowDepth1;
                y2 = flowDepth2;

                a1 = getAreaC(y1, TLyFull);
                r1 = xsect_getRofYC(y1, TLyFull);
                a2 = getAreaC(y2, TLyFull);
                yMid = 0.5 * (y1 + y2);
                aMid = getAreaC(yMid,TLyFull);
                rMid = xsect_getRofYC(yMid, TLyFull);
                if ( y1 >= TLyFull && y2 >= TLyFull) isFull = TRUE;
                if (Tlink.flowClass(ii) == DRY ||
                    Tlink.flowClass(ii) == UP_DRY ||
                    Tlink.flowClass(ii) == DN_DRY ||
                    isClosed || aMid <= FUDGE)
                {
                    // printf("aMid is : %f\n",aMid);
                    Tconduit.a1(k) = 0.5 * (a1 + a2);
                    Tconduit.q1(k) = 0.0;
                    Tlink.dqdh(ii)  = GRAV * dt * aMid / length;
                    Tlink.froude(ii) = 0.0;
                    Tlink.newDepth(ii) = min(yMid, TLyFull);
                    Tlink.newVolume(ii) = Tconduit.a1(k) * TClength;
                    Tlink.newFlow(ii) = 0.0;
                    return;
                }

            v = qLast / aMid;
           
            if (fabs(v) > MAXVELOCITY)
                v = MAXVELOCITY * sgn(qLast);
            Tlink.froude(ii) = link_getFroudeC(v,yMid, TLtypee, TLyFull, TLwMax);
            if ( Tlink.flowClass(ii) == SUBCRITICAL &&
                Tlink.froude(ii) > 1.0 ) Tlink.flowClass(ii) = SUPCRITICAL;
            if      ( Tlink.froude(ii) <= 0.5 ) sigma = 1.0;
            else if ( Tlink.froude(ii) >= 1.0 ) sigma = 0.0;
            else    sigma = 2.0 * (1.0 - Tlink.froude(ii));
            rho = 1.0;
            if ( !isFull && qLast > 0.0 && h1 >= h2 ) rho = sigma;
            aWtd = a1 + (aMid - a1) * rho;
            rWtd = r1 + (rMid - r1) * rho;
            sigma =1.0;
            if ( isFull) sigma = 0.0;
            dq1 = dt * Tconduit.roughFactor(k) / pow(rWtd, 1.33333) * fabs(v);
            dq2 = dt *GRAV * aWtd * (h2 - h1) / length;
            dq3 = 0.0;
            dq4 = 0.0;
            if ( sigma > 0.0 )
            {
                dq3 = 2.0 * v * (aMid - aOld) * sigma;
                dq4 = dt * v * v * (a2 - a1) / length * sigma;
            }
            dq5 = 0.0;
            dq6 = 0.0;
            denom = 1.0 + dq1 + dq5;
            q = (qOld - dq2 + dq3 + dq4 + dq6) / denom;
            Tlink.dqdh(ii) = 1.0 / denom  * GRAV * dt * aWtd / length;
            Tlink.inletControl(ii) = FALSE;
            Tlink.normalFlow(ii) = FALSE;
            if ( q > 0.0 )
            {
                if (y1 < TLyFull &&
                        ( Tlink.flowClass(ii) == SUBCRITICAL || 
                        Tlink.flowClass(ii) == SUPCRITICAL ))
                    q = checkNormalFlowC(q, y1, y2, a1, r1,n1,n2,TNtypee1, TNtypee2, TCbeta, TLnormalFlow,TLtypee, TLyFull, TLwMax);
            }
            if ( Steps > 0 )
            {
                q = (1.0 - omega) * qLast + omega * q;
                if ( q * qLast < 0.0 ) q = 0.001 * sgn(q);
            }

            if ( Tlink.qLimit(ii) > 0.0 )
            {
                if ( fabs(q) > Tlink.qLimit(ii) ) q = sgn(q) * Tlink.qLimit(ii);
            }

            if( q >  FUDGE && Tnode.newDepth(n1) <= FUDGE ) q =  FUDGE;
            if( q < -FUDGE && Tnode.newDepth(n2) <= FUDGE ) q = -FUDGE;
            // //znEI理想算例设置，其他算例需删除
            // if (ii==0){
            //       q = 1;
            //   }
            // //
            Tconduit.a1(k) = aMid;
            Tconduit.q1(k) = q;
            Tlink.newDepth(ii)  = min(yMid, Tlink.yFull(ii));
            aMid = (a1 + a2) / 2.0;
            Tconduit.fullState(k) = link_getFullStateC(a1, a2, TLaFull);
            Tlink.newVolume(ii) = aMid * TClength;
            Tlink.newFlow(ii) = q;
            // printf("the new flow is : %f\n",Tlink.newFlow(ii));
        
            });
            Kokkos::fence();
        

    // for (int ii = 0; ii < Nobjects[LINK]; ii++){
    //             int n1 = Tlink.node1(ii);
    //             int n2 = Tlink.node2(ii);
    //             if (Tlink.newFlow(ii) >= 0.0)
    //                 {
    //                     Tnode.outflow(n1) += Tlink.newFlow(ii);
    //                     Tnode.inflow(n2) += Tlink.newFlow(ii);
    //                 }
    //             else
    //                 {
    //                     Tnode.inflow(n1) -= Tlink.newFlow(ii);
    //                     Tnode.outflow(n2) -= Tlink.newFlow(ii);
    //                 }

    //                 Tnode.newSurfArea(n1) += Tlink.surfArea1(ii);
    //                 Tnode.newSurfArea(n2) += Tlink.surfArea2(ii);

    //                 Tnode.sumdqdh(n1) += Tlink.dqdh(ii);
    //                 Tnode.sumdqdh(n2) += Tlink.dqdh(ii);
    // }
            {
            const int nLink = Nobjects[LINK];
            Kokkos::parallel_for(
                nLink,
                KOKKOS_LAMBDA(int ii) {
                    const int n1 = Tlink.node1(ii);
                    const int n2 = Tlink.node2(ii);
                    const double qf = Tlink.newFlow(ii);
                    const double sa1 = Tlink.surfArea1(ii);
                    const double sa2 = Tlink.surfArea2(ii);
                    const double dqdhL = Tlink.dqdh(ii);
                    if (qf >= 0.0) {
                        Kokkos::atomic_add(&Tnode.outflow(n1), qf);
                        Kokkos::atomic_add(&Tnode.inflow(n2), qf);
                    } else {
                        Kokkos::atomic_add(&Tnode.inflow(n1), -qf);
                        Kokkos::atomic_add(&Tnode.outflow(n2), -qf);
                    }
                    Kokkos::atomic_add(&Tnode.newSurfArea(n1), sa1);
                    Kokkos::atomic_add(&Tnode.newSurfArea(n2), sa2);
                    Kokkos::atomic_add(&Tnode.sumdqdh(n1), dqdhL);
                    Kokkos::atomic_add(&Tnode.sumdqdh(n2), dqdhL);
                });
            Kokkos::fence();
        }
}

KOKKOS_INLINE_FUNCTION char link_getFullStateC(double a1, double a2, double aFull)
{
    if ( a1 >= aFull )
    {
        if ( a2 >= aFull ) return ALL_FULL;
        else return UP_FULL;
    }
    if ( a2 >= aFull ) return DN_FULL;
    return 0;

}
KOKKOS_INLINE_FUNCTION int link_setFlapGateC(double q, double TLhasFlapGate, double TLdirection )
{
    if ( TLhasFlapGate) {
        if ( q * TLdirection < 0.0 ) return TRUE;
    }
    return FALSE;
}


KOKKOS_INLINE_FUNCTION double checkNormalFlowC(double q, double y1, double y2, double a1, double r1, 
double n1, double n2, int TNtypee1, int TNtypee2, double TCbeta, double TLnormalFlow,int TLtypee, double TLyFull,double TLwMax)
{
        int    check  = FALSE;
        int    hasOutfall = (TNtypee1 == OUTFALL || TNtypee2 == OUTFALL);
        double qNorm;
        double f1;
        if ( y1 < y2) check = TRUE;
        if ( !check &&
            !hasOutfall )
        {
            if ( y1 > FUDGE && y2 > FUDGE )
            {
                f1 = link_getFroudeC(q/a1, y1, TLtypee, TLyFull,TLwMax);
                if ( f1 >= 1.0 ) check = TRUE;
            }
        }
        if ( check )
        {
            qNorm = TCbeta * a1 * pow(r1, 2./3.);
            if ( qNorm < q )
            {
                TLnormalFlow = TRUE;
                return qNorm;
            }
        }
        return q;
}

KOKKOS_INLINE_FUNCTION double link_getFroudeC(double v, double y, int TLtypee, double TLyFull,double TLwMax)
{
if ( y <= FUDGE ) return 0.0;
if (TLyFull - y <= FUDGE ) return 0.0;
y =getAreaC(y, TLyFull) / xsect_getWofYC(y, TLyFull, TLwMax);
return fabs(v) / sqrt(GRAV * y);
}  

KOKKOS_INLINE_FUNCTION double xsect_getRofYC(double y, double TLyFull)
{
    double yNorm = y / TLyFull; 
    double area = getAreaC(y, TLyFull);
    double r = 0.5*TLyFull;
    if (y >= TLyFull) return area / (2.0*PI*r);
    if (y <= 0.0) return 0.0;
    if (y < r)
    {
        double alpha = acos((r - y) / r);
        return area / (2*alpha*r);
    };
    if (y > r)
    {
        double alpha = acos((y - r) / r);
        return area / (2*PI*r-2*alpha*r);
    };

    return area / (PI*r);
    
}

KOKKOS_INLINE_FUNCTION double getAreaC(double y, double TLyFull)
{

        double r = 0.5*TLyFull;
        if (y >= TLyFull) return PI*r*r;
        if (y <= 0.0) return 0.0;
        if (y < r)
        {
            double alpha = acos((r - y) / r);
            return alpha * r * r - (r - y) * r * sin(alpha);
        }
        if (y > r){
            double alpha = acos((y - r) / r);
            return PI * r * r - alpha * r * r + (y - r) * r * sin(alpha);          
        }

        
        return PI * r * r / 2.0;
        
}


KOKKOS_INLINE_FUNCTION double getWidthC(double y, double TLyFull,double TLwMax  )
{
   
    if (y / TLyFull >= CrownCutoff)
            y = CrownCutoff * TLyFull;
        return xsect_getWofYC(y, TLyFull, TLwMax);  
}

KOKKOS_INLINE_FUNCTION int getFlowClassC(double* fasnh, double qLast,double h1,double h2,double y1,double y2,double *yC, double *yN, int n1, int n2,double TNtypee1,double TNtypee2,double TLflowClass,
double TLyFull, double TLaFull, double TLqc,double TLwMax, double TCqMax,double TCbeta, double TLsFull, double TLsMax, double TLtypee, double TNinvertElev1,double TNinvertElev2,double TNnewDepth1,double TNnewDepth2)
{
        double z1, z2, flowClass;
        double ycMin, ycMax;
        z1 = 0;
        z2 = 0;
        if (TNtypee1 == OUTFALL)
        {
            z1 = max(0.0, (z1 - TNnewDepth1));
        }
            
        if (TNtypee2 == OUTFALL)
        {
            z2 = max(0.0, (z2 - TNnewDepth2));
        }
        flowClass = SUBCRITICAL;
        *fasnh = 1.0;
        if (y1 > FUDGE && y2 > FUDGE)
        {
            if (qLast < 0.0)
            {
                if (z1 > 0.0)
                {
                    *yN = link_getYnormC(fabs(qLast), TLtypee, TCqMax, TCbeta, TLsFull, TLsMax, TLaFull, TLyFull);
                    *yC = xsect_getYcritC(fabs(qLast), TLyFull, TLaFull, TLqc, TLwMax);
                    ycMin = min(*yN, *yC);
                    if (y1 < ycMin)
                        flowClass = UP_CRITICAL;
                }
            }
            else
            {
                if (z2 > 0.0)
                {
                    *yN = link_getYnormC(fabs(qLast), TLtypee, TCqMax, TCbeta, TLsFull, TLsMax, TLaFull, TLyFull);
                    *yC = xsect_getYcritC(fabs(qLast), TLyFull, TLaFull, TLqc, TLwMax);
                    ycMin = min(*yN, *yC);
                    ycMax = min(*yN, *yC);
                    if (y2 < ycMin)
                        flowClass = DN_CRITICAL;
                    else if (y2 < ycMax)
                    {
                        if (ycMax - ycMin < FUDGE)
                            *fasnh = 0.0;
                        else
                            *fasnh = (ycMax - y2) / (ycMax - ycMin);
                    }
                }
            }
        }
        else if (y1 <= FUDGE && y2 <= FUDGE)
            flowClass = DRY;
        else if (y2 > FUDGE)
        {
            if (h2 < TNinvertElev1)
                flowClass = UP_DRY;
            else if (z1 > 0.0)
            {
                *yN = link_getYnormC(fabs(qLast), TLtypee, TCqMax, TCbeta, TLsFull, TLsMax, TLaFull, TLyFull);
                *yC = xsect_getYcritC(fabs(qLast), TLyFull, TLaFull, TLqc, TLwMax);
                flowClass = UP_CRITICAL;
            }
        }
        else
        {
            if (h1 < TNinvertElev2)
                flowClass = DN_DRY;

            else if (z2 > 0.0)
            {
                *yN = link_getYnormC(fabs(qLast), TLtypee, TCqMax, TCbeta, TLsFull, TLsMax, TLaFull, TLyFull);
                *yC = xsect_getYcritC(fabs(qLast), TLyFull, TLaFull, TLqc, TLwMax);
                flowClass = DN_CRITICAL;
            }
        }
        return flowClass;
}

KOKKOS_INLINE_FUNCTION double xsect_getYcritC (double qLast, double TLyFull, double TLaFull, double TLqc, double TLwMax)
{
    double q2g = pow(qLast,2) / GRAV;
    double y, r;
    if (q2g ==0.0)
        return 0.0;
    y = 1.01 * pow(q2g / TLyFull, 1./4.);
    if (y >=TLyFull) y = 0.97 * TLyFull;
    r = TLaFull / (PI / 4.0 * pow(TLyFull, 2));

    if ( r >= 0.5 && r <= 2.0 )
        y = getYcritEnumC(qLast, y, TLyFull, TLaFull, TLqc, TLwMax);
    else y = getYcritRidderC(qLast, y, TLyFull, TLaFull, TLqc, TLwMax);
    return min(y, TLyFull);
}

KOKKOS_INLINE_FUNCTION double getYcritRidderC(double qLast, double y0, double TLyFull, double TLqc, double TLaFull, double TLwMax)
{
        double  y1 = 0.0;
        double  y2 = 0.99 * TLyFull;
        double  yc;
        double q0, q1, q2;
        TLqc = 0.0;
        q2 = getQcriticalC(y2, TLqc,TLyFull,TLaFull,TLwMax);
        if (q2 < qLast ) return TLyFull;
        q0 = getQcriticalC(y0, TLqc,TLyFull,TLaFull,TLwMax);
        q1 = getQcriticalC(0.5*TLyFull, TLqc,TLyFull,TLaFull,TLwMax);

        if ( q0 > qLast )
        {
            y2 = y0;
            if ( q1 < qLast ) y1 = 0.5*TLyFull;
        }else{
            y1 = y0;
            if ( q1 > qLast ) y2 = 0.5*TLyFull;
        }

        TLqc = qLast;
        yc = findroot_RidderC(y1, y2, 0.001, TLqc,TLyFull,TLaFull,TLwMax);
        return yc;
}

KOKKOS_INLINE_FUNCTION double findroot_RidderC(double x1, double x2, double xacc,double TLqc,double TLyFull,double TLaFull,double TLwMax)
{
     int j;
        double ans, fhi, flo, fm, fnew, s, xhi, xlo, xm, xnew;

        flo = getQcriticalC(x1, TLqc,TLyFull,TLaFull,TLwMax);
        fhi = getQcriticalC(x2, TLqc,TLyFull,TLaFull,TLwMax);
        if ( flo == 0.0 ) return x1;
        if ( fhi == 0.0 ) return x2;
        ans = 0.5*(x1+x2);
        if ( (flo > 0.0 && fhi < 0.0) || (flo < 0.0 && fhi > 0.0) )
        {
            xlo = x1;
            xhi = x2;
            for (j=1; j<=MAXIT; j++) {
                xm = 0.5*(xlo + xhi);
                fm = getQcriticalC(xm, TLqc,TLyFull,TLaFull,TLwMax);
                s = sqrt( fm*fm - flo*fhi );
                if (s == 0.0) return ans;
                xnew = xm + (xm-xlo)*( (flo >= fhi ? 1.0 : -1.0)*fm/s );
                if ( fabs(xnew - ans) <= xacc ) break;
                ans = xnew;
                fnew = getQcriticalC(ans, TLqc,TLyFull,TLaFull,TLwMax);
                if ( sign(fm, fnew) != fm)
                {
                    xlo = xm;
                    flo = fm;
                    xhi = ans;
                    fhi = fnew;
                }
                else if ( sign(flo, fnew) != flo )
                {
                    xhi = ans;
                    fhi = fnew;
                }
                else if ( sign(fhi, fnew) != fhi)
                {
                    xlo = ans;
                    flo = fnew;
                }
                else return ans;
                if ( fabs(xhi - xlo) <= xacc ) return ans;
            }
            return ans;
        }
        return -1.e20;
}

KOKKOS_INLINE_FUNCTION double getYcritEnumC (double qLast, double y, double TLyFull, double TLaFull, double TLqc, double TLwMax)
{
    double     q0, dy, qc, yc;
    int        i1, i;
    dy = TLyFull / 25.;
    i1 = (int)(y / dy);
    TLqc = 0.0;
    q0 = getQcriticalC(i1*dy, TLqc,TLyFull,TLaFull,TLwMax);
    if ( q0 < qLast )
    {
        yc = TLyFull;
        for ( i = i1+1; i <= 25; i++)
        {
            qc = getQcriticalC(i*dy, TLqc,TLyFull,TLaFull,TLwMax);
            if ( qc >= qLast )
            {
                yc = ( (qLast-q0) / (qc - q0) + ((double)i-1) ) * dy;
                break;
            }
            q0 = qc;
        }
    }
    else
    {
        yc = 0.0;
        for ( i = i1-1; i >= 0; i--)
        {
            qc = getQcriticalC(i*dy, TLqc,TLyFull,TLaFull,TLwMax);
            if ( qc < qLast )
            {
                yc = ( (qLast-qc) / (q0-qc) + (double)i ) * dy;
                break;
            }
            q0 = qc;
        }
    }
    return yc;
}

KOKKOS_INLINE_FUNCTION double getQcriticalC (double yc, double TLqc, double TLyFull, double TLaFull, double TLwMax)
{
        double a, w, qc;
        a = getAreaC(yc, TLyFull);
        w = xsect_getWofYC(yc, TLyFull, TLwMax);
        qc = -TLqc;
        if ( w > 0.0 )  qc = a * sqrt(GRAV * a / w) - TLqc;
        return qc;
}

KOKKOS_INLINE_FUNCTION double  xsect_getWofYC(double yc, double TLyFull, double TLwMax)
{
    double yNorm = yc / TLyFull;
    double r = 0.5*TLyFull;
    if (yc <= 0.0) return 0.0;
    if (yc >= TLyFull) return TLwMax;
    if (yc < r)
    {
        double alpha = acos((r - yc) / r);
        return 2.0 * r * sin(alpha);
    }
    if (yc > r)
    {
        double alpha = acos((yc - r) / r);
        return 2.0 * r * sin(alpha);
    }
 
    return 2.0 * r;

}



KOKKOS_INLINE_FUNCTION double link_getYnormC (double qLast, int TLtypee, double TCqMax, double TCbeta, double TLsFull, double TLsMax, double TLaFull, double TLyFull)
{
    
    int    k;
    double s, a, y;

    if ( qLast > TCqMax ) qLast = TCqMax;
    if ( qLast <= 0.0 ) return 0.0;
    s = qLast / TCbeta;
    a = xsect_getAofSC(s, TLsFull, TLaFull, TLsMax);
    y = circ_getYofAC(a,TLyFull);
    return y;
}

KOKKOS_INLINE_FUNCTION double circ_getYofAC(double a, double TLyFull )
{
     double r = 0.5*TLyFull;
     double y_low, y_mid, y_high, area_mid;
     int max_iter = 1000;
     double tol = 1.0e-3;

     y_low = 0.0;            
     y_high = 2 * r; 

    if (a <= 0)
        return 0.0;
    if (a == PI * r*r)
        return 2.0*r;

    for (int i = 0; i < max_iter; i++)
    {
        y_mid = (y_low + y_high) / 2;
        area_mid = getAreaC(y_mid, TLyFull);

        if (abs(area_mid - a) < tol)
            return y_mid;
        else if (area_mid < a)
            y_low = y_mid;
        else y_high = y_mid;
    }
    return 0.0;
}

KOKKOS_INLINE_FUNCTION double getYcircularC(double alpha)
{
    double theta;
    if ( alpha >= 1.0 ) return 1.0;
    if ( alpha <= 0.0 ) return 0.0;
    if ( alpha <= 1.0e-5 )
    {
        theta = pow(37.6911*alpha, 1./3.);
        return theta * theta / 16.0;
    }
    theta = getThetaOfAlphaC(alpha);
    return (1.0 - cos(theta/2.)) / 2.0;
}

KOKKOS_INLINE_FUNCTION double getThetaOfAlphaC(double alpha)
{
      int    k;
    double theta, theta1, ap, d;

    if ( alpha > 0.04 ) theta = 1.2 + 5.08 * (alpha - 0.04) / 0.96;
    else theta = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha);
    theta1 = theta;
    ap  = (2.0*PI) * alpha;
    for (k = 1; k <= 40; k++ )
    {
        d = - (ap - theta + sin(theta)) / (1.0 - cos(theta));
        if ( d > 1.0 ) d = sign( 1.0, d );
        theta = theta - d;
        if ( fabs(d) <= 0.0001 ) return theta;
    }
    return theta1;  
}

KOKKOS_INLINE_FUNCTION double xsect_getAofSC (double s, double TLsFull, double TLaFull, double TLsMax )
{
      double psi = s / TLsFull;
    if ( s <= 0.0 ) return 0.0;
    if ( s > TLsMax ) s = TLsMax;

    return circ_getAofSC(s, TLsFull, TLaFull);

}

KOKKOS_INLINE_FUNCTION double circ_getAofSC(double s, double TLsFull, double TLaFull)
{
    double psi = s / TLsFull;
    double x_low, x_high, x_mid, area_mid;
    int max_iter = 1000;
    double tol = 1.0e-2;
    if (psi == 0.0) return 0.0;
    if (psi >= 1.0) return TLaFull;

    if (psi <= 0.015) return TLaFull * getAcircularC(psi);

     x_low = 0.0;            
     x_high = 1; 

    for (int i = 0; i < max_iter; i++)
    {
        x_mid = (x_low + x_high) / 2;
        area_mid = -0.89906279*x_mid*x_mid*x_mid + 1.58653555*x_mid*x_mid + 0.4134494*x_mid - 0.00288671;

        if (abs(area_mid - psi) < tol)
            return TLaFull * x_mid;
        else if (area_mid < psi)
            x_low = x_mid;
        else x_high = x_mid;
    }
    return 0.0;
}

    inline void findSurfArea(int j, double q, double length, double h1, double h2,
                  double y1, double y2, Node& Tnode, Link& Tlink, Conduit& Tconduit)
    {
        int     n1, n2;                    
        double  flowDepth1;                
        double  flowDepth2;                
        double  flowDepthMid;              
        double  width1;                    
        double  width2;                   
        double  widthMid;                  
        double  surfArea1 = 0.0;          
        double  surfArea2 = 0.0;           
        double  criticalDepth;             
        double  normalDepth;             
        double  fullDepth;                 
        double  fasnh = 1.0;               

        n1 = Tlink.node1(j);
        n2 = Tlink.node2(j);
        flowDepth1 = y1;
        flowDepth2 = y2;
        normalDepth = (flowDepth1 + flowDepth2) / 2.0;
        criticalDepth = normalDepth;

        fullDepth = Tlink.yFull(j);
        if (flowDepth1 >= fullDepth && flowDepth2 >= fullDepth)
        {
            Tlink.flowClass(j) = SUBCRITICAL;
        }
        
        else Tlink.flowClass(j) = getFlowClass(j, q, h1, h2, y1, y2,
                                              criticalDepth, normalDepth, fasnh, Tnode, Tlink, Tconduit);
        switch ( Tlink.flowClass(j) )
        {
            
        case SUBCRITICAL:
            flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
            if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
            width1 =   getWidth(j, flowDepth1, Tlink);
            width2 =   getWidth(j, flowDepth2, Tlink);
            widthMid = getWidth(j, flowDepthMid, Tlink);
            surfArea1 = (width1 + widthMid) * length / 4.;
            surfArea2 = (widthMid + width2) * length / 4. * fasnh;
            break;

        case UP_CRITICAL:
            flowDepth1 = criticalDepth;
            if ( normalDepth < criticalDepth ) flowDepth1 = normalDepth;
            flowDepth1 = max(flowDepth1, FUDGE);
            h1 = Tnode.invertElev(n1) + Tlink.offset1(j) + flowDepth1;
            flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
            if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
            width2   = getWidth(j, flowDepth2, Tlink);
            widthMid = getWidth(j, flowDepthMid, Tlink);
            surfArea2 = (widthMid + width2) * length * 0.5;
            break;

        case DN_CRITICAL:
            flowDepth2 = criticalDepth;
            if ( normalDepth < criticalDepth ) flowDepth2 = normalDepth;
            flowDepth2 = max(flowDepth2, FUDGE);
            h2 = Tnode.invertElev(n2) + Tlink.offset2(j) + flowDepth2;
            width1 = getWidth(j, flowDepth1, Tlink);
            flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
            if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
            widthMid = getWidth(j, flowDepthMid, Tlink);
            surfArea1 = (width1 + widthMid) * length * 0.5;
            break;

        case UP_DRY:
            flowDepth1 = FUDGE;
            flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
            if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
            width1 = getWidth(j, flowDepth1, Tlink);
            width2 = getWidth(j, flowDepth2, Tlink);
            widthMid = getWidth(j, flowDepthMid, Tlink);
            surfArea2 = (widthMid + width2) * length / 4.;
            if ( Tlink.offset1(j) <= 0.0 )
            {
                surfArea1 = (width1 + widthMid) * length / 4.;
            }
            break;

        case DN_DRY:
            flowDepth2 = FUDGE;
            flowDepthMid = 0.5 * (flowDepth1 + flowDepth2);
            if ( flowDepthMid < FUDGE ) flowDepthMid = FUDGE;
            width1 = getWidth(j, flowDepth1, Tlink);
            width2 = getWidth(j, flowDepth2, Tlink);
            widthMid = getWidth(j, flowDepthMid, Tlink);
            surfArea1 = (widthMid + width1) * length / 4.;
            if ( Tlink.offset2(j) <= 0.0 )
            {
                surfArea2 = (width2 + widthMid) * length / 4.;
            }
            break;

        case DRY:
            surfArea1 = FUDGE * length / 2.0;
            surfArea2 = surfArea1;
            break;
        }
        Tlink.surfArea1(j) = surfArea1;
        Tlink.surfArea2(j) = surfArea2;
        y1 = flowDepth1;
        y2 = flowDepth2;
    }

    inline int getFlowClass(int j, double q, double h1, double h2, double y1, double y2,
    double yC, double yN, double fasnh, Node& Tnode, Link& Tlink, Conduit& Tconduit)
    
    {
        int    n1, n2;                    
        int    flowClass;                 
        double ycMin, ycMax;               
        double z1, z2;                   
        n1 = Tlink.node1(j);
        n2 = Tlink.node2(j);
        z1 = Tlink.offset1(j);
        z2 = Tlink.offset2(j);
        if ( Tnode.typee(n1) == OUTFALL ) z1 = max(0.0, (z1 - Tnode.newDepth(n1)));
        if ( Tnode.typee(n2) == OUTFALL ) z2 = max(0.0, (z2 - Tnode.newDepth(n2)));

        flowClass = SUBCRITICAL;
        fasnh = 1.0;
        if ( y1 > FUDGE && y2 > FUDGE )
        {
            if ( q < 0.0 )
            {
                if ( z1 > 0.0 )
                {
                    yN = link_getYnorm(j, fabs(q), Tlink, Tconduit);
                    yC = xsect_getYcrit(j, fabs(q), Tlink);
                    ycMin = min(yN, yC);
                    if ( y1 < ycMin ) flowClass = UP_CRITICAL;
                }
            }

            else
            {
                if ( z2 > 0.0 )
                {
                    yN = link_getYnorm(j, fabs(q), Tlink, Tconduit);
                    yC = xsect_getYcrit(j, fabs(q), Tlink);
                    ycMin = min(yN, yC);
                    ycMax = max(yN, yC);
                    if ( y2 < ycMin ) flowClass = DN_CRITICAL;
                    else if ( y2 < ycMax )
                    {
                        if ( ycMax - ycMin < FUDGE ) fasnh = 0.0;
                        else fasnh = (ycMax - y2) / (ycMax - ycMin);
                    }
                }
            }
        }

        else if ( y1 <= FUDGE && y2 <= FUDGE ) flowClass = DRY;
        else if ( y2 > FUDGE )
        {
            if ( h2 < Tnode.invertElev(n1) + Tlink.offset1(j) ) flowClass = UP_DRY;
            else if ( z1 > 0.0 )
            {
                yN = link_getYnorm(j, fabs(q), Tlink, Tconduit);
                yC = xsect_getYcrit(j, fabs(q), Tlink);
                flowClass = UP_CRITICAL;
            }
        }
        else
        {
            if ( h1 < Tnode.invertElev(n2) + Tlink.offset2(j) ) flowClass = DN_DRY;
            else if ( z2 > 0.0 )
            {
                yN = link_getYnorm(j, fabs(q), Tlink, Tconduit);
                yC = xsect_getYcrit(j, fabs(q), Tlink);
                flowClass = DN_CRITICAL;
            }
        }
        return flowClass;
    }

    inline void updateNodeFlows(int i, double dt, Node& Tnode, Link& Tlink)

    {
        int    n1 = Tlink.node1(i);
        int    n2 = Tlink.node2(i);
        double q = Tlink.newFlow(i);
        if ( q >= 0.0 )
        {
            Tnode.outflow(n1) += q;
            Tnode.inflow(n2)  += q;
        }
        else
        {   
            Tnode.inflow(n1)   -= q;
            Tnode.outflow(n2)  -= q;
        }
        Tnode.newSurfArea(n1) += Tlink.surfArea1(i);
        Tnode.newSurfArea(n2) += Tlink.surfArea2(i);
        Tnode.sumdqdh(n1) += Tlink.dqdh(i);
        Tnode.sumdqdh(n2) += Tlink.dqdh(i);
    }

inline double  link_getYnorm(int j, double q, Link& Tlink, Conduit &Tconduit)

{
    int    k;
    double s, a, y;

    q = fabs(q);
    k = Tlink.subIndex(j);
    if ( q > Tconduit.qMax(k) ) q = Tconduit.qMax(k);
    if ( q <= 0.0 ) return 0.0;
    s = q / Tconduit.beta(k);
    a = xxsect.xsect_getAofS(j, s, Tlink);
    y = xxsect.circ_getYofA(j, a, Tlink);
    return y;
}

    inline double getWidth(int j, double y, Link& Tlink)

    {
        if (y / Tlink.yFull(j) >= CrownCutoff && !xxsect.xsect_isOpen(Tlink.typee(j)))
            y = CrownCutoff * Tlink.yFull(j);
        return xxsect.xsect_getWofY(j, y, Tlink);
    }

     
    inline double getArea(int j, double y, double wSlot, Link& Tlink)

    {
        if ( y >= Tlink.yFull(j) ) return Tlink.aFull(j) + (y - Tlink.yFull(j)) * wSlot;
        return xxsect.xsect_getAofY(j, y, Tlink);
    }

    inline double xsect_getYcrit(int j, double q, Link& Tlink)
    {
        double q2g = pow(q,2) / GRAV;
        double y, r;
        if (q2g ==0.0)
            return 0.0;
        y = 1.01 * pow(q2g / Tlink.yFull(j), 1./4.);
        if (y >=Tlink.yFull(j)) y = 0.97 * Tlink.yFull(j);
        r = Tlink.aFull(j) / (PI / 4.0 * pow(Tlink.yFull(j),2));
        if ( r >= 0.5 && r <= 2.0 )
            y = getYcritEnum(j, q, y, Tlink);
        else y = getYcritRidder(j, q, y, Tlink);
        return min(y, Tlink.yFull(j));
        
    }

    inline double getYcritEnum(int j, double q, double y0, Link& Tlink)
{
    double     q0, dy, qc, yc;
    int        i1, i;
    dy = Tlink.yFull(j) / 25.;
    i1 = (int)(y0 / dy);
    Tlink.qc(j) = 0.0;
    q0 = getQcritical(i1*dy, j, Tlink);
    if ( q0 < q )
    {
        yc = Tlink.yFull(j);
        for ( i = i1+1; i <= 25; i++)
        {
            qc = getQcritical(i*dy, j, Tlink);
            if ( qc >= q )
            {
                yc = ( (q-q0) / (qc - q0) + ((double)i-1) ) * dy;
                break;
            }
            q0 = qc;
        }
    }

    else
    {
        yc = 0.0;
        for ( i = i1-1; i >= 0; i--)
        {
            qc = getQcritical(i*dy, j, Tlink);
            if ( qc < q )
            {
                yc = ( (q-qc) / (q0-qc) + (double)i ) * dy;
                break;
            }
            q0 = qc;
        }
    }
    return yc;
}

    inline double getQcritical(double yc, int j, Link& Tlink)

    {
        double a, w, qc;
 
        a = xxsect.xsect_getAofY(j, yc, Tlink);
        w = xxsect.xsect_getWofY(j, yc, Tlink);
        qc = -Tlink.qc(j);
        if ( w > 0.0 )  qc = a * sqrt(GRAV * a / w) - Tlink.qc(j);
        return qc;
    }

    inline double getHydRad(int j, double y, Link& Tlink)

    {
        if (y >= Tlink.yFull(j)) return Tlink.rFull(j);
        return xxsect.xsect_getRofY(j, y, Tlink);
    }

    inline double link_getLength(int j, Link& Tlink, Conduit& Tconduit)

    {
        if ( Tlink.typee(j) == CONDUIT )
            {int k = Tlink.subIndex(j);
            return Tconduit.length(k);}
        return 0.0;
    }

    
    inline   double link_getFroude(int j, double v, double y, Link& Tlink)

    {
        if ( Tlink.typee(j) != CONDUIT ) return 0.0;
        if ( y <= FUDGE ) return 0.0;
        if ( !xxsect.xsect_isOpen(Tlink.typee(j)) &&
            Tlink.yFull(j) - y <= FUDGE ) return 0.0;
        y = xxsect.xsect_getAofY(j, y, Tlink) / xxsect.xsect_getWofY(j, y, Tlink);
        return fabs(v) / sqrt(GRAV * y);
    }


  double checkNormalFlow(int j, double q, double y1, double y2, double a1,
                        double r1, Node& Tnode, Link& Tlink, Conduit& Tconduit)

    {
        int    check  = FALSE;
        int    k = Tlink.subIndex(j);
        int    n1 = Tlink.node1(j);
        int    n2 = Tlink.node2(j);
        int    hasOutfall = (Tnode.typee(n1) == OUTFALL || Tnode.typee(n2) == OUTFALL);
        double qNorm;
        double f1;

        if ( NormalFlowLtd == SLOPE || NormalFlowLtd == BOTH || hasOutfall )
        {
            if ( y1 < y2) check = TRUE;
        }

        if ( !check && (NormalFlowLtd == FROUDE || NormalFlowLtd == BOTH) &&
            !hasOutfall )
        {
            if ( y1 > FUDGE && y2 > FUDGE )
            {
                f1 = link_getFroude(j, q/a1, y1, Tlink);
                if ( f1 >= 1.0 ) check = TRUE;
            }
        }

        if ( check )
        {
            qNorm = Tconduit.beta(k) * a1 * pow(r1, 2./3.);
            if ( qNorm < q )
            {
                Tlink.normalFlow(j) = TRUE;
                return qNorm;
            }
        }
        return q;
    }

    int link_setFlapGate(int j, int n1, int n2, double q, Node& Tnode, Link& Tlink, Outfall& Toutfall)

    {
        if ( Tlink.hasFlapGate(j) )
        {
            if ( q * (double)Tlink.direction(j) < 0.0 ) return TRUE;
        }

        return FALSE;
    }

    char  link_getFullState(double a1, double a2, double aFull)
 
{
    if ( a1 >= aFull )
    {
        if ( a2 >= aFull ) return ALL_FULL;
        else return UP_FULL;
    }
    if ( a2 >= aFull ) return DN_FULL;
    return 0;
}

    inline double getYcritRidder(int j, double q, double y0, Link& Tlink)

    {
        double  y1 = 0.0;
        double  y2 = 0.99 * Tlink.yFull(j);
        double  yc;
        double q0, q1, q2;

        Tlink.qc(j) = 0.0;
        q2 = getQcritical(y2, j, Tlink);
        if (q2 < q ) return Tlink.yFull(j);
        q0 = getQcritical(y0, j, Tlink);
        q1 = getQcritical(0.5*Tlink.yFull(j), j, Tlink);
        if ( q0 > q )
        {
            y2 = y0;
            if ( q1 < q ) y1 = 0.5*Tlink.yFull(j);
        }else{
            y1 = y0;
            if ( q1 > q ) y2 = 0.5*Tlink.yFull(j);
        }
        Tlink.qc(j) = q;
        yc = findroot_Ridder(y1, y2, 0.001, j, Tlink);
        return yc;
    }

    inline double findroot_Ridder(double x1, double x2, double xacc, int jj, Link& Tlink)
    {
        int j;
        double ans, fhi, flo, fm, fnew, s, xhi, xlo, xm, xnew;

        flo = getQcritical(x1, jj, Tlink);
        fhi = getQcritical(x2, jj, Tlink);
        if ( flo == 0.0 ) return x1;
        if ( fhi == 0.0 ) return x2;
        ans = 0.5*(x1+x2);
        if ( (flo > 0.0 && fhi < 0.0) || (flo < 0.0 && fhi > 0.0) )
        {
            xlo = x1;
            xhi = x2;
            for (j=1; j<=MAXIT; j++) {
                xm = 0.5*(xlo + xhi);
                fm = getQcritical(xm, jj, Tlink);
                s = sqrt( fm*fm - flo*fhi );
                if (s == 0.0) return ans;
                xnew = xm + (xm-xlo)*( (flo >= fhi ? 1.0 : -1.0)*fm/s );
                if ( fabs(xnew - ans) <= xacc ) break;
                ans = xnew;
                fnew = getQcritical(ans, jj, Tlink);
                if ( sign(fm, fnew) != fm)
                {
                    xlo = xm;
                    flo = fm;
                    xhi = ans;
                    fhi = fnew;
                }
                else if ( sign(flo, fnew) != flo )
                {
                    xhi = ans;
                    fhi = fnew;
                }
                else if ( sign(fhi, fnew) != fhi)
                {
                    xlo = ans;
                    flo = fnew;
                }
                else return ans;
                if ( fabs(xhi - xlo) <= xacc ) return ans;
            }
            return ans;
        }
        return -1.e20;
    }

 
};

#endif

