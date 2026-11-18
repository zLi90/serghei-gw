#ifndef DRAINAGE_PUMP_H
#define DRAINAGE_PUMP_H

#include "../define.h"
#include "DrainageState.h"
#include "enums.h"

namespace DrainagePump {

/** Linear interpolation on a single pump curve (PUMP3: x=head, y=flow). */
KOKKOS_INLINE_FUNCTION
real curveLookupLinear(const PumpCurves& curves, int curveId, real x)
{
    if (curveId < 0 || curveId >= curves.nCurves || curves.nPts <= 0)
        return 0.0;

    const int start = curves.startIdx(curveId);
    const int npt = curves.nPtsPerCurve(curveId);
    if (npt <= 0 || start < 0 || start > curves.nPts - npt) return 0.0;

    const int end = start + npt - 1;
    real x1 = curves.x(start);
    real y1 = curves.y(start);
    if (x <= x1) return y1;

    real x2 = curves.x(end);
    real y2 = curves.y(end);
    if (x >= x2) return y2;

    for (int i = start; i < end; ++i)
    {
        x2 = curves.x(i + 1);
        if (x <= x2)
        {
            x1 = curves.x(i);
            y1 = curves.y(i);
            y2 = curves.y(i + 1);
            if (x2 <= x1) return y1;
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
        }
    }
    return y2;
}

/** Match SWMM table_getSlope, including its endpoint behavior. */
KOKKOS_INLINE_FUNCTION
real curveSlope(const PumpCurves& curves, int curveId, real x)
{
    if (curveId < 0 || curveId >= curves.nCurves) return 0.0;
    const int start = curves.startIdx(curveId);
    const int npt = curves.nPtsPerCurve(curveId);
    if (npt < 2 || start < 0 || start > curves.nPts - npt) return 0.0;
    for (int i = start; i < start + npt - 1; ++i) {
        if (x <= curves.x(i + 1)) {
            const real dx = curves.x(i + 1) - curves.x(i);
            return dx > 0.0 ? (curves.y(i + 1) - curves.y(i)) / dx : 0.0;
        }
    }
    return 0.0;
}

KOKKOS_INLINE_FUNCTION
real pumpHead(int linkIdx, const Node& Tnode, const Link& Tlink)
{
    const int n1 = Tlink.node1(linkIdx);
    const int n2 = Tlink.node2(linkIdx);
    return (Tnode.invertElev(n2) + Tnode.newDepth(n2))
         - (Tnode.invertElev(n1) + Tnode.newDepth(n1));
}

KOKKOS_INLINE_FUNCTION
real pumpIdealFlow(int linkIdx, const Node& Tnode, const Link& Tlink)
{
    const int n1 = Tlink.node1(linkIdx);
    const real setting = Tlink.targetSetting(linkIdx);
    if (setting <= 0.0) return 0.0;

    real qIn = Tnode.inflow(n1) + Tnode.overflow(n1);
    if (qIn < 0.0) qIn = 0.0;
    return qIn * setting;
}

/** TYPE3 pump flow from current node heads (A=node1, B=node2). */
KOKKOS_INLINE_FUNCTION
real pumpType3Flow(int linkIdx, const Node& Tnode, const Link& Tlink, const Pump& Tpump,
                   const PumpCurves& curves)
{
    const int k = Tlink.subIndex(linkIdx);
    const real setting = Tlink.targetSetting(linkIdx);
    if (setting <= 0.0) return 0.0;

    real head = pumpHead(linkIdx, Tnode, Tlink);
    if (head < 0.0) head = 0.0;

    const int cid = Tpump.pumpCurve(k);
    real qIn = curveLookupLinear(curves, cid, head);
    if (qIn < 0.0) qIn = 0.0;
    return qIn * setting;
}

KOKKOS_INLINE_FUNCTION
real pumpLinkFlow(int linkIdx, const Node& Tnode, const Link& Tlink, const Pump& Tpump,
                  const PumpCurves& curves)
{
    const int k = Tlink.subIndex(linkIdx);
    if (Tpump.type(k) == IDEAL_PUMP)
        return pumpIdealFlow(linkIdx, Tnode, Tlink);
    return pumpType3Flow(linkIdx, Tnode, Tlink, Tpump, curves);
}

/** Prevent pumping more than available at inlet node (wet well). */
KOKKOS_INLINE_FUNCTION
real pumpModFlow(int pumpType, int n1, real q, real dt, const Node& Tnode)
{
    if (q <= 0.0 || dt <= 0.0) return 0.0;
    const real remainingInflow = Tnode.inflow(n1) - Tnode.outflow(n1);
    if (Tnode.typee(n1) == STORAGE) {
        if (Tnode.fullVolume(n1) > 0.0) {
            const real qMax = remainingInflow + Tnode.oldVolume(n1) / dt;
            if (q > qMax) q = qMax;
        }
    } else if (pumpType == TYPE3_PUMP) {
        const real area = Tnode.newSurfArea(n1) > MINSURFAREA
                        ? Tnode.newSurfArea(n1) : MINSURFAREA;
        const real dV = 0.5 * (Tnode.oldNetInflow(n1) + remainingInflow - q) * dt;
        if (Tnode.oldDepth(n1) + dV / area <= 0.0 && q > remainingInflow)
            q = remainingInflow;
    } else {
        // Unlike independent IDEAL links, parallel pumps must share the inlet flow.
        const real qMax = remainingInflow + Tnode.overflow(n1);
        if (q > qMax) q = qMax;
    }
    return q > 0.0 ? q : 0.0;
}

/** Run after conduit flows have been accumulated. One device thread preserves
 *  SWMM's link-order dependency without host access to device Views. */
inline void routePumps(int nLinks, real dt, Node& nodes, Link& links,
                       Pump& pumps, const PumpCurves& curves)
{
    Kokkos::parallel_for("drainage_pumps", 1, KOKKOS_LAMBDA(int) {
        for (int j = 0; j < nLinks; ++j) {
            if (links.typee(j) != PUMP) continue;
            const int k = links.subIndex(j);
            const int n1 = links.node1(j);
            const int n2 = links.node2(j);
            links.setting(j) = links.targetSetting(j);
            const real requested = pumpLinkFlow(j, nodes, links, pumps, curves);
            const real q = pumpModFlow(pumps.type(k), n1, requested, dt, nodes);
            real dqdh = 0.0;
            pumps.offCurve(k) = 0;
            if (pumps.type(k) == TYPE3_PUMP && links.setting(j) > 0.0) {
                const real h = pumpHead(j, nodes, links);
                const real head = h > 0.0 ? h : 0.0;
                dqdh = -curveSlope(curves, pumps.pumpCurve(k), head);
                if (head < pumps.xMin(k)) pumps.offCurve(k) = -1;
                if (head > pumps.xMax(k)) pumps.offCurve(k) = 1;
            }
            pumps.limited(k) = q < requested;
            links.newFlow(j) = q;
            links.newDepth(j) = 0.0;
            links.newVolume(j) = 0.0;
            links.dqdh(j) = dqdh;
            links.froude(j) = 0.0;
            links.flowClass(j) = DRY;
            links.surfArea1(j) = 0.0;
            links.surfArea2(j) = 0.0;
            nodes.outflow(n1) += q;
            nodes.inflow(n2) += q;
            nodes.sumdqdh(n1) += dqdh;
            nodes.sumdqdh(n2) += dqdh;
        }
    });
}

inline void updateStatistics(int nLinks, real dt, const Node& nodes,
                              const Link& links, Pump& pumps)
{
    Kokkos::parallel_for("pump_statistics", nLinks, KOKKOS_LAMBDA(int j) {
        if (links.typee(j) != PUMP) return;
        const int k = links.subIndex(j);
        const real q = links.newFlow(j);
        pumps.volume(k) += q * dt;
        if (q > FUDGE) {
            pumps.runTime(k) += dt;
            if (links.oldFlow(j) <= FUDGE) pumps.startUps(k) += 1;
            const real h = pumpHead(j, nodes, links);
            if (h > 0.0) pumps.energy(k) += RHOW * GRAV * q * h * dt / 3.6e6;
            if (pumps.offCurve(k) < 0) pumps.offCurveLow(k) += dt;
            if (pumps.offCurve(k) > 0) pumps.offCurveHigh(k) += dt;
        }
    });
}

/** Startup / shutoff depth limits on pump inlet before a routing step. */
KOKKOS_INLINE_FUNCTION
void pumpSetTargetSetting(int j, const Node& Tnode, const Link& Tlink, const Pump& Tpump)
{
    const int k = Tlink.subIndex(j);
    const int n1 = Tlink.node1(j);
    Tlink.targetSetting(j) = Tlink.setting(j);
    if (Tpump.yOff(k) > 0.0 &&
        Tlink.setting(j) > 0.0 &&
        Tnode.newDepth(n1) < Tpump.yOff(k))
        Tlink.targetSetting(j) = 0.0;
    if (Tpump.yOn(k) > 0.0 &&
        Tlink.setting(j) == 0.0 &&
        Tnode.newDepth(n1) > Tpump.yOn(k))
        Tlink.targetSetting(j) = 1.0;
}

/** Interpolate river stage from a time series (device-callable). */
KOKKOS_INLINE_FUNCTION
real riverStageLookup(const RiverStages& rs, int seriesId, real simTime)
{
    if (seriesId < 0 || seriesId >= rs.nSeries || rs.nPts <= 0)
        return 0.0;

    const int start = rs.startIdx(seriesId);
    const int npt = rs.nPtsPerSeries(seriesId);
    if (npt <= 0) return 0.0;

    const int end = start + npt - 1;
    if (simTime <= rs.time(start)) return rs.stage(start);
    if (simTime >= rs.time(end)) return rs.stage(end);

    for (int i = start; i < end; ++i)
    {
        real t2 = rs.time(i + 1);
        if (simTime <= t2)
        {
            real t1 = rs.time(i);
            real s1 = rs.stage(i);
            real s2 = rs.stage(i + 1);
            if (t2 <= t1) return s1;
            return s1 + (s2 - s1) * (simTime - t1) / (t2 - t1);
        }
    }
    return rs.stage(end);
}

} // namespace DrainagePump

#endif
