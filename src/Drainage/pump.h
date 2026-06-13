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
    if (npt <= 0) return 0.0;

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
    const int n1 = Tlink.node1(linkIdx);
    const int n2 = Tlink.node2(linkIdx);
    (void)n1;

    const real setting = Tlink.targetSetting(linkIdx);
    if (setting <= 0.0) return 0.0;

    real head = (Tnode.invertElev(n2) + Tnode.newDepth(n2))
              - (Tnode.invertElev(n1) + Tnode.newDepth(n1));
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
    if (q <= 0.0) return q;

    if (pumpType == IDEAL_PUMP && Tnode.typee(n1) != STORAGE) return q;
    real qMax = Tnode.inflow(n1);
    if (Tnode.fullVolume(n1) > 0.0)
        qMax += Tnode.oldVolume(n1) / dt;
    if (q > qMax) q = qMax;
    return q > 0.0 ? q : 0.0;
}

/** Startup / shutoff depth limits on pump inlet (host, before routing step). */
inline void pumpSetTargetSetting(int j, Node& Tnode, Link& Tlink, Pump& Tpump)
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
