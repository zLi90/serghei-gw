#ifndef DRAINAGE_STORAGE_H
#define DRAINAGE_STORAGE_H

#include "../define.h"
#include "DrainageState.h"
#include "enums.h"
#include <cmath>

namespace DrainageStorage {

/** Surface area A(d) for storage unit k at depth d (m). */
KOKKOS_INLINE_FUNCTION
real storage_getSurfArea(int k, real d, const Storage& Tstorage)
{
    if (k < 0) return 0.0;
    const int shape = Tstorage.shape(k);
    const real a0 = Tstorage.a0(k);
    const real a1 = Tstorage.a1(k);
    const real a2 = Tstorage.a2(k);

    switch (shape) {
        case FUNCTIONAL_SHAPE:
            return a0 + a1 * pow(d, a2);
        case CYLINDRICAL_SHAPE:
            return a0 + d * (a1 + d * a2);
        default:
            return 0.0;
    }
}

/** Stored volume V(d) for storage unit k at depth d (m). */
KOKKOS_INLINE_FUNCTION
real storage_getVolume(int k, real d, const Storage& Tstorage)
{
    if (k < 0 || d <= 0.0) return 0.0;
    const int shape = Tstorage.shape(k);
    const real a0 = Tstorage.a0(k);
    const real a1 = Tstorage.a1(k);
    const real a2 = Tstorage.a2(k);

    switch (shape) {
        case FUNCTIONAL_SHAPE: {
            const real n = a2 + 1.0;
            return a0 * d + a1 / n * pow(d, n);
        }
        case CYLINDRICAL_SHAPE:
            return d * (a0 + d * (a1 / 2.0 + d * a2 / 3.0));
        default:
            return 0.0;
    }
}

/** Depth from stored volume v (m^3); Newton solve for general FUNCTIONAL. */
KOKKOS_INLINE_FUNCTION
real storage_getDepth(int k, real v, real fullDepth, real fullVolume, const Storage& Tstorage)
{
    if (k < 0 || v <= 0.0) return 0.0;
    if (fullVolume > 0.0 && v >= fullVolume) return fullDepth;

    const int shape = Tstorage.shape(k);
    const real a0 = Tstorage.a0(k);
    const real a1 = Tstorage.a1(k);
    const real a2 = Tstorage.a2(k);

    if (shape == CYLINDRICAL_SHAPE || (shape == FUNCTIONAL_SHAPE && a1 == 0.0)) {
        if (a0 > 0.0) return v / a0;
        return 0.0;
    }

    real y = (fullDepth > 0.0) ? 0.5 * fullDepth : 1.0;
    for (int iter = 0; iter < 20; ++iter) {
        const real f = storage_getVolume(k, y, Tstorage) - v;
        const real df = storage_getSurfArea(k, y, Tstorage);
        if (fabs(df) < TOL12) break;
        const real dy = f / df;
        y -= dy;
        if (y < 0.0) y = 0.0;
        if (fabs(dy) < TOL8) break;
    }
    return y;
}

KOKKOS_INLINE_FUNCTION
real node_getSurfArea(int j, real d, const Node& Tnode, const Storage& Tstorage)
{
    if (Tnode.typee(j) != STORAGE) return 0.0;
    return storage_getSurfArea(Tnode.subIndex(j), d, Tstorage);
}

KOKKOS_INLINE_FUNCTION
real node_getVolume(int j, real d, const Node& Tnode, const Storage& Tstorage)
{
    if (Tnode.typee(j) == STORAGE) {
        if (d >= Tnode.fullDepth(j) && Tnode.fullVolume(j) > 0.0)
            return Tnode.fullVolume(j);
        return storage_getVolume(Tnode.subIndex(j), d, Tstorage);
    }
    if (Tnode.fullDepth(j) > 0.0)
        return Tnode.fullVolume(j) * (d / Tnode.fullDepth(j));
    return 0.0;
}

} // namespace DrainageStorage

#endif
