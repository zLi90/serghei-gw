/**
 * @file Afgen.h
 * @brief AFGEN (Arbitrary Function GENerator) linear interpolation utility.
 *
 * Provides piecewise-linear interpolation over tabulated parameter data.
 * In the WOFOST crop growth model, many crop parameters (e.g., partitioning
 * fractions, temperature response curves) are specified as ordered pairs
 * (x, y) in lookup tables. This class performs efficient linear interpolation
 * on such tables and is marked KOKKOS_INLINE_FUNCTION for use on both
 * host and device (GPU) code paths.
 *
 * Table format: a flat 1D array of alternating x,y values:
 *   [x0, y0, x1, y1, x2, y2, ...]
 * where x0 < x1 < x2 < ... (monotonically increasing abscissae).
 *
 * Boundary behaviour:
 *   - If x <= x0, returns y0 (clamp to first value).
 *   - If x >= x_{n-1}, returns y_{n-1} (clamp to last value).
 *
 * Original reference: WOFOST 7.1 documentation, "AFGEN" subroutine.
 *
 * @see CropState.h for table storage (std::map of Kokkos Views).
 * @see CropInit.h for table parsing from input files.
 */

#ifndef _AFGEN_H_
#define _AFGEN_H_

#include "../define.h"

class Afgen {
public:
    /**
     * @brief Perform linear interpolation on a tabulated parameter function.
     *
     * Given a table of monotonically-ordered (x, y) pairs stored as a flat
     * Kokkos 1D View [x0, y0, x1, y1, ...], return the interpolated y value
     * at the given x using piecewise-linear interpolation.
     *
     * @param table  Kokkos View containing alternating x,y pairs.
     *               Must have even extent >= 2 (i.e., at least one x,y pair).
     * @param x      The abscissa value at which to evaluate the function.
     * @return       Interpolated ordinate value y(x). Returns the boundary
     *               y-value if x falls outside the table range.
     *
     * @note Time complexity is O(n/2) where n = table extent.
     *       Tables in WOFOST typically have 3--10 pairs, so linear search is
     *       efficient and GPU-friendly (no binary search overhead).
     */
    KOKKOS_INLINE_FUNCTION
    static real lookup(const realArr& table, real x) {
        int n = table.extent(0);

        /* Guard: empty or single-element table cannot be interpolated */
        if (n < 2) return 0.0;

        /* --- Boundary clamping --- */
        if (x <= table(0))     return table(1);      /* below first x -> return first y */
        if (x >= table(n - 2)) return table(n - 1);  /* above last  x -> return last  y */

        /* --- Linear search for the enclosing interval [x_i, x_{i+1}) --- */
        for (int i = 0; i < n - 2; i += 2) {
            if (x >= table(i) && x < table(i + 2)) {
                real x1 = table(i);
                real y1 = table(i + 1);
                real x2 = table(i + 2);
                real y2 = table(i + 3);

                real slope = (y2 - y1) / (x2 - x1);
                return y1 + slope * (x - x1);
            }
        }

        /* Fallback: should not be reached if table is well-formed */
        return table(n - 1);
    }
};

#endif /* _AFGEN_H_ */
