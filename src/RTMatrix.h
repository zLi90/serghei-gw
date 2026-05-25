/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file RTMatrix.h
 * @brief Sparse matrix data structure and assembly utilities for groundwater
 *        reactive transport simulations (3D subsurface domain).
 *
 * This header defines the RTMatrix class, which stores a sparse linear system
 * A * x = b in Compressed Row Storage (CRS / CSR) format.  The sparsity
 * pattern is derived from a 19-point finite-difference stencil that couples
 * each groundwater cell to its 6 face-neighbours and 12 edge-neighbours in a
 * structured Cartesian grid.  The class handles:
 *
 *   - Dimensioning the matrix (rows, columns, non-zero count) for 1D-z,
 *     2D-xz, and full 3D-xyz simulations.
 *   - Allocating Kokkos views for the CRS arrays (ptr, ind, val), the
 *     diagonal, the right-hand side, the solution vector, and auxiliary
 *     arrays used by iterative solvers (PCG / BiCGSTAB / GMRES).
 *   - Building the row-pointer array @c rt_ptr so that @c rt_ptr(row)
 *     gives the starting index into @c rt_ind / @c rt_val for that row.
 *
 * Stencil connectivity (19-point, 3D):
 *   - Corner nodes: up to  7 non-zero entries
 *   - Edge nodes:   up to 10 non-zero entries
 *   - Face nodes:   up to 14 non-zero entries
 *   - Interior nodes: 19 non-zero entries (self + 18 neighbours)
 *
 * Units:
 *   - Concentration c  : mg/L
 *   - Matrix values A  : L/s  (transport coefficients)
 *   - RHS b            : mg/s (mass loading)
 *   - Solution x       : mg/L (concentration)
 *
 * @see RTMatrixSW.h   Surface-water counterpart (2D grid, 9-point stencil)
 * @see RTSolverGW.h   Linear solvers that operate on RTMatrix
 */
#ifndef _RT_MATRIX_H_
#define _RT_MATRIX_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"

/**
 * @class RTMatrix
 * @brief Sparse matrix in CRS format for the groundwater reactive-transport
 *        linear system.
 *
 * The matrix is square with dimension nCell (number of active groundwater
 * cells).  Non-zero structure follows a 19-point stencil on a structured
 * Cartesian grid.  Solver work vectors (r, z, p, q, w) are owned here so
 * that they can be accessed by external solver implementations without
 * additional allocation.
 */
class RTMatrix {

public:
    /* ------------------------------------------------------------------ */
    /*  Matrix dimensions                                                  */
    /* ------------------------------------------------------------------ */
    int rt_nrow;  ///< Number of rows    (= nCell, the total number of groundwater cells)
    int rt_ncol;  ///< Number of columns (= nCell, the matrix is square)
    int rt_nnz;   ///< Total number of non-zero entries in the sparse matrix
    int rt_nx;    ///< Grid extent in the x-direction [cells]
    int rt_ny;    ///< Grid extent in the y-direction [cells]
    int rt_nz;    ///< Grid extent in the z-direction [cells]
    int rt_cg_iter; ///< Number of CG iterations performed in the last solve

    /* ------------------------------------------------------------------ */
    /*  CRS (Compressed Row Storage) arrays                                */
    /* ------------------------------------------------------------------ */
    intArr rt_ptr;   ///< Row-pointer array [nrow+1].  rt_ptr(i) gives the
                     ///< starting index into rt_ind / rt_val for row i.
                     ///< rt_ptr(nrow) = nnz (one-past-end sentinel).
    intArr rt_ind;   ///< Column-index array [nnz].  For each non-zero entry,
                     ///< stores the column index of that entry.
    realArr rt_val;  ///< Value array [nnz].  Stores the numerical value of
                     ///< each non-zero matrix entry (transport coefficients).
    realArr rt_diag; ///< Diagonal values [nrow].  Inverse of the diagonal
                     ///< entries, used as the Jacobi preconditioner.

    /* ------------------------------------------------------------------ */
    /*  Right-hand side and solution vectors                               */
    /* ------------------------------------------------------------------ */
    realArr rt_rhs;  ///< Right-hand side vector b [nrow]  [mg/s]
    realArr rt_x;    ///< Solution vector x (concentration) [nrow]  [mg/L]

    /* ------------------------------------------------------------------ */
    /*  Transpose CRS arrays (used by some preconditioners)                */
    /* ------------------------------------------------------------------ */
    intArr rt_ptrT;  ///< Row-pointer array of the transpose matrix
    intArr rt_indT;  ///< Column-index array of the transpose matrix
    realArr rt_valT; ///< Value array of the transpose matrix
    realArr rt_lt;   ///< Lower-triangular factor (ILU-type preconditioners)
    realArr rt_ut;   ///< Upper-triangular factor (ILU-type preconditioners)

    /* ------------------------------------------------------------------ */
    /*  Work vectors for iterative solvers (PCG / BiCGSTAB / GMRES)        */
    /* ------------------------------------------------------------------ */
    realArr rt_r;  ///< Residual vector  r = b - A*x  [nrow]
    realArr rt_z;  ///< Preconditioned residual  z = M^{-1} * r  [nrow]
    realArr rt_p;  ///< Search direction vector (CG)  [nrow]
    realArr rt_q;  ///< Matrix-vector product  q = A * p  [nrow]
    realArr rt_w;  ///< Auxiliary work vector (GMRES Arnoldi)  [nrow]


    /* ================================================================== */
    /*  Initialization                                                     */
    /* ================================================================== */

    /**
     * @brief Initialize the sparse matrix structure for the given groundwater
     *        domain.
     *
     * This method:
     *   1. Reads the grid dimensions (nx, ny, nz) from the domain.
     *   2. Computes the total number of non-zero entries based on the stencil
     *      connectivity for 1D-z, 2D-xz, or 3D-xyz configurations.
     *   3. Allocates all Kokkos views (CRS arrays, diagonal, RHS, solution,
     *      solver work vectors).
     *   4. Fills the row-pointer array @c rt_ptr by calling get_irow() for
     *      each grid cell.
     *
     * @param gdom  Groundwater domain object providing grid geometry.
     *
     * @note The actual matrix values (rt_val, rt_rhs) are NOT filled here;
     *       they are assembled elsewhere during the transport time step.
     */
    void init(GwDomain &gdom) {
        int ii, jj, kk, ndom = gdom.nCell;
        rt_nx = gdom.nx;   rt_ny = gdom.ny;   rt_nz = gdom.nz;
        rt_nrow = ndom;  // Matrix dimension equals the number of active cells
        rt_ncol = ndom;  // Square system: nrow == ncol

        // ------------------------------------------------------------------
        // Compute the total number of non-zero entries based on the
        // stencil connectivity for each spatial dimensionality.
        // ------------------------------------------------------------------
        // Number of non-zeros for 1D-z, 2D-xz, 3D-xyz simulations
        if (rt_nx == 1 & rt_ny == 1 & rt_nz > 1) {
            // 1D-z: vertical column only.  Interior cells have 3 entries
            // (self + 1 above + 1 below); boundary cells have 2 entries.
            rt_nnz = rt_nrow + (rt_nz-2)*2 + 2;
        }
        else if (rt_nx > 1 & rt_ny == 1 & rt_nz > 1) {
            // 2D-xz: extended stencil on a vertical cross-section.
            // Corner cells: 4 entries; edge cells: 6 entries;
            // interior cells: 9 entries.
            rt_nnz = rt_nrow + 4*3 + ((rt_nx-2)*2 + (rt_nz-2)*2)*5 + (rt_nx-2)*(rt_nz-2)*8;
        }
        else if (rt_nx > 1 & rt_ny > 1 & rt_nz > 1) {
            // 3D-xyz: full 19-point stencil
            //   8 corner  cells: 7 non-zeros each
            //   edge      cells: 10 non-zeros each (on 12 edges)
            //   face      cells: 14 non-zeros each (on 6 faces)
            //   interior  cells: 19 non-zeros each
            rt_nnz = rt_nrow
             + 8*6  // 8 corner points, each with 7 non-zero entries (self + 6 neighbours)
             + ((rt_nx-2)*4 + (rt_ny-2)*4 + (rt_nz-2)*4)*9  // edge points: each edge has (nx-2) points, each with 10 non-zero entries
             + ((rt_nx-2)*(rt_ny-2)*2 + (rt_nx-2)*(rt_nz-2)*2 + (rt_ny-2)*(rt_nz-2)*2)*13  // face points: each face has (nx-2)*(ny-2) points, each with 14 non-zero entries
             + (rt_nx-2)*(rt_ny-2)*(rt_nz-2)*18;  // interior points: each with 19 non-zero entries (19-point finite-difference stencil)
        }
        else {
            printf("ERROR : RTDomain must be 1D-z, 2D-xz or 3D-xyz!\n");
        }

        // ------------------------------------------------------------------
        // Allocate Kokkos views for the CRS arrays, RHS, solution, and
        // solver work vectors.
        // ------------------------------------------------------------------
        rt_ptr = intArr("rt_ptr", rt_nrow+1);      rt_ind = intArr("rt_ind", rt_nnz);
        rt_val = realArr("rt_val", rt_nnz);          rt_diag = realArr("rt_diag", rt_nrow);
        rt_rhs = realArr("rt_rhs", rt_nrow);         rt_x = realArr("rt_x", rt_nrow);
        rt_r = realArr("rt_r", rt_nrow);             rt_z = realArr("rt_z", rt_nrow);
        rt_p = realArr("rt_p", rt_nrow);             rt_q = realArr("rt_q", rt_nrow);
        rt_w = realArr("rt_w", rt_nrow);

        // ------------------------------------------------------------------
        // Build the row-pointer array rt_ptr.
        // For each row idx, rt_ptr(idx) stores the cumulative count of
        // non-zero entries in all preceding rows, which serves as the
        // starting index into rt_ind / rt_val for that row.
        // The sentinel rt_ptr(nrow) = nnz.
        // ------------------------------------------------------------------
        for (int idx = 0; idx < rt_nrow+1; idx++) {
            gdom.unpackIndices(idx, kk, jj, ii);
            if (idx == rt_nrow)
            {
                // Sentinel: one-past-end index equals total non-zero count
                rt_ptr(idx) = rt_nnz;
            }
            else
            {
                // Compute the starting non-zero index for this row based
                // on the grid position (i, j, k)
                rt_ptr(idx) = get_irow(ii, jj, kk, rt_nx, rt_ny, rt_nz);
            }
        }
    }


    /* ================================================================== */
    /*  Row-pointer computation                                            */
    /* ================================================================== */

    /**
     * @brief Compute the starting non-zero index (into rt_ind / rt_val) for
     *        the matrix row corresponding to grid cell (i, j, k).
     *
     * The calculation accounts for the position-dependent stencil size:
     *   - Cells on the domain boundary have fewer connections.
     *   - Interior cells have the full 19-point stencil.
     *
     * The row index is obtained by accumulating the non-zero counts of all
     * cells that precede the target cell in a layer-by-layer, row-by-row
     * ordering (k-major, then j, then i).
     *
     * Stencil constants (number of non-zero entries per cell type):
     *   - JDY = 7   : corner cells in 3D (on a boundary in all 3 axes)
     *   - BDY = 10  : edge cells in 3D   (on a boundary in 2 of 3 axes)
     *   - MDY = 14  : face cells in 3D   (on a boundary in 1 of 3 axes)
     *   - NDY = 19  : interior cells      (not on any boundary)
     *
     * @param i   Cell index in the x-direction [0, nx-1]
     * @param j   Cell index in the y-direction [0, ny-1]
     * @param k   Cell index in the z-direction [0, nz-1]
     * @param nx  Total number of cells in x
     * @param ny  Total number of cells in y
     * @param nz  Total number of cells in z
     * @return    The starting index into the non-zero arrays for this row.
     */
    inline int get_irow(int i, int j, int k, int nx, int ny, int nz) {
    int rt_nrowi = 0, rt_nlayer1 = 0, rt_nlayerk = 0, rt_ncol1 = 0, rt_ncolj = 0, rt_irow;

    // Stencil connectivity constants: non-zero entries per cell type
    int JDY = 7;   ///< Corner cells: 7 non-zeros (self + 6 neighbours)
    int BDY = 10;  ///< Edge cells:   10 non-zeros
    int MDY = 14;  ///< Face cells:   14 non-zeros
    int NDY = 19;  ///< Interior cells: 19 non-zeros (full stencil)

    // ==================================================================
    //  First layer (k == 0): cells on the bottom of the domain
    // ==================================================================
    if (k == 0)
    {
        // --------------------------------------------------------------
        //  First layer, first column (j == 0)
        // --------------------------------------------------------------
        if (j == 0)
        {
            // 3D case (ny > 1): (i,0,0) is a corner cell
            if (ny > 1)
            {
                if (i > 0)
                {
                    // Interior cells along the i-direction on the first row
                    rt_nrowi = JDY + (i - 1) * BDY;
                }
            }
            // 2D x-z case (ny == 1): only x and z dimensions
            else
            {
                // Cells in the first row (not the origin cell)
                if (i > 0)
                {
                    rt_nrowi = 4 + (i - 1) * 6;
                }
            }
        }
        // --------------------------------------------------------------
        //  First layer, non-first column (j > 0)
        // --------------------------------------------------------------
        else
        {
            // Total non-zeros contributed by the first column of this layer
            rt_ncol1 = JDY*2 + BDY*(nx-2);
            // Total non-zeros contributed by previous columns in this layer
            rt_ncolj = (j - 1) * BDY*2 + (j - 1) * MDY*(nx-2);

            if (j == ny - 1)
            {
                // Last column (far boundary in y)
                if (i > 0)
                {
                    rt_nrowi = JDY + (i - 1) * BDY;
                }
            }
            else
            {
                // Interior column
                if (i > 0)
                {
                    rt_nrowi = BDY + (i - 1) * MDY;
                }
            }
        }
    }
    // ==================================================================
    //  Interior and top layers (k > 0)
    // ==================================================================
    else
    {
        if (nx > 1)
        {
            // 3D case
            if (ny > 1)
            {
                // Non-zero count contributed by the entire first layer (k=0)
                rt_nlayer1 = 4*JDY + BDY*2*(nx-2+ny-2) + MDY*(nx-2)*(ny-2);
                // Non-zero count contributed by all interior layers above k=0
                rt_nlayerk = (k - 1) * (JDY*4 + MDY*2*(nx-2+ny-2) + NDY*(nx-2)*(ny-2));
            }
            // 2D x-z case
            else
            {
                rt_nlayer1 = 4 * 2 + 6 * (nx - 2);
                rt_nlayerk = (k - 1) * (6 * 2 + 9 * (nx - 2));
            }
        }
        // 1D-z case (nx == 1)
        else
        {
            rt_nlayer1 = 2;
            rt_nlayerk = (k - 1) * 3;
        }

        // --------------------------------------------------------------
        //  Top layer (k == nz-1): cells on the top of the domain
        // --------------------------------------------------------------
        if (k == nz - 1)
        {
            // Top layer, first column (j == 0)
            if (j == 0)
            {
                // 3D
                if (ny > 1)
                {
                    if (i > 0)
                    {
                        rt_nrowi = 4 + (i - 1) * 5;
                    }
                }
                // 2D x-z
                else
                {
                    if (i > 0)
                    {
                        rt_nrowi = 4 + (i - 1) * 6;
                    }
                }
            }
            // Top layer, non-first column (j > 0)
            else
            {
                // Total non-zeros from the first column of the top layer
                rt_ncol1 = JDY*2 + BDY*(nx-2);
                // Total non-zeros from previous columns in the top layer
                rt_ncolj = (j - 1) * BDY*2 + (j - 1) * MDY*(nx-2);
                if (j == ny - 1)
                {
                    // Last column (far y-boundary)
                    if (i > 0)
                    {
                        rt_nrowi = JDY + (i - 1) * BDY;
                    }
                }
                else
                {
                    // Interior column in the top layer
                    if (i > 0)
                    {
                        rt_nrowi = BDY + (i - 1) * MDY;
                    }
                }
            }
        }
        // --------------------------------------------------------------
        //  Interior layers (0 < k < nz-1): neither bottom nor top
        // --------------------------------------------------------------
        else
        {
            // Interior layer, first column (j == 0)
            if (j == 0)
            {
                // 3D
                if (ny > 1)
                {
                    if (i > 0)
                    {
                        rt_nrowi = BDY + (i - 1) * MDY;
                    }
                }
                // 2D x-z
                else
                {
                    if (i > 0)
                    {
                        rt_nrowi = 6 + (i - 1) * 9;
                    }
                }
            }
            // Interior layer, non-first column (j > 0)
            else
            {
                // Total non-zeros from the first column of this interior layer
                rt_ncol1 = BDY*2 + MDY*(nx-2);
                // Total non-zeros from previous columns in this interior layer
                rt_ncolj = (j - 1) * MDY*2 + (j - 1) * NDY*(nx-2);
                if (j == ny - 1)
                {
                    // Last column in an interior layer
                    if (i > 0)
                    {
                        rt_nrowi = BDY + (i - 1) * MDY;
                    }
                }
                else
                {
                    // Interior column in an interior layer
                    if (i > 0)
                    {
                        rt_nrowi = MDY + (i - 1) * NDY;
                    }
                }
            }
        }
    }
    // Accumulate all contributions to get the final row-pointer offset
    rt_irow = rt_nlayer1 + rt_nlayerk + rt_ncol1 + rt_ncolj + rt_nrowi;
    return rt_irow;

    }


};


#endif
