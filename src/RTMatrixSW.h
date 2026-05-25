/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file RTMatrixSW.h
 * @brief Sparse matrix data structure for surface water reactive transport
 *        simulations (2D horizontal domain).
 *
 * This header defines the RTMatrixSW class, which stores a sparse linear
 * system A * x = b in Compressed Row Storage (CRS / CSR) format for the
 * surface water transport layer.  The sparsity pattern follows an extended
 * 9-point finite-difference stencil on a 2D Cartesian grid (nx x ny),
 * coupling each surface cell to its face-neighbours (N, S, E, W) and
 * edge-neighbours (NE, NW, SE, SW).
 *
 * Compared to the groundwater counterpart (RTMatrix.h), this class operates
 * on a 2D horizontal grid without a vertical dimension, resulting in a
 * simpler stencil and fewer non-zero entries per row.
 *
 * Stencil connectivity (extended 9-point, 2D):
 *   - Corner cells: 3 non-zero entries  (JDY = 3)
 *   - Edge cells:   4 non-zero entries  (BDY = 4)
 *   - Face cells:   5 non-zero entries  (MDY = 5)
 *   - Interior cells: 9 non-zero entries (self + 8 neighbours)
 *
 * Units:
 *   - Concentration c  : mg/L
 *   - Matrix values A  : L/s  (transport coefficients)
 *   - RHS b            : mg/s (mass loading)
 *   - Solution x       : mg/L (concentration)
 *
 * @see RTMatrix.h      Groundwater counterpart (3D grid, 19-point stencil)
 * @see RTSolver-kkpcg.h  Linear solvers that operate on RTMatrix/RTMatrixSW
 */
#ifndef _RT_MATRIX_SW_H_
#define _RT_MATRIX_SW_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "Domain.h"

/**
 * @class RTMatrixSW
 * @brief Sparse matrix in CRS format for the surface-water reactive-transport
 *        linear system.
 *
 * The matrix is square with dimension nCell (number of active surface cells).
 * Non-zero structure follows an extended 9-point stencil on a structured 2D
 * Cartesian grid.  Solver work vectors (r, z, p, q, w) are owned here so
 * that they can be accessed by external solver implementations without
 * additional allocation.
 */
class RTMatrixSW {

public:
    /* ------------------------------------------------------------------ */
    /*  Matrix dimensions                                                  */
    /* ------------------------------------------------------------------ */
    int rt_nrow;    ///< Number of rows    (= nCell, total surface water cells)
    int rt_ncol;    ///< Number of columns (= nCell, the matrix is square)
    int rt_nnz;     ///< Total number of non-zero entries in the sparse matrix
    int rt_nx;      ///< Grid extent in the x-direction [cells]
    int rt_ny;      ///< Grid extent in the y-direction [cells]
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
     * @brief Initialize the sparse matrix structure for the given surface
     *        water domain.
     *
     * This method:
     *   1. Reads the 2D grid dimensions (nx, ny) from the domain.
     *   2. Computes the total number of non-zero entries based on the
     *      extended 9-point stencil connectivity for 1D-x, 1D-y, or 2D-xy
     *      configurations.
     *   3. Allocates all Kokkos views (CRS arrays, diagonal, RHS, solution,
     *      solver work vectors).
     *   4. Fills the row-pointer array @c rt_ptr by calling get_irow() for
     *      each grid cell.
     *
     * @param dom  Surface water domain object providing grid geometry.
     *
     * @note The actual matrix values (rt_val, rt_rhs) are NOT filled here;
     *       they are assembled elsewhere during the transport time step.
     */
    void init(Domain &dom) {
        int ii, jj, ndom = dom.nCell;
        rt_nx = dom.nx;   rt_ny = dom.ny;
        rt_nrow = ndom;  // Matrix dimension equals the number of active cells
        rt_ncol = ndom;  // Square system: nrow == ncol

        // ------------------------------------------------------------------
        // Compute the total number of non-zero entries based on the
        // stencil connectivity for each spatial dimensionality.
        // ------------------------------------------------------------------
        // Number of non-zeros for 1D-x, 1D-y, 2D-xy simulations

        // 1D-x: horizontal line in x.  Interior cells have 3 entries
        // (self + 2 neighbours); boundary cells have 2 entries.
        if (rt_nx > 1 & rt_ny == 1 ) {
            rt_nnz = rt_nrow + (rt_nx-2)*2 + 2;
        }
        // 1D-y: horizontal line in y (symmetric to 1D-x)
        else if (rt_nx == 1 & rt_ny > 1 ) {
            rt_nnz = rt_nrow + (rt_ny-2)*2 + 2;
        }
        // 2D-xy: extended 9-point stencil.
        // Corner cells: 4 entries; edge cells: 6 entries;
        // interior cells: 9 entries.
        else if (rt_nx > 1 & rt_ny > 1 ) {
            rt_nnz = rt_nrow + 4*3 + ((rt_nx-2)*2 + (rt_ny-2)*2)*5 + (rt_nx-2)*(rt_ny-2)*8;
        }
        else {
            printf("ERROR : RTDomainSW must be 1D-x, 1D-y or 2D-xy!\n");
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
            dom.unpackIndices(idx, jj, ii);

            if (idx == rt_nrow)
            {
                // Sentinel: one-past-end index equals total non-zero count
                rt_ptr(idx) = rt_nnz;
            }
            else
            {
                // Compute the starting non-zero index for this row based
                // on the 2D grid position (i, j)
                rt_ptr(idx) = get_irow(ii, jj, rt_nx, rt_ny);
            }
        }
    }


    /* ================================================================== */
    /*  Row-pointer computation                                            */
    /* ================================================================== */

    /**
     * @brief Compute the starting non-zero index (into rt_ind / rt_val) for
     *        the matrix row corresponding to surface water grid cell (i, j).
     *
     * The calculation accounts for the position-dependent stencil size:
     *   - Cells on the domain boundary have fewer connections.
     *   - Interior cells have the full 9-point stencil.
     *
     * The row index is obtained by accumulating the non-zero counts of all
     * cells that precede the target cell in a row-by-row ordering
     * (j-major, then i).
     *
     * Stencil constants (number of non-zero entries per cell type):
     *   - JDY = 3   : corner cells (on a boundary in both x and y)
     *   - BDY = 4   : edge cells   (on a boundary in one direction)
     *   - MDY = 5   : cells adjacent to the boundary (interior face)
     *   - NDY = 5   : (same as MDY, used for consistency)
     *
     * @param i   Cell index in the x-direction [0, nx-1]
     * @param j   Cell index in the y-direction [0, ny-1]
     * @param nx  Total number of cells in x
     * @param ny  Total number of cells in y
     * @return    The starting index into the non-zero arrays for this row.
     */
    inline int get_irow(int i, int j, int nx, int ny) {
    int rt_nrowi = 0, rt_nlayer1 = 0, rt_nlayerk = 0, rt_ncol1 = 0, rt_ncolj = 0, rt_irow;

    // Stencil connectivity constants: non-zero entries per cell type
    int JDY = 3;   ///< Corner cells: 3 non-zeros
    int BDY = 4;   ///< Edge cells:   4 non-zeros
    int MDY = 5;   ///< Face cells:   5 non-zeros
    int NDY = 5;   ///< Interior/face cells: 5 non-zeros (same as MDY for 2D)

        // First column (j == 0): western boundary of the domain
        if (j == 0)
        {
            if ( nx > 1)
            {
                if (ny > 1) {
                    // 2D-xy: cell (i,0) is on the western edge
                    if (i > 0)
                    {
                        // Interior cells along the i-direction on the first row
                        rt_nrowi = JDY + (i - 1) * BDY;
                    }
                }
                else {
                    // 1D-y: only the y dimension; cells at j=0
                    rt_nrowi = (JDY-1) + (nx -1)*(BDY-1);
                }
            }
            else {
                // 1D-x: only the x dimension; single column
                rt_nrowi = JDY -1;
            }
        }
        // Non-first column (j > 0)
        else
        {
            if (nx>1){
                // Accumulate non-zeros from the first column
                rt_nlayer1 = JDY * 2 + BDY * (nx - 2);
                // Accumulate non-zeros from intermediate columns up to j-1
                rt_nlayerk = (i - 1) * BDY * 2 + (nx - 2)*(j - 1) * NDY;

                if (j == ny-1) {
                    // Last column (eastern boundary)
                    if (i>0){
                        rt_nrowi = JDY + (i - 1) * BDY;
                    }
                }
                else {
                    // Interior column
                    rt_nrowi = BDY + (i - 1)*NDY;
                }
            }
            else {
                // 1D-x: non-first column, single cell per column
                rt_nrowi = (JDY - 1)+(j-1)*(BDY-1);
            }
        }

    // Accumulate all contributions to get the final row-pointer offset
    rt_irow = rt_nlayer1 + rt_nlayerk + rt_ncol1 + rt_ncolj + rt_nrowi;
    return rt_irow;

    }


};

#endif
