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

    /* ------------------------------------------------------------------ */
    /*  [P1-3] Frozen CRS graph gather map (filled once by fillGraph())    */
    /* ------------------------------------------------------------------ */
    intArr rt_row_of;   ///< [nnz] Row (cell idom) owning each non-zero entry.
                        ///< Geometry-only, permanent for the whole run.
    intArr rt_slot_of;  ///< [nnz] RTcoef column that feeds each non-zero
                        ///< entry, so each assembly only needs
                        ///< rt_val(j) = RTcoef(rt_row_of(j), rt_slot_of(j)).
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
     *   4. Fills the row-pointer array @c rt_ptr by cumulatively counting
     *      the exact stencil entries inserted by RTlinear_system().
     *
     * @param gdom  Groundwater domain object providing grid geometry.
     *
     * @note The actual matrix values (rt_val, rt_rhs) are NOT filled here;
     *       they are assembled elsewhere during the transport time step.
     */
    KOKKOS_INLINE_FUNCTION
    int count_stencil_entries(int i, int j, int k, int nx, int ny, int nz) const
    {
        int count = 1; // Diagonal entry.

        if (k > 0 && j > 0) count++;
        if (i > 0 && k > 0) count++;
        if (k > 0) count++;
        if (i < nx - 1 && k > 0) count++;
        if (k > 0 && j < ny - 1) count++;

        if (j > 0 && i > 0) count++;
        if (j > 0) count++;
        if (j > 0 && i < nx - 1) count++;

        if (i > 0) count++;
        if (i < nx - 1) count++;

        if (j < ny - 1 && i > 0) count++;
        if (j < ny - 1) count++;
        if (j < ny - 1 && i < nx - 1) count++;

        if (k < nz - 1 && j > 0) count++;
        if (i > 0 && k < nz - 1) count++;
        if (k < nz - 1) count++;
        if (i < nx - 1 && k < nz - 1) count++;
        if (k < nz - 1 && j < ny - 1) count++;

        return count;
    }

    void init(GwDomain &gdom) {
        int ii, jj, kk, ndom = gdom.nCell;
        rt_nx = gdom.nx;   rt_ny = gdom.ny;   rt_nz = gdom.nz;
        rt_nrow = ndom;  // Matrix dimension equals the number of active cells
        rt_ncol = ndom;  // Square system: nrow == ncol

        // ------------------------------------------------------------------
        // Compute the total number of non-zero entries by summing the exact
        // boundary-checked 19-point stencil used in RTlinear_system().
        // This avoids any mismatch between the CRS row pointer and matrix
        // fill order at 3-D boundaries.
        // ------------------------------------------------------------------
        rt_nnz = 0;
        for (int idx = 0; idx < rt_nrow; idx++) {
            gdom.unpackIndices(idx, kk, jj, ii);
            rt_nnz += count_stencil_entries(ii, jj, kk, rt_nx, rt_ny, rt_nz);
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
        // Build the row-pointer array by cumulative summation of the same
        // stencil counts.  rt_ptr(row+1)-rt_ptr(row) now exactly matches the
        // number of entries written in RTlinear_system() for that row.
        // ------------------------------------------------------------------
        rt_ptr(0) = 0;
        for (int idx = 0; idx < rt_nrow; idx++) {
            gdom.unpackIndices(idx, kk, jj, ii);
            rt_ptr(idx + 1) = rt_ptr(idx) + count_stencil_entries(ii, jj, kk, rt_nx, rt_ny, rt_nz);
        }
    }


    /* ================================================================== */
    /*  [P1-3] One-time CRS graph fill                                     */
    /* ================================================================== */

    /**
     * @brief One-time fill of the CRS graph and the (row, slot) gather map.
     *
     * The stencil connectivity written by the CRS fill stage of
     * RTlinear_system() is geometry-only: the column indices never change
     * during a simulation.  This method replicates the exact boundary-checked
     * 19-entry insertion order of the assembly to fill, for every non-zero j:
     *   - rt_ind(j)    : column index of non-zero j (permanent),
     *   - rt_row_of(j) : row (cell idom) owning non-zero j,
     *   - rt_slot_of(j): RTcoef column that feeds non-zero j.
     *
     * After this call each assembly only gathers values:
     *   rt_val(j) = RTcoef(rt_row_of(j), rt_slot_of(j)),
     * which writes the same value to the same position as the previous
     * branch-based fill (bitwise identical).
     *
     * The insertion order (and thus the slot order) is:
     *   25, 19, 6, 17, 23, 13, 4, 11, 2, 0, 1, 12, 3, 9, 24, 18, 5, 15, 22
     * matching count_stencil_entries().
     *
     * @param gdom  Groundwater domain object providing grid geometry.
     */
	    void fillGraph(GwDomain &gdom) {
	        rt_row_of  = intArr("rt_row_of",  rt_nnz);
	        rt_slot_of = intArr("rt_slot_of", rt_nnz);

	        const int nx = gdom.nx;
	        const int ny = gdom.ny;
	        const int nz = gdom.nz;
	        const int nxy = nx * ny;
	        const int nrow = rt_nrow;
	        auto rt_ptr_l = rt_ptr;
	        auto rt_ind_l = rt_ind;
	        auto rt_row_of_l = rt_row_of;
	        auto rt_slot_of_l = rt_slot_of;

	        Kokkos::parallel_for("RT_Fill_CRS_Graph", nrow, KOKKOS_LAMBDA(int idom) {
	            int ii = idom % nx;
	            int jj = (idom / nx) % ny;
	            int kk = idom / nxy;
	            int irow = rt_ptr_l(idom);

	            if (kk > 0 && jj > 0) {
	                rt_ind_l(irow) = idom - nxy - nx;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 25;  irow++;   // C_{j-1,k-1}
	            }
	            if (ii > 0 && kk > 0) {
	                rt_ind_l(irow) = idom - 1 - nxy;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 19;  irow++;   // C_{i-1,k-1}
	            }
	            if (kk > 0) {
	                rt_ind_l(irow) = idom - nxy;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 6;   irow++;   // C_{k-1}
	            }
	            if (ii < nx-1 && kk > 0) {
	                rt_ind_l(irow) = idom + 1 - nxy;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 17;  irow++;   // C_{i+1,k-1}
	            }
	            if (kk > 0 && jj < ny-1) {
	                rt_ind_l(irow) = idom - nxy + nx;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 23;  irow++;   // C_{k-1,j+1}
	            }
	            if (jj > 0 && ii > 0) {
	                rt_ind_l(irow) = idom - nx - 1;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 13;  irow++;   // C_{j-1,i-1}
	            }
	            if (jj > 0) {
	                rt_ind_l(irow) = idom - nx;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 4;   irow++;   // C_{j-1,i}
	            }
	            if (jj > 0 && ii < nx-1) {
	                rt_ind_l(irow) = idom - nx + 1;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 11;  irow++;   // C_{j-1,i+1}
	            }
	            if (ii > 0) {
	                rt_ind_l(irow) = idom - 1;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 2;   irow++;   // C_{i-1}
	            }
	            // Diagonal: always present
	            rt_ind_l(irow) = idom;
	            rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 0;   irow++;       // C_{i,j,k}
	            if (ii < nx-1) {
	                rt_ind_l(irow) = idom + 1;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 1;   irow++;   // C_{i+1}
	            }
	            if (jj < ny-1 && ii > 0) {
	                rt_ind_l(irow) = idom + nx - 1;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 12;  irow++;   // C_{j+1,i-1}
	            }
	            if (jj < ny-1) {
	                rt_ind_l(irow) = idom + nx;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 3;   irow++;   // C_{j+1,i}
	            }
	            if (jj < ny-1 && ii < nx-1) {
	                rt_ind_l(irow) = idom + nx + 1;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 9;   irow++;   // C_{j+1,i+1}
	            }
	            if (kk < nz-1 && jj > 0) {
	                rt_ind_l(irow) = idom + nxy - nx;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 24;  irow++;   // C_{k+1,j-1}
	            }
	            if (ii > 0 && kk < nz-1) {
	                rt_ind_l(irow) = idom - 1 + nxy;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 18;  irow++;   // C_{i-1,k+1}
	            }
	            if (kk < nz-1) {
	                rt_ind_l(irow) = idom + nxy;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 5;   irow++;   // C_{k+1}
	            }
	            if (ii < nx-1 && kk < nz-1) {
	                rt_ind_l(irow) = idom + 1 + nxy;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 15;  irow++;   // C_{i+1,k+1}
	            }
	            if (kk < nz-1 && jj < ny-1) {
	                rt_ind_l(irow) = idom + nxy + nx;
	                rt_row_of_l(irow) = idom;  rt_slot_of_l(irow) = 22;  irow++;   // C_{k+1,j+1}
	            } });
	        Kokkos::fence();
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
