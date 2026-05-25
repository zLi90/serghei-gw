/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file RTSolver-kkpcg.h
 * @brief KokkosKernels-based Preconditioned Conjugate Gradient (PCG) solver
 *        and manual CG / Jacobi solver for the groundwater reactive-transport
 *        sparse linear system.
 *
 * This header provides the RTSolverGW template class with two solver
 * backends:
 *
 *   1. kkpcg() : Delegates to KokkosKernels' pcgsolve() with a Gauss-Seidel
 *                preconditioner.  This is the high-performance path that
 *                leverages GPU-optimised sparse linear algebra kernels.
 *
 *   2. cg()    : A hand-coded Preconditioned Conjugate Gradient solver with
 *                Jacobi (diagonal) preconditioning.  This serves as a
 *                reference implementation and a fallback when KokkosKernels
 *                is not available.
 *
 * Both solvers operate on an RTMatrix object stored in Compressed Row Storage
 * (CRS / CSR) format.
 *
 * Physical units:
 *   - Matrix A entries : transport coefficients  [1/s]
 *   - RHS b            : mass loading rate       [mg/s]
 *   - Solution x       : solute concentration    [mg/L]
 *
 * @see RTMatrix.h       Sparse matrix storage (CRS format)
 * @see RTSolverGW.h     Alternative solver with GMRES support
 */
#ifndef _RT_SOLVER_H_
#define _RT_SOLVER_H_

#include "KokkosKernels_config.h"
#include "KokkosSparse_pcg.hpp"
#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_default_types.hpp"
#include <iostream>
#include "GwDomain.h"
#include "RTMatrix.h"
#include "RTMatrix.h"
#include <vector>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosKernels_Handle.hpp"
#include <KokkosSparse_spiluk.hpp>
#include <KokkosSparse_sptrsv.hpp"
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas.hpp>
#include <KokkosSparse_gauss_seidel.hpp>
#include <KokkosBlas1_dot.hpp>
#include <cmath>
#include <memory>
#include <fstream>

#include "KokkosSparse_bicgstab.hpp"

/**
 * @brief View type aliases for KokkosKernels CRS matrix interfacing.
 *
 * These aliases match the view types expected by the KokkosKernels
 * solver routines.  They are defined at file scope because the
 * KokkosKernels API expects them as template arguments.
 */
using lno_view_t = Kokkos::View<int *>;
using lno_nnz_view_t = Kokkos::View<int *>;
using scalar_view_t = Kokkos::View<double *>;

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosKernels::Experimental::Example;

/**
 * @class RTSolverGW
 * @brief Template class providing the KokkosKernels PCG solver and a
 *        manual CG solver for the groundwater reactive-transport system.
 *
 * @tparam execution_space  The Kokkos execution space (e.g. Serial, OpenMP,
 *                          Cuda, HIP).  All parallel kernels are dispatched
 *                          on this space.
 */
template <typename execution_space>
class RTSolverGW
{
    // ==================================================================
    //  Type aliases derived from the KokkosSparse CrsMatrix template.
    //  These define the view types used by KokkosKernels solver backends.
    // ==================================================================
    typedef typename KokkosSparse::CrsMatrix<default_scalar, default_lno_t, execution_space, void, default_size_type> crsMat_t;
    typedef typename crsMat_t::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
    typedef typename scalar_view_t::value_type scalar_t;
    typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type::value_type size_type;
    typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type::value_type lno_t;

public:
    /* ------------------------------------------------------------------ */
    /*  Solver control parameters                                          */
    /* ------------------------------------------------------------------ */
    int iter;       ///< Iteration count of the last solver run
    int iter_max;   ///< Maximum allowed iterations (default: 400 * 20 = 8000)
    int gsteps;     ///< Number of Gauss-Seidel / restart groups (default: 20)
    int nrow;       ///< Matrix dimension (number of rows), set during init
    int nnz;        ///< Total number of non-zero entries, set during init
    real eps;       ///< Achieved residual norm from the last solve
    real eps_min;   ///< Convergence tolerance for the residual norm (default: 1e-7)

    /* ------------------------------------------------------------------ */
    /*  KokkosKernels views for interfacing with the library PCG solver    */
    /* ------------------------------------------------------------------ */
    lno_view_t rt_ptr;       ///< Row-pointer view (KokkosKernels format)
    lno_nnz_view_t rt_ind;   ///< Column-index view (KokkosKernels format)
    scalar_view_t rt_val;    ///< Matrix values view (KokkosKernels format)
    scalar_view_t rt_rhs;    ///< Right-hand side view (KokkosKernels format)
    scalar_view_t rt_vecx;   ///< Solution vector view (KokkosKernels format)
    scalar_view_t rt_diag;   ///< Diagonal values view (KokkosKernels format)

    /* ================================================================== */
    /*  Initialization                                                     */
    /* ================================================================== */

    /**
     * @brief Initialize the solver with matrix dimensions and copy the
     *        row-pointer array from the RTMatrix object.
     *
     * Sets default solver parameters (iter_max = 8000, eps_min = 1e-7) and,
     * if SERGHEI_KOKKOSKERNELS_SOLVER is enabled, allocates KokkosKernels
     * views and deep-copies the row-pointer array.
     *
     * @param rtA   The sparse matrix object (provides dimensions).
     * @param gdom  The groundwater domain (provides grid info).
     */
    inline void init(RTMatrix &rtA, GwDomain &gdom)
    {
        nrow = rtA.rt_nrow;
        nnz = rtA.rt_nnz;
        iter_max = 400;
        gsteps = 20;
        iter_max = iter_max * gsteps;
        eps_min = 1e-7;
#if SERGHEI_KOKKOSKERNELS_SOLVER
        rt_ptr = lno_view_t("rt_ptr", rtA.rt_nrow + 1);      // Starting index of non-zeros for each row
        rt_ind = lno_nnz_view_t("rt_ind", rtA.rt_nnz);        // Column indices of non-zero entries
        rt_val = scalar_view_t("rt_val", rtA.rt_nnz);          // Values of non-zero entries
        rt_rhs = scalar_view_t("rt_rhs", rtA.rt_nrow);         // Right-hand side vector b
        rt_vecx = scalar_view_t("rt_vecx", rtA.rt_nrow);       // Approximate solution vector x
        rt_diag = scalar_view_t("rt_diag", rtA.rt_nrow);       // Diagonal entries (for preconditioning)
        Kokkos::deep_copy(rt_ptr, rtA.rt_ptr);
#endif
    }

    /* ================================================================== */
    /*  Diagnostic utilities                                               */
    /* ================================================================== */

    /**
     * @brief Print the full dense representation of the sparse matrix along
     *        with the RHS and solution vectors.
     *
     * This method converts the CRS matrix to dense format and prints all
     * entries (zero and non-zero).  Intended for small-scale debugging only.
     *
     * @param rtA  The sparse matrix to print.
     */
    void printRTMatrix(const RTMatrix &rtA)
    {
        int n = rtA.rt_nrow; // Number of matrix rows
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                double val = 0.0; // Default value for absent (zero) entries
                for (int k = rtA.rt_ptr[i]; k < rtA.rt_ptr[i + 1]; ++k)
                {
                    if (rtA.rt_ind(k) == j)
                    {
                        val = rtA.rt_val(k); // Found a non-zero entry
                        break;               // No need to search further in this row
                    }
                }
                std::cout << "A[" << i << "][" << j << "] = " << val << " ";
            }
            std::cout << " | b[" << i << "] = " << rtA.rt_rhs[i] << " | x[" << i << "] = " << rtA.rt_x[i] << std::endl;
        }
    }

    /* ================================================================== */
    /*  KokkosKernels PCG solver                                           */
    /* ================================================================== */

    /**
     * @brief Solve the linear system using KokkosKernels' Preconditioned
     *        Conjugate Gradient (PCG) with Gauss-Seidel preconditioning.
     *
     * This method:
     *   1. Extracts the diagonal of A (used as the preconditioner basis).
     *   2. Deep-copies the CRS arrays (indices, values, RHS, diagonal)
     *      into KokkosKernels-compatible views.
     *   3. Constructs a KokkosSparse::CrsMatrix object.
     *   4. Invokes KokkosKernels::pcgsolve() with Gauss-Seidel preconditioning.
     *   5. Copies the solution back to the RTMatrix object.
     *
     * The Gauss-Seidel preconditioner is configured via a KernelHandle.
     *
     * @param rtA  The sparse matrix containing the system to solve.
     *             On entry: rt_x holds the initial guess, rt_rhs holds b,
     *                       rt_val / rt_ind / rt_ptr hold the CRS matrix.
     *             On exit:  rt_x holds the computed solution,
     *                       rt_cg_iter holds the iteration count.
     */
    void kkpcg(RTMatrix &rtA)
    {
        bool usePreconditioner = 1;
        const unsigned cg_iteration_limit = iter_max;
        const double cg_iteration_tolerance = eps_min;

        // Extract the diagonal of A for the preconditioner
        decompose(rtA);

        // Copy CRS data into KokkosKernels views
        Kokkos::deep_copy(rt_ind, rtA.rt_ind);
        Kokkos::deep_copy(rt_val, rtA.rt_val);
        Kokkos::deep_copy(rt_rhs, rtA.rt_rhs);
        Kokkos::deep_copy(rt_diag, rtA.rt_diag);

        // Construct the KokkosSparse CRS matrix from the views
        crsMat_t matA = crsMat_t("matA", rtA.rt_nrow, rtA.rt_nrow, rtA.rt_nnz, rt_val, rt_ptr, rt_ind);

        // Initialize KokkosKernels solver handle and result struct
        KokkosKernels::Experimental::Example::CGSolveResult cg_result;
        typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, execution_space, execution_space, execution_space> KernelHandle;
        KernelHandle kh;
        kh.create_gs_handle();

        // Run the Preconditioned CG solver
        KokkosKernels::Experimental::Example::pcgsolve(kh, matA, rt_rhs, rt_vecx, rt_diag, cg_iteration_limit, cg_iteration_tolerance, &cg_result, usePreconditioner);
        Kokkos::fence();

        // Store the iteration count and clean up the handle
        rtA.rt_cg_iter = cg_result.iteration;
        kh.destroy_gs_handle();

        // Copy solution back to the RTMatrix object
        Kokkos::deep_copy(rtA.rt_x, rt_vecx);
    }

    /* ================================================================== */
    /*  Manual CG solver with Jacobi preconditioning                       */
    /* ================================================================== */

    /**
     * @brief Solve the linear system using a hand-coded Preconditioned
     *        Conjugate Gradient (PCG) algorithm with Jacobi preconditioning.
     *
     * Algorithm outline:
     *   1. Extract the diagonal of A and compute the initial residual
     *      r = b - A * x_0.
     *   2. Apply the Jacobi preconditioner: z = diag(A)^{-1} * r.
     *   3. Compute rho = r^T * z.
     *   4. Update the search direction: p = z + beta * p_old.
     *   5. Compute the matrix-vector product: q = A * p.
     *   6. Compute the step length: alpha = rho / (p^T * q).
     *   7. Update the solution: x = x + alpha * p.
     *   8. Update the residual: r = r - alpha * q.
     *   9. Check convergence: ||r|| < eps_min.
     *
     * @param rtA   The sparse matrix containing the system to solve.
     * @param gdom  The groundwater domain (unused but kept for interface
     *              compatibility).
     * @return      Number of CG iterations performed.
     */
    int cg(RTMatrix &rtA, GwDomain &gdom)
    {
        real rho, rhoOld, alpha, beta;
        rhoOld = 0.0;

        // Extract diagonal entries of A for Jacobi preconditioning
        decompose(rtA);

        // Compute the initial residual r = b - A * x_0
        get_residual(rtA);
        iter = 0;
        eps = 1.0;

        while (iter<iter_max & eps> eps_min)
        {
            // Apply Jacobi preconditioner: z(i) = r(i) * diag_inv(i)
            precJACO(rtA);

            // rho = r^T * z  (preconditioned inner product)
            rho = dot(rtA.rt_r, rtA.rt_z, rtA.rt_nrow);

            if (iter == 0)
            {
                // First iteration: set p = z (no previous direction)
                beta = 0.0;
                Kokkos::deep_copy(rtA.rt_p, rtA.rt_z);
            }
            else
            {
                // Subsequent iterations: p = z + beta * p_old
                beta = rho / rhoOld;
                kxpy(rtA.rt_p, beta, rtA.rt_p, rtA.rt_z, rtA.rt_nrow);
            }

            // Matrix-vector product: q = A * p
            mul_MV(rtA.rt_q, rtA, rtA.rt_p);

            // Step length: alpha = rho / (p^T * q)
            alpha = rho / dot(rtA.rt_q, rtA.rt_p, rtA.rt_nrow);

            // Update solution: x = x + alpha * p
            update_X(rtA, alpha);

            // Update residual: r = r - alpha * q
            kxpy(rtA.rt_r, -alpha, rtA.rt_q, rtA.rt_r, rtA.rt_nrow);

            rhoOld = rho;

            // Compute residual norm for convergence check
            eps = pow(dot(rtA.rt_r, rtA.rt_r, rtA.rt_nrow), 0.5);
            iter += 1;
        }
        return iter;
    }

    /* ================================================================== */
    /*  Low-level linear algebra kernels                                   */
    /* ================================================================== */

    /**
     * @brief Apply the Jacobi preconditioner:  z(i) = r(i) / A(i,i).
     *
     * Uses the pre-computed inverse diagonal stored in rt_diag.
     *
     * @param rtA  Sparse matrix with rt_diag containing 1/A(i,i).
     */
    void precJACO(RTMatrix rtA)
    {
        Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idom) { rtA.rt_z(idom) = rtA.rt_r(idom) * rtA.rt_diag(idom); });
    }

    /**
     * @brief Reset a vector to zero.
     *
     * @param x  The vector to reset.
     * @param n  Length of the vector.
     */
    void reset(realArr x, int n)
    {
        Kokkos::parallel_for(n, KOKKOS_LAMBDA(int idom) { x(idom) = 0.0; });
    }

    /**
     * @brief Compute the residual:  r = b - A * x.
     *
     * @param rtA  Sparse matrix; uses rt_ptr, rt_ind, rt_val, rt_rhs, rt_x.
     */
    void get_residual(RTMatrix &rtA)
    {
        Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idom) {
            int icol;
            rtA.rt_r(idom) = rtA.rt_rhs(idom);
            for (icol = rtA.rt_ptr(idom); icol < rtA.rt_ptr(idom+1); icol++) {
                rtA.rt_r(idom) -= rtA.rt_val(icol) * rtA.rt_x(rtA.rt_ind(icol));
            } });
    }

    /**
     * @brief Compute the dot product of two vectors.
     *
     * @param v1  First input vector  [n]
     * @param v2  Second input vector [n]
     * @param n   Vector length
     * @return    The scalar dot product  sum_i(v1[i] * v2[i])
     */
    real dot(realArr v1, realArr v2, int n)
    {
        real out;
        Kokkos::parallel_reduce(n, KOKKOS_LAMBDA(int idx, real &prod) { prod += v1(idx) * v2(idx); }, out);
        return out;
    }

    /**
     * @brief Scaled vector addition:  out = k * x + y.
     *
     * Used for updating search directions and residuals in the CG algorithm.
     *
     * @param out  Output vector [n]
     * @param k    Scalar multiplier applied to x
     * @param x    Input vector to be scaled  [n]
     * @param y    Input vector added as-is    [n]
     * @param n    Vector length
     */
    void kxpy(realArr out, real k, realArr x, realArr y, int n)
    {
        Kokkos::parallel_for(n, KOKKOS_LAMBDA(int idx) { out(idx) = k * x(idx) + y(idx); });
    }

    /**
     * @brief Sparse matrix-vector multiplication:  out = A * x.
     *
     * Iterates over each row of the CRS matrix and computes the dot product
     * of the row with the input vector.
     *
     * @param out  Output vector [nrow]
     * @param rtA  Sparse matrix in CRS format
     * @param x    Input vector  [nrow]
     */
    void mul_MV(realArr out, RTMatrix rtA, realArr x)
    {
        Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
            int icol;
            out(idx) = 0.0;
            for (icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx+1); icol++) {
                out(idx) += rtA.rt_val(icol) * x(rtA.rt_ind(icol));
            } });
    }

    /**
     * @brief Update the solution vector:  x = x + alpha * p.
     *
     * @param rtA    Sparse matrix whose solution vector rt_x is updated.
     * @param alpha  Step length scalar.
     */
    void update_X(RTMatrix rtA, real alpha)
    {
        Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idx) { rtA.rt_x(idx) += alpha * rtA.rt_p(idx); });
    }

    /**
     * @brief Extract and invert the diagonal of the sparse matrix for
     *        Jacobi preconditioning.
     *
     * For each row i, scans the non-zero entries to find A(i,i) and stores
     * 1.0 / A(i,i) in rt_diag(i).  This pre-computation allows the Jacobi
     * preconditioner to be applied as a simple element-wise multiplication.
     *
     * @param rtA  Sparse matrix whose diagonal is extracted.
     */
    void decompose(RTMatrix rtA)
    {
        Kokkos::parallel_for(rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
            int icol;
            rtA.rt_diag(idx) = 0.0;
            for (icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx+1); icol++) {
                if (rtA.rt_ind(icol) == idx) {rtA.rt_diag(idx) = 1.0 / rtA.rt_val(icol);}
            } });
    }
};
#endif
