/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file RTSolverGW.h
 * @brief Linear solvers for the groundwater reactive-transport sparse system.
 *
 * This header provides the RTSolverGW template class, which implements
 * several iterative solvers for the linear system A * x = b arising from
 * the implicit discretization of the groundwater transport equation:
 *
 *   - CG        : Conjugate Gradient with Jacobi preconditioning
 *   - kkpcg     : KokkosKernels-backed Preconditioned CG (offloaded to
 *                 KokkosKernels library routines)
 *   - Gauss-Seidel : Serial Gauss-Seidel smoother (for debugging / small
 *                    problems)
 *   - GMRES     : Generalised Minimum Residual with Jacobi preconditioning
 *                 and restarted Arnoldi (m_restart = 30)
 *
 * The class operates on an RTMatrix object (CRS-format sparse matrix) and
 * a GwDomain object (grid geometry).  All Kokkos-parallel operations are
 * parameterised by the execution_space template argument, enabling
 * GPU offloading when Kokkos is configured for CUDA / HIP / SYCL.
 *
 * Physical units:
 *   - Matrix A entries : transport coefficients  [1/s]
 *   - RHS b            : mass loading rate       [mg/s]
 *   - Solution x       : solute concentration    [mg/L]
 *
 * @see RTMatrix.h       Sparse matrix storage (CRS format)
 * @see RTSolver-kkpcg.h  KokkosKernels PCG solver wrapper
 */
#ifndef _RT_SOLVERGW_H_
#define _RT_SOLVERGW_H_

#include "KokkosKernels_config.h"
#include "KokkosSparse_pcg.hpp"
#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_default_types.hpp"
#include <iostream>
#include "GwDomain.h"
#include "RTMatrix.h"
#include <vector>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosKernels_Handle.hpp>
#include <KokkosSparse_spiluk.hpp>
#include <KokkosSparse_sptrsv.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas.hpp>
#include <KokkosSparse_gauss_seidel.hpp>
#include <KokkosBlas1_dot.hpp>
#include <cmath>
#include <memory>
#include <fstream>

// Note: File-level using-declarations are intentionally NOT placed here to
// avoid confusion with the identically-named typedefs inside the class.
// The KokkosKernels namespaces are imported below for convenience.

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosKernels::Experimental::Example;

/**
 * @brief Force LayoutRight for 2D Kokkos views used as Krylov basis matrices
 *        in GMRES.  With LayoutRight, subview(v, k, ALL) yields a contiguous
 *        1D view, which is both cache-friendly on CPUs and efficient for
 *        GPU memory access patterns.
 */
using realArr2_LR = Kokkos::View<real **, Kokkos::LayoutRight,
                                 Kokkos::Device<Kokkos::DefaultExecutionSpace,
                                                Kokkos::SharedSpace>>;

/**
 * @class RTSolverGW
 * @brief Template class providing iterative linear solvers for the
 *        groundwater reactive-transport system.
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
    //  These define the view types for row pointers, column indices,
    //  and matrix values used by the KokkosKernels solver backends.
    // ==================================================================
    typedef typename KokkosSparse::CrsMatrix<default_scalar, default_lno_t,
                                             execution_space, void, default_size_type>
        crsMat_t;
    typedef typename crsMat_t::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
    typedef typename scalar_view_t::value_type scalar_t;

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
     * @brief Scan the RHS and solution vectors for abnormal numerical values.
     *
     * Checks each entry of rt_rhs and rt_x for NaN, infinity, or absolute
     * values exceeding a threshold (1e10).  Prints diagnostic information
     * for any abnormal values found.
     *
     * @param rtA  The sparse matrix whose RHS and solution vectors are checked.
     */
    void printRTMatrix(const RTMatrix &rtA)
    {
        const double THRESHOLD = 1e10;
        int n = rtA.rt_nrow;
        bool foundAbnormal = false;
        for (int i = 0; i < n; ++i)
        {
            // Check RHS vector for NaN, Inf, or excessively large values
            if (std::isnan(rtA.rt_rhs(i)) || std::isinf(rtA.rt_rhs(i)) || std::abs(rtA.rt_rhs(i)) > THRESHOLD)
            {
                std::cout << "Abnormal value - b[" << i << "] = " << rtA.rt_rhs(i) << std::endl;
                foundAbnormal = true;
            }
            // Check solution vector for NaN, Inf, or excessively large values
            if (std::isnan(rtA.rt_x(i)) || std::isinf(rtA.rt_x(i)) || std::abs(rtA.rt_x(i)) > THRESHOLD)
            {
                std::cout << "Abnormal value - x[" << i << "] = " << rtA.rt_x(i) << std::endl;
                foundAbnormal = true;
            }
        }
        if (!foundAbnormal)
            std::cout << "No abnormal numerical values detected" << std::endl;
    }

    /* ================================================================== */
    /*  Gauss-Seidel solver                                                */
    /* ================================================================== */

    /**
     * @brief Solve the linear system using a serial Gauss-Seidel iteration.
     *
     * This solver runs entirely on the host (CPU) using std::vector.  It is
     * intended for debugging and for small-scale problems where the overhead
     * of Kokkos parallelisation is not justified.
     *
     * The method iterates until either:
     *   - The L2-norm of the update (x - x_old) falls below eps_min, or
     *   - The maximum iteration count iter_max is reached.
     *
     * If any NaN values appear in the solution, printRTMatrix() is called
     * automatically for diagnostics.
     *
     * @param rtA  The sparse matrix containing the system to solve.
     *             On entry: rt_x holds the initial guess, rt_rhs holds b.
     *             On exit:  rt_x holds the computed solution.
     */
    void Gauss_Seidel(RTMatrix &rtA)
    {
        int n = rtA.rt_nrow;
        std::vector<double> x(n), b(n);
        // Copy Kokkos views to host std::vectors for serial access
        for (int i = 0; i < n; ++i)
        {
            b[i] = rtA.rt_rhs(i);
            x[i] = rtA.rt_x(i);
        }

        for (int it = 0; it < iter_max; ++it)
        {
            std::vector<double> x_old = x;
            for (int i = 0; i < n; ++i)
            {
                double sum = 0.0, diag = 0.0;
                for (int j = rtA.rt_ptr(i); j < rtA.rt_ptr(i + 1); ++j)
                {
                    int col = rtA.rt_ind(j);
                    if (col == i)
                        diag = rtA.rt_val(j);       // Diagonal entry
                    else
                        sum += rtA.rt_val(j) * x[col];  // Off-diagonal contribution
                }
                if (diag != 0)
                    x[i] = (b[i] - sum) / diag;
            }
            // Compute L2-norm of the update
            double error = 0.0;
            for (int i = 0; i < n; ++i)
                error += std::pow(x[i] - x_old[i], 2);
            error = std::sqrt(error);
            if (error < eps_min)
                break;
        }
        // Copy solution back to Kokkos view; diagnose NaN if present
        for (int i = 0; i < n; ++i)
        {
            rtA.rt_x(i) = x[i];
            if (std::isnan(rtA.rt_x(i)))
                printRTMatrix(rtA);
        }
    }

    /* ================================================================== */
    /*  Low-level linear algebra kernels (template-based for flexibility)  */
    /* ================================================================== */

    /**
     * @brief Compute the dot product of two vectors.
     *
     * Templated to accept views of any Kokkos Layout (LayoutLeft,
     * LayoutRight, LayoutStride), which is necessary for GMRES where
     * subviews of the Krylov basis matrix may have LayoutStride.
     *
     * @param v1  First input vector  [n]
     * @param v2  Second input vector [n]
     * @param n   Vector length
     * @return    The scalar dot product  sum_i(v1[i] * v2[i])
     */
    template <typename V1, typename V2>
    real dot(const V1 &v1, const V2 &v2, int n)
    {
        real out = 0.0;
        Kokkos::parallel_reduce("rt_dot", n, KOKKOS_LAMBDA(int idx, real &prod) { prod += v1(idx) * v2(idx); }, out);
        return out;
    }

    /**
     * @brief Scaled vector addition:  out = k * x + y.
     *
     * Used in CG for updating the search direction (beta * p + z) and
     * residual update (-alpha * q + r).  Templated for layout flexibility.
     *
     * @param out  Output vector [n]
     * @param k    Scalar multiplier applied to x
     * @param x    Input vector to be scaled  [n]
     * @param y    Input vector added as-is    [n]
     * @param n    Vector length
     */
    template <typename VOut, typename VX, typename VY>
    void kxpy(const VOut &out, real k, const VX &x, const VY &y, int n)
    {
        Kokkos::parallel_for("rt_kxpy", n, KOKKOS_LAMBDA(int idx) { out(idx) = k * x(idx) + y(idx); });
    }

    /**
     * @brief Scaled vector plus scalar:  out = k * x + y.
     *
     * Similar to kxpy but with a scalar y instead of a vector.  Used in
     * GMRES for normalising the Krylov basis vectors.
     *
     * @param out  Output vector [n]
     * @param k    Scalar multiplier applied to x
     * @param x    Input vector to be scaled  [n]
     * @param y    Scalar added to each element
     * @param n    Vector length
     */
    template <typename VOut, typename VX>
    void kxpy2(const VOut &out, real k, const VX &x, double y, int n)
    {
        Kokkos::parallel_for("rt_kxpy2", n, KOKKOS_LAMBDA(int idx) { out(idx) = k * x(idx) + y; });
    }

    /**
     * @brief Sparse matrix-vector multiplication:  out = A * x.
     *
     * Iterates over each row of the CRS matrix and computes the dot product
     * of the row with the input vector.  Templated to accept views of any
     * Kokkos Layout.
     *
     * @param out  Output vector [nrow]
     * @param rtA  Sparse matrix in CRS format
     * @param x    Input vector  [nrow]
     */
    template <typename VOut, typename VX>
    void mul_MV(const VOut &out, RTMatrix &rtA, const VX &x)
    {
        Kokkos::parallel_for("rt_mul_MV", rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
            real acc = 0.0;
            for (int icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx + 1); ++icol) {
                acc += rtA.rt_val(icol) * x(rtA.rt_ind(icol));
            }
            out(idx) = acc; });
    }

    /**
     * @brief Compute the residual:  r = b - A * x.
     *
     * The residual is stored in rtA.rt_r and is used to monitor convergence
     * and as input to preconditioners.
     *
     * @param rtA  Sparse matrix; uses rt_ptr, rt_ind, rt_val, rt_rhs, rt_x.
     */
    void get_residual(RTMatrix &rtA)
    {
        Kokkos::parallel_for("rt_residual", rtA.rt_nrow, KOKKOS_LAMBDA(int idx) {
            rtA.rt_r(idx) = rtA.rt_rhs(idx);
            for (int icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx + 1); ++icol) {
                rtA.rt_r(idx) -= rtA.rt_val(icol) * rtA.rt_x(rtA.rt_ind(icol));
            } });
    }

    /**
     * @brief Apply the Jacobi preconditioner:  z = M^{-1} * r.
     *
     * The Jacobi preconditioner uses the inverse of the diagonal entries:
     *   z(i) = r(i) / A(i,i)
     * If the diagonal entry is zero, z(i) is set to zero.
     *
     * @param rtA  Sparse matrix (used to find diagonal entries)
     * @param r    Input residual vector
     * @param z    Output preconditioned residual vector
     */
    template <typename VR, typename VZ>
    void Jacobi_Precondition(RTMatrix &rtA, const VR &r, const VZ &z)
    {
        int n = rtA.rt_nrow;
        Kokkos::parallel_for("Jacobi_Precondition", n, KOKKOS_LAMBDA(int i) {
            double diag_val = 0.0;
            for (int j = rtA.rt_ptr(i); j < rtA.rt_ptr(i + 1); ++j) {
                if (rtA.rt_ind(j) == i) { diag_val = rtA.rt_val(j); break; }
            }
            z(i) = (diag_val != 0.0) ? r(i) / diag_val : 0.0; });
    }

    /* ================================================================== */
    /*  GMRES solver with Jacobi preconditioning                           */
    /* ================================================================== */

    /**
     * @brief Solve the linear system using restarted GMRES with Jacobi
     *        preconditioning (GMRES(m)).
     *
     * The Generalised Minimum Residual method is suitable for non-symmetric
     * linear systems.  This implementation:
     *   1. Computes the initial preconditioned residual.
     *   2. Builds a Krylov subspace of dimension m_restart using the
     *      Modified Gram-Schmidt (MGS) orthogonalisation.
     *   3. Solves the resulting upper Hessenberg least-squares problem via
     *      Givens rotations.
     *   4. Restarts if convergence is not achieved within m_restart steps.
     *
     * Convergence criterion:
     *   |g_{k+1}| < eps_min_local  (1e-6 by default)
     *
     * @param rtA  The sparse matrix containing the system to solve.
     *             On entry: rt_x holds the initial guess, rt_rhs holds b.
     *             On exit:  rt_x holds the computed solution.
     * @return     Total number of GMRES iterations performed.
     */
    int Jacobi_GMRES_Solve(RTMatrix &rtA)
    {

        int m_restart = 30;            ///< Krylov subspace dimension before restart
        const int n = rtA.rt_nrow;
        const real eps_min_local = 1e-6;  ///< Convergence tolerance for GMRES
        int iter_total = 0;
        real eps_local = 1.0;

        // References to work vectors in the matrix object
        realArr &rt_z = rtA.rt_z;
        realArr &rt_w = rtA.rt_w;

        // Krylov basis matrix with LayoutRight for contiguous subview access
        realArr2_LR v("v_LR", m_restart + 1, n);

        // Givens rotation coefficients and projected RHS
        realArr c_arr("c", m_restart), s_arr("s", m_restart);
        realArr g_arr("g", m_restart + 1);

        // Upper Hessenberg matrix (projected system)
        realArr2 H("H", m_restart + 1, m_restart);

        // Main GMRES iteration loop (with restarts)
        while (iter_total < iter_max && eps_local > eps_min_local)
        {

            // Compute initial preconditioned residual: z = M^{-1}(b - A*x)
            get_residual(rtA);
            Jacobi_Precondition(rtA, rtA.rt_r, rt_z);

            // Compute the norm of the preconditioned residual (beta)
            real beta = std::sqrt(dot(rt_z, rt_z, n));
            if (beta < eps_min_local)
                break;

            // Normalise z to obtain the first Krylov basis vector v[0]
            Kokkos::parallel_for("rt_norm_z", n, KOKKOS_LAMBDA(const int i) { rt_z(i) /= beta; });

            // v[0] = z / ||z||
            auto v0 = Kokkos::subview(v, 0, Kokkos::ALL());
            Kokkos::deep_copy(v0, rt_z);

            // Initialise the projected RHS: g = ||beta|| * e_1
            g_arr(0) = beta;
            for (int i = 1; i <= m_restart; ++i)
                g_arr(i) = 0.0;

            // Arnoldi iteration: build the Krylov basis
            int k;
            for (k = 0; k < m_restart; ++k)
            {
                auto vk = Kokkos::subview(v, k, Kokkos::ALL());

                // Compute w = M^{-1} * A * v[k]  (preconditioned matrix-vector product)
                mul_MV(rt_w, rtA, vk);
                Jacobi_Precondition(rtA, rt_w, rt_w);

                // Modified Gram-Schmidt orthogonalisation against all previous
                // basis vectors v[0], ..., v[k]
                for (int i = 0; i <= k; ++i)
                {
                    auto vi = Kokkos::subview(v, i, Kokkos::ALL());
                    H(i, k) = dot(vi, rt_w, n);
                    kxpy(rt_w, -H(i, k), vi, rt_w, n);
                }
                // Norm of the orthogonalised vector
                H(k + 1, k) = std::sqrt(dot(rt_w, rt_w, n));

                // Normalise to obtain the next basis vector v[k+1]
                if (H(k + 1, k) != 0.0)
                {
                    auto vkp1 = Kokkos::subview(v, k + 1, Kokkos::ALL());
                    kxpy2(vkp1, 1.0 / H(k + 1, k), rt_w, 0.0, n);
                }

                // Apply all previous Givens rotations to the new column of H
                for (int i = 0; i < k; ++i)
                {
                    real temp = c_arr(i) * H(i, k) + s_arr(i) * H(i + 1, k);
                    H(i + 1, k) = -s_arr(i) * H(i, k) + c_arr(i) * H(i + 1, k);
                    H(i, k) = temp;
                }

                // Compute new Givens rotation to eliminate H(k+1, k)
                real gamma = std::sqrt(H(k, k) * H(k, k) + H(k + 1, k) * H(k + 1, k));
                if (gamma == 0)
                    gamma = 1e-16;
                c_arr(k) = H(k, k) / gamma;     // Cosine of rotation angle
                s_arr(k) = H(k + 1, k) / gamma;  // Sine of rotation angle
                H(k, k) = gamma;

                // Apply the new rotation to the projected RHS g
                real gk = c_arr(k) * g_arr(k) + s_arr(k) * g_arr(k + 1);
                g_arr(k + 1) = -s_arr(k) * g_arr(k) + c_arr(k) * g_arr(k + 1);
                g_arr(k) = gk;

                // Check convergence: |g(k+1)| is the residual norm
                eps_local = std::fabs(g_arr(k + 1));
                if (eps_local < eps_min_local)
                {
                    k++;
                    break;
                }
            }

            // Back-substitution to solve the upper-triangular system
            // H(0:k, 0:k) * y = g(0:k)
            realArr y("y", k);
            for (int i = k - 1; i >= 0; --i)
            {
                y(i) = g_arr(i);
                for (int j = i + 1; j < k; ++j)
                    y(i) -= H(i, j) * y(j);
                y(i) /= H(i, i);
            }

            // Update the solution: x = x + V(:,0:k) * y
            Kokkos::parallel_for("rt_update_x", n, KOKKOS_LAMBDA(int idx) {
                real delta = 0.0;
                for (int j = 0; j < k; ++j) delta += v(j, idx) * y(j);
                rtA.rt_x(idx) += delta; });

            iter_total += k;
            eps_local = std::fabs(g_arr(k));
        }
        return iter_total;
    }
};

#endif
