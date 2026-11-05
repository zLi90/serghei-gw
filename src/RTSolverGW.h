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
#include <cstdlib>
#include <memory>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

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

#if SERGHEI_RTGW_SPMV_KOKKOS_MODE == 1
    // CrsMatrix wrapper typed to alias RTMatrix's own SharedSpace
    // views directly (verified: this constructor accepts rtA.rt_val /
    // rt_ptr / rt_ind and the values view aliases rtA.rt_val, not a copy),
    // so the wrapper tracks every re-assembly with zero data movement.
    typedef Kokkos::Device<execution_space, Kokkos::SharedSpace> rt_spmv_device_t;
    typedef KokkosSparse::CrsMatrix<real, int, rt_spmv_device_t, void, int> rtSpmvMat_t;
#endif

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
    int gmres_m_restart; ///< GMRES restart dimension

    /* ------------------------------------------------------------------ */
    /*  KokkosKernels views for interfacing with the library PCG solver    */
    /* ------------------------------------------------------------------ */
    lno_view_t rt_ptr;       ///< Row-pointer view (KokkosKernels format)
    lno_nnz_view_t rt_ind;   ///< Column-index view (KokkosKernels format)
    scalar_view_t rt_val;    ///< Matrix values view (KokkosKernels format)
    scalar_view_t rt_rhs;    ///< Right-hand side view (KokkosKernels format)
    scalar_view_t rt_vecx;   ///< Solution vector view (KokkosKernels format)
    scalar_view_t rt_diag;   ///< Diagonal values view (KokkosKernels format)
    realArr2_LR gmres_v;     ///< Cached Krylov basis for GMRES
    realArr2_LR gmres_H;     ///< Cached Hessenberg matrix for GMRES
    realArr gmres_c;         ///< Cached Givens cosine coefficients
    realArr gmres_s;         ///< Cached Givens sine coefficients
    realArr gmres_g;         ///< Cached projected RHS
    realArr gmres_y;         ///< Cached least-squares solution
#if SERGHEI_RTGW_SPMV_KOKKOS_MODE == 1
    rtSpmvMat_t rt_crs_A;    ///<  Aliased CRS wrapper for KokkosSparse::spmv
#endif

    /* ------------------------------------------------------------------ */
    /*   Upgradeable preconditioner state (SERGHEI_RTGW_SOLVER_PREC)         */
    /* ------------------------------------------------------------------ */
#if SERGHEI_RTGW_SOLVER_PREC != 0
    intArr prec_diag_j;       ///< [nrow] nnz position of each row's diagonal entry in the frozen graph
    bool prec_ready = false;  ///< one-time structural setup succeeded
    bool prec_active = false; ///< current factorisation is healthy (else Jacobi fallback)
#endif
#if SERGHEI_RTGW_SOLVER_PREC == 1
    realArr prec_lu;          ///< [nnz] ILU(0) factor values stored in A's own pattern
    realArr prec_uinv;        ///< [nrow] inverse U diagonal per row
    intArr prec_iw;           ///< [nrow] scatter work array: column -> position in the row being factorised
#endif
#if SERGHEI_RTGW_SOLVER_PREC == 2
    realArr prec_bj_l;        ///< [nrow] Thomas elimination factor of the z-column blocks
    realArr prec_bj_u;        ///< [nrow] Thomas modified pivot of the z-column blocks
    intArr prec_sub_j;        ///< [nrow] nnz position of the z-neighbour below (col == i - plane), -1 if absent
    intArr prec_sup_j;        ///< [nrow] nnz position of the z-neighbour above (col == i + plane), -1 if absent
    int prec_plane = 0;       ///< nx * ny (row stride between consecutive z layers)
    int prec_nz = 0;          ///< number of z layers (block size)
    int prec_ncol = 0;        ///< number of vertical column blocks (= nx * ny)
#endif

    /* ------------------------------------------------------------------ */
    /*  [P3-4] Linear-solver statistics & phase profiling                  */
    /*  (counters always on; phase timers only with SERGHEI_RTGW_SOLVER_PROFILE=1)  */
    /* ------------------------------------------------------------------ */
    /// Per-species accumulation bucket for the solver counters and, when
    /// profiling is enabled (SERGHEI_RTGW_SOLVER_PROFILE macro or env), the phase timers.
    struct LSvSpecStats
    {
        long n_solves = 0;         ///< linear solves performed
        long n_iters = 0;          ///< total GMRES (Arnoldi) iterations
        long instant_exits = 0;    ///< solves converging at the initial residual check (0 iterations)
        long restarts = 0;         ///< restart cycles entered
        long nonconv = 0;          ///< solves leaving the loop on iter_max
        double t_total = 0.0;      ///< whole-solve time (cross-check vs rtgw.linsol)
        double eps_min = 1e300;    ///< smallest final preconditioned residual
        double eps_max = 0.0;      ///< largest final preconditioned residual
        double t_precond = 0.0;    ///< inverse-diagonal update + preconditioner factorisation
        double t_init_res = 0.0;   ///< restart-start residual + preconditioner apply + beta dot
        double t_spmv_apply = 0.0; ///< Arnoldi SpMV + preconditioner apply
        double t_ortho = 0.0;      ///< MGS dot/axpy chain + normalisation
        double t_givens = 0.0;     ///< Givens rotations + convergence check (host scalars)
        double t_backsub = 0.0;    ///< back-substitution + solution update kernel
        std::vector<int> iters;    ///< per-solve iteration counts (median / p90)
        std::vector<double> eps;   ///< per-solve final residuals (median)
    };
    std::vector<LSvSpecStats> lsv;   ///< per-species buckets, grown on demand
    long lsv_half_solves[2] = {0, 0}; ///< solves per transport half-step (Step1 / Step3)
    long lsv_half_iters[2] = {0, 0};
    double lsv_half_time[2] = {0.0, 0.0};
    long lsv_prec_fallbacks = 0;     ///< factorisations degraded to Jacobi (SERGHEI_RTGW_SOLVER_PREC 1/2)
    bool lsv_profile = false;        ///< phase timers active (SERGHEI_RTGW_SOLVER_PROFILE macro or env)
    long lsv_nrow = 0;               ///< matrix dimension (informational)
    long lsv_nnz = 0;                ///< matrix non-zeros (informational)

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
        iter = 0;
        eps = 1.0;
        gmres_m_restart = 30;
        gmres_v = realArr2_LR("rt_gmres_v", gmres_m_restart + 1, nrow);
        gmres_H = realArr2_LR("rt_gmres_H", gmres_m_restart + 1, gmres_m_restart);
        gmres_c = realArr("rt_gmres_c", gmres_m_restart);
        gmres_s = realArr("rt_gmres_s", gmres_m_restart);
        gmres_g = realArr("rt_gmres_g", gmres_m_restart + 1);
        gmres_y = realArr("rt_gmres_y", gmres_m_restart);
#if SERGHEI_KOKKOSKERNELS_SOLVER
        rt_ptr = lno_view_t("rt_ptr", rtA.rt_nrow + 1);      // Starting index of non-zeros for each row
        rt_ind = lno_nnz_view_t("rt_ind", rtA.rt_nnz);        // Column indices of non-zero entries
        rt_val = scalar_view_t("rt_val", rtA.rt_nnz);          // Values of non-zero entries
        rt_rhs = scalar_view_t("rt_rhs", rtA.rt_nrow);         // Right-hand side vector b
        rt_vecx = scalar_view_t("rt_vecx", rtA.rt_nrow);       // Approximate solution vector x
        rt_diag = scalar_view_t("rt_diag", rtA.rt_nrow);       // Diagonal entries (for preconditioning)
        Kokkos::deep_copy(rt_ptr, rtA.rt_ptr);
#endif
#if SERGHEI_RTGW_SPMV_KOKKOS_MODE == 1
        // Wrap rtA's graph and values once.  The values view aliases
        // rtA.rt_val, so every re-assembly is visible without any copy; the
        // graph (rt_ptr / rt_ind) is frozen by RTMatrix::init / fillGraph
        // at start-up.
        static_assert(std::is_same<execution_space, Kokkos::DefaultExecutionSpace>::value,
                      "SERGHEI_RTGW_SPMV_KOKKOS_MODE=1 requires the default execution space "
                      "(RTMatrix views live in DefaultExecutionSpace/SharedSpace)");
        rt_crs_A = rtSpmvMat_t("rt_crs_A", rtA.rt_nrow, rtA.rt_nrow,
                               (typename rtSpmvMat_t::size_type)rtA.rt_nnz,
                               rtA.rt_val, rtA.rt_ptr, rtA.rt_ind);
#endif
#if SERGHEI_RTGW_SOLVER_PREC != 0
        // One-time preconditioner structural setup.  Safe here:
        // serghei.h runs rtA.fillGraph() (freezing the CRS graph) before
        // this init, so rt_ptr / rt_ind hold their final values.
        RT_Preconditioner_Setup(rtA, gdom);
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
            z(i) = r(i) * rtA.rt_diag(i); });
    }

    /**
     * @brief Cache inverse diagonal values for Jacobi preconditioning.
     *
     * The matrix is reassembled before each transport solve.  GMRES applies
     * the preconditioner many times per solve, so scanning each CRS row once
     * here is substantially cheaper than scanning it inside every
     * preconditioner call.
     */
    void Update_Jacobi_Inverse_Diagonal(RTMatrix &rtA)
    {
        int n = rtA.rt_nrow;
        Kokkos::parallel_for("RT_Update_Jacobi_InvDiag", n, KOKKOS_LAMBDA(int i) {
            real diag_val = 0.0;
            for (int j = rtA.rt_ptr(i); j < rtA.rt_ptr(i + 1); ++j) {
                if (rtA.rt_ind(j) == i) { diag_val = rtA.rt_val(j); break; }
            }
            rtA.rt_diag(i) = (diag_val != 0.0) ? 1.0 / diag_val : 0.0; });
    }

    /* ================================================================== */
    /*   Upgradeable preconditioner (SERGHEI_RTGW_SOLVER_PREC)                */
    /*                                                                     */
    /*  0 = Jacobi (default, code paths above unchanged),                  */
    /*  1 = ILU(0) on the frozen 19-point CRS pattern,                     */
    /*  2 = block-Jacobi with vertical (z) soil-column blocks.             */
    /*                                                                     */
    /*  The factorisation kernels are self-contained rather than calling   */
    /*  KokkosKernels spiluk / sptrsv, for two reasons read from the       */
    /*  installed KokkosKernels 4.2.00 sources:                            */
    /*   (a) spiluk's level-scheduled numeric phase launches one kernel    */
    /*       per level, and with the natural row ordering of this stencil  */
    /*       (every row i depends on row i-1 through L) the schedule       */
    /*       degenerates to ~nrow single-row levels, i.e. ~nrow kernel     */
    /*       launches per factorisation;                                   */
    /*   (b) sptrsv divides by the stored diagonal and expects an explicit */
    /*       L diagonal, while spiluk emits a unit-lower L without a       */
    /*       stored diagonal and an inverted U diagonal, so the two        */
    /*       library entry points are not drop-in compatible.              */
    /* ================================================================== */

    /**
     * @brief One-time structural setup for the upgraded preconditioner.
     *
     * Locates, in the frozen CRS graph built by RTMatrix::fillGraph(), the
     * nnz position of every row's diagonal entry (all modes) and, for
     * mode 2, the nnz positions of the vertical (z) neighbours.  For
     * mode 1 it also validates the layout required by the in-place IKJ
     * ILU(0): entries strictly ascending by column, i.e. rows laid out as
     * [lower ... | diagonal | ... upper].  fillGraph() inserts entries in
     * exactly this order, so the check is a fail-fast guard against future
     * graph changes, not a per-run cost.
     *
     * On any structural violation prec_ready stays false and every solve
     * silently falls back to Jacobi (with a one-time warning).
     */
    void RT_Preconditioner_Setup(RTMatrix &rtA, GwDomain &gdom)
    {
#if SERGHEI_RTGW_SOLVER_PREC != 0
        const int n = rtA.rt_nrow;
        prec_diag_j = intArr("prec_diag_j", n);

        auto ptr_l = rtA.rt_ptr;
        auto ind_l = rtA.rt_ind;
        auto diagj_l = prec_diag_j;
        Kokkos::parallel_for("RT_Prec_LocateDiag", n, KOKKOS_LAMBDA(int i) {
            int dj = -1;
            for (int k = ptr_l(i); k < ptr_l(i + 1); ++k)
                if (ind_l(k) == i) { dj = k; break; }
            diagj_l(i) = dj; });
        Kokkos::fence();

        // Host-side fail-fast validation on mirrors of the frozen graph.
        bool ok = true;
        {
            auto h_ptr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rtA.rt_ptr);
            auto h_ind = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rtA.rt_ind);
            auto h_dj = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), prec_diag_j);
            for (int i = 0; i < n && ok; ++i)
            {
                if (h_dj(i) < 0) ok = false; // diagonal must exist in every row
#if SERGHEI_RTGW_SOLVER_PREC == 1
                for (int k = h_ptr(i) + 1; ok && k < h_ptr(i + 1); ++k)
                    if (h_ind(k) <= h_ind(k - 1)) ok = false; // must be column-ascending
#endif
            }
        }
        if (!ok)
        {
            std::cerr << "[RTGW_PREC] WARNING: CRS graph violates the structural "
                         "invariant required by SERGHEI_RTGW_SOLVER_PREC=" << SERGHEI_RTGW_SOLVER_PREC
                      << " (missing diagonal or non-ascending row).  Falling back to "
                         "Jacobi preconditioning for the whole run." << std::endl;
            return;
        }

#if SERGHEI_RTGW_SOLVER_PREC == 1
        prec_lu = realArr("prec_lu", rtA.rt_nnz);
        prec_uinv = realArr("prec_uinv", n);
        prec_iw = intArr("prec_iw", n);
        Kokkos::deep_copy(prec_iw, -1);
#endif
#if SERGHEI_RTGW_SOLVER_PREC == 2
        prec_bj_l = realArr("prec_bj_l", n);
        prec_bj_u = realArr("prec_bj_u", n);
        prec_sub_j = intArr("prec_sub_j", n);
        prec_sup_j = intArr("prec_sup_j", n);
        prec_plane = gdom.nx * gdom.ny;
        prec_nz = gdom.nz;
        prec_ncol = prec_plane;
        auto sub_l = prec_sub_j;
        auto sup_l = prec_sup_j;
        const int plane_l = prec_plane;
        Kokkos::parallel_for("RT_Prec_LocateZ", n, KOKKOS_LAMBDA(int i) {
            int sj = -1, pj = -1;
            for (int k = ptr_l(i); k < ptr_l(i + 1); ++k)
            {
                const int c = ind_l(k);
                if (c == i - plane_l) sj = k;
                else if (c == i + plane_l) pj = k;
            }
            sub_l(i) = sj;
            sup_l(i) = pj; });
        Kokkos::fence();
#endif
        prec_ready = true;
#endif // SERGHEI_RTGW_SOLVER_PREC != 0
    }

    /**
     * @brief Factorise the upgraded preconditioner for the current matrix.
     *
     * Called at the top of Jacobi_GMRES_Solve, i.e. once per freshly
     * assembled system (per species per Picard iteration), which matches
     * the matrix-update cadence ILU-type factorisations require.  The
     * cost lands inside the rtgw.linsol timer, so A/B benchmarks measure
     * factorisation + iteration reduction together, as they should.
     *
     * Mode 1 (ILU(0)): factorises in place on a copy of rt_val in A's own
     * pattern.  Row layout [lower | diag | upper] with ascending columns
     * means the k-loop from the row start to the diagonal visits L entries
     * in increasing column order, as the IKJ recurrence requires.  The
     * elimination chain (row i reads final factors of rows j < i) is
     * strictly serial, so this is a single-threaded kernel; on the CPU
     * target it costs a few ms per factorisation.
     *
     * Mode 2 (block-Jacobi): one Thomas factorisation per vertical soil
     * column (rows col, col+plane, ..., col+(nz-1)*plane), nx*ny columns
     * factorised in parallel.
     *
     * A health check (finite values, positive pivots) guards the solve;
     * on failure prec_active is false and every apply of this solve falls
     * back to Jacobi.
     */
    void RT_Preconditioner_Factorize(RTMatrix &rtA)
    {
#if SERGHEI_RTGW_SOLVER_PREC != 0
        prec_active = false;
        if (!prec_ready) { lsv_prec_fallbacks++; return; } // [P3-4] stats: Jacobi fallback
        const int n = rtA.rt_nrow;
        auto diagj_l = prec_diag_j;
#if SERGHEI_RTGW_SOLVER_PREC == 1
        auto lu = prec_lu;
        auto uinv = prec_uinv;
        auto iw = prec_iw;
        auto aval = rtA.rt_val;
        auto jacinv = rtA.rt_diag;
        auto ptr_l = rtA.rt_ptr;
        auto ind_l = rtA.rt_ind;
        Kokkos::deep_copy(lu, aval);
        Kokkos::fence();
        Kokkos::parallel_for("RT_ILU0_Factorize",
            Kokkos::RangePolicy<execution_space>(0, 1),
            KOKKOS_LAMBDA(int) {
                for (int i = 0; i < n; ++i)
                {
                    const int dj = diagj_l(i);
                    for (int k = ptr_l(i); k < ptr_l(i + 1); ++k) iw(ind_l(k)) = k;
                    for (int k = ptr_l(i); k < dj; ++k) // L entries, ascending columns
                    {
                        const int j = ind_l(k);
                        const real fact = lu(k) * uinv(j); // a_ij / u_jj
                        lu(k) = fact;
                        for (int kk = diagj_l(j) + 1; kk < ptr_l(j + 1); ++kk)
                        {
                            const int pos = iw(ind_l(kk)); // final U row of row j
                            if (pos != -1) lu(pos) -= lu(kk) * fact;
                        }
                    }
                    const real d = lu(dj);
                    uinv(i) = (d != 0.0) ? 1.0 / d : jacinv(i); // Jacobi fallback on zero pivot
                    for (int k = ptr_l(i); k < ptr_l(i + 1); ++k) iw(ind_l(k)) = -1;
                } });
        Kokkos::fence();
        double sumabs = 0.0;
        Kokkos::parallel_reduce("RT_ILU0_Check", rtA.rt_nnz,
            KOKKOS_LAMBDA(int k, double &s) { s += fabs(lu(k)); }, sumabs);
        real minu = 1e300;
        Kokkos::parallel_reduce("RT_ILU0_MinPivot", n,
            KOKKOS_LAMBDA(int i, real &m) { if (uinv(i) < m) m = uinv(i); },
            Kokkos::Min<real>(minu));
        prec_active = std::isfinite(sumabs) && std::isfinite((double)minu) && (minu > 0.0);
        if (!prec_active) lsv_prec_fallbacks++; // [P3-4] stats: Jacobi fallback this solve
#elif SERGHEI_RTGW_SOLVER_PREC == 2
        auto l = prec_bj_l;
        auto u = prec_bj_u;
        auto sub_l = prec_sub_j;
        auto sup_l = prec_sup_j;
        auto aval = rtA.rt_val;
        const int plane = prec_plane;
        const int nz = prec_nz;
        const int ncol = prec_ncol;
        Kokkos::parallel_for("RT_BJ_Factorize", ncol, KOKKOS_LAMBDA(int col) {
            {
                const int i = col;
                const real b = aval(diagj_l(i));
                l(i) = 0.0;
                u(i) = (b != 0.0) ? b : 1.0;
            }
            for (int kk = 1; kk < nz; ++kk)
            {
                const int i = col + kk * plane;
                const int ip = i - plane;
                const real a = (sub_l(i) >= 0) ? aval(sub_l(i)) : 0.0;
                const real b = aval(diagj_l(i));
                const real cp = (sup_l(ip) >= 0) ? aval(sup_l(ip)) : 0.0;
                const real up = u(ip);
                if (up != 0.0)
                {
                    l(i) = a / up;
                    u(i) = b - l(i) * cp;
                }
                else
                {
                    l(i) = 0.0; // degrade this row to a plain Jacobi row
                    u(i) = (b != 0.0) ? b : 1.0;
                }
            } });
        Kokkos::fence();
        double sumabs = 0.0;
        Kokkos::parallel_reduce("RT_BJ_Check", n,
            KOKKOS_LAMBDA(int i, double &s) { s += fabs(l(i)) + fabs(u(i)); }, sumabs);
        real minu = 1e300;
        Kokkos::parallel_reduce("RT_BJ_MinPivot", n,
            KOKKOS_LAMBDA(int i, real &m) { if (u(i) < m) m = u(i); },
            Kokkos::Min<real>(minu));
        prec_active = std::isfinite(sumabs) && std::isfinite((double)minu) && (minu > 0.0);
        if (!prec_active) lsv_prec_fallbacks++; // [P3-4] stats: Jacobi fallback this solve
#endif
#endif // SERGHEI_RTGW_SOLVER_PREC != 0
    }

    /**
     * @brief Apply the active preconditioner:  z = M^{-1} r.
     *
     * tmp is a scratch vector distinct from r and z (the caller passes the
     * GMRES work vector that is currently free).  In-place r == z is safe:
     * the forward sweep writes only tmp, and the backward sweep consumes
     * tmp while overwriting z in reverse row order, so entries of r are
     * never read after being overwritten.
     *
     * Mode 1 applies the ILU(0) factors: forward solve with the unit lower
     * triangle (no stored diagonal), backward solve with the upper triangle
     * whose diagonal is stored already inverted (multiplication).  Both
     * sweeps are strictly serial for the same dependency reason as the
     * factorisation; each costs a fraction of a millisecond at production
     * grid sizes.
     *
     * Mode 2 applies the Thomas factors of each vertical column block;
     * columns are independent so the sweep parallelises over nx*ny blocks.
     *
     * If the factorisation is unhealthy or SERGHEI_RTGW_SOLVER_PREC == 0, this
     * reduces to the plain Jacobi apply.
     */
    void RT_Preconditioner_Apply(RTMatrix &rtA, const realArr &r, realArr &z, realArr &tmp)
    {
#if SERGHEI_RTGW_SOLVER_PREC == 1
        if (prec_active)
        {
            const int n = rtA.rt_nrow;
            auto lu = prec_lu;
            auto uinv = prec_uinv;
            auto diagj_l = prec_diag_j;
            auto ptr_l = rtA.rt_ptr;
            auto ind_l = rtA.rt_ind;
            Kokkos::parallel_for("RT_ILU0_Lsolve",
                Kokkos::RangePolicy<execution_space>(0, 1),
                KOKKOS_LAMBDA(int) {
                    for (int i = 0; i < n; ++i)
                    {
                        real s = r(i);
                        for (int k = ptr_l(i); k < diagj_l(i); ++k)
                            s -= lu(k) * tmp(ind_l(k));
                        tmp(i) = s;
                    } });
            Kokkos::parallel_for("RT_ILU0_Usolve",
                Kokkos::RangePolicy<execution_space>(0, 1),
                KOKKOS_LAMBDA(int) {
                    for (int i = n - 1; i >= 0; --i)
                    {
                        real s = tmp(i);
                        for (int k = diagj_l(i) + 1; k < ptr_l(i + 1); ++k)
                            s -= lu(k) * z(ind_l(k));
                        z(i) = s * uinv(i);
                    } });
            return;
        }
#elif SERGHEI_RTGW_SOLVER_PREC == 2
        if (prec_active)
        {
            auto l = prec_bj_l;
            auto u = prec_bj_u;
            auto sup_l = prec_sup_j;
            auto aval = rtA.rt_val;
            const int plane = prec_plane;
            const int nz = prec_nz;
            const int ncol = prec_ncol;
            Kokkos::parallel_for("RT_BJ_Apply", ncol, KOKKOS_LAMBDA(int col) {
                tmp(col) = r(col);
                for (int kk = 1; kk < nz; ++kk)
                {
                    const int i = col + kk * plane;
                    tmp(i) = r(i) - l(i) * tmp(i - plane);
                }
                {
                    const int i = col + (nz - 1) * plane;
                    z(i) = tmp(i) / u(i);
                }
                for (int kk = nz - 2; kk >= 0; --kk)
                {
                    const int i = col + kk * plane;
                    const real c = (sup_l(i) >= 0) ? aval(sup_l(i)) : 0.0;
                    z(i) = (tmp(i) - c * z(i + plane)) / u(i);
                } });
            return;
        }
#endif
        Jacobi_Precondition(rtA, r, z);
    }

    /**
     * @brief Assemble the RTGW linear-solver statistics report block.
     *
     * Called once from the end-of-run timer report (serghei.h passes the
     * string into writeLogFile) so the block lands in log.out right after
     * the LOAD BALANCE section.  The counters are always accumulated
     * (host-side increments only, no effect on the numerical path); the
     * phase timers carry data only when profiling is enabled (compile-time
     * SERGHEI_RTGW_SOLVER_PROFILE or the exported env variable), and profiled runs
     * must not be used for final performance numbers.
     *
     * t_total deliberately measures the same span as the rtgw.linsol
     * timer around each solve call, so the reported t_total is a direct
     * cross-check of that timer.
     */
    std::string solverStatsReport()
    {
        std::ostringstream os;
        os << "-----RTGW LINEAR SOLVER STATISTICS-----" << std::endl;
        os << "rtgw.lsv.prec_mode\t" << SERGHEI_RTGW_SOLVER_PREC << std::endl;
        os << "rtgw.lsv.spmv_mode\t" << SERGHEI_RTGW_SPMV_KOKKOS_MODE << std::endl;
        os << "rtgw.lsv.matrix\tnrow\t" << lsv_nrow << "\tnnz\t" << lsv_nnz << std::endl;

        // Aggregate over species
        long ns = 0, ni = 0, nrst = 0, ninst = 0, nncv = 0;
        double ttot = 0.0;
        std::vector<int> allIters;
        std::vector<double> allEps;
        double emin = 1e300, emax = 0.0;
        for (const LSvSpecStats &s : lsv)
        {
            ns += s.n_solves; ni += s.n_iters; nrst += s.restarts;
            ninst += s.instant_exits; nncv += s.nonconv; ttot += s.t_total;
            emin = (s.eps_min < emin) ? s.eps_min : emin;
            emax = (s.eps_max > emax) ? s.eps_max : emax;
            allIters.insert(allIters.end(), s.iters.begin(), s.iters.end());
            allEps.insert(allEps.end(), s.eps.begin(), s.eps.end());
        }
        os << "rtgw.lsv.n_solves\t" << ns << std::endl;
        os << "rtgw.lsv.iters_total\t" << ni << std::endl;
        if (ns > 0)
        {
            std::sort(allIters.begin(), allIters.end());
            std::sort(allEps.begin(), allEps.end());
            os << "rtgw.lsv.iters_per_solve\tmedian\t" << allIters[ns / 2]
               << "\tmean\t" << (double)ni / (double)ns
               << "\tp90\t" << allIters[(ns * 9) / 10]
               << "\tmax\t" << allIters[ns - 1] << std::endl;
            os << "rtgw.lsv.instant_exits\t" << ninst << "\t(" << 100.0 * (double)ninst / (double)ns << "%)" << std::endl;
            os << "rtgw.lsv.restarts\t" << nrst << std::endl;
            os << "rtgw.lsv.nonconv\t" << nncv << std::endl;
            os << "rtgw.lsv.eps\tmedian\t" << allEps[ns / 2]
               << "\tmin\t" << emin << "\tmax\t" << emax << std::endl;
        }
        os << "rtgw.lsv.t_total\t" << ttot << "\t(cross-check vs rtgw.linsol)" << std::endl;
        os << "rtgw.lsv.prec_fallbacks\t" << lsv_prec_fallbacks << std::endl;
        os << "rtgw.lsv.half_step\tStep1 solves\t" << lsv_half_solves[0] << "\titers\t" << lsv_half_iters[0]
           << "\tt\t" << lsv_half_time[0]
           << "\t| Step3 solves\t" << lsv_half_solves[1] << "\titers\t" << lsv_half_iters[1]
           << "\tt\t" << lsv_half_time[1] << std::endl;

        // Phase breakdown (meaningful only when profiled)
        if (lsv_profile)
        {
            double pre = 0, inr = 0, spmv = 0, ort = 0, giv = 0, bsub = 0;
            for (const LSvSpecStats &s : lsv)
            { pre += s.t_precond; inr += s.t_init_res; spmv += s.t_spmv_apply;
              ort += s.t_ortho; giv += s.t_givens; bsub += s.t_backsub; }
            const double psum = pre + inr + spmv + ort + giv + bsub;
            os << "rtgw.lsv.phases\tt_precond\t" << pre << "\t(" << 100.0 * pre / psum << "%)"
               << "\tt_init_res\t" << inr << "\t(" << 100.0 * inr / psum << "%)" << std::endl;
            os << "rtgw.lsv.phases\tt_spmv_apply\t" << spmv << "\t(" << 100.0 * spmv / psum << "%)"
               << "\tt_ortho\t" << ort << "\t(" << 100.0 * ort / psum << "%)" << std::endl;
            os << "rtgw.lsv.phases\tt_givens\t" << giv << "\t(" << 100.0 * giv / psum << "%)"
               << "\tt_backsub\t" << bsub << "\t(" << 100.0 * bsub / psum << "%)" << std::endl;
        }
        else
        {
            os << "rtgw.lsv.phases\tnot profiled (build with -DSERGHEI_RT_PROFILE=1 or export SERGHEI_RTGW_SOLVER_PROFILE=1)" << std::endl;
        }

        // Per-species lines
        for (size_t is = 0; is < lsv.size(); ++is)
        {
            const LSvSpecStats &s = lsv[is];
            if (s.n_solves == 0) continue;
            os << "rtgw.lsv.spec" << is << "\tsolves\t" << s.n_solves
               << "\titers\t" << s.n_iters
               << "\tmean\t" << (double)s.n_iters / (double)s.n_solves
               << "\tinstant\t" << s.instant_exits
               << "\tnonconv\t" << s.nonconv
               << "\tt\t" << s.t_total << std::endl;
        }
        return os.str();
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
    int Jacobi_GMRES_Solve(RTMatrix &rtA, bool print_diag = false, const char *diag_half = "",
                           int diag_iSpec = -1, int diag_picard_iter = -1,
                           real diag_etime = 0.0, real diag_gdom_dt = 0.0, real diag_rt_dt = 0.0)
    {

        int m_restart = gmres_m_restart; ///< Krylov subspace dimension before restart
        const int n = rtA.rt_nrow;
        const real eps_min_local = eps_min;  ///< Convergence tolerance for GMRES
        int iter_total = 0;
        real eps_local = 1.0;
        int restart_count = 0;
        Kokkos::Timer gmres_timer;

        // References to work vectors in the matrix object
        realArr &rt_z = rtA.rt_z;
        realArr &rt_w = rtA.rt_w;

        // Opt-in per-solve diagnostics for preconditioner A/B tests:
        // enabled either at compile time (-DSERGHEI_RT_GMRES_DIAG=1, const.h)
        // or per run without recompiling by exporting SERGHEI_RTGW_SOLVER_GMRES_DIAG=1,
        // turning on the existing [RTGW_GMRES_DIAG] channel (iteration counts,
        // residual, timing).
        static const bool env_diag = (SERGHEI_RTGW_SOLVER_GMRES_DIAG != 0) ||
                                     (std::getenv("SERGHEI_RTGW_SOLVER_GMRES_DIAG") != nullptr);
        if (!print_diag && env_diag) print_diag = true;

        // [P3-4] Solver statistics: resolve the per-species bucket and the
        // phase-profiling switch (compile-time -DSERGHEI_RT_PROFILE=1 in
        // const.h, or the exported SERGHEI_RTGW_SOLVER_PROFILE environment variable;
        // either one enables).  Counters below are host-side increments
        // only (always on); the phase timers only read Kokkos::Timer when
        // profiling is on, so the default path is untouched.
        static const bool env_prof = (SERGHEI_RTGW_SOLVER_PROFILE != 0) ||
                                     (std::getenv("SERGHEI_RTGW_SOLVER_PROFILE") != nullptr);
        lsv_profile = env_prof;
        lsv_nrow = n;
        lsv_nnz = rtA.rt_nnz;
        int lsv_spec = (diag_iSpec >= 0) ? diag_iSpec : 0;
        if ((int)lsv.size() <= lsv_spec) lsv.resize(lsv_spec + 1);
        LSvSpecStats &lsv_s = lsv[lsv_spec];
        const int lsv_half = (diag_half && diag_half[4] == '3') ? 1 : 0; // "Step1" / "Step3"
        bool lsv_instant = false;
        Kokkos::Timer ph_timer;

        ph_timer.reset();
        Update_Jacobi_Inverse_Diagonal(rtA);
#if SERGHEI_RTGW_SOLVER_PREC != 0
        // Factorise the upgraded preconditioner for this freshly
        // assembled system (per species / Picard iteration).  On failure
        // every apply below degrades to Jacobi for this solve.
        RT_Preconditioner_Factorize(rtA);
#endif
        if (lsv_profile) lsv_s.t_precond += ph_timer.seconds();
        if (print_diag)
            std::cerr << "[RTGW_GMRES_DIAG]"
                      << " event=GMRES_BEGIN"
                      << " half=" << diag_half
                      << " etime=" << diag_etime
                      << " gdom_dt=" << diag_gdom_dt
                      << " rt_dt=" << diag_rt_dt
                      << " iSpec=" << diag_iSpec
                      << " picard_iter=" << diag_picard_iter
                      << " nrow=" << n
                      << " nnz=" << rtA.rt_nnz
                      << " iter_max=" << iter_max
                      << " restart_m=" << m_restart
                      << std::endl;

        auto &v = gmres_v;
        auto &H = gmres_H;
        auto &c_arr = gmres_c;
        auto &s_arr = gmres_s;
        auto &g_arr = gmres_g;
        auto &y = gmres_y;

        // Main GMRES iteration loop (with restarts)
        while (iter_total < iter_max && eps_local > eps_min_local)
        {
            restart_count++;
            ph_timer.reset();

            // Compute initial preconditioned residual: z = M^{-1}(b - A*x)
#if SERGHEI_RTGW_SOLVER_PREC == 0
            // [P1-5a] Fused residual + Jacobi preconditioner in one kernel.
            // Per element the operation sequence is identical to the previous
            // get_residual() followed by Jacobi_Precondition() (the residual
            // is simply kept in a register instead of round-tripping through
            // rt_r), so the result is bitwise identical.  rt_r is no longer
            // written here; it has no other reader inside GMRES.
            Kokkos::parallel_for("RT_Precond_Residual", n, KOKKOS_LAMBDA(int idx) {
                real r = rtA.rt_rhs(idx);
                for (int icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx + 1); ++icol) {
                    r -= rtA.rt_val(icol) * rtA.rt_x(rtA.rt_ind(icol));
                }
                rt_z(idx) = r * rtA.rt_diag(idx); });
#else
            // Split residual and preconditioner apply so the upgraded
            // preconditioner can be used.  rt_w is free here as the apply's
            // scratch vector.
            Kokkos::parallel_for("RT_Residual", n, KOKKOS_LAMBDA(int idx) {
                real r = rtA.rt_rhs(idx);
                for (int icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx + 1); ++icol) {
                    r -= rtA.rt_val(icol) * rtA.rt_x(rtA.rt_ind(icol));
                }
                rt_z(idx) = r; });
            RT_Preconditioner_Apply(rtA, rt_z, rt_z, rt_w);
#endif

            // Compute the norm of the preconditioned residual (beta)
            real beta = std::sqrt(dot(rt_z, rt_z, n));
            if (lsv_profile) lsv_s.t_init_res += ph_timer.seconds();
            eps_local = beta;
            if (print_diag)
                std::cerr << "[RTGW_GMRES_DIAG]"
                          << " event=GMRES_RESTART_BEGIN"
                          << " half=" << diag_half
                          << " etime=" << diag_etime
                          << " iSpec=" << diag_iSpec
                          << " picard_iter=" << diag_picard_iter
                          << " restart=" << restart_count
                          << " iter_total=" << iter_total
                          << " beta=" << beta
                          << " elapsed_s=" << gmres_timer.seconds()
                          << std::endl;
            if (beta < eps_min_local)
            {
                lsv_instant = (restart_count == 1); // [P3-4] converged at the restart-start check
                break;
            }

            // Normalise z to obtain the first Krylov basis vector v[0]
            // [P1-5c] Fused normalise + copy: v[0] = z / beta in a single
            // kernel (division kept, so bitwise identical to the previous
            // "rt_z /= beta" followed by deep_copy(v0, rt_z)).  rt_z itself
            // is no longer normalised in place; it is recomputed at the top
            // of the next restart and not read in between.
            auto v0 = Kokkos::subview(v, 0, Kokkos::ALL());
            Kokkos::parallel_for("RT_Init_V0", n, KOKKOS_LAMBDA(const int i) { v0(i) = rt_z(i) / beta; });

            // Initialise the projected RHS: g = ||beta|| * e_1
            g_arr(0) = beta;
            for (int i = 1; i <= m_restart; ++i)
                g_arr(i) = 0.0;

            // Arnoldi iteration: build the Krylov basis
            int k;
            for (k = 0; k < m_restart; ++k)
            {
                auto vk = Kokkos::subview(v, k, Kokkos::ALL());
                ph_timer.reset();

                // Compute w = M^{-1} * A * v[k]  (preconditioned matrix-vector product)
#if SERGHEI_RTGW_SOLVER_PREC == 0
#if SERGHEI_RTGW_SPMV_KOKKOS_MODE == 1
                // Optional KokkosSparse::spmv backend.  The library
                // kernel uses a different reduction order than the fused
                // kernel below (L2, ~1 ulp class differences); keep only if
                // the A/B benchmark shows a clear rtgw.linsol win.
                KokkosSparse::spmv("N", 1.0, rt_crs_A, vk, 0.0, rt_w);
                Jacobi_Precondition(rtA, rt_w, rt_w);
#else
                // [P1-5b] Fused SpMV + Jacobi preconditioner in one kernel.
                // Per element this performs exactly acc = sum(A(i,:) * v);
                // w(i) = acc * diag(i), i.e. the same operations in the same
                // order as mul_MV() followed by Jacobi_Precondition(), so the
                // result is bitwise identical while saving one full pass
                // (read + write) over rt_w per Arnoldi step.
                Kokkos::parallel_for("RT_SpMV_Jacobi", n, KOKKOS_LAMBDA(int idx) {
                    real acc = 0.0;
                    for (int icol = rtA.rt_ptr(idx); icol < rtA.rt_ptr(idx + 1); ++icol) {
                        acc += rtA.rt_val(icol) * vk(rtA.rt_ind(icol));
                    }
                    rt_w(idx) = acc * rtA.rt_diag(idx); });
#endif
#else
                // Unpreconditioned SpMV followed by the upgraded
                // preconditioner apply; rt_z is free here as scratch.
#if SERGHEI_RTGW_SPMV_KOKKOS_MODE == 1
                KokkosSparse::spmv("N", 1.0, rt_crs_A, vk, 0.0, rt_w);
#else
                mul_MV(rt_w, rtA, vk);
#endif
                RT_Preconditioner_Apply(rtA, rt_w, rt_w, rt_z);
#endif
                if (lsv_profile) lsv_s.t_spmv_apply += ph_timer.seconds();

                // Modified Gram-Schmidt orthogonalisation against all previous
                // basis vectors v[0], ..., v[k]
                ph_timer.reset();
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
                if (lsv_profile) lsv_s.t_ortho += ph_timer.seconds();

                // Apply all previous Givens rotations to the new column of H
                ph_timer.reset();
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
                if (lsv_profile) lsv_s.t_givens += ph_timer.seconds();
                if (print_diag && ((k + 1) % 10 == 0 || eps_local < eps_min_local))
                    std::cerr << "[RTGW_GMRES_DIAG]"
                              << " event=GMRES_PROGRESS"
                              << " half=" << diag_half
                              << " etime=" << diag_etime
                              << " iSpec=" << diag_iSpec
                              << " picard_iter=" << diag_picard_iter
                              << " restart=" << restart_count
                              << " inner_iter=" << (k + 1)
                              << " iter_total_est=" << (iter_total + k + 1)
                              << " eps=" << eps_local
                              << " elapsed_s=" << gmres_timer.seconds()
                              << std::endl;
                if (eps_local < eps_min_local)
                {
                    k++;
                    break;
                }
            }

            // Back-substitution to solve the upper-triangular system
            // H(0:k, 0:k) * y = g(0:k)
            ph_timer.reset();
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
            if (lsv_profile) lsv_s.t_backsub += ph_timer.seconds();

            iter_total += k;
            eps_local = std::fabs(g_arr(k));
            if (print_diag)
                std::cerr << "[RTGW_GMRES_DIAG]"
                          << " event=GMRES_RESTART_END"
                          << " half=" << diag_half
                          << " etime=" << diag_etime
                          << " iSpec=" << diag_iSpec
                          << " picard_iter=" << diag_picard_iter
                          << " restart=" << restart_count
                          << " iter_total=" << iter_total
                          << " eps=" << eps_local
                          << " elapsed_s=" << gmres_timer.seconds()
                          << std::endl;
        }
        // [P3-4] Per-solve counter accumulation (host-side only).
        const double lsv_t = gmres_timer.seconds();
        lsv_s.n_solves++;
        lsv_s.n_iters += iter_total;
        lsv_s.restarts += restart_count;
        lsv_s.iters.push_back(iter_total);
        lsv_s.eps.push_back((double)eps_local);
        lsv_s.t_total += lsv_t;
        if (eps_local < lsv_s.eps_min) lsv_s.eps_min = eps_local;
        if (eps_local > lsv_s.eps_max) lsv_s.eps_max = eps_local;
        if (lsv_instant) lsv_s.instant_exits++;
        if (eps_local > eps_min_local) lsv_s.nonconv++; // left loop on iter_max
        lsv_half_solves[lsv_half]++;
        lsv_half_iters[lsv_half] += iter_total;
        lsv_half_time[lsv_half] += lsv_t;

        iter = iter_total;
        eps = eps_local;
        if (print_diag)
            std::cerr << "[RTGW_GMRES_DIAG]"
                      << " event=GMRES_END"
                      << " half=" << diag_half
                      << " etime=" << diag_etime
                      << " iSpec=" << diag_iSpec
                      << " picard_iter=" << diag_picard_iter
                      << " restart_count=" << restart_count
                      << " iter_total=" << iter_total
                      << " eps=" << eps_local
                      << " elapsed_s=" << gmres_timer.seconds()
                      << std::endl;
        return iter_total;
    }
};

#endif
