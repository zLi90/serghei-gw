/**
 * @file RTFunctionSW.h
 * @brief Surface water reactive transport solver using Strang operator splitting.
 *
 * This file implements the RTFunctionSW class, which provides a complete
 * reactive transport framework for dissolved species in surface water flows.
 * The solver uses a Strang splitting approach to separate advection-dispersion
 * (transport) from biogeochemical reactions within each global timestep:
 *
 *   1. Transport half-step  (dt/2)
 *   2. Full reaction step    (dt)   -- optional, nitrogen cycle only
 *   3. Transport half-step  (dt/2)
 *
 * == Transport Solver ==
 * Conservative finite-volume method (FVM) on a uniform 2D grid with:
 *   - Riemann-interface water fluxes from the hydrodynamic solver
 *   - TVD (Total Variation Diminishing) concentration reconstruction at
 *     cell interfaces for monotone, second-order advection
 *   - Optional molecular + mechanical dispersion (longitudinal / transverse)
 *   - Rainfall source terms and external boundary conditions
 *
 * == Reaction Solver (Nitrogen Cycle, 5 species) ==
 *   Species indices: 0 = NH4 (ammonium), 1 = NO3 (nitrate), 2 = DO (dissolved oxygen),
 *                    3 = DOC (dissolved organic carbon), 4 = DON (dissolved organic nitrogen)
 *
 *   Monod kinetics with temperature correction for:
 *     1. Mineralization   (DON -> NH4, consuming DO)
 *     2. Nitrification    (NH4 -> NO3, consuming DO)
 *     3. Denitrification  (NO3 -> N2, consuming DOC, inhibited by DO)
 *     4. Aerobic respiration (DOC consumption, consuming DO)
 *     5. Surface-water specific sources: atmospheric reaeration (DO),
 *        DOC/DON release from sediments
 *
 *   Fully implicit integration via Newton-Raphson iteration with an
 *   analytically computed Jacobian and LU-decomposition linear solver.
 *
 * == Surface-Subsurface Exchange ==
 *   Computes solute mass exchange between surface water and groundwater
 *   based on the volumetric exchange flux qss from the coupled hydrologic model.
 *
 * == MPI Parallelism ==
 *   Halo exchange for solute concentrations is handled by RTExchangeSW,
 *   called at each stage boundary of the Strang splitting to ensure
 *   consistency across MPI domain decompositions.
 */

#ifndef _RT_FUNCTION_SW_H_
#define _RT_FUNCTION_SW_H_

#if SERGHEI_SURFACE_TRANSPORT

#include <Kokkos_Core.hpp>
#include "define.h"

#include "const.h"
#include "BC.h"
#include "Domain.h"
#include "State.h"
#include "SArray.h"
#include "Indexing.h"
#include "RTStateSW.h"
#include "RTExchangeSW.h"
#include "Parallel.h"
#include <set>
#include <math.h>

class RTFunctionSW
{

private:
    RTExchangeSW rtexch; ///< MPI halo exchanger for surface water solute fields

    /* -----------------------------------------------------------------------
       Constants for the 5-species nitrogen cycle Newton solver
       Species ordering: 0:NH4, 1:NO3, 2:DO, 3:DOC, 4:DON
       ----------------------------------------------------------------------- */
    static const int MAX_SPEC = 5;                 ///< Number of chemical species
    static const int MAX_NEWTON_ITER = 10;         ///< Maximum Newton-Raphson iterations
    static constexpr real NEWTON_TOLERANCE = 1e-6; ///< Convergence tolerance on residual norm
    static constexpr real SMALL_CONC = 1e-15;      ///< Floor value to avoid division by zero in Monod terms

    /* -----------------------------------------------------------------------
       calculate_residual
       -----------------------------------------------------------------------
       Computes the implicit residual vector for the 5-species nitrogen cycle
       reaction system using the trapezoidal (Crank-Nicolson) rule:

         Resid[k] = C_new[k] - C_old[k] - 0.5 * dt * (dCdt_new[k] + dCdt_old[k])

       The residual is the function F(C_new) = 0 that the Newton solver drives
       to zero.  The trapezoidal rule gives second-order accuracy in time.

       @param C_new      [in]  Candidate new concentrations at the end of the timestep (5 elements)
       @param C_old      [in]  Known concentrations at the start of the timestep (5 elements)
       @param dt         [in]  Reaction sub-timestep length [s]
       @param Resid      [out] Residual vector R = C_new - C_old - 0.5*dt*(f(C_new)+f(C_old))
       @param rt         [in]  RTStateSW object containing all kinetic and stoichiometric parameters
       @param Temp_Corr  [in]  Temperature correction factor (theta^(T_opt - T_env))
    ----------------------------------------------------------------------- */
    KOKKOS_INLINE_FUNCTION
    void calculate_residual(const real C_new[], const real C_old[], real dt, real Resid[],
                            const RTStateSW &rt, real Temp_Corr)
    {
        // --- Surface-water physical parameters ---
        const real Ka_sw = rt.Ka_sw;         ///< Surface reaeration rate coefficient [1/s]
        const real DO_sat = rt.DO_sat;       ///< DO saturation concentration [mg/L]
        const real DOC_eq_sw = rt.DOC_eq_sw; ///< Surface DOC equilibrium concentration [mg/L]
        const real DON_eq_sw = rt.DON_eq_sw; ///< Surface DON equilibrium concentration [mg/L]
        const real K_rel_sw = rt.K_rel_sw;   ///< Surface DOC/DON release rate coefficient [1/s]

        // --- Stoichiometric yield coefficients (mg product / mg reactant) ---
        // These are derived from molecular weight ratios:
        //   NH4 (18.04 g/mol), NO3 (62.00 g/mol), O2 (32.00 g/mol),
        //   N in NH4 (14.01 g/mol), C in DOC (12.01 g/mol)
        const real Y_NH4_DON = 18.04 / 14.01;         ///< mg NH4 produced per mg DON mineralized
        const real Y_NO3_NH4 = 62.00 / 18.04;         ///< mg NO3 produced per mg NH4 nitrified
        const real Y_O2_NH4 = 64.00 / 18.04;          ///< mg O2 consumed per mg NH4 nitrified
        const real Y_O2_DON = 32.00 / 14.01;          ///< mg O2 consumed per mg DON mineralized
        const real Y_O2_DOC = 32.00 / 12.01;          ///< mg O2 consumed per mg DOC respired aerobically
        const real Y_NO3_DOC = (0.8 * 62.00) / 12.01; ///< mg NO3 consumed per mg DOC during denitrification (0.8 efficiency factor)

        // Lambda: compute the reaction rate vector dCdt from Monod kinetics
        auto calc_dCdt = [&](const real C[], real dCdt[])
        {
            // Clamp concentrations to non-negative to avoid spurious kinetics
            real NH4 = fmax(0.0, C[0]);
            real NO3 = fmax(0.0, C[1]);
            real DO = fmax(0.0, C[2]);
            real DOC = fmax(0.0, C[3]);
            real DON = fmax(0.0, C[4]);

            // --- Reaction 1: Mineralization (DON -> NH4) ---
            // Ammonification of dissolved organic nitrogen, consuming dissolved oxygen.
            // Rate depends on DON and DO via dual Monod limitation.
            real R_min = rt.Rate_Max_Min * Temp_Corr * (DON / (rt.K_Monod_DON + DON + SMALL_CONC)) * (DO / (rt.K_Monod_DO_Min + DO + SMALL_CONC));

            // --- Reaction 2: Nitrification (NH4 -> NO3) ---
            // Oxidation of ammonium to nitrate, consuming dissolved oxygen.
            real R_nit = rt.Rate_Max_Nit * Temp_Corr * (NH4 / (rt.K_Monod_NH4 + NH4 + SMALL_CONC)) * (DO / (rt.K_Monod_DO + DO + SMALL_CONC));

            // --- Reaction 3: Denitrification (NO3 -> N2) ---
            // Reduction of nitrate under anoxic conditions, consuming DOC.
            // Inhibited by dissolved oxygen via a non-competitive inhibition term.
            real R_denit = rt.Rate_Max_Denit * Temp_Corr * (NO3 / (rt.K_Monod_NO3 + NO3 + SMALL_CONC)) * (DOC / (rt.K_Monod_DOC + DOC + SMALL_CONC)) * (rt.Ki_Inhib_DO / (rt.Ki_Inhib_DO + DO + SMALL_CONC));

            // --- Reaction 4: Aerobic respiration (DOC consumption) ---
            // Heterotrophic consumption of dissolved organic carbon, consuming DO.
            real R_aer = rt.Rate_Max_Hetero * Temp_Corr * (DOC / (rt.K_Monod_DOC_Aerobic + DOC + SMALL_CONC)) * (DO / (rt.K_Monod_DO_Aerobic + DO + SMALL_CONC));

            // --- Reaction 5: Surface-water specific source terms ---
            // Note: Unlike the subsurface module, the surface module does not multiply
            // reaeration by a moisture/saturation factor (always water-covered).
            real R_air = Ka_sw * fmax(0.0, DO_sat - DO);            ///< Atmospheric reaeration of DO
            real R_rel_DOC = K_rel_sw * fmax(0.0, DOC_eq_sw - DOC); ///< DOC release from sediments
            real R_rel_DON = K_rel_sw * fmax(0.0, DON_eq_sw - DON); ///< DON release from sediments

            // --- Assemble species rate equations (mass balance) ---
            dCdt[0] = Y_NH4_DON * R_min - R_nit;                                      ///< NH4: gains from mineralization, loses to nitrification
            dCdt[1] = Y_NO3_NH4 * R_nit - Y_NO3_DOC * R_denit;                        ///< NO3: gains from nitrification, loses to denitrification
            dCdt[2] = R_air - Y_O2_DON * R_min - Y_O2_NH4 * R_nit - Y_O2_DOC * R_aer; ///< DO: reaeration, consumed by mineralization, nitrification, respiration
            dCdt[3] = R_rel_DOC - R_denit - R_aer;                                    ///< DOC: sediment release, consumed by denitrification and respiration
            dCdt[4] = R_rel_DON - R_min;                                              ///< DON: sediment release, consumed by mineralization
        };

        // Evaluate reaction rates at both old and new concentration states
        real dCdt_old[5] = {0}, dCdt_new[5] = {0};
        calc_dCdt(C_old, dCdt_old);
        calc_dCdt(C_new, dCdt_new);

        // Trapezoidal (Crank-Nicolson) residual: R = C_new - C_old - 0.5*dt*(f_new + f_old)
        for (int k = 0; k < 5; k++)
        {
            Resid[k] = C_new[k] - C_old[k] - 0.5 * dt * (dCdt_new[k] + dCdt_old[k]);
        }
    }

    /* -----------------------------------------------------------------------
       calculate_jacobian
       -----------------------------------------------------------------------
       Computes the 5x5 analytical Jacobian matrix J = dResid/dC_new for the
       reaction residual.  Only the partial derivatives with respect to C_new
       are needed because C_old is treated as constant in the Newton iteration.

       The Jacobian structure is:
         J[i][j] = delta_ij - 0.5 * dt * (sum of stoichiometric-weighted rate derivatives)

       Each Monod term d(C/(K+C))/dC = K / (K+C)^2 is computed analytically
       via the deriv_monod helper lambda.

       @param C_new      [in]  Current candidate concentrations (5 elements)
       @param dt         [in]  Reaction sub-timestep [s]
       @param Jacobian   [out] 5x5 Jacobian matrix J_ij = dResid_i / dC_new_j
       @param rt         [in]  RTStateSW with kinetic parameters
       @param Temp_Corr  [in]  Temperature correction factor
    ----------------------------------------------------------------------- */
    KOKKOS_INLINE_FUNCTION
    void calculate_jacobian(const real C_new[], real dt, real Jacobian[MAX_SPEC][MAX_SPEC],
                            const RTStateSW &rt, real Temp_Corr)
    {
        const real Ka_sw = rt.Ka_sw;
        const real K_rel_sw = rt.K_rel_sw;

        // Stoichiometric yields (same as in calculate_residual)
        const real Y_NH4_DON = 18.04 / 14.01;
        const real Y_NO3_NH4 = 62.00 / 18.04;
        const real Y_O2_NH4 = 64.00 / 18.04;
        const real Y_O2_DON = 32.00 / 14.01;
        const real Y_O2_DOC = 32.00 / 12.01;
        const real Y_NO3_DOC = (0.8 * 62.00) / 12.01;

        // Clamp concentrations to non-negative
        real NH4 = fmax(0.0, C_new[0]);
        real NO3 = fmax(0.0, C_new[1]);
        real DO = fmax(0.0, C_new[2]);
        real DOC = fmax(0.0, C_new[3]);
        real DON = fmax(0.0, C_new[4]);

        // Helper: derivative of Monod term C/(K+C) with respect to C is K/(K+C)^2
        auto deriv_monod = [](real C, real K)
        { return K / pow(K + C + SMALL_CONC, 2); };

        // --- Partial derivatives of reaction rates with respect to each species ---

        // Mineralization rate derivatives
        real dRmin_dDON = rt.Rate_Max_Min * Temp_Corr * deriv_monod(DON, rt.K_Monod_DON) * (DO / (rt.K_Monod_DO_Min + DO + SMALL_CONC));
        real dRmin_dDO = rt.Rate_Max_Min * Temp_Corr * (DON / (rt.K_Monod_DON + DON + SMALL_CONC)) * deriv_monod(DO, rt.K_Monod_DO_Min);

        // Nitrification rate derivatives
        real dRnit_dNH4 = rt.Rate_Max_Nit * Temp_Corr * deriv_monod(NH4, rt.K_Monod_NH4) * (DO / (rt.K_Monod_DO + DO + SMALL_CONC));
        real dRnit_dDO = rt.Rate_Max_Nit * Temp_Corr * (NH4 / (rt.K_Monod_NH4 + NH4 + SMALL_CONC)) * deriv_monod(DO, rt.K_Monod_DO);

        // Denitrification rate derivatives (includes DO inhibition term I_DO)
        real I_DO = rt.Ki_Inhib_DO / (rt.Ki_Inhib_DO + DO + SMALL_CONC);
        real dI_DO_dDO = -rt.Ki_Inhib_DO / pow(rt.Ki_Inhib_DO + DO + SMALL_CONC, 2);
        real dRdenit_dNO3 = rt.Rate_Max_Denit * Temp_Corr * deriv_monod(NO3, rt.K_Monod_NO3) * (DOC / (rt.K_Monod_DOC + DOC + SMALL_CONC)) * I_DO;
        real dRdenit_dDOC = rt.Rate_Max_Denit * Temp_Corr * (NO3 / (rt.K_Monod_NO3 + NO3 + SMALL_CONC)) * deriv_monod(DOC, rt.K_Monod_DOC) * I_DO;
        real dRdenit_dDO = rt.Rate_Max_Denit * Temp_Corr * (NO3 / (rt.K_Monod_NO3 + NO3 + SMALL_CONC)) * (DOC / (rt.K_Monod_DOC + DOC + SMALL_CONC)) * dI_DO_dDO;

        // Aerobic respiration rate derivatives
        real dRaer_dDOC = rt.Rate_Max_Hetero * Temp_Corr * deriv_monod(DOC, rt.K_Monod_DOC_Aerobic) * (DO / (rt.K_Monod_DO_Aerobic + DO + SMALL_CONC));
        real dRaer_dDO = rt.Rate_Max_Hetero * Temp_Corr * (DOC / (rt.K_Monod_DOC_Aerobic + DOC + SMALL_CONC)) * deriv_monod(DO, rt.K_Monod_DO_Aerobic);

        // Atmospheric reaeration derivative: d(Ka*(DO_sat - DO))/dDO = -Ka
        real dRair_dDO = -Ka_sw;

        // Initialize Jacobian as identity matrix (from the C_new - C_old term in residual)
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                Jacobian[i][j] = (i == j) ? 1.0 : 0.0;

        // --- Fill in the reaction Jacobian entries ---
        // Each entry: J[i][j] -= 0.5 * dt * (chain-rule expansion of d(dCdt_i)/dC_j)

        // Row 0: NH4 equation
        Jacobian[0][0] -= 0.5 * dt * (-dRnit_dNH4);                       ///< dR_NH4/dNH4: nitrification loses NH4
        Jacobian[0][2] -= 0.5 * dt * (Y_NH4_DON * dRmin_dDO - dRnit_dDO); ///< dR_NH4/dDO: mineralization gains, nitrification loses
        Jacobian[0][4] -= 0.5 * dt * (Y_NH4_DON * dRmin_dDON);            ///< dR_NH4/dDON: mineralization gains NH4

        // Row 1: NO3 equation
        Jacobian[1][0] -= 0.5 * dt * (Y_NO3_NH4 * dRnit_dNH4);                          ///< dR_NO3/dNH4: nitrification gains NO3
        Jacobian[1][1] -= 0.5 * dt * (-Y_NO3_DOC * dRdenit_dNO3);                       ///< dR_NO3/dNO3: denitrification loses NO3
        Jacobian[1][2] -= 0.5 * dt * (Y_NO3_NH4 * dRnit_dDO - Y_NO3_DOC * dRdenit_dDO); ///< dR_NO3/dDO: nitrification gains, denitrification inhibition
        Jacobian[1][3] -= 0.5 * dt * (-Y_NO3_DOC * dRdenit_dDOC);                       ///< dR_NO3/dDOC: denitrification loses NO3

        // Row 2: DO equation
        Jacobian[2][0] -= 0.5 * dt * (-Y_O2_NH4 * dRnit_dNH4);                                                         ///< dR_DO/dNH4: nitrification consumes DO
        Jacobian[2][2] -= 0.5 * dt * (dRair_dDO - Y_O2_DON * dRmin_dDO - Y_O2_NH4 * dRnit_dDO - Y_O2_DOC * dRaer_dDO); ///< dR_DO/dDO: reaeration, mineralization, nitrification, respiration
        Jacobian[2][3] -= 0.5 * dt * (-Y_O2_DOC * dRaer_dDOC);                                                         ///< dR_DO/dDOC: respiration consumes DO
        Jacobian[2][4] -= 0.5 * dt * (-Y_O2_DON * dRmin_dDON);                                                         ///< dR_DO/dDON: mineralization consumes DO

        // Row 3: DOC equation
        Jacobian[3][1] -= 0.5 * dt * (-dRdenit_dNO3);                         ///< dR_DOC/dNO3: denitrification consumes DOC
        Jacobian[3][2] -= 0.5 * dt * (-dRdenit_dDO - dRaer_dDO);              ///< dR_DOC/dDO: denitrification inhibition, respiration
        Jacobian[3][3] -= 0.5 * dt * (-K_rel_sw - dRdenit_dDOC - dRaer_dDOC); ///< dR_DOC/dDOC: sediment release, denitrification, respiration

        // Row 4: DON equation
        Jacobian[4][2] -= 0.5 * dt * (-dRmin_dDO);             ///< dR_DON/dDO: mineralization consumes DON
        Jacobian[4][4] -= 0.5 * dt * (-K_rel_sw - dRmin_dDON); ///< dR_DON/dDON: sediment release, mineralization
    }

    /* -----------------------------------------------------------------------
       solve_small_system_NxN
       -----------------------------------------------------------------------
       Solves a small dense linear system  J * Delta_C = -Resid  using
       LU decomposition without pivoting (Doolittle algorithm), followed by
       forward and backward substitution.

       The system size is at most 5x5, so pivoting is not critical for
       performance.  The matrix is overwritten with its LU factors in-place.

       Algorithm:
         1. LU decomposition: decompose A into lower (L) and upper (U) triangular
            matrices stored in-place.  L has unit diagonal.
         2. Forward substitution: solve  L * y = -Resid
         3. Backward substitution: solve  U * Delta_C = y

       @param Jacobian  [in/out]  The Jacobian matrix (n_spec x n_spec), overwritten with LU factors
       @param Resid     [in]      The residual vector (n_spec)
       @param Delta_C   [out]     The Newton update vector (n_spec)
       @param n_spec    [in]      Number of species (system dimension)
    ----------------------------------------------------------------------- */
    KOKKOS_INLINE_FUNCTION
    void solve_small_system_NxN(real Jacobian[MAX_SPEC][MAX_SPEC], real Resid[],
                                real Delta_C[], int n_spec)
    {
        // Copy Jacobian into working array to avoid modifying the original
        real A[MAX_SPEC][MAX_SPEC];
        for (int i = 0; i < n_spec; i++)
            for (int j = 0; j < n_spec; j++)
                A[i][j] = Jacobian[i][j];

        // --- Step 1: LU decomposition (Doolittle, no pivoting) ---
        // After this loop, A stores L (below diagonal, unit diagonal) and U (on and above diagonal)
        for (int k = 0; k < n_spec; k++)
        {
            // Compute U[k][j] for j >= k
            for (int j = k; j < n_spec; j++)
            {
                real sum = 0.0;
                for (int p = 0; p < k; p++)
                    sum += A[k][p] * A[p][j];
                A[k][j] -= sum;
            }
            // Compute L[i][k] for i > k
            for (int i = k + 1; i < n_spec; i++)
            {
                real sum = 0.0;
                for (int p = 0; p < k; p++)
                    sum += A[i][p] * A[p][k];
                A[i][k] = (A[i][k] - sum) / A[k][k];
            }
        }

        // --- Step 2: Forward substitution  L * y = -Resid ---
        real y[MAX_SPEC];
        for (int i = 0; i < n_spec; i++)
        {
            real sum = 0.0;
            for (int j = 0; j < i; j++)
                sum += A[i][j] * y[j];
            y[i] = (-Resid[i] - sum); // L has unit diagonal, so no division needed
        }

        // --- Step 3: Backward substitution  U * Delta_C = y ---
        for (int i = n_spec - 1; i >= 0; i--)
        {
            real sum = 0.0;
            for (int j = i + 1; j < n_spec; j++)
                sum += A[i][j] * Delta_C[j];
            Delta_C[i] = (y[i] - sum) / A[i][i];
        }
    }

public:
    /* -----------------------------------------------------------------------
       initializeExchange
       -----------------------------------------------------------------------
       Allocates MPI halo exchange buffers for the surface water solute fields.
       Must be called once during initialization before any transport steps.

       @param rtsw  [in/out]  Surface water reactive transport state
       @param dom   [in]      Computational domain (grid and partition info)
    ----------------------------------------------------------------------- */
    void initializeExchange(RTStateSW &rtsw, Domain &dom)
    {
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering initializeExchange() n_mass=" << rtsw.n_mass << std::endl;
#endif
        rtexch.allocate(dom, rtsw.n_mass);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> initializeExchange() completed OK" << std::endl;
#endif
    }

    /* -----------------------------------------------------------------------
       exchangeMPI
       -----------------------------------------------------------------------
       Performs MPI halo exchange for all solute species in the surface water
       reactive transport state.  Delegates to RTExchangeSW.

       @param rtsw  [in/out]  Surface water reactive transport state
       @param dom   [in]      Computational domain
       @param par   [in]      Parallel decomposition information
    ----------------------------------------------------------------------- */
    void exchangeMPI(RTStateSW &rtsw, Domain &dom, Parallel &par)
    {
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering exchangeMPI()" << std::endl;
#endif
        rtexch.exchangeMPI(rtsw, dom, par, 1);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> exchangeMPI() completed OK" << std::endl;
#endif
    }

    /* -----------------------------------------------------------------------
       calculate_dt
       -----------------------------------------------------------------------
       Computes the maximum allowable transport timestep based on the CFL
       (Courant-Friedrichs-Lewy) condition for advection:

         dt = CFL * dx / max(|velocity|)

       where CFL is set to 0.5 for stability, and the maximum velocity is
       determined by a parallel reduction over all grid cells.

       @param state  [in]  Hydrodynamic state (hu, hv, h for velocity computation)
       @param dom    [in]  Domain with grid spacing dx
       @return            CFL-limited timestep [s]
    ----------------------------------------------------------------------- */
    real calculate_dt(State &state, Domain &dom)
    {
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering calculate_dt()" << std::endl;
#endif
        real max_vel = 0;
        // Parallel reduction to find the maximum velocity magnitude across all cells
        Kokkos::parallel_reduce("max_vel", dom.nCellMem, KOKKOS_LAMBDA(const int iGlob, real &maxv) {
            // Compute depth-averaged velocity components from momentum (hu, hv) and depth (h)
            real u = state.hu4rtsw(iGlob,0)/state.h4rtsw(iGlob,0);

            real v = state.hv4rtsw(iGlob,0)/state.h4rtsw(iGlob,0);
            real vel = sqrt(u*u + v*v);
            if (vel > maxv) maxv = vel; }, Kokkos::Max<real>(max_vel));

        const real CFL = 0.5;                        ///< CFL number for stability
        real dt = CFL * dom.dx() / (max_vel + TOL6); // TOL6 prevents division by zero
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> calculate_dt() completed: max_vel=" << max_vel << " dt=" << dt << std::endl;
#endif
        return dt;
    }

public:
    // =========================================================================
    //  Conservative finite-volume advection-dispersion solver (FVM)
    //
    //  Solves the transport equation for each dissolved species using:
    //    - Riemann water fluxes at cell interfaces (from hydrodynamic solver)
    //    - TVD concentration reconstruction for second-order advection
    //    - Hydrodynamic dispersion (longitudinal + transverse) and molecular diffusion
    //    - Rainfall source terms and external boundary conditions
    //
    //  The kernel is "fused": all species are processed inside a single
    //  parallel_for over grid cells, which maximizes data locality for the
    //  shared hydrodynamic quantities (h, hu, hv, interface fluxes).
    //
    //  @param rtsw      [in/out] Surface water reactive transport state (concentrations updated)
    //  @param state     [in]     Hydrodynamic state (h, hu, hv, interface fluxes)
    //  @param dom       [in]     Computational domain
    //  @param ebc       [in]     External boundary conditions for reactive transport
    //  @param ss        [in]     Source/sink data (rainfall, etc.)
    //  @param step_part [in]     Strang splitting part indicator:
    //                              1 = first half-step (dt/2): h_start = h_old, h_end = h_mid
    //                              2 = second half-step (dt/2): h_start = h_mid, h_end = h_new
    // =========================================================================
    void transport(RTStateSW &rtsw, State &state, Domain &dom, ExternalBoundaries &ebc, const SourceSinkData &ss, int step_part)
    {
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering transport() step_part=" << step_part << " n_mass=" << rtsw.n_mass << " dt=" << rtsw.dt << std::endl;
#endif
        const real dx = dom.dx();
        const real dy = dom.dx();
        const real dt = rtsw.dt;
        const real hmin = state.hmin;

        // --- Fused transport kernel ---
        // All species are processed per cell to reuse hydrodynamic data (h, Q, etc.)
        // from cache rather than re-reading from global memory for each species.
        Kokkos::parallel_for("Fused_SW_Transport", dom.nCellMem, KOKKOS_LAMBDA(const int iGlob) {
            int i, j, nxhc, nyhc;
            nxhc = dom.nx + 2 * dom.hc;
            nyhc = dom.ny + 2 * dom.hc;
            unpackIndicesUniformGrid(iGlob, nyhc, nxhc, j, i);

            // --- Step 1: Pre-read hydrodynamic data shared by all species ---
            // Time slot 1 = previous state, time slot 0 = current state
            real h_old = state.h4rtsw(iGlob, 1);
            real h_new = state.h4rtsw(iGlob, 0);
            // Interpolate water depth at start and end of this transport sub-step
            // based on the Strang splitting position (step_part 1 or 2)
            real h_start = (step_part == 1) ? h_old : 0.5 * (h_old + h_new);
            real h_end   = (step_part == 1) ? 0.5 * (h_old + h_new) : h_new;

            // Compute neighbor cell indices (or -1 at domain boundaries)
            int east  = (i < nxhc - 1) ? iGlob + 1 : -1;
            int west  = (i > 0)        ? iGlob - 1 : -1;
            int south = (j < nyhc - 1) ? iGlob + nxhc : -1;
            int north = (j > 0)        ? iGlob - nxhc : -1;

            // Pre-read Riemann water fluxes at the four cell interfaces
            real Q_E = (east != -1) ? state.InterfaceFlux_x(iGlob) : 0.0;   ///< East face flux (positive = eastward)
            real Q_W = (west != -1) ? state.InterfaceFlux_x(west)  : 0.0;   ///< West face flux (positive = eastward)
            real Q_S = (south != -1) ? state.InterfaceFlux_y(iGlob) : 0.0;  ///< South face flux (positive = southward)
            real Q_N = (north != -1) ? state.InterfaceFlux_y(north) : 0.0;  ///< North face flux (positive = southward)

            // Determine flow direction at each interface for upwinding
            bool flow_E = (Q_E >= 0.0);
            bool flow_W = (Q_W >= 0.0);
            bool flow_S = (Q_S >= 0.0);
            bool flow_N = (Q_N >= 0.0);

            // Pre-read neighbor water depths (for dispersion face-averaging)
            real h_E = (east != -1) ? state.h4rtsw(east, 1) : 0.0;
            real h_W = (west != -1) ? state.h4rtsw(west, 1) : 0.0;
            real h_S = (south != -1) ? state.h4rtsw(south, 1) : 0.0;
            real h_N = (north != -1) ? state.h4rtsw(north, 1) : 0.0;

            // --- Step 2: Loop over all species within this cell ---
            // This inner loop maximizes cache reuse for the hydrodynamic data read above.
            for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++) {

                // Initialize solute mass from start-of-step concentration and water depth
                real mass = 0.0;
                if (h_start > hmin) {
                    mass = h_start * rtsw.c(iSpec, iGlob, 0);
                }

                // ================== X-direction fluxes ==================

                // East face flux (from current cell iGlob to east neighbor)
                if (east != -1 && (h_old > hmin || h_E > hmin)) {
                    // TVD reconstruction of concentration at the east face
                    real C_face = reconstruct_TVD_concentration(rtsw, iSpec, iGlob, east, flow_E, dom, i, j, 1);
                    real adv_flux = Q_E * C_face;
                    // Dispersion flux: D * h * dC/dx (only if face has water)
                    real diff_flux = 0.0;
                    real h_face = 0.5 * (h_old + h_E);
                    if (h_face > hmin) {
                        // Compute face velocity components for dispersion tensor
                        real u_f = Q_E / h_face;
                        real v_f = 0.5 * (state.hv4rtsw(iGlob,1) + state.hv4rtsw(east,1)) / h_face;
                        real vel = sqrt(u_f*u_f + v_f*v_f);
                        // Dispersion tensor: D = D_mol + alpha_T * |v| + (alpha_L - alpha_T) * u^2 / |v|
                        real Dxx = rtsw.diffusion_molecular(iSpec);
                        if (vel > TOL6) Dxx += rtsw.alpha_T(iSpec) * vel + (rtsw.alpha_L(iSpec) - rtsw.alpha_T(iSpec)) * u_f * u_f / vel;
                        diff_flux = - h_face * Dxx * (rtsw.c(iSpec, east, 0) - rtsw.c(iSpec, iGlob, 0)) / dx;
                    }
                    mass -= (adv_flux + diff_flux) * dt / dx; // Outflux from this cell
                }

                // West face flux (from west neighbor into current cell)
                if (west != -1 && (h_W > hmin || h_old > hmin)) {
                    // TVD reconstruction at the west face (west->iGlob interface)
                    real C_face = reconstruct_TVD_concentration(rtsw, iSpec, west, iGlob, flow_W, dom, i-1, j, 1);
                    real adv_flux = Q_W * C_face;
                    real diff_flux = 0.0;
                    real h_face = 0.5 * (h_W + h_old);
                    if (h_face > hmin) {
                        real u_f = Q_W / h_face;
                        real v_f = 0.5 * (state.hv4rtsw(west,1) + state.hv4rtsw(iGlob,1)) / h_face;
                        real vel = sqrt(u_f*u_f + v_f*v_f);
                        real Dxx = rtsw.diffusion_molecular(iSpec);
                        if (vel > TOL6) Dxx += rtsw.alpha_T(iSpec) * vel + (rtsw.alpha_L(iSpec) - rtsw.alpha_T(iSpec)) * u_f * u_f / vel;
                        diff_flux = - h_face * Dxx * (rtsw.c(iSpec, iGlob, 0) - rtsw.c(iSpec, west, 0)) / dx;
                    }
                    mass += (adv_flux + diff_flux) * dt / dx; // Influx to this cell
                }

                // ================== Y-direction fluxes ==================

                // South face flux (from current cell to south neighbor)
                if (south != -1 && (h_old > hmin || h_S > hmin)) {
                    real C_face = reconstruct_TVD_concentration(rtsw, iSpec, iGlob, south, flow_S, dom, i, j, 2);
                    real adv_flux = Q_S * C_face;
                    real diff_flux = 0.0;
                    real h_face = 0.5 * (h_old + h_S);
                    if (h_face > hmin) {
                        real v_f = Q_S / h_face;
                        real u_f = 0.5 * (state.hu4rtsw(iGlob,1) + state.hu4rtsw(south,1)) / h_face;
                        real vel = sqrt(u_f*u_f + v_f*v_f);
                        // Dispersion tensor in Y: use v_f for longitudinal component
                        real Dyy = rtsw.diffusion_molecular(iSpec);
                        if (vel > TOL6) Dyy += rtsw.alpha_T(iSpec) * vel + (rtsw.alpha_L(iSpec) - rtsw.alpha_T(iSpec)) * v_f * v_f / vel;
                        diff_flux = - h_face * Dyy * (rtsw.c(iSpec, south, 0) - rtsw.c(iSpec, iGlob, 0)) / dy;
                    }
                    mass -= (adv_flux + diff_flux) * dt / dy;
                }

                // North face flux (from north neighbor into current cell)
                if (north != -1 && (h_N > hmin || h_old > hmin)) {
                    real C_face = reconstruct_TVD_concentration(rtsw, iSpec, north, iGlob, flow_N, dom, i, j-1, 2);
                    real adv_flux = Q_N * C_face;
                    real diff_flux = 0.0;
                    real h_face = 0.5 * (h_N + h_old);
                    if (h_face > hmin) {
                        real v_f = Q_N / h_face;
                        real u_f = 0.5 * (state.hu4rtsw(north,1) + state.hu4rtsw(iGlob,1)) / h_face;
                        real vel = sqrt(u_f*u_f + v_f*v_f);
                        real Dyy = rtsw.diffusion_molecular(iSpec);
                        if (vel > TOL6) Dyy += rtsw.alpha_T(iSpec) * vel + (rtsw.alpha_L(iSpec) - rtsw.alpha_T(iSpec)) * v_f * v_f / vel;
                        diff_flux = - h_face * Dyy * (rtsw.c(iSpec, iGlob, 0) - rtsw.c(iSpec, north, 0)) / dy;
                    }
                    mass += (adv_flux + diff_flux) * dt / dy;
                }

                // ================== Source/sink: rainfall ==================
                // Rain adds solute mass at the rain concentration for the given species
                if (dom.isRain && h_end > hmin) {
                    mass += ss.rainRate(iGlob) * rtsw.RainCon_arr(iSpec) * dt;
                }

                const real MIN_PHYSICAL_MASS = TOL15; ///< Truncation threshold for mass floor

                // --- Update concentration from computed mass ---
                if (h_end > state.hmin) {
                    if (mass < MIN_PHYSICAL_MASS) {
                        // Mass is effectively zero; set concentration to zero
                        rtsw.c(iSpec, iGlob, 1) = 0.0;
                        rtsw.solute_mass(iSpec, iGlob, 1) = 0.0;
                    } else {
                        real calculated_c = mass / h_end;
                        rtsw.c(iSpec, iGlob, 1) = calculated_c;
                        // Back-calculate mass from concentration to ensure strict conservation
                        rtsw.solute_mass(iSpec, iGlob, 1) = calculated_c * h_end;
                    }
                } else {
                    // Water depth below threshold: no solute in this dry cell
                    rtsw.c(iSpec, iGlob, 1) = 0.0;
                    rtsw.solute_mass(iSpec, iGlob, 1) = 0.0;
                }
            } });
        Kokkos::fence();
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> transport() Fused_SW_Transport kernel completed" << std::endl;
#endif

        // --- External source/sink terms (separate from the fused kernel) ---
        // Processed per-species; includes boundary inflows, point sources, etc.
        for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
        {
            RTSwSS::ComputeRTSWSourceSink(rtsw, state, dom, ss, iSpec);
        }
        Kokkos::fence();
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> transport() source/sink computed for " << rtsw.n_mass << " species" << std::endl;
#endif

        // --- Apply reactive transport boundary conditions ---
        for (int k = 0; k < ebc.rtbc.size(); k++)
        {
            ebc.rtbc[k].applyrtbc(rtsw, dom, state);
        }
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> transport() BC applied (" << ebc.rtbc.size() << " BC segments)" << std::endl;
#endif

        // --- Synchronize time slots: copy slot 1 (new) to slot 0 (current) ---
        // This prepares the state for the next sub-step or the reaction solver.
        Kokkos::parallel_for("update_state", dom.nCellMem, KOKKOS_LAMBDA(const int iGlob) {
            for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++) {
                rtsw.c(iSpec, iGlob, 0) = rtsw.c(iSpec, iGlob, 1);
                rtsw.solute_mass(iSpec, iGlob, 0) = rtsw.solute_mass(iSpec, iGlob, 1);
                rtsw.csw4gw(iSpec, iGlob) = rtsw.c(iSpec, iGlob, 1); ///< Also update concentration visible to groundwater coupling
            } });
        Kokkos::fence();
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> transport() step_part=" << step_part << " completed OK" << std::endl;
#endif
    }

    /* -----------------------------------------------------------------------
       computeGwExchange_For_RT
       -----------------------------------------------------------------------
       Computes the solute mass exchange between surface water and groundwater
       based on the volumetric water exchange flux qss from the coupled
       hydrologic model.

       For groundwater -> surface water (q_vol > 0):
         - If surface water is dry (h_mid < hmin): concentration is set directly
           from the influx concentration
         - If surface water is wet: implicit update using mixing equation:
             C_new = (C_old + dt * flux_mass / h_mid) / (1 + dt * q_vol / h_mid)

       For surface water -> groundwater (q_vol < 0):
         - Surface water concentration remains unchanged (loss is handled by
           the groundwater module)

       @param state  [in]     Hydrodynamic state (qss = volumetric exchange flux)
       @param rtsw   [in/out] Surface water RT state (concentrations updated)
       @param dom    [in]     Computational domain
    ----------------------------------------------------------------------- */
    inline void computeGwExchange_For_RT(State &state, RTStateSW &rtsw, const Domain &dom)
    {
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering computeGwExchange_For_RT() n_mass=" << rtsw.n_mass << std::endl;
#endif
        Kokkos::parallel_for(dom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii = dom.getIndex(idom);
            real q_vol = state.qss(idom); ///< Volumetric water exchange flux [m/s] (positive = GW -> SW)
            real h = state.h4rtsw(ii, 0);
            // Compute mid-step water depth for the exchange calculation
            real h_old = state.h4rtsw(ii, 1);
            real h_new = state.h4rtsw(ii, 0);
            real h_mid = 0.5 * (h_old + h_new);
            const real h_min = state.hmin;

            for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
            {
                real flux_mass = rtsw.ConQss(iSpec, idom); ///< Solute mass flux concentration from groundwater [mg/L]

                if (q_vol > 0.0) // Groundwater inflow to surface water
                {
                    real C_old = rtsw.c(iSpec, ii, 0);
                    if (h_mid < h_min) {
                        // Dry surface cell: concentration is set by the GW inflow
                        rtsw.c(iSpec, ii, 0) = flux_mass / q_vol;
                    } else {
                        // Wet surface cell: implicit mixing update
                        real numerator = C_old + (rtsw.dt * flux_mass / h_mid);
                        real denominator = 1.0 + (rtsw.dt * q_vol / h_mid);
                        rtsw.c(iSpec, ii, 0) = numerator / denominator;
                    }
                }
                else if (q_vol < 0.0)
                {
                    // Surface water outflow to groundwater: SW concentration unchanged
                }

                // Safety: clamp to non-negative and guard against NaN
                if (rtsw.c(iSpec, ii, 0) < 0.0) rtsw.c(iSpec, ii, 0) = 0.0;
                if (isnan(rtsw.c(iSpec, ii, 0))) rtsw.c(iSpec, ii, 0) = 0.0;

                // Synchronize both time slots after the exchange
                rtsw.c(iSpec, ii, 1) = rtsw.c(iSpec, ii, 0);

                // Immediately recalculate solute_mass from the updated concentration
                // to maintain consistency for the subsequent second transport half-step.
                if (h_mid > h_min) {
                    rtsw.solute_mass(iSpec, ii, 0) = h_mid * rtsw.c(iSpec, ii, 0);
                    rtsw.solute_mass(iSpec, ii, 1) = rtsw.solute_mass(iSpec, ii, 0);
                } else {
                    rtsw.solute_mass(iSpec, ii, 0) = 0.0;
                    rtsw.solute_mass(iSpec, ii, 1) = 0.0;
                }
            } });
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> computeGwExchange_For_RT() completed OK" << std::endl;
#endif
    }

    /* -----------------------------------------------------------------------
       reaction_integration_kernel
       -----------------------------------------------------------------------
       Performs fully implicit reaction integration for the 5-species nitrogen
       cycle using Newton-Raphson iteration.

       Algorithm:
         1. Skip dry cells (h < hmin) and cells with negligible concentrations
         2. Initialize C_new = C_old
         3. Newton iteration loop (up to MAX_NEWTON_ITER):
            a. Compute residual F(C_new)
            b. Check convergence: ||F|| < NEWTON_TOLERANCE
            c. Compute analytical Jacobian J = dF/dC_new
            d. Solve linear system J * Delta = -F via LU decomposition
            e. Update C_new += Delta, clamping negatives to small positive value
         4. If converged, write updated concentrations and masses back

       The mid-step water depth h_mid (from the Strang splitting) is used to
       compute solute mass, ensuring mass consistency with the transport steps.

       @param rtsw  [in/out] Surface water reactive transport state (concentrations updated)
       @param state [in]     Hydrodynamic state (water depth)
       @param dom   [in]     Computational domain
    ----------------------------------------------------------------------- */
    void reaction_integration_kernel(RTStateSW &rtsw, State &state, Domain &dom)
    {
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering reaction_integration_kernel() n_mass=" << rtsw.n_mass << " NCS=" << rtsw.Nitrogen_Cycle_Simulation << std::endl;
#endif

        Kokkos::parallel_for("SW_Reaction", dom.nCellMem, KOKKOS_LAMBDA(int iGlob) {

            real h = state.h4rtsw(iGlob, 0);
            // Compute mid-step water depth for mass consistency in Strang splitting
            real h_old = state.h4rtsw(iGlob, 1);
            real h_new = state.h4rtsw(iGlob, 0);
            real h_mid = 0.5 * (h_old + h_new); // Mid-point depth from the Strang splitting
            // Skip dry cells
            if (h_mid < state.hmin) return;

            real dt = rtsw.dt;
            real T_env = 20.0;  // Reference environmental temperature [deg C]
            // Temperature correction factor: theta^(T_opt - T_env)
            real Temp_Corr = pow(rtsw.Temp_Coeff_Theta, T_env - rtsw.Opt_Temp);

            const int n_spec = rtsw.n_mass;
            // Require at least 5 species for the nitrogen cycle reaction system
            if (n_spec < 5) return;

            // --- Early exit check: skip reaction if all concentrations are negligible ---
            bool has_substance = false;
            for (int k = 0; k < n_spec; k++) {
                if (rtsw.c(k, iGlob, 1) > TOL12) { // Physical zero threshold
                    has_substance = true;
                    break;
                }
            }
            if (!has_substance) {
                // Force all concentrations to exact zero to prevent floating-point noise
                for (int k = 0; k < n_spec; k++) rtsw.c(k, iGlob, 1) = 0.0;
                return;
            }

            // --- Allocate local arrays for Newton solver ---
            real C_old[MAX_SPEC];
            real C_new[MAX_SPEC];
            real Resid[MAX_SPEC];
            real Jacobian[MAX_SPEC][MAX_SPEC];
            real Delta_C[MAX_SPEC];

            // Initialize: C_new starts from C_old (the Newton initial guess)
            for (int k = 0; k < n_spec; k++)
            {
                C_old[k] = rtsw.c(k, iGlob, 1);
                C_new[k] = C_old[k];
            }

            // --- Newton-Raphson iteration ---
            bool newton_converged = false;
            for (int iter = 0; iter < MAX_NEWTON_ITER; iter++)
            {
                // Step (a): Compute the residual F(C_new)
                calculate_residual(C_new, C_old, dt, Resid, rtsw, Temp_Corr);

                // Step (b): Check convergence using L2 norm of residual
                real norm_resid = 0.0;
                for (int k = 0; k < n_spec; k++)
                    norm_resid += Resid[k] * Resid[k];
                if (sqrt(norm_resid) < NEWTON_TOLERANCE)
                {
                    newton_converged = true;
                    break;
                }

                // Step (c): Compute the analytical Jacobian J = dF/dC_new
                calculate_jacobian(C_new, dt, Jacobian, rtsw, Temp_Corr);
                // Step (d): Solve the linear system J * Delta = -Resid via LU decomposition
                solve_small_system_NxN(Jacobian, Resid, Delta_C, n_spec);

                // Step (e): Update C_new with damped/clamped correction
                for (int k = 0; k < n_spec; k++)
                {
                    real next_val = C_new[k] + Delta_C[k];
                    // Clamp negative concentrations to a small positive value
                    if (next_val < 0.0)
                        next_val = 1e-15;
                    C_new[k] = next_val;
                }
            }

            // If Newton iteration failed to converge, report error and skip this cell
            if (!newton_converged)
            {
                printf("[ERROR] RTGW Newton iteration failed to converge at cell iGlob=%d after %d iterations. Residual norm too large. Skipping cell.\n", iGlob, MAX_NEWTON_ITER);
                return;
            }

            // --- Write converged concentrations and masses back ---
            for (int k = 0; k < n_spec; k++)
            {
                rtsw.c(k, iGlob, 1) = C_new[k];
                rtsw.c(k, iGlob, 0) = C_new[k];

                // Use the mid-step water depth h_mid to compute solute mass.
                // This ensures that the mass passed to the second transport
                // half-step is consistent with the Strang splitting.
                rtsw.solute_mass(k, iGlob, 1) = h_mid * C_new[k];
                rtsw.solute_mass(k, iGlob, 0) = h_mid * C_new[k];
            } });
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> reaction_integration_kernel() completed OK" << std::endl;
#endif
    }

    /* -----------------------------------------------------------------------
       rt_solve_sw
       -----------------------------------------------------------------------
       Top-level Strang splitting solver for surface water reactive transport.

       Implements the Strang operator-splitting scheme for the coupled
       advection-dispersion-reaction equations:

         Step 1: Transport half-step  (dt/2)  -- advection + dispersion
         Step 2: Full reaction step    (dt)   -- nitrogen cycle biogeochemistry
         Step 3: Transport half-step  (dt/2)  -- advection + dispersion

       MPI halo exchanges are performed at every stage boundary to ensure
       that each sub-solver sees up-to-date concentration data at domain
       boundaries.

       Timing is recorded per stage (transport, reaction, communication)
       into the domain timer structure for performance profiling.

       @param rtsw    [in/out] Surface water reactive transport state
       @param state   [in]     Hydrodynamic state
       @param dom     [in/out] Computational domain (timers updated)
       @param ebc     [in]     External boundary conditions for reactive transport
       @param ss      [in]     Source/sink data
       @param rtintsw [in/out] RT integrator for flux bookkeeping (mass before/after reaction)
       @param par     [in]     Parallel decomposition information
    ----------------------------------------------------------------------- */
    void rt_solve_sw(RTStateSW &rtsw, State &state, Domain &dom, ExternalBoundaries &ebc,
                     const SourceSinkData &ss, RTIntegratorSW &rtintsw, Parallel &par)
    {
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rt_solve_sw() dt_global=" << dom.dt << " n_mass=" << rtsw.n_mass << std::endl;
#endif
        real dt_global = dom.dt;
        Kokkos::Timer timer;

        // Set sub-timestep to half the global timestep for the first transport step
        rtsw.dt = dt_global * 0.5;

        // --- Strang Step 1: First transport half-step (dt/2) ---
        // Pre-transport MPI halo exchange
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 1] Pre-transport MPI exchange" << std::endl;
#endif
        timer.reset();
        rtexch.exchangeMPI(rtsw, dom, par, 0);
        rtexch.exchangeMPI(rtsw, dom, par, 1);
        Kokkos::fence();
        dom.timers.rtsw.comm += timer.seconds();

#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 1] Entering transport() half-step 1" << std::endl;
#endif
        timer.reset();
        transport(rtsw, state, dom, ebc, ss, 1); // step_part = 1
        Kokkos::fence();
        dom.timers.rtsw.transport += timer.seconds();
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 1] transport() half-step 1 completed" << std::endl;
#endif

        // Post-transport halo exchange: ensure reaction solver reads up-to-date concentrations
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 1] Post-transport MPI exchange" << std::endl;
#endif
        timer.reset();
        rtexch.exchangeMPI(rtsw, dom, par, 0);
        rtexch.exchangeMPI(rtsw, dom, par, 1);
        Kokkos::fence();
        dom.timers.rtsw.comm += timer.seconds();

        // --- Strang Step 2: Full reaction step (dt) ---
        // Only executed if the nitrogen cycle simulation flag is enabled
        if (rtsw.Nitrogen_Cycle_Simulation == 1)
        {
            rtsw.dt = dt_global; // Full timestep for the reaction solver
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
            std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 2] Entering reaction_integration_kernel()" << std::endl;
#endif
            // Record pre-reaction mass for flux bookkeeping
            rtintsw.recordMassBefore(rtsw, state, dom);

            timer.reset();
            reaction_integration_kernel(rtsw, state, dom);
            Kokkos::fence();
            dom.timers.rtsw.react += timer.seconds();
            // Compute the reaction-induced flux from pre/post mass difference
            rtintsw.finalizeReactionFlux(rtsw, state, dom, dt_global);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
            std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 2] reaction_integration_kernel() completed" << std::endl;
#endif

            // Post-reaction halo exchange
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
            std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 2] Post-reaction MPI exchange" << std::endl;
#endif
            timer.reset();
            rtexch.exchangeMPI(rtsw, dom, par, 0);
            rtexch.exchangeMPI(rtsw, dom, par, 1);
            Kokkos::fence();
            dom.timers.rtsw.comm += timer.seconds();
        }

        // --- Strang Step 3: Second transport half-step (dt/2) ---
        rtsw.dt = dt_global * 0.5;

#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 3] Entering transport() half-step 2" << std::endl;
#endif
        timer.reset();
        transport(rtsw, state, dom, ebc, ss, 2); // step_part = 2
        Kokkos::fence();
        dom.timers.rtsw.transport += timer.seconds();
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> [Strang Step 3] transport() half-step 2 completed" << std::endl;
#endif

        // Final halo exchange after the complete Strang cycle
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Final MPI exchange" << std::endl;
#endif
        timer.reset();
        rtexch.exchangeMPI(rtsw, dom, par, 0);
        rtexch.exchangeMPI(rtsw, dom, par, 1);
        Kokkos::fence();
        dom.timers.rtsw.comm += timer.seconds();

        // Restore full timestep
        rtsw.dt = dt_global;

        // Integrator hook (currently unused but reserved for future extensions)
        timer.reset();
        // rtintsw.integrate(rtsw, state, dom, ss, ebc);
        Kokkos::fence();
        dom.timers.rtsw.integrate += timer.seconds();

#if DEBUG_SERGHEI_SURFACE_TRANSPORT
        std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve_sw() completed OK" << std::endl;
#endif
    }

    // =========================================================================
    //  compute_flux
    //
    //  Generalized interface flux computation with dynamic scheme switching.
    //  Supports three advection schemes:
    //    scheme_type = 1: First-order upwind
    //    scheme_type = 2: TVD with Van Leer limiter
    //    scheme_type = 3: TVD with Superbee limiter
    //
    //  For TVD schemes, the reconstructed face concentration is:
    //    C_face = C_upwind + 0.5 * phi * (1 - Co) * (C_downwind - C_upwind)
    //  where Co is the local Courant number and phi is the flux limiter.
    //
    //  @param c3           [in]  Three-cell stencil: c3(0)=upwind-of-upwind,
    //                             c3(1)=upwind, c3(2)=downwind
    //  @param u            [in]  Interface velocity (signed)
    //  @param dt           [in]  Timestep
    //  @param dl           [in]  Grid spacing in the flux direction
    //  @param scheme_type  [in]  Advection scheme selector (1, 2, or 3)
    //  @return                  Reconstructed concentration at the cell interface
    // =========================================================================
    KOKKOS_INLINE_FUNCTION real compute_flux(const SArray<real, 3> &c3, real u, real dt, real dl, int scheme_type)
    {
        real c0 = c3(0); // Upwind of upwind
        real c1 = c3(1); // Upwind cell
        real c2 = c3(2); // Downwind cell

        // First-order upwind: simply return the upwind cell concentration
        if (scheme_type == 1)
        {
            return c1;
        }

        // Compute local Courant number, clamped to < 1 for stability
        real co = fabs(u) * dt / dl;
        if (co >= 1.0)
            co = 0.999; // Safety clamp to prevent instability

        // Compute the upwind ratio r = (C_upwind - C_up-upwind) / (C_downwind - C_upwind)
        real r = 0.0;
        if (fabs(c2 - c1) > 1e-12)
        {
            r = (c1 - c0) / (c2 - c1);
        }

        real phi = 0.0;

        if (scheme_type == 2)
        {
            // --- TVD Van Leer limiter ---
            // phi(r) = 2r / (1 + r) for r > 0, else 0
            if (r > 0.0)
            {
                phi = 2.0 * r / (1.0 + r);
            }
            phi = fmin(phi, 2.0);
            phi = fmax(phi, 0.0);
        }
        else if (scheme_type == 3)
        {
            // --- TVD Superbee limiter ---
            // phi(r) = max(min(2r,1), min(r,2))
            real r1 = (2.0 * r < 1.0) ? 2.0 * r : 1.0;
            real r2 = (r < 2.0) ? r : 2.0;
            phi = (r1 > r2) ? r1 : r2;
            phi = fmin(phi, 2.0);
            if (r <= 0.0)
                phi = 0.0;
        }

        // Lax-Wendroff flux form: C_face = C_upwind + 0.5 * phi * (1 - Co) * (C_down - C_up)
        // The sign convention is unified because the upwind cell is already selected
        // by the flow-direction logic in the calling transport routine.
        return c1 + 0.5 * phi * (1.0 - co) * (c2 - c1);
    }

    // =========================================================================
    //  reconstruct_TVD_concentration
    //
    //  Reconstructs the concentration at a cell interface using the TVD scheme
    //  specified by rtsw.Advection_Scheme.  This function identifies the
    //  three-cell stencil (upwind-of-upwind, upwind, downwind) based on the
    //  flow direction, then applies the selected flux limiter.
    //
    //  The upwind-of-upwind cell may not exist at domain boundaries; in that
    //  case the reconstruction degrades gracefully to first-order upwind by
    //  setting c0 = c1 (the upwind cell value).
    //
    //  @param rtsw            [in]  Surface water RT state (concentrations, scheme type)
    //  @param iSpec           [in]  Species index
    //  @param left_id         [in]  Global index of the left (lower) cell at the interface
    //  @param right_id        [in]  Global index of the right (upper) cell at the interface
    //  @param flow_positive   [in]  True if flow is in the positive direction (left -> right)
    //  @param dom             [in]  Computational domain
    //  @param i               [in]  X-index of the interface (used for boundary detection)
    //  @param j               [in]  Y-index of the interface (used for boundary detection)
    //  @param dir             [in]  Direction: 1 = X (east face), 2 = Y (south face)
    //  @return                      Reconstructed concentration at the cell interface
    // =========================================================================
    KOKKOS_INLINE_FUNCTION
    real reconstruct_TVD_concentration(const RTStateSW &rtsw, int iSpec, int left_id, int right_id, bool flow_positive, const Domain &dom, int i, int j, int dir)
    {
        int upwind_id, downwind_id, up_upwind_id = -1;
        int nxhc = dom.nx + 2 * dom.hc;

        // Identify the three-cell stencil based on flow direction
        if (flow_positive)
        {
            // Flow goes left -> right: upwind = left, downwind = right
            upwind_id = left_id;
            downwind_id = right_id;
            // Upwind-of-upwind is one cell further upwind
            if (dir == 1 && i > 0)
                up_upwind_id = left_id - 1; // X-direction: one cell to the west
            if (dir == 2 && j > 0)
                up_upwind_id = left_id - nxhc; // Y-direction: one cell to the north
        }
        else
        {
            // Flow goes right -> left: upwind = right, downwind = left
            upwind_id = right_id;
            downwind_id = left_id;
            // Upwind-of-upwind is one cell further upwind
            if (dir == 1 && i < nxhc - 2)
                up_upwind_id = right_id + 1; // X-direction: one cell to the east
            if (dir == 2 && j < dom.ny + 2 * dom.hc - 2)
                up_upwind_id = right_id + nxhc; // Y-direction: one cell to the south
        }

        // Read concentrations from the three-cell stencil
        real c1 = rtsw.c(iSpec, upwind_id, 0);   ///< Upwind cell concentration
        real c2 = rtsw.c(iSpec, downwind_id, 0); ///< Downwind cell concentration
        // At domain boundaries, degrade to first-order by duplicating the upwind value
        real c0 = (up_upwind_id != -1) ? rtsw.c(iSpec, up_upwind_id, 0) : c1;

        // First-order upwind: return upwind concentration directly
        if (rtsw.Advection_Scheme == 1)
            return c1;

        // Compute the upwind ratio r = (C_upwind - C_up-upwind) / (C_downwind - C_upwind)
        real r = 0.0;
        if (fabs(c2 - c1) > 1e-12)
        {
            r = (c1 - c0) / (c2 - c1);
        }

        real phi = 0.0;

        // TVD Van Leer limiter: phi(r) = 2r / (1+r)
        if (rtsw.Advection_Scheme == 2)
        {
            if (r > 0.0)
                phi = 2.0 * r / (1.0 + r);
            phi = fmin(phi, 2.0);
        }
        // TVD Superbee limiter: phi(r) = max(min(2r,1), min(r,2))
        else if (rtsw.Advection_Scheme == 3)
        {
            real r1 = (2.0 * r < 1.0) ? 2.0 * r : 1.0;
            real r2 = (r < 2.0) ? r : 2.0;
            phi = (r1 > r2) ? r1 : r2;
            if (r <= 0.0)
                phi = 0.0;
        }

        // Reconstructed interface concentration: C_face = C_upwind + 0.5 * phi * (C_downwind - C_upwind)
        return c1 + 0.5 * phi * (c2 - c1);
    }
};

#endif
#endif