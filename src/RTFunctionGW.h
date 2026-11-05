/**
 * @file RTFunctionGW.h
 * @brief Groundwater Reactive Transport Solver — function collection.
 *
 * This header implements the core numerical methods for simulating reactive
 * solute transport in variably-saturated groundwater systems.  It covers:
 *
 *   1. **Flux limiter functions** — Minmod, Superbee, and Van Leer limiters
 *      used by the TVD (Total Variation Diminishing) high-resolution scheme.
 *   2. **Dispersion tensor computation** — velocity coordinate transformation,
 *      inter-cell velocity interpolation, and the full 3-D hydrodynamic
 *      dispersion tensor (longitudinal + transverse + molecular diffusion).
 *   3. **Linear system assembly** — implicit advection-dispersion equation
 *      discretisation on a structured grid with upwind weighting, cross-
 *      dispersion terms, sorption (linear / Freundlich / Langmuir / kinetic),
 *      first-order decay, and TVD flux correction (Picard iteration).
 *   4. **Monod-based nitrogen cycle reactions** — mineralisation, nitrification,
 *      denitrification, aerobic respiration, soil re-aeration, and organic
 *      matter release, solved via a cell-local fully-implicit Newton-Raphson
 *      (Crank-Nicolson) integrator with analytical Jacobian.
 *   5. **Strang-splitting driver** — the main rt_solve() function performs
 *      Transport(dt/2) -> Reaction(dt) -> Transport(dt/2) with per-species
 *      Picard iterations, MPI halo exchange, and detailed performance timers.
 *
 * Key design choices:
 *   - All performance-critical loops use Kokkos parallel_for / parallel_reduce
 *     for portability across GPU (CUDA/HIP) and multi-core CPU backends.
 *   - The 27-point stencil supports full 3-D cross-dispersion on non-uniform
 *     vertical grids (variable dz).
 *   - MPI domain decomposition is handled via halo exchange after each major
 *     computational phase.
 *
 * @note Nitrogen-cycle species are indexed as:
 *       0 = NH4 (ammonium), 1 = NO3 (nitrate), 2 = DO (dissolved oxygen),
 *       3 = DOC (dissolved organic carbon), 4 = DON (dissolved organic nitrogen).
 *       NH4 and NO3 use molecular mass; DON and DOC use elemental mass.
 *
 * @copyright Part of the SERGHEI framework.
 */

#ifndef _RT_FUNCTIONGW_H_
#define _RT_FUNCTIONGW_H_

#include "const.h"
#include <sstream>
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
#include "MPIUtils.h"
#include "GwState.h"
#include "RTMatrix.h"
#include "RTSolverGW.h"
#include "RTStateGW.h"
#include "State.h"
#include "RTIntegrator.h"
#include <set>
#include <limits>
#include <math.h>
#include <cmath>

/**
 * @class RTFunctionGW
 * @brief Collects all functions for groundwater reactive transport computation.
 *
 * The class contains both private helper methods (flux limiters, gradient
 * ratio, flux correction) and public methods that drive the full transport
 * and reaction workflow.
 */
class RTFunctionGW
{

private:
	Kokkos::Timer timer;  ///< General-purpose Kokkos timer (reserved)
	Kokkos::Timer timer2; ///< Secondary Kokkos timer (reserved)
	realArr2 orig_RTcoef; ///< Unmodified physical stencil and RHS before RT boundary constraints.

	/// [P1-4] Interior row indices (idom) of the boundary cells whose RTcoef
	/// must be backed up into orig_RTcoef.  Only the residual-based branch of
	/// compute_boundary_mass_fluxes() reads orig_RTcoef, and only for these
	/// rows (columns 0..25), so the backup can skip all interior rows.
	intArr bres_rows;
	int nbres_rows = 0;

	/// [P2-4] One-time pure-geometry reciprocals (dx, dy and the dz grid
	/// never change during a simulation).  The face-coefficient precompute
	/// and the assembly kernel multiply these instead of dividing (~1 ulp
	/// difference per replaced term, L2).
	realArr inv_dz;       ///< [nCellMem] 1/dz(iGlob)
	realArr inv_dz_denom; ///< [nCellMem] 1 / (0.5*(dz(i)+dz(i+z))) with the same clamp logic as the old per-species precompute
	realArr inv_2dxdz;    ///< [nCellMem] 1/(2 dx dz(iGlob))
	realArr inv_2dydz;    ///< [nCellMem] 1/(2 dy dz(iGlob))
	realArr inv_dxdz_span; ///< [nCellMem] 1/(dx dz_span(iGlob))
	realArr inv_dydz_span; ///< [nCellMem] 1/(dy dz_span(iGlob))
	real inv_dx = 0.0, inv_dy = 0.0;     ///< 1/dx, 1/dy
	real inv_dxdx = 0.0, inv_dydy = 0.0; ///< 1/dx^2, 1/dy^2
	real inv_2dxdy = 0.0;                ///< 1/(2 dx dy)

	/// [P2-5] Per-species staging buffer for the 19 Picard-invariant stencil
	/// columns (all RTcoef columns except 0 and 7), written once per species
	/// per half-step by RTprepare_species_coefficients() and copied by every
	/// assembly iteration.  Columns 0 (sorption mass diagonal) and 7 (mass
	/// RHS + TVD correction) depend on rt.c(..,1) and stay in the assembly
	/// kernel.  A single buffer is sufficient because the species loop fully
	/// consumes it (all Picard iterations) before the next species overwrites it.
	realArr2 RTcoef_inv;

	/* =====================================================================
	 *  Section: Flux Limiter Functions (TVD High-Resolution Scheme)
	 *
	 *  These limiters are applied to the ratio of consecutive gradients r
	 *  to suppress spurious oscillations near sharp concentration fronts
	 *  while maintaining higher-than-first-order accuracy in smooth regions.
	 * ===================================================================== */

	/**
	 * @brief Minmod flux limiter.
	 *
	 * The most diffusive (most stable) limiter.  Selects the smaller of
	 * unity and the gradient ratio, ensuring TVD compliance.
	 *
	 *   psi(r) = max( 0, min(1, r) )
	 *
	 * @param r   Gradient ratio: (C_upstream - C_upstream_of_upstream) /
	 *            (C_downstream - C_upstream).  Dimensionless.
	 * @return    Limiter value in [0, 1].  Dimensionless.
	 */
	KOKKOS_INLINE_FUNCTION
	real minmod_limiter(real r) const
	{
		if (r <= 0.0)
			return 0.0;
		if (r >= 1.0)
			return 1.0;
		return r;
	}

	/**
	 * @brief Superbee flux limiter.
	 *
	 * A compressive limiter that can sharpen fronts but may produce
	 * slight over-shoots in highly advective regimes.
	 *
	 *   psi(r) = max( 0, min(1, 2r), min(2, r) )
	 *
	 * @param r   Gradient ratio (dimensionless).
	 * @return    Limiter value in [0, 2].  Dimensionless.
	 */
	KOKKOS_INLINE_FUNCTION
	real superbee_limiter(real r) const
	{
		if (r <= 0.0)
			return 0.0;
		real val1 = (r <= 0.5) ? 2.0 * r : 1.0;
		real val2 = (r <= 2.0) ? r : 2.0;
		return fmax(0.0, fmin(1.0, val1));
	}

	/**
	 * @brief Van Leer flux limiter.
	 *
	 * A smooth (differentiable) limiter that offers a good balance between
	 * compressiveness and stability.
	 *
	 *   psi(r) = (r + |r|) / (1 + |r|)
	 *
	 * @param r   Gradient ratio (dimensionless).
	 * @return    Limiter value in [0, 2].  Dimensionless.
	 */
	KOKKOS_INLINE_FUNCTION
	real vanleer_limiter(real r) const
	{
		if (r <= 0.0)
			return 0.0;
		return (r + fabs(r)) / (1.0 + fabs(r));
	}

	/**
	 * @brief Compute the gradient ratio r for flux limiters.
	 *
	 * The ratio compares the upwind gradient to the local gradient:
	 *   r = (C_upstream - C_upstream_of_upstream) /
	 *       (C_downstream - C_upstream)
	 *
	 * A sign-preserving tiny value is used for the denominator when it
	 * approaches zero.  This prevents flipping the sign of negative
	 * gradients, which would inject erroneous flux and cause the solution
	 * to blow up.
	 *
	 * @param C_upstream              Concentration at the upstream cell [mg/L].
	 * @param C_downstream            Concentration at the downstream cell [mg/L].
	 * @param C_upstream_upstream     Concentration one cell further upstream [mg/L].
	 * @return                        Gradient ratio r (dimensionless).
	 */
	KOKKOS_INLINE_FUNCTION
	real compute_gradient_ratio(real C_upstream, real C_downstream, real C_upstream_upstream) const
	{
		real denominator = C_downstream - C_upstream;
		// Preserve the sign of the tiny denominator!  Otherwise a negative
		// gradient would be forcibly flipped, injecting wrong flux and
		// causing numerical blow-up.
		if (fabs(denominator) < 1.0e-20)
		{
			denominator = (denominator >= 0.0) ? 1.0e-20 : -1.0e-20;
		}
		return (C_upstream - C_upstream_upstream) / denominator;
	}

	/**
	 * @brief Compute the high-order flux correction (anti-diffusive flux).
	 *
	 * The correction is the difference between the second-order central
	 * flux and the first-order upwind flux, modulated by the chosen
	 * flux limiter:
	 *
	 *   F_correction = psi(r) * (F_high - F_low)
	 *
	 * where:
	 *   F_low  = v * C_upstream              (first-order upwind)
	 *   F_high = v * 0.5 * (C_up + C_down)   (second-order central)
	 *
	 * This correction is added explicitly to the RHS during Picard
	 * iterations to upgrade the baseline first-order scheme to a
	 * high-resolution TVD scheme.
	 *
	 * @param C_upstream              Upstream concentration [mg/L].
	 * @param C_downstream            Downstream concentration [mg/L].
	 * @param C_upstream_upstream     Second upstream concentration [mg/L].
	 * @param velocity                Inter-cell velocity [m/s].
	 * @param dt                      Time step [s].
	 * @param dx                      Cell spacing [m].
	 * @param limiter_type            1 = Minmod, 2 = Superbee, 3 = Van Leer.
	 * @return                        Flux correction [mg/(L*s)].
	 */
	KOKKOS_INLINE_FUNCTION
	real compute_flux_correction(real C_upstream, real C_downstream, real C_upstream_upstream,
								 real velocity, real dt, real dx, int limiter_type) const
	{
		// First-order upwind flux [mg/(L*s)]
		real flux_low = velocity * C_upstream;

		// Second-order central flux (high-order) [mg/(L*s)]
		real flux_high = velocity * 0.5 * (C_upstream + C_downstream);

		// Compute the gradient ratio
		real r = compute_gradient_ratio(C_upstream, C_downstream, C_upstream_upstream);

		// Apply the selected flux limiter
		real psi = 0.0;
		switch (limiter_type)
		{
		case 1: // Minmod (most stable)
			psi = minmod_limiter(r);
			break;
		case 2: // Superbee (most compressive)
			psi = superbee_limiter(r);
			break;
		case 3: // Van Leer (balanced)
			psi = vanleer_limiter(r);
			break;
		default:
			psi = minmod_limiter(r);
			break;
		}

		// Return the anti-diffusive flux correction
		return psi * (flux_high - flux_low);
	}

	/**
	 * @brief Equilibrium sorption concentration S(C) for one species.
	 *
	 * Linear:      S = Kd * C
	 * Freundlich:  S = Kf * C^Nf
	 * Langmuir:    S = (Alpha_D * Beta_D * C) / (1 + Alpha_D * C)
	 */
	KOKKOS_INLINE_FUNCTION
	real sorption_equilibrium_concentration(const RTStateGW &rt, int iSpec, real C) const
	{
		int sorption_type = rt.Sorption_Type(iSpec);
		if (sorption_type == 1)
		{
			return rt.Kd(iSpec) * C;
		}
		if (sorption_type == 2)
		{
			if (C <= 0.0)
				return 0.0;
			return rt.Kf(iSpec) * pow(C, rt.Nf(iSpec));
		}
		if (sorption_type == 3)
		{
			real alpha = rt.Alpha_D(iSpec);
			real beta_L = rt.Beta_D(iSpec);
			return (alpha * beta_L * C) / (1.0 + alpha * C);
		}
		return 0.0;
	}

	/**
	 * @brief Tangent slope dS/dC for the equilibrium sorption isotherm.
	 */
	KOKKOS_INLINE_FUNCTION
	real sorption_equilibrium_derivative(const RTStateGW &rt, int iSpec, real C) const
	{
		int sorption_type = rt.Sorption_Type(iSpec);
		if (sorption_type == 1)
		{
			return rt.Kd(iSpec);
		}
		if (sorption_type == 2)
		{
			if (C <= 0.0)
				return 0.0;
			real Nf = rt.Nf(iSpec);
			return rt.Kf(iSpec) * Nf * pow(C, Nf - 1.0);
		}
		if (sorption_type == 3)
		{
			real alpha = rt.Alpha_D(iSpec);
			real beta_L = rt.Beta_D(iSpec);
			real denom = 1.0 + alpha * C;
			return (alpha * beta_L) / (denom * denom);
		}
		return 0.0;
	}

public:
	/* -----------------------------------------------------------------
	 *  [P3-4] Picard iteration statistics (always-on counters, reported
	 *  to log.out together with the linear-solver statistics block).
	 * ----------------------------------------------------------------- */
	std::vector<long> rtPicardEpisodes; ///< [n_mass] Picard episodes entered (species half-steps)
	std::vector<long> rtPicardIters;    ///< [n_mass] total Picard iterations executed
	std::vector<long> rtPicardConv1;    ///< [n_mass] episodes converged at iteration 1
	std::vector<long> rtPicardConv2;    ///< [n_mass] episodes converged at iteration 2
	std::vector<long> rtPicardConv3;    ///< [n_mass] episodes converged at iteration 3

	inline void rtPicardEnsure(int n_mass)
	{
		if ((int)rtPicardIters.size() < n_mass)
		{
			rtPicardEpisodes.assign(n_mass, 0);
			rtPicardIters.assign(n_mass, 0);
			rtPicardConv1.assign(n_mass, 0);
			rtPicardConv2.assign(n_mass, 0);
			rtPicardConv3.assign(n_mass, 0);
		}
	}

	/// Episodes that exhausted picard_iter_max without meeting the
	/// tolerance (derived: episodes minus the converged buckets).
	inline long rtPicardMaxed(int iSpec) const
	{
		return rtPicardEpisodes[iSpec] - rtPicardConv1[iSpec] - rtPicardConv2[iSpec] - rtPicardConv3[iSpec];
	}

	std::string picardStatsReport()
	{
		std::ostringstream os;
		os << "-----RTGW PICARD STATISTICS-----" << std::endl;
		long ep = 0, it = 0;
		for (size_t sp = 0; sp < rtPicardIters.size(); ++sp) { ep += rtPicardEpisodes[sp]; it += rtPicardIters[sp]; }
		os << "rtgw.pic.episodes\t" << ep << "\titers\t" << it;
		if (ep > 0) os << "\tmean\t" << (double)it / (double)ep;
		os << std::endl;
		for (size_t sp = 0; sp < rtPicardIters.size(); ++sp)
		{
			if (rtPicardEpisodes[sp] == 0) continue;
			os << "rtgw.pic.spec" << sp
			   << "\tepisodes\t" << rtPicardEpisodes[sp]
			   << "\titers\t" << rtPicardIters[sp]
			   << "\tmean\t" << (double)rtPicardIters[sp] / (double)rtPicardEpisodes[sp]
			   << "\tconv@1\t" << rtPicardConv1[sp]
			   << "\tconv@2\t" << rtPicardConv2[sp]
			   << "\tconv@3\t" << rtPicardConv3[sp]
			   << "\tmaxed\t" << rtPicardMaxed(sp) << std::endl;
		}
		return os.str();
	}

	/* =====================================================================
	 *  Section: Dispersion Tensor Computation
	 *
	 *  Computes the physical seepage velocity vector and the full 3-D
	 *  hydrodynamic dispersion tensor at cell interfaces (i+1/2, j+1/2,
	 *  k+1/2).  Supports multiple chemical species with species-specific
	 *  longitudinal / transverse dispersivities and molecular diffusion
	 *  coefficients.
	 * ===================================================================== */

	/**
	 * @brief Compute physical velocity vectors and the hydrodynamic
	 *        dispersion tensor for all species.
	 *
	 * The method performs three stages:
	 *   1. **Velocity transformation** — converts Darcy fluxes from the
	 *      native coordinate convention (q(0): positive leftward, q(1):
	 *      positive inward, q(2): positive upward) to the grid convention
	 *      (i: rightward, j: outward, k: downward).  The result is stored
	 *      in rt.aveVB(:,0..2) with the speed magnitude in (:,3).
	 *   2. **Inter-cell velocity interpolation** — averages the velocity
	 *      components at each cell face (i+1/2, j+1/2, k+1/2) and stores
	 *      them in rt.q_dispersion(:,0..11).
	 *   3. **Dispersion tensor assembly** — for each species, computes the
	 *      9 components of the symmetric dispersion tensor (Dxx, Dyy, Dzz,
	 *      Dxy, Dxz, Dyx, Dyz, Dzx, Dzy) from the longitudinal and
	 *      transverse dispersivities and the face-centred velocities.
	 *      Molecular diffusion (weighted by water content) is added to the
	 *      diagonal.  NaN values arising from zero-velocity faces are
	 *      replaced with zero.  Dimensionality is enforced (1-D, 2-D, or
	 *      3-D) by zeroing inappropriate components.
	 *
	 * @param rt    Reactive transport state (velocity, dispersivities,
	 *              diffusion coefficients, output tensor).
	 * @param gw    Groundwater state (Darcy fluxes, water content).
	 * @param gdom  Groundwater domain (grid dimensions, spacing).
	 * @param gmpi  MPI communication handler (halo exchange).
	 * @param par   Parallel execution context (rank info).
	 */
	inline void dispersion_tensor(RTStateGW &rt, GwState &gw, GwDomain &gdom, GwMPI &gmpi, Parallel &par)
	{
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering dispersion_tensor()" << std::endl;
#endif
		// -----------------------------------------------------------------
		// Stage 1: Velocity vector coordinate transformation
		//
		// Native convention: q(0) positive leftward,  q(1) positive inward,
		//                    q(2) positive upward.
		// Grid convention:   i increases rightward, j increases outward,
		//                    k increases downward.
		// Therefore:  v_x = -q(0),  v_y = -q(1),  v_z = -q(2)
		// -----------------------------------------------------------------
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			rt.aveVB(iGlob, 0) = -(gw.q_new(iGlob, 0)); // vx [m/s]
			rt.aveVB(iGlob, 1) = -(gw.q_new(iGlob, 1)); // vy [m/s]
			rt.aveVB(iGlob, 2) = -(gw.q_new(iGlob, 2)); // vz [m/s]

			// Velocity magnitude |v| [m/s]
			// [P1-6] x*x instead of pow(x,2) (same value, cheaper on GPU)
			const real vx = rt.aveVB(iGlob, 0);
			const real vy = rt.aveVB(iGlob, 1);
			const real vz = rt.aveVB(iGlob, 2);
			rt.aveVB(iGlob, 3) = sqrt(vx * vx + vy * vy + vz * vz); });
		Kokkos::fence();

		// [Critical Fix 1]: Exchange velocity vectors across MPI boundaries.
		// Without this, halo cells have zero velocity, which completely
		// blocks advection across domain boundaries.
		gmpi.mpi_sendrecv(rt.aveVB, gdom, par);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> dispersion_tensor(): velocity transform & MPI exchange done" << std::endl;
#endif

		// -----------------------------------------------------------------
		// Stage 2: Inter-cell velocity interpolation
		//
		// Compute velocity components at each of the three face types
		// (i+1/2, j+1/2, k+1/2) using four-point bilinear averaging
		// of the cell-centred values from the surrounding cells.
		// -----------------------------------------------------------------
		Kokkos::parallel_for("Calculation_dispersion_coefficient", gdom.nCellMem, KOKKOS_LAMBDA(int idom) {

            int iGlob = idom;
            // Decode (i,j,k) indices for safe neighbour velocity access
            int kk = idom / (gdom.nxhc * gdom.nyhc);
            int rem = idom % (gdom.nxhc * gdom.nyhc);
            int jj = rem / gdom.nxhc;
            int ii = rem % gdom.nxhc;

            // [Critical Fix 2]: Lambda with clamped index access, replacing the
            // blunt "if(ii<1) return" approach so that boundary cells still
            // transition smoothly without early-exit artefacts.
            auto V = [&](int c_k, int c_j, int c_i, int comp) {
                if (c_k < 0) c_k = 0; else if (c_k >= gdom.nzhc) c_k = gdom.nzhc - 1;
                if (c_j < 0) c_j = 0; else if (c_j >= gdom.nyhc) c_j = gdom.nyhc - 1;
                if (c_i < 0) c_i = 0; else if (c_i >= gdom.nxhc) c_i = gdom.nxhc - 1;
                int safe_idx = c_k * gdom.nxhc * gdom.nyhc + c_j * gdom.nxhc + c_i;
                return rt.aveVB(safe_idx, comp);
            };

            auto Dz = [&](int c_k, int c_j, int c_i) {
                if (c_k < 0) c_k = 0; else if (c_k >= gdom.nzhc) c_k = gdom.nzhc - 1;
                if (c_j < 0) c_j = 0; else if (c_j >= gdom.nyhc) c_j = gdom.nyhc - 1;
                if (c_i < 0) c_i = 0; else if (c_i >= gdom.nxhc) c_i = gdom.nxhc - 1;
                int safe_idx = c_k * gdom.nxhc * gdom.nyhc + c_j * gdom.nxhc + c_i;
                return gdom.dz(safe_idx);
            };

			// --- Face i+1/2 velocity components ---
			// vx at i+1/2 (direct cell-centre value) [m/s]
			rt.q_dispersion(iGlob, 0) = V(kk, jj, ii, 0);
			// vy at i+1/2 (four-point average) [m/s]
			rt.q_dispersion(iGlob, 1) = 0.5 * ((V(kk, jj - 1, ii, 1) + V(kk, jj, ii, 1)) * 0.5 + (V(kk, jj - 1, ii + 1, 1) + V(kk, jj, ii + 1, 1)) * 0.5);
			// vz at i+1/2 (four-point average) [m/s]
			rt.q_dispersion(iGlob, 2) = 0.5 * ((V(kk - 1, jj, ii, 2) + V(kk, jj, ii, 2)) * 0.5 + (V(kk - 1, jj, ii + 1, 2) + V(kk, jj, ii + 1, 2)) * 0.5);
			// Speed at i+1/2 [m/s]
			// [P1-6] speed at face i+1/2 with cached products (x*x, not pow)
			// [P2-1] Species-invariant face geometry: the velocity products,
			// the speed and its reciprocal are staged here once per half-step
			// so that Stage 3 no longer recomputes them for every species.
			const real fx0 = rt.q_dispersion(iGlob, 0), fx1 = rt.q_dispersion(iGlob, 1), fx2 = rt.q_dispersion(iGlob, 2);
			const real spd_fx = sqrt(fx0 * fx0 + fx1 * fx1 + fx2 * fx2);
			rt.q_dispersion(iGlob, 3) = spd_fx;
			rt.q_face_geo(iGlob, 0) = fx0 * fx0;
			rt.q_face_geo(iGlob, 1) = fx1 * fx1;
			rt.q_face_geo(iGlob, 2) = fx2 * fx2;
			rt.q_face_geo(iGlob, 3) = fx0 * fx1;
			rt.q_face_geo(iGlob, 4) = fx0 * fx2;
			rt.q_face_geo(iGlob, 5) = fx1 * fx2;
			rt.q_face_geo(iGlob, 6) = spd_fx;
			rt.q_face_geo(iGlob, 7) = (spd_fx > 0.0) ? 1.0 / spd_fx : 0.0;

			// --- Face j+1/2 velocity components ---
			// vx at j+1/2 (four-point average) [m/s]
			//! rt.q_dispersion(iGlob, 4) = 0.5 * ((rt.aveVB(iGlob, 1) + rt.aveVB(iGlob - 1, 1)) * 0.5 + (rt.aveVB(iGlob + gdom.nxhc, 1) + rt.aveVB(iGlob + gdom.nxhc - 1, 1)) * 0.5);
			rt.q_dispersion(iGlob, 4) = 0.5 * ((V(kk, jj, ii, 0) + V(kk, jj, ii - 1, 0)) * 0.5 + (V(kk, jj + 1, ii, 0) + V(kk, jj + 1, ii - 1, 0)) * 0.5);
			// vy at j+1/2 (direct cell-centre value) [m/s]
			rt.q_dispersion(iGlob, 5) = V(kk, jj, ii, 1);
			// vz at j+1/2 (four-point average) [m/s]
			rt.q_dispersion(iGlob, 6) = 0.5 * ((V(kk - 1, jj, ii, 2) + V(kk, jj, ii, 2)) * 0.5 + (V(kk - 1, jj + 1, ii, 2) + V(kk, jj + 1, ii, 2)) * 0.5);
			// Speed at j+1/2 [m/s]
			// [P2-1] Species-invariant face geometry (see face i+1/2 above).
			const real fy4 = rt.q_dispersion(iGlob, 4), fy5 = rt.q_dispersion(iGlob, 5), fy6 = rt.q_dispersion(iGlob, 6);
			const real spd_fy = sqrt(fy4 * fy4 + fy5 * fy5 + fy6 * fy6);
			rt.q_dispersion(iGlob, 7) = spd_fy;
			rt.q_face_geo(iGlob, 8) = fy4 * fy4;
			rt.q_face_geo(iGlob, 9) = fy5 * fy5;
			rt.q_face_geo(iGlob, 10) = fy6 * fy6;
			rt.q_face_geo(iGlob, 11) = fy4 * fy5;
			rt.q_face_geo(iGlob, 12) = fy4 * fy6;
			rt.q_face_geo(iGlob, 13) = fy5 * fy6;
			rt.q_face_geo(iGlob, 14) = spd_fy;
			rt.q_face_geo(iGlob, 15) = (spd_fy > 0.0) ? 1.0 / spd_fy : 0.0;

// --- Face k+1/2 velocity components ---
			// Dynamic calculation of interpolation weights based on vertical grid thickness
			real dz_k = Dz(kk, jj, ii);
			real dz_kp1 = Dz(kk + 1, jj, ii);
			real weight_k = dz_kp1 / (dz_k + dz_kp1);
			real weight_kp1 = dz_k / (dz_k + dz_kp1);

			// vx at k+1/2 (Performing non-uniform vertical interpolation on cell-centred VX) [m/s]
			real vx_k = 0.5 * (V(kk, jj, ii, 0) + V(kk, jj, ii - 1, 0));
			real vx_kp1 = 0.5 * (V(kk + 1, jj, ii, 0) + V(kk + 1, jj, ii - 1, 0));
			rt.q_dispersion(iGlob, 8) = weight_k * vx_k + weight_kp1 * vx_kp1;

			// vy at k+1/2 (Performing non-uniform vertical interpolation on cell-centred VY) [m/s]
			real vy_k = 0.5 * (V(kk, jj, ii, 1) + V(kk, jj - 1, ii, 1));
			real vy_kp1 = 0.5 * (V(kk + 1, jj, ii, 1) + V(kk + 1, jj - 1, ii, 1));
			rt.q_dispersion(iGlob, 9) = weight_k * vy_k + weight_kp1 * vy_kp1;
			rt.q_dispersion(iGlob, 10) = V(kk, jj, ii, 2);
			// Speed at k+1/2 [m/s]
			// [P2-1] Species-invariant face geometry (see face i+1/2 above).
			const real fz8 = rt.q_dispersion(iGlob, 8), fz9 = rt.q_dispersion(iGlob, 9), fz10 = rt.q_dispersion(iGlob, 10);
			const real spd_fz = sqrt(fz8 * fz8 + fz9 * fz9 + fz10 * fz10);
			rt.q_dispersion(iGlob, 11) = spd_fz;
			rt.q_face_geo(iGlob, 16) = fz8 * fz8;
			rt.q_face_geo(iGlob, 17) = fz9 * fz9;
			rt.q_face_geo(iGlob, 18) = fz10 * fz10;
			rt.q_face_geo(iGlob, 19) = fz8 * fz9;
			rt.q_face_geo(iGlob, 20) = fz8 * fz10;
			rt.q_face_geo(iGlob, 21) = fz9 * fz10;
			rt.q_face_geo(iGlob, 22) = spd_fz;
			rt.q_face_geo(iGlob, 23) = (spd_fz > 0.0) ? 1.0 / spd_fz : 0.0; });

		// -----------------------------------------------------------------
		// Stage 3: Dispersion tensor assembly (per species)
		//
		// Computes the 9 components of the hydrodynamic dispersion tensor
		// using the Scheidegger-Bear formulation:
		//   D_ij = alpha_T * |v| * delta_ij + (alpha_L - alpha_T) * v_i * v_j / |v|
		// then adds molecular diffusion (weighted by water content) to the
		// diagonal.
		// -----------------------------------------------------------------
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> dispersion_tensor(): computing dispersion tensor for " << rt.n_mass << " species" << std::endl;
#endif
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {

			// [P2-1] All face velocity products, speeds and reciprocal speeds
			// are species-invariant and were staged once per half-step by
			// Stage 2 into q_face_geo (per face f, base = 8*f:
			// [vx^2, vy^2, vz^2, vx*vy, vx*vz, vy*vz, |v|, 1/|v|]).
			// This kernel only forms the species-specific alpha_L / alpha_T
			// linear combinations.
			const real aL = rt.alpha_L(iSpec), aT = rt.alpha_T(iSpec);

			// --- Dxx, Dxy, Dxz at face i+1/2 [m^2/s] ---
			// [P2-1b] Division by |v| is replaced by multiplication with the
			// Stage-2 reciprocal (~1 ulp difference, L2).  Zero-velocity faces
			// (and pathological NaN/Inf speeds, where the reciprocal is 0)
			// take the explicit zero branch: this yields exactly the same
			// values as the former 0/0 = NaN -> isnan sweep, without ever
			// producing a NaN.
			const real inv_fx = rt.q_face_geo(iGlob, 7);
			if (inv_fx > 0.0)
			{
				rt.dcal(iSpec, iGlob, 0) = (aT * (rt.q_face_geo(iGlob, 1) + rt.q_face_geo(iGlob, 2)) + aL * rt.q_face_geo(iGlob, 0)) * inv_fx;
				rt.dcal(iSpec, iGlob, 3) = ((aL - aT) * rt.q_face_geo(iGlob, 3)) * inv_fx;
				rt.dcal(iSpec, iGlob, 4) = ((aL - aT) * rt.q_face_geo(iGlob, 4)) * inv_fx;
			}
			else
			{
				rt.dcal(iSpec, iGlob, 0) = 0.0;
				rt.dcal(iSpec, iGlob, 3) = 0.0;
				rt.dcal(iSpec, iGlob, 4) = 0.0;
			}

			//! --- Dyy, Dyx, Dyz at face j+1/2 [m^2/s] ---
			const real inv_fy = rt.q_face_geo(iGlob, 15);
			if (inv_fy > 0.0)
			{
				rt.dcal(iSpec, iGlob, 1) = (aT * (rt.q_face_geo(iGlob, 8) + rt.q_face_geo(iGlob, 10)) + aL * rt.q_face_geo(iGlob, 9)) * inv_fy;
				rt.dcal(iSpec, iGlob, 5) = ((aL - aT) * rt.q_face_geo(iGlob, 11)) * inv_fy;
				rt.dcal(iSpec, iGlob, 6) = ((aL - aT) * rt.q_face_geo(iGlob, 13)) * inv_fy;
			}
			else
			{
				rt.dcal(iSpec, iGlob, 1) = 0.0;
				rt.dcal(iSpec, iGlob, 5) = 0.0;
				rt.dcal(iSpec, iGlob, 6) = 0.0;
			}

			//! --- Dzz, Dzx, Dzy at face k+1/2 [m^2/s] ---
			const real inv_fz = rt.q_face_geo(iGlob, 23);
			if (inv_fz > 0.0)
			{
				rt.dcal(iSpec, iGlob, 2) = (aT * (rt.q_face_geo(iGlob, 16) + rt.q_face_geo(iGlob, 17)) + aL * rt.q_face_geo(iGlob, 18)) * inv_fz;
				rt.dcal(iSpec, iGlob, 7) = ((aL - aT) * rt.q_face_geo(iGlob, 20)) * inv_fz;
				rt.dcal(iSpec, iGlob, 8) = ((aL - aT) * rt.q_face_geo(iGlob, 21)) * inv_fz;
			}
			else
			{
				rt.dcal(iSpec, iGlob, 2) = 0.0;
				rt.dcal(iSpec, iGlob, 7) = 0.0;
				rt.dcal(iSpec, iGlob, 8) = 0.0;
			}

			// Add molecular diffusion (weighted by water content) to diagonal
			// D_total = D_hydrodynamic + theta * D_molecular  [m^2/s]
			rt.dcal(iSpec, iGlob, 0) += gw.wc(iGlob, 0) * rt.diffusion_molecular(iSpec);
			rt.dcal(iSpec, iGlob, 1) += gw.wc(iGlob, 0) * rt.diffusion_molecular(iSpec);
			rt.dcal(iSpec, iGlob, 2) += gw.wc(iGlob, 0) * rt.diffusion_molecular(iSpec);
			// Enforce dimensionality: zero out dispersion components that
			// correspond to inactive spatial dimensions.
			if (gdom.nx == 1) // 2-D YZ plane: qx = 0, so Dxx = Dxy = Dyx = Dxz = Dzx = 0
			{
				rt.dcal(iSpec, iGlob, 0) = 0;
				rt.dcal(iSpec, iGlob, 3) = 0;
				rt.dcal(iSpec, iGlob, 4) = 0;
				rt.dcal(iSpec, iGlob, 5) = 0;
				rt.dcal(iSpec, iGlob, 7) = 0;
			}
			else if (gdom.ny == 1) // 2-D XZ plane: qy = 0
			{
				rt.dcal(iSpec, iGlob, 1) = 0;
				rt.dcal(iSpec, iGlob, 3) = 0;
				rt.dcal(iSpec, iGlob, 5) = 0;
				rt.dcal(iSpec, iGlob, 6) = 0;
				rt.dcal(iSpec, iGlob, 8) = 0;
			}
			else if (gdom.nz == 1) // 2-D XY plane: qz = 0
			{
				rt.dcal(iSpec, iGlob, 2) = 0;
				rt.dcal(iSpec, iGlob, 4) = 0;
				rt.dcal(iSpec, iGlob, 6) = 0;
				rt.dcal(iSpec, iGlob, 7) = 0;
				rt.dcal(iSpec, iGlob, 8) = 0;
			} });
		}
	}

	/**
	 * @brief One-time precompute of pure-geometry reciprocals and wz.
	 *
	 * dx, dy and the dz grid never change during a simulation, so every
	 * geometric divisor used by the face-coefficient precompute
	 * ([P2-4]: replaces divisions with reciprocal multiplications, ~1 ulp,
	 * L2) and wz plus the c_difzz denominator (pure geometry, previously
	 * recomputed for every species; values bitwise identical) are computed
	 * once here.  Called lazily from rt_solve().
	 *
	 * @param gdom  Domain geometry (grid dimensions, spacing).
	 * @param rt    Reactive transport state (wz written in place).
	 */
	inline void RTprecompute_geometry(GwDomain &gdom, RTStateGW &rt)
	{
		inv_dx = 1.0 / gdom.dx;
		inv_dy = 1.0 / gdom.dy;
		inv_dxdx = 1.0 / (gdom.dx * gdom.dx);
		inv_dydy = 1.0 / (gdom.dy * gdom.dy);
		inv_2dxdy = 1.0 / (2.0 * gdom.dx * gdom.dy);

		inv_dz = realArr("inv_dz", gdom.nCellMem);
		inv_dz_denom = realArr("inv_dz_denom", gdom.nCellMem);
		inv_2dxdz = realArr("inv_2dxdz", gdom.nCellMem);
		inv_2dydz = realArr("inv_2dydz", gdom.nCellMem);
		inv_dxdz_span = realArr("inv_dxdz_span", gdom.nCellMem);
		inv_dydz_span = realArr("inv_dydz_span", gdom.nCellMem);

		const int plane = gdom.nxhc * gdom.nyhc;
		const int nCellMem_l = gdom.nCellMem;
		const real dx_l = gdom.dx, dy_l = gdom.dy;
		auto wz_l = rt.wz;
		auto inv_dz_l = inv_dz;
		auto inv_dz_denom_l = inv_dz_denom;
		auto inv_2dxdz_l = inv_2dxdz;
		auto inv_2dydz_l = inv_2dydz;
		auto inv_dxdz_span_l = inv_dxdz_span;
		auto inv_dydz_span_l = inv_dydz_span;

		Kokkos::parallel_for("RT_Geometry_Reciprocal", gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			// wz and the c_difzz denominator replicate the clamp logic of the
			// former per-species precompute kernel (identical values, now
			// computed once).
			const int iGlob_zp = (iGlob + plane < nCellMem_l) ? iGlob + plane : iGlob;
			const real dz_here = gdom.dz(iGlob);
			const real dz_zp = gdom.dz(iGlob_zp);
			const real dz_sum = dz_here + dz_zp;
			wz_l(iGlob) = (dz_sum > 0.0) ? dz_zp / dz_sum : 0.5;
			const real dz_denom = (dz_sum > 0.0) ? (0.5 * dz_sum) : 1.0;
			inv_dz_denom_l(iGlob) = 1.0 / dz_denom;
			inv_dz_l(iGlob) = 1.0 / dz_here;

			// dz_span reciprocals.  Only interior cells ever read these
			// entries; the guarded neighbour access below only affects
			// never-read halo slots.
			const int iGlob_zm = iGlob - plane;
			const int iGlob_zp_raw = iGlob + plane;
			const real dz_im = (iGlob_zm >= 0) ? gdom.dz(iGlob_zm) : dz_here;
			const real dz_ip = (iGlob_zp_raw < nCellMem_l) ? gdom.dz(iGlob_zp_raw) : dz_here;
			const real dz_span = 0.5 * dz_im + dz_here + 0.5 * dz_ip;
			inv_dxdz_span_l(iGlob) = 1.0 / (dx_l * dz_span);
			inv_dydz_span_l(iGlob) = 1.0 / (dy_l * dz_span);
			inv_2dxdz_l(iGlob) = 1.0 / (2.0 * dx_l * dz_here);
			inv_2dydz_l(iGlob) = 1.0 / (2.0 * dy_l * dz_here); });
		Kokkos::fence();
	}

	/* =====================================================================
	 *  Section: Picard-invariant precompute (performance, Level-1)
	 *
	 *  [P1-1]/[P1-2] The face coefficients and upwind weights below depend
	 *  only on quantities frozen for the whole transport half-step:
	 *  aveVB (written once per half-step by dispersion_tensor()),
	 *  dcal(iSpec,..) (once per species per half-step), rt.dt and the grid
	 *  geometry.  They used to be recomputed inside every Picard iteration
	 *  of RTlinear_system(); hoisting them out is bitwise identical.
	 * ===================================================================== */

	/**
	 * @brief Compute the Picard-invariant face coefficients and the invariant
	 *        stencil columns for one species.
	 *
	 * Every coefficient computed here depends only on the half-step-invariant
	 * face velocities (aveVB), the species dispersion tensor (dcal), the
	 * transport sub-step rt.dt, and the grid geometry.  Must be called once
	 * per species before the Picard loop (and after dispersion_tensor()).
	 *
	 * The face-coefficient loop covers the full halo memory (nCellMem)
	 * because the assembly kernel for interior cells reads negative-side
	 * face coefficients such as c_advxx(iGlob - 1) and c_difxx(iGlob - 1),
	 * which live on halo/face slots at physical X-/Y-/Z- faces.
	 *
	 * [P2-4] Divisors are pure geometry and have been inverted once by
	 * RTprecompute_geometry(); the expressions below multiply those
	 * reciprocals instead of dividing (~1 ulp per term, L2).  wz itself is
	 * pure geometry and is now written by RTprecompute_geometry() only.
	 *
	 * [P2-5] The second kernel additionally stages the 19 Picard-invariant
	 * stencil columns (RTcoef columns 1-6 and 8-26; the expressions are
	 * identical to the former assembly-kernel lines, term order preserved)
	 * into RTcoef_inv, so every assembly iteration only has to recompute the
	 * two Picard-dependent columns 0 (sorption mass diagonal) and 7 (mass
	 * RHS + TVD correction).
	 *
	 * @param gdom   Domain geometry (grid dimensions, spacing).
	 * @param rt     Reactive transport state (coefficients written in place).
	 * @param iSpec  Species index (0-based).
	 */
	inline void RTprepare_species_coefficients(GwDomain &gdom, RTStateGW &rt, int iSpec)
	{
		if (RTcoef_inv.extent(0) != (size_t)gdom.nCell || RTcoef_inv.extent(1) != 27)
		{
			RTcoef_inv = realArr2("RTcoef_inv", gdom.nCell, 27);
		}
		auto coef_inv = RTcoef_inv;
		const real inv_dx_l = inv_dx, inv_dy_l = inv_dy;
		const real inv_dxdx_l = inv_dxdx, inv_dydy_l = inv_dydy, inv_2dxdy_l = inv_2dxdy;
		auto inv_dz_denom_l = inv_dz_denom;
		auto inv_dz_l = inv_dz;
		auto inv_2dxdz_l = inv_2dxdz;
		auto inv_2dydz_l = inv_2dydz;
		auto inv_dxdz_span_l = inv_dxdz_span;
		auto inv_dydz_span_l = inv_dydz_span;

		Kokkos::parallel_for("Precompute_RT_Face_Coefficients", gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			rt.c_advxx(iGlob) = rt.aveVB(iGlob, 0) * rt.dt * inv_dx_l;
			rt.c_advyy(iGlob) = rt.aveVB(iGlob, 1) * rt.dt * inv_dy_l;
			rt.c_advzz(iGlob) = rt.aveVB(iGlob, 2) * rt.dt;

			rt.c_difxx(iGlob) = rt.dcal(iSpec, iGlob, 0) * rt.dt * inv_dxdx_l;
			rt.c_difyy(iGlob) = rt.dcal(iSpec, iGlob, 1) * rt.dt * inv_dydy_l;
			rt.c_difzz(iGlob) = rt.dcal(iSpec, iGlob, 2) * rt.dt * inv_dz_denom_l(iGlob);

			rt.c_difxy(iGlob) = rt.dcal(iSpec, iGlob, 3) * rt.dt * inv_2dxdy_l;
			rt.c_difxz(iGlob) = rt.dcal(iSpec, iGlob, 4) * rt.dt;
			rt.c_difyx(iGlob) = rt.dcal(iSpec, iGlob, 5) * rt.dt * inv_2dxdy_l;
			rt.c_difyz(iGlob) = rt.dcal(iSpec, iGlob, 6) * rt.dt;
			rt.c_difzx(iGlob) = rt.dcal(iSpec, iGlob, 7) * rt.dt;
			rt.c_difzy(iGlob) = rt.dcal(iSpec, iGlob, 8) * rt.dt; });

		// [P2-5] Picard-invariant stencil columns.  Every expression is
		// verbatim from the former Build_Matrix_Fused section (the divisions
		// by geometric products are the [P2-4] reciprocal-multiply form).
		// Runs after the face-coefficient kernel above on the same execution
		// space, so the c_* / Up_Weighting arrays are already final.
		Kokkos::parallel_for("Precompute_Invariant_Coef", gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

			const int plane = gdom.nxhc * gdom.nyhc;
			const real inv_dz_i = inv_dz_l(iGlob);
			const real inv_2dxdz_i = inv_2dxdz_l(iGlob);
			const real inv_2dydz_i = inv_2dydz_l(iGlob);
			const real inv_dxdz_span_i = inv_dxdz_span_l(iGlob);
			const real inv_dydz_span_i = inv_dydz_span_l(iGlob);

			real local_coef[27] = {0.0};

			// Face neighbour coefficients (moved to LHS, note signs)
			// local_coef[1]: C_{i+1}
			local_coef[1] = -(rt.c_difxx(iGlob) - (1.0 - rt.Up_Weighting_xp(iGlob)) * rt.c_advxx(iGlob));
			// local_coef[2]: C_{i-1}
			local_coef[2] = -(rt.c_difxx(iGlob - 1) + rt.Up_Weighting_xm(iGlob) * rt.c_advxx(iGlob - 1));
			// local_coef[3]: C_{j+1}
			local_coef[3] = -(rt.c_difyy(iGlob) - (1.0 - rt.Up_Weighting_yp(iGlob)) * rt.c_advyy(iGlob));
			// local_coef[4]: C_{j-1}
			local_coef[4] = -(rt.c_difyy(iGlob - gdom.nxhc) + rt.Up_Weighting_ym(iGlob) * rt.c_advyy(iGlob - gdom.nxhc));
			// local_coef[5]: C_{k+1}
			local_coef[5] = -(rt.c_difzz(iGlob) * inv_dz_i - (1.0 - rt.Up_Weighting_zp(iGlob)) * rt.c_advzz(iGlob) * inv_dz_i);
			// local_coef[6]: C_{k-1}
			local_coef[6] = -(rt.c_difzz(iGlob - plane) * inv_dz_i + rt.Up_Weighting_zm(iGlob) * rt.c_advzz(iGlob - plane) * inv_dz_i);

			// =================================================================
			// Cross-dispersion stencil entries (edge and corner neighbours of
			// the 27-point stencil, from the cross-derivative terms).
			// C_{i,j+1}: accumulated from XY and ZY cross-terms
			local_coef[8] = -0.5 * rt.c_difxy(iGlob) + 0.5 * rt.c_difxy(iGlob - 1) - rt.wz(iGlob) * rt.c_difzy(iGlob) * inv_2dydz_i + (1 - rt.wz(iGlob - plane)) * rt.c_difzy(iGlob - plane) * inv_2dydz_i;
			local_coef[3] += local_coef[8];
			// C_{i+1,j+1}
			local_coef[9] = -0.5 * rt.c_difxy(iGlob) - 0.5 * rt.c_difyx(iGlob);
			// Correction to C_{j-1} from XY cross-dispersion
			local_coef[10] = 0.5 * rt.c_difxy(iGlob) - 0.5 * rt.c_difxy(iGlob - 1) + rt.wz(iGlob) * rt.c_difzy(iGlob) * inv_2dydz_i - (1 - rt.wz(iGlob - plane)) * rt.c_difzy(iGlob - plane) * inv_2dydz_i;
			local_coef[4] += local_coef[10];
			// C_{i+1,j-1}
			local_coef[11] = 0.5 * rt.c_difxy(iGlob) + 0.5 * rt.c_difyx(iGlob - gdom.nxhc);
			// C_{i-1,j+1}
			local_coef[12] = 0.5 * rt.c_difxy(iGlob - 1) + 0.5 * rt.c_difyx(iGlob);
			// C_{i-1,j-1}
			local_coef[13] = -0.5 * rt.c_difxy(iGlob - 1) - 0.5 * rt.c_difyx(iGlob - gdom.nxhc);

			// --- Z-direction cross-dispersion contributions ---
			// Correction to C_{k+1} from XZ and YZ cross-terms
			local_coef[14] = -0.5 * rt.c_difxz(iGlob) * inv_dxdz_span_i + 0.5 * rt.c_difxz(iGlob - 1) * inv_dxdz_span_i - 0.5 * rt.c_difyz(iGlob) * inv_dydz_span_i + 0.5 * rt.c_difyz(iGlob - gdom.nxhc) * inv_dydz_span_i;
			local_coef[5] += local_coef[14];
			//! C_{i+1,k+1}
			local_coef[15] = -0.5 * rt.c_difxz(iGlob) * inv_dxdz_span_i - (1 - rt.wz(iGlob)) * rt.c_difzx(iGlob) * inv_2dxdz_i;
			// Correction to C_{k-1} from XZ and YZ cross-terms
			local_coef[16] = 0.5 * rt.c_difxz(iGlob) * inv_dxdz_span_i - 0.5 * rt.c_difxz(iGlob - 1) * inv_dxdz_span_i + 0.5 * rt.c_difyz(iGlob) * inv_dydz_span_i - 0.5 * rt.c_difyz(iGlob - gdom.nxhc) * inv_dydz_span_i;
			local_coef[6] += local_coef[16];
			//! C_{i+1,k-1}
			local_coef[17] = 0.5 * rt.c_difxz(iGlob) * inv_dxdz_span_i + rt.wz(iGlob - plane) * rt.c_difzx(iGlob - plane) * inv_2dxdz_i;
			// C_{i-1,k+1}
			local_coef[18] = 0.5 * rt.c_difxz(iGlob - 1) * inv_dxdz_span_i + (1 - rt.wz(iGlob)) * rt.c_difzx(iGlob) * inv_2dxdz_i;
			//! C_{i-1,k-1}
			local_coef[19] = -0.5 * rt.c_difxz(iGlob - 1) * inv_dxdz_span_i - rt.wz(iGlob - plane) * rt.c_difzx(iGlob - plane) * inv_2dxdz_i;
			// Correction to C_{i+1} from YX and ZX cross-terms
			local_coef[20] = -0.5 * rt.c_difyx(iGlob) + 0.5 * rt.c_difyx(iGlob - gdom.nxhc) - rt.wz(iGlob) * rt.c_difzx(iGlob) * inv_2dxdz_i + (1 - rt.wz(iGlob - plane)) * rt.c_difzx(iGlob - plane) * inv_2dxdz_i;
			local_coef[1] += local_coef[20];
			// Correction to C_{i-1} from YX and ZX cross-terms
			local_coef[21] = 0.5 * rt.c_difyx(iGlob) - 0.5 * rt.c_difyx(iGlob - gdom.nxhc) + rt.wz(iGlob) * rt.c_difzx(iGlob) * inv_2dxdz_i - (1 - rt.wz(iGlob - plane)) * rt.c_difzx(iGlob - plane) * inv_2dxdz_i;
			local_coef[2] += local_coef[21];
			// C_{j+1,k+1}
			local_coef[22] = -0.5 * rt.c_difyz(iGlob) * inv_dydz_span_i - (1 - rt.wz(iGlob)) * rt.c_difzy(iGlob) * inv_2dydz_i;
			// C_{j+1,k-1}
			local_coef[23] = 0.5 * rt.c_difyz(iGlob) * inv_dydz_span_i + rt.wz(iGlob - plane) * rt.c_difzy(iGlob - plane) * inv_2dydz_i;
			// C_{j-1,k+1}
			local_coef[24] = 0.5 * rt.c_difyz(iGlob - gdom.nxhc) * inv_dydz_span_i + (1 - rt.wz(iGlob)) * rt.c_difzy(iGlob) * inv_2dydz_i;
			// C_{j-1,k-1}
			local_coef[25] = -0.5 * rt.c_difyz(iGlob - gdom.nxhc) * inv_dydz_span_i - rt.wz(iGlob - plane) * rt.c_difzy(iGlob - plane) * inv_2dydz_i;

			// Stage everything except the Picard-dependent columns 0 and 7
			// (column 26 is structurally zero).
			for (int k = 0; k < 27; ++k)
				if (k != 0 && k != 7)
					coef_inv(idom, k) = local_coef[k]; });

		Kokkos::fence();
	}

	/**
	 * @brief Compute the upwind weighting factors for one transport half-step.
	 *
	 * The weights depend only on the sign of the half-step-invariant face
	 * velocities (aveVB): they are neither species- nor Picard-dependent.
	 * Called once per half-step right after dispersion_tensor(); the
	 * assembly kernel and the RT boundary conditions (applyRTMatBC) read
	 * exactly the same values as before the hoisting.
	 *
	 * @param gdom  Domain geometry (grid dimensions).
	 * @param rt    Reactive transport state (weights written in place).
	 */
	inline void RTcompute_upwind_weights(GwDomain &gdom, RTStateGW &rt)
	{
		Kokkos::parallel_for("RT_Upwind_Weights", gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

			// --- Upwind weighting factors (based on corrected velocity) ---
			// v > 0 (flow in positive direction): weight = 1.0 (use C_i)
			// v < 0 (flow in negative direction): weight = 0.0 (use C_{i+1})
			// v = 0: weight = 0.5 (central average)
			rt.Up_Weighting_xp(iGlob) = (rt.aveVB(iGlob, 0) > 0) ? rt.Up_Weighting_vplus : (rt.aveVB(iGlob, 0) < 0) ? rt.Up_Weighting_vminus
																													: 0.5;
			rt.Up_Weighting_xm(iGlob) = (rt.aveVB(iGlob - 1, 0) > 0) ? rt.Up_Weighting_vplus : (rt.aveVB(iGlob - 1, 0) < 0) ? rt.Up_Weighting_vminus
																																: 0.5;
			rt.Up_Weighting_yp(iGlob) = (rt.aveVB(iGlob, 1) > 0) ? rt.Up_Weighting_vplus : (rt.aveVB(iGlob, 1) < 0) ? rt.Up_Weighting_vminus
																													: 0.5;
			rt.Up_Weighting_ym(iGlob) = (rt.aveVB(iGlob - gdom.nxhc, 1) > 0) ? rt.Up_Weighting_vplus : (rt.aveVB(iGlob - gdom.nxhc, 1) < 0) ? rt.Up_Weighting_vminus
																																				: 0.5;
			rt.Up_Weighting_zp(iGlob) = (rt.aveVB(iGlob, 2) > 0) ? rt.Up_Weighting_vplus : (rt.aveVB(iGlob, 2) < 0) ? rt.Up_Weighting_vminus
																													: 0.5;
			rt.Up_Weighting_zm(iGlob) = (rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) > 0) ? rt.Up_Weighting_vplus : (rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) < 0) ? rt.Up_Weighting_vminus
																																											: 0.5; });
	}

	/* =====================================================================
	 *  Section: Linear System Assembly
	 *
	 *  Assembles the fully-implicit advection-dispersion equation for a
	 *  single chemical species on the 3-D structured grid.  The resulting
	 *  27-point stencil is stored in a compressed row storage (CRS) sparse
	 *  matrix format.  Sorption, decay, and TVD flux corrections are
	 *  incorporated during Picard iteration.
	 * ===================================================================== */

	/**
	 * @brief Assemble the linear system for one species and one Picard
	 *        iteration.
	 *
	 * The discretisation is fully implicit in time and uses a 27-point
	 * stencil in 3-D (6 face neighbours + 12 edge neighbours + 8 corner
	 * neighbours + the diagonal).  The advection term uses upwind weighting
	 * with the corrected velocity, while the dispersion term uses central
	 * differences.  Cross-dispersion terms (Dxy, Dxz, etc.) couple diagonal
	 * stencil entries.
	 *
	 * The right-hand side includes the mass storage term (old concentration
	 * times retardation factor) and, during Picard iterations (iter > 0),
	 * the TVD flux correction from the high-resolution scheme.
	 *
	 * After the interior assembly, boundary conditions are applied, root-
	 * uptake / source-sink terms are added, and MPI halo-boundary
	 * contributions are moved to the RHS.
	 *
	 * @param gw           Groundwater state (water content, etc.).
	 * @param gdom         Domain geometry (grid dimensions, spacing).
	 * @param rtA          Sparse matrix in CRS format (output).
	 * @param par          Parallel execution context.
	 * @param rt           Reactive transport state.
	 * @param rtgbc        Vector of reactive transport boundary condition objects.
	 * @param ss           Source/sink data (e.g. root uptake).
	 * @param iSpec        Species index (0-based).
	 * @param picard_iter  Current Picard iteration number (0 = first iteration,
	 *                     no TVD correction).
	 */
	inline void RTlinear_system(GwState &gw, GwDomain &gdom, RTMatrix &rtA, Parallel &par, RTStateGW &rt, std::vector<RTBCGW> &rtgbc, const SourceSink &ss, int iSpec = 0, int picard_iter = 0)
	{
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering RTlinear_system() for iSpec=" << iSpec << " picard_iter=" << picard_iter << std::endl;
#endif
		// -----------------------------------------------------------------
		// [P1-1] The face-coefficient precompute (wz, c_adv**, c_dif**) is
		// Picard-invariant and now runs once per species per half-step in
		// RTprepare_species_coefficients(), called by rt_solve() before the
		// Picard loop.  [P2-5] That function also stages the 19 Picard-
		// invariant stencil columns into RTcoef_inv; this kernel copies them
		// and only recomputes the truly Picard-dependent work (sorption mass
		// terms and TVD flux corrections, both functions of rt.c(..,1)).
		// -----------------------------------------------------------------
		auto coef_inv_l = RTcoef_inv;
		auto inv_dz_l = inv_dz;
		const real inv_dx_l = inv_dx, inv_dy_l = inv_dy;
		Kokkos::parallel_for("Build_Matrix_Fused", gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

			// [P1-7]/[P2-4] Cached geometry and one-time reciprocals: the
			// former divisions by dz are reciprocal multiplications (~1 ulp,
			// L2); dz_i is kept for the TVD flux-correction calls below.
			const int plane = gdom.nxhc * gdom.nyhc;
			const real dz_i = gdom.dz(iGlob);
			const real inv_dz_i = inv_dz_l(iGlob);
			const real dt_ = rt.dt;

			// Allocate local coefficient array in registers / fast cache.
			// [P2-5] Columns 1-6 and 8-26 are Picard-invariant: staged once
			// per species per half-step in RTcoef_inv by
			// RTprepare_species_coefficients() and copied here.  Buffer
			// columns 0 and 7 are never written and are immediately
			// overwritten below.
			real local_coef[27];
			for (int k = 0; k < 27; ++k)
				local_coef[k] = coef_inv_l(idom, k);

			// [P1-2] Upwind weighting factors depend only on the sign of the
			// half-step-invariant face velocities; they are now precomputed
			// once per half-step by RTcompute_upwind_weights() and are only
			// read here (applyRTMatBC reads the same values).

			// =================================================================
			// Mass conservation storage term with sorption and decay
			// =================================================================
			real intermediate_variable = 0.0; // Intermediate variable for kinetic sorption
			int sorption_type = rt.Sorption_Type(iSpec);
			real Kd_eff = 0.0; // Dynamic equivalent partition coefficient [L/kg]
			real sorption_conc = 0.0;
			real mass_lhs_diag = 0.0;
			real mass_rhs = 0.0;

				// During Picard iteration, use the concentration from the
				// previous iteration to evaluate nonlinear terms.
				real C_iter = fmax(0.0, rt.c(iSpec, iGlob, 1));

				Kd_eff = sorption_equilibrium_derivative(rt, iSpec, C_iter);
				sorption_conc = sorption_equilibrium_concentration(rt, iSpec, C_iter);
			if (sorption_type == 0) 
			{
				// No-sorption mode: disconnect the liquid equation from c_solid
				// so stale solid concentrations cannot re-enter the mass balance.
				mass_lhs_diag = gw.wc_new(iGlob) + dt_ * rt.Lambda_1(iSpec) * gw.wc_new(iGlob);
				mass_rhs = gw.wc_old(iGlob) * rt.c(iSpec, iGlob, 0);
			}

			else if (sorption_type == 4) // Non-equilibrium (kinetic) sorption
			{
				
				real beta_k = rt.Beta(iSpec);
				real Kd = rt.Kd(iSpec);
				real intermediate_variable = 0.0;
				if (Kd > 0.0) {
					intermediate_variable = beta_k / (Kd * (rt.rho_b / dt_ + beta_k / Kd + rt.Lambda_2(iSpec) * rt.rho_b));
				}
				// [P1-6] pow(Beta,2.0) -> beta_k * beta_k (same value, avoids pow)
				mass_lhs_diag = gw.wc_new(iGlob) + dt_ * (beta_k + rt.Lambda_1(iSpec) * gw.wc_new(iGlob) - intermediate_variable * beta_k * beta_k);
				mass_rhs = gw.wc_old(iGlob) * rt.c(iSpec, iGlob, 0) + intermediate_variable * rt.rho_b * rt.c_solid(iSpec, iGlob, 0);
			}
			else 
			{
				// Equilibrium adsorption (linear/Freundlich/Langmuir) and non-adsorption models
				// LHS = wc_new + rho_b * dS/dC + dt * Lambda_1 * wc_new (Natural decay in the liquid phase)
				mass_lhs_diag = gw.wc_new(iGlob) + rt.rho_b * Kd_eff + dt_ * rt.Lambda_1(iSpec) * gw.wc_new(iGlob);
				
				// RHS includes the Picard/Newton tangent intercept so nonlinear sorption remains conservative.
				real sorption_old = sorption_equilibrium_concentration(rt, iSpec, rt.c(iSpec, iGlob, 0));
				real tangent_intercept = rt.rho_b * (Kd_eff * C_iter - sorption_conc);
				mass_rhs = gw.wc_old(iGlob) * rt.c(iSpec, iGlob, 0) + rt.rho_b * sorption_old + tangent_intercept;
			}

			

			// --- Main diagonal (face-neighbour columns 1-6 and cross-dispersion
			// columns 8-26 are final from the [P2-5] buffer copy above) ---
			// local_coef[0]: diagonal (i,j,k)
			local_coef[0] = mass_lhs_diag + rt.c_difxx(iGlob) + rt.c_difxx(iGlob - 1) + rt.Up_Weighting_xp(iGlob) * rt.c_advxx(iGlob) - (1.0 - rt.Up_Weighting_xm(iGlob)) * rt.c_advxx(iGlob - 1) + rt.c_difyy(iGlob) + rt.c_difyy(iGlob - gdom.nxhc) + rt.Up_Weighting_yp(iGlob) * rt.c_advyy(iGlob) - (1.0 - rt.Up_Weighting_ym(iGlob)) * rt.c_advyy(iGlob - gdom.nxhc) + rt.c_difzz(iGlob) * inv_dz_i + rt.c_difzz(iGlob - plane) * inv_dz_i + rt.Up_Weighting_zp(iGlob) * rt.c_advzz(iGlob) * inv_dz_i - (1.0 - rt.Up_Weighting_zm(iGlob)) * rt.c_advzz(iGlob - plane) * inv_dz_i;

			// local_coef[7]: RHS (previous time-step mass + source/sink terms)
			local_coef[7] = mass_rhs;

			// -----------------------------------------------------------------
			// High-order TVD flux correction (added only during Picard
			// iterations iter > 0).  Strategy: keep the LHS as first-order
			// upwind; add the high-order correction explicitly to the RHS.
			// -----------------------------------------------------------------
			if (picard_iter > 0)
			{
				// Limiter type: 1=Minmod (most stable), 2=Superbee (may
				// overshoot), 3=Van Leer (recommended balance)
				int limiter_type = 1; // Default: Minmod limiter

				// === X-direction high-order flux correction ===

				// Face i+1/2 flux correction
				if (ii < gdom.nx - 1)
				{
					int iGlob_plus = iGlob + 1;
					int iGlob_minus = iGlob - 1;
					real C_upstream, C_downstream, C_upstream_upstream;

					// Determine upstream/downstream based on velocity direction
					if (rt.aveVB(iGlob, 0) >= 0) // v > 0: upstream = i, downstream = i+1
					{
						C_upstream = rt.c(iSpec, iGlob, 1);
						C_downstream = rt.c(iSpec, iGlob_plus, 1);
						C_upstream_upstream = (ii > 0) ? rt.c(iSpec, iGlob_minus, 1) : rt.c(iSpec, iGlob, 1);
					}
					else // v < 0: upstream = i+1, downstream = i
					{
						C_upstream = rt.c(iSpec, iGlob_plus, 1);
						C_downstream = rt.c(iSpec, iGlob, 1);
						C_upstream_upstream = (ii < gdom.nx - 2) ? rt.c(iSpec, iGlob_plus + 1, 1) : rt.c(iSpec, iGlob_plus, 1);
					}

					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
															rt.aveVB(iGlob, 0), dt_, gdom.dx, limiter_type);
					// Add divergence of flux correction to RHS:
					//   dt/dx * (F_corr_{i+1/2} - F_corr_{i-1/2})
					local_coef[7] -= dt_ * inv_dx_l * flux_corr;
				}

				// Face i-1/2 flux correction
				if (ii > 0)
				{
					int iGlob_minus = iGlob - 1;
					int iGlob_minus2 = iGlob - 2;
					real C_upstream, C_downstream, C_upstream_upstream;

					// Determine upstream/downstream based on velocity direction
					if (rt.aveVB(iGlob_minus, 0) >= 0) // v > 0: upstream = i-1, downstream = i
					{
						C_upstream = rt.c(iSpec, iGlob_minus, 1);
						C_downstream = rt.c(iSpec, iGlob, 1);
						C_upstream_upstream = (ii > 1) ? rt.c(iSpec, iGlob_minus2, 1) : rt.c(iSpec, iGlob_minus, 1);
					}
					else // v < 0: upstream = i, downstream = i-1
					{
						C_upstream = rt.c(iSpec, iGlob, 1);
						C_downstream = rt.c(iSpec, iGlob_minus, 1);
						C_upstream_upstream = (ii < gdom.nx - 1) ? rt.c(iSpec, iGlob + 1, 1) : rt.c(iSpec, iGlob, 1);
					}

					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
															rt.aveVB(iGlob_minus, 0), dt_, gdom.dx, limiter_type);
					// Subtract from RHS (flux divergence at i-1/2)
					local_coef[7] += dt_ * inv_dx_l * flux_corr;
				}

				// === Y-direction high-order flux correction ===
				if (jj < gdom.ny - 1)
				{
					int iGlob_plus = iGlob + gdom.nxhc;
					int iGlob_minus = iGlob - gdom.nxhc;
					real C_upstream, C_downstream, C_upstream_upstream;

					if (rt.aveVB(iGlob, 1) >= 0)
					{
						C_upstream = rt.c(iSpec, iGlob, 1);
						C_downstream = rt.c(iSpec, iGlob_plus, 1);
						C_upstream_upstream = (jj > 0) ? rt.c(iSpec, iGlob_minus, 1) : rt.c(iSpec, iGlob, 1);
					}
					else
					{
						C_upstream = rt.c(iSpec, iGlob_plus, 1);
						C_downstream = rt.c(iSpec, iGlob, 1);
						C_upstream_upstream = (jj < gdom.ny - 2) ? rt.c(iSpec, iGlob_plus + gdom.nxhc, 1) : rt.c(iSpec, iGlob_plus, 1);
					}

					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
															rt.aveVB(iGlob, 1), dt_, gdom.dy, limiter_type);
					local_coef[7] -= dt_ * inv_dy_l * flux_corr;
				}

				if (jj > 0)
				{
					int iGlob_minus = iGlob - gdom.nxhc;
					int iGlob_minus2 = iGlob - 2 * gdom.nxhc;
					real C_upstream, C_downstream, C_upstream_upstream;

					if (rt.aveVB(iGlob_minus, 1) >= 0)
					{
						C_upstream = rt.c(iSpec, iGlob_minus, 1);
						C_downstream = rt.c(iSpec, iGlob, 1);
						C_upstream_upstream = (jj > 1) ? rt.c(iSpec, iGlob_minus2, 1) : rt.c(iSpec, iGlob_minus, 1);
					}
					else
					{
						C_upstream = rt.c(iSpec, iGlob, 1);
						C_downstream = rt.c(iSpec, iGlob_minus, 1);
						C_upstream_upstream = (jj < gdom.ny - 1) ? rt.c(iSpec, iGlob + gdom.nxhc, 1) : rt.c(iSpec, iGlob, 1);
					}

					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
															rt.aveVB(iGlob_minus, 1), dt_, gdom.dy, limiter_type);
					local_coef[7] += dt_ * inv_dy_l * flux_corr;
				}

				// === Z-direction high-order flux correction ===
				if (kk < gdom.nz - 1)
				{
					int iGlob_plus = iGlob + gdom.nxhc * gdom.nyhc;
					int iGlob_minus = iGlob - gdom.nxhc * gdom.nyhc;
					real C_upstream, C_downstream, C_upstream_upstream;

					if (rt.aveVB(iGlob, 2) >= 0)
					{
						C_upstream = rt.c(iSpec, iGlob, 1);
						C_downstream = rt.c(iSpec, iGlob_plus, 1);
						C_upstream_upstream = (kk > 0) ? rt.c(iSpec, iGlob_minus, 1) : rt.c(iSpec, iGlob, 1);
					}
					else
					{
						C_upstream = rt.c(iSpec, iGlob_plus, 1);
						C_downstream = rt.c(iSpec, iGlob, 1);
						C_upstream_upstream = (kk < gdom.nz - 2) ? rt.c(iSpec, iGlob_plus + gdom.nxhc * gdom.nyhc, 1) : rt.c(iSpec, iGlob_plus, 1);
					}

					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
															rt.aveVB(iGlob, 2), dt_, dz_i, limiter_type);
					local_coef[7] -= dt_ * inv_dz_i * flux_corr;
				}

				if (kk > 0)
				{
					int iGlob_minus = iGlob - gdom.nxhc * gdom.nyhc;
					int iGlob_minus2 = iGlob - 2 * gdom.nxhc * gdom.nyhc;
					real C_upstream, C_downstream, C_upstream_upstream;

					if (rt.aveVB(iGlob_minus, 2) >= 0)
					{
						C_upstream = rt.c(iSpec, iGlob_minus, 1);
						C_downstream = rt.c(iSpec, iGlob, 1);
						C_upstream_upstream = (kk > 1) ? rt.c(iSpec, iGlob_minus2, 1) : rt.c(iSpec, iGlob_minus, 1);
					}
					else
					{
						C_upstream = rt.c(iSpec, iGlob, 1);
						C_downstream = rt.c(iSpec, iGlob_minus, 1);
						C_upstream_upstream = (kk < gdom.nz - 1) ? rt.c(iSpec, iGlob + gdom.nxhc * gdom.nyhc, 1) : rt.c(iSpec, iGlob, 1);
					}

					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
															rt.aveVB(iGlob_minus, 2), dt_, dz_i, limiter_type);
					local_coef[7] += dt_ * inv_dz_i * flux_corr;
				}
			}


			// =================================================================
			// Cross-dispersion stencil entries (local_coef[8..26] and their
			// face-column corrections to local_coef[1..6]) are Picard-
			// invariant and now come from the [P2-5] RTcoef_inv buffer; see
			// RTprepare_species_coefficients().
			// =================================================================

			// Final step: write all local coefficients to global memory.
			// [P1-7] The rt_rhs zero-fill that used to sit here is dead: the
			// final rt_rhs is unconditionally (over)written by the CRS value
			// gather at the end of this function, and neither applyRTMatBC()
			// nor ComputeRTGWSourceSink() touches rtA.rt_rhs.
			for(int k = 0; k < 27; k++) {
				rt.RTcoef(idom, k) = local_coef[k];
			} });

		// [P1-4] Backup of the pre-BC physical stencil for the residual-based
		// branch of compute_boundary_mass_fluxes().  Only the boundary rows
		// are ever read there (columns 0..25), so the previous full
		// nCell x 30 deep_copy is replaced by a boundary-row copy.  The
		// backed-up values are bitwise identical to the full copy.
		const int nrows_bak = nbres_rows;
		auto rows_local = bres_rows;
		auto orig_local = orig_RTcoef;
		Kokkos::parallel_for("Backup_Boundary_RTcoef", nrows_bak, KOKKOS_LAMBDA(int i) {
			const int idom = rows_local(i);
			for (int k = 0; k < 27; ++k)
				orig_local(idom, k) = rt.RTcoef(idom, k); });

		// -----------------------------------------------------------------
		// Stage 4: Apply boundary conditions to the matrix
		// -----------------------------------------------------------------
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> RTlinear_system(): applying RT boundary conditions" << std::endl;
#endif
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyRTMatBC(rt, gw, gdom, par, iSpec);
		}

		// Handle root water uptake and groundwater RT source/sink terms
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> RTlinear_system(): computing RT GW source/sink" << std::endl;
#endif
		RTGwSS::ComputeRTGWSourceSink(rt, gw, gdom, ss, iSpec);

		// =================================================================
		// [Critical Fix 4]: Move MPI halo-boundary stencil contributions
		// to the RHS.  For cells at the domain boundary whose 27-point
		// stencil extends beyond the local partition, the corresponding
		// neighbour coefficients are multiplied by the known halo
		// concentration and subtracted from the RHS, effectively making
		// them explicit contributions.
		// =================================================================
		Kokkos::parallel_for("MPI_RHS_Update", gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

            // X, Y, Z face-neighbour contributions
            if (!(ii < gdom.nx - 1)) rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 1) * rt.c(iSpec, iGlob + 1, 1);
            if (!(ii > 0))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 2) * rt.c(iSpec, iGlob - 1, 1);
            if (!(jj < gdom.ny - 1)) rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 3) * rt.c(iSpec, iGlob + gdom.nxhc, 1);
            if (!(jj > 0))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 4) * rt.c(iSpec, iGlob - gdom.nxhc, 1);
            if (!(kk < gdom.nz - 1)) rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 5) * rt.c(iSpec, iGlob + gdom.nxhc * gdom.nyhc, 1);
            if (!(kk > 0))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 6) * rt.c(iSpec, iGlob - gdom.nxhc * gdom.nyhc, 1);

            // XY / YX cross-dispersion contributions
            if (!(ii < gdom.nx - 1 && jj < gdom.ny - 1)) rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 9) * rt.c(iSpec, iGlob + 1 + gdom.nxhc, 1);
            if (!(ii < gdom.nx - 1 && jj > 0))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 11) * rt.c(iSpec, iGlob + 1 - gdom.nxhc, 1);
            if (!(ii > 0 && jj < gdom.ny - 1))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 12) * rt.c(iSpec, iGlob - 1 + gdom.nxhc, 1);
            if (!(ii > 0 && jj > 0))                     rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 13) * rt.c(iSpec, iGlob - 1 - gdom.nxhc, 1);

            // XZ / ZX cross-dispersion contributions
            if (!(ii < gdom.nx - 1 && kk < gdom.nz - 1)) rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 15) * rt.c(iSpec, iGlob + 1 + gdom.nxhc * gdom.nyhc, 1);
            if (!(ii < gdom.nx - 1 && kk > 0))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 17) * rt.c(iSpec, iGlob + 1 - gdom.nxhc * gdom.nyhc, 1);
            if (!(ii > 0 && kk < gdom.nz - 1))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 18) * rt.c(iSpec, iGlob - 1 + gdom.nxhc * gdom.nyhc, 1);
            if (!(ii > 0 && kk > 0))                     rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 19) * rt.c(iSpec, iGlob - 1 - gdom.nxhc * gdom.nyhc, 1);

            // YZ / ZY cross-dispersion contributions
            if (!(jj < gdom.ny - 1 && kk < gdom.nz - 1)) rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 22) * rt.c(iSpec, iGlob + gdom.nxhc + gdom.nxhc * gdom.nyhc, 1);
            if (!(jj < gdom.ny - 1 && kk > 0))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 23) * rt.c(iSpec, iGlob + gdom.nxhc - gdom.nxhc * gdom.nyhc, 1);
            if (!(jj > 0 && kk < gdom.nz - 1))           rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 24) * rt.c(iSpec, iGlob - gdom.nxhc + gdom.nxhc * gdom.nyhc, 1);
            if (!(jj > 0 && kk > 0))                     rt.RTcoef(idom, 7) -= rt.RTcoef(idom, 25) * rt.c(iSpec, iGlob - gdom.nxhc - gdom.nxhc * gdom.nyhc, 1); });
		// [P1-8] No fence needed here: the CRS gather kernels below run on
		// the same execution space (launch-order sequenced), and the caller
		// fences before reading the linsys timer.
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> RTlinear_system(): MPI RHS update completed" << std::endl;
#endif

		// -----------------------------------------------------------------
		// Stage 5: Fill the sparse matrix values in CRS (Compressed Row
		//          Storage) format.
		//
		// [P1-3] The CRS graph (rt_ind) and the per-entry (row, RTcoef-slot)
		// gather map are geometry-only: they were filled once by
		// RTMatrix::fillGraph() at initialisation, replicating the exact
		// boundary-checked 19-entry insertion order that used to live here.
		// Each assembly therefore only gathers the (possibly BC-modified)
		// coefficient values into rt_val and copies column 7 into rt_rhs;
		// the written values and their positions are identical to the
		// previous branch-based fill.
		// -----------------------------------------------------------------
		Kokkos::parallel_for("CRS_Gather_Values", rtA.rt_nnz, KOKKOS_LAMBDA(int j) {
			rtA.rt_val(j) = rt.RTcoef(rtA.rt_row_of(j), rtA.rt_slot_of(j)); });
		Kokkos::parallel_for("CRS_Set_RHS", rtA.rt_nrow, KOKKOS_LAMBDA(int i) {
			rtA.rt_rhs(i) = rt.RTcoef(i, 7); });
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> RTlinear_system() completed OK for iSpec=" << iSpec << std::endl;
#endif
	}

	inline void compute_boundary_mass_fluxes(RTStateGW &rt, GwDomain &gdom, std::vector<RTBCGW> &rtbc, int iSpec)
	{
		for (size_t k = 0; k < rtbc.size(); k++)
		{
			if (rtbc[k].ncellsBC <= 0)
				continue;

			int current_bctype = rtbc[k].rtbctype;
			if (iSpec < (int)rtbc[k].spec_bctype.size())
				current_bctype = rtbc[k].spec_bctype[iSpec];

			if (rtbc[k].spec_mass_in_rate.size() != rt.n_mass)
			{
				rtbc[k].spec_mass_in_rate.assign(rt.n_mass, 0.0);
				rtbc[k].spec_mass_out_rate.assign(rt.n_mass, 0.0);
			}

			if (current_bctype == SUB_RT_BC_SWE)
			{
				rtbc[k].spec_mass_in_rate[iSpec] = 0.0;
				rtbc[k].spec_mass_out_rate[iSpec] = 0.0;
				continue;
			}

			real q_in_local = 0.0;
			real q_out_local = 0.0;
			int ncells = rtbc[k].ncellsBC;
			auto bcells_local = rtbc[k].bcells;
			auto gcells_local = rtbc[k].gcells;
			int dir = rtbc[k].direction;
			int plane = gdom.nxhc * gdom.nyhc;

			if (current_bctype == SUB_RT_BC_FD)
			{
				Kokkos::parallel_reduce("RTGW_FD_Advective_Flux", ncells, KOKKOS_LAMBDA(int ibc, real &tmp_in, real &tmp_out) {
					int iGlob = bcells_local(ibc);
					int iGhost = gcells_local(ibc);

					real vn_out = 0.0;
					real area = 0.0;

					if (dir == XPLUS)
					{
						vn_out = rt.aveVB(iGlob, 0);
						area = gdom.dy * gdom.dz(iGlob);
					}
					else if (dir == XMINUS)
					{
						vn_out = -rt.aveVB(iGlob - 1, 0);
						area = gdom.dy * gdom.dz(iGlob);
					}
					else if (dir == YPLUS)
					{
						vn_out = rt.aveVB(iGlob, 1);
						area = gdom.dx * gdom.dz(iGlob);
					}
					else if (dir == YMINUS)
					{
						vn_out = -rt.aveVB(iGlob - gdom.nxhc, 1);
						area = gdom.dx * gdom.dz(iGlob);
					}
					else if (dir == ZPLUS)
					{
						vn_out = rt.aveVB(iGlob, 2);
						area = gdom.dx * gdom.dy;
					}
					else if (dir == ZMINUS)
					{
						vn_out = -rt.aveVB(iGlob - plane, 2);
						area = gdom.dx * gdom.dy;
					}
					else
					{
						return;
					}

					real c_inside = fmax(0.0, rt.c(iSpec, iGlob, 1));
					real c_external = fmax(0.0, rt.c(iSpec, iGhost, 1));
					real c_upwind = (vn_out >= 0.0) ? c_inside : c_external;
					real q_mass_out = vn_out * c_upwind * area * 1000.0;

					if (q_mass_out >= 0.0)
						tmp_out += q_mass_out;
					else
						tmp_in += -q_mass_out; }, Kokkos::Sum<real>(q_in_local), Kokkos::Sum<real>(q_out_local));
				// [P1-8] No fence needed: parallel_reduce already synchronised.

				rtbc[k].spec_mass_in_rate[iSpec] = q_in_local;
				rtbc[k].spec_mass_out_rate[iSpec] = q_out_local;
				rtbc[k].QMassInflow = q_in_local;
				rtbc[k].QMassOutflow = q_out_local;
				rtbc[k].QMassTot = q_out_local - q_in_local;
				continue;
			}

			auto orig_RTcoef_local = orig_RTcoef;

			Kokkos::parallel_reduce("RTGW_Boundary_Residual_Flux", ncells, KOKKOS_LAMBDA(int ibc, real &tmp_in, real &tmp_out) {
				int iGlob = bcells_local(ibc);
				int iGhost = gcells_local(ibc);
				int ii, jj, kk;
				gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
				int idom = (kk - gdom.hc) * gdom.nx * gdom.ny + (jj - gdom.hc) * gdom.nx + (ii - gdom.hc);

				const int nxhc = gdom.nxhc;
				const int plane = gdom.nxhc * gdom.nyhc;

				real c0 = rt.c(iSpec, iGlob, 1);
				real c1 = rt.c(iSpec, iGlob + 1, 1);
				real c2 = rt.c(iSpec, iGlob - 1, 1);
				real c3 = rt.c(iSpec, iGlob + nxhc, 1);
				real c4 = rt.c(iSpec, iGlob - nxhc, 1);
				real c5 = rt.c(iSpec, iGlob + plane, 1);
				real c6 = rt.c(iSpec, iGlob - plane, 1);

				real c9  = rt.c(iSpec, iGlob + 1 + nxhc, 1);
				real c11 = rt.c(iSpec, iGlob + 1 - nxhc, 1);
				real c12 = rt.c(iSpec, iGlob - 1 + nxhc, 1);
				real c13 = rt.c(iSpec, iGlob - 1 - nxhc, 1);

				real c15 = rt.c(iSpec, iGlob + 1 + plane, 1);
				real c17 = rt.c(iSpec, iGlob + 1 - plane, 1);
				real c18 = rt.c(iSpec, iGlob - 1 + plane, 1);
				real c19 = rt.c(iSpec, iGlob - 1 - plane, 1);

				real c22 = rt.c(iSpec, iGlob + nxhc + plane, 1);
				real c23 = rt.c(iSpec, iGlob + nxhc - plane, 1);
				real c24 = rt.c(iSpec, iGlob - nxhc + plane, 1);
				real c25 = rt.c(iSpec, iGlob - nxhc - plane, 1);

				// A Dirichlet boundary constrains the boundary cell itself.  For
				// residual-based mass-flux recovery, extend that boundary state
				// constantly in the outward normal direction.  This removes only
				// the artificial exterior dispersive flux; the constraint reaction
				// needed to maintain C_bc is still recovered from the physical row.
				// Apply the same mapping to every face and all associated cross terms.
				if (current_bctype == SUB_RT_BC_Dirichlet_CONST ||
				    current_bctype == SUB_RT_BC_Dirichlet_T)
				{
					if (dir == XPLUS)
					{
						c1 = c0; c9 = c3; c11 = c4; c15 = c5; c17 = c6;
					}
					else if (dir == XMINUS)
					{
						c2 = c0; c12 = c3; c13 = c4; c18 = c5; c19 = c6;
					}
					else if (dir == YPLUS)
					{
						c3 = c0; c9 = c1; c12 = c2; c22 = c5; c23 = c6;
					}
					else if (dir == YMINUS)
					{
						c4 = c0; c11 = c1; c13 = c2; c24 = c5; c25 = c6;
					}
					else if (dir == ZPLUS)
					{
						c5 = c0; c15 = c1; c18 = c2; c22 = c3; c24 = c4;
					}
					else if (dir == ZMINUS)
					{
						c6 = c0; c17 = c1; c19 = c2; c23 = c3; c25 = c4;
					}
				}

				real lhs_val = 0.0;
				lhs_val += orig_RTcoef_local(idom, 0) * c0;
				lhs_val += orig_RTcoef_local(idom, 1) * c1;
				lhs_val += orig_RTcoef_local(idom, 2) * c2;
				lhs_val += orig_RTcoef_local(idom, 3) * c3;
				lhs_val += orig_RTcoef_local(idom, 4) * c4;
				lhs_val += orig_RTcoef_local(idom, 5) * c5;
				lhs_val += orig_RTcoef_local(idom, 6) * c6;
				lhs_val += orig_RTcoef_local(idom, 9) * c9;
				lhs_val += orig_RTcoef_local(idom, 11) * c11;
				lhs_val += orig_RTcoef_local(idom, 12) * c12;
				lhs_val += orig_RTcoef_local(idom, 13) * c13;
				lhs_val += orig_RTcoef_local(idom, 15) * c15;
				lhs_val += orig_RTcoef_local(idom, 17) * c17;
				lhs_val += orig_RTcoef_local(idom, 18) * c18;
				lhs_val += orig_RTcoef_local(idom, 19) * c19;
				lhs_val += orig_RTcoef_local(idom, 22) * c22;
				lhs_val += orig_RTcoef_local(idom, 23) * c23;
				lhs_val += orig_RTcoef_local(idom, 24) * c24;
				lhs_val += orig_RTcoef_local(idom, 25) * c25;

				real residual = orig_RTcoef_local(idom, 7) - lhs_val;
				real volume = gdom.dx * gdom.dy * gdom.dz(iGlob);
				real rate = residual * volume * 1000.0 / rt.dt;

				bool positive_rate_is_inflow = (dir == 1 || dir == 3 || dir == 5);
				bool is_inflow = positive_rate_is_inflow ? (rate > 0.0) : (rate < 0.0);

				if (is_inflow)
					tmp_in += fabs(rate);
				else
					tmp_out += fabs(rate);

				// The mapped physical-row residual recovers the Dirichlet
				// constraint/dispersive contribution after suppressing the
				// artificial exterior concentration jump.  The prescribed water
				// flux through the physical face is a separate advective mass
				// transfer and must be included explicitly.  Use one outward-normal
				// sign convention for all six faces and the same staggered face
				// velocity used by the transport matrix.
				if (current_bctype == SUB_RT_BC_Dirichlet_CONST ||
				    current_bctype == SUB_RT_BC_Dirichlet_T)
				{
					real vn_out = 0.0;
					real area = 0.0;

					if (dir == XPLUS)
					{
						vn_out = rt.aveVB(iGlob, 0);
						area = gdom.dy * gdom.dz(iGlob);
					}
					else if (dir == XMINUS)
					{
						vn_out = -rt.aveVB(iGlob - 1, 0);
						area = gdom.dy * gdom.dz(iGlob);
					}
					else if (dir == YPLUS)
					{
						vn_out = rt.aveVB(iGlob, 1);
						area = gdom.dx * gdom.dz(iGlob);
					}
					else if (dir == YMINUS)
					{
						vn_out = -rt.aveVB(iGlob - nxhc, 1);
						area = gdom.dx * gdom.dz(iGlob);
					}
					else if (dir == ZPLUS)
					{
						vn_out = rt.aveVB(iGlob, 2);
						area = gdom.dx * gdom.dy;
					}
					else if (dir == ZMINUS)
					{
						vn_out = -rt.aveVB(iGlob - plane, 2);
						area = gdom.dx * gdom.dy;
					}

					const real c_inside = fmax(0.0, c0);
					const real c_external = fmax(0.0, rt.c(iSpec, iGhost, 1));
					const real c_upwind = (vn_out >= 0.0) ? c_inside : c_external;
					const real advective_mass_out = vn_out * c_upwind * area * 1000.0;

					if (advective_mass_out >= 0.0)
						tmp_out += advective_mass_out;
					else
						tmp_in += -advective_mass_out;
				} }, Kokkos::Sum<real>(q_in_local), Kokkos::Sum<real>(q_out_local));
			// [P1-8] No fence needed: parallel_reduce already synchronised.

			rtbc[k].spec_mass_in_rate[iSpec] = q_in_local;
			rtbc[k].spec_mass_out_rate[iSpec] = q_out_local;
			rtbc[k].QMassInflow = q_in_local;
			rtbc[k].QMassOutflow = q_out_local;
			rtbc[k].QMassTot = q_out_local - q_in_local;
		}
	}

	/* =====================================================================
	 *  Section: Monod Nitrogen Cycle Reaction Integration
	 *
	 *  Implements a fully-implicit, multi-species Monod kinetic model for
	 *  the nitrogen cycle in groundwater.  The five coupled species are:
	 *
	 *    0 = NH4  (ammonium)         — molecular mass: 18.04 g/mol
	 *    1 = NO3  (nitrate)          — molecular mass: 62.00 g/mol
	 *    2 = DO   (dissolved oxygen) — molecular mass: 32.00 g/mol
	 *    3 = DOC  (dissolved organic carbon)  — elemental mass: 12.01 g/mol (as C)
	 *    4 = DON  (dissolved organic nitrogen) — elemental mass: 14.01 g/mol (as N)
	 *
	 *  Concentration mass basis is species-specific, not a universal
	 *  nitrogen basis: NH4 and NO3 are compound masses, DO is O2 mass,
	 *  DOC is elemental carbon mass, and DON is elemental nitrogen mass.
	 *  Input concentrations are not automatically converted between
	 *  compound-mass and elemental-mass conventions.
	 *
	 *  The reaction network includes:
	 *    - Mineralisation  (DON -> NH4, aerobic)
	 *    - Nitrification   (NH4 -> NO3, aerobic)
	 *    - Denitrification (NO3 + DOC -> N2, anaerobic)
	 *    - Aerobic respiration (DOC + DO -> CO2)
	 *    - Soil re-aeration and organic matter release
	 *
	 *  Solved via Newton-Raphson with Crank-Nicolson time integration
	 *  and an analytical Jacobian.
	 * ===================================================================== */

	/**
	 * @brief Maximum number of species in the coupled reaction system.
	 *
	 * GPU constraint: NVIDIA CUDA allocates only 1024 bytes of stack
	 * space per thread by default, so this is kept at 5.
	 */
	static const int MAX_SPEC = 5;
	static const int MAX_NEWTON_ITER = 10;		   ///< Maximum Newton-Raphson iterations per cell
	static constexpr real NEWTON_TOLERANCE = 1e-6; ///< Newton convergence tolerance (L2 norm of residual)
	static constexpr real SMALL_CONC = 1e-15;	   ///< Minimum concentration to avoid log/0

	/**
	 * @brief Compute the residual vector for the Crank-Nicolson implicit
	 *        reaction system.
	 *
	 * Evaluates the Monod kinetic rates at both the old and new
	 * concentrations, then forms the Crank-Nicolson residual:
	 *
	 *   R_k = C_new_k - C_old_k - 0.5 * dt * (dCdt_new_k + dCdt_old_k) / Rf_k
	 *
	 * The Monod kinetics include:
	 *   1. Mineralisation:    R_min = mu_max_min * f(T) * [DON/(K_DON+DON)] * [DO/(K_DO_min+DO)]
	 *   2. Nitrification:     R_nit = mu_max_nit * f(T) * [NH4/(K_NH4+NH4)] * [DO/(K_DO+DO)]
	 *   3. Denitrification:   R_denit = mu_max_denit * f(T) * [NO3/(K_NO3+NO3)] * [DOC/(K_DOC+DOC)] * [Ki/(Ki+DO)]
	 *   4. Aerobic respiration: R_aer = mu_max_hetero * f(T) * [DOC/(K_DOC_aer+DOC)] * [DO/(K_DO_aer+DO)]
	 *   5. Groundwater sources: re-aeration + DOC/DON release
	 *
	 * @param C_new      New concentrations for the five species on their
	 *                   species-specific mass bases [mg/L or mg element/L].
	 * @param C_old      Old concentrations on the same species-specific bases.
	 * @param dt         Time step [s].
	 * @param Resid      Output residual vector (length 5).
	 * @param rt         Reactive transport state (kinetic parameters).
	 * @param Temp_Corr  Temperature correction factor (dimensionless).
	 * @param Rf         Retardation factors for each species [dimensionless].
	 * @param theta      Volumetric water content [dimensionless, m^3/m^3].
	 */
	KOKKOS_INLINE_FUNCTION
	void calculate_residual(const real C_new[], const real C_old[], real dt, real Resid[],
							const RTStateGW &rt, real Temp_Corr, const real Rf[], real theta, real theta_s)
	{
		if (theta < 1e-5)
			theta = 1e-5;
		const real SMALL_CONC = 1e-15;

		// --- Physical constants and environmental parameters ---
		// const real theta_s = 0.51;             // Saturated water content [m^3/m^3]
		const real Ka_DO = rt.Ka_DO;		 // Re-aeration rate constant [1/s]
		const real DO_sat = rt.DO_sat;		 // DO saturation concentration [mg/L]
		const real DOC_eq = rt.DOC_eq;		 // DOC release equilibrium concentration [mg C/L]
		const real DON_eq = rt.DON_eq;		 // DON release equilibrium concentration [mg N/L]
		const real K_rel_DON = rt.K_rel_DON; // DON release rate [1/s]
		const real K_rel_DOC = rt.K_rel_DOC; // DOC release rate [1/s]

		// =========================================================
		// Stoichiometric mass ratios. The base rates below use the
		// concentration basis of the directly consumed species:
		// R_min [mg N/L/s], R_nit [mg NH4/L/s], and R_denit/R_aer
		// [mg C/L/s]. These coefficients convert each base rate to the
		// basis of the affected species; they are not additional inputs.
		// =========================================================
		const real Y_NH4_DON = 18.04 / 14.01;		  // DON -> NH4  (molecular mass ratio)
		const real Y_NO3_NH4 = 62.00 / 18.04;		  // NH4 -> NO3  (molecular mass ratio)
		const real Y_O2_NH4 = 64.00 / 18.04;		  // Nitrification oxygen demand (2 mol O2 / 1 mol NH4)
		const real Y_O2_DON = 32.00 / 14.01;		  // Mineralisation oxygen demand (partial oxidation)
		const real Y_O2_DOC = 32.00 / 12.01;		  // Aerobic respiration (1 mol O2 / 1 mol C)
		const real Y_NO3_DOC = (0.8 * 62.00) / 12.01; // Denitrification C/N ratio (0.8 mol NO3 / 1 mol C)

		// Local lambda to compute concentration rate of change dC/dt
		auto calc_dCdt = [&](const real C[], real dCdt[])
		{
			real NH4 = fmax(0.0, C[0]);
			real NO3 = fmax(0.0, C[1]);
			real DO = fmax(0.0, C[2]);
			real DOC = fmax(0.0, C[3]);
			real DON = fmax(0.0, C[4]);

			// Base reaction rates. Units are species-specific; see the
			// rate comments below rather than interpreting all rates as mg/L/s.

			// 1. Mineralisation rate: consumes DON-N, [mg N/L/s].
			real R_min = rt.Rate_Max_Min * Temp_Corr * (DON / (rt.K_Monod_DON + DON + SMALL_CONC)) * (DO / (rt.K_Monod_DO_Min + DO + SMALL_CONC));

			// 2. Nitrification rate: consumes NH4 compound mass, [mg NH4/L/s].
			real R_nit = rt.Rate_Max_Nit * Temp_Corr * (NH4 / (rt.K_Monod_NH4 + NH4 + SMALL_CONC)) * (DO / (rt.K_Monod_DO + DO + SMALL_CONC));

			// 3. Denitrification rate: consumes DOC-C, [mg C/L/s].
			real R_denit = rt.Rate_Max_Denit * Temp_Corr * (NO3 / (rt.K_Monod_NO3 + NO3 + SMALL_CONC)) * (DOC / (rt.K_Monod_DOC + DOC + SMALL_CONC)) * (rt.Ki_Inhib_DO / (rt.Ki_Inhib_DO + DO + SMALL_CONC));

			// 4. Aerobic respiration rate: consumes DOC-C, [mg C/L/s].
			real R_aer = rt.Rate_Max_Hetero * Temp_Corr * (DOC / (rt.K_Monod_DOC_Aerobic + DOC + SMALL_CONC)) * (DO / (rt.K_Monod_DO_Aerobic + DO + SMALL_CONC));

			// 5. Groundwater source terms. R_air is [mg O2/L/s],
			// R_rel_DOC is [mg C/L/s], and R_rel_DON is [mg N/L/s].
			real air_filled_porosity = fmax(0.0, theta_s - theta);
			real R_air = Ka_DO * air_filled_porosity * fmax(0.0, DO_sat - DO);
			real R_rel_DOC = K_rel_DOC * fmax(0.0, DOC_eq - DOC);
			real R_rel_DON = K_rel_DON * fmax(0.0, DON_eq - DON);

			// Assemble source/sink terms on each species' concentration basis.
			dCdt[0] = Y_NH4_DON * R_min - R_nit;									  // NH4: mineralisation gain - nitrification loss
			dCdt[1] = Y_NO3_NH4 * R_nit - Y_NO3_DOC * R_denit;						  // NO3: nitrification gain - denitrification loss
			dCdt[2] = R_air - Y_O2_DON * R_min - Y_O2_NH4 * R_nit - Y_O2_DOC * R_aer; // DO: re-aeration - mineralisation/nitrification/respiration O2 demand
			dCdt[3] = R_rel_DOC - R_denit - R_aer;									  // DOC: release - denitrification consumption - respiration
			dCdt[4] = R_rel_DON - R_min;											  // DON: release - mineralisation
		};

		real dCdt_old[5] = {0}, dCdt_new[5] = {0};
		calc_dCdt(C_old, dCdt_old);
		calc_dCdt(C_new, dCdt_new);

		// Crank-Nicolson residual:
		//   R_k = C_new_k - C_old_k - 0.5 * dt * (dCdt_new_k + dCdt_old_k) / Rf_k
		for (int k = 0; k < 5; k++)
		{
			Resid[k] = C_new[k] - C_old[k] - 0.5 * dt * (dCdt_new[k] + dCdt_old[k]) / Rf[k];
		}
	}

	/**
	 * @brief Compute the 5x5 analytical Jacobian matrix for the coupled
	 *        Monod reaction system.
	 *
	 * The Jacobian J[i][j] = dR_i / dC_new_j is derived by applying the
	 * chain rule to the Monod rate expressions.  It is used by the
	 * Newton-Raphson solver.
	 *
	 * @param C_new      Current guess concentrations (5 species) [mg/L].
	 * @param dt         Time step [s].
	 * @param Jacobian   Output 5x5 Jacobian matrix.
	 * @param rt         Reactive transport state (kinetic parameters).
	 * @param Temp_Corr  Temperature correction factor (dimensionless).
	 * @param Rf         Retardation factors for each species [dimensionless].
	 * @param theta      Volumetric water content [m^3/m^3].
	 */
	KOKKOS_INLINE_FUNCTION
	void calculate_jacobian(const real C_new[], real dt, real Jacobian[MAX_SPEC][MAX_SPEC],
							const RTStateGW &rt, real Temp_Corr, const real Rf[], real theta, real theta_s)
	{
		const real SMALL_CONC = 1e-15;

		// const real theta_s = 0.51;
		const real Ka_DO = rt.Ka_DO;		 // Re-aeration rate constant [1/s]
		const real K_rel_DOC = rt.K_rel_DOC; // DOC release rate [1/s]
		const real K_rel_DON = rt.K_rel_DON; // DON release rate [1/s]
		const real DO_sat = rt.DO_sat;

		// Stoichiometric mass ratios, expressed as mg product per mg base rate.
		const real Y_NH4_DON = 18.04 / 14.01;
		const real Y_O2_DON = 32.00 / 14.01;
		const real Y_NO3_NH4 = 62.00 / 18.04;
		const real Y_O2_NH4 = 64.00 / 18.04;
		const real Y_NO3_DOC = (0.8 * 62.00) / 12.01;
		const real Y_O2_DOC = 32.00 / 12.01;

		real NH4 = fmax(0.0, C_new[0]);
		real NO3 = fmax(0.0, C_new[1]);
		real DO = fmax(0.0, C_new[2]);
		real DOC = fmax(0.0, C_new[3]);
		real DON = fmax(0.0, C_new[4]);

		// =========================================================
		// Partial derivatives (chain rule on Monod expressions)
		// Using dynamic parameters read from the rt object.
		// =========================================================

		// --- Monod term derivative helper: d/dC [ C / (K + C) ] = K / (K + C)^2 ---
		auto deriv_monod = [=](real C, real K)
		{
			return K / pow(K + C + SMALL_CONC, 2);
		};

		// 1. Mineralisation partial derivatives
		real dRmin_dDON = rt.Rate_Max_Min * Temp_Corr * deriv_monod(DON, rt.K_Monod_DON) * (DO / (rt.K_Monod_DO_Min + DO + SMALL_CONC));
		real dRmin_dDO = rt.Rate_Max_Min * Temp_Corr * (DON / (rt.K_Monod_DON + DON + SMALL_CONC)) * deriv_monod(DO, rt.K_Monod_DO_Min);

		// 2. Nitrification partial derivatives
		real dRnit_dNH4 = rt.Rate_Max_Nit * Temp_Corr * deriv_monod(NH4, rt.K_Monod_NH4) * (DO / (rt.K_Monod_DO + DO + SMALL_CONC));
		real dRnit_dDO = rt.Rate_Max_Nit * Temp_Corr * (NH4 / (rt.K_Monod_NH4 + NH4 + SMALL_CONC)) * deriv_monod(DO, rt.K_Monod_DO);

		// 3. Denitrification partial derivatives (includes DO inhibition)
		real I_DO = rt.Ki_Inhib_DO / (rt.Ki_Inhib_DO + DO + SMALL_CONC);
		real dI_DO_dDO = -rt.Ki_Inhib_DO / pow(rt.Ki_Inhib_DO + DO + SMALL_CONC, 2);
		real dRdenit_dNO3 = rt.Rate_Max_Denit * Temp_Corr * deriv_monod(NO3, rt.K_Monod_NO3) * (DOC / (rt.K_Monod_DOC + DOC + SMALL_CONC)) * I_DO;
		real dRdenit_dDOC = rt.Rate_Max_Denit * Temp_Corr * (NO3 / (rt.K_Monod_NO3 + NO3 + SMALL_CONC)) * deriv_monod(DOC, rt.K_Monod_DOC) * I_DO;
		real dRdenit_dDO = rt.Rate_Max_Denit * Temp_Corr * (NO3 / (rt.K_Monod_NO3 + NO3 + SMALL_CONC)) * (DOC / (rt.K_Monod_DOC + DOC + SMALL_CONC)) * dI_DO_dDO;

		// 4. Aerobic respiration partial derivatives
		real dRaer_dDOC = rt.Rate_Max_Hetero * Temp_Corr * deriv_monod(DOC, rt.K_Monod_DOC_Aerobic) * (DO / (rt.K_Monod_DO_Aerobic + DO + SMALL_CONC));
		real dRaer_dDO = rt.Rate_Max_Hetero * Temp_Corr * (DOC / (rt.K_Monod_DOC_Aerobic + DOC + SMALL_CONC)) * deriv_monod(DO, rt.K_Monod_DO_Aerobic);

		// 5. Groundwater source-term partial derivatives
		real air_filled_porosity = fmax(0.0, theta_s - theta);
		real dRair_dDO = (DO < DO_sat) ? (-Ka_DO * air_filled_porosity) : 0.0;

		// --- Assemble the Jacobian matrix ---
		// Start with identity (from the C_new - C_old term)
		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				Jacobian[i][j] = (i == j) ? 1.0 : 0.0;

		// Fill the 5x5 Jacobian with Crank-Nicolson time-integrated derivatives
		// Eq 0: NH4
		Jacobian[0][0] -= 0.5 * dt * (-dRnit_dNH4) / Rf[0];
		Jacobian[0][2] -= 0.5 * dt * (Y_NH4_DON * dRmin_dDO - dRnit_dDO) / Rf[0];
		Jacobian[0][4] -= 0.5 * dt * (Y_NH4_DON * dRmin_dDON) / Rf[0];

		// Eq 1: NO3
		Jacobian[1][0] -= 0.5 * dt * (Y_NO3_NH4 * dRnit_dNH4) / Rf[1];
		Jacobian[1][1] -= 0.5 * dt * (-Y_NO3_DOC * dRdenit_dNO3) / Rf[1];
		Jacobian[1][2] -= 0.5 * dt * (Y_NO3_NH4 * dRnit_dDO - Y_NO3_DOC * dRdenit_dDO) / Rf[1];
		Jacobian[1][3] -= 0.5 * dt * (-Y_NO3_DOC * dRdenit_dDOC) / Rf[1];

		// Eq 2: DO
		Jacobian[2][0] -= 0.5 * dt * (-Y_O2_NH4 * dRnit_dNH4) / Rf[2];
		Jacobian[2][2] -= 0.5 * dt * (dRair_dDO - Y_O2_DON * dRmin_dDO - Y_O2_NH4 * dRnit_dDO - Y_O2_DOC * dRaer_dDO) / Rf[2];
		Jacobian[2][3] -= 0.5 * dt * (-Y_O2_DOC * dRaer_dDOC) / Rf[2];
		Jacobian[2][4] -= 0.5 * dt * (-Y_O2_DON * dRmin_dDON) / Rf[2];

		// Eq 3: DOC
		Jacobian[3][1] -= 0.5 * dt * (-dRdenit_dNO3) / Rf[3];
		Jacobian[3][2] -= 0.5 * dt * (-dRdenit_dDO - dRaer_dDO) / Rf[3];
		Jacobian[3][3] -= 0.5 * dt * (-K_rel_DOC - dRdenit_dDOC - dRaer_dDOC) / Rf[3];

		// Eq 4: DON
		Jacobian[4][2] -= 0.5 * dt * (-dRmin_dDO) / Rf[4];
		Jacobian[4][4] -= 0.5 * dt * (-K_rel_DON - dRmin_dDON) / Rf[4];
	}

	/**
	 * @brief Solve a small dense linear system using LU decomposition
	 *        (Doolittle algorithm) with forward/back substitution.
	 *
	 * Used internally by the Newton-Raphson solver to compute the
	 * concentration update:
	 *   J * Delta_C = -Residual
	 *
	 * @param Jacobian   Input 5x5 matrix (overwritten with LU factors).
	 * @param Resid      Residual vector (length n_spec).
	 * @param Delta_C    Output solution vector (length n_spec).
	 * @param n_spec     Number of species (typically 5).
	 */
	KOKKOS_INLINE_FUNCTION
	void solve_small_system_NxN(real Jacobian[MAX_SPEC][MAX_SPEC], real Resid[],
								real Delta_C[], int n_spec)
	{
		// Copy Jacobian to local working array (LU decomposition is in-place)
		real A[MAX_SPEC][MAX_SPEC];
		for (int i = 0; i < n_spec; i++)
		{
			for (int j = 0; j < n_spec; j++)
			{
				A[i][j] = Jacobian[i][j];
			}
		}

		// LU decomposition (Doolittle, in-place, no pivoting)
		for (int k = 0; k < n_spec; k++)
		{
			// Upper triangular part U[k][j]
			for (int j = k; j < n_spec; j++)
			{
				real sum = 0.0;
				for (int p = 0; p < k; p++)
				{
					sum += A[k][p] * A[p][j];
				}
				A[k][j] -= sum;
			}
			// Lower triangular part L[i][k]
			for (int i = k + 1; i < n_spec; i++)
			{
				real sum = 0.0;
				for (int p = 0; p < k; p++)
				{
					sum += A[i][p] * A[p][k];
				}
				A[i][k] = (A[i][k] - sum) / A[k][k];
			}
		}

		// Forward substitution: L * y = -Residual
		real y[MAX_SPEC];
		for (int i = 0; i < n_spec; i++)
		{
			real sum = 0.0;
			for (int j = 0; j < i; j++)
			{
				sum += A[i][j] * y[j];
			}
			y[i] = (-Resid[i] - sum);
		}

		// Back substitution: U * Delta_C = y
		for (int i = n_spec - 1; i >= 0; i--)
		{
			real sum = 0.0;
			for (int j = i + 1; j < n_spec; j++)
			{
				sum += A[i][j] * Delta_C[j];
			}
			Delta_C[i] = (y[i] - sum) / A[i][i];
		}
	}

	/**
	 * @brief Compute the L2 norm of a vector.
	 *
	 * @param vec     Input vector.
	 * @param n_spec  Vector length.
	 * @return        L2 norm (Euclidean norm).
	 */
	KOKKOS_INLINE_FUNCTION
	real calculate_norm(const real vec[], int n_spec) const
	{
		real sum = 0.0;
		for (int i = 0; i < n_spec; i++)
		{
			sum += vec[i] * vec[i];
		}
		return sqrt(sum);
	}

	/* =====================================================================
	 *  Section: Reaction Integration Kernel
	 *
	 *  Performs cell-local fully-implicit Newton-Raphson integration of
	 *  the coupled Monod nitrogen cycle reactions for all grid cells.
	 * ===================================================================== */

	/**
	 * @brief Execute the reaction integration kernel for all cells.
	 *
	 * For each cell, the method:
	 *   1. Checks if concentrations are above a physical-zero threshold.
	 *   2. Computes retardation factors for each species based on
	 *      sorption type.
	 *   3. Runs a Newton-Raphson loop (up to MAX_NEWTON_ITER) to solve
	 *      the fully-implicit Crank-Nicolson reaction system.
	 *   4. Writes the updated concentrations back to the global arrays.
	 *
	 * Only active when rt.Nitrogen_Cycle_Simulation == 1.
	 *
	 * @param rt    Reactive transport state.
	 * @param gw    Groundwater state (water content).
	 * @param gdom  Domain geometry.
	 * @param par   Parallel execution context.
	 */
	inline void reaction_integration_kernel(RTStateGW &rt, GwState &gw, GwDomain &gdom, Parallel &par)
	{
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering reaction_integration_kernel()" << std::endl;
#endif
		// ============================================================
		// Multi-species coupled reaction processing (Monod nitrogen
		// cycle) — no longer includes separate sorption and decay
		// modules (those are handled in the transport step).
		// ============================================================
		if (rt.Nitrogen_Cycle_Simulation == 1)
		{
			Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
            real dt = rt.dt;

            real theta = fmax(1e-5, gw.wc(iGlob, 0));  // Volumetric water content [m^3/m^3]
			int ivg = gw.soilID(iGlob) * NVG;
			real theta_s = gw.vgTable(ivg + 2); // wcsat from van Genuchten table for the soil type

			// TODO: Couple this kernel to the thermal module. Until then,
			// groundwater reaction rates use a fixed correction factor of 1.0.
			real Temp_Corr = 1.0;

            const int n_spec = 5; // Fixed at 5 for the nitrogen cycle species
            if (rt.n_mass < 5) return;

			// 1. Fast check: skip reaction computation if all species
			//    concentrations are below the physical-zero threshold.
            bool has_substance = false;
            for (int k = 0; k < n_spec; k++) {
                if (rt.c(k, iGlob, 1) > TOL12) { // Physical zero threshold
                    has_substance = true;
                    break;
                }
            }
            if (!has_substance) {
                // Force concentrations to zero to prevent tiny values like 1e-166
                for (int k = 0; k < n_spec; k++) rt.c(k, iGlob, 1) = 0.0;
                return;
            }

            // 2. Extract retardation factors Rf for each species
			real Rf[MAX_SPEC];
			for (int k = 0; k < n_spec; k++) {
				Rf[k] = 1.0;
				real C_liquid = rt.c(k, iGlob, 1);
				real dSdC = sorption_equilibrium_derivative(rt, k, C_liquid);
				if (rt.Sorption_Type(k) != 0) {
					Rf[k] = 1.0 + (rt.rho_b * dSdC) / theta;
				}
            }

            // 3. Newton-Raphson iteration for the fully-implicit reaction system
            real C_old[MAX_SPEC];
            real C_new[MAX_SPEC];
            real Resid[MAX_SPEC];
            real Jacobian[MAX_SPEC][MAX_SPEC];
            real Delta_C[MAX_SPEC];

            for (int k = 0; k < n_spec; k++) {
                C_old[k] = rt.c(k, iGlob, 1);
                C_new[k] = C_old[k];
            }

            for (int iter = 0; iter < MAX_NEWTON_ITER; iter++) {
                calculate_residual(C_new, C_old, dt, Resid, rt, Temp_Corr, Rf, theta, theta_s);

                // Check convergence (L2 norm of residual)
                real norm_resid = 0.0;
                for (int k = 0; k < n_spec; k++) norm_resid += Resid[k] * Resid[k];
                if (sqrt(norm_resid) < NEWTON_TOLERANCE) break;

                // Compute Jacobian and solve the linear system
                calculate_jacobian(C_new, dt, Jacobian, rt, Temp_Corr, Rf, theta, theta_s);

                solve_small_system_NxN(Jacobian, Resid, Delta_C, n_spec);

                // Update concentrations with positivity enforcement
                for (int k = 0; k < n_spec; k++) {
                    real next_val = C_new[k] + Delta_C[k];
                    if (next_val < 0) next_val = 1e-15;
                    C_new[k] = next_val;
                }
            }

            // Write updated concentrations back to global arrays
            for (int k = 0; k < n_spec; k++) rt.c(k, iGlob, 1) = C_new[k]; });

			Kokkos::fence();
		}
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> reaction_integration_kernel() completed OK (Nitrogen_Cycle_Simulation=" << rt.Nitrogen_Cycle_Simulation << ")" << std::endl;
#endif
	}

	/* =====================================================================
	 *  Section: Main Solver — Strang Splitting Driver
	 *
	 *  The overall algorithm uses Strang (operator) splitting to separate
	 *  the physical transport from the chemical reactions:
	 *
	 *    Step 1: Transport(dt/2)   — implicit advection-dispersion
	 *    Step 2: Reaction(dt)      — cell-local Monod kinetics
	 *    Step 3: Transport(dt/2)   — implicit advection-dispersion
	 *
	 *  Each transport step involves:
	 *    - MPI halo exchange for concentrations
	 *    - Dispersion tensor computation
	 *    - Per-species Picard iteration (linear system assembly + solve)
	 *    - Concentration and solid-phase update with sorption
	 *    - Convergence check via MPI global reduction
	 *
	 *  Detailed timing is recorded for each computational phase.
	 * ===================================================================== */

	/**
	 * @brief Main reactive transport solver using Strang splitting.
	 *
	 * @tparam execution_space  Kokkos execution space (e.g. Kokkos::DefaultExecutionSpace).
	 * @tparam type_solver      Linear solver type (e.g. RTSolverGW with
	 *                          Gauss-Seidel or GMRES).
	 *
	 * @param rt           Reactive transport state.
	 * @param rtA          Sparse matrix (CRS format) for the linear system.
	 * @param gw           Groundwater state.
	 * @param gdom         Domain geometry.
	 * @param rtgbc        Boundary condition objects.
	 * @param gmpi         MPI communication handler.
	 * @param par          Parallel execution context.
	 * @param rtint        Reactive transport integrator (mass tracking).
	 * @param RTSolverGW   Linear solver instance.
	 * @param ss           Source/sink data (e.g. root uptake).
	 */
	template <typename execution_space, typename type_solver>
	inline void rt_solve(RTStateGW &rt, RTMatrix &rtA, GwState &gw, GwDomain &gdom, std::vector<RTBCGW> &rtgbc, GwMPI &gmpi, GwMPI &rtgmpi, Parallel &par, RTIntegrator &rtint, type_solver &RTSolverGW, const SourceSink &ss)
	{
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rt_solve() (Strang Splitting) n_mass=" << rt.n_mass << " dt=" << gdom.dt << std::endl;
#endif
		int picard_iter_max = 3; // Maximum Picard iterations per species
		real picard_tol = 1e-6;	 // Picard convergence tolerance
		real dt_original = gdom.dt;

		// Store the original time step
		rt.dt = dt_original;
		if (orig_RTcoef.extent(0) != (size_t)gdom.nCell || orig_RTcoef.extent(1) != 30)
		{
			orig_RTcoef = realArr2("orig_RTcoef", gdom.nCell, 30);
			// [P1-4] Collect the interior rows of the boundary cells whose
			// RTcoef must be backed up.  Only the residual-based branch of
			// compute_boundary_mass_fluxes() reads orig_RTcoef (columns
			// 0..25); FD/SWE-type boundaries never touch it.  Built once per
			// allocation (rtgbc is fixed after initialisation).
			std::set<int> rows;
			for (size_t k = 0; k < rtgbc.size(); k++)
			{
				if (rtgbc[k].ncellsBC <= 0)
					continue;
				bool needs_backup =
					(rtgbc[k].rtbctype != SUB_RT_BC_FD && rtgbc[k].rtbctype != SUB_RT_BC_SWE);
				for (int s = 0; s < (int)rtgbc[k].spec_bctype.size(); ++s)
					if (rtgbc[k].spec_bctype[s] != SUB_RT_BC_FD &&
						rtgbc[k].spec_bctype[s] != SUB_RT_BC_SWE)
						needs_backup = true; // any species on the residual branch
				if (!needs_backup)
					continue;
				auto h_bcells = Kokkos::create_mirror_view(rtgbc[k].bcells);
				Kokkos::deep_copy(h_bcells, rtgbc[k].bcells);
				for (int i = 0; i < rtgbc[k].ncellsBC; ++i)
				{
					int ii, jj, kk;
					gdom.unpackIndicesHalo(h_bcells(i), kk, jj, ii);
					rows.insert((kk - gdom.hc) * gdom.nx * gdom.ny
								+ (jj - gdom.hc) * gdom.nx + (ii - gdom.hc));
				}
			}
			nbres_rows = (int)rows.size();
			bres_rows = intArr("bres_rows", nbres_rows);
			auto h_rows = Kokkos::create_mirror_view(bres_rows);
			int c = 0;
			for (int r : rows)
				h_rows(c++) = r;
			Kokkos::deep_copy(bres_rows, h_rows);
		}

		// [P2-4] One-time pure-geometry reciprocals and wz (the grid never
		// changes during a simulation).
		if (inv_dz.extent(0) != (size_t)gdom.nCellMem)
		{
			RTprecompute_geometry(gdom, rt);
		}

		// Performance timers for each computational phase
		Kokkos::Timer timer_bc, timer_comm, timer_disp, timer_linsys, timer_linsol;
		Kokkos::Timer timer_react, timer_update, timer_audit, timer_prepare, timer_postprocess;
		Kokkos::Timer timer_total;
		timer_total.reset();

		// [Critical Fix 3]: Synchronise water content across MPI
		// boundaries so that the storage matrix diagonal is seamless.
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve(): MPI sync water content" << std::endl;
#endif
		timer_comm.reset();
		gmpi.mpi_sendrecv1(gw.wc_old, gdom, par);
		gmpi.mpi_sendrecv1(gw.wc_new, gdom, par);
		Kokkos::fence();
		gdom.timers.rtgw.comm += timer_comm.seconds();
		// Compute mid-point water content for the Strang splitting
		timer_prepare.reset();
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
		    rt.wc_original_old(iGlob) = gw.wc_old(iGlob);
		    rt.wc_original_new(iGlob) = gw.wc_new(iGlob);
		    rt.wc_mid(iGlob) = 0.5 * (gw.wc_old(iGlob) + gw.wc_new(iGlob)); });
		Kokkos::fence();
		gdom.timers.rtgw.prepare += timer_prepare.seconds();

		// ============================================================
		// Step 1: Physical Transport — dt/2
		// ============================================================
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve(): === Step 1: Transport (dt/2) ===" << std::endl;
#endif
		rt.dt = dt_original / 2.0;
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
		timer_audit.reset();
		const std::vector<real> mass_before_half_1 =
			rtint.computeLocalTotalMasses(rt, gdom, gw.wc_old);
		Kokkos::fence();
		gdom.timers.rtgw.massAudit += timer_audit.seconds();
#endif

		// Use mid-point water content for the first half-step
		timer_prepare.reset();
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { gw.wc_new(iGlob) = rt.wc_mid(iGlob); });
		Kokkos::fence();
		gdom.timers.rtgw.prepare += timer_prepare.seconds();

		// Apply boundary conditions
		timer_bc.reset();
		// Start a complete boundary refresh with an empty exchange field.
		// applyConcentrationBC() is called for multiple boundary objects; it
		// must not clear ConQss itself or a later non-SWE boundary will erase
		// the species-0 exchange flux written by the top SWE boundary.
		Kokkos::deep_copy(rt.ConQss, 0.0);
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec, false);
			}
		}
		Kokkos::fence(); // Ensure BC application completes
		gdom.timers.rtgw.bc += timer_bc.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Step1: BC application completed" << std::endl;
#endif
		// MPI concentration halo exchange
		timer_comm.reset();
		rtgmpi.mpi_sendrecv_rt(rt.c, gdom, par, rt.n_mass, 1);
		gdom.timers.rtgw.comm += timer_comm.seconds();

		// Compute dispersion tensor
		timer_disp.reset();
		dispersion_tensor(rt, gw, gdom, gmpi, par);
		Kokkos::fence(); // Ensure dispersion tensor computation completes
		gdom.timers.rtgw.disp += timer_disp.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Step1: dispersion tensor computed, starting Picard iteration" << std::endl;
#endif
		// [P1-2] Upwind weights depend only on the (now final) aveVB signs:
		// compute once per half-step instead of inside every assembly kernel.
		timer_prepare.reset();
		RTcompute_upwind_weights(gdom, rt);
		Kokkos::fence();
		gdom.timers.rtgw.prepare += timer_prepare.seconds();
		// Back up current concentrations as the "old" time level for Picard
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
			{
				rt.c(iSpec, iGlob,  0) = rt.c(iSpec, iGlob, 1);
				rt.c_solid(iSpec, iGlob,  0) = rt.c_solid(iSpec, iGlob, 1); // Synchronise solid-phase initial concentration
			} });

		// Per-species Picard iteration loop
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			// [P1-1] Picard-invariant face coefficients: computed once per
			// species per half-step instead of once per Picard iteration.
			timer_prepare.reset();
			RTprepare_species_coefficients(gdom, rt, iSpec);
			Kokkos::fence();
			gdom.timers.rtgw.prepare += timer_prepare.seconds();

			rtPicardEnsure(rt.n_mass); // [P3-4] Picard statistics
			rtPicardEpisodes[iSpec]++;
			for (int picard_iter = 0; picard_iter < picard_iter_max; picard_iter++)
			{
				rtPicardIters[iSpec]++; // [P3-4]
				// Initialise the linear solver with the current guess
				Kokkos::parallel_for("Init_Solver_X", gdom.nCell, KOKKOS_LAMBDA(int idom) {
	            int ii, jj, kk, iGlob;
	            gdom.unpackIndices(idom, kk, jj, ii);
	            iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

	            rtA.rt_x(idom) = rt.c(iSpec, iGlob, 1); });
				Kokkos::fence();

				// Assemble the linear system
				timer_linsys.reset();
				RTlinear_system(gw, gdom, rtA, par, rt, rtgbc, ss, iSpec, picard_iter);
				Kokkos::fence();
				gdom.timers.rtgw.linsys += timer_linsys.seconds();

				// Solve the linear system
				timer_linsol.reset();
				if (rt.rt_scheme == 1)
				{
					RTSolverGW.Gauss_Seidel(rtA);
				}
				else if (rt.rt_scheme == 2)
				{
					RTSolverGW.Jacobi_GMRES_Solve(rtA, false, "Step1", iSpec, picard_iter, 0.0, 0.0, rt.dt);
				}
				Kokkos::fence();
				gdom.timers.rtgw.linsol += timer_linsol.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
				if (par.masterproc)
					std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Step1: linear solve done for iSpec=" << iSpec << " picard_iter=" << picard_iter << " scheme=" << rt.rt_scheme << std::endl;
#endif
				// Update concentrations and solid-phase sorption
				timer_update.reset();
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
				real clipped_mass_local = 0.0;
				Kokkos::parallel_reduce("RTGW_Step1_NegativeClipMass", gdom.nCell, KOKKOS_LAMBDA(int idom, real &tmp) {
						if (rtA.rt_x(idom) < 0.0)
						{
							int ii, jj, kk;
							gdom.unpackIndices(idom, kk, jj, ii);
							const int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc
							                + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
							const real volume = gdom.dx * gdom.dy * gdom.dz(iGlob);
							tmp += -rtA.rt_x(idom) * gw.wc_new(iGlob) * volume * 1000.0;
						} }, Kokkos::Sum<real>(clipped_mass_local));
				rtint.AccumNegativeClippedMass_spec(iSpec) += clipped_mass_local;
#endif
				real err_max_local = 0.0;
				Kokkos::parallel_reduce("Update_Concentration", gdom.nCell, KOKKOS_LAMBDA(int idom, real &lmax) {
					int ii, jj, kk, iGlob;
					gdom.unpackIndices(idom, kk, jj, ii);
					iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

					// Update liquid-phase concentration
					real C_old_iter = rt.c(iSpec, iGlob, 1);

					rt.c(iSpec, iGlob, 1) = rtA.rt_x(idom);
					if (rt.c(iSpec, iGlob, 1) < 0.0)
						rt.c(iSpec, iGlob, 1) = 0.0;

					// Update solid-phase concentration based on the new
					// liquid-phase value (diagnostic sorption update)
					int sorption_type = rt.Sorption_Type(iSpec);
					real C_liquid = rt.c(iSpec, iGlob, 1);

					if (sorption_type == 1 || sorption_type == 2 || sorption_type == 3)
					{
						rt.c_solid(iSpec, iGlob, 1) = sorption_equilibrium_concentration(rt, iSpec, C_liquid);
					}
					else if (sorption_type == 4)  // Kinetic sorption: implicit update
					{
						real dt = rt.dt;
						real Kd = rt.Kd(iSpec);
						real beta_k = rt.Beta(iSpec);
						real lambda_2 = rt.Lambda_2(iSpec);
						real rho_b = rt.rho_b;
						real S_old = rt.c_solid(iSpec, iGlob, 0);
						if (Kd > 0.0)
						{
							real denom = rho_b / dt + beta_k / Kd + lambda_2 * rho_b;
							rt.c_solid(iSpec, iGlob, 1) = (rho_b / dt * S_old + beta_k * C_liquid) / denom;
						}
					}

					real diff = fabs(rt.c(iSpec, iGlob, 1) - C_old_iter);
					if (diff > lmax)
						lmax = diff; }, Kokkos::Max<real>(err_max_local));
				// [P1-8] No fence needed: parallel_reduce has already
				// synchronised before delivering err_max_local to the host.

				// MPI global reduction for Picard convergence check
				real err_max_global = 0.0;
				gdom.timers.rtgw.update += timer_update.seconds();
				timer_comm.reset();
				SERGHEI_MPI_CHECK(MPI_Allreduce(&err_max_local, &err_max_global, 1, SERGHEI_MPI_REAL, MPI_MAX, MPI_COMM_WORLD));
				gdom.timers.rtgw.convergence += timer_comm.seconds();

				if (err_max_global < picard_tol)
				{
					// [P3-4] record the iteration at which convergence occurred
					if (picard_iter == 0) rtPicardConv1[iSpec]++;
					else if (picard_iter == 1) rtPicardConv2[iSpec]++;
					else rtPicardConv3[iSpec]++;
					break; // Converged: exit Picard loop for this species
				}

				// MPI concentration exchange and BC re-application for
				// the next Picard iteration
				timer_comm.reset();
				rtgmpi.mpi_sendrecv_rt_species(rt.c, gdom, par, iSpec, 1);
				Kokkos::fence();
				gdom.timers.rtgw.comm += timer_comm.seconds();

				timer_bc.reset();
				for (int k = 0; k < rtgbc.size(); k++)
				{
					rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec, false);
				}
				Kokkos::fence();
				gdom.timers.rtgw.bc += timer_bc.seconds();
			}
			timer_comm.reset();
			rtgmpi.mpi_sendrecv_rt_species(rt.c, gdom, par, iSpec, 1);
			Kokkos::fence();
			gdom.timers.rtgw.comm += timer_comm.seconds();

			timer_bc.reset();
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec, false);
			}
			Kokkos::fence();
			gdom.timers.rtgw.bc += timer_bc.seconds();

			timer_audit.reset();
			compute_boundary_mass_fluxes(rt, gdom, rtgbc, iSpec);
			Kokkos::fence();
			gdom.timers.rtgw.boundaryAudit += timer_audit.seconds();
		}
		// Refresh boundary flux diagnostics with the converged concentrations
		// before integrating the first transport half-step.
		timer_bc.reset();
		Kokkos::deep_copy(rt.ConQss, 0.0);
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec, false);
			}
		}
		Kokkos::fence();
		gdom.timers.rtgw.bc += timer_bc.seconds();
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
		timer_audit.reset();
		const std::vector<real> mass_after_half_1 =
			rtint.computeLocalTotalMasses(rt, gdom, gw.wc_new);
		rtint.accumulateTransportFluxes(gw, rt, gdom, rtgbc, dt_original / 2.0,
										&mass_before_half_1, &mass_after_half_1, 1);
#else
		timer_audit.reset();
		rtint.accumulateTransportFluxes(gw, rt, gdom, rtgbc, dt_original / 2.0);
#endif
		Kokkos::fence();
		gdom.timers.rtgw.massAudit += timer_audit.seconds();

		// ============================================================
		// Step 2: Chemical Reaction — dt
		// ============================================================
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve(): === Step 2: Reaction (dt) ===" << std::endl;
#endif
		rt.dt = dt_original;

		// Record liquid mass before reaction (for mass balance diagnostics)
		timer_audit.reset();
		if (rtint.hasReaction())
			rtint.recordLiquidMassBefore(rt, gw, gdom);
		Kokkos::fence();
		gdom.timers.rtgw.reactAudit += timer_audit.seconds();

		// Execute reaction integration kernel
		timer_react.reset();
		reaction_integration_kernel(rt, gw, gdom, par);
		Kokkos::fence();
		gdom.timers.rtgw.reactKernel += timer_react.seconds();

		// Immediately after the reaction step ends, update the solid-phase concentration in real time based on the latest liquid-phase concentration!<\/h>
		// Otherwise, the retained solid-phase mass will desorb again when the next Transport step begins, causing ‘mass revival’
		timer_update.reset();
		Kokkos::parallel_for("Update_Solid_After_Reaction", gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			for (int iSpec = 0; iSpec < rt.n_mass; iSpec++) {
				int sorption_type = rt.Sorption_Type(iSpec);
				real C_liquid = rt.c(iSpec, iGlob, 1);

				if (sorption_type == 1 || sorption_type == 2 || sorption_type == 3)
				{
					rt.c_solid(iSpec, iGlob, 1) = sorption_equilibrium_concentration(rt, iSpec, C_liquid);
				}
			} });
		Kokkos::fence();
		gdom.timers.rtgw.solidUpdate += timer_update.seconds();

		// Compute reaction loss/gain rates (mass balance diagnostics)
		timer_audit.reset();
		if (rtint.hasReaction())
			rtint.finalizeReactionFlux(rt, gw, gdom, dt_original);
		Kokkos::fence();
		gdom.timers.rtgw.reactAudit += timer_audit.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Step2: reaction completed" << std::endl;
#endif
		// Update the "old" concentration arrays after reaction
		timer_update.reset();
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
    	for (int iSpec = 0; iSpec < rt.n_mass; iSpec++) {
            rt.c(iSpec, iGlob,  0) = rt.c(iSpec, iGlob, 1);
            rt.c_solid(iSpec, iGlob,  0) = rt.c_solid(iSpec, iGlob, 1);
        } });
		Kokkos::fence();
		gdom.timers.rtgw.update += timer_update.seconds();

		// ============================================================
		// Step 3: Physical Transport — dt/2
		// ============================================================
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve(): === Step 3: Transport (dt/2) ===" << std::endl;
#endif
		rt.dt = dt_original / 2.0;

		// Use mid-point as old water content and original new as target
		timer_prepare.reset();
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			gw.wc_old(iGlob) = rt.wc_mid(iGlob);
			gw.wc_new(iGlob) = rt.wc_original_new(iGlob); });
		Kokkos::fence();
		gdom.timers.rtgw.prepare += timer_prepare.seconds();
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
		timer_audit.reset();
		const std::vector<real> mass_before_half_3 =
			rtint.computeLocalTotalMasses(rt, gdom, gw.wc_old);
		Kokkos::fence();
		gdom.timers.rtgw.massAudit += timer_audit.seconds();
#endif

		// Apply boundary conditions
		timer_bc.reset();
		Kokkos::deep_copy(rt.ConQss, 0.0);
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec, false);
			}
		}
		Kokkos::fence();
		gdom.timers.rtgw.bc += timer_bc.seconds();

		// MPI concentration halo exchange
		timer_comm.reset();
		rtgmpi.mpi_sendrecv_rt(rt.c, gdom, par, rt.n_mass, 1);
		gdom.timers.rtgw.comm += timer_comm.seconds();

		// Compute dispersion tensor
		timer_disp.reset();
		dispersion_tensor(rt, gw, gdom, gmpi, par);
		Kokkos::fence();
		gdom.timers.rtgw.disp += timer_disp.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Step3: dispersion tensor computed, starting Picard iteration" << std::endl;
#endif
		// [P1-2] Upwind weights: once per half-step (see Step 1).
		timer_prepare.reset();
		RTcompute_upwind_weights(gdom, rt);
		Kokkos::fence();
		gdom.timers.rtgw.prepare += timer_prepare.seconds();
		// Per-species Picard iteration loop (Step 3)
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			// [P1-1] Picard-invariant face coefficients: once per species.
			timer_prepare.reset();
			RTprepare_species_coefficients(gdom, rt, iSpec);
			Kokkos::fence();
			gdom.timers.rtgw.prepare += timer_prepare.seconds();

			rtPicardEnsure(rt.n_mass); // [P3-4] Picard statistics
			rtPicardEpisodes[iSpec]++;
			for (int picard_iter = 0; picard_iter < picard_iter_max; picard_iter++)
			{
				rtPicardIters[iSpec]++; // [P3-4]
				// Initialise linear solver guess
				Kokkos::parallel_for("Init_Solver_X", gdom.nCell, KOKKOS_LAMBDA(int idom) {
	            int ii, jj, kk, iGlob;
	            gdom.unpackIndices(idom, kk, jj, ii);
	            iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

	            rtA.rt_x(idom) = rt.c(iSpec, iGlob, 1); });
				Kokkos::fence();

				// Assemble linear system
				timer_linsys.reset();
				RTlinear_system(gw, gdom, rtA, par, rt, rtgbc, ss, iSpec, picard_iter);
				Kokkos::fence();
				gdom.timers.rtgw.linsys += timer_linsys.seconds();

				// Solve linear system
				timer_linsol.reset();
				if (rt.rt_scheme == 1)
				{
					RTSolverGW.Gauss_Seidel(rtA);
				}
				else if (rt.rt_scheme == 2)
				{
					RTSolverGW.Jacobi_GMRES_Solve(rtA, false, "Step3", iSpec, picard_iter, 0.0, 0.0, rt.dt);
				}
				Kokkos::fence();
				gdom.timers.rtgw.linsol += timer_linsol.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
				if (par.masterproc)
					std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Step3: linear solve done for iSpec=" << iSpec << " picard_iter=" << picard_iter << std::endl;
#endif
				// Update concentrations and solid-phase
				timer_update.reset();
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
				real clipped_mass_local = 0.0;
				Kokkos::parallel_reduce("RTGW_Step3_NegativeClipMass", gdom.nCell, KOKKOS_LAMBDA(int idom, real &tmp) {
						if (rtA.rt_x(idom) < 0.0)
						{
							int ii, jj, kk;
							gdom.unpackIndices(idom, kk, jj, ii);
							const int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc
							                + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
							const real volume = gdom.dx * gdom.dy * gdom.dz(iGlob);
							tmp += -rtA.rt_x(idom) * gw.wc_new(iGlob) * volume * 1000.0;
						} }, Kokkos::Sum<real>(clipped_mass_local));
				rtint.AccumNegativeClippedMass_spec(iSpec) += clipped_mass_local;
#endif

				real err_max_local = 0.0;
				Kokkos::parallel_reduce("Update_Concentration", gdom.nCell, KOKKOS_LAMBDA(int idom, real &lmax) {
						int ii, jj, kk, iGlob;
						gdom.unpackIndices(idom, kk, jj, ii);
						iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

						real C_old_iter = rt.c(iSpec, iGlob, 1);
						rt.c(iSpec, iGlob, 1) = rtA.rt_x(idom);
						if (rt.c(iSpec, iGlob, 1) < 0.0)
							rt.c(iSpec, iGlob, 1) = 0.0;

						// Update solid-phase concentration (diagnostic)
						int sorption_type = rt.Sorption_Type(iSpec);
						real C_liquid = rt.c(iSpec, iGlob, 1);

						if (sorption_type == 1 || sorption_type == 2 || sorption_type == 3)
						{
							rt.c_solid(iSpec, iGlob, 1) = sorption_equilibrium_concentration(rt, iSpec, C_liquid);
						}
						else if (sorption_type == 4)
						{
							real dt = rt.dt;
							real Kd = rt.Kd(iSpec);
							real beta_k = rt.Beta(iSpec);
							real lambda_2 = rt.Lambda_2(iSpec);
							real rho_b = rt.rho_b;
							real S_old = rt.c_solid(iSpec, iGlob, 0);
							if (Kd > 0.0)
							{
								real denom = rho_b / dt + beta_k / Kd + lambda_2 * rho_b;
								rt.c_solid(iSpec, iGlob, 1) = (rho_b / dt * S_old + beta_k * C_liquid) / denom;
							}
						}

						real diff = fabs(rt.c(iSpec, iGlob, 1) - C_old_iter);
						if (diff > lmax)
							lmax = diff; }, Kokkos::Max<real>(err_max_local));
				Kokkos::fence();

				// MPI global reduction for Picard convergence
				real err_max_global = 0.0;
				gdom.timers.rtgw.update += timer_update.seconds();
				timer_comm.reset();
				SERGHEI_MPI_CHECK(MPI_Allreduce(&err_max_local, &err_max_global, 1, SERGHEI_MPI_REAL, MPI_MAX, MPI_COMM_WORLD));
				gdom.timers.rtgw.convergence += timer_comm.seconds();

				if (err_max_global < picard_tol)
				{
					// [P3-4] record the iteration at which convergence occurred
					if (picard_iter == 0) rtPicardConv1[iSpec]++;
					else if (picard_iter == 1) rtPicardConv2[iSpec]++;
					else rtPicardConv3[iSpec]++;
					break; // Converged: exit Picard loop
				}

				// MPI exchange and BC re-application for next iteration
				timer_comm.reset();
				rtgmpi.mpi_sendrecv_rt_species(rt.c, gdom, par, iSpec, 1);
				Kokkos::fence();
				gdom.timers.rtgw.comm += timer_comm.seconds();

				timer_bc.reset();
				for (int k = 0; k < rtgbc.size(); k++)
				{
					rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec, false);
				}
				Kokkos::fence();
				gdom.timers.rtgw.bc += timer_bc.seconds();
			}
			timer_comm.reset();
			rtgmpi.mpi_sendrecv_rt_species(rt.c, gdom, par, iSpec, 1);
			Kokkos::fence();
			gdom.timers.rtgw.comm += timer_comm.seconds();

			timer_bc.reset();
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec, false);
			}
			Kokkos::fence();
			gdom.timers.rtgw.bc += timer_bc.seconds();

			timer_audit.reset();
			compute_boundary_mass_fluxes(rt, gdom, rtgbc, iSpec);
			Kokkos::fence();
			gdom.timers.rtgw.boundaryAudit += timer_audit.seconds();
		}
		// Refresh boundary flux diagnostics with the converged concentrations
		// before integrating the second transport half-step.
		timer_bc.reset();
		Kokkos::deep_copy(rt.ConQss, 0.0);
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec, false);
			}
		}
		Kokkos::fence();
		gdom.timers.rtgw.bc += timer_bc.seconds();
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
		timer_audit.reset();
		const std::vector<real> mass_after_half_3 =
			rtint.computeLocalTotalMasses(rt, gdom, gw.wc_new);
		rtint.accumulateTransportFluxes(gw, rt, gdom, rtgbc, dt_original / 2.0,
										&mass_before_half_3, &mass_after_half_3, 3);
#else
		timer_audit.reset();
		rtint.accumulateTransportFluxes(gw, rt, gdom, rtgbc, dt_original / 2.0);
#endif
		Kokkos::fence();
		gdom.timers.rtgw.massAudit += timer_audit.seconds();

		// -----------------------------------------------------------------
		// Post-processing: restore time step and synchronise old/new arrays
		// -----------------------------------------------------------------
		rt.dt = dt_original;
		timer_comm.reset();
		real dt_tmp = gdom.dt;
		SERGHEI_MPI_CHECK(MPI_Allreduce(&dt_tmp, &gdom.dt, 1, SERGHEI_MPI_REAL, MPI_MIN, MPI_COMM_WORLD));
		gdom.timers.rtgw.comm += timer_comm.seconds();

		// Final synchronisation: copy new concentrations to old arrays
		timer_postprocess.reset();
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
			{
				rt.c(iSpec, iGlob,  0) = rt.c(iSpec, iGlob, 1);
				rt.c_solid(iSpec, iGlob,  0) = rt.c_solid(iSpec, iGlob, 1);
			} });

		Kokkos::fence();
		gdom.timers.rtgw.postprocess += timer_postprocess.seconds();

		gdom.timers.rtgw.total += timer_total.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve() completed OK. Total timer=" << timer_total.seconds() << "s" << std::endl;
#endif
	}
};
#endif
