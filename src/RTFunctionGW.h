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
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
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

public:
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
			rt.aveVB(iGlob, 3) = sqrt(pow(rt.aveVB(iGlob, 0), 2) +
										  pow(rt.aveVB(iGlob, 1), 2) +
										  pow(rt.aveVB(iGlob, 2), 2)); });
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

			// --- Face i+1/2 velocity components ---
			// vx at i+1/2 (direct cell-centre value) [m/s]
			rt.q_dispersion(iGlob, 0) = rt.aveVB(iGlob, 0);
			// vy at i+1/2 (four-point average) [m/s]
			rt.q_dispersion(iGlob, 1) = 0.5 * ((rt.aveVB(iGlob - gdom.nxhc, 1) + rt.aveVB(iGlob, 1)) * 0.5 + (rt.aveVB(iGlob + 1 - gdom.nxhc, 1) + rt.aveVB(iGlob + 1, 1)) * 0.5);
			// vz at i+1/2 (four-point average) [m/s]
			rt.q_dispersion(iGlob, 2) = 0.5 * ((rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) + rt.aveVB(iGlob, 2)) * 0.5 + (rt.aveVB(iGlob + 1 - gdom.nxhc * gdom.nyhc, 2) + rt.aveVB(iGlob + 1, 2)) * 0.5);
			// Speed at i+1/2 [m/s]
			rt.q_dispersion(iGlob, 3) = sqrt(pow(rt.q_dispersion(iGlob, 0), 2) + pow(rt.q_dispersion(iGlob, 1), 2) + pow(rt.q_dispersion(iGlob, 2), 2));

			// --- Face j+1/2 velocity components ---
			// vx at j+1/2 (four-point average) [m/s]
			//! rt.q_dispersion(iGlob, 4) = 0.5 * ((rt.aveVB(iGlob, 1) + rt.aveVB(iGlob - 1, 1)) * 0.5 + (rt.aveVB(iGlob + gdom.nxhc, 1) + rt.aveVB(iGlob + gdom.nxhc - 1, 1)) * 0.5);
			rt.q_dispersion(iGlob, 4) = 0.5 * ((rt.aveVB(iGlob, 0) + rt.aveVB(iGlob - 1, 0)) * 0.5 + (rt.aveVB(iGlob + gdom.nxhc, 0) + rt.aveVB(iGlob + gdom.nxhc - 1, 0)) * 0.5);
			// vy at j+1/2 (direct cell-centre value) [m/s]
			rt.q_dispersion(iGlob, 5) = rt.aveVB(iGlob, 1);
			// vz at j+1/2 (four-point average) [m/s]
			rt.q_dispersion(iGlob, 6) = 0.5 * ((rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) + rt.aveVB(iGlob, 2)) * 0.5 + (rt.aveVB(iGlob + gdom.nxhc - gdom.nxhc * gdom.nyhc, 2) + rt.aveVB(iGlob + gdom.nxhc, 2)) * 0.5);
			// Speed at j+1/2 [m/s]
			rt.q_dispersion(iGlob, 7) = sqrt(pow(rt.q_dispersion(iGlob, 4), 2) + pow(rt.q_dispersion(iGlob, 5), 2) + pow(rt.q_dispersion(iGlob, 6), 2));

			// --- Face k+1/2 velocity components ---
			// vx at k+1/2 (four-point average) [m/s]
			rt.q_dispersion(iGlob, 8) = 0.5 * ((rt.aveVB(iGlob, 0) + rt.aveVB(iGlob + gdom.nxhc * gdom.nyhc, 0)) * 0.5 + (rt.aveVB(iGlob - 1, 0) + rt.aveVB(iGlob + gdom.nxhc * gdom.nyhc - 1, 0)) * 0.5);
			// vy at k+1/2 (four-point average) [m/s]
			rt.q_dispersion(iGlob, 9) = 0.5 * ((rt.aveVB(iGlob, 1) + rt.aveVB(iGlob - gdom.nxhc, 1)) * 0.5 + (rt.aveVB(iGlob - gdom.nxhc + gdom.nxhc * gdom.nyhc, 1) + rt.aveVB(iGlob + gdom.nxhc * gdom.nyhc, 1)) * 0.5);
			// vz at k+1/2 (direct cell-centre value) [m/s]
			rt.q_dispersion(iGlob, 10) = rt.aveVB(iGlob, 2);
			// Speed at k+1/2 [m/s]
			rt.q_dispersion(iGlob, 11) = sqrt(pow(rt.q_dispersion(iGlob, 8), 2) + pow(rt.q_dispersion(iGlob, 9), 2) + pow(rt.q_dispersion(iGlob, 10), 2)); });

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

			// --- Dxx at face i+1/2 [m^2/s] ---
			// Dxx = (alpha_T * (vy^2 + vz^2) + alpha_L * vx^2) / |v|
			rt.dcal(iSpec, iGlob, 0) = (rt.alpha_T(iSpec) * (pow(rt.q_dispersion(iGlob, 1), 2) + pow(rt.q_dispersion(iGlob, 2), 2)) + rt.alpha_L(iSpec) * pow(rt.q_dispersion(iGlob, 0), 2)) / rt.q_dispersion(iGlob, 3);

			// Dxy at face i+1/2 [m^2/s]
			rt.dcal(iSpec, iGlob, 3) = (rt.alpha_L(iSpec) - rt.alpha_T(iSpec)) * rt.q_dispersion(iGlob, 0) * rt.q_dispersion(iGlob, 1) / rt.q_dispersion(iGlob, 3);
			// Dxz at face i+1/2 [m^2/s]
			rt.dcal(iSpec, iGlob, 4) = (rt.alpha_L(iSpec) - rt.alpha_T(iSpec)) * rt.q_dispersion(iGlob, 0) * rt.q_dispersion(iGlob, 2) / rt.q_dispersion(iGlob, 3);

			//! --- Dyy at face j+1/2 [m^2/s] ---
			rt.dcal(iSpec, iGlob, 1) = (rt.alpha_T(iSpec) * (pow(rt.q_dispersion(iGlob, 4), 2) + pow(rt.q_dispersion(iGlob, 6), 2)) + rt.alpha_L(iSpec) * pow(rt.q_dispersion(iGlob, 5), 2)) / rt.q_dispersion(iGlob, 7);
			// Dyx at face j+1/2 [m^2/s]
			rt.dcal(iSpec, iGlob, 5) = (rt.alpha_L(iSpec) - rt.alpha_T(iSpec)) * rt.q_dispersion(iGlob, 5) * rt.q_dispersion(iGlob, 4) / rt.q_dispersion(iGlob, 7);
			// Dyz at face j+1/2 [m^2/s]
			rt.dcal(iSpec, iGlob, 6) = (rt.alpha_L(iSpec) - rt.alpha_T(iSpec)) * rt.q_dispersion(iGlob, 5) * rt.q_dispersion(iGlob, 6) / rt.q_dispersion(iGlob, 7);

			//! --- Dzz at face k+1/2 [m^2/s] ---
			rt.dcal(iSpec, iGlob, 2) = (rt.alpha_T(iSpec) * (pow(rt.q_dispersion(iGlob, 8), 2) + pow(rt.q_dispersion(iGlob, 9), 2)) + rt.alpha_L(iSpec) * pow(rt.q_dispersion(iGlob, 10), 2)) / rt.q_dispersion(iGlob, 11);
			// Dzx at face k+1/2 [m^2/s]
			rt.dcal(iSpec, iGlob, 7) = (rt.alpha_L(iSpec) - rt.alpha_T(iSpec)) * rt.q_dispersion(iGlob, 10) * rt.q_dispersion(iGlob, 8) / rt.q_dispersion(iGlob, 11);
			// Dzy at face k+1/2 [m^2/s]
			rt.dcal(iSpec, iGlob, 8) = (rt.alpha_L(iSpec) - rt.alpha_T(iSpec)) * rt.q_dispersion(iGlob, 10) * rt.q_dispersion(iGlob, 9) / rt.q_dispersion(iGlob, 11);


			// If at a boundary, the dispersion coefficients should be set to 0
			// to naturally enforce the boundary condition.
			// If the qz velocity is zero at the bottom boundary (hydrodynamic
			// no-flux BC), then Dzz = 0.
			// TODO: Consider combining with boundary condition modifications to
			// the dispersion tensor.

			if (rt.q_dispersion(iGlob, 10) == 0) // qz,i,j,k+1/2 == 0
			{
				rt.dcal(iSpec, iGlob, 2) = 0;
			}

			if (rt.q_dispersion(iGlob, 5) == 0) // qy,i,j+1/2,k == 0
			{
				rt.dcal(iSpec, iGlob, 1) = 0;
			}



			// Replace any NaN values in the dispersion tensor with zero
			// (can arise from 0/0 at zero-velocity faces)
			for (int i = 0; i < 9; ++i)
			{
				if (std::isnan(rt.dcal(iSpec, iGlob, i)))
				{
					rt.dcal(iSpec, iGlob, i) = 0;
				}
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
		// Fused kernel: clear coefficients + main stencil assembly +
		//               cross-dispersion terms  =>  single parallel_for
		// -----------------------------------------------------------------
		Kokkos::parallel_for("Build_Matrix_Fused", gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

			// Allocate local coefficient array in registers / fast cache;
			// initialisation to zero also clears any stale data.
			real local_coef[27] = {0.0};

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
																																										: 0.5;

			// Vertical geometric weight (used for cross-dispersion z-terms)
			rt.wz(iGlob) = 0.5;

			// --- Flux coefficients ---
			// Advection coefficients: v * dt / dx  [dimensionless Courant number]
			// (note: aveVB is already corrected to follow the coordinate axes)
			rt.c_advxx(iSpec, iGlob) = rt.aveVB(iGlob, 0) * rt.dt / gdom.dx;
			rt.c_advyy(iSpec, iGlob) = rt.aveVB(iGlob, 1) * rt.dt / gdom.dy;
			rt.c_advzz(iSpec, iGlob) = rt.aveVB(iGlob, 2) * rt.dt; // divided by dz later

			// Dispersion coefficients: D * dt / dx^2  [dimensionless]
			rt.c_difxx(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 0) * rt.dt / (gdom.dx * gdom.dx);
			rt.c_difyy(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 1) * rt.dt / (gdom.dy * gdom.dy);
			rt.c_difzz(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 2) * rt.dt / (0.5 * gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc));

			// Cross-dispersion coefficients [dimensionless]
			rt.c_difxy(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 3) * rt.dt / (2 * gdom.dx * gdom.dy);
			rt.c_difxz(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 4) * rt.dt;
			rt.c_difyx(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 5) * rt.dt / (2 * gdom.dx * gdom.dy);
			rt.c_difyz(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 6) * rt.dt;
			rt.c_difzx(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 7) * rt.dt;
			rt.c_difzy(iSpec, iGlob) = rt.dcal(iSpec, iGlob, 8) * rt.dt;

			// =================================================================
			// Mass conservation storage term with sorption and decay
			// =================================================================
			real min_conc = 1e-6;
			real Rf = 1.0;                  // Retardation factor [dimensionless]
			real intermediate_variable = 0.0; // Intermediate variable for kinetic sorption
			int sorption_type = rt.Sorption_Type(iSpec);

			// During Picard iteration, use the concentration from the
			// previous iteration to evaluate nonlinear terms.
			real C_iter = fmax(0.0, rt.c(iSpec, iGlob, 1));

			if (sorption_type == 1) // Linear equilibrium sorption
			{
				Rf = 1.0 + rt.rho_b * rt.Kd(iSpec) / gw.wc_old(iGlob);
			}
			else if (sorption_type == 2) // Freundlich isotherm sorption
			{
				if (C_iter >= min_conc) {
					Rf = 1.0 + (rt.rho_b * rt.Kf(iSpec) * rt.Nf(iSpec) * pow(C_iter, rt.Nf(iSpec) - 1.0)) / gw.wc_old(iGlob);
				}
			}
			else if (sorption_type == 3) // Langmuir isotherm sorption
			{
				real alpha = rt.Alpha_D(iSpec);
				real beta_L = rt.Beta_D(iSpec);
				Rf = 1.0 + (rt.rho_b / gw.wc_old(iGlob)) * ((alpha * beta_L) / pow(1.0 + alpha * C_iter, 2.0));
			}
			else if (sorption_type == 4) // Non-equilibrium (kinetic) sorption
			{
				real beta_k = rt.Beta(iSpec);
				real Kd = rt.Kd(iSpec);
				if (Kd > 0.0) {
					intermediate_variable = beta_k / (Kd * (rt.rho_b / rt.dt + beta_k / Kd + rt.Lambda_2(iSpec) * rt.rho_b));
				}
				Rf = 1.0;
			}

			// Apply retardation factor to the storage term (LHS diagonal and RHS)
			real mass_lhs_diag = gw.wc_new(iGlob) * Rf;
			real mass_rhs = gw.wc_old(iGlob) * rt.c(iSpec, iGlob, 0) * Rf;

			// Additional reaction source/sink terms (decay and kinetic terms)
			if (sorption_type == 4)
			{
				mass_lhs_diag += rt.dt * (rt.Beta(iSpec) + rt.Lambda_1(iSpec) * gw.wc_new(iGlob) - intermediate_variable * pow(rt.Beta(iSpec), 2.0));
				mass_rhs += intermediate_variable * rt.rho_b * rt.c_solid(iSpec, iGlob, 0); // previous solid-phase concentration
			}
			else
			{
				// Equilibrium sorption and no-sorption modes: apply only
				// first-order liquid-phase natural decay
				mass_lhs_diag += rt.dt * rt.Lambda_1(iSpec) * gw.wc_new(iGlob);
			}

			// --- Assemble main diagonal and face-neighbour coefficients ---
			// local_coef[0]: diagonal (i,j,k)
			local_coef[0] = mass_lhs_diag + rt.c_difxx(iSpec, iGlob) + rt.c_difxx(iSpec, iGlob - 1) + rt.Up_Weighting_xp(iGlob) * rt.c_advxx(iSpec, iGlob) - (1.0 - rt.Up_Weighting_xm(iGlob)) * rt.c_advxx(iSpec, iGlob - 1) + rt.c_difyy(iSpec, iGlob) + rt.c_difyy(iSpec, iGlob - gdom.nxhc) + rt.Up_Weighting_yp(iGlob) * rt.c_advyy(iSpec, iGlob) - (1.0 - rt.Up_Weighting_ym(iGlob)) * rt.c_advyy(iSpec, iGlob - gdom.nxhc) + rt.c_difzz(iSpec, iGlob) / gdom.dz(iGlob) + rt.c_difzz(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / gdom.dz(iGlob) + rt.Up_Weighting_zp(iGlob) * rt.c_advzz(iSpec, iGlob) / gdom.dz(iGlob) - (1.0 - rt.Up_Weighting_zm(iGlob)) * rt.c_advzz(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / gdom.dz(iGlob);

			// Face neighbour coefficients (moved to LHS, note signs)
			// local_coef[1]: C_{i+1}
			local_coef[1] = -(rt.c_difxx(iSpec, iGlob) - (1.0 - rt.Up_Weighting_xp(iGlob)) * rt.c_advxx(iSpec, iGlob));
			// local_coef[2]: C_{i-1}
			local_coef[2] = -(rt.c_difxx(iSpec, iGlob - 1) + rt.Up_Weighting_xm(iGlob) * rt.c_advxx(iSpec, iGlob - 1));
			// local_coef[3]: C_{j+1}
			local_coef[3] = -(rt.c_difyy(iSpec, iGlob) - (1.0 - rt.Up_Weighting_yp(iGlob)) * rt.c_advyy(iSpec, iGlob));
			// local_coef[4]: C_{j-1}
			local_coef[4] = -(rt.c_difyy(iSpec, iGlob - gdom.nxhc) + rt.Up_Weighting_ym(iGlob) * rt.c_advyy(iSpec, iGlob - gdom.nxhc));
			// local_coef[5]: C_{k+1}
			local_coef[5] = -(rt.c_difzz(iSpec, iGlob) / gdom.dz(iGlob) - (1.0 - rt.Up_Weighting_zp(iGlob)) * rt.c_advzz(iSpec, iGlob) / gdom.dz(iGlob));
			// local_coef[6]: C_{k-1}
			local_coef[6] = -(rt.c_difzz(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / gdom.dz(iGlob) + rt.Up_Weighting_zm(iGlob) * rt.c_advzz(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / gdom.dz(iGlob));

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
															rt.aveVB(iGlob, 0), rt.dt, gdom.dx, limiter_type);
					// Add divergence of flux correction to RHS:
					//   dt/dx * (F_corr_{i+1/2} - F_corr_{i-1/2})
					local_coef[7] -= rt.dt / gdom.dx * flux_corr;
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
															rt.aveVB(iGlob_minus, 0), rt.dt, gdom.dx, limiter_type);
					// Subtract from RHS (flux divergence at i-1/2)
					local_coef[7] += rt.dt / gdom.dx * flux_corr;
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
															rt.aveVB(iGlob, 1), rt.dt, gdom.dy, limiter_type);
					local_coef[7] -= rt.dt / gdom.dy * flux_corr;
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
															rt.aveVB(iGlob_minus, 1), rt.dt, gdom.dy, limiter_type);
					local_coef[7] += rt.dt / gdom.dy * flux_corr;
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
															rt.aveVB(iGlob, 2), rt.dt, gdom.dz(iGlob), limiter_type);
					local_coef[7] -= rt.dt / gdom.dz(iGlob) * flux_corr;
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
															rt.aveVB(iGlob_minus, 2), rt.dt, gdom.dz(iGlob), limiter_type);
					local_coef[7] += rt.dt / gdom.dz(iGlob) * flux_corr;
				}
			}

			// =================================================================
			// Cross-dispersion stencil entries
			//
			// The following coefficients correspond to edge and corner
			// neighbours in the 27-point stencil.  They arise from the
			// cross-derivative terms in the dispersion tensor.
			// =================================================================

			// Correction to C_{j+1} from XY cross-dispersion
			local_coef[8] = -0.5 * rt.c_difxy(iSpec, iGlob) + 0.5 * rt.c_difxy(iSpec, iGlob - 1);
			local_coef[3] += local_coef[8];
			// C_{i,j+1}: accumulated from XY and ZY cross-terms
			local_coef[8] = -0.5 * rt.c_difxy(iSpec, iGlob) + 0.5 * rt.c_difxy(iSpec, iGlob - 1) - rt.wz(iGlob) * rt.c_difzy(iSpec, iGlob) / (2 * gdom.dy * gdom.dz(iGlob)) + (1 - rt.wz(iGlob - gdom.nxhc * gdom.nyhc)) * rt.c_difzy(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dy * gdom.dz(iGlob));
			local_coef[3] += local_coef[8];
			// C_{i+1,j+1}
			local_coef[9] = -0.5 * rt.c_difxy(iSpec, iGlob) - 0.5 * rt.c_difyx(iSpec, iGlob);
			// Correction to C_{j-1} from XY cross-dispersion
			local_coef[10] = 0.5 * rt.c_difxy(iSpec, iGlob) - 0.5 * rt.c_difxy(iSpec, iGlob - 1) + rt.wz(iGlob) * rt.c_difzy(iSpec, iGlob) / (2 * gdom.dy * gdom.dz(iGlob)) - (1 - rt.wz(iGlob - gdom.nxhc * gdom.nyhc)) * rt.c_difzy(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dy * gdom.dz(iGlob));
			local_coef[4] += local_coef[10];
			// C_{i+1,j-1}
			local_coef[11] = 0.5 * rt.c_difxy(iSpec, iGlob) + 0.5 * rt.c_difyx(iSpec, iGlob - gdom.nxhc);
			// C_{i-1,j+1}
			local_coef[12] = 0.5 * rt.c_difxy(iSpec, iGlob - 1) + 0.5 * rt.c_difyx(iSpec, iGlob);
			// C_{i-1,j-1}
			local_coef[13] = -0.5 * rt.c_difxy(iSpec, iGlob - 1) - 0.5 * rt.c_difyx(iSpec, iGlob - gdom.nxhc);

			// --- Z-direction cross-dispersion contributions ---
			// Correction to C_{k+1} from XZ and YZ cross-terms
			local_coef[14] = -0.5 * rt.c_difxz(iSpec, iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + 0.5 * rt.c_difxz(iSpec, iGlob - 1) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - 0.5 * rt.c_difyz(iSpec, iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + 0.5 * rt.c_difyz(iSpec, iGlob - gdom.nxhc) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc)));
			local_coef[5] += local_coef[14];
			//! C_{i+1,k+1}
			local_coef[15] = -0.5 * rt.c_difxz(iSpec, iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - (1 - rt.wz(iGlob)) * rt.c_difzx(iSpec, iGlob) / (2 * gdom.dx * gdom.dz(iGlob));
			// Correction to C_{k-1} from XZ and YZ cross-terms
			local_coef[16] = 0.5 * rt.c_difxz(iSpec, iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - 0.5 * rt.c_difxz(iSpec, iGlob - 1) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + 0.5 * rt.c_difyz(iSpec, iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - 0.5 * rt.c_difyz(iSpec, iGlob - gdom.nxhc) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc)));
			local_coef[6] += local_coef[16];
			//! C_{i+1,k-1}
			local_coef[17] = 0.5 * rt.c_difxz(iSpec, iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + rt.wz(iGlob - gdom.nxhc * gdom.nyhc) * rt.c_difzx(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dx * gdom.dz(iGlob));
			// C_{i-1,k+1}
			local_coef[18] = 0.5 * rt.c_difxz(iSpec, iGlob - 1) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + (1 - rt.wz(iGlob)) * rt.c_difzx(iSpec, iGlob) / (2 * gdom.dx * gdom.dz(iGlob));
			//! C_{i-1,k-1}
			local_coef[19] = -0.5 * rt.c_difxz(iSpec, iGlob - 1) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - rt.wz(iGlob - gdom.nxhc * gdom.nyhc) * rt.c_difzx(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dx * gdom.dz(iGlob));
			// Correction to C_{i+1} from YX and ZX cross-terms
			local_coef[20] = -0.5 * rt.c_difyx(iSpec, iGlob) + 0.5 * rt.c_difyx(iSpec, iGlob - gdom.nxhc) - rt.wz(iGlob) * rt.c_difzx(iSpec, iGlob) / (2 * gdom.dx * gdom.dz(iGlob)) + (1 - rt.wz(iGlob - gdom.nxhc * gdom.nyhc)) * rt.c_difzx(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dx * gdom.dz(iGlob));
			local_coef[1] += local_coef[20];
			// Correction to C_{i-1} from YX and ZX cross-terms
			local_coef[21] = 0.5 * rt.c_difyx(iSpec, iGlob) - 0.5 * rt.c_difyx(iSpec, iGlob - gdom.nxhc) + rt.wz(iGlob) * rt.c_difzx(iSpec, iGlob) / (2 * gdom.dx * gdom.dz(iGlob)) - (1 - rt.wz(iGlob - gdom.nxhc * gdom.nyhc)) * rt.c_difzx(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dx * gdom.dz(iGlob));
			local_coef[2] += local_coef[21];
			// C_{j+1,k+1}
			local_coef[22] = -0.5 * rt.c_difyz(iSpec, iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - (1 - rt.wz(iGlob)) * rt.c_difzy(iSpec, iGlob) / (2 * gdom.dy * gdom.dz(iGlob));
			// C_{j+1,k-1}
			local_coef[23] = 0.5 * rt.c_difyz(iSpec, iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + rt.wz(iGlob - gdom.nxhc * gdom.nyhc) * rt.c_difzy(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dy * gdom.dz(iGlob));
			// C_{j-1,k+1}
			local_coef[24] = 0.5 * rt.c_difyz(iSpec, iGlob - gdom.nxhc) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + (1 - rt.wz(iGlob)) * rt.c_difzy(iSpec, iGlob) / (2 * gdom.dy * gdom.dz(iGlob));
			// C_{j-1,k-1}
			local_coef[25] = -0.5 * rt.c_difyz(iSpec, iGlob - gdom.nxhc) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - rt.wz(iGlob - gdom.nxhc * gdom.nyhc) * rt.c_difzy(iSpec, iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dy * gdom.dz(iGlob));

			// Final step: write all local coefficients to global memory
			for(int k = 0; k < 27; k++) {
				rt.RTcoef(idom, k) = local_coef[k];
			}
			rtA.rt_rhs(idom) = 0.0; });

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
		Kokkos::fence();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> RTlinear_system(): MPI RHS update completed" << std::endl;
#endif

		// -----------------------------------------------------------------
		// Stage 5: Fill the sparse matrix in CRS (Compressed Row Storage)
		//          format.
		//
		// Each row corresponds to one cell (idom).  Non-zero entries are
		// inserted only for neighbours that exist (boundary check).
		// The ordering follows a lexicographic sweep so that the
		// CRS row pointer (rt_ptr) advances correctly.
		// -----------------------------------------------------------------
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk;
			int irow = rtA.rt_ptr(idom);
			gdom.unpackIndices(idom, kk, jj, ii);

			// Boundary-checked insertion of non-zero stencil entries

			if (kk > 0 && jj > 0){
				rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny - gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 25);  // C_{j-1,k-1}
				irow++;
			}

	// C_{i-1,k-1}
			if (ii > 0 && kk > 0){
				rtA.rt_ind(irow) = idom -1 - gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 19);
				irow++;
			}
	// C_{k-1}
		    if (kk > 0)	{
			rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny;
			rtA.rt_val(irow) = rt.RTcoef(idom,6);
			irow++;
			}
	// C_{i+1,k-1}
			if (ii < gdom.nx-1 && kk > 0){
				rtA.rt_ind(irow) = idom + 1 -gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 17);
				irow++;
			}

	//! C_{k-1,j+1}
			if (kk > 0 && jj < gdom.ny-1){
				rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny + gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 23);
				irow++;
			}

	//! C_{j-1,i-1}
			if (jj > 0 && ii > 0){
				rtA.rt_ind(irow) = idom - gdom.nx - 1;
				rtA.rt_val(irow) = rt.RTcoef(idom, 13);
				irow++;
			}
	//! C_{j-1,i}
			if (jj > 0)	{
				rtA.rt_ind(irow) = idom - gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 4);
				irow++;
			}
	//! C_{j-1,i+1}
			if (jj > 0 && ii < gdom.nx-1){
				rtA.rt_ind(irow) = idom - gdom.nx + 1;
				rtA.rt_val(irow) = rt.RTcoef(idom, 11);
				irow++;
			}

	// C_{i-1}
			if (ii > 0)	{
				rtA.rt_ind(irow) = idom - 1;
				rtA.rt_val(irow) = rt.RTcoef(idom,2);
				irow++;
			}
	// C_{i,j,k}  (diagonal)
			rtA.rt_ind(irow) = idom;
			rtA.rt_val(irow) = rt.RTcoef(idom,0);
			irow++;
	// C_{i+1}
			if (ii < gdom.nx-1)	{
				rtA.rt_ind(irow) = idom + 1;
				rtA.rt_val(irow) = rt.RTcoef(idom,1);
				irow++;
			}

	//! C_{j+1,i-1}
			if (jj < gdom.ny-1 && ii > 0){
				rtA.rt_ind(irow) = idom + gdom.nx - 1;
				rtA.rt_val(irow) = rt.RTcoef(idom, 12);
				irow++;
			}
	//! C_{j+1,i}
			if (jj < gdom.ny-1)	{
				rtA.rt_ind(irow) = idom + gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom,3);
				irow++;
			}
	//! C_{j+1,i+1}
			if (jj < gdom.ny-1 && ii < gdom.nx-1){
				rtA.rt_ind(irow) = idom + gdom.nx + 1;
				rtA.rt_val(irow) = rt.RTcoef(idom, 9);
				irow++;
			}

	//! C_{k+1,j-1}
			if (kk < gdom.nz-1 && jj > 0){
				rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny - gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 24);
				irow++;
			}

	// C_{i-1,k+1}
			if (ii > 0 && kk < gdom.nz-1){
				rtA.rt_ind(irow) = idom -1 + gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 18);
				irow++;
			}

	// C_{k+1}
			if (kk < gdom.nz-1)	{
				rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom,5);
				irow++;
			}
	// C_{i+1,k+1}
			if (ii < gdom.nx-1 && kk < gdom.nz-1){
				rtA.rt_ind(irow) = idom + 1 + gdom.nx * gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 15);
				irow++;
			}
	//! C_{k+1,j+1}
			if (kk < gdom.nz-1 && jj < gdom.ny-1){
				rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny + gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 22);
				irow++;
			}

			rtA.rt_rhs(idom) = rt.RTcoef(idom,7); });
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> RTlinear_system() completed OK for iSpec=" << iSpec << std::endl;
#endif
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
	 *  NH4 and NO3 use molecular mass; DON and DOC use elemental mass.
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
	 *   5. Soil source:       re-aeration + DOC/DON release
	 *
	 * @param C_new      New (guess) concentrations for 5 species [mg/L].
	 * @param C_old      Old (previous time-step) concentrations [mg/L].
	 * @param dt         Time step [s].
	 * @param Resid      Output residual vector (length 5).
	 * @param rt         Reactive transport state (kinetic parameters).
	 * @param Temp_Corr  Temperature correction factor (dimensionless).
	 * @param Rf         Retardation factors for each species [dimensionless].
	 * @param theta      Volumetric water content [dimensionless, m^3/m^3].
	 */
	KOKKOS_INLINE_FUNCTION
	void calculate_residual(const real C_new[], const real C_old[], real dt, real Resid[],
							const RTStateGW &rt, real Temp_Corr, const real Rf[], real theta)
	{
		if (theta < 1e-5)
			theta = 1e-5;
		const real SMALL_CONC = 1e-15;

		// --- Physical constants and environmental parameters ---
		const real theta_s = 0.51;			 // Saturated water content [m^3/m^3]
		const real Ka_DO = rt.Ka_DO;		 // Re-aeration rate constant [1/s]
		const real DO_sat = rt.DO_sat;		 // DO saturation concentration [mg/L]
		const real DOC_eq = rt.DOC_eq;		 // DOC release equilibrium concentration [mg C/L]
		const real DON_eq = rt.DON_eq;		 // DON release equilibrium concentration [mg N/L]
		const real K_rel_DON = rt.K_rel_DON; // DON release rate [1/s]
		const real K_rel_DOC = rt.K_rel_DOC; // DOC release rate [1/s]

		// =========================================================
		// Stoichiometric yield coefficients (mass ratios, mg/mg)
		// These are pure chemical/mass-conservation quantities,
		// hence kept as hardcoded constants.
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

			// =========================================================
			// Mass reaction rates [mg/(L*s)]
			// =========================================================

			// 1. Mineralisation rate (DON -> NH4)
			real R_min = rt.Rate_Max_Min * Temp_Corr * (DON / (rt.K_Monod_DON + DON + SMALL_CONC)) * (DO / (rt.K_Monod_DO_Min + DO + SMALL_CONC));

			// 2. Nitrification rate (NH4 -> NO3)
			real R_nit = rt.Rate_Max_Nit * Temp_Corr * (NH4 / (rt.K_Monod_NH4 + NH4 + SMALL_CONC)) * (DO / (rt.K_Monod_DO + DO + SMALL_CONC));

			// 3. Denitrification rate (NO3 + DOC -> N2)
			real R_denit = rt.Rate_Max_Denit * Temp_Corr * (NO3 / (rt.K_Monod_NO3 + NO3 + SMALL_CONC)) * (DOC / (rt.K_Monod_DOC + DOC + SMALL_CONC)) * (rt.Ki_Inhib_DO / (rt.Ki_Inhib_DO + DO + SMALL_CONC));

			// 4. Aerobic respiration rate (DOC + DO -> CO2)
			real R_aer = rt.Rate_Max_Hetero * Temp_Corr * (DOC / (rt.K_Monod_DOC_Aerobic + DOC + SMALL_CONC)) * (DO / (rt.K_Monod_DO_Aerobic + DO + SMALL_CONC));

			// 5. Soil source terms (re-aeration and organic matter release)
			real R_air = Ka_DO * fmax(0.0, (theta - theta_s)) * fmax(0.0, DO_sat - DO);
			real R_rel_DOC = K_rel_DOC * fmax(0.0, DOC_eq - DOC);
			real R_rel_DON = K_rel_DON * fmax(0.0, DON_eq - DON);

			// =========================================================
			// Source/sink assembly for each species [mg/(L*s)]
			// =========================================================
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
							const RTStateGW &rt, real Temp_Corr, const real Rf[], real theta)
	{
		const real SMALL_CONC = 1e-15;

		const real theta_s = 0.51;
		const real Ka_DO = rt.Ka_DO;		 // Re-aeration rate constant [1/s]
		const real K_rel_DOC = rt.K_rel_DOC; // DOC release rate [1/s]
		const real K_rel_DON = rt.K_rel_DON; // DON release rate [1/s]

		// Stoichiometric yield coefficients (mass ratios, mg/mg)
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

		// 5. Soil source term partial derivatives
		real dRair_dDO = -Ka_DO * fmax(0.0, (theta - theta_s));

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

			// TODO: Couple with thermal reaction results in the future;
			// currently uses a fixed environmental temperature for correction.
            real T_env = 20.0;  // Future: connect to thermal module temperature [deg C]
            // real Temp_Corr = pow(rt.Temp_Coeff_Theta, T_env - rt.Opt_Temp);
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
                int sorption_type = rt.Sorption_Type(k);
				real C_liquid = rt.c(k, iGlob, 1);

                if (sorption_type == 1) {           // Linear sorption
                    Rf[k] = 1.0 + (rt.rho_b * rt.Kd(k)) / theta;
                } else if (sorption_type == 2 && C_liquid >= 1e-6) { // Freundlich
                    Rf[k] = 1.0 + (rt.rho_b * rt.Kf(k) * rt.Nf(k) * pow(C_liquid, rt.Nf(k) - 1.0)) / theta;
                } else if (sorption_type == 3) {    // Langmuir
                    real alpha = rt.Alpha_D(k);
                    real beta_L = rt.Beta_D(k);
                    Rf[k] = 1.0 + (rt.rho_b / theta) * ((alpha * beta_L) / pow(1.0 + alpha * C_liquid, 2.0));
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
                calculate_residual(C_new, C_old, dt, Resid, rt, Temp_Corr, Rf, theta);

                // Check convergence (L2 norm of residual)
                real norm_resid = 0.0;
                for (int k = 0; k < n_spec; k++) norm_resid += Resid[k] * Resid[k];
                if (sqrt(norm_resid) < NEWTON_TOLERANCE) break;

                // Compute Jacobian and solve the linear system
                calculate_jacobian(C_new, dt, Jacobian, rt, Temp_Corr, Rf, theta);

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
	inline void rt_solve(RTStateGW &rt, RTMatrix &rtA, GwState &gw, GwDomain &gdom, std::vector<RTBCGW> &rtgbc, GwMPI &gmpi, Parallel &par, RTIntegrator &rtint, type_solver &RTSolverGW, const SourceSink &ss)
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

		// Performance timers for each computational phase
		Kokkos::Timer timer_bc, timer_comm, timer_disp, timer_linsys, timer_linsol, timer_react, timer_update, timer_integrate;
		Kokkos::Timer timer_total;
		timer_total.reset();

		// [Critical Fix 3]: Synchronise water content across MPI
		// boundaries so that the storage matrix diagonal is seamless.
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve(): MPI sync water content" << std::endl;
#endif
		gmpi.mpi_sendrecv1(gw.wc_old, gdom, par);
		gmpi.mpi_sendrecv1(gw.wc_new, gdom, par);
		// Compute mid-point water content for the Strang splitting
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
		    rt.wc_original_old(iGlob) = gw.wc_old(iGlob);
		    rt.wc_original_new(iGlob) = gw.wc_new(iGlob);
		    rt.wc_mid(iGlob) = 0.5 * (gw.wc_old(iGlob) + gw.wc_new(iGlob)); });
		Kokkos::fence();

		// ============================================================
		// Step 1: Physical Transport — dt/2
		// ============================================================
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve(): === Step 1: Transport (dt/2) ===" << std::endl;
#endif
		rt.dt = dt_original / 2.0;

		// Use mid-point water content for the first half-step
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { gw.wc_new(iGlob) = rt.wc_mid(iGlob); });
		Kokkos::fence();

		// Apply boundary conditions
		timer_bc.reset();
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec);
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
		gmpi.mpi_sendrecv_rt(rt.c, gdom, par, rt.n_mass, 1);
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
			for (int picard_iter = 0; picard_iter < picard_iter_max; picard_iter++)
			{
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
					RTSolverGW.Jacobi_GMRES_Solve(rtA);
				}
				Kokkos::fence();
				gdom.timers.rtgw.linsol += timer_linsol.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
				if (par.masterproc)
					std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Step1: linear solve done for iSpec=" << iSpec << " picard_iter=" << picard_iter << " scheme=" << rt.rt_scheme << std::endl;
#endif
				// Update concentrations and solid-phase sorption
				timer_update.reset();
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

					if (sorption_type == 1)       // Linear sorption: S = Kd * C
					{
						rt.c_solid(iSpec, iGlob, 1) = rt.Kd(iSpec) * C_liquid;
					}
					else if (sorption_type == 2)  // Freundlich: S = Kf * C^Nf
					{
						rt.c_solid(iSpec, iGlob, 1) = rt.Kf(iSpec) * pow(C_liquid, rt.Nf(iSpec));
					}
					else if (sorption_type == 3)  // Langmuir: S = (alpha*beta*C) / (1 + alpha*C)
					{
						real alpha = rt.Alpha_D(iSpec);
						real beta_L = rt.Beta_D(iSpec);
						rt.c_solid(iSpec, iGlob, 1) = (alpha * beta_L * C_liquid) / (1.0 + alpha * C_liquid);
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

						real diff = fabs(rt.c(iSpec, iGlob, 1) - C_old_iter);
						if (diff > lmax)
							lmax = diff;
					} }, Kokkos::Max<real>(err_max_local));
				Kokkos::fence();

				// MPI global reduction for Picard convergence check
				real err_max_global = 0.0;
				MPI_Allreduce(&err_max_local, &err_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

				if (err_max_global < picard_tol)
				{
					break; // Converged: exit Picard loop for this species
				}

				gdom.timers.rtgw.update += timer_update.seconds();

				// MPI concentration exchange and BC re-application for
				// the next Picard iteration
				timer_comm.reset();
				gmpi.mpi_sendrecv_rt(rt.c, gdom, par, rt.n_mass, 1);
				Kokkos::fence();
				gdom.timers.rtgw.comm += timer_comm.seconds();

				timer_bc.reset();
				for (int k = 0; k < rtgbc.size(); k++)
				{
					rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec);
				}
				Kokkos::fence();
				gdom.timers.rtgw.bc += timer_bc.seconds();
			}
		}

		// ============================================================
		// Step 2: Chemical Reaction — dt
		// ============================================================
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve(): === Step 2: Reaction (dt) ===" << std::endl;
#endif
		rt.dt = dt_original;

		// Record liquid mass before reaction (for mass balance diagnostics)
		rtint.recordLiquidMassBefore(rt, gw, gdom);

		// Execute reaction integration kernel
		timer_react.reset();
		reaction_integration_kernel(rt, gw, gdom, par);
		Kokkos::fence();
		gdom.timers.rtgw.react += timer_react.seconds();

		// Compute reaction loss/gain rates (mass balance diagnostics)
		rtint.finalizeReactionFlux(rt, gw, gdom, dt_original);
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
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			gw.wc_old(iGlob) = rt.wc_mid(iGlob);
			gw.wc_new(iGlob) = rt.wc_original_new(iGlob); });
		Kokkos::fence();

		// Apply boundary conditions
		timer_bc.reset();
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec);
			}
		}
		Kokkos::fence();
		gdom.timers.rtgw.bc += timer_bc.seconds();

		// MPI concentration halo exchange
		timer_comm.reset();
		gmpi.mpi_sendrecv_rt(rt.c, gdom, par, rt.n_mass, 1);
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
		// Per-species Picard iteration loop (Step 3)
		for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
		{
			for (int picard_iter = 0; picard_iter < picard_iter_max; picard_iter++)
			{
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
					RTSolverGW.Jacobi_GMRES_Solve(rtA);
				}
				Kokkos::fence();
				gdom.timers.rtgw.linsol += timer_linsol.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
				if (par.masterproc)
					std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Step3: linear solve done for iSpec=" << iSpec << " picard_iter=" << picard_iter << std::endl;
#endif
				// Update concentrations and solid-phase
				timer_update.reset();

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

						if (sorption_type == 1)
						{
							rt.c_solid(iSpec, iGlob, 1) = rt.Kd(iSpec) * C_liquid;
						}
						else if (sorption_type == 2)
						{
							rt.c_solid(iSpec, iGlob, 1) = rt.Kf(iSpec) * pow(C_liquid, rt.Nf(iSpec));
						}
						else if (sorption_type == 3)
						{
							real alpha = rt.Alpha_D(iSpec);
							real beta_L = rt.Beta_D(iSpec);
							rt.c_solid(iSpec, iGlob, 1) = (alpha * beta_L * C_liquid) / (1.0 + alpha * C_liquid);
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

							real diff = fabs(rt.c(iSpec, iGlob, 1) - C_old_iter);
							if (diff > lmax)
								lmax = diff;
						} }, Kokkos::Max<real>(err_max_local));
				Kokkos::fence();

				// MPI global reduction for Picard convergence
				real err_max_global = 0.0;
				MPI_Allreduce(&err_max_local, &err_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

				if (err_max_global < picard_tol)
				{
					break; // Converged: exit Picard loop
				}

				gdom.timers.rtgw.update += timer_update.seconds();
				// MPI exchange and BC re-application for next iteration
				timer_comm.reset();
				gmpi.mpi_sendrecv_rt(rt.c, gdom, par, rt.n_mass, 1);
				Kokkos::fence();
				gdom.timers.rtgw.comm += timer_comm.seconds();

				timer_bc.reset();
				for (int k = 0; k < rtgbc.size(); k++)
				{
					rtgbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec);
				}
				Kokkos::fence();
				gdom.timers.rtgw.bc += timer_bc.seconds();
			}
		}

		// -----------------------------------------------------------------
		// Post-processing: restore time step and synchronise old/new arrays
		// -----------------------------------------------------------------
		rt.dt = dt_original;
		real dt_tmp = gdom.dt;
		int ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		// Final synchronisation: copy new concentrations to old arrays
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
			{
				rt.c(iSpec, iGlob,  0) = rt.c(iSpec, iGlob, 1);
				rt.c_solid(iSpec, iGlob,  0) = rt.c_solid(iSpec, iGlob, 1);
			} });

		timer_integrate.reset();
		// rtint.integrate(gw, rt, gdom, rtgbc);  // Reserved for future use
		Kokkos::fence();
		gdom.timers.rtgw.integrate += timer_integrate.seconds();

		gdom.timers.rtgw.total += timer_total.seconds();
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rt_solve() completed OK. Total timer=" << timer_total.seconds() << "s" << std::endl;
#endif
	}
};
#endif