
#ifndef _EDGES_H_
#define _EDGES_H_

#include "define.h"
#include "Parallel.h"
#include "SArray.h"
#include "Domain.h"
#include "Exchange.h"
#include "Indexing.h"
#include "Solvers.h"

class Edges
{

	Kokkos::Timer timer, timerdt;

public:
	inline void computeDeltaStateSW(State &state, Domain &dom, Exchange &exch, Parallel &par)
	{
		timer.reset();
		solve(state, dom, exch, par);
		computeTimeStepReduction(dom, state);
		dom.timers.swe.flux.total += timer.seconds();
	}

	inline void solve(State &state, Domain &dom, Exchange &exch, Parallel &par)
	{

		// Riemann solver X
		computeDeltaFluxXRoe(state, dom, par);
		// Riemann solver Y
		computeDeltaFluxYRoe(state, dom, par);
	}

	inline void computeDeltaFluxXRoe(State &state, Domain const &dom, Parallel &par)
	{

		auto ncells = dom.nCellMem;

		Kokkos::parallel_for("computeDeltaFluxXRoe", ncells, KOKKOS_LAMBDA(int iGlob) {
			int i, j;
			int id1, id2;
			unpackIndicesUniformGrid(iGlob, dom.ny + 2 * dom.hc, dom.nx + 2 * dom.hc, j, i);
			if (i > dom.hc - 2 && i < dom.nx + dom.hc && j > dom.hc - 1 && j < dom.ny + dom.hc)
			{							// note the hc-2 (first valid halo-inner wall)
				SArray<real, 3> upwM, upwP;
				SArray<real, 5> s1, s2; // 3 sw variables plus z and roughness

				id1 = iGlob; // j*(dom.nx+2*hc)+i
				id2 = j * (dom.nx + 2 * dom.hc) + i + 1;

				s1(idH) = state.h(id1);
				s2(idH) = state.h(id2);

				bool nodata = state.isnodata(id1) || state.isnodata(id2);

				if ((s1(idH) > 0. || s2(idH) > 0.) && // only wet-wet
					!(nodata) &&					  // no data-nodata
					!(dom.iW && i == dom.hc - 1) &&	  // no data - outer halo
					!(dom.iE && i == dom.nx + dom.hc - 1))
				{ // no data - outer halo

					s1(idHU) = state.hu(id1);
					s2(idHU) = state.hu(id2);
					s1(idHV) = state.hv(id1);
					s2(idHV) = state.hv(id2);
					s1(idZ) = state.z(id1);
					s2(idZ) = state.z(id2);
					s1(idR) = state.roughness(id1);
					s2(idR) = state.roughness(id2);

					Solver solver;
					solver.roe(s1, s2, dom.dt, dom.dx(), 1, 0);
					
					#if SERGHEI_SURFACE_TRANSPORT
					real interface_flux = 0.0;
					solver.roeSolver(s1, s2, upwM, upwP, interface_flux, dom.dt, dom.dx(), 1, 0);
					// Save the x-direction interface discharge for scalar-transport coupling
					// Stored at id1: the flux across the face between id1 (left) and id2 (right)
                    state.InterfaceFlux_x(id1) = interface_flux;
					#endif

					state.dsw0(id1) = solver.upwM(0);
					state.dsw0(id1 + ncells) = solver.upwM(1);
					state.dsw0(id1 + 2 * ncells) = solver.upwM(2);

					state.dsw1(id2) = solver.upwP(0);
					state.dsw1(id2 + ncells) = solver.upwP(1);
					state.dsw1(id2 + 2 * ncells) = solver.upwP(2);

#if SERGHEI_SCALAR_TRANSPORT
					state.ade.upwinding(id1, id2, s1(idH), s2(idH), solver.numFlux);
#if SERGHEI_SCALAR_DIFFUSION
					state.ade.edgeDiffusion(id1, id2, s1(idH), s2(idH), solver, dom.dt, dom.dx(), 1, 0);
#endif
#endif

#if SERGHEI_SEDIMENT_TRANSPORT && SERGHEI_UPWIND_BED
					state.sediment.upwinding(id1, id2, s1(idZ), s2(idZ), dom.dt, dom.dx(), state.ade);
#endif
				}
				else
				{
					#if SERGHEI_SURFACE_TRANSPORT
					// Safety: dry or boundary cells get zero interface flux
                    state.InterfaceFlux_x(id1) = 0.0;
					#endif

				}
			} });
	}

	inline void computeDeltaFluxYRoe(State &state, Domain const &dom, Parallel &par)
	{

		auto ncells = dom.nCellMem;

		Kokkos::parallel_for("computeDeltaFluxXRoe", ncells, KOKKOS_LAMBDA(int iGlob) {
			int i, j;
			int id1, id2;
			unpackIndicesUniformGrid(iGlob, dom.ny + 2 * dom.hc, dom.nx + 2 * dom.hc, j, i);
			if (i > dom.hc - 1 && i < dom.nx + dom.hc && j > dom.hc - 2 && j < dom.ny + dom.hc)
			{							// note the hc-2 (first valid halo-inner wall)
				SArray<real, 3> upwM, upwP;
				SArray<real, 5> s1, s2; // 3 sw variables plus z and roughness

				id1 = iGlob; // j*(dom.nx+2*hc)+i
				id2 = (j + 1) * (dom.nx + 2 * dom.hc) + i;

				s1(idH) = state.h(id1);
				s2(idH) = state.h(id2);

				bool nodata = state.isnodata(id1) || state.isnodata(id2);

				if ((s1(idH) > 0. || s2(idH) > 0.) && // only wet-wet
					!(nodata) &&					  // no data-nodata
					!(dom.iN && j == dom.hc - 1) &&	  // no data - outer halo
					!(dom.iS && j == dom.ny + dom.hc - 1))
				{ // no data - outer halo

					s1(idHU) = state.hu(id1);
					s2(idHU) = state.hu(id2);
					s1(idHV) = state.hv(id1);
					s2(idHV) = state.hv(id2);
					s1(idZ) = state.z(id1);
					s2(idZ) = state.z(id2);
					s1(idR) = state.roughness(id1);
					s2(idR) = state.roughness(id2);

					Solver solver;
					solver.roe(s1, s2, dom.dt, dom.dx(), 0, -1);

					#if SERGHEI_SURFACE_TRANSPORT
					real interface_flux = 0.0;
					solver.roeSolver(s1, s2, upwM, upwP, interface_flux, dom.dt, dom.dx(), 0, -1);
					// Save the y-direction interface discharge for scalar-transport coupling
					// Stored at id1: the flux across the face between id1 (upper) and id2 (lower)
					state.InterfaceFlux_y(id1) = interface_flux;
					#endif

					// note that we have sum to not overwrite the x-contributions
					state.dsw0(id1) += solver.upwM(0);
					state.dsw0(id1 + ncells) += solver.upwM(1);
					state.dsw0(id1 + 2 * ncells) += solver.upwM(2);

					state.dsw1(id2) += solver.upwP(0);
					state.dsw1(id2 + ncells) += solver.upwP(1);
					state.dsw1(id2 + 2 * ncells) += solver.upwP(2);

#if SERGHEI_SCALAR_TRANSPORT
					state.ade.upwinding(id1, id2, s1(idH), s2(idH), solver.numFlux);
#if SERGHEI_SCALAR_DIFFUSION
					state.ade.edgeDiffusion(id1, id2, s1(idH), s2(idH), solver, dom.dt, dom.dx(), 0, -1);
#endif
#endif

#if SERGHEI_SEDIMENT_TRANSPORT && SERGHEI_UPWIND_BED
					state.sediment.upwinding(id1, id2, s1(idZ), s2(idZ), dom.dt, dom.dx(), state.ade);
#endif
				}
				else{
					#if SERGHEI_SURFACE_TRANSPORT
					// Safety: dry or boundary cells get zero interface flux
                    state.InterfaceFlux_y(id1) = 0.0;
					#endif
				}
			} });
	}

	inline void computeTimeStepReduction(Domain &dom, State &state)
	{
		real dtloc = dom.dt;
		Kokkos::fence();
		timerdt.reset();
		Kokkos::parallel_reduce("computeTimeStepReduction", dom.nCell, KOKKOS_LAMBDA(int iGlob, real &dt) {
			int ii = dom.getIndex(iGlob);

			dt = min(dt, dom.dt);
			real h = state.h(ii);
			real dh = state.dsw0(ii) + state.dsw1(ii);
			if (dh > TOLDRY)
			{ // only positive dh can make negative water depths
				dt = min(dt, (h + TOLDRY) * dom.dx() / dh);
			} }, Kokkos::Min<real>(dtloc));

		Kokkos::fence();

		dom.timers.swe.flux.dt += timerdt.seconds();
		dom.dt = dtloc;
		if (dom.nsubdom > 1)
		{
			timerdt.reset();
			SERGHEI_MPI_CHECK(MPI_Allreduce(&dtloc, &dom.dt, 1, SERGHEI_MPI_REAL, MPI_MIN, MPI_COMM_WORLD));
			dom.timers.swe.flux.mpi += timerdt.seconds();
		}
#if SERGHEI_DEBUG_DT
		std::cout << "time = " << dom.etime << "\tdt_neg = " << dom.dt << std::endl;
#endif
	}
};

#endif
