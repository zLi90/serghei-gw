#ifndef _GW_INTEGRATOR_H_
#define _GW_INTEGRATOR_H_

#if SERGHEI_RE_MODEL // [Code2] Macro already renamed from SERGHEI_SUBSURFACE_MODEL to SERGHEI_RE_MODEL

#include "GwState.h"
#include "SourceSink.h"
#include "GwDomain.h"
#include "GwBC.h"
#include "Indexing.h"

class GwIntegrator
{

	Kokkos::Timer timer;

public:
	real Vtot, Vtot_glob;
	real Vexch, Vexch_glob;

	SourceSink *ss;
	GwState *gw;
	GwDomain *gdom;
	std::vector<GwBC> *gwbc;
	std::vector<GwSS> *gwss;

	int ncellsBC, ncellsBC_glob, ncellsIT, ncellsIT_glob;
	real QinBC, QoutBC, QinSS, QoutSS;
	real QinBC_glob, QoutBC_glob, QinSS_glob, QoutSS_glob;

	// [Code1] Extended source/sink tracking for soil evaporation, root transpiration, and tile drainage
	real QoutSS_SoilEvap, QoutSS_RootTransp, QoutSS_TileDrainage;
	real QoutSS_SoilEvap_glob, QoutSS_RootTransp_glob, QoutSS_TileDrainage_glob;

	void initialize(GwState &gw_, GwDomain &gdom_, std::vector<GwBC> &gwbc_, std::vector<GwSS> &gwss_)
	{
		gw = &gw_;
		gdom = &gdom_;
		gwss = &gwss_;
		gwbc = &gwbc_;
	}

	void integrate(GwState const &gw, GwDomain const &gdom, std::vector<GwBC> &gwbc, std::vector<GwSS> &gwss)
	{
		int ierr = 0;
		// get total volume
		// [Code2] Uses gdom.hc (properly scoped) instead of bare hc
		Vtot = 0;
		Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
	  			int ii, jj, kk, iGlob;
            	gdom.unpackIndices(idx, kk, jj, ii);
            	iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
            	tmp += gw.wc(iGlob,1) * gdom.dx * gdom.dx * gdom.dz(iGlob); }, Kokkos::Sum<real>(Vtot));
		Kokkos::fence();
		Vtot_glob = 0.0;
		ierr = MPI_Allreduce(&Vtot, &Vtot_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// [Code1] Root zone water content calculation for crop growth model coupling
		// Computes average water content in the root zone for each source/sink domain
		// #if CROP_GROWTH_MODEL
		for (int k = 0; k < gwss.size(); k++)
		{
			real wc_sum = 0.0;	 // sum of water content in root zone
			int ncells_root = 0; // number of root zone cells

			// get the layer index corresponding to maximum root depth
			int k_root_max = gwss[k].k_max_root;

			// iterate over cells from layer 0 to k_root_max, computing root zone average wc
			for (int kk = 0; kk <= k_root_max; kk++)
			{
				for (int jj = 0; jj < gdom.ny; jj++)
				{
					for (int ii = 0; ii < gdom.nx; ii++)
					{
						int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
						wc_sum += gw.wc(iGlob, 0); // accumulate water content
						ncells_root++;			   // count root zone cells
					}
				}
			}

			// compute root zone average water content
			if (ncells_root > 0)
			{
				gwss[k].wc_root_zone = wc_sum / ncells_root;
			}
			else
			{
				gwss[k].wc_root_zone = 0.0; // no root zone cells, average wc = 0
			}
		}
		// #endif
		// [Code1] End root zone water content calculation

		// get surface-subsurface exchange rate [m3/s]
		// [Code2] Base structure preserved
		Vexch = 0;
#if SERGHEI_SWE_MODEL
		Kokkos::parallel_reduce(gdom.ny * gdom.nx, KOKKOS_LAMBDA(int idx, real &tmp) {
		        int ii, jj, iGlob;
		        unpackIndicesUniformGrid(idx, gdom.ny, gdom.nx, jj, ii);
		        iGlob = jj*gdom.nx + ii;
		        tmp += gw.qss(iGlob) * gdom.dx * gdom.dx; }, Kokkos::Sum<real>(Vexch));
		Kokkos::fence();
#endif
		Vexch_glob = 0.0;
		ierr = MPI_Allreduce(&Vexch, &Vexch_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// boundary flow
		// [Code2] Base boundary flow logic preserved
		QinBC = 0.0;
		QoutBC = 0.0;
		for (int k = 0; k < gwbc.size(); k++)
		{
			int _ncellsBC;
			MPI_Allreduce(&(gwbc[k].ncellsBC), &_ncellsBC, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			ncellsBC += _ncellsBC;
			QinBC += gwbc[k].Qinflow;
			QoutBC += gwbc[k].Qoutflow;
		}
		QinBC_glob = 0.0;
		QoutBC_glob = 0.0;
		ncellsBC_glob = 0;
		ierr = MPI_Allreduce(&QinBC, &QinBC_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		ierr = MPI_Allreduce(&QoutBC, &QoutBC_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		ierr = MPI_Allreduce(&ncellsBC, &ncellsBC_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// source/sink terms
		// [Merge] Base structure from Code2, extended with Code1's soil evaporation / root transpiration / tile drainage tracking
		QinSS = 0.0;
		QoutSS = 0.0;
		// [Code1] Initialize extended source/sink tracking variables
		QoutSS_SoilEvap = 0.0;
		QoutSS_RootTransp = 0.0;
		QoutSS_TileDrainage = 0.0;
		for (int k = 0; k < gwss.size(); k++)
		{
			int _ncellsIT;
			MPI_Allreduce(&(gwss[k].ncellsIT), &_ncellsIT, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			ncellsIT += _ncellsIT;
			QinSS += gwss[k].Qinflow;
			QoutSS += gwss[k].Qoutflow;
			// [Code1] Accumulate soil evaporation, root transpiration, and tile drainage outflows
			QoutSS_SoilEvap += gwss[k].Qoutflow_SoilEvap;
			QoutSS_RootTransp += gwss[k].Qoutflow_RootTransp;
			QoutSS_TileDrainage += gwss[k].Qoutflow_TileDrainage;
		}
		QinSS_glob = 0.0;
		QoutSS_glob = 0.0;
		ncellsIT_glob = 0;
		// [Code1] Initialize extended global source/sink tracking variables
		QoutSS_SoilEvap_glob = 0.0;
		QoutSS_RootTransp_glob = 0.0;
		QoutSS_TileDrainage_glob = 0.0;
		ierr = MPI_Allreduce(&QinSS, &QinSS_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		ierr = MPI_Allreduce(&QoutSS, &QoutSS_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		// [Code1] MPI_Allreduce for extended source/sink terms
		ierr = MPI_Allreduce(&QoutSS_SoilEvap, &QoutSS_SoilEvap_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		ierr = MPI_Allreduce(&QoutSS_RootTransp, &QoutSS_RootTransp_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		ierr = MPI_Allreduce(&QoutSS_TileDrainage, &QoutSS_TileDrainage_glob, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
		ierr = MPI_Allreduce(&ncellsIT, &ncellsIT_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
};

#endif

#endif