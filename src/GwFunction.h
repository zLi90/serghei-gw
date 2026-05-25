#ifndef _GW_FUNCTION_H_
#define _GW_FUNCTION_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwBC.h"
#include "GwDomain.h"
#include "GwMatrix.h"
#include "GwMPI.h"
#include "GwSolver.h"
#include "GwState.h"
#include "GwIntegrator.h"
#include "State.h"
#include "SourceSink.h"
#include "FileIO.h"
#include <set>
#include <math.h>
#if SERGHEI_WAVE_MODEL
#include "WaveBoundary.h"
#endif
// [CODE1] Crop growth model (WOFOST) coupling support
#if CROP_GROWTH_MODEL
#include "./cropsrc/Wofost72.h"
#endif

class GwFunction   {

private:
	Kokkos::Timer timer;
	Kokkos::Timer timer2;

public:
    /* --------------------------------------------------
        Top-level PCA solver
    -------------------------------------------------- */
    // [CODE1] CROP_GROWTH_MODEL overload: adds Wofost72 parameter for crop coupling
#if CROP_GROWTH_MODEL
	template <typename execution_space, typename type_solver>
    inline void pca_solve(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc,
            GwMatrix &A, type_solver &gsolver, std::vector<GwSS> &gss, GwMPI &gmpi, GwIntegrator &gint, Parallel &par, Wofost72 &wofost
#if SERGHEI_WAVE_MODEL
                , WaveBoundaryModel *waveBoundary = nullptr, WaveBoundaryState *waveState = nullptr, SourceSinkData *swss = nullptr, Domain const *dom = nullptr, FileIO const *io = nullptr
#endif
            )  {
#else
	template <typename execution_space, typename type_solver>
    inline void pca_solve(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc,
            GwMatrix &A, type_solver &gsolver, std::vector<GwSS> &gss, GwMPI &gmpi, GwIntegrator &gint, Parallel &par
#if SERGHEI_WAVE_MODEL
                , WaveBoundaryModel *waveBoundary = nullptr, WaveBoundaryState *waveState = nullptr, SourceSinkData *swss = nullptr, Domain const *dom = nullptr, FileIO const *io = nullptr
#endif
            )  {
#endif
        int iter, ierr=1;
        real dt_tmp;

		timer.reset();
		timer2.reset();
#if SERGHEI_WAVE_MODEL
        if (waveBoundary != nullptr && waveState != nullptr && swss != nullptr && dom != nullptr && io != nullptr && dom->isWave)
        {
            waveBoundary->update_wave_boundary(*waveState, *swss, *dom, gdom, *io);
            waveBoundary->apply_wave_head_to_bc(*waveState, gbc);
        }
#endif

// printf(" START: h = %f\n", gw.h(200,1));
		for (int k = 0; k < gbc.size(); k++) {
            gbc[k].applyHBC(gw, gdom, par);
        }
		// [CODE2] Hierarchical timer system (re namespace)
		gdom.timers.re.gwBC += timer.seconds();

        face_conductivity(gw, gdom, gbc, gmpi, par);

		timer.reset();
        // [CODE1] CROP_GROWTH_MODEL: linear_system with Wofost72 parameter
#if CROP_GROWTH_MODEL
        linear_system(gw, gdom, gbc, gss, A, par, wofost);
#else
        linear_system(gw, gdom, gbc, gss, A, par);
#endif
		gdom.timers.re.gwlinsys += timer.seconds();

        timer.reset();
        #if SERGHEI_KOKKOSKERNELS_SOLVER
        gsolver.kkpcg(A);
        #else
        iter = gsolver.cg(A, gdom);
        #endif
        gdom.timers.re.gwlinsol += timer.seconds();

        // [CODE2] Use gdom.hc instead of bare hc
        Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
            gw.h(iGlob,1) = A.x(idom);
        });

        // printf(" SOLVE: h = %f\n", gw.h(200,1));

		timer.reset();
        gmpi.mpi_sendrecv(gw.h, gdom, par);
		gdom.timers.re.gwMPI += timer.seconds();

		timer.reset();
        for (int k = 0; k < gbc.size(); k++) {
            gbc[k].applyHBC(gw, gdom, par);
        }
		gdom.timers.re.gwBC += timer.seconds();

        face_conductivity(gw, gdom, gbc, gmpi, par);

        // [CODE2] Original face_flux
        face_flux(gw, gdom, gbc, gmpi, par);
        // [CODE1] face_flux_new for solute transport (gw.q_new used by RTFunctionGW)
        face_flux_new(gw, gdom, gbc, gmpi, par);

		timer.reset();
        // [CODE1] CROP_GROWTH_MODEL: update_wc with Wofost72 parameter
#if CROP_GROWTH_MODEL
        update_wc(gw, gdom, gss, wofost);
#else
        update_wc(gw, gdom, gss);
#endif
		gdom.timers.re.gwUpdateWC += timer.seconds();

		timer.reset();
        gmpi.mpi_sendrecv(gw.h, gdom, par);
        gmpi.mpi_sendrecv(gw.wc, gdom, par);
		gdom.timers.re.gwMPI += timer.seconds();

		timer.reset();
        for (int k = 0; k < gbc.size(); k++) {
            gbc[k].applyHBC(gw, gdom, par);
        }
		gdom.timers.re.gwBC += timer.seconds();

        dt_waco(gw, gdom);
        dt_tmp = gdom.dt;
        ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        // [MERGED] Combine Code2's gdom.hc with Code1's wc_new/wc_old tracking for RT equation
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            // [CODE1] Track wc_new/wc_old for reactive transport (solute transport) equation
            gw.wc_new(iGlob) = gw.wc(iGlob, 1);  // wc(n+1) for RT equation
            gw.wc_old(iGlob) = gw.wc(iGlob, 0);  // wc(n) for RT equation
            // [CODE2] Standard time level update
            gw.h(iGlob,0) = gw.h(iGlob,1);
            gw.wc(iGlob,0) = gw.wc(iGlob,1);
        });

		timer.reset();
        gint.integrate(gw, gdom, gbc, gss);
		gdom.timers.re.gwIntegrate += timer.seconds();

		gdom.timers.re.gw += timer2.seconds();

        // printf(" FINAL: h = %f\n", gw.h(200,1));
        // printf(" ----- \n\n");
    }

    /* --------------------------------------------------
        Top-level Picard solver
    -------------------------------------------------- */
    // [CODE1] CROP_GROWTH_MODEL overload: adds Wofost72 parameter for crop coupling
#if CROP_GROWTH_MODEL
	template <typename execution_space, typename type_solver>
    inline void picard_solve(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc,
            GwMatrix &A, type_solver &gsolver, std::vector<GwSS> &gss, GwMPI &gmpi, GwIntegrator &gint, Parallel &par, Wofost72 &wofost
#if SERGHEI_WAVE_MODEL
            , WaveBoundaryModel *waveBoundary = nullptr, WaveBoundaryState *waveState = nullptr, SourceSinkData *swss = nullptr, Domain const *dom = nullptr, FileIO const *io = nullptr
#endif
            )  {
#else
	template <typename execution_space, typename type_solver>
    inline void picard_solve(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc,
            GwMatrix &A, type_solver &gsolver, std::vector<GwSS> &gss, GwMPI &gmpi, GwIntegrator &gint, Parallel &par
#if SERGHEI_WAVE_MODEL
            , WaveBoundaryModel *waveBoundary = nullptr, WaveBoundaryState *waveState = nullptr, SourceSinkData *swss = nullptr, Domain const *dom = nullptr, FileIO const *io = nullptr
#endif
            )  {
#endif
        int iter, iter_cg, iter_max = 50, ierr=1;
        real eps_diff = 1.0, eps_old = 1.0, eps = 1.0, eps_diff_tmp = 1.0, eps_old_emp = 1.0, eps_tmp = 1.0;
		real eps_min = 1e-5, dt_tmp;
		timer.reset();
		timer2.reset();
#if SERGHEI_WAVE_MODEL
        if (waveBoundary != nullptr && waveState != nullptr && swss != nullptr && dom != nullptr && io != nullptr && dom->isWave)
        {
            waveBoundary->update_wave_boundary(*waveState, *swss, *dom, gdom, *io);
            waveBoundary->apply_wave_head_to_bc(*waveState, gbc);
        }
#endif
        for (int k = 0; k < gbc.size(); k++) {
            gbc[k].applyHBC(gw, gdom, par);
        }
		// [CODE2] Hierarchical timer system
		gdom.timers.re.gwBC += timer.seconds();

        face_conductivity(gw, gdom, gbc, gmpi, par);

        iter = 0;
        while (iter < iter_max && eps_diff/eps_old > eps_min && eps > eps_min) {

			timer.reset();
            // [CODE1] CROP_GROWTH_MODEL: linear_system with Wofost72 parameter
#if CROP_GROWTH_MODEL
            linear_system(gw, gdom, gbc, gss, A, par, wofost);
#else
            linear_system(gw, gdom, gbc, gss, A, par);
#endif
			gdom.timers.re.gwlinsys += timer.seconds();

            timer.reset();
            #if SERGHEI_KOKKOSKERNELS_SOLVER
            gsolver.kkpcg(A);
            #else
            iter = gsolver.cg(A, gdom);
            #endif
            gdom.timers.re.gwlinsol += timer.seconds();

            // [CODE2] Use gdom.hc
            Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                gw.h(iGlob,0) = gw.h(iGlob,1);
                gw.h(iGlob,1) = A.x(idom);
            });

            face_conductivity(gw, gdom, gbc, gmpi, par);

            // [CODE2] Original face_flux
            face_flux(gw, gdom, gbc, gmpi, par);
            // [MERGED] Add face_flux_new call for solute transport consistency (PCA solver already has this)
            face_flux_new(gw, gdom, gbc, gmpi, par);

			eps_old = eps_tmp;
			eps_tmp = get_eps(gw, gdom);
			eps_diff_tmp = myfabs(eps_old - eps_tmp);
			MPI_Allreduce(&eps_tmp, &eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&eps_diff_tmp, &eps_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			timer.reset();
            // [CODE1] CROP_GROWTH_MODEL: update_wc with Wofost72 parameter
#if CROP_GROWTH_MODEL
            update_wc(gw, gdom, gss, wofost);
#else
            update_wc(gw, gdom, gss);
#endif
			gdom.timers.re.gwUpdateWC += timer.seconds();

			// printf("    > RANK -%d- : Picard loop %d completed : epsOLD=%f, eps=%f, epsDIFF=%f\n",par.myrank,iter,eps_old,eps,eps_diff);
            iter += 1;
        }
        // printf("    > Picard loop converges in %d iterations with eps = %f, %f\n",iter,eps,eps_diff);

        dt_iter(gw, gdom, iter);
        dt_tmp = gdom.dt;
        ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        // [MERGED] Combine Code2's gdom.hc with Code1's wc_new/wc_old tracking for RT equation
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            // [CODE1] Track wc_new/wc_old for reactive transport (solute transport) equation
            // This was a TODO in Code1's original Picard solver, now completed in merge
            gw.wc_new(iGlob) = gw.wc(iGlob, 1);  // wc(n+1) for RT equation
            gw.wc_old(iGlob) = gw.wc(iGlob, 0);  // wc(n) for RT equation
            // [CODE2] Standard time level update
            gw.h(iGlob,0) = gw.h(iGlob,1);
            gw.wc(iGlob,0) = gw.wc(iGlob,1);
        });

		timer.reset();
        gint.integrate(gw, gdom, gbc, gss);
		gdom.timers.re.gwIntegrate += timer.seconds();

		gdom.timers.re.gw += timer2.seconds();
    }

    /* --------------------------------------------------
        Get face conductivity
    -------------------------------------------------- */
    // [MERGED] Identical logic in both branches, using Code2's gdom.hc and timer system
    inline void face_conductivity(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc, GwMPI &gmpi, Parallel &par)	{
        // Initialize K to zero (this also set boundary K=0 by default)
		timer.reset();
        Kokkos::parallel_for( gdom.nCellMem , KOKKOS_LAMBDA(int iGlob) {
            gw.k(iGlob,0) = 0.0; gw.k(iGlob,1) = 0.0;   gw.k(iGlob,2) = 0.0;
        });
        // Get relatively permeability at cell centers
        Kokkos::parallel_for( gdom.nCellMem , KOKKOS_LAMBDA(int iGlob) {
            int ii, jj, kk, ivg;
            real s, alpha, n, m, wcm, wcr, wcs, nume, deno;
            gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
            ivg = gw.soilID(iGlob) * NVG;
            wcs = gw.vgTable(ivg+2);
            wcr = gw.vgTable(ivg+3);
            n = gw.vgTable(ivg+4);
            alpha = gw.vgTable(ivg+6);
            m = 1.0 - 1.0 / n;
            wcm = wcr + (wcs-wcr)*mypow((1.0 + mypow(myfabs(gdom.aev)*alpha,n)), m);
            s = mypow(1.0 + mypow(myfabs(alpha*gw.h(iGlob,1)), n), -m);
            nume = 1.0-mypow(1.0-mypow(s*(wcs-wcr)/(wcm-wcr),1.0/m),m);
            deno = 1.0-mypow(1.0-mypow((wcs-wcr)/(wcm-wcr),1.0/m),m);
            if (deno == 0.0)    {gw.k(iGlob,3) = 1.0;}
            else {gw.k(iGlob,3) = mypow(s,0.5) * mypow(nume/deno, 2.0);}
        	gw.k(iGlob,3) = mypow(s,0.5) * mypow(1-mypow(1-mypow(s,1.0/m),m), 2.0);
            if (gw.k(iGlob,3) > 1.0 | gw.h(iGlob,1) >= gdom.aev)	{gw.k(iGlob,3) = 1.0;}
            // set no data cells impermeable
            if (gdom.isnodata(iGlob) == 1)  {gw.k(iGlob,3) = 0.0;}
        });
        // Get K on interior cell faces
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob, ivg, ivgx, ivgy, ivgz, ivgback;
            real ks, ksx, ksy, ksz, ksback;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
            ivg = gw.soilID(iGlob) * NVG;                       ks = gw.vgTable(ivg);
            ivgx = gw.soilID(iGlob+1) * NVG;                    ksx = gw.vgTable(ivgx);
            ivgy = gw.soilID(iGlob+gdom.nxhc) * NVG;            ksy = gw.vgTable(ivgy);
            ivgz = gw.soilID(iGlob+gdom.nxhc*gdom.nyhc) * NVG;  ksz = gw.vgTable(ivgz);
            // Kx
			if (ii == 0)	{
				gw.k(iGlob,0) = 0.5 * (ks * gw.k(iGlob,3) + ksx * gw.k(iGlob+1,3));
				if (ks * ksx == 0.0)    {gw.k(iGlob,0) = 0.0;}
				if (par.px > 0)	{
					ivgback = gw.soilID(iGlob-1) * NVG;     ksback = gw.vgTable(ivgback);
	                gw.k(iGlob-1,0) = 0.5 * (ks * gw.k(iGlob,3) + ksback * gw.k(iGlob-1,3));
	                if (ks * ksback == 0.0)  {gw.k(iGlob-1,0) = 0.0;}
				}
			}
			else if (ii == gdom.nx-1)	{
				if (par.px < par.nproc_x-1)	{
					gw.k(iGlob,0) = 0.5 * (ks * gw.k(iGlob,3) + ksx * gw.k(iGlob+1,3));
					if (ks * ksx == 0.0)    {gw.k(iGlob,0) = 0.0;}
				}
			}
			else {
				gw.k(iGlob,0) = 0.5 * (ks * gw.k(iGlob,3) + ksx * gw.k(iGlob+1,3));
				if (ks * ksx == 0.0)    {gw.k(iGlob,0) = 0.0;}
			}
            // Ky
			if (jj == 0)	{
				gw.k(iGlob,1) = 0.5 * (ks * gw.k(iGlob,3) + ksy * gw.k(iGlob+gdom.nxhc,3));
                if (ks * ksy == 0.0)    {gw.k(iGlob,1) = 0.0;}
				if (par.py > 0)	{
					ivgback = gw.soilID(iGlob-gdom.nxhc) * NVG;     ksback = gw.vgTable(ivgback);
	                gw.k(iGlob-gdom.nxhc,1) = 0.5 * (ks * gw.k(iGlob,3) + ksback * gw.k(iGlob-gdom.nxhc,3));
	                if (ks * ksback == 0.0)  {gw.k(iGlob-gdom.nxhc,1) = 0.0;}
				}
			}
			else if (jj == gdom.ny-1)	{
				if (par.py < par.nproc_y-1)	{
					gw.k(iGlob,1) = 0.5 * (ks * gw.k(iGlob,3) + ksy * gw.k(iGlob+gdom.nxhc,3));
	                if (ks * ksy == 0.0)    {gw.k(iGlob,1) = 0.0;}
				}
			}
			else {
				gw.k(iGlob,1) = 0.5 * (ks * gw.k(iGlob,3) + ksy * gw.k(iGlob+gdom.nxhc,3));
                if (ks * ksy == 0.0)    {gw.k(iGlob,1) = 0.0;}
			}
            // Kz
            if (kk < gdom.nz-1) {
                gw.k(iGlob,2) = 0.5 * (ks * gw.k(iGlob,3) + ksz * gw.k(iGlob+gdom.nxhc*gdom.nyhc,3));
                if (ks * ksz == 0.0)    {gw.k(iGlob,2) = 0.0;}
            }
        });
		gdom.timers.re.gwUpdateK += timer.seconds();

        // MPI exchange of K
		timer.reset();
        gmpi.mpi_sendrecv(gw.k, gdom, par);
		gdom.timers.re.gwMPI += timer.seconds();

        // Apply boundary conditions
		timer.reset();
        for (int k = 0; k < gbc.size(); k++) {
            gbc[k].applyKBC(gw, gdom, par);
        }
		gdom.timers.re.gwBC += timer.seconds();

        // Zero K for reduced dimension simulation
        if (gdom.nx == 1)   {
            Kokkos::parallel_for( gdom.nCellMem , KOKKOS_LAMBDA(int iGlob) {gw.k(iGlob,0) = 0.0;});
        }
        if (gdom.ny == 1)   {
            Kokkos::parallel_for( gdom.nCellMem , KOKKOS_LAMBDA(int iGlob) {gw.k(iGlob,1) = 0.0;});
        }
	}

    // /* --------------------------------------------------
    //     End of conductivity block
    // -------------------------------------------------- */


    // /* --------------------------------------------------
    //     Get face flux (original - used for GW solver internally)
    // -------------------------------------------------- */
    // [CODE2] Original face_flux function (identical logic in both branches, using Code2 style)
    inline void face_flux(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc, GwMPI &gmpi, Parallel &par)	{
        // Initialize Q to zero (this also set boundary Q=0 by default)
		timer.reset();
        Kokkos::parallel_for( gdom.nCellMem , KOKKOS_LAMBDA(int iGlob) {
            gw.q(iGlob,0) = 0.0; gw.q(iGlob,1) = 0.0;   gw.q(iGlob,2) = 0.0;
        });
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
            gw.q(iGlob,0) = gw.k(iGlob,0) * gdom.cosx(iGlob) * (gw.h(iGlob+1,1) - gw.h(iGlob,1)) / gdom.dx
                + gw.k(iGlob,0) * gdom.sinx(iGlob);
            gw.q(iGlob,1) = gw.k(iGlob,1) * gdom.cosy(iGlob) * (gw.h(iGlob+gdom.nxhc,1) - gw.h(iGlob,1)) / gdom.dy
                + gw.k(iGlob,1) * gdom.siny(iGlob);
            // gw.q(iGlob,2) = gw.k(iGlob,2) * (gw.h(iGlob+gdom.nxhc*gdom.nyhc,1) - gw.h(iGlob,1)) / gdom.dz(iGlob)
            //     - gw.k(iGlob,2);
			gw.q(iGlob,2) = gw.k(iGlob,2) * (gw.h(iGlob+gdom.nxhc*gdom.nyhc,1) - gw.h(iGlob,1)) / (0.5*(gdom.dz(iGlob)+gdom.dz(iGlob+gdom.nxhc*gdom.nyhc)))
                - gw.k(iGlob,2);
        });
		gdom.timers.re.gwUpdateQ += timer.seconds();

        // MPI exchange of flux
		timer.reset();
        gmpi.mpi_sendrecv(gw.q, gdom, par);
		gdom.timers.re.gwMPI += timer.seconds();

        // Apply boundary conditions
		timer.reset();
        for (int k = 0; k < gbc.size(); k++) {
            gbc[k].applyQBC(gw, gdom, par);
        }
		gdom.timers.re.gwBC += timer.seconds();
	}

    // /* --------------------------------------------------
    //     Get face flux NEW (for solute transport coupling)
    // -------------------------------------------------- */
    // [CODE1] face_flux_new: computes gw.q_new which is used by RTFunctionGW for solute transport
    // This function stores fluxes in a separate array (q_new) to avoid interfering with the
    // solver's internal q array. RTFunctionGW uses q_new for velocity computation in the
    // reactive transport equation.
    inline void face_flux_new(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc, GwMPI &gmpi, Parallel &par)
    {
        // Initialize q_new to zero
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            gw.q_new(iGlob,0) = 0.0; gw.q_new(iGlob,1) = 0.0; gw.q_new(iGlob,2) = 0.0;
        });
        // [CODE1] Compute q_new using same flux formulation as face_flux but storing in q_new
        Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
            gw.q_new(iGlob, 0) = gw.k(iGlob, 0) * gdom.cosx(iGlob) * (gw.h(iGlob + 1, 1) - gw.h(iGlob, 1)) / gdom.dx + gw.k(iGlob, 0) * gdom.sinx(iGlob);
            gw.q_new(iGlob, 1) = gw.k(iGlob, 1) * gdom.cosy(iGlob) * (gw.h(iGlob + gdom.nxhc, 1) - gw.h(iGlob, 1)) / gdom.dy + gw.k(iGlob, 1) * gdom.siny(iGlob);
            gw.q_new(iGlob, 2) = gw.k(iGlob, 2) * (gw.h(iGlob + gdom.nxhc * gdom.nyhc, 1) - gw.h(iGlob, 1)) / (0.5 * (gdom.dz(iGlob) + gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - gw.k(iGlob, 2);
        });
        // MPI exchange of q_new
        gmpi.mpi_sendrecv(gw.q_new, gdom, par);
        // Apply boundary conditions (using same BC as face_flux)
        for (int k = 0; k < gbc.size(); k++)
        {
            gbc[k].applyQBC(gw, gdom, par);
        }
    }

    // /* --------------------------------------------------
    //     End of flux block
    // -------------------------------------------------- */


    // /* --------------------------------------------------
    //     Get matrix coefficients
    // -------------------------------------------------- */
    // [CODE1] CROP_GROWTH_MODEL overload: adds Wofost72 parameter for crop coupling in source/sink
#if CROP_GROWTH_MODEL
    inline void linear_system(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc, std::vector<GwSS> &gss, GwMatrix &A, Parallel &par, Wofost72 &wofost){
#else
    inline void linear_system(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc, std::vector<GwSS> &gss, GwMatrix &A, Parallel &par)	{
#endif
        // Calculate matrix coefficients
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, ivg, iGlob, iGlobSW;
            real wcs, wcr, wcm, n, m, alpha, nume, deno, ch = 0.0, ss = 1e-5;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
            iGlobSW = (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
            ivg = gw.soilID(iGlob) * NVG;
            wcs = gw.vgTable(ivg+2);     wcr = gw.vgTable(ivg+3);
            n = gw.vgTable(ivg+4);       alpha = gw.vgTable(ivg+6);
            m = 1.0 - 1.0 / n;
            if (gw.h(iGlob,1) < gdom.aev)	{
                wcm = wcr + (wcs-wcr)*mypow((1.0 + mypow(myfabs(gdom.aev)*alpha,n)), m);
                nume = alpha*n*m*(wcm - wcr)*mypow(myfabs(alpha*gw.h(iGlob,1)),n-1.0);
                deno = mypow((1.0 + mypow(myfabs(alpha*gw.h(iGlob,1)),n)), m+1);
                ch = nume / deno;
            }
            gw.coef(idom,0) = ch + ss*gw.wc(iGlob,1)/wcs;
            gw.coef(idom,1) = - gdom.dt * gw.k(iGlob,0) * gdom.cosx(iGlob) / mypow(gdom.dx, 2.0);
            gw.coef(idom,2) = - gdom.dt * gw.k(iGlob-1,0) * gdom.cosx(iGlob-1) / mypow(gdom.dx, 2.0);
            gw.coef(idom,3) = - gdom.dt * gw.k(iGlob,1) * gdom.cosy(iGlob) / mypow(gdom.dy, 2.0);
            gw.coef(idom,4) = - gdom.dt * gw.k(iGlob-gdom.nxhc,1) * gdom.cosy(iGlob-gdom.nxhc) / mypow(gdom.dy, 2.0);
            // gw.coef(idom,5) = - gdom.dt * gw.k(iGlob,2) / mypow(gdom.dz(iGlob), 2.0);
            // gw.coef(idom,6) = - gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / mypow(gdom.dz(iGlob), 2.0);
			gw.coef(idom,5) = - gdom.dt * gw.k(iGlob,2) / gdom.dz(iGlob) / (0.5*(gdom.dz(iGlob)+gdom.dz(iGlob+gdom.nxhc*gdom.nyhc)));
            gw.coef(idom,6) = - gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / gdom.dz(iGlob) / (0.5*(gdom.dz(iGlob)+gdom.dz(iGlob-gdom.nxhc*gdom.nyhc)));
            gw.coef(idom,7) = (ch + ss*gw.wc(iGlob,1)/wcs)*gw.h(iGlob,1)
                - gdom.dt*(gw.k(iGlob,2) - gw.k(iGlob-gdom.nxhc*gdom.nyhc,2)) / gdom.dz(iGlob)
                + gdom.dt*(gw.k(iGlob,0) * gdom.sinx(iGlob) - gw.k(iGlob-1,0) * gdom.sinx(iGlob-1))/gdom.dx
                + gdom.dt*(gw.k(iGlob,1) * gdom.siny(iGlob) - gw.k(iGlob-gdom.nxhc,1) * gdom.siny(iGlob-gdom.nxhc))/gdom.dy;
            if (gdom.gw_scheme != 1) {
                gw.coef(idom,7) -= (gw.wc(iGlob,1) - gw.wc(iGlob,0));
            }
            // Apply internal boundary conditions (needed when MPI is used)
            if (ii == 0 && par.px > 0)    {gw.coef(idom,7) -= gw.coef(idom,2) * gw.h(iGlob-1,1);}
            else if (ii == gdom.nx-1 && par.px < par.nproc_x-1)   {gw.coef(idom,7) -= gw.coef(idom,1) * gw.h(iGlob+1,1);}
            if (jj == 0 && par.py > 0)    {gw.coef(idom,7) -= gw.coef(idom,4) * gw.h(iGlob-gdom.nxhc,1);}
            else if (jj == gdom.ny-1 && par.py < par.nproc_y-1)   {gw.coef(idom,7) -= gw.coef(idom,3) * gw.h(iGlob+gdom.nxhc,1);}
            
        });
        // Apply outer boundary conditions
        for (int k = 0; k < gbc.size(); k++) {
            gbc[k].applyMatBC(gw, gdom, par);
        }

        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            gw.coef(idom,0) -= (gw.coef(idom,1)+gw.coef(idom,2)+gw.coef(idom,3)+gw.coef(idom,4)+gw.coef(idom,5)+gw.coef(idom,6));
        });
        
        // [CODE2] Remove dependencies on internal NODATA cells
        // This prevents the solver from coupling NODATA cells with active cells
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
            // no data cells 
            if (gdom.isnodata(iGlob) == 1)  {
            	gw.coef(idom,0) = 1e12; gw.coef(idom,7) = 1e12; gw.coef(idom,5) = 0.0;	gw.coef(idom,6) = 0.0;
            	gw.coef(idom,1) = 0.0; gw.coef(idom,2) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,4) = 0.0;	
            }
            else {
            	if (gdom.isnodata(iGlob+1) == 1)	{gw.coef(idom,1) = 0.0;}
            	if (gdom.isnodata(iGlob-1) == 1)	{gw.coef(idom,2) = 0.0;}
            	if (gdom.isnodata(iGlob+gdom.nxhc) == 1)	{gw.coef(idom,3) = 0.0;}
            	if (gdom.isnodata(iGlob-gdom.nxhc) == 1)	{gw.coef(idom,4) = 0.0;}
            }
        });

        // Apply internal source/sink terms
        // [CODE1] CROP_GROWTH_MODEL: source/sink with Wofost72 for root water uptake
		for (int k = 0; k < gss.size(); k++) {
#if CROP_GROWTH_MODEL
            gss[k].applyMatSS(gw, gdom, wofost);
#else
            gss[k].applyMatSS(gw, gdom);
#endif
        }

        // Insert coefficients into Matrix A
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob, irow = A.ptr(idom);
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
        	if (kk > 0)	{A.ind(irow) = idom - gdom.nx*gdom.ny;	A.val(irow) = gw.coef(idom,6);  irow++;}
        	if (jj > 0)	{
        		A.ind(irow) = idom - gdom.nx;	        A.val(irow) = gw.coef(idom,4);  irow++;
        	}
        	if (ii > 0)	{
        		A.ind(irow) = idom - 1;		        A.val(irow) = gw.coef(idom,2);  irow++;
        	}
        	A.ind(irow) = idom;	A.val(irow) = gw.coef(idom,0);	irow++;
        	if (ii < gdom.nx-1)	{
        		A.ind(irow) = idom + 1;		        A.val(irow) = gw.coef(idom,1);  irow++;
        	}
        	if (jj < gdom.ny-1)	{
        		A.ind(irow) = idom + gdom.nx;	        A.val(irow) = gw.coef(idom,3);  irow++;
        	}
        	if (kk < gdom.nz-1)	{A.ind(irow) = idom + gdom.nx*gdom.ny;	A.val(irow) = gw.coef(idom,5);  irow++;}
        	A.rhs(idom) = gw.coef(idom,7);
        });

	}
    // /* --------------------------------------------------
    //     End matrix coefficients
    // -------------------------------------------------- */

    // /* --------------------------------------------------
    //     Get water content
    // -------------------------------------------------- */
    // [CODE1] CROP_GROWTH_MODEL overload: adds Wofost72 parameter for crop coupling in source/sink
#if CROP_GROWTH_MODEL
    inline void update_wc(GwState &gw, GwDomain &gdom, std::vector<GwSS> &gss, Wofost72 &wofost)
#else
    inline void update_wc(GwState &gw, GwDomain &gdom, std::vector<GwSS> &gss)
#endif
    {
        // Update wc with explicit scheme
        if (gdom.gw_scheme == 1)    {
            // [MERGED] Combine Code1's CROP_GROWTH_MODEL evaporation with Code2's NODATA checks
            Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob, ivg, iGlobSW;
                real coef, qqx, qqy, qqz, wcs, ss = 1e-5;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
				iGlobSW = (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                ivg = gw.soilID(iGlob) * NVG;
                wcs = gw.vgTable(ivg+2);
                coef = 1.0 + ss*(gw.h(iGlob,1) - gw.h(iGlob,0)) / wcs;
                qqx = (gw.q(iGlob,0) - gw.q(iGlob-1,0)) / gdom.dx;
                qqy = (gw.q(iGlob,1) - gw.q(iGlob-gdom.nxhc,1)) / gdom.dy;
                qqz = (gw.q(iGlob,2) - gw.q(iGlob-gdom.nxhc*gdom.nyhc,2)) / gdom.dz(iGlob);
                // [CODE2] Only update for non-NODATA cells
                if (gdom.isnodata(iGlob) == 0)  {
                    gw.wc(iGlob,1) = (gw.wc(iGlob,0) + gdom.dt * (qqx + qqy + qqz)) / coef;
    				// [CODE1] Evaporation: surface water evaporation coupling
    				// Applied only on top layer (kk==0) when no ponded water above
    				if (kk == 0 && gw.h(iGlob-gdom.nxhc*gdom.nyhc,1) <= 0 && gdom.isEvap == 1)	{
    					gw.wc(iGlob,1) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);
    				}
                    gw.wc(iGlob,2) = 0.0;
                }
                
            });
			// source/sink terms
			// [CODE1] CROP_GROWTH_MODEL: source/sink with Wofost72 for root water uptake / transpiration
			for (int k = 0; k < gss.size(); k++) {
#if CROP_GROWTH_MODEL
                gss[k].applyWCSS(gw, gdom, wofost);
#else
			    gss[k].applyWCSS(gw, gdom);
#endif
			}
            // Choose h or wc at the interface
            // [MERGED] Combine Code1's logic with Code2's NODATA checks
            Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob, iGlobSW, ivg, flag;
        		real wcs, wcr, wcm, n, m, alpha, sbar;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                iGlobSW = (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                ivg = gw.soilID(iGlob) * NVG;
                wcs = gw.vgTable(ivg+2);     wcr = gw.vgTable(ivg+3);
                n = gw.vgTable(ivg+4);       alpha = gw.vgTable(ivg+6);
                m = 1.0 - 1.0 / n;
                wcm = wcr + (wcs-wcr)*mypow((1.0 + mypow(myfabs(gdom.aev)*alpha,n)), m);
				// make choice
                if (gw.wc(iGlob,1) - wcs > TOL8NEG)	{
                    gw.wc(iGlob,2) = gw.wc(iGlob,1)-wcs;    gw.wc(iGlob,1) = wcs;
                }
                else    {
                    flag = 0;
                    if (gw.wc(iGlob+1,1) - wcs > TOL8NEG & gw.k(iGlob,0) > 0.0)    {flag = 1;}
                    else if (gw.wc(iGlob-1,1) - wcs > TOL8NEG & gw.k(iGlob-1,0) > 0.0) {flag = 1;}
                    else if (gw.wc(iGlob+gdom.nxhc,1) - wcs > TOL8NEG & gw.k(iGlob,1) > 0.0)  {flag = 1;}
                    else if (gw.wc(iGlob-gdom.nxhc,1) - wcs > TOL8NEG & gw.k(iGlob-gdom.nxhc,1) > 0.0)  {flag = 1;}
                    else if (gw.wc(iGlob+gdom.nxhc*gdom.nyhc,1) - wcs > TOL8NEG & gw.k(iGlob,2) > 0.0)  {flag = 1;}
                    else if (gw.wc(iGlob-gdom.nxhc*gdom.nyhc,1) - wcs > TOL8NEG & gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) > 0.0)  {flag = 1;}

                    // Use head form for the top layer
                    // Not sure if this works for impermeable top boundary ?
                    if (ii == 0 && gw.h(iGlob-1,0) >= 0.0)  {flag = 1;}
                    if (ii == gdom.nx-1 && gw.h(iGlob+1,0) >= 0.0)  {flag = 1;}

                    if (kk == 0)    {
						// HERE: If a fixed pressure h=0 is enforced on the top boundary, there are 2 situations:
						// 1 : h=0 means a thin-layer of inundated water
						// 2 : h=0 means atmosphere pressure, no water
						// In Case2, flag should be 0. In Case1, flag should be 1.
						// To avoid confusing, in Case1, we recommend to set h=0.0001 to represent the thin-layer of water
                        if (gw.wc(iGlob,1) < wcs-TOL8NEG && gw.h(iGlob-gdom.nxhc*gdom.nyhc,1) == 0.0)   {flag = 0;}
                        else {flag = 1;}
                    }
                    // [CODE2] Only apply h/wc switch for non-NODATA cells
                    if (gdom.isnodata(iGlob) == 0)  {
                        if (flag == 1)  {
                            real tmp = gw.wc(iGlob,1);
                            sbar = mypow(1.0 + mypow(myfabs(alpha*gw.h(iGlob,1)), n), -m);
                            if (gw.h(iGlob,1) > gdom.aev)   {gw.wc(iGlob,1) = wcs;}
                            else {gw.wc(iGlob,1) = sbar * (wcm - wcr) + wcr;}
                            gw.wc(iGlob,2) = tmp - gw.wc(iGlob,1);
                        }
                        else    {
                            if (gw.wc(iGlob,1) < wcs)   {
                                if (gw.wc(iGlob,1) < wcr)   {gw.wc(iGlob,1) = wcr + 1e-5;}
                                gw.h(iGlob,1) = -(1.0/alpha) * (mypow(mypow((wcm-wcr)/(gw.wc(iGlob,1)-wcr),(1/m)) - 1.0, 1/n));
                            }
                            else {gw.h(iGlob,1) = 0.0;}
                        }
                    }
                }
            });
        }
        else {
            // [MERGED] Implicit scheme: Code2's NODATA checks + Code1's logic
            Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob, ivg, flag;
        		real wcs, wcr, wcm, n, m, alpha, sbar;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                ivg = gw.soilID(iGlob) * NVG;
                wcs = gw.vgTable(ivg+2);     wcr = gw.vgTable(ivg+3);
                n = gw.vgTable(ivg+4);       alpha = gw.vgTable(ivg+6);
                m = 1.0 - 1.0 / n;
                wcm = wcr + (wcs-wcr)*mypow((1.0 + mypow(myfabs(gdom.aev)*alpha,n)), m);
                sbar = mypow(1.0 + mypow(myfabs(alpha*gw.h(iGlob,1)), n), -m);
                // [CODE2] Only update for non-NODATA cells
                if (gdom.isnodata(iGlob) == 0)  {
                    if (gw.h(iGlob,1) > gdom.aev)   {gw.wc(iGlob,1) = wcs;}
                    else {gw.wc(iGlob,1) = sbar * (wcm - wcr) + wcr;}
                    if (gw.wc(iGlob,1) > wcs)	{gw.wc(iGlob,2) += (gw.wc(iGlob,1)-wcs); gw.wc(iGlob,1) = wcs;}
            		else if (gw.wc(iGlob,1) < wcr+1e-5)	{gw.wc(iGlob,1) = wcr+1e-5;}
                }
            });
        }

    }
    // /* --------------------------------------------------
    //     End water content block
    // -------------------------------------------------- */



    // /* --------------------------------------------------
    //     Update dt based on either water content or iteration
    // -------------------------------------------------- */
    // [MERGED] Identical in both branches, using Code2's gdom.hc style
    inline void dt_waco(GwState &gw, GwDomain &gdom)	{
    	real dwc_max, dt_old;
    	dt_old = gdom.dt;
        Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA (int idx, real &tmp) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idx, kk, jj, ii);
            // gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
            real dwc = myfabs(gw.wc(iGlob,1) - gw.wc(iGlob,0));
			tmp = (dwc > tmp) ? dwc : tmp;
		} , Kokkos::Max<real>(dwc_max) );
    	if (dwc_max > 0.02)	{gdom.dt = gdom.dt * 0.9;}
    	else if (dwc_max >= 0.0 & dwc_max < 0.01)	{gdom.dt = gdom.dt * 1.1;}
    	if (gdom.dt > gdom.dt_max)	{gdom.dt = gdom.dt_max;}
    	else if (gdom.dt < gdom.dt_init)	{gdom.dt = gdom.dt_init;}
    }

    inline void dt_iter(GwState &gw, GwDomain &gdom, int iter)	{
    	real dt_old;
    	dt_old = gdom.dt;
        if (iter < 7)   {gdom.dt = gdom.dt * 1.1;}
        else if (iter > 11)  {gdom.dt = gdom.dt * 0.9;}
        if (gdom.dt > gdom.dt_max)	{gdom.dt = gdom.dt_max;}
    	else if (gdom.dt < gdom.dt_init)	{gdom.dt = gdom.dt_init;}
    }

    // /* --------------------------------------------------
    //     End dt update
    // -------------------------------------------------- */



    // /* --------------------------------------------------
    //     Get the maximum difference between iterations
    // -------------------------------------------------- */
    inline real get_eps(GwState &gw, GwDomain &gdom)	{
    	real eps;
        Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA (int idx, real &tmp) {
            real dwc = myfabs(gw.h(idx,1) - gw.h(idx,0));
			tmp = (dwc > tmp) ? dwc : tmp;
		} , Kokkos::Max<real>(eps) );
        return eps;
    }

    // /* --------------------------------------------------
    //     End eps computation
    // -------------------------------------------------- */
};

#endif