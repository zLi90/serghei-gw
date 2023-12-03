#ifndef _GW_FUNCTION_H_
#define _GW_FUNCTION_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMatrix.h"
#include "GwMPI.h"
#include "GwSolver.h"
#include "GwState.h"
#include "State.h"
#include <set>
#include <math.h>


class GwFunction   {

public:
    /* --------------------------------------------------
        Top-level PCA solver
    -------------------------------------------------- */
    inline void pca_solve(GwState &gw, State &state, GwDomain &gdom, GwBC &gbc,
            GwMatrix &A, GwSolver &gsolver, SourceSinkData &ss, GwMPI &gmpi, Parallel &par)  {
        int iter, ierr=1;
        real dt_tmp;
        enforce_swe_bc(gw, state, gdom, gbc, ss);
        face_conductivity(gw, gdom, gbc, gmpi, par);
        linear_system(gw, gdom, gbc, A, par);
        iter = gsolver.cg(A, gdom);
        Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            gw.h(iGlob,1) = A.x(idom);
        });
        gmpi.mpi_sendrecv(gw.h, gdom, par);
        enforce_lateral_bc(gw, gdom, gbc, par);

        face_conductivity(gw, gdom, gbc, gmpi, par);

        face_flux(gw, state, gdom, gbc, gmpi, par);

        update_wc(gw, state, gdom);
        gmpi.mpi_sendrecv(gw.h, gdom, par);
        gmpi.mpi_sendrecv(gw.wc, gdom, par);
        enforce_lateral_bc(gw, gdom, gbc, par);

        dt_waco(gw, gdom);
        dt_tmp = gdom.dt;
        ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            gw.h(iGlob,0) = gw.h(iGlob,1);  gw.wc(iGlob,0) = gw.wc(iGlob,1);
        });

        gw.Vtot = integrate(gw, gdom);
        gw.Vexch = get_Vexchange(state, gdom);

    }

    /* --------------------------------------------------
        Top-level Picard solver
    -------------------------------------------------- */
    inline void picard_solve(GwState &gw, State &state, GwDomain &gdom, GwBC &gbc,
            GwMatrix &A, GwSolver &gsolver, SourceSinkData &ss, GwMPI &gmpi, Parallel &par)  {
        int iter, iter_cg, iter_max = 100, ierr=1;
        real eps_diff = 1.0, eps_tmp, eps_old = 1.0, eps = 1.0, eps_min = 5e-6, dt_tmp;
        enforce_swe_bc(gw, state, gdom, gbc, ss);
        face_conductivity(gw, gdom, gbc, gmpi, par);

        iter = 0;
        while (iter < iter_max && eps_diff/eps_old > eps_min && eps > eps_min) {
            linear_system(gw, gdom, gbc, A, par);
            iter_cg = gsolver.cg(A, gdom);
            Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                gw.h(iGlob,0) = gw.h(iGlob,1);
                gw.h(iGlob,1) = A.x(idom);
            });
            gmpi.mpi_sendrecv(gw.h, gdom, par);
            enforce_lateral_bc(gw, gdom, gbc, par);

            face_conductivity(gw, gdom, gbc, gmpi, par);

            face_flux(gw, state, gdom, gbc, gmpi, par);

            eps_old = eps;
            eps = get_eps(gw, gdom);
            eps_tmp = fabs(eps_old - eps);
            ierr = MPI_Allreduce(&eps_tmp, &eps_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            update_wc(gw, state, gdom);
            gmpi.mpi_sendrecv(gw.wc, gdom, par);

            iter += 1;
        }
        // printf("    > Picard loop converges in %d iterations with eps = %f, %f\n",iter,eps,eps_diff);

        dt_iter(gw, gdom, iter);
        dt_tmp = gdom.dt;
        ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            gw.h(iGlob,0) = gw.h(iGlob,1);  gw.wc(iGlob,0) = gw.wc(iGlob,1);
        });

        gw.Vtot = integrate(gw, gdom);
        gw.Vexch = get_Vexchange(state, gdom);
    }

    // /* --------------------------------------------------
    //     Get pressure BC from SWE module
    // -------------------------------------------------- */
    inline void enforce_swe_bc(GwState &gw, State &state, GwDomain &gdom, GwBC &gbc, SourceSinkData &ss)	{
        Kokkos::parallel_for("enforce_swe_BC", gdom.nhalo , KOKKOS_LAMBDA(int idx) {
            int ii, jj, kk, iGlobSW, ivg, iGlob = gdom.hpair(idx,1);
            real ks;
            gdom.unpackIndicesGw(iGlob, gdom.nzhc, gdom.nyhc, gdom.nxhc, kk, jj, ii);
            iGlobSW = jj*gdom.nxhc + ii;
            ivg = gw.soilID(iGlob) * NVG;
            ks = gw.vgTable(ivg);
            if (gdom.hpair(idx,2) == -3)  {
                if (gbc.bctypeZM == SUB_BC_H_SWE)    {
                    gw.h(gdom.hpair(idx,0),1) = state.h(iGlobSW);
                    if (gdom.isRain == 1)    {
                        gdom.qrain(iGlobSW) = ss.rainRate(iGlobSW);
                        if (state.h(iGlobSW) <= 0.0)  {gw.wc(gdom.hpair(idx,0),1) = gw.wc(gdom.hpair(idx,0),1)/2.0;}
                    }
                    // assign BC type for the top layer
                    // Surface Ponding
                    if (gw.h(gdom.hpair(idx,0),1) > 0.0)    {
                        real q_infilt = 2.0 * ks * (gw.h(gdom.hpair(idx,1),1) - gw.h(gdom.hpair(idx,0),1)) / gdom.dz(iGlob) - ks;
                        if (-q_infilt * gdom.dt <= gw.h(gdom.hpair(idx,0),1))    {
                            gbc.topBC(iGlobSW) = 1;
                        }
                        else    {gbc.topBC(iGlobSW) = 2;}
                    }
                    else {
                        // Exfiltration
                        if (gw.h(gdom.hpair(idx,1),1) > gw.h(gdom.hpair(idx,0),1) + 0.5*gdom.dz(iGlob)) {
                            gbc.topBC(iGlobSW) = 1;
                        }
                        // No flow
                        else {gbc.topBC(iGlobSW) = 0;}
                    }
                }
            }
        });
	}

    inline void enforce_lateral_bc(GwState &gw, GwDomain &gdom, GwBC &gbc, Parallel &par)    {
        if (gbc.bctypeXM != SUB_BC_NOFLOW || gbc.bctypeXP != SUB_BC_NOFLOW ||
            gbc.bctypeYM != SUB_BC_NOFLOW || gbc.bctypeYP != SUB_BC_NOFLOW ) {
            Kokkos::parallel_for("enforce_side_BC",  gdom.nCell , KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, ii2, ii3, iGlob, ivg;
                real wcs, wcr, alpha, n;
                gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                ii2 = kk*gdom.ny_glob + (par.j_beg+jj);
                ii3 = kk*gdom.nx_glob + (par.i_beg+ii);
                ivg = gw.soilID(iGlob) * gw.nVGparam;
                wcs = gw.vgTable(ivg+2);    wcr = gw.vgTable(ivg+3);
                n = gw.vgTable(ivg+4);  alpha = gw.vgTable(ivg+6);
                if (ii == 0 && par.px == 0 && gbc.bctypeXM != SUB_BC_NOFLOW)    {
                    gw.h(iGlob-1,1) = gw.hbcX(ii2,0);
                    gw.wc(iGlob-1,1) = h2wc(gw.h(iGlob-1,1), alpha, n, wcs, wcr);
                }
                else if (ii == gdom.nx-1 && par.px == par.nproc_x-1 && gbc.bctypeXP != SUB_BC_NOFLOW)    {
                    gw.h(iGlob+1,1) = gw.hbcX(ii2,1);
                    gw.wc(iGlob+1,1) = h2wc(gw.h(iGlob+1,1), alpha, n, wcs, wcr);
                }
                if (jj == 0 && par.py == 0 && gbc.bctypeYM != SUB_BC_NOFLOW)    {
                    gw.h(iGlob-gdom.nxhc,1) = gw.hbcY(ii3,0);
                    gw.wc(iGlob-gdom.nxhc,1) = h2wc(gw.h(iGlob-gdom.nxhc,1), alpha, n, wcs, wcr);
                }
                else if (jj == gdom.ny-1 && par.py == par.nproc_y-1 && gbc.bctypeYP != SUB_BC_NOFLOW)    {
                    gw.h(iGlob+gdom.nxhc,1) = gw.hbcY(ii3,1);
                    gw.wc(iGlob+gdom.nxhc,1) = h2wc(gw.h(iGlob+gdom.nxhc,1), alpha, n, wcs, wcr);
                }
            });
        }
    }

    //
    // /* --------------------------------------------------
    //     End of pressure BC block
    // -------------------------------------------------- */
    
    
    // /* --------------------------------------------------
    //     Get rainfall BC from SWE module
    // -------------------------------------------------- */
    inline void enforce_rainfall_bc(State &state, GwDomain &gdom, Domain &dom, SourceSinkData &ss)	{
    	Kokkos::parallel_for( dom.nCellMem , KOKKOS_LAMBDA(int iGlob) {gdom.qrain(iGlob) = ss.rainRate(iGlob); });
        Kokkos::parallel_for( dom.nCell , KOKKOS_LAMBDA (int idom) {
            int iGlob = dom.getIndex(idom);
            state.h(iGlob) += ss.rainRate(iGlob)*dom.dt;
        });
    }
    
    // /* --------------------------------------------------
    //     End of rainfall BC block
    // -------------------------------------------------- */


    /* --------------------------------------------------
        Get face conductivity
    -------------------------------------------------- */
    inline void face_conductivity(GwState &gw, GwDomain &gdom, GwBC &gbc, GwMPI &gmpi, Parallel &par)	{
        // Get relatively permeability at cell centers
        
        Kokkos::parallel_for("relative_permeability",  gdom.nCellMem , KOKKOS_LAMBDA(int iGlob) {
            int ii, jj, kk, ivg;
            real s, alpha, n, m, wcm, wcr, wcs, nume, deno;
            gdom.unpackIndicesGw(iGlob, gdom.nzhc, gdom.nyhc, gdom.nxhc, kk, jj, ii);
            ivg = gw.soilID(iGlob) * NVG;
            wcs = gw.vgTable(ivg+2);
            wcr = gw.vgTable(ivg+3);
            n = gw.vgTable(ivg+4);
            alpha = gw.vgTable(ivg+6);
            m = 1.0 - 1.0 / n;
            wcm = wcr + (wcs-wcr)*pow((1.0 + pow(fabs(gdom.aev)*alpha,n)), m);
            s = pow(1.0 + pow(fabs(alpha*gw.h(iGlob,1)), n), -m);
            nume = 1.0-pow(1.0-pow(s*(wcs-wcr)/(wcm-wcr),1.0/m),m);
            deno = 1.0-pow(1.0-pow((wcs-wcr)/(wcm-wcr),1.0/m),m);
            if (deno == 0.0)    {gw.k(iGlob,3) = 1.0;}
            else {gw.k(iGlob,3) = pow(s,0.5) * pow(nume/deno, 2.0);}
        	gw.k(iGlob,3) = pow(s,0.5) * pow(1-pow(1-pow(s,1.0/m),m), 2.0);
            if (gw.k(iGlob,3) > 1.0 | gw.h(iGlob,1) >= gdom.aev)	{gw.k(iGlob,3) = 1.0;}
            // set no data cells impermeable
            if (gdom.isnodata(iGlob) == 1)  {gw.k(iGlob,3) = 0.0;}
        });
        // Get K on interior cell faces
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob, ivg;
            real ks;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            ivg = gw.soilID(iGlob) * NVG;
            ks = gw.vgTable(ivg);
            gw.k(iGlob,0) = 0.5 * ks * (gw.k(iGlob,3) + gw.k(iGlob+1,3));
            gw.k(iGlob,1) = 0.5 * ks * (gw.k(iGlob,3) + gw.k(iGlob+gdom.nxhc,3));
            gw.k(iGlob,2) = 0.5 * ks * (gw.k(iGlob,3) + gw.k(iGlob+gdom.nxhc*gdom.nyhc,3));
        });
        // Set K=0 for impervious layers
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob, ivg, ivgx, ivgy, ivgz;
            real ks, ksx, ksy, ksz;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            ivg = gw.soilID(iGlob) * NVG;
            ivgx = gw.soilID(iGlob+1) * NVG;
            ivgy = gw.soilID(iGlob+gdom.nxhc) * NVG;
            ivgz = gw.soilID(iGlob+gdom.nxhc*gdom.nyhc) * NVG;
            ks = gw.vgTable(ivg);
            ksx = gw.vgTable(ivgx);
            ksy = gw.vgTable(ivgy);
            ksz = gw.vgTable(ivgz);
            if (ks * ksx == 0.0) {gw.k(iGlob,0) = 0.0;}
            if (ks * ksy == 0.0) {gw.k(iGlob,1) = 0.0;}
            if (ks * ksz == 0.0) {gw.k(iGlob,2) = 0.0;}
        });
        // MPI exchange of K
        gmpi.mpi_sendrecv(gw.k, gdom, par);
        // Apply boundary conditions
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob, iGlobSW, ivg;
            real ks, ksx, ksy;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            iGlobSW = (hc+jj)*gdom.nxhc + ii + hc;
            ivg = gw.soilID(iGlob) * NVG;
            ks = gw.vgTable(ivg);
            // Kx boundary
            if (ii == 0 && par.px == 0) {
                if (gbc.bctypeXM == SUB_BC_NOFLOW)  {gw.k(iGlob-1,0) = 0.0;}
                else {
                    gw.k(iGlob-1,0) = 0.5 * ks * (gw.k(iGlob,3) + gw.k(iGlob-1,3));
                }
            }
            else if (ii == gdom.nx-1 && par.px == par.nproc_x-1)    {
                if (gbc.bctypeXP == SUB_BC_NOFLOW)  {gw.k(iGlob,0) = 0.0;}
                else {
                    gw.k(iGlob,0) = 0.5 * ks * (gw.k(iGlob,3) + gw.k(iGlob+1,3));
                }
            }
            if (gdom.nx == 1)  {gw.k(iGlob-1,0) = 0.0; gw.k(iGlob,0) = 0.0;}
            // Ky boundary
            if (jj == 0 && par.py == 0) {
                if (gbc.bctypeYM == SUB_BC_NOFLOW)  {gw.k(iGlob-gdom.nxhc,1) = 0.0;}
                else {
                    gw.k(iGlob-gdom.nxhc,1) = 0.5 * ks * (gw.k(iGlob,3) + gw.k(iGlob-gdom.nxhc,3));
                }
            }
            else if (jj == gdom.ny-1 && par.py == par.nproc_y-1)    {
                if (gbc.bctypeYP == SUB_BC_NOFLOW)  {gw.k(iGlob,1) = 0.0;}
                else {
                    gw.k(iGlob,1) = 0.5 * ks * (gw.k(iGlob,3) + gw.k(iGlob+gdom.nxhc,3));
                }
            }
            if (gdom.ny == 1)  {gw.k(iGlob-gdom.nxhc,1) = 0.0; gw.k(iGlob,1) = 0.0;}
            // Kz boundary
            if (kk == 0)    {
                if (gbc.bctypeZM == SUB_BC_NOFLOW)   {gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) = 0.0;}
                else if (gbc.bctypeZM == SUB_BC_H_CONST) {gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) = ks;}
                else if (gbc.bctypeZM == SUB_BC_H_SWE)  {
                    if (gbc.topBC(iGlobSW) == 0)    {
                        gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) = 0.0;
                    }
                    else {
                        gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) = ks;
                    }
                }
                else {
                    gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) = 0.5 * ks * (gw.k(iGlob,3) + gw.k(iGlob-gdom.nxhc*gdom.nyhc,3));
                }
            }
            else if (kk == gdom.nz-1)   {
                // For now only support impervious bottom
                gw.k(iGlob,2) = 0.0;
                gw.k(iGlob+gdom.nxhc*gdom.nyhc,2) = 0.0;
            }
        });
	}
    // /* --------------------------------------------------
    //     End of conductivity block
    // -------------------------------------------------- */





    // /* --------------------------------------------------
    //     Get face flux
    // -------------------------------------------------- */
    inline void face_flux(GwState &gw, State &state, GwDomain &gdom, GwBC &gbc, GwMPI &gmpi, Parallel &par)	{
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            gw.q(iGlob,0) = gw.k(iGlob,0) * gdom.cosx(iGlob) * (gw.h(iGlob+1,1) - gw.h(iGlob,1)) / gdom.dx
                + gw.k(iGlob,0) * gdom.sinx(iGlob);
            gw.q(iGlob,1) = gw.k(iGlob,1) * gdom.cosy(iGlob) * (gw.h(iGlob+gdom.nxhc,1) - gw.h(iGlob,1)) / gdom.dy
                + gw.k(iGlob,1) * gdom.siny(iGlob);
            gw.q(iGlob,2) = gw.k(iGlob,2) * (gw.h(iGlob+gdom.nxhc*gdom.nyhc,1) - gw.h(iGlob,1)) / gdom.dz(iGlob)
                - gw.k(iGlob,2);
        });
        // MPI exchange of flux
        gmpi.mpi_sendrecv(gw.q, gdom, par);
        // Apply boundary conditions
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, ivg, iGlob, iGlobSW;
            real q_max, wcs;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            ivg = gw.soilID(iGlob) * NVG;
            wcs = gw.vgTable(ivg+2);
            // Get qx
            if (ii == 0 && par.px == 0) {
                if (gbc.bctypeXM == SUB_BC_NOFLOW)   {gw.q(iGlob-1,0) = 0.0;}
                else if (gbc.bctypeXM == SUB_BC_Q_CONST) {gw.q(iGlob-1,0) = gbc.qbcXM;}
                else    {
                    gw.q(iGlob-1,0) = 2.0 * gw.k(iGlob-1,0) * (gw.h(iGlob,1) - gw.h(iGlob-1,1)) / gdom.dx;
                }
            }
            else if (ii == gdom.nx-1 && par.px == par.nproc_x-1)    {
                if (gbc.bctypeXP == SUB_BC_NOFLOW)   {gw.q(iGlob,0) = 0.0;}
                else if (gbc.bctypeXP == SUB_BC_Q_CONST) {gw.q(iGlob,0) = gbc.qbcXP;}
                else    {
                    gw.q(iGlob,0) = 2.0 * gw.k(iGlob,0) * (gw.h(iGlob+1,1) - gw.h(iGlob,1)) / gdom.dx;
                }
            }
            // Get qy
            if (jj == 0 && par.py == 0) {
                if (gbc.bctypeYM == SUB_BC_NOFLOW)   {gw.q(iGlob-gdom.nxhc,1) = 0.0;}
                else if (gbc.bctypeYM == SUB_BC_Q_CONST) {gw.q(iGlob-gdom.nxhc,1) = gbc.qbcYM;}
                else    {
                    gw.q(iGlob-gdom.nxhc,1) = 2.0 * gw.k(iGlob-gdom.nxhc,1) * (gw.h(iGlob,1) - gw.h(iGlob-gdom.nxhc,1)) / gdom.dy;
                }
            }
            else if (jj == gdom.ny-1 && par.py == par.nproc_y-1)    {
                if (gbc.bctypeYP == SUB_BC_NOFLOW)   {gw.q(iGlob,1) = 0.0;}
                else if (gbc.bctypeYP == SUB_BC_Q_CONST) {gw.q(iGlob,1) = gbc.qbcYP;}
                else    {
                    gw.q(iGlob,1) = 2.0 * gw.k(iGlob,1) * (gw.h(iGlob+gdom.nxhc,1) - gw.h(iGlob,1)) / gdom.dy;
                }
            }
            // Get qz
            iGlobSW = packIndicesUniformGrid(gdom.nyhc, gdom.nxhc, jj+hc, ii+hc);
            if (kk == 0)    {
                if (gbc.bctypeZM == SUB_BC_NOFLOW)   {gw.q(iGlob-gdom.nxhc*gdom.nyhc,2) = 0.0;}
                else if (gbc.bctypeZM == SUB_BC_Q_CONST) {
                    if (gdom.isRain == 1) {gw.q(iGlob-gdom.nxhc*gdom.nyhc,2) = -gdom.qrain(iGlobSW);}
                    else {gw.q(iGlob-gdom.nxhc*gdom.nyhc,2) = gbc.qbcZM;}
                }
                else    {
                    gw.q(iGlob-gdom.nxhc*gdom.nyhc,2) = 2.0 * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) *
                        (gw.h(iGlob,1) - gw.h(iGlob-gdom.nxhc*gdom.nyhc,1)) / gdom.dz(iGlob) - gw.k(iGlob-gdom.nxhc*gdom.nyhc,2);
                    // Surface-subsurface exchange flux
                    if (gbc.bctypeZM == SUB_BC_H_SWE)    {
                        if (gbc.topBC(iGlobSW) == 0)    {
                            gw.q(iGlob-gdom.nxhc*gdom.nyhc,2) = 0.0;
                        }
                        else if (gbc.topBC(iGlobSW) == 2)   {
                            gw.q(iGlob-gdom.nxhc*gdom.nyhc,2) = -gw.h(iGlob-gdom.nxhc*gdom.nyhc,1) / gdom.dt;
                        }
                        // Get exchange flux
                        state.qss(iGlobSW) = gw.q(iGlob-gdom.nxhc*gdom.nyhc,2);
                        // Remove exchange flux when subsurface is fully saturated
                        if (state.qss(iGlobSW) < 0.0 & gw.wc(iGlob,1) >= wcs)   {state.qss(iGlobSW) = 0.0;}
                    }
                }
            }
            else if (kk == gdom.nz-1)   {
                if (gbc.bctypeZP == SUB_BC_NOFLOW)   {gw.q(iGlob,2) = 0.0;}
                else if (gbc.bctypeZP == SUB_BC_Q_CONST) {gw.q(iGlob,2) = gbc.qbcZP;}
                else    {gw.q(iGlob,2) = 2.0 * gw.q(iGlob,2);}
            }
        });
	}
    // /* --------------------------------------------------
    //     End of flux block
    // -------------------------------------------------- */




    // /* --------------------------------------------------
    //     Get matrix coefficients
    // -------------------------------------------------- */
    inline void linear_system(GwState &gw, GwDomain &gdom, GwBC &gbc, GwMatrix &A, Parallel &par)	{
        // Calculate matrix coefficients
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, ivg, iGlob, iGlobSW;
            real wcs, wcr, wcm, n, m, alpha, nume, deno, ch = 0.0, ss = 1e-5;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            iGlobSW = packIndicesUniformGrid(gdom.nyhc, gdom.nxhc, jj+hc, ii+hc);
            ivg = gw.soilID(iGlob) * NVG;
            wcs = gw.vgTable(ivg+2);     wcr = gw.vgTable(ivg+3);
            n = gw.vgTable(ivg+4);       alpha = gw.vgTable(ivg+6);
            m = 1.0 - 1.0 / n;
            if (gw.h(iGlob,1) < gdom.aev)	{
                wcm = wcr + (wcs-wcr)*pow((1.0 + pow(fabs(gdom.aev)*alpha,n)), m);
                nume = alpha*n*m*(wcm - wcr)*pow(fabs(alpha*gw.h(iGlob,1)),n-1.0);
                deno = pow((1.0 + pow(fabs(alpha*gw.h(iGlob,1)),n)), m+1);
                ch = nume / deno;
            }
            gw.coef(idom,0) = ch + ss*gw.wc(iGlob,1)/wcs;
            gw.coef(idom,1) = - gdom.dt * gw.k(iGlob,0) * gdom.cosx(iGlob) / pow(gdom.dx, 2.0);
            gw.coef(idom,2) = - gdom.dt * gw.k(iGlob-1,0) * gdom.cosx(iGlob-1) / pow(gdom.dx, 2.0);
            gw.coef(idom,3) = - gdom.dt * gw.k(iGlob,1) * gdom.cosy(iGlob) / pow(gdom.dy, 2.0);
            gw.coef(idom,4) = - gdom.dt * gw.k(iGlob-gdom.nxhc,1) * gdom.cosy(iGlob-gdom.nxhc) / pow(gdom.dy, 2.0);
            gw.coef(idom,5) = - gdom.dt * gw.k(iGlob,2) / pow(gdom.dz(iGlob), 2.0);
            gw.coef(idom,6) = - gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / pow(gdom.dz(iGlob), 2.0);
            gw.coef(idom,7) = (ch + ss*gw.wc(iGlob,1)/wcs)*gw.h(iGlob,1)
                - gdom.dt*(gw.k(iGlob,2) - gw.k(iGlob-gdom.nxhc*gdom.nyhc,2)) / gdom.dz(iGlob)
                + gdom.dt*(gw.k(iGlob,0) * gdom.sinx(iGlob) - gw.k(iGlob-1,0) * gdom.sinx(iGlob-1))/gdom.dx
                + gdom.dt*(gw.k(iGlob,1) * gdom.siny(iGlob) - gw.k(iGlob-gdom.nxhc,1) * gdom.siny(iGlob-gdom.nxhc))/gdom.dy;
            if (gdom.gw_scheme != 1) {
                gw.coef(idom,7) -= (gw.wc(iGlob,1) - gw.wc(iGlob,0));
            }
            // Adjust coefficients on boundaries
            // As of 202302, lateral bc must be Dirichlet type
            if (ii == 0)    {
                if (par.px == 0)    {gw.coef(idom,2) = gw.coef(idom,2) * 2.0;}
                gw.coef(idom,7) -= gw.coef(idom,2) * gw.h(iGlob-1,1);
            }
            else if (ii == gdom.nx-1)   {
                if (par.px == par.nproc_x-1)    {gw.coef(idom,1) = gw.coef(idom,1) * 2.0;}
                gw.coef(idom,7) -= gw.coef(idom,1) * gw.h(iGlob+1,1);
            }
            if (jj == 0)    {
                if (par.py == 0)    {gw.coef(idom,4) = gw.coef(idom,4) * 2.0;}
                gw.coef(idom,7) -= gw.coef(idom,4) * gw.h(iGlob-gdom.nxhc,1);
            }
            else if (jj == gdom.ny-1)   {
                if (par.py == par.nproc_y-1)    {gw.coef(idom,3) = gw.coef(idom,3) * 2.0;}
                gw.coef(idom,7) -= gw.coef(idom,3) * gw.h(iGlob+gdom.nxhc,1);
            }
            if (kk == 0)    {
                if (gbc.bctypeZM == SUB_BC_H_SWE)    {
                    if (gbc.topBC(iGlobSW) == 2)    {
                        real q_infilt = gw.h(iGlob-gdom.nxhc*gdom.nyhc,1) / gdom.dt;
                        gw.coef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / gdom.dz(iGlob);
                        gw.coef(idom,7) += gdom.dt * q_infilt / gdom.dz(iGlob);
                        gw.coef(idom,6) = 0.0;
                    }
                    else {
                        gw.coef(idom,6) = gw.coef(idom,6) * 2.0;
                        gw.coef(idom,7) += gw.coef(idom,6) * gw.h(iGlob-gdom.nxhc*gdom.nyhc,1);
                    }
                }
                else if (gbc.bctypeZM == SUB_BC_Q_CONST)    {
                    gw.coef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / gdom.dz(iGlob);
                    if (gdom.isRain == 1) {gw.coef(idom,7) -= gdom.dt * gdom.qrain(iGlobSW) / gdom.dz(iGlob);}
                    else {gw.coef(idom,7) += gdom.dt * gbc.qbcZM / gdom.dz(iGlob);}
                    gw.coef(idom,6) = 0.0;
                }
                else    {
                    gw.coef(idom,6) = gw.coef(idom,6) * 2.0;
                    gw.coef(idom,7) += gw.coef(idom,6) * gw.h(iGlob-gdom.nxhc*gdom.nyhc,1);
                }
            }
            else if (kk == gdom.nz-1)   {
                gw.coef(idom,5) = gw.coef(idom,5) * 2.0;
                gw.coef(idom,7) += gw.coef(idom,5) * gw.h(iGlob+gdom.nxhc*gdom.nyhc,1);
            }
            gw.coef(idom,0) -= (gw.coef(idom,1)+gw.coef(idom,2)+gw.coef(idom,3)+gw.coef(idom,4)+gw.coef(idom,5)+gw.coef(idom,6));

            // real fac = 1e5;
            // printf(" -%d,%d- : (%f, %f, %f, %f, %f) --> %f --> %f\n",ii,kk,fac*gw.coef(idom,2),
            //     fac*gw.coef(idom,6),fac*gw.coef(idom,0),fac*gw.coef(idom,5),fac*gw.coef(idom,1),fac*gw.coef(idom,7),
            //     gw.h(iGlob,1));
            // if (ii == gdom.nx-1)    {printf(" ----- \n");}

        });
        // printf(" ----- \n\n");
        // Insert coefficients into Matrix A
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, irow = A.ptr(idom);
			gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
        	if (kk > 0)	{A.ind(irow) = idom - gdom.nx*gdom.ny;	A.val(irow) = gw.coef(idom,6);  irow++;}
        	if (jj > 0)	{A.ind(irow) = idom - gdom.nx;	        A.val(irow) = gw.coef(idom,4);  irow++;}
        	if (ii > 0)	{A.ind(irow) = idom - 1;		        A.val(irow) = gw.coef(idom,2);  irow++;}
        	A.ind(irow) = idom;	A.val(irow) = gw.coef(idom,0);	irow++;
        	if (ii < gdom.nx-1)	{A.ind(irow) = idom + 1;		        A.val(irow) = gw.coef(idom,1);  irow++;}
        	if (jj < gdom.ny-1)	{A.ind(irow) = idom + gdom.nx;	        A.val(irow) = gw.coef(idom,3);  irow++;}
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
    inline void update_wc(GwState &gw, State &state, GwDomain &gdom)	{
        // Update wc with explicit scheme
        if (gdom.gw_scheme == 1)    {
            Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob, ivg, iGlobSW;
                real coef, qqx, qqy, qqz, wcs, ss = 1e-5;
                gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                ivg = gw.soilID(iGlob) * NVG;
                wcs = gw.vgTable(ivg+2);
                coef = 1.0 + ss*(gw.h(iGlob,1) - gw.h(iGlob,0)) / wcs;
                qqx = (gw.q(iGlob,0) - gw.q(iGlob-1,0)) / gdom.dx;
                qqy = (gw.q(iGlob,1) - gw.q(iGlob-gdom.nxhc,1)) / gdom.dy;
                qqz = (gw.q(iGlob,2) - gw.q(iGlob-gdom.nxhc*gdom.nyhc,2)) / gdom.dz(iGlob);
                gw.wc(iGlob,1) = (gw.wc(iGlob,0) + gdom.dt * (qqx + qqy + qqz)) / coef;
                gw.wc(iGlob,2) = 0.0;
            });
            // Choose h or wc at the interface
            Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob, iGlobSW, ivg, flag;
        		real wcs, wcr, wcm, n, m, alpha, sbar;
                gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                iGlobSW = packIndicesUniformGrid(gdom.nyhc, gdom.nxhc, jj+hc, ii+hc);
                ivg = gw.soilID(iGlob) * NVG;
                wcs = gw.vgTable(ivg+2);     wcr = gw.vgTable(ivg+3);
                n = gw.vgTable(ivg+4);       alpha = gw.vgTable(ivg+6);
                m = 1.0 - 1.0 / n;
                wcm = wcr + (wcs-wcr)*pow((1.0 + pow(fabs(gdom.aev)*alpha,n)), m);

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
                    // if (kk == 0 && gw.h(iGlob-gdom.nxhc*gdom.nyhc,1) >= 0.0)  {flag = 1;}

                    if (kk == 0)    {
                        if (gw.wc(iGlob,1) < wcs+TOL8NEG && gw.h(iGlob-gdom.nxhc*gdom.nyhc,1) >= 0.0)   {flag = 0;}
                        else {flag = 1;}
                    }

                    if (flag == 1)  {
                        real tmp = gw.wc(iGlob,1);
                        sbar = pow(1.0 + pow(fabs(alpha*gw.h(iGlob,1)), n), -m);
                        if (gw.h(iGlob,1) > gdom.aev)   {gw.wc(iGlob,1) = wcs;}
                        else {gw.wc(iGlob,1) = sbar * (wcm - wcr) + wcr;}
                        gw.wc(iGlob,2) = tmp - gw.wc(iGlob,1);
                    }
                    else    {
                        if (gw.wc(iGlob,1) < wcs)   {
                            if (gw.wc(iGlob,1) < wcr)   {gw.wc(iGlob,1) = wcr;}
                            else {
                                gw.h(iGlob,1) = -(1.0/alpha) * (pow(pow((wcm-wcr)/(gw.wc(iGlob,1)-wcr),(1/m)) - 1.0, 1/n));
                            }
                        }
                        else {gw.h(iGlob,1) = 0.0;}
                    }
                }
            });
        }
        else {
            Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob, ivg, flag;
        		real wcs, wcr, wcm, n, m, alpha, sbar;
                gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                ivg = gw.soilID(iGlob) * NVG;
                wcs = gw.vgTable(ivg+2);     wcr = gw.vgTable(ivg+3);
                n = gw.vgTable(ivg+4);       alpha = gw.vgTable(ivg+6);
                m = 1.0 - 1.0 / n;
                wcm = wcr + (wcs-wcr)*pow((1.0 + pow(fabs(gdom.aev)*alpha,n)), m);
                sbar = pow(1.0 + pow(fabs(alpha*gw.h(iGlob,1)), n), -m);
                if (gw.h(iGlob,1) > gdom.aev)   {gw.wc(iGlob,1) = wcs;}
                else {gw.wc(iGlob,1) = sbar * (wcm - wcr) + wcr;}
                if (gw.wc(iGlob,1) > wcs)	{gw.wc(iGlob,2) += (gw.wc(iGlob,1)-wcs); gw.wc(iGlob,1) = wcs;}
        		// else if (gw.wc(iGlob,1) < wcr+0.001)	{gw.wc(iGlob,1) = wcr+0.001;}

            });
        }

    }
    // /* --------------------------------------------------
    //     End water content block
    // -------------------------------------------------- */





    // /* --------------------------------------------------
    //     Update dt based on either water content or iteration
    // -------------------------------------------------- */
    inline void dt_waco(GwState &gw, GwDomain &gdom)	{
    	real dwc_max, dt_old;
    	dt_old = gdom.dt;
        Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA (int idx, real &tmp) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            real dwc = fabs(gw.wc(iGlob,1) - gw.wc(iGlob,0));
			tmp = (dwc > tmp) ? dwc : tmp;
		} , Kokkos::Max<real>(dwc_max) );
    	if (dwc_max > 0.02)	{gdom.dt = gdom.dt * 0.9;}
    	else if (dwc_max >= 0.0 & dwc_max < 0.01)	{gdom.dt = gdom.dt * 1.1;}
    	if (gdom.dt > gdom.dt_max)	{gdom.dt = gdom.dt_max;}
    	else if (gdom.dt < gdom.dt_init)	{gdom.dt = gdom.dt_init;}
    	// if (gdom.dt != dt_old)	{
        //    std::cerr<< GOK "  dt is changed to" << gdom.dt << " (Max. d_wc = "<< dwc_max <<")\n";
    	// }
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
            real dwc = fabs(gw.h(idx,1) - gw.h(idx,0));
			tmp = (dwc > tmp) ? dwc : tmp;
		} , Kokkos::Max<real>(eps) );
        return eps;
    }

    // /* --------------------------------------------------
    //     End eps computation
    // -------------------------------------------------- */

    // /* --------------------------------------------------
    //     Integrator
    // -------------------------------------------------- */
    inline real integrate(GwState &gw, GwDomain &gdom)	{
    	real V_tot;
        Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA (int idx, real &tmp) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            tmp += gw.wc(iGlob,1) * gdom.dx * gdom.dx * gdom.dz(iGlob);
		} , Kokkos::Sum<real>(V_tot) );
        return V_tot;
    }
    inline real get_Vexchange(State &state, GwDomain &gdom)	{
    	real V_exch;
        Kokkos::parallel_reduce(gdom.ny*gdom.nx, KOKKOS_LAMBDA (int idx, real &tmp) {
            int ii, jj, iGlob;
            // printf(" nx, ny = %d, %d : idx = %d \n",gdom.nx, gdom.ny, idx);
            unpackIndicesUniformGrid(idx, gdom.ny, gdom.nx, jj, ii);
            // gdom.unpackIndices(idx, jj, ii);
            iGlob = (hc+jj)*gdom.nxhc + ii + hc;
            tmp += state.qss(iGlob) * gdom.dx * gdom.dx * gdom.dt;
		} , Kokkos::Sum<real>(V_exch) );
        return V_exch;
    }

    // /* --------------------------------------------------
    //     End integrating computation
    // -------------------------------------------------- */




};

#endif
