#ifndef _RT_INTEGRATOR_H_
#define _RT_INTEGRATOR_H_

#if SERGHEI_SUBSURFACE_TRANSPORT

#include "RTState.h"
#include "SourceSink.h"
#include "GwDomain.h"
#include "RTBC.h"
#include "Indexing.h"

class RTIntegrator	{

    Kokkos::Timer timer;

    public:
        real LiquidMassTot, LiquidMassTot_glob;
        real SolidMassTot, SolidMassTot_glob;
        real MassExch, MassExch_glob;

        SourceSink *ss;
        RTState *rt;
        GwDomain *gdom;
        std::vector<RTBC>* rtbc;

        int ncellsBC, ncellsBC_glob, ncellsIT, ncellsIT_glob;
        real QMassInBC, QMassOutBC, QMassInSS, QMassOutSS;
        real QMassInBC_glob, QMassOutBC_glob, QMassInSS_glob, QMassOutSS_glob;
        real QMassReaction, QMassReaction_glob;

        void initialize(RTState &rt_, GwDomain &gdom_, std::vector<RTBC> &rtbc_){
            rt = &rt_;
            gdom = &gdom_;
            rtbc = &rtbc_;
            // ss = &rtss_;
        }

        void integrate(GwState const &gw, RTState const &rt, GwDomain const &gdom, std::vector<RTBC> &rtbc)	{
        	int ierr=0;

        	// get total mass of liquid phase solute in the system [mg / L * m3]
        	LiquidMassTot = 0;
            
        	Kokkos::parallel_reduce(gdom.nCell , KOKKOS_LAMBDA (int idx, real &tmp) {
                
        		int ii, jj, kk, iGlob;
            	gdom.unpackIndices(idx, kk, jj, ii);
            	iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                real dxi = gdom.dx;
			    real dyi = gdom.dy;
			    real dzi = gdom.dz(iGlob);
			    real Vol = dxi * dyi * dzi;
                //浓度单位为mg/L，体积单位为m3，质量单位为mg，temp单位转换为mg
            	tmp += gw.wc_new(iGlob) * rt.c(iGlob, 1) * Vol * 1000.0;
            }, Kokkos::Sum<real>(LiquidMassTot) );
            Kokkos::fence();

            LiquidMassTot_glob = 0.0;
            ierr=MPI_Allreduce(&LiquidMassTot, &LiquidMassTot_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            //debug:输出总质量
            // printf("Total solute mass in the system: %e \n",LiquidMassTot_glob);
            
    // get total mass of solid phase solute in the system 
            SolidMassTot = 0;
            Kokkos::parallel_reduce(gdom.nCell , KOKKOS_LAMBDA (int idx, real &tmp) {

                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                real dxi = gdom.dx;
                real dyi = gdom.dy;
                real dzi = gdom.dz(iGlob);
                real Vol = dxi * dyi * dzi;
                //浓度单位为mg/L，体积单位为m3，质量单位为mg，temp单位转换为mg
                tmp += rt.c_solid(iGlob, 1) * Vol*1000.0;
            
            }, Kokkos::Sum<real>(SolidMassTot) );
            Kokkos::fence();
            SolidMassTot_glob = 0.0;
            ierr=MPI_Allreduce(&SolidMassTot, &SolidMassTot_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

    // get surface-subsurface solute mass exchange rate

            MassExch = 0;
           #if SERGHEI_SURFACE_TRANSPORT
            Kokkos::parallel_reduce(gdom.ny*gdom.nx, KOKKOS_LAMBDA (int idx, real &tmp) {
                int ii, jj, iGlob;
                unpackIndicesUniformGrid(idx, gdom.ny, gdom.nx, jj, ii);
                iGlob = jj*gdom.nx + ii;
                tmp += rt.ConQss(iGlob) * gdom.dx * gdom.dx;
            } , Kokkos::Sum<real>(MassExch) );
            Kokkos::fence();
            #endif
            MassExch_glob = 0.0;
            ierr=MPI_Allreduce(&MassExch, &MassExch_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            // boundary mass flow
            QMassInBC = 0.0;
            QMassOutBC = 0.0;
            for (int k = 0; k < rtbc.size(); k++) {
                int _ncellsBC;
                MPI_Allreduce(&(rtbc[k].ncellsBC), &_ncellsBC, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                ncellsBC += _ncellsBC;

                //质量通量，单位为 mg/s
                QMassInBC += rtbc[k].QMassInflow*1000;
                QMassOutBC += rtbc[k].QMassOutflow*1000;
            }
            QMassInBC_glob=0.0;
            QMassOutBC_glob=0.0;
            ncellsBC_glob = 0;
            ierr=MPI_Allreduce(&QMassInBC, &QMassInBC_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
            ierr=MPI_Allreduce(&QMassOutBC, &QMassOutBC_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
            ierr=MPI_Allreduce(&ncellsBC, &ncellsBC_glob, 1, MPI_INT , MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            //debug:输出边界质量流量
            // printf("Total solute mass flow rate at the boundaries: Inflow = %e , Outflow = %e \n",QMassInBC_glob, QMassOutBC_glob);
        
    //计算反应项质量消耗
            QMassReaction = 0.0;
            Kokkos::parallel_reduce(gdom.nCell , KOKKOS_LAMBDA (int idx, real &tmp) {

                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                real dxi = gdom.dx;
                real dyi = gdom.dy;
                real dzi = gdom.dz(iGlob);
                real Vol = dxi * dyi * dzi;
                // 计算本时间步内消耗的质量 
                // !(mg)
                tmp += rt.lambda * gw.wc_new(iGlob) * rt.c(iGlob, 1) * Vol * rt.dt * 1000.0;

            }, Kokkos::Sum<real>(QMassReaction) );
            Kokkos::fence();
            QMassReaction_glob = 0.0;
            ierr=MPI_Allreduce(&QMassReaction, &QMassReaction_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            //debug:输出反应消耗质量
            // printf("Total solute mass reacted in the system: %e \n",QMassReaction_glob);    
        
        }

};
#endif

#endif