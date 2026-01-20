#ifndef _GW_INTEGRATOR_H_
#define _GW_INTEGRATOR_H_

#if SERGHEI_SUBSURFACE_MODEL

#include "GwState.h"
#include "SourceSink.h"
#include "GwDomain.h"
#include "GwBC.h"
#include "Indexing.h"

class GwIntegrator	{

	Kokkos::Timer timer;

	public:
		real Vtot, Vtot_glob;
		real Vexch, Vexch_glob;

		SourceSink *ss;
		GwState *gw;
		GwDomain *gdom;
		std::vector<GwBC>* gwbc;
		std::vector<GwSS>* gwss;

		int ncellsBC, ncellsBC_glob, ncellsIT, ncellsIT_glob;
		real QinBC, QoutBC, QinSS, QoutSS, QoutSS_qe;
		real QinBC_glob, QoutBC_glob, QinSS_glob, QoutSS_glob, QoutSS_qe_glob;


		void initialize(GwState &gw_, GwDomain &gdom_, std::vector<GwBC> &gwbc_, std::vector<GwSS> &gwss_){
			gw = &gw_;
			gdom = &gdom_;
			gwss = &gwss_;
			gwbc = &gwbc_;
		}

		void integrate(GwState const &gw, GwDomain const &gdom, std::vector<GwBC> &gwbc, std::vector<GwSS> &gwss)	{
	 		int ierr=0;
	 		// get total volume
	 		Vtot = 0;
	  		Kokkos::parallel_reduce(gdom.nCell , KOKKOS_LAMBDA (int idx, real &tmp) {
				
	  			int ii, jj, kk, iGlob;
            	gdom.unpackIndices(idx, kk, jj, ii);
            	iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            	tmp += gw.wc(iGlob,1) * gdom.dx * gdom.dx * gdom.dz(iGlob);
    		} , Kokkos::Sum<real>(Vtot) );
			Kokkos::fence();
			Vtot_glob = 0.0;
			ierr=MPI_Allreduce(&Vtot, &Vtot_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

//! zzb 20241213 根系区含水率均值计算函数
// #if CROP_GROWTH_MODEL
    // 计算根区平均含水率
    for (int k = 0; k < gwss.size(); k++) {
        real wc_sum = 0.0; // 根区含水率总和
        int ncells_root = 0; // 根区单元数量

        // 获取最大根深所在层数
        int k_root_max = gwss[k].k_max_root;

        // 遍历 0 到 k_root_max 层的网格单元，计算根区平均含水率
        for (int kk = 0; kk <= k_root_max; kk++) {
            for (int jj = 0; jj < gdom.ny; jj++) {
                for (int ii = 0; ii < gdom.nx; ii++) {
                    int iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
                    wc_sum += gw.wc(iGlob, 0); // 累加含水率
                    ncells_root++; // 统计根区单元数量
                }
            }
        }

        // 计算根区平均含水率
        if (ncells_root > 0) {
            gwss[k].wc_root_zone = wc_sum / ncells_root;
			// printf("wc_sum: %f\n, ncells_root: %d\n, wc_root_zone: %f\n", wc_sum, ncells_root, gwss[k].wc_root_zone);
        } else {
            gwss[k].wc_root_zone = 0.0; // 如果没有根区单元，平均含水率为 0
        }
		// printf("ncells_root: %d\n, wc_root_zone: %f\n", ncells_root, gwss[k].wc_root_zone);
    }	
// #endif
//! zzb 20241213 根系区含水率均值计算函数

			// get surface-subsurface exchange rate [m3/s]
			Vexch = 0;
			#if SERGHEI_SWE_MODEL
			Kokkos::parallel_reduce(gdom.ny*gdom.nx, KOKKOS_LAMBDA (int idx, real &tmp) {
		        int ii, jj, iGlob;
		        unpackIndicesUniformGrid(idx, gdom.ny, gdom.nx, jj, ii);
		        iGlob = jj*gdom.nx + ii;
		        tmp += gw.qss(iGlob) * gdom.dx * gdom.dx;
			} , Kokkos::Sum<real>(Vexch) );
			Kokkos::fence();
			#endif
			Vexch_glob = 0.0;
			ierr=MPI_Allreduce(&Vexch, &Vexch_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			// boundary flow
			QinBC = 0.0;
			QoutBC = 0.0;
			for (int k = 0; k < gwbc.size(); k++) {
				int _ncellsBC;
				MPI_Allreduce(&(gwbc[k].ncellsBC), &_ncellsBC, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				ncellsBC += _ncellsBC;
				QinBC += gwbc[k].Qinflow;
				QoutBC += gwbc[k].Qoutflow;
			}
			QinBC_glob=0.0;
			QoutBC_glob=0.0;
        	ncellsBC_glob = 0;
			ierr=MPI_Allreduce(&QinBC, &QinBC_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
			ierr=MPI_Allreduce(&QoutBC, &QoutBC_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
		    ierr=MPI_Allreduce(&ncellsBC, &ncellsBC_glob, 1, MPI_INT , MPI_SUM, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			// source/sink terms
			QinSS = 0.0;
			QoutSS = 0.0;
			QoutSS_qe = 0.0;//!zzb 20250402 实际蒸发量
			for (int k = 0; k < gwss.size(); k++) {
				int _ncellsIT;
				MPI_Allreduce(&(gwss[k].ncellsIT), &_ncellsIT, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				ncellsIT += _ncellsIT;
				QinSS += gwss[k].Qinflow;
				QoutSS += gwss[k].Qoutflow;
				QoutSS_qe += gwss[k].Qoutflow_qe;//!zzb 20250402 实际蒸发量
				

			}
			QinSS_glob=0.0;
			QoutSS_glob=0.0;
			QoutSS_qe_glob=0.0;//!zzb 20250402 实际蒸发量
        	ncellsIT_glob = 0;
			ierr=MPI_Allreduce(&QinSS, &QinSS_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
			ierr=MPI_Allreduce(&QoutSS, &QoutSS_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
			ierr=MPI_Allreduce(&QoutSS_qe, &QoutSS_qe_glob, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
			ierr=MPI_Allreduce(&ncellsIT, &ncellsIT_glob, 1, MPI_INT , MPI_SUM, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

		}


};

#endif

#endif
