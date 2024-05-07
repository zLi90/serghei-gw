#ifndef _RT_FUNCTION_H_
#define _RT_FUNCTION_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
#include "GwState.h"
#include "RTMatrix.h"
#include "RTSolver.h"
#include "State.h"
#include <set>

class RTFunction
{

private:
	Kokkos::Timer timer;

public:
	/* --------------------------------------------------
		Top-level PCA solver
	-------------------------------------------------- */
	template <typename execution_space, typename type_solver>
	inline void rt_pca_solve(RTState &rt, RTMatrix &rtA, GwState &gw, GwDomain &gdom, GwMPI &gmpi, Parallel &par, type_solver &rtsolver)
	{
		int iter, ierr=1; 
		real dt_tmp;

		//扩散系数计算
		dispersion_tensor(rt, gw, gdom);

		RTlinear_system(gw, gdom,  rtA, par, rt);

        	timer.reset();
		#if SERGHEI_KOKKOSKERNELS_SOLVER
		rtsolver.kkpcg(rtA);
		#else
		iter = rtsolver.cg(rtA, gdom);
		#endif
		gdom.timers.solver += timer.seconds();	

		// Update concentration	
        Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            rt.c(iGlob,1) = rtA.rt_x(idom);
        });
        gmpi.mpi_sendrecv(rt.c, gdom, par);

	   // Update time step
        dt_tmp = gdom.dt;

        ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            rt.c(iGlob,0) = rt.c(iGlob,1);
        });
	    gmpi.mpi_sendrecv(rt.c, gdom, par);

	}

    /* --------------------------------------------------
        Top-level Picard solver
    -------------------------------------------------- */
	template <typename execution_space, typename type_solver>
    inline void rt_picard_solve(RTState &rt, RTMatrix &rtA,  GwState &gw, GwDomain &gdom, 
             type_solver &rtsolver,  GwMPI &gmpi, GwIntegrator &gint, Parallel &par)  {
        int iter, iter_cg, iter_max = 50, ierr=1;
        real eps_diff = 1.0, eps_tmp, eps_old = 1.0, eps = 1.0, eps_min = 5e-6, dt_tmp;
	//施加边界条件
     //    for (int k = 0; k < gbc.size(); k++) {
     //        gbc[k].applyHBC(gw, gdom, par);
     //    }
		// for (int i = 225; i <= 250; ++i) {
    			// rt.c(i, 1) = 15;
		// 	}
		// rt.c(130, 1) = 15;
            Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
		// 	std::cout << "-------gdom.nxhc----- " <<gdom.nxhc << std::endl;
		// 	std::cout << "-------gdom.nyhc----- " <<gdom.nyhc << std::endl;
		// std::cout << "-------hc----- " <<hc << std::endl;
			
            });
		
		// rt.c(55, 1) = 10;
		dispersion_tensor(rt, gw, gdom);
        iter = 0;
        while (iter < iter_max && eps_diff/eps_old > eps_min && eps > eps_min) {
		
		  
            RTlinear_system(gw, gdom,  rtA, par, rt);
		//   std::cout << "---------RTdt:----------- " << std::endl;
            timer.reset();
            #if SERGHEI_KOKKOSKERNELS_SOLVER
            rtsolver.kkpcg(rtA);
            #else
            iter = rtsolver.cg(rtA, gdom);
            #endif
            gdom.timers.solver += timer.seconds();
            Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                rt.c(iGlob,0) = rt.c(iGlob,1);
                rt.c(iGlob,1) = rtA.rt_x(idom);
            });
            gmpi.mpi_sendrecv(rt.c, gdom, par);



            eps_old = eps;
            eps = get_eps(rt, gdom);
            eps_tmp = fabs(eps_old - eps);
            ierr = MPI_Allreduce(&eps_tmp, &eps_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


            gmpi.mpi_sendrecv(rt.c, gdom, par);

            iter += 1;
		          // std::cout << "---------RTdt:----------- " <<rt.dt << std::endl;  
				// std::cout << "---------GWdt:----------- " << gdom.dt << std::endl;    
		//   std::cout << "---------RTiter2222:----------- " << iter << std::endl;
        }
        // printf("    > Picard loop converges in %d iterations with eps = %f, %f\n",iter,eps,eps_diff);

        dt_iter(rt, gw, gdom, iter);
        dt_tmp = rt.dt;
        ierr = MPI_Allreduce(&dt_tmp, &rt.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            rt.c(iGlob,0) = rt.c(iGlob,1);  
        });
	// std::cout << "---------rt.dt:----------- " << rt.dt << std::endl;
     //    gint.integrate(gw, gdom, gbc, gss);
    }
	
	// Calculation of the hydrodynamic dispersion coefficient tensor
	inline void dispersion_tensor(RTState &rt, GwState &gw, GwDomain &gdom)
	{
			Kokkos::parallel_for(
		    gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			    // rt.c(iGlob, 0) 表示下一时间步长的浓度值
			    // rt.c(iGlob, 1) 表示当前步长的浓度值
				real phi = 0.45;//孔隙度
				//平均孔隙流速计算,x,y,z和平均流速方向
				rt.aveV(iGlob, 0) = gw.q(iGlob,0)/phi;
				rt.aveV(iGlob, 1) = gw.q(iGlob,1)/phi;
				rt.aveV(iGlob, 2) = gw.q(iGlob,2)/phi;
				rt.aveV(iGlob, 3) = sqrt(pow(rt.aveV(iGlob, 0), 2) + pow(rt.aveV(iGlob, 1), 2) + pow(rt.aveV(iGlob, 2), 2));
			
				// std::cout << "-------rt.aveV(iGlob, 3)----- " <<rt.aveV(iGlob, 3) << std::endl;
				// 机械弥散系数计算
				//Dxx=DT(qy+qz)+DL*qx
				rt.dcal(iGlob,0) = (0.3*(pow(rt.aveV(iGlob, 1), 2)+pow(rt.aveV(iGlob, 2), 2))+6*pow(rt.aveV(iGlob, 0), 2))/rt.aveV(iGlob,3);
				rt.dcal(iGlob,1) = (0.3*(pow(rt.aveV(iGlob, 0), 2)+pow(rt.aveV(iGlob, 2), 2))+6*pow(rt.aveV(iGlob, 1), 2))/rt.aveV(iGlob,3);
				rt.dcal(iGlob,2) = (0.3*(pow(rt.aveV(iGlob, 0), 2)+pow(rt.aveV(iGlob, 1), 2))+6*pow(rt.aveV(iGlob, 2), 2))/rt.aveV(iGlob,3);
				// std::cout << "-------rt.d(iGlob,0)----- " <<rt.d(iGlob,0) << std::endl;
				//将rt.d中nan值替换
				Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
					for (int i = 0; i < 3; ++i) {
						if (std::isnan(rt.dcal(iGlob,i))) {
							rt.dcal(iGlob,i) = 0;
							
						}
					}
					
					});
				//对gw.wc进行处理,将0值替换为ThetaR  : 0.029 
				// Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
				// 	// for (int i = 0; i < 3; ++i) {
				// 		if (std::isnan(gw.wc(iGlob,0)) || gw.wc(iGlob,0) == 0) {
				// 			gw.wc(iGlob,0) = 0.029;
							
				// 		}
				// 	// }
				// 	});
					// std::cout << "-------gw.wc(iGlob,0)----- " <<gw.wc(iGlob,0) << std::endl;
				//弥散系数=机械弥散+分子弥散
				//1e-5为分子弥散系数
				// rt.dcal(iGlob,0) += 1e-3*(pow(gw.wc(iGlob,0), 7.0/3.0) / pow(phi, 2.0));
				// rt.dcal(iGlob,1) += 1e-3*(pow(gw.wc(iGlob,0), 7.0/3.0) / pow(phi, 2.0));
				// rt.dcal(iGlob,2) += 1e-3*(pow(gw.wc(iGlob,0), 7.0/3.0) / pow(phi, 2.0));
				
				rt.dcal(iGlob,0) += 1e-7;
				rt.dcal(iGlob,1) += 1e-7;
				rt.dcal(iGlob,2) += 1e-7;
				// std::cout << "-------rt.dcal(iGlob,0)----- " <<rt.dcal(iGlob,2) << std::endl;	
			});
	}
	    /* --------------------------------------------------
        Get RT matrix coefficients
    -------------------------------------------------- */
    inline void RTlinear_system(GwState &gw, GwDomain &gdom, RTMatrix &rtA, Parallel &par, RTState &rt)	{
        // Calculate matrix coefficients
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, ivg, iGlob, iGlobSW;
		  
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            iGlobSW = (hc+jj)*gdom.nxhc + ii + hc;

		rt.RTcoef(idom, 0) = -(1 + 2 * (rt.dcal(iGlob,0) * rt.dt / (gdom.dx * gdom.dx) + rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob))));
		rt.RTcoef(idom, 1) = rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ) - rt.aveV(iGlob, 0) * rt.dt / (2 * gdom.dx );
          rt.RTcoef(idom, 2) = rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ) + rt.aveV(iGlob, 0) * rt.dt / (2 * gdom.dx );
          rt.RTcoef(idom, 3) = rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) - rt.aveV(iGlob, 1) * rt.dt / (2 * gdom.dy);
          rt.RTcoef(idom, 4) = rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.aveV(iGlob, 1) * rt.dt / (2 * gdom.dy);
          rt.RTcoef(idom, 5) = rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)) - rt.aveV(iGlob, 2) * rt.dt / (2 * gdom.dz(iGlob));
          rt.RTcoef(idom, 6) = rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)) + rt.aveV(iGlob, 2) * rt.dt / (2 * gdom.dz(iGlob));
		rt.RTcoef(idom, 7) = -rt.c(iGlob, 1);//pca迭代方法中的右端项为已知浓度值Cn
		
		// rt.RTcoef(idom, 0) = -gw.wc(iGlob,0)*((1 + 2 * (rt.dcal(iGlob,0) * rt.dt / (gdom.dx * gdom.dx) + rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)))));
		// rt.RTcoef(idom, 1) = gw.wc(iGlob,0)*(rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ) - rt.aveV(iGlob, 0) * rt.dt / (2 * gdom.dx ));
          // rt.RTcoef(idom, 2) = gw.wc(iGlob,0)*(rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ) + rt.aveV(iGlob, 0) * rt.dt / (2 * gdom.dx ));
          // rt.RTcoef(idom, 3) = gw.wc(iGlob,0)*(rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) - rt.aveV(iGlob, 1) * rt.dt / (2 * gdom.dy));
          // rt.RTcoef(idom, 4) = gw.wc(iGlob,0)*(rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.aveV(iGlob, 1) * rt.dt / (2 * gdom.dy));
          // rt.RTcoef(idom, 5) = gw.wc(iGlob,0)*(rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)) - rt.aveV(iGlob, 2) * rt.dt / (2 * gdom.dz(iGlob)));
          // rt.RTcoef(idom, 6) = gw.wc(iGlob,0)*(rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)) + rt.aveV(iGlob, 2) * rt.dt / (2 * gdom.dz(iGlob)));
		// rt.RTcoef(idom, 7) = -gw.wc(iGlob,0)*rt.c(iGlob, 1);//pca迭代方法中的右端项为已知浓度值Cn


		if (gdom.gw_scheme != 1)	//picard迭代
		{
			// rt.RTcoef(idom, 7) 需要修改
			// rt.RTcoef(idom, 7) = -gw.wc(iGlob,0)*rt.c(iGlob, 1);
			rt.RTcoef(idom, 7) = -rt.c(iGlob, 1);	
		}
		
	}   );
		 /*--------------后续完善部分
		  // Apply internal boundary conditions (needed when MPI is used)

		// Apply outer boundary conditions
			for (int k = 0; k < gbc.size(); k++) {
				gbc[k].applyMatBC(gw, gdom, par);
			}
			-------------后续完善部分*/
    
		// Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
          //   rt.RTcoef(idom,0) = rt.RTcoef(idom,0);
       	//  });//矩阵A对角线项H(i,j,k)(n+1)的系数
		
		/*--------------后续完善部分
		// Apply internal source/sink terms
		for (int k = 0; k < gss.size(); k++) {
            gss[k].applyMatSS(gw, gdom);
        	}
		-------------后续完善部分*/

					 			
		  // Insert coefficients into Matrix rtA.

	   Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
			
			int ii, jj, kk;
			// idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;

			int irow = rtA.rt_ptr(idom); 
			
			
			gdom.unpackIndices(idom, kk, jj, ii);

		     if (kk > 0)	{rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny;	rtA.rt_val(irow) = rt.RTcoef(idom,6);  irow++;}
        		if (jj > 0)	{rtA.rt_ind(irow) = idom - gdom.nx;	        rtA.rt_val(irow) = rt.RTcoef(idom,4);  irow++;}
			if (ii > 0)	{rtA.rt_ind(irow) = idom - 1;		        rtA.rt_val(irow) = rt.RTcoef(idom,2);  irow++;}
			rtA.rt_ind(irow) = idom;	rtA.rt_val(irow) = rt.RTcoef(idom,0);	irow++;
			if (ii < gdom.nx-1)	{rtA.rt_ind(irow) = idom + 1;		        rtA.rt_val(irow) = rt.RTcoef(idom,1);  irow++;}
			if (jj < gdom.ny-1)	{rtA.rt_ind(irow) = idom + gdom.nx;	        rtA.rt_val(irow) = rt.RTcoef(idom,3);  irow++;}
			if (kk < gdom.nz-1)	{rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny;	rtA.rt_val(irow) = rt.RTcoef(idom,5);  irow++;}
		
			rtA.rt_rhs(idom) = rt.RTcoef(idom,7);

		   });
	}

/*----------------20240426修改-----------*///非迭代方法时间步长控制
    inline void dt_waco(RTState &rt, GwState &gw, GwDomain &gdom)	{
    	real dwc_max, dt_old;
    	dt_old = rt.dt;//在serghei.h中初始化rt.dt=gdom.dt_init;
        Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA (int idx, real &tmp) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idx, kk, jj, ii);
            // gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            real dwc = fabs(gw.wc(iGlob,1) - gw.wc(iGlob,0));
			tmp = (dwc > tmp) ? dwc : tmp;
		} , Kokkos::Max<real>(dwc_max) );
    	if (dwc_max > 0.02)	{gdom.dt = gdom.dt * 0.9;}
    	else if (dwc_max >= 0.0 & dwc_max < 0.01)	{gdom.dt = gdom.dt * 1.1;}
    	if (gdom.dt > gdom.dt_max)	{gdom.dt = gdom.dt_max;}
    	else if (gdom.dt < gdom.dt_init)	{gdom.dt = gdom.dt_init;}
    }
/*----------------20240426修改-----------*///非迭代方法时间步长控制

    inline void dt_iter(RTState &rt, GwState &gw, GwDomain &gdom, int iter)	{
    	real dt_old;
    	dt_old = rt.dt;
        if (iter < 7)   {rt.dt = rt.dt * 1.1;}
        else if (iter > 11)  {rt.dt = rt.dt * 0.9;}
        if (rt.dt > gdom.dt_max)	{rt.dt = gdom.dt_max;}
    	else if (rt.dt < gdom.dt_init)	{rt.dt = gdom.dt_init;}
    }

    inline real get_eps(RTState &rt, GwDomain &gdom)	{
    	real eps;
        Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA (int idx, real &tmp) {
            real dwc = fabs(rt.c(idx,1) - rt.c(idx,0));
			tmp = (dwc > tmp) ? dwc : tmp;
		} , Kokkos::Max<real>(eps) );
        return eps;
    }


};
#endif
