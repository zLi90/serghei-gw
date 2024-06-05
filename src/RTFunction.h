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
#include <limits> 

class RTFunction
{

private:
	Kokkos::Timer timer;

public:
	/* --------------------------------------------------
		Top-level PCA solver
	-------------------------------------------------- */
	template <typename execution_space, typename type_solver>
	inline void rt_pca_solve(RTState &rt, RTMatrix &rtA, GwState &gw, GwDomain &gdom, std::vector<RTBC> &rtgbc, GwMPI &gmpi, Parallel &par, type_solver &rtsolver)
	{
		int iter, ierr=1; 
		real dt_tmp;

		//施加边界条件
        for (int k = 0; k < rtgbc.size(); k++) {
		  rtgbc[k].applyConcentrationBC(rt, gw, gdom, par); 		  
        }
		//扩散系数计算
		dispersion_tensor(rt, gw, gdom);

		RTlinear_system(gw, gdom,  rtA, par, rt, rtgbc);

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
		//施加边界条件
        for (int k = 0; k < rtgbc.size(); k++) {
		  rtgbc[k].applyConcentrationBC(rt, gw, gdom, par); 		  
        }
        gmpi.mpi_sendrecv(rt.c, gdom, par);
		//施加边界条件
        for (int k = 0; k < rtgbc.size(); k++) {
		  rtgbc[k].applyConcentrationBC(rt, gw, gdom, par); 		  
        }

        dt_con(rt, gw, gdom);//计算浓度残差值，并更新
	   // Update time step
        dt_tmp = gdom.dt;

        ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            rt.c(iGlob,0) = rt.c(iGlob,1);
        });
	//     gmpi.mpi_sendrecv(rt.c, gdom, par);

// // 调用函数找到最大值及其索引
// 	std::pair<real, int> max_value_pair = find_max_value(rt);
// 	real max_value = max_value_pair.first;
// 	int max_index = max_value_pair.second;
// 	std::cout << "最大浓度值: " << max_value << std::endl;
// 	std::cout << "对应的索引: " << max_index << std::endl;
// //调用函数找到最大Peclet
// 	std::pair<real, int> peclet_pair = find_peclet(rt,gdom);
// 	real peclet_value = peclet_pair.first;
// 	int peclet_index = peclet_pair.second;
// 	std::cout << "最大Peclet: " << peclet_value << std::endl;
// 	std::cout << "对应的索引: " << peclet_index << std::endl;  

// //调用函数找到最大孔隙均速
// 	std::pair<real, int> aveV_pair = find_aveV(rt);
// 	real aveV_value = aveV_pair.first;
// 	int aveV_index = aveV_pair.second;
// 	std::cout << "最大孔隙均速: " << aveV_value << std::endl;
// 	std::cout << "对应的索引: " << aveV_index << std::endl;  

	}

    /* --------------------------------------------------
        Top-level Picard solver
    -------------------------------------------------- */
	template <typename execution_space, typename type_solver>
    inline void rt_picard_solve(RTState &rt, RTMatrix &rtA,  GwState &gw, GwDomain &gdom, std::vector<RTBC> &rtgbc,
             type_solver &rtsolver, GwMPI &gmpi, GwIntegrator &gint, Parallel &par)  {


        int iter, iter_cg, iter_max = 100, ierr=1;
        real eps_diff = 1.0, eps_tmp, eps_old = 1.0, eps = 1.0,  dt_tmp;
	real eps_min = 5e-5;//eps_min = 5e-6,20240510修改了收敛残差标准
		//施加边界条件
        for (int k = 0; k < rtgbc.size(); k++) {
		  rtgbc[k].applyConcentrationBC(rt, gw, gdom, par); 		  
        }
            Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
		  });			
        iter = 0;
     //计算扩散系数
	dispersion_tensor(rt, gw, gdom); 	

        while (iter < iter_max && eps_diff/eps_old > eps_min && eps > eps_min) {

            RTlinear_system(gw, gdom,  rtA, par, rt, rtgbc);
		  
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

		dispersion_tensor(rt, gw, gdom); 

            eps_old = eps;
            eps = get_eps(rt, gdom);
            eps_tmp = fabs(eps_old - eps);
            ierr = MPI_Allreduce(&eps_tmp, &eps_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            iter += 1;
        }
        // printf("    > Picard loop converges in %d iterations with eps = %f, %f\n",iter,eps,eps_diff);

        dt_iter(rt, gw, gdom, iter);
        dt_tmp = rt.dt;
        ierr = MPI_Allreduce(&dt_tmp, &rt.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            rt.c(iGlob,0) = rt.c(iGlob,1); 
        });
     //    gint.integrate(gw, gdom, gbc, gss);

/*
// 调用函数找到最大值及其索引
	std::pair<real, int> max_value_pair = find_max_value(rt);
	real max_value = max_value_pair.first;
	int max_index = max_value_pair.second;
	std::cout << "最大浓度值: " << max_value << std::endl;
	std::cout << "对应的索引: " << max_index << std::endl;
//调用函数找到最大Peclet
	std::pair<real, int> peclet_pair = find_peclet(rt,gdom);
	real peclet_value = peclet_pair.first;
	int peclet_index = peclet_pair.second;
	std::cout << "最大Peclet: " << peclet_value << std::endl;
	std::cout << "对应的索引: " << peclet_index << std::endl;  

//调用函数找到最大孔隙均速
	std::pair<real, int> aveV_pair = find_aveV(rt);
	real aveV_value = aveV_pair.first;
	int aveV_index = aveV_pair.second;
	std::cout << "最大孔隙均速: " << aveV_value << std::endl;
	std::cout << "对应的索引: " << aveV_index << std::endl; 
	 */
}


// 寻找 rt.c(iGlob, 1) 中的最大值及其索引
	std::pair<real, int> find_max_value(const RTState &rt) {
	real max_value = 0;
	int max_index = 0;
	for (int i = 1; i < rt.c.extent(0); ++i) {
		real current_value = rt.c(i, 1);
		if (current_value > max_value) {
			max_value = current_value;
			max_index = i;
		}
	}
	return std::make_pair(max_value, max_index);
	}

	//计算Peclet
	std::pair<real, int> find_peclet(const RTState &rt, GwDomain &gdom) {
	real max_value = 0;
	int max_index = 0;
	for (int i = 1; i < rt.c.extent(0); ++i) {
		real current_value = fabs(gdom.dx*rt.aveV(i, 3)/rt.dcal(i,0));
		if (current_value > max_value) {
			max_value = current_value;
			max_index = i;
		}
	}
	return std::make_pair(max_value, max_index);
	}

	//计算最大孔隙均速
	std::pair<real, int> find_aveV(const RTState &rt) {
	real max_value = 0;
	int max_index = 0;
	for (int i = 1; i < rt.c.extent(0); ++i) {
		real current_value = fabs(rt.aveV(i, 3));
		if (current_value > max_value) {
			max_value = current_value;
			max_index = i;
		}
	}
	return std::make_pair(max_value, max_index);
	}

	
	// Calculation of the hydrodynamic dispersion coefficient tensor
	inline void dispersion_tensor(RTState &rt, GwState &gw, GwDomain &gdom)
	{
				Kokkos::parallel_for(
		    		gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			    // rt.c(iGlob, 0) 表示下一时间步长的浓度值
			    // rt.c(iGlob, 1) 表示当前步长的浓度值
				
//平均孔隙流速计算，平均孔隙流速计算,x,y,z和平均流速方向			
				int ii, jj, kk, ivg;
				real  wcs;
				gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
				ivg = gw.soilID(iGlob) * NVG;
				wcs = gw.vgTable(ivg+2);            			  
				//aveVB为边界流速				
				rt.aveVB(iGlob, 0) = -(gw.q(iGlob,0)/wcs);				
				rt.aveVB(iGlob, 1) = (gw.q(iGlob,1)/wcs);
				rt.aveVB(iGlob, 2) = -(gw.q(iGlob,2)/wcs);
				rt.aveVB(iGlob, 3) = sqrt(pow(rt.aveVB(iGlob, 0), 2) + pow(rt.aveVB(iGlob, 1), 2) + pow(rt.aveVB(iGlob, 2), 2));
				
				//aveV为节点流速
				rt.aveV(iGlob, 0) = (rt.aveVB(iGlob, 0) + rt.aveVB(iGlob-1, 0))/2;			
				rt.aveV(iGlob, 1) = (rt.aveVB(iGlob, 1) + rt.aveVB(iGlob-gdom.nxhc, 1))/2;				
				rt.aveV(iGlob, 2) = (rt.aveVB(iGlob, 2) + rt.aveVB(iGlob-gdom.nxhc*gdom.nyhc, 2))/2;
				rt.aveV(iGlob, 3) = sqrt(pow(rt.aveV(iGlob, 0), 2) + pow(rt.aveV(iGlob, 1), 2) + pow(rt.aveV(iGlob, 2), 2));
				
		
// 机械弥散系数计算
				
				//Dxx
				rt.dcal(iGlob,0) = (rt.alpha_T*(pow(rt.aveV(iGlob, 1), 2)+pow(rt.aveV(iGlob, 2), 2))+rt.alpha_L*pow(rt.aveV(iGlob, 0), 2))/rt.aveV(iGlob,3);
				//DYY
				rt.dcal(iGlob,1) = (rt.alpha_T*(pow(rt.aveV(iGlob, 0), 2)+pow(rt.aveVB(iGlob, 2), 2))+rt.alpha_L*pow(rt.aveV(iGlob, 1), 2))/rt.aveV(iGlob,3);
				// DZZ
				rt.dcal(iGlob,2) = (rt.alpha_T*(pow(rt.aveV(iGlob, 0), 2)+pow(rt.aveV(iGlob, 1), 2))+rt.alpha_L*pow(rt.aveV(iGlob, 2), 2))/rt.aveV(iGlob,3);				
				//Dxy
				rt.dcal(iGlob,3) = fabs((rt.alpha_L - rt.alpha_T)*rt.aveV(iGlob, 0)*rt.aveV(iGlob, 1)/rt.aveV(iGlob,3));
				//Dxz
				rt.dcal(iGlob,4) =fabs ((rt.alpha_L - rt.alpha_T)*rt.aveV(iGlob, 0)*rt.aveV(iGlob, 2)/rt.aveV(iGlob,3));
				//Dyz
				rt.dcal(iGlob,5) = fabs ((rt.alpha_L - rt.alpha_T)*rt.aveV(iGlob, 1)*rt.aveV(iGlob, 2)/rt.aveV(iGlob,3));


				
//将rt.d中nan值替换
				for (int i = 0; i < 6; ++i) {
					if (std::isnan(rt.dcal(iGlob,i))) {
							rt.dcal(iGlob,i) = 0;}};								
																										
				/*------弥散系数=机械弥散+分子弥散----------*/

				rt.dcal(iGlob,0) += gw.wc(iGlob,0)*rt.d_base;
				rt.dcal(iGlob,1) += gw.wc(iGlob,0)*rt.d_base;
				rt.dcal(iGlob,2) += gw.wc(iGlob,0)*rt.d_base;

//根据研究区域维数对水动力弥散系数进行调整
				if (gdom.nx == 1) 
				{
					rt.dcal(iGlob,0) = 0;
					rt.dcal(iGlob,4) = 0; 
					rt.dcal(iGlob,5) = 0;
					}				
				else if (gdom.ny == 1) 
				{
					rt.dcal(iGlob,1) = 0; 
					rt.dcal(iGlob,3) = 0; 
					rt.dcal(iGlob,5) = 0;
					}
				else if (gdom.nz == 1) 
				{
					rt.dcal(iGlob,2) = 0;
					rt.dcal(iGlob,3) = 0; 
					rt.dcal(iGlob,4) = 0;
					}
			});
	}
	    /* --------------------------------------------------
        Get RT matrix coefficients
    -------------------------------------------------- */
    inline void RTlinear_system(GwState &gw, GwDomain &gdom, RTMatrix &rtA, Parallel &par, RTState &rt ,std::vector<RTBC> &rtgbc)	{
        // Calculate matrix coefficients
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, ivg, iGlob, iGlobSW;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            iGlobSW = (hc+jj)*gdom.nxhc + ii + hc;



//添加判断条件，判断速度方向，如果速度方向为正，则取上游加权，否则取下游加权
//x方向速度方向判断
double Up_Weighting_x,Up_Weighting_y,Up_Weighting_z;
		if (rt.aveV(iGlob, 0) > 0)
		{
			Up_Weighting_x = rt.Up_Weighting_vplus;
		}
		else
		{
			Up_Weighting_x = rt.Up_Weighting_vminus;
		}

//y方向速度方向判断
		if (rt.aveV(iGlob, 1) > 0)
		{
			Up_Weighting_y = rt.Up_Weighting_vplus;
		}
		else
		{
			Up_Weighting_y = rt.Up_Weighting_vminus;
		}
//z方向速度方向判断
		if (rt.aveV(iGlob, 2) > 0)
		{
			Up_Weighting_z = rt.Up_Weighting_vplus;
		}
		else
		{
			Up_Weighting_z = rt.Up_Weighting_vminus;
		}

// 弥散项和对流项系数，不考虑含水率
		// rt.c_advxx(iGlob) = rt.aveV(iGlob,0) * rt.dt / gdom.dx;
		// rt.c_advyy(iGlob) = rt.aveV(iGlob,1) * rt.dt / gdom.dy;
		// // rt.c_advyy(iGlob) = 0;
		// rt.c_advzz(iGlob) = rt.aveV(iGlob,2) * rt.dt / gdom.dz(iGlob);
		// rt.c_difxx(iGlob) = rt.dcal(iGlob,0) * rt.dt / gdom.dx / gdom.dx;
		// rt.c_difyy(iGlob) = rt.dcal(iGlob,1) * rt.dt / gdom.dy / gdom.dy;
		// // rt.c_difyy(iGlob) = 0;
		// rt.c_difzz(iGlob) = rt.dcal(iGlob,2) * rt.dt / gdom.dz(iGlob) / gdom.dz(iGlob);

		// std::cout << "-------rt.c_advyy(iGlob)----- " <<rt.c_advyy(iGlob) << std::endl;
		// std::cout << "-------rt.c_difyy(iGlob)----- " <<rt.c_difyy(iGlob) << std::endl;

// 弥散项和对流项系数，考虑含水率
		rt.c_advxx(iGlob) = rt.aveV(iGlob,0) * rt.dt / gdom.dx;
		rt.c_advyy(iGlob) = rt.aveV(iGlob,1) * rt.dt / gdom.dy;		
		rt.c_advzz(iGlob) = rt.aveV(iGlob,2) * rt.dt / gdom.dz(iGlob);
		rt.c_difxx(iGlob) = gw.wc(iGlob,0)*rt.dcal(iGlob,0) * rt.dt / gdom.dx / gdom.dx;
		rt.c_difyy(iGlob) = gw.wc(iGlob,0)*rt.dcal(iGlob,1) * rt.dt / gdom.dy / gdom.dy;	
		rt.c_difzz(iGlob) = gw.wc(iGlob,0)*rt.dcal(iGlob,2) * rt.dt / gdom.dz(iGlob) / gdom.dz(iGlob);
		//交叉项弥散项系数
		rt.c_difxy(iGlob) = gw.wc(iGlob,0)*rt.dcal(iGlob,3) * rt.dt / gdom.dx / gdom.dy;
		rt.c_difxz(iGlob) = gw.wc(iGlob,0)*rt.dcal(iGlob,4) * rt.dt / gdom.dx / gdom.dz(iGlob);
		rt.c_difyz(iGlob) = gw.wc(iGlob,0)*rt.dcal(iGlob,5) * rt.dt / gdom.dy / gdom.dz(iGlob);
// A矩阵系数				
//           rt.RTcoef(idom, 0) = 1.0 + rt.c_advxx(iGlob) * (2 * Up_Weighting_x- 1) + rt.c_advyy(iGlob) * (2 * Up_Weighting_y- 1) + rt.c_advzz(iGlob) * (2 * Up_Weighting_z - 1) + 2.0 * (rt.c_difxx(iGlob)+ rt.c_difyy(iGlob) + rt.c_difzz(iGlob));
// 		rt.RTcoef(idom, 1) = - (rt.c_difxx(iGlob) - (1 - Up_Weighting_x) * rt.c_advxx(iGlob)) ;//Ci+1
// 		rt.RTcoef(idom, 2) = - (Up_Weighting_x * rt.c_advxx(iGlob) + rt.c_difxx(iGlob)); // # Ci-1
//           rt.RTcoef(idom, 3) =  - (rt.c_difyy(iGlob) - (1 - Up_Weighting_y) * rt.c_advyy(iGlob)) ; //Cj+1
//           rt.RTcoef(idom, 4) = - (Up_Weighting_y * rt.c_advyy(iGlob) + rt.c_difyy(iGlob)); //Cj-1	
//           rt.RTcoef(idom, 5) = - (rt.c_difzz(iGlob) - (1 - Up_Weighting_z) * rt.c_advzz(iGlob));// Ck+1 
// 		rt.RTcoef(idom, 6) = - (Up_Weighting_z * rt.c_advzz(iGlob) + rt.c_difzz(iGlob));//Ck-1		          
//           rt.RTcoef(idom, 7) = rt.c(iGlob, 1);
//A矩阵系数，考虑交叉项弥散系数				
          rt.RTcoef(idom, 0) = 1.0 + rt.c_advxx(iGlob) * (2 * Up_Weighting_x- 1) + rt.c_advyy(iGlob) * (2 * Up_Weighting_y- 1) + rt.c_advzz(iGlob) * (2 * Up_Weighting_z - 1) 
		+ 2.0 * (rt.c_difxx(iGlob)+ rt.c_difyy(iGlob) + rt.c_difzz(iGlob)+ rt.c_difxz(iGlob) + rt.c_difyz(iGlob)+ rt.c_difxy(iGlob));
		rt.RTcoef(idom, 1) = - (rt.c_difxx(iGlob) + rt.c_difxy(iGlob) + rt.c_difxz(iGlob) - (1 - Up_Weighting_x) * rt.c_advxx(iGlob)) ;//Ci+1
		rt.RTcoef(idom, 2) = - (Up_Weighting_x * rt.c_advxx(iGlob) + rt.c_difxx(iGlob)+ rt.c_difxy(iGlob) + rt.c_difxz(iGlob)); // # Ci-1
          rt.RTcoef(idom, 3) =  - (rt.c_difyy(iGlob)+ rt.c_difxy(iGlob) + rt.c_difyz(iGlob) - (1 - Up_Weighting_y) * rt.c_advyy(iGlob)) ; //Cj+1
          rt.RTcoef(idom, 4) = - (Up_Weighting_y * rt.c_advyy(iGlob) + rt.c_difyy(iGlob)+ rt.c_difxy(iGlob) + rt.c_difyz(iGlob)); //Cj-1	
          rt.RTcoef(idom, 5) = - (rt.c_difzz(iGlob)+ rt.c_difxz(iGlob) + rt.c_difyz(iGlob) - (1 - Up_Weighting_z) * rt.c_advzz(iGlob));// Ck+1 
		rt.RTcoef(idom, 6) = - (Up_Weighting_z * rt.c_advzz(iGlob) + rt.c_difzz(iGlob) + rt.c_difxz(iGlob) + rt.c_difyz(iGlob));//Ck-1		          
          
		rt.RTcoef(idom, 7) = rt.c(iGlob, 1);


//上游加权形式(原始)
		// rt.RTcoef(idom, 0)	= 1 + 2 *(rt.c_difxx(iGlob)+rt.c_difyy(iGlob)+rt.c_difzz(iGlob)) + (2*0.7-1)*rt.c_advxx(iGlob) + (2*0.7-1)*rt.c_advyy(iGlob) + (2*0.7-1)*rt.c_advzz(iGlob);
		// rt.RTcoef(idom, 1)	= -(rt.c_difxx(iGlob) - (1-0.7)*rt.c_advxx(iGlob));//ci+1
		// rt.RTcoef(idom, 2)	= -(rt.c_difxx(iGlob) + 0.7*rt.c_advxx(iGlob));//ci-1
		// // rt.RTcoef(idom, 3)	= -(rt.c_difyy(iGlob) - (1-Up_Weighting_y)*rt.c_advyy(iGlob));//cj+1
		// rt.RTcoef(idom, 3) = 0;
		// // rt.RTcoef(idom, 4)	= -(rt.c_difyy(iGlob) + Up_Weighting_y*rt.c_advyy(iGlob));//cj-1
		// rt.RTcoef(idom, 4) = 0;
		// rt.RTcoef(idom, 5)	= -(rt.c_difzz(iGlob) - (1-0.7)*rt.c_advzz(iGlob));//ck+1
		// rt.RTcoef(idom, 6)	= -(rt.c_difzz(iGlob) + 0.7*rt.c_advzz(iGlob));//ck-1
		// rt.RTcoef(idom, 7)	= rt.c(iGlob, 1);//右端项

//上游加权形式
		// rt.RTcoef(idom, 0) = 
		// 	gw.wc(iGlob,0) *1 + 2 * gw.wc(iGlob,0) * (rt.dcal(iGlob,0) * rt.dt / (gdom.dx * gdom.dx) + rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)))
		// +rt.aveV(iGlob, 0) * rt.dt *(2*Up_Weighting_x-1)/ gdom.dx + rt.aveV(iGlob, 1) * rt.dt *(2*Up_Weighting_y-1)/ gdom.dy + rt.aveV(iGlob, 2) * rt.dt *(2*Up_Weighting_z-1) / gdom.dz(iGlob);
		// rt.RTcoef(idom, 1) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ))-(rt.aveV(iGlob, 0) *(1-Up_Weighting_x)* rt.dt / gdom.dx) );
          // rt.RTcoef(idom, 2) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx )) + (rt.aveV(iGlob, 0) *Up_Weighting_x * rt.dt / gdom.dx)) ;
          // rt.RTcoef(idom, 3) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy))-(rt.aveV(iGlob, 1) * (1-Up_Weighting_y)* rt.dt / gdom.dy)) ;
          // rt.RTcoef(idom, 4) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy)) + (rt.aveV(iGlob, 1)  * Up_Weighting_y * rt.dt / gdom.dy));
          // rt.RTcoef(idom, 5) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)))-(rt.aveV(iGlob, 2) *(1-Up_Weighting_z)* rt.dt / gdom.dz(iGlob)) );
          // rt.RTcoef(idom, 6) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob))) + (rt.aveV(iGlob, 2) * Up_Weighting_z * rt.dt /  gdom.dz(iGlob)));
		// rt.RTcoef(idom, 7) = gw.wc(iGlob,0)*rt.c(iGlob, 1);
//上游加权形式,DXX,i+1/2,j,k形式,弥散项中心差分
		// rt.RTcoef(idom, 0) = 
		// 	gw.wc(iGlob,0) *1 
		// +rt.aveV(iGlob, 0) * rt.dt *(2*Up_Weighting_x-1)/ gdom.dx + rt.aveV(iGlob, 1) * rt.dt *(2*Up_Weighting_y-1)/ gdom.dy + rt.aveV(iGlob, 2) * rt.dt *(2*Up_Weighting_z-1) / gdom.dz(iGlob);
		// rt.RTcoef(idom, 1) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob,0) - rt.dcal(iGlob-1,0)) * rt.dt / 2*(gdom.dx  * gdom.dx ))-(rt.aveV(iGlob, 0) *(1-Up_Weighting_x)* rt.dt / gdom.dx) );
          // rt.RTcoef(idom, 2) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob-1,0) - rt.dcal(iGlob,0)) * rt.dt / 2*(gdom.dx  * gdom.dx )) + (rt.aveV(iGlob, 0) *Up_Weighting_x * rt.dt / gdom.dx)) ;
          // rt.RTcoef(idom, 3) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob,1) - rt.dcal(iGlob-gdom.nxhc,1)) * rt.dt / 2*(gdom.dy * gdom.dy))-(rt.aveV(iGlob, 1) * (1-Up_Weighting_y)* rt.dt / gdom.dy)) ;
          // rt.RTcoef(idom, 4) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob-gdom.nxhc,1) - rt.dcal(iGlob,1)) * rt.dt / 2*(gdom.dy * gdom.dy)) + (rt.aveV(iGlob, 1)  * Up_Weighting_y * rt.dt / gdom.dy));
          // rt.RTcoef(idom, 5) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob,2) - rt.dcal(iGlob-gdom.nxhc*gdom.nyhc,2)) * rt.dt / 2*(gdom.dz(iGlob) * gdom.dz(iGlob)))-(rt.aveV(iGlob, 2) *(1-Up_Weighting_z)* rt.dt / gdom.dz(iGlob)) );
          // rt.RTcoef(idom, 6) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob-gdom.nxhc*gdom.nyhc,2) - rt.dcal(iGlob,2)) * rt.dt / 2*(gdom.dz(iGlob) * gdom.dz(iGlob)))+(rt.aveV(iGlob, 2) * Up_Weighting_z * rt.dt /  gdom.dz(iGlob)));
		// rt.RTcoef(idom, 7) = gw.wc(iGlob,0)*rt.c(iGlob, 1);
//上游加权形式,DXX,i+1/2,j,k形式,弥散项中心差分,考虑Dxy,DXZ,DYZ
		// rt.RTcoef(idom, 0) = 
		// 	gw.wc(iGlob,0) *1 
		// +rt.aveV(iGlob, 0) * rt.dt *(2*Up_Weighting_x-1)/ gdom.dx + rt.aveV(iGlob, 1) * rt.dt *(2*Up_Weighting_y-1)/ gdom.dy + rt.aveV(iGlob, 2) * rt.dt *(2*Up_Weighting_z-1) / gdom.dz(iGlob);
		// rt.RTcoef(idom, 1) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob,0) - rt.dcal(iGlob-1,0)) * rt.dt / 2*(gdom.dx  * gdom.dx ))-(rt.aveV(iGlob, 0) *(1-Up_Weighting_x)* rt.dt / gdom.dx) );
          // rt.RTcoef(idom, 2) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob-1,0) - rt.dcal(iGlob,0)) * rt.dt / 2*(gdom.dx  * gdom.dx )) + (rt.aveV(iGlob, 0) *Up_Weighting_x * rt.dt / gdom.dx)) ;
          // rt.RTcoef(idom, 3) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob,1) - rt.dcal(iGlob-gdom.nxhc,1)) * rt.dt / 2*(gdom.dy * gdom.dy))-(rt.aveV(iGlob, 1) * (1-Up_Weighting_y)* rt.dt / gdom.dy)) ;
          // rt.RTcoef(idom, 4) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob-gdom.nxhc,1) - rt.dcal(iGlob,1)) * rt.dt / 2*(gdom.dy * gdom.dy)) + (rt.aveV(iGlob, 1)  * Up_Weighting_y * rt.dt / gdom.dy));
          // rt.RTcoef(idom, 5) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob,2) - rt.dcal(iGlob-gdom.nxhc*gdom.nyhc,2)) * rt.dt / 2*(gdom.dz(iGlob) * gdom.dz(iGlob)))-(rt.aveV(iGlob, 2) *(1-Up_Weighting_z)* rt.dt / gdom.dz(iGlob)) );         
		// rt.RTcoef(idom, 6) = -(gw.wc(iGlob,0)*((rt.dcal(iGlob-gdom.nxhc*gdom.nyhc,2) - rt.dcal(iGlob,2)) * rt.dt / 2*(gdom.dz(iGlob) * gdom.dz(iGlob)))+(rt.aveV(iGlob, 2) * Up_Weighting_z * rt.dt /  gdom.dz(iGlob)));
		// // Ci+1,j+1
		// rt.RTcoef(idom, 8) = -(gw.wc(iGlob,0) * rt.dcal(iGlob,3) * rt.dt / 2*(gdom.dx  * gdom.dy ));
		// //Ci+1,j-1
		// rt.RTcoef(idom, 9) = (gw.wc(iGlob,0) * rt.dcal(iGlob,3) * rt.dt / 2*(gdom.dx  * gdom.dy ));
		// //Ci-1,j+1
		// rt.RTcoef(idom, 10) = (gw.wc(iGlob,0) * rt.dcal(iGlob,3) * rt.dt / 2*(gdom.dx  * gdom.dy ));
		// //Ci-1,j-1
		// rt.RTcoef(idom, 11) = -(gw.wc(iGlob,0) * rt.dcal(iGlob,3) * rt.dt / 2*(gdom.dx  * gdom.dy ));
		// //Ci+1,k+1
		// rt.RTcoef(idom, 12) = -(gw.wc(iGlob,0) * rt.dcal(iGlob,4) * rt.dt / 2*(gdom.dx  * gdom.dz(iGlob) ));
		// //Ci+1,k-1
		// rt.RTcoef(idom, 13) = (gw.wc(iGlob,0) * rt.dcal(iGlob,4) * rt.dt / 2*(gdom.dx  * gdom.dz(iGlob) ));
		// //Ci-1,k+1
		// rt.RTcoef(idom, 14) = (gw.wc(iGlob,0) * rt.dcal(iGlob,4) * rt.dt / 2*(gdom.dx  * gdom.dz(iGlob) ));
		// //Ci-1,k-1
		// rt.RTcoef(idom, 15) = -(gw.wc(iGlob,0) * rt.dcal(iGlob,4) * rt.dt / 2*(gdom.dx  * gdom.dz(iGlob) ));
		// //Cj+1,k+1
		// rt.RTcoef(idom, 16) = -(gw.wc(iGlob,0) * rt.dcal(iGlob,5) * rt.dt / 2*(gdom.dy  * gdom.dz(iGlob) ));
		// //Cj+1,k-1
		// rt.RTcoef(idom, 17) = (gw.wc(iGlob,0) * rt.dcal(iGlob,5) * rt.dt / 2*(gdom.dy  * gdom.dz(iGlob) ));
		// //Cj-1,k+1
		// rt.RTcoef(idom, 18) = (gw.wc(iGlob,0) * rt.dcal(iGlob,5) * rt.dt / 2*(gdom.dy  * gdom.dz(iGlob) ));
		// //Cj-1,k-1
		// rt.RTcoef(idom, 19) = -(gw.wc(iGlob,0) * rt.dcal(iGlob,5) * rt.dt / 2*(gdom.dy  * gdom.dz(iGlob) ));		
		
		// rt.RTcoef(idom, 7) = gw.wc(iGlob,0)*rt.c(iGlob, 1);


//上游加权形式不考虑含水量
		// rt.RTcoef(idom, 0) = 
		// 	  1 + 2 *   (rt.dcal(iGlob,0) * rt.dt / (gdom.dx * gdom.dx) + rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)))
		// +rt.aveV(iGlob, 0) * rt.dt *(2*Up_Weighting_x-1)/ gdom.dx + rt.aveV(iGlob, 1) * rt.dt *(2*Up_Weighting_y-1)/ gdom.dy + rt.aveV(iGlob, 2) * rt.dt *(2*Up_Weighting_z-1) / gdom.dz(iGlob);
		// //c(i+1,j,k)
		// rt.RTcoef(idom, 1) = -( (rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ))-(rt.aveV(iGlob, 0) *(1-Up_Weighting_x)* rt.dt / gdom.dx) );
          // //c(i-1,j,k)
		// rt.RTcoef(idom, 2) = -( (rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx )) + (rt.aveV(iGlob, 0) *Up_Weighting_x * rt.dt / gdom.dx)) ;
          // //c(i,j+1,k)
	     // // rt.RTcoef(idom, 3) = -( (rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy))-(rt.aveV(iGlob, 1) * (1-Up_Weighting_y)* rt.dt / gdom.dy)) ;
          // rt.RTcoef(idom, 3) = 0;
		// //c(i,j-1,k)
		// // rt.RTcoef(idom, 4) = -( (rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy)) + (rt.aveV(iGlob, 1)  * Up_Weighting_y * rt.dt / gdom.dy));
          // rt.RTcoef(idom, 3) = 0;
		// //c(i,j,k+1)
		// rt.RTcoef(idom, 5) = -( (rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)))-(rt.aveV(iGlob, 2) *(1-Up_Weighting_z)* rt.dt / gdom.dz(iGlob)) );
          
		// //c(i,j,k-1)
		// rt.RTcoef(idom, 6) = -( (rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob))) + (rt.aveV(iGlob, 2) * Up_Weighting_z * rt.dt /  gdom.dz(iGlob)));
		
		// //c(i,j,k)
		// rt.RTcoef(idom, 7) =  rt.c(iGlob, 1);


//参照python代码		
		// double dx = 0.25;
		// double dy = 0.25;
		
		// double u = 0;
		// double v = 1e-5;
		
		// double D = 1e-6;
		// double c_advx = u * rt.dt / gdom.dx;
		// double c_advz = v * rt.dt / gdom.dy;
		// double c_difx = D * rt.dt / gdom.dx / gdom.dx;
		// double c_difz = D * rt.dt / gdom.dy / gdom.dy;		
		
		// rt.c_advxx(iGlob) = rt.aveV(iGlob,0) * rt.dt / gdom.dx;
		// rt.c_advyy(iGlob) = rt.aveV(iGlob,1) * rt.dt / gdom.dy;
		// rt.c_advzz(iGlob) = rt.aveV(iGlob,2) * rt.dt / gdom.dz(iGlob);
		// rt.c_difxx(iGlob) = rt.dcal(iGlob,0) * rt.dt / gdom.dx / gdom.dx;
		// rt.c_difyy(iGlob) = rt.dcal(iGlob,1) * rt.dt / gdom.dy / gdom.dy;
		// rt.c_difzz(iGlob) = rt.dcal(iGlob,2) * rt.dt / gdom.dz(iGlob) / gdom.dz(iGlob);

          // rt.RTcoef(idom, 2) = - (0.7 * rt.c_advxx(iGlob) + rt.c_difxx(iGlob)); // # Ci-1
          // rt.RTcoef(idom, 1) = - (rt.c_difxx(iGlob) - (1 - 0.7) * rt.c_advxx(iGlob)) ;//Ci+1 
          // rt.RTcoef(idom, 6) = - (0.7 * rt.c_advzz(iGlob) + rt.c_difzz(iGlob));//Cj-1
          // rt.RTcoef(idom, 5) = - (rt.c_difzz(iGlob) - (1 - 0.7) * rt.c_advzz(iGlob));// Cj+1 
		// rt.RTcoef(idom, 3) = 0;
		// rt.RTcoef(idom, 4) = 0;
          // rt.RTcoef(idom, 0) = 1.0 + rt.c_advxx(iGlob) * (2 * 0.7- 1) + rt.c_advzz(iGlob) * (2 * 0.7 - 1) + 2.0 * (rt.c_difxx(iGlob) + rt.c_difzz(iGlob));
          // rt.RTcoef(idom, 7) = rt.c(iGlob, 1);

		

//上游加权形式考虑含水量,包括对流项
		// rt.RTcoef(idom, 0) = 
		// 	gw.wc(iGlob,0) *1 + 2 * gw.wc(iGlob,0) * (rt.dcal(iGlob,0) * rt.dt / (gdom.dx * gdom.dx) + rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)))
		// +gw.wc(iGlob,0)*rt.aveV(iGlob, 0) * rt.dt *(2*Up_Weighting_x-1)/ gdom.dx + gw.wc(iGlob,0)*rt.aveV(iGlob, 1) * rt.dt *(2*Up_Weighting_y-1)/ gdom.dy + gw.wc(iGlob,0)*rt.aveV(iGlob, 2) * rt.dt *(2*Up_Weighting_z-1) / gdom.dz(iGlob);
		// rt.RTcoef(idom, 1) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ))-(gw.wc(iGlob,0)*rt.aveV(iGlob, 0) *(1-Up_Weighting_x)* rt.dt / gdom.dx) );
          // rt.RTcoef(idom, 2) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx )) + (gw.wc(iGlob,0)*rt.aveV(iGlob, 0) *Up_Weighting_x * rt.dt / gdom.dx)) ;
          // rt.RTcoef(idom, 3) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy))-(gw.wc(iGlob,0)*rt.aveV(iGlob, 1) * (1-Up_Weighting_y)* rt.dt / gdom.dy)) ;
          // rt.RTcoef(idom, 4) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy)) + (gw.wc(iGlob,0)*rt.aveV(iGlob, 1)  * Up_Weighting_y * rt.dt / gdom.dy));
          // rt.RTcoef(idom, 5) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)))-(gw.wc(iGlob,0)*rt.aveV(iGlob, 2) *(1-Up_Weighting_z)* rt.dt / gdom.dz(iGlob)) );
          // rt.RTcoef(idom, 6) = -(gw.wc(iGlob,0)*(rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob))) + (gw.wc(iGlob,0)*rt.aveV(iGlob, 2) * Up_Weighting_z * rt.dt /  gdom.dz(iGlob)));
		// rt.RTcoef(idom, 7) = gw.wc(iGlob,0)*rt.c(iGlob, 1);		
		// // std::cout << "-------* rt.Up_Weighting----- " <<rt.Up_Weighting << std::endl;

		if (rt.rt_scheme != 1)	//picard迭代
		{
		
			rt.RTcoef(idom,7) -= (gw.wc(iGlob,1) - gw.wc(iGlob,0));//为了消除迭代中的线性误差
			// rt.RTcoef(idom,7) = rt.RTcoef(idom,7);
		}
		
	}   );
		 /*--------------后续完善部分*/
		  // Apply internal boundary conditions (needed when MPI is used)

		// Apply outer boundary conditions
			for (int k = 0; k < rtgbc.size(); k++) {
				rtgbc[k].applyRTMatBC(rt, gw, gdom, par);
			}
			
    
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
//弥散项中心差分
			// //Ci+1,j+1
			// if (ii < gdom.nx-1 && jj< gdom.ny-1){
			// 	rtA.rt_ind(irow) = idom + gdom.nx + 1;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,8);  
			// 	irow++;
			// }
			// //Ci+1,j-1
			// if (ii < gdom.nx-1 && jj > 0){
			// 	rtA.rt_ind(irow) = idom + 1 - gdom.nx;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,9);  
			// 	irow++;				
			// }
			// //ci-1,j+1
			// if (ii > 0 && jj < gdom.ny - 1){
			// 	rtA.rt_ind(irow) = idom -1+gdom.nx;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,10);  
			// 	irow++;				
			// }
			// //ci-1,j-1
			// if (ii > 0 && jj > 0){
			// 	rtA.rt_ind(irow) = idom - gdom.nx - 1;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,11);  
			// 	irow++;				
			// }
			// //ci+1,k+1
			// if (ii < gdom.nx-1 && kk < gdom.nz-1){
			// 	rtA.rt_ind(irow) = idom + 1 + gdom.nx * gdom.ny;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,12);  
			// 	irow++;				
			// }
			// //ci+1,k-1
			// if (ii < gdom.nx-1 && kk > 0){
			// 	rtA.rt_ind(irow) = idom + 1 -gdom.nx*gdom.ny;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,13);  
			// 	irow++;				
			// }
			// //ci-1,k+1
			// if (ii > 0 && kk < gdom.nz-1){
			// 	rtA.rt_ind(irow) = idom -1 + gdom.nx*gdom.ny;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,14);  
			// 	irow++;				
			// }
			// //ci-1,k-1
			// if (ii > 0 && kk > 0){
			// 	rtA.rt_ind(irow) = idom -1 - gdom.nx*gdom.ny;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,15);  
			// 	irow++;				
			// }
			// //cj+1,k+1
			// if (jj < gdom.ny - 1 && kk < gdom.nz - 1){
			// 	rtA.rt_ind(irow) = idom + 1 + gdom.nx*gdom.ny;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,16);  
			// 	irow++;				
			// }
			// //cj+1,k-1
			// if (jj < gdom.ny -1 && kk > 0){
			// 	rtA.rt_ind(irow) = idom + gdom.nx - gdom.nx*gdom.ny;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,17);  
			// 	irow++;				
			// }
			// //cj-1,k+1
			// if (jj > 0 && kk < gdom.nz -1){
			// 	rtA.rt_ind(irow) = idom -gdom.nx + gdom.nx*gdom.ny;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,18);  
			// 	irow++;				
			// }
			// //cj-1,k-1
			// if (jj > 0 && kk > 0){
			// 	rtA.rt_ind(irow) = idom -gdom.nx - gdom.nx*gdom.ny;
			// 	rtA.rt_val(irow) = rt.RTcoef(idom,19);  
			// 	irow++;				
			// }
//弥散项中心差分

			// //添加Ci+1,j+1等
			// if (ii < gdom.nx-1 && jj < gdom.ny-1)	{rtA.rt_ind(irow) = idom + gdom.nx + 1;		        rtA.rt_val(irow) = rt.RTcoef(idom,8);  irow++;}
			// //Ci+1,k+1
			// if (ii < gdom.nx-1 && kk < gdom.nz-1)	{rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny +1;	        rtA.rt_val(irow) = rt.RTcoef(idom,9);  irow++;}
			// //Cj+1,k+1
			// if (jj < gdom.ny-1 && kk < gdom.nz-1)	{rtA.rt_ind(irow) = idom + gdom.nx + gdom.nx*gdom.ny;	rtA.rt_val(irow) = rt.RTcoef(idom,10);  irow++;}
			rtA.rt_rhs(idom) = rt.RTcoef(idom,7);
		   });

			

	}

/*----------------20240426修改-----------*///非迭代方法时间步长控制
    inline void dt_con(RTState &rt, GwState &gw, GwDomain &gdom)	{
    	real dc_max, dt_old;
    	dt_old = rt.dt;//在serghei.h中初始化rt.dt=gdom.dt_init;
        Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA (int idx, real &tmp) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idx, kk, jj, ii);
            // gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            real dc = fabs(rt.c(iGlob,1) - rt.c(iGlob,0));
			tmp = (dc > tmp) ? dc : tmp;
		} , Kokkos::Max<real>(dc_max) );
    	if (dc_max > 0.02)	{rt.dt = rt.dt * 0.9;}
    	else if (dc_max >= 0.0 & dc_max < 0.01)	{rt.dt = rt.dt * 1.1;}
    	if (rt.dt > gdom.dt_max)	{rt.dt = gdom.dt_max;}
    	else if (rt.dt < gdom.dt_init)	{rt.dt = gdom.dt_init;}
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
