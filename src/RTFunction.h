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

// Type of adsorption reaction
#define Equilibrium_Linear_Model 1
#define Equilibrium_Freundlich_Model 2
#define Equilibrium_Langmuir_Model 3
#define Nonequilibrium_Model 4


class RTFunction
{

private:
	Kokkos::Timer timer;

public:
	/* --------------------------------------------------
		rt_Gauss_Seidel_solve
	-------------------------------------------------- */
	template <typename execution_space, typename type_solver>
	inline void rt_Gauss_Seidel_solve(RTState &rt, RTMatrix &rtA, GwState &gw, GwDomain &gdom, std::vector<RTBC> &rtgbc, GwMPI &gmpi, Parallel &par, type_solver &rtsolver)
	{
		int iter, ierr = 1;
		real dt_tmp;

		// 施加边界条件
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
		}
		// 扩散系数计算
		dispersion_tensor(rt, gw, gdom);

		RTlinear_system(gw, gdom, rtA, par, rt, rtgbc);

		timer.reset();

		//选择求解方法
		if (rt.rt_scheme == 3){
				rtsolver.Eigen(rtA);
				}
		else if (rt.rt_scheme == 4){
				rtsolver.Gauss_Seidel(rtA);
				}		
		gdom.timers.rtlinsol += timer.seconds();

		// Update concentration
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            rt.c(iGlob,1) = rtA.rt_x(idom); //0为n+1时刻，1为n时刻
			
			//20240814添加,固相浓度更新
			//2:Equilibrium-Nonlinear sorption,Freundlich;
			if (rt.AdsorptionDesorptionModel == Equilibrium_Freundlich_Model){			
				// rt.c_solid(iGlob,0) = (rt.beta * rtA.rt_x(idom) + rt.rho_b * rt.c_solid(iGlob, 1) / rt.dt ) / (rt.rho_b / rt.dt + rt.beta / rt.Kd + rt.lambda_2 * rt.rho_b) ;
				rt.c_solid(iGlob,0) = rt.Kf * pow(rtA.rt_x(idom), rt.Nf);
				rt.c_solid(iGlob,1) = rt.c_solid(iGlob,0);
			}
			if (rt.AdsorptionDesorptionModel == Nonequilibrium_Model){			
				rt.c_solid(iGlob,0) = (rt.beta * rtA.rt_x(idom) + rt.rho_b * rt.c_solid(iGlob, 1) / rt.dt ) / (rt.rho_b / rt.dt + rt.beta / rt.Kd + rt.lambda_2 * rt.rho_b) ;
				rt.c_solid(iGlob,1) = rt.c_solid(iGlob,0);
			}			
			
			});						
		gmpi.mpi_sendrecv(rt.c, gdom, par);
		// 施加边界条件
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
		}
		gmpi.mpi_sendrecv(rt.c, gdom, par);
		// 施加边界条件
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
		}

		dt_con(rt, gw, gdom); // 计算浓度残差值，并更新
		// Update time step
		dt_tmp = gdom.dt;

		ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { rt.c(iGlob, 0) = rt.c(iGlob, 1); });
	}
	/* --------------------------------------------------
		Top-level PCA solver
	-------------------------------------------------- */
	template <typename execution_space, typename type_solver>
	inline void rt_pca_solve(RTState &rt, RTMatrix &rtA, GwState &gw, GwDomain &gdom, std::vector<RTBC> &rtgbc, GwMPI &gmpi, Parallel &par, type_solver &rtsolver)
	{
		int iter, ierr = 1;
		real dt_tmp;

		// 施加边界条件
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
		}
		// 扩散系数计算
		dispersion_tensor(rt, gw, gdom);

		RTlinear_system(gw, gdom, rtA, par, rt, rtgbc);

		timer.reset();
#if SERGHEI_KOKKOSKERNELS_SOLVER
		rtsolver.kkpcg(rtA);
#else
		iter = rtsolver.cg(rtA, gdom);
#endif
		gdom.timers.rtlinsol += timer.seconds();

		// Update concentration
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            rt.c(iGlob,1) = rtA.rt_x(idom); });
		gmpi.mpi_sendrecv(rt.c, gdom, par);
		// 施加边界条件
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
		}
		gmpi.mpi_sendrecv(rt.c, gdom, par);
		// 施加边界条件
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
		}

		dt_con(rt, gw, gdom); // 计算浓度残差值，并更新
		// Update time step
		dt_tmp = gdom.dt;

		ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { rt.c(iGlob, 0) = rt.c(iGlob, 1); });
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
		// 	std::cout << "对应的索引: "

		
	}

	/* --------------------------------------------------
		Top-level Picard solver
	-------------------------------------------------- */
	template <typename execution_space, typename type_solver>
	inline void rt_picard_solve(RTState &rt, RTMatrix &rtA, GwState &gw, GwDomain &gdom, std::vector<RTBC> &rtgbc,
								type_solver &rtsolver, GwMPI &gmpi, GwIntegrator &gint, Parallel &par)
	{

		int iter, iter_cg, iter_max = 100, ierr = 1;
		real eps_diff = 1.0, eps_tmp, eps_old = 1.0, eps = 1.0, dt_tmp;
		real eps_min = 5e-5; // eps_min = 5e-6,20240510修改了收敛残差标准
		// 施加边界条件
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
		}
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc; });
		iter = 0;
		// 计算扩散系数
		dispersion_tensor(rt, gw, gdom);

		while (iter < iter_max && eps_diff / eps_old > eps_min && eps > eps_min)
		{

			RTlinear_system(gw, gdom, rtA, par, rt, rtgbc);

			timer.reset();
#if SERGHEI_KOKKOSKERNELS_SOLVER
			rtsolver.kkpcg(rtA);
#else
			iter = rtsolver.cg(rtA, gdom);
#endif
			gdom.timers.rtlinsol += timer.seconds();
			Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                rt.c(iGlob,0) = rt.c(iGlob,1);
                rt.c(iGlob,1) = rtA.rt_x(idom); });
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
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { rt.c(iGlob, 0) = rt.c(iGlob, 1); });
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
	std::pair<real, int> find_max_value(const RTState &rt)
	{
		real max_value = 0;
		int max_index = 0;
		for (int i = 1; i < rt.c.extent(0); ++i)
		{
			real current_value = rt.c(i, 1);
			if (current_value > max_value)
			{
				max_value = current_value;
				max_index = i;
			}
		}
		return std::make_pair(max_value, max_index);
	}

	// 计算Peclet
	std::pair<real, int> find_peclet(const RTState &rt, GwDomain &gdom)
	{
		real max_value = 0;
		int max_index = 0;
		for (int i = 1; i < rt.c.extent(0); ++i)
		{
			real current_value = fabs(gdom.dx * rt.aveV(i, 3) / rt.dcal(i, 0));
			if (current_value > max_value)
			{
				max_value = current_value;
				max_index = i;
			}
		}
		return std::make_pair(max_value, max_index);
	}

	// 计算最大孔隙均速
	std::pair<real, int> find_aveV(const RTState &rt)
	{
		real max_value = 0;
		int max_index = 0;
		for (int i = 1; i < rt.c.extent(0); ++i)
		{
			real current_value = fabs(rt.aveV(i, 3));
			if (current_value > max_value)
			{
				max_value = current_value;
				max_index = i;
			}
		}
		return std::make_pair(max_value, max_index);
	}


// 计算n+1时刻q值

	// Calculation of the hydrodynamic dispersion coefficient tensor
	inline void dispersion_tensor(RTState &rt, GwState &gw, GwDomain &gdom)
	{
		Kokkos::parallel_for(
			gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
				// rt.c(iGlob, 0) 表示下一时间步长的浓度值
				// rt.c(iGlob, 1) 表示当前步长的浓度值
				// aveVB为边界流速
				// rt.aveVB(iGlob, 0) = -(gw.q(iGlob, 0));
				// rt.aveVB(iGlob, 1) = (gw.q(iGlob, 1));
				// rt.aveVB(iGlob, 2) = -(gw.q(iGlob, 2));
				// rt.aveVB(iGlob, 3) = sqrt(pow(rt.aveVB(iGlob, 0), 2) + pow(rt.aveVB(iGlob, 1), 2) + pow(rt.aveVB(iGlob, 2), 2));

				//!gw.q_new的计算中采用的水头h时间为0时刻的水头（n+1)时刻的水头
				rt.aveVB(iGlob, 0) = -(gw.q_new(iGlob, 0));
				rt.aveVB(iGlob, 1) = (gw.q_new(iGlob, 1));
				rt.aveVB(iGlob, 2) = -(gw.q_new(iGlob, 2));	
				// rt.aveVB(iGlob, 0) = -(gw.q(iGlob, 0));
				// rt.aveVB(iGlob, 1) = (gw.q(iGlob, 1));
				// rt.aveVB(iGlob, 2) = -(gw.q(iGlob, 2));		
				rt.aveVB(iGlob, 3) = sqrt(pow(rt.aveVB(iGlob, 0), 2) + pow(rt.aveVB(iGlob, 1), 2) + pow(rt.aveVB(iGlob, 2), 2));

				// aveV为节点流速
				rt.aveV(iGlob, 0) = (rt.aveVB(iGlob, 0) + rt.aveVB(iGlob - 1, 0)) / 2;
				rt.aveV(iGlob, 1) = (rt.aveVB(iGlob, 1) + rt.aveVB(iGlob - gdom.nxhc, 1)) / 2;
				rt.aveV(iGlob, 2) = (rt.aveVB(iGlob, 2) + rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2)) / 2;
				rt.aveV(iGlob, 3) = sqrt(pow(rt.aveV(iGlob, 0), 2) + pow(rt.aveV(iGlob, 1), 2) + pow(rt.aveV(iGlob, 2), 2));
				
				
// 机械弥散系数计算
				// // Dxx
				// rt.dcal(iGlob, 0) = (rt.alpha_T * (pow(rt.aveVB(iGlob, 1), 2) + pow(rt.aveVB(iGlob, 2), 2)) + rt.alpha_L * pow(rt.aveVB(iGlob, 0), 2)) / rt.aveVB(iGlob, 3);

				// // DYY
				// rt.dcal(iGlob, 1) = (rt.alpha_T * (pow(rt.aveVB(iGlob, 0), 2) + pow(rt.aveVB(iGlob, 2), 2)) + rt.alpha_L * pow(rt.aveVB(iGlob, 1), 2)) / rt.aveVB(iGlob, 3);
				// // DZZ
				// rt.dcal(iGlob, 2) = (rt.alpha_T * (pow(rt.aveVB(iGlob, 0), 2) + pow(rt.aveVB(iGlob, 1), 2)) + rt.alpha_L * pow(rt.aveVB(iGlob, 2), 2)) / rt.aveVB(iGlob, 3);

				// // Dxy
				// rt.dcal(iGlob, 3) = (rt.alpha_L - rt.alpha_T) * rt.aveVB(iGlob, 0) * rt.aveVB(iGlob, 1) / rt.aveVB(iGlob, 3);
				// // Dxz
				// rt.dcal(iGlob, 4) = (rt.alpha_L - rt.alpha_T) * rt.aveVB(iGlob, 0) * rt.aveVB(iGlob, 2) / rt.aveVB(iGlob, 3);
				// // Dyz
				// rt.dcal(iGlob, 5) = (rt.alpha_L - rt.alpha_T) * rt.aveVB(iGlob, 1) * rt.aveVB(iGlob, 2) / rt.aveVB(iGlob, 3);
//机械弥散系数计算20240921
//! Dxx
				//qx,i+1/2,j,k
				rt.q_dispersion(iGlob, 0) = rt.aveVB(iGlob, 0);
				//qy,i+1/2,j,k
				rt.q_dispersion(iGlob, 1) = 0.5*((rt.aveVB(iGlob-gdom.nxhc,1)+rt.aveVB(iGlob,1))*0.5+(rt.aveVB(iGlob+1-gdom.nxhc,1)+rt.aveVB(iGlob+1,1))*0.5);
				//!qz,i+1/2,j,k
				rt.q_dispersion(iGlob, 2) = 0.5*((rt.aveVB(iGlob-gdom.nxhc*gdom.nyhc,2)+rt.aveVB(iGlob,2))*0.5+(rt.aveVB(iGlob+1-gdom.nxhc*gdom.nyhc,2)+rt.aveVB(iGlob+1,2))*0.5);
											

				//q,i+1/2,j,k
				rt.q_dispersion(iGlob, 3) = sqrt(pow(rt.q_dispersion(iGlob, 0), 2) + pow(rt.q_dispersion(iGlob, 1), 2) + pow(rt.q_dispersion(iGlob, 2), 2));
				// Dxxi+1/2,j,k
				rt.dcal(iGlob, 0) = (rt.alpha_T * (pow(rt.q_dispersion(iGlob, 1), 2) + pow(rt.q_dispersion(iGlob, 2), 2)) + rt.alpha_L * pow(rt.q_dispersion(iGlob, 0), 2)) / rt.q_dispersion(iGlob, 3);
				
				//Dxyi+1/2,j,k
				rt.dcal(iGlob, 3) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 0) * rt.q_dispersion(iGlob, 1) / rt.q_dispersion(iGlob, 3);
				//!Dxzi+1/2,j,k
				rt.dcal(iGlob, 4) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 0) * rt.q_dispersion(iGlob, 2) / rt.q_dispersion(iGlob, 3);

//! DYY
				//qx,i,j+1/2,k
				rt.q_dispersion(iGlob, 4) = 0.5*((rt.aveVB(iGlob,1)+rt.aveVB(iGlob-1,1))*0.5+(rt.aveVB(iGlob+gdom.nxhc,1)+rt.aveVB(iGlob+gdom.nxhc-1,1))*0.5);
				//qy,i,j+1/2,k
				rt.q_dispersion(iGlob, 5) = rt.aveVB(iGlob, 1);
				//qz,i,j+1/2,k
				rt.q_dispersion(iGlob, 6) = 0.5*((rt.aveVB(iGlob-gdom.nxhc*gdom.nyhc,2)+rt.aveVB(iGlob,2))*0.5+(rt.aveVB(iGlob+gdom.nxhc-gdom.nxhc*gdom.nyhc,2)+rt.aveVB(iGlob+gdom.nxhc,2))*0.5);
				//q,i,j+1/2,k
				rt.q_dispersion(iGlob, 7) = sqrt(pow(rt.q_dispersion(iGlob, 4), 2) + pow(rt.q_dispersion(iGlob, 5), 2) + pow(rt.q_dispersion(iGlob, 6), 2));
				// DYYi,j+1/2,k
				rt.dcal(iGlob, 1) = (rt.alpha_T * (pow(rt.q_dispersion(iGlob, 4), 2) + pow(rt.q_dispersion(iGlob, 6), 2)) + rt.alpha_L * pow(rt.q_dispersion(iGlob, 5), 2)) / rt.q_dispersion(iGlob, 7);
				//Dyx,i,j+1/2,k
				rt.dcal(iGlob, 5) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 5) * rt.q_dispersion(iGlob, 4) / rt.q_dispersion(iGlob, 7);
				//Dyzi,j+1/2,k
				rt.dcal(iGlob, 6) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 5) * rt.q_dispersion(iGlob, 6) / rt.q_dispersion(iGlob, 7);
//! DZZ	
				//qx,i,j,k+1/2
				rt.q_dispersion(iGlob, 8) = 0.5*((rt.aveVB(iGlob,0)+rt.aveVB(iGlob+gdom.nxhc*gdom.nyhc,0))*0.5+(rt.aveVB(iGlob-1,0)+rt.aveVB(iGlob+gdom.nxhc*gdom.nyhc-1,0))*0.5);
				//qy,i,j,k+1/2
				rt.q_dispersion(iGlob, 9) = 0.5*((rt.aveVB(iGlob,1)+rt.aveVB(iGlob-gdom.nxhc,1))*0.5+(rt.aveVB(iGlob-gdom.nxhc+gdom.nxhc*gdom.nyhc,1)+rt.aveVB(iGlob+gdom.nxhc*gdom.nyhc,1))*0.5);
				//qz,i,j,k+1/2
				rt.q_dispersion(iGlob, 10) = rt.aveVB(iGlob, 2);
				//q,i,j,k+1/2
				rt.q_dispersion(iGlob, 11) = sqrt(pow(rt.q_dispersion(iGlob, 8), 2) + pow(rt.q_dispersion(iGlob, 9), 2) + pow(rt.q_dispersion(iGlob, 10), 2));
				// Dzz,i,j,k+1/2											
				rt.dcal(iGlob, 2) = (rt.alpha_T * (pow(rt.q_dispersion(iGlob, 8), 2) + pow(rt.q_dispersion(iGlob, 9), 2)) + rt.alpha_L * pow(rt.q_dispersion(iGlob, 10), 2)) / rt.q_dispersion(iGlob, 11);
				//!Dzx,i,j,k+1/2
				rt.dcal(iGlob, 7) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 10) * rt.q_dispersion(iGlob, 8) / rt.q_dispersion(iGlob, 11);
				// rt.dcal(iGlob, 7) = rt.dcal(iGlob, 4);
				// rt.dcal(iGlob, 4) = -rt.dcal(iGlob, 7);
				
				// rt.dcal(iGlob, 7) = 0;
				// rt.dcal(iGlob, 4) = 0;

				//Dzy,i,j,k+1/2
				rt.dcal(iGlob, 8) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 10) * rt.q_dispersion(iGlob, 9) / rt.q_dispersion(iGlob, 11);
				
				//!如果网格qz速度为0，则Dzz=0，即当底边界为水动力无通量边界时
				//todo 后续考虑是不是要结合边界条件对弥散系数张量进行修改
				
				if (rt.q_dispersion(iGlob, 10)== 0)//qz,i,j,k+1/2
				{
					rt.dcal(iGlob, 2) = 0;									
				}

//!调试	
				// if (iGlob ==1320)
				// {
				
				// std::cout<<"222rt.q_dispersion(iGlob, 8)\n"<<rt.q_dispersion(iGlob, 8)<<std::endl;	
				// // // std::cout<<"333rt.q_dispersion(iGlob, 9)\n"<<rt.q_dispersion(iGlob, 9)<<std::endl;
				// std::cout<<"222rt.q_dispersion(iGlob, 10)\n"<<rt.q_dispersion(iGlob, 10)<<std::endl;
				// std::cout<<"222rt.q_dispersion(iGlob, 11)\n"<<rt.q_dispersion(iGlob, 11)<<std::endl;
				// std::cout<<"222rt.dcal(iGlob, 7)\n"<<rt.dcal(iGlob, 7)<<std::endl;
				// std::cout<<"2222--------------------------------------------\n"<<std::endl;

				// std::cout<<"333rt.dcal(iGlob, 4)\n"<<rt.dcal(iGlob, 4)<<std::endl;
				// std::cout<<"333rt.q_dispersion(iGlob, 0)\n"<<rt.q_dispersion(iGlob, 0)<<std::endl;
				// std::cout<<"333rt.q_dispersion(iGlob, 2)\n"<<rt.q_dispersion(iGlob, 2)<<std::endl;
				// std::cout<<"333rt.q_dispersion(iGlob, 3)\n"<<rt.q_dispersion(iGlob, 3)<<std::endl;
				// // // std::cout<<"2222--------------------------------------------\n"<<std::endl;				
				// std::cout<<"333rt.aveVB(iGlob, 0)\n"<<rt.aveVB(iGlob, 0)<<std::endl;
				// std::cout<<"333rt.aveVB(iGlob, 2)\n"<<rt.aveVB(iGlob, 2)<<std::endl;
				// std::cout<<"333gw.q_new(iGlob, 2)\n"<<gw.q_new(iGlob, 2)<<std::endl;
				// std::cout<<"333gw.q(iGlob, 2)\n"<<gw.q(iGlob, 2)<<std::endl;
				// }
//!调试					



				// 将rt.d中nan值替换
				for (int i = 0; i < 9; ++i)
				{
					if (std::isnan(rt.dcal(iGlob, i)))
					{
						rt.dcal(iGlob, i) = 0;
					}
				};
				// std::cout<<"Dxxi+1/2,j,k Dxxi+1/2,j,k rt.dcal(iGlob, 0)\n"<<rt.dcal(iGlob, 0)<<std::endl;
				// std::cout<<"DZZ,i,j,k+1/2 rt.dcal(iGlob, 2)\n"<<rt.dcal(iGlob, 2)<<std::endl;
				/*------弥散系数=机械弥散+分子弥散----------*/

				rt.dcal(iGlob, 0) += gw.wc(iGlob, 0) * rt.diffusion_molecular;
				rt.dcal(iGlob, 1) += gw.wc(iGlob, 0) * rt.diffusion_molecular;
				rt.dcal(iGlob, 2) += gw.wc(iGlob, 0) * rt.diffusion_molecular;


				// 根据研究区域维数对水动力弥散系数进行调整
				if (gdom.nx == 1)//yz纬度，qx=0，Dxx=0,Dxy=DYX=DXZ=DZX=0
				{
					rt.dcal(iGlob, 0) = 0;
					rt.dcal(iGlob, 3) = 0;
					rt.dcal(iGlob, 4) = 0;
					rt.dcal(iGlob, 5) = 0;
					rt.dcal(iGlob, 7) = 0;
				}
				else if (gdom.ny == 1)
				{
					rt.dcal(iGlob, 1) = 0;
					rt.dcal(iGlob, 3) = 0;
					rt.dcal(iGlob, 5) = 0;
					rt.dcal(iGlob, 6) = 0;
					rt.dcal(iGlob, 8) = 0;
				}
				else if (gdom.nz == 1)
				{
					rt.dcal(iGlob, 2) = 0;
					rt.dcal(iGlob, 4) = 0;
					rt.dcal(iGlob, 6) = 0;
					rt.dcal(iGlob, 7) = 0;
					rt.dcal(iGlob, 8) = 0;
				}



			});
	}
	/* --------------------------------------------------
	Get RT matrix coefficients
-------------------------------------------------- */
	inline void RTlinear_system(GwState &gw, GwDomain &gdom, RTMatrix &rtA, Parallel &par, RTState &rt, std::vector<RTBC> &rtgbc)
	{
		// Calculate matrix coefficients
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, ivg, iGlob, iGlobSW;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
			iGlobSW = (hc + jj) * gdom.nxhc + ii + hc;

			// 添加判断条件，判断速度方向，如果速度方向为正，则取上游加权，否则取下游加权					
			//对流项加权系数
			// double Up_Weighting_xp = (rt.aveVB(iGlob, 0) > 0) ? rt.Up_Weighting_vplus : 
			// 						(rt.aveVB(iGlob, 0) < 0) ? rt.Up_Weighting_vminus : 0.5;
			// double Up_Weighting_xm = (rt.aveVB(iGlob - 1, 0) > 0) ? rt.Up_Weighting_vplus : 
			// 						(rt.aveVB(iGlob - 1, 0) < 0) ? rt.Up_Weighting_vminus : 0.5;
			// double Up_Weighting_yp = (rt.aveVB(iGlob, 1) > 0) ? rt.Up_Weighting_vplus : 
			// 						(rt.aveVB(iGlob, 1) < 0) ? rt.Up_Weighting_vminus : 0.5;
			// double Up_Weighting_ym = (rt.aveVB(iGlob - gdom.nxhc, 1) > 0) ? rt.Up_Weighting_vplus : 
			// 						(rt.aveVB(iGlob - gdom.nxhc, 1) < 0) ? rt.Up_Weighting_vminus : 0.5;
			// double Up_Weighting_zp = (rt.aveVB(iGlob, 2) > 0) ? rt.Up_Weighting_vplus : 
			// 						(rt.aveVB(iGlob, 2) < 0) ? rt.Up_Weighting_vminus : 0.5;
			// double Up_Weighting_zm = (rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) > 0) ? rt.Up_Weighting_vplus : 
			// 						(rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) < 0) ? rt.Up_Weighting_vminus : 0.5;

			rt.Up_Weighting_xp(iGlob) = (rt.aveVB(iGlob, 0) > 0) ? rt.Up_Weighting_vplus : 
									(rt.aveVB(iGlob, 0) < 0) ? rt.Up_Weighting_vminus : 0.5;
			rt.Up_Weighting_xm(iGlob) = (rt.aveVB(iGlob - 1, 0) > 0) ? rt.Up_Weighting_vplus : 
									(rt.aveVB(iGlob - 1, 0) < 0) ? rt.Up_Weighting_vminus : 0.5;
			rt.Up_Weighting_yp(iGlob) = (rt.aveVB(iGlob, 1) > 0) ? rt.Up_Weighting_vplus : 
									(rt.aveVB(iGlob, 1) < 0) ? rt.Up_Weighting_vminus : 0.5;
			rt.Up_Weighting_ym(iGlob) = (rt.aveVB(iGlob - gdom.nxhc, 1) > 0) ? rt.Up_Weighting_vplus : 
									(rt.aveVB(iGlob - gdom.nxhc, 1) < 0) ? rt.Up_Weighting_vminus : 0.5;
			rt.Up_Weighting_zp(iGlob) = (rt.aveVB(iGlob, 2) > 0) ? rt.Up_Weighting_vplus : 
									(rt.aveVB(iGlob, 2) < 0) ? rt.Up_Weighting_vminus : 0.5;
			rt.Up_Weighting_zm(iGlob) = (rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) > 0) ? rt.Up_Weighting_vplus : 
									(rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) < 0) ? rt.Up_Weighting_vminus : 0.5;
			
			//弥散项加权系数
			// rt.wz(iGlob) = gdom.dz(iGlob + gdom.nxhc*gdom.nyhc) 
			// 				/ (gdom.dz(iGlob) + gdom.dz(iGlob + gdom.nxhc*gdom.nyhc));
			


//!20240821修改
			// 弥散项和对流项系数
			rt.c_advxx(iGlob) = rt.aveVB(iGlob, 0) * rt.dt / gdom.dx;
			rt.c_advyy(iGlob) = rt.aveVB(iGlob, 1) * rt.dt / gdom.dy;
			rt.c_advzz(iGlob) = rt.aveVB(iGlob, 2) * rt.dt ;
			rt.c_difxx(iGlob) = rt.dcal(iGlob, 0) * rt.dt / gdom.dx / gdom.dx;
			// rt.c_difxx(iGlob) = 0;
			rt.c_difyy(iGlob) = rt.dcal(iGlob, 1) * rt.dt / gdom.dy / gdom.dy;
			// rt.c_difyy(iGlob) = 0;
			rt.c_difzz(iGlob) = rt.dcal(iGlob, 2) * rt.dt / (0.5*gdom.dz(iGlob) + 0.5*gdom.dz(iGlob + gdom.nxhc * gdom.nyhc));
			// rt.c_difzz(iGlob) = 0;
			rt.c_difxy(iGlob) = rt.dcal(iGlob, 3) * rt.dt / (2 * gdom.dx * gdom.dy);
			// rt.c_difxy(iGlob) = 0;
			rt.c_difxz(iGlob) = rt.dcal(iGlob, 4) * rt.dt ;
			// rt.c_difxz(iGlob) = 0;
			rt.c_difyx(iGlob) = rt.dcal(iGlob, 5) * rt.dt / (2 * gdom.dx * gdom.dy);
			// rt.c_difyx(iGlob) = 0;
			rt.c_difyz(iGlob) = rt.dcal(iGlob, 6) * rt.dt ;		
			// rt.c_difyz(iGlob) = 0;	
			rt.c_difzx(iGlob) = rt.dcal(iGlob, 7) * rt.dt ;
			// rt.c_difzx(iGlob) = 0;
			rt.c_difzy(iGlob) = rt.dcal(iGlob, 8) * rt.dt ;
			// rt.c_difzy(iGlob) = 0;
			

//todo
//*
//!-----------	输出结果显示，gw.wc(iGlob, 0)-gw.wc(iGlob, 1)为0-----------------//		
//! std::cout<<"gw.wc_new(iGlob)-gw.wc_old(iGlob)\n"<<gw.wc_new(iGlob)-gw.wc_old(iGlob)<<std::endl;
//!-----------	输出结果显示，gw.wc(iGlob, 0)-gw.wc(iGlob, 1)为0-----------------//
        // if (iGlob == 315){
        //     std::cout << "-------gw.wc_new(iGlob)----- " <<gw.wc_new(iGlob) << std::endl;
        //     std::cout << "-------gw.wc_old(iGlob)----- " <<gw.wc_old(iGlob) << std::endl;
        //     std::cout << "-------gw.wc_new(iGlob)-gw.wc_old(iGlob)------ " <<gw.wc_new(iGlob)-gw.wc_old(iGlob) << std::endl;  
		// 	std::cout << "-------gw.wc(iGlob, 0)----- " <<gw.wc(iGlob, 0) << std::endl;  
		// 	std::cout << "-------gw.wc(iGlob, 1)----- " <<gw.wc(iGlob, 1) << std::endl;        
        // }  

//!矩阵系数（修改使用new和old）
// gw.wc_new(iGlob) = gw.wc(iGlob, 1);
// gw.wc_old(iGlob) = gw.wc(iGlob, 0);
			rt.RTcoef(idom, 0) = 1.0 * gw.wc_new(iGlob) + rt.c_difxx(iGlob) + rt.c_difxx(iGlob - 1) + rt.Up_Weighting_xp(iGlob) * rt.c_advxx(iGlob) - (1 - rt.Up_Weighting_xm(iGlob)) * rt.c_advxx(iGlob - 1) 
			+ rt.c_difyy(iGlob) + rt.c_difyy(iGlob - gdom.nxhc) + rt.Up_Weighting_yp(iGlob) * rt.c_advyy(iGlob) - (1 - rt.Up_Weighting_ym(iGlob)) * rt.c_advyy(iGlob - gdom.nxhc) 
			+ rt.c_difzz(iGlob) /gdom.dz(iGlob) + rt.c_difzz(iGlob - gdom.nxhc * gdom.nyhc)/gdom.dz(iGlob) + rt.Up_Weighting_zp(iGlob) * rt.c_advzz(iGlob)/gdom.dz(iGlob) - (1 - rt.Up_Weighting_zm(iGlob)) * rt.c_advzz(iGlob - gdom.nxhc * gdom.nyhc)/gdom.dz(iGlob);
			// + (gw.wc_new(iGlob) - gw.wc_old(iGlob));
			// Ci+1
			rt.RTcoef(idom, 1) = -(rt.c_difxx(iGlob) - (1 - rt.Up_Weighting_xp(iGlob)) * rt.c_advxx(iGlob)); // Ci+1
			//  Ci-1
			rt.RTcoef(idom, 2) = -(rt.Up_Weighting_xm(iGlob) * rt.c_advxx(iGlob - 1) + rt.c_difxx(iGlob - 1)); // # Ci-1
			// Cj+1
			rt.RTcoef(idom, 3) = -(rt.c_difyy(iGlob) - (1 - rt.Up_Weighting_yp(iGlob)) * rt.c_advyy(iGlob)); // Cj+1
			// Cj-1
			rt.RTcoef(idom, 4) = -(rt.Up_Weighting_ym(iGlob) * rt.c_advyy(iGlob - gdom.nxhc) + rt.c_difyy(iGlob - gdom.nxhc)); // Cj-1
			// Ck+1
			rt.RTcoef(idom, 5) = -(rt.c_difzz(iGlob)/gdom.dz(iGlob) - (1 - rt.Up_Weighting_zp(iGlob)) * rt.c_advzz(iGlob)/gdom.dz(iGlob)); // Ck+1
			// Ck-1
			rt.RTcoef(idom, 6) = -(rt.Up_Weighting_zm(iGlob) * rt.c_advzz(iGlob - gdom.nxhc * gdom.nyhc)/gdom.dz(iGlob) + rt.c_difzz(iGlob - gdom.nxhc * gdom.nyhc)/gdom.dz(iGlob)); // Ck-1
			//rhs
			rt.RTcoef(idom, 7) = rt.c(iGlob, 1) * gw.wc_old(iGlob);

			rt.wz(iGlob) = 0.5;

//!矩阵系数，print矩阵模拟用
			// //矩阵系数（修改使用new和old）
			// rt.RTcoef(idom, 0) = 100;
			// // Ci+1
			// rt.RTcoef(idom, 1) = 1; // Ci+1
			// //  Ci-1
			// rt.RTcoef(idom, 2) = 2; // # Ci-1
			// // Cj+1
			// rt.RTcoef(idom, 3) = 3; // Cj+1
			// // Cj-1
			// rt.RTcoef(idom, 4) = 4; // Cj-1
			// // Ck+1
			// rt.RTcoef(idom, 5) = 5; // Ck+1
			// rt.RTcoef(idom, 6) = 6; // Ck-1
			// //rhs
			// rt.RTcoef(idom, 7) = 7;
//! A矩阵系数，交叉弥散项
			//ci,j+1//添加到上面
			rt.RTcoef(idom, 8) = -0.5*rt.c_difxy(iGlob) + 0.5*rt.c_difxy(iGlob-1)-rt.wz(iGlob)*rt.c_difzy(iGlob)/ (2 * gdom.dy *  gdom.dz(iGlob))+(1-rt.wz(iGlob-gdom.nxhc*gdom.nyhc))*rt.c_difzy(iGlob- gdom.nxhc * gdom.nyhc)/(2 * gdom.dy *  gdom.dz(iGlob)); // Ci,j+1
			rt.RTcoef(idom, 3) += rt.RTcoef(idom, 8); 

			//Ci+1,j+1
			rt.RTcoef(idom, 9) = -0.5*rt.c_difxy(iGlob) - 0.5*rt.c_difyx(iGlob); // Ci+1,j+1
			
			//Cj-1//添加到上面
			rt.RTcoef(idom, 10) = 0.5*rt.c_difxy(iGlob) - 0.5*rt.c_difxy(iGlob-1) + rt.wz(iGlob)*rt.c_difzy(iGlob) / (2 * gdom.dy *  gdom.dz(iGlob))- (1-rt.wz(iGlob-gdom.nxhc*gdom.nyhc))*rt.c_difzy(iGlob-gdom.nxhc * gdom.nyhc)/(2 * gdom.dy *  gdom.dz(iGlob)); // Cj-1
			
			rt.RTcoef(idom, 4) += rt.RTcoef(idom, 10);
			//Ci+1,j-1
			rt.RTcoef(idom, 11) = 0.5*rt.c_difxy(iGlob) + 0.5*rt.c_difyx(iGlob-gdom.nxhc); // Ci+1,j-1
			//Ci-1,j+1
			rt.RTcoef(idom, 12) = 0.5*rt.c_difxy(iGlob-1) + 0.5*rt.c_difyx(iGlob); // Ci-1,j+1
			//Ci-1,j-1
			rt.RTcoef(idom, 13) = -0.5*rt.c_difxy(iGlob-1) - 0.5*rt.c_difyx(iGlob-1); // Ci-1,j-1
			//Ck+1//添加到上面
			rt.RTcoef(idom, 14) = -0.5*rt.c_difxz(iGlob)/ (gdom.dx * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))) + 0.5*rt.c_difxz(iGlob-1) / (gdom.dx * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))) - 0.5*rt.c_difyz(iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc)))+ 0.5 *rt.c_difyz(iGlob-gdom.nxhc)/ (gdom.dy * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))); // Ck+1
			rt.RTcoef(idom, 5) += rt.RTcoef(idom, 14);
			//!Ci+1,k+1
			rt.RTcoef(idom, 15) = -0.5*rt.c_difxz(iGlob)/ (gdom.dx * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))) - (1-rt.wz(iGlob))*rt.c_difzx(iGlob)/ (2 * gdom.dx *  gdom.dz(iGlob)); // Ci+1,k+1
			//Ck-1//添加到上面
			rt.RTcoef(idom, 16) = 0.5*rt.c_difxz(iGlob)/ (gdom.dx * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc)))  - 0.5*rt.c_difxz(iGlob-1) / (gdom.dx * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))) + 0.5*rt.c_difyz(iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc)))- 0.5*rt.c_difyz(iGlob-gdom.nxhc)/ (gdom.dy * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))); // Ck-1
			rt.RTcoef(idom, 6) += rt.RTcoef(idom, 16);
			//!Ci+1,k-1
			rt.RTcoef(idom, 17) = 0.5*rt.c_difxz(iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))) + rt.wz(iGlob-gdom.nxhc*gdom.nyhc) * rt.c_difzx(iGlob-gdom.nxhc*gdom.nyhc)/ (2 * gdom.dx *  gdom.dz(iGlob)); // Ci+1,k-1
			//Ci-1,k+1
			rt.RTcoef(idom, 18) = 0.5*rt.c_difxz(iGlob-1) / (gdom.dx * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))) + (1-rt.wz(iGlob))*rt.c_difzx(iGlob)/ (2 * gdom.dx *  gdom.dz(iGlob)); // Ci-1,k+1
			//!Ci-1,k-1
			rt.RTcoef(idom, 19) = -0.5*rt.c_difxz(iGlob-1) / (gdom.dx * (0.5 * gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))) - rt.wz(iGlob-gdom.nxhc*gdom.nyhc)*rt.c_difzx(iGlob-gdom.nxhc*gdom.nyhc)/ (2 * gdom.dx *  gdom.dz(iGlob)); // Ci-1,k-1
			//Ci+1///添加到上面
			rt.RTcoef(idom, 20) = -0.5*rt.c_difyx(iGlob) + 0.5*rt.c_difyx(iGlob-gdom.nxhc) - rt.wz(iGlob)*rt.c_difzx(iGlob) / (2 * gdom.dx *  gdom.dz(iGlob))+ (1-rt.wz(iGlob-gdom.nxhc*gdom.nyhc))*rt.c_difzx(iGlob-gdom.nxhc*gdom.nyhc)/ (2 * gdom.dx *  gdom.dz(iGlob)); // Ci+1
			rt.RTcoef(idom, 1) += rt.RTcoef(idom, 20);
			//Ci-1//添加到上面
			rt.RTcoef(idom, 21) = 0.5*rt.c_difyx(iGlob) - 0.5*rt.c_difyx(iGlob-gdom.nxhc) + rt.wz(iGlob)*rt.c_difzx(iGlob) / (2 * gdom.dx *  gdom.dz(iGlob))- (1-rt.wz(iGlob-gdom.nxhc*gdom.nyhc))*rt.c_difzx(iGlob-gdom.nxhc*gdom.nyhc)/ (2 * gdom.dx *  gdom.dz(iGlob)); // Ci-1
			rt.RTcoef(idom, 2) += rt.RTcoef(idom, 21);
			//Cj+1,k+1
			rt.RTcoef(idom, 22) = -0.5*rt.c_difyz(iGlob) / ( gdom.dy * (0.5*gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc)))- (1-rt.wz(iGlob))*rt.c_difzy(iGlob)/ (2 * gdom.dy *  gdom.dz(iGlob)); // Cj+1,k+1
			//Cj+1,k-1
			rt.RTcoef(idom, 23) = 0.5*rt.c_difyz(iGlob) / (gdom.dy * (0.5*gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc)))+ rt.wz(iGlob-gdom.nxhc*gdom.nyhc)*rt.c_difzy(iGlob-gdom.nxhc*gdom.nyhc)/ (2 * gdom.dy *  gdom.dz(iGlob)); // Cj+1,k-1
			//Cj-1,k+1
			rt.RTcoef(idom, 24) = 0.5*rt.c_difyz(iGlob-gdom.nxhc) / (gdom.dy * (0.5*gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc)))+ (1-rt.wz(iGlob))*rt.c_difzy(iGlob)/ (2 * gdom.dy *  gdom.dz(iGlob)); // Cj-1,k+1
			//Cj-1,k-1
			rt.RTcoef(idom, 25) = -0.5*rt.c_difyz(iGlob-gdom.nxhc) / ( gdom.dy * (0.5*gdom.dz(iGlob-gdom.nxhc*gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob+gdom.nxhc*gdom.nyhc))) - rt.wz(iGlob-gdom.nxhc*gdom.nyhc)*rt.c_difzy(iGlob-gdom.nxhc*gdom.nyhc)/ (2 * gdom.dy *  gdom.dz(iGlob)); // Cj-1,k-1
//!交叉弥散项模拟值，print矩阵用

// !考虑反应源汇项

			double intermediate_variable = 0;//定义非平衡吸附中间变量
			double min_conc = 1e-6;//定义最小浓度值


			if (rt.ReactionModule == 1)
			{
				// 阻滞因子计算,根据前一次迭代计算浓度进行更新
				if (rt.AdsorptionDesorptionModel == Equilibrium_Linear_Model) // 表面平衡-线性吸附模式;
				{
					rt.Rf(iGlob) = 1 + rt.rho_b * rt.Kd / gw.wc_old(iGlob);
				}
				else if (rt.AdsorptionDesorptionModel == Equilibrium_Freundlich_Model) // 2:Freundlich等温吸附模式;
				{
					
					rt.Rf(iGlob) = 1 + (rt.rho_b * rt.Kf * rt.Nf * pow(rt.c(iGlob, 1), (rt.Nf - 1))) / gw.wc_old(iGlob);
					//如果浓度值低于1e-6，则Rf=1，防止Rf值过大影响矩阵求解
					if (rt.c(iGlob, 1) < min_conc){
						rt.Rf(iGlob) = 1;
					}
				}
				else if (rt.AdsorptionDesorptionModel == Equilibrium_Langmuir_Model) // 3:Langmuir等温吸附模式（固体表面吸附位有限
				{
					rt.Rf(iGlob) = 1 + (rt.rho_b  / gw.wc_old(iGlob)) * ((rt.Kl ) / (pow((1 + rt.eta * rt.c(iGlob, 1)), 2)) );

				}
				// else
				// {
				// 	std::cerr << RERROR "Check parameter AdsorptionDesorptionModel \n";
				// 	exit(-1);
				// }
				
				else if (rt.AdsorptionDesorptionModel == Nonequilibrium_Model)//4:非平衡吸附模式
				{					
					intermediate_variable = rt.beta / (rt.Kd * (rt.rho_b / rt.dt + rt.beta / rt.Kd + rt.lambda_2 * rt.rho_b) );
					rt.Rf(iGlob) = 1;
				}
				else if (rt.AdsorptionDesorptionModel != 1 && rt.AdsorptionDesorptionModel != 2 && rt.AdsorptionDesorptionModel != 3 && rt.AdsorptionDesorptionModel != 4)
				{
					std::cout << "parameter AdsorptionDesorptionModel must be 1 or 2 or 3 or 4" << std::endl;
				}
//更新吸附反应后的矩阵系数
				rt.RTcoef(idom, 0) = rt.Rf(iGlob) * gw.wc_new(iGlob) + rt.lambda * gw.wc_new(iGlob) 
				+ rt.c_difxx(iGlob) + rt.c_difxx(iGlob - 1) + rt.Up_Weighting_xp(iGlob) * rt.c_advxx(iGlob) - (1 - rt.Up_Weighting_xm(iGlob)) * rt.c_advxx(iGlob - 1) 
				+ rt.c_difyy(iGlob) + rt.c_difyy(iGlob - gdom.nxhc) + rt.Up_Weighting_yp(iGlob) * rt.c_advyy(iGlob) - (1 - rt.Up_Weighting_ym(iGlob)) * rt.c_advyy(iGlob - gdom.nxhc) 
				+ rt.c_difzz(iGlob) /gdom.dz(iGlob) + rt.c_difzz(iGlob - gdom.nxhc * gdom.nyhc)/gdom.dz(iGlob) + rt.Up_Weighting_zp(iGlob) * rt.c_advzz(iGlob)/gdom.dz(iGlob) - (1 - rt.Up_Weighting_zm(iGlob)) * rt.c_advzz(iGlob - gdom.nxhc * gdom.nyhc)/gdom.dz(iGlob)
				+ (gw.wc_new(iGlob) - gw.wc_old(iGlob))
				+ rt.dt * (rt.beta + (rt.lambda_1 * gw.wc_new(iGlob)) - (intermediate_variable * pow (rt.beta,2)));
			

				rt.RTcoef(idom, 7) = rt.c(iGlob, 1) * gw.wc_old(iGlob)* rt.Rf(iGlob)
				+ intermediate_variable * rt.rho_b * rt.c_solid(iGlob,1);

			}
			if (rt.ReactionModule != 1 && rt.ReactionModule != 0)
			{
				std::cout << "parameter ReactionModule must be 0 or 1" << std::endl;
			}

			

			if (rt.rt_scheme == 2) // picard迭代
			{

				rt.RTcoef(idom, 7) -= (gw.wc(iGlob, 1) - gw.wc(iGlob, 0)); // 为了消除迭代中的线性误差
				// rt.RTcoef(idom,7) = rt.RTcoef(idom,7);
			} });
		/*--------------后续完善部分*/
		// Apply internal boundary conditions (needed when MPI is used)

		// Apply outer boundary conditions
		for (int k = 0; k < rtgbc.size(); k++)
		{
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

		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			
			int ii, jj, kk;
			// idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
			
			//获取当前单元格idom在rtA矩阵中的行起始位置
			//rtA.rt_ptr数组，存储稀疏矩阵中每一行非零元素的起始索引
			//irow指向当前处理的单元格在rtA.rt_ind和rtA.rt_val数组中的起始位置
			int irow = rtA.rt_ptr(idom); 				
			gdom.unpackIndices(idom, kk, jj, ii);//将idom转换为三维坐标
/*
!重新填充，矩阵系数填充矩阵考虑Dxy修改后
*/ 
//ci-1,k-1
			if (ii > 0 && kk > 0){
				rtA.rt_ind(irow) = idom -1 - gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 19);
				// rtA.rt_val(irow) = 19;				
				irow++;				
			}
//Ck-1
		    if (kk > 0)	{
			rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny;// rtA.rt_ind非零元素的列索引
			rtA.rt_val(irow) = rt.RTcoef(idom,6); // rtA.rt_val非零元素的值 
			irow++;
			}	
//ci+1,k-1
			if (ii < gdom.nx-1 && kk > 0){
				rtA.rt_ind(irow) = idom + 1 -gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 17); 
				// rtA.rt_val(irow) = 17;				
				irow++;				
			}	
//ci-1
			if (ii > 0)	{
				rtA.rt_ind(irow) = idom - 1;		        
				rtA.rt_val(irow) = rt.RTcoef(idom,2);  
				irow++;
				}	
//Ci,j,k
			rtA.rt_ind(irow) = idom;	
			rtA.rt_val(irow) = rt.RTcoef(idom,0);	
			irow++;
//Ci+1
			if (ii < gdom.nx-1)	{
				rtA.rt_ind(irow) = idom + 1;		        
				rtA.rt_val(irow) = rt.RTcoef(idom,1);  
				irow++;			
				}			
//ci-1,k+1
			if (ii > 0 && kk < gdom.nz-1){
				rtA.rt_ind(irow) = idom -1 + gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 18);
				// rtA.rt_val(irow) = 18;				 
				irow++;				
			}				

//Ck+1
			if (kk < gdom.nz-1)	{
				rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny;	
				rtA.rt_val(irow) = rt.RTcoef(idom,5);  
				irow++;
				}
//ci+1,k+1
			if (ii < gdom.nx-1 && kk < gdom.nz-1){
				rtA.rt_ind(irow) = idom + 1 + gdom.nx * gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 15); 
				// rtA.rt_val(irow) = 15;				
				irow++;				
			}	

			rtA.rt_rhs(idom) = rt.RTcoef(idom,7); 
			});
	}

	/*----------------20240426修改-----------*/ // 非迭代方法时间步长控制
	inline void dt_con(RTState &rt, GwState &gw, GwDomain &gdom)
	{
		real dc_max, dt_old;
		dt_old = rt.dt; // 在serghei.h中初始化rt.dt=gdom.dt_init;
		Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idx, kk, jj, ii);
            // gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            real dc = fabs(rt.c(iGlob,1) - rt.c(iGlob,0));
			tmp = (dc > tmp) ? dc : tmp; }, Kokkos::Max<real>(dc_max));
		if (dc_max > 0.02)
		{
			rt.dt = rt.dt * 0.9;
		}
		else if (dc_max >= 0.0 & dc_max < 0.01)
		{
			rt.dt = rt.dt * 1.1;
		}
		if (rt.dt > gdom.dt_max)
		{
			rt.dt = gdom.dt_max;
		}
		else if (rt.dt < gdom.dt_init)
		{
			rt.dt = gdom.dt_init;
		}
	}
	/*----------------20240426修改-----------*/ // 非迭代方法时间步长控制

	inline void dt_iter(RTState &rt, GwState &gw, GwDomain &gdom, int iter)
	{
		real dt_old;
		dt_old = rt.dt;
		if (iter < 7)
		{
			rt.dt = rt.dt * 1.1;
		}
		else if (iter > 11)
		{
			rt.dt = rt.dt * 0.9;
		}
		if (rt.dt > gdom.dt_max)
		{
			rt.dt = gdom.dt_max;
		}
		else if (rt.dt < gdom.dt_init)
		{
			rt.dt = gdom.dt_init;
		}
	}

	inline real get_eps(RTState &rt, GwDomain &gdom)
	{
		real eps;
		Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
            real dwc = fabs(rt.c(idx,1) - rt.c(idx,0));
			tmp = (dwc > tmp) ? dwc : tmp; }, Kokkos::Max<real>(eps));
		return eps;
	}
};
#endif