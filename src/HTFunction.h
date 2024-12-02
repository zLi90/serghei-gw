#ifndef _HT_FUNCTION_H_
#define _HT_FUNCTION_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
#include "GwState.h"
#include "HTMatrix.h"
#include "HTSolver.h"
#include "HTState.h"
#include "State.h"
#include <set>
#include <limits>

class HTFunction
{

private:
	Kokkos::Timer timer;

public:
	/* --------------------------------------------------
		ht_Gauss_Seidel_solve
	-------------------------------------------------- */
	template <typename execution_space, typename type_solver>
	inline void ht_Gauss_Seidel_solve(HTState &ht, RTState &rt, HTMatrix &htA, GwState &gw, GwDomain &gdom, std::vector<HTBC> &htgbc, GwMPI &gmpi, Parallel &par, type_solver &htsolver)
	{
		int iter, ierr = 1;
		real dt_tmp;

		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}
		// 扩散系数计算
		//dispersion_tensor(ht, gw, gdom);

		HTlinear_system(gw, gdom, htA, par, ht, htgbc);

		timer.reset();

		htsolver.Gauss_Seidel(htA);//求解含水层温度
		gdom.timers.solver += timer.seconds();

		// Update temperature
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            ht.T(iGlob,1) = htA.ht_x(idom); //0为n+1时刻，1为n时刻
								
			});						
		gmpi.mpi_sendrecv(ht.T, gdom, par);
		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}
		gmpi.mpi_sendrecv(ht.T, gdom, par);
		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}

		dt_tem(ht, gw, gdom, rt); // 计算温度残差值，并更新
		// Update time step
		dt_tmp = gdom.dt;

		ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { ht.T(iGlob, 0) = ht.T(iGlob, 1); });
	}
	/* --------------------------------------------------
		Top-level PCA solver
	-------------------------------------------------- */
	template <typename execution_space, typename type_solver>
	inline void ht_pca_solve(HTState &ht,RTState &rt, HTMatrix &htA, GwState &gw, GwDomain &gdom, std::vector<HTBC> &htgbc, GwMPI &gmpi, Parallel &par, type_solver &htsolver)
	{
		int iter, ierr = 1;
		real dt_tmp;

		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}
		// 扩散系数计算
		//dispersion_tensor(rt, gw, gdom);

		HTlinear_system(gw, gdom, htA, par, ht, htgbc);

		timer.reset();
#if SERGHEI_KOKKOSKERNELS_SOLVER
		htsolver.kkpcg(htA);
#else
		iter = htsolver.cg(htA, gdom);
#endif
		gdom.timers.solver += timer.seconds();

		// Update temperature
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            ht.T(iGlob,1) = htA.ht_x(idom); });
		gmpi.mpi_sendrecv(ht.T, gdom, par);
		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}
		gmpi.mpi_sendrecv(ht.T, gdom, par);
		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}

		dt_tem(ht, gw, gdom, rt); // 计算温度残差值，并更新
		// Update time step
		dt_tmp = gdom.dt;

		ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { ht.T(iGlob, 0) = ht.T(iGlob, 1); });
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
	inline void ht_picard_solve(HTState &ht, RTState &rt, HTMatrix &htA, GwState &gw, GwDomain &gdom, std::vector<HTBC> &htgbc,
								type_solver &htsolver, GwMPI &gmpi, GwIntegrator &gint, Parallel &par)
	{

		int iter, iter_cg, iter_max = 100, ierr = 1;
		real eps_diff = 1.0, eps_tmp, eps_old = 1.0, eps = 1.0, dt_tmp;
		real eps_min = 5e-5; // eps_min = 5e-6,20240510修改了收敛残差标准
		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc; });
		iter = 0;
		// 计算扩散系数
		//dispersion_tensor(rt, gw, gdom);

		while (iter < iter_max && eps_diff / eps_old > eps_min && eps > eps_min)
		{

			HTlinear_system(gw, gdom, htA, par, ht, htgbc);

			timer.reset();
#if SERGHEI_KOKKOSKERNELS_SOLVER
			htsolver.kkpcg(htA);
#else
			iter = htsolver.cg(htA, gdom);
#endif
			gdom.timers.solver += timer.seconds();
			Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idom, kk, jj, ii);
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                ht.T(iGlob,0) = ht.T(iGlob,1);
                ht.T(iGlob,1) = htA.ht_x(idom); });
			gmpi.mpi_sendrecv(ht.T, gdom, par);

			//dispersion_tensor(rt, gw, gdom);

			eps_old = eps;
			eps = get_eps(ht, gdom);
			eps_tmp = fabs(eps_old - eps);
			ierr = MPI_Allreduce(&eps_tmp, &eps_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			iter += 1;
		}
		// printf("    > Picard loop converges in %d iterations with eps = %f, %f\n",iter,eps,eps_diff);

		dt_iter(ht, rt, gw, gdom, iter);
		dt_tmp = ht.dt;
		ierr = MPI_Allreduce(&dt_tmp, &ht.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { ht.T(iGlob, 0) = ht.T(iGlob, 1); });
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

	//20241105----------------//
	template <typename execution_space, typename type_solver>
	inline void ht_solve(HTState &ht, RTState &rt, HTMatrix &htA, GwState &gw, GwDomain &gdom, std::vector<HTBC> &htgbc, GwMPI &gmpi, Parallel &par, type_solver &htsolver)
	{
		int iter, ierr = 1;
		real dt_tmp;

		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}
		// 扩散系数计算
		//dispersion_tensor(ht, gw, gdom);

		//advection_dispersion(ht, gw, gdom, dt_tmp);

		timer.reset();

		//htsolver.Gauss_Seidel(htA);//求解含水层温度
		gdom.timers.solver += timer.seconds();
		iter = htsolver.cg(htA, gdom);

		// Update temperature
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
			ht.Cp(iGlob, 1) = (1.92* ht.theta_n + 4.18 * gw.wc(iGlob, 1)) * 1e6;
			ht.lambda(iGlob, 1) = (ht.b1 + ht.b2 * gw.wc(iGlob, 1) + ht.b3 * pow(gw.wc(iGlob, 1), 0.5))/864 ;
            ht.T(iGlob,1) = ht.T(iGlob, 1) + ht.dt * ht.lambda(iGlob, 1) / ht.Cp(iGlob, 1) * ( (ht.T(iGlob + 1, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - 1, 1)) / (gdom.dx * gdom.dx) 
				            + (ht.T(iGlob + gdom.nxhc, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - gdom.nxhc, 1)) / (gdom.dy * gdom.dy) 
							+ (ht.T(iGlob + gdom.nxhc * gdom.nyhc, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - gdom.nxhc * gdom.nyhc, 1)) / (gdom.dz(iGlob) * gdom.dz(iGlob)) )
							- ht.dt * ht.Cw / ht.Cp(iGlob, 1) * ( gw.q(iGlob, 0) * (ht.T(iGlob, 1) - ht.T(iGlob - 1, 1)) / gdom.dx + gw.q(iGlob, 1) * (ht.T(iGlob, 1) - ht.T(iGlob - gdom.nxhc, 1)) / gdom.dy 
							+ gw.q(iGlob, 2) * (ht.T(iGlob, 1) - ht.T(iGlob - gdom.nxhc * gdom.nyhc, 1)) / gdom.dz(iGlob)) ; //0为n+1时刻，1为n时刻
								
			});						
		gmpi.mpi_sendrecv(ht.T, gdom, par);
		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}
		gmpi.mpi_sendrecv(ht.T, gdom, par);
		// 施加边界条件
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyTemperatureBC(ht, gw, gdom, par);
		}

		dt_tem(ht, gw, gdom, rt); // 计算温度残差值，并更新
		// Update time step
		dt_tmp = gdom.dt;

		ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { ht.T(iGlob, 0) = ht.T(iGlob, 1); });
	}


	/*template <typename execution_space, typename type_solver>
	inline void ht_solve(HTState &ht, GwState &gw, GwDomain &gdom, GwMPI &gmpi, Parallel &par, type_solver &htsolver)
	{
		int ierr;//CFL number
        real dt_tmp;
        real cfl = 0.01; // CFL number
		real tolerance = 1e-3;
	    // real dt = dt_cfl; // 初始化时间步长为 CFL 时间步长
		real dt_cfl = 0.001;
		dt_tmp = dt_cfl;
		
		while (ierr < 100)
		{
			ht.dt = dt_tmp;
			advection_dispersion(ht, gw, gdom, dt_tmp);
			
			if (check_convergence(ht, gw, gdom, tolerance))
			{
				ht.dt =dt_tmp;
				std::cout << "-------debugzzb-rt.dt00000----- " << ht.dt << std::endl;
				
						// Update concentration arrays
				Kokkos::parallel_for(
		    			gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			    		ht.T(iGlob, 1) = ht.T(iGlob, 0);//收敛，则将更新后的Cn值传递给Cn+1
		   			 });
		    		gmpi.mpi_sendrecv(ht.T, gdom, par);

				
				break;
			}
			else
			{
				dt_tmp *=0.5;
				ierr++;
				std::cout << "-------debugzzb-ierr----- " << ierr << std::endl;
			
			if (ierr >= 100 || dt_tmp < 1e-3) 
			{
			std::cout << "-------达到最大迭代次数，认为收敛----- " << std::endl;
			// 更新温度数组
                Kokkos::parallel_for(
                    gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
                        ht.T(iGlob, 1) = ht.T(iGlob, 0); // 如果达到最大迭代次数，也将更新后的 Cn 值传递给 Cn+1
                    });
                gmpi.mpi_sendrecv(ht.T, gdom, par);

                break;
			}
			else
			{
				std::cout << "-------未达到最大迭代次数----- " << std::endl;

			// std::cout << "-------dt_tmp1111----- " << dt_tmp << std::endl;
			//查看未收敛情况下，重新计算前的误差
			Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
				ht.T(iGlob, 1) = ht.T(iGlob, 0);
			});
			}
		}
		}
	}*/

	real find_max_residual(const HTState &ht)
	{
		real max_residual = 0.0;
		// 遍历所有单元格
	    for (int i = 0; i < ht.residual.extent(0); ++i) 
	    {
		// 计算当前单元格的残差值
		real residual_value = ht.residual(i);
		
		// 更新最大残差值
		max_residual = fmax(max_residual, residual_value);
	    }

	return max_residual;
	}

	inline void advection_dispersion(HTState &ht, GwState &gw, GwDomain &gdom, real dt_tmp)
	{
		Kokkos::parallel_for(
		    gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
				ht.Cp(iGlob, 1) = (1.92* ht.theta_n + 4.18 * gw.wc(iGlob, 1)) * 1e6;
			    ht.lambda(iGlob, 1) = ht.b1 + ht.b2 * gw.wc(iGlob, 1) + ht.b3 * pow(gw.wc(iGlob, 1), 0.5) ;

				ht.T(iGlob, 1) = ht.T(iGlob, 1) + ht.dt * ht.lambda(iGlob, 1) / ht.Cp(iGlob, 1) * ( (ht.T(iGlob + 1, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - 1, 1)) / (gdom.dx * gdom.dx) 
				            + (ht.T(iGlob + gdom.nxhc, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - gdom.nxhc, 1)) / (gdom.dy * gdom.dy) 
							+ (ht.T(iGlob + gdom.nxhc * gdom.nyhc, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - gdom.nxhc * gdom.nyhc, 1)) / (gdom.dz(iGlob) * gdom.dz(iGlob)) )
							- ht.dt * ht.Cw / ht.Cp(iGlob, 1) * ( gw.q(iGlob, 0) * (ht.T(iGlob, 1) - ht.T(iGlob - 1, 1)) / gdom.dx + gw.q(iGlob, 1) * (ht.T(iGlob, 1) - ht.T(iGlob - gdom.nxhc, 1)) / gdom.dy 
							+ gw.q(iGlob, 2) * (ht.T(iGlob, 1) - ht.T(iGlob - gdom.nxhc * gdom.nyhc, 1)) / gdom.dz(iGlob)) ;

			});
	}

	bool check_convergence(HTState &ht, GwState &gw, GwDomain &gdom, real tolerance)
	{
		
		bool converged = true;
		Kokkos::parallel_reduce(
		    "CheckConvergence", ht.T.extent(0), KOKKOS_LAMBDA(int iGlob, bool &converged) {
			    // Compute concentration change

			    ht.residual(iGlob) = fabs((ht.T(iGlob, 0) - ht.T(iGlob, 1)) / ht.dt - ( ht.lambda(iGlob, 1) / ht.Cp(iGlob, 1) * ( (ht.T(iGlob + 1, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - 1, 1)) / (gdom.dx * gdom.dx) 
				            + (ht.T(iGlob + gdom.nxhc, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - gdom.nxhc, 1)) / (gdom.dy * gdom.dy) 
							+ (ht.T(iGlob + gdom.nxhc * gdom.nyhc, 1) - 2 * ht.T(iGlob, 1) + ht.T(iGlob - gdom.nxhc * gdom.nyhc, 1)) / (gdom.dz(iGlob) * gdom.dz(iGlob)) )
							-  ht.Cw / ht.Cp(iGlob, 1) * ( gw.q(iGlob, 0) * (ht.T(iGlob, 1) - ht.T(iGlob - 1, 1)) / gdom.dx + gw.q(iGlob, 1) * (ht.T(iGlob, 1) - ht.T(iGlob - gdom.nxhc, 1)) / gdom.dy 
							+ gw.q(iGlob, 2) * (ht.T(iGlob, 1) - ht.T(iGlob - gdom.nxhc * gdom.nyhc, 1)) / gdom.dz(iGlob))));
			    if (ht.residual(iGlob) > tolerance)
			    {
				    converged = true;
				   
			    }
			    
		    },
		    Kokkos::LAnd<bool>(converged)
		    );
		
		
		return converged;
					    
	}
	//20241105-----------------//


	// 寻找 ht.T(iGlob, 1) 中的最大值及其索引
	std::pair<real, int> find_max_value(const HTState &ht)
	{
		real max_value = 0;
		int max_index = 0;
		for (int i = 1; i < ht.T.extent(0); ++i)
		{
			real current_value = ht.T(i, 1);
			if (current_value > max_value)
			{
				max_value = current_value;
				max_index = i;
			}
		}
		return std::make_pair(max_value, max_index);
	}

	// 计算Peclet
	/*std::pair<real, int> find_peclet(const HTState &ht, GwDomain &gdom)
	{
		real max_value = 0;
		int max_index = 0;
		for (int i = 1; i < ht.T.extent(0); ++i)
		{
			real current_value = fabs(gdom.dx * ht.aveV(i, 3) / ht.dcal(i, 0));
			if (current_value > max_value)
			{
				max_value = current_value;
				max_index = i;
			}
		}
		return std::make_pair(max_value, max_index);
	}*/

	// 计算最大孔隙均速
	/*std::pair<real, int> find_aveV(const RTState &rt)
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
	}*/

	// Calculation of the hydrodynamic dispersion coefficient tensor
	/*inline void dispersion_tensor(RTState &rt, GwState &gw, GwDomain &gdom)
	{
		Kokkos::parallel_for(
			gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
				// rt.c(iGlob, 0) 表示下一时间步长的浓度值
				// rt.c(iGlob, 1) 表示当前步长的浓度值

				// 平均孔隙流速计算，平均孔隙流速计算,x,y,z和平均流速方向
				//  int ii, jj, kk, ivg;
				//  real  wcs;
				//  gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
				//  ivg = gw.soilID(iGlob) * NVG;
				//  wcs = gw.vgTable(ivg+2);
				// aveVB为边界流速
				rt.aveVB(iGlob, 0) = -(gw.q(iGlob, 0));
				rt.aveVB(iGlob, 1) = (gw.q(iGlob, 1));
				rt.aveVB(iGlob, 2) = -(gw.q(iGlob, 2));
				// rt.aveVB(iGlob, 2) = 5.787e-7;

				rt.aveVB(iGlob, 3) = sqrt(pow(rt.aveVB(iGlob, 0), 2) + pow(rt.aveVB(iGlob, 1), 2) + pow(rt.aveVB(iGlob, 2), 2));

				// aveV为节点流速
				rt.aveV(iGlob, 0) = (rt.aveVB(iGlob, 0) + rt.aveVB(iGlob - 1, 0)) / 2;
				rt.aveV(iGlob, 1) = (rt.aveVB(iGlob, 1) + rt.aveVB(iGlob - gdom.nxhc, 1)) / 2;
				rt.aveV(iGlob, 2) = (rt.aveVB(iGlob, 2) + rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2)) / 2;
				rt.aveV(iGlob, 3) = sqrt(pow(rt.aveV(iGlob, 0), 2) + pow(rt.aveV(iGlob, 1), 2) + pow(rt.aveV(iGlob, 2), 2));

				// 机械弥散系数计算

				// Dxx
				rt.dcal(iGlob, 0) = (rt.alpha_T * (pow(rt.aveVB(iGlob, 1), 2) + pow(rt.aveVB(iGlob, 2), 2)) + rt.alpha_L * pow(rt.aveVB(iGlob, 0), 2)) / rt.aveVB(iGlob, 3);
				// DYY
				rt.dcal(iGlob, 1) = (rt.alpha_T * (pow(rt.aveVB(iGlob, 0), 2) + pow(rt.aveVB(iGlob, 2), 2)) + rt.alpha_L * pow(rt.aveVB(iGlob, 1), 2)) / rt.aveVB(iGlob, 3);
				// DZZ
				rt.dcal(iGlob, 2) = (rt.alpha_T * (pow(rt.aveVB(iGlob, 0), 2) + pow(rt.aveVB(iGlob, 1), 2)) + rt.alpha_L * pow(rt.aveVB(iGlob, 2), 2)) / rt.aveVB(iGlob, 3);
				// Dxy
				rt.dcal(iGlob, 3) = fabs((rt.alpha_L - rt.alpha_T) * rt.aveVB(iGlob, 0) * rt.aveVB(iGlob, 1) / rt.aveVB(iGlob, 3));
				// Dxz
				rt.dcal(iGlob, 4) = fabs((rt.alpha_L - rt.alpha_T) * rt.aveVB(iGlob, 0) * rt.aveVB(iGlob, 2) / rt.aveVB(iGlob, 3));
				// Dyz
				rt.dcal(iGlob, 5) = fabs((rt.alpha_L - rt.alpha_T) * rt.aveV(iGlob, 1) * rt.aveV(iGlob, 2) / rt.aveV(iGlob, 3));

				// rt.dcal(iGlob,0) = 1*1e-7;
				// rt.dcal(iGlob,1) = 0;
				// rt.dcal(iGlob,2) = 0.5*1e-7;
				// 将rt.d中nan值替换
				for (int i = 0; i < 6; ++i)
				{
					if (std::isnan(rt.dcal(iGlob, i)))
					{
						rt.dcal(iGlob, i) = 0;
					}
				};

				//------弥散系数=机械弥散+分子弥散----------

				rt.dcal(iGlob, 0) += gw.wc(iGlob, 0) * rt.diffusion_molecular;
				rt.dcal(iGlob, 1) += gw.wc(iGlob, 0) * rt.diffusion_molecular;
				rt.dcal(iGlob, 2) += gw.wc(iGlob, 0) * rt.diffusion_molecular;

				// 根据研究区域维数对水动力弥散系数进行调整
				if (gdom.nx == 1)
				{
					rt.dcal(iGlob, 0) = 0;
					rt.dcal(iGlob, 4) = 0;
					rt.dcal(iGlob, 5) = 0;
				}
				else if (gdom.ny == 1)
				{
					rt.dcal(iGlob, 1) = 0;
					rt.dcal(iGlob, 3) = 0;
					rt.dcal(iGlob, 5) = 0;
				}
				else if (gdom.nz == 1)
				{
					rt.dcal(iGlob, 2) = 0;
					rt.dcal(iGlob, 3) = 0;
					rt.dcal(iGlob, 4) = 0;
				}
			});
	}*/
	/* --------------------------------------------------
	Get HT matrix coefficients
-------------------------------------------------- */
	inline void HTlinear_system(GwState &gw, GwDomain &gdom, HTMatrix &htA, Parallel &par, HTState &ht, std::vector<HTBC> &htgbc)
	{
		// Calculate matrix coefficients
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, ivg, iGlob, iGlobSW;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
			iGlobSW = (hc + jj) * gdom.nxhc + ii + hc;

			// 添加判断条件，判断速度方向，如果速度方向为正，则取上游加权，否则取下游加权
			// x方向速度方向判断
			/*double Up_Weighting_x, Up_Weighting_y, Up_Weighting_z;
			if (rt.aveV(iGlob, 0) > 0)
			{
				Up_Weighting_x = rt.Up_Weighting_vplus;
			}
			else
			{
				Up_Weighting_x = rt.Up_Weighting_vminus;
			}

			// y方向速度方向判断
			if (rt.aveV(iGlob, 1) > 0)
			{
				Up_Weighting_y = rt.Up_Weighting_vplus;
			}
			else
			{
				Up_Weighting_y = rt.Up_Weighting_vminus;
			}
			// z方向速度方向判断
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

			rt.c_advxx(iGlob) = rt.aveVB(iGlob, 0) * rt.dt / gdom.dx;
			rt.c_advyy(iGlob) = rt.aveVB(iGlob, 1) * rt.dt / gdom.dy;
			rt.c_advzz(iGlob) = rt.aveVB(iGlob, 2) * rt.dt / gdom.dz(iGlob);
			rt.c_difxx(iGlob) = rt.dcal(iGlob, 0) * rt.dt / gdom.dx / gdom.dx;
			rt.c_difyy(iGlob) = rt.dcal(iGlob, 1) * rt.dt / gdom.dy / gdom.dy;
			rt.c_difzz(iGlob) = rt.dcal(iGlob, 2) * rt.dt / gdom.dz(iGlob) / gdom.dz(iGlob);*/


			ht.Cp(iGlob, 1) = (1.92* ht.theta_n + 4.18 * gw.wc(iGlob, 1)) * 1e6;
			ht.lambda(iGlob, 1) = (ht.b1 + ht.b2 * gw.wc(iGlob, 1) + ht.b3 * pow(gw.wc(iGlob, 1), 0.5)) / 864 ;

			// A矩阵系数
			ht.HTcoef(idom, 0) = ht.Cp(iGlob, 1) + 2 * ht.lambda(iGlob, 1) * ht.dt * ( 1 / (gdom.dx * gdom.dx) + 1 / (gdom.dy * gdom.dy) + 1 / (gdom.dz(iGlob) * gdom.dz(iGlob)))
		             + ht.Cw * ht.dt * ( gw.q(iGlob, 0) / gdom.dx + gw.q(iGlob, 1) / gdom.dy + gw.q(iGlob, 2) / gdom.dz(iGlob));

		    ht.HTcoef(idom, 1) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dx  * gdom.dx )); //Ti+1

            ht.HTcoef(idom, 2) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dx  * gdom.dx ) + ht.Cw * gw.q(iGlob, 0) * ht.dt / gdom.dx); //Ti-1

            ht.HTcoef(idom, 3) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dy  * gdom.dy )) ; //Tj+1

            ht.HTcoef(idom, 4) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dy  * gdom.dy ) + ht.Cw * gw.q(iGlob, 1) * ht.dt / gdom.dy); //Tj-1

            ht.HTcoef(idom, 5) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dz(iGlob) * gdom.dz(iGlob))); //Tk+1

            ht.HTcoef(idom, 6) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)) + ht.Cw * gw.q(iGlob, 2) * ht.dt /  gdom.dz(iGlob)); //Tk-1

		    ht.HTcoef(idom, 7) = ht.Cp(iGlob, 1) * ht.T(iGlob, 1);


			/*ht.HTcoef(idom, 0) = ht.Cp(iGlob, 1) + 2 * ht.lambda(iGlob, 1) * ht.dt * ( 1 / (gdom.dx * gdom.dx) + 1 / (gdom.dy * gdom.dy) + 1 / (gdom.dz(iGlob) * gdom.dz(iGlob)));

		    ht.HTcoef(idom, 1) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dx  * gdom.dx )); //Ti+1

            ht.HTcoef(idom, 2) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dx  * gdom.dx )); //Ti-1

            ht.HTcoef(idom, 3) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dy  * gdom.dy )) ; //Tj+1

            ht.HTcoef(idom, 4) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dy  * gdom.dy )); //Tj-1

            ht.HTcoef(idom, 5) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dz(iGlob) * gdom.dz(iGlob))); //Tk+1

            ht.HTcoef(idom, 6) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dz(iGlob) * gdom.dz(iGlob))); //Tk-1

		    ht.HTcoef(idom, 7) = ht.Cp(iGlob, 1) * ht.T(iGlob, 1);*/


			/*ht.HTcoef(idom, 0) = ht.Cp(iGlob, 1) + 2 * ht.lambda(iGlob, 1) * ht.dt * ( 1 / (gdom.dx * gdom.dx) + 1 / (gdom.dy * gdom.dy) + 1 / (gdom.dz(iGlob) * gdom.dz(iGlob)))
		             - ht.Cw * ht.dt * ( gw.q(iGlob, 0) / gdom.dx + gw.q(iGlob, 1) / gdom.dy + gw.q(iGlob, 2) / gdom.dz(iGlob));

		    ht.HTcoef(idom, 1) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dx  * gdom.dx ) - ht.Cw * gw.q(iGlob, 0) * ht.dt / gdom.dx); //Ti+1

            ht.HTcoef(idom, 2) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dx  * gdom.dx )); //Ti-1

            ht.HTcoef(idom, 3) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dy  * gdom.dy ) - ht.Cw * gw.q(iGlob, 1) * ht.dt / gdom.dy) ; //Tj+1

            ht.HTcoef(idom, 4) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dy  * gdom.dy )); //Tj-1

            ht.HTcoef(idom, 5) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)) - ht.Cw * gw.q(iGlob, 2) * ht.dt /  gdom.dz(iGlob)); //Tk+1

            ht.HTcoef(idom, 6) = -( ht.lambda(iGlob, 1) * ht.dt / (gdom.dz(iGlob) * gdom.dz(iGlob))); //Tk-1

		    ht.HTcoef(idom, 7) = ht.Cp(iGlob, 1) * ht.T(iGlob, 1);*/
			

			// 考虑反应源汇项

			//double intermediate_variable;//定义非平衡吸附中间变量


	

			if (ht.ht_scheme == 2) // picard迭代
			{

				ht.HTcoef(idom, 7) -= (gw.wc(iGlob, 1) - gw.wc(iGlob, 0)); // 为了消除迭代中的线性误差
				// rt.RTcoef(idom,7) = rt.RTcoef(idom,7);
			} });
		/*--------------后续完善部分*/
		// Apply internal boundary conditions (needed when MPI is used)

		// Apply outer boundary conditions
		for (int k = 0; k < htgbc.size(); k++)
		{
			htgbc[k].applyHTMatBC(ht, gw, gdom, par);
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

		// Insert coefficients into Matrix htA.

		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			
			int ii, jj, kk;
			// idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;

			int irow = htA.ht_ptr(idom); 
			
			
			gdom.unpackIndices(idom, kk, jj, ii);

		     if (kk > 0)	{htA.ht_ind(irow) = idom - gdom.nx*gdom.ny;	htA.ht_val(irow) = ht.HTcoef(idom,6);  irow++;}
        		if (jj > 0)	{htA.ht_ind(irow) = idom - gdom.nx;	        htA.ht_val(irow) = ht.HTcoef(idom,4);  irow++;}
			if (ii > 0)	{htA.ht_ind(irow) = idom - 1;		        htA.ht_val(irow) = ht.HTcoef(idom,2);  irow++;}
			htA.ht_ind(irow) = idom;	htA.ht_val(irow) = ht.HTcoef(idom,0);	irow++;
			if (ii < gdom.nx-1)	{htA.ht_ind(irow) = idom + 1;		        htA.ht_val(irow) = ht.HTcoef(idom,1);  irow++;}
			if (jj < gdom.ny-1)	{htA.ht_ind(irow) = idom + gdom.nx;	        htA.ht_val(irow) = ht.HTcoef(idom,3);  irow++;}
			if (kk < gdom.nz-1)	{htA.ht_ind(irow) = idom + gdom.nx*gdom.ny;	htA.ht_val(irow) = ht.HTcoef(idom,5);  irow++;}

			htA.ht_rhs(idom) = ht.HTcoef(idom,7); });
	}

	/*----------------20240426修改-----------*/ // 非迭代方法时间步长控制
	inline void dt_tem(HTState &ht, GwState &gw, GwDomain &gdom, RTState &rt)
	{
		real dc_max, dt_old;
		dt_old = ht.dt; // 在serghei.h中初始化ht.dt=gdom.dt_init;
		Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndices(idx, kk, jj, ii);
            // gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            real dc = fabs(ht.T(iGlob,1) - ht.T(iGlob,0));
			tmp = (dc > tmp) ? dc : tmp; }, Kokkos::Max<real>(dc_max));
		if (dc_max > 0.02)
		{
			ht.dt = ht.dt * 0.9;
		}
		else if (dc_max >= 0.0 & dc_max < 0.01)
		{
			ht.dt = ht.dt * 1.1;
		}
		if (ht.dt > rt.dt_max)
		{
			ht.dt = rt.dt_max;
		}
		else if (ht.dt < rt.dt_init)
		{
			ht.dt = rt.dt_init;
		}
	}
	/*----------------20240426修改-----------*/ // 非迭代方法时间步长控制

	inline void dt_iter(HTState &ht, RTState &rt, GwState &gw, GwDomain &gdom, int iter)
	{
		real dt_old;
		dt_old = ht.dt;
		if (iter < 7)
		{
			ht.dt = ht.dt * 1.1;
		}
		else if (iter > 11)
		{
			ht.dt = ht.dt * 0.9;
		}
		if (ht.dt > rt.dt_max)
		{
			ht.dt = rt.dt_max;
		}
		else if (ht.dt < rt.dt_init)
		{
			ht.dt = rt.dt_init;
		}
	}

	inline real get_eps(HTState &ht, GwDomain &gdom)
	{
		real eps;
		Kokkos::parallel_reduce(gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
            real dwc = fabs(ht.T(idx,1) - ht.T(idx,0));
			tmp = (dwc > tmp) ? dwc : tmp; }, Kokkos::Max<real>(eps));
		return eps;
	}
};
#endif