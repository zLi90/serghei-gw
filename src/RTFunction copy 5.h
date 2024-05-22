#ifndef _RT_FUNCTION_H_
#define _RT_FUNCTION_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
#include "GwState.h"
#include "State.h"
#include <set>

class RTFunction
{

private:
	Kokkos::Timer timer;

public:
	/* --------------------------------------------------
		Top-level RT solver
	-------------------------------------------------- */
	inline void rt_solve(RTState &rt, GwState &gw, GwDomain &gdom, GwMPI &gmpi, Parallel &par)
	{
		int ierr; // CFL number
		real dt_tmp;
		real cfl = 0.01; // CFL number
		real tolerance = 1e-3;
	// real dt = dt_cfl; // 初始化时间步长为 CFL 时间步长
		real dt_cfl = 0.001;
		dt_tmp = dt_cfl;

	// 	//施加边界条件
     //    for (int k = 0; k < gbc.size(); k++) {
     //        gbc[k].applyHBC(gw, gdom, par);
     //    }
		//施加边界条件
		// rt.c(617, 1) = 15;


		// rt.dt = dt;
		// while (dt_tmp > 1e-3) {
		while (ierr < 100){
		// while (tolerance < 0.01) {
			// std::cout << "-------debugzzb-dt1111----- " << dt_tmp << std::endl;
			rt.dt = dt_tmp;
			advection_dispersion(rt, gw, gdom, dt_tmp);

		// Check for convergence
			if (check_convergence(rt, tolerance)) {
			// 如果收敛，更新当前时间步长并退出循环
				rt.dt = dt_tmp;
				std::cout << "-------debugzzb-rt.dt00000----- " << rt.dt << std::endl;
				
						// Update concentration arrays
				Kokkos::parallel_for(
		    			gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			    		rt.c(iGlob, 0) = rt.c(iGlob, 1);//收敛，则将更新后的Cn值传递给Cn+1
		   			 });
		    		gmpi.mpi_sendrecv(rt.c, gdom, par);

				real mass = calculate_total_concentration_mass(rt);
				real max_residual_out = find_max_residual(rt);
				std::cout << "-------total_mass000----- " << mass << std::endl;
				std::cout << "-------max_residual_out000----- " << max_residual_out << std::endl;

				break;
		} 
			else {

			// 缩小时间步长
			dt_tmp *= 0.5;
			ierr++;
			std::cout << "-------debugzzb-ierr----- " << ierr << std::endl;
			
			if (ierr >= 100 || dt_tmp < 1e-3) {
			std::cout << "-------达到最大迭代次数，认为收敛----- " << std::endl;
			// 更新浓度数组
                Kokkos::parallel_for(
                    gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
                        rt.c(iGlob, 0) = rt.c(iGlob, 1); // 如果达到最大迭代次数，也将更新后的 Cn 值传递给 Cn+1
                    });
                gmpi.mpi_sendrecv(rt.c, gdom, par);

                // 计算总浓度质量和最大残差
                real mass111 = calculate_total_concentration_mass(rt);
                real max_residual_out111 = find_max_residual(rt);
                std::cout << "-------total_mass111----- " << mass111 << std::endl;
                std::cout << "-------max_residual_out111----- " << max_residual_out111 << std::endl;

                break;
			}
			else {
			std::cout << "-------未达到最大迭代次数----- " << std::endl;

			// std::cout << "-------dt_tmp1111----- " << dt_tmp << std::endl;
			//查看未收敛情况下，重新计算前的误差
			Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
				rt.c(iGlob, 1) = rt.c(iGlob, 0);
			});
			real mass222 = calculate_total_concentration_mass(rt);
			real max_residual_out222 = find_max_residual(rt);
			std::cout << "-------total_mass222----- " << mass222 << std::endl;
			std::cout << "-------max_residual_out222----- " << max_residual_out222 << std::endl;
			std::cout << "-------debugzzb-rt.dt222----- " << rt.dt << std::endl;
			
			}
			// advection_dispersion(rt, gw, gdom, dt_tmp);
			// real mass = calculate_total_concentration_mass(rt);
			// real max_residual_out = find_max_residual(rt);
			// std::cout << "-------total_mass222----- " << mass << std::endl;
			// std::cout << "-------max_residual_out222----- " << max_residual_out << std::endl;
		}

			
		}
	}
/*----------------------------**********************************-*/
	// 计算dt
/*
	inline real calculate_cfl_dt(GwDomain &gdom, real cfl, GwState &gw)
	{
		real max_velocity = 0.0;
		Kokkos::parallel_reduce(
		    "MaxVelocity", gdom.nCellMem, KOKKOS_LAMBDA(int iGlob, real &max_vel) {
			    real u = gw.q(iGlob, 0);
			    real v = gw.q(iGlob, 1);
			    real w = gw.q(iGlob, 2);
			    real vel = sqrt(u * u + v * v + w * w);
			    max_vel = fmax(max_vel, vel);
		    },
		    Kokkos::Max<real>(max_velocity));

		// Ensure denominator is not zero
		max_velocity += 1e-10;

		real min_dz = std::numeric_limits<real>::max(); // 初始化为最大值

		// 使用普通的循环来寻找最小值
		for (int iGlob = 0; iGlob < gdom.nCellMem; ++iGlob) {
		// 更新最小值
		min_dz = fmin(min_dz, gdom.dz(iGlob));
		}


		// Calculate CFL timestep
		real dt_cfl = cfl * fmin(fmin(gdom.dx, gdom.dy), min_dz) / max_velocity;

		return dt_cfl;
	}*/	


	//全网格总质量浓度计算
	real calculate_total_concentration_mass(const RTState &rt) {
	real total_mass = 0.0;

	// 遍历所有网格
	for (int i = 0; i < rt.c.extent(0); ++i) {
		// 将当前网格的浓度质量累加到总质量中
		total_mass += rt.c(i, 0);
	}
	
	return total_mass;
	}

//寻找最大残差值
	real find_max_residual(const RTState &rt) {
	real max_residual = 0.0;

	// 遍历所有单元格
	for (int i = 0; i < rt.residual.extent(0); ++i) {
		// 计算当前单元格的残差值
		real residual_value = rt.residual(i);
		
		// 更新最大残差值
		max_residual = fmax(max_residual, residual_value);
	}

	return max_residual;
	}

	// 计算对流扩散项
	inline void advection_dispersion(RTState &rt, GwState &gw, GwDomain &gdom, real dt_tmp)
	// inline void advection_dispersion(RTState &rt, GwState &gw, GwDomain &gdom)
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
				//对流项
			    rt.advection(iGlob) = -rt.aveV(iGlob, 0) * (rt.c(iGlob, 1) - rt.c(iGlob - 1, 1)) / gdom.dx - rt.aveV(iGlob, 1) * (rt.c(iGlob, 1) - rt.c(iGlob - gdom.nxhc, 1)) / gdom.dy - rt.aveV(iGlob, 2) * (rt.c(iGlob, 1) - rt.c(iGlob - gdom.nxhc * gdom.nyhc, 1)) / gdom.dz(iGlob);				

				
				// std::cout << "-------rt.aveV(iGlob, 3)----- " <<rt.aveV(iGlob, 3) << std::endl;
				// 机械弥散系数计算
				rt.dcal(iGlob,0) = (0.3*(pow(rt.aveV(iGlob, 1), 2)+pow(rt.aveV(iGlob, 2), 2))+6*pow(rt.aveV(iGlob, 0), 2))/rt.aveV(iGlob,3);
				rt.dcal(iGlob,1) = (0.3*(pow(rt.aveV(iGlob, 0), 2)+pow(rt.aveV(iGlob, 2), 2))+6*pow(rt.aveV(iGlob, 1), 2))/rt.aveV(iGlob,3);
				rt.dcal(iGlob,2) = (0.3*(pow(rt.aveV(iGlob, 0), 2)+pow(rt.aveV(iGlob, 1), 2))+6*pow(rt.aveV(iGlob, 2), 2))/rt.aveV(iGlob,3);
				// std::cout << "-------rt.d(iGlob,0)----- " <<rt.d(iGlob,0) << std::endl;
				//将rt.d中nan值替换
				Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
					for (int i = 0; i < 3; ++i) {
						if (std::isnan(rt.dcal(iGlob,i))) {
							rt.dcal(iGlob,i) = 1e-6;
						}
					}
					});
				//弥散系数=机械弥散+分子弥散
				//1e-5为分子弥散系数
				rt.dcal(iGlob,0) += 1e-3*(pow(gw.wc(iGlob,0), 7.0/3.0) / pow(phi, 2.0));
				rt.dcal(iGlob,1) += 1e-3*(pow(gw.wc(iGlob,0), 7.0/3.0) / pow(phi, 2.0));
				rt.dcal(iGlob,2) += 1e-3*(pow(gw.wc(iGlob,0), 7.0/3.0) / pow(phi, 2.0));
				// std::cout << "-------rt.dcal(iGlob,0)----- " <<rt.dcal(iGlob,0) << std::endl;	

				//弥散项计算
			    rt.dispersion(iGlob,0) = gw.wc(iGlob,0)*rt.dcal(iGlob,0) * (rt.c(iGlob + 1, 1) - 2 * rt.c(iGlob, 1) + rt.c(iGlob - 1, 1)) / (gdom.dx * gdom.dx);
			    rt.dispersion(iGlob,1) =  gw.wc(iGlob,0)*rt.dcal(iGlob,1) * (rt.c(iGlob + gdom.nxhc, 1) - 2 * rt.c(iGlob, 1) + rt.c(iGlob - gdom.nxhc, 1)) / (gdom.dy * gdom.dy);
			    rt.dispersion(iGlob,2) =gw.wc(iGlob,0)* rt.dcal(iGlob,2)* (rt.c(iGlob + gdom.nxhc * gdom.nyhc, 1) - 2 * rt.c(iGlob, 1) + rt.c(iGlob - gdom.nxhc * gdom.nyhc, 1)) / (gdom.dz(iGlob) * gdom.dz(iGlob));
			
			    rt.dispersion(iGlob,3) = rt.dispersion(iGlob,0) + rt.dispersion(iGlob,1) + rt.dispersion(iGlob,2);
				// std::cout << "-------gw.wc(iGlob,0)----- " << gw.wc(iGlob,0) << std::endl;
				// std::cout << "-------gw.wc(iGlob,1)----- " << gw.wc(iGlob,1) << std::endl;
			    // Update concentration1
			    rt.c(iGlob, 1) = rt.c(iGlob, 1) + rt.dt * (rt.advection(iGlob) + rt.dispersion(iGlob,3))/gw.wc(iGlob,0);
			//     std::cout << "-------debugzzb-rt1111111111.dt----- " << rt.dt << std::endl;
			// rt.residual(iGlob) = fabs((rt.c(iGlob, 0) - rt.c(iGlob, 1)) / rt.dt - (rt.advection(iGlob) + rt.dispersion(iGlob,3)));
 			// std::cout << "-------rt.residual(iGlob)----- " << rt.residual(iGlob) << std::endl;
				
			
			//     std::cout << "-------debugzzb-rt.rt.c(iGlob, 1)55555----- " << rt.c(iGlob, 1) << std::endl;
			//     std::cout << "-------debugzzb-rt.rt.c(iGlob, 0)55555----- " << rt.c(iGlob, 0) << std::endl;
		    		// Update concentration2
				// rt.residual(iGlob) = fabs((rt.c(iGlob, 0) - rt.c(iGlob, 1)) / rt.dt - (rt.advection(iGlob,0) + rt.dispersion(iGlob,0)));
				
				// rt.c(iGlob, 1) -= rt.residual(iGlob) / ((rt.c(iGlob, 0) - rt.c(iGlob, 1)) / rt.dt - (gw.q(iGlob, 0) + gw.q(iGlob, 1) + gw.q(iGlob, 2)) * rt.c(iGlob, 1) - (rt.dispersion_x(iGlob,0) +rt. dispersion_y(iGlob,0) + rt.dispersion_z(iGlob,0))) * rt.dt;
		    });	
	}

	/*--------检查数值解的收敛性方法二------------*/
	bool check_convergence(RTState &rt, real tolerance)
	{
		
		bool converged = true;
		Kokkos::parallel_reduce(
		    "CheckConvergence", rt.c.extent(0), KOKKOS_LAMBDA(int iGlob, bool &converged) {
			    // Compute concentration change

			    rt.residual(iGlob) = fabs((rt.c(iGlob, 0) - rt.c(iGlob, 1)) / rt.dt - (rt.advection(iGlob) + rt.dispersion(iGlob,3)));
			//     rt.residual(iGlob) = fabs((rt.c(iGlob, 0) - rt.c(iGlob, 1)));
			//     std::cout << "-------debugzzb-rt.residual(iGlob)----- " << rt.residual(iGlob) << std::endl;			    
				// std::cout << "-------debugzzb- tolerance----- " <<  tolerance << std::endl;			    
			    // Check convergence
			    if (rt.residual(iGlob) > tolerance)
			    {
				    converged = true;
				   
			    }
			    
		    },
		    Kokkos::LAnd<bool>(converged)
		    );
		
		
		return converged;
					    
	}


	    /* --------------------------------------------------
        Get RT matrix coefficients
    
    inline void RTlinear_system(GwState &gw, GwDomain &gdom, std::vector<GwBC> &gbc, std::vector<GwSS> &gss, RTMatrix &rtA, Parallel &par, RTState &rt,)	{
        // Calculate matrix coefficients
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, ivg, iGlob, iGlobSW;
            gdom.unpackIndices(idom, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            iGlobSW = (hc+jj)*gdom.nxhc + ii + hc;
		  
		  rt.RTcoef(idom, 0) = 1 + 2 * (rt.dcal(iGlob,0) * rt.dt / (gdom.dx * gdom.dx) + rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)));;
		  rt.RTcoef(idom, 1) = -rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ) - rt.aveV(iGlob, 0) * rt.dt / (2 * gdom.dx );
            rt.RTcoef(idom, 2) = -rt.dcal(iGlob,0) * rt.dt / (gdom.dx  * gdom.dx ) + rt.aveV(iGlob, 0) * rt.dt / (2 * gdom.dx );
            rt.RTcoef(idom, 3) = -rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) - rt.aveV(iGlob, 1) * rt.dt / (2 * gdom.dy);
            rt.RTcoef(idom, 4) = -rt.dcal(iGlob,1) * rt.dt / (gdom.dy * gdom.dy) + rt.aveV(iGlob, 1) * rt.dt / (2 * gdom.dy);
            rt.RTcoef(idom, 5) = -rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)) - rt.aveV(iGlob, 2) * rt.dt / (2 * gdom.dz(iGlob));
            rt.RTcoef(idom, 6) = -rt.dcal(iGlob,2) * rt.dt / (gdom.dz(iGlob) * gdom.dz(iGlob)) + rt.aveV(iGlob, 2) * rt.dt / (2 * gdom.dz(iGlob));
		  rt.RTcoef(idom, 7) = -rt.c(iGlob, 1);
		 }); 
		  // Insert coefficients into Matrix rtA
		   Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, irow = rtA.ptr(idom);
			gdom.unpackIndices(idom, kk, jj, ii);
		     if (kk > 0)	{rtA.ind(irow) = idom - gdom.nx*gdom.ny;	rtA.val(irow) = gw.coef(idom,6);  irow++;}
        		if (jj > 0)	{rtA.ind(irow) = idom - gdom.nx;	        rtA.val(irow) = gw.coef(idom,4);  irow++;}
			if (ii > 0)	{rtA.ind(irow) = idom - 1;		        rtA.val(irow) = gw.coef(idom,2);  irow++;}
			rtA.ind(irow) = idom;	rtA.val(irow) = gw.coef(idom,0);	irow++;
			if (ii < gdom.nx-1)	{rtA.ind(irow) = idom + 1;		        rtA.val(irow) = gw.coef(idom,1);  irow++;}
			if (jj < gdom.ny-1)	{rtA.ind(irow) = idom + gdom.nx;	        rtA.val(irow) = gw.coef(idom,3);  irow++;}
			if (kk < gdom.nz-1)	{rtA.ind(irow) = idom + gdom.nx*gdom.ny;	rtA.val(irow) = gw.coef(idom,5);  irow++;}
		
			rtA.rhs(idom) = rt.RTcoef(idom,7);



		   });

-------------------------------------------------- */
		  
						
};

#endif
