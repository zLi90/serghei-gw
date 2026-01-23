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
#include "RTState.h"
#include "State.h"
#include "RTIntegrator.h"
#include <set>
#include <limits>
#include <math.h>
#include <cmath>

class RTFunction
{

private:
	Kokkos::Timer timer;
	Kokkos::Timer timer2;

	/* --------------------------------------------------
		Minmod 通量限制器
		ψ(r) = max(0, min(1, r))
	-------------------------------------------------- */
	KOKKOS_INLINE_FUNCTION
	real minmod_limiter(real r) const
	{
		if (r <= 0.0)
			return 0.0;
		if (r >= 1.0)
			return 1.0;
		return r;
	}

	/* --------------------------------------------------
		Superbee 通量限制器
		ψ(r) = max(0, min(1, 2r), min(2, r))
	-------------------------------------------------- */
	KOKKOS_INLINE_FUNCTION
	real superbee_limiter(real r) const
	{
		if (r <= 0.0)
			return 0.0;
		real val1 = (r <= 0.5) ? 2.0 * r : 1.0;
		real val2 = (r <= 2.0) ? r : 2.0;
		return fmax(0.0, fmin(1.0, val1));
	}

	/* --------------------------------------------------
		Van Leer 通量限制器
		ψ(r) = (r + |r|) / (1 + |r|)
	-------------------------------------------------- */
	KOKKOS_INLINE_FUNCTION
	real vanleer_limiter(real r) const
	{
		if (r <= 0.0)
			return 0.0;
		return (r + fabs(r)) / (1.0 + fabs(r));
	}

	/* --------------------------------------------------
		计算梯度比 r (用于通量限制器)
		r = (C_upstream - C_upstream_upstream) / (C_downstream - C_upstream)
	-------------------------------------------------- */
	KOKKOS_INLINE_FUNCTION
	real compute_gradient_ratio(real C_upstream, real C_downstream, real C_upstream_upstream) const
	{
		real denominator = C_downstream - C_upstream;
		if (fabs(denominator) < 1e-12)
			return 0.0;
		return (C_upstream - C_upstream_upstream) / denominator;
	}

	/* --------------------------------------------------
		计算高阶通量修正项 (反扩散通量)
		Flux_correction = ψ(r) * (Flux_high - Flux_low)
	-------------------------------------------------- */
	KOKKOS_INLINE_FUNCTION
	real compute_flux_correction(real C_upstream, real C_downstream, real C_upstream_upstream,
								 real velocity, real dt, real dx, int limiter_type) const
	{
		// 一阶迎风通量
		real flux_low = velocity * C_upstream;

		// 二阶中心通量 (高阶)
		real flux_high = velocity * 0.5 * (C_upstream + C_downstream);

		// 计算梯度比
		real r = compute_gradient_ratio(C_upstream, C_downstream, C_upstream_upstream);

		// 应用通量限制器
		real psi = 0.0;
		switch (limiter_type)
		{
		case 1: // Minmod
			psi = minmod_limiter(r);
			break;
		case 2: // Superbee
			psi = superbee_limiter(r);
			break;
		case 3: // Van Leer
			psi = vanleer_limiter(r);
			break;
		default:
			psi = minmod_limiter(r);
			break;
		}

		// 返回通量修正项 (反扩散通量)
		return psi * (flux_high - flux_low);
	}

public:
	/* --------------------------------------------------
		rt_solve: 主求解函数
	-------------------------------------------------- */
	template <typename execution_space, typename type_solver>
	inline void rt_solve(RTState &rt, RTMatrix &rtA, GwState &gw, GwDomain &gdom, std::vector<RTBC> &rtgbc, GwMPI &gmpi, Parallel &par, RTIntegrator &rtint, type_solver &rtsolver)
	{
		int iter, ierr = 1;
		real dt_tmp;
		int picard_iter_max = 3; // Picard迭代次数，建议3-5次
		real picard_tol = 1e-6;	 // Picard迭代收敛容差

		rt.dt = gdom.dt;

		// 1. 施加边界条件 (浓度)
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
		}

		// 确保在计算通量前，Ghost Cell的浓度和水动力变量是最新的
		gmpi.mpi_sendrecv(rt.c, gdom, par);

		// 2. 计算弥散张量和物理速度
		// 注意：此时会处理坐标系转换
		dispersion_tensor(rt, gw, gdom, gmpi, par);

		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			rt.c(iGlob, 0) = rt.c(iGlob, 1); // 备份 Cn 到 index 0
		});

		// 3. Picard迭代循环 (用于高阶通量修正)
		for (int picard_iter = 0; picard_iter < picard_iter_max; picard_iter++)
		{
			// 3.1 组装线性系统矩阵 (LHS保持一阶迎风，RHS添加高阶修正)
			RTlinear_system(gw, gdom, rtA, par, rt, rtgbc, picard_iter);

			// 3.2 选择求解方法
			if (rt.rt_scheme == 1)
			{
				rtsolver.Gauss_Seidel(rtA);
			}
			else if (rt.rt_scheme == 2)
			{
				// rtsolver.Jacobi_BICGSTAB_Solve(rtA);
				rtsolver.Jacobi_GMRES_Solve(rtA);
			}

			// 3.3 更新浓度
			Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
				int ii, jj, kk, iGlob;
				gdom.unpackIndices(idom, kk, jj, ii);
				iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;

				// 更新当前时刻浓度 (0为n+1时刻，1为n时刻)
				rt.c(iGlob, 1) = rtA.rt_x(idom); 

				// 物理约束：浓度非负
				if (rt.c(iGlob, 1) < 0.0) rt.c(iGlob, 1) = 0.0; });

			// 3.4 通信同步更新后的浓度
			gmpi.mpi_sendrecv(rt.c, gdom, par);

			// 3.5 施加边界条件以确保边界值正确
			for (int k = 0; k < rtgbc.size(); k++)
			{
				rtgbc[k].applyConcentrationBC(rt, gw, gdom, par);
			}
			gmpi.mpi_sendrecv(rt.c, gdom, par);

			// 3.6 检查Picard迭代收敛 (可选)
			// 如果需要收敛性检查，可以在这里添加
		}

		// 4. 更新固相浓度 (用于下一次迭代或输出)
		// 注意：这里仅做显式更新用于记录，矩阵求解中已经隐式包含
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;

			if (rt.AdsorptionDesorptionModel == Equilibrium_Freundlich_Model)
			{
				rt.c_solid(iGlob, 0) = rt.Kf * pow(rt.c(iGlob, 1), rt.Nf);
				rt.c_solid(iGlob, 1) = rt.c_solid(iGlob, 0);
			}
			else if (rt.AdsorptionDesorptionModel == Equilibrium_Linear_Model)
			{
				rt.c_solid(iGlob, 0) = rt.Kd * rt.c(iGlob, 1);
				rt.c_solid(iGlob, 1) = rt.c_solid(iGlob, 0);
			}
			else if (rt.AdsorptionDesorptionModel == Nonequilibrium_Model)
			{
				// 非平衡吸附需要单独的时间积分
				rt.c_solid(iGlob, 0) = (rt.beta * rt.c(iGlob, 1) + rt.rho_b * rt.c_solid(iGlob, 1) / rt.dt) / (rt.rho_b / rt.dt + rt.beta / rt.Kd + rt.lambda_2 * rt.rho_b);
				rt.c_solid(iGlob, 1) = rt.c_solid(iGlob, 0);
			} });

		// 5. 时间步长控制与同步
		dt_tmp = gdom.dt;
		ierr = MPI_Allreduce(&dt_tmp, &gdom.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		// 准备下一时间步：将当前新值作为旧值 (如果非迭代模式)
		// 注意：通常在主时间循环外进行 swap，这里保留原逻辑
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { rt.c(iGlob, 0) = rt.c(iGlob, 1); });

		rtint.integrate(gw, rt, gdom, rtgbc);
	}

	/* --------------------------------------------------
		dispersion_tensor: 计算物理速度与水动力弥散系数张量
	-------------------------------------------------- */
	inline void dispersion_tensor(RTState &rt, GwState &gw, GwDomain &gdom, GwMPI &gmpi, Parallel &par)
	{
		// 1. 速度矢量转换 (Coordinate Transformation)
		// 原始定义：q(0)向左为正，q(1)向内为正，q(2)向上为正
		// 网格定义：i向右增，j向外增，k向下增
		// 因此：v_x = -q(0), v_y = -q(1), v_z = -q(2)
		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			rt.aveVB(iGlob, 0) = -(gw.q_new(iGlob, 0)); // vx
			rt.aveVB(iGlob, 1) = -(gw.q_new(iGlob, 1)); // vy
			rt.aveVB(iGlob, 2) = -(gw.q_new(iGlob, 2)); // vz
			
			// 计算速度模长 |v|
			rt.aveVB(iGlob, 3) = sqrt(pow(rt.aveVB(iGlob, 0), 2) + 
									  pow(rt.aveVB(iGlob, 1), 2) + 
									  pow(rt.aveVB(iGlob, 2), 2)); });

		// 2. 计算界面处的水动力弥散系数 (i+1/2, j+1/2, k+1/2)
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
			// qx,i+1/2,j,k
			rt.q_dispersion(iGlob, 0) = rt.aveVB(iGlob, 0);
			// qy,i+1/2,j,k
			rt.q_dispersion(iGlob, 1) = 0.5 * ((rt.aveVB(iGlob - gdom.nxhc, 1) + rt.aveVB(iGlob, 1)) * 0.5 + (rt.aveVB(iGlob + 1 - gdom.nxhc, 1) + rt.aveVB(iGlob + 1, 1)) * 0.5);
			// qz,i+1/2,j,k
			rt.q_dispersion(iGlob, 2) = 0.5 * ((rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) + rt.aveVB(iGlob, 2)) * 0.5 + (rt.aveVB(iGlob + 1 - gdom.nxhc * gdom.nyhc, 2) + rt.aveVB(iGlob + 1, 2)) * 0.5);
			// q,i+1/2,j,k
			rt.q_dispersion(iGlob, 3) = sqrt(pow(rt.q_dispersion(iGlob, 0), 2) + pow(rt.q_dispersion(iGlob, 1), 2) + pow(rt.q_dispersion(iGlob, 2), 2));

			// qx,i,j+1/2,k
			//! rt.q_dispersion(iGlob, 4) = 0.5 * ((rt.aveVB(iGlob, 1) + rt.aveVB(iGlob - 1, 1)) * 0.5 + (rt.aveVB(iGlob + gdom.nxhc, 1) + rt.aveVB(iGlob + gdom.nxhc - 1, 1)) * 0.5);
			rt.q_dispersion(iGlob, 4) = 0.5 * ((rt.aveVB(iGlob, 0) + rt.aveVB(iGlob - 1, 0)) * 0.5 + (rt.aveVB(iGlob + gdom.nxhc, 0) + rt.aveVB(iGlob + gdom.nxhc - 1, 0)) * 0.5);
			// qy,i,j+1/2,k
			rt.q_dispersion(iGlob, 5) = rt.aveVB(iGlob, 1);
			// qz,i,j+1/2,k
			rt.q_dispersion(iGlob, 6) = 0.5 * ((rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2) + rt.aveVB(iGlob, 2)) * 0.5 + (rt.aveVB(iGlob + gdom.nxhc - gdom.nxhc * gdom.nyhc, 2) + rt.aveVB(iGlob + gdom.nxhc, 2)) * 0.5);
			// q,i,j+1/2,k
			rt.q_dispersion(iGlob, 7) = sqrt(pow(rt.q_dispersion(iGlob, 4), 2) + pow(rt.q_dispersion(iGlob, 5), 2) + pow(rt.q_dispersion(iGlob, 6), 2));

			// qx,i,j,k+1/2
			rt.q_dispersion(iGlob, 8) = 0.5 * ((rt.aveVB(iGlob, 0) + rt.aveVB(iGlob + gdom.nxhc * gdom.nyhc, 0)) * 0.5 + (rt.aveVB(iGlob - 1, 0) + rt.aveVB(iGlob + gdom.nxhc * gdom.nyhc - 1, 0)) * 0.5);
			// qy,i,j,k+1/2
			rt.q_dispersion(iGlob, 9) = 0.5 * ((rt.aveVB(iGlob, 1) + rt.aveVB(iGlob - gdom.nxhc, 1)) * 0.5 + (rt.aveVB(iGlob - gdom.nxhc + gdom.nxhc * gdom.nyhc, 1) + rt.aveVB(iGlob + gdom.nxhc * gdom.nyhc, 1)) * 0.5);
			// qz,i,j,k+1/2
			rt.q_dispersion(iGlob, 10) = rt.aveVB(iGlob, 2);
			// q,i,j,k+1/2
			rt.q_dispersion(iGlob, 11) = sqrt(pow(rt.q_dispersion(iGlob, 8), 2) + pow(rt.q_dispersion(iGlob, 9), 2) + pow(rt.q_dispersion(iGlob, 10), 2)); });

		// 3. 计算弥散张量分量
		// 并添加分子扩散 (使用界面平均含水率)

		Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {		
		
		// ! DXX
		// Dxxi+1/2,j,k
		rt.dcal(iGlob, 0) = (rt.alpha_T * (pow(rt.q_dispersion(iGlob, 1), 2) + pow(rt.q_dispersion(iGlob, 2), 2)) + rt.alpha_L * pow(rt.q_dispersion(iGlob, 0), 2)) / rt.q_dispersion(iGlob, 3);

		// Dxyi+1/2,j,k
		rt.dcal(iGlob, 3) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 0) * rt.q_dispersion(iGlob, 1) / rt.q_dispersion(iGlob, 3);
		// Dxzi+1/2,j,k
		rt.dcal(iGlob, 4) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 0) * rt.q_dispersion(iGlob, 2) / rt.q_dispersion(iGlob, 3);

		//! DYY
		// DYYi,j+1/2,k
		rt.dcal(iGlob, 1) = (rt.alpha_T * (pow(rt.q_dispersion(iGlob, 4), 2) + pow(rt.q_dispersion(iGlob, 6), 2)) + rt.alpha_L * pow(rt.q_dispersion(iGlob, 5), 2)) / rt.q_dispersion(iGlob, 7);
		// Dyx,i,j+1/2,k
		rt.dcal(iGlob, 5) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 5) * rt.q_dispersion(iGlob, 4) / rt.q_dispersion(iGlob, 7);
		// Dyzi,j+1/2,k
		rt.dcal(iGlob, 6) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 5) * rt.q_dispersion(iGlob, 6) / rt.q_dispersion(iGlob, 7);
		//! DZZ
		// Dzz,i,j,k+1/2
		rt.dcal(iGlob, 2) = (rt.alpha_T * (pow(rt.q_dispersion(iGlob, 8), 2) + pow(rt.q_dispersion(iGlob, 9), 2)) + rt.alpha_L * pow(rt.q_dispersion(iGlob, 10), 2)) / rt.q_dispersion(iGlob, 11);
		// Dzx,i,j,k+1/2
		rt.dcal(iGlob, 7) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 10) * rt.q_dispersion(iGlob, 8) / rt.q_dispersion(iGlob, 11);
		// Dzy,i,j,k+1/2
		rt.dcal(iGlob, 8) = (rt.alpha_L - rt.alpha_T) * rt.q_dispersion(iGlob, 10) * rt.q_dispersion(iGlob, 9) / rt.q_dispersion(iGlob, 11); 
		
	
		// !若处于边界，c_difxy等系数应在dispersion_tensor中被设为0，从而自然处理边界
		//! 如果网格qz速度为0，则Dzz=0，即当底边界为水动力无通量边界时
		// todo 后续考虑是不是要结合边界条件对弥散系数张量进行修改
		
		if (rt.q_dispersion(iGlob, 10) == 0) // qz,i,j,k+1/2
		{
			rt.dcal(iGlob, 2) = 0;
		}

		if (rt.q_dispersion(iGlob, 5) == 0) // qy,i,j+1/2,k
		{
			rt.dcal(iGlob, 1) = 0;
		}



		// 将rt.d中nan值替换
		for (int i = 0; i < 9; ++i)
		{
			if (std::isnan(rt.dcal(iGlob, i)))
			{
				rt.dcal(iGlob, i) = 0;
			}
		} 
				
		rt.dcal(iGlob, 0) += gw.wc(iGlob, 0) * rt.diffusion_molecular;
		rt.dcal(iGlob, 1) += gw.wc(iGlob, 0) * rt.diffusion_molecular;
		rt.dcal(iGlob, 2) += gw.wc(iGlob, 0) * rt.diffusion_molecular; 
		
		// Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {	
		// 根据研究区域维数对水动力弥散系数进行调整
		if (gdom.nx == 1) // yz纬度，qx=0，Dxx=0,Dxy=DYX=DXZ=DZX=0
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
		} });
	}

	/* --------------------------------------------------
	   RTlinear_system: 组装线性系统 (质量守恒修正版 + 高阶通量修正)
	-------------------------------------------------- */
	inline void RTlinear_system(GwState &gw, GwDomain &gdom, RTMatrix &rtA, Parallel &par, RTState &rt, std::vector<RTBC> &rtgbc, int picard_iter = 0)
	{
		// 1. 重要：初始化系数矩阵，防止累加错误
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			for (int k = 0; k < 27; ++k) rt.RTcoef(idom, k) = 0.0; });

		// 2. 计算矩阵系数
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;

			// --- 迎风格式权重 (基于校正后的正向速度) ---
			// v > 0 (向右): weight=1.0 (取Ci); v < 0 (向左): weight=0.0 (取Ci+1)
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
			
			// 垂直方向几何权重
			rt.wz(iGlob) = 0.5;

			// --- 通量系数计算 (Flux Coefficients) ---
			// 对流项系数: v * dt / dx (注意 aveVB 已经修正为沿坐标轴正向)
			rt.c_advxx(iGlob) = rt.aveVB(iGlob, 0) * rt.dt / gdom.dx;
			rt.c_advyy(iGlob) = rt.aveVB(iGlob, 1) * rt.dt / gdom.dy;
			rt.c_advzz(iGlob) = rt.aveVB(iGlob, 2) * rt.dt; // 后面再除dz

			// 弥散项系数: D * dt / dx^2
			rt.c_difxx(iGlob) = rt.dcal(iGlob, 0) * rt.dt / (gdom.dx * gdom.dx);
			rt.c_difyy(iGlob) = rt.dcal(iGlob, 1) * rt.dt / (gdom.dy * gdom.dy);
			rt.c_difzz(iGlob) = rt.dcal(iGlob, 2) * rt.dt / (0.5 * gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc));
			
			// 交叉项系数
			rt.c_difxy(iGlob) = rt.dcal(iGlob, 3) * rt.dt / (2 * gdom.dx * gdom.dy);
			rt.c_difxz(iGlob) = rt.dcal(iGlob, 4) * rt.dt; 
			rt.c_difyx(iGlob) = rt.dcal(iGlob, 5) * rt.dt / (2 * gdom.dx * gdom.dy);
			rt.c_difyz(iGlob) = rt.dcal(iGlob, 6) * rt.dt;		
			rt.c_difzx(iGlob) = rt.dcal(iGlob, 7) * rt.dt;
			rt.c_difzy(iGlob) = rt.dcal(iGlob, 8) * rt.dt;

			// --- 质量守恒累积项 (Mass Accumulation) ---
			// 方程: d(theta * C + rho * S)/dt = ...
			// 离散: (theta_new * C_new - theta_old * C_old)/dt + (rho * S_new - rho * S_old)/dt
			
			real mass_lhs_diag = gw.wc_new(iGlob);     // theta_new
			real mass_rhs      = gw.wc_old(iGlob) * rt.c(iGlob, 0); // theta_old * C_old

			// 吸附反应处理
			if (rt.ReactionModule == 1)
			{
				if (rt.AdsorptionDesorptionModel == Equilibrium_Linear_Model)
				{
					// S = Kd * C -> d(rho*Kd*C)/dt
					mass_lhs_diag += rt.rho_b * rt.Kd;
					mass_rhs      += rt.rho_b * rt.Kd * rt.c(iGlob, 0);
				}
				else if (rt.AdsorptionDesorptionModel == Equilibrium_Freundlich_Model)
				{
					// S = Kf * C^Nf -> 线性化: S_new ~ S_old + dS/dC * (C_new - C_old)
					real C_curr = rt.c(iGlob, 1);
					if(C_curr < 1e-9) C_curr = 1e-9;
					
					real capacity = rt.Kf * rt.Nf * pow(C_curr, rt.Nf - 1.0);
					mass_lhs_diag += rt.rho_b * capacity;
					// RHS 包含 S_old 的贡献
					mass_rhs += rt.rho_b * (rt.Kf * pow(C_curr, rt.Nf)) - rt.rho_b * capacity * C_curr + rt.rho_b * capacity * rt.c(iGlob, 0);
					// 简化为 C_new * Capacity = C_old * Capacity + (S_old - S(C_old))? 
					// 标准Picard: LHS = Capacity, RHS = Capacity * C_old + Residual. 
					// 此处简化为直接使用 Capacity 修正质量矩阵
				}
				else if (rt.AdsorptionDesorptionModel == Nonequilibrium_Model)
				{
					// 复杂的非平衡源汇项，此处保持原有参数逻辑但修正位置
					double intermediate = rt.beta / (rt.Kd * (rt.rho_b / rt.dt + rt.beta / rt.Kd + rt.lambda_2 * rt.rho_b));
					// ... 需要根据具体方程推导，此处暂略，建议仔细核对非平衡方程离散格式
				}
			}

			// --- 组装主对角线及相邻项 ---
			// 1. 对角线 (i,j,k)
			rt.RTcoef(idom, 0) = mass_lhs_diag 
				+ rt.c_difxx(iGlob) + rt.c_difxx(iGlob - 1)
				+ rt.Up_Weighting_xp(iGlob) * rt.c_advxx(iGlob) - (1.0 - rt.Up_Weighting_xm(iGlob)) * rt.c_advxx(iGlob - 1)
				+ rt.c_difyy(iGlob) + rt.c_difyy(iGlob - gdom.nxhc)
				+ rt.Up_Weighting_yp(iGlob) * rt.c_advyy(iGlob) - (1.0 - rt.Up_Weighting_ym(iGlob)) * rt.c_advyy(iGlob - gdom.nxhc)
				+ rt.c_difzz(iGlob) / gdom.dz(iGlob) + rt.c_difzz(iGlob - gdom.nxhc * gdom.nyhc) / gdom.dz(iGlob)
				+ rt.Up_Weighting_zp(iGlob) * rt.c_advzz(iGlob) / gdom.dz(iGlob) - (1.0 - rt.Up_Weighting_zm(iGlob)) * rt.c_advzz(iGlob - gdom.nxhc * gdom.nyhc) / gdom.dz(iGlob);

			// 添加衰减项 lambda * theta * C
			if(rt.ReactionModule == 1) {
				rt.RTcoef(idom, 0) += rt.lambda * gw.wc_new(iGlob) * rt.dt; // 若方程未同除dt，则需注意量纲
			}

			// 2. 邻居项 (移项到LHS，注意符号)
			// Ci+1
			rt.RTcoef(idom, 1) = -(rt.c_difxx(iGlob) - (1.0 - rt.Up_Weighting_xp(iGlob)) * rt.c_advxx(iGlob));
			// Ci-1
			rt.RTcoef(idom, 2) = -(rt.c_difxx(iGlob - 1) + rt.Up_Weighting_xm(iGlob) * rt.c_advxx(iGlob - 1));
			// Cj+1
			rt.RTcoef(idom, 3) = -(rt.c_difyy(iGlob) - (1.0 - rt.Up_Weighting_yp(iGlob)) * rt.c_advyy(iGlob));
			// Cj-1
			rt.RTcoef(idom, 4) = -(rt.c_difyy(iGlob - gdom.nxhc) + rt.Up_Weighting_ym(iGlob) * rt.c_advyy(iGlob - gdom.nxhc));
			// Ck+1
			rt.RTcoef(idom, 5) = -(rt.c_difzz(iGlob)/gdom.dz(iGlob) - (1.0 - rt.Up_Weighting_zp(iGlob)) * rt.c_advzz(iGlob)/gdom.dz(iGlob));
			// Ck-1
			rt.RTcoef(idom, 6) = -(rt.c_difzz(iGlob - gdom.nxhc * gdom.nyhc)/gdom.dz(iGlob) + rt.Up_Weighting_zm(iGlob) * rt.c_advzz(iGlob - gdom.nxhc * gdom.nyhc)/gdom.dz(iGlob));

			// 3. RHS (包含上一时刻质量 + 源汇项)
			rt.RTcoef(idom, 7) = mass_rhs;

			// 4. 高阶通量修正项 (仅在Picard迭代中添加)
			// 策略：LHS保持一阶迎风，RHS显式添加高阶修正
			if (picard_iter > 0)
			{
				// 通量限制器类型：1=Minmod, 2=Superbee, 3=Van Leer
				int limiter_type = 2;  // 默认使用Superbee限制器

				// X方向高阶通量修正
				// i+1/2 界面通量修正
				if (ii < gdom.nx - 1)
				{
					int iGlob_plus = iGlob + 1;
					int iGlob_minus = iGlob - 1;
					real C_upstream, C_downstream, C_upstream_upstream;
					
					// 根据速度方向确定上下游
					if (rt.aveVB(iGlob, 0) >= 0)  // v > 0: 上游是i，下游是i+1
					{
						C_upstream = rt.c(iGlob, 1);
						C_downstream = rt.c(iGlob_plus, 1);
						C_upstream_upstream = (ii > 0) ? rt.c(iGlob_minus, 1) : rt.c(iGlob, 1);
					}
					else  // v < 0: 上游是i+1，下游是i
					{
						C_upstream = rt.c(iGlob_plus, 1);
						C_downstream = rt.c(iGlob, 1);
						C_upstream_upstream = (ii < gdom.nx - 2) ? rt.c(iGlob_plus + 1, 1) : rt.c(iGlob_plus, 1);
					}
					
					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
					                                           rt.aveVB(iGlob, 0), rt.dt, gdom.dx, limiter_type);
					// 添加到RHS：通量修正项的散度 (dt/dx * (F_corr_i+1/2 - F_corr_i-1/2))
					rt.RTcoef(idom, 7) -= rt.dt / gdom.dx * flux_corr;
				}

				// i-1/2 界面通量修正
				if (ii > 0)
				{
					int iGlob_minus = iGlob - 1;
					int iGlob_minus2 = iGlob - 2;
					real C_upstream, C_downstream, C_upstream_upstream;
					
					// 根据速度方向确定上下游
					if (rt.aveVB(iGlob_minus, 0) >= 0)  // v > 0: 上游是i-1，下游是i
					{
						C_upstream = rt.c(iGlob_minus, 1);
						C_downstream = rt.c(iGlob, 1);
						C_upstream_upstream = (ii > 1) ? rt.c(iGlob_minus2, 1) : rt.c(iGlob_minus, 1);
					}
					else  // v < 0: 上游是i，下游是i-1
					{
						C_upstream = rt.c(iGlob, 1);
						C_downstream = rt.c(iGlob_minus, 1);
						C_upstream_upstream = (ii < gdom.nx - 1) ? rt.c(iGlob + 1, 1) : rt.c(iGlob, 1);
					}
					
					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
					                                           rt.aveVB(iGlob_minus, 0), rt.dt, gdom.dx, limiter_type);
					// 从RHS减去：通量修正项的散度 (dt/dx * (F_corr_i+1/2 - F_corr_i-1/2))
					rt.RTcoef(idom, 7) += rt.dt / gdom.dx * flux_corr;
				}

				// Y方向高阶通量修正
				if (jj < gdom.ny - 1)
				{
					int iGlob_plus = iGlob + gdom.nxhc;
					int iGlob_minus = iGlob - gdom.nxhc;
					real C_upstream, C_downstream, C_upstream_upstream;
					
					if (rt.aveVB(iGlob, 1) >= 0)
					{
						C_upstream = rt.c(iGlob, 1);
						C_downstream = rt.c(iGlob_plus, 1);
						C_upstream_upstream = (jj > 0) ? rt.c(iGlob_minus, 1) : rt.c(iGlob, 1);
					}
					else
					{
						C_upstream = rt.c(iGlob_plus, 1);
						C_downstream = rt.c(iGlob, 1);
						C_upstream_upstream = (jj < gdom.ny - 2) ? rt.c(iGlob_plus + gdom.nxhc, 1) : rt.c(iGlob_plus, 1);
					}
					
					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
					                                           rt.aveVB(iGlob, 1), rt.dt, gdom.dy, limiter_type);
					rt.RTcoef(idom, 7) -= rt.dt / gdom.dy * flux_corr;
				}

				if (jj > 0)
				{
					int iGlob_minus = iGlob - gdom.nxhc;
					int iGlob_minus2 = iGlob - 2 * gdom.nxhc;
					real C_upstream, C_downstream, C_upstream_upstream;
					
					if (rt.aveVB(iGlob_minus, 1) >= 0)
					{
						C_upstream = rt.c(iGlob_minus, 1);
						C_downstream = rt.c(iGlob, 1);
						C_upstream_upstream = (jj > 1) ? rt.c(iGlob_minus2, 1) : rt.c(iGlob_minus, 1);
					}
					else
					{
						C_upstream = rt.c(iGlob, 1);
						C_downstream = rt.c(iGlob_minus, 1);
						C_upstream_upstream = (jj < gdom.ny - 1) ? rt.c(iGlob + gdom.nxhc, 1) : rt.c(iGlob, 1);
					}
					
					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
					                                           rt.aveVB(iGlob_minus, 1), rt.dt, gdom.dy, limiter_type);
					rt.RTcoef(idom, 7) += rt.dt / gdom.dy * flux_corr;
				}

				// Z方向高阶通量修正
				if (kk < gdom.nz - 1)
				{
					int iGlob_plus = iGlob + gdom.nxhc * gdom.nyhc;
					int iGlob_minus = iGlob - gdom.nxhc * gdom.nyhc;
					real C_upstream, C_downstream, C_upstream_upstream;
					
					if (rt.aveVB(iGlob, 2) >= 0)
					{
						C_upstream = rt.c(iGlob, 1);
						C_downstream = rt.c(iGlob_plus, 1);
						C_upstream_upstream = (kk > 0) ? rt.c(iGlob_minus, 1) : rt.c(iGlob, 1);
					}
					else
					{
						C_upstream = rt.c(iGlob_plus, 1);
						C_downstream = rt.c(iGlob, 1);
						C_upstream_upstream = (kk < gdom.nz - 2) ? rt.c(iGlob_plus + gdom.nxhc * gdom.nyhc, 1) : rt.c(iGlob_plus, 1);
					}
					
					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
					                                           rt.aveVB(iGlob, 2), rt.dt, gdom.dz(iGlob), limiter_type);
					rt.RTcoef(idom, 7) -= rt.dt / gdom.dz(iGlob) * flux_corr;
				}

				if (kk > 0)
				{
					int iGlob_minus = iGlob - gdom.nxhc * gdom.nyhc;
					int iGlob_minus2 = iGlob - 2 * gdom.nxhc * gdom.nyhc;
					real C_upstream, C_downstream, C_upstream_upstream;
					
					if (rt.aveVB(iGlob_minus, 2) >= 0)
					{
						C_upstream = rt.c(iGlob_minus, 1);
						C_downstream = rt.c(iGlob, 1);
						C_upstream_upstream = (kk > 1) ? rt.c(iGlob_minus2, 1) : rt.c(iGlob_minus, 1);
					}
					else
					{
						C_upstream = rt.c(iGlob, 1);
						C_downstream = rt.c(iGlob_minus, 1);
						C_upstream_upstream = (kk < gdom.nz - 1) ? rt.c(iGlob + gdom.nxhc * gdom.nyhc, 1) : rt.c(iGlob, 1);
					}
					
					real flux_corr = compute_flux_correction(C_upstream, C_downstream, C_upstream_upstream,
					                                           rt.aveVB(iGlob_minus, 2), rt.dt, gdom.dz(iGlob), limiter_type);
					rt.RTcoef(idom, 7) += rt.dt / gdom.dz(iGlob) * flux_corr;
				}
			} });

		// 3. 交叉弥散项 (Cross-Dispersion Terms) - 保持原有结构但修正索引
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk, iGlob;
			gdom.unpackIndices(idom, kk, jj, ii);
			iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;

			// 修正后的 Cj+1 项累加
			rt.RTcoef(idom, 8) = -0.5 * rt.c_difxy(iGlob) + 0.5 * rt.c_difxy(iGlob - 1); // 简化示例
			rt.RTcoef(idom, 3) += rt.RTcoef(idom, 8);
			// ci,j+1//添加到上面
			rt.RTcoef(idom, 8) = -0.5 * rt.c_difxy(iGlob) + 0.5 * rt.c_difxy(iGlob - 1) - rt.wz(iGlob) * rt.c_difzy(iGlob) / (2 * gdom.dy * gdom.dz(iGlob)) + (1 - rt.wz(iGlob - gdom.nxhc * gdom.nyhc)) * rt.c_difzy(iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dy * gdom.dz(iGlob)); // Ci,j+1
			rt.RTcoef(idom, 3) += rt.RTcoef(idom, 8);
			// Ci+1,j+1
			rt.RTcoef(idom, 9) = -0.5 * rt.c_difxy(iGlob) - 0.5 * rt.c_difyx(iGlob); // Ci+1,j+1
			// Cj-1//添加到上面
			rt.RTcoef(idom, 10) = 0.5 * rt.c_difxy(iGlob) - 0.5 * rt.c_difxy(iGlob - 1) + rt.wz(iGlob) * rt.c_difzy(iGlob) / (2 * gdom.dy * gdom.dz(iGlob)) - (1 - rt.wz(iGlob - gdom.nxhc * gdom.nyhc)) * rt.c_difzy(iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dy * gdom.dz(iGlob)); // Cj-1
			rt.RTcoef(idom, 4) += rt.RTcoef(idom, 10);
			// Ci+1,j-1
			rt.RTcoef(idom, 11) = 0.5 * rt.c_difxy(iGlob) + 0.5 * rt.c_difyx(iGlob - gdom.nxhc); // Ci+1,j-1
			// Ci-1,j+1
			rt.RTcoef(idom, 12) = 0.5 * rt.c_difxy(iGlob - 1) + 0.5 * rt.c_difyx(iGlob); // Ci-1,j+1
			// Ci-1,j-1
			rt.RTcoef(idom, 13) = -0.5 * rt.c_difxy(iGlob - 1) - 0.5 * rt.c_difyx(iGlob - 1); // Ci-1,j-1
			// Ck+1//添加到上面
			rt.RTcoef(idom, 14) = -0.5 * rt.c_difxz(iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + 0.5 * rt.c_difxz(iGlob - 1) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - 0.5 * rt.c_difyz(iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + 0.5 * rt.c_difyz(iGlob - gdom.nxhc) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))); // Ck+1
			rt.RTcoef(idom, 5) += rt.RTcoef(idom, 14);
			//! Ci+1,k+1
			rt.RTcoef(idom, 15) = -0.5 * rt.c_difxz(iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - (1 - rt.wz(iGlob)) * rt.c_difzx(iGlob) / (2 * gdom.dx * gdom.dz(iGlob)); // Ci+1,k+1
			// Ck-1//添加到上面
			rt.RTcoef(idom, 16) = 0.5 * rt.c_difxz(iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - 0.5 * rt.c_difxz(iGlob - 1) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + 0.5 * rt.c_difyz(iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - 0.5 * rt.c_difyz(iGlob - gdom.nxhc) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))); // Ck-1
			rt.RTcoef(idom, 6) += rt.RTcoef(idom, 16);
			//! Ci+1,k-1
			rt.RTcoef(idom, 17) = 0.5 * rt.c_difxz(iGlob) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + rt.wz(iGlob - gdom.nxhc * gdom.nyhc) * rt.c_difzx(iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dx * gdom.dz(iGlob)); // Ci+1,k-1
			// Ci-1,k+1
			rt.RTcoef(idom, 18) = 0.5 * rt.c_difxz(iGlob - 1) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + (1 - rt.wz(iGlob)) * rt.c_difzx(iGlob) / (2 * gdom.dx * gdom.dz(iGlob)); // Ci-1,k+1
			//! Ci-1,k-1
			rt.RTcoef(idom, 19) = -0.5 * rt.c_difxz(iGlob - 1) / (gdom.dx * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - rt.wz(iGlob - gdom.nxhc * gdom.nyhc) * rt.c_difzx(iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dx * gdom.dz(iGlob)); // Ci-1,k-1
			// Ci+1///添加到上面
			rt.RTcoef(idom, 20) = -0.5 * rt.c_difyx(iGlob) + 0.5 * rt.c_difyx(iGlob - gdom.nxhc) - rt.wz(iGlob) * rt.c_difzx(iGlob) / (2 * gdom.dx * gdom.dz(iGlob)) + (1 - rt.wz(iGlob - gdom.nxhc * gdom.nyhc)) * rt.c_difzx(iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dx * gdom.dz(iGlob)); // Ci+1
			rt.RTcoef(idom, 1) += rt.RTcoef(idom, 20);
			// Ci-1//添加到上面
			rt.RTcoef(idom, 21) = 0.5 * rt.c_difyx(iGlob) - 0.5 * rt.c_difyx(iGlob - gdom.nxhc) + rt.wz(iGlob) * rt.c_difzx(iGlob) / (2 * gdom.dx * gdom.dz(iGlob)) - (1 - rt.wz(iGlob - gdom.nxhc * gdom.nyhc)) * rt.c_difzx(iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dx * gdom.dz(iGlob)); // Ci-1
			rt.RTcoef(idom, 2) += rt.RTcoef(idom, 21);
			// Cj+1,k+1
			rt.RTcoef(idom, 22) = -0.5 * rt.c_difyz(iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - (1 - rt.wz(iGlob)) * rt.c_difzy(iGlob) / (2 * gdom.dy * gdom.dz(iGlob)); // Cj+1,k+1
			// Cj+1,k-1
			rt.RTcoef(idom, 23) = 0.5 * rt.c_difyz(iGlob) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + rt.wz(iGlob - gdom.nxhc * gdom.nyhc) * rt.c_difzy(iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dy * gdom.dz(iGlob)); // Cj+1,k-1
			// Cj-1,k+1
			rt.RTcoef(idom, 24) = 0.5 * rt.c_difyz(iGlob - gdom.nxhc) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) + (1 - rt.wz(iGlob)) * rt.c_difzy(iGlob) / (2 * gdom.dy * gdom.dz(iGlob)); // Cj-1,k+1
			// Cj-1,k-1
			rt.RTcoef(idom, 25) = -0.5 * rt.c_difyz(iGlob - gdom.nxhc) / (gdom.dy * (0.5 * gdom.dz(iGlob - gdom.nxhc * gdom.nyhc) + gdom.dz(iGlob) + 0.5 * gdom.dz(iGlob + gdom.nxhc * gdom.nyhc))) - rt.wz(iGlob - gdom.nxhc * gdom.nyhc) * rt.c_difzy(iGlob - gdom.nxhc * gdom.nyhc) / (2 * gdom.dy * gdom.dz(iGlob)); // Cj-1,k-1
		});

		// 4. 应用边界条件对矩阵的修正
		for (int k = 0; k < rtgbc.size(); k++)
		{
			rtgbc[k].applyRTMatBC(rt, gw, gdom, par);
		}

		// 5. 将系数填充到稀疏矩阵 (CRS格式)
		Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
			int ii, jj, kk;			
			int irow = rtA.rt_ptr(idom); 				
			gdom.unpackIndices(idom, kk, jj, ii);

			// 填充逻辑需严格进行边界检查 (Boundary Checks)
			// 只有当邻居存在时，才填充非零元素

			if (kk > 0 && jj > 0){
				// rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny - 1;
				rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny - gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 25);
				// rtA.rt_val(irow) = 0;			
				irow++;				
			}

//ci-1,k-1
			if (ii > 0 && kk > 0){
				rtA.rt_ind(irow) = idom -1 - gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 19);
				// rtA.rt_val(irow) = 0;
							
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
				// rtA.rt_val(irow) = 0;
								
				irow++;				
			}	

//!ck-1,j+1
			if (kk > 0 && jj < gdom.ny-1){
				// rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny + 1;
				rtA.rt_ind(irow) = idom - gdom.nx*gdom.ny + gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 23);
				// rtA.rt_val(irow) = 0;
								
				irow++;				
			}

//!cj-1,i-1
			if (jj > 0 && ii > 0){
				rtA.rt_ind(irow) = idom - gdom.nx - 1;
				rtA.rt_val(irow) = rt.RTcoef(idom, 13);
				// rtA.rt_val(irow) = 0;
								
				irow++;				
			}
//!cj-1,i
			if (jj > 0)	{
				rtA.rt_ind(irow) = idom - gdom.nx;	
				rtA.rt_val(irow) = rt.RTcoef(idom, 4); 
				
				irow++;
				}
//!cj-1,i+1
			if (jj > 0 && ii < gdom.nx-1){
				rtA.rt_ind(irow) = idom - gdom.nx + 1;
				rtA.rt_val(irow) = rt.RTcoef(idom, 11); 
				// rtA.rt_val(irow) = 0;
								
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

//!cj+1,i-1
			if (jj < gdom.ny-1 && ii > 0){
				rtA.rt_ind(irow) = idom + gdom.nx - 1;
				rtA.rt_val(irow) = rt.RTcoef(idom, 12);
				// rtA.rt_val(irow) = 0;
							
				irow++;				
			}
//!cj+1,i
			if (jj < gdom.ny-1)	{
				rtA.rt_ind(irow) = idom + gdom.nx;	
				rtA.rt_val(irow) = rt.RTcoef(idom,3); 
				 
				irow++;
				}
//!cj+1,i+1
			if (jj < gdom.ny-1 && ii < gdom.nx-1){
				rtA.rt_ind(irow) = idom + gdom.nx + 1;
				rtA.rt_val(irow) = rt.RTcoef(idom, 9); 
				// rtA.rt_val(irow) = 0;
								
				irow++;				
			}

//!ck+1,j-1
			if (kk < gdom.nz-1 && jj > 0){
				// rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny - 1;
				rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny - gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 24);
				// rtA.rt_val(irow) = 0;
							
				irow++;				
			}

//ci-1,k+1
			if (ii > 0 && kk < gdom.nz-1){
				rtA.rt_ind(irow) = idom -1 + gdom.nx*gdom.ny;
				rtA.rt_val(irow) = rt.RTcoef(idom, 18);
				// rtA.rt_val(irow) = 0;
								 
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
				// rtA.rt_val(irow) = 0;
							
				irow++;				
			}	
//!ck+1,j+1
			if (kk < gdom.nz-1 && jj < gdom.ny-1){
				// rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny + 1;
				rtA.rt_ind(irow) = idom + gdom.nx*gdom.ny + gdom.nx;
				rtA.rt_val(irow) = rt.RTcoef(idom, 22);
				// rtA.rt_val(irow) = 0;
								
				irow++;				
			}

			rtA.rt_rhs(idom) = rt.RTcoef(idom,7); });
	}

	
};
#endif