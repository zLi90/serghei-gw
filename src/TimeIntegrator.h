/* -*- mode: c++ -*- */

#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include <stdlib.h>

#include "define.h"
#include "Parallel.h"
#include "Domain.h"
#include "State.h"
#include "BC.h"
#include "Edges.h"
#include "Indexing.h"
#include "FileIO.h"
#include "SourceSink.h"

#include <Kokkos_Core.hpp>

class TimeIntegrator
{

	Edges edge;
	Kokkos::Timer timer;

public:
	inline void stepForward(State &state, SourceSinkData &ss, std::vector<ExtBC> &extbc, Domain &dom, Exchange &exch, Parallel &par, FileIO &io)
	{

#if !SERGHEI_SWE_GW
		computeDt(state, dom, io);
#endif

		edge.computeDeltaStateSW(state, dom, exch, par);

		ss.ComputeSWSourceSink(state, dom);

		computeNewState(state, dom, ss);

		for (int k = 0; k < extbc.size(); k++)
		{ // should be done before the exchange (water depth might be modified).
			extbc[k].apply(state, dom);
		}

		exch.exchangeMPIh(state, dom, exch, par); // only neccesary to exchange the h (for wet-dry) but for the moment we exchange everything

		wetDryCorrection(state, dom);

		exch.exchangeMPIhuhv(state, dom, exch, par); // neccesary to exchange again because of the wet/dry correction

		for (int k = 0; k < extbc.size(); k++)
		{ // after getting the final values, the discharge is integrated at every BC. The reason for not doing this before is because the previous kernels could eventually modify the boundary cell values.
			extbc[k].integrate(state, dom);
		}
	}

	inline void computeGwExchange(State &state, const Domain &dom)
	{
		Kokkos::parallel_for(dom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii = dom.getIndex(idom);
            state.h(ii) += state.qss(idom) * dom.dt;
            // if(state.h(ii)<TOL12) {state.h(ii)=0.0;} 

			//!zzb 修改，为了避免农田下渗或者蒸发导致地表水深为0，无法施加流量边界
			// 将水深限制在一个非常小的正值以上
		if (state.h(ii)<TOL12){state.h(ii)=state.hmin;} });
	}

	inline void computeNewState(State &state, Domain &dom, SourceSinkData &ss)
	{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		timer.reset();

		// 网格参数（常量）
		const real cellArea = dom.dx() * dom.dx();
		const real dt = dom.dt;
		const real hmin = state.hmin;

		// 蒸发通量归约变量
		real evapFluxSum = 0.0;

		// 合并为一个完整的并行计算
		Kokkos::parallel_reduce("computeNewState", dom.nCell, KOKKOS_LAMBDA(int iGlob, real &localEvapSum) {
			int ii = dom.getIndex(iGlob);

			// 获取当前状态
			real z = state.z(ii);
			real hold = state.h(ii);
			real huold = state.hu(ii);
			real hvold = state.hv(ii);
			bool nodata = state.isnodata(ii);

			real hf, huf, hvf;
			int ncells = dom.nCellMem;

			// ==== 第一步：通量计算 ====
			hf = hold - dom.dt * (state.dsw0(ii) + state.dsw1(ii)) / dom.dx();

			int i, j;
			unpackIndicesUniformGrid(iGlob, dom.ny + 2 * hc, dom.nx + 2 * hc, j, i);

			// ==== 第二步：降雨计算 ====
			if (dom.isRain)
			{
				// if (j < dom.ny + hc - 1 && j > hc + 1) // 不包括东西边界
				// {

				hf += ss.rainRate(ii) * dom.dt;

				// 负降雨（蒸发）处理
				if (hf < 0.0)
				{
					hf = 0.0;
				}
				// }
			}

// ==== 第三步：蒸发计算 ====
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
			real evapRateActual = 0.0;

			// 1. 获取基础蒸发率 (确保这里获取的是潜在蒸发能力，而非受土壤限制后的蒸发)
			real evapRate = ss.evapRate(ii);

			// 2. 计算潜在蒸发深度 (确保 evapRate 单位是 m/s, dt 是 s, 结果是 m)
			// 既然 evapRate 是负值(代表流出)，这里乘以 -1 变为正的深度值便于比较
			real evapPotentialDepth = -1.0 * evapRate * dt;

			// 调试建议：如果在 Debug 模式，打印关键值
			// printf("Cell: %d, hf: %f, hmin: %f, Rate: %e, PotDepth: %e\n", ii, hf, hmin, evapRate, evapPotentialDepth);

			// 3. 逻辑分支处理
			// 必须先检查 hf 是否大于 hmin，否则没有水可供蒸发
			if (hf > hmin)
			{
				// 针对 evapRate 可能为正(凝结)或负(蒸发)的防御性编程
				// 如果是蒸发(Rate < 0 -> PotDepth > 0)
				if (evapPotentialDepth > 0.0)
				{
					if (evapPotentialDepth < (hf - hmin)) // 剩余水深足够蒸发，且蒸发后保留 hmin 底水
					{
						// 情况1：水足够
						hf -= evapPotentialDepth;
						evapRateActual = evapRate; // 直接使用原速率，减少浮点运算误差
					}
					else
					{
						// 情况2：水不够蒸发了 (会蒸干到 hmin)
						real availableDepth = hf - hmin;
						evapRateActual = -1.0 * availableDepth / dt; // 反算实际能发生的蒸发率(负值)
						hf = hmin;
					}
				}
				else
				{
					// 情况4：凝结 (evapRate >= 0)，直接增加水深
					// 假设 evapRate 为正表示水进入
					hf -= evapPotentialDepth; // PotDepth 是负的，减去负值等于加
					evapRateActual = evapRate;
				}
			}
			else
			{
				// 情况3：水深不足 (hf <= hmin)
				evapRateActual = 0.0;
				// 可以在这里强制 hf = hmin 保持数值稳定，或者不做处理
			}

			localEvapSum += evapRateActual * cellArea;
			state.surfaceEvaporation(ii) = evapRateActual;
#endif

			// ==== 第四步：下渗计算 ====
			if (ss.inf.model)
			{
				ss.inf.rate(ii) = min(ss.inf.rate(ii), hf / dom.dt);
				ss.inf.rate(ii) = max(ss.inf.rate(ii), 0.0);
				hf -= ss.inf.rate(ii) * dom.dt;
			}

			// ==== 第五步：干湿处理 ====
			if (hf < TOL_MACHINE_ACCURACY || nodata)
			{
				hf = 0.0;
			}

			// ==== 第六步：动量计算 ====
			if (hf < state.hmin)
			{
				huf = 0.0;
				hvf = 0.0;
			}
			else
			{
				real mx = huold - (state.dsw0(ii + ncells) + state.dsw1(ii + ncells)) * dom.dt / dom.dx();
				real my = hvold - (state.dsw0(ii + 2 * ncells) + state.dsw1(ii + 2 * ncells)) * dom.dt / dom.dx();

#if SERGHEI_POINTWISE_FRICTION
				real nsq = state.roughness(ii) * state.roughness(ii);
				real modM = sqrt(mx * mx / hold / hold + my * my / hold / hold);
				if (nsq > TOL12 && hold >= state.hmin && modM > TOL12)
				{
					real tt = dom.dt * GRAV * nsq * modM / (hold * cbrt(hold));
					real ff = sqrt(1.0 + 4.0 * tt);
					huf = -0.5 * (mx - mx * ff) / tt;
					hvf = -0.5 * (my - my * ff) / tt;
				}
				else
				{
					huf = mx;
					hvf = my;
				}
#else
				huf = mx;
				hvf = my;
#endif
			}

			// 动量归零处理
			if (fabs(huf) < TOL_ZERO_MOMENTUM)
				huf = 0.0;
			if (fabs(hvf) < TOL_ZERO_MOMENTUM)
				hvf = 0.0;

			// ==== 第七步：更新状态 ====
			state.h(ii) = hf;
			state.hu(ii) = huf;
			state.hv(ii) = hvf;

#if SERGHEI_SURFACE_TRANSPORT
			state.h4rtsw(ii, 0) = state.h(ii);
			state.hu4rtsw(ii, 0) = state.hu(ii);
			state.hv4rtsw(ii, 0) = state.hv(ii);
#endif

			// ==== 第八步：重置通量贡献 ====
			state.dsw0(ii) = 0.0;
			state.dsw0(ii + ncells) = 0.0;
			state.dsw0(ii + 2 * ncells) = 0.0;
			state.dsw1(ii) = 0.0;
			state.dsw1(ii + ncells) = 0.0;
			state.dsw1(ii + 2 * ncells) = 0.0;

#if SERGHEI_MAXFLOOD
			if (hf > state.hMax(ii))
			{
				state.hMax(ii) = hf;
				state.time_hMax(ii) = dom.etime;
			}
			real mom = mysqrt(huf * huf + hvf * hvf);
			if (mom > state.momentumMax(ii))
				state.momentumMax(ii) = mom;
#endif
		},
								Kokkos::Sum<real>(evapFluxSum));

// 保存蒸发通量
// if (dom.isEvap)
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
		// {
		ss.evapFluxActual = evapFluxSum;
		// }
#endif
		dom.timers.swe += timer.seconds();
	}

	inline void dtMatchOutput(Domain &dom, const FileIO &io)
	{
		// correction to match output times
		// std::cout << GGD << dom.etime << "\t" << dom.etime+dom.dt << "\t" << io.numOut << "\t" << io.outFreq*io.numOut << "\t" << io.numOut*io.outFreq + dom.startTime << std::endl;
		if (dom.etime + dom.dt > dom.startTime + io.numOut * io.outFreq)
			dom.dt = io.numOut * io.outFreq + dom.startTime - dom.etime;
		if (dom.etime + dom.dt > dom.endTime)
		{
			dom.dt = dom.endTime - dom.etime;
		}
#if SERGHEI_DEBUG_DT
		std::cout << "time = " << dom.etime << "\tdt_cor = " << dom.dt << std::endl;
#endif
	}

	inline void computeDt(State &state, Domain &dom, FileIO &io)
	{
		timer.reset();

		dom.dt = 1.e7;

		Kokkos::parallel_reduce("reduceDt", dom.nCell, KOKKOS_LAMBDA(int iGlob, real &dt) {
    int ii = dom.getIndex(iGlob);
	 	real h=state.h(ii);
	 	real hu=state.hu(ii);
	 	real hv=state.hv(ii);
		dt=min(dt,1.e6);
		if(h>TOL12){
			dt=min(dt,dom.dx()/(fabs(hu/h)+sqrt(GRAV*h)));
			dt=min(dt,dom.dx()/(fabs(hv/h)+sqrt(GRAV*h)));
		} }, Kokkos::Min<real>(dom.dt));

		Kokkos::fence();

		real dtloc = dom.dt;
		int ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, SERGHEI_MPI_REAL, MPI_MIN, MPI_COMM_WORLD);

		dom.dt *= dom.cfl;

#if SERGHEI_DEBUG_DT
		std::cout << "time = " << dom.etime << "\tdt_cfl = " << dom.dt << "\tdom.cfl = " << dom.cfl << std::endl;
#endif
		if (dom.dt > 1e5)
		{ // it means that evertyhing is dry.
			if (dom.isRain)
			{
				// if there is rain, we impose a time step equivalent for h=1
				// this is to make sure we capture the start of the rain
				// TODO: improve this using the known rainfall signal
				dom.dt = dom.dx() / (1 + sqrt(GRAV)); // dx() purposely used here to fail upon compilation when moving to adaptive mesh
#if SERGHEI_DEBUG_DT
				std::cout << "time = " << dom.etime << "\tdt_rain = " << dom.dt << std::endl;
#endif
			}
		}
		// dom.dt = 10;//!手动控制dt时间
		dtMatchOutput(dom, io);
		dom.timers.dt = timer.seconds();
	}

	inline void wetDryCorrection(State &state, Domain &dom)
	{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		timer.reset();

		Kokkos::parallel_for("wetDryCorrection", dom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
			int i, j;
			unpackIndicesUniformGrid(iGlob, dom.ny + 2 * hc, dom.nx + 2 * hc, j, i);
			if (i > hc - 1 && i < dom.nx + hc && j > hc - 1 && j < dom.ny + hc)
			{ // computational domain without halo cells
				real hij = state.h(iGlob);
				real zij = state.z(iGlob);
				int isB = state.isBound(iGlob);
				if (hij >= state.hmin)
				{
					if (((hij + zij < state.z(iGlob + 1)) && state.h(iGlob + 1) < TOL_WETDRY) || ((hij + zij < state.z(iGlob - 1)) && state.h(iGlob - 1) < TOL_WETDRY) || (isB == 0 && state.isnodata(iGlob + 1)) || (isB == 0 && state.isnodata(iGlob - 1)))
					{
						state.hu(iGlob) = 0.0;
					}
					if (((hij + zij < state.z(iGlob + dom.nx + 2 * hc)) && state.h(iGlob + dom.nx + 2 * hc) < TOL_WETDRY) || ((hij + zij < state.z(iGlob - (dom.nx + 2 * hc))) && state.h(iGlob - (dom.nx + 2 * hc)) < TOL_WETDRY) || (isB == 0 && state.isnodata(iGlob + dom.nx + 2 * hc)) || (isB == 0 && state.isnodata(iGlob - (dom.nx + 2 * hc))))
					{
						state.hv(iGlob) = 0.0;
					}
				}
			} });

		dom.timers.swe += timer.seconds();
	}

}; // end of TimeIntegrator class

#endif
