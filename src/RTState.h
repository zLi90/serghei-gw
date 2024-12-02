/* -*- mode: c++; c-default-style: "linux" -*- */

/**
 * Reactive transport
 **/

#ifndef _RT_STATE_H_
#define _RT_STATE_H_

#if SERGHEI_SUBSURFACE_TRANSPORT

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include <set>
#include <math.h>

class RTState
{

public:
	// transport time step size
	real dt_init, dt_max, dt, dtOld;
	// 上游加权系数
	real Up_Weighting_vplus, Up_Weighting_vminus;
	// number of substances
	int n_mass;
	// molecular diffusivity
	real diffusion_molecular;
	// 横向弥散度纵向弥散度
	real alpha_T, alpha_L;
	//Peclet
	real Pe;

	// 反应相关参数
	real ReactionModule, AdsorptionDesorptionModel; // 反应模块开关，吸附解吸类型
	real rho_b;										// 多孔介质的密度
	
	real k1;                  //first-order reaction rate constant
	real k2;

	//Equilibrium or nonlinear sorption
	real Kf;										// Freundlich等温吸附常数
	real Nf;										// Nf为常数
	real alpha_D;									// Langmuir吸附常数
	real beta_D;									// 固体所能吸附的最大容量
	real Kd;										// 固液分配系数，L/g
	//Nonequilibrium sorption
	// real NonequilibriumSorptionModel;				// Nonequilibrium sorption model type
	real f;                                         //平衡吸附交换点位占总吸附点位的比例
	real beta;									// Nonequilibrium sorption model parameter,first-order mass transfer coefficient between the dissolved and sorbed phases, T^-1
	real lambda_1;									// first-order reaction rate for the dissolved phase, T-1
	real lambda_2;									// first-order reaction rate for the sorbed (solid) phase, T-1
	//real intermediate_variable;                      // 中间变量
	
	//Radioactive decay or biodegradation
	real lambda;									// 污染物的衰变系数量纲为[T-1]

	realArr Rf;										// 各单元阻滞因子
    
	// Liquid phase concentration
	realArr2 c;
	//Solid phase concentration
	realArr2 c_solid;

	// dispersion tensor
	// realArr2 d;
	realArr2 dcal;

	// 地下水浓度初始值模式
	int RT_Aquifer_initialMode, RT_Solid_initialMode, rt_scheme; // 地下水浓度初始值模式，固相初始值模式，反应传输计算方案

	real total_mass;
	realArr2 aveV;
	realArr2 aveVB;
	realArr2 RTcoef;
	// real phizzb;
	realArr tau;

	// python测试参数
	realArr c_advxx;
	realArr c_advyy;
	realArr c_advzz;
	realArr c_difxx;
	realArr c_difyy;
	realArr c_difzz;
	realArr c_difxy;
	realArr c_difxz;
	realArr c_difyz;

	// Allocate state variables for groundwater
	inline void
	allocate(GwDomain &gdom)
	{
		// concentration is (n by 2): c_now, c_old
		c = realArr2("c", gdom.nCellMem, 2);
		c_solid = realArr2("c_solid", gdom.nCellMem, 2);
		// dispersion is (n by 3): d_xx, d_yy, d_zz
		// Note that we ignore d_xy, d_xz, etc. for now,
		// but they should be implemented later
		// d = realArr2("d", gdom.nCellMem, 3);
		dcal = realArr2("dcal", gdom.nCellMem, 6);

		// advection = realArr("advection", gdom.nCellMem);

		aveV = realArr2("aveV", gdom.nCellMem, 4);
		// 边界流速
		aveVB = realArr2("aveVB", gdom.nCellMem, 4);
		// dispersion = realArr2("dispersion", gdom.nCellMem,4);
		// residual = realArr("residual", gdom.nCellMem);
		tau = realArr("tau", gdom.nCellMem);

		RTcoef = realArr2("RTcoef", gdom.nCell, 20);

		c_advxx = realArr("c_advxx", gdom.nCellMem);
		c_advyy = realArr("c_advyy", gdom.nCellMem);
		c_advzz = realArr("c_advzz", gdom.nCellMem);
		c_difxx = realArr("c_difxx", gdom.nCellMem);
		c_difyy = realArr("c_difyy", gdom.nCellMem);
		c_difzz = realArr("c_difzz", gdom.nCellMem);
		c_difxy = realArr("c_difxy", gdom.nCellMem);
		c_difxz = realArr("c_difxz", gdom.nCellMem);
		c_difyz = realArr("c_difyz", gdom.nCellMem);

		Rf = realArr("Rf", gdom.nCellMem);
	}
};

#endif

#endif
