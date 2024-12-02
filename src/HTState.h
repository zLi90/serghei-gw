/* -*- mode: c++; c-default-style: "linux" -*- */

/**
 * Reactive transport
 **/

#ifndef _HT_STATE_H_
#define _HT_STATE_H_

#if SERGHEI_SUBSURFACE_HEAT

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include <set>
#include <math.h>

class HTState
{

public:
	//heat transport time step size
	real ht_init, ht_max, dt, dtOld;
	// 上游加权系数
	real Up_Weighting_vplus, Up_Weighting_vminus;
	// number of substances
	int n_mass;
	
	//the volumetric soil heat capacity
	realArr2 Cp;
	//the volumetric heat capacity of the liquid phase
	real Cw;
	//thermal conductivity :λ
	realArr2 lambda;

	real b1, b2, b3;
	real theta_n;
	// temperature
	realArr2 T;
	realArr2 HTcoef;

	// 地下水温度初始值模式
	int HT_Aquifer_initialMode, ht_scheme; // 地下水温度初始值模式，反应传输计算方案

	real total_mass;
	realArr residual;
	
	
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
		// concentration is (n by 2): T_now, T_old
		T = realArr2("T", gdom.nCellMem, 2);
		
		
		Cp = realArr2("Cp", gdom.nCellMem, 3);
		lambda = realArr2("lambda", gdom.nCellMem, 3);
        HTcoef = realArr2("HTcoef", gdom.nCell, 20);

		

		// dispersion = realArr2("dispersion", gdom.nCellMem,4);
		// residual = realArr("residual", gdom.nCellMem);
		tau = realArr("tau", gdom.nCellMem);

		

		c_advxx = realArr("c_advxx", gdom.nCellMem);
		c_advyy = realArr("c_advyy", gdom.nCellMem);
		c_advzz = realArr("c_advzz", gdom.nCellMem);
		c_difxx = realArr("c_difxx", gdom.nCellMem);
		c_difyy = realArr("c_difyy", gdom.nCellMem);
		c_difzz = realArr("c_difzz", gdom.nCellMem);
		c_difxy = realArr("c_difxy", gdom.nCellMem);
		c_difxz = realArr("c_difxz", gdom.nCellMem);
		c_difyz = realArr("c_difyz", gdom.nCellMem);

		
	}
};

#endif

#endif
