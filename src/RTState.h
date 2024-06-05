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
	//上游加权系数
	real Up_Weighting_vplus, Up_Weighting_vminus;
	// number of substances
	int n_mass;
	// molecular diffusivity
	real d_base;
	//横向弥散度
	real alpha_T;
	//纵向弥散度
	real alpha_L;
	// scalar concentration
	realArr2 c;
	// dispersion tensor
	// realArr2 d;
	realArr2 dcal;
	//地下水浓度初始值模式
	int RTinitialMode, rt_scheme;

	real total_mass;
	realArr2 aveV;
	realArr2 aveVB;
	realArr2 RTcoef;
	// real phizzb;
	realArr tau; 

	//python测试参数
	realArr c_advxx;
	realArr c_advyy ;
	realArr c_advzz ;
	realArr c_difxx ;
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
		// dispersion is (n by 3): d_xx, d_yy, d_zz
		// Note that we ignore d_xy, d_xz, etc. for now,
		// but they should be implemented later
		// d = realArr2("d", gdom.nCellMem, 3);
		dcal = realArr2("dcal", gdom.nCellMem, 6);

		


		// advection = realArr("advection", gdom.nCellMem);

		aveV = realArr2("aveV", gdom.nCellMem,4);
		//边界流速
		aveVB = realArr2("aveVB", gdom.nCellMem,4);
		// dispersion = realArr2("dispersion", gdom.nCellMem,4);
		// residual = realArr("residual", gdom.nCellMem);
		tau =realArr("tau", gdom.nCellMem);

		RTcoef = realArr2("RTcoef", gdom.nCell, 20);

		//python测试参数
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
