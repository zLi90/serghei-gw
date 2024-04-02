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


class RTState   {

public:
	// transport time step size 
	real dt;
	// number of substances
	int n_mass;
	// molecular diffusivity 
	real d_base;
    // scalar concentration
    realArr2 c;
    // dispersion tensor
    realArr2 d;

    // Allocate state variables for groundwater
    inline void allocate (GwDomain &gdom) {
		// concentration is (n by 2): c_now, c_old
        c = realArr2("c", gdom.nCellMem, 2);
		// dispersion is (n by 3): d_xx, d_yy, d_zz
		// Note that we ignore d_xy, d_xz, etc. for now,
		// but they should be implemented later
        d = realArr2("d", gdom.nCellMem, 3);
    }

};

#endif

#endif
