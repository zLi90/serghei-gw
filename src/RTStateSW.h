/* -*- mode: c++; c-default-style: "linux" -*- */

/**
 * Reactive transport for surface flow
 **/

#ifndef _RT_STATE_SW_H_
#define _RT_STATE_SW_H_

#if SERGHEI_SURFACE_TRANSPORT

#include "const.h"
#include "define.h"
#include "Domain.h"
#include "SArray.h"
#include "Indexing.h"
#include "RTState.h"
#include <set>
#include <math.h>


class RTStateSW 
{

public:
	
	real dt; // transport time step size
	realArr2 c, solute_mass;//, solute mass
	realArr csw4gw;//!用于地表地下水溶质交互时，地表水溶质浓度为边界条件，将sw溶质浓度复制至gw中
	realArr ConQss;//地表地下溶质交换通量q*c
	realArr Pe; // Peclet number
	realArr h; // water depth
	int RtInitialModeSW, InitialCon, n_mass;
	int east, west, north, south; // indices for neighboring cells
	real c0; // initial concentration for reactive transport

	// Allocate state variables for groundwater
	inline void
	allocate(Domain &dom)
	{

		Pe = realArr("Pe", dom.nCellMem);
		solute_mass = realArr2("solute_mass", dom.nCellMem, 2);

		c = realArr2("c", dom.nCellMem, 3);

		h = realArr("h", dom.nCellMem);

		csw4gw = realArr("csw4gw", dom.nCellMem);
		ConQss = realArr("ConQss", dom.nCell);


	}
};

#endif

#endif
