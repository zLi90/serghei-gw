#ifndef _GW_DOMAIN_H_
#define _GW_DOMAIN_H_

#include "define.h"
#include "Domain.h"

class GwDomain : public Domain {

public:
    real dt_init, dt_max, dt, dtOld;
    real thickH, topZ, xll, yll, zll, dx, dy, aev, hmin;
    int nz_glob, nhalo, nx, ny, nz, nSoilID, gw_scheme, isRain;
    int nxhc, nyhc, nzhc, n_substep;
    realArr z, dz, sinx, cosx, siny, cosy, qrain;
    intArr2 hpair;
    intArr isnodata;
};



class GwBC {

public:
	int ncellsBC = 0; //number of bcells
	intArr bcells; //array of indexes of boundary cells
	real normalx, normaly; //direction set by user for inflow/outflow
	int location; //1->west, 2->north, 3->east, 4-> south
	int bctypeXP, bctypeXM, bctypeYP, bctypeYM, bctypeZP, bctypeZM;
    intArr topBC; // cell-by-cell BC type for the top layer 1->head 2->flux
	real qbcXP, qbcXM, qbcYP, qbcYM, qbcZP, qbcZM;
	real hbcXP, hbcXM, hbcYP, hbcYM, hbcZP, hbcZM;
	int isInDomain;
  	realArr bcvals;

};

#endif
