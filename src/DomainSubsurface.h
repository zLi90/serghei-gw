#ifndef _DOMAIN_SUBSURFACE_H_
#define _DOMAIN_SUBSURFACE_H_

#include "define.h"

class DomainSubsurface {

 public:

  real cfl;
  real simLength;
  real etime;
  real dt;

  // raster variables

  int nx_glob;
  int ny_glob;
  int nz_glob;
  int nx;
  int ny;
  int nz;
  int ncells;
  int nCellDomain;

  boolArr isnodata;

  real xll;
  real yll;
  real zll;

  real dx;

  // bottom elevation of the subsurface domain
  real bottomZ;

  // physics variables

  int nsoilID; // number of soil types; used to initialize soilID

  // other variables

  int BCtype;
  int nIter;
  int countIterDt;

};

#endif
