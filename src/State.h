
#ifndef _STATE_H_
#define _STATE_H_

#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include <set>

#ifndef SERGHEI_MAXFLOOD
#define SERGHEI_MAXFLOOD 0
#endif

typedef struct{
  real h=0;
  real hu=0;
  real hv=0;
  real z=0;
} swState;

// https://stackoverflow.com/questions/14712837/is-mpi-allreduce-on-a-structure-with-fields-of-the-same-type-portable
// define an MPI structure
/*
void defineMPIswState(MPI_Datatype *tstype) {
    const int count = 4;
    int          blocklens[count];
    MPI_Datatype types[count];
    MPI_Aint     disps[count];

    for (int i=0; i < count; i++) {
        types[i] = SERGHEI_MPI_REAL;
        blocklens[i] = 1;
    }

    disps[0] = offsetof(swState,h);
    disps[1] = offsetof(swState,hu);
    disps[2] = offsetof(swState,hv);
    disps[3] = offsetof(swState,z);

    MPI_Type_create_struct(count, blocklens, disps, types, tstype);
    MPI_Type_commit(tstype);
}

// define a reduction operation
void MPISUM_swState(void *in, void *inout, int *len, MPI_Datatype *type){
    // ignore type, just trust that it's our struct type

    swState *invals    = (swState*) in;
    swState *inoutvals = (swState*) inout;

    for (int i=0; i<*len; i++) {
      inoutvals[i].h  += invals[i].h;
      inoutvals[i].hu += invals[i].hu;
      inoutvals[i].hv += invals[i].hv;
      inoutvals[i].z  += invals[i].hv;
    }
    return;
}
*/



class State {

public:

  // SW variables
  realArr h;
  realArr hu;
  realArr hv;

  //elevation
  realArr z;

  //roughness
  realArr roughness;
  real hmin;

  //deltaFluxes
  realArr dsw0; //3 variables (h,hu,hv). left and south contribs
  realArr dsw1; //3 variables (h,hu,hv). right and north contribs

  boolArr isnodata; //contains 0 if is a regular cell, 1 if is nodata cell
  intArr isBound; //positive values for inlet boundaries, negative values for outlet bvoundaries, 0 for inner cells

  #if SERGHEI_MAXFLOOD
    realArr hMax;
    realArr momentumMax;
    realArr time_hMax;
  #endif


  inline void allocate(Domain &dom){
    h 				= realArr( "h" , dom.nCellMem );
    hu 			= realArr( "hu" , dom.nCellMem );
    hv 			= realArr( "hv" , dom.nCellMem );
    z 				= realArr( "z" , dom.nCellMem );
    roughness 	= realArr( "roughness" , dom.nCellMem );
    isnodata 	= boolArr( "isnodata" , dom.nCellMem );
	 isBound 	= intArr( "isBound" , dom.nCellMem );
    dsw0 			= realArr( "dsw0" , 3*dom.nCellMem );
    dsw1 			= realArr( "dsw1" , 3*dom.nCellMem );

    #if SERGHEI_MAXFLOOD
      hMax = realArr("hMax",dom.nCellMem);
      momentumMax = realArr("momMax",dom.nCellMem);
      time_hMax = realArr("timehMax",dom.nCellMem);

      Kokkos::parallel_for("initialise_maxflood",dom.nCell,KOKKOS_CLASS_LAMBDA(int iGlob) {
        int ii = dom.getIndex(iGlob);
        hMax(ii) = 0;
        momentumMax(ii)=0;
        time_hMax(ii)=0;
      });
    #endif
  }

  inline void filterDomain(const Domain &dom){
    Kokkos::parallel_for("filter_domain",dom.nCell,KOKKOS_CLASS_LAMBDA(int iGlob) {
      int ii = dom.getIndex(iGlob);
      if(isnodata(ii)){
        h(ii) = SERGHEI_NAN;
        hu(ii) = SERGHEI_NAN;
        hv(ii) = SERGHEI_NAN;
      }
    });
  }
};


class ShallowWater{
public:
  std::string initialMode;
  std::string frictionModel;
  std::string roughnessInput;
  std::set<std::string> initialModes = {"dry","h","h+z","file","netcdf"};
  std::set<std::string> frictionModels = {"none","manning","darcyweisbach","chezy"};
  real roughness = 0 ;
  real initialValue = 0;
  real hmin = -1;
};

#endif
