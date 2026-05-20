
#ifndef _CONST_H_
#define _CONST_H_

#include "math.h"
#include "const.h"
#include "mpi.h"
#include <Kokkos_Core.hpp>

// SWE boundary conditions
#define BC_PERIODIC 1
#define BC_REFLECTIVE 2
#define BC_TRANSMISSIVE 3

// Subsurface initial mode
#define IC_SAT 1
#define IC_H 2
#define IC_WC 3
#define IC_WT 4

#define SERGHEI_FLOAT 1
#define SERGHEI_DOUBLE 2

#ifndef SERGHEI_REAL
  #define SERGHEI_REAL SERGHEI_DOUBLE
#endif

#ifndef SERGHEI_MPI_REAL
  #if SERGHEI_REAL == SERGHEI_DOUBLE
    #define SERGHEI_MPI_REAL MPI_DOUBLE
  #elif SERGHEI_REAL == SERGHEI_FLOAT
    #define SERGHEI_MPI_REAL MPI_FLOAT
  #endif
#endif

#if SERGHEI_REAL == SERGHEI_DOUBLE
  #define SERGHEI_VTK_REAL "DOUBLE"
#elif SERGHEI_REAL == SERGHEI_FLOAT
  #define SERGHEI_VTK_REAL "FLOAT"
#endif

#if SERGHEI_REAL == SERGHEI_DOUBLE
  typedef double         real;
  #define TOL_MASS_ERROR TOL8
  #define TOLDRY TOL12
  #define TOL_ZERO_MOMENTUM TOL12
  #define TOL_WETDRY TOL12
  #define TOL_MACHINE_ACCURACY TOL12
#endif
#if SERGHEI_REAL == SERGHEI_FLOAT
  typedef float         real;
  #define TOL_MASS_ERROR TOL5
  #define TOLDRY TOL6
  #define TOL_ZERO_MOMENTUM TOL12
  #define TOL_WETDRY TOL12
  #define TOL_MACHINE_ACCURACY TOL6
#endif
typedef unsigned long ulong;
typedef unsigned int  uint;
#if SERGHEI_DRAINAGE_MODEL
int Steps = 0;
const double OMEGA  = 0.5;
const  double MAXVELOCITY =  50.; 

//drainage parameters
const double Cw = 0.6;
const double w = 1.6;
const double Co = 0.62;
const double Amh = 0.09;
#define  N_TRANSECT_TBL  51 
char     MemPoolAllocated;    
int  Ntokens;                   
int  Mobjects[MAX_OBJ_TYPES];   
int  Mnodes[MAX_NODE_TYPES];    
int  Mlinks[MAX_LINK_TYPES];    
int  Mevents;                  
// const int MAXERRS = 100;       
const double MIN_DELTA_Z = 0.001;

#endif


#ifdef __NVCC__
  typedef Kokkos::View<real*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> realArr;
  typedef Kokkos::View<real**     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> realArr2;
  typedef Kokkos::View<int*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> intArr;
  typedef Kokkos::View<int**     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> intArr2;
  typedef Kokkos::View<bool*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> boolArr;
  typedef Kokkos::View<double*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> doubleArr;
#else
  typedef Kokkos::View<real*     ,Kokkos::LayoutRight> realArr;
  typedef Kokkos::View<real**     ,Kokkos::LayoutRight> realArr2;
  typedef Kokkos::View<int*     ,Kokkos::LayoutRight> intArr;
  typedef Kokkos::View<int**     ,Kokkos::LayoutRight> intArr2;
  typedef Kokkos::View<bool*     ,Kokkos::LayoutRight> boolArr;
  typedef Kokkos::View<double*     ,Kokkos::LayoutRight> doubleArr;
#endif

#ifdef __NVCC__
#define _HOSTDEV __host__ __device__
#else
#define _HOSTDEV
#endif

KOKKOS_INLINE_FUNCTION real operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

KOKKOS_INLINE_FUNCTION double mypow ( double const x , double const p ) { return pow (x,p); }
KOKKOS_INLINE_FUNCTION float  mypow ( float  const x , float  const p ) { return powf(x,p); }
KOKKOS_INLINE_FUNCTION double mysqrt( double const x ) { return sqrt (x); }
KOKKOS_INLINE_FUNCTION float  mysqrt( float  const x ) { return sqrtf(x); }
KOKKOS_INLINE_FUNCTION double myfabs( double const x ) { return fabs (x); }
KOKKOS_INLINE_FUNCTION float  myfabs( float  const x ) { return fabsf(x); }
KOKKOS_INLINE_FUNCTION int  myfabs( int  const x ) { return abs(x); }

/*
template <class T> KOKKOS_INLINE_FUNCTION T min( T const v1 , T const v2 ) {
  if (v1 < v2) { return v1; }
  else         { return v2; }
}

template <class T> KOKKOS_INLINE_FUNCTION T max( T const v1 , T const v2 ) {
  if (v1 > v2) { return v1; }
  else         { return v2; }
}
*/
template <class T1, class T2> KOKKOS_INLINE_FUNCTION T1 min (T1 const v1, T2 const v2){
  if (v1 < v2) { return (T1) v1; }
  else         { return (T1) v2; }
}

template <class T1, class T2> KOKKOS_INLINE_FUNCTION T1 max (T1 const v1, T2 const v2){
  if (v1 > v2) { return (T1) v1; }
  else         { return (T1) v2; }
}

template <class T> KOKKOS_INLINE_FUNCTION int sgn(T const val) {
    return (T(0) < val) - (val < T(0));
}
template <class T> KOKKOS_INLINE_FUNCTION T sign(T const x, T const y) {
    return (y >= T(0)) ? std::fabs(x) : -std::fabs(x);
}

/*
// a bit hacky, perhaps better to include in template
KOKKOS_INLINE_FUNCTION float min(float v1, double v2){
  return(min(v1,(float) v2));
}
KOKKOS_INLINE_FUNCTION float max(float v1, double v2){
  return(max(v1,(float) v2));
}
*/


class SergheiTimers{
public:
  real total=0;
  real init=0;
  real out=0;
  real swe=0;
  real raininf=0;
  real sweflux=0;
  real exchange=0;
  real integrate=0;
  real integrateMPI=0;
  real dt=0;
  real sweBC=0;
  real halo=0;
  #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
    real update_particles=0;
    real initialization_particles=0;
    real output_particles=0;
  #endif
  // timers for the subsurface solver
  real gw = 0;
  real gwBC = 0;
  real gwlinsys = 0;
  real gwlinsol = 0;
  real gwUpdateK = 0;
  real gwUpdateQ = 0;
  real gwUpdateWC = 0;
  real gwMPI = 0;
  real gwIntegrate = 0;
      // timers for the drainage model
  real drainage = 0;
  real drainageExchange = 0;
  real drainageLinkFlows = 0;
  real drainageNodeDepths = 0;
};

#endif
