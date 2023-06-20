
#ifndef _CONST_H_
#define _CONST_H_

#include "math.h"
#include "const.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

typedef double         real;
typedef unsigned long ulong;
typedef unsigned int  uint;

#ifdef __NVCC__
  typedef Kokkos::View<real*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> realArr;
  typedef Kokkos::View<real**     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> realArr2;
  typedef Kokkos::View<int*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> intArr;
  typedef Kokkos::View<int**     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> intArr2;
  typedef Kokkos::View<bool*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> boolArr;

#else
  typedef Kokkos::View<real*     ,Kokkos::LayoutRight> realArr;
  typedef Kokkos::View<real**     ,Kokkos::LayoutRight> realArr2;
  typedef Kokkos::View<int*     ,Kokkos::LayoutRight> intArr;
  typedef Kokkos::View<int**     ,Kokkos::LayoutRight> intArr2;
  typedef Kokkos::View<bool*     ,Kokkos::LayoutRight> boolArr;
#endif
// dual views
typedef Kokkos::DualView<int**> dualInt;
typedef Kokkos::DualView<double**> dualDbl;
typedef dualDbl::execution_space dspace;

#ifdef __NVCC__
#define _HOSTDEV __host__ __device__
#else
#define _HOSTDEV
#endif

inline _HOSTDEV real operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

inline _HOSTDEV double mypow ( double const x , double const p ) { return pow (x,p); }
inline _HOSTDEV float  mypow ( float  const x , float  const p ) { return powf(x,p); }
inline _HOSTDEV double mysqrt( double const x ) { return sqrt (x); }
inline _HOSTDEV float  mysqrt( float  const x ) { return sqrtf(x); }
inline _HOSTDEV double myfabs( double const x ) { return fabs (x); }
inline _HOSTDEV float  myfabs( float  const x ) { return fabsf(x); }

template <class T> inline _HOSTDEV T minn( T const v1 , T const v2 ) {
  if (v1 < v2) { return v1; }
  else         { return v2; }
}
template <class T> inline _HOSTDEV T maxx( T const v1 , T const v2 ) {
  if (v1 > v2) { return v1; }
  else         { return v2; }
}

template <class T> inline _HOSTDEV int sgn(T const val) {
    return (T(0) < val) - (val < T(0));
}

class SergheiTimers{
public:
  Kokkos::Timer serghei;
  Kokkos::Timer out;
  Kokkos::Timer swe;
  real total=0;
  real Tinit=0;
  real Tout=0;
  real Tswe=0;
  real Traininf=0;
  real Tsweflux=0;
  real Texchange=0;
  real Tintegrate=0;
};

#endif
