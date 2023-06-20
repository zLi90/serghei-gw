#ifndef _GEOMETRY_
#define _GEOMETRY_

#include "SArray.h"

#define _X 0
#define _Y 1
#define _Z 2

#define _SOUTH 0
#define _EAST 1
#define _NORTH 2
#define _WEST 3
#define N_SPATIAL_DIM 2
#define N_CELL_EDGES 4


namespace geometry{
  typedef SArray<real,N_SPATIAL_DIM> point;


  inline real distance(const point &p, const point &q){
    point d;
    real ss=0.;
    for(int i=0; i<N_SPATIAL_DIM; i++){
      d(i) = p(i) - q(i);
      ss += d(i)*d(i);
    }
    return(mysqrt(ss));
  };

  inline point crossProduct(const point &a, const point &b){
    point c;
    c(0) = a(1)*b(2) - a(2)*b(1);
    c(1) = a(2)*b(0) - a(0)*b(2);
    c(2) = a(0)*b(1) - a(1)*b(0);
    return(c);
  };

  inline geometry::point normalToLineInPlane(const point &a, const point &b){
    point c;
    c(0) = a(1) - b(1);
    c(1) = b(0) - a(0);
    return(c);
  };
};
#endif
