#ifndef _PARTICLETRACKING_H_
#define _PARTICLETRACKING_H_

#ifndef SERGHEI_PARTICLE_TRACKING
#define SERGHEI_PARTICLE_TRACKING 0
#endif

#include "define.h"
#include "tools.h"
#include "Indexing.h"
#include "Domain.h"

#if SERGHEI_PARTICLE_TRACKING

#define ZERO_VELOCITY 1E-6
#define ZERO_DEPTH 1E-12
#define RANDOM_COORDINATES 100
#define RANDOM_CELLS 101


typedef struct{
  geometry::point x;
  real t0;
  real lifetime;
} particleType;

typedef struct{
  geometry::point v;
  real dt;
  int edge;
  real dEdge;
} particleSupport;

#ifdef __NVCC__
  typedef Kokkos::View<particleType*,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> particleArr;
#else
  typedef Kokkos::View<particleType*,Kokkos::LayoutRight> particleArr;
#endif

class ParticleTracker{
  public:
    int N;  // number of particles
    particleArr particles;  // particles
    int initialDist;
    int seed = (int) time(NULL);

  inline void initalise(){
    // N needs to have been previously set
    particles = particleArr("particles",N);
    Kokkos::parallel_for("initialiseParticles",N, KOKKOS_LAMBDA(int i){
      particles(i).t0 = particles(i).lifetime = 0.;
    });
  };

  inline void initialDistribution(const Domain &dom){ //, Geometry *geom){   // TODO, define a geometry class
    switch(initialDist){
	    srand(seed);
      case RANDOM_CELLS:
        int i,j;
        // TODO collect i and j for the 1D cells array
        for(int ip=0; ip<N; ip++){
          int r = rand() % dom.nCell;
          particles(ip).x(_X) = dom.xll + i * dom.dx;
          particles(ip).x(_Y) = dom.yll + j * dom.dx;
        }
      break;
      /*
      case RANDOM_COORDINATES:
        for(int ip=0; ip<N; ip++){
          for(int kk=0; kk<N_SPATIAL_DIM; kk++) particles(ip).x(kk) = rand() % (geom.max[kk]-geom.min[kk]) + geom.min[kk];
        }
        break;
        */
    }
  }

private:
  inline particleSupport getParticleSupport(const particleType &p, const Domain &dom, const State &state){
    particleSupport ps;

    int iGlob = dom.getCellForPoint(p.x);
    int ii = dom.getIndex(iGlob);
    // get velocity field
    ps.v(_X) = ps.v(_Y) = 0;
    if(state.h(ii) > ZERO_DEPTH){
      ps.v(_X) = state.hu(ii)/state.h(ii);
      ps.v(_Y) = state.hv(ii)/state.h(ii);
    }

    // velocity unit vector
    real vmag = sqrt(ps.v(_X)*ps.v(_X) + ps.v(_Y)+ps.v(_Y));

    // if velocity is zero, early exit
    if(vmag < ZERO_VELOCITY) return(ps);

    // unit vector of velocity
    geometry::point unitv = ps.v / vmag;

    // get southwest corner
    geometry::point SW;
    int i,j;
    unpackIndices(iGlob,dom.ny,dom.nx,j,i);
    SW(_X) = i*dom.dx;
    SW(_Y) = j*dom.dx; // VERIFY THIS - is it ordered like a raster or not?

    real dEast,dWest,dNorth,dSouth;
    int easting=0,northing=0;
    int dEdge[N_CELL_EDGES];
    for(int kk=0; kk<N_CELL_EDGES; kk++) dEdge[kk] = dom.dx*1000;
    if(unitv(_X) != ZERO_VELOCITY){
      // evaluate distance to intersection with east and west edges
      dEdge[_EAST] = (SW(_X) + dom.dx - p.x(_X))/unitv(_X);
      dEdge[_WEST] = (SW(_X) - p.x(_X))/unitv(_X);
    }
    if(unitv(_Y) != ZERO_VELOCITY){
      // evaluate distance to intersection with north and south edges
      dEdge[_NORTH] = (SW(_Y) + dom.dx - p.x(_Y))/unitv(_Y);
      dEdge[_SOUTH] = (SW(_Y) - p.x(_Y))/unitv(_Y);
    }

    // the shortest positive intersection distance is the actual intersection with cell edges
    // negative distances are antiparallel to the velocity vector, and should be discarded
    ps.dEdge = dom.dx*1000;
    for(int kk=0; kk<N_CELL_EDGES; kk++){
      if(fabs(dEdge[kk]) < ps.dEdge || dEdge[kk]>0){
        ps.edge = kk;
        ps.dEdge = dEdge[kk];
      }
    }


    // time step to enforce particle to stay in cell
    ps.dt = ps.dEdge / vmag;
    // check that the future time is not beyond the future eulerian time
    ps.dt = min(ps.dt, dom.etime-ps.dt);
    return(ps);
  };

  public:
    inline void update(const Domain &dom, const State &state){
      // TODO domain decomposition?
      Kokkos::parallel_for("updateParticles",N, KOKKOS_LAMBDA(int i){
        while(particles(i).lifetime - particles(i).t0 < dom.etime){
          particleSupport ps = getParticleSupport(particles(i),dom,state);
          for(int kk=0; kk<N_SPATIAL_DIM; kk++) particles(i).x(kk)  = particles(i).x(kk) + ps.dt * ps.v(kk);
          particles(i).lifetime += ps.dt;
        }
      });
    };




};

#endif
#endif
