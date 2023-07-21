/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "define.h"
#include "Exchange.h"
#include "TimeIntegrator.h"
#include "Indexing.h"
#include "Parser.h"
#include "SWSourceSink.h"


class Initializer{

public:


  void initializeMPI( int *argc , char ***argv , Parallel &par ) {
    int ierr = MPI_Init( argc , argv );
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);

    //Determine if I'm the master process
    if (par.myrank == 0) {
      par.masterproc = 1;
    } else {
      par.masterproc = 0;
    }
  }

  int initialize(State &state, SourceSinkData &ss, ExternalBoundaries &ebc, Domain &dom, Parallel &par, TimeIntegrator &tint, surfaceIntegrator &sint, boundaryIntegrator &bint, Parser &parser, Exchange &exch, FileIO &io, std::string inFolder, std::string outFolder) {

    if(!parser.readDimensions(inFolder, dom, state, par, io)){
      return 0;
    }

    //buildDomainDecomposition(state, dom, par);
    dom.buildDomainDecomposition(par);

    dom.initialise();
    state.allocate(dom);

    if(!parser.readInputFiles(inFolder, dom, state, ss, ebc, par, io)){
      return 0;
    }

    dom.getStatistics();
    ss.allocate(dom);
    sint.initialize(state,dom,ss);
    bint.initialize(ebc.extbc);

    if(par.masterproc){
      if(!parser.createOutputDir(outFolder)){
	return 0;
      }
    }

    exch.exchangeIniMPI(state,dom,exch,par);

    //initialize boundaries
    boundaryIni(state,dom,par,ebc.extbc);

    // Output the initial model state
    io.outputIni(state, dom, ss, par,outFolder);

    return 1;

  }



	inline void boundaryIni(State &state, Domain &dom, Parallel &par, std::vector<ExtBC> &extbc){

		//impose boundary conditions in the outer (full) domain in the case of reflective boundary conditions. Periodic and transmissive are default in exchangeIniMPI

		if(dom.BCtype==BC_REFLECTIVE){
			//boundary conditions (halo) for outer domain (periodic/transmissive by default with the exchange in exchangeIniMPI)
			if(dom.iE){
				boundaryReflectiveE(state,dom);
			}
			if(dom.iW){
				boundaryReflectiveW(state,dom);
			}
			if(dom.iS){
				boundaryReflectiveS(state,dom);
			}
			if(dom.iN){
				boundaryReflectiveN(state,dom);
			}

		}

		//numerical boundaries. Remove the high walls
		for (int k = 0; k < extbc.size(); k ++) {
			removeElevationNumericalBoundaries(state, extbc[k], dom,par);
		}


	}

	inline void removeElevationNumericalBoundaries(State &state, ExtBC &extbc, Domain &dom,  Parallel &par){

		Kokkos::parallel_for("remove_elevation_numerical_boundaries", extbc.ncellsBC, KOKKOS_LAMBDA (int iGlob) {
			int ii = extbc.bcells[iGlob];
			int i, j;
			unpackIndices(ii,dom.ny+2*haloc,dom.nx+2*haloc,j,i);

			if(i==dom.nx+haloc-1&&dom.iE){
				state.z(ii+1)=state.z(ii);
			}
			if(i==haloc&&dom.iW) {
				state.z(ii-1)=state.z(ii);
			}
			if(j==dom.ny+haloc-1&&dom.iS){
				state.z(ii+dom.nx+2*haloc)=state.z(ii);
			}
			if(j==haloc&&dom.iN){
				state.z(ii-(dom.nx+2*haloc))=state.z(ii);
			}

		});

	}


	inline void boundaryReflectiveW(State &state, Domain &dom){

		Kokkos::parallel_for("boundaryReflectiveW", dom.ny*haloc , KOKKOS_LAMBDA (int iGlob) {
			int rx, ry;
			unpackIndices(iGlob,dom.ny,haloc,ry,rx);
			int ii=(haloc+ry)*(dom.nx+2*haloc)+haloc;
			//west boundary
			state.z(ii-rx-1)	=1e4; //10000m high
			state.h(ii-rx-1)	=0.0;
			state.hu(ii-rx-1)	=0.0;
			state.hv(ii-rx-1)	=0.0;
		 });
	}

	inline void boundaryReflectiveE(State &state,Domain &dom){

		Kokkos::parallel_for("boundaryReflectiveE", dom.ny*haloc , KOKKOS_LAMBDA (int iGlob) {
			int rx, ry;
			unpackIndices(iGlob,dom.ny,haloc,ry,rx);
			int ii=(haloc+ry)*(dom.nx+2*haloc)+haloc+dom.nx-1;
			//east boundary
			state.z(ii+rx+1) 	= 1e4; //10000m high
			state.h(ii+rx+1)	=0.0;
			state.hu(ii+rx+1)	=0.0;
			state.hv(ii+rx+1)	=0.0;

		 });
	}

	inline void boundaryReflectiveN(State &state,Domain &dom){

		Kokkos::parallel_for("boundaryReflectiveN", haloc*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      	int rx, ry;
      	unpackIndices(iGlob,haloc,dom.nx,ry,rx);
			int ii=haloc*(dom.nx+2*haloc)+haloc+rx;
			//north boundary
			state.z(ii-(ry+1)*(dom.nx+2*haloc))	=1e4; //10000m high
			state.h(ii-(ry+1)*(dom.nx+2*haloc))	=0.0;
			state.hu(ii-(ry+1)*(dom.nx+2*haloc))=0.0;
			state.hv(ii-(ry+1)*(dom.nx+2*haloc))=0.0;

		 });
	}

	inline void boundaryReflectiveS(State &state,Domain &dom){

		Kokkos::parallel_for("boundaryReflectiveS", haloc*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      	int rx, ry;
      	unpackIndices(iGlob,haloc,dom.nx,ry,rx);
			int ii=(haloc+dom.ny-1)*(dom.nx+2*haloc)+haloc+rx;
			//south boundary
			state.z(ii+(ry+1)*(dom.nx+2*haloc))	=1e4; //10000m high
			state.h(ii+(ry+1)*(dom.nx+2*haloc))	=0.0;
			state.hu(ii+(ry+1)*(dom.nx+2*haloc))=0.0;
			state.hv(ii+(ry+1)*(dom.nx+2*haloc))=0.0;

		 });
	}

};

#endif
