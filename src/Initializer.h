/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "define.h"
#include "Exchange.h"
#include "TimeIntegrator.h"
#include "Indexing.h"
#include "Parser.h"
#include "SWSourceSink.h"
#include "Subsurface.h"


class Initializer{

public:

	bool read = 1;

  void initializeMPI( int *argc , char ***argv , Parallel &par ) {
  	#if SERGHEI_DEBUG_WORKFLOW
	    std::cerr << GGD "Initialising MPI " << std::endl;
	  #endif
    int ierr = MPI_Init( argc , argv );
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);

    //Determine if I'm the master process
    if (par.myrank == SERGHEI_MASTERPROC) {
      par.masterproc = 1;
    } else {
      par.masterproc = 0;
    }
  }

  int initialize(State &state, SourceSinkData &ss, ExternalBoundaries &ebc, Domain &dom, Parallel &par, TimeIntegrator &tint, surfaceIntegrator &sint, boundaryIntegrator &bint, Parser &parser, Exchange &exch, FileIO &io, std::string inFolder, std::string outFolder){

    if(read) if(!parser.readDimensions(inFolder, dom, state, par, io)) return 0;
    
    dom.buildDomainDecomposition(par);
    
    dom.initialise();
    state.allocate(dom);
    
		if(read) if(!parser.readInputFiles(inFolder, dom, state, ss, ebc, par, io)) return 0;

    dom.getStatistics();
    ss.allocate(dom);
    sint.initialize(state,dom,ss);
    bint.initialize(ebc.extbc);
    //state.filterDomain(dom);

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
#if SERGHEI_SUBSURFACE_MODEL
    // initialize subsurface domain, ZhiLi20210219
  int initializeSubsurface(SubsurfaceState &statesub, State &state, DomainSubsurface &domsub, Domain &dom, Parallel &par, TimeIntegrator &tint, Parser &parser, FileIO &io, std::string inFolder, std::string outFolder)
    {
      int flag = -1;

      // read the subsurface domain information
      if(!parser.readSubsurfaceDimensions(inFolder, domsub, par))
	{
	  std::cerr << GOK << " Reading in subsurface dimensions failed." << std::endl;
	  flag = 0;
	  return flag;
	}

      // As of 2020-03-01 : Assume nz = nz_glob (no vertical decomposition)
      domsub.nz = domsub.nz_glob;
      domsub.nx = dom.nx;
      domsub.ny = dom.ny;
      // Assume serial
      domsub.nx_glob = domsub.nx;
      domsub.ny_glob = domsub.ny;
      domsub.dx = dom.dx;

      // allocate subsurface domain
      allocateDomainSubsurface(statesub, state, domsub, tint);

      /* now the state is allocated and array dimensions are known */
      /* we can proceed to read the van Genuchten soil parameters */
      if(!parser.readVGParameters(inFolder, statesub, domsub, par))
	{
	  std::cerr << RERROR << " Reading in van Genuchten parameters failed." << std::endl;
	  flag = 0;
	  return flag;
	}
    // read input / initial conditions for the subsurface domain
    if(!parser.readSubsurfaceState(inFolder, domsub, statesub, state, par, io)) return 0;

      // ignore source/sink and domain integrator for now


      // ignore BC initialization for now


      flag = 1;
      return flag;
    }
  #endif
    #if SERGHEI_SUBSURFACE_MODEL
    void allocateDomainSubsurface(SubsurfaceState &statesub, State &state, DomainSubsurface &dom, TimeIntegrator &tint)
    {
        // Initialize the time
        dom.nIter = 0;
        dom.countIterDt=0;
        dom.nCellMem=(dom.ny+2*hc)*(dom.nx+2*hc)*(dom.nz+2*hc); //cells plus halo cells, hc is number of overlapping halo cells defined in const.h
        dom.nCell = dom.nx*dom.ny*dom.nz; // WARNING UCM

	/* initialize edge values for Riemann solver */
        // tint.initialize(dom);

        //allocate the variables in the subsurface domain
        statesub.allocate(dom);
        // calculate dz
        Kokkos::parallel_for( dom.nCellMem , KOKKOS_LAMBDA (int iGlob) {
          int i, j, k, iGlobSW;

          unpackIndices(iGlob, dom.nz+2*hc, dom.ny+2*hc, dom.nx+2*hc, k, j, i);
          iGlobSW = packIndices(dom.ny+2*hc, dom.nx+2*hc, j, i);

          // only calculate dz for interior cells (not halo cells)
          if(i>hc-1 && i<dom.nx+hc && j>hc-1 && j<dom.ny+hc && k>hc-1 && k<dom.nz+hc)
          {
              if (state.z(iGlobSW) <= dom.bottomZ)
              {
                  std::cerr<< RERROR "Bottom of subsurface domain must be lower than DEM!\n";
              }
              // As of 20200301, dz is assumed uniform in vertical directions
              // For implementing non-uniform dz in the future, these lines need
              // to be changed.
              statesub.dz(iGlob) = (state.z(iGlobSW) - dom.bottomZ) / dom.nz_glob;
              statesub.z(iGlob) = state.z(iGlobSW) - k*statesub.dz(iGlob);
          }
    	});

    }
    #endif


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
			unpackIndicesUniformGrid(ii,dom.ny+2*hc,dom.nx+2*hc,j,i);

			if(i==dom.nx+hc-1&&dom.iE){
				state.z(ii+1)=state.z(ii);
			}
			if(i==hc&&dom.iW) {
				state.z(ii-1)=state.z(ii);
			}
			if(j==dom.ny+hc-1&&dom.iS){
				state.z(ii+dom.nx+2*hc)=state.z(ii);
			}
			if(j==hc&&dom.iN){
				state.z(ii-(dom.nx+2*hc))=state.z(ii);
			}

		});

	}


	inline void boundaryReflectiveW(State &state, Domain &dom){

		Kokkos::parallel_for("boundaryReflectiveW", dom.ny*hc , KOKKOS_LAMBDA (int iGlob) {
			int rx, ry;
			unpackIndicesUniformGrid(iGlob,dom.ny,hc,ry,rx);
			int ii=(hc+ry)*(dom.nx+2*hc)+hc;
			//west boundary
			state.z(ii-rx-1)	=1e4; //10000m high
			state.h(ii-rx-1)	=0.0;
			state.hu(ii-rx-1)	=0.0;
			state.hv(ii-rx-1)	=0.0;
		 });
	}

	inline void boundaryReflectiveE(State &state,Domain &dom){

		Kokkos::parallel_for("boundaryReflectiveE", dom.ny*hc , KOKKOS_LAMBDA (int iGlob) {
			int rx, ry;
			unpackIndicesUniformGrid(iGlob,dom.ny,hc,ry,rx);
			int ii=(hc+ry)*(dom.nx+2*hc)+hc+dom.nx-1;
			//east boundary
			state.z(ii+rx+1) 	= 1e4; //10000m high
			state.h(ii+rx+1)	=0.0;
			state.hu(ii+rx+1)	=0.0;
			state.hv(ii+rx+1)	=0.0;

		 });
	}

	inline void boundaryReflectiveN(State &state,Domain &dom){

		Kokkos::parallel_for("boundaryReflectiveN", hc*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      	int rx, ry;
      	unpackIndicesUniformGrid(iGlob,hc,dom.nx,ry,rx);
			int ii=hc*(dom.nx+2*hc)+hc+rx;
			//north boundary
			state.z(ii-(ry+1)*(dom.nx+2*hc))	=1e4; //10000m high
			state.h(ii-(ry+1)*(dom.nx+2*hc))	=0.0;
			state.hu(ii-(ry+1)*(dom.nx+2*hc))=0.0;
			state.hv(ii-(ry+1)*(dom.nx+2*hc))=0.0;

		 });
	}

	inline void boundaryReflectiveS(State &state,Domain &dom){

		Kokkos::parallel_for("boundaryReflectiveS", hc*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      	int rx, ry;
      	unpackIndicesUniformGrid(iGlob,hc,dom.nx,ry,rx);
			int ii=(hc+dom.ny-1)*(dom.nx+2*hc)+hc+rx;
			//south boundary
			state.z(ii+(ry+1)*(dom.nx+2*hc))	=1e4; //10000m high
			state.h(ii+(ry+1)*(dom.nx+2*hc))	=0.0;
			state.hu(ii+(ry+1)*(dom.nx+2*hc))=0.0;
			state.hv(ii+(ry+1)*(dom.nx+2*hc))=0.0;

		 });
	}

};

#endif
