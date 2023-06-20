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

  int initialize(State &state, SourceSinkData &ss, std::vector<ExtBC> &extbc, Domain &dom, Parallel &par, TimeIntegrator &tint, surfaceIntegrator &sint, boundaryIntegrator &bint, Parser &parser, Exchange &exch, FileIO &io, std::string inFolder, std::string outFolder) {

    if(!parser.readDimensions(inFolder, dom, state, par, io)){
      return 0;
    }

    //buildDomainDecomposition(state, dom, par);
    dom.buildDomainDecomposition(par);

    dom.initialise();
    state.allocate(dom);

    if(!parser.readInputFiles(inFolder, dom, state, ss, extbc, par, io)){
      return 0;
    }

    ss.allocate(dom);
    sint.initialize(state,dom,ss);
    bint.initialize(extbc);

    if(par.masterproc){
      if(!parser.createOutputDir(outFolder)){
	return 0;
      }
    }

    exch.exchangeIniMPI(state,dom,exch,par);

    //initialize boundaries
    boundaryIni(state,dom,par,extbc);

    // Output the initial model state
    io.outputIni(state, dom, ss, par,outFolder);

    return 1;

  }
// #if SERGHEI_SUBSURFACE_MODEL
//     // initialize subsurface domain, ZhiLi20210219
//   int initializeSubsurface(SubsurfaceState &statesub, State &state, SubBC &subbc, DomainSubsurface &domsub, Domain &dom, Parallel &par, TimeIntegrator &tint, Parser &parser, FileIO &io, std::string inFolder, std::string outFolder)
//     {
//       int flag = -1;
//
//       // read the subsurface domain information
//       if(!parser.readSubsurfaceDimensions(inFolder, domsub, par))
// 	{
// 	  std::cerr << GOK << " Reading in subsurface dimensions failed." << std::endl;
// 	  flag = 0;
// 	  return flag;
// 	}
//
//       // As of 2020-03-01 : Assume nz = nz_glob (no vertical decomposition)
//       domsub.nz = domsub.nz_glob;
//       domsub.nx = dom.nx;
//       domsub.ny = dom.ny;
//       domsub.xll = dom.xll;
//       domsub.yll = dom.yll;
//       domsub.zll = domsub.maxZ;
//       // Assume serial
//       domsub.nx_glob = domsub.nx;
//       domsub.ny_glob = domsub.ny;
//       domsub.dx = dom.dx;
//
//       if (dom.isRain)   {domsub.isRain = 1;}
//
//       // allocate subsurface domain
//       allocateDomainSubsurface(statesub, state, domsub, tint);
//       /* now the state is allocated and array dimensions are known */
//       /* we can proceed to read the van Genuchten soil parameters */
//       if(!parser.readVGParameters(inFolder, statesub, domsub, par))
// 	{
// 	  std::cerr << RERROR << " Reading in van Genuchten parameters failed." << std::endl;
// 	  flag = 0;
// 	  return flag;
// 	}
//     // read boundary conditions
//     if(!parser.readSubBCFile(inFolder, domsub, subbc, par, statesub))    return 0;
//     // read input / initial conditions for the subsurface domain
//     if(!parser.readSubsurfaceState(inFolder, domsub, statesub, state, par, io, subbc)) return 0;
//
//       // ignore source/sink and domain integrator for now
//
//       flag = 1;
//       #if SERGHEI_SUBSURFACE_MODEL
//       io.outputSubsurface(statesub, domsub, par,outFolder);
//       #endif
//
//       return flag;
//     }
//   #endif
//
//     #if SERGHEI_SUBSURFACE_MODEL
//     void allocateDomainSubsurface(SubsurfaceState &statesub, State &state, DomainSubsurface &dom, TimeIntegrator &tint)
//     {
//         // Initialize the time
//         dom.etime = 0;
//         dom.nIter = 0;
//         dom.n_iter_gmres = 0;
//         dom.countIterDt=0;
//         dom.ncells=(dom.ny+2*haloc)*(dom.nx+2*haloc)*(dom.nz+2*haloc); //cells plus halo cells, haloc is number of overlapping halo cells defined in const.h
//         dom.nCellDomain = dom.nx*dom.ny*dom.nz; // WARNING UCM
//         dom.isnodata = boolArr ("isnodata", dom.ncells);
//
// 	    // We need the global max elevation to calculate vertical discretization
//         // Not sure if there are existing functions to do this.
//         // int dom.maxZ = 2.0;
//
//         //allocate the variables in the subsurface domain
//         statesub.allocate(dom);
//         // calculate dz
//         Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA (int iGlob) {
//           int i, j, k, iGlobSW, nxhalo, nyhalo;
//           real dist;
//           nxhalo = dom.nx+2*haloc;
//           nyhalo = dom.ny+2*haloc;
//           unpackIndices(iGlob, dom.nz+2*haloc, nyhalo, nxhalo, k, j, i);
//           iGlobSW = packIndices(nyhalo, nxhalo, j, i);
//           dom.isnodata(iGlob) = 0;
//           // first calculate dz for interior cells (not halo cells)
//           if(i>haloc-1 && i<dom.nx+haloc && j>haloc-1 && j<dom.ny+haloc && k>haloc-1 && k<dom.nz+haloc)
//           {
//               if (state.z(iGlobSW) <= dom.bottomZ) {
//                   std::cerr<< RERROR "Bottom of subsurface domain must be lower than DEM!\n";
//               }
//               #if FOLLOW_TERRAIN
//               // ----- Terrain-following mesh by Zhi Li 20220307 -----
//               //          Uniform dz within 1 column for now
//               dom.isnodata(iGlob) = 0;
//               if (k == 1)   {
//                   state.qss(iGlobSW) = 0.0;
//                   statesub.ktop(iGlobSW) = 1;
//               }
//               // get dz
//               statesub.dz(iGlob) = (state.z(iGlobSW) - dom.bottomZ) / dom.nz_glob;
//               // get center of the grid
//               statesub.z(iGlob) = state.z(iGlobSW) - (k-0.5)*statesub.dz(iGlob);
//               #else
//               // Old vertical discretization following Frehg
//               statesub.dz(iGlob) = (dom.maxZ - dom.bottomZ) / dom.nz_glob;
//               statesub.z(iGlob) = dom.maxZ - k*statesub.dz(iGlob);
//               if (statesub.z(iGlob) < state.z(iGlobSW)-1e-7){
//                   dom.isnodata(iGlob) = 0;
//                   if (statesub.z(iGlob) >= state.z(iGlobSW)-statesub.dz(iGlob)-1e-7) {
//                       dist = state.z(iGlobSW) - statesub.z(iGlob);
//                       state.qss(iGlobSW) = 0.0;
//                       statesub.ktop(iGlobSW) = k;
//                       statesub.dz(iGlob) = dist;
//                   }
//               }
//               #endif
//           }
//     	});
//
//         // get the angles
//         #if FOLLOW_TERRAIN
//         Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA (int iGlob) {
//           int i, j, k, nxhalo;
//           real hdiff, dist;
//           nxhalo = dom.nx+2*haloc;
//           unpackIndices(iGlob, dom.nz+2*haloc, dom.ny+2*haloc, dom.nx+2*haloc, k, j, i);
//           if (j>haloc-1 && j<dom.ny+haloc && k>haloc-1 && k<dom.nz+haloc) {
//               // x
//               if (i == haloc-1 | i == dom.nx) {
//                   statesub.sinx(iGlob) = 0.0;
//                   statesub.cosx(iGlob) = 1.0;
//               }
//               else  {
//                   hdiff = myfabs(statesub.z(iGlob+1) - statesub.z(iGlob));
//                   dist = mysqrt(mypow(hdiff, 2.0) + mypow(dom.dx, 2.0));
//                   statesub.sinx(iGlob) = hdiff / dist;
//                   statesub.cosx(iGlob) = dom.dx / dist;
//               }
//               // y
//               if (j == haloc-1 | j == dom.ny) {
//                   statesub.siny(iGlob) = 0.0;
//                   statesub.cosy(iGlob) = 1.0;
//               }
//               else  {
//                   hdiff = myfabs(statesub.z(iGlob+nxhalo) - statesub.z(iGlob));
//                   dist = mysqrt(mypow(hdiff, 2.0) + mypow(dom.dx, 2.0));
//                   statesub.siny(iGlob) = hdiff / dist;
//                   statesub.cosy(iGlob) = dom.dx / dist;
//               }
//           }
//         });
//         #else
//         Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA (int iGlob) {
//           int i, j, k;
//           unpackIndices(iGlob, dom.nz+2*haloc, dom.ny+2*haloc, dom.nx+2*haloc, k, j, i);
//           // if (j>haloc-1 && j<dom.ny+haloc && k>haloc-1 && k<dom.nz+haloc) {
//           //     statesub.sinx(iGlob) = 0.0;
//           //     statesub.siny(iGlob) = 0.0;
//           //     statesub.cosx(iGlob) = 1.0;
//           //     statesub.cosy(iGlob) = 1.0;
//           // }
//           statesub.sinx(iGlob) = 0.0;
//           statesub.siny(iGlob) = 0.0;
//           statesub.cosx(iGlob) = 1.0;
//           statesub.cosy(iGlob) = 1.0;
//         });
//         #endif
//
//
//         // then calculate dz for halo cells
//         Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA (int iGlob) {
//           int i, j, k;
//           unpackIndices(iGlob, dom.nz+2*haloc, dom.ny+2*haloc, dom.nx+2*haloc, k, j, i);
//           if (j>haloc-1 && j<dom.ny+haloc && k>haloc-1 && k<dom.nz+haloc) {
//               if (i == haloc-1) {
//                   statesub.dz(iGlob) = statesub.dz(iGlob+1);
//                   statesub.z(iGlob) = statesub.z(iGlob+1);
//               }
//               else if (i == dom.nx+haloc) {
//                   statesub.dz(iGlob) = statesub.dz(iGlob-1);
//                   statesub.z(iGlob) = statesub.z(iGlob-1);
//               }
//           }
//         });
//
//         Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA (int iGlob) {
//           int i, j, k;
//           unpackIndices(iGlob, dom.nz+2*haloc, dom.ny+2*haloc, dom.nx+2*haloc, k, j, i);
//           if (k>haloc-1 && k<dom.nz+haloc) {
//               if (j == haloc-1) {
//                   statesub.dz(iGlob) = statesub.dz(iGlob+dom.nx+2*haloc);
//                   statesub.z(iGlob) = statesub.z(iGlob+dom.nx+2*haloc);
//               }
//               else if (j == dom.ny+haloc) {
//                   statesub.dz(iGlob) = statesub.dz(iGlob-dom.nx-2*haloc);
//                   statesub.z(iGlob) = statesub.z(iGlob-dom.nx-2*haloc);
//               }
//           }
//         });
//
//         Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA (int iGlob) {
//           int i, j, k, nxny, iGlobSW;
//           nxny = (dom.nx+2*haloc) * (dom.ny+2*haloc);
//           unpackIndices(iGlob, dom.nz+2*haloc, dom.ny+2*haloc, dom.nx+2*haloc, k, j, i);
//           if (k == haloc-1) {
//               statesub.dz(iGlob) = statesub.dz(iGlob+nxny);
//               statesub.z(iGlob) = statesub.z(iGlob+nxny) + statesub.dz(iGlob);
//           }
//           else if (k == dom.nz+haloc) {
//               statesub.dz(iGlob) = statesub.dz(iGlob-nxny);
//               statesub.z(iGlob) = statesub.z(iGlob-nxny) - statesub.dz(iGlob);
//               iGlobSW = packIndices(dom.ny+2*haloc, dom.nx+2*haloc, j, i);
//               if (i == haloc-1) {
//                   statesub.ktop(iGlobSW) = statesub.ktop(iGlobSW+1);
//               }
//               else if (i == dom.nx+haloc) {
//                   statesub.ktop(iGlobSW) = statesub.ktop(iGlobSW-1);
//               }
//               if (j == haloc-1) {
//                   statesub.ktop(iGlobSW) = statesub.ktop(iGlobSW+dom.nx+2*haloc);
//               }
//               else if (j == dom.ny+haloc) {
//                   statesub.ktop(iGlobSW) = statesub.ktop(iGlobSW-dom.nx-2*haloc);
//               }
//           }
//         });
//
//     }
//     #endif


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

		Kokkos::parallel_for( extbc.ncellsBC, KOKKOS_LAMBDA (int iGlob) {
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

		Kokkos::parallel_for( dom.ny*haloc , KOKKOS_LAMBDA (int iGlob) {
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

		Kokkos::parallel_for( dom.ny*haloc , KOKKOS_LAMBDA (int iGlob) {
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

		Kokkos::parallel_for( haloc*dom.nx , KOKKOS_LAMBDA (int iGlob) {
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

		Kokkos::parallel_for( haloc*dom.nx , KOKKOS_LAMBDA (int iGlob) {
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
