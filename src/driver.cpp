/* -*- mode: c++; c-default-style: "linux" -*- */

#include "stdlib.h"
#include <iostream>
#include <string>
#include "define.h"
#include "Domain.h"
#include "Parallel.h"
#include "Parser.h"
#include "State.h"
#include "BC.h"
#include "Initializer.h"
#include "TimeIntegrator.h"
#include "FileIO.h"
#include "Exchange.h"
#include "SWSourceSink.h"
#include "DomainIntegrator.h"
#include "Vegetation.h"
#include "ParticleTracking.h"
#include "tools.h"

#if SERGHEI_SUBSURFACE_MODEL
#include "GwDomain.h"
#include "GwFunction.h"
#include "GwInit.h"
#include "GwMPI.h"
#include "GwState.h"
#endif

int main(int argc, char** argv) {
	#if SERGHEI_DEBUG_WORKFLOW
	std::cerr << GGD "Start" << std::endl;
	#endif


	// InitArguments deprecated in Kokkos v4.0.0
	// Kokkos::InitArguments args;
	Kokkos::InitializationSettings kokkosSettings;

	//these variables are needed to call  initializeMPI
	Parallel       par;
	Initializer    init;

	//Read input and create output
	if (argc != 4){
		std::cerr << RERROR "The program is run as: ./nprogram inputFolder/ outputFolder/ Nthreads" << "\n";
		return 0;
	}
	std::string inFolder = argv[1];
	std::string outFolder = argv[2];
	par.nthreads=atoi(argv[3]);

	#if SERGHEI_DEBUG_WORKFLOW
	std::cerr << GGD "Initialising MPI AND CUDA" << std::endl;
	#endif

	// Initialize MPI According to the documentation, MPI_Init should be called before Kokkos::initialize
	init.initializeMPI( &argc , &argv , par );

	#ifdef __NVCC__
		kokkosSettings.set_device_id(par.myrank%par.nthreads);
	#else
		if(par.nthreads!=0){
			kokkosSettings.set_num_threads(par.nthreads);
		}
	#endif
	#if SERGHEI_DEBUG_KOKKOS_SETUP
		printKokkosInitArguments(args,par);
		#if __NVCC__
			printKokkosCuda(args,par);
		#endif
	#endif

	#if SERGHEI_DEBUG_WORKFLOW
	std::cerr << GGD "Initialising Kokoks - rank " << par.myrank << std::endl;
	#endif
	// INITIALIZE KOKKOS
	Kokkos::initialize(kokkosSettings);


	#if SERGHEI_DEBUG_WORKFLOW
		std::cerr << GGD "Program instantiated, creating objects - rank " << par.myrank << std::endl;
	#endif

	{ //these scope guards are needed to avoid annoying warnings

		// Create the model objects, these are created one for each MPI rank (or subdomain)
		State          state;
		SourceSinkData ss;
		ExternalBoundaries ebc;
		Domain         dom;
		Parser         parser;
		FileIO         io;
		Exchange       exch;
		TimeIntegrator tint;
		surfaceIntegrator sint;
		boundaryIntegrator bint;
		#if SERGHEI_PARTICLE_TRACKING
		ParticleTracker parTrack;
		#endif
		#if SERGHEI_TOOLS
		Observations obs;
		#endif
		#if SERGHEI_SUBSURFACE_MODEL
		GwState gw;
		GwFunction gwf;
		GwBC gbc;
		GwDomain gdom;
		GwInit ginit;
		GwMPI gmpi;
		#endif
		double oldVolume,newVolume, diffVolume;
		double accumDt=0.0;
    	Kokkos::Timer timer; // auxiliary timer object

		// Initialize the model
		if(!init.initialize(state, ss, ebc, dom, par, tint, sint, bint, parser, exch, io, inFolder, outFolder)){
			std::cerr << RERROR "Unable to start the simulation" << "\n";
			return 0;
		};

		// Initialize subsurface model if activated
		#if SERGHEI_SUBSURFACE_MODEL
		if (!ginit.initialize_gw(gw, gdom, state, dom, gbc, gmpi, par, io, ss, inFolder, outFolder)) {
			std::cerr << RERROR "Unable to initialize the subsurface domain" << "\n"; return 0;
		};
		GwMatrix A(gdom);
		GwSolver gsolver;
		gsolver.init(A, gdom);
		if( par.masterproc){std::cerr << GOK "Subsurface Solver has been initialized! " << std::endl;}
		#endif

		#if SERGHEI_TOOLS
  		if(!obs.readInputFiles(inFolder,par)) return 0;
		if(!obs.configure(dom,outFolder)) return 0;	// observations for surface domain
		//obs.printGauges(dom);
		obs.update(state,par,dom);

		if( par.masterproc){
			obs.writeLinesSamplingCoordinates(outFolder);
			obs.writeGauges(dom.etime);
			obs.writeLines(dom.etime);
		}
		#endif

		#if SERGHEI_DEBUG_WORKFLOW

		for (int k = 0; k < ebc.extbc.size(); k ++) {
		  std::cout << GGD << GRAY << __FILE__ << ":" << __LINE__ << RESET << "\tExtBC[" << k << "]: " << ebc.extbc[k].bcvals(0) << ", " << ebc.extbc[k].bcvals(1) << ", " << ebc.extbc[k].bcvals(2) << std::endl;
		}

		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "Initialisation finished, starting to run main loop" << std::endl;
		for(int i = 0; i < ebc.extbc.size(); i ++) {
		  std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "ncellsBC for segment " << i << ": " << ebc.extbc[i].ncellsBC << "\n";
		}
		#if SERGHEI_DEBUG_BOUNDARY
		bint.integrate(ebc.extbc,dom,1);
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "ncellsBC (integrated) " << bint.ncellsBC << "\n";
		std::cout << GGD<< GRAY << __PRETTY_FUNCTION__ << RESET <<  "outflow discharge (integrated) " << bint.outflowDischarge << "\n";
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow accumulated (integrated) " << bint.outflowAccumulated << std::endl;
		#endif
		#endif

		//integrator at the beginning or the simulation
		sint.integrate(state,dom,ss);

		// Write initial time series data
		io.writeTimeSeriesIni(state,dom,par,ss,sint,bint,ebc.extbc,outFolder);
		#if SERGHEI_SUBSURFACE_MODEL
		io.writeSubTimeSeriesIni(gw,gdom,dom,par,outFolder);
		#endif

		// capture initialisation time
		dom.timers.init = timer.seconds();

		if (par.masterproc) std::cerr << "\n" << GOK "SIMULATION STARTS\n";

		#if SERGHEI_SUBSURFACE_MODEL
		gdom.dt = gdom.dt_init;
		gdom.dtOld = gdom.dt_init;
		dom.dt = gdom.dt;
		#else
		tint.computeDt(state,dom,io);
		#endif

		// The Main Time Loop
		while (dom.etime < dom.simLength) {
			//previous mass
			oldVolume=sint.surfaceVolumeG;

			bint.integrate(ebc.extbc,dom,1);//has to be called here (previous time step) with mode==1 (boundary flows)

			tint.stepForward(state, ss, ebc.extbc, dom, exch, par, io);

			// run subsurface model
			#if SERGHEI_SUBSURFACE_MODEL
			// rainfall
	        if (dom.isRain) {
	            Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA(int idom) {gdom.qrain(idom) = ss.rainRate(idom); });
				// apply rainfall to the surface domain
		        Kokkos::parallel_for( dom.nCellDomain , KOKKOS_LAMBDA (int idom) {
		            int iGlob = dom.getIndex(idom);
		            state.h(iGlob) += ss.rainRate(iGlob)*dom.dt;
		        });
	        }
			// Asynchronous coupling
			if (gdom.async)	{
				if (gdom.etime + gdom.dt < dom.etime)	{
					std::cerr << "     Asynchrnous coupling, execute GW at dt = " << gdom.dt <<"\n";
					gdom.etime += gdom.dt;
					if (gdom.gw_scheme == 1)	{gwf.pca_solve(gw, state, gdom, gbc, A, gsolver, ss, gmpi, par);}
					else {gwf.picard_solve(gw, state, gdom, gbc, A, gsolver, ss, gmpi, par);}
				}
			}
			else {
				gdom.etime = dom.etime;
				if (gdom.gw_scheme == 1)	{gwf.pca_solve(gw, state, gdom, gbc, A, gsolver, ss, gmpi, par);}
				else {gwf.picard_solve(gw, state, gdom, gbc, A, gsolver, ss, gmpi, par);}
			}
			tint.computeGwExchange(state , dom);
			#endif

			// Unify dt
			tint.computeDt(state,dom,io);
			#if SERGHEI_SUBSURFACE_MODEL
			if (!gdom.async)	{
				if (dom.dt < gdom.dt)	{gdom.dt = dom.dt;}
				else {dom.dt = gdom.dt;}
			}
			else if (dom.dt > gdom.dt)	{dom.dt = gdom.dt;}
			#endif

			oldVolume+=(bint.inflowDischargeG - bint.outflowDischargeG)*dom.dt; //Boundary fluxes with the new dt

			bint.integrate(ebc.extbc,dom,0);//called here with mode==0 (adjusted volume)

			oldVolume+=bint.adjustedVolumeG; //Some mass changes can occur through the boundaries

			sint.integrate(state,dom,ss); //new mass after the new time step integration
			oldVolume+= (sint.rainFluxG-sint.infFluxG)*dom.dt; //after integrate, we have to sum the rain and inf mass


			//new mass
			newVolume=sint.surfaceVolumeG;

			if(fabs(oldVolume)>TOL12){
				diffVolume=(newVolume-oldVolume)/oldVolume*100.;
			}else{
				diffVolume=0.0;
			}
			dom.etime += dom.dt;
			dom.nIter++;
			dom.countIterDt++;
			accumDt+=dom.dt;

			if (dom.nIter%io.nScreen==0 || fabs(dom.etime - io.numOut*io.outFreq) <= 0.5*dom.dt) {
			//if (fabs(dom.etime - io.numOut*io.outFreq) <= 0.5*dom.dt) {
				if (par.masterproc) {
					std::cerr << std::fixed;
					std::cerr << GSTAR "TIME: " << dom.etime << " average dt: " << accumDt/dom.countIterDt <<"\n";
					std::cerr.precision(9);
					std::cerr << std::scientific;
					std::cerr << "     Diff Volume: " << diffVolume <<"\n";
					std::cerr << std::fixed;
					std::cerr.precision(12);
					std::cerr << "     Inflow Discharge: " << bint.inflowDischargeG <<"\n";
					std::cerr << "     Outflow Discharge: " << bint.outflowDischargeG <<"\n";
					if(fabs(diffVolume)>TOL_MASS_ERROR){
						std::cerr << YEXC "   Rain Volume:\t" << sint.rainFluxG*dom.dt <<"\n";
						#if SERGHEI_DEBUG_MASS_CONS > 1
                            getchar();
                        #endif
					}
				}
				// write output
				if (fabs(dom.etime - io.numOut*io.outFreq) <= 0.5*dom.dt) {
					io.output(state, dom, ss, par,outFolder);
					#if SERGHEI_SUBSURFACE_MODEL
					io.outputSubsurface(gw, gdom, par,outFolder);
					#endif
					if(par.masterproc) std::cerr << GIO "File " << io.numOut-1 << " written" <<"\n"; //io.numOut already updated
				}

				// real nout = io.numOut;
				// int ierr = MPI_Bcast(&io.numOut, 1, MPI_INT, 0, MPI_COMM_WORLD);
				// if(par.masterproc) dom.timers.Tout+=dom.timers.out.seconds();

				if(par.masterproc) std::cerr << "-------------------------------------------------\n";
				dom.countIterDt=0;
				accumDt=0.0;

			}

			#if SERGHEI_PARTICLE_TRACKING
			parTrack.update(dom,state);
			#endif

			if (dom.etime >= io.numObs*io.obsFreq) {
				#if SERGHEI_TOOLS
				obs.update(state,par,dom);
				#endif
				io.writeTimeSeries(state,dom,par,sint,bint);
				#if SERGHEI_SUBSURFACE_MODEL
				io.writeSubsurfaceTimeSeries(gw,gdom,dom,par);
				#endif
			  if (par.masterproc){
					#if SERGHEI_TOOLS
          obs.write(dom);
					#endif
			  }

			}
		} 		// end of time loop

		dom.timers.total = timer.seconds();
		if (par.masterproc){
			std::cerr << GOK "SIMULATION FINISHED\n";
			std::cerr << GOK "Time elapsed: " << dom.timers.total << "\n";
		}
		io.writeLogFile(dom,par,outFolder);

	io.closeOutputStreams();
	#if SERGHEI_TOOLS
	if(par.masterproc) obs.closeOutputStreams();
	#endif
	}


	Kokkos::finalize();

	MPI_Finalize();

}
