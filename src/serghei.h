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

class SERGHEI{
public:
	Parallel       par;
	Initializer    init;

	State               state;
	Domain              dom;
	FileIO              io;

 private:
	SourceSinkData      ss;
	ExternalBoundaries  ebc;
	Parser              parser;
	Exchange            exch;
	TimeIntegrator      tint;
	surfaceIntegrator   sint;
	boundaryIntegrator  bint;
  #if SERGHEI_TOOLS
	  Observations obs;
	#endif

  #if SERGHEI_PARTICLE_TRACKING
		ParticleTracker parTrack;
	#endif

	double oldVolume,newVolume, diffVolume;
	double accumDt=0.0;

  // Kokkos objects
  Kokkos::Timer timer;
  Kokkos::InitializationSettings kokkosSettings;

public:
	std::string inFolder, outFolder;

////////////// METHODS ///////////////
public:

  int start(int argc, char **argv){
    #if SERGHEI_DEBUG_WORKFLOW
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	  #endif


    init.initializeMPI( &argc , &argv , par );

    #ifdef __NVCC__
		  kokkosSettings.set_device_id(par.myrank%par.nthreads);
	  #else
		  if(par.nthreads!=0) kokkosSettings.set_num_threads(par.nthreads);
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

    Kokkos::initialize(kokkosSettings);

	  #if SERGHEI_DEBUG_WORKFLOW
		  std::cerr << GGD "Program instantiated, creating objects - rank " << par.myrank << std::endl;
	  #endif

    // Initialize the model
		if(!init.initialize(state, ss, ebc, dom, par, tint, sint, bint, parser, exch, io, inFolder, outFolder)){
			std::cerr << RERROR "Unable to start the simulation" << "\n";
			return 0;
		};


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

		// capture initialisation time
		dom.timers.init = timer.seconds();

		std::cout << GOK << "Initialisation complete. Initialisation time: " << dom.timers.init << " [s]" << std::endl;
    return 1;
  }

  int compute(){
    #if SERGHEI_DEBUG_WORKFLOW
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	  #endif
		if (par.masterproc){
			std::cout << std::endl << GOK "SIMULATION STARTS" << std::endl;
			std::cout << BDASH << "Start time: " << dom.startTime << std::endl;
			std::cout << BDASH << "End time: " << dom.endTime << std::endl;
		}

		while (dom.etime < dom.endTime) {

			//previous mass
			oldVolume=sint.surfaceVolumeG;

			bint.integrate(ebc.extbc,dom,1);//has to be called here (previous time step) with mode==1 (boundary flows)

			tint.stepForward(state, ss, ebc.extbc, dom, exch, par, io);

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

			if (dom.nIter%io.nScreen==0 || fabs(dom.etime - dom.startTime - io.numOut*io.outFreq) < TOL12) {
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
				if(fabs(dom.etime - dom.startTime - io.numOut*io.outFreq) < TOL12){
					io.output(state, dom, ss, par,outFolder);
					if(par.masterproc) std::cerr << GIO "File " << io.numOut-1 << " written" << std::endl; //io.numOut already updated
				}
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

				io.writeTimeSeries(state,dom,par,sint,bint,ebc.extbc);
			  if (par.masterproc){
					#if SERGHEI_TOOLS
          obs.write(dom);
					#endif
			  }

			}

		} 		// end of time loop

    return 1;
  }

  int finalise(){
    #if SERGHEI_DEBUG_WORKFLOW
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	  #endif

    dom.timers.total = timer.seconds();
		if (par.masterproc){
			std::cerr << GOK "SIMULATION FINISHED\n";
			std::cerr << GOK "Time elapsed: " << dom.timers.total << std::endl;
		}
		io.writeLogFile(dom,par,outFolder);

	  io.closeOutputStreams();
	  #if SERGHEI_TOOLS
	    if(par.masterproc) obs.closeOutputStreams();
	  #endif


    return 1;

  }


};
