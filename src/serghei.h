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
#include "SourceSink.h"
#include "DomainIntegrator.h"
#include "Vegetation.h"
#include "ParticleTracking.h"
#include "tools.h"

#if SERGHEI_SUBSURFACE_MODEL
#include "GwDomain.h"
#include "GwFunction.h"
#include "GwInit.h"
#include "GwMPI.h"
#include "GwMatrix.h"
#include "GwState.h"
#include "GwSolver.h"
#include "GwIntegrator.h"
#endif
#if SERGHEI_DRAINAGE_MODEL
#define EXTT extern

#include "Drainage/globals.h"                  
#include "Drainage/funcs.h"          
#include "Drainage/DrainageDywave.h"
#endif

class SERGHEI{
public:
	Parallel       par;
	Initializer    init;
	FileIO				io;

	State               state;
	Domain              dom;
	#if SERGHEI_SUBSURFACE_MODEL
	GwState gw;
	GwDomain gdom;
	GwInit ginit;
	GwMPI gmpi;
	GwFunction gwf;
	SubsurfaceBoundaries	gbc;
	GwMatrix A;
	GwIntegrator gint;
	#ifdef __NVCC__
	GwSolver<Kokkos::Cuda> gsolver;
	#else
	GwSolver<Kokkos::OpenMP> gsolver;
	#endif
	#endif
	#if SERGHEI_DRAINAGE_MODEL
	ReadFileFuncs RFF;
	DrainageDywave Ddyw;
	Node Tnode;
	Link Tlink;
	Outfall Toutfall;
	Storage Tstorage;
	Conduit Tconduit;
	Pump Tpump;
	PumpCurves TpumpCurves;
	RiverStages TRiver;
	#endif

	SourceSink      ss;
	//SourceSinkData      ss;

  Exchange            exch;
	ExternalBoundaries  ebc;

 private:
	Parser              parser;
	TimeIntegrator      tint;
	surfaceIntegrator   sint;
	boundaryIntegrator  bint;
	#if SERGHEI_TOOLS
		Observations obs;
	#endif
		
  #if SERGHEI_LPT
	ParticleTracker parTrack;
  #endif

	double oldVolume,newVolume, diffVolume;
	double accumDt=0.0;
  
  // Kokkos objects
  Kokkos::Timer timer;
  Kokkos::Timer timer_particles; 
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
		  printKokkosInitArguments(par);
		  #if __NVCC__
			  printKokkosCuda(args,par);
		  #endif
	  #endif

		#if SERGHEI_DEBUG_WORKFLOW
			std::cerr << GGD "Initialising Kokkos - rank " << par.myrank << std::endl;
		#endif

		Kokkos::initialize(kokkosSettings);

		#if SERGHEI_DEBUG_WORKFLOW
			std::cerr << GGD "Program instantiated, creating objects - rank " << par.myrank << std::endl;
		#endif

		// Initialize the model
		if(!init.initialize(state, ss.swss, ebc, dom, par, tint, sint, bint, parser, exch, io, inFolder, outFolder)){
			std::cerr << RERROR "Unable to start the simulation" << "\n";
			return 0;
		};
		#if SERGHEI_LPT
		  if(!parser.readParticles(inFolder,par,&parTrack)){
			std::cerr << RERROR "Unable to start the simulation because of LPT initialization files" << "\n";
		    return 0;
		  }
	      parTrack.initialiseParticles(dom, state);

	      if(io.outFormat==OUT_VTK){
			#if SERGHEI_PARTICLE_NO_OUTPUT
			#else
	                 io.outputIniParticle(dom, parTrack, outFolder);
			#endif
		  }else{
			io.outputInitParticlesNETCDF(parTrack,dom,par,outFolder);
		  }
		  std::cerr << GOK "LPT module has been initialized! Number of particles: " << parTrack.N_par << "\n" << std::endl;
		#endif


		// Initialize subsurface model if activated
		#if SERGHEI_SUBSURFACE_MODEL
		if (!ginit.initialize_gw(gw, gdom, state, dom, gbc, gmpi, gint, par, io, ss, inFolder, outFolder)) {
			std::cerr << RERROR "Unable to initialize the subsurface domain" << "\n"; return 0;
		};
		A.init(gdom);
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
		sint.integrate(state,dom,ss.swss);
		// Write initial time series data
		io.writeTimeSeriesIni(state,dom,par,ss.swss,sint,bint,ebc.extbc,outFolder);
		#if SERGHEI_SUBSURFACE_MODEL
		io.writeSubTimeSeriesIni(gdom, gint, par, outFolder);
		io.outputSubsurface(gw, gdom, par, outFolder);
		#endif
		// capture initialisation time
		dom.timers.init = timer.seconds();

		if(par.masterproc) std::cout << GOK << "Initialisation complete. Initialisation time: " << dom.timers.init << " [s]" << std::endl;
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
		#if SERGHEI_SUBSURFACE_MODEL
		gdom.dt = gdom.dt_init;
		gdom.dtOld = gdom.dt_init;
		dom.dt = gdom.dt;
		gdom.cg_iter = 0;
		#else
		tint.computeDt(state,dom,io);
		#endif
		#if SERGHEI_DRAINAGE_MODEL
		Ddyw.routingStep = InitRouteStep;
		// Read drainage parameter file
		std::string drainageParamFile = inFolder + "Drainage-parameter.input";
		if(!parser.readDrainageParameterFile(drainageParamFile, Ddyw, par)) return 0;
		RFF.project_readInput(inFolder, Tnode, Tlink, Tconduit, Toutfall, Tstorage, Tpump, TpumpCurves, TRiver); 
		
		RFF.project_validate(Tnode, Tlink, Tconduit, Toutfall, Tstorage, Tpump, TpumpCurves);
		Ddyw.flowrout_init(Tnode, Tlink, Tconduit, TRiver, dom.etime);
		io.writeDrainageTimeSeriesIni(par, dom, sint, Ddyw, outFolder, RFF, Tlink, Tnode);
		std::cerr << GOK "The DRAINAGE input information has been read.\n";
		#endif

		// Main Time Loop
		while (dom.etime < dom.endTime) {
			//previous mass
			oldVolume=sint.surfaceVolumeG;
			bint.integrate(ebc.extbc,dom,1);//has to be called here (previous time step) with mode==1 (boundary flows)
			// run surface model
			#if SERGHEI_SWE_MODEL
			tint.stepForward(state, ss.swss, ebc.extbc, dom, exch, par, io);
			//std::cerr << GSTAR "TIME: " << dom.etime << " SW dt: " << dom.dt <<"\n";
			#if SERGHEI_DRAINAGE_MODEL
			RFF.simTime = dom.etime;
			Ddyw.drainageSimTime = dom.etime;

			Kokkos::Timer drainageExchangeTimer;
			Ddyw.DrainageExchange_solve(state, dom, Tnode); 
			dom.timers.drainageExchange += drainageExchangeTimer.seconds();
			Ddyw.snapshotExchangeFluxes(Tnode);

			Kokkos::Timer routingTimer;
			Ddyw.routing_execute(Ddyw.routingStep, Tnode, Tlink, Tconduit, Toutfall, Tstorage,
			    Tpump, TpumpCurves, TRiver, dom.timers);
    		dom.timers.drainage += routingTimer.seconds();
			#endif
			#else
			ss.swss.ComputeSWSourceSink(state, dom);
			#endif
			// run subsurface model
			#if SERGHEI_SUBSURFACE_MODEL
				#if SERGHEI_SWE_MODEL
				// If both surface and subsurface modules are on
				// Rainfall is first read by the surface module, then copy to the subsurface
				if (gdom.isRain) {Kokkos::deep_copy(gdom.rainRate, ss.swss.rainRate);}
				Kokkos::deep_copy(gw.hs, state.h);
				#endif
			// Asynchronous coupling
			// timer.reset();
			if (gdom.async)	{
				// printf("asy:gdom.dt = %f\n",gdom.dt);
				// printf("asy:dom.dt = %f\n",dom.dt);
				if (gdom.etime + gdom.dt < dom.etime)	{
					gdom.etime += gdom.dt;
					#ifdef __NVCC__
					//	PC scheme
					if (gdom.gw_scheme == 1)	{
						gwf.pca_solve<Kokkos::Cuda>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
					}
					//	Modified Picard scheme
					else {
						gwf.picard_solve<Kokkos::Cuda>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
					}
					#else
					if (gdom.gw_scheme == 1)	{
						gwf.pca_solve<Kokkos::OpenMP>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
					}
					else {
						gwf.picard_solve<Kokkos::OpenMP>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
					}
					#endif
					//std::cerr << GSTAR "TIME: " << dom.etime << " GW dt: " << gdom.dt <<"\n\n";
				}
			}
			else {
				// printf("syn:gdom.dt = %f\n",gdom.dt);
				// printf("syn:dom.dt = %f\n",dom.dt);
				gdom.etime = dom.etime;
				#ifdef __NVCC__
				if (gdom.gw_scheme == 1)	{
					gwf.pca_solve<Kokkos::Cuda>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
				}
				else {
					gwf.picard_solve<Kokkos::Cuda>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
				}
				#else
				if (gdom.gw_scheme == 1)	{
					gwf.pca_solve<Kokkos::OpenMP>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
				}
				else {
					gwf.picard_solve<Kokkos::OpenMP>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
				}
				#endif
			}
			gdom.cg_iter += A.cg_iter;
			// gdom.timers.gw += timer.seconds();
				// surface-subsurface exchange
				#if SERGHEI_SWE_MODEL
					Kokkos::deep_copy(state.qss, gw.qss);
					tint.computeGwExchange(state , dom);
				#endif
			#endif

			oldVolume+=(bint.inflowDischargeG - bint.outflowDischargeG)*dom.dt; //Boundary fluxes with the new dt
			bint.integrate(ebc.extbc,dom,0);//called here with mode==0 (adjusted volume)
			oldVolume+=bint.adjustedVolumeG; //Some mass changes can occur through the boundaries
			sint.integrate(state,dom,ss.swss); //new mass after the new time step integration
			oldVolume+= (sint.rainFluxG-sint.infFluxG)*dom.dt; //after integrate, we have to sum the rain and inf mass
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

			#if SERGHEI_LPT
			  parTrack.update(dom,state);
			#endif

			if (dom.nIter%io.nScreen==0 || fabs(dom.etime - dom.startTime - io.numOut*io.outFreq) <= TOL12) {
				if (par.masterproc) {
					std::cerr << std::fixed;
					std::cerr << GSTAR "TIME: " << dom.etime << " average dt: " << accumDt/dom.countIterDt <<"\n";
					std::cerr.precision(9);
					std::cerr << std::scientific;
					std::cerr << std::fixed;
					std::cerr.precision(12);
					#if SERGHEI_SWE_MODEL
					std::cerr << "     Surface Volume:\t" << newVolume <<"\n";
					std::cerr << "     Surface inflow: " << bint.inflowDischargeG <<"\n";
					std::cerr << "     Surface outflow Volume: " << bint.outflowDischargeG*dom.dt <<"\n";
					#if SERGHEI_SUBSURFACE_MODEL
					std::cerr << "     Exchange Volume: " << gint.Vexch_glob <<"\n";
					#endif
					#endif

					if(fabs(diffVolume)>TOL_MASS_ERROR){
						// std::cerr << YEXC "   Old Volume:\t" << oldVolume <<"\n";
						// std::cerr << YEXC "   New Volume:\t" << newVolume <<"\n";
						// std::cerr << YEXC "   Diff Volume:\t" << newVolume-oldVolume <<"\n";
						// std::cerr << YEXC "   Inflow Volume:\t" << bint.inflowDischargeG*dom.dt <<"\n";
						// std::cerr << YEXC "   Outflow Volume:\t" << bint.outflowDischargeG*dom.dt <<"\n";
						// std::cerr << YEXC "   Adjusted Volume:\t" << bint.adjustedVolumeG <<"\n";
						std::cerr << YEXC "   Rain Volume:\t" << sint.rainFluxG*dom.dt <<"\n";
						// std::cerr << YEXC "   Inf Volume:\t" << sint.infFluxG*dom.dt <<"\n";
						#if SERGHEI_DEBUG_MASS_CONS > 1
                            getchar();
                        #endif
					}
				}
				if(fabs(dom.etime - dom.startTime - io.numOut*io.outFreq) <= TOL12){
					
					#if SERGHEI_LPT==0
					  #if SERGHEI_SWE_MODEL
					  io.output(state, dom, ss.swss, par,outFolder);
					  #endif
					  #if SERGHEI_SUBSURFACE_MODEL
					  io.outputSubsurface(gw, gdom, par,outFolder);
					  #endif
					#else
					  io.output(state, dom, ss.swss, par,outFolder, parTrack);
					#endif

					if(par.masterproc) std::cerr << GIO "File " << io.numOut-1 << " written" << std::endl; //io.numOut already updated
				}
				if(par.masterproc) std::cerr << "-------------------------------------------------\n";
				dom.countIterDt=0;
				accumDt=0.0;
			}

			if (dom.etime >= io.numObs*io.obsFreq) {
				#if SERGHEI_TOOLS
				obs.update(state,par,dom);
				#endif
				io.writeTimeSeries(state,dom,par,sint,bint,ebc.extbc);
				#if SERGHEI_SUBSURFACE_MODEL
				io.writeSubsurfaceTimeSeries(gdom, gint);
				#endif
				#if SERGHEI_DRAINAGE_MODEL
				io.writeDrainageOutputTimeSeries (par, dom, sint, Ddyw, Tlink, Tnode);
				#endif
				if (par.masterproc){
					#if SERGHEI_TOOLS
					obs.write(dom);
					#endif
				}
			}

			// Unify dt for coupled simulations
				#if SERGHEI_SWE_MODEL
				tint.computeDt(state,dom,io);

				#if SERGHEI_DRAINAGE_MODEL
				printf("the swe dt and drainage dt are: %f, %f\n", dom.dt, Ddyw.routingStep);
				if (dom.dt < Ddyw.routingStep)
				{
					Ddyw.routingStep = dom.dt;
				}
				else{dom.dt = Ddyw.routingStep;}
				#endif
			#endif 
			#if SERGHEI_SWE_GW
				tint.computeDt(state,dom,io);
				#if SERGHEI_DRAINAGE_MODEL
				if (dom.dt < Ddyw.routingStep)
				{
					Ddyw.routingStep = dom.dt;
				}
				else{dom.dt = Ddyw.routingStep;}
				#endif
				#if SERGHEI_SUBSURFACE_MODEL
				if (!gdom.async)	{
					if (dom.dt < gdom.dt)	{gdom.dt = dom.dt;}
					else {dom.dt = gdom.dt;}
				}
				else if (dom.dt > gdom.dt)	{dom.dt = gdom.dt;}
				#endif
			#elif SERGHEI_SUBSURFACE_MODEL
				dom.dt = gdom.dt;
				tint.dtMatchOutput(dom,io);
				gdom.dt = dom.dt;
			#elif !SERGHEI_SWE_MODEL
        std::cout << RERROR << "Impossible configuration without SWE nor GW model" << std::endl;
				return 0;
			#endif
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
		#if SERGHEI_SUBSURFACE_MODEL
		dom.timers.gw = gdom.timers.gw;
		dom.timers.gwBC = gdom.timers.gwBC;
		dom.timers.gwlinsys = gdom.timers.gwlinsys;
		dom.timers.gwlinsol = gdom.timers.gwlinsol;
		dom.timers.gwUpdateK = gdom.timers.gwUpdateK;
		dom.timers.gwUpdateQ = gdom.timers.gwUpdateQ;
		dom.timers.gwUpdateWC = gdom.timers.gwUpdateWC;
		dom.timers.gwIntegrate = gdom.timers.gwIntegrate;
		dom.timers.gwMPI = gdom.timers.gwMPI;
		dom.timers.out += gdom.timers.out;
		dom.cg_iter = gdom.cg_iter;
		#endif
		io.writeLogFile(dom, par, outFolder);
		#if SERGHEI_LPT
		  io.writeParticleFile(dom,parTrack,par,outFolder);
		#endif
		io.closeOutputStreams();
		#if SERGHEI_TOOLS
			if(par.masterproc) obs.closeOutputStreams();
		#endif
		return 1;
	}


};
