/* -*- mode: c++; c-default-style: "linux" -*- */
// ============================================================================
// SERGHEI - Merged serghei.h
// ============================================================================
// This file merges two development branches:
//   - CODE1 (src copy/): Root water uptake, evapotranspiration, solute transport,
//     crop growth model coupling, surface-subsurface solute exchange
//   - CODE2 (src/): Granular timer system, scalar transport, sediment,
//     LPT, improved memory management, Kokkos::DefaultExecutionSpace
//
// Merge rules:
//   - Macro renames: SERGHEI_SWE_GW → SERGHEI_SWE_RE,
//                    SERGHEI_SUBSURFACE_MODEL → SERGHEI_RE_MODEL
//   - [CODE1] marks code from Code 1 (user's branch)
//   - [CODE2] marks code from Code 2 (other developer's branch)
//   - [MERGE] marks fusion/conflict resolution code
// ============================================================================

#include "stdlib.h"
#include <iostream>
#include <string>
// [CODE2] Code 2 uses io.inFolder / io.outFolder instead of separate members
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

// [CODE2] Subsurface model uses SERGHEI_RE_MODEL (renamed from SERGHEI_SUBSURFACE_MODEL)
#if SERGHEI_RE_MODEL
#include "GwDomain.h"
#include "GwFunction.h"
#include "GwInit.h"
#include "GwMPI.h"
#include "GwMatrix.h"
#include "GwState.h"
#include "GwSolver.h"
#include "GwIntegrator.h"
#if SERGHEI_WAVE_MODEL
#include "WaveBoundaryInit.h"
#include "WaveBoundary.h"
#endif
#endif

// [CODE1] Subsurface transport module
#if SERGHEI_SUBSURFACE_TRANSPORT
#include "RTStateGW.h"
#include "RTInitGW.h"
#include "RTFunctionGW.h"
#include "RTMatrix.h"
#include "RTSolverGW.h"
#include "RTIntegrator.h"
#include "GwMPI.h"
#endif

// [CODE1] Surface transport module
#if SERGHEI_SURFACE_TRANSPORT
#include "RTStateSW.h"
#include "RTInitSW.h"
#include "RTFunctionSW.h"
#endif

// [CODE1] Crop growth model (WOFOST)
#if CROP_GROWTH_MODEL
#include "./cropsrc/Wofost72.h"
#include "./cropsrc/CropState.h"
#include "./cropsrc/MeteoState.h"
#include "./cropsrc/MeteoInit.h"
#include "./cropsrc/CropInit.h"
#include "./cropsrc/WofostInit.h"
#include "./cropsrc/EvapotranspirationState.h"
#include "SourceSinkCrop.h"
#endif

class SERGHEI
{
public:
	Parallel par;
	Initializer init;
	FileIO io;
	State state;
	Domain dom;

// [CODE2] Subsurface model members (using Kokkos::DefaultExecutionSpace)
#if SERGHEI_RE_MODEL
	GwState gw;
	GwDomain gdom;
	GwInit ginit;
	GwMPI gmpi;
	GwFunction gwf;
	SubsurfaceBoundaries gbc;
	GwMatrix A;
	GwIntegrator gint;
	GwSolver<Kokkos::DefaultExecutionSpace> gsolver;
#if SERGHEI_WAVE_MODEL
	WaveBoundaryState waveBC;
	WaveBoundaryInit waveInit;
        WaveBoundaryModel waveBoundary;
#endif
#endif

// [CODE1] Subsurface transport members
#if SERGHEI_SUBSURFACE_TRANSPORT
	RTStateGW rtgw;
	RTInitgw rtinitgw;
	RTFunctionGW rtfgw;
	RTSubsurfaceBoundaries rtgbc;
	RTMatrix rtA;
	RTIntegrator rtint;
	RTSolverGW<Kokkos::DefaultExecutionSpace> rtsolvergw;
#endif

// [CODE1] Surface transport members
#if SERGHEI_SURFACE_TRANSPORT
	RTStateSW rtsw;
	RTInitSW rtinitsw;
	RTFunctionSW rtfsw;
	RTIntegratorSW rtintsw;
#endif

	SourceSink ss;
	Exchange exch;
	ExternalBoundaries ebc;

// [CODE1] Crop growth model members
#if CROP_GROWTH_MODEL
	Wofost72 wofost;
	CropState cropParam;
	MeteoState meteo;
	MeteoInit meteoInit;
	CropInit cropInit;
	EvapotranspirationState ets;
	WofostInit wofostInit;
	double time_since_last_crop_step = 0.0;
	const double CROP_STEP_SIZE = 86400.0; // [CODE1] 1 day in seconds
	bool wofost_initialized = false;
#endif

private:
	Parser parser;
	TimeIntegrator tint;
	surfaceIntegrator sint;
	boundaryIntegrator bint;

#if SERGHEI_TOOLS
	Observations obs;
#endif

// [CODE2] Particle tracking (renamed macro: SERGHEI_PARTICLE_TRACKING → SERGHEI_LPT)
#if SERGHEI_LPT
	ParticleTracker parTrack;
#endif

	double oldVolume, newVolume, diffVolume;
// [CODE2] Suspended sediment tracking
#if SERGHEI_SUSPENDED_SEDIMENT
	real oldSolidVolume, newSolidVolume, diffSolidVolume;
#endif
	double accumDt = 0.0;

	// [CODE2] Kokkos objects
	Kokkos::Timer timer;
	Kokkos::Timer timer_particles; // [CODE2] Separate timer for particles
	Kokkos::InitializationSettings kokkosSettings;

	////////////// METHODS ///////////////
public:
	// [MERGE] start() method: primarily Code 2 structure with Code 1 additions
	int start(int argc, char **argv)
	{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif

		init.initializeMPI(&argc, &argv, par);

// [CODE2] Use KOKKOS_ENABLE_CUDA instead of __NVCC__
#ifdef KOKKOS_ENABLE_CUDA
		kokkosSettings.set_device_id(par.myrank % par.nthreads);
#else
		if (par.nthreads != 0)
			kokkosSettings.set_num_threads(par.nthreads);
#endif

#if SERGHEI_DEBUG_KOKKOS_SETUP
		printKokkosInitArguments(par);
#if KOKKOS_ENABLE_CUDA
		printKokkosCuda(args, par);
#endif
#endif

#if SERGHEI_DEBUG_WORKFLOW
		std::cerr << GGD "Initialising Kokkos - rank " << par.myrank << std::endl;
#endif

		Kokkos::initialize(kokkosSettings);

#if SERGHEI_DEBUG_WORKFLOW
		std::cerr << GGD "Program instantiated, creating objects - rank " << par.myrank << std::endl;
#endif

		// [CODE2] Initialize the surface model
		// [MERGE] Code 2 uses io.inFolder/io.outFolder; Code 1 passed inFolder/outFolder separately
		if (!init.initialize(state, ss.swss, ebc, dom, par, tint, sint, bint, parser, exch, io))
		{
			if (par.masterproc)
				std::cerr << RERROR "Unable to start the simulation" << std::endl;
			return 0;
		};

#if SERGHEI_WAVE_MODEL
#if !SERGHEI_RE_MODEL
		if (dom.isWave)
		{
			if (par.masterproc)
				std::cerr << RERROR << "Wave module requires SERGHEI_RE_MODEL to be enabled." << std::endl;
			return 0;
		}
#endif
#endif

// [CODE2] LPT initialization
#if SERGHEI_LPT
		if (!parser.readParticles(io.inFolder, par, &parTrack))
		{
			if (par.masterproc)
				std::cerr << RERROR "Unable to start the simulation because of LPT initialization files" << std::endl;
			return 0;
		}
		parTrack.initialiseParticles(dom, state);

		if (io.outFormat == OUT_VTK)
		{
#if SERGHEI_PARTICLE_NO_OUTPUT
#else
			io.outputIniParticle(dom, parTrack, io.outFolder);
#endif
		}
		else
		{
			io.outputInitParticlesNETCDF(parTrack, dom, par, io.outFolder);
		}
		if (par.masterproc)
			std::cerr << GOK "LPT module has been initialized! Number of particles: " << parTrack.N_par << std::endl
					  << std::endl;
#endif

// [CODE2] Initialize subsurface model if activated (using SERGHEI_RE_MODEL)
#if SERGHEI_RE_MODEL
		if (!ginit.initialize_gw(gw, gdom, state, dom, gbc, gmpi, gint, par, io, ss, io.inFolder, io.outFolder))
		{
			std::cerr << RERROR "Unable to initialize the subsurface domain" << "\n";
			return 0;
		};
#if SERGHEI_WAVE_MODEL
		if (dom.isWave)
		{
			if (!waveInit.initialize(waveBC, gdom, gbc.gwbc, par, io))
			{
				std::cerr << RERROR << "Unable to initialize wave boundary preprocessing" << "\n";
				return 0;
			}
		}
#endif
		A.init(gdom);
		gsolver.init(A, gdom);
		if (par.masterproc)
		{
			std::cerr << GOK "Subsurface Solver has been initialized! " << std::endl;
		}
#endif

// [CODE1] Initialize subsurface transport module
#if SERGHEI_SUBSURFACE_TRANSPORT
		DEBUG_SUBSURFACE_TRANSPORT_PRINT("Initializing subsurface transport module\n");
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtinitgw.initialize_rt()" << std::endl;
#endif
		if (!rtinitgw.initialize_rt(rtgw, gw, gdom, rtgbc, gmpi, par, io, ss, io.inFolder, io.outFolder))
		{
			std::cerr << RERROR "Unable to initialize the transport module" << "\n";
			return 0;
		}
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtinitgw.initialize_rt() completed OK" << std::endl;
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtA.init(gdom)" << std::endl;
#endif
		rtA.init(gdom);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtA.init(gdom) completed OK" << std::endl;
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtsolvergw.init(rtA, gdom)" << std::endl;
#endif
		rtsolvergw.init(rtA, gdom);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtsolvergw.init(rtA, gdom) completed OK" << std::endl;
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtint.initialize()" << std::endl;
#endif
		rtint.initialize(rtgw, gdom, rtgbc.rtgwbc, ss);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtint.initialize() completed OK" << std::endl;
#endif
		DEBUG_SUBSURFACE_TRANSPORT_PRINT("Subsurface transport module initialized\n");
#endif

// [CODE1] Initialize surface transport module
#if SERGHEI_SURFACE_TRANSPORT
		DEBUG_SURFACE_TRANSPORT_PRINT("Initializing surface transport module\n");
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtinitsw.initialize_rt()" << std::endl;
#endif
		if (!rtinitsw.initialize_rt(rtsw, state, dom, ss, ebc, io, par, io.inFolder, io.outFolder))
		{
			std::cerr << RERROR "Unable to initialize the SW transport module" << "\n";
			return 0;
		}
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtinitsw.initialize_rt() completed OK" << std::endl;
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtintsw.initialize()" << std::endl;
#endif
		rtintsw.initialize(rtsw, state, dom, ss.swss, ebc);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtintsw.initialize() completed OK" << std::endl;
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtfsw.initializeExchange()" << std::endl;
#endif
		rtfsw.initializeExchange(rtsw, dom);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtfsw.initializeExchange() completed OK" << std::endl;
#endif
		if (par.masterproc)
		{
			std::cerr << GOK "Surface transport MPI exchange has been initialized! " << std::endl;
		}
		DEBUG_SURFACE_TRANSPORT_PRINT("Surface transport module initialized\n");
#endif

// [CODE1] Initialize WOFOST crop growth model
#if CROP_GROWTH_MODEL
		DEBUG_CROP_PRINT("Initializing WOFOST crop model\n");
		if (par.masterproc)
			std::cout << GOK "Initializing WOFOST Crop Model..." << std::endl;

		if (!wofostInit.initialize_wofost(wofost, cropParam, meteo,
										  cropInit, meteoInit,
										  gw, gdom, dom,
										  par, io, io.inFolder, io.outFolder))
		{
			std::cerr << RERROR "Failed to initialize WOFOST" << std::endl;
			return 0;
		}
		wofost_initialized = true;

		for (int k = 0; k < ss.gwss.size(); k++)
		{
			if (ss.gwss[k].use_realtime_data)
			{
				ss.gwss[k].allocateWofostCoupling(dom.nCell);
			}
		}
		if (par.masterproc)
		{
			std::cout << GOK "WOFOST initialized with " << dom.nCell << " crop cells" << std::endl;
		}
		DEBUG_CROP_PRINT("WOFOST crop model initialized\n");
#endif

// [CODE2] SERGHEI_TOOLS observations
#if SERGHEI_TOOLS
		if (!obs.readInputFiles(io.inFolder, par))
			return 0;
		if (!obs.configure(dom, io.outFolder))
			return 0;
		obs.update(state, par, dom);
		if (par.masterproc)
		{
			obs.writeLinesSamplingCoordinates(io.outFolder);
			obs.writeGauges(dom.etime);
			obs.writeLines(dom.etime);
		}
#endif

#if SERGHEI_DEBUG_WORKFLOW
		for (int k = 0; k < ebc.extbc.size(); k++)
		{
			std::cout << GGD << GRAY << __FILE__ << ":" << __LINE__ << RESET << "\tExtBC[" << k << "]: " << ebc.extbc[k].bcvals(0) << ", " << ebc.extbc[k].bcvals(1) << ", " << ebc.extbc[k].bcvals(2) << std::endl;
		}
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "Initialisation finished, starting to run main loop" << std::endl;
		for (int i = 0; i < ebc.extbc.size(); i++)
		{
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "ncellsBC for segment " << i << ": " << ebc.extbc[i].ncellsBC << "\n";
		}
#if SERGHEI_DEBUG_BOUNDARY
		bint.integrate(ebc.extbc, dom, 1);
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "ncellsBC (integrated) " << bint.ncellsBC << "\n";
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow discharge (integrated) " << bint.outflowDischarge << "\n";
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow accumulated (integrated) " << bint.outflowAccumulated << std::endl;
#endif
#endif

		// [CODE2] Integrator at the beginning of the simulation
		sint.integrate(state, dom, ss.swss);
		// [CODE2] Write initial time series data
		io.writeTimeSeriesIni(state, dom, par, ss.swss, sint, bint, ebc.extbc, io.outFolder);

// [CODE1] Surface transport initial time series
#if SERGHEI_SURFACE_TRANSPORT
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtintsw.integrate() (initial)" << std::endl;
#endif
		rtintsw.integrate(rtsw, state, dom, ss.swss, ebc);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtintsw.integrate() (initial) completed OK" << std::endl;
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering io.writeRTSurfaceTimeSeriesIni()" << std::endl;
#endif
		io.writeRTSurfaceTimeSeriesIni(dom, rtsw, rtintsw, par, io.outFolder);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> io.writeRTSurfaceTimeSeriesIni() completed OK" << std::endl;
#endif
#endif

// [CODE2/MERGE] Subsurface initial output (Code 2 signature with gdom first)
#if SERGHEI_RE_MODEL
		io.writeSubTimeSeriesIni(gdom, gint, par, io.outFolder);
		io.outputSubsurface(gw, gdom, par, io.outFolder);
#endif

// [CODE1] Subsurface transport initial output
#if SERGHEI_SUBSURFACE_TRANSPORT
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtint.integrate() (initial)" << std::endl;
#endif
		rtint.integrate(gw, rtgw, gdom, rtgbc.rtgwbc);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtint.integrate() completed OK" << std::endl;
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering io.writeRTSubsurfaceTimeSeriesIni()" << std::endl;
#endif
		io.writeRTSubsurfaceTimeSeriesIni(dom, rtgw, rtint, par, io.outFolder);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> io.writeRTSubsurfaceTimeSeriesIni() completed OK" << std::endl;
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering io.outputTransport() (initial)" << std::endl;
#endif
		io.outputTransport(rtgw, gdom, par, io.outFolder);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> io.outputTransport() completed OK" << std::endl;
#endif
#endif

// [CODE1] Crop growth initial output
#if CROP_GROWTH_MODEL
		io.writeRootZoneTimeSeriesIni(dom, par, io.outFolder, ss.gwss);
		io.outputCrop(wofost, dom, par, io.outFolder);
#endif

		// [CODE2] Capture initialisation time using granular timer
		dom.timers.swe.init.total = timer.seconds();

		if (par.masterproc)
			std::cout << GOK << "Initialisation complete. Initialisation time: " << dom.timers.swe.init.total << " [s]" << std::endl;
		return 1;
	}

	// ====================================================================
	// compute() - Main simulation loop
	// [MERGE] Based on Code 2 structure with Code 1's transport, crop, ET features
	// ====================================================================
	int compute()
	{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		if (par.masterproc)
		{
			std::cout << std::endl
					  << GOK "SIMULATION STARTS" << std::endl;
			std::cout << BDASH << "Start time: " << dom.startTime << std::endl;
			std::cout << BDASH << "End time: " << dom.endTime << std::endl;
		}

// [CODE2] Initialize timestep
#if SERGHEI_RE_MODEL
		gdom.dt = gdom.dt_init;
		gdom.dtOld = gdom.dt_init;
		dom.dt = gdom.dt;
		gdom.cg_iter = 0;
#else
		tint.computeDt(state, dom, io);
#endif

// [CODE1] Initialize transport timesteps
#if SERGHEI_SUBSURFACE_TRANSPORT
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Initializing transport timesteps: rtgw.dt=" << gdom.dt_init << std::endl;
#endif
		rtgw.dt = gdom.dt_init;
		rtgw.dtOld = gdom.dt_init;
// [MERGE] transport dt sync: only set if RE_MODEL is not active
// (when RE_MODEL is active, dt is set in the RE_MODEL block above)
#if !SERGHEI_RE_MODEL
		dom.dt = rtgw.dt;
#endif
#endif

		// Main Time Loop
		while (dom.etime < dom.endTime)
		{
			// printf("Model Time = %f\n", dom.etime);
			// previous mass
			oldVolume = sint.surfaceVolumeG;
#if SERGHEI_SUSPENDED_SEDIMENT // [CODE2]
			oldSolidVolume = sint.surfaceSolidVolumeG;
#endif

			// [CODE2] Boundary integration (mode=1: boundary flows)
			bint.integrate(ebc.extbc, dom, 1);

// [CODE2] Run surface model
#if SERGHEI_SWE_MODEL
			tint.stepForward(state, ss.swss, ebc.extbc, dom, exch, par, io);
#else
			ss.swss.ComputeSWSourceSink(state, dom);
#endif

// ============================================================
// [CODE1] Subsurface model block with transport, crop, ET coupling
// [CODE2] Uses Kokkos::DefaultExecutionSpace and SERGHEI_RE_MODEL
// ============================================================
#if SERGHEI_RE_MODEL
#if SERGHEI_SWE_MODEL
			// [CODE2] If both surface and subsurface modules are on
			// Rainfall is first read by the surface module, then copy to the subsurface
			{
				Kokkos::Timer timer_coupling;
				timer_coupling.reset();
				if (gdom.isRain)
				{
					Kokkos::deep_copy(gdom.rainRate, ss.swss.rainRate);
				}
				Kokkos::deep_copy(gw.hs, state.h);
// [CODE1] Evapotranspiration: copy surface evaporation data to groundwater
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
				Kokkos::deep_copy(gw.surfaceWaterEvapActual, state.surfaceEvaporation);
#endif
				dom.timers.coupling.total += timer_coupling.seconds();
			}
#endif

// [CODE1] WOFOST crop growth model step
#if CROP_GROWTH_MODEL
			time_since_last_crop_step += dom.dt;
			if (time_since_last_crop_step >= CROP_STEP_SIZE && wofost_initialized)
			{
				int days_passed = (int)(dom.etime / 86400.0);
				int current_doy = ((meteo.start_doy + days_passed - 1) % 365) + 1;

				DEBUG_CROP_PRINT("Running WOFOST step at t=%f, DOY=%d\n", dom.etime, current_doy);
				if (par.masterproc)
				{
					std::cout << GOK "Running WOFOST Step at t=" << dom.etime
							  << " (Meteo DOY: " << current_doy << ")" << std::endl;
				}

				Kokkos::Timer timer_crop_total;
				timer_crop_total.reset();

				// [CODE1] Compute root zone moisture using WOFOST root depth
				if (par.masterproc)
				{
					std::cout << GOK "About to compute root zone moisture with WOFOST RD" << std::endl;
				}
				{
					Kokkos::Timer timer_crop;
					timer_crop.reset();
					SourceSinkCrop::compute_root_zone_moisture(gw, gdom, wofost.s.RD, wofost.ets.RZSM);
					// [MERGE] Accumulate crop timer (Code 2 timers structure)
					dom.timers.crop.rzm += timer_crop.seconds();
				}

				if (par.masterproc)
				{
					std::cout << GOK "Root zone moisture computed using WOFOST RD" << std::endl;
					std::cout << GOK "Sample wofost.ets.RZSM[0]=" << wofost.ets.RZSM(0) << std::endl;
				}

				// [CODE1] WOFOST growth calculation
				if (par.masterproc)
				{
					std::cout << GOK "About to run WOFOST calc_rates_and_integrate" << std::endl;
				}
				{
					Kokkos::Timer timer_crop;
					timer_crop.reset();
					wofost.calc_rates_and_integrate(cropParam, meteo, dom.etime,
													wofost.ets.RZSM, current_doy, 1.0);
					// [MERGE] Accumulate crop timer
					dom.timers.crop.calc += timer_crop.seconds();
				}

				if (par.masterproc)
				{
					std::cout << GOK "WOFOST calc_rates_and_integrate completed" << std::endl;
				}

				// [CODE1] Synchronize WOFOST data (LAI, RD) back to source/sink
				int nCellSwSurface = gdom.nx * gdom.ny;
				{
					Kokkos::Timer timer_crop;
					timer_crop.reset();
					for (int k = 0; k < ss.gwss.size(); k++)
					{
						if (ss.gwss[k].use_realtime_data)
						{
							auto src_lai = wofost.s.LAI;
							auto src_rd = wofost.s.RD;
							auto dst_lai = ss.gwss[k].wofost_lai;
							auto dst_rd = ss.gwss[k].wofost_rd;
							Kokkos::parallel_for("copy_wofost_data", nCellSwSurface, KOKKOS_LAMBDA(int i) {
								dst_lai(i) = src_lai(i);
								dst_rd(i)  = src_rd(i); });
							Kokkos::fence();
						}
					}
					// [MERGE] Accumulate crop timer
					dom.timers.crop.sync += timer_crop.seconds();
				}

				dom.timers.crop.total += timer_crop_total.seconds();

				time_since_last_crop_step = 0.0;

				if (par.masterproc)
				{
					double lai_val = wofost.s.LAI(0);
					double rd_val = wofost.s.RD(0);
					std::cout << GOK "WOFOST step completed. LAI[0]=" << lai_val
							  << ", RD[0]=" << rd_val << ", TAGP[0]=" << wofost.TAGP(0) << std::endl;
				}
				DEBUG_CROP_PRINT("WOFOST step completed. LAI[0]=%f, RD[0]=%f\n", wofost.s.LAI(0), wofost.s.RD(0));

// [CODE1] Output crop data
#if CROP_GROWTH_MODEL
				io.outputCrop(wofost, dom, par, io.outFolder);
#endif
			}
#endif // CROP_GROWTH_MODEL

// [CODE1] Copy surface transport concentration to subsurface transport
#if SERGHEI_SUBSURFACE_TRANSPORT
#if SERGHEI_SURFACE_TRANSPORT
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
			if (par.masterproc)
				std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Copying surface concentration to subsurface (csw, RainCon) at t=" << dom.etime << std::endl;
#endif
			Kokkos::deep_copy(rtgw.csw, rtsw.csw4gw);
			Kokkos::deep_copy(rtgw.RainCon_arr, rtsw.RainCon_arr);
#endif
#endif

			// [CODE2] Asynchronous coupling for subsurface
			if (gdom.async)
			{
				if (gdom.etime + gdom.dt < dom.etime)
				{
					gdom.etime += gdom.dt;
					// [CODE2] PC scheme / Modified Picard with DefaultExecutionSpace
					// [CODE1] Pass wofost to solver when CROP_GROWTH_MODEL is active
					if (gdom.gw_scheme == 1)
					{
#if CROP_GROWTH_MODEL
#if SERGHEI_WAVE_MODEL
						gwf.pca_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, wofost, &waveBoundary, &waveBC, &ss.swss, &dom, &io);
#else
						gwf.pca_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, wofost);
#endif
#else
#if SERGHEI_WAVE_MODEL
						gwf.pca_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, &waveBoundary, &waveBC, &ss.swss, &dom, &io);
#else
						gwf.pca_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
#endif
#endif
					}
					else
					{
#if CROP_GROWTH_MODEL
#if SERGHEI_WAVE_MODEL
						gwf.picard_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, wofost, &waveBoundary, &waveBC, &ss.swss, &dom, &io);
#else
						gwf.picard_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, wofost);
#endif
#else
#if SERGHEI_WAVE_MODEL
						gwf.picard_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, &waveBoundary, &waveBC, &ss.swss, &dom, &io);
#else
						gwf.picard_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
#endif
#endif
					}
// [CODE1] Subsurface transport solve (async)
#if SERGHEI_SUBSURFACE_TRANSPORT
					rtgw.dt = gdom.dt;
					DEBUG_SUBSURFACE_TRANSPORT_PRINT("Calling rt_solve for subsurface transport (async)\n");
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtfgw.rt_solve() [ASYNC] at gdom.etime=" << gdom.etime << " dt=" << gdom.dt << std::endl;
#endif
					rtfgw.rt_solve<Kokkos::DefaultExecutionSpace>(rtgw, rtA, gw, gdom, rtgbc.rtgwbc, gmpi, par, rtint, rtsolvergw, ss);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtfgw.rt_solve() [ASYNC] completed OK" << std::endl;
#endif
#endif
				}
			}
			else
			{
				gdom.etime = dom.etime;
				if (gdom.gw_scheme == 1)
				{
#if CROP_GROWTH_MODEL
#if SERGHEI_WAVE_MODEL
					gwf.pca_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, wofost, &waveBoundary, &waveBC, &ss.swss, &dom, &io);
#else
					gwf.pca_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, wofost);
#endif
#else
#if SERGHEI_WAVE_MODEL
					gwf.pca_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, &waveBoundary, &waveBC, &ss.swss, &dom, &io);
#else
					gwf.pca_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
#endif
#endif
				}
				else
				{
#if CROP_GROWTH_MODEL
#if SERGHEI_WAVE_MODEL
					gwf.picard_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, wofost, &waveBoundary, &waveBC, &ss.swss, &dom, &io);
#else
					gwf.picard_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, wofost);
#endif
#else
#if SERGHEI_WAVE_MODEL
					gwf.picard_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par, &waveBoundary, &waveBC, &ss.swss, &dom, &io);
#else
					gwf.picard_solve<Kokkos::DefaultExecutionSpace>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
#endif
#endif
				}
// [CODE1] Subsurface transport solve (sync)
#if SERGHEI_SUBSURFACE_TRANSPORT
				rtgw.dt = gdom.dt;
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
				if (par.masterproc)
					std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtfgw.rt_solve() [SYNC] at dom.etime=" << dom.etime << " dt=" << gdom.dt << std::endl;
#endif
				rtfgw.rt_solve<Kokkos::DefaultExecutionSpace>(rtgw, rtA, gw, gdom, rtgbc.rtgwbc, gmpi, par, rtint, rtsolvergw, ss);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
				if (par.masterproc)
					std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtfgw.rt_solve() [SYNC] completed OK" << std::endl;
#endif
#endif
			}
			gdom.cg_iter += A.cg_iter;

// [CODE2] Surface-subsurface exchange
#if SERGHEI_SWE_MODEL
			{
				Kokkos::Timer timer_coupling;
				timer_coupling.reset();
				Kokkos::deep_copy(state.qss, gw.qss);
				tint.computeGwExchange(state, dom);
				dom.timers.coupling.total += timer_coupling.seconds();
			}
#endif
#endif // SERGHEI_RE_MODEL

// [CODE1] Surface transport step (all RTSW operations)
#if SERGHEI_SURFACE_TRANSPORT
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
			if (par.masterproc)
				std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> === Surface transport step at t=" << dom.etime << " dt=" << dom.dt << " ===" << std::endl;
#endif
			Kokkos::Timer timer_rtsw_total;
			timer_rtsw_total.reset();

			rtsw.dt = dom.dt;
// [CODE1] Copy subsurface solute exchange to surface transport
#if SERGHEI_SUBSURFACE_TRANSPORT
			DEBUG_SURFACE_TRANSPORT_PRINT("Copying ConQss from subsurface to surface transport\n");
			Kokkos::deep_copy(rtsw.ConQss, rtgw.ConQss);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
			if (par.masterproc)
				std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Copied ConQss from subsurface to surface" << std::endl;
#endif
#endif
			// [CODE1] Compute GW exchange for surface transport (solute exchange)
			{
				Kokkos::Timer timer_rtswexch;
				DEBUG_SURFACE_TRANSPORT_PRINT("Computing GW exchange for surface transport\n");
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
				if (par.masterproc)
					std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtfsw.computeGwExchange_For_RT()" << std::endl;
#endif
				rtfsw.computeGwExchange_For_RT(state, rtsw, dom);
				Kokkos::fence();
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
				if (par.masterproc)
					std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtfsw.computeGwExchange_For_RT() completed OK" << std::endl;
#endif
				dom.timers.rtsw.gwExch += timer_rtswexch.seconds();
			}

			DEBUG_SURFACE_TRANSPORT_PRINT("Running surface transport solver\n");
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
			if (par.masterproc)
				std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtfsw.rt_solve_sw()" << std::endl;
#endif
			rtfsw.rt_solve_sw(rtsw, state, dom, ebc, ss.swss, rtintsw, par);
			Kokkos::fence();
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
			if (par.masterproc)
				std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtfsw.rt_solve_sw() completed OK" << std::endl;
#endif

			dom.timers.rtsw.total += timer_rtsw_total.seconds();
#endif

			// [CODE2] Mass balance bookkeeping
			oldVolume += (bint.inflowDischargeG - bint.outflowDischargeG) * dom.dt;
#if SERGHEI_SUSPENDED_SEDIMENT
			oldSolidVolume += (bint.inflowSolidDischargeG - bint.outflowSolidDischargeG) * dom.dt;
#endif

			bint.integrate(ebc.extbc, dom, 0); // [CODE2] mode=0: adjusted volume
			oldVolume += bint.adjustedVolumeG;
			sint.integrate(state, dom, ss.swss);
			oldVolume += (sint.rainFluxG - sint.infFluxG) * dom.dt;
#if SERGHEI_SUSPENDED_SEDIMENT // [CODE2]
			oldVolume += sint.BedExchangeVolG;
			oldSolidVolume += sint.BedExchangeSolidG;
#endif

			// new mass
			newVolume = sint.surfaceVolumeG;
			if (fabs(oldVolume) > TOL12)
			{
				diffVolume = (newVolume - oldVolume) / oldVolume * 100.;
			}
			else
			{
				diffVolume = 0.0;
			}
#if SERGHEI_SUSPENDED_SEDIMENT // [CODE2]
			newSolidVolume = sint.surfaceSolidVolumeG;
			if (fabs(oldSolidVolume) > TOL12)
			{
				diffSolidVolume = (newSolidVolume - oldSolidVolume) / oldSolidVolume * 100.;
			}
			else
			{
				diffSolidVolume = 0.0;
			}
#endif

			// [CODE2] Advance time
			dom.etime += dom.dt;
			dom.nTimeSteps++;
			dom.countIterDt++;
			accumDt += dom.dt;

// [CODE2] Particle tracking update
#if SERGHEI_LPT
			parTrack.update(dom, state);
#endif

			// ============================================================
			// [CODE2] Screen output and file output section
			// ============================================================
			if (dom.nTimeSteps % io.nScreen == 0 || fabs(dom.etime - dom.startTime - io.numOut * io.outFreq) <= TOL12)
			{
				if (par.masterproc)
				{
					std::cerr << std::fixed;
					std::cerr << GSTAR "TIME: " << dom.etime << " av_dt: " << accumDt / dom.countIterDt << " iter: " << dom.nTimeSteps << "\n";
					std::cerr.precision(9);
					std::cerr << std::scientific;
					std::cerr << "     Inner Volume: " << newVolume << "\n";
					std::cerr << "     Diff Volume: " << diffVolume << "\n";
					std::cerr << std::fixed;
					std::cerr.precision(12);
#if SERGHEI_SWE_MODEL
					std::cerr << "     Surface Volume:\t" << newVolume << "\n";
					std::cerr << "     Surface inflow: " << bint.inflowDischargeG << "\n";
					std::cerr << "     Surface outflow Volume: " << bint.outflowDischargeG * dom.dt << "\n";
#if SERGHEI_RE_MODEL
					std::cerr << "     Exchange Volume: " << gint.Vexch_glob << "\n";
#endif
#endif
#if SERGHEI_RE_MODEL
					std::cerr << "     Subsurface Volume: " << gint.Vtot_glob << "\n";
#endif

#if SERGHEI_SURFACE_TRANSPORT
					std::cerr.precision(9);
					std::cerr << std::scientific;
					for (int iSpec = 0; iSpec < rtsw.n_mass; ++iSpec)
					{
						std::cerr << "     SW Solute[" << iSpec << "] Total Mass [mg]: " << rtintsw.MassTot_spec_glob(iSpec) << "\n";
					}
#endif
#if SERGHEI_SUBSURFACE_TRANSPORT
					std::cerr.precision(9);
					std::cerr << std::scientific;
					for (int iSpec = 0; iSpec < rtgw.n_mass; ++iSpec)
					{
						std::cerr << "     GW Solute[" << iSpec << "] Total Mass [mg]: " << rtint.LiquidMassTot_spec_glob(iSpec) << "\n";
					}
					std::cerr << std::fixed;
					std::cerr.precision(12);
#endif

					if (fabs(diffVolume) > TOL_MASS_ERROR)
					{
						std::cerr << YEXC "   Rain Volume:\t" << sint.rainFluxG * dom.dt << "\n";
#if SERGHEI_DEBUG_MASS_CONS > 1
						getchar();
#endif
					}

#if SERGHEI_SUSPENDED_SEDIMENT // [CODE2]
					std::cerr.precision(9);
					std::cerr << std::scientific;
					std::cerr << "     Inner Solid Volume: " << newSolidVolume << "\n";
					std::cerr << "     Diff Solid Volume: " << diffSolidVolume << "\n";
					std::cerr << std::fixed;
					std::cerr.precision(12);
					std::cerr << "     Inflow Solid Discharge: " << bint.inflowSolidDischargeG << "\n";
					std::cerr << "     Outflow Solid Discharge: " << bint.outflowSolidDischargeG << "\n";
					if (fabs(diffSolidVolume) > TOL8)
					{
						std::cerr << YEXC "   Old Solid Volume:\t" << oldSolidVolume << "\n";
						std::cerr << YEXC "   New Solid Volume:\t" << newSolidVolume << "\n";
						std::cerr << YEXC "   Diff Solid Volume:\t" << newSolidVolume - oldSolidVolume << "\n";
						std::cerr << YEXC "   Inflow Solid Volume:\t" << bint.inflowSolidDischargeG * dom.dt << "\n";
						std::cerr << YEXC "   Outflow Solid Volume:\t" << bint.outflowSolidDischargeG * dom.dt << "\n";
#if SERGHEI_DEBUG_MASS_CONS > 1
						getchar();
#endif
					}
#endif
				}

				if (fabs(dom.etime - dom.startTime - io.numOut * io.outFreq) <= TOL12)
				{
					Kokkos::Timer timer_io;
					timer_io.reset();
// [CODE2] File output with LPT support
#if SERGHEI_LPT == 0
#if SERGHEI_SWE_MODEL
					io.output(state, dom, ss.swss, par);
#endif
#if SERGHEI_RE_MODEL
					io.outputSubsurface(gw, gdom, par, io.outFolder);
#endif
#else
					io.output(state, dom, ss.swss, par, io.outFolder, parTrack);
#endif

// [CODE1] Surface transport output
#if SERGHEI_SURFACE_TRANSPORT
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering io.outputTransportSW() at t=" << dom.etime << std::endl;
#endif
					io.outputTransportSW(rtsw, state, dom, par, io.outFolder);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> io.outputTransportSW() completed OK" << std::endl;
#endif
#endif

// [CODE1] Subsurface transport output
#if SERGHEI_SUBSURFACE_TRANSPORT
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering io.outputTransport() at t=" << dom.etime << std::endl;
#endif
					io.outputTransport(rtgw, gdom, par, io.outFolder);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> io.outputTransport() completed OK" << std::endl;
#endif
#endif

					if (par.masterproc)
						std::cerr << GIO "File " << io.numOut - 1 << " written" << std::endl;
					dom.timers.io_loop.output += timer_io.seconds();
				}
				if (par.masterproc)
					std::cerr << "-------------------------------------------------\n";
				dom.countIterDt = 0;
				accumDt = 0.0;
			}

// [CODE2] Particle tracking (kept for backward compatibility)
#if SERGHEI_PARTICLE_TRACKING
			parTrack.update(dom, state);
#endif

			// ============================================================
			// [MERGE] Time series and observation output section
			// ============================================================
			if (dom.etime >= io.numObs * io.obsFreq)
			{
#if SERGHEI_TOOLS // [CODE2]
				obs.update(state, par, dom);
#endif

// [CODE1] Surface transport integration for time series
#if SERGHEI_SURFACE_TRANSPORT
				{
					Kokkos::Timer timer_rtsw_int;
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtintsw.integrate() at t=" << dom.etime << std::endl;
#endif
					rtintsw.integrate(rtsw, state, dom, ss.swss, ebc);
					Kokkos::fence();
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtintsw.integrate() completed OK" << std::endl;
#endif
					dom.timers.rtsw.integrate += timer_rtsw_int.seconds();
				}
#endif

				// [MERGE] writeTimeSeries must be called by all MPI ranks (internal MPI comm)
				io.writeTimeSeries(state, dom, par, sint, bint, ebc.extbc);

				if (par.masterproc)
				{
// [CODE1] Surface transport time series
#if SERGHEI_SURFACE_TRANSPORT
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering io.writeRTSurfaceTimeSeries() at t=" << dom.etime << std::endl;
#endif
					io.writeRTSurfaceTimeSeries(dom, rtintsw);
#if DEBUG_SERGHEI_SURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[SW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> io.writeRTSurfaceTimeSeries() completed OK" << std::endl;
#endif
#endif
// [CODE2] Subsurface time series (Code 2 signature)
#if SERGHEI_RE_MODEL
					io.writeSubsurfaceTimeSeries(gdom, gint);
#endif

// [CODE1] Subsurface transport time series
#if SERGHEI_SUBSURFACE_TRANSPORT

#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering rtint.integrate() at t=" << dom.etime << std::endl;
#endif

					Kokkos::Timer timer_rtgw_integrator;
					timer_rtgw_integrator.reset();
					rtint.integrate(gw, rtgw, gdom, rtgbc.rtgwbc);
					Kokkos::fence();
					gdom.timers.rtgw.integrate += timer_rtgw_integrator.seconds();

#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> rtint.integrate() completed OK" << std::endl;
#endif

#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Entering io.writeRTSubsurfaceTimeSeries() at t=" << dom.etime << std::endl;
#endif
					io.writeRTSubsurfaceTimeSeries(dom, rtint);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
					if (par.masterproc)
						std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> io.writeRTSubsurfaceTimeSeries() completed OK" << std::endl;
#endif
#endif
// [CODE1] Crop root zone water content time series
#if CROP_GROWTH_MODEL
					io.writeRootZoneWaterContentTimeSeries(dom, ss.gwss);
#endif
				}

				if (par.masterproc)
				{
#if SERGHEI_TOOLS // [CODE2]
					obs.write(dom);
#endif
				}
			}

			// ============================================================
			// [MERGE] Timestep synchronization section
			// SERGHEI_SWE_RE replaces SERGHEI_SWE_GW per merge instructions
			// ============================================================
#if SERGHEI_SWE_RE
			tint.computeDt(state, dom, io);
#if SERGHEI_RE_MODEL
			if (!gdom.async)
			{
				if (dom.dt < gdom.dt)
				{
					gdom.dt = dom.dt;
				}
				else
				{
					dom.dt = gdom.dt;
				}
			}
			else if (dom.dt > gdom.dt)
			{
				dom.dt = gdom.dt;
			}
// [CODE1] Sync transport timestep
#if SERGHEI_SUBSURFACE_TRANSPORT
			rtgw.dt = gdom.dt;
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
			if (par.masterproc)
				std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Sync transport dt: rtgw.dt=" << rtgw.dt << " gdom.dt=" << gdom.dt << " dom.dt=" << dom.dt << std::endl;
#endif
#endif
#endif
#elif SERGHEI_RE_MODEL
			dom.dt = gdom.dt;
			tint.dtMatchOutput(dom, io);
			if (dom.dt <= 0.0)
			{
				if (par.masterproc)
				{
					std::cerr << RERROR << "Non-positive dt detected after dtMatchOutput. "
					          << "etime=" << dom.etime << ", dt=" << dom.dt
					          << ", numOut=" << io.numOut << ", outFreq=" << io.outFreq
					          << ". This usually indicates inconsistent output counters."
					          << std::endl;
				}
				return 0;
			}
			gdom.dt = dom.dt;
#elif !SERGHEI_SWE_MODEL
			std::cout << RERROR << "Impossible configuration without SWE nor GW model" << std::endl;
			return 0;
#endif
		} // end of time loop
		return 1;
	}

	// ====================================================================
	// finalise() - Simulation cleanup and timer reporting
	// [MERGE] Primarily Code 2 structure with Code 1's transport/crop timers
	// ====================================================================
	int finalise()
	{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif

// [CODE2] LPT final file output
#if SERGHEI_LPT
		io.writeParticleFile(dom, parTrack, par, io.outFolder);
#endif

		dom.timers.total = timer.seconds();

// [CODE2] Gather RE timers from gdom
#if SERGHEI_RE_MODEL
		dom.timers.re = gdom.timers.re;
		dom.cg_iter = gdom.cg_iter;
#endif

// [CODE1] Gather subsurface transport timers (total independently timed in rt_solve)
#if SERGHEI_SUBSURFACE_TRANSPORT
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
		if (par.masterproc)
			std::cerr << GGD << "[GW_TRANSPORT] " << __FILE__ << ":" << __LINE__ << " -> Gathering subsurface transport timers in finalise()" << std::endl;
#endif
		dom.timers.rtgw = gdom.timers.rtgw;
#endif

// [CODE1] Surface transport timers (total independently timed in compute loop)
#if SERGHEI_SURFACE_TRANSPORT
#endif

// [CODE1] Crop timers already accumulated in dom.timers.crop during compute()
#if CROP_GROWTH_MODEL
#endif

		if (par.masterproc)
		{
			std::cerr << GOK "SIMULATION FINISHED\n";
			std::cerr << GOK "Time elapsed: " << dom.timers.total << std::endl;
		}

		// [CODE2] Timer closure, gathering, and relative computation
		dom.timers.closure();
		dom.timers.gather(par);
		dom.relative = dom.timers;
		dom.relative.computeRelative(dom.timers);
		dom.relative.gather(par);

		io.writeLogFile(dom, par, io.outFolder);
		io.closeOutputStreams();

#if SERGHEI_TOOLS
		if (par.masterproc)
			obs.closeOutputStreams();
#endif
		return 1;
	}

}; // end of class SERGHEI
