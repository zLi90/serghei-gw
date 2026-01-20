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
#include <thread> //!zzb

#include <cstdio>
#include <cstdlib>

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

#if SERGHEI_SUBSURFACE_TRANSPORT
#include "RTState.h"
#include "RTInit.h"
#include "RTFunction.h"
#include "RTMatrix.h"
#include "RTSolver.h"
#include "RTIntegrator.h"
#include "GwMPI.h"
#endif

#if SERGHEI_SURFACE_TRANSPORT
#include "RTStateSW.h"
#include "RTInitSW.h"
#include "RTFunctionSW.h"
#endif

#if CROP_GROWTH_MODEL
#include "./cropsrc/Wofost72.h" // 包含所有 WOFOST 相关头文件
#include "./cropsrc/CropState.h"
#include "./cropsrc/MeteoState.h"
#include "./cropsrc/MeteoInit.h"
#include "./cropsrc/CropInit.h"
#include "./cropsrc/WofostInit.h"
#include "./cropsrc/EvapotranspirationState.h" // 包含蒸散发状态类
#include "SourceSinkCrop.h"					   // 包含作物源汇项处理类
#endif

class SERGHEI
{
public:
	Parallel par;
	Initializer init;
	FileIO io;

	State state;
	Domain dom;
#if SERGHEI_SUBSURFACE_MODEL
	GwState gw;
	GwDomain gdom;
	GwInit ginit;
	GwMPI gmpi;
	GwFunction gwf;
	SubsurfaceBoundaries gbc;
	GwMatrix A;
	GwIntegrator gint;
#ifdef __NVCC__
	GwSolver<Kokkos::Cuda> gsolver;
#else
	GwSolver<Kokkos::OpenMP> gsolver;
#endif
#endif

#if SERGHEI_SUBSURFACE_TRANSPORT
	RTState rt;
	RTInit rtinit;
	RTFunction rtf;
	RTSubsurfaceBoundaries rtgbc;
	RTMatrix rtA;
	RTIntegrator rtint;
#ifdef __NVCC__
	RTSolver<Kokkos::Cuda> rtsolver;
#else
	RTSolver<Kokkos::OpenMP> rtsolver;
#endif
#endif

#if SERGHEI_SURFACE_TRANSPORT
	RTStateSW rtsw;
	RTInitSW rtinitsw;
	RTFunctionSW rtfsw;
#endif

	SourceSink ss;
	// SourceSinkData      ss;

	Exchange exch;
	ExternalBoundaries ebc;

#if CROP_GROWTH_MODEL
	Wofost72 wofost;
	CropState cropParam;					// 静态参数
	MeteoState meteo;						// 气象数据
	MeteoInit meteoInit;					// 气象读取器
	CropInit cropInit;						// 参数读取器
	EvapotranspirationState ets;			// 蒸散发状态
	WofostInit wofostInit;					// WOFOST 初始化器
	realArr current_root_zone_SM;			// 当前根区含水率 (每个地表单元的平均根区含水率)
	int n_crop_cells;						// 作物单元数量 (地表网格单元数)
	double time_since_last_crop_step = 0.0; // 作物模型计时器
	const double CROP_STEP_SIZE = 86400.0;	// 1天 (秒)
	bool wofost_initialized = false;		// WOFOST 初始化标志
#endif

private:
	Parser parser;
	TimeIntegrator tint;
	surfaceIntegrator sint;
	boundaryIntegrator bint;
#if SERGHEI_TOOLS
	Observations obs;
#endif
#if SERGHEI_PARTICLE_TRACKING
	ParticleTracker parTrack;
#endif

	double oldVolume, newVolume, diffVolume;
	double accumDt = 0.0;

	// Kokkos objects
	Kokkos::Timer timer;
	Kokkos::InitializationSettings kokkosSettings;

public:
	std::string inFolder, outFolder;

	////////////// METHODS ///////////////
public:
	int start(int argc, char **argv)
	{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif

		init.initializeMPI(&argc, &argv, par);

#ifdef __NVCC__
		kokkosSettings.set_device_id(par.myrank % par.nthreads);
#else
		if (par.nthreads != 0)
			kokkosSettings.set_num_threads(par.nthreads);
#endif
#if SERGHEI_DEBUG_KOKKOS_SETUP
		printKokkosInitArguments(par);
#if __NVCC__
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

		// Initialize the model
		if (!init.initialize(state, ss.swss, ebc, dom, par, tint, sint, bint, parser, exch, io, inFolder, outFolder))
		{
			std::cerr << RERROR "Unable to start the simulation" << "\n";
			return 0;
		};

// Initialize subsurface model if activated
#if SERGHEI_SUBSURFACE_MODEL
		if (!ginit.initialize_gw(gw, gdom, state, dom, gbc, gmpi, gint, par, io, ss, inFolder, outFolder))
		{
			std::cerr << RERROR "Unable to initialize the subsurface domain" << "\n";
			return 0;
		};

		A.init(gdom);
		gsolver.init(A, gdom);
		if (par.masterproc)
		{
			std::cerr << GOK "Subsurface Solver has been initialized! " << std::endl;
		}
#endif

//! zzb
// Initialize subsurface transport module if activated
#if SERGHEI_SUBSURFACE_TRANSPORT
		if (!rtinit.initialize_rt(rt, gw, gdom, rtgbc, gmpi, par, io, ss, inFolder, outFolder))
		{
			std::cerr << RERROR "Unable to initialize the transport module"
					  << "\n";
			return 0;
		};

		rtA.init(gdom);
		rtsolver.init(rtA, gdom);
#endif
		//! zzb

// Initialize surface transport module if activated
#if SERGHEI_SURFACE_TRANSPORT

		if (!rtinitsw.initialize_rt(rtsw, state, dom, ss.swss, ebc, io, par, inFolder, outFolder))
		{
			std::cerr << RERROR "Unable to initialize the SW transport module" << "\n";
			return 0;
		}
#endif

		// --- WOFOST Initialization ---
#if CROP_GROWTH_MODEL
		if (par.masterproc)
			std::cout << GOK "Initializing WOFOST Crop Model..." << std::endl;

		// 1. Determine number of crop cells (surface grid cells)
		n_crop_cells = dom.nCellMem; // Use surface domain size

		// 2. Allocate memory for root zone moisture array
		current_root_zone_SM = realArr("RootZoneSM", n_crop_cells);

		// 3. Initialize WOFOST model
		if (!wofostInit.initialize_wofost(wofost, cropParam, meteo,
										  cropInit, meteoInit,
										  gw, gdom, dom,
										  par, io, inFolder, outFolder))
		{
			std::cerr << RERROR "Failed to initialize WOFOST" << std::endl;
			return 0;
		}

		// 4. Set initialization flag
		wofost_initialized = true;

		if (par.masterproc)
		{
			std::cout << GOK "WOFOST initialized with " << n_crop_cells << " crop cells" << std::endl;
		}
#endif

#if SERGHEI_TOOLS
		if (!obs.readInputFiles(inFolder, par))
			return 0;
		if (!obs.configure(dom, outFolder))
			return 0; // observations for surface domain
		// obs.printGauges(dom);
		obs.update(state, par, dom);
		if (par.masterproc)
		{
			obs.writeLinesSamplingCoordinates(outFolder);
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
		// integrator at the beginning or the simulation
		sint.integrate(state, dom, ss.swss);
		// Write initial time series data
		io.writeTimeSeriesIni(state, dom, par, ss.swss, sint, bint, ebc.extbc, outFolder);
#if SERGHEI_SUBSURFACE_MODEL
		io.writeSubTimeSeriesIni(gdom, gint, par, outFolder);
		io.outputSubsurface(gw, gdom, par, outFolder);
#endif

#if SERGHEI_SUBSURFACE_TRANSPORT
		io.writeRTSubsurfaceTimeSeriesIni(gdom, rt, rtint, par, outFolder);
		io.outputTransport(rt, gdom, par, outFolder);
#endif

//! zzb
#if CROP_GROWTH_MODEL
		io.writeRootZoneTimeSeriesIni(gdom, gint, par, outFolder, ss.gwss);
		io.outputCrop(wofost, dom, par, outFolder);
#endif
		//! zzb
		// capture initialisation time
		dom.timers.init = timer.seconds();

		if (par.masterproc)
			std::cout << GOK << "Initialisation complete. Initialisation time: " << dom.timers.init << " [s]" << std::endl;
		return 1;
	}

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
#if SERGHEI_SUBSURFACE_MODEL
		gdom.dt = gdom.dt_init;
		gdom.dtOld = gdom.dt_init;
		dom.dt = gdom.dt;
		gdom.cg_iter = 0;

#else
		tint.computeDt(state, dom, io);

#endif
//! zzb
#if SERGHEI_SUBSURFACE_TRANSPORT
		rt.dt = gdom.dt_init;
		rt.dtOld = gdom.dt_init;
		dom.dt = rt.dt;

		// rt.dt = gdom.dt;

#endif

		// Main Time Loop
		while (dom.etime < dom.endTime)
		{

			// previous mass
			oldVolume = sint.surfaceVolumeG;
			bint.integrate(ebc.extbc, dom, 1); // has to be called here (previous time step) with mode==1 (boundary flows)
											   // run surface model

#if SERGHEI_SWE_MODEL
			tint.stepForward(state, ss.swss, ebc.extbc, dom, exch, par, io);
#if SERGHEI_SUBSURFACE_MODEL
			if (dom.isEvap)
			{
				gdom.dt4swEvap = dom.dt;
			}
#endif
#else

			// run source-sink model
			ss.swss.ComputeSWSourceSink(state, dom);

#endif
			//! run surface transport model

#if SERGHEI_SURFACE_TRANSPORT

			rtfsw.transport(rtsw, state, dom, ebc, ss.swss);

#endif

			// #endif

			// run subsurface model
#if SERGHEI_SUBSURFACE_MODEL
#if SERGHEI_SWE_MODEL
			// If both surface and subsurface modules are on
			// Rainfall is first read by the surface module, then copy to the subsurface
			if (gdom.isRain)
			{
				Kokkos::deep_copy(gdom.rainRate, ss.swss.rainRate);
			}
			Kokkos::deep_copy(gw.hs, state.h);
// 将地表水模块中实际水面蒸发量复制到地下水模块中
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
			Kokkos::deep_copy(gw.surfaceWaterEvapActual, state.surfaceEvaporation);
#endif
#endif

#if CROP_GROWTH_MODEL
			// Update the crop timer
			time_since_last_crop_step += dom.dt;

			// Check if it's time to run WOFOST (Daily)
			if (time_since_last_crop_step >= CROP_STEP_SIZE && wofost_initialized)
			{
				// ================== FIX START ==================
				// Calculate Current DOY based on Meteo Start DOY + Simulation Time
				// 1. Calculate how many days have passed since simulation start
				int days_passed = (int)(dom.etime / 86400.0);

				// 2. Calculate current DOY (handling 365 day rollover simply for now)
				// If start is 95, at t=0, doy=95. At t=86400, doy=96.
				int current_doy = ((meteo.start_doy + days_passed - 1) % 365) + 1;
				// ================== FIX END ==================

				if (par.masterproc)
				{
					std::cout << GOK "Running WOFOST Step at t=" << dom.etime
							  << " (Meteo DOY: " << current_doy << ")" << std::endl;
				}

				// 1. Check if RD is valid before computing root zone moisture
				// For the first step, RD might be zero or uninitialized
				bool rd_valid = false;
				double rd_sum = 0.0;

				// Check a few sample cells to see if RD has valid values
				Kokkos::parallel_reduce("CheckRD", 5, KOKKOS_LAMBDA(const int i, double &sum) {
                    if (i < wofost.s.RD.extent(0)) {
                        sum += wofost.s.RD(i);
                    } }, rd_sum);

				if (rd_sum > 1e-6)
				{
					rd_valid = true;
				}

				if (rd_valid)
				{
					// 1. Calculate average root zone moisture from SERGHEI state
					if (par.masterproc)
					{
						std::cout << GOK "About to compute root zone moisture with WOFOST RD" << std::endl;
					}
					SourceSinkCrop::compute_root_zone_moisture(gw, gdom, wofost.s.RD, current_root_zone_SM);

					if (par.masterproc)
					{
						std::cout << GOK "Root zone moisture computed using WOFOST RD" << std::endl;
					}
				}
				else
				{
					// First step: use initial root depth from crop parameters
					if (par.masterproc)
					{
						std::cout << GOK "First WOFOST step: using initial root depth RDI=" << cropParam.p.RDI << " cm" << std::endl;
					}

					// Create temporary RD array with initial root depth
					realArr temp_RD("TempRD", gdom.nCellSwMem);
					Kokkos::parallel_for("InitRD", gdom.nCellSwMem, KOKKOS_LAMBDA(const int i) {
						temp_RD(i) = cropParam.p.RDI; // cm
					});

					if (par.masterproc)
					{
						std::cout << GOK "About to compute root zone moisture with initial RD" << std::endl;
					}
					SourceSinkCrop::compute_root_zone_moisture(gw, gdom, temp_RD, current_root_zone_SM);
					if (par.masterproc)
					{
						std::cout << GOK "Root zone moisture computed with initial RD" << std::endl;
					}
				}

				// 2. Prepare time variables
				// int current_doy = (int)(dom.etime / 86400.0) % 365 + 1;
				if (par.masterproc)
				{
					std::cout << GOK "About to run WOFOST calc_rates_and_integrate" << std::endl;
				}

				// 3. Run WOFOST daily calculation
				wofost.calc_rates_and_integrate(cropParam, meteo, dom.etime,
												current_root_zone_SM, current_doy, 1.0);

				if (par.masterproc)
				{
					std::cout << GOK "WOFOST calc_rates_and_integrate completed" << std::endl;
				}

				// 4. Reset timer for next daily step
				time_since_last_crop_step = 0.0;

				if (par.masterproc)
				{
					// 输出第一个单元的 LAI 和 RD 值（因为它们是数组）
					double lai_val = wofost.s.LAI(0);
					double rd_val = wofost.s.RD(0);
					double TAGP_val = wofost.TAGP(0);
					std::cout << GOK "WOFOST step completed. LAI[0]=" << lai_val
							  << ", RD[0]=" << rd_val << ", TAGP[0]=" << wofost.TAGP(0) << std::endl;
				}

#if CROP_GROWTH_MODEL
				io.outputCrop(wofost, dom, par, outFolder);
#endif
			}
#endif

			//! 20250923
#if SERGHEI_SUBSURFACE_TRANSPORT
#if SERGHEI_SURFACE_TRANSPORT
			Kokkos::deep_copy(rt.csw, rtsw.csw4gw); //! 将sw网格中溶质浓度复制到gw溶质class中
#endif
#endif
				//! 20250923
			// Asynchronous coupling
			// timer.reset();
			if (gdom.async)
			{
				if (gdom.etime + gdom.dt < dom.etime)
				{
					gdom.etime += gdom.dt;
#ifdef __NVCC__
					//	PC scheme
					if (gdom.gw_scheme == 1)
					{
						gwf.pca_solve<Kokkos::Cuda>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
					}
					//	Modified Picard scheme
					else
					{
						gwf.picard_solve<Kokkos::Cuda>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
					}

					//! 20250923
#if SERGHEI_SUBSURFACE_TRANSPORT
					rtf.rt_solve<Kokkos::Cuda>(rt, rtA, gw, gdom, rtgbc.rtgwbc, gmpi, par, rtint, rtsolver);
#endif
					//! 20250923
#else
					if (gdom.gw_scheme == 1)
					{
						gwf.pca_solve<Kokkos::OpenMP>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
					}
					else
					{
						gwf.picard_solve<Kokkos::OpenMP>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
					}
					//! 20250923
#if SERGHEI_SUBSURFACE_TRANSPORT
					rtf.rt_solve<Kokkos::OpenMP>(rt, rtA, gw, gdom, rtgbc.rtgwbc, gmpi, par, rtint, rtsolver);
					// rtf.rt_solve<Kokkos::OpenMP>(rt, rtA, gw, gdom, rtgbc.rtgwbc, gmpi, par, rtsolver);
#endif

					//! 20250923
#endif
				}
			}
			else
			{ //! zzb 同步情况
				gdom.etime = dom.etime;
#ifdef __NVCC__
				if (gdom.gw_scheme == 1)
				{
					gwf.pca_solve<Kokkos::Cuda>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
				}
				else
				{
					gwf.picard_solve<Kokkos::Cuda>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
				}

				//! 20250923
#if SERGHEI_SUBSURFACE_TRANSPORT
				rtf.rt_solve<Kokkos::Cuda>(rt, rtA, gw, gdom, rtgbc.rtgwbc, gmpi, par, rtint, rtsolver);
#endif
				//! 20250923

#else //! zzb cpu计算

				//! zzb 当dom.etime=86400整数倍时，运行ginit.initialize_gw函数，为了每天更新gwss中input文件
				// #if CROP_GROWTH_MODEL

				// 				if (
				// 					// std::fabs(std::fmod(dom.etime, 86400.0)) == 0.0
				// 					// ||
				// 					// std::fabs(std::fmod(dom.etime, 86400.0) - 86400.0) <= 50.0)
				// 					fabs(dom.etime - dom.startTime - io.numOut * io.outFreq) <= TOL12)
				// 				{

				// 					//! zzb 1.调用python脚本运行crop_growth_model
				// 					// 定义要执行的 Python 脚本的命令
				// 					const char *cmd = "/bin/python3 /home/zhangzhibo/serghei/serghei-20240910-RTM/src/wofost/scripts/replace_values_in_etinput.py";

				// 					// 使用 popen 打开子进程
				// 					FILE *pipe = popen(cmd, "r");
				// 					if (!pipe)
				// 					{
				// 						perror("popen failed");
				// 						return EXIT_FAILURE;
				// 					}

				// 					// 读取 Python 脚本的输出
				// 					const size_t buffer_size = 256;
				// 					char buffer[buffer_size];
				// 					while (fgets(buffer, buffer_size, pipe) != nullptr)
				// 					{
				// 						// 处理输出，例如打印到控制台
				// 						std::printf("%s", buffer);
				// 					}

				// 					// 关闭子进程
				// 					int status = pclose(pipe);
				// 					if (status == -1)
				// 					{
				// 						perror("pclose failed");
				// 						return EXIT_FAILURE;
				// 					}

				// 					// 检查子进程的退出状态
				// 					if (WIFEXITED(status))
				// 					{
				// 						std::printf("Python script exited with status %d\n", WEXITSTATUS(status));
				// 					}
				// 					else if (WIFSIGNALED(status))
				// 					{
				// 						std::printf("Python script killed by signal %d\n", WTERMSIG(status));
				// 					}
				// 					// return EXIT_SUCCESS;

				// 					//! zzb 2.逐日读取更新gwss中input文件
				// 					std::string fNameIn = inFolder + "gwss.input";
				// 					if (!ginit.readGwSSFile(fNameIn, gw, gdom, ss, par))
				// 					{
				// 						std::cerr << RERROR "Unable to run the ginit.readGwSSFile when crop_growth_model on" << "\n";
				// 						return 0;
				// 					}
				// 				}

				// #endif
				if (gdom.gw_scheme == 1)
				{
					gwf.pca_solve<Kokkos::OpenMP>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
				}
				else
				{
					gwf.picard_solve<Kokkos::OpenMP>(gw, gdom, gbc.gwbc, A, gsolver, ss.gwss, gmpi, gint, par);
				}
				//! 20250923
#if SERGHEI_SUBSURFACE_TRANSPORT
				rtf.rt_solve<Kokkos::OpenMP>(rt, rtA, gw, gdom, rtgbc.rtgwbc, gmpi, par, rtint, rtsolver);
// rtf.rt_solve<Kokkos::OpenMP>(rt, rtA, gw, gdom, rtgbc.rtgwbc, gmpi, par, rtsolver);
#endif
				//! 20250923
#endif
			}
			gdom.cg_iter += A.cg_iter;
			// gdom.timers.gw += timer.seconds();
			// surface-subsurface exchange
#if SERGHEI_SWE_MODEL
			Kokkos::deep_copy(state.qss, gw.qss);
			tint.computeGwExchange(state, dom);
#endif

#if SERGHEI_SURFACE_TRANSPORT
			Kokkos::deep_copy(rtsw.ConQss, rt.ConQss);
			rtfsw.computeGwExchange_For_RT(state, rtsw, dom);
#endif

#endif

			// printf("11111\n");

			oldVolume += (bint.inflowDischargeG - bint.outflowDischargeG) * dom.dt; // Boundary fluxes with the new dt
			bint.integrate(ebc.extbc, dom, 0);										// called here with mode==0 (adjusted volume)
			oldVolume += bint.adjustedVolumeG;										// Some mass changes can occur through the boundaries
			sint.integrate(state, dom, ss.swss);									// new mass after the new time step integration
			oldVolume += (sint.rainFluxG - sint.infFluxG) * dom.dt;					// after integrate, we have to sum the rain and inf mass
			newVolume = sint.surfaceVolumeG;

			if (fabs(oldVolume) > TOL12)
			{
				diffVolume = (newVolume - oldVolume) / oldVolume * 100.;
			}
			else
			{
				diffVolume = 0.0;
			}
			dom.etime += dom.dt;
			dom.nIter++;
			dom.countIterDt++;
			accumDt += dom.dt;

			if (dom.nIter % io.nScreen == 0 || fabs(dom.etime - dom.startTime - io.numOut * io.outFreq) <= TOL12)
			{
				if (par.masterproc)
				{
					std::cerr << std::fixed;
					std::cerr << GSTAR "TIME: " << dom.etime << " average dt: " << accumDt / dom.countIterDt << "\n";
					std::cerr.precision(9);
					std::cerr << std::scientific;
					std::cerr << std::fixed;
					std::cerr.precision(12);
#if SERGHEI_SWE_MODEL
					std::cerr << "     Surface Volume:\t" << newVolume << "\n";
					std::cerr << "     Surface inflow: " << bint.inflowDischargeG << "\n";
					std::cerr << "     Surface outflow Volume: " << bint.outflowDischargeG * dom.dt << "\n";
#if SERGHEI_SUBSURFACE_MODEL
					std::cerr << "     Exchange Volume: " << gint.Vexch_glob << "\n";
#endif
#endif

					if (fabs(diffVolume) > TOL_MASS_ERROR)
					{
						std::cerr << YEXC "   Old Volume:\t" << oldVolume << "\n";
						std::cerr << YEXC "   New Volume:\t" << newVolume << "\n";
						std::cerr << YEXC "   Diff Volume:\t" << newVolume - oldVolume << "\n";
						std::cerr << YEXC "   Inflow Volume:\t" << bint.inflowDischargeG * dom.dt << "\n";
						std::cerr << YEXC "   Outflow Volume:\t" << bint.outflowDischargeG * dom.dt << "\n";
						std::cerr << YEXC "   Adjusted Volume:\t" << bint.adjustedVolumeG << "\n";
						std::cerr << YEXC "   Rain Volume:\t" << sint.rainFluxG * dom.dt << "\n";
						std::cerr << YEXC "   Inf Volume:\t" << sint.infFluxG * dom.dt << "\n";
#if SERGHEI_DEBUG_MASS_CONS > 1
						getchar();
#endif
					}
				}
				if (fabs(dom.etime - dom.startTime - io.numOut * io.outFreq) <= TOL12)
				{
#if SERGHEI_SWE_MODEL
					io.output(state, dom, ss.swss, par, outFolder);
#if SERGHEI_SURFACE_TRANSPORT
					io.outputTransportSW(rtsw, dom, par, outFolder); // todo
#endif
#endif
#if SERGHEI_SUBSURFACE_MODEL
					io.outputSubsurface(gw, gdom, par, outFolder);
#endif

					// #if CROP_GROWTH_MODEL
					// 					io.outputCrop(wofost, dom, par, outFolder);
					// #endif

#if SERGHEI_SUBSURFACE_TRANSPORT
					io.outputTransport(rt, gdom, par, outFolder);
#endif

					//! zzb
					if (par.masterproc)
						std::cerr << GIO "File " << io.numOut - 1 << " written" << std::endl; // io.numOut already updated
				}
				if (par.masterproc)
					std::cerr << "-------------------------------------------------\n";
				dom.countIterDt = 0;
				accumDt = 0.0;
			}

#if SERGHEI_PARTICLE_TRACKING
			parTrack.update(dom, state);
#endif

			if (dom.etime >= io.numObs * io.obsFreq)
			{
#if SERGHEI_TOOLS
				obs.update(state, par, dom);
#endif
				io.writeTimeSeries(state, dom, par, sint, bint, ebc.extbc);
#if SERGHEI_SUBSURFACE_MODEL
				io.writeSubsurfaceTimeSeries(gdom, gint);
#endif
//! zzb
#if SERGHEI_SUBSURFACE_TRANSPORT
				io.writeRTSubsurfaceTimeSeries(gdom, rtint);
#endif
#if CROP_GROWTH_MODEL
				io.writeRootZoneWaterContentTimeSeries(gdom, gint, ss.gwss);
#endif
				//! zzb
				if (par.masterproc)
				{
#if SERGHEI_TOOLS
					obs.write(dom);
#endif
				}
			}

// Unify dt for coupled simulations
#if SERGHEI_SWE_GW

			tint.computeDt(state, dom, io);
#if SERGHEI_SUBSURFACE_MODEL
			if (!gdom.async)
			{ // 同步
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
#endif

#elif SERGHEI_SUBSURFACE_MODEL
			dom.dt = gdom.dt;
			tint.dtMatchOutput(dom, io);
			gdom.dt = dom.dt;
#elif !SERGHEI_SWE_MODEL
			std::cout << RERROR << "Impossible configuration without SWE nor GW model" << std::endl;
			return 0;

#if SERGHEI_SUBSURFACE_TRANSPORT

			// 选取gdom.dt和rt.dt中较小的作为rt.dt
			//  if (rt.dt < dom.dt)
			//  {
			//  	dom.dt = rt.dt;
			//  }

			// else
			// {
			// 	rt.dt = dom.dt;
			// }
			//! zzb 直接规定rt.dt值（debug）
			// rt.dt = 20.0;

			rt.dt = gdom.dt;

#endif
#endif
			// printf("dom.etime=%f, dom.dt=%f, gdom.dt=%f, rt.dt=%f, rtsw.dt=%f\n",dom.etime, dom.dt, gdom.dt, rt.dt, rtsw.dt);

		} // end of time loop

		return 1;
	}

	int finalise()
	{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		dom.timers.total = timer.seconds();
		if (par.masterproc)
		{
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
		io.closeOutputStreams();
#if SERGHEI_TOOLS
		if (par.masterproc)
			obs.closeOutputStreams();
#endif
		return 1;
	}
};
