#ifndef _FILEIO_H_
#define _FILEIO_H_

#include <chrono>
#include <ctime>
#include <unistd.h>
#ifndef KOKKOS_ENABLE_CUDA
#include <cpuid.h>
#endif
#include "const.h"
#include <iostream>
#include "define.h"
#include "State.h"
#include "SourceSink.h"
#include "BC.h"
#include "mpi.h"
#include "Indexing.h"
#include "DomainIntegrator.h"
#include "tools.h"
#include "ParticleTracking.h"
#include "netcdfIO.h"
#include "units.h"

#include "GwState.h"
#include "GwDomain.h"
#include "GwIntegrator.h"

// [CODE 1] 核心模块的引入：溶质运移与作物模型
#include "RTStateGW.h"
#include "RTStateSW.h"
#include "RTIntegrator.h"
#include "RTIntegratorSW.h"

#if CROP_GROWTH_MODEL
#include "./cropsrc/Wofost72.h"
#endif

class Config
{
public:
#if SERGHEI_SCALAR_TRANSPORT
  ScalarTransport st;
#endif
#if SERGHEI_SEDIMENT_TRANSPORT
  SedimentTransport sedt;
#endif
};

#ifndef SERGHEI_NC_MODE
#define SERGHEI_NC_MODE NC_CLOBBER
#endif

#ifndef SERGHEI_WRITE_HZ
#define SERGHEI_WRITE_HZ 0
#endif
#ifndef SERGHEI_WRITE_SUBDOMS
#define SERGHEI_WRITE_SUBDOMS 0
#endif

#ifndef SERGHEI_NC_REAL
#if SERGHEI_REAL == SERGHEI_DOUBLE
#define SERGHEI_NC_REAL NC_DOUBLE
#elif SERGHEI_REAL == SERGHEI_FLOAT
#define SERGHEI_NC_REAL NC_FLOAT
#endif
#endif

#if SERGHEI_NC_REAL == NC_DOUBLE
#define ncmpi_put_vara_real_all ncmpi_put_vara_double_all
#define ncmpi_put_att_real ncmpi_put_att_double
#endif
#if SERGHEI_NC_REAL == NC_FLOAT
#define ncmpi_put_vara_real_all ncmpi_put_vara_float_all
#define ncmpi_put_att_real ncmpi_put_att_float
#endif

#define SERGHEI_NC_FILL_VALUE_KEY "_FillValue"
#ifndef SERGHEI_NC_ENABLE_NAN
#define SERGHEI_NC_ENABLE_NAN 1
#endif

#ifndef SERGHEI_INPUT_NETCDF
#define SERGHEI_INPUT_NETCDF 0
#endif

#if SERGHEI_NC_REAL == NC_FLOAT && SERGHEI_REAL == SERGHEI_DOUBLE
typedef Kokkos::View<float *, Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::SharedSpace>> ncArr;
#else
typedef Kokkos::View<real *, Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::SharedSpace>> ncArr;
#endif

#if SERGHEI_NC_REAL == NC_FLOAT
typedef float ncreal;
#elif SERGHEI_NC_REAL == NC_DOUBLE
typedef double ncreal;
#endif

class FileIO
{
public:
  std::string inFolder, outFolder;
  ncStream ncin;
  ncStreamStrings ncinS;
  bool writeRain = 0;
  bool writeRainA = 0;
#if SERGHEI_NETCDF_FORCING
  ncStreamStrings aforcingS;
#endif
  bool readLandUse = 0;
  bool writeLandUse = 0;
  bool writeSoilMap = 0;
  bool writeRoughness = 0;
  bool writeInfParameters = 0;
  ShallowWater sw;

protected:
  int ncid;
  int tDim, xDim, yDim, zDim, nDim;
  int tVar, xVar, yVar, hVar, hzVar, uVar, vVar, zVar, z3Var, hdVar, wcVar, qVar;
#if SERGHEI_SEDIMENT_TRANSPORT
  int intTauDtVar, intTau32DtVar, dzVar, bedExVar;
#endif
  int infVar, infVolVar;
  int rainVar, rainAVar;
  int roughnessVar;
  int landuseVar, soilmapVar;
  int cinfCapVar;

  // [CODE 1] 新增功能变量支持
  //! 地下水溶质运移与水文属性
  int nVar, cVar, c_solidVar, qxVar, qyVar, qzVar;
  //! 地表水面蒸发
  int evapVar;
  //! 地下水蒸发与植物蒸腾
  int evapGWVar, transpGWVar;
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
  int c_difxxVar, c_difzzVar, c_difxzVar, c_difzxVar;
#endif
  // 地表水溶质运移
  int cSWVar;
  //! 作物输出
  int cropTimeVar, laiVar, rdVar, wsoVar, tagpVar, dvsVar, wlvVar, wstVar, wrtVar, hiVar, totalVar, rftraVar, gassVar, asrcVar, rzsmVar;

  std::ofstream domainOutputFile;
  std::ofstream SubsurfaceOutputFile;
  std::ofstream RootZoneWaterContentOutputFile; // [CODE 1]
  std::ofstream RTSubsurfaceOutputFile;         // [CODE 1]
  std::ofstream RTSurfaceOutputFile;            // [CODE 1]
  std::ofstream logFile;
  int nOut;           // expected number of spatial outputs
  int numOutCrop = 0; //[CODE 1]

#if SERGHEI_MAXFLOOD > 0
  int hMaxVar, momMaxVar, timehMaxVar;
#endif
#if SERGHEI_SCALAR_TRANSPORT
  int *phiVar;
#endif
#if SERGHEI_WRITE_SUBDOMS
  int subdomVar;
#endif

#if SERGHEI_LPT
  std::string lpt_filename;
  int tpDim, pDim;
  int ncidP, pVar, tpVar, xpVar, ypVar, zpVar;
  std::ofstream particlesFile;
#endif

private:
  Kokkos::Timer timer;
#if SERGHEI_LPT
  Kokkos::Timer timer_particles;
#endif

public:
  std::string infMapFile;
  std::string infFile;
  std::string ncInputFilename;
  bool allowIni = 1;
  real outFreq;
  int outFormat;
  int numOut = 0; // spatial output counter

  real obsFreq;
  int numObs = 0;
  int nScreen = 100; // if this value is not specified, the information is displayed every 1000 iterations

  Config conf;

  // Writes spatial fields for initial state
  void outputIni(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par)
  {
    nOut = floor(dom.simLength / outFreq);
    numOut = 0;
    if (outFormat == OUT_NETCDF)
    {
      outputInitNETCDF(state, dom, ss, par, outFolder);
    }
    if (outFormat == OUT_VTK)
    {
      outputVTK(state, dom, ss, par, outFolder);
    }
    if (outFormat == OUT_BIN)
    {
      initBIN(state, dom, par, outFolder);
      outputBIN(state, dom, par, outFolder);
    }
    numOut++;
  }

  // [CODE 1] Surface transport
#if SERGHEI_SURFACE_TRANSPORT
  void outputIniRTSW(const RTStateSW &rtsw, const State &state, Domain const &dom, Parallel const &par, std::string dir)
  {
    numOut = 0;
    nOut = floor(dom.simLength / outFreq);
    outputInitNETCDFRTSW(rtsw, state, dom, par, dir);
    numOut++;
  }

  void outputTransportSW(const RTStateSW &rtsw, const State &state, Domain const &dom, Parallel const &par, std::string dir)
  {
    timer.reset();
    numOut--;
    outputNETCDFTransportSW(rtsw, state, dom, par, dir);
    numOut++;
    dom.timers.swe.io.out += timer.seconds(); // Adapting to Code 2's timer system
  }
#endif

  //[CODE 1] Subsurface transport
#if SERGHEI_SUBSURFACE_TRANSPORT
  void outputIniRT(const RTStateGW &rt, GwDomain const &gdom, Parallel const &par, std::string dir)
  {
    nOut = floor(gdom.simLength / outFreq);
    outputInitNETCDFRT(rt, gdom, par, dir);
  }
  void outputTransport(const RTStateGW &rt, GwDomain const &gdom, Parallel const &par, std::string dir)
  {
    timer.reset();
    numOut--;
    outputNETCDFTransport(rt, gdom, par, dir);
    numOut++;
#if SERGHEI_RE_MODEL
    gdom.timers.re.out += timer.seconds();
#endif
  }
#endif

  // [CODE 1] Crop Growth
#if CROP_GROWTH_MODEL
  void outputIniCrop(const Wofost72 &wofost, Domain const &dom, Parallel const &par, std::string dir)
  {
    nOut = floor(dom.simLength / outFreq);
    numOutCrop = 0;
    if (outFormat == OUT_NETCDF)
    {
      outputInitNETCDFCrop(wofost, dom, par, dir);
    }
    numOutCrop++;
  }

  void outputCrop(const Wofost72 &wofost, Domain const &dom, Parallel const &par, std::string dir)
  {
    timer.reset();
    if (outFormat == OUT_NETCDF)
    {
      outputNETCDFCrop(wofost, dom, par, dir);
    }
    numOutCrop++;
    dom.timers.swe.io.out += timer.seconds();
  }
#endif

#if SERGHEI_RE_MODEL
  void outputIniSub(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir)
  {
    nOut = floor(gdom.simLength / outFreq);
    outputInitNETCDFSub(gw, gdom, par, dir);
  }

  void outputSubsurface(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir)
  {
    timer.reset();
#if SERGHEI_SWE_MODEL
    numOut--;
#endif
    outputNETCDFSubsurface(gw, gdom, par, dir);
    numOut++;
    gdom.timers.re.out += timer.seconds();
  }
#endif

// Writes spatial fields
#if SERGHEI_LPT == 0
  void output(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par)
  {
#else
  void output(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir, ParticleTracker const &particles)
  {
#endif
    timer.reset();

    if (outFormat == OUT_NETCDF)
    {
      outputNETCDF(state, dom, ss, par, outFolder);
#if SERGHEI_LPT
      timer_particles.reset();
#if SERGHEI_PARTICLE_NO_OUTPUT
#else
      outputNETCDF_particles(particles, dom, par, dir); // Particle Tracking
#endif
      dom.timers.lpt.out += timer_particles.seconds();
#endif
    }
    if (outFormat == OUT_VTK)
    {
      outputVTK(state, dom, ss, par, outFolder);

#if SERGHEI_LPT
      timer_particles.reset();
#if SERGHEI_PARTICLE_NO_OUTPUT
#else
      outputParticle(dom, particles, dir); // Particle Tracking
#endif
      dom.timers.lpt.out += timer_particles.seconds();
#endif
    }
    if (outFormat == OUT_BIN)
    {
      outputBIN(state, dom, par, outFolder);
    }

#if SERGHEI_LPT
    timer_particles.reset();
#if SERGHEI_PARTICLE_NO_OUTPUT
#else
#if SERGHEI_LPT
    writeParticle_TemporalFile(dom, particles, dir);
#endif
#endif
    dom.timers.lpt.out += timer_particles.seconds();
#endif

    numOut++;
    dom.timers.swe.io.out += timer.seconds();
#if SERGHEI_LPT
    // This is not a nice solution
    dom.timers.swe.io.out -= dom.timers.lpt.out;
#endif
  }

  template <typename T, typename TNC>
  void writeNetCDFfield(const Domain &dom, int ncid, int ncvar, MPI_Offset *st, MPI_Offset *ct, const boolArr &mask, const T &myview, TNC &data)
  {

    constexpr auto nodata = []
    {
      if constexpr (std::is_same<typename T::value_type, real>::value)
        return static_cast<real>(SERGHEI_NAN);
      else if constexpr (std::is_same<typename T::value_type, float>::value)
        return static_cast<float>(SERGHEI_NAN);
      else if constexpr (std::is_same<typename T::value_type, unsigned short>::value)
        return static_cast<ushort>(NC_FILL_USHORT);
      else
        static_assert(!std::is_same<T, T>::value, "nodata is not defined for this type");
    }();

    Kokkos::parallel_for("writeNCDFfield", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
      int ii = dom.getIndex(iGlob);
      int jj = iGlob;
      // check if the view and the data buffer have the same size. If they don't, the view is dom.nCellMem and the indices need mapping.
      if (data.extent(0) != myview.extent(0))
        jj = ii;
      data(iGlob) = myview(jj);
#if SERGHEI_NC_ENABLE_NAN
      if (mask(ii))
        data(iGlob) = nodata; // mask is of size dom.nCellMem, so, the indices need mapping for sure.
#endif
    });
    Kokkos::fence();

    if constexpr (std::is_same<T, realArr>::value)
      ncwrap(ncmpi_put_vara_real_all(ncid, ncvar, st, ct, data.data()), __LINE__);
    if constexpr (std::is_same<T, ushortArr>::value)
      ncwrap(ncmpi_put_vara_ushort_all(ncid, ncvar, st, ct, data.data()), __LINE__);
  };

  // NetCDF initialiser and writer
  void outputInitNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir)
  {
    int dimids[3];
    MPI_Offset st[3], ct[3];
    doubleArr xCoord = doubleArr("xCoord", dom.nx);
    doubleArr yCoord = doubleArr("yCoord", dom.ny);
    ncArr data = ncArr("data", dom.nCell);
    static char timeUnits[] = "seconds";
    std::string filename;
    std::string longname;
    std::string units;

    filename = dir + "output.nc";

#if SERGHEI_SCALAR_TRANSPORT
    phiVar = (int *)malloc(conf.st.nScalar * sizeof(int));
#endif

    // Create the file
    ncwrap(ncmpi_create(MPI_COMM_WORLD, filename.c_str(), SERGHEI_NC_MODE, MPI_INFO_NULL, &ncid), __LINE__, par.myrank);

    // Create the dimensions
    ncwrap(ncmpi_def_dim(ncid, "t", (MPI_Offset)NC_UNLIMITED, &tDim), __LINE__, par.myrank);
    ncwrap(ncmpi_def_dim(ncid, "x", (MPI_Offset)dom.nx_glob, &xDim), __LINE__, par.myrank);
    ncwrap(ncmpi_def_dim(ncid, "y", (MPI_Offset)dom.ny_glob, &yDim), __LINE__, par.myrank);
    // Create the variables
    dimids[0] = tDim;
    ncwrap(ncmpi_def_var(ncid, "t", NC_DOUBLE, 1, dimids, &tVar), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, tVar, "units", strlen(timeUnits), timeUnits), __LINE__, par.myrank);
    dimids[0] = xDim;
    ncwrap(ncmpi_def_var(ncid, "x", NC_DOUBLE, 1, dimids, &xVar), __LINE__, par.myrank);
    dimids[0] = yDim;
    ncwrap(ncmpi_def_var(ncid, "y", NC_DOUBLE, 1, dimids, &yVar), __LINE__, par.myrank);

    // time dependent variables
    dimids[0] = tDim;
    dimids[1] = yDim;
    dimids[2] = xDim;
    ncwrap(ncmpi_def_var(ncid, "h", SERGHEI_NC_REAL, 3, dimids, &hVar), __LINE__, par.myrank);
#if SERGHEI_WRITE_HZ
    ncwrap(ncmpi_def_var(ncid, "h+z", SERGHEI_NC_REAL, 3, dimids, &hzVar), __LINE__, par.myrank);
#endif
    ncwrap(ncmpi_def_var(ncid, "u", SERGHEI_NC_REAL, 3, dimids, &uVar), __LINE__, par.myrank);
    ncwrap(ncmpi_def_var(ncid, "v", SERGHEI_NC_REAL, 3, dimids, &vVar), __LINE__, par.myrank);
    if (ss.inf.model)
    {
      ncwrap(ncmpi_def_var(ncid, "inf", SERGHEI_NC_REAL, 3, dimids, &infVar), __LINE__, par.myrank);
      ncwrap(ncmpi_def_var(ncid, "infVol", SERGHEI_NC_REAL, 3, dimids, &infVolVar), __LINE__, par.myrank);
    }

    // [CODE 1] 表面蒸发
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
  
      ncwrap(ncmpi_def_var(ncid, "surfaceEvaporation", SERGHEI_NC_REAL, 3, dimids, &evapVar), __LINE__, par.myrank);
    
#endif

#if SERGHEI_RE_MODEL
    ncwrap(ncmpi_def_var(ncid, "qss", SERGHEI_NC_REAL, 3, dimids, &qVar), __LINE__);
#endif

    if (writeRain)
      ncwrap(ncmpi_def_var(ncid, "rain", SERGHEI_NC_REAL, 3, dimids, &rainVar), __LINE__, par.myrank);
#if SERGHEI_NETCDF_FORCING
    if (writeRainA)
      ncwrap(ncmpi_def_var(ncid, "rainAccum", SERGHEI_NC_REAL, 3, dimids, &rainAVar), __LINE__, par.myrank);
#endif

    // variables constant in time
    int dimids2[2];
    dimids2[0] = yDim;
    dimids2[1] = xDim;
    if (writeRoughness)
    {
      ncwrap(ncmpi_def_var(ncid, "roughness", SERGHEI_NC_REAL, 2, dimids2, &roughnessVar), __LINE__, par.myrank);
      longname.assign("Hydraulic roughness");
      ncwrap(ncmpi_put_att_text(ncid, roughnessVar, "long_name", longname.length(), longname.c_str()), __LINE__);
      units.assign("s m^(-1/3)");
      ncwrap(ncmpi_put_att_text(ncid, roughnessVar, "units", units.length(), units.c_str()), __LINE__);
    }

    if constexpr (SERGHEI_MAPCLASS)
    {
      ushort nanvalue = NC_FILL_USHORT;
      if (writeLandUse && state.landuse.use)
      {
        if constexpr (SERGHEI_DEBUG_MAPCLASS)
          std::cout << GGD << __PRETTY_FUNCTION__ << "::" << __LINE__ << ". Writing land use map to netcdf" << std::endl;
        ncwrap(ncmpi_def_var(ncid, "landuse", NC_USHORT, 2, dimids2, &landuseVar), __LINE__, par.myrank);
        longname.assign("Land use class");
        ncwrap(ncmpi_put_att_text(ncid, landuseVar, "long_name", longname.length(), longname.c_str()), __LINE__);
        units.assign("-");
        ncwrap(ncmpi_put_att_ushort(ncid, landuseVar, SERGHEI_NC_FILL_VALUE_KEY, NC_USHORT, 1, &nanvalue), __LINE__);
      }

      if (writeSoilMap && ss.inf.soilmap.use)
      {
        if constexpr (SERGHEI_DEBUG_MAPCLASS)
          std::cout << GGD << __PRETTY_FUNCTION__ << "::" << __LINE__ << ". Writing soil map to netcdf" << std::endl;

        ncwrap(ncmpi_def_var(ncid, "soilmap", NC_USHORT, 2, dimids2, &soilmapVar), __LINE__, par.myrank);
        longname.assign("Soil class");
        ncwrap(ncmpi_put_att_text(ncid, soilmapVar, "long_name", longname.length(), longname.c_str()), __LINE__);
        ncwrap(ncmpi_put_att_ushort(ncid, soilmapVar, SERGHEI_NC_FILL_VALUE_KEY, NC_USHORT, 1, &nanvalue), __LINE__);
      }
    }

    if (writeInfParameters)
    {
      if (ss.inf.model == INF_CONSTANT)
      {
        ncwrap(ncmpi_def_var(ncid, "cinfCap", SERGHEI_NC_REAL, 2, dimids2, &cinfCapVar), __LINE__, par.myrank);
        longname.assign("Constant infiltration capacity");
        ncwrap(ncmpi_put_att_text(ncid, cinfCapVar, "long_name", longname.length(), longname.c_str()), __LINE__);
        units.assign("mm/s");
        ncwrap(ncmpi_put_att_text(ncid, cinfCapVar, "units", units.length(), units.c_str()), __LINE__);
      }
      else
      {
        std::cout << YEXC << "NetCDF output for infiltration parameters other than constant infiltration has not yet been implemented." << std::endl;
      }
    }

#if SERGHEI_MAXFLOOD > 0
    int nc_ndims = 2 + SERGHEI_MAXFLOOD - 1;
    if (nc_ndims == 2)
    {
      dimids[0] = yDim;
      dimids[1] = xDim;
    }
    ncwrap(ncmpi_def_var(ncid, "hMax", SERGHEI_NC_REAL, nc_ndims, dimids, &hMaxVar), __LINE__, par.myrank);
    ncwrap(ncmpi_def_var(ncid, "momMax", SERGHEI_NC_REAL, nc_ndims, dimids, &momMaxVar), __LINE__, par.myrank);
    ncwrap(ncmpi_def_var(ncid, "timehMax", SERGHEI_NC_REAL, nc_ndims, dimids, &timehMaxVar), __LINE__, par.myrank);

    longname.assign("Maximum water depth");
    ncwrap(ncmpi_put_att_text(ncid, hMaxVar, "long_name", longname.length(), longname.c_str()), __LINE__);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncid, hMaxVar, "units", units.length(), units.c_str()), __LINE__);

    longname.assign("Maximum momentum");
    ncwrap(ncmpi_put_att_text(ncid, momMaxVar, "long_name", longname.length(), longname.c_str()), __LINE__);
    units.assign("m^2/s");
    ncwrap(ncmpi_put_att_text(ncid, momMaxVar, "units", units.length(), units.c_str()), __LINE__);

    longname.assign("Time to maximum depth");
    ncwrap(ncmpi_put_att_text(ncid, timehMaxVar, "long_name", longname.length(), longname.c_str()), __LINE__);
    units.assign("s");
    ncwrap(ncmpi_put_att_text(ncid, timehMaxVar, "units", units.length(), units.c_str()), __LINE__);
#endif

#if SERGHEI_WRITE_SUBDOMS
    dimids[0] = yDim;
    dimids[1] = xDim;
    ncwrap(ncmpi_def_var(ncid, "subdom", NC_INT, 2, dimids, &subdomVar), __LINE__);
#endif

#if SERGHEI_EROSIVE_SHEAR
    dimids[0] = tDim;
    dimids[1] = yDim;
    dimids[2] = xDim;
    ncwrap(ncmpi_def_var(ncid, "intTauDt", SERGHEI_NC_REAL, 3, dimids, &intTauDtVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "phiAcc", SERGHEI_NC_REAL, 3, dimids, &intTau32DtVar), __LINE__);
#endif

#if SERGHEI_SCALAR_TRANSPORT
    dimids[0] = tDim;
    dimids[1] = yDim;
    dimids[2] = xDim;
    for (int iphi = 0; iphi < conf.st.nScalar; iphi++)
    {
      ncwrap(ncmpi_def_var(ncid, conf.st.scalarNames[iphi].c_str(), SERGHEI_NC_REAL, 3, dimids, &phiVar[iphi]), __LINE__);
    }
#endif

#if SERGHEI_SEDIMENT_TRANSPORT
    dimids[0] = tDim;
    dimids[1] = yDim;
    dimids[2] = xDim;
    ncwrap(ncmpi_def_var(ncid, "z", SERGHEI_NC_REAL, 3, dimids, &zVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "SoilLoss", SERGHEI_NC_REAL, 3, dimids, &dzVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "bedEx", SERGHEI_NC_REAL, 3, dimids, &bedExVar), __LINE__);
#else
    dimids[0] = yDim;
    dimids[1] = xDim;
    ncwrap(ncmpi_def_var(ncid, "z", SERGHEI_NC_REAL, 2, dimids, &zVar), __LINE__, par.myrank);
#endif

    // define global attributes
    static char title[] = "SERGHEI simulation";
    ncwrap(ncmpi_put_att(ncid, NC_GLOBAL, "title", NC_CHAR, strlen(title) + 1, title), __LINE__, par.myrank);

    std::string source = std::string("SERGHEI ");
    ncwrap(ncmpi_put_att(ncid, NC_GLOBAL, "source", NC_CHAR, source.length() + 1, source.c_str()), __LINE__, par.myrank);

    auto nowtime = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(nowtime);
    source = std::string("Simulation started on ") + std::ctime(&now_time);
    ncwrap(ncmpi_put_att(ncid, NC_GLOBAL, "history", NC_CHAR, source.length() - 1, source.c_str()), __LINE__, par.myrank);

    longname.assign("projection_x_coordinate");
    ncwrap(ncmpi_put_att_text(ncid, xVar, "standard_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    longname.assign("x coordinate of projection");
    ncwrap(ncmpi_put_att_text(ncid, xVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncid, xVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    longname.assign("projection_y_coordinate");
    ncwrap(ncmpi_put_att_text(ncid, yVar, "standard_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    longname.assign("y coordinate of projection");
    ncwrap(ncmpi_put_att_text(ncid, yVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncid, yVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    longname.assign("Bed elevation");
    ncwrap(ncmpi_put_att_text(ncid, zVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncid, zVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    longname.assign("Water depth");
    ncwrap(ncmpi_put_att_text(ncid, hVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncid, hVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

#if SERGHEI_WRITE_HZ
    longname.assign("Water elevation");
    ncwrap(ncmpi_put_att_text(ncid, hzVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncid, hzVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);
#endif

    longname.assign("Velocity x-component");
    ncwrap(ncmpi_put_att_text(ncid, uVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    units.assign("m/s");
    ncwrap(ncmpi_put_att_text(ncid, uVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    longname.assign("Velocity y-component");
    ncwrap(ncmpi_put_att_text(ncid, vVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    units.assign("m/s");
    ncwrap(ncmpi_put_att_text(ncid, vVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    if (ss.inf.model)
    {
      longname.assign("Infiltration rate");
      ncwrap(ncmpi_put_att_text(ncid, infVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
      units.assign("m/s");
      ncwrap(ncmpi_put_att_text(ncid, infVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

      longname.assign("Accumulated infiltration volume");
      ncwrap(ncmpi_put_att_text(ncid, infVolVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
      units.assign("m^3");
      ncwrap(ncmpi_put_att_text(ncid, infVolVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);
    }

    // [CODE 1] 表面蒸发
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL

      longname.assign("Surface evaporation rate");
      ncwrap(ncmpi_put_att_text(ncid, evapVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
      units.assign("m/s");
      ncwrap(ncmpi_put_att_text(ncid, evapVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

#endif

    if (writeRain)
    {
      longname.assign("Rainfall rate");
      ncwrap(ncmpi_put_att_text(ncid, rainVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
      units.assign("mm/h");
      ncwrap(ncmpi_put_att_text(ncid, rainVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);
    }

#if SERGHEI_NETCDF_FORCING
    if (writeRainA)
    {
      longname.assign("Accumulated rainfall");
      ncwrap(ncmpi_put_att_text(ncid, rainAVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
      units.assign("m");
      ncwrap(ncmpi_put_att_text(ncid, rainAVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);
    }
#endif

#if SERGHEI_WRITE_SUBDOMS
    longname.assign("Subdomain index");
    ncwrap(ncmpi_put_att_text(ncid, subdomVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    units.assign("NA");
    ncwrap(ncmpi_put_att_text(ncid, subdomVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);
    static int range[] = {0, par.nranks - 1};
    ncwrap(ncmpi_put_att_int(ncid, subdomVar, "valid_range", NC_INT, 2, range), __LINE__, par.myrank);
#endif

    // End "define" mode
    ncwrap(ncmpi_enddef(ncid), __LINE__, par.myrank);

#if SERGHEI_DEBUG_OUTPUT
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " NetCDF define mode completed" << std::endl;
    ;
#endif

    // Compute x, y coordinates
    Kokkos::parallel_for("compute_grid_coord_x", dom.nx, KOKKOS_LAMBDA(int i) { xCoord(i) = dom.xll + (par.i_beg + i + 0.5) * dom.dxConst; });
    Kokkos::parallel_for("compute_grid_coord_x", dom.ny, KOKKOS_LAMBDA(int j) { yCoord(j) = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j + 0.5) * dom.dxConst; });
    Kokkos::fence();

    // Write out x, y coordinates
    st[0] = par.i_beg;
    ct[0] = dom.nx;
    ncwrap(ncmpi_put_vara_double_all(ncid, xVar, st, ct, xCoord.data()), __LINE__, par.myrank);
    st[0] = par.j_beg;
    ct[0] = dom.ny;
    ncwrap(ncmpi_put_vara_double_all(ncid, yVar, st, ct, yCoord.data()), __LINE__, par.myrank);

    st[0] = par.j_beg;
    st[1] = par.i_beg;
    ct[0] = dom.ny;
    ct[1] = dom.nx;

#if SERGHEI_NC_ENABLE_MISSING_VALUE
    real missing_value = dom.MISSING_VALUE;
    ncwrap(nc_put_att_float(ncid, zVar, SERGHEI_NC_MISSING_VALUE, SERGHEI_NC_REAL, 1, &missing_value));
#endif

// write fixed bed elevation
#if !SERGHEI_SEDIMENT_TRANSPORT
    writeNetCDFfield(dom, ncid, zVar, st, ct, state.isnodata, state.z, data);
#endif

#if SERGHEI_WRITE_SUBDOMS
    dimids[0] = yDim;
    dimids[1] = xDim;
    realArr subdom = realArr("subdom", dom.nCellMem);
    Kokkos::parallel_for("subdom", dom.nCellMem, KOKKOS_LAMBDA(int iGlob) { subdom(iGlob) = dom.id; });
    Kokkos::fence();
    writeNetCDFfield(dom, ncid, subdomVar, st, ct, state.isnodata, subdom, data);
#endif

    if constexpr (SERGHEI_MAPCLASS)
    {
      ushortArr usdata = ushortArr("usdata", dom.nCell);
      if (writeLandUse && state.landuse.use)
        writeNetCDFfield(dom, ncid, landuseVar, st, ct, state.isnodata, state.landuse.mapClass, usdata);
      if (writeSoilMap && ss.inf.soilmap.use)
        writeNetCDFfield(dom, ncid, soilmapVar, st, ct, state.isnodata, ss.inf.soilmap.mapClass, usdata);
    }

    if (writeRoughness)
      writeNetCDFfield(dom, ncid, roughnessVar, st, ct, state.isnodata, state.roughness, data);

    if (writeInfParameters)
    {
      if (ss.inf.model == INF_CONSTANT)
      {
        if (ss.inf.spatial == ss.inf.spatialNetCDF || ss.inf.spatial == ss.inf.spatialRaster)
        {
          if (ss.inf.spatial == ss.inf.spatialNetCDF)
          {
            real factor = Units::rainFactor("m/s", "mm/s");
            Kokkos::parallel_for("factorInfConstParam", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
                        int ii = dom.getIndex(iGlob);
                        ss.inf.capacity(ii) *= factor; });
            Kokkos::fence();
            writeNetCDFfield(dom, ncid, cinfCapVar, st, ct, state.isnodata, ss.inf.capacity, data);
            Kokkos::parallel_for("factorInfConstParam", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
                        int ii = dom.getIndex(iGlob);
                        ss.inf.capacity(ii) /= factor; });
          }
        }
        else
        {
          Kokkos::parallel_for("writeInfConstParam", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
                    int ii = dom.getIndex(iGlob);
                    int id = ss.inf.infLabel(ii);
                    data(iGlob) = ss.inf.constCap(id); });
          Kokkos::fence();
          writeNetCDFfield(dom, ncid, cinfCapVar, st, ct, state.isnodata, data, data);
        }
      }
    }

    writeStateNETCDF(state, dom, ss, par);

    ncwrap(ncmpi_close(ncid), __LINE__, par.myrank);
  }

  // NetCDF writer - requires the initialisation
  void outputNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir)
  {

    std::string filename;
    filename = dir + "output.nc";

    // Create the file
    ncwrap(ncmpi_open(MPI_COMM_WORLD, filename.c_str(), NC_WRITE, MPI_INFO_NULL, &ncid), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "h", &hVar), __LINE__, par.myrank);
#if SERGHEI_WRITE_HZ
    ncwrap(ncmpi_inq_varid(ncid, "h+z", &hzVar), __LINE__, par.myrank);
#endif
    ncwrap(ncmpi_inq_varid(ncid, "u", &uVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "v", &vVar), __LINE__, par.myrank);
    if (ss.inf.model)
    {
      ncwrap(ncmpi_inq_varid(ncid, "inf", &infVar), __LINE__, par.myrank);
      ncwrap(ncmpi_inq_varid(ncid, "infVol", &infVolVar), __LINE__, par.myrank);
    }

    //[CODE 1] 表面蒸发
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL

      ncwrap(ncmpi_inq_varid(ncid, "surfaceEvaporation", &evapVar), __LINE__, par.myrank);

#endif

    if (writeRain)
      ncwrap(ncmpi_inq_varid(ncid, "rain", &rainVar), __LINE__, par.myrank);
#if SERGHEI_NETCDF_FORCING
    if (writeRainA)
      ncwrap(ncmpi_inq_varid(ncid, "rainAccum", &rainAVar), __LINE__, par.myrank);
#endif
#if SERGHEI_RE_MODEL
    ncwrap(ncmpi_inq_varid(ncid, "qss", &qVar), __LINE__);
#endif

#if SERGHEI_EROSIVE_SHEAR
    ncwrap(ncmpi_inq_varid(ncid, "intTauDt", &intTauDtVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "phiAcc", &intTau32DtVar), __LINE__);
#endif
#if SERGHEI_SCALAR_TRANSPORT
    for (int iphi = 0; iphi < conf.st.nScalar; iphi++)
    {
      ncwrap(ncmpi_inq_varid(ncid, conf.st.scalarNames[iphi].c_str(), &phiVar[iphi]), __LINE__);
    }
#endif
#if SERGHEI_SEDIMENT_TRANSPORT
    ncwrap(ncmpi_inq_varid(ncid, "z", &zVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "SoilLoss", &dzVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "bedEx", &bedExVar), __LINE__);
#endif
#if SERGHEI_MAXFLOOD > 0
    ncwrap(ncmpi_inq_varid(ncid, "hMax", &hMaxVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "momMax", &momMaxVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "timehMax", &timehMaxVar), __LINE__);
#endif

    writeStateNETCDF(state, dom, ss, par);

    ncwrap(ncmpi_close(ncid), __LINE__, par.myrank);
  }

  // Builds the NetCDF dataset for the state
  void writeStateNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par)
  {
#if SERGHEI_DEBUG_OUTPUT
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
    ;
#endif
    ncArr data = ncArr("data", dom.nCell);

    MPI_Offset st[3], ct[3];
    double timeIter[numOut + 1];

    // write t. As the first one is written in the first iteration we should add +1
    for (int i = 0; i < numOut + 1; i++)
    {
      timeIter[i] = i * 1.0 * outFreq + dom.startTime;
    }

    st[0] = 0;
    ct[0] = numOut + 1;

    ncwrap(ncmpi_put_vara_double_all(ncid, tVar, st, ct, timeIter), __LINE__, par.myrank);

    st[0] = numOut;
    st[1] = par.j_beg;
    st[2] = par.i_beg;
    ct[0] = 1;
    ct[1] = dom.ny;
    ct[2] = dom.nx;

    // Write out depth
    writeNetCDFfield(dom, ncid, hVar, st, ct, state.isnodata, state.h, data);

#if SERGHEI_WRITE_HZ
    Kokkos::parallel_for("ncwrap_h+z", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
      int i, j;
      dom.unpackIndices(iGlob, j, i);
      int ii = (dom.hc + j) * (dom.nx + 2 * dom.hc) + dom.hc + i; // index with the extended domain (including halo cells)
      data(iGlob) = state.h(ii) + state.z(ii);
#if SERGHEI_NC_ENABLE_NAN
      if (state.isnodata(ii))
        data(iGlob) = SERGHEI_NAN;
#endif
    });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, hzVar, st, ct, data.data()), __LINE__, par.myrank);
#endif

    // [MERGE FIX] 消除多余重复的 qss 写出块，完美合并
#if SERGHEI_RE_MODEL
    Kokkos::parallel_for("ncwrap_qss", dom.ny * dom.nx, KOKKOS_LAMBDA(int iGlob) {
			int i,j;
			dom.unpackIndices(iGlob,j,i);
			int ii=(dom.hc+j)*(dom.nx+2*dom.hc) + dom.hc + i;
		    data(iGlob) = state.qss(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, qVar, st, ct, data.data()), __LINE__);
#endif

    // Write out x-velocity
    Kokkos::parallel_for("ncwrap_u", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
      int ii = dom.getIndex(iGlob);
      if (state.h(ii) > TOL12)
      {
        data(iGlob) = state.hu(ii) / state.h(ii);
      }
      else
      {
        data(iGlob) = 0.0;
      }
#if SERGHEI_NC_ENABLE_NAN
      if (state.isnodata(ii))
      {
        data(iGlob) = SERGHEI_NAN;
      }
#endif
    });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, uVar, st, ct, data.data()), __LINE__, par.myrank);

    // Write out y-velocity
    Kokkos::parallel_for("ncwrap_v", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
      int ii = dom.getIndex(iGlob);
      if (state.h(ii) > TOL12)
      {
        data(iGlob) = state.hv(ii) / state.h(ii);
      }
      else
      {
        data(iGlob) = 0.0;
      }
#if SERGHEI_NC_ENABLE_NAN
      if (state.isnodata(ii))
      {
        data(iGlob) = SERGHEI_NAN;
      }
#endif
    });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, vVar, st, ct, data.data()), __LINE__, par.myrank);

    if (ss.inf.model)
    {
      // infiltration rates
      writeNetCDFfield(dom, ncid, infVar, st, ct, state.isnodata, ss.inf.rate, data);
      // accumulated infiltration
      writeNetCDFfield(dom, ncid, infVolVar, st, ct, state.isnodata, ss.inf.infVol, data);
    }

    // [CODE 1] 表面蒸发输出
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL

      writeNetCDFfield(dom, ncid, evapVar, st, ct, state.isnodata, state.surfaceEvaporation, data);
    
#endif

    if (writeRain)
    {
      real rainFactor = Units::rainFactor("m/s", "mm/h");
      Kokkos::parallel_for("writeRain", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
					int ii = dom.getIndex(iGlob);
			ss.rainRate(ii) *= rainFactor; });
      writeNetCDFfield(dom, ncid, rainVar, st, ct, state.isnodata, ss.rainRate, data);
    }
#if SERGHEI_NETCDF_FORCING
    if (writeRainA)
      writeNetCDFfield(dom, ncid, rainAVar, st, ct, state.isnodata, ss.rainAccum, data);
#endif

#if SERGHEI_MAXFLOOD > 0
#if SERGHEI_MAXFLOOD == 1 // write only at the end, otherwise writes at every time step
    if (numOut == nOut)
    {
      st[0] = par.j_beg;
      st[1] = par.i_beg;
      ct[0] = dom.ny;
      ct[1] = dom.nx;
#endif
      writeNetCDFfield(dom, ncid, hMaxVar, st, ct, state.isnodata, state.hMax, data);
      writeNetCDFfield(dom, ncid, momMaxVar, st, ct, state.isnodata, state.momentumMax, data);
      writeNetCDFfield(dom, ncid, timehMaxVar, st, ct, state.isnodata, state.time_hMax, data);
#if SERGHEI_MAXFLOOD == 1
    }
#endif
#endif

// write scalars
#if SERGHEI_SCALAR_TRANSPORT
    if (conf.st.mode)
    {
      for (int iphi = 0; iphi < conf.st.nScalar; iphi++)
      {
        Kokkos::parallel_for("ncwrap_phi", dom.ny * dom.nx, KOKKOS_LAMBDA(int iGlob) {
          int ii = dom.getIndex(iGlob);
          if (state.h(ii) > TOL12)
          {
            data(iGlob) = state.ade.hphi(ii, iphi) / state.h(ii);
          }
          else
          {
            data(iGlob) = 0.0;
          }
#if SERGHEI_NC_ENABLE_NAN
          if (state.isnodata(ii))
          {
            data(iGlob) = SERGHEI_NAN;
          }
#endif
        });
        Kokkos::fence();
        ncwrap(ncmpi_put_vara_real_all(ncid, phiVar[iphi], st, ct, data.data()), __LINE__);
      }
    }
#endif

// write integrated Tau*Dt
#if SERGHEI_EROSIVE_SHEAR
    writeNetCDFfield(dom, ncid, intTauDtVar, st, ct, state.isnodata, state.shearAccum, data);
    writeNetCDFfield(dom, ncid, intTau32DtVar, st, ct, state.isnodata, state.phiTot, data);
#endif

// write bed changes
#if SERGHEI_SEDIMENT_TRANSPORT
    // bed elevation z
    Kokkos::parallel_for("ncwrap_z", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
      int ii = dom.getIndex(iGlob);
      data(iGlob) = state.z(ii);
#if SERGHEI_NC_ENABLE_NAN
      if (state.isnodata(ii))
      {
        data(iGlob) = SERGHEI_NAN;
      }
#endif
    });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, zVar, st, ct, data.data()), __LINE__);

    // erosion-deposition dz
    Kokkos::parallel_for("ncwrap_dz", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
      int ii = dom.getIndex(iGlob);
      data(iGlob) = state.z(ii) - state.zini(ii);
#if SERGHEI_NC_ENABLE_NAN
      if (state.isnodata(ii))
      {
        data(iGlob) = SERGHEI_NAN;
      }
#endif
    });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, dzVar, st, ct, data.data()), __LINE__);

    // entraintment
    Kokkos::parallel_for("ncwrap_bedEx", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
      int ii = dom.getIndex(iGlob);
      data(iGlob) = state.sediment.bedExchangeVol(ii);
#if SERGHEI_NC_ENABLE_NAN
      if (state.isnodata(ii))
      {
        data(iGlob) = SERGHEI_NAN;
      }
#endif
    });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, bedExVar, st, ct, data.data()), __LINE__);
#endif
  }

  // Error reporting routine for the PNetCDF I/O
  void ncwrap(int ierr, int line, int rank)
  {
    if (ierr != NC_NOERR)
    {
      std::cerr << RERROR "NetCDF reports error from rank " << rank << " at " << __FILE__ << ":" << line << std::endl
                << RERROR << ncmpi_strerror(ierr) << std::endl;
      exit(-1);
    }
  }
  void ncwrap(int ierr, int line)
  {
    if (ierr != NC_NOERR)
    {
      std::cerr << RERROR "NetCDF reports error at " << __FILE__ << ":" << line << std::endl
                << RERROR << ncmpi_strerror(ierr) << std::endl;
      exit(-1);
    }
  }

  // Write DAT file
  void outputDAT(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir)
  {

    std::string filename_h;
    std::string filename_u;
    std::string filename_v;
    std::string header_DAT;
    std::string el_h;
    std::string el_u;
    std::string el_v;
    filename_h = dir + "h_" + std::to_string(numOut) + ".input";
    filename_u = dir + "u_" + std::to_string(numOut) + ".input";
    filename_v = dir + "v_" + std::to_string(numOut) + ".input";

    int ncells = dom.nCell;

#if SERGHEI_EROSIVE_SHEAR
    std::string el_shear;
    std::string filename_shear;
    std::string el_shear32;
    std::string filename_shear32;
    filename_shear = dir + "shear_" + std::to_string(numOut) + ".input";
    filename_shear32 = dir + "shear32_" + std::to_string(numOut) + ".input";
#endif

    // Create a header for the raster
    header_DAT = "ncols          " + std::to_string(dom.nx) + "\n";
    header_DAT += "nrows          " + std::to_string(dom.ny) + "\n";
    header_DAT += "xllcorner      " + std::to_string(dom.xll) + "\n";
    header_DAT += "yllcorner      " + std::to_string(dom.yll) + "\n";
    header_DAT += "cellsize       " + std::to_string(dom.dx()) + "\n";
    header_DAT += "NODATA_value    -9999\n";

    // saving index
    int ii;
    int iGlob, i, j;

    // Open file and write header for h
    std::ofstream fOutStream_h(filename_h, std::ios::out);
    std::ofstream fOutStream_u(filename_u, std::ios::out);
    std::ofstream fOutStream_v(filename_v, std::ios::out);
#if SERGHEI_EROSIVE_SHEAR
    std::ofstream fOutStream_shear(filename_shear, std::ios::out);
    std::ofstream fOutStream_shear32(filename_shear, std::ios::out);
#endif

    // Write headers
    fOutStream_h << header_DAT;
    fOutStream_u << header_DAT;
    fOutStream_v << header_DAT;
#if SERGHEI_EROSIVE_SHEAR
    fOutStream_shear << header_DAT;
    fOutStream_shear32 << header_DAT;
#endif

    for (iGlob = 0; iGlob <= ncells; iGlob++)
    {
      unpackIndicesUniformGrid(iGlob, dom.ny, dom.nx, j, i);
      ii = (dom.hc + j) * (dom.nx + 2 * dom.hc) + dom.hc + i;
      el_h = std::to_string(state.h[ii]);
      el_u = std::to_string(state.hu[ii]);
      el_v = std::to_string(state.hv[ii]);
#if SERGHEI_EROSIVE_SHEAR
      el_shear = std::to_string(state.shearAccum[ii]);
      el_shear32 = std::to_string(state.phiTot[ii]);
#endif
      if (i % dom.nx == 0 && i > 0)
      {
        fOutStream_h << el_h << "\n";
        fOutStream_u << el_u << "\n";
        fOutStream_v << el_v << "\n";
#if SERGHEI_EROSIVE_SHEAR
        fOutStream_shear << el_shear << "\n";
        fOutStream_shear32 << el_shear32 << "\n";
#endif
      }
      else
      {
        fOutStream_h << el_h << " ";
        fOutStream_u << el_u << " ";
        fOutStream_v << el_v << " ";
#if SERGHEI_EROSIVE_SHEAR
        fOutStream_shear << el_shear << " ";
        fOutStream_shear32 << el_shear32 << " ";
#endif
      }
    }
    fOutStream_h.close();
    fOutStream_u.close();
    fOutStream_v.close();
#if SERGHEI_EROSIVE_SHEAR
    fOutStream_shear.close();
    fOutStream_shear32.close();
#endif
  }

#if SERGHEI_LPT == 1

  void outputParticle(Domain const &dom, ParticleTracker const &particles, std::string dir)
  {

    std::string filename;
    filename = dir + "particles" + std::to_string(numOut) + ".vtk";

    std::ofstream particleFile(filename);

    if (particleFile.is_open())
    {
      particleFile << "# vtk DataFile Version 2.0\n";
      particleFile << "Particles\n";
      particleFile << "ASCII\n";
      particleFile << "DATASET UNSTRUCTURED_GRID\n";

      int num_points;
      int noactive = 0;
      num_points = particles.N_par;
      for (int i = 0; i < num_points; i++)
      {
        if (particles.particles[i].active == 0)
        {
          noactive += 1;
        }
      }
      particleFile << "POINTS " << num_points - noactive << " float" << std::endl;
      for (int i = 0; i < num_points; i++)
      {
        if (particles.particles[i].active == 1)
        {
          particleFile << std::setprecision(9) << particles.particles[i].x(0) << " " << particles.particles[i].x(1) << " " << particles.particles[i].z << std::endl;
        }
      }

      int num_cells = num_points - noactive;
      particleFile << "CELLS " << num_cells << " " << 2 * num_cells << std::endl;
      for (int i = 0; i < num_cells; i++)
      {
        particleFile << "1 " << i << std::endl;
      }

      particleFile << "CELL_TYPES " << num_cells << std::endl;
      for (int i = 0; i < num_cells; i++)
      {
        particleFile << "1" << std::endl;
      }
      particleFile.close();
    }
  }

  void outputIniParticle(Domain const &dom, ParticleTracker const &particles, std::string dir)
  {

    std::string filename;
    int numOutIni;
    numOutIni = 0;
    filename = dir + "particles" + std::to_string(numOutIni) + ".vtk";

    std::ofstream particleFile(filename);

    if (particleFile.is_open())
    {
      particleFile << "# vtk DataFile Version 2.0\n";
      particleFile << "Particles\n";
      particleFile << "ASCII\n";
      particleFile << "DATASET UNSTRUCTURED_GRID\n";

      int num_points;
      num_points = particles.N_par;
      int noactive = 0;
      num_points = particles.N_par;
      for (int i = 0; i < num_points; i++)
      {
        if (particles.particles[i].active == 0)
        {
          noactive += 1;
        }
      }
      particleFile << "POINTS " << num_points - noactive << " float" << std::endl;
      for (int i = 0; i < num_points; i++)
      {
        if (particles.particles[i].active == 1)
        {
          particleFile << std::setprecision(9) << particles.particles[i].x(0) << " " << particles.particles[i].x(1) << " " << particles.particles[i].z << std::endl;
        }
      }

      int num_cells = num_points;
      particleFile << "CELLS " << num_cells - noactive << " " << 2 * (num_cells - noactive) << std::endl;
      for (int i = 0; i < num_cells - noactive; i++)
      {
        particleFile << "1 " << i << std::endl;
      }

      particleFile << "CELL_TYPES " << num_cells - noactive << std::endl;
      for (int i = 0; i < num_cells - noactive; i++)
      {
        particleFile << "1" << std::endl;
      }
      particleFile.close();
    }
  }

  void writeParticleFile(Domain const &dom, ParticleTracker const &particles, Parallel const &par, std::string dir)
  {
    std::string filename = dir + "particles.out";
    particlesFile.open(filename);
    particlesFile << "ID Particle \t CoveredDistance \t TotalTime \t TravelTime \t StopTime" << std::endl;
    int num_points;
    num_points = particles.N_par;
    for (int i = 0; i < num_points; i++)
    {
      particlesFile << std::setprecision(6) << i << " " << particles.particles[i].CoveredDistance << " " << particles.particles[i].TotalLife << " " << particles.particles[i].TravelTime << " " << particles.particles[i].StopTime << std::endl;
    }
    particlesFile.close();
  }

  void writeParticle_TemporalFile(Domain const &dom, ParticleTracker const &particles, std::string dir)
  {
    std::string filename = dir + "particles_temporal.out";
    particlesFile.open(filename, std::ios::app);
    if (particlesFile.is_open())
    {
      int num_points;
      real total_distance;
      num_points = particles.N_par;
      total_distance = 0.0;
      for (int i = 0; i < num_points; i++)
      {
        total_distance += particles.particles[i].CoveredDistance;
      }
      particlesFile << dom.etime << " " << total_distance << " " << std::endl;
      particlesFile.close();
    }
  }

  void outputInitParticlesNETCDF(ParticleTracker const &particles, Domain const &dom, Parallel const &par, std::string dir)
  {

    int dimids[2];
    MPI_Offset st[3], ct[3];
    ncArr data = ncArr("data", particles.N_par);
    static char timeUnits[] = "seconds";
    std::string filename;
    std::string longname;
    std::string units;

    lpt_filename = dir + "lpt.nc";

    // Create the file
    ncwrap(ncmpi_create(MPI_COMM_WORLD, lpt_filename.c_str(), SERGHEI_NC_MODE, MPI_INFO_NULL, &ncidP), __LINE__);

    ncwrap(ncmpi_def_dim(ncidP, "time", (MPI_Offset)NC_UNLIMITED, &tpDim), __LINE__);
    ncwrap(ncmpi_def_dim(ncidP, "particle", (MPI_Offset)particles.N_par, &pDim), __LINE__);
    dimids[0] = tpDim;
    ncwrap(ncmpi_def_var(ncidP, "time", NC_DOUBLE, 1, dimids, &tpVar), __LINE__);
    ncwrap(ncmpi_put_att_text(ncidP, tpVar, "units", strlen(timeUnits), timeUnits), __LINE__);
    dimids[0] = pDim;
    ncwrap(ncmpi_def_var(ncidP, "particle", NC_INT, 1, dimids, &pVar), __LINE__);
    static char dimlessUnits[] = "1";
    ncwrap(ncmpi_put_att_text(ncidP, pVar, "units", strlen(dimlessUnits), dimlessUnits), __LINE__);
    dimids[0] = tpDim;
    dimids[1] = pDim;
    ncwrap(ncmpi_def_var(ncidP, "x", SERGHEI_NC_REAL, 2, dimids, &xpVar), __LINE__);
    ncwrap(ncmpi_def_var(ncidP, "y", SERGHEI_NC_REAL, 2, dimids, &ypVar), __LINE__);
    ncwrap(ncmpi_def_var(ncidP, "z", SERGHEI_NC_REAL, 2, dimids, &zpVar), __LINE__);

    // define global attributes
    static char title[] = "SERGHEI-LPT simulation";
    ncwrap(ncmpi_put_att(ncidP, NC_GLOBAL, "title", NC_CHAR, strlen(title) + 1, title), __LINE__);
    static char conventions[] = "CF-1.6";
    ncwrap(ncmpi_put_att(ncidP, NC_GLOBAL, "Conventions", NC_CHAR, strlen(conventions), conventions), __LINE__);

    longname.assign("x-coordinate");
    ncwrap(ncmpi_put_att_text(ncidP, xpVar, "long_name", longname.length(), longname.c_str()), __LINE__);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncidP, xpVar, "units", units.length(), units.c_str()), __LINE__);

    longname.assign("y-coordinate");
    ncwrap(ncmpi_put_att_text(ncidP, ypVar, "long_name", longname.length(), longname.c_str()), __LINE__);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncidP, ypVar, "units", units.length(), units.c_str()), __LINE__);

    longname.assign("z-coordinate");
    ncwrap(ncmpi_put_att_text(ncidP, zpVar, "long_name", longname.length(), longname.c_str()), __LINE__);
    units.assign("m");
    ncwrap(ncmpi_put_att_text(ncidP, zpVar, "units", units.length(), units.c_str()), __LINE__);

    // End "define" mode
    ncwrap(ncmpi_enddef(ncidP), __LINE__);

#if SERGHEI_DEBUG_OUTPUT
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " NetCDF define mode completed" << std::endl;
    ;
#endif

    st[0] = 0;
    ct[0] = particles.N_par;

    intArr pID = intArr("particleID", particles.N_par);

    Kokkos::parallel_for("generate_particle_id", particles.N_par, KOKKOS_LAMBDA(int ii) { pID(ii) = ii; });
    ncwrap(ncmpi_put_vara_int_all(ncidP, pVar, st, ct, pID.data()), __LINE__);

    int savedNumOut = numOut;
    numOut = 0;
    writeStateNETCDF_particles(particles, dom, par);
    numOut = savedNumOut;

    ncwrap(ncmpi_close(ncidP), __LINE__);
  }

  void writeStateNETCDF_particles(ParticleTracker const &particles, Domain const &dom, Parallel const &par)
  {
#if SERGHEI_DEBUG_OUTPUT
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
    ;
#endif

    MPI_Offset st[2], ct[2];
    double timeIter[numOut + 1];

    for (int i = 0; i < numOut + 1; i++)
    {
      timeIter[i] = i * 1.0 * outFreq + dom.startTime;
    }
    st[0] = 0;
    ct[0] = numOut + 1;

    ncwrap(ncmpi_put_vara_double_all(ncidP, tpVar, st, ct, timeIter), __LINE__);

    st[0] = numOut;
    st[1] = 0;
    ct[0] = 1;
    ct[1] = particles.N_par;

    int inactive = 0;
    Kokkos::parallel_reduce("particle_count_active", particles.N_par, KOKKOS_LAMBDA(int ii, int inactive_l) {
			if(particles.particles(ii).active!=1) inactive_l++; }, Kokkos::Sum<int>(inactive));

    int Np;
    Np = particles.N_par;
    ncArr data = ncArr("data", Np);
    const auto nan_val = static_cast<ncreal>(SERGHEI_NAN);

    // Write out x-coordinate
    Kokkos::parallel_for("ncwrap_x_par", Np, KOKKOS_LAMBDA(int ii) { data(ii) = (particles.particles(ii).active == 1) ? particles.particles(ii).x(0) : nan_val; });
    Kokkos::fence();

    ncwrap(ncmpi_put_vara_real_all(ncidP, xpVar, st, ct, data.data()), __LINE__);

    // Write out y-coordinate
    Kokkos::parallel_for("ncwrap_y_par", Np, KOKKOS_LAMBDA(int ii) { data(ii) = (particles.particles(ii).active == 1) ? particles.particles(ii).x(1) : nan_val; });
    Kokkos::fence();

    ncwrap(ncmpi_put_vara_real_all(ncidP, ypVar, st, ct, data.data()), __LINE__);

    // Write out z-coordinate
    Kokkos::parallel_for("ncwrap_z_par", Np, KOKKOS_LAMBDA(int ii) { data(ii) = (particles.particles(ii).active == 1) ? particles.particles(ii).z : nan_val; });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncidP, zpVar, st, ct, data.data()), __LINE__);
  }

  void outputNETCDF_particles(ParticleTracker const &particles, Domain const &dom, Parallel const &par, std::string dir)
  {
    lpt_filename = dir + "lpt.nc";
    ncwrap(ncmpi_open(MPI_COMM_WORLD, lpt_filename.c_str(), NC_WRITE, MPI_INFO_NULL, &ncidP), __LINE__);
    ncwrap(ncmpi_inq_varid(ncidP, "time", &tpVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncidP, "x", &xpVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncidP, "y", &ypVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncidP, "z", &zpVar), __LINE__);
    writeStateNETCDF_particles(particles, dom, par);
    ncwrap(ncmpi_close(ncidP), __LINE__);
  }
#endif

  // Writes VTK file
  void outputVTK(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir)
  {

    std::string filename;
    filename = dir + "result" + std::to_string(numOut) + ".vtk";
    real *xCoord;
    real *yCoord;
    real *data_cpu;

    int ncells = dom.nCell;

    int nVars = 4; // z, h, (u,v)
#if SERGHEI_WRITE_HZ
    nVars++;
#endif
#if SERGHEI_EROSIVE_SHEAR
    nVars++; // intTauDt
#endif
#if SERGHEI_SEDIMENT_TRANSPORT
    nVars++; // dz
#endif
#if SERGHEI_DEBUG_BOUNDARY
    nVars++;
#endif
#if SERGHEI_LPT
    nVars++;
#endif
#if SERGHEI_VERTICAL_VELOCITY
    nVars++;
#endif
    if (ss.inf.model)
      nVars = nVars + 2; // 2 more variables: inf, infVol

    realArr data = realArr("data", nVars * ncells);
#ifdef KOKKOS_ENABLE_CUDA
    cudaMallocHost(&data_cpu, nVars * ncells * sizeof(real));
#else
    data_cpu = data.data();
#endif

    xCoord = (real *)malloc((dom.nx + 1) * sizeof(real));
    yCoord = (real *)malloc((dom.ny + 1) * sizeof(real));

#if SERGHEI_DEBUG_OUTPUT
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " VTK data allocated" << std::endl;
    ;
#endif

    // Compute x, y coordinates
    for (int i = 0; i < dom.nx + 1; i++)
    {
      xCoord[i] = dom.xll + (par.i_beg + i) * dom.dxConst;
    };

    for (int j = 0; j < dom.ny + 1; j++)
    {
      yCoord[j] = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j) * dom.dxConst;
    };

#if SERGHEI_DEBUG_OUTPUT
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " VTK coordinate data generated" << std::endl;
    ;
#endif
    int nnodes = (dom.nx + 1) * (dom.ny + 1);
    int of_sw = 3;
    int offset_z = 0;
    int offset_h = ncells;
    int offset_hu = 2 * ncells;
    int offset_hv = 3 * ncells;
#if SERGHEI_EROSIVE_SHEAR
    int of_sa = of_sw;
    of_sa++;
    of_sw = of_sa;
    int offset_sa = of_sa * ncells;
#endif
#if SERGHEI_SEDIMENT_TRANSPORT
    int of_dz = of_sw;
    of_dz++;
    of_sw = of_dz;
    int offset_dz = of_dz * ncells;
#endif
#if SERGHEI_LPT && SERGHEI_VERTICAL_VELOCITY == 0
    int offset_particles = 4 * ncells;
#elif SERGHEI_LPT && SERGHEI_VERTICAL_VELOCITY
    int offset_particles = 4 * ncells;
    int offset_vertical_velocity = 5 * ncells;
#elif SERGHEI_LPT == 0 && SERGHEI_VERTICAL_VELOCITY == 1
  int offset_vertical_velocity = 4 * ncells;
#endif
    int of_bc = of_sw;
#if SERGHEI_DEBUG_BOUNDARY
    of_bc++;
    int offset_bc = of_bc * ncells;
#endif
    int of_inf = of_bc + 1;
    int offset_infrate = of_inf * ncells;
    int offset_infVol = (of_inf + 1) * ncells;

    Kokkos::parallel_for("vtkwrap_all", ncells, KOKKOS_LAMBDA(int iGlob) {
      int ii = dom.getIndex(iGlob);
      data(iGlob + offset_z) = state.z(ii);
      data(iGlob + offset_h) = state.h(ii);
      data(iGlob + offset_hu) = state.hu(ii);
      data(iGlob + offset_hv) = state.hv(ii);
      if (ss.inf.model)
      {
        data(iGlob + offset_infrate) = ss.inf.rate(ii);
        data(iGlob + offset_infVol) = ss.inf.infVol(ii);
      }
#if SERGHEI_LPT && SERGHEI_VERTICAL_VELOCITY
      if (dom.etime == 0.0)
      {
        data(iGlob + offset_particles) = state.tparticles(ii);
      }
      else
      {
        data(iGlob + offset_particles) = state.tparticles(ii) / (dom.etime);
      }
      data(iGlob + offset_vertical_velocity) = state.w(ii);
#elif SERGHEI_LPT == 0 && SERGHEI_VERTICAL_VELOCITY == 1
      data(iGlob + offset_vertical_velocity) = state.w(ii);
#endif
#if SERGHEI_DEBUG_BOUNDARY
      data(iGlob + offset_bc) = state.isBound(ii);
#endif
    });
    Kokkos::fence();

#ifdef KOKKOS_ENABLE_CUDA
    cudaMemcpyAsync(data_cpu, data.data(), nVars * ncells * sizeof(real), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
#endif

    int i, j, iGlob;

    std::ofstream fOutStream(filename);
    if (fOutStream.is_open())
    {
      fOutStream << "# vtk DataFile Version 3.0.\nOutputfile\nASCII\nDATASET UNSTRUCTURED_GRID\n";

      fOutStream << "POINTS " << nnodes << SERGHEI_VTK_REAL << "\n";
      for (j = 0; j <= dom.ny; j++)
      {
        for (i = 0; i <= dom.nx; i++)
        {
          fOutStream << std::setprecision(9) << xCoord[i] << " " << yCoord[j] << " 0.0\n";
        }
      }
      fOutStream << "CELLS " << ncells << " " << 5 * ncells << "\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        int ii, jj;
        dom.unpackIndices(iGlob, jj, ii);
        fOutStream << "4 " << jj * (dom.nx + 1) + ii << " " << (jj + 1) * (dom.nx + 1) + ii << " " << (jj + 1) * (dom.nx + 1) + ii + 1 << " " << jj * (dom.nx + 1) + ii + 1 << "\n";
      }

      fOutStream << "CELL_TYPES " << ncells << "\n";
      for (i = 0; i < ncells; i++)
      {
        fOutStream << "9\n";
      }

      fOutStream << "CELL_DATA " << ncells << "\n";
      fOutStream << "SCALARS z " << SERGHEI_VTK_REAL << "\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob] << "\n";
      }

      fOutStream << "SCALARS h " << SERGHEI_VTK_REAL << "\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob + offset_h] << "\n";
      }

#if SERGHEI_WRITE_HZ
      fOutStream << "SCALARS h+z " << SERGHEI_VTK_REAL << "\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob] + data_cpu[iGlob + offset_h] << "\n";
      }
#endif

      fOutStream << "VECTORS velocity " << SERGHEI_VTK_REAL << "\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        if (data_cpu[iGlob + ncells] > TOL12)
        {
          fOutStream << std::setprecision(9) << data_cpu[iGlob + offset_hu] / data_cpu[iGlob + offset_h] << " " << data_cpu[iGlob + offset_hv] / data_cpu[iGlob + offset_h] << " 0.0\n";
        }
        else
        {
          fOutStream << "0.0 0.0 0.0\n";
        }
      }

      if (ss.inf.model)
      {
        fOutStream << "SCALARS inf " << SERGHEI_VTK_REAL << "\n";
        fOutStream << "LOOKUP_TABLE default\n";
        for (iGlob = 0; iGlob < ncells; iGlob++)
        {
          fOutStream << std::setprecision(9) << data_cpu[iGlob + offset_infrate] << "\n";
        }

        fOutStream << "SCALARS infVol " << SERGHEI_VTK_REAL << "\n";
        fOutStream << "LOOKUP_TABLE default\n";
        for (iGlob = 0; iGlob < ncells; iGlob++)
        {
          fOutStream << std::setprecision(9) << data_cpu[iGlob + offset_infVol] << "\n";
        }
      }
#if SERGHEI_LPT
      fOutStream << "SCALARS t_particles " << SERGHEI_VTK_REAL << "\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob + offset_particles] << "\n";
      }
#endif
#if SERGHEI_VERTICAL_VELOCITY
      fOutStream << "SCALARS w_velocity " << SERGHEI_VTK_REAL << "\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob + offset_vertical_velocity] << "\n";
      }
#endif

#if SERGHEI_EROSIVE_SHEAR
      fOutStream << "SCALARS shearAccum double\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob + offset_sa] << "\n";
      }
#endif

#if SERGHEI_SEDIMENT_TRANSPORT
      fOutStream << "SCALARS SoilLoss double\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob + offset_dz] << "\n";
      }
#endif

#if SERGHEI_DEBUG_BOUNDARY
      fOutStream << "SCALARS bc int\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << data_cpu[iGlob + offset_bc] << "\n";
      }
#endif

#if SERGHEI_SCALAR_TRANSPORT
      if (conf.st.mode)
      {
        double phi = 0;
        for (int iphi = 0; iphi < conf.st.nScalar; iphi++)
        {
          fOutStream << "SCALARS " << conf.st.scalarNames[iphi] << " double\n";
          fOutStream << "LOOKUP_TABLE default\n";
          for (iGlob = 0; iGlob < ncells; iGlob++)
          {
            int ii = dom.getIndex(iGlob);
            if (state.h(ii) > TOL12)
            {
              phi = state.ade.hphi(ii, iphi) / state.h(ii);
            }
            else
            {
              phi = 0.;
            }
            fOutStream << std::setprecision(9) << phi << "\n";
          }
        }
      }
#endif
      fOutStream.close();
    }
  }

  // Writes a binary raster file - initialisation
  void initBIN(const State &state, Domain const &dom, Parallel const &par, std::string dir)
  {

    std::string filename;
    realArr data = realArr("data", dom.nCell);

    // header
    filename = dir + std::to_string(par.myrank) + "result.hdr";
    std::ofstream fOutStream1(filename);
    if (fOutStream1.is_open())
    {
      fOutStream1 << "ncols " << dom.nx << std::endl;
      fOutStream1 << "nrows " << dom.ny << std::endl;
      fOutStream1 << "xllcorner " << dom.xll + par.i_beg * dom.dxConst << std::endl;
      fOutStream1 << "yllcorner " << dom.yll + dom.ny_glob * dom.dxConst - par.j_beg * dom.dxConst << std::endl;
      fOutStream1 << "cellsize " << dom.dxConst << std::endl;
      fOutStream1 << "nodata_value -9999" << std::endl;
      fOutStream1 << "byteorder msbfirst" << std::endl;
      fOutStream1.close();
    }

    Kokkos::fence();

    // z
    Kokkos::parallel_for("binwrap_init_z", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
				int ii = dom.getIndex(iGlob);
      	data(iGlob) = state.z(ii); });

    Kokkos::fence();

    filename = dir + std::to_string(par.myrank) + "elevation.bin";
    std::ofstream fOutStream2(filename.c_str(), std::ios::binary);
    fOutStream2.write((char *)&data[0], dom.nCell * sizeof(real));
    fOutStream2.close();

    Kokkos::fence();

    // h
    Kokkos::parallel_for("binwrap_init_h", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
				int ii = dom.getIndex(iGlob);
      	data(iGlob) = state.h(ii); });
    Kokkos::fence();

    filename = dir + std::to_string(par.myrank) + "result" + std::to_string(numOut) + ".bin";
    std::ofstream fOutStream3(filename.c_str(), std::ios::binary);
    fOutStream3.write((char *)&data[0], dom.nCell * sizeof(real));
    fOutStream3.close();

    Kokkos::fence();

    if (par.nranks > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // Writes a binary raster file
  void outputBIN(const State &state, Domain const &dom, Parallel const &par, std::string dir)
  {

    std::string filename;
    realArr data = realArr("data", dom.nCell);

    Kokkos::parallel_for("binwrap_h", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
            int ii = dom.getIndex(iGlob);
            data(iGlob) = state.h(ii); });

    Kokkos::fence();

    filename = dir + std::to_string(par.myrank) + "result" + std::to_string(numOut) + ".bin";
    std::ofstream fOutStream(filename.c_str(), std::ios::binary);
    fOutStream.write((char *)&data[0], dom.nCell * sizeof(real));
    fOutStream.close();

    if (par.nranks > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // Write header and initial state for time series files
  int writeTimeSeriesIni(const State &state, Domain const &dom, Parallel const &par,
                         SourceSinkData &ss, surfaceIntegrator &sint, boundaryIntegrator &bint,
                         std::vector<ExtBC> &extbc, std::string dir)
  {

    numObs = 0;
    std::string filename = dir + "domainTimeSeries.out";
    domainOutputFile.open(filename);

    if (domainOutputFile.is_open())
    {

      // Write the header
      domainOutputFile << "Time ";
      domainOutputFile << "SurfaceVolume ";

      int _ncellsBC = 0;
      bint.integrate(extbc, dom, 1);
      _ncellsBC = bint.ncellsBC;

      if (_ncellsBC)
      {
        // global inflow/outflow
        domainOutputFile << "BoundaryInflow ";
        domainOutputFile << "BoundaryInflowAccum ";
        domainOutputFile << "BoundaryOutflow ";
        domainOutputFile << "BoundaryOutflowAccum ";
        // boundary flow per BC
        for (int i = 0; i < extbc.size(); i++)
        {
          domainOutputFile << "BoundaryFlow_" << i << " ";
          domainOutputFile << "BoundaryAccum_" << i << " ";
        }
      }

      if (dom.isRain)
      {
        domainOutputFile << "RainFlux ";
        domainOutputFile << "RainAccum ";
      }

      //[CODE 1] 增加表面蒸发
 
        domainOutputFile << "SurfaceEvapFlux ";
        domainOutputFile << "SurfaceEvapAccum ";
      

      if (ss.inf.model)
      {
        domainOutputFile << "InfFlux ";
        domainOutputFile << "InfAccum ";
      }

#if SERGHEI_SUSPENDED_SEDIMENT
      domainOutputFile << "SurfaceSolidVolume ";
#endif

      domainOutputFile << std::endl;
    }
    else
    {
      std::cerr << RERROR "Could not create domainTimeSeries.out file" << std::endl;
      return 0;
    }

    // write the data
    writeTimeSeries(state, dom, par, sint, bint, *bint.extbc);

    numObs++;

    return 1;
  }

  // Writes time series files
  void writeTimeSeries(const State &state, Domain const &dom, Parallel const &par, surfaceIntegrator &sint, boundaryIntegrator &bint, std::vector<ExtBC> &extbc)
  {
    timer.reset();
    for (int i = 0; i < extbc.size(); i++)
      extbc[i].reduce(par);
    if (par.masterproc)
      writeDomainTimeSeries(state, dom, par, sint, bint, extbc);
    numObs++;
    dom.timers.swe.io.out += timer.seconds();
  }

  void writeDomainTimeSeries(State const &state, Domain const &dom,
                             Parallel const &par, surfaceIntegrator &sint,
                             boundaryIntegrator &bint, const std::vector<ExtBC> &extbc)
  {

    // Write the data
    std::cout.precision(OUTPUT_PRECISION);
    domainOutputFile << std::scientific << dom.etime << " ";
    domainOutputFile << std::scientific << sint.surfaceVolumeG << " ";

    if (bint.ncellsBC)
    {
      domainOutputFile << std::scientific << bint.inflowDischargeG << " ";
      domainOutputFile << std::scientific << bint.inflowAccumulatedG << " ";
      domainOutputFile << std::scientific << bint.outflowDischargeG << " ";
      domainOutputFile << std::scientific << bint.outflowAccumulatedG << " ";
      for (int i = 0; i < extbc.size(); i++)
        domainOutputFile << extbc[i].netQ << " ";
      for (int i = 0; i < extbc.size(); i++)
        domainOutputFile << extbc[i].netVol << " ";
    }

    if (dom.isRain)
    {
      domainOutputFile << std::scientific << sint.rainFluxG << " ";
      domainOutputFile << std::scientific << sint.rainAccumG << " ";
    }

    //[CODE 1]
 
      domainOutputFile << std::scientific << sint.evapFluxG << " ";
      domainOutputFile << std::scientific << sint.evapAccumG << " ";
   

    if (sint.ss->inf.model)
    {
      domainOutputFile << std::scientific << sint.infFluxG << " ";
      domainOutputFile << std::scientific << sint.infAccumG << " ";
    }

#if SERGHEI_SUSPENDED_SEDIMENT
    domainOutputFile << std::scientific << sint.surfaceSolidVolumeG << " ";
#endif

    domainOutputFile << std::endl;
  }

  void closeOutputStreams()
  {
    domainOutputFile.close();
  }

  void writeTimerRank(const Parallel &par, const real &time, const std::string &name)
  {
    real *logdata;
    logdata = (real *)malloc(par.nranks * sizeof(real));

    MPI_Gather(&time, 1, MPI_DOUBLE, logdata, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (par.masterproc)
      logFile << name << "\t\t";
    for (int ii = 0; ii < par.nranks; ii++)
      logFile << " :\t" << logdata[ii];
    logFile << std::endl;
    real avg = 0;
    for (int ii = 0; ii < par.nranks; ii++)
      avg += logdata[ii];
    avg = avg / par.nranks;
    if (par.masterproc)
      logFile << name << "_lb\t\t";
    for (int ii = 0; ii < par.nranks; ii++)
      logFile << " :\t" << logdata[ii] / avg;
    logFile << std::endl;
  }

  void writeLogFile(Domain const &dom, Parallel const &par, std::string dir)
  {
    if (par.masterproc)
    {
      std::string filename = dir + "log.out";
      logFile.open(filename);
      if (logFile.is_open())
      {
        logFile << "DomainArea [m2]: " << dom.areaGlobal << std::endl;
        logFile << "nCell(SW) : " << dom.nCellGlobal << std::endl;
        logFile << "nCellValid(SW) : " << dom.nCellValidGlobal << std::endl;

        logFile << "-----SERGHEI TIMERS-----" << std::endl;
        dom.timers.report(logFile, par.nranks);
        logFile << std::endl;
        logFile << "-----RELATIVE TIMERS-----" << std::endl;
        dom.relative.report(logFile, par.nranks);
        logFile << std::endl;
        logFile << "-----LOAD BALANCE-----" << std::endl;
        dom.timers.reportLoadBalance(logFile, par.nranks);
        logFile << "-------------------------------- " << std::endl;
        auto nowtime = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(nowtime);
        logFile << std::endl
                << "DateTime : " << std::ctime(&now_time) << std::endl;
        char hostbuffer[256];
        int hostname = gethostname(hostbuffer, sizeof(hostbuffer));
        logFile << "Machine : " << hostbuffer << std::endl;
#if defined(KOKKOS_ENABLE_CUDA)
        cudaDeviceProp deviceProp;
        cudaError_t result = cudaGetDeviceProperties(&deviceProp, 0);
        logFile << "GPU : " << deviceProp.name << std::endl;
#elif defined(KOKKOS_ENABLE_HIP)
        hipDeviceProp_t deviceProp;
        hipError_t result = hipGetDeviceProperties(&deviceProp, 0);
        logFile << "GPU : " << deviceProp.name << std::endl;
#elif defined(KOKKOS_ENABLE_SYCL)
      sycl::queue q{sycl::gpu_selector_v};
      auto device = q.get_device();
      logFile << "GPU : " << device.get_info<sycl::info::device::name>() << std::endl;
#else
#if !(defined(__arm__) || defined(__aarch64__))
      std::string CPUBrandString;
      CPUBrandString.resize(49);
      uint *CPUInfo = reinterpret_cast<uint *>(CPUBrandString.data());
      for (uint i = 0; i < 3; i++)
        __cpuid(0x80000002 + i, CPUInfo[i * 4 + 0], CPUInfo[i * 4 + 1], CPUInfo[i * 4 + 2], CPUInfo[i * 4 + 3]);
      CPUBrandString.assign(CPUBrandString.data()); // correct null terminator
      logFile << "CPU: " << CPUBrandString << std::endl;
#endif
#endif
        logFile << "nTasks : " << par.nranks << std::endl;

        // write compilation setup
        logFile << "\nSERGHEI_GIT_VERSION : " << SERGHEI_GIT_VERSION << std::endl;
        logFile << "\n-------------------------\nMODEL COMPONENT SETUP" << std::endl;
        logFile << "SERGHEI_TOOLS : " << SERGHEI_TOOLS << std::endl;
        logFile << "SERGHEI_RE_MODEL : " << SERGHEI_RE_MODEL << std::endl;
        logFile << "SERGHEI_VEGETATION_MODEL : " << SERGHEI_VEGETATION_MODEL << std::endl;
        logFile << "SERGHEI_SCALAR_TRANSPORT : " << SERGHEI_SCALAR_TRANSPORT << std::endl;
        logFile << "SERGHEI_SCALAR_DIFFUSION : " << SERGHEI_SCALAR_DIFFUSION << std::endl;
        logFile << "SERGHEI_SEDIMENT_TRANSPORT : " << SERGHEI_SEDIMENT_TRANSPORT << std::endl;
        logFile << "SERGHEI_SUSPENDED_SEDIMENT : " << SERGHEI_SUSPENDED_SEDIMENT << std::endl;
        // [CODE 1] Added components
#ifdef CROP_GROWTH_MODEL
        logFile << "CROP_GROWTH_MODEL : " << CROP_GROWTH_MODEL << std::endl;
#endif
#ifdef SW_GW_EVAPORATION_TRANSPIRATION_MODEL
        logFile << "SW_GW_EVAPORATION_TRANSPIRATION_MODEL : " << SW_GW_EVAPORATION_TRANSPIRATION_MODEL << std::endl;
#endif
#ifdef SERGHEI_SUBSURFACE_TRANSPORT
        logFile << "SERGHEI_SUBSURFACE_TRANSPORT : " << SERGHEI_SUBSURFACE_TRANSPORT << std::endl;
#endif
#ifdef SERGHEI_SURFACE_TRANSPORT
        logFile << "SERGHEI_SURFACE_TRANSPORT : " << SERGHEI_SURFACE_TRANSPORT << std::endl;
#endif
        logFile << "SERGHEI_REBALANCE_SOLVER_CONTRIBUTIONS : " << SERGHEI_REBALANCE_SOLVER_CONTRIBUTIONS << std::endl;
        logFile << "SERGHEI_UPWIND_BED : " << SERGHEI_UPWIND_BED << std::endl;
        logFile << "SERGHEI_BEDLOAD_SEDIMENT : " << SERGHEI_BEDLOAD_SEDIMENT << std::endl;
        logFile << "SERGHEI_PARTICLE_TRACKING : " << SERGHEI_LPT << std::endl;
        logFile << "SERGHEI_FRICTION_MODEL : " << SERGHEI_FRICTION_MODEL << std::endl;
        logFile << "SERGHEI_POINTWISE_FRICTION : " << SERGHEI_POINTWISE_FRICTION << std::endl;
        logFile << "SERGHEI_EROSIVE_SHEAR : " << SERGHEI_EROSIVE_SHEAR << std::endl;
        logFile << "SERGHEI_RAINFALL_POLYGONS : " << SERGHEI_RAINFALL_POLYGONS << std::endl;
        logFile << "SERGHEI_SWE_DRY_RUNOFF_START_DT : " << SERGHEI_SWE_DRY_RUNOFF_START_DT << std::endl;
        logFile << "SERGHEI_MAXFLOOD : " << SERGHEI_MAXFLOOD << std::endl;
        logFile << "SERGHEI_REAL : " << SERGHEI_REAL << std::endl;
        logFile << "SERGHEI_NC_MODE : " << SERGHEI_NC_MODE << std::endl;
        logFile << "SERGHEI_NC_REAL : ";
        if (SERGHEI_NC_REAL == NC_FLOAT)
          logFile << "NC_FLOAT";
        if (SERGHEI_NC_REAL == NC_DOUBLE)
          logFile << "NC_DOUBLE";
        logFile << std::endl;
        logFile << "SERGHEI_WRITE_HZ : " << SERGHEI_WRITE_HZ << std::endl;
        logFile << "SERGHEI_WRITE_SUBDOMS : " << SERGHEI_WRITE_SUBDOMS << std::endl;
        logFile << "SERGHEI_INPUT_NETCDF : " << SERGHEI_INPUT_NETCDF << std::endl;
        logFile << "SERGHEI_NC_ENABLE_NAN : " << SERGHEI_NC_ENABLE_NAN << std::endl;
        logFile << "SERGHEI_NETCDF_FORCING : " << SERGHEI_NETCDF_FORCING << std::endl;

        logFile << "\n-------------------------\nDEBUG FLAGS" << std::endl;
        logFile << "SERGHEI_DEBUG_PARALLEL_DECOMPOSITION : " << SERGHEI_DEBUG_PARALLEL_DECOMPOSITION << std::endl;
        logFile << "SERGHEI_DEBUG_WORKFLOW : " << SERGHEI_DEBUG_WORKFLOW << std::endl;
        logFile << "SERGHEI_DEBUG_KOKKOS_SETUP : " << SERGHEI_DEBUG_KOKKOS_SETUP << std::endl;
        logFile << "SERGHEI_DEBUG_BOUNDARY : " << SERGHEI_DEBUG_BOUNDARY << std::endl;
        logFile << "SERGHEI_DEBUG_DT : " << SERGHEI_DEBUG_DT << std::endl;
        logFile << "SERGHEI_DEBUG_TOOLS : " << SERGHEI_DEBUG_TOOLS << std::endl;
        logFile << "SERGHEI_DEBUG_MASS_CONS : " << SERGHEI_DEBUG_MASS_CONS << std::endl;
        logFile << "SERGHEI_DEBUG_INFILTRATION : " << SERGHEI_DEBUG_INFILTRATION << std::endl;
        logFile << "SERGHEI_DEBUG_MPI : " << SERGHEI_DEBUG_MPI << std::endl;
        logFile << "SERGHEI_DEBUG_SCALAR_TRANSPORT : " << SERGHEI_DEBUG_SCALAR_TRANSPORT << std::endl;
        logFile << "SERGHEI_DEBUG_SEDIMENT : " << SERGHEI_DEBUG_SEDIMENT << std::endl;
        logFile << "SERGHEI_DEBUG_OUTPUT : " << SERGHEI_DEBUG_OUTPUT << std::endl;
        logFile << "SERGHEI_DEBUG_INPUT_NETCDF" << SERGHEI_DEBUG_INPUT_NETCDF << std::endl;
        logFile << "SERGHEI_DEBUG_RAINFALL" << SERGHEI_DEBUG_RAINFALL << std::endl;
      }
      std::cout << GOK << "Log file written" << std::endl;
    }
  }

  // Writes VTK file for subsurface
#if SERGHEI_RE_MODEL

  void outputInitNETCDFSub(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir)
  {
    int dimids[4], nxhalo, nyhalo;
    MPI_Offset st[3], ct[3];
    realArr xCoord = realArr("xCoord", gdom.nx);
    realArr yCoord = realArr("yCoord", gdom.ny);
    realArr zCoord = realArr("zCoord", gdom.nz);
    realArr data = realArr("data", gdom.nCell);
    static char title[] = "seconds";
    std::string filename;

    filename = dir + "output_subsurface.nc";

    // Create the file
    ncwrap(ncmpi_create(MPI_COMM_WORLD, filename.c_str(), NC_CLOBBER, MPI_INFO_NULL, &ncid), __LINE__);

    // Create the dimensions
    ncwrap(ncmpi_def_dim(ncid, "t", (MPI_Offset)NC_UNLIMITED, &tDim), __LINE__);
    ncwrap(ncmpi_def_dim(ncid, "x", (MPI_Offset)gdom.nx_glob, &xDim), __LINE__);
    ncwrap(ncmpi_def_dim(ncid, "y", (MPI_Offset)gdom.ny_glob, &yDim), __LINE__);
    ncwrap(ncmpi_def_dim(ncid, "z", (MPI_Offset)gdom.nz_glob, &zDim), __LINE__);
    // Create the variables
    dimids[0] = tDim;
    ncwrap(ncmpi_def_var(ncid, "t", NC_DOUBLE, 1, dimids, &tVar), __LINE__);
    ncwrap(ncmpi_put_att_text(ncid, tVar, "units", strlen(title), title), __LINE__);
    dimids[0] = xDim;
    ncwrap(ncmpi_def_var(ncid, "x", NC_DOUBLE, 1, dimids, &xVar), __LINE__);
    dimids[0] = yDim;
    ncwrap(ncmpi_def_var(ncid, "y", NC_DOUBLE, 1, dimids, &yVar), __LINE__);
    dimids[0] = zDim;
    ncwrap(ncmpi_def_var(ncid, "z", NC_DOUBLE, 1, dimids, &zVar), __LINE__);

    dimids[0] = zDim;
    dimids[1] = yDim;
    dimids[2] = xDim;
    ncwrap(ncmpi_def_var(ncid, "z3d", NC_DOUBLE, 3, dimids, &z3Var), __LINE__);
    dimids[0] = tDim;
    dimids[1] = zDim;
    dimids[2] = yDim;
    dimids[3] = xDim;
    ncwrap(ncmpi_def_var(ncid, "hd", NC_DOUBLE, 4, dimids, &hdVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "wc", NC_DOUBLE, 4, dimids, &wcVar), __LINE__);

    //! zzb 20241011 [CODE 1] 新增 qx,qy,qz, evap
    ncwrap(ncmpi_def_var(ncid, "qx", NC_DOUBLE, 4, dimids, &qxVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "qy", NC_DOUBLE, 4, dimids, &qyVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "qz", NC_DOUBLE, 4, dimids, &qzVar), __LINE__);

    if (gdom.hasET == 1)
    {
      ncwrap(ncmpi_def_var(ncid, "evapGW", NC_DOUBLE, 4, dimids, &evapGWVar), __LINE__);
      ncwrap(ncmpi_def_var(ncid, "transpGW", NC_DOUBLE, 4, dimids, &transpGWVar), __LINE__);
    }

    // End "define" mode
    ncwrap(ncmpi_enddef(ncid), __LINE__);

    // Compute x, y, z coordinates
    Kokkos::parallel_for(gdom.nx, KOKKOS_LAMBDA(int i) { xCoord(i) = gdom.xll + (par.i_beg + i + 0.5) * gdom.dx; });
    Kokkos::parallel_for(gdom.ny, KOKKOS_LAMBDA(int j) { yCoord(j) = gdom.yll + gdom.ny_glob * gdom.dx - (par.j_beg + j + 0.5) * gdom.dx; });
    Kokkos::parallel_for(gdom.nz, KOKKOS_LAMBDA(int k) { zCoord(k) = -(k + 0.5) * gdom.dz(k); });
    Kokkos::fence();

    // Write out x, y coordinates
    st[0] = par.i_beg;
    ct[0] = gdom.nx;
    ncwrap(ncmpi_put_vara_double_all(ncid, xVar, st, ct, xCoord.data()), __LINE__);
    st[0] = par.j_beg;
    ct[0] = gdom.ny;
    ncwrap(ncmpi_put_vara_double_all(ncid, yVar, st, ct, yCoord.data()), __LINE__);
    st[0] = 0;
    ct[0] = gdom.nz;
    ncwrap(ncmpi_put_vara_double_all(ncid, zVar, st, ct, zCoord.data()), __LINE__);

    // Write z for the 3D domain
    nxhalo = gdom.nx + 2 * gdom.hc;
    nyhalo = gdom.ny + 2 * gdom.hc;
    st[0] = 0;
    st[1] = par.j_beg;
    st[2] = par.i_beg;
    ct[0] = gdom.nz;
    ct[1] = gdom.ny;
    ct[2] = gdom.nx;
    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
            data(idom) = gdom.z(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, z3Var, st, ct, data.data()), __LINE__);

    writeGwNETCDF(gw, gdom, par);
    ncwrap(ncmpi_close(ncid), __LINE__);
  }

  void outputNETCDFSubsurface(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir)
  {
    std::string filename;
    filename = dir + "output_subsurface.nc";
    // Create the file
    ncwrap(ncmpi_open(MPI_COMM_WORLD, filename.c_str(), NC_WRITE, MPI_INFO_NULL, &ncid), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "hd", &hdVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "wc", &wcVar), __LINE__);

    //! zzb 20241011 [CODE 1]
    ncwrap(ncmpi_inq_varid(ncid, "qx", &qxVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "qy", &qyVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "qz", &qzVar), __LINE__);
    if (gdom.hasET == 1)
    {
      ncwrap(ncmpi_inq_varid(ncid, "evapGW", &evapGWVar), __LINE__);
      ncwrap(ncmpi_inq_varid(ncid, "transpGW", &transpGWVar), __LINE__);
    }

    writeGwNETCDF(gw, gdom, par);

    ncwrap(ncmpi_close(ncid), __LINE__);
  }

  void writeGwNETCDF(const GwState &gw, GwDomain const &gdom, Parallel const &par)
  {
    realArr data = realArr("data", gdom.nCell);
    MPI_Offset st[4], ct[4];
    double timeIter[numOut + 1];
    int nxhalo = gdom.nx + 2 * gdom.hc, nyhalo = gdom.ny + 2 * gdom.hc;
    // write t. As the first one is written in the first iteration we should add +1
    for (int i = 0; i < numOut + 1; i++)
    {
      timeIter[i] = i * 1.0;
    }
    st[0] = 0;
    ct[0] = numOut + 1;
    ncwrap(ncmpi_put_vara_double_all(ncid, tVar, st, ct, timeIter), __LINE__);

    st[0] = numOut;
    st[1] = 0;
    st[2] = par.j_beg;
    st[3] = par.i_beg;
    ct[0] = 1;
    ct[1] = gdom.nz;
    ct[2] = gdom.ny;
    ct[3] = gdom.nx;

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
            data(idom) = gw.h(iGlob,1); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, hdVar, st, ct, data.data()), __LINE__);

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
            data(idom) = gw.wc(iGlob,1); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, wcVar, st, ct, data.data()), __LINE__);

    //! zzb 20241011 [CODE 1] 新增 qx, qy, qz 输出
    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
            data(idom) = gw.q(iGlob,0); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, qxVar, st, ct, data.data()), __LINE__);

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
            data(idom) = gw.q(iGlob,1); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, qyVar, st, ct, data.data()), __LINE__);

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
            data(idom) = gw.q(iGlob,2); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, qzVar, st, ct, data.data()), __LINE__);

    //! zzb subsurface evaporation and transpiration
    if (gdom.hasET == 1)
    {
      Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
              int ii, jj, kk, iGlob;
              gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
              iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
              data(idom) = gw.evapGW(iGlob); });
      Kokkos::fence();
      ncwrap(ncmpi_put_vara_double_all(ncid, evapGWVar, st, ct, data.data()), __LINE__);
    }

    if (gdom.hasET == 1)
    {
      Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
              int ii, jj, kk, iGlob;
              gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
              iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
              data(idom) = gw.transpGW(iGlob); });
      Kokkos::fence();
      ncwrap(ncmpi_put_vara_double_all(ncid, transpGWVar, st, ct, data.data()), __LINE__);
    }
  }

  // VTK writer for the subsurface
  void outputVTKsubsurface(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir)
  {

    std::string filename;
    filename = dir + "result_subsurface" + std::to_string(numOut) + ".vtk";
    real *xCoord;
    real *yCoord;
    real *zCoord;
    real *data_cpu;

    int ncells = gdom.ny * gdom.nx * gdom.nz;

    int nVars = 2; // 2 variables to write: h, wc

    realArr data = realArr("data", nVars * ncells);
#ifdef KOKKOS_ENABLE_CUDA_
    cudaMallocHost(&data_cpu, nVars * ncells * sizeof(real));
#else
    data_cpu = data.data();
#endif

    xCoord = (real *)malloc((gdom.nx + 1) * sizeof(real));
    yCoord = (real *)malloc((gdom.ny + 1) * sizeof(real));
    zCoord = (real *)malloc((gdom.nz + 1) * sizeof(real));

    // Compute x, y coordinates
    for (int i = 0; i < gdom.nx + 1; i++)
    {
      xCoord[i] = gdom.xll + (par.i_beg + i) * gdom.dx;
    };

    for (int j = 0; j < gdom.ny + 1; j++)
    {
      yCoord[j] = gdom.yll + gdom.ny_glob * gdom.dx - (par.j_beg + j) * gdom.dy;
    };

    int nnodes = (gdom.nx + 1) * (gdom.ny + 1) * (gdom.nz + 1);

    Kokkos::parallel_for(ncells, KOKKOS_LAMBDA(int idom) {
          int i,j,k;
          int nxhalo = gdom.nx+2*gdom.hc;
          int nyhalo = gdom.ny+2*gdom.hc;
          gdom.unpackIndicesGw(idom,gdom.nz,gdom.ny,gdom.nx,k,j,i);
          int iGlob = (gdom.hc+k)*nxhalo*nyhalo + (gdom.hc+j)*nxhalo + i + gdom.hc;
          data(idom) = gw.h(iGlob,1);
          data(idom+ncells) = gw.wc(iGlob,1);
          // get z coordinate
          if (i == 0 & j == 0) {zCoord[k] = gdom.z(iGlob);} });
    Kokkos::fence();

#ifdef KOKKOS_ENABLE_CUDA
    cudaMemcpyAsync(data_cpu, data.data(), nVars * ncells * sizeof(real), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
#endif

    int i, j, k, iGlob;

    std::ofstream fOutStream(filename);
    if (fOutStream.is_open())
    {
      fOutStream << "# vtk DataFile Version 3.0.\nOutputfile\nASCII\nDATASET UNSTRUCTURED_GRID\n";

      fOutStream << "POINTS " << nnodes << " double\n";
      for (k = 0; k <= gdom.nz; k++)
      {
        for (j = 0; j <= gdom.ny; j++)
        {
          for (i = 0; i <= gdom.nx; i++)
          {
            fOutStream << std::setprecision(9) << xCoord[i] << " " << yCoord[j] << " " << zCoord[k] << " 0.0\n";
          }
        }
      }
      fOutStream << "CELLS " << ncells << " " << 2 * ncells << "\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        int ii, jj, kk;
        gdom.unpackIndicesGw(iGlob, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
        fOutStream << "4 " << jj * (gdom.nx + 1) + ii << " " << (jj + 1) * (gdom.nx + 1) + ii << " " << (jj + 1) * (gdom.nx + 1) + ii + 1 << " " << jj * (gdom.nx + 1) + ii + 1 << "\n";
      }

      fOutStream << "CELL_TYPES " << ncells << "\n";
      for (i = 0; i < ncells; i++)
      {
        fOutStream << "9\n";
      }

      fOutStream << "CELL_DATA " << ncells << "\n";
      fOutStream << "SCALARS h double\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob] << "\n";
      }

      fOutStream << "SCALARS wc double\n";
      fOutStream << "LOOKUP_TABLE default\n";
      for (iGlob = 0; iGlob < ncells; iGlob++)
      {
        fOutStream << std::setprecision(9) << data_cpu[iGlob + ncells] << "\n";
      }

      fOutStream.close();
    }
  }

  // Write header and initial state for time series files
  int writeSubTimeSeriesIni(GwDomain const &gdom, GwIntegrator const &gint, Parallel const &par, std::string dir)
  {

    numObs = 0;
    std::string filename = dir + "SubsurfaceTimeSeries.out";
    SubsurfaceOutputFile.open(filename);

    if (SubsurfaceOutputFile.is_open())
    {
      // Write the header
      SubsurfaceOutputFile << "Time ";
      SubsurfaceOutputFile << "SubSurfaceVolume [m3] ";
      SubsurfaceOutputFile << "ExchangeRate [m3/s] ";
      SubsurfaceOutputFile << "BoundaryInflow ";
      SubsurfaceOutputFile << "BoundaryOutflow ";
      SubsurfaceOutputFile << "Source/SinkInflow [m3/s] ";
      SubsurfaceOutputFile << "Source/SinkOutflow [m3/s] ";
      // [CODE 1] 新增根系蒸腾与土壤蒸发边界流量输出
      SubsurfaceOutputFile << "Source/SinkOutflow_RootTransp[m3/s] ";
      SubsurfaceOutputFile << "Source/SinkOutflow_SoilEvap[m3/s] ";
      SubsurfaceOutputFile << std::endl;
    }
    else
    {
      std::cerr << RERROR "Could not create domainTimeSeries.out file" << std::endl;
      return 0;
    }

    if (par.masterproc)
    {
      writeSubsurfaceTimeSeries(gdom, gint);
    }
    numObs++;
    return 1;
  }

  void writeSubsurfaceTimeSeries(GwDomain const &gdom, GwIntegrator const &gint)
  {
    // Write the data
    std::cout.precision(OUTPUT_PRECISION);
    SubsurfaceOutputFile << std::scientific << gdom.etime << " ";
    SubsurfaceOutputFile << std::scientific << gint.Vtot_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.Vexch_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QinBC_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QoutBC_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QinSS_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QoutSS_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QoutSS_RootTransp_glob << " "; // [CODE 1]
    SubsurfaceOutputFile << std::scientific << gint.QoutSS_SoilEvap_glob << " ";   // [CODE 1]
    SubsurfaceOutputFile << std::endl;
  }

  // =========================================================================
  // 【修复 4：根系含水率输出】[CODE 1]
  // =========================================================================
  int writeRootZoneTimeSeriesIni(Domain const &dom, Parallel const &par, std::string dir, std::vector<GwSS> &gwss)
  {
    numObs = 0;
    if (par.masterproc)
    {
      std::string filename = dir + "RootZoneWaterContentTimeSeries.out";
      RootZoneWaterContentOutputFile.open(filename);

      if (RootZoneWaterContentOutputFile.is_open())
      {
        RootZoneWaterContentOutputFile << "Time RootZoneWaterContent_Mean [-] \n";
        writeRootZoneWaterContentTimeSeries(dom, gwss);
      }
      else
      {
        std::cerr << RERROR "Could not create RootZoneWaterContentOutputFile.out file" << std::endl;
        return 0;
      }
    }
    numObs++;
    return 1;
  }

  void writeRootZoneWaterContentTimeSeries(Domain const &dom, std::vector<GwSS> &gwss)
  {
    std::cout.precision(OUTPUT_PRECISION);
    RootZoneWaterContentOutputFile << std::scientific << dom.etime << " ";
    for (int k = 0; k < gwss.size(); k++)
    {
      RootZoneWaterContentOutputFile << std::scientific << gwss[k].wc_root_zone << " ";
    }
    RootZoneWaterContentOutputFile << std::endl;
  }

#endif

  // =========================================================================
  // [CODE 1] 地下水多组分溶质运移核心 IO 代码区
  // =========================================================================
#if SERGHEI_SUBSURFACE_TRANSPORT

  int writeRTSubsurfaceTimeSeriesIni(Domain const &dom, RTStateGW const &rt, RTIntegrator const &rtint, Parallel const &par, std::string dir)
  {
    numObs = 0;
    if (par.masterproc)
    {
      std::string filename = dir + "RTSubsurfaceTimeSeries.out";
      RTSubsurfaceOutputFile.open(filename);

      if (RTSubsurfaceOutputFile.is_open())
      {
        RTSubsurfaceOutputFile << "Time[s] ";
        for (int is = 0; is < rt.n_mass; is++)
        {
          std::string specName = "spec" + std::to_string(is);
          RTSubsurfaceOutputFile << specName << "_LiquidMass[mg] ";
          RTSubsurfaceOutputFile << specName << "_SolidMass[mg] ";
          RTSubsurfaceOutputFile << specName << "_MassExch[mg/s] ";
          RTSubsurfaceOutputFile << specName << "_BCInflow[mg/s] ";
          RTSubsurfaceOutputFile << specName << "_BCOutflow[mg/s] ";
          RTSubsurfaceOutputFile << specName << "_SSInflow[mg/s] ";
          RTSubsurfaceOutputFile << specName << "_SSOutflow[mg/s] ";
          RTSubsurfaceOutputFile << specName << "_Reaction[mg/s] ";
        }
        RTSubsurfaceOutputFile << std::endl;
        writeRTSubsurfaceTimeSeries(dom, rtint);
      }
      else
      {
        std::cerr << RERROR "Could not create RTSubsurfaceTimeSeries.out file" << std::endl;
        return 0;
      }
    }
    numObs++;
    return 1;
  }

  void writeRTSubsurfaceTimeSeries(Domain const &dom, RTIntegrator const &rtint)
  {
    std::cout.precision(OUTPUT_PRECISION);
    RTSubsurfaceOutputFile << std::scientific << dom.etime << " ";

    for (int is = 0; is < rtint.n_mass; is++)
    {
      RTSubsurfaceOutputFile << std::scientific << rtint.LiquidMassTot_spec_glob(is) << " ";
      RTSubsurfaceOutputFile << std::scientific << rtint.SolidMassTot_spec_glob(is) << " ";
      RTSubsurfaceOutputFile << std::scientific << rtint.MassExch_spec_glob(is) << " ";
      RTSubsurfaceOutputFile << std::scientific << rtint.QMassInBC_spec_glob(is) << " ";
      RTSubsurfaceOutputFile << std::scientific << rtint.QMassOutBC_spec_glob(is) << " ";
      RTSubsurfaceOutputFile << std::scientific << rtint.QMassInSS_spec_glob(is) << " ";
      RTSubsurfaceOutputFile << std::scientific << rtint.QMassOutSS_spec_glob(is) << " ";
      RTSubsurfaceOutputFile << std::scientific << rtint.QMassReaction_spec_glob(is) << " ";
    }
    RTSubsurfaceOutputFile << std::endl;
  }

  void outputInitNETCDFRT(const RTStateGW &rt, GwDomain const &gdom, Parallel const &par, std::string dir)
  {
    int dimids[5], nxhalo, nyhalo;
    MPI_Offset st[3], ct[3];
    realArr xCoord = realArr("xCoord", gdom.nx);
    realArr yCoord = realArr("yCoord", gdom.ny);
    realArr zCoord = realArr("zCoord", gdom.nz);
    realArr data = realArr("data", gdom.nCell);
    realArr2 data2 = realArr2("data2", rt.n_mass, gdom.nCell);
    static char title[] = "seconds";
    std::string filename;

    filename = dir + "output_transport.nc";
    ncwrap(ncmpi_create(MPI_COMM_WORLD, filename.c_str(), NC_CLOBBER, MPI_INFO_NULL, &ncid), __LINE__);
    ncwrap(ncmpi_def_dim(ncid, "t", (MPI_Offset)NC_UNLIMITED, &tDim), __LINE__);
    ncwrap(ncmpi_def_dim(ncid, "n", (MPI_Offset)rt.n_mass, &nDim), __LINE__);
    ncwrap(ncmpi_def_dim(ncid, "x", (MPI_Offset)gdom.nx_glob, &xDim), __LINE__);
    ncwrap(ncmpi_def_dim(ncid, "y", (MPI_Offset)gdom.ny_glob, &yDim), __LINE__);
    ncwrap(ncmpi_def_dim(ncid, "z", (MPI_Offset)gdom.nz_glob, &zDim), __LINE__);

    dimids[0] = tDim;
    ncwrap(ncmpi_def_var(ncid, "t", NC_DOUBLE, 1, dimids, &tVar), __LINE__);
    ncwrap(ncmpi_put_att_text(ncid, tVar, "units", strlen(title), title), __LINE__);
    dimids[0] = nDim;
    ncwrap(ncmpi_def_var(ncid, "n", NC_DOUBLE, 1, dimids, &nVar), __LINE__);
    dimids[0] = xDim;
    ncwrap(ncmpi_def_var(ncid, "x", NC_DOUBLE, 1, dimids, &xVar), __LINE__);
    dimids[0] = yDim;
    ncwrap(ncmpi_def_var(ncid, "y", NC_DOUBLE, 1, dimids, &yVar), __LINE__);
    dimids[0] = zDim;
    ncwrap(ncmpi_def_var(ncid, "z", NC_DOUBLE, 1, dimids, &zVar), __LINE__);

    dimids[0] = zDim;
    dimids[1] = yDim;
    dimids[2] = xDim;
    ncwrap(ncmpi_def_var(ncid, "z3d", NC_DOUBLE, 3, dimids, &z3Var), __LINE__);
    dimids[0] = tDim;
    dimids[1] = nDim;
    dimids[2] = zDim;
    dimids[3] = yDim;
    dimids[4] = xDim;
    ncwrap(ncmpi_def_var(ncid, "c", NC_DOUBLE, 5, dimids, &cVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "c_solid", NC_DOUBLE, 5, dimids, &c_solidVar), __LINE__);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
    dimids[0] = tDim;
    dimids[1] = nDim;
    dimids[2] = zDim;
    dimids[3] = yDim;
    dimids[4] = xDim;
    ncwrap(ncmpi_def_var(ncid, "c_difxz", NC_DOUBLE, 5, dimids, &c_difxzVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "c_difzx", NC_DOUBLE, 5, dimids, &c_difzxVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "c_difxx", NC_DOUBLE, 5, dimids, &c_difxxVar), __LINE__);
    ncwrap(ncmpi_def_var(ncid, "c_difzz", NC_DOUBLE, 5, dimids, &c_difzzVar), __LINE__);
#endif

    ncwrap(ncmpi_enddef(ncid), __LINE__);

    Kokkos::parallel_for(gdom.nx, KOKKOS_LAMBDA(int i) { xCoord(i) = gdom.xll + (par.i_beg + i + 0.5) * gdom.dx; });
    Kokkos::parallel_for(gdom.ny, KOKKOS_LAMBDA(int j) { yCoord(j) = gdom.yll + gdom.ny_glob * gdom.dx - (par.j_beg + j + 0.5) * gdom.dx; });
    Kokkos::parallel_for(gdom.nz, KOKKOS_LAMBDA(int k) { zCoord(k) = -(k + 0.5) * gdom.dz(k); });
    Kokkos::fence();

    st[0] = par.i_beg;
    ct[0] = gdom.nx;
    ncwrap(ncmpi_put_vara_double_all(ncid, xVar, st, ct, xCoord.data()), __LINE__);
    st[0] = par.j_beg;
    ct[0] = gdom.ny;
    ncwrap(ncmpi_put_vara_double_all(ncid, yVar, st, ct, yCoord.data()), __LINE__);
    st[0] = 0;
    ct[0] = gdom.nz;
    ncwrap(ncmpi_put_vara_double_all(ncid, zVar, st, ct, zCoord.data()), __LINE__);

    nxhalo = gdom.nx + 2 * gdom.hc;
    nyhalo = gdom.ny + 2 * gdom.hc;
    st[0] = 0;
    st[1] = par.j_beg;
    st[2] = par.i_beg;
    ct[0] = gdom.nz;
    ct[1] = gdom.ny;
    ct[2] = gdom.nx;
    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
        int ii, jj, kk, iGlob;
        gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
        iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
        data(idom) = gdom.z(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, z3Var, st, ct, data.data()), __LINE__);

    writeRTNETCDF(rt, gdom, par);
    ncwrap(ncmpi_close(ncid), __LINE__);
  }

  void outputNETCDFTransport(const RTStateGW &rt, GwDomain const &gdom, Parallel const &par, std::string dir)
  {
    std::string filename;
    filename = dir + "output_transport.nc";
    ncwrap(ncmpi_open(MPI_COMM_WORLD, filename.c_str(), NC_WRITE, MPI_INFO_NULL, &ncid), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "c", &cVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "c_solid", &c_solidVar), __LINE__);
#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
    ncwrap(ncmpi_inq_varid(ncid, "c_difxz", &c_difxzVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "c_difzx", &c_difzxVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "c_difxx", &c_difxxVar), __LINE__);
    ncwrap(ncmpi_inq_varid(ncid, "c_difzz", &c_difzzVar), __LINE__);
#endif
    writeRTNETCDF(rt, gdom, par);
    ncwrap(ncmpi_close(ncid), __LINE__);
  }

  void writeRTNETCDF(const RTStateGW &rt, GwDomain const &gdom, Parallel const &par)
  {
    realArr data = realArr("data", gdom.nCell);
    realArr2 data2 = realArr2("data2", rt.n_mass, gdom.nCell);
    MPI_Offset st[5], ct[5];
    double timeIter[numOut + 1];
    double n_massIter[rt.n_mass];
    int nxhalo = gdom.nx + 2 * gdom.hc, nyhalo = gdom.ny + 2 * gdom.hc;

    for (int i = 0; i < numOut + 1; i++)
    {
      timeIter[i] = i * 1.0;
    }
    st[0] = 0;
    ct[0] = numOut + 1;
    ncwrap(ncmpi_put_vara_double_all(ncid, tVar, st, ct, timeIter), __LINE__);
    for (int i = 0; i < rt.n_mass; i++)
    {
      n_massIter[i] = i * 1.0;
    }
    st[0] = 0;
    ct[0] = rt.n_mass;
    ncwrap(ncmpi_put_vara_double_all(ncid, nVar, st, ct, n_massIter), __LINE__);

    st[0] = numOut;
    st[1] = 0;
    st[2] = 0;
    st[3] = par.j_beg;
    st[4] = par.i_beg;
    ct[0] = 1;
    ct[1] = rt.n_mass;
    ct[2] = gdom.nz;
    ct[3] = gdom.ny;
    ct[4] = gdom.nx;

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
        int ii, jj, kk, iGlob;
        gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
        iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
        for (int s = 0; s < rt.n_mass; ++s)
        {
          data2(s, idom) = rt.c(s, iGlob,1);
        } });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, cVar, st, ct, data2.data()), __LINE__);

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
        int ii, jj, kk, iGlob;
        gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
        iGlob = (gdom.hc+kk)*nxhalo*nyhalo + (gdom.hc+jj)*nxhalo + ii + gdom.hc;
                for (int s = 0; s < rt.n_mass; ++s)
        {
          data2(s, idom) = rt.c_solid(s, iGlob, 1);
        } });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, c_solidVar, st, ct, data2.data()), __LINE__);

#if DEBUG_SERGHEI_SUBSURFACE_TRANSPORT
    st[0] = numOut;
    st[1] = 0;
    st[2] = 0;
    st[3] = par.j_beg;
    st[4] = par.i_beg;
    ct[0] = 1;
    ct[1] = rt.n_mass;
    ct[2] = gdom.nz;
    ct[3] = gdom.ny;
    ct[4] = gdom.nx;
    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
      int ii, jj, kk, iGlob;
      gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
      iGlob = (gdom.hc + kk) * nxhalo * nyhalo + (gdom.hc + jj) * nxhalo + ii + gdom.hc;
      for (int s = 0; s < rt.n_mass; ++s)
      {
        data2(s, idom) = rt.dcal(s, iGlob, 4); // dcalxz
      } });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, c_difxzVar, st, ct, data2.data()), __LINE__);

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
      int ii, jj, kk, iGlob;
      gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
      iGlob = (gdom.hc + kk) * nxhalo * nyhalo + (gdom.hc + jj) * nxhalo + ii + gdom.hc;
      for (int s = 0; s < rt.n_mass; ++s)
      {
        data2(s, idom) = rt.dcal(s, iGlob, 7); // dcalzx
      } });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, c_difzxVar, st, ct, data2.data()), __LINE__);

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
      int ii, jj, kk, iGlob;
      gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
      iGlob = (gdom.hc + kk) * nxhalo * nyhalo + (gdom.hc + jj) * nxhalo + ii + gdom.hc;
      for (int s = 0; s < rt.n_mass; ++s)
      {
        data2(s, idom) = rt.dcal(s, iGlob, 0); // dcalxx
      } });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, c_difxxVar, st, ct, data2.data()), __LINE__);

    Kokkos::parallel_for(gdom.nCell, KOKKOS_LAMBDA(int idom) {
      int ii, jj, kk, iGlob;
      gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
      iGlob = (gdom.hc + kk) * nxhalo * nyhalo + (gdom.hc + jj) * nxhalo + ii + gdom.hc;
      for (int s = 0; s < rt.n_mass; ++s)
      {
        data2(s, idom) = rt.dcal(s, iGlob, 2); // dcalzz
      } });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, c_difzzVar, st, ct, data2.data()), __LINE__);
#endif
  }
#endif

  // =========================================================================
  // [CODE 1] 地表水溶质运移核心 IO 代码区
  // =========================================================================
#if SERGHEI_SURFACE_TRANSPORT

  int writeRTSurfaceTimeSeriesIni(Domain const &dom, RTStateSW const &rtsw, RTIntegratorSW const &rtintsw, Parallel const &par, std::string dir)
  {
    numObs = 0;
    if (par.masterproc)
    {
      std::string filename = dir + "RTSurfaceTimeSeries.out";
      RTSurfaceOutputFile.open(filename);

      if (RTSurfaceOutputFile.is_open())
      {
        RTSurfaceOutputFile << "Time[s] ";
        for (int is = 0; is < rtsw.n_mass; is++)
        {
          std::string specName = "spec" + std::to_string(is);
          RTSurfaceOutputFile << specName << "_SurfaceMass[mg] ";
          RTSurfaceOutputFile << specName << "_MassExch[mg/s] ";
          RTSurfaceOutputFile << specName << "_BCInflow[mg/s] ";
          RTSurfaceOutputFile << specName << "_BCOutflow[mg/s] ";
          RTSurfaceOutputFile << specName << "_SSInflow[mg/s] ";
          RTSurfaceOutputFile << specName << "_SSOutflow[mg/s] ";
          RTSurfaceOutputFile << specName << "_Reaction[mg/s] ";
          RTSurfaceOutputFile << specName << "_RainInflow[mg/s] ";
        }
        RTSurfaceOutputFile << std::endl;
        writeRTSurfaceTimeSeries(dom, rtintsw);
      }
      else
      {
        std::cerr << RERROR "Could not create RTSurfaceTimeSeries.out file" << std::endl;
        return 0;
      }
    }
    numObs++;
    return 1;
  }

  void writeRTSurfaceTimeSeries(Domain const &dom, RTIntegratorSW const &rtintsw)
  {
    std::cout.precision(OUTPUT_PRECISION);
    RTSurfaceOutputFile << std::scientific << dom.etime << " ";

    for (int is = 0; is < rtintsw.n_mass; is++)
    {
      RTSurfaceOutputFile << std::scientific << rtintsw.MassTot_spec_glob(is) << " ";
      RTSurfaceOutputFile << std::scientific << rtintsw.MassExch_spec_glob(is) << " ";
      RTSurfaceOutputFile << std::scientific << rtintsw.QMassInBC_spec_glob(is) << " ";
      RTSurfaceOutputFile << std::scientific << rtintsw.QMassOutBC_spec_glob(is) << " ";
      RTSurfaceOutputFile << std::scientific << rtintsw.QMassInSS_spec_glob(is) << " ";
      RTSurfaceOutputFile << std::scientific << rtintsw.QMassOutSS_spec_glob(is) << " ";
      RTSurfaceOutputFile << std::scientific << rtintsw.QMassReaction_spec_glob(is) << " ";
      RTSurfaceOutputFile << std::scientific << rtintsw.QMassRain_spec_glob(is) << " ";
    }
    RTSurfaceOutputFile << std::endl;
  }

  void outputInitNETCDFRTSW(const RTStateSW &rtsw, const State &state, Domain const &dom, Parallel const &par, std::string dir)
  {
    int dimids[4];
    MPI_Offset st[3], ct[3];
    doubleArr xCoord = doubleArr("xCoord", dom.nx);
    doubleArr yCoord = doubleArr("yCoord", dom.ny);
    ncArr data = ncArr("data", dom.nCell);
    static char timeUnits[] = "seconds";
    std::string filename;
    std::string longname;
    std::string units;

    filename = dir + "output_transportSW.nc";
    // Create the file
    ncwrap(ncmpi_create(MPI_COMM_WORLD, filename.c_str(), SERGHEI_NC_MODE, MPI_INFO_NULL, &ncid), __LINE__, par.myrank);

    // Create the dimensions
    ncwrap(ncmpi_def_dim(ncid, "t", (MPI_Offset)NC_UNLIMITED, &tDim), __LINE__, par.myrank);
    ncwrap(ncmpi_def_dim(ncid, "n", (MPI_Offset)rtsw.n_mass, &nDim), __LINE__, par.myrank);
    ncwrap(ncmpi_def_dim(ncid, "x", (MPI_Offset)dom.nx_glob, &xDim), __LINE__, par.myrank);
    ncwrap(ncmpi_def_dim(ncid, "y", (MPI_Offset)dom.ny_glob, &yDim), __LINE__, par.myrank);

    // Create the variables
    dimids[0] = tDim;
    ncwrap(ncmpi_def_var(ncid, "t", NC_DOUBLE, 1, dimids, &tVar), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, tVar, "units", strlen(timeUnits), timeUnits), __LINE__, par.myrank);
    dimids[0] = nDim;
    ncwrap(ncmpi_def_var(ncid, "n", NC_DOUBLE, 1, dimids, &nVar), __LINE__, par.myrank);
    dimids[0] = xDim;
    ncwrap(ncmpi_def_var(ncid, "x", NC_DOUBLE, 1, dimids, &xVar), __LINE__, par.myrank);
    dimids[0] = yDim;
    ncwrap(ncmpi_def_var(ncid, "y", NC_DOUBLE, 1, dimids, &yVar), __LINE__, par.myrank);

    dimids[0] = tDim;
    dimids[1] = nDim;
    dimids[2] = yDim;
    dimids[3] = xDim;

    ncwrap(ncmpi_def_var(ncid, "csw", NC_DOUBLE, 4, dimids, &cSWVar), __LINE__, par.myrank);
#if SERGHEI_NC_ENABLE_NAN
    double fill_val = SERGHEI_NAN;
    ncwrap(ncmpi_put_att_double(ncid, cSWVar, "_FillValue", NC_DOUBLE, 1, &fill_val), __LINE__, par.myrank);
#endif
    // End "define" mode
    ncwrap(ncmpi_enddef(ncid), __LINE__, par.myrank);

    Kokkos::parallel_for("compute_grid_coord_x", dom.nx, KOKKOS_LAMBDA(int i) { xCoord(i) = dom.xll + (par.i_beg + i + 0.5) * dom.dxConst; });
    Kokkos::parallel_for("compute_grid_coord_y", dom.ny, KOKKOS_LAMBDA(int j) { yCoord(j) = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j + 0.5) * dom.dxConst; });
    Kokkos::fence();

    st[0] = par.i_beg;
    ct[0] = dom.nx;
    ncwrap(ncmpi_put_vara_double_all(ncid, xVar, st, ct, xCoord.data()), __LINE__, par.myrank);
    st[0] = par.j_beg;
    ct[0] = dom.ny;
    ncwrap(ncmpi_put_vara_double_all(ncid, yVar, st, ct, yCoord.data()), __LINE__, par.myrank);

    st[0] = par.j_beg;
    st[1] = par.i_beg;
    ct[0] = dom.ny;
    ct[1] = dom.nx;

#if SERGHEI_NC_ENABLE_MISSING_VALUE
    real missing_value = dom.MISSING_VALUE;
    ncwrap(nc_put_att_float(ncid, zVar, SERGHEI_NC_MISSING_VALUE, SERGHEI_NC_REAL, 1, &missing_value));
#endif

#if SERGHEI_WRITE_SUBDOMS
    dimids[0] = yDim;
    dimids[1] = xDim;
    realArr subdom = realArr("subdom", dom.nCellMem);
    Kokkos::parallel_for("subdom", dom.nCellMem, KOKKOS_LAMBDA(int iGlob) { subdom(iGlob) = dom.id; });
    Kokkos::fence();
    writeNetCDFfield(dom, ncid, subdomVar, st, ct, state.isnodata, subdom, data);
#endif

    writeRTNETCDFSW(rtsw, state, dom, par);

    ncwrap(ncmpi_close(ncid), __LINE__, par.myrank);
  }

  void outputNETCDFTransportSW(const RTStateSW &rtsw, const State &state, Domain const &dom, Parallel const &par, std::string dir)
  {
    std::string filename;
    filename = dir + "output_transportSW.nc";
    // Create the file
    ncwrap(ncmpi_open(MPI_COMM_WORLD, filename.c_str(), NC_WRITE, MPI_INFO_NULL, &ncid), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "csw", &cSWVar), __LINE__, par.myrank);

    writeRTNETCDFSW(rtsw, state, dom, par);

    ncwrap(ncmpi_close(ncid), __LINE__, par.myrank);
  }

  void writeRTNETCDFSW(const RTStateSW &rtsw, const State &state, Domain const &dom, Parallel const &par)
  {
    realArr data = realArr("data", dom.nCell);
    realArr2 data2 = realArr2("data2", rtsw.n_mass, dom.nCell);
    MPI_Offset st[4], ct[4];
    double timeIter[numOut + 1];
    double n_massIter[rtsw.n_mass];

    for (int i = 0; i < numOut + 1; i++)
    {
      timeIter[i] = i * 1.0 + dom.startTime;
    }

    st[0] = 0;
    ct[0] = numOut + 1;
    ncwrap(ncmpi_put_vara_double_all(ncid, tVar, st, ct, timeIter), __LINE__, par.myrank);

    for (int i = 0; i < rtsw.n_mass; i++)
    {
      n_massIter[i] = i * 1.0;
    }
    st[0] = 0;
    ct[0] = rtsw.n_mass;
    ncwrap(ncmpi_put_vara_double_all(ncid, nVar, st, ct, n_massIter), __LINE__, par.myrank);

    st[0] = numOut;
    st[1] = 0;
    st[2] = par.j_beg;
    st[3] = par.i_beg;
    ct[0] = 1;
    ct[1] = rtsw.n_mass;
    ct[2] = dom.ny;
    ct[3] = dom.nx;

    Kokkos::parallel_for(dom.nCell, KOKKOS_LAMBDA(int iGlob) {
        int i, j;
        dom.unpackIndices(iGlob, j, i);
        int ii = (dom.hc + j) * (dom.nx + 2 * dom.hc) + dom.hc + i; 
        for (int s = 0; s < rtsw.n_mass; ++s)
        {
            if (state.isnodata(ii)) {
                data2(s, iGlob) = SERGHEI_NAN; 
            } else {
                data2(s, iGlob) = rtsw.c(s, ii, 1); 
            }
        } });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_double_all(ncid, cSWVar, st, ct, data2.data()), __LINE__, par.myrank);
  }
#endif

  // =========================================================================
  // [CODE 1] WOFOST 作物生长耦合模型 IO 代码区
  // =========================================================================
#if CROP_GROWTH_MODEL
  void outputInitNETCDFCrop(const Wofost72 &wofost, Domain const &dom, Parallel const &par, std::string dir)
  {
    int dimids[3];
    MPI_Offset st[3], ct[3];
    doubleArr xCoord = doubleArr("xCoord", dom.nx);
    doubleArr yCoord = doubleArr("yCoord", dom.ny);

    static char timeUnits[] = "seconds";
    std::string filename;
    std::string longname;
    std::string units;

    filename = dir + "crop.nc";

    ncwrap(ncmpi_create(MPI_COMM_WORLD, filename.c_str(), SERGHEI_NC_MODE, MPI_INFO_NULL, &ncid), __LINE__, par.myrank);
    ncwrap(ncmpi_def_dim(ncid, "t", (MPI_Offset)NC_UNLIMITED, &tDim), __LINE__, par.myrank);
    ncwrap(ncmpi_def_dim(ncid, "x", (MPI_Offset)dom.nx_glob, &xDim), __LINE__, par.myrank);
    ncwrap(ncmpi_def_dim(ncid, "y", (MPI_Offset)dom.ny_glob, &yDim), __LINE__, par.myrank);

    dimids[0] = tDim;
    ncwrap(ncmpi_def_var(ncid, "t", NC_DOUBLE, 1, dimids, &cropTimeVar), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, cropTimeVar, "units", strlen(timeUnits), timeUnits), __LINE__, par.myrank);

    dimids[0] = xDim;
    ncwrap(ncmpi_def_var(ncid, "x", NC_DOUBLE, 1, dimids, &xVar), __LINE__, par.myrank);

    dimids[0] = yDim;
    ncwrap(ncmpi_def_var(ncid, "y", NC_DOUBLE, 1, dimids, &yVar), __LINE__, par.myrank);

    dimids[0] = tDim;
    dimids[1] = yDim;
    dimids[2] = xDim;

    ncwrap(ncmpi_def_var(ncid, "DVS", SERGHEI_NC_REAL, 3, dimids, &dvsVar), __LINE__, par.myrank);
    longname.assign("Development stage");
    units.assign("-");
    ncwrap(ncmpi_put_att_text(ncid, dvsVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, dvsVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "LAI", SERGHEI_NC_REAL, 3, dimids, &laiVar), __LINE__, par.myrank);
    longname.assign("Leaf area index");
    units.assign("m^2/m^2");
    ncwrap(ncmpi_put_att_text(ncid, laiVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, laiVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "RD", SERGHEI_NC_REAL, 3, dimids, &rdVar), __LINE__, par.myrank);
    longname.assign("Rooting depth");
    units.assign("cm");
    ncwrap(ncmpi_put_att_text(ncid, rdVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, rdVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "TAGP", SERGHEI_NC_REAL, 3, dimids, &tagpVar), __LINE__, par.myrank);
    longname.assign("Total above ground biomass");
    units.assign("kg/ha");
    ncwrap(ncmpi_put_att_text(ncid, tagpVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, tagpVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "TWSO", SERGHEI_NC_REAL, 3, dimids, &wsoVar), __LINE__, par.myrank);
    longname.assign("Total storage organ biomass");
    units.assign("kg/ha");
    ncwrap(ncmpi_put_att_text(ncid, wsoVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, wsoVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "TWLV", SERGHEI_NC_REAL, 3, dimids, &wlvVar), __LINE__, par.myrank);
    longname.assign("Total leaf biomass");
    units.assign("kg/ha");
    ncwrap(ncmpi_put_att_text(ncid, wlvVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, wlvVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "TWST", SERGHEI_NC_REAL, 3, dimids, &wstVar), __LINE__, par.myrank);
    longname.assign("Total stem biomass");
    units.assign("kg/ha");
    ncwrap(ncmpi_put_att_text(ncid, wstVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, wstVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "TWRT", SERGHEI_NC_REAL, 3, dimids, &wrtVar), __LINE__, par.myrank);
    longname.assign("Total root biomass");
    units.assign("kg/ha");
    ncwrap(ncmpi_put_att_text(ncid, wrtVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, wrtVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "RFTRA", SERGHEI_NC_REAL, 3, dimids, &rftraVar), __LINE__, par.myrank);
    longname.assign("Transpiration reduction factor");
    units.assign("-");
    ncwrap(ncmpi_put_att_text(ncid, rftraVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, rftraVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    ncwrap(ncmpi_def_var(ncid, "RZSM", SERGHEI_NC_REAL, 3, dimids, &rzsmVar), __LINE__, par.myrank);
    longname.assign("Root zone soil moisture content");
    units.assign("cm3/cm3");
    ncwrap(ncmpi_put_att_text(ncid, rzsmVar, "long_name", longname.length(), longname.c_str()), __LINE__, par.myrank);
    ncwrap(ncmpi_put_att_text(ncid, rzsmVar, "units", units.length(), units.c_str()), __LINE__, par.myrank);

    static char title[] = "SERGHEI-WOFOST Crop Simulation Output";
    ncwrap(ncmpi_put_att(ncid, NC_GLOBAL, "title", NC_CHAR, strlen(title) + 1, title), __LINE__, par.myrank);

    ncwrap(ncmpi_enddef(ncid), __LINE__, par.myrank);

    Kokkos::parallel_for("compute_grid_coord_x", dom.nx, KOKKOS_LAMBDA(int i) { xCoord(i) = dom.xll + (par.i_beg + i + 0.5) * dom.dxConst; });
    Kokkos::parallel_for("compute_grid_coord_y", dom.ny, KOKKOS_LAMBDA(int j) { yCoord(j) = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j + 0.5) * dom.dxConst; });
    Kokkos::fence();

    st[0] = par.i_beg;
    ct[0] = dom.nx;
    ncwrap(ncmpi_put_vara_double_all(ncid, xVar, st, ct, xCoord.data()), __LINE__, par.myrank);
    st[0] = par.j_beg;
    ct[0] = dom.ny;
    ncwrap(ncmpi_put_vara_double_all(ncid, yVar, st, ct, yCoord.data()), __LINE__, par.myrank);

    writeCropNETCDF(wofost, dom, par);

    ncwrap(ncmpi_close(ncid), __LINE__, par.myrank);
  }

  void outputNETCDFCrop(const Wofost72 &wofost, Domain const &dom, Parallel const &par, std::string dir)
  {
    std::string filename = dir + "crop.nc";

    ncwrap(ncmpi_open(MPI_COMM_WORLD, filename.c_str(), NC_WRITE, MPI_INFO_NULL, &ncid), __LINE__, par.myrank);

    ncwrap(ncmpi_inq_varid(ncid, "t", &cropTimeVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "DVS", &dvsVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "LAI", &laiVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "RD", &rdVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "TAGP", &tagpVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "TWSO", &wsoVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "TWLV", &wlvVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "TWST", &wstVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "TWRT", &wrtVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "RFTRA", &rftraVar), __LINE__, par.myrank);
    ncwrap(ncmpi_inq_varid(ncid, "RZSM", &rzsmVar), __LINE__, par.myrank);

    writeCropNETCDF(wofost, dom, par);

    ncwrap(ncmpi_close(ncid), __LINE__, par.myrank);
  }

  void writeCropNETCDF(const Wofost72 &wofost, Domain const &dom, Parallel const &par)
  {
    ncArr data = ncArr("data", dom.nCell);
    MPI_Offset st[3], ct[3];

    std::vector<double> timeIter(numOutCrop + 1);

    for (int i = 0; i < numOutCrop + 1; i++)
    {
      timeIter[i] = i * 1.0 + dom.startTime;
    }

    st[0] = 0;
    ct[0] = numOutCrop + 1;
    ncwrap(ncmpi_put_vara_double_all(ncid, cropTimeVar, st, ct, timeIter.data()), __LINE__, par.myrank);

    st[0] = numOutCrop;
    st[1] = par.j_beg;
    st[2] = par.i_beg;
    ct[0] = 1;
    ct[1] = dom.ny;
    ct[2] = dom.nx;

    auto view_dvs = wofost.s.DVS;
    auto view_lai = wofost.s.LAI;
    auto view_rd = wofost.s.RD;
    auto view_tagp = wofost.TAGP;
    auto view_twso = wofost.sods.TWSO;
    auto view_twlv = wofost.lds.TWLV;
    auto view_twst = wofost.sds.TWST;
    auto view_twrt = wofost.rds.TWRT;
    auto view_rftra = wofost.ets.RFTRA;
    auto view_rzsm = wofost.ets.RZSM;

    Kokkos::parallel_for("write_crop_dvs", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_dvs(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, dvsVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_lai", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_lai(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, laiVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_rd", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_rd(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, rdVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_tagp", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_tagp(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, tagpVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_twso", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_twso(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, wsoVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_twlv", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_twlv(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, wlvVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_twst", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_twst(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, wstVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_twrt", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_twrt(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, wrtVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_rftra", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_rftra(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, rftraVar, st, ct, data.data()), __LINE__, par.myrank);

    Kokkos::parallel_for("write_crop_rzsm", dom.nCell, KOKKOS_LAMBDA(int iGlob) { data(iGlob) = view_rzsm(iGlob); });
    Kokkos::fence();
    ncwrap(ncmpi_put_vara_real_all(ncid, rzsmVar, st, ct, data.data()), __LINE__, par.myrank);
  }
#endif

  // parallel netcdf input functionality

  // generic error handling function
  static void handle_error(Parallel &par, int status, int lineno)
  {

    if (par.masterproc)
      std::cerr << RERROR << "Error at line: " << lineno << ": " << ncmpi_strerror(status) << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
};
#endif