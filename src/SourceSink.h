#ifndef _SOURCESINK_H_
#define _SOURCESINK_H_

#include "define.h"
#include "SArray.h"
#include "State.h"
#include "GwDomain.h"
#include "GwState.h"

// Atmospheric forcing support (NetCDF, wind, etc.)
#include "atmForcing.h"

// WOFOST crop growth model coupling
#if CROP_GROWTH_MODEL
#include "./cropsrc/Wofost72.h"
#endif

#define INF_NONE 0
#define INF_CONSTANT 1
#define INF_HORTON 2
#define INF_GREENAMPT 3

// Rainfall polygon support
#ifndef SERGHEI_RAINFALL_POLYGONS
#define SERGHEI_RAINFALL_POLYGONS 0
#endif

// ============================================================
// Rainfall polygon time interpolation functions
// ============================================================
#if SERGHEI_RAINFALL_POLYGONS
KOKKOS_INLINE_FUNCTION int findTimeBlockByPolygons(const realArr raintime, const int idx0, const int idx1, real const &t)
{
  int ti = idx0;
  if (t >= raintime(idx1))
  {
    ti = idx1;
  }
  else
  {
    for (int ii = idx0; ii < idx1; ii++)
    {
      if (t >= raintime(ii) && t < raintime(ii + 1))
        ti = ii;
    }
  }
  return (ti);
};

KOKKOS_INLINE_FUNCTION real interpolatePiecewiseByPolygons(const realArr rainvalue, int idx)
{
  return (rainvalue(idx));
};
#endif

// class State;	// forward declaration
class SourceSinkData; // forward declaration

// ============================================================
// [COMMON] Time interpolation utility functions
// ============================================================
KOKKOS_INLINE_FUNCTION void findTimeBlock(TimeSeries &ts, real const &t)
{
  // Optional forcing series may be empty.  Do not dereference an empty view.
  if (ts.np <= 0 || ts.time.data() == nullptr)
  {
    ts.timeIndex = 0;
    return;
  }
  if (ts.timeIndex < 0)
    ts.timeIndex = 0;
  if (ts.timeIndex >= ts.np)
    ts.timeIndex = ts.np - 1;
  if (t >= ts.time(ts.np - 1))
  {
    ts.timeIndex = ts.np - 1;
  }
  else
  {
    if (t >= ts.time(ts.timeIndex + 1))
      ts.timeIndex++;
  }
};

KOKKOS_INLINE_FUNCTION real interpolatePiecewise(TimeSeries const &ts, real const &t, int const spaceIndex)
{
  return (ts.value(ts.np * spaceIndex + ts.timeIndex));
};

KOKKOS_INLINE_FUNCTION real interpolateLinear(TimeSeries &ts, real const &t)
{
#if SERGHEI_DEBUG_WORKFLOW > 1
  std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
  int ii, jj;
  findTimeBlock(ts, t);
  ii = ts.timeIndex;
  jj = ii + 1;
  if (ii == ts.np - 1)
    return ts.value(ii);
  if (ts.time(jj) == ts.time(ii))
    return ts.value(ii);
  real v = ts.value(ii) + (ts.value(jj) - ts.value(ii)) / (ts.time(jj) - ts.time(ii)) * (t - ts.time(ii));
  return (v);
};

KOKKOS_INLINE_FUNCTION real interpolateValues(TimeSeries &ts, real const &t, int icol)
{
  int ii, jj;
  findTimeBlock(ts, t);
  ii = ts.timeIndex;
  jj = ii + 1;
  if (ii == ts.np - 1)
    jj = ii;
  real v = ts.values(ii, icol) + (ts.values(jj, icol) - ts.values(ii, icol)) / (ts.time(jj) - ts.time(ii)) * (t - ts.time(ii));
  return (v);
};

// ============================================================
//[COMMON] Constant infiltration helper class
// ============================================================
class ConstantInfiltration
{
private:
  real _constCap;

public:
  ConstantInfiltration(real constCap) : _constCap(constCap) {}
  real operator()()
  {
    return _constCap;
  }
};

// ============================================================
// InfiltrationModel class
// Base from CODE2 (spatial infiltration, color types, capacity array, checkModel)
// + backward compatibility with CODE1's assignModel interface
// ============================================================
class InfiltrationModel
{
public:
  // Spatial infiltration support
  static constexpr int spatialNone = 0;
  static constexpr int spatialRaster = 1;
  static constexpr int spatialSoilMap = 2;
  static constexpr int spatialLanduse = 3;
  static constexpr int spatialNetCDF = 4;
  int spatial = spatialRaster;
  MapClass soilmap; // soil map for infiltration

private:
  realArr infTime;
  KOKKOS_INLINE_FUNCTION real greenAmpt(const int ii, const real) const
  {
    real infCap = 0.;
    return (infCap);
  }

public:
  int model = -999;
  color nLabels = 0; // uses color (ushort) for memory efficiency
  realArr constCap;
  // Horton
  realArr k;
  realArr fc;
  realArr f0;
  // Green-Ampt
  real ks = -999;
  real psi = -999;
  real dtheta = -999;

  real infDry = 1E-8;

  // State variables for output and GreenAmpt model
  realArr infVol;    // accumulated infiltration volume
  realArr rate;      // infiltration rate, m/s
  realArr capacity;  // prescribed infiltration capacity, m/s
  colorArr infLabel; // uses colorArr for labels

  void allocate(const Domain &dom)
  {
    if (model)
    {
      rate = realArr("rate", dom.nCellMem);
      infVol = realArr("infVol", dom.nCellMem);
      if (model == INF_HORTON)
        infTime = realArr("infTime", dom.nCellMem);
      // Allocate capacity array for raster/NetCDF spatial infiltration
      if (spatial == spatialRaster || spatial == spatialNetCDF)
        capacity = realArr("capacity", dom.nCellMem);
    }
  }

  // Renamed from assignModel to checkModel with improved error handling
  int checkModel(Parallel &par)
  {
    int error = 0;
#if SERGHEI_DEBUG_INFILTRATION
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "model: " << model << std::endl;
#endif
    switch (model)
    {
    case INF_NONE:
      if (par.masterproc)
        std::cerr << BDASH << "No infiltration capacity" << std::endl;
      break;
    case INF_CONSTANT:
      if (par.masterproc)
        std::cerr << BDASH << "Constant infiltration capacity" << std::endl;
      for (int id = 1; id < nLabels; id++)
      {
        if (constCap(id) < 0)
        {
          if (par.masterproc)
            std::cerr << RERROR << "Infiltration rate not found for constant infiltration model" << std::endl;
          error++;
        }
      }
      break;
    case INF_HORTON:
      if (par.masterproc)
        std::cerr << BDASH << "Horton infiltration capacity" << std::endl;
      for (int ii = 0; ii < nLabels; ii++)
      {
        if (k(ii) < 0)
        {
          if (par.masterproc)
            std::cerr << RERROR << "Shape factor not found for Horton infiltration model" << std::endl;
          error++;
        }
        if (f0(ii) < 0)
        {
          if (par.masterproc)
            std::cerr << RERROR << "Initial infiltration capacity not found for Horton infiltration model" << std::endl;
          error++;
        }
        if (fc(ii) < 0)
        {
          if (par.masterproc)
            std::cerr << RERROR << "Asymptotic infiltration capacity not found for Horton infiltration model" << std::endl;
          error++;
        }
      }
      break;
    case INF_GREENAMPT:
      if (par.masterproc)
      {
        std::cerr << BDASH << "Green-Ampt infiltration capacity" << std::endl;
        std::cerr << RERROR << "Not enabled yet" << std::endl;
      }
      error++;
      if (ks < 0)
      {
        if (par.masterproc)
          std::cerr << RERROR << "Saturated hydraulic conductivity not found for Green-Ampt infiltration model" << std::endl;
        error++;
      }
      if (psi < 0)
      {
        if (par.masterproc)
          std::cerr << RERROR << "Average suction head not found for Green-Ampt infiltration model" << std::endl;
        error++;
      }
      if (dtheta < 0)
      {
        if (par.masterproc)
          std::cerr << RERROR << "Water content difference not found for Green-Ampt infiltration model" << std::endl;
        error++;
      }
      break;
    default:
      if (par.masterproc)
        std::cerr << RERROR << "Error processing data in infiltration.input using infiltration model " << model << "." << std::endl;
      error++;
      break;
    }
    if (error > 0)
      return 0;
    return 1;
  }

  // Backward-compatible alias
  inline int assignModel(Parallel &par) { return checkModel(par); }

  // ComputeInfiltrationCapacity
  // CODE2's enhanced version with spatial infiltration support + KOKKOS_CLASS_LAMBDA
  inline void ComputeInfiltrationCapacity(const Domain &dom)
  {
#if SERGHEI_DEBUG_WORKFLOW
    std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
#endif
    if (model)
    {
      switch (model)
      {
      case INF_CONSTANT:
        // Support raster/NetCDF-based spatial infiltration
        if (spatial == spatialNetCDF || spatial == spatialRaster)
        {
          Kokkos::deep_copy(rate, capacity);
        }
        else
        {
          Kokkos::parallel_for("inf_constant", dom.nCell, KOKKOS_CLASS_LAMBDA(int iGlob) {
                  int ii = dom.getIndex(iGlob);
                  int id = infLabel(ii);
                  rate(ii) = constCap(id); });
        }
        break;
      case INF_HORTON:
        // Uses KOKKOS_CLASS_LAMBDA for direct member access
        Kokkos::parallel_for("inf_horton", dom.nCell, KOKKOS_CLASS_LAMBDA(int iGlob) {
              int ii = dom.getIndex(iGlob);
              int id = infLabel(ii);
              real t = infTime(ii) + dom.dt;
              infTime(ii) = t;
              rate(ii) = fc(id) + (f0(id)-fc(id))*exp(-k(id) * t); });
        break;
      }
    }
  }
};

// ============================================================
// Surface water reactive transport source/sink class
// Handles fertilizer application (dissolution model) and other
// surface solute source/sink terms
// ============================================================
#if SERGHEI_SURFACE_TRANSPORT
class RTSwSS
{
public:
  std::string id;
  int rtsstype; // 3: areal mass flux application (broadcast fertilisation) [mg/m2/s]
  int ncellsIT = 0;
  intArr icells;

  std::vector<int> spec_sstype;
  std::vector<real> spec_ssval_const;
  std::vector<std::string> spec_ssfile;
  std::vector<int> has_ssfile;
  std::vector<TimeSeries> spec_ts;
  std::vector<realArr> undissolved_mass; // undissolved solid fertilizer mass on the surface [mg/m2]
  realArr dissolved_mass_accum;          // mass actually dissolved into the surface liquid phase this hydrodynamic step [mg]

  MPI_Comm comm;

  inline int find_icells_sw(const State &state, std::string &id_in, const Domain &dom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly)
  {
    int foundInSubdom = -1;
    std::vector<int> tmpicells;
    std::vector<int> subdomains;

    for (int iGlob = 0; iGlob < dom.nCell; iGlob++)
    {
      int i, j;
      dom.unpackIndices(iGlob, j, i);
      int ii = dom.getHaloExtension(i, j);

      if (!state.isnodata(ii))
      {
        real xCoord = dom.xll + (par.i_beg + i + 0.5) * dom.dxConst;
        real yCoord = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j + 0.5) * dom.dxConst;

        if (geometry::isInsidePoly(nPoly, xPoly, yPoly, xCoord, yCoord))
        {
          tmpicells.push_back(ii);
        }
      }
    }

    ncellsIT = int(tmpicells.size());
    if (ncellsIT > 0)
      foundInSubdom = par.myrank;

    int ncells_all;
    SERGHEI_MPI_CHECK(MPI_Allreduce(&ncellsIT, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    int *subdoms = (int *)malloc(par.nranks * sizeof(int));
    SERGHEI_MPI_CHECK(MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD));
#if SERGHEI_DEBUG_BOUNDARY
    if (par.masterproc)
      std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "ncellsIT " << ncells_all << std::endl;
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "RTSWSS subdomains: ";
    for (int i = 0; i < par.nranks; i++)
    {
      std::cout << GGD << " ";
      if (subdoms[i] == par.myrank)
        std::cout << RED;
      std::cout << subdoms[i] << "\t" << RESET;
    }
    std::cout << std::endl;
#endif

    for (int i = 0; i < par.nranks; i++)
    {
      if (subdoms[i] >= 0)
        subdomains.push_back(subdoms[i]);
    }
    free(subdoms);

    MPI_Group group, subgroup;
    SERGHEI_MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &group));
    SERGHEI_MPI_CHECK(MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup));
    SERGHEI_MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm));

    if (ncells_all > 0)
    {
      icells = intArr("icells_rtsw", ncellsIT);
      undissolved_mass.resize(spec_sstype.size());
      for (size_t i = 0; i < spec_sstype.size(); i++)
      {
        undissolved_mass[i] = realArr("undissolved_mass", ncellsIT);
        Kokkos::deep_copy(undissolved_mass[i], 0.0);
      }
      dissolved_mass_accum = realArr("dissolved_mass_accum", spec_sstype.size());
      Kokkos::deep_copy(dissolved_mass_accum, 0.0);
#ifdef __NVCC__
      cudaMemcpyAsync(icells.data(), tmpicells.data(), ncellsIT * sizeof(int), cudaMemcpyHostToDevice);
      cudaDeviceSynchronize();
#else
      std::memcpy(icells.data(), tmpicells.data(), ncellsIT * sizeof(int));
#endif
    }
    else
    {
      if (par.masterproc)
        std::cerr << YEXC << "SWRTM: No internal cells found for SW RT Source/Sink id '" << id_in << "'" << std::endl;
      return 0;
    }
    return 1;
  }

  // Forward declaration of core computation function (implemented after SourceSink class)
  template <typename RTStateSWType>
  static void ComputeRTSWSourceSink(RTStateSWType &rtsw, const State &state, const Domain &dom, const class SourceSinkData &ss, int iSpec, int step_part);
};
#endif // SERGHEI_SURFACE_TRANSPORT

// ============================================================
// SourceSinkData class
// CODE2's enhanced infrastructure (initializeSW, wind, NetCDF, polygon rain, granular timer)
// + CODE1's surface reactive transport source/sink (RTSwSS) integration
// + CODE1's custom evaporation computation
// ============================================================
class SourceSinkData
{
public:
// NetCDF atmospheric forcing support
#if SERGHEI_NETCDF_FORCING
  AtmosphericForcing aforcing;
  realArr rainAccum;
#endif

  TimeSeries rain, evap, wind; // added wind
  // Surface ET partition forcing. In crop-coupled runs, ET0 is combined
  // with the current WOFOST LAI; evap remains the static-LAI fallback.
  TimeSeries et0;
  real etExtinctionCoefficient = 0.5;
  bool hasEtPartitionForcing = false;
#if CROP_GROWTH_MODEL
  bool invalidDynamicLaiWarningIssued = false;
#endif
  InfiltrationModel inf;
  realArr rainRate, evapRate, windspd, winddir; // added wind arrays

  // Actual surface water evaporation volume for SW-GW coupling
  real evapFluxActual;

  // Rainfall polygon support
  intArr rainPol;

#if SERGHEI_RAINFALL_POLYGONS
  int nrainpol = 0;
  int nrainpoints = 0;
  std::vector<TimeSeries> rainSeries;

  intArr startIdxRain;
  intArr npRain;
  intArr timeIdx;
  realArr raintime;
  realArr rainvalue;
#endif

  // Surface reactive transport source/sink
#if SERGHEI_SURFACE_TRANSPORT
  std::vector<RTSwSS> rtswss;
#endif

  // Initialization function
  // CODE2's comprehensive initializeSW with NetCDF, polygon, wind support
  void initializeSW(Domain const &dom, State const &state, Parallel const &par)
  {
    if (dom.isRain)
    {
#if SERGHEI_NETCDF_FORCING
      rainAccum = realArr("rainAccum", dom.nCellMem);
#else
      rainRate = realArr("rainRate", dom.nCellMem);
      Kokkos::deep_copy(rainRate, 0.0);
#if SERGHEI_RAINFALL_POLYGONS
      rainPol = intArr("rainPol", dom.nCellMem);
      Kokkos::deep_copy(rainPol, -1);
      for (int pp = 0; pp < nrainpol; pp++)
      {
        TimeSeries &rain = rainSeries[pp];
        find_raincells(dom, state, par, rain, pp);
      }
      startIdxRain = intArr("startIdxRain", nrainpol);
      npRain = intArr("npRain", nrainpol);
      timeIdx = intArr("timeIdx", nrainpol);
      Kokkos::parallel_for(nrainpol, KOKKOS_CLASS_LAMBDA(int iGlob) {
              startIdxRain(iGlob) = 0;
              npRain(iGlob) = 0;
              timeIdx(iGlob) = 0; });
      raintime = realArr("raintime", nrainpoints);
      rainvalue = realArr("rainvalue", nrainpoints);
      Kokkos::parallel_for(nrainpoints, KOKKOS_CLASS_LAMBDA(int iGlob) {
              raintime(iGlob) = 0.0;
              rainvalue(iGlob) = 0.0; });
      int startposition = 0;
      for (int pp = 0; pp < nrainpol; pp++)
      {
        TimeSeries &rain = rainSeries[pp];
        store_raintime(dom, rain, startposition, pp);
      }
#endif
#endif
    }
    if (dom.isEvap)
    {
      evapRate = realArr("evapRate", dom.nCellMem);
    }
    if (dom.isWind)
    {
      windspd = realArr("windspd", dom.nCellMem);
      winddir = realArr("winddir", dom.nCellMem);
    }
    if (inf.model)
    {
      inf.allocate(dom);
    }
  }

  // Backward-compatible alias
  void allocateSW(Domain const &dom)
  {
    if (dom.isRain)
    {
      rainRate = realArr("rainRate", dom.nCellMem);
    }
    if (dom.isEvap)
    {
      evapRate = realArr("evapRate", dom.nCellMem);
    }
    if (inf.model)
    {
      inf.allocate(dom);
    }
  }

// Rainfall polygon cell finding
#if SERGHEI_RAINFALL_POLYGONS
  inline int find_raincells(Domain const &dom, State const &state, Parallel const &par, TimeSeries &rain, int pp)
  {
    int nPoly = rain.nver;
    realArr &xPoly = rain.xPoly;
    realArr &yPoly = rain.yPoly;
    int count = 0;
    for (int iGlob = 0; iGlob < dom.nCell; iGlob++)
    {
      int i, j;
      dom.unpackIndices(iGlob, j, i);
      int ii = dom.getHaloExtension(i, j);
      if (!state.isnodata(ii))
      {
        real xCoord = dom.xll + (par.i_beg + i + 0.5) * dom.dxConst;
        real yCoord = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j + 0.5) * dom.dxConst;
        if (geometry::isInsidePoly(nPoly, xPoly, yPoly, xCoord, yCoord))
        {
          rainPol(ii) = pp;
          count++;
        }
      }
    }
    return (1);
  }

  inline int store_raintime(Domain const &dom, TimeSeries const &rain, int &startposition, int pp)
  {
    startIdxRain(pp) = startposition;
    npRain(pp) = rain.np;
    timeIdx(pp) = startposition;
    int idx;
    for (int ii = 0; ii < rain.np; ii++)
    {
      idx = startposition + ii;
      raintime(idx) = rain.time(ii);
      rainvalue(idx) = rain.value(ii);
    }
    startposition += rain.np;
    return (1);
  }

  inline void ComputeRainByPolygons(const State &state, const Domain &dom)
  {
    if (dom.isRain)
    {
      Kokkos::parallel_for(nrainpol, KOKKOS_CLASS_LAMBDA(int pp) {
        int idx0=timeIdx(pp);
        int idx1=startIdxRain(pp)+npRain(pp)-1;
        timeIdx(pp)=findTimeBlockByPolygons(raintime,idx0,idx1,dom.etime); });
      Kokkos::parallel_for("rain_interpolation", dom.nCell, KOKKOS_CLASS_LAMBDA(int iGlob) {
        int ii = dom.getIndex(iGlob);
        real Irain = 0.0;
        if(!state.isnodata(ii)){
          int pp=rainPol(ii);
          if(pp>=0)
            Irain = interpolatePiecewiseByPolygons(rainvalue, timeIdx(pp));
        }
        rainRate(ii) = Irain; });
    }
  }
#endif // SERGHEI_RAINFALL_POLYGONS

  // ComputeRain - CODE2's version with CODE1's alternative indexing comment
  inline void ComputeRain(const Domain &dom)
  {
    if (dom.isRain)
    {
      int rainx = rain.nx;
      int rainy = rain.ny;
      int intervalx = dom.nx / rainx;
      int intervaly = dom.ny / rainy;

      realArr &rr_p = rainRate;
      findTimeBlock(rain, dom.etime);
      TimeSeries rrain = rain;

      Kokkos::parallel_for("rain_interpolation", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
        int ix;
        int iy;
        dom.unpackIndices(iGlob, iy, ix);
        int ii = dom.getHaloExtension(ix, iy);
        int _x = ix / intervalx;
        int _y = iy / intervaly;
        int rain_glob = _x + _y * rainx;
        real rainValue = interpolatePiecewise(rrain, dom.etime, rain_glob);
        rr_p(ii) = rainValue;
#if SERGHEI_DEBUG_RAINFALL
        Kokkos::printf("%s[DEBUG] %s%s _x: %lf _y: %lf ix: %d iy: %d, ~> rainfall: %lf \train_glob: %lf\n", GRAY, __PRETTY_FUNCTION__, RESET, _x, _y, ix, iy, rr_p(iGlob), rain_glob);
#endif
      });
    }
  }

  // ComputeEvap
  // Uses CODE2's spatial interpolation (interpolatePiecewise) as default
  // CODE1's uniform linear interpolation available as alternative (commented)
  inline void ComputeEvap(const Domain &dom)
  {
    if (dom.isEvap && evap.np > 0)
    {
      realArr &rr_e = evapRate;
      findTimeBlock(evap, dom.etime);
      TimeSeries revap = evap;
      // Spatially varying evaporation via piecewise interpolation
      Kokkos::parallel_for("evap_interpolation", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
              int ix;
              int iy;
              dom.unpackIndices (iGlob, iy, ix);
              int ii = dom.getHaloExtension(ix,iy);
              int evap_glob = 0;
              real evapValue = interpolatePiecewise(revap, dom.etime, evap_glob);
              rr_e(ii) = evapValue; });
    }
  }

#if CROP_GROWTH_MODEL
  // Compute potential surface-water evaporation from ET0 and the current
  // spatial WOFOST LAI. The actual water-limited flux is still imposed by
  // TimeIntegrator::evaporation().
  inline void ComputeEvap(const Domain &dom, const realArr &dynamicLai)
  {
    if (!dom.isEvap)
      return;

    // Preserve the legacy/static forcing path when this run has no ET
    // partition data supplied by et.input, or when the crop field cannot be
    // mapped one-to-one to the local physical surface cells. The extent check
    // prevents an out-of-bounds device access if crop initialization failed.
    if (!hasEtPartitionForcing || et0.np <= 0)
    {
      ComputeEvap(dom);
      return;
    }
    if (dynamicLai.data() == nullptr || dynamicLai.extent(0) != dom.nCell)
    {
      if (!invalidDynamicLaiWarningIssued)
      {
        std::cerr << RERROR
                  << "Surface evaporation cannot use dynamic WOFOST LAI: "
                  << "expected " << dom.nCell << " cells, received "
                  << dynamicLai.extent(0)
                  << ". Falling back to the static-LAI et.input series."
                  << std::endl;
        invalidDynamicLaiWarningIssued = true;
      }
      ComputeEvap(dom);
      return;
    }

    real currentEt0 = interpolateLinear(et0, dom.etime);
    real extinctionCoefficient = etExtinctionCoefficient;
    realArr &rr_e = evapRate;

    Kokkos::parallel_for("evap_dynamic_wofost_lai", dom.nCell,
      KOKKOS_LAMBDA(int iGlob) {
        int ix, iy;
        dom.unpackIndices(iGlob, iy, ix);
        int ii = dom.getHaloExtension(ix, iy);
        real lai = dynamicLai(iGlob);
        // The negated comparison also catches NaN without relying on a
        // backend-specific isfinite implementation.
        if (!(lai >= 0.0))
          lai = 0.0;
        rr_e(ii) = currentEt0 * exp(-extinctionCoefficient * lai);
      });
  }
#endif

  // Wind computation
  inline void ComputeWind(const Domain &dom)
  {
    if (dom.isWind)
    {
      realArr &rr_w = windspd;
      realArr &rr_d = winddir;
      findTimeBlock(wind, dom.etime);
      TimeSeries rwind = wind;
      Kokkos::parallel_for("wind_interpolation", dom.nCell, KOKKOS_LAMBDA(int idom) {
              int ix, iy;
              dom.unpackIndices (idom, iy, ix);
              int ii = dom.getHaloExtension(ix,iy);
              real spdValue = interpolatePiecewise(rwind, dom.etime, 0);
              rr_w(ii) = spdValue;
              real dirValue = interpolatePiecewise(rwind, dom.etime, 1);
              rr_d(ii) = dirValue; });
    }
  }

  // ComputeSWSourceSink
  // CODE2's structure: Kokkos::fence(), granular timer, polygon/NetCDF rain, wind
  // CODE1's timer approach integrated into CODE2's granular timer structure
  inline void ComputeSWSourceSink(const State &state, const Domain &dom)
  {
    Kokkos::fence(); // Ensure previous GPU operations complete
    Kokkos::Timer timer;
#if SERGHEI_DEBUG_WORKFLOW
    std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
#endif

// Multiple rainfall computation paths
#if SERGHEI_RAINFALL_POLYGONS
    ComputeRainByPolygons(state, dom);
#elif SERGHEI_NETCDF_FORCING
    aforcing.computeRain(dom, state, rainRate);
#else
    ComputeRain(dom);
#endif
    ComputeEvap(dom);
    ComputeWind(dom); // Wind computation
    inf.ComputeInfiltrationCapacity(dom);

    Kokkos::fence(); // Ensure all GPU operations complete before timing
    // Granular timer: uses swe.ss.raininf for detailed timing breakdown
    dom.timers.swe.ss.raininf += timer.seconds();
  }

#if CROP_GROWTH_MODEL
  inline void ComputeSWSourceSink(const State &state, const Domain &dom,
                                  const realArr &dynamicLai)
  {
    Kokkos::fence();
    Kokkos::Timer timer;
#if SERGHEI_DEBUG_WORKFLOW
    std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
#endif

#if SERGHEI_RAINFALL_POLYGONS
    ComputeRainByPolygons(state, dom);
#elif SERGHEI_NETCDF_FORCING
    aforcing.computeRain(dom, state, rainRate);
#else
    ComputeRain(dom);
#endif
    ComputeEvap(dom, dynamicLai);
    ComputeWind(dom);
    inf.ComputeInfiltrationCapacity(dom);

    Kokkos::fence();
    dom.timers.swe.ss.raininf += timer.seconds();
  }
#endif

}; // end class SourceSinkData

// ============================================================
// Subsurface source/sink class (GwSS)
// #if SERGHEI_RE_MODEL (renamed from CODE1's SERGHEI_SUBSURFACE_MODEL)
//
// CODE2's base structure + CODE1's advanced features:
//   - WOFOST crop growth coupling
//   - Beer's Law ET partitioning (potential transpiration/evaporation)
//   - Feddes water stress function
//   - Dynamic root density distribution (Vogeland model with numerical integration)
//   - Tile drainage
//   - SW-GW evaporation-transpiration coupling
//   - Root zone cell finding and allocation
//   - ET0 / LAI / RD time series
// ============================================================
#if SERGHEI_RE_MODEL
class GwSS
{
public:
// subsurface ss directions
#define XPLUS 1
#define XMINUS 2
#define YPLUS 3
#define YMINUS 4
#define ZPLUS 5
#define ZMINUS 6

  // The type of source/sink
  //  0 : Evapotranspiration (from PM equation)
  //  1 : Flux (constant or time series)
  //  2 : Head for drainage (constant or time series)
  //  3 : Tile drainage

  std::string id; // identifier string
  int sstype;
  int direction;
  int ndepth;
  int ncellsIT = 0;

  // Root zone cells
  int ncel_ROOT = 0; // number of root zone cells
  int k_max_root;    // maximum root depth layer index

  intArr icells;
  intArr root_icells; // array of indexes of root zone cells
  real wc_root_zone;  // average water content in root zone
  realArr ssvals, ssdata;

  real Qinflow, Qoutflow;
  // Detailed outflow breakdown
  real Qoutflow_SoilEvap, Qoutflow_RootTransp, Qoutflow_TileDrainage;
  real Cpipe;

  TimeSeries ts;
  // ET-related time series
  TimeSeries et0_series, evap, tran;
  TimeSeries lai_series, rd_series;

  // Water stress and root distribution parameters
  real f;              // extinction coefficient for Beer's Law
  real h1, h2, h3, h4; // Feddes water stress thresholds
  real px, py, pz, xs, ys, zs, xm, ym;
  realArr coef_wat, coef_root;
  realArr wc_root; // root zone water content array

  // Tile drain parameters
  real drain_de;         // effective diameter [L]
  real drain_red_factor; // reduction factor [-]
  real drain_Cd = 1.0;   // dimensionless correction factor [-]

  // WOFOST coupling parameters
  bool use_realtime_data = false;
  realArr wofost_lai; // WOFOST LAI per cell
  realArr wofost_rd;  // WOFOST root depth per cell

  void allocateWofostCoupling(int nCells)
  {
    if (use_realtime_data)
    {
      // WOFOST initialization may already have allocated halo-sized arrays.
      // Preserve that layout; only allocate the compact physical-cell form
      // when no coupling storage exists yet.
      if (wofost_lai.data() == nullptr)
        wofost_lai = realArr("wofost_lai_ss", nCells);
      if (wofost_rd.data() == nullptr)
        wofost_rd = realArr("wofost_rd_ss", nCells);
      Kokkos::deep_copy(wofost_lai, 0.0);
      Kokkos::deep_copy(wofost_rd, 0.0);
    }
  }

  MPI_Comm comm;

  void allocateGW(GwDomain const &gdom)
  {
    ssdata = realArr("ssdata", gdom.nCellMem);
    for (int idx = 0; idx < gdom.nCellMem; idx++)
    {
      ssdata(idx) = 0.0;
    }
    if (sstype == 0)
    {
      coef_wat = realArr("wat", ncellsIT);
      coef_root = realArr("root", ncellsIT);
    }
  }

  // Allocate root zone coefficient arrays
  void allocateRootCoef(GwDomain const &gdom)
  {
    if (sstype == 0)
    {
      wc_root = realArr("wc_root", ncel_ROOT);
    }
  }

  // Find root zone internal cells
  inline int find_root_icells(GwState &gw, std::string &id, GwDomain &gdom, Parallel &par, real zm)
  {
    int foundInSubdom;
    std::vector<int> tmpicells;
    std::vector<int> subdomains;

    int kmax = 0, kmin = gdom.nz;
    for (int kk = 0; kk < gdom.nz; kk++)
    {
      for (int jj = 0; jj < gdom.ny; jj++)
      {
        for (int ii = 0; ii < gdom.nx; ii++)
        {
          int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
          foundInSubdom = -1;
          real zCoord = -gdom.depth(iGlob);

          if (zCoord <= zm)
          {
            tmpicells.push_back(iGlob);
            if (kk > kmax)
              kmax = kk;
            if (kk < kmin)
              kmin = kk;
          }
        }
      }
    }

    if (kmax > kmin)
      ndepth = kmax - kmin;
    else
      ndepth = 1;

    int ncel_ROOT = int(tmpicells.size());
    if (ncel_ROOT > 0)
      foundInSubdom = par.myrank;

    int ncel_ROOT_all;
    SERGHEI_MPI_CHECK(MPI_Allreduce(&ncel_ROOT, &ncel_ROOT_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    int *subdoms = (int *)malloc(par.nranks * sizeof(int));
    SERGHEI_MPI_CHECK(MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD));
    for (int i = 0; i < par.nranks; i++)
    {
      if (subdoms[i] >= 0)
        subdomains.push_back(subdoms[i]);
    }
    free(subdoms);

    MPI_Group group, subgroup;
    SERGHEI_MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &group));
    SERGHEI_MPI_CHECK(MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup));
    SERGHEI_MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm));

    if (ncel_ROOT_all > 0)
    {
      root_icells = intArr("root_icells", ncel_ROOT);
#ifdef __NVCC__
      cudaMemcpyAsync(root_icells.data(), tmpicells.data(), ncel_ROOT * sizeof(int), cudaMemcpyHostToDevice);
      cudaDeviceSynchronize();
#else
      std::memcpy(root_icells.data(), tmpicells.data(), ncel_ROOT * sizeof(int));
#endif
    }
    else
    {
      if (par.masterproc)
      {
        std::cerr << RERROR << "No root zone cells found for subsurface boundary with id '" << id << "'" << std::endl;
      }
      return 0;
    }
    return 1;
  }

  // Find maximum root depth layer index
  inline int find_max_root_depth_layer(GwDomain &gdom, real zm)
  {
    int k_max_root = -1;
    for (int kk = 0; kk < gdom.nz; kk++)
    {
      for (int jj = 0; jj < gdom.ny; jj++)
      {
        for (int ii = 0; ii < gdom.nx; ii++)
        {
          int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
          real zCoord = -1 * gdom.z(iGlob);
          if (zCoord <= zm)
            k_max_root = kk;
        }
      }
    }
    if (k_max_root == -1)
    {
      std::cerr << RERROR << "No cells found in the root zone for zm = " << zm << std::endl;
    }
    return k_max_root;
  }

  // Find internal cells for applying source/sink conditions
  // CODE2's structure with CODE1's depth-based z coordinate calculation
  inline int find_icells(GwState &gw, std::string &id, GwDomain &gdom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly, realArr &zPoly)
  {
    int foundInSubdom;
    std::vector<int> tmpicells;
    std::vector<int> tmpgcells;
    std::vector<int> subdomains;
    int kmax = 0, kmin = gdom.nz, idx = 0;
    for (int kk = 0; kk < gdom.nz; kk++)
    {
      for (int jj = 0; jj < gdom.ny; jj++)
      {
        for (int ii = 0; ii < gdom.nx; ii++)
        {
          int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
          foundInSubdom = -1;
          real xCoord = gdom.xll + (par.i_beg + ii + 0.5) * gdom.dx;
          real yCoord = gdom.yll + gdom.ny_glob * gdom.dx - (par.j_beg + jj + 0.5) * gdom.dx;
          // Use depth-based z coordinate (positive downward, negated)
          real zCoord = -gdom.depth(iGlob);
          if (geometry::isInsidePoly3D(nPoly, xPoly, yPoly, zPoly, xCoord, yCoord, zCoord))
          {
            tmpicells.push_back(iGlob);
            if (kk > kmax)
            {
              kmax = kk;
            }
            if (kk < kmin)
            {
              kmin = kk;
            }
          }
        }
      }
    }
    if (kmax > kmin)
    {
      ndepth = kmax - kmin;
    }
    else
    {
      ndepth = 1;
    }

    ncellsIT = int(tmpicells.size());
    if (ncellsIT > 0)
      foundInSubdom = par.myrank;

    int ncells_all;
    SERGHEI_MPI_CHECK(MPI_Allreduce(&ncellsIT, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    int *subdoms;
    subdoms = (int *)malloc(par.nranks * sizeof(int));
    SERGHEI_MPI_CHECK(MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD));
    for (int i = 0; i < par.nranks; i++)
    {
      if (subdoms[i] >= 0)
      {
        subdomains.push_back(subdoms[i]);
      }
    }
    free(subdoms);
    MPI_Group group, subgroup;
    SERGHEI_MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &group));
    SERGHEI_MPI_CHECK(MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup));
    SERGHEI_MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm));
    if (ncells_all > 0)
    {
      icells = intArr("icells", ncellsIT);
#ifdef __NVCC__
      cudaMemcpyAsync(icells.data(), tmpicells.data(), ncellsIT * sizeof(int), cudaMemcpyHostToDevice);
      cudaDeviceSynchronize();
#else
      std::memcpy(icells.data(), tmpicells.data(), ncellsIT * sizeof(int));
#endif
    }
    else
    {
      if (par.masterproc)
      {
        std::cerr << RERROR << "No internal cells found for subsurface boundary with id '" << id << "'" << std::endl;
      }
      return 0;
    }
    return 1;
  }

  // =========================================================================
  // Extracts the 3D interior cell block delimited by 8 vertices (4 top, 4 bottom).
  // Used for tile drains, deep local injection and other selections located in
  // the interior of the study area
  // =========================================================================
  inline int find_icells_tile_drainage(GwState &gw, std::string &id, GwDomain &gdom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly, realArr &zPoly)
  {
    if (nPoly != 8)
    {
      if (par.masterproc)
        std::cerr << RERROR << "find_icells_tile_drainage requires exactly 8 points (4 for top, 4 for bottom)!" << std::endl;
      return 0;
    }

    int foundInSubdom = -1;
    std::vector<int> tmpicells;
    std::vector<int> subdomains;

    int kmax = 0, kmin = gdom.nz;
    for (int kk = 0; kk < gdom.nz; kk++)
    {
      for (int jj = 0; jj < gdom.ny; jj++)
      {
        for (int ii = 0; ii < gdom.nx; ii++)
        {
          int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
          real xCoord = gdom.xll + (par.i_beg + ii + 0.5) * gdom.dx;
          real yCoord = gdom.yll + gdom.ny_glob * gdom.dx - (par.j_beg + jj + 0.5) * gdom.dx;
          real zCoord = -gdom.depth(iGlob);

          // Step 1: XY 2D polygon test (using first 4 vertices)
          int count = 0;
          for (int i = 0; i < 4; i++)
          {
            real x1 = xPoly(i), y1 = yPoly(i);
            real x2 = xPoly((i + 1) % 4), y2 = yPoly((i + 1) % 4);
            if ((yCoord > fmin(y1, y2)) && (yCoord <= fmax(y1, y2)) && (xCoord <= fmax(x1, x2)))
            {
              if (y1 != y2)
              {
                double xInter = x1 + (yCoord - y1) * (x2 - x1) / (y2 - y1);
                if (x1 == x2 || xCoord <= xInter)
                  count++;
              }
            }
          }

          if (count % 2 != 0)
          {
            // Step 2: IDW interpolation for top/bottom elevations
            double sum_z_top = 0.0, sum_w_top = 0.0;
            double sum_z_bot = 0.0, sum_w_bot = 0.0;
            double z_top = 0.0, z_bottom = 0.0;
            int exact_match_top = 0, exact_match_bot = 0;

            // 1. Top face elevation (vertices 0-3)
            for (int i = 0; i < 4; i++)
            {
              double dx = xCoord - xPoly(i);
              double dy = yCoord - yPoly(i);
              double dist2 = dx * dx + dy * dy;
              if (dist2 < 1e-12)
              {
                z_top = zPoly(i);
                exact_match_top = 1;
                break;
              }
              double w = 1.0 / dist2;
              sum_w_top += w;
              sum_z_top += zPoly(i) * w;
            }
            if (!exact_match_top)
              z_top = sum_z_top / sum_w_top;

            // 2. Bottom face elevation (vertices 4-7)
            for (int i = 4; i < 8; i++)
            {
              double dx = xCoord - xPoly(i);
              double dy = yCoord - yPoly(i);
              double dist2 = dx * dx + dy * dy;
              if (dist2 < 1e-12)
              {
                z_bottom = zPoly(i);
                exact_match_bot = 1;
                break;
              }
              double w = 1.0 / dist2;
              sum_w_bot += w;
              sum_z_bot += zPoly(i) * w;
            }
            if (!exact_match_bot)
              z_bottom = sum_z_bot / sum_w_bot;

            // Ensure the Z ordering (deeper elevations are smaller, e.g. -0.3 < -0.1, so z_top should exceed z_bottom)
            if (z_bottom > z_top)
            {
              double temp = z_top;
              z_top = z_bottom;
              z_bottom = temp;
            }

            // Step 3: 3D elevation test: the cell Z must lie between Top and Bottom
            if (zCoord <= z_top && zCoord >= z_bottom)
            {
              tmpicells.push_back(iGlob);
              if (kk > kmax)
                kmax = kk;
              if (kk < kmin)
                kmin = kk;
            }
          }
        }
      }
    }

    if (kmax > kmin)
      ndepth = kmax - kmin;
    else
      ndepth = 1;

    ncellsIT = int(tmpicells.size());
    if (ncellsIT > 0)
      foundInSubdom = par.myrank;

    int ncells_all;
    SERGHEI_MPI_CHECK(MPI_Allreduce(&ncellsIT, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    int *subdoms = (int *)malloc(par.nranks * sizeof(int));
    SERGHEI_MPI_CHECK(MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD));
    for (int i = 0; i < par.nranks; i++)
    {
      if (subdoms[i] >= 0)
        subdomains.push_back(subdoms[i]);
    }
    free(subdoms);

    MPI_Group group, subgroup;
    SERGHEI_MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &group));
    SERGHEI_MPI_CHECK(MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup));
    SERGHEI_MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm));

    if (ncells_all > 0)
    {
      icells = intArr("icells", ncellsIT);
#ifdef __NVCC__
      cudaMemcpyAsync(icells.data(), tmpicells.data(), ncellsIT * sizeof(int), cudaMemcpyHostToDevice);
      cudaDeviceSynchronize();
#else
      std::memcpy(icells.data(), tmpicells.data(), ncellsIT * sizeof(int));
#endif
    }
    else
    {
      if (par.masterproc)
        std::cerr << YEXC << "No internal cells found for tile drainage block with id '" << id << "'" << std::endl;
      return 0;
    }
    return 1;
  }

  // ================================================================
  // Beer's Law ET partitioning
  // Splits ET0 into potential transpiration and potential soil evaporation
  // ================================================================
  KOKKOS_INLINE_FUNCTION
  void compute_potential_fluxes(real et0, real lai, real k_ext,
                                real &pot_transp, real &pot_evap) const
  {
    real fraction = exp(-k_ext * lai);
    pot_transp = et0 * (1.0 - fraction);
    pot_evap = et0 * fraction;
  }

  // Numerical integration of root distribution function (trapezoidal rule)
  KOKKOS_INLINE_FUNCTION
  double integrateCRoot(double zs, double zm, double pz) const
  {
    if (zm < 1e-10)
      return 1e-10;
    constexpr int N = 1000;
    const double dz = zm / (N - 1);
    double integral = 0.0;
    for (int i = 0; i < N; ++i)
    {
      const double z = i * dz;
      const double expo = (pz / zm) * Kokkos::fabs(zs - z);
      const double c_root = (1.0 - z / zm) * Kokkos::exp(-expo);
      integral += c_root * dz;
    }
    if (integral < 1e-10)
      return 1e-10;
    return integral;
  }

  // Raw root density shape. The final allocation is normalized per column
  // using each layer's actual overlap with the root zone.
  KOKKOS_INLINE_FUNCTION
  real compute_root_density_raw(real z_depth, real root_depth) const
  {
    real zs = 0.2 * root_depth; // 20% of root depth
    real pz = 1.0;
    if (root_depth < 1e-6)
      return 0.0;
    else if (z_depth <= root_depth)
    {
      real expo = (pz / root_depth) * myfabs(zs - z_depth);
      real c_root = (1.0 - z_depth / root_depth) * exp(-expo);
      return fmax(0.0, c_root);
    }
    else
    {
      return 0.0;
    }
  }

  // Dimensionless fraction of a column's root uptake assigned to one layer.
  // This remains normalized for nonuniform dz and arbitrary vertical grids.
  KOKKOS_INLINE_FUNCTION
  real compute_column_root_fraction(int ii, int jj, int kk, real root_depth,
                                    GwDomain const &gdom) const
  {
    if (root_depth <= 1e-6)
      return 0.0;

    real denominator = 0.0;
    real numerator = 0.0;
    for (int kk_col = gdom.hc; kk_col < gdom.nz + gdom.hc; ++kk_col)
    {
      const int iGlob_col = kk_col * gdom.nxhc * gdom.nyhc +
                            jj * gdom.nxhc + ii;
      const real dz = gdom.dz(iGlob_col);
      const real center = gdom.depth(iGlob_col);
      const real top = center - 0.5 * dz;
      const real bottom = center + 0.5 * dz;
      const real overlap = fmin(root_depth, bottom) - fmax(0.0, top);
      if (overlap <= 0.0)
      {
        if (top >= root_depth)
          break;
        continue;
      }

      const real weight = compute_root_density_raw(center, root_depth) * overlap;
      denominator += weight;
      if (kk_col == kk)
        numerator = weight;
    }

    if (denominator <= 1e-30)
      return 0.0;
    return numerator / denominator;
  }

  // Feddes water stress reduction function
  KOKKOS_INLINE_FUNCTION
  real compute_root_water_uptake_declining(real h_val) const
  {
    real alpha_w = 0.0;
    if (h_val > h4 && h_val < h1)
    {
      if (h_val >= h2)
        alpha_w = (h_val - h1) / (h2 - h1);
      else if (h_val > h3)
        alpha_w = 1.0;
      else
        alpha_w = (h_val - h4) / (h3 - h4);
    }
    return alpha_w;
  }

  // Helper for potential ET from LAI
  KOKKOS_INLINE_FUNCTION
  void compute_potential_et_from_lai(real lai, real et0, real &pot_transp, real &pot_evap, real k_ext = 0.5) const
  {
    real exp_term = exp(-k_ext * lai);
    pot_transp = et0 * (1.0 - exp_term);
    pot_evap = et0 * exp_term;
  }

  // ================================================================
  // applyMatSS - Matrix assembly for source/sink
  // CODE1's advanced ET computation (Beer's Law, Feddes, root density, WOFOST coupling)
  // + CODE2's sstype==1 and sstype==2 handling with gdom.hc
  // ================================================================
#if CROP_GROWTH_MODEL
  inline void applyMatSS(GwState &gw, GwDomain &gdom, Wofost72 &wofost)
#else
  inline void applyMatSS(GwState &gw, GwDomain &gdom)
#endif
  {
    if (ncellsIT > 0)
    {
      // ET (Penman-Monteith)
      if (sstype == 0)
      {
        // Advanced ET computation with ET0, LAI, RD time series
        real current_et0 = interpolateLinear(et0_series, gdom.etime);
        real ts_lai = 0.0, ts_rd = 0.0;
#if CROP_GROWTH_MODEL
        // WOFOST coupled mode: LAI and RD come from WOFOST per-cell data
#else
        {
          ts_lai = interpolateLinear(lai_series, gdom.etime);
          ts_rd = interpolateLinear(rd_series, gdom.etime);
        }
#endif
        // Capture member variables for lambda
        real f_local = f;
        real h1_local = h1, h2_local = h2, h3_local = h3, h4_local = h4;

        Kokkos::parallel_for("gw_et_mat", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
            int ii, jj, kk;
            int iGlob = icells[idx];
            gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
            int idom = (kk - gdom.hc) * gdom.nx * gdom.ny + (jj - gdom.hc) * gdom.nx + ii - gdom.hc;

            // Index computation for SW-GW coupling
            int iGlobSW_Halo = jj * gdom.nxhc + ii;
            int iGlobSW_Internal = (jj - gdom.hc) * gdom.nx + (ii - gdom.hc);
            bool is_internal_cell = (ii >= gdom.hc && ii < gdom.nx + gdom.hc &&
                                     jj >= gdom.hc && jj < gdom.ny + gdom.hc);

            // Get LAI and RD (WOFOST or time series)
            real local_lai = 0.0;
            real local_rd  = 0.0;
#if CROP_GROWTH_MODEL
            if (is_internal_cell) {
              local_lai = wofost_lai(iGlobSW_Internal);
              local_rd  = wofost_rd(iGlobSW_Internal) * 0.01; // cm -> m
            }
#else
            local_lai = ts_lai;
            local_rd  = ts_rd;
#endif
            // Beer's Law: split ET0 into T_pot and E_pot
            real T_pot, E_pot;
            compute_potential_fluxes(current_et0, local_lai, f_local, T_pot, E_pot);
            // Feddes water stress coefficient
            real h_val   = gw.h(iGlob, 1);
            real alpha_w = compute_root_water_uptake_declining(h_val);

            // Root density distribution
            real beta = compute_column_root_fraction(ii, jj, kk, local_rd, gdom);

            // Transpiration sink term
            real transp_sink = 1.0 * T_pot * alpha_w * beta;

            // Update RHS
            gw.coef(idom, 7) += (gdom.dt * transp_sink / gdom.dz(iGlob));
            gw.transpGW(iGlob) = transp_sink;

            // Soil evaporation (top layer only)
            if (kk == gdom.hc)
            {
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
              // Coupled SW-GW evaporation model
              real qeTotalRate = E_pot;
              real swEvapRate = gw.surfaceWaterEvapActual(iGlobSW_Halo);
              real qeSoilWaterActual = qeTotalRate - swEvapRate;

              // Physical constraint correction
              if (qeTotalRate < 0.0){
                if (qeSoilWaterActual > 0.0) qeSoilWaterActual = 0.0;
                if (qeSoilWaterActual < qeTotalRate) qeSoilWaterActual = qeTotalRate;
              }
#else
              // Non-coupled mode: all potential evaporation acts on soil
              real qeSoilWaterActual = E_pot;
#endif
              gw.coef(idom, 7) += gdom.dt * qeSoilWaterActual / gdom.dz(iGlob);
              gw.evapGW(iGlob) = qeSoilWaterActual;

              // Upper boundary head limit
              if (gw.h(iGlob, 1) < -1000.0){
                gw.h(iGlob, 1) = -1000.0;
              }
            } });
      }
      // Flux Source/Sink (CODE2's structure with gdom.hc)
      else if (sstype == 1)
      {
        real qbc = interpolateLinear(ts, gdom.etime);
        Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                  int ii, jj, kk, idom, iGlob = icells[idx];
                  gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                  idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
                  gw.coef(idom,7) += gdom.dt * qbc; });
      }
      // Internal Drainage with Fixed Head
      else if (sstype == 2)
      {
        real hbc = ssvals[0];
        if (direction == XPLUS || direction == XMINUS)
        {
          Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                  int ii, jj, kk, idom, iGlob = icells[idx];
                  gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                  idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
                  int jp = idom+1, jm = idom-1, kp = idom+gdom.nx*gdom.ny, km = idom-gdom.nx*gdom.ny;
                  gw.coef(jp,7) -= Cpipe * gw.coef(idom,3) * hbc;
                  gw.coef(jm,7) -= Cpipe * gw.coef(idom,4) * hbc;
                  gw.coef(kp,7) -= Cpipe * gw.coef(idom,5) * hbc;
                  gw.coef(km,7) -= Cpipe * gw.coef(idom,6) * hbc;
                  gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                  gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                  gw.coef(jp,4) = 0.0;  gw.coef(jm,3) = 0.0;
                  gw.coef(kp,6) = 0.0;  gw.coef(km,5) = 0.0; });
        }
        else if (direction == YPLUS || direction == YMINUS)
        {
          Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                  int ii, jj, kk, idom, iGlob = icells[idx];
                  gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                  idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
                  int ip = idom+1, im = idom-1, kp = idom+gdom.nx*gdom.ny, km = idom-gdom.nx*gdom.ny;
                  gw.coef(ip,7) -= Cpipe * gw.coef(idom,1) * hbc;
                  gw.coef(im,7) -= Cpipe * gw.coef(idom,2) * hbc;
                  gw.coef(kp,7) -= Cpipe * gw.coef(idom,5) * hbc;
                  gw.coef(km,7) -= Cpipe * gw.coef(idom,6) * hbc;
                  gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                  gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                  gw.coef(ip,2) = 0.0;  gw.coef(im,1) = 0.0;
                  gw.coef(kp,6) = 0.0;  gw.coef(km,5) = 0.0; });
        }
        else
        {
          Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                  int ii, jj, kk, idom, iGlob = icells[idx];
                  gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                  idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
                  int ip = idom+1, im = idom-1, jp = idom+1, jm = idom-1;
                  gw.coef(jp,7) -= Cpipe * gw.coef(idom,3) * hbc;
                  gw.coef(jm,7) -= Cpipe * gw.coef(idom,4) * hbc;
                  gw.coef(ip,7) -= Cpipe * gw.coef(idom,1) * hbc;
                  gw.coef(im,7) -= Cpipe * gw.coef(idom,2) * hbc;
                  gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                  gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                  gw.coef(jp,4) = 0.0;  gw.coef(jm,3) = 0.0;
                  gw.coef(ip,2) = 0.0;  gw.coef(im,1) = 0.0; });
        }
      }
      // Tile Drains: implicit source/sink assembly
      else if (sstype == 3)
      {
        Kokkos::parallel_for("gw_ss_drain_mat", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
              int ii, jj, kk, iGlob = icells[idx];
              gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
              int idom = (kk - gdom.hc) * gdom.nx * gdom.ny + (jj - gdom.hc) * gdom.nx + ii - gdom.hc;
              real h_val = gw.h(iGlob, 1);
              if (h_val > 0.0){
                  int ivg = gw.soilID(iGlob) * NVG;
                  real ks = gw.vgTable(ivg);
                  real Gamma = 4.0 * drain_Cd * ks / (gdom.dx * gdom.dy);
                  gw.coef(idom, 0) += gdom.dt * Gamma;
              } });
      }
    }
  }

  // ================================================================
  // applyWCSS - Water content update for source/sink
  // CODE1's advanced ET computation with detailed outflow tracking
  // + CODE2's sstype==1 and sstype==2 handling
  // ================================================================
#if CROP_GROWTH_MODEL
  inline void applyWCSS(GwState &gw, GwDomain &gdom, Wofost72 &wofost)
#else
  inline void applyWCSS(GwState &gw, GwDomain &gdom)
#endif
  {
    if (ncellsIT > 0)
    {
      Qoutflow = 0.0;
      // Detailed outflow breakdown
      Qoutflow_SoilEvap = 0.0;
      Qoutflow_RootTransp = 0.0;
      Qoutflow_TileDrainage = 0.0;
      Qinflow = 0.0;

      // ET (Penman-Monteith)
      if (sstype == 0)
      {
        // ET0 from time series
        real current_et0 = interpolateLinear(et0_series, gdom.etime);
        real ts_lai = 0.0, ts_rd = 0.0;
#if CROP_GROWTH_MODEL
        // WOFOST coupled mode
#else
        {
          ts_lai = interpolateLinear(lai_series, gdom.etime);
          ts_rd = interpolateLinear(rd_series, gdom.etime);
        }
#endif
        // Capture for lambda
        real f_local = f;

        // parallel_reduce with dual accumulation (transpiration + soil evap)
        Kokkos::parallel_reduce("gw_et_wcss", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp, real &tmp_qe) {
                int ii, jj, kk;
                int iGlob = icells[idx];
                gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                int idom = (kk - gdom.hc) * gdom.nx * gdom.ny + (jj - gdom.hc) * gdom.nx + ii - gdom.hc;

                int iGlobSW_Halo = jj * gdom.nxhc + ii;
                int iGlobSW_Internal = (jj - gdom.hc) * gdom.nx + (ii - gdom.hc);
                bool is_internal_cell = (ii >= gdom.hc && ii < gdom.nx + gdom.hc &&
                                         jj >= gdom.hc && jj < gdom.ny + gdom.hc);

                real local_lai = 0.0;
                real local_rd = 0.0;
#if CROP_GROWTH_MODEL
                if (is_internal_cell){
                  local_lai = wofost_lai(iGlobSW_Internal);
                  local_rd = wofost_rd(iGlobSW_Internal) * 0.01; // cm -> m
                }
#else
                local_lai = ts_lai;
                local_rd = ts_rd;
#endif

                // Beer's Law + Feddes + Root Distribution
                real T_pot, E_pot;
                compute_potential_fluxes(current_et0, local_lai, f_local, T_pot, E_pot);
                // T_pot = current_et0;
                // E_pot = 0.0;

                real h_val = gw.h(iGlob, 1);
                real alpha_w = compute_root_water_uptake_declining(h_val);
                real beta = compute_column_root_fraction(ii, jj, kk, local_rd, gdom);

                real transp_sink = 1.0 * T_pot * alpha_w * beta;
                gw.wc(iGlob, 1) += gdom.dt * transp_sink / gdom.dz(iGlob);
                tmp += transp_sink * gdom.dx * gdom.dy;

                if (kk == gdom.hc){
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
                  real qeTotalRate = E_pot;
                  real swEvapRate = gw.surfaceWaterEvapActual(iGlobSW_Halo);
                  real qeSoilWaterActual = qeTotalRate - swEvapRate;
                  if (qeTotalRate < 0.0){
                    if (qeSoilWaterActual > 0.0) qeSoilWaterActual = 0.0;
                    if (qeSoilWaterActual < qeTotalRate) qeSoilWaterActual = qeTotalRate;
                  }
#else
                  real qeSoilWaterActual = E_pot;
#endif
                  gw.wc(iGlob, 1) += gdom.dt * qeSoilWaterActual / gdom.dz(iGlob);
                  tmp_qe += qeSoilWaterActual * gdom.dx * gdom.dy;

                  if (gw.h(iGlob, 1) < -1000) {
                    gw.h(iGlob, 1) = -1000;
                  }
                } }, Kokkos::Sum<real>(Qoutflow_RootTransp), Kokkos::Sum<real>(Qoutflow_SoilEvap));

#if CROP_GROWTH_MODEL
        if (wofost.ets.GW_TPOT_ACC.extent(0) > 0 && wofost.ets.GW_TACT_ACC.extent(0) > 0)
        {
          auto gw_tpot_acc = wofost.ets.GW_TPOT_ACC;
          auto gw_tact_acc = wofost.ets.GW_TACT_ACC;
          Kokkos::parallel_for("gw_wofost_feedback_accumulate", ncellsIT,
                               KOKKOS_CLASS_LAMBDA(int idx) {
            int ii, jj, kk;
            int iGlob = icells[idx];
            gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
            if (ii < gdom.hc || ii >= gdom.nx + gdom.hc ||
                jj < gdom.hc || jj >= gdom.ny + gdom.hc)
              return;

            int iSurf = (jj - gdom.hc) * gdom.nx + (ii - gdom.hc);
            real lai = wofost_lai(iSurf);
            real rd = wofost_rd(iSurf) * 0.01;
            real T_pot_local, E_pot_local;
            compute_potential_fluxes(current_et0, lai, f_local, T_pot_local, E_pot_local);

            real beta_local = compute_column_root_fraction(ii, jj, kk, rd, gdom);
            if (beta_local <= 0.0)
              return;

            real alpha_local = compute_root_water_uptake_declining(gw.h(iGlob, 1));
            real cell_area = gdom.dx * gdom.dy;
            real dt = gdom.dt;
            Kokkos::atomic_add(&gw_tpot_acc(iSurf), Kokkos::fabs(T_pot_local * beta_local) * cell_area * dt);
            Kokkos::atomic_add(&gw_tact_acc(iSurf), Kokkos::fabs(T_pot_local * alpha_local * beta_local) * cell_area * dt);
          });
          Kokkos::fence();
        }
#endif
      }
      // Flux Source/Sink
      else if (sstype == 1)
      {
        real qbc = interpolateLinear(ts, gdom.etime);
        Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp) {
                              int ii, jj, kk, idom, iGlob = icells[idx];
                              gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                              idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
                              gw.wc(iGlob,1) += gdom.dt * qbc;
                              if (qbc > 0)	{tmp += qbc * gdom.dt;}
                              else {tmp -= qbc * gdom.dt;} }, Kokkos::Sum<real>(Qinflow));
      }
      // Internal Drainage with Fixed Head
      else if (sstype == 2)
      {
        real hbc = ssvals[0], flux;
        if (direction == XPLUS || direction == XMINUS)
        {
          Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp) {
                              int ii, jj, kk, idom, iGlob = icells[idx];
                              gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                              idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
                              int jp = iGlob+gdom.nxhc, jm = iGlob-gdom.nxhc, kp = iGlob+gdom.nxhc*gdom.nyhc, km = iGlob-gdom.nxhc*gdom.nyhc;
                              if (hbc < gw.h(jp,1))	{
                                  tmp += Cpipe * gdom.dx * gdom.dz(iGlob) * gw.k(iGlob,1) * (hbc - gw.h(jp,1)) / gdom.dy;
                              }
                              if (hbc < gw.h(jm,1))	{
                                  tmp += Cpipe * gdom.dx * gdom.dz(iGlob) * gw.k(jm,1) * (hbc - gw.h(jm,1)) / gdom.dy;
                              }
                              if (hbc < gw.h(kp,1))	{
                                  tmp += Cpipe * gdom.dx * gdom.dy * gw.k(iGlob,2) * (hbc - gw.h(kp,1)) / gdom.dz(kp);
                              }
                              if (hbc < gw.h(km,1))	{
                                  tmp += Cpipe * gdom.dx * gdom.dy * gw.k(km,2) * (hbc - gw.h(km,1)) / gdom.dz(km);
                              } }, Kokkos::Sum<real>(flux));
        }
        else if (direction == YPLUS || direction == YMINUS)
        {
          Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp) {
                              int ii, jj, kk, idom, iGlob = icells[idx];
                              gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                              idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
                              int ip = iGlob+1, im = iGlob-1, kp = iGlob+gdom.nxhc*gdom.nyhc, km = iGlob-gdom.nxhc*gdom.nyhc;
                              if (hbc < gw.h(ip,1))	{
                                  tmp += Cpipe * gdom.dy * gdom.dz(iGlob) * gw.k(iGlob,0) * (hbc - gw.h(ip,1)) / gdom.dx;
                              }
                              if (hbc < gw.h(im,1))	{
                                  tmp += Cpipe * gdom.dy * gdom.dz(iGlob) * gw.k(im,0) * (hbc - gw.h(im,1)) / gdom.dx;
                              }
                              if (hbc < gw.h(kp,1))	{
                                  tmp += Cpipe * gdom.dx * gdom.dy * gw.k(iGlob,2) * (hbc - gw.h(kp,1)) / gdom.dz(kp);
                              }
                              if (hbc < gw.h(km,1))	{
                                  tmp += Cpipe * gdom.dx * gdom.dy * gw.k(km,2) * (hbc - gw.h(km,1)) / gdom.dz(km);
                              } }, Kokkos::Sum<real>(flux));
        }
        else
        {
          Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp) {
                              int ii, jj, kk, idom, iGlob = icells[idx];
                              gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                              idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
                              int ip = iGlob+1, im = iGlob-1, jp = iGlob+gdom.nxhc, jm = iGlob-gdom.nxhc;
                              if (hbc < gw.h(ip,1))	{
                                  tmp += Cpipe * gdom.dy * gdom.dz(iGlob) * gw.k(iGlob,0) * (hbc - gw.h(ip,1)) / gdom.dx;
                              }
                              if (hbc < gw.h(im,1))	{
                                  tmp += Cpipe * gdom.dy * gdom.dz(iGlob) * gw.k(im,0) * (hbc - gw.h(im,1)) / gdom.dx;
                              }
                              if (hbc < gw.h(jp,1))	{
                                  tmp += Cpipe * gdom.dx * gdom.dz(iGlob) * gw.k(iGlob,1) * (hbc - gw.h(jp,1)) / gdom.dy;
                              }
                              if (hbc < gw.h(jm,1))	{
                                  tmp += Cpipe * gdom.dx * gdom.dz(iGlob) * gw.k(jm,1) * (hbc - gw.h(jm,1)) / gdom.dy;
                              } }, Kokkos::Sum<real>(flux));
        }
        Qoutflow = -flux;
      }
      // Tile drains: water content update and outflow statistics
      else if (sstype == 3)
      {
        real flux = 0.0;
        Kokkos::parallel_reduce("gw_ss_drain_wc", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp) {
                  int ii, jj, kk, iGlob = icells[idx];
                  gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                  
                  real h_val = gw.h(iGlob, 1);
                  if (h_val > 0.0) {
                      int ivg = gw.soilID(iGlob) * NVG;
                      real ks = gw.vgTable(ivg);
                      real Gamma = 4.0 * drain_Cd * ks / (gdom.dx * gdom.dy);
                      real q_drain = Gamma * h_val; 
                      tmp += q_drain * gdom.dx * gdom.dy * gdom.dz(iGlob);
                  } }, Kokkos::Sum<real>(flux));

        Qoutflow_TileDrainage = flux;
      }
    }
  }
};
#endif // SERGHEI_RE_MODEL

// ==========================================
// Groundwater reactive transport source/sink class RTGwSS
// ==========================================
#if SERGHEI_RE_MODEL
#define DEFAULT_TSCF 0.8 // transpiration stream concentration factor (0.0 = pure-water transpiration, 1.0 = full solute uptake)

class SourceSink; // forward declaration, used by the computation functions below
#if SERGHEI_SUBSURFACE_TRANSPORT
class RTGwSS
{
public:
  std::string id;
  // 1: volumetric mass flux injection [M/L3/T]
  // 2: point-source concentration injection carried by the flow [M/L3]
  // 3: areal mass flux application (broadcast fertilisation) [M/L2/T],
  //    automatically distributed over the vertical grid thickness
  int rtsstype;
  int ncellsIT = 0;
  intArr icells;
  realArr col_depth; // total vertical thickness of each selected water column (H_col)

  std::vector<int> spec_sstype; // 0: constant, 1: linearly interpolated time series, 2: step-constant time series (instantaneous fertilisation)
  std::vector<real> spec_ssval_const;
  std::vector<std::string> spec_ssfile;
  std::vector<int> has_ssfile;
  std::vector<TimeSeries> spec_ts;
  std::vector<realArr> spec_ssvals;

  MPI_Comm comm;

  // Find the interior cells of the region and precompute the vertical thickness
  inline int find_icells(GwState &gw, std::string &id_in, GwDomain &gdom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly, realArr &zPoly)
  {
    int foundInSubdom = -1;
    std::vector<int> tmpicells;
    std::vector<int> subdomains;

    // Traverse the interior cells only (excluding halo layers)
    for (int kk = 0; kk < gdom.nz; kk++)
    {
      for (int jj = 0; jj < gdom.ny; jj++)
      {
        for (int ii = 0; ii < gdom.nx; ii++)
        {
          // Use the correct layer offset when computing the global index
          int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
          real xCoord = gdom.xll + (par.i_beg + ii + 0.5) * gdom.dx;
          real yCoord = gdom.yll + gdom.ny_glob * gdom.dx - (par.j_beg + jj + 0.5) * gdom.dx;
          real zCoord = -gdom.depth(iGlob);
          if (geometry::isInsidePoly3D(nPoly, xPoly, yPoly, zPoly, xCoord, yCoord, zCoord))
          {
            tmpicells.push_back(iGlob);
          }
        }
      }
    }

    ncellsIT = int(tmpicells.size());
    if (ncellsIT > 0)
      foundInSubdom = par.myrank;

    int ncells_all;
    SERGHEI_MPI_CHECK(MPI_Allreduce(&ncellsIT, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

    int *subdoms = (int *)malloc(par.nranks * sizeof(int));
    SERGHEI_MPI_CHECK(MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD));
    for (int i = 0; i < par.nranks; i++)
    {
      if (subdoms[i] >= 0)
        subdomains.push_back(subdoms[i]);
    }
    free(subdoms);

    MPI_Group group, subgroup;
    SERGHEI_MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &group));
    SERGHEI_MPI_CHECK(MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup));
    SERGHEI_MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm));

    if (ncells_all > 0)
    {
      icells = intArr("icells_rt", ncellsIT);

      auto icells_host = Kokkos::create_mirror_view(icells);
      for (int i = 0; i < ncellsIT; i++)
      {
        icells_host(i) = tmpicells[i];
      }
      Kokkos::deep_copy(icells, icells_host);

      // ========================================================
      // Mass-conservation core: precompute H_col, the total thickness of the
      // 2D water column over which the polygon-selected cells are distributed
      // ========================================================
      col_depth = realArr("col_depth", gdom.nCellSwMem);
      Kokkos::deep_copy(col_depth, 0.0);

      intArr icells_local = this->icells;
      auto dz_local = gdom.dz;
      int nx_h = gdom.nxhc;
      int ny_h = gdom.nyhc;
      auto col_depth_local = this->col_depth;

      int ncellmem_check = gdom.nCellMem;

      realArr dz_temp = realArr("dz_temp", ncellmem_check);
      Kokkos::deep_copy(dz_temp, dz_local);

      Kokkos::parallel_for("Calc_Col_Depth", ncellsIT, KOKKOS_LAMBDA(int idx) {
                int iGlob = icells_local(idx);
                if(iGlob >= 0 && iGlob < ncellmem_check) {
                  int rem = iGlob % (nx_h * ny_h);
                  int jj_h = rem / nx_h;
                  int ii_h = rem % nx_h;
                  int iGlobSW = jj_h * nx_h + ii_h;
                  if(iGlobSW >= 0 && iGlobSW < nx_h * ny_h) {
                    Kokkos::atomic_add(&col_depth_local(iGlobSW), dz_temp(iGlob));
                  }
                } });
      Kokkos::fence();
    }
    else
    {
      if (par.masterproc)
        std::cerr << RERROR << "No internal cells found for RT Source/Sink id '" << id_in << "'" << std::endl;
      return 0;
    }
    return 1;
  }

  template <typename RTStateGWType>
  static void ComputeRTGWSourceSink(RTStateGWType &rtgw, const GwState &gw, const GwDomain &gdom, const SourceSink &ss, int iSpec);
};
#endif
#endif

// ============================================================
// SourceSink class
// Merges CODE2 structure with CODE1's solute transport vectors
// ============================================================
class SourceSink
{
public:
  std::vector<std::string> id;
  SourceSinkData swss;
#if SERGHEI_RE_MODEL
  std::vector<GwSS> gwss;
#endif

#if SERGHEI_SUBSURFACE_TRANSPORT
  std::vector<RTGwSS> rtgwss;
#endif
};

// ==========================================
// Implementation of the RTGwSS core computation functions (must follow the SourceSink declaration)
// ==========================================
#if SERGHEI_SUBSURFACE_TRANSPORT
template <typename RTStateGWType>
inline void RTGwSS::ComputeRTGWSourceSink(RTStateGWType &rtgw, const GwState &gw, const GwDomain &gdom, const SourceSink &ss, int iSpec)
{
  // Read the species-specific TSCF coefficient from the species state object
  real tscf = rtgw.tscf(iSpec);

  // 1. Transpiration (root uptake) & soil evaporation
  Kokkos::parallel_for("GW_RT_Transpiration", gdom.nCell, KOKKOS_LAMBDA(int idom) {
        int ii, jj, kk, iGlob;
        gdom.unpackIndices(idom, kk, jj, ii);
        iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

        real transp_rate = gw.transpGW(iGlob); 
        real evap_rate = gw.evapGW(iGlob);
        real C_curr = rtgw.c(iSpec, iGlob, 0); 
        
        // 1.1 Root uptake. RTcoef stores the volumetric storage equation
        // wc*C, while transp_rate is an areal water flux [m/s]. Divide by
        // the cell thickness to convert the sink to a volumetric increment,
        // matching the groundwater water-balance discretisation.
        if (transp_rate < 0.0 && tscf > 0.0) {
            const real dz = gdom.dz(iGlob);
            if (dz > 0.0) {
                const real removed_storage_transp =
                    (-transp_rate) * C_curr * tscf * rtgw.dt / dz;
                rtgw.RTcoef(idom, 7) -= removed_storage_transp;
            }
        }

        // 1.2 Soil evaporation (assumed pure-water evaporation: salts stay in
        // the soil, so the concentration effect arises automatically)
        if (evap_rate < 0.0) {
            real removed_mass_evap = -1 * evap_rate * C_curr * 0.0 * rtgw.dt;
            rtgw.RTcoef(idom, 7) -= removed_mass_evap;
        } });

  // 2. Pure solute cell source/sink terms (broadcast fertilisation, deep injection, etc.)
  for (size_t k = 0; k < ss.rtgwss.size(); k++)
  {
    const RTGwSS &rt_ss_obj = ss.rtgwss[k];
    if (rt_ss_obj.ncellsIT <= 0)
      continue;

    // Obtain the injection concentration or mass flux rate at the current time
    real C_inject = 0.0;
    if (rt_ss_obj.spec_sstype[iSpec] == 0)
    {
      C_inject = rt_ss_obj.spec_ssval_const[iSpec];
    }
    else if (rt_ss_obj.spec_sstype[iSpec] == 1)
    {
      TimeSeries ts = rt_ss_obj.spec_ts[iSpec];
      real t = gdom.etime;
      int ii_ts = 0;
      while (ii_ts < ts.np - 1 && t >= ts.time(ii_ts + 1))
        ii_ts++;
      int jj_ts = (ii_ts == ts.np - 1) ? ii_ts : ii_ts + 1;
      if (ii_ts == jj_ts)
        C_inject = ts.value(ii_ts);
      else
        C_inject = ts.value(ii_ts) + (ts.value(jj_ts) - ts.value(ii_ts)) / (ts.time(jj_ts) - ts.time(ii_ts)) * (t - ts.time(ii_ts));
    }
    else if (rt_ss_obj.spec_sstype[iSpec] == 2)
    {
      TimeSeries ts = rt_ss_obj.spec_ts[iSpec];
      real t = gdom.etime;
      int ii_ts = 0;
      while (ii_ts < ts.np - 1 && t >= ts.time(ii_ts + 1))
        ii_ts++;
      C_inject = ts.value(ii_ts);
    }

    int rtsstype = rt_ss_obj.rtsstype;
    auto icells_local = rt_ss_obj.icells;
    auto col_depth_local = rt_ss_obj.col_depth;
    int ncellsIT_local = rt_ss_obj.ncellsIT;

    Kokkos::parallel_for("GW_RT_Custom_Wells", rt_ss_obj.ncellsIT, KOKKOS_LAMBDA(int idx) {
            int iGlob = icells_local(idx);
            int nx_h = gdom.nxhc; int ny_h = gdom.nyhc;
            int kk_h = iGlob / (nx_h * ny_h);
            int rem = iGlob % (nx_h * ny_h);
            int jj_h = rem / nx_h;  int ii_h = rem % nx_h;
            int ii = ii_h - gdom.hc; int jj = jj_h - gdom.hc; int kk = kk_h - gdom.hc;

            if (ii >= 0 && ii < gdom.nx && jj >= 0 && jj < gdom.ny && kk >= 0 && kk < gdom.nz) {
                int idom = kk * gdom.nx * gdom.ny + jj * gdom.nx + ii;
                
                if (rtsstype == 1) { 
                    rtgw.RTcoef(idom, 7) += C_inject * rtgw.dt / 1000.0;
                } 
                else if (rtsstype == 3) { 
                    int iGlobSW = jj_h * nx_h + ii_h;
                    real H_col = col_depth_local(iGlobSW);
                    if (H_col > 1e-6) {
                        rtgw.RTcoef(idom, 7) += (C_inject / H_col) * (rtgw.dt) / 1000.0;
                    }
                }
            } });
  }
}
#endif

// ==========================================
// Implementation of the RTSwSS core computation functions
// ==========================================
#if SERGHEI_SURFACE_TRANSPORT
template <typename RTStateSWType>
inline void RTSwSS::ComputeRTSWSourceSink(RTStateSWType &rtsw, const State &state, const Domain &dom, const SourceSinkData &ss, int iSpec, int step_part)
{
  for (size_t k = 0; k < ss.rtswss.size(); k++)
  {
    const RTSwSS &sw_ss_obj = ss.rtswss[k];
    if (sw_ss_obj.ncellsIT <= 0)
      continue;

    // 1. Extract the external source mass flux rate at the current time [mg/m2/s]
    real C_inject = 0.0;
    if (sw_ss_obj.spec_sstype[iSpec] == 0)
    {
      C_inject = sw_ss_obj.spec_ssval_const[iSpec];
    }
    else if (sw_ss_obj.spec_sstype[iSpec] == 1 || sw_ss_obj.spec_sstype[iSpec] == 2)
    {
      TimeSeries ts = sw_ss_obj.spec_ts[iSpec];
      real t = dom.etime;
      int ii_ts = 0;
      while (ii_ts < ts.np - 1 && t >= ts.time(ii_ts + 1))
        ii_ts++;

      if (sw_ss_obj.spec_sstype[iSpec] == 2)
      {
        C_inject = ts.value(ii_ts);
      }
      else
      {
        int jj_ts = (ii_ts == ts.np - 1) ? ii_ts : ii_ts + 1;
        if (ii_ts == jj_ts)
          C_inject = ts.value(ii_ts);
        else
          C_inject = ts.value(ii_ts) + (ts.value(jj_ts) - ts.value(ii_ts)) / (ts.time(jj_ts) - ts.time(ii_ts)) * (t - ts.time(ii_ts));
      }
    }

    // 2. Dissolution logic in parallel
    int rtsstype = sw_ss_obj.rtsstype;
    intArr icells_local = sw_ss_obj.icells;
    auto undissolved_local = sw_ss_obj.undissolved_mass[iSpec];
    auto dissolved_mass_accum_local = sw_ss_obj.dissolved_mass_accum;
    int ncells = sw_ss_obj.ncellsIT;
    real dt_local = rtsw.dt;
    real cell_area = dom.cellArea();

    // --- Core physical parameters ---
    const real k_diss = 1.0e-3;         // dissolution rate constant [1/s]
    const real h_min_dissolve = 1.0e-3; // minimum water depth for dissolution: 1 mm
    // Solubility limit [mg/L]:
    // Although urea is highly soluble, the saturation concentration is capped
    // at 500 mg/L for numerical stability and as a common field-application
    // upper bound. This prevents unbounded concentration spikes.
    const real Solubility_Limit = 500000.0;

    Kokkos::parallel_for("SW_RT_Safe_Dissolution", ncells, KOKKOS_LAMBDA(int idx) {
        int hydro_idx = icells_local(idx);
        int iGlob = rtsw.hydroToRtIndex(hydro_idx, dom.hc);
        real h_old = state.h4rtsw(hydro_idx, 1);
        // Match transport(): state.h contains the current-step SW-GW exchange,
        // while h4rtsw(..., 0) is the pre-exchange hydrodynamic snapshot.
        real h_new = state.h(hydro_idx);
        real h = (step_part == 1) ? 0.5 * (h_old + h_new) : h_new;

        if (rtsstype == 3) { 
            // A. Mass conservation: first add all new fertilizer mass to the solid pool [mg/m2]
            undissolved_local(idx) += C_inject * dt_local;

            real M_solid = undissolved_local(idx);

            // B. Dissolution condition: sufficient water depth and solid fertilizer present
            if (h > h_min_dissolve && M_solid > 1.0e-8) {
                real C_curr = rtsw.c(iSpec, iGlob, 1);
                
                // Stop dissolving once the concentration reaches the solubility limit
                if (C_curr < Solubility_Limit) {
                    
                    // 1. Kinetics-limited amount the solid wants to dissolve [mg/m2]
                    real M_kin = M_solid * (1.0 - exp(-k_diss * dt_local));

                    // 2. Solubility-limited amount the water column can host [mg/m2]
                    // capacity = (saturation concentration - current concentration) * depth * 1000
                    real M_sat = (Solubility_Limit - C_curr) * h * 1000.0;
                    if (M_sat < 0.0) M_sat = 0.0;

                    // 3. The actual mass dissolved this step is the minimum of the two
                    real M_to_dissolve = fmin(M_kin, M_sat);
                    
                    // C. Update the liquid-phase concentration and mass
                    real delta_C = M_to_dissolve / (h * 1000.0);
                    rtsw.c(iSpec, iGlob, 1) += delta_C;
                    
                    // Important: keep the mass variable in sync for the transport half-step
                    rtsw.solute_mass(iSpec, iGlob, 1) = rtsw.c(iSpec, iGlob, 1) * h;
                    
                    // D. Update the solid pool: subtract the dissolved mass
                    undissolved_local(idx) -= M_to_dissolve;

                    // Mass-conservation bookkeeping counts only mass actually dissolved into the liquid phase
                    Kokkos::atomic_add(&dissolved_mass_accum_local(iSpec),
                                       M_to_dissolve * cell_area);
                }
            }
            // If the water depth is insufficient (h < 1 mm), M_to_dissolve = 0 and all mass stays in the solid pool
        } });
  }
}
#endif

#endif
