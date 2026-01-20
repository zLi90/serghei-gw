#ifndef _SOURCESINK_H_
#define _SOURCESINK_H_

#include "define.h"
#include "SArray.h"
#include "State.h"

#include "GwDomain.h"
#include "GwState.h"

#define INF_NONE 0
#define INF_CONSTANT 1
#define INF_HORTON 2
#define INF_GREENAMPT 3

// class State;	// forward declaration

/*
  TimeSeries provides a construct/class to store time series.
 */
class TimeSeries
{

public:
  int np;     // number of points in time
  int nc;     // number of grid cells with different time series values
  int nx = 1; // number of partitions in x direction
  int ny = 1; // number of partitions in y direction

  realArr time;
  realArr value;
  realArr2 values;
  //! zzb添加
  realArr LAI_value;
  realArr zm_value;
  int timeIndex = 0;

  /*
    // WARNING valid only for piece-wise constant time data
    inline real interpolate (real const &t, int spaceIndex){
      if(t >= time (np - 1)){
        timeIndex = np - 1;
      }
      else{
         if (t >= time (timeIndex + 1)) timeIndex++;
      }
      return (value (np * spaceIndex + timeIndex));
    }
  */
  void initialise(int n)
  {
    np = n;
    time = realArr("time", np);
    value = realArr("value", np);
    //! zzb添加
    LAI_value = realArr("LAI_value", np);
    zm_value = realArr("zm_value", np);
  };
};

KOKKOS_INLINE_FUNCTION void findTimeBlock(TimeSeries &ts, real const &t)
{
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
    jj = ii;
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

class ConstantInfiltration
{

private:
  real _constCap;

public:
  ConstantInfiltration(real constCap) : _constCap(constCap) {} // constructor
  real operator()()
  {
    return _constCap;
  }
};

class InfiltrationModel
{

private:
  /*
    KOKKOS_INLINE_FUNCTION real horton(const int ii,const real dt) const
    {
      real t = infTime(ii) + dt;
      infTime(ii) = t;
      real infCap = fc + (f0-fc)*exp(-k * t);
      return(infCap);
    }
  */
  // TODO need to program GreenAmpt model
  KOKKOS_INLINE_FUNCTION real greenAmpt(const int ii, const real) const
  {
    real infCap = 0.;
    return (infCap);
  }

  realArr infTime;

public:
  int model = -999;
  int nLabels = 0;
  realArr constCap;
  // Horton
  realArr k;
  realArr fc;
  realArr f0;
  // Green-Ampt
  real ks = -999;
  real psi = -999;
  real dtheta = -999;

  real infDry = 1E-8; // [L] threshold to consider dry for infiltration purposes

  // This is not a state variable
  // which is why it is here and not in class State.
  // It is necessary for output
  // and because it is a variable in the GreenAmpt model
  realArr infVol;  // accumulated infiltration volume
  realArr rate;    // infiltration rate
  intArr infLabel; // labels for heterogeneous infiltration

  // this is a function pointer which allows to redirect
  // to the specific infiltration capacity function.
  // the goal is to avoid evaluating which model to use every time step
  // real (InfiltrationModel::*capacity)(const int ii, const real dt) const;

  // Define infiltration capacity models
  /*
  real constant(const int ii, const real t) const{
    return(constCap);
  }
  */

  void allocate(const Domain &dom)
  {
    if (model)
    {
      rate = realArr("rate", dom.nCellMem);
      infVol = realArr("infVol", dom.nCellMem);
      if (model == INF_HORTON)
        infTime = realArr("infTime", dom.nCellMem);
    }
  }

  int assignModel(Parallel &par)
  {
    int error = 0;

#if SERGHEI_DEBUG_INFILTRATION
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "model: " << model << std::endl;
#endif
    switch (model)
    {
    case INF_NONE:
      if (par.masterproc)
      {
        std::cerr << BDASH << "No infiltration capacity" << std::endl;
      }
      break;
    case INF_CONSTANT:
      // capacity = &InfiltrationModel::constant;
      std::cerr << BDASH << "Constant infiltration capacity" << std::endl;
      for (int id = 1; id < nLabels; id++)
      {
        if (constCap(id) < 0)
        {
          if (par.masterproc)
          {
            std::cerr << RERROR << "Infiltration rate not found for constant infiltration model" << std::endl;
          }
          error++;
        }
      }
      break;
    case INF_HORTON:
      // capacity = &InfiltrationModel::horton;
      if (par.masterproc)
      {
        std::cerr << BDASH << "Horton infiltration capacity" << std::endl;
      }
      for (int ii = 0; ii < nLabels; ii++)
      {
        // std::cerr << ii << "\t" << k(ii) << "\t" << fc(ii) << "\t" << f0(ii) << std::endl;
        if (k(ii) < 0)
        {
          if (par.masterproc)
          {
            std::cerr << RERROR << "Shape factor not found for Horton infiltration model" << std::endl;
          }
          error++;
        }
        if (f0(ii) < 0)
        {
          if (par.masterproc)
          {
            std::cerr << RERROR << "Initial infiltration capacity not found for Horton infiltration model" << std::endl;
          }
          error++;
        }
        if (fc(ii) < 0)
        {
          if (par.masterproc)
          {
            std::cerr << RERROR << "Asymptotic infiltration capacity not found for Horton infiltration model" << std::endl;
          }
          error++;
        }
      }
      break;
    case INF_GREENAMPT:
      //			capacity = &InfiltrationModel::greenAmpt;
      if (par.masterproc)
      {
        std::cerr << BDASH << "Green-Ampt infiltration capacity" << std::endl;
        std::cerr << RERROR << "Not enabled yet" << std::endl;
      }
      error++;
      if (ks < 0)
      {
        if (par.masterproc)
        {
          std::cerr << RERROR << "Saturated hydraulic conductivity not found for Green-Ampt infiltration model" << std::endl;
        }
        error++;
      }
      if (psi < 0)
      {
        if (par.masterproc)
        {
          std::cerr << RERROR << "Average suction head not found for Green-Ampt infiltration model" << std::endl;
        }
        error++;
      }
      if (dtheta < 0)
      {
        if (par.masterproc)
        {
          std::cerr << RERROR << "Water content difference not found for Green-Ampt infiltration model" << std::endl;
        }
        error++;
      }
      break;
    default:
      if (par.masterproc)
      {
        std::cerr << RERROR << "Error processing data in infiltration.input using infiltration model " << model << "." << std::endl;
      }
      error++;
      break;
    }
    if (error > 0)
      return 0;
    return 1;
  }

  inline void ComputeInfiltrationCapacity(const Domain &dom)
  {
#if SERGHEI_DEBUG_WORKFLOW
    std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
#endif
    if (model)
    {
      realArr &inf_p = rate;
      intArr infLabel = this->infLabel;
      realArr constCap = this->constCap;

      switch (model)
      {
      case INF_CONSTANT:
        Kokkos::parallel_for("inf_constant", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
                            int ii = dom.getIndex(iGlob);
                            int id = infLabel(ii);
                            inf_p(ii) = constCap(id); });
        break;
      case INF_HORTON:
        realArr fc = this->fc;
        realArr f0 = this->f0;
        realArr k = this->k;
        realArr &infTime_p = infTime;
        Kokkos::parallel_for("inf_horton", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
                            int ii = dom.getIndex(iGlob);
                            int id = infLabel(ii);
                            real t = infTime_p(ii) + dom.dt;
                            infTime_p(ii) = t;
                            inf_p(ii) = fc(id) + (f0(id)-fc(id))*exp(-k(id) * t); });
        break;
      }
    }
  }
};

class SourceSinkData
{

public:
  TimeSeries rain, evap;

  InfiltrationModel inf;
  realArr rainRate, evapRate;
  real evapFluxActual; // 实际地表水蒸发体积

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

  inline void ComputeRain(const Domain &dom)
  {
    if (dom.isRain)
    {

      // ----------------------------------------------------------------------
      // get global values to map to the correct rain subdomain
      // ----------------------------------------------------------------------
      // int nx = dom.nx_glob; // computational cell number in x direction
      // int ny = dom.ny_glob; // computational cell number in y direction

      int rainx = rain.nx; // rain subdomain number in x direction
      int rainy = rain.ny; // rain subdomain number in y direction

      int intervalx = dom.nx / rainx; // approximate number of cells in
                                      // a subdomain in x direction
      int intervaly = dom.ny / rainy; // approximate number of cells in
                                      // a subdomain in y direction
      // ----------------------------------------------------------------------

      realArr &rr_p = rainRate;

      findTimeBlock(rain, dom.etime);
      TimeSeries rrain = rain;

      Kokkos::parallel_for("rain_interpolation", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
        int ix;
        int iy;

        // dom.unpackIndices (iGlob, iy, ix);
        // int ii = dom.getHaloExtension(ix,iy);
        unpackIndicesUniformGrid(iGlob, dom.ny, dom.nx, iy, ix);
        int ii = (hc + iy) * (dom.nx + 2 * hc) + hc + ix;

        int _x = ix / intervalx;
        int _y = iy / intervaly;

        int rain_glob = _x + _y * rainx;

        real rainValue = interpolatePiecewise(rrain, dom.etime, rain_glob);

        // printf("rain_glob: %d, _x: %d, _y: %d, ix: %d, iy: %d, ii: %d, rainValue: %e\n", rain_glob, _x, _y, ix, iy, ii, rainValue);

        rr_p(ii) = rainValue;

#if SERGHEI_DEBUG_RAINFALL
        std::cerr << GGD "_x : " << _x << " _j: " << _y << " ix: " << ix << ", iy: " << iy << " ~> rainfall " << rr_p(iGlob) << std::endl;
        std::cerr << GGD "rain_glob " << rain_glob << std::endl;
#endif
      });

      /** basically the same as above but parallel for-ized. this
          works with MPI. we may think about a switch that uses this
          portion of code when compiled for CPU.

      Kokkos::parallel_for (dom.nCell, KOKKOS_LAMBDA (int iGlob)
                {

            int ix; // global x coordinate
            int iy; // global y coordinate
            unpackIndices (iGlob, ny, nx, iy, ix);

            int _x = ix / intervalx;
            int _y = iy / intervaly;

            int rain_glob = _x + _y * rainx;

            real rainValue = rain.interpolate (dom.etime, rain_glob);

            int ii = getIndex (iGlob, dom);
            rr_p (ii) = rainValue;

                });
      **/
    }
#if SERGHEI_DEBUG_RAINFALL
    std::cerr << GGD "-----------" << std::endl;
    ;
#endif
  }

  inline void ComputeEvap(const Domain &dom)
  {

    if (dom.isEvap)
    {

      // printf("1111111ComputeEvap dom.etime: %f\n", dom.etime);

      int intervalx = dom.nx; // approximate number of cells in
      int intervaly = dom.ny; // approximate number of cells in
      realArr &rr_e = evapRate;
      // findTimeBlock(evap,dom.etime);
      TimeSeries revap = evap;

      // for (int i = 0; i < 2; i++) {

      //     printf("i:%d, evap.time: %f, evap.value: %e\n", i, evap.time[i], evap.value[i]);
      //   }
      real evapValue = interpolateLinear(evap, dom.etime);
      Kokkos::parallel_for("evap_interpolation", dom.nCell, KOKKOS_LAMBDA(int iGlob) {
        int ix;
        int iy;
        dom.unpackIndices(iGlob, iy, ix);
        int ii = dom.getHaloExtension(ix, iy);
        int evap_glob = 0;
        //! zzb 20250926 修改，原始空间分布插值
        // real evapValue = interpolatePiecewise(revap, dom.etime, evap_glob);

        //! zzb 20240912 修改，采用gwss中线性插值
        // real evapValue = interpolateLinear(evap, dom.etime);
        // real evapValue =  0.00001;

        rr_e(ii) = evapValue;

        // printf("ii: %d, ix: %d, iy: %d, evapValue: %e\n", ii, ix, iy, evapValue);
      });
    }
  }

  inline void ComputeSWSourceSink(const State &state, const Domain &dom)
  {
    Kokkos::Timer timer;
#if SERGHEI_DEBUG_WORKFLOW
    std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
#endif
    ComputeRain(dom);
    ComputeEvap(dom);
    inf.ComputeInfiltrationCapacity(dom);
    // no rate correction is necessary here beacuse the rate correction is done in ComputeNewState, according to the new water depth
    // timerRainInf += timer.seconds();
    dom.timers.raininf += timer.seconds();
  }
};

// Subsurface source/sinks
#if SERGHEI_SUBSURFACE_MODEL
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

  int sstype;
  int direction;     // direction of source/sink, only needed for drainage ss
  int ndepth;        // number of cells in the vertical direction of the polygon
  int ncellsIT = 0;  // number of internal source/sink cells
  int ncel_ROOT = 0; //! zzb  number of root zone cells

  intArr icells;      // array of indexes of boundary cells
  intArr root_icells; //! zzb array of indexes of root zone cells
  int k_max_root;     //! zzb 最大根深所在层数
  real wc_root_zone;  //! zzb 根区的含水率均值
  realArr ssvals, ssdata;
  real Qinflow, Qoutflow, Qoutflow_qe, Cpipe;
  TimeSeries ts;
  TimeSeries evap, tran;
  // TimeSeries zm_series;//!zzb zm time series
  TimeSeries lai_series, zm_series;

  // Real-time crop parameters for dynamic coupling with WOFOST
  realArr lai_realtime;
  realArr rd_realtime;
  bool use_realtime_data = false;

  // water stress and root distribution function
  // real lai;
  real f;
  real h1, h2, h3, h4;
  real px, py, pz, xs, ys, zs, xm, ym;
  realArr coef_wat, coef_root;
  realArr wc_root;

  MPI_Comm comm; // communicator for ranks associated to the BC
                 //! zzb 20241008修改
                 // private:
                 //     double* coef_root_h = nullptr;
                 //! zzb 20241008修改

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

  //! zzb 分配根区
  void allocateRootCoef(GwDomain const &gdom)
  {
    if (sstype == 0)
    {
      wc_root = realArr("wc_root", ncel_ROOT);
    }
  }
  //! zzb 查找最大zm根区网格函数
  inline int find_root_icells(GwState &gw, std::string &id, GwDomain &gdom, Parallel &par, real zm)
  {
    int foundInSubdom;           // 用于跟踪哪些子域与根区相关
    std::vector<int> tmpicells;  // 根区内部单元的索引数组
    std::vector<int> subdomains; // 记录哪些子域与根区相关

    // 遍历整个域，查找根区的内部单元
    int kmax = 0, kmin = gdom.nz, idx = 0;
    for (int kk = 0; kk < gdom.nz; kk++)
    {
      for (int jj = 0; jj < gdom.ny; jj++)
      {
        for (int ii = 0; ii < gdom.nx; ii++)
        {
          int iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
          foundInSubdom = -1;
          // real zCoord = gdom.z(iGlob);
          real zCoord = -gdom.depth(iGlob);

          // 判断单元是否在根区（zCoord <= zm）
          if (zCoord <= zm)
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

    // 计算根区的深度
    if (kmax > kmin)
    {
      ndepth = kmax - kmin;
    }
    else
    {
      ndepth = 1;
    }

    // 计算根区的单元数量
    int ncel_ROOT = int(tmpicells.size());
    if (ncel_ROOT > 0)
    {
      foundInSubdom = par.myrank;
    } // 如果当前子域中有至少一个单元在根区，标记为找到

    // 全局统计根区的单元数量
    int ncel_ROOT_all;
    MPI_Allreduce(&ncel_ROOT, &ncel_ROOT_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // 收集包含根区单元的子域信息
    int *subdoms;
    subdoms = (int *)malloc(par.nranks * sizeof(int));
    MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD);
    for (int i = 0; i < par.nranks; i++)
    {
      if (subdoms[i] >= 0)
      {
        subdomains.push_back(subdoms[i]);
      }
    }

    // 创建 MPI 通信子组
    MPI_Group group, subgroup;
    MPI_Comm_group(MPI_COMM_WORLD, &group);
    MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup);
    MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm);

    // 如果没有找到根区单元，输出错误信息
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

  //! zzb 查找最大根深处对应层数
  inline int find_max_root_depth_layer(GwDomain &gdom, real zm)
  {
    int k_max_root = -1; // 初始化最大根深所在层数

    // 遍历整个域，查找最大根深所在的层数
    for (int kk = 0; kk < gdom.nz; kk++)
    {
      for (int jj = 0; jj < gdom.ny; jj++)
      {
        for (int ii = 0; ii < gdom.nx; ii++)
        {
          int iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
          real zCoord = -1 * gdom.z(iGlob);

          // 判断当前单元是否在根区（zCoord <= zm）
          if (zCoord <= zm)
          {
            // 更新最大根深所在层数
            k_max_root = kk;
          }
        }
      }
    }

    // 如果没有找到任何单元，返回 -1 表示错误
    if (k_max_root == -1)
    {
      std::cerr << RERROR << "No cells found in the root zone for zm = " << zm << std::endl;
    }

    return k_max_root;
  }

  // find internal cells for applying source/sink conditions
  inline int find_icells(GwState &gw, std::string &id, GwDomain &gdom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly, realArr &zPoly)
  {
    int foundInSubdom;           // to keep track of which subdomains are associated to this boundary
    std::vector<int> tmpicells;  // array of indexes of internal cells
    std::vector<int> tmpgcells;  // array of indexes of ghost cells
    std::vector<int> subdomains; // keeps track of which subdomains are associated to the BC
                                 // Loop over the entire domain to find internal source/sink cells
    int kmax = 0, kmin = gdom.nz, idx = 0;
    for (int kk = 0; kk < gdom.nz; kk++)
    {
      for (int jj = 0; jj < gdom.ny; jj++)
      {
        for (int ii = 0; ii < gdom.nx; ii++)
        {
          int iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
          foundInSubdom = -1;
          real xCoord = gdom.xll + (par.i_beg + ii + 0.5) * gdom.dx;
          real yCoord = gdom.yll + gdom.ny_glob * gdom.dx - (par.j_beg + jj + 0.5) * gdom.dx;
          // real zCoord = gdom.z(iGlob);
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
      foundInSubdom = par.myrank; // if at least one cell in this subdomain (rank) is in the BC, tag as found

    int ncells_all;
    MPI_Allreduce(&ncellsIT, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int *subdoms;
    subdoms = (int *)malloc(par.nranks * sizeof(int));
    MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD);
    for (int i = 0; i < par.nranks; i++)
    {
      if (subdoms[i] >= 0)
      {
        subdomains.push_back(subdoms[i]);
      }
    }
    MPI_Group group, subgroup;
    MPI_Comm_group(MPI_COMM_WORLD, &group);
    MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup);
    MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm);
    // we need the total internal cells detected by all subdomain to launch an error otherwise
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

  //! zzb 对c_root函数进行积分，数值积分法，梯形法则
  KOKKOS_INLINE_FUNCTION
  double integrateCRoot(double zs, double zm, double pz) const
  {

    constexpr int N = 1000; // 积分点数
    double integral = 0.0;
    const double dz = zm / (N - 1); // 自动计算步长

    // Kokkos 并行归约
    Kokkos::parallel_reduce(
        "integrateCRoot",
        N,
        KOKKOS_LAMBDA(int i, double &local_sum) {
          const double z = i * dz;
          const double expo = (pz / zm) * Kokkos::fabs(zs - z);
          const double c_root = (1.0 - z / zm) * Kokkos::exp(-expo);
          local_sum += c_root * dz;
        },
        integral);

    return integral;
  }

  // get coefficients for root water uptake declining and root distribution

  void rootCoef(const GwState &gw, const GwDomain &gdom)
  {

    //! zzb 根据根长时间序列插值计算当前时刻根长值
    real zm = interpolateLinear(zm_series, gdom.etime);

    real zs = 0.2 * zm; //! zzb 20240912修改,zs取zm的20%
    real pz = 1.0;      //! pz,zs经验系数,pz=1.0

    real c_root_integral = integrateCRoot(zs, zm, pz);

    // find_root_icells(gw, gdom, par, zm);

    //! zzb 对c_root函数进行积分，数值积分法，梯形法则
    // 定义积分步长
    // real dz = 0.001;

    // // 计算 c_root 的积分值（归一化因子）
    // real c_root_integral = 0.0;
    // for (real z = 0; z <= zm; z += dz) {
    //     real zs = 0.2*zm; // zs 取 zm 的 20%
    //     real pz = 1.0;  // pz 经验系数
    //     real expo = (pz / zm) * myfabs(zs - z);
    //     real c_root = (1.0 - z / zm) * exp(-expo);
    //     c_root_integral += c_root * dz;
    // }

    Kokkos::parallel_for("root", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
      real c_wat = 1.0, c_root = 0.0, expo, x, y, z;
      int iGlob = icells[idx];
      real h = gw.h(iGlob, 1);

      // get coef_wat
      //!(1)Feddes reduction function
      real h1 = 0, h2 = -0.01, h3 = -5, h4 = -160;//小麦fedds
      //real h1 = 0, h2 = -0.01, h3 = -3, h4 = -160;//grass fedds
      // real h1 = -0.01, h2 = -0.25, h3 = -8, h4 = -80;//!zzb 现在是代码中直接定值
     // real h1 = 1, h2 = 0.55, h3 = -2.5, h4 = -150; // 水稻fedds
      if (h <= h1 && h > h2)
      {
        c_wat = (h - h1) / (h2 - h1);
      }
      else if (h <= h3 && h > h4)
      {
        c_wat = (h - h4) / (h3 - h4);
      }
      else if (h > h1 || h <= h4)
      {
        c_wat = 0.0;
      }

      coef_wat(idx) = c_wat;

      // printf("idx:%d, h: %e, c_wat: %e\n", idx, h, c_wat);

      //!(2)S-shaped reduction function(ignore the effect of salinity on root water uptake)
      //! p is a fitting parameter(-),usually p=3
      //! 芦苇h50 = -24.56，杏树h50=-0.533
      // real h50 = -0.533, p1 = 3, h_50 = -65, p2 = 2.3;
      // c_wat = 1 / (1 + pow((h / h50), p1)) ;
      // c_wat = 1;

      // get coordinates x, y, z
      x = gdom.x(iGlob);
      y = gdom.y(iGlob);
      // z = gdom.z(iGlob);
      z = gdom.depth(iGlob);

      // get coef_root
      // Three-dimensional root water uptake declining function
      // expo = px/(xm*myfabs(xs-x)) + py/(ym*myfabs(ys-y)) + pz/(zm*myfabs(zs-z));
      // c_root = (1.0-x/xm)*(1.0-y/ym)*(1.0-z/zm)*exp(-expo);

      // real zs = 0.2*zm;//!zzb 20240912修改,zs取zm的20%
      // real pz = 1.0;//!pz,zs经验系数,pz=1.0

      // one-dimensional root water uptake declining function
      //  x*, y*, and z* are indicated as Depth of Maximum Intensity or Radius of Maximum Intensity;
      //  px, py, and pz are assumed to be equal to one for x> x*, y> y*, z> z*,

      // printf("idx:%d, z: %e, zm: %f\n", idx, z, zm);
      if (z <= zm)
      { //! zzb 将c_root计算范围限制在zm内

        expo = (pz / zm) * myfabs(zs - z);
        c_root = (1.0 - z / zm) * exp(-expo);

        coef_root(idx) = c_root / c_root_integral * (gdom.thickH / gdom.nz_glob);
        // printf("idx:%d, z: %e, c_root: %e, c_root_integral:%e, coef_root(idx): %e\n", idx, z, c_root, c_root_integral, coef_root(idx));
      }

      // real c_root_integral = integrateCRoot(zs, zm, pz);

      // coef_root(idx) = c_root / c_root_integral * (gdom.thickH / gdom.nz_glob);

      // printf("555555idx:%d, z: %e, c_root: %e, c_root_integral:%e, coef_root(idx): %e\n", idx, z, c_root, c_root_integral, coef_root(idx));

      // for (int i = 0; i < ncellsIT; ++i) {
      // printf("coef_root(idx): %f\n", coef_root(idx));

      // }
      // }

      //! 根系密度函数 b求和

      double sum = 0.0;
      for (int i = 0; i < ncellsIT; ++i)
      {
        sum += coef_root(i);
      }
      // std::cout << "ncellsIT \n" <<ncellsIT<< std::endl;
      // std::cout << "Sum of coef_root: \n" << sum << std::endl;
    });
  }

  // //!zzb 对c_root函数进行积分，数值积分法，梯形法则
  //   double integrateCRoot(double zs, double zm, double pz) const {

  //       real integral = 0.0;
  //       real dz = 0.001; // 积分步长

  //       for (real z = 0; z <= zm; z += dz) {
  //           real expo = (pz / zm) * myfabs(zs - z);
  //           real c_root = (1.0 - z / zm) * exp(-expo);
  //           integral += c_root * dz;
  //       }
  //       return integral;
  //   }

  // Helper function to compute potential transpiration and evaporation from LAI
  KOKKOS_INLINE_FUNCTION
  void compute_potential_et_from_lai(real lai, real et0, real& pot_transp, real& pot_evap, real k_ext = 0.5) const {
      real exp_term = exp(-k_ext * lai);
      pot_transp = et0 * (1.0 - exp_term);
      pot_evap = et0 * exp_term;
  }

  inline void applyMatSS(GwState &gw, GwDomain &gdom)
  {

    if (ncellsIT > 0)
    {
      // ET (Penman-Monteith)

      if (sstype == 0)
      {
        real qt, qe, zm;
        
        // Use real-time data if available, otherwise use time series interpolation
        if (use_realtime_data && lai_realtime.data() != nullptr && rd_realtime.data() != nullptr) {
            // Get ET0 from time series (base reference evapotranspiration)
            real et0 = interpolateLinear(tran, gdom.etime); // tran contains ET0
            
            // For each cell, compute qt and qe based on real-time LAI
            // Note: This requires a different approach since we need per-cell values
            // We'll handle this in the parallel loop below
            zm = 0.0; // Will be set per cell from rd_realtime
        } else {
            // Fallback to original time series interpolation
            qt = interpolateLinear(tran, gdom.etime);
            zm = interpolateLinear(zm_series, gdom.etime);
            qe = interpolateLinear(evap, gdom.etime);
        }

        // Handle real-time data vs time series interpolation
        if (use_realtime_data && lai_realtime.data() != nullptr && rd_realtime.data() != nullptr) {
            // For real-time coupling, we need to handle per-cell calculations in the parallel loop
            // Get base ET0 from time series
            real et0_base = interpolateLinear(tran, gdom.etime);
            
            // Update root coefficients using real-time root depth
            // We'll modify rootCoef to accept real-time RD if available
            rootCoef(gw, gdom);
            
            Kokkos::parallel_for("gw_et_realtime", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                idom = (kk - hc) * gdom.nx * gdom.ny + (jj - hc) * gdom.nx + ii - hc;
                int iGlobSW = jj * gdom.nxhc + ii;
                
                // Get real-time LAI and root depth for this surface cell
                real lai_val = lai_realtime(iGlobSW);
                real rd_val = rd_realtime(iGlobSW) * 0.01; // cm to m
                
                // Compute potential transpiration and evaporation from LAI
                real pot_transp, pot_evap;
                compute_potential_et_from_lai(lai_val, et0_base, pot_transp, pot_evap);
                
                // Apply transpiration with root distribution and water stress
                gw.coef(idom, 7) += gdom.dt * pot_transp * coef_wat(idx) * coef_root(idx) / gdom.dz(iGlob);
                gw.transpGW(iGlob) = pot_transp * coef_wat(idx) * coef_root(idx);
                
                // Handle evaporation for top layer
                if (kk == 1) {
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
                    real swEvapRate = gw.surfaceWaterEvapActual(iGlobSW);
                    real qeSoilWaterActual = pot_evap - swEvapRate;
                    
                    if (pot_evap < 0.0) {
                        if (qeSoilWaterActual > 0.0) qeSoilWaterActual = 0.0;
                        if (qeSoilWaterActual < pot_evap) qeSoilWaterActual = pot_evap;
                    }
#else
                    real qeSoilWaterActual = pot_evap;
#endif
                    
                    gw.coef(idom, 7) += gdom.dt * qeSoilWaterActual / gdom.dz(iGlob);
                    gw.evapGW(iGlob) = qeSoilWaterActual;
                    
                    if (gw.h(iGlob, 1) < -1000.0) {
                        gw.h(iGlob, 1) = -1000.0;
                    }
                }
            });
        } else {
            // Original time series interpolation approach
            real qt = interpolateLinear(tran, gdom.etime);
            real zm = interpolateLinear(zm_series, gdom.etime);
            real qe = interpolateLinear(evap, gdom.etime);
            
            k_max_root = find_max_root_depth_layer(gdom, zm);
            rootCoef(gw, gdom);
            
            Kokkos::parallel_for("gw_et_timeseries", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                idom = (kk - hc) * gdom.nx * gdom.ny + (jj - hc) * gdom.nx + ii - hc;
                int iGlobSW = jj * gdom.nxhc + ii;
                
                gw.coef(idom, 7) += gdom.dt * qt * coef_wat(idx) * coef_root(idx) / gdom.dz(iGlob);
                gw.transpGW(iGlob) = qt * coef_wat(idx) * coef_root(idx);
                
                if (kk == 1) {
#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL
                    real swEvapRate = gw.surfaceWaterEvapActual(iGlobSW);
                    real qeSoilWaterActual = qe - swEvapRate;
                    
                    if (qe < 0.0) {
                        if (qeSoilWaterActual > 0.0) qeSoilWaterActual = 0.0;
                        if (qeSoilWaterActual < qe) qeSoilWaterActual = qe;
                    }
#else
                    real qeSoilWaterActual = qe;
#endif
                    
                    gw.coef(idom, 7) += gdom.dt * qeSoilWaterActual / gdom.dz(iGlob);
                    gw.evapGW(iGlob) = qeSoilWaterActual;
                    
                    if (gw.h(iGlob, 1) < -1000.0) {
                        gw.h(iGlob, 1) = -1000.0;
                    }
                }
            });
        }
      }

      // Flux Source/Sink
      else if (sstype == 1)
      {
        real qbc = interpolateLinear(ts, gdom.etime);
        Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;

                        gw.coef(idom,7) += gdom.dt * qbc; });
      }
      // Internal Drainage with Fixed Head
      // Note: For now, this only supports draining in the saturated zone
      else if (sstype == 2)
      {
        // real hbc = interpolateLinear(ts, gdom.etime);
        real hbc = ssvals[0];
        if (direction == XPLUS || direction == XMINUS)
        {
          Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        real flux = 0.0;
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                  		// get index of cells next to the source/sink
                        int jp = idom+1, jm = idom-1, kp = idom+gdom.nx*gdom.ny, km = idom-gdom.nx*gdom.ny;
                        // treat the source/sink as a pressure boundary
                        gw.coef(jp,7) -= Cpipe * gw.coef(idom,3) * hbc;
                        gw.coef(jm,7) -= Cpipe * gw.coef(idom,4) * hbc;
                        gw.coef(kp,7) -= Cpipe * gw.coef(idom,5) * hbc;
                        gw.coef(km,7) -= Cpipe * gw.coef(idom,6) * hbc;
                        // exclude the source/sink cell from linear system
                        gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                        gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                        gw.coef(jp,4) = 0.0;	gw.coef(jm,3) = 0.0;
                        gw.coef(kp,6) = 0.0;	gw.coef(km,5) = 0.0; });
        }
        else if (direction == YPLUS || direction == YMINUS)
        {
          Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        int ip = idom+1, im = idom-1, kp = idom+gdom.nx*gdom.ny, km = idom-gdom.nx*gdom.ny;
                        // Seepage as a fixed H condition
                        gw.coef(ip,7) -= Cpipe * gw.coef(idom,1) * hbc;
                        gw.coef(im,7) -= Cpipe * gw.coef(idom,2) * hbc;
                        gw.coef(kp,7) -= Cpipe * gw.coef(idom,5) * hbc;
                        gw.coef(km,7) -= Cpipe * gw.coef(idom,6) * hbc;

                        gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                        gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                        gw.coef(ip,2) = 0.0;	gw.coef(im,1) = 0.0;
                        gw.coef(kp,6) = 0.0;	gw.coef(km,5) = 0.0; });
        }
        else
        {
          Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx) {
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        real flux = 0.0;
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        // get index of cells next to the source/sink
                        int ip = idom+1, im = idom-1, jp = idom+1, jm = idom-1;
                        // treat the source/sink as a pressure boundary
                        gw.coef(jp,7) -= Cpipe * gw.coef(idom,3) * hbc;
                        gw.coef(jm,7) -= Cpipe * gw.coef(idom,4) * hbc;
                        gw.coef(ip,7) -= Cpipe * gw.coef(idom,1) * hbc;
                        gw.coef(im,7) -= Cpipe * gw.coef(idom,2) * hbc;
                        // exclude the source/sink cell from linear system
                        gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                        gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                        gw.coef(jp,4) = 0.0;	gw.coef(jm,3) = 0.0;
                        gw.coef(ip,2) = 0.0;	gw.coef(im,1) = 0.0; });
        }
      }
    }
  }

  inline void applyWCSS(GwState &gw, GwDomain &gdom)
  {
    if (ncellsIT > 0)
    {
      Qoutflow = 0.0;
      Qoutflow_qe = 0.0;

      Qinflow = 0.0;
      // ET (Penman-Monteith)
      if (sstype == 0)
      {
        real qt = interpolateLinear(tran, gdom.etime);

        real qe = interpolateLinear(evap, gdom.etime);

        // real qt = -1.157e-8;//!zzb 0.01m/d =1.157e-7m/s

        Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp, real &tmp_qe) {
          int ii, jj, kk, ivg, idom, iGlob = icells[idx];
          gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
          idom = (kk - hc) * gdom.nx * gdom.ny + (jj - hc) * gdom.nx + ii - hc;
          int iGlobSW = jj * gdom.nxhc + ii;

          gw.wc(iGlob, 1) += gdom.dt * qt * coef_wat(idx) * coef_root(idx) / gdom.dz(iGlob);
          // tmp += qt * gdom.dt * coef_wat(idx) * coef_root(idx)  / gdom.dz(iGlob);

          tmp += qt * coef_wat(idx) * coef_root(idx);
          // if (coef_root(idx)!=0){
          // printf("idx:%d, qt: %e,coef_wat(idx):%e, coef_root(idx):%e\n", idx, qt,coef_wat(idx), coef_root(idx));
          // printf("tmp:%e\n", tmp);
          // }
          // 假设这是一个全局变量或者在合适的范围内定义的静态变量
          // static double qe_potential = -1.157e-6; // 初始假设 qe 为 -1
          // double qe = 0;
          //! zzb 三段式蒸发模型
          // todo thetaf田间持水量0.3455对应饱和含水量0.36
          // double thetaf =  0.3455;
          // if (kk == 1) {
          // 	if (gw.wc(iGlob,1) >= 0.65*thetaf) {
          // 		qe = qe_potential;
          // 		// std::cout << "111bcvals(ibc)"<<bcvals(ibc) << "\n";
          // 	}
          // 	else if (gw.wc(iGlob,1) < 0.65*thetaf && gw.wc(iGlob,1) >= 0.07) {
          // 		qe = qe_potential*((gw.wc(iGlob,1)-0.07)/(0.65*thetaf-0.07));
          // 		// std::cout << "222bcvals(ibc)"<<bcvals(ibc) << "\n";
          // 	}
          // 	else {
          // 		qe = 0;
          // 		// std::cout << "333bcvals(ibc)"<<bcvals(ibc) << "\n";
          // 	}

          // }
          // printf("11111111kk:%d, qt: %e , qe: %e \n",kk,  qt, qe);
          if (kk == 1)
          {

#if SW_GW_EVAPORATION_TRANSPIRATION_MODEL

            // 1. 获取基础变量
            // qe 通常为负值 (例如 -10 mm/d)，代表通量流出
            real qeTotalRate = qe; 
            
            // 地表水实际蒸发强度 (已在地表水模块计算完成，应为负值或0)
            // 单位与 qe 一致 (m/s)
            real swEvapRate = gw.surfaceWaterEvapActual(iGlobSW); 

            // 2. 计算土壤蒸发强度 (剩余能量法)
            // 逻辑：土壤只能蒸发“地表水没用完”的那部分能量
            // 公式：Total = Surface + Soil  =>  Soil = Total - Surface
            // 示例：总能力(-10) - 地表实际(-2) = 土壤(-8)
            real qeSoilWaterActual = qeTotalRate - swEvapRate;

            // 3. 物理约束修正 (防守性编程)
            // 修正1：如果计算出的土壤蒸发比总蒸发还大（符号为负，绝对值更大），说明数值异常
            if (qeTotalRate < 0.0) // 正常蒸发情况
            {
                // 如果算出来是正数（变成了降水），或比总能力还负（过度蒸发），需要截断
                if (qeSoilWaterActual > 0.0) qeSoilWaterActual = 0.0;
                if (qeSoilWaterActual < qeTotalRate) qeSoilWaterActual = qeTotalRate;
            }
            
            // 4. (可选) 蒸发类型标记，仅用于调试或统计
            int SoilWaterEvapType = 999;
            real epsilon = 1e-10; // 浮点数容差
            
            if (fabs(swEvapRate - qeTotalRate) < epsilon) {
                SoilWaterEvapType = 0; // 类型0：完全由地表水满足 (Soil = 0)
            } else if (fabs(swEvapRate) < epsilon) {
                SoilWaterEvapType = 2; // 类型2：完全由土壤水满足 (Surface = 0)
            } else {
                SoilWaterEvapType = 1; // 类型1：混合蒸发 (地表水蒸干了，剩下蒸土)
            }

            // printf("Type: %d, Total: %e, Surf: %e, Soil: %e\n", SoilWaterEvapType, qeTotalRate, swEvapRate, qeSoilWaterActual);

#else
            // 非耦合模式，全部潜在蒸发都作用于土壤
            real qeSoilWaterActual = qe;
#endif
            gw.wc(iGlob, 1) += gdom.dt * qeSoilWaterActual / gdom.dz(iGlob);


            tmp_qe += qeSoilWaterActual;

            // printf("tmp_qe:%e, qeSoilWaterActual: %e\n", tmp_qe, qeSoilWaterActual);
            //! zzb gwbc形式模拟蒸发比例，限制上边界水头为-1000

            if (gw.h(iGlob, 1) < -1000)
            {
              gw.h(iGlob, 1) = -1000;
            }
          } }, Kokkos::Sum<real>(Qoutflow), Kokkos::Sum<real>(Qoutflow_qe));
      }
      // Flux Source/Sink
      else if (sstype == 1)
      {
        real qbc = interpolateLinear(ts, gdom.etime);
        Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp) {
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        gw.wc(iGlob,1) += gdom.dt * qbc;
                        if (qbc > 0)	{tmp += qbc * gdom.dt;}
                        else {tmp -= qbc * gdom.dt;} }, Kokkos::Sum<real>(Qinflow));
      }
      else if (sstype == 2)
      {
        // real hbc = interpolateLinear(ts, gdom.etime);
        real hbc = ssvals[0], flux;
        if (direction == XPLUS || direction == XMINUS)
        {
          Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA(int idx, real &tmp) {
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        // Calculate the cumulative outflow
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
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        // Only consider drainage under fully saturated condition
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
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        // Only consider drainage under fully saturated condition
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
    }
  }
};
#endif

class SourceSink
{
public:
  std::vector<std::string> id;
  SourceSinkData swss;
#if SERGHEI_SUBSURFACE_MODEL
  std::vector<GwSS> gwss;
#endif
};

#endif
