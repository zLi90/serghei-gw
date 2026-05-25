/* -*- mode: c++; c-default-style: "linux" -*- */
#pragma once

#if SERGHEI_WAVE_MODEL

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>

#include "define.h"
#include "GwBC.h"
#include "GwDomain.h"
#include "FileIO.h"
#include "Parallel.h"

class WaveBoundaryState
{
public:
  int enabled = 0;
  int nActive = 0;
  intArr gwbcIndex;
  intArr bcLocalIndex;
  intArr iGlobCell;
  realArr fetchLength;
  realArr xCoord;
  realArr yCoord;
  realArr xProj;
  realArr meanLakeLevel;
  realArr Hs;
  realArr Ts;
  realArr kWave;
  realArr eta;
  realArr hBound;

  void allocate(int n)
  {
    nActive = n;
    gwbcIndex = intArr("wave_gwbc_index", nActive);
    bcLocalIndex = intArr("wave_bc_local_index", nActive);
    iGlobCell = intArr("wave_iglob_cell", nActive);
    fetchLength = realArr("wave_fetch", nActive);
    xCoord = realArr("wave_xcoord", nActive);
    yCoord = realArr("wave_ycoord", nActive);
    xProj = realArr("wave_xproj", nActive);
    meanLakeLevel = realArr("wave_mean_lake_level", nActive);
    Hs = realArr("wave_hs", nActive);
    Ts = realArr("wave_ts", nActive);
    kWave = realArr("wave_kwave", nActive);
    eta = realArr("wave_eta", nActive);
    hBound = realArr("wave_hbound", nActive);
  }
};

class WaveBoundaryInit
{
private:
  bool readFetchRaster(
      const std::string &fname,
      int expectedNx,
      int expectedNy,
      std::vector<real> &fetchGrid,
      Parallel const &par) const
  {
    std::ifstream in(fname);
    if (!in.is_open())
      return false;

    std::string key;
    int nx = -1;
    int ny = -1;
    real tmp = 0.0;
    for (int i = 0; i < 6; ++i)
    {
      if (!(in >> key >> tmp))
        return false;
      if (i == 0)
        nx = static_cast<int>(tmp);
      if (i == 1)
        ny = static_cast<int>(tmp);
    }

    if (nx != expectedNx || ny != expectedNy)
    {
      if (par.masterproc)
      {
        std::cerr << RERROR << "fetch raster dimensions do not match domain size (" << nx << "," << ny
                  << ") vs (" << expectedNx << "," << expectedNy << ")" << std::endl;
      }
      return false;
    }

    fetchGrid.resize(nx * ny, SERGHEI_NAN);
    for (int j = 0; j < ny; ++j)
    {
      for (int i = 0; i < nx; ++i)
      {
        if (!(in >> fetchGrid[j * nx + i]))
        {
          return false;
        }
      }
    }
    return true;
  }

  bool readFetchList(const std::string &fname, int expectedCount, std::vector<real> &values) const
  {
    std::ifstream in(fname);
    if (!in.is_open())
      return false;

    values.clear();
    real v = 0.0;
    while (in >> v)
      values.push_back(v);

    return static_cast<int>(values.size()) == expectedCount;
  }

public:
  bool initialize(
      WaveBoundaryState &waveState,
      GwDomain const &gdom,
      std::vector<GwBC> const &gwbc,
      Parallel const &par,
      FileIO const &io)
  {
    waveState.enabled = 0;
    waveState.nActive = 0;

    if (!io.wave.enabled)
      return true;

    std::vector<int> hostGwbc;
    std::vector<int> hostLocal;
    std::vector<int> hostIGlob;
    std::vector<real> hostX;
    std::vector<real> hostY;

    for (int k = 0; k < static_cast<int>(gwbc.size()); ++k)
    {
      if (gwbc[k].direction != ZMINUS)
        continue;

      bool isDirichletHead =
          (gwbc[k].bctype == SUB_BC_H_CONST) ||
          (gwbc[k].bctype == SUB_BC_H_T) ||
          (gwbc[k].bctype == SUB_BC_WT_CONST) ||
          (gwbc[k].bctype == SUB_BC_WT_T);

      if (!isDirichletHead)
        continue;

      for (int ibc = 0; ibc < gwbc[k].ncellsBC; ++ibc)
      {
        int iGlob = gwbc[k].bcells(ibc);
        hostGwbc.push_back(k);
        hostLocal.push_back(ibc);
        hostIGlob.push_back(iGlob);
        hostX.push_back(gdom.x(iGlob));
        hostY.push_back(gdom.y(iGlob));
      }
    }

    if (hostIGlob.empty())
    {
      if (par.masterproc)
        std::cerr << YEXC << "Wave module enabled, but no active top Dirichlet RE boundary cells were found." << std::endl;
      return false;
    }

    waveState.allocate(static_cast<int>(hostIGlob.size()));
    waveState.enabled = 1;

    for (int i = 0; i < waveState.nActive; ++i)
    {
      waveState.gwbcIndex(i) = hostGwbc[i];
      waveState.bcLocalIndex(i) = hostLocal[i];
      waveState.iGlobCell(i) = hostIGlob[i];
      waveState.xCoord(i) = hostX[i];
      waveState.yCoord(i) = hostY[i];
      waveState.xProj(i) = 0.0;
      waveState.meanLakeLevel(i) = io.wave.meanLakeLevel;
      waveState.fetchLength(i) = SERGHEI_NAN;
      waveState.Hs(i) = 0.0;
      waveState.Ts(i) = 1.0;
      waveState.kWave(i) = 0.0;
      waveState.eta(i) = 0.0;
      waveState.hBound(i) = io.wave.meanLakeLevel;
    }

    std::string fetchPath = io.inFolder + io.wave.fetchFile;
    if (io.wave.fetchFile.size() > 0 && io.wave.fetchFile[0] == '/')
      fetchPath = io.wave.fetchFile;

    std::vector<real> fetchValues;
    bool ok = false;

    if (!io.wave.fetchSource.compare("precomputed"))
    {
      std::vector<real> fetchGrid;
      ok = readFetchRaster(fetchPath, gdom.nx_glob, gdom.ny_glob, fetchGrid, par);
      if (ok)
      {
        for (int i = 0; i < waveState.nActive; ++i)
        {
          int ii, jj, kk;
          gdom.unpackIndicesHalo(waveState.iGlobCell(i), kk, jj, ii);
          int ig = par.i_beg + (ii - gdom.hc);
          int jg = par.j_beg + (jj - gdom.hc);
          int idx = jg * gdom.nx_glob + ig;
          waveState.fetchLength(i) = fetchGrid[idx];
        }
      }
    }
    else
    {
      ok = readFetchList(fetchPath, waveState.nActive, fetchValues);
      if (ok)
      {
        for (int i = 0; i < waveState.nActive; ++i)
          waveState.fetchLength(i) = fetchValues[i];
      }
    }

    if (!ok)
    {
      if (par.masterproc)
      {
        std::cerr << RERROR << "Could not read fetch values from " << fetchPath << std::endl;
      }
      return false;
    }

    for (int i = 0; i < waveState.nActive; ++i)
    {
      if (!(waveState.fetchLength(i) > 0.0))
      {
        if (par.masterproc)
        {
          std::cerr << RERROR << "Invalid fetch length at mapped wave boundary cell index " << i << std::endl;
        }
        return false;
      }
    }

    if (par.masterproc)
    {
      std::cout << GOK << "Wave boundary preprocessing completed" << std::endl;
      std::cout << BDASH << "active RE top boundary cells: " << waveState.nActive << std::endl;
      std::cout << BDASH << "fetch source file: " << fetchPath << std::endl;
    }

    return true;
  }
};

#endif

