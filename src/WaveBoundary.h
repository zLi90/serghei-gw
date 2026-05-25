/* -*- mode: c++; c-default-style: "linux" -*- */
#pragma once

#if SERGHEI_WAVE_MODEL

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "define.h"
#include "FileIO.h"
#include "GwBC.h"
#include "GwDomain.h"
#include "SourceSink.h"
#include "WaveBoundaryInit.h"

class WaveBoundaryModel
{
private:
  KOKKOS_INLINE_FUNCTION
  real clipPositive(real v, real vmin) const
  {
    return (v > vmin) ? v : vmin;
  }

  int getWaveMethodId(std::string const &waveMethod) const
  {
    if (!waveMethod.compare("jonswap_simple"))
      return 2;
    return 1; // default SMB
  }

  int getPhaseModeId(std::string const &phaseMode) const
  {
    if (!phaseMode.compare("xcoordinate"))
      return 2;
    return 1; // default fetchprojection
  }

  real interpolateSeriesAt(TimeSeries const &ts, real t, int icol) const
  {
    if (ts.np <= 0)
      return 0.0;
    if (ts.np == 1)
      return ts.value(icol * ts.np);

    if (t <= ts.time(0))
      return ts.value(icol * ts.np);
    if (t >= ts.time(ts.np - 1))
      return ts.value(icol * ts.np + ts.np - 1);

    int i0 = 0;
    for (int i = 0; i < ts.np - 1; ++i)
    {
      if (t >= ts.time(i) && t < ts.time(i + 1))
      {
        i0 = i;
        break;
      }
    }
    int i1 = i0 + 1;
    real t0 = ts.time(i0);
    real t1 = ts.time(i1);
    real v0 = ts.value(icol * ts.np + i0);
    real v1 = ts.value(icol * ts.np + i1);
    real w = (t - t0) / (t1 - t0);
    return v0 + (v1 - v0) * w;
  }

public:
  void update_wave_boundary(
      WaveBoundaryState &waveState,
      SourceSinkData &swss,
      Domain const &dom,
      GwDomain const &gdom,
      FileIO const &io)
  {
    if (!waveState.enabled || waveState.nActive <= 0)
      return;

    real U10 = 0.0;
    real windDirDeg = 0.0;

    if (!io.wave.windSource.compare("constant"))
    {
      U10 = dom.waveWindSpeed;
      windDirDeg = dom.waveWindDirection;
    }
    else
    {
      // Use robust interpolation independent of mutable timeIndex.
      U10 = interpolateSeriesAt(swss.wind, gdom.etime, 0);
      windDirDeg = interpolateSeriesAt(swss.wind, gdom.etime, 1);
    }

    U10 = (U10 > 0.0) ? U10 : 0.0;
    const real g = GRAV;
    const real omegaMin = 2.0 * PI / 120.0; // lower frequency bound (T <= 120 s)
    const real TsMin = 0.5;                 // avoid division instability

    const real windDirRad = windDirDeg * PI / 180.0;
    const real windCx = cos(windDirRad);
    const real windCy = sin(windDirRad);
    const int methodId = getWaveMethodId(io.wave.waveMethod);
    const int phaseModeId = getPhaseModeId(io.wave.phaseMode);
    const real t = gdom.etime;

    // Kernel A: compute Hs, Ts, k and phase coordinate.
    Kokkos::parallel_for("compute_wave_params", waveState.nActive, KOKKOS_CLASS_LAMBDA(int i) {
      const real F = waveState.fetchLength(i);
      const real x = waveState.xCoord(i);
      const real y = waveState.yCoord(i);

      real HsLocal = 0.0;
      real TsLocal = TsMin;

      if (U10 > 0.0 && F > 0.0)
      {
        const real X = g * F / (U10 * U10);

        if (methodId == 2)
        {
          // Simplified JONSWAP fetch-limited proxy
          const real hsNd = 0.0016 * mypow(X, 0.5);
          const real tpNd = 0.2857 * mypow(X, 0.33);
          HsLocal = hsNd * (U10 * U10 / g);
          TsLocal = tpNd * (U10 / g);
        }
        else
        {
          // SMB fetch-limited growth curves
          const real hsNd = 0.283 * tanh(0.0125 * mypow(X, 0.42));
          const real tpNd = 7.54 * tanh(0.077 * mypow(X, 0.25));
          HsLocal = hsNd * (U10 * U10 / g);
          TsLocal = tpNd * (U10 / g);
        }
      }

      HsLocal = clipPositive(HsLocal, 0.0);
      TsLocal = clipPositive(TsLocal, TsMin);

      const real omega = max(2.0 * PI / TsLocal, omegaMin);
      const real kLocal = (omega * omega) / g; // deep-water approximation

      real xProjLocal = x;
      if (phaseModeId == 1)
      {
        xProjLocal = x * windCx + y * windCy;
      }

      waveState.Hs(i) = HsLocal;
      waveState.Ts(i) = TsLocal;
      waveState.kWave(i) = kLocal;
      waveState.xProj(i) = xProjLocal;
    });

    // Kernel B: transient elevation and dynamic boundary head.
    Kokkos::parallel_for("compute_wave_head", waveState.nActive, KOKKOS_CLASS_LAMBDA(int i) {
      const real HsLocal = waveState.Hs(i);
      const real TsLocal = waveState.Ts(i);
      const real kLocal = waveState.kWave(i);
      const real xProjLocal = waveState.xProj(i);
      const real Hmean = waveState.meanLakeLevel(i);

      const real phase = 2.0 * PI * t / TsLocal - kLocal * xProjLocal;
      const real etaLocal = 0.5 * HsLocal * cos(phase);
      const real hBoundLocal = Hmean + etaLocal;

      waveState.eta(i) = etaLocal;
      waveState.hBound(i) = hBoundLocal;
    });
  }

  void apply_wave_head_to_bc(
      WaveBoundaryState const &waveState,
      std::vector<GwBC> &gbc) const
  {
    if (!waveState.enabled || waveState.nActive <= 0)
      return;

    Kokkos::fence();
    for (int i = 0; i < waveState.nActive; ++i)
    {
      int k = waveState.gwbcIndex(i);
      int ibc = waveState.bcLocalIndex(i);
      gbc[k].bcvals(ibc) = waveState.hBound(i);
    }
  }
};

#endif

