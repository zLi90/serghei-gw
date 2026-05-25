#!/usr/bin/env python3
"""
Visualize wave-generated top-boundary head perturbations from SERGHEI inputs.

This script reads:
  - wave.input
  - fetch.input (or boundary-list fetch file)
  - wind.input (when wave windSource = windfile)

It mirrors the wave formulas implemented in src/WaveBoundary.h:
  - SMB fetch-limited growth (default)
  - jonswap_simple fetch-limited proxy
  - phase modes: fetchprojection or xcoordinate
  - eta = 0.5 * Hs * cos(2*pi*t/Ts - k*x_proj)
  - hBound = meanLakeLevel + eta

Outputs:
  - fetch map
  - Hs/Ts maps at reference time
  - eta/hBound snapshot maps at selected times
  - spatiotemporal eta heatmap
  - point eta/hBound time-series
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# User options (hard-coded; edit here)
# ---------------------------------------------------------------------------
# Folder containing wave.input/fetch.input and optional wind.input
CASE_DIR = Path("../t1-lake2d/input-closedwave")

# Output folder for figures. Set to None to use CASE_DIR / "fig-wave".
OUTDIR = Path("../t1-lake2d/fig-wave-test")

# Time axis settings used mainly for constant-wind cases.
T_END = 86400.0  # [s]
DT = 300.0       # [s], used only when NT <= 1 for windfile mode
NT = 240         # number of sampled time points

# Number of eta/hBound snapshot maps to save.
N_SNAPSHOTS = 4

# Fixed colorbar ranges for time-snapshot maps.
# Set to None to use matplotlib automatic scaling.
# These ranges are applied consistently to ALL map_eta_t*.png figures.
ETA_VMIN = -0.03
ETA_VMAX = 0.03

# These ranges are applied consistently to ALL map_hbound_t*.png figures.
HBOUND_VMIN = 0.47
HBOUND_VMAX = 0.53


G = 9.81
PI = math.pi
TS_MIN = 0.5
OMEGA_MIN = 2.0 * PI / 120.0


@dataclass
class Raster:
    nx: int
    ny: int
    xll: float
    yll: float
    dx: float
    nodata: float
    values: np.ndarray  # shape (ny, nx)


def _clean_line(line: str) -> str:
    # Match parser behavior: allow comments after //
    line = line.split("//", 1)[0]
    line = line.split("#", 1)[0]
    return line.strip()


def parse_kv_file(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = _clean_line(raw)
        if not line:
            continue
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        out[k.strip().lower()] = v.strip()
    return out


def parse_raster(path: Path) -> Raster:
    with path.open() as f:
        header = {}
        for _ in range(6):
            line = _clean_line(f.readline())
            if not line:
                raise ValueError(f"Malformed raster header in {path}")
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Malformed raster header line '{line}' in {path}")
            key = parts[0].lower()
            val = float(parts[1])
            header[key] = val
        nx = int(header["ncols"])
        ny = int(header["nrows"])
        xll = float(header["xllcorner"])
        yll = float(header["yllcorner"])
        dx = float(header["cellsize"])
        nodata = float(header["nodata_value"])

        vals: List[float] = []
        for raw in f:
            line = _clean_line(raw)
            if not line:
                continue
            vals.extend(float(x) for x in line.split())
    if len(vals) != nx * ny:
        raise ValueError(f"Raster data size mismatch in {path}: expected {nx*ny}, got {len(vals)}")
    arr = np.array(vals, dtype=float).reshape(ny, nx)
    return Raster(nx=nx, ny=ny, xll=xll, yll=yll, dx=dx, nodata=nodata, values=arr)


def parse_wind_file(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Parser.h readWindFile format:
    # np <int>
    # CwT <real>
    # hwmin <real>
    # t0 U0 dir0
    # ...
    lines = [_clean_line(ln) for ln in path.read_text().splitlines()]
    lines = [ln for ln in lines if ln]
    if len(lines) < 4:
        raise ValueError(f"wind.input appears too short: {path}")

    nparts = lines[0].split()
    if len(nparts) < 2 or nparts[0].lower() != "np":
        raise ValueError(f"First line must be 'np <N>' in {path}")
    n = int(float(nparts[1]))

    t = []
    speed = []
    direction = []
    for ln in lines[3:]:
        p = ln.split()
        if len(p) < 3:
            continue
        t.append(float(p[0]))
        speed.append(float(p[1]))
        direction.append(float(p[2]))
    if len(t) != n:
        raise ValueError(f"wind.input record mismatch in {path}: np={n}, rows={len(t)}")
    return np.array(t), np.array(speed), np.array(direction)


def interp_series(times: np.ndarray, values: np.ndarray, t: np.ndarray) -> np.ndarray:
    # Clamp outside range, linear interpolate inside (matches C++ logic).
    return np.interp(t, times, values, left=values[0], right=values[-1])


def compute_wave(
    fetch: np.ndarray,
    xproj: np.ndarray,
    mean_lake_level: float,
    t: np.ndarray,
    wind_speed_t: np.ndarray,
    wave_method: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns:
      Hs_tn (nt, ncell)
      Ts_tn (nt, ncell)
      eta_tn (nt, ncell)
      hbound_tn (nt, ncell)
    """
    nt = t.size
    ncell = fetch.size
    fetch2 = fetch[None, :]  # (1, ncell)
    U = np.maximum(wind_speed_t[:, None], 0.0)  # (nt, 1)

    Hs = np.zeros((nt, ncell), dtype=float)
    Ts = np.full((nt, ncell), TS_MIN, dtype=float)

    active = (U > 0.0) & (fetch2 > 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        X = G * fetch2 / (U * U)

    if wave_method == "jonswap_simple":
        hs_nd = 0.0016 * np.power(np.maximum(X, 0.0), 0.5)
        tp_nd = 0.2857 * np.power(np.maximum(X, 0.0), 0.33)
    else:
        hs_nd = 0.283 * np.tanh(0.0125 * np.power(np.maximum(X, 0.0), 0.42))
        tp_nd = 7.54 * np.tanh(0.077 * np.power(np.maximum(X, 0.0), 0.25))

    with np.errstate(divide="ignore", invalid="ignore"):
        Hs_calc = hs_nd * (U * U / G)
        Ts_calc = tp_nd * (U / G)
    Hs[active] = Hs_calc[active]
    Ts[active] = Ts_calc[active]
    Hs = np.maximum(Hs, 0.0)
    Ts = np.maximum(Ts, TS_MIN)

    omega = np.maximum(2.0 * PI / Ts, OMEGA_MIN)
    k_wave = omega * omega / G
    phase = 2.0 * PI * t[:, None] / Ts - k_wave * xproj[None, :]

    eta = 0.5 * Hs * np.cos(phase)
    hbound = mean_lake_level + eta
    return Hs, Ts, eta, hbound


def build_time_axis(
    wave_cfg: Dict[str, str],
    case_dir: Path,
    t_end: float,
    dt: float,
    nt: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    wind_source = wave_cfg.get("windsource", "windfile").lower()
    if wind_source == "constant":
        ws = float(wave_cfg["windspeed"])
        wd = float(wave_cfg["winddirection"])
        times = np.linspace(0.0, t_end, nt) if nt > 1 else np.array([0.0])
        speed = np.full_like(times, ws, dtype=float)
        direction = np.full_like(times, wd, dtype=float)
        return times, speed, direction

    wind_path = case_dir / "wind.input"
    if not wind_path.exists():
        raise FileNotFoundError(f"windSource=windfile but missing {wind_path}")
    t_w, s_w, d_w = parse_wind_file(wind_path)

    if nt > 1:
        times = np.linspace(float(t_w[0]), float(t_w[-1]), nt)
    else:
        if dt > 0.0:
            tmax = max(float(t_w[-1]), t_end)
            n = int(math.floor((tmax - t_w[0]) / dt)) + 1
            times = t_w[0] + np.arange(n, dtype=float) * dt
        else:
            times = t_w.copy()
    speed = interp_series(t_w, s_w, times)
    direction = interp_series(t_w, d_w, times)
    return times, speed, direction


def ensure_outdir(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)


def save_map(
    arr2d: np.ndarray,
    title: str,
    cbar_label: str,
    outpath: Path,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    plt.figure(figsize=(9, 3))
    im = plt.imshow(arr2d, cmap="viridis", aspect="auto", vmin=vmin, vmax=vmax)
    plt.title(title)
    plt.xlabel("x index")
    plt.ylabel("y index")
    cbar = plt.colorbar(im)
    cbar.set_label(cbar_label)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    if CASE_DIR.is_absolute():
        case_dir = CASE_DIR
    else:
        case_dir = (script_dir / CASE_DIR).resolve()

    wave_path = case_dir / "wave.input"
    if not wave_path.exists():
        raise FileNotFoundError(f"Missing {wave_path}")
    wave_cfg = parse_kv_file(wave_path)

    if int(float(wave_cfg.get("waveenable", "0"))) == 0:
        raise RuntimeError("waveEnable is 0 in wave.input. No wave to visualize.")

    mean_lake = float(wave_cfg["meanlakelevel"])
    wave_method = wave_cfg.get("wavemethod", "smb").lower()
    if wave_method not in {"smb", "jonswap_simple"}:
        raise ValueError(f"Unsupported wave method '{wave_method}'")

    fetch_source = wave_cfg.get("fetchsource", "precomputed").lower()
    fetch_file = wave_cfg.get("fetchfile", "fetch.input")
    fetch_path = Path(fetch_file)
    if not fetch_path.is_absolute():
        fetch_path = case_dir / fetch_file
    if not fetch_path.exists():
        raise FileNotFoundError(f"Missing fetch file: {fetch_path}")

    if OUTDIR is None:
        outdir = case_dir / "fig-wave"
    else:
        if OUTDIR.is_absolute():
            outdir = OUTDIR
        else:
            outdir = (script_dir / OUTDIR).resolve()
    ensure_outdir(outdir)

    if fetch_source == "precomputed":
        rast = parse_raster(fetch_path)
        fetch2d = rast.values.copy()
        valid = (fetch2d != rast.nodata) & np.isfinite(fetch2d) & (fetch2d > 0.0)
        if not np.any(valid):
            raise RuntimeError("No valid positive fetch cells found in fetch raster.")

        jj, ii = np.indices((rast.ny, rast.nx))
        x = rast.xll + (ii + 0.5) * rast.dx
        # Match common ASC row order (north -> south), consistent with model convention usage.
        y = rast.yll + (rast.ny - (jj + 0.5)) * rast.dx
        x_flat = x[valid]
        y_flat = y[valid]
        fetch_flat = fetch2d[valid]
        shape = fetch2d.shape
    elif fetch_source == "boundarylist":
        vals = np.loadtxt(fetch_path, dtype=float)
        vals = np.atleast_1d(vals).astype(float)
        valid = np.isfinite(vals) & (vals > 0.0)
        if not np.any(valid):
            raise RuntimeError("No valid positive fetch entries found in boundary list.")
        fetch_flat = vals[valid]
        x_flat = np.arange(fetch_flat.size, dtype=float)
        y_flat = np.zeros(fetch_flat.size, dtype=float)
        shape = (1, fetch_flat.size)
        fetch2d = np.full(shape, np.nan)
        fetch2d[0, :] = fetch_flat
    else:
        raise ValueError(f"Unsupported fetchSource '{fetch_source}'")

    times, wind_speed, wind_dir_deg = build_time_axis(
        wave_cfg=wave_cfg,
        case_dir=case_dir,
        t_end=T_END,
        dt=DT,
        nt=max(2, NT),
    )

    phase_mode = wave_cfg.get("phasemode", "fetchprojection").lower()
    if phase_mode == "fetchprojection":
        theta = np.deg2rad(wind_dir_deg)
        xproj_tn = x_flat[None, :] * np.cos(theta)[:, None] + y_flat[None, :] * np.sin(theta)[:, None]
    elif phase_mode == "xcoordinate":
        xproj_tn = np.broadcast_to(x_flat[None, :], (times.size, x_flat.size))
    else:
        raise ValueError(f"Unsupported phaseMode '{phase_mode}'")

    # Compute wave fields; when xproj depends on time (fetchprojection + variable wind),
    # evaluate each time step with the corresponding projected coordinate.
    Hs_tn = np.zeros((times.size, x_flat.size), dtype=float)
    Ts_tn = np.zeros((times.size, x_flat.size), dtype=float)
    eta_tn = np.zeros((times.size, x_flat.size), dtype=float)
    hb_tn = np.zeros((times.size, x_flat.size), dtype=float)

    for it in range(times.size):
        Hs, Ts, eta, hb = compute_wave(
            fetch=fetch_flat,
            xproj=xproj_tn[it, :],
            mean_lake_level=mean_lake,
            t=np.array([times[it]], dtype=float),
            wind_speed_t=np.array([wind_speed[it]], dtype=float),
            wave_method=wave_method,
        )
        Hs_tn[it, :] = Hs[0, :]
        Ts_tn[it, :] = Ts[0, :]
        eta_tn[it, :] = eta[0, :]
        hb_tn[it, :] = hb[0, :]

    # Helper to scatter flattened valid cells back to 2D map.
    def back_to_map(vflat: np.ndarray) -> np.ndarray:
        out = np.full(shape, np.nan, dtype=float)
        if fetch_source == "precomputed":
            out[valid] = vflat
        else:
            out[0, : vflat.size] = vflat
        return out

    # 1) Base maps
    save_map(fetch2d, "Fetch length", "Fetch [m]", outdir / "map_fetch.png")
    save_map(back_to_map(Hs_tn[0]), f"Hs at t={times[0]:.1f}s", "Hs [m]", outdir / "map_hs_t0.png")
    save_map(back_to_map(Ts_tn[0]), f"Ts at t={times[0]:.1f}s", "Ts [s]", outdir / "map_ts_t0.png")

    # 2) Snapshot maps
    snap_idx = np.linspace(0, times.size - 1, max(1, N_SNAPSHOTS), dtype=int)
    for i, it in enumerate(snap_idx, start=1):
        save_map(
            back_to_map(eta_tn[it]),
            f"eta at t={times[it]:.1f}s",
            "eta [m]",
            outdir / f"map_eta_t{i:02d}.png",
            vmin=ETA_VMIN,
            vmax=ETA_VMAX,
        )
        save_map(
            back_to_map(hb_tn[it]),
            f"hBound at t={times[it]:.1f}s",
            "hBound [m]",
            outdir / f"map_hbound_t{i:02d}.png",
            vmin=HBOUND_VMIN,
            vmax=HBOUND_VMAX,
        )

    # 3) Spatiotemporal heatmap (flattened active-cell axis).
    plt.figure(figsize=(10, 4))
    im = plt.imshow(
        eta_tn.T,
        origin="lower",
        aspect="auto",
        extent=[times[0], times[-1], 0, eta_tn.shape[1] - 1],
        cmap="coolwarm",
    )
    plt.xlabel("Time [s]")
    plt.ylabel("Active boundary cell index")
    plt.title("Spatiotemporal wave elevation eta")
    cbar = plt.colorbar(im)
    cbar.set_label("eta [m]")
    plt.tight_layout()
    plt.savefig(outdir / "heatmap_eta_time_cell.png", dpi=180)
    plt.close()

    # 4) Point time series (min/max fetch and center cell by index).
    i_max = int(np.argmax(fetch_flat))
    i_min = int(np.argmin(fetch_flat))
    i_mid = int(fetch_flat.size // 2)
    idxs = [("max_fetch", i_max), ("min_fetch", i_min), ("mid_index", i_mid)]

    plt.figure(figsize=(10, 4))
    for name, idx in idxs:
        plt.plot(times, eta_tn[:, idx], label=f"eta ({name}, idx={idx})")
    plt.axhline(0.0, color="k", lw=0.8, alpha=0.5)
    plt.xlabel("Time [s]")
    plt.ylabel("eta [m]")
    plt.title("Wave elevation time series at selected boundary cells")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "timeseries_eta_selected_cells.png", dpi=180)
    plt.close()

    plt.figure(figsize=(10, 4))
    for name, idx in idxs:
        plt.plot(times, hb_tn[:, idx], label=f"hBound ({name}, idx={idx})")
    plt.axhline(mean_lake, color="k", lw=0.8, alpha=0.5, label="meanLakeLevel")
    plt.xlabel("Time [s]")
    plt.ylabel("hBound [m]")
    plt.title("Top boundary head (meanLakeLevel + eta)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "timeseries_hbound_selected_cells.png", dpi=180)
    plt.close()

    print("Wave visualization complete.")
    print(f"  Case dir: {case_dir}")
    print(f"  Figure dir: {outdir}")
    print(f"  Time samples: {times.size}, active cells: {fetch_flat.size}")
    print(f"  Wave method: {wave_method}, phase mode: {phase_mode}")


if __name__ == "__main__":
    main()
