#!/usr/bin/env python3
"""
Post-process SERGHEI lake-groundwater exchange outputs.

This script focuses on:
1) Domain-integrated exchange time series (water + solute)
2) Spatial distribution maps from NetCDF outputs (if available)
3) Printing all NetCDF keys (variables + dimensions) to screen

Expected files in OUTDIR:
  - SubsurfaceTimeSeries.out
  - RTSubsurfaceTimeSeries.out (optional, for solute)
  - output_subsurface.nc (optional, for spatial water exchange map)
  - output_transport.nc (optional, for spatial solute map/proxy exchange)
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import xarray as xr
except Exception:  # pragma: no cover
    xr = None

try:
    from netCDF4 import Dataset
except Exception:  # pragma: no cover
    Dataset = None

# ---------------------------------------------------------------------
# User options (edit here)
# ---------------------------------------------------------------------
# Output folder produced by SERGHEI (contains *.out and optionally *.nc).
OUTDIR = Path("../t1-lake2d/out-base")

# Where plots/CSVs will be saved.
# - None: automatically use OUTDIR / "figures"
# - Path(...): use a custom location
FIGDIR: Optional[Path] = None  # None -> OUTDIR / "figures"

# Solute species index to post-process from RTSubsurfaceTimeSeries.out
# and output_transport.nc.
SPECIES_INDEX = 0

# Time index used for spatial snapshots from NetCDF:
# -1 means the last saved time step.
TIME_INDEX = -1

# Vertical layer index used for spatial maps:
# for this setup, 0 corresponds to the top groundwater layer.
LAYER_INDEX = 1

# Optional colorbar ranges for spatial figures (None = matplotlib auto-range).
# qz map range [m/s]
QZ_VMIN: Optional[float] = -1e-7
QZ_VMAX: Optional[float] = 1e-7
# concentration map range [mg/L]
C_VMIN: Optional[float] = 0.0
C_VMAX: Optional[float] = 1.0
# proxy solute exchange map range [mg/(m2*s)]
SOLUTE_EXCHANGE_VMIN: Optional[float] = 0.0
SOLUTE_EXCHANGE_VMAX: Optional[float] = 5e-4


def _read_numeric_table(path: Path) -> np.ndarray:
    """
    Read whitespace-delimited numeric table while ignoring first header line.

    The SERGHEI *.out files are written as:
      line 1: header with names/units
      line 2+: numeric values
    """
    rows = []
    with path.open("r", encoding="utf-8") as f:
        _ = f.readline()  # header
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append([float(x) for x in line.split()])
    if not rows:
        raise ValueError(f"No numeric rows found in {path}")
    return np.asarray(rows, dtype=float)


def _read_header_tokens(path: Path) -> list[str]:
    """Read first-line header tokens from a SERGHEI table file."""
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip()
    if not header:
        raise ValueError(f"Empty header in {path}")
    raw = header.split()
    tokens: list[str] = []
    i = 0
    while i < len(raw):
        name = raw[i]
        # Merge split unit token, e.g. "SubSurfaceVolume [m3]" -> "SubSurfaceVolume[m3]".
        if i + 1 < len(raw) and raw[i + 1].startswith("["):
            name = f"{name}{raw[i + 1]}"
            i += 1
        tokens.append(name)
        i += 1
    return tokens


def _build_water_df(path: Path) -> pd.DataFrame:
    """
    Parse SubsurfaceTimeSeries.out into a labeled dataframe.

    Important interpretation used here:
      - bc_inflow_m3s: mostly lateral inflow (model-dependent).
      - bc_outflow_m3s: interpreted by the user as top-boundary exchange term.
      - exchange_rate_reported_m3s: file column "ExchangeRate", which can be
        zero in RE-only runs depending on compile/runtime coupling path.
    """
    header_tokens = _read_header_tokens(path)
    arr = _read_numeric_table(path)
    if arr.shape[1] < 2:
        raise ValueError(f"Unexpected column count in {path}: {arr.shape[1]}")

    # New format: BoundaryFluxNet_i and BoundaryFluxAbs_i for each boundary.
    net_cols = [i for i, name in enumerate(header_tokens) if name.startswith("BoundaryFluxNet_")]
    abs_cols = [i for i, name in enumerate(header_tokens) if name.startswith("BoundaryFluxAbs_")]

    if net_cols and abs_cols:
        # Aggregate over all listed boundaries so one domain-level net/abs curve is produced.
        net_flux = arr[:, net_cols].sum(axis=1)
        abs_flux = arr[:, abs_cols].sum(axis=1)
    else:
        # Backward compatibility with old files that only had one boundary flux column.
        net_flux = arr[:, -1]
        abs_flux = np.abs(net_flux)

    df = pd.DataFrame({"time_s": arr[:, 0], "boundary_flux_net_m3s": net_flux, "boundary_flux_abs_m3s": abs_flux})
    return df


def _build_solute_df(path: Path, species: int) -> pd.DataFrame:
    """
    Parse RTSubsurfaceTimeSeries.out for one species.

    The file stores 8 columns per species after the time column:
      LiquidMass, SolidMass, MassExch, BCInflow, BCOutflow,
      SSInflow, SSOutflow, Reaction
    """
    header_tokens = _read_header_tokens(path)
    arr = _read_numeric_table(path)
    n_base = 1
    n_per_species = 8
    start = n_base + species * n_per_species
    stop = start + n_per_species
    if arr.shape[1] < stop:
        raise ValueError(
            f"Requested species {species}, but file has only {arr.shape[1]} columns."
        )
    sp = arr[:, start:stop]
    df = pd.DataFrame(
        {
            "time_s": arr[:, 0],
            "liquid_mass_mg": sp[:, 0],
            "solid_mass_mg": sp[:, 1],
            "mass_exchange_reported_mgs": sp[:, 2],
            "bc_inflow_mgs": sp[:, 3],
            "bc_outflow_mgs": sp[:, 4],
            "ss_inflow_mgs": sp[:, 5],
            "ss_outflow_mgs": sp[:, 6],
            "reaction_mgs": sp[:, 7],
        }
    )
    df["net_boundary_exchange_mgs"] = df["bc_inflow_mgs"] - df["bc_outflow_mgs"]
    # User-requested top-boundary solute exchange proxy.
    df["top_boundary_exchange_mgs"] = df["bc_outflow_mgs"]

    # New RT format: per-species per-boundary net/abs solute flux columns.
    # Example names: spec0_BoundaryFluxNet_0[mg/s], spec0_BoundaryFluxAbs_0[mg/s]
    spec_prefix = f"spec{species}_"
    net_cols = [
        i
        for i, name in enumerate(header_tokens)
        if name.startswith(spec_prefix + "BoundaryFluxNet_")
    ]
    abs_cols = [
        i
        for i, name in enumerate(header_tokens)
        if name.startswith(spec_prefix + "BoundaryFluxAbs_")
    ]
    if net_cols and abs_cols:
        df["boundary_flux_net_mgs"] = arr[:, net_cols].sum(axis=1)
        df["boundary_flux_abs_mgs"] = arr[:, abs_cols].sum(axis=1)
    else:
        # Backward compatibility for older RTSubsurfaceTimeSeries files:
        # approximate net from BCInflow-BCOutflow and abs intensity from BCInflow+BCOutflow.
        df["boundary_flux_net_mgs"] = df["bc_inflow_mgs"] - df["bc_outflow_mgs"]
        df["boundary_flux_abs_mgs"] = df["bc_inflow_mgs"] + df["bc_outflow_mgs"]
    return df


def _plot_timeseries(
    water_df: pd.DataFrame,
    solute_df: Optional[pd.DataFrame],
    fig_dir: Path,
) -> None:
    """
    Create time-series plots requested by the user.

    Water:
      - plot both boundary-flux net and absolute exchange intensity

    Solute:
      - same structure using BCInflow/BCOutflow and top-boundary proxy
    """
    fig_dir.mkdir(parents=True, exist_ok=True)
    t_hr = water_df["time_s"] / 3600.0

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t_hr, water_df["boundary_flux_net_m3s"], label="Boundary flux net", lw=2)
    ax.plot(t_hr, water_df["boundary_flux_abs_m3s"], label="Boundary flux abs", lw=2)
    ax.axhline(0.0, color="k", ls="--", lw=1)
    ax.set_xlabel("Time [h]")
    ax.set_ylabel("Water exchange rate [m3/s]")
    ax.set_title("Boundary Flux Net vs Abs")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / "timeseries_water_exchange.png", dpi=200)
    plt.close(fig)

    if solute_df is not None:
        t_hr_sol = solute_df["time_s"] / 3600.0
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(
            t_hr_sol,
            solute_df["bc_inflow_mgs"],
            label="BCInflow (lateral-dominant)",
            lw=2,
        )
        ax.plot(
            t_hr_sol,
            solute_df["bc_outflow_mgs"],
            label="BCOutflow (top-boundary solute exchange)",
            lw=2,
        )
        ax.plot(
            t_hr_sol,
            solute_df["top_boundary_exchange_mgs"],
            label="Top-boundary solute exchange",
            lw=2,
            ls="--",
        )
        ax.plot(
            t_hr_sol,
            solute_df["mass_exchange_reported_mgs"],
            label="MassExch column (reported)",
            lw=2,
            alpha=0.8,
        )
        ax.axhline(0.0, color="k", ls="--", lw=1)
        ax.set_xlabel("Time [h]")
        ax.set_ylabel("Solute exchange rate [mg/s]")
        ax.set_title("Solute Boundary Flux Components and Top-Boundary Exchange")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(fig_dir / "timeseries_solute_exchange.png", dpi=200)
        plt.close(fig)

        # Second solute figure: net vs absolute boundary exchange intensity.
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(t_hr_sol, solute_df["boundary_flux_net_mgs"], label="Boundary flux net", lw=2)
        ax.plot(t_hr_sol, solute_df["boundary_flux_abs_mgs"], label="Boundary flux abs", lw=2)
        ax.axhline(0.0, color="k", ls="--", lw=1)
        ax.set_xlabel("Time [h]")
        ax.set_ylabel("Solute exchange rate [mg/s]")
        ax.set_title("Solute Boundary Flux Net vs Abs")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(fig_dir / "timeseries_solute_boundary_flux_net_abs.png", dpi=200)
        plt.close(fig)


def _extract_2d(data_array: "xr.DataArray", time_index: int, layer_index: int) -> np.ndarray:
    """Extract a y-x map from arrays shaped either [t,z,y,x] or [t,y,x]."""
    arr = data_array
    if "t" in arr.dims:
        arr = arr.isel(t=time_index)
    if "z" in arr.dims:
        arr = arr.isel(z=layer_index)
    # Remaining dims expected to be y, x
    return np.asarray(arr)


def _print_netcdf_keys(nc_path: Path) -> None:
    """
    Print all available keys in a NetCDF file.

    "Keys" here include:
      - dimensions
      - variables

    We prefer xarray when available; otherwise fallback to netCDF4.
    """
    if not nc_path.exists():
        print(f"[NetCDF keys] File not found: {nc_path}")
        return

    print(f"\n[NetCDF keys] {nc_path.name}")
    if xr is not None:
        with xr.open_dataset(nc_path) as ds:
            print(f"  dimensions: {list(ds.dims.keys())}")
            print(f"  variables:  {list(ds.variables.keys())}")
        return

    if Dataset is not None:
        with Dataset(nc_path, mode="r") as ds:
            print(f"  dimensions: {list(ds.dimensions.keys())}")
            print(f"  variables:  {list(ds.variables.keys())}")
        return

    print("  Cannot inspect keys (xarray/netCDF4 not available).")


def _plot_spatial_maps(
    outdir: Path,
    fig_dir: Path,
    species: int,
    time_index: int,
    layer_index: int,
) -> None:
    if xr is None:
        print("xarray is not available; skipping NetCDF-based spatial plots.")
        return

    fig_dir.mkdir(parents=True, exist_ok=True)
    subs_nc = outdir / "output_subsurface.nc"
    tr_nc = outdir / "output_transport.nc"

    qz_map = None
    c_map = None

    if subs_nc.exists():
        with xr.open_dataset(subs_nc) as ds_sub:
            if "qz" in ds_sub:
                qz_map = _extract_2d(ds_sub["qz"], time_index, layer_index)
                fig, ax = plt.subplots(figsize=(8, 4))
                im = ax.imshow(
                    qz_map,
                    origin="upper",
                    aspect="auto",
                    cmap="RdBu_r",
                    vmin=QZ_VMIN,
                    vmax=QZ_VMAX,
                )
                ax.set_title("Vertical Flux qz (top layer)")
                ax.set_xlabel("x index")
                ax.set_ylabel("y index")
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label("qz [m/s]")
                fig.tight_layout()
                fig.savefig(fig_dir / "map_qz_toplayer.png", dpi=220)
                plt.show()
                plt.close(fig)
            else:
                print(f"{subs_nc.name} does not contain variable 'qz'.")
    else:
        print(f"{subs_nc} not found; skipping water exchange spatial map.")

    if tr_nc.exists():
        with xr.open_dataset(tr_nc) as ds_tr:
            if "c" in ds_tr:
                da = ds_tr["c"]
                if "n" not in da.dims:
                    raise ValueError("Variable 'c' exists but no 'n' species dimension found.")
                c_map = _extract_2d(da.isel(n=species), time_index, layer_index)
                fig, ax = plt.subplots(figsize=(8, 4))
                im = ax.imshow(
                    c_map,
                    origin="upper",
                    aspect="auto",
                    cmap="viridis",
                    vmin=C_VMIN,
                    vmax=C_VMAX,
                )
                ax.set_title(f"Solute Concentration c (species={species}, top layer)")
                ax.set_xlabel("x index")
                ax.set_ylabel("y index")
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label("c [mg/L]")
                fig.tight_layout()
                fig.savefig(fig_dir / f"map_c_species{species}_toplayer.png", dpi=220)
                plt.show()
                plt.close(fig)
            else:
                print(f"{tr_nc.name} does not contain variable 'c'.")
    else:
        print(f"{tr_nc} not found; skipping concentration/solute maps.")

    # Proxy solute exchange flux map: qz * c_top * 1000 [mg/(m2*s)]
    # Sign convention follows qz.
    if qz_map is not None and c_map is not None:
        solute_exchange = qz_map * c_map * 1000.0
        fig, ax = plt.subplots(figsize=(8, 4))
        im = ax.imshow(
            solute_exchange,
            origin="upper",
            aspect="auto",
            cmap="RdBu_r",
            vmin=SOLUTE_EXCHANGE_VMIN,
            vmax=SOLUTE_EXCHANGE_VMAX,
        )
        ax.set_title(f"Proxy Solute Exchange Flux (species={species})")
        ax.set_xlabel("x index")
        ax.set_ylabel("y index")
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("qz*c*1000 [mg/(m2*s)]")
        fig.tight_layout()
        fig.savefig(fig_dir / f"map_solute_exchange_proxy_species{species}.png", dpi=220)
        plt.show()
        plt.close(fig)


def main() -> None:
    """Main workflow: read -> save CSV -> plot -> report."""
    outdir = OUTDIR
    if not outdir.exists():
        raise FileNotFoundError(f"Output directory not found: {outdir}")

    fig_dir = FIGDIR if FIGDIR is not None else outdir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    sub_ts = outdir / "SubsurfaceTimeSeries.out"
    if not sub_ts.exists():
        raise FileNotFoundError(f"Missing required file: {sub_ts}")
    water_df = _build_water_df(sub_ts)
    water_df.to_csv(fig_dir / "water_exchange_timeseries.csv", index=False)

    rt_ts = outdir / "RTSubsurfaceTimeSeries.out"
    solute_df = None
    if rt_ts.exists():
        solute_df = _build_solute_df(rt_ts, species=SPECIES_INDEX)
        solute_df.to_csv(fig_dir / f"solute_exchange_species{SPECIES_INDEX}.csv", index=False)
    else:
        print(f"{rt_ts} not found. Solute time-series plots will be skipped.")

    # Always print all keys available in expected NetCDF outputs.
    _print_netcdf_keys(outdir / "output_subsurface.nc")
    _print_netcdf_keys(outdir / "output_transport.nc")

    _plot_timeseries(water_df, solute_df, fig_dir)
    _plot_spatial_maps(
        outdir=outdir,
        fig_dir=fig_dir,
        species=SPECIES_INDEX,
        time_index=TIME_INDEX,
        layer_index=LAYER_INDEX,
    )

    print(f"Post-processing complete. Figures and CSVs saved to: {fig_dir}")


if __name__ == "__main__":
    main()

