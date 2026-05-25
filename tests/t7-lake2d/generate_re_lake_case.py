#!/usr/bin/env python3
"""
Generate a SERGHEI RE-only lake-groundwater test case.

This script is intentionally configured with hard-coded parameters (no CLI arguments),
so users can edit values directly in the CONFIG block below.

Key features:
- 2D x-z style setup (ny = 1)
- Flat central lake bed and two linearly sloped banks
- Partially inundated top boundary (lake level above bed, below bank crest)
- Fixed lateral groundwater table boundaries
- Homogeneous soil by default
- Optional wind-wave boundary inputs (wave/fetch/wind files) for RE top BC forcing
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple


@dataclass(frozen=True)
class Config:
    # -------------------------------------------------------------------------
    # Output location
    # -------------------------------------------------------------------------
    outdir: Path = Path("../t1-lake2d/input-waveconst")

    # -------------------------------------------------------------------------
    # Grid/domain (requested defaults)
    # -------------------------------------------------------------------------
    nx: int = 120
    ny: int = 1
    nz: int = 50
    dx: float = 0.2
    dz: float = 0.1
    xll: float = 0.0
    yll: float = 0.0
    nodata_value: float = -9999.0

    # -------------------------------------------------------------------------
    # Lake/bank geometry
    # z(x):
    #   left bank (linear down), flat bed (center), right bank (linear up)
    # -------------------------------------------------------------------------
    bed_elev: float = 0.0
    bank_width: float = 5.0
    bank_slope: float = 0.2  # dz/dx [m/m]

    # -------------------------------------------------------------------------
    # Hydrologic levels (must produce partial inundation)
    # -------------------------------------------------------------------------
    lake_level_elev: float = 0.5       # top BC over lake polygon (WT const)
    initial_wt_elev: float = 0.5      # wt.input initial condition
    left_wt_elev: float = 0.6         # left lateral BC
    right_wt_elev: float = 0.6        # right lateral BC
    boundary_side_offset: float = 0.01  # extra x-buffer for left/right boundary polygons

    # -------------------------------------------------------------------------
    # Global controls
    # -------------------------------------------------------------------------
    sim_length: float = 86400.0
    cfl: float = 0.8
    out_freq: float = 3600.0
    obs_freq: float = 3600.0
    n_screen: int = 1000
    par_nx: int = 1
    par_ny: int = 1
    out_format: str = "NETCDF"  # BIN / VTK / NETCDF (BIN recommended)

    # -------------------------------------------------------------------------
    # RE solver controls
    # -------------------------------------------------------------------------
    dz_multiplier: float = 1.0
    dt_init: float = 0.01
    dt_max: float = 30.0
    gw_scheme: int = 1  # 1=PCA, 2=Picard
    cg_iter: int = 1000
    cg_tol: float = 1e-8
    async_coupling: int = 0
    aev: float = 0.0

    # -------------------------------------------------------------------------
    # Soil / van Genuchten (homogeneous by default)
    # Use one value for homogeneous case; semicolon-join logic handled below.
    # -------------------------------------------------------------------------
    n_soil: int = 1
    vg_ks: Tuple[float, ...] = (2.9e-6,)
    vg_phi: Tuple[float, ...] = (0.44,)
    vg_thetas: Tuple[float, ...] = (0.44,)
    vg_thetar: Tuple[float, ...] = (0.1,)
    vg_n: Tuple[float, ...] = (1.38,)
    vg_alpha: Tuple[float, ...] = (0.93,)

    # -------------------------------------------------------------------------
    # Groundwater solute transport (single species)
    # -------------------------------------------------------------------------
    enable_rt: bool = True
    rt_n_mass: int = 1
    rt_initial_aq_conc: float = 1.0  # uniform in RE domain
    rt_diffusion_molecular: float = 1.0e-9
    rt_alpha_L: float = 0.1
    rt_alpha_T: float = 0.01
    rt_scheme: int = 2
    rt_up_weighting_vplus: float = 1.0
    rt_up_weighting_vminus: float = 0.0
    rt_reaction_module: int = 0
    rt_top_bc_conc: float = 0.0      # zero concentration in lake water/top boundary

    # -------------------------------------------------------------------------
    # Optional wave boundary forcing for RE top boundary
    #
    # Goal:
    #   Generate wave.input (+ fetch.input / wind.input) so SERGHEI can
    #   activate the RE wave boundary module (SERGHEI_WAVE_MODEL).
    #
    # Typical setup for a "t1-lake2d-like with wind-wave" test:
    #   - keep geometry and hydrology identical
    #   - set enable_wave = True
    #   - keep wave_mean_lake_level close to lake_level_elev
    #   - use constant wind first, then switch to windfile if needed
    # -------------------------------------------------------------------------
    enable_wave: bool = True

    # waveMethod in wave.input:
    #   "smb"      : SMB-type wave estimate
    #   "fetchlaw" : fetch-law variant (if supported in your branch)
    wave_method: str = "smb"  # smb | fetchlaw

    # phaseMode in wave.input:
    #   "fetchprojection" : project along wind/fetch direction
    #   "none"            : no additional phase handling
    wave_phase_mode: str = "fetchprojection"  # fetchprojection | none

    # meanLakeLevel in wave.input [m]:
    #   Reference lake level used by wave boundary calculations.
    #   A practical default is to match lake_level_elev.
    wave_mean_lake_level: float = 0.5

    # windSource in wave.input:
    #   "constant" : use wave_wind_speed + wave_wind_direction below
    #   "windfile" : read time-varying wind from wind.input generated by script
    wave_wind_source: str = "constant"  # constant | windfile

    # Constant wind parameters (used only when wave_wind_source == "constant"):
    #   wave_wind_speed     [m/s]
    #   wave_wind_direction [deg] (model convention, typically azimuth-like)
    wave_wind_speed: float = 5.0
    wave_wind_direction: float = 90.0

    # fetchSource in wave.input:
    #   "precomputed" : read 2D fetch raster (fetch.input)
    #   "boundarylist": read boundary-cell list values (still written to fetch file
    #                   path; format/usage depends on solver branch)
    wave_fetch_source: str = "precomputed"  # precomputed | boundarylist

    # Name of fetch file referenced by wave.input.
    # Default and recommended: "fetch.input".
    wave_fetch_file: str = "fetch.input"

    # Uniform fetch distance used when writing fetch.input raster [m].
    # Increase for longer effective fetch over open-water stretches.
    wave_fetch_length: float = 10.0

    # wind.input fields (used only when wave_wind_source == "windfile"):
    #   np    : number of wind time records
    #   CwT   : wind coupling coefficient in wind.input header
    #   hwmin : minimum water depth threshold in wind.input header [m]
    # The script writes a simple piecewise-constant (or 2-point linear) series.
    wind_np: int = 2
    wind_cwt: float = 0.0013
    wind_hwmin: float = 0.01


CFG = Config()


def ensure_list_length(values: Tuple[float, ...], n: int, name: str) -> List[float]:
    if len(values) == 1:
        return [values[0]] * n
    if len(values) != n:
        raise ValueError(f"{name} length must be 1 or n_soil={n}, got {len(values)}")
    return list(values)


def semicolon(values: List[float]) -> str:
    return ";".join(f"{v:.10g}" for v in values)


def bed_profile(cfg: Config) -> List[float]:
    length = cfg.nx * cfg.dx
    if cfg.bank_width <= 0.0:
        return [cfg.bed_elev] * cfg.nx
    if 2.0 * cfg.bank_width > length:
        raise ValueError("bank_width is too large: 2 * bank_width must be <= domain length")

    crest = cfg.bed_elev + cfg.bank_slope * cfg.bank_width
    z = []
    for i in range(cfg.nx):
        x = cfg.xll + (i + 0.5) * cfg.dx
        xr = x - cfg.xll
        if xr < cfg.bank_width:
            frac = xr / cfg.bank_width
            zi = crest - (crest - cfg.bed_elev) * frac
        elif xr <= (length - cfg.bank_width):
            zi = cfg.bed_elev
        else:
            frac = (xr - (length - cfg.bank_width)) / cfg.bank_width
            zi = cfg.bed_elev + (crest - cfg.bed_elev) * frac
        z.append(zi)
    return z


def find_lake_span(bed: List[float], lake_level: float) -> Tuple[int, int]:
    wet = [i for i, z in enumerate(bed) if z < lake_level]
    if not wet:
        raise ValueError("No inundated cells: increase lake_level_elev above bed elevation.")
    i0, i1 = wet[0], wet[-1]
    if i0 == 0 and i1 == len(bed) - 1:
        raise ValueError(
            "Domain fully inundated. Set lake_level_elev lower than highest bank crest."
        )
    return i0, i1


def write_polygon(path: Path, pts: List[Tuple[float, float]]) -> None:
    lines = [f"NPOINTS {len(pts)}"] + [f"{x:.10g} {y:.10g}" for x, y in pts]
    path.write_text("\n".join(lines) + "\n")


def build_top_polygons(cfg: Config, lake_i0: int, lake_i1: int) -> dict:
    """Create top polygons for lake and dry-bank segments."""
    width = cfg.ny * cfg.dx
    m = 0.25 * cfg.dx
    x0 = cfg.xll
    x1 = cfg.xll + cfg.nx * cfg.dx
    xl0 = cfg.xll + lake_i0 * cfg.dx
    xl1 = cfg.xll + (lake_i1 + 1) * cfg.dx

    polys = {}
    polys["lake"] = [
        (xl0 - m, cfg.yll - m),
        (xl0 - m, cfg.yll + width + m),
        (xl1 + m, cfg.yll + width + m),
        (xl1 + m, cfg.yll - m),
    ]

    if lake_i0 > 0:
        polys["top_left_dry"] = [
            (x0 - m, cfg.yll - m),
            (x0 - m, cfg.yll + width + m),
            (xl0 + m, cfg.yll + width + m),
            (xl0 + m, cfg.yll - m),
        ]
    if lake_i1 < cfg.nx - 1:
        polys["top_right_dry"] = [
            (xl1 - m, cfg.yll - m),
            (xl1 - m, cfg.yll + width + m),
            (x1 + m, cfg.yll + width + m),
            (x1 + m, cfg.yll - m),
        ]

    # side polygons
    side_offset = cfg.boundary_side_offset
    polys["left_side"] = [
        (x0 - m, cfg.yll - m),
        (x0 - m, cfg.yll + width + m),
        (x0 + m + side_offset, cfg.yll + width + m),
        (x0 + m + side_offset, cfg.yll - m),
    ]
    polys["right_side"] = [
        (x1 - m - side_offset, cfg.yll - m),
        (x1 - m - side_offset, cfg.yll + width + m),
        (x1 + m, cfg.yll + width + m),
        (x1 + m, cfg.yll - m),
    ]
    polys["full_top"] = [
        (x0 - m, cfg.yll - m),
        (x0 - m, cfg.yll + width + m),
        (x1 + m, cfg.yll + width + m),
        (x1 + m, cfg.yll - m),
    ]
    return polys


def write_parameters(path: Path, cfg: Config) -> None:
    text = (
        "// Global run control\n"
        f"simLength : {cfg.sim_length}\n"
        f"cfl : {cfg.cfl}\n"
        f"outFreq : {cfg.out_freq}\n"
        f"obsFreq : {cfg.obs_freq}\n"
        f"nScreen : {cfg.n_screen}\n"
        f"parNx : {cfg.par_nx}\n"
        f"parNy : {cfg.par_ny}\n"
        f"outFormat : {cfg.out_format}\n"
        "BCtype : REFLECTIVE\n"
    )
    path.write_text(text)


def write_dem(path: Path, cfg: Config, bed: List[float]) -> None:
    lines = [
        f"ncols {cfg.nx}",
        f"nrows {cfg.ny}",
        f"xllcorner {cfg.xll}",
        f"yllcorner {cfg.yll}",
        f"cellsize {cfg.dx}",
        f"nodata_value {cfg.nodata_value}",
        " ".join(f"{z:.10g}" for z in bed),
    ]
    path.write_text("\n".join(lines) + "\n")


def write_subsurface(path: Path, cfg: Config) -> None:
    text = (
        "// RE solver / grid control\n"
        f"height : {cfg.nz * cfg.dz}\n"
        f"dz_multiplier : {cfg.dz_multiplier}\n"
        f"dz_base : {cfg.dz}\n"
        f"ndepth : {cfg.nz}\n"
        f"nSoilID : {cfg.n_soil}\n"
        "initialMode : 4\n"
        f"dt_init : {cfg.dt_init}\n"
        f"dt_max : {cfg.dt_max}\n"
        f"gw_scheme : {cfg.gw_scheme}\n"
        f"cg_iter : {cfg.cg_iter}\n"
        f"cg_tol : {cfg.cg_tol}\n"
        f"async : {cfg.async_coupling}\n"
        f"aev : {cfg.aev}\n"
    )
    path.write_text(text)


def write_soil_id(path: Path, cfg: Config) -> None:
    ndata = cfg.nx * cfg.ny * cfg.nz
    if cfg.n_soil == 1:
        path.write_text("n_soil 1\n0\n")
    else:
        path.write_text(f"n_soil {cfg.n_soil}\n" + " ".join(["0"] * ndata) + "\n")


def write_vg(path: Path, cfg: Config) -> None:
    ks = ensure_list_length(cfg.vg_ks, cfg.n_soil, "vg_ks")
    phi = ensure_list_length(cfg.vg_phi, cfg.n_soil, "vg_phi")
    ts = ensure_list_length(cfg.vg_thetas, cfg.n_soil, "vg_thetas")
    tr = ensure_list_length(cfg.vg_thetar, cfg.n_soil, "vg_thetar")
    n = ensure_list_length(cfg.vg_n, cfg.n_soil, "vg_n")
    a = ensure_list_length(cfg.vg_alpha, cfg.n_soil, "vg_alpha")

    text = (
        "// van Genuchten parameters\n"
        f"alpha : {semicolon(a)}\n"
        f"n : {semicolon(n)}\n"
        f"Ks : {semicolon(ks)}\n"
        f"Phi : {semicolon(phi)}\n"
        f"ThetaR : {semicolon(tr)}\n"
        f"ThetaS : {semicolon(ts)}\n"
    )
    path.write_text(text)


def write_wt(path: Path, cfg: Config) -> None:
    vals = " ".join([f"{cfg.initial_wt_elev:.10g}"] * (cfg.nx * cfg.ny))
    lines = [
        f"ncols {cfg.nx}",
        f"nrows {cfg.ny}",
        f"nodata_value {cfg.nodata_value}",
        vals,
    ]
    path.write_text("\n".join(lines) + "\n")


def write_gwbc(path: Path, cfg: Config, has_left_dry: bool, has_right_dry: bool) -> None:
    blocks = []

    # lake segment at top (partially inundated)
    blocks.append(
        (
            "lake_top",
            6,
            4,  # SUB_BC_WT_CONST
            cfg.lake_level_elev,
            "polygonlake.input",
        )
    )

    # lateral fixed groundwater table boundaries
    blocks.append(("left", 2, 4, cfg.left_wt_elev, "polygonleft.input"))
    blocks.append(("right", 1, 4, cfg.right_wt_elev, "polygonright.input"))

    lines = [f"bccount : {len(blocks)}", ""]
    for bid, direction, bctype, bcval, poly in blocks:
        lines.append(f"id : {bid}")
        lines.append(f"direction : {direction}")
        lines.append(f"bctype : {bctype}")
        if bcval is not None:
            lines.append(f"bcvals : {bcval}")
        lines.append(f"polygon : {poly}")
        lines.append("")
    path.write_text("\n".join(lines))


def write_readme(path: Path, cfg: Config, crest: float, lake_i0: int, lake_i1: int) -> None:
    length = cfg.nx * cfg.dx
    lake_x0 = cfg.xll + lake_i0 * cfg.dx
    lake_x1 = cfg.xll + (lake_i1 + 1) * cfg.dx
    text = (
        "Generated SERGHEI RE-only lake-groundwater case (partially inundated)\n"
        "---------------------------------------------------------------------\n"
        f"nx={cfg.nx}, ny={cfg.ny}, nz={cfg.nz}\n"
        f"dx={cfg.dx} m, dz={cfg.dz} m, domain length={length} m\n"
        f"bed elevation={cfg.bed_elev} m, bank crest={crest} m\n"
        f"lake level={cfg.lake_level_elev} m (must satisfy bed < lake < crest)\n"
        f"lake top span in x: [{lake_x0}, {lake_x1}] m\n"
        f"left WT={cfg.left_wt_elev} m, right WT={cfg.right_wt_elev} m\n"
        f"initial WT={cfg.initial_wt_elev} m\n"
        f"left/right side polygon x-offset={cfg.boundary_side_offset} m\n"
        f"RT enabled={cfg.enable_rt}, RT initial concentration={cfg.rt_initial_aq_conc}\n"
        f"RT top boundary concentration (lake)={cfg.rt_top_bc_conc}\n"
        f"Wave enabled={cfg.enable_wave}, windSource={cfg.wave_wind_source}, "
        f"fetchSource={cfg.wave_fetch_source}\n"
        "\n"
        "Run example:\n"
        "  <serghei_bin> <this_folder/> <new_output_folder/> <n_threads>\n"
    )
    path.write_text(text)


def write_rttransportgw(path: Path, cfg: Config) -> None:
    text = (
        "// Groundwater reactive transport setup (single species)\n"
        f"n_mass : {cfg.rt_n_mass}\n"
        "id : 0\n"
        "aq_mode : 0\n"
        f"aq_val : {cfg.rt_initial_aq_conc}\n"
        f"diffusion_molecular : {cfg.rt_diffusion_molecular}\n"
        f"alpha_L : {cfg.rt_alpha_L}\n"
        f"alpha_T : {cfg.rt_alpha_T}\n"
        f"rt_scheme : {cfg.rt_scheme}\n"
        f"Up_Weighting_vplus : {cfg.rt_up_weighting_vplus}\n"
        f"Up_Weighting_vminus : {cfg.rt_up_weighting_vminus}\n"
        f"ReactionModule : {cfg.rt_reaction_module}\n"
    )
    path.write_text(text)


def write_rtgwbc(path: Path, cfg: Config) -> None:
    # left/right are fixed (Dirichlet) concentration boundaries equal to
    # the initial groundwater concentration; top (lake) is fixed zero concentration.
    text = (
        "bccount : 3\n"
        "\n"
        "id : rtleft\n"
        "direction : 2\n"
        "rtbctype : 2\n"
        f"bcvals : {cfg.rt_initial_aq_conc}\n"
        "polygon : polygonleft.input\n"
        "spec_0_bctype : 2\n"
        f"spec_0_bcval : {cfg.rt_initial_aq_conc}\n"
        "\n"
        "id : rtright\n"
        "direction : 1\n"
        "rtbctype : 2\n"
        f"bcvals : {cfg.rt_initial_aq_conc}\n"
        "polygon : polygonright.input\n"
        "spec_0_bctype : 2\n"
        f"spec_0_bcval : {cfg.rt_initial_aq_conc}\n"
        "\n"
        "id : rtlake_top\n"
        "direction : 6\n"
        "rtbctype : 2\n"
        f"bcvals : {cfg.rt_top_bc_conc}\n"
        "polygon : polygonlake.input\n"
        "spec_0_bctype : 2\n"
        f"spec_0_bcval : {cfg.rt_top_bc_conc}\n"
    )
    path.write_text(text)


def write_wave(path: Path, cfg: Config) -> None:
    # Write wave.input consumed by Parser::readWaveFile().
    lines = [
        "// Wave boundary module control",
        f"waveEnable : {1 if cfg.enable_wave else 0}",
        f"waveMethod : {cfg.wave_method}",
        f"windSource : {cfg.wave_wind_source}",
        f"fetchSource : {cfg.wave_fetch_source}",
        f"fetchFile : {cfg.wave_fetch_file}",
        f"phaseMode : {cfg.wave_phase_mode}",
        f"meanLakeLevel : {cfg.wave_mean_lake_level}",
    ]
    if cfg.wave_wind_source == "constant":
        lines.append(f"windSpeed : {cfg.wave_wind_speed}")
        lines.append(f"windDirection : {cfg.wave_wind_direction}")
    path.write_text("\n".join(lines) + "\n")


def write_fetch(path: Path, cfg: Config) -> None:
    # Write fetch.input as a DEM-like raster:
    # one fetch value per surface cell (here uniform for simplicity).
    vals = " ".join([f"{cfg.wave_fetch_length:.10g}"] * (cfg.nx * cfg.ny))
    lines = [
        f"ncols {cfg.nx}",
        f"nrows {cfg.ny}",
        f"xllcorner {cfg.xll}",
        f"yllcorner {cfg.yll}",
        f"cellsize {cfg.dx}",
        f"nodata_value {cfg.nodata_value}",
        vals,
    ]
    path.write_text("\n".join(lines) + "\n")


def write_wind(path: Path, cfg: Config) -> None:
    # Write wind.input:
    # time [s], wind speed [m/s], wind direction [deg].
    # For wind_np > 1, times are evenly spaced over sim_length.
    lines = [
        f"np {cfg.wind_np}",
        f"CwT {cfg.wind_cwt}",
        f"hwmin {cfg.wind_hwmin}",
    ]
    if cfg.wind_np <= 1:
        lines.append(f"0 {cfg.wave_wind_speed} {cfg.wave_wind_direction}")
    else:
        dt = cfg.sim_length / float(cfg.wind_np - 1)
        for i in range(cfg.wind_np):
            t = i * dt
            lines.append(f"{t:.10g} {cfg.wave_wind_speed:.10g} {cfg.wave_wind_direction:.10g}")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    cfg = CFG
    if cfg.nx <= 0 or cfg.ny <= 0 or cfg.nz <= 0:
        raise ValueError("nx, ny, nz must be positive")
    if cfg.dx <= 0.0 or cfg.dz <= 0.0:
        raise ValueError("dx and dz must be positive")
    if cfg.boundary_side_offset < 0.0:
        raise ValueError("boundary_side_offset must be >= 0")
    if cfg.wave_wind_source not in {"constant", "windfile"}:
        raise ValueError("wave_wind_source must be 'constant' or 'windfile'")
    if cfg.wave_fetch_source not in {"precomputed", "boundarylist"}:
        raise ValueError("wave_fetch_source must be 'precomputed' or 'boundarylist'")

    bed = bed_profile(cfg)
    crest = cfg.bed_elev + cfg.bank_slope * cfg.bank_width
    bed_min = min(bed)
    if not (bed_min < cfg.lake_level_elev < crest):
        raise ValueError(
            f"Need partial inundation: bed_min({bed_min}) < lake_level({cfg.lake_level_elev}) < crest({crest})."
        )
    lake_i0, lake_i1 = find_lake_span(bed, cfg.lake_level_elev)

    outdir = cfg.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # Clean up stale files from older script versions
    for stale in ["polygontop_left_dry.input", "polygontop_right_dry.input"]:
        p = outdir / stale
        if p.exists():
            p.unlink()

    polys = build_top_polygons(cfg, lake_i0, lake_i1)
    write_polygon(outdir / "polygonlake.input", polys["lake"])
    write_polygon(outdir / "polygonleft.input", polys["left_side"])
    write_polygon(outdir / "polygonright.input", polys["right_side"])
    write_polygon(outdir / "polygontop.input", polys["full_top"])

    write_parameters(outdir / "parameters.input", cfg)
    write_dem(outdir / "dem.input", cfg, bed)
    write_subsurface(outdir / "subsurface.input", cfg)
    write_soil_id(outdir / "soilID.input", cfg)
    write_vg(outdir / "vg.input", cfg)
    write_wt(outdir / "wt.input", cfg)
    write_gwbc(
        outdir / "gwbc.input",
        cfg,
        has_left_dry=("top_left_dry" in polys),
        has_right_dry=("top_right_dry" in polys),
    )
    if cfg.enable_rt:
        write_rttransportgw(outdir / "rttransportgw.input", cfg)
        write_rtgwbc(outdir / "rtgwbc.input", cfg)
    if cfg.enable_wave:
        # Generate wave forcing files only when requested.
        # This keeps a baseline RE-only case uncluttered by optional inputs.
        write_wave(outdir / "wave.input", cfg)
        write_fetch(outdir / cfg.wave_fetch_file, cfg)
        if cfg.wave_wind_source == "windfile":
            write_wind(outdir / "wind.input", cfg)
    write_readme(outdir / "README_generated_case.txt", cfg, crest, lake_i0, lake_i1)

    print(f"Generated RE-only partially inundated case in: {outdir}")
    print("Edit hard-coded values in CONFIG class if needed.")


if __name__ == "__main__":
    main()
