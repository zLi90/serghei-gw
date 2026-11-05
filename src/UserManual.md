# SERGHEI User Manual — Coupled Surface–Subsurface Hydrodynamics, Solute Transport & WOFOST Crop Growth

> **Scope**: This manual covers the **surface-water hydrodynamics** module ([Section 7](#7-surface-water-hydrodynamics-swe)), the **groundwater hydrodynamics** (Richards-equation subsurface flow) module ([Section 8](#8-groundwater-hydrodynamics-subsurface-flow)), the **solute (reactive) transport** module (surface water and groundwater) and the **WOFOST 7.2 crop growth** module within SERGHEI, together with the **triggered irrigation / triggered drainage** surface boundary types (`bctype` 13/14 in `extbc.input`).
>
> **Reference test cases**: `tests/t9-gaoyou-Ncycle-WOFOST/` (paddy rice with nitrogen cycle), `tests/t15-gaoyou-wofost-triggerIrrigation/` (rule-based irrigation and drainage management, depth triggers), `tests/t16-gaoyou-wofost-triggerIrrigation/` (W1 shallow-wet-dry schedule with a root-zone soil-moisture trigger stage)

---

## Table of Contents

- [1. Surface Water Solute Transport](#1-surface-water-solute-transport)
  - [1.1 Transport Parameters — rttransportsw.input](#11-transport-parameters--rttransportswinput)
  - [1.2 Boundary Conditions — rtswbc.input](#12-boundary-conditions--rtswbcinput)
  - [1.3 Source/Sink Terms — rtswss.input](#13-sourcesink-terms--rtswssinput)
  - [1.4 Reaction Parameters — reactionsw.input](#14-reaction-parameters--reactionswinput)
- [2. Groundwater Solute Transport](#2-groundwater-solute-transport)
  - [2.1 Transport Parameters — rttransportgw.input](#21-transport-parameters--rttransportgwinput)
  - [2.2 Boundary Conditions — rtgwbc.input](#22-boundary-conditions--rtgwbcinput)
  - [2.3 Reaction Parameters — reactiongw.input](#23-reaction-parameters--reactiongwinput)
  - [2.4 Groundwater Source/Sink — rtgwss.input](#24-groundwater-sourcesink--rtgwssinput)
  - [2.5 Initial Condition Files](#25-initial-condition-files)
- [3. WOFOST Crop Growth Model](#3-wofost-crop-growth-model)
  - [3.1 Crop Parameters — cropparameter.input](#31-crop-parameters--cropparameterinput)
  - [3.2 Meteorological Data — cropmeteo.input](#32-meteorological-data--cropmeteoinput)
  - [3.3 Agricultural Management — agro.input](#33-agricultural-management--agroinput)
- [3.4 WOFOST–SERGHEI Coupling Interface](#34-wofostserghei-coupling-interface)
- [4. Triggered Irrigation and Drainage Boundaries](#4-triggered-irrigation-and-drainage-boundaries)
  - [4.1 Overview](#41-overview)
  - [4.2 Boundary Definition — extbc.input](#42-boundary-definition--extbcinput)
  - [4.3 Irrigation Strategy File](#43-irrigation-strategy-file)
  - [4.4 Drainage Strategy File](#44-drainage-strategy-file)
  - [4.5 Trigger Polygons](#45-trigger-polygons)
  - [4.6 Runtime Behaviour and Output](#46-runtime-behaviour-and-output)
- [5. Timeseries File Format](#5-timeseries-file-format)
- [6. Output Files and Headers](#6-output-files-and-headers)
- [7. Surface Water Hydrodynamics (SWE)](#7-surface-water-hydrodynamics-swe)
  - [7.1 Simulation Parameters — parameters.input](#71-simulation-parameters--parametersinput)
  - [7.2 Terrain — dem.input](#72-terrain--deminput)
  - [7.3 Initial State and Friction — sw.input](#73-initial-state-and-friction--swinput)
  - [7.4 External Boundaries — extbc.input](#74-external-boundaries--extbcinput)
  - [7.5 Rainfall — rainfall.input (optional)](#75-rainfall--rainfallinput-optional)
  - [7.6 Infiltration — infiltration.input (optional)](#76-infiltration--infiltrationinput-optional)
  - [7.7 Wind — wind.input (optional)](#77-wind--windinput-optional)
- [8. Groundwater Hydrodynamics (Subsurface Flow)](#8-groundwater-hydrodynamics-subsurface-flow)
  - [8.1 Domain and Solver — subsurface.input](#81-domain-and-solver--subsurfaceinput)
  - [8.2 Soil Hydraulic Parameters — vg.input](#82-soil-hydraulic-parameters--vginput)
  - [8.3 Soil Map — soilID.input](#83-soil-map--soilidinput)
  - [8.4 Boundary Conditions — gwbc.input](#84-boundary-conditions--gwbcinput)
  - [8.5 Source/Sink Zones — gwss.input](#85-sourcesink-zones--gwssinput)
  - [8.6 Initial Conditions](#86-initial-conditions)
  - [8.7 Surface–Subsurface Coupling (summary)](#87-surfacesubsurface-coupling-summary)
- [Appendix A: Species Ordering Convention](#appendix-a-species-ordering-convention)
- [Appendix B: Input File Quick-Reference Checklist](#appendix-b-input-file-quick-reference-checklist)

---

## 1. Surface Water Solute Transport

### 1.1 Transport Parameters — `rttransportsw.input`

**File path**: `<input_folder>/rttransportsw.input`
**Source code**: `src/RTInitSW.h` → `readRTFileSW()`

This file defines the number of solute species, per-species transport properties, initial conditions, and global solver settings for surface water reactive transport.

#### File Format

Colon-separated key–value pairs. Per-species parameters are grouped under an `id` selector. Lines starting with `//` are comments.

#### Global Parameters

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `n_mass` | int | **Yes** | Total number of solute species (chemical components). Must be the first species-related key in the file. |
| `ReactionModule` | int | **Yes** | Reaction module switch: `0` = disabled, `1` = enabled (requires `reactionsw.input`). |
| `Advection_Scheme` | int | **Yes** | Numerical scheme for advective transport: |

| `Advection_Scheme` Value | Name | Description |
|:---:|------|-------------|
| 1 | 1st-order Upwind | First-order upstream weighting. Most diffusive, most stable. |
| 2 | TVD Van Leer | Total Variation Diminishing with Van Leer flux limiter. Recommended for most applications. |
| 3 | TVD Superbee | Total Variation Diminishing with Superbee flux limiter. Sharpest fronts, may produce oscillations near discontinuities. |

#### Per-Species Parameters

Each species block begins with `id : <species_index>` (0-based). All subsequent keys apply to that species until the next `id` key.

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `aq_mode` | int | **Yes** | Initial condition mode for aqueous phase: `0` = constant value (`aq_val`), `1` = spatial file (`aq_file`). |
| `aq_val` | real | If `aq_mode=0` | Uniform initial concentration [mg/L]. The parser does not convert from another concentration unit. |
| `aq_file` | string | If `aq_mode=1` | Filename for spatially distributed initial concentration (ESRI ASCII Grid or simplified format; see [Section 2.5](#25-initial-condition-files)). |
| `diffusion_molecular` | real | **Yes** | Molecular diffusion coefficient [m²/s]. Typical range: 1.0×10⁻¹⁰ to 1.0×10⁻⁸. |
| `alpha_L` | real | **Yes** | Longitudinal dispersivity [m]. Typical range: 0.01 to 100. |
| `alpha_T` | real | **Yes** | Transverse dispersivity [m]. Typical range: 0.001 to 10. Often 1/10 of `alpha_L`. |
| `RainCon` | real | **Yes** | Solute concentration in rainfall for this species [mg/L]. Set to `0.0` if rainfall carries no solute. |

#### Example

```text
// --- Reactive Transport Parameters ---
n_mass : 5        // Total number of species (NH4, NO3, DO, DOC, DON)

// ========== Species 0: NH4 ==========
id : 0
aq_mode : 1
aq_file : rtinitialswNH4.input
RainCon : 0.5                           // [mg/L]
diffusion_molecular : 1.76e-9           // [m^2/s]
alpha_L : 0.5                           // [m]
alpha_T : 0.05                          // [m]

// ========== Species 1: NO3 ==========
id : 1
aq_mode : 1
aq_file : rtinitialswNO3.input
RainCon : 0.4                           // [mg/L]
diffusion_molecular : 1.76e-9           // [m^2/s]
alpha_L : 0.5                           // [m]
alpha_T : 0.05                          // [m]

// ========== Species 2: DO ==========
id : 2
aq_mode : 0
aq_val : 5.0                            // [mg/L]
RainCon : 4.5                           // [mg/L]
diffusion_molecular : 1.76e-9           // [m^2/s]
alpha_L : 0.5                           // [m]
alpha_T : 0.05                          // [m]

// ... species 3 (DOC), 4 (DON) follow the same pattern ...

ReactionModule : 1
Advection_Scheme : 2                    // TVD Van Leer
```

---

### 1.2 Boundary Conditions — `rtswbc.input`

**File path**: `<input_folder>/rtswbc.input`
**Source code**: `src/RTInitSW.h` → `readRTBCFileSW()` (parsing and pairing), `src/BC.h` → `applyrtbc()` (gate resolution and application), `src/RTFunctionSW.h` (per-transport-half-step sync of the irrigation/drainage on/off state)

Defines solute boundary conditions at external boundaries for each species. Boundary geometry (cell indices) is inherited from the hydrodynamic boundary conditions (`extbc.input`).

> **Pairing rule**: boundary blocks in `rtswbc.input` are paired with the boundaries in `extbc.input` **by order** (k-th `id:` block ↔ k-th hydrodynamic boundary); the name after `id:` is not parsed. Gated types 5/6 rely on this pairing to read the hydrodynamic boundary type and its runtime on/off state — **a mismatched block order silently breaks the gating**.

#### Boundary Condition Types

| Value | Constant Name | Description |
|:---:|------|-------------|
| 1 | `RTSW_BC_NONE` | No reactive transport BC applied at this boundary. |
| 2 | `RTSW_BC_CON_CONST` | **Constant concentration** (Dirichlet). The boundary cells are held at a fixed concentration specified by `spec_X_bcval`. |
| 3 | `RTSW_BC_ZEROGRAD` | **Zero-gradient** (Neumann, free outflow). The concentration at the boundary is extrapolated from the interior (zero concentration gradient normal to the boundary). Used for free-discharge outlets. |
| 4 | `RTSW_BC_CON_T` | **Time-varying concentration** (Dirichlet time series). The boundary concentration varies in time according to an external timeseries file specified by `spec_X_bcfile`. See [Section 5](#5-timeseries-file-format) for timeseries format. |
| 5 | `RTSW_BC_TRIGGERED_IRRIGATION` | **Trigger-gated Dirichlet**, pairs with hydraulic `bctype : 13`. While irrigation is *active*, the boundary behaves as type 2/4 (Dirichlet at `spec_X_bcval`, or at the `spec_X_bcfile` timeseries if given — `bcfile` takes precedence); while *inactive*, it falls back to zero-gradient so the boundary concentration follows the interior water and no solute is forced in or out. The advective inflow term of the mass budget already captures the delivered solute while active (the transport kernel advects the Dirichlet value in through the ghost cells), so no extra booking is needed. |
| 6 | `RTSW_BC_TRIGGERED_DRAINAGE` | **Zero-gradient outflow**, explicit type for the triggered drainage boundary (`bctype : 14`). Identical to type 3: solute only leaves while drainage actually discharges; nothing is forced while the outlet is closed. |

#### File Format

```text
bccount : <number_of_boundaries>

id : <boundary_name>            // name not parsed; k-th id block matches the k-th boundary in extbc.input
direction : <nx> <ny>           // Normal vector (will be normalized to unit length)
rtbctype : <type>               // Boundary-level default type (optional; inherited by species without spec_X_bctype)

// Per-species BC (X = species index, 0-based)
spec_X_bctype : <type>          // BC type (1–6, see table above)
spec_X_bcval  : <value>         // Constant concentration [mg/L] (for types 2 and 5)
spec_X_bcfile : <filename>      // Timeseries file (for types 4 and, optionally, 5)
```

#### Cross-Validation Rules

- **Type 2** (constant): Use `spec_X_bcval`. Do **not** specify `spec_X_bcfile`.
- **Type 4** (time-varying): Use `spec_X_bcfile`. Do **not** specify `spec_X_bcval`.
- **Type 5** (triggered irrigation): either `spec_X_bcval` (constant while active) or `spec_X_bcfile` (timeseries while active; when both are given the timeseries takes precedence and `bcval` is silently ignored — provide one or the other). For a non-triggered hydraulic boundary the flag stays off and type 5 degenerates to zero-gradient.
- **Timeseries loading condition**: `spec_X_bcfile` is loaded **only when `spec_X_bctype` is 4 or 5**. If only the boundary-level `rtbctype : 5` is given and the `spec_X_bctype` line is omitted, the timeseries is silently ignored and the irrigation degenerates to clear water (Dirichlet 0) without any error. Always write `spec_X_bctype : 5` explicitly for species with a `bcfile`.
- **Avoid bare type 4 on triggered hydraulic boundaries**: type 4 is not gated and keeps imposing Dirichlet while irrigation is off (a concentration source/sink with no water inflow, not captured by the boundary-flux accounting). Use type 5 for fertilizer-loaded irrigation water.
- The parser will reject mismatched combinations.

#### Example — Inlet (time-varying) + Outlet (free outflow)

```text
bccount : 2

// --- Boundary 0: irrigation inlet ---
id : rtirrigation
direction : 0.0 1.0
spec_0_bctype : 4                        // NH4: time-varying
spec_0_bcfile : rtirrigationNH4.input
spec_1_bctype : 4                        // NO3: time-varying
spec_1_bcfile : rtirrigationNO3.input
spec_2_bctype : 4                        // DO: time-varying
spec_2_bcfile : rtirrigationDO.input
spec_3_bctype : 4                        // DOC: time-varying
spec_3_bcfile : rtirrigationDOC.input
spec_4_bctype : 4                        // DON: time-varying
spec_4_bcfile : rtirrigationDON.input

// --- Boundary 1: outlet ---
id : rtoutlet
direction : 0.0 1.0
spec_0_bctype : 3                        // NH4: free outflow
spec_1_bctype : 3                        // NO3: free outflow
spec_2_bctype : 3                        // DO: free outflow
spec_3_bctype : 3                        // DOC: free outflow
spec_4_bctype : 3                        // DON: free outflow
```

#### Example — Gated solute BCs for triggered irrigation/drainage (types 5/6, paired with hydraulic `bctype : 13/14`)

When pairing with the triggered hydraulic boundaries ([Section 4](#4-triggered-irrigation-and-drainage-boundaries)), the recommended solute types are **5 (irrigation)** and **6 (drainage)**. The semantics of the timeseries under gating deserve emphasis:

> **Meaning of the `spec_X_bcfile` timeseries under type 5**: it is the **concentration of the irrigation water itself as a function of calendar time** (seconds since simulation start) — *not* a schedule of when irrigation/fertilization happens. When and how much water enters is decided dynamically by the hydraulic trigger at every time step; the timeseries only answers "*if* water were flowing at time t, what concentration would it carry". While irrigation is **active**, the boundary imposes a Dirichlet at the concentration interpolated at the current time, and the actual inflow flux = Q(t)·c(t) — the delivered mass of each real trigger episode emerges from integrating over its actual duration. While irrigation is **inactive**, the series is not read at all and the boundary falls back to zero-gradient. **Trigger times therefore do not need to be known in advance**: write the series as "fertilization windows" — an episode triggered inside a window carries the window's concentration; an episode triggered outside carries clear water (series value 0).

**Scenario A — fertilizer concentration varying in time (timeseries)**

`rtswbc.input` (block order matching `extbc.input`: block 1 pairs with `bctype : 13`, block 2 with `bctype : 14`):

```text
// --- Boundary 0: triggered irrigation (hydraulic bctype : 13) ---
id : rtirrigation
direction : 0.0 1.0
rtbctype : 5                 // gated: watering = Dirichlet (from series), no watering = zero-gradient

spec_0_bctype : 5            // NH4: must be written explicitly, otherwise the bcfile below is not loaded
spec_0_bcfile : rtirrigationNH4.input
spec_1_bctype : 5            // NO3
spec_1_bcfile : rtirrigationNO3.input
spec_2_bctype : 5            // DO
spec_2_bcfile : rtirrigationDO.input
spec_3_bctype : 5            // DOC
spec_3_bcfile : rtirrigationDOC.input
spec_4_bctype : 5            // DON: clear water (no bcval/bcfile -> Dirichlet 0 while active)

// --- Boundary 1: triggered drainage (hydraulic bctype : 14) ---
id : rtoutlet
direction : 0.0 1.0
rtbctype : 6                 // always zero-gradient outflow; one key suffices, all species inherit (3 is fully equivalent)
```

`rtirrigationNH4.input` (format in [Section 5](#5-timeseries-file-format); **no comments allowed in data files**):

```text
ndata 14
0.0         0.0
489480.0    0.0
489600.0    0.768
497040.0    0.768
497200.0    0.0
835080.0    0.0
835200.0    0.639
842400.0    0.639
842600.0    0.0
1270440.0   0.0
1270740.0   1.348
1290540.0   1.348
1290600.0   0.0
10886400.0  0.0
```

This defines three "fertilization windows" (0.768 / 0.639 / 1.348 mg/L at ~DAT 5.7 d, 9.7 d, 14.7 d; 126 d = 10 886 400 s). An irrigation episode triggered inside a window delivers water at that window's concentration; outside the windows, clear water; while irrigation is off, zero-gradient with no solute input at all.

**Scenario B — constant irrigation-water concentration**

```text
// --- Boundary 0: triggered irrigation (hydraulic bctype : 13) ---
id : rtirrigation
direction : 0.0 1.0
rtbctype : 5

spec_0_bctype : 5
spec_0_bcval : 0.768     // NH4  mg/L (fill with measured fertigation-water quality)
spec_1_bctype : 5
spec_1_bcval : 0.10      // NO3
spec_2_bctype : 5
spec_2_bcval : 8.0       // DO (near-saturated dissolved oxygen)
spec_3_bctype : 5
spec_3_bcval : 2.0       // DOC
spec_4_bctype : 5        // DON: clear water

// --- Boundary 1: triggered drainage (hydraulic bctype : 14) ---
id : rtoutlet
direction : 0.0 1.0
rtbctype : 6
```

Every irrigation episode over the whole season delivers solute at the same concentration; the delivered mass of each episode is still set by that episode's discharge and duration.

**Notes**

1. Do not give both `bcfile` and `bcval` for the same species: the timeseries wins and `bcval` is silently ignored.
2. Species with a `bcfile` must carry an explicit `spec_X_bctype : 5` (loading condition, see cross-validation rules).
3. Start the series with a `t = 0` row and end it at the simulation length — beyond the last point the value is **held constant**, so a nonzero last value would stay in effect until the end of the run.
4. Verification: after a short run, inspect `outRTSurfaceTimeSeries.out` in the input folder — per-species `_BCInflow[mg/s]` should be 0 while irrigation is off, `_BCInflow[mg/s] ≈ Q·c` while on, and `_BalanceResidual[mg/s]` should stay closed throughout.

---

### 1.3 Source/Sink Terms — `rtswss.input`

**File path**: `<input_folder>/rtswss.input`
**Source code**: `src/RTInitSW.h` → `readRTSSFileSW()`, `src/SourceSink.h`

Defines spatial source/sink zones where solute mass is injected into (or extracted from) the surface water domain. Each zone has a polygon defining its spatial extent.

#### Global Source/Sink Type (`rtsstype`)

| `rtsstype` Value | Description |
|:---:|-------------|
| 3 | Area mass flux injection (surface spreading, e.g., fertilizer) [mg/m²/s]. The code automatically distributes the flux over grid cell areas. |

#### Per-Species Source/Sink Type (`spec_X_sstype`)

| Value | Description |
|:---:|-------------|
| 0 | **Constant**. The source/sink rate is a fixed value given by `spec_X_ssval`. |
| 1 | **Linear interpolation timeseries**. The rate varies over time according to a timeseries file (`spec_X_ssfile`). Linear interpolation is applied between data points. |
| 2 | **Step-function timeseries**. Same as 1 but uses step-wise (piecewise constant) interpolation between data points. Used for instantaneous events like fertilizer application. |

#### File Format

```text
sscount : <number_of_sources>

id : <source_name>              // Source identifier
rtsstype : <global_type>        // Global source/sink type (3 = area flux)
polygon : <polygon_file>        // Polygon defining the spatial zone

// Per-species (X = species index)
spec_X_sstype : <type>          // 0=constant, 1=linear timeseries
spec_X_ssval  : <value>         // Area mass flux [mg/m^2/s] (if sstype=0)
spec_X_ssfile : <filename>      // Timeseries file (if sstype=1)
```

#### Example — Fertilizer Application

```text
sscount : 1

id : fertilizer
rtsstype : 3                            // Area mass flux [mg/m^2/s]
polygon : polygonfer.input              // Fertilizer application area

spec_0_sstype : 1                       // NH4: linear timeseries
spec_0_ssfile : rtfertNH4-DON-NH4.input

spec_1_sstype : 0                       // NO3: constant (no NO3 fertilizer)
spec_1_ssval  : 0

spec_2_sstype : 0                       // DO: no source
spec_2_ssval  : 0

spec_3_sstype : 0                       // DOC: no source
spec_3_ssval  : 0

spec_4_sstype : 1                       // DON: linear timeseries
spec_4_ssfile : rtfertDON-DON-NH4.input
```

> **Note**: If `rtswss.input` is not found, the model proceeds without surface water solute sources/sinks (a warning is printed).

---

### 1.4 Reaction Parameters — `reactionsw.input`

**File path**: `<input_folder>/reactionsw.input`
**Source code**: `src/RTInitSW.h` → `readRTReactionFile()`, `src/RTStateSW.h`

This file is **only required** when `ReactionModule : 1` is set in `rttransportsw.input`. It configures the multi-component Monod kinetic reaction model for surface water, including nitrogen cycle processes, reaeration, and environmental correction factors.

#### Reaction Mass Basis and Fixed Species Order

The active reaction kernel is a five-species system with fixed zero-based indices:

| Index | Name | Concentration basis |
|---:|---|---|
| 0 | NH4 | Compound mass, `mg NH4/L` (molecular mass about 18.04 g/mol) |
| 1 | NO3 | Compound mass, `mg NO3/L` (molecular mass about 62.00 g/mol) |
| 2 | DO | Oxygen mass, `mg O2/L` |
| 3 | DOC | Elemental carbon mass, `mg C/L` |
| 4 | DON | Elemental nitrogen mass, `mg N/L` |

This is not a universal `mg N/L` model. The input reader does not convert `mg NH4-N/L` to `mg NH4/L`, or `mg NO3-N/L` to `mg NO3/L`; use the bases above for initial conditions, boundaries, source files, and reaction parameters. `SPECIES_ID` and `SPECIES_Name` are informational and do not remap the concentration array. The fixed reaction kernel requires `n_mass >= 5`; when fewer than five species are configured, the reaction step returns without solving the nitrogen-cycle equations.

The base rates are `R_min [mg N/L/s]` (DON mineralization), `R_nit [mg NH4/L/s]` (NH4 nitrification), and `R_denit`/`R_aer [mg C/L/s]` (DOC consumption). The source terms apply these hard-coded mass ratios: `18.04/14.01 = 1.2877` for NH4 per mineralized N, `62.00/18.04 = 3.4368` for NO3 per oxidized NH4, `64.00/18.04 = 3.5477` for O2 per nitrified NH4, `32.00/14.01 = 2.2841` for O2 per mineralized N, `32.00/12.01 = 2.6644` for O2 per respired C, and `(0.8*62.00)/12.01 ≈ 4.130` for NO3 per denitrified C. The `0.8` factor is a hard-coded denitrification efficiency factor, not an input parameter.

**Validation status.** The parser checks that required keys exist using a `-999` sentinel, but generally does not enforce physical minimum/maximum values. Ranges below are recommended calibration ranges, not program-enforced limits.

#### Global Configuration

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `ReactModel_Type` | int | **Yes** | Reaction model type: `0` = no reaction, `1` = multi-component Monod kinetics. |
| `Nitrogen_Cycle_Simulation` | int | **Yes** | Nitrogen cycle switch: `0` = off, `1` = on. Enables NH4→NO3→N2 pathway. |
| `SPECIES_ID` | int list | No | Informational. Species index list (semicolon-separated). Example: `0; 1; 2; 3; 4`. |
| `SPECIES_Name` | string list | No | Informational. Species name list (semicolon-separated). Example: `NH4; NO3; DO; DOC; DON`. |

#### First-Order Decay Rates

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `Lambda_1` | real list | [1/s] | First-order decay rate in **liquid phase** for each species (semicolon-separated, one value per species). Typical range: 0 to 1×10⁻⁵. Set to `0.0` for non-decaying species. |
| `Lambda_2` | real list | [1/s] | First-order decay rate in **solid phase** for each species. Same format as `Lambda_1`. Set to `0.0` for surface water (no solid phase). |

#### Nitrification (NH4⁺ → NO3⁻)

Biological oxidation of ammonium to nitrate by nitrifying bacteria (Nitrosomonas).

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Nit` | real | [mg NH4/L/s] | 1×10⁻⁷ – 1×10⁻⁴ | Maximum nitrification rate on the NH4 compound-mass basis. |
| `K_Monod_NH4` | real | [mg NH4/L] | 0.01 – 5.0 | NH4 half-saturation concentration on the compound-mass basis. |
| `K_Monod_DO` | real | [mg O2/L] | 0.01 – 5.0 | DO half-saturation concentration. |

#### Denitrification (NO3⁻ → N2)

Reduction of nitrate to nitrogen gas by heterotrophic bacteria under anaerobic conditions. This is the dominant nitrogen loss pathway in flooded rice paddies.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Denit` | real | [mg C/L/s] | 1×10⁻⁷ – 3×10⁻⁴ | Maximum denitrification rate on the DOC-C basis; nitrate consumption uses the hard-coded 0.8 factor. |
| `K_Monod_NO3` | real | [mg NO3/L] | 0.01 – 5.0 | NO3 half-saturation concentration on the compound-mass basis. |
| `K_Monod_DOC` | real | [mg C/L] | 0.01 – 20.0 | DOC half-saturation concentration. |
| `Ki_Inhib_DO` | real | [mg O2/L] | 0.01 – 5.0 | DO inhibition concentration; larger DO reduces denitrification. |

#### Aerobic Respiration (DOC + O2 → CO2)

Heterotrophic consumption of dissolved organic carbon with oxygen as electron acceptor.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Hetero` | real | [mg C/L/s] | 1×10⁻⁷ – 3×10⁻⁴ | Maximum aerobic respiration rate on the DOC-C basis. |
| `K_Monod_DOC_Aerobic` | real | [mg C/L] | 0.01 – 20.0 | DOC half-saturation concentration on the elemental-carbon basis. |
| `K_Monod_DO_Aerobic` | real | [mg O2/L] | 0.1 – 2.0 | Half-saturation constant for DO. |

#### Mineralization (DON → NH4⁺)

Decomposition of dissolved organic nitrogen to ammonium.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Min` | real | [mg N/L/s] | 1×10⁻⁸ – 1×10⁻⁴ | Maximum mineralization rate on the DON-N basis. |
| `K_Monod_DON` | real | [mg N/L] | 0.01 – 10.0 | DON half-saturation concentration on the elemental-N basis. |
| `K_Monod_DO_Min` | real | [mg O2/L] | 0.01 – 5.0 | DO half-saturation concentration for mineralization. |

#### Surface-Water-Specific Source/Sink Parameters

These parameters govern atmospheric reaeration and sediment-water solute exchange processes unique to surface water.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Ka_sw` | real | [1/s] | 0 – 1×10⁻⁴ | Atmospheric reaeration rate constant. Set to `0.0` to disable reaeration for validation cases. Typical: 5×10⁻⁵. |
| `DO_sat` | real | [mg O2/L] | 0 – 15 | Saturated dissolved oxygen concentration used by reaeration. |
| `DOC_eq_sw` | real | [mg C/L] | 2 – 20 | Equilibrium DOC concentration in surface water (sediment release endpoint). Fertile paddy fields: ~10 mg/L. |
| `DON_eq_sw` | real | [mg N/L] | 0.1 – 5.0 | Equilibrium DON concentration in surface water (sediment release endpoint). |
| `K_rel_sw` | real | [1/s] | 0 – 1×10⁻⁵ | Release rate constant for sediment-water solute exchange. Set to `0.0` to disable. |

#### Environmental Correction Factors

Temperature correction for surface-water biological rates. `Opt_pH` is currently retained for input compatibility only.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Temp_Coeff_Theta` | real | [-] | 1.00 – 1.10 | Base of the current surface-water multiplier `theta^(20-Opt_Temp)`. |
| `Opt_Temp` | real | [°C] | 0 – 40 | Reference temperature in the surface-water multiplier; current `T_env` is fixed at 20°C. |
| `Opt_pH` | real | [-] | 0 – 14 | Parsed and required, but not used in the current reaction equations. |

#### Example

```text
# === Global ===
ReactModel_Type           : 1                  // Multispecies Monod kinetics
Nitrogen_Cycle_Simulation : 1                  // Enable nitrogen cycle
SPECIES_ID   : 0; 1; 2; 3; 4
SPECIES_Name : NH4; NO3; DO; DOC; DON

# === Decay ===
Lambda_1 : 0.0; 0.0; 0.0; 0.0; 0.0
Lambda_2 : 0.0; 0.0; 0.0; 0.0; 0.0

# === Nitrification ===
Rate_Max_Nit  : 1.0e-5        // [mg NH4/L/s]
K_Monod_NH4   : 0.5           // [mg NH4/L]
K_Monod_DO    : 1.5           // [mg O2/L]

# === Denitrification ===
Rate_Max_Denit : 5.0e-6       // [mg C/L/s]
K_Monod_NO3    : 0.5          // [mg NO3/L]
K_Monod_DOC    : 5.0          // [mg C/L]
Ki_Inhib_DO    : 1.0          // [mg O2/L]

# === Aerobic Respiration ===
Rate_Max_Hetero       : 5.0e-5  // [mg C/L/s]
K_Monod_DOC_Aerobic   : 5.0     // [mg C/L]
K_Monod_DO_Aerobic    : 1.0     // [mg O2/L]

# === Mineralization ===
Rate_Max_Min   : 5.0e-6        // [mg N/L/s]
K_Monod_DON    : 0.21015       // [mg N/L]
K_Monod_DO_Min : 0.48          // [mg O2/L]

# === Surface Water Source/Sink ===
Ka_sw      : 0.0               // [1/s] (disabled for validation)
DO_sat     : 8.1               // [mg/L]
DOC_eq_sw  : 10.0              // [mg C/L]
DON_eq_sw  : 0.5               // [mg N/L]
K_rel_sw   : 0.0               // [1/s]

# === Environmental Correction ===
Temp_Coeff_Theta : 1.047
Opt_Temp         : 25.0        // [°C]
Opt_pH           : 7.0
```

---

## 2. Groundwater Solute Transport

### 2.1 Transport Parameters — `rttransportgw.input`

**File path**: `<input_folder>/rttransportgw.input`
**Source code**: `src/RTInitGW.h` → `readRTFile()`

Same format as `rttransportsw.input` with additional keys for groundwater-specific settings. Key differences:

All aqueous concentrations consumed by the transport and reaction code use **mg/L**. The parser does not convert values supplied in kg/m3 or other concentration units. Molecular diffusion and dispersivities are entered in SI units as shown below.

#### Additional Global Parameters

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `rt_scheme` | int | **Yes** | Numerical solver for the transport matrix: `1` = Gauss-Seidel, `2` = Kokkos BiCGSTAB (recommended for large 3D domains). |
| `Up_Weighting_vplus` | real | **Yes** | Upstream weighting coefficient for positive velocity direction [-]. Range: [0.5, 1.0]. `1.0` = full upstream, `0.5` = central difference. |
| `Up_Weighting_vminus` | real | **Yes** | Upstream weighting coefficient for negative velocity direction [-]. Range: [0.0, 0.5]. `0.0` = full upstream, `0.5` = central difference. |

#### Additional Per-Species Parameters

Groundwater transport includes a **solid phase** (adsorbed concentration) in addition to the aqueous phase.

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `solid_mode` | int | If `ReactionModule=1` | Solid-phase initial condition mode: `0` = constant value (`solid_val`), `1` = spatial file (`solid_file`). It is not required when reactions are disabled. |
| `solid_val` | real | If `solid_mode=0` | Uniform initial solid-phase (adsorbed) concentration [mg/kg]. |
| `solid_file` | string | If `solid_mode=1` | Filename for spatially distributed solid-phase initial concentration. |
| `aq_phase` | int | No | Not parsed by the current groundwater reader. Use `aq_mode`; do not rely on `aq_phase` as an alias. |

#### Example

```text
n_mass : 5

// === Species 0: NH4 ===
id : 0
aq_mode : 1
aq_file : rtinitialgwNH4.input
solid_mode : 1
solid_file : rtinitialgwsolidNH4.input
diffusion_molecular : 1.76e-9
alpha_L : 0.5
alpha_T : 0.05

// === Species 2: DO ===
id : 2
aq_mode : 0
aq_val : 1.0                            // [mg/L]
solid_mode : 0
solid_val : 0.0                         // [mg/kg] (DO not adsorbed)
diffusion_molecular : 1.76e-9
alpha_L : 0.5
alpha_T : 0.05

// ... other species ...

ReactionModule : 1
rt_scheme : 2                            // BiCGSTAB solver
Up_Weighting_vplus : 1.0                 // Full upstream (positive)
Up_Weighting_vminus : 0.0               // Full upstream (negative)
```

---

### 2.2 Boundary Conditions — `rtgwbc.input`

**File path**: `<input_folder>/rtgwbc.input`
**Source code**: `src/RTBCGW.h` (lines 47–55), `src/RTInitGW.h`

Defines solute boundary conditions for the 3D groundwater domain. Boundaries are specified by polygon files and face directions.

#### Face Direction Codes

| Direction | Face | Description |
|:---:|------|-------------|
| 1 | +x | East face |
| 2 | −x | West face |
| 3 | +y | North face |
| 4 | −y | South face |
| 5 | +z | **Top face** (land surface, upward) |
| 6 | −z | **Bottom face** (aquifer base, downward) |

#### Boundary Condition Types

| Value | Constant Name | Description |
|:---:|------|-------------|
| 1 | `SUB_RT_BC_NOFLOW` | **No-flow**. Neither water nor solute crosses the boundary. |
| 2 | `SUB_RT_BC_Dirichlet_CONST` | **Constant concentration** (Dirichlet). Fixed concentration specified by `spec_X_bcval`. |
| 3 | `SUB_RT_BC_Neumann_CONST` | **Constant flux** (Neumann). Fixed solute flux [mg/m²/s]. |
| 4 | `SUB_RT_BC_Cauchy_CONST` | **Mixed (Cauchy)**. Combined concentration/flux condition. |
| 5 | `SUB_RT_BC_Dirichlet_T` | **Time-varying Dirichlet**. Concentration from timeseries file. |
| 6 | `SUB_RT_BC_Neumann_T` | **Time-varying Neumann**. Flux from timeseries file. |
| 7 | `SUB_RT_BC_Cauchy_T` | **Time-varying Cauchy**. Mixed condition from timeseries file. |
| 8 | `SUB_RT_BC_SWE` | **Surface-water exchange**. Couples groundwater solute transport with overland flow concentrations. Used at the land surface to exchange solute between surface and subsurface domains. |
| 9 | `SUB_RT_BC_FD` | **Free drainage** (zero-gradient outflow). Solute leaves the domain at the same concentration as the interior cell. Used at the aquifer base. |

#### Example — Top (SWE coupling) + Bottom (free drainage)

```text
bccount : 2

// --- Top boundary: coupled with surface water ---
id : rttop
direction : 6                            // ZMINUS: upper/surface face in the groundwater convention
polygon : polygontop.input
spec_0_bctype : 8                        // NH4: surface-water exchange
spec_1_bctype : 8                        // NO3: surface-water exchange
spec_2_bctype : 8                        // DO: surface-water exchange
spec_3_bctype : 8                        // DOC: surface-water exchange
spec_4_bctype : 8                        // DON: surface-water exchange

// --- Bottom boundary: free drainage ---
id : rtbottom
direction : 5                            // ZPLUS: lower/bottom face in the groundwater convention
polygon : polygontop.input
spec_0_bctype : 9                        // NH4: free drainage
spec_1_bctype : 9
spec_2_bctype : 9
spec_3_bctype : 9
spec_4_bctype : 9
```

---

### 2.3 Reaction Parameters — `reactiongw.input`

**File path**: `<input_folder>/reactiongw.input`
**Source code**: `src/RTInitGW.h` → `readRTReactionFile()`, `src/RTStateGW.h`

Contains all surface water reaction parameters plus **sorption**, **bulk density**, and **groundwater-specific** reaeration/release parameters. Only required when `ReactionModule : 1`.

The reaction implementation uses aqueous concentrations in mg/L, solid-phase concentrations in mg/kg, distribution coefficients in L/kg, and bulk density internally in kg/L. `Bulk_Density` is the only density value in this file: enter it as g/cm3. Numerically, 1 g/cm3 equals 1 kg/L, so no numerical conversion is required after parsing. Values above 10 are rejected as likely kg/m3 input.

#### Reaction Mass Basis and Solver Preconditions

The groundwater reaction kernel uses the fixed five-species order and mixed mass basis: index 0 `NH4` is `mg NH4/L`, index 1 `NO3` is `mg NO3/L`, index 2 `DO` is `mg O2/L`, index 3 `DOC` is `mg C/L`, and index 4 `DON` is `mg N/L`. It is not a universal nitrogen-mass model, and the input reader does not convert `mg NH4-N/L` or `mg NO3-N/L` into compound mass. `SPECIES_ID` and `SPECIES_Name` do not reorder the arrays. The fixed-index reaction solver requires `n_mass >= 5`; otherwise it returns without solving the nitrogen-cycle reactions.

The groundwater base rates are `R_min [mg N/L/s]`, `R_nit [mg NH4/L/s]`, and `R_denit`/`R_aer [mg C/L/s]`. The source terms use the hard-coded mass ratios `18.04/14.01`, `62.00/18.04`, `64.00/18.04`, `32.00/14.01`, `32.00/12.01`, and `(0.8*62.00)/12.01`. Groundwater additionally divides each reaction source term by its species retardation factor `Rf`. Reaeration is `Ka_DO * max(0, theta_s-theta) * max(0, DO_sat-DO)`, so it is reduced when the pore space is water-filled.

**Validation status.** Required reaction keys are checked for presence using `-999`; most physical ranges are not parser-enforced. The ranges below are recommended calibration ranges. `Bulk_Density` is the exception: it must be supplied in `g/cm³`, and values above `10` are rejected as likely `kg/m³` input. For an active five-species run, provide five values for every per-species list and use `Sorption_Type` values `0` through `4` only.

#### Additional Global Parameters

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Bulk_Density` | real | [g/cm³] | 1.1 – 1.4 | Bulk density of the porous medium (ρb). Input must be in g/cm³; the model stores and uses it internally as kg/L. Typical rice plow layer: 1.30 g/cm³. |

#### Sorption Parameters (per-species, semicolon-separated)

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `Sorption_Type` | int list | [-] | Sorption model for each species: `0` = none, `1` = linear (`Kd`), `2` = Freundlich (`Kf`, `Nf`), `3` = Langmuir (`Alpha_D`, `Beta_D`), `4` = kinetic/non-equilibrium (`Kd`, `Beta`). |
| `Kd` | real list | [L/kg] | Distribution coefficient for **Linear** sorption. With `C` in mg/L and `S` in mg/kg, use L/kg. Typical NH4 in paddy: 2.0–5.0. Set `0.0` for non-adsorbing species. |
| `Kf` | real list | [(mg/kg)/(mg/L)^Nf] | Freundlich constant. Set `0.0` if not using Freundlich. |
| `Nf` | real list | [-] | Freundlich nonlinearity exponent. Typical: 0.5–1.2. Set `1.0` if not used. |
| `Alpha_D` | real list | [L/mg] | Langmuir affinity constant in `S = Alpha_D*Beta_D*C/(1+Alpha_D*C)`. Set `0.0` if not using Langmuir. |
| `Beta_D` | real list | [mg/kg] | Langmuir maximum sorption capacity. Set `0.0` if not using Langmuir. |
| `Beta` | real list | [1/s] | Kinetic sorption rate constant used by `Sorption_Type=4`. It is not a Langmuir parameter. |

#### Sorption Model Summary

| `Sorption_Type` | Model | Required Parameters | Isotherm Equation |
|:---:|-------|---------------------|-------------------|
| 0 | None (no sorption) | — | — |
| 1 | Linear | `Kd` | S = Kd × C |
| 2 | Freundlich | `Kf`, `Nf` | S = Kf × C^Nf |
| 3 | Langmuir | `Alpha_D`, `Beta_D` | S = (Alpha_D × Beta_D × C) / (1 + Alpha_D × C) |
| 4 | Kinetic/non-equilibrium | `Kd`, `Beta` | Implemented in the transport assembly; `Beta` is a first-order mass-transfer rate [1/s] |

where S = adsorbed concentration [mg/kg], C = dissolved concentration [mg/L].

#### Root Solute Uptake — `TSCF`

`TSCF` is a required per-species list in `reactiongw.input`. The name means **Transpiration Stream Concentration Factor**. It controls the fraction of the dissolved groundwater concentration that is carried out of the groundwater domain with root transpiration. It is applied independently to each species and uses the fixed species order `NH4; NO3; DO; DOC; DON`.

| Key | Type | Units | Recommended physical range | Parser status | Implementation meaning |
|---|---|---|---:|---|---|
| `TSCF` | real list | [-] | `0.0–1.0` | Only presence is checked; numerical bounds are not enforced | Species-specific fraction multiplying the root-transpiration solute flux. |

For a cell and species, the active source/sink calculation is conceptually:

```text
solute_outflow_root = |transpGW| * cell_area * C_dissolved * 1000 * TSCF   [mg/s]
```

where `transpGW < 0` denotes upward root transpiration, `C_dissolved` is the current dissolved concentration in `[mg/L]`, the factor `1000` converts `L` to `m³`, and `TSCF` is dimensionless. The same solute mass is removed from the local reactive-transport equation. Root uptake is therefore a groundwater source/sink process, not an additional Monod reaction term.

Recommended interpretation:

- `TSCF = 0.0`: transpiration removes water but no dissolved mass of that species. This is the appropriate neutral value for DO, DOC, and DON unless plant uptake or transport data justify otherwise.
- `0.0 < TSCF < 1.0`: only a fraction of the dissolved concentration is transported with the transpiration stream. Use this for partial solute exclusion or species-dependent uptake.
- `TSCF = 1.0`: the transpiring water carries solute at the local dissolved concentration; this is the upper end of the usual passive-advection interpretation.
- `TSCF > 1.0`: mathematically accepted by the current parser and produces concentration enrichment in the transpiration stream, but it should only be used with explicit experimental evidence and documented as an enrichment factor.
- `TSCF < 0.0`: the parser does not reject it, but the active calculation treats it as disabled because uptake is applied only when `TSCF > 0`; do not use negative values.

For the five-species nitrogen-cycle example, a conservative starting list is:

```text
TSCF : 0.0; 0.0; 0.0; 0.0; 0.0   // NH4, NO3, DO, DOC, DON
```

If plant nitrogen uptake is explicitly represented, a sensitivity case may assign nonzero values to NH4 and NO3, for example `TSCF : 0.3; 0.3; 0.0; 0.0; 0.0`. This is a calibration assumption, not a universal rice value. Do not infer a nonzero TSCF merely because WOFOST computes root water uptake: water uptake and solute uptake are separate, and soil evaporation is implemented as pure-water loss (`TSCF` has no effect on evaporation).

The current code does not enforce a physical `0–1` range, so users must validate the list themselves. Provide exactly one value per configured species; for the active five-species reaction kernel, provide five values. The relevant implementation is in `src/RTInitGW.h`, `src/RTStateGW.h`, `src/SourceSink.h`, and `src/RTIntegrator.h`.

#### Groundwater-Specific Source/Sink Parameters

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Ka_DO` | real | [1/s] | 0 – 1×10⁻⁵ | Reaeration rate constant for groundwater. Typical: 1×10⁻⁵ (≈ 0.5/3600). |
| `DOC_eq` | real | [mg C/L] | 2 – 20 | DOC release equilibrium concentration from sediment. |
| `DON_eq` | real | [mg N/L] | 0.1 – 5.0 | DON release equilibrium concentration from sediment. |
| `K_rel_DOC` | real | [1/s] | 0 – 1×10⁻⁵ | DOC release rate constant. Set `0.0` to disable. |
| `K_rel_DON` | real | [1/s] | 0 – 1×10⁻⁵ | DON release rate constant. Set `0.0` to disable. |

#### Groundwater Reaction Parameters

The following parameters are read as scalar values and are required when `ReactModel_Type : 1`. The recommended ranges are calibration guidance only; the current parser checks key presence rather than these numerical limits.

| Key | Type | Units | Recommended range | Implementation meaning |
|---|---|---|---:|---|
| `Rate_Max_Nit` | real | [mg NH4/L/s] | 1×10⁻⁷ to 1×10⁻⁴ | Maximum nitrification base rate; directly removes NH4 compound mass. |
| `K_Monod_NH4` | real | [mg NH4/L] | 0.01 to 5.0 | NH4 half-saturation concentration. |
| `K_Monod_DO` | real | [mg O2/L] | 0.01 to 5.0 | DO half-saturation concentration for nitrification. |
| `Rate_Max_Denit` | real | [mg C/L/s] | 1×10⁻⁷ to 3×10⁻⁴ | Maximum denitrification base rate; directly removes DOC-C. |
| `K_Monod_NO3` | real | [mg NO3/L] | 0.01 to 5.0 | NO3 half-saturation concentration on the compound-mass basis. |
| `K_Monod_DOC` | real | [mg C/L] | 0.01 to 20.0 | DOC half-saturation concentration for denitrification. |
| `Ki_Inhib_DO` | real | [mg O2/L] | 0.01 to 5.0 | DO inhibition concentration in `Ki/(Ki+DO)`. |
| `Rate_Max_Hetero` | real | [mg C/L/s] | 1×10⁻⁷ to 3×10⁻⁴ | Maximum aerobic DOC respiration rate. |
| `K_Monod_DOC_Aerobic` | real | [mg C/L] | 0.01 to 20.0 | DOC half-saturation concentration for aerobic respiration. |
| `K_Monod_DO_Aerobic` | real | [mg O2/L] | 0.01 to 5.0 | DO half-saturation concentration for aerobic respiration. |
| `Rate_Max_Min` | real | [mg N/L/s] | 1×10⁻⁸ to 1×10⁻⁴ | Maximum mineralization base rate; directly removes DON-N. |
| `K_Monod_DON` | real | [mg N/L] | 0.01 to 10.0 | DON half-saturation concentration on the elemental-N basis. |
| `K_Monod_DO_Min` | real | [mg O2/L] | 0.01 to 5.0 | DO half-saturation concentration for mineralization. |
| `Ka_DO` | real | [1/s] | 0 to 1×10⁻⁴ | Groundwater reaeration coefficient. The source is additionally multiplied by air-filled porosity. |
| `DO_sat` | real | [mg O2/L] | 0 to 15 | DO saturation concentration used by reaeration. |
| `DOC_eq` | real | [mg C/L] | 0 to 50 | Equilibrium DOC concentration for release toward the liquid phase. |
| `DON_eq` | real | [mg N/L] | 0 to 10 | Equilibrium DON concentration for release toward the liquid phase. |
| `K_rel_DOC` | real | [1/s] | 0 to 1×10⁻⁴ | DOC release coefficient; zero disables this release source. |
| `K_rel_DON` | real | [1/s] | 0 to 1×10⁻⁴ | DON release coefficient; zero disables this release source. |
| `Temp_Coeff_Theta` | real | [-] | 1.00 to 1.10 | Parsed environmental temperature coefficient, currently not effective in the groundwater kernel because `Temp_Corr` is fixed at 1.0. |
| `Opt_Temp` | real | [°C] | 0 to 40 | Parsed reference temperature, currently not effective in the groundwater kernel. |
| `Opt_pH` | real | [-] | 0 to 14 | Parsed and required, but not used in the current groundwater reaction equations. |

`Lambda_1` and `Lambda_2` are per-species first-order decay coefficients in `[1/s]`, normally recommended from `0` to `1×10⁻⁵`; they are applied in the transport assembly rather than in the five-species reaction source equations. `Lambda_1` acts on dissolved concentration and `Lambda_2` on sorbed concentration. A zero value disables that decay pathway.

#### Example

```text
# === Global ===
ReactModel_Type           : 1
Bulk_Density              : 1.30                // input [g/cm^3], internal kg/L
Nitrogen_Cycle_Simulation : 1
SPECIES_ID   : 0; 1; 2; 3; 4
SPECIES_Name : NH4; NO3; DO; DOC; DON

# === Sorption ===
Sorption_Type : 1; 0; 0; 0; 1          // NH4: linear, DON: linear
Kd            : 3.0; 0.0; 0.0; 0.0; 0.5  // [L/kg]
Kf            : 0.0; 0.0; 0.0; 0.0; 0.0
Nf            : 1.0; 1.0; 1.0; 1.0; 1.0
Alpha_D       : 0.0; 0.0; 0.0; 0.0; 0.0
Beta_D        : 0.0; 0.0; 0.0; 0.0; 0.0
Beta          : 0.0; 0.0; 0.0; 0.0; 0.0

# === Decay ===
Lambda_1 : 0.0; 0.0; 0.0; 0.0; 0.0
Lambda_2 : 0.0; 0.0; 0.0; 0.0; 0.0

# === Nitrification ===
Rate_Max_Nit  : 5.0e-6                 // [mg NH4/L/s]
K_Monod_NH4   : 0.5                    // [mg NH4/L]
K_Monod_DO    : 0.5                    // [mg O2/L]

# === Denitrification ===
Rate_Max_Denit : 5.0e-5                // [mg C/L/s]
K_Monod_NO3    : 0.5                   // [mg NO3/L]
K_Monod_DOC    : 2.0                   // [mg C/L]
Ki_Inhib_DO    : 0.2                   // [mg O2/L]

# === Aerobic Respiration ===
Rate_Max_Hetero       : 1.0e-4         // [mg C/L/s]
K_Monod_DOC_Aerobic   : 5.0            // [mg C/L]
K_Monod_DO_Aerobic    : 0.5            // [mg O2/L]

# === Mineralization ===
Rate_Max_Min  : 3.0e-6                // [mg N/L/s]
K_Monod_DON   : 0.21015               // [mg N/L]
K_Monod_DO_Min: 0.48                  // [mg O2/L]

# === Reaeration / Release ===
Ka_DO      : 1.0e-5                    // [1/s]
DO_sat     : 8.1                       // [mg/L]
DOC_eq     : 20.0                      // [mg C/L]
DON_eq     : 2.0                       // [mg N/L]
K_rel_DOC  : 0.0                       // [1/s]
K_rel_DON  : 0.0                       // [1/s]

# === Environmental Correction ===
Temp_Coeff_Theta : 1.047
Opt_Temp         : 25.0
Opt_pH           : 7.0
```

---

### 2.4 Groundwater Source/Sink — `rtgwss.input`

**File path**: `<input_folder>/rtgwss.input`
**Source code**: `src/RTInitGW.h`

Defines solute source/sink zones in the groundwater domain. The format mirrors `rtgwbc.input` with polygon-based spatial selection.

#### File Format

```text
sscount : <number_of_sources>

id : <source_name>
rtsstype : <type>                // Global source/sink type
polygon : <polygon_file>         // Spatial extent polygon

// Per-species values or time series
spec_X_sstype : <type>
spec_X_ssval  : <value>
spec_X_ssfile : <timeseries_file>
```

#### Source/Sink Types (`sstype`)

| Value | Description |
|:---:|-------------|
| `rtsstype=1` | Volumetric mass injection [mg/L/s] in the selected groundwater cells. |
| `rtsstype=3` | Area mass flux injection [mg/m²/s], distributed over the selected vertical column thickness. |

`rtsstype=0` belongs to the separate WOFOST water source/sink structure, not the reactive-solute `RTGwSS` calculation. `rtsstype=2` is not handled by the current reactive-solute injection branch and should not be used here.

#### Example — Root Zone for WOFOST Coupling

```text
sscount : 1
id : root1
rtsstype : 1                     // Reactive-solute volumetric injection
polygon : polygonroot.input      // Root zone polygon
spec_0_sstype : 0
spec_0_ssval : 0.0                // [mg/L/s]
```

Root water uptake and soil evaporation are supplied through the groundwater water-balance state and WOFOST coupling; they are not configured by `rtsstype=0` in this reactive-solute file.

---

### 2.5 Initial Condition Files

**Source code**: `src/RTInitSW.h` → `readRtICFileSW()`, `src/RTInitGW.h` → `readRtICFileGW()`

Initial concentration fields can be specified as spatial grids. Surface water accepts the flexible ESRI-style header parser; groundwater uses the three-record labelled header described below.

#### File Naming Convention

| Domain | Phase | Pattern | Example |
|--------|-------|---------|---------|
| Surface water | Liquid (aqueous) | `rtinitialsw<Name>.input` | `rtinitialswNH4.input` |
| Groundwater | Liquid (aqueous) | `rtinitialgw<Name>.input` | `rtinitialgwNO3.input` |
| Groundwater | Solid (adsorbed) | `rtinitialgwsolid<Name>.input` | `rtinitialgwsolidNH4.input` |

#### Groundwater format: three header records plus a 3-D field

```text
ncols <ncols>
nrows <nrows>
nodata_value <value>
<data values...>     // ncols × nrows × nz values, layer-major/row-major
```

The groundwater reader consumes the first three header records as a label followed by a number. Use the labelled form above. `ncols` and `nrows` must equal the global horizontal domain dimensions; the data section must contain `nx_glob * ny_glob * nz_glob` values. There is no `nz` header field.

#### Surface-water format: simplified or ESRI-style 2-D grid

```text
NCOLS <ncols>
NROWS <nrows>
XLLCORNER <x>         // Ignored by the model, logged only
YLLCORNER <y>         // Ignored by the model, logged only
CELLSIZE <size>       // Ignored by the model, logged only
NODATA_value <value>
<data values...>      // Row-major, nrows × ncols values
```

#### Requirements

- Grid dimensions must **exactly match** the global domain: `NCOLS = nx_glob`, `NROWS = ny_glob`.
- Data values are in row-major order.
- Nodata values are read as numeric concentration values; the current groundwater reader does not mask them.
- Units: liquid phase [mg/L], solid phase [mg/kg].

---

## 3. WOFOST Crop Growth Model

### 3.1 Crop Parameters — `cropparameter.input`

**File path**: `<input_folder>/cropparameter.input`
**Source code**: `src/cropsrc/CropInit.h` → `readCropParameter()`, `src/cropsrc/CropState.h`

This file contains all WOFOST crop growth parameters plus soil physical properties for the root zone. The file uses colon-separated key–value pairs. Tabulated parameters are specified as space-separated x,y pairs.

#### Meta Data

| Key | Type | Description |
|-----|------|-------------|
| `CRPNAM` | string | Crop name (informational). |
`CROP_NO` is not parsed by the current `CropInit` reader. Use `CROPNAM` only for the crop name metadata; all numerical crop parameters use the keys listed below.

#### Emergence Parameters

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `TBASEM` | real | [°C] | Lower threshold temperature for emergence. Below this, no emergence occurs. |
| `TEFFMX` | real | [°C] | Maximum effective temperature for emergence. Above this, no additional thermal time accumulates. |
| `TSUMEM` | real | [°C·d] | Temperature sum from sowing to emergence. |

#### Phenology Parameters

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `IDSL` | int | [-] | Development trigger mode: `0` = temperature only, `1` = daylength only, `2` = temperature + daylength. |
| `DLO` | real | [hr] | Optimum daylength for development. Set to `-99` if `IDSL = 0`. |
| `DLC` | real | [hr] | Critical daylength (lower threshold). Set to `-99` if `IDSL = 0`. |
| `TSUM1` | real | [°C·d] | Temperature sum from emergence to anthesis. |
| `TSUM2` | real | [°C·d] | Temperature sum from anthesis to maturity. |
| `DVSI` | real | [-] | Initial development stage (DVS) at simulation start. `0.0` = emergence, `1.0` = anthesis, `2.0` = maturity. |
| `DVSEND` | real | [-] | Development stage at physiological maturity/harvest (typically `2.0`). |

#### Initial Conditions

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `TDWI` | real | [kg/ha] | 10 – 200 | Initial total crop dry weight at emergence. |
| `LAIEM` | real | [ha/ha] | 0.0001 – 0.01 | Leaf area index at emergence. |
| `RGRLAI` | real | [ha/ha/d] | 0.001 – 0.05 | Maximum relative increase in LAI during exponential growth phase. |

#### Green Area Parameters

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `SPAN` | real | [d] | Life span of leaves growing at 35°C. Controls leaf senescence rate. |
| `TBASE` | real | [°C] | Lower threshold temperature for ageing of leaves. |
| `SPA` | real | [ha/kg] | Specific pod area. Set to `0.0` for crops without pods. |

#### Tabulated Parameters

All tables are specified as space-separated x,y pairs. The x-variable is typically DVS (development stage) or temperature. The model uses linear interpolation between points.

| Table | x-variable | y-variable | Description |
|-------|-----------|-----------|-------------|
| `DTSMTB` | Mean temp [°C] | Daily increase in temp. sum [°C/d] | Defines effective temperature function for phenology. |
| `SLATB` | DVS [-] | Specific leaf area [ha/kg] | Leaf area per unit leaf biomass. Controls LAI growth. |
| `SSATB` | DVS [-] | Specific stem area [ha/kg] | Stem area per unit stem biomass. Often set to 0. |
| `KDIFTB` | DVS [-] | Extinction coefficient for diffuse visible light [-] | Controls light interception by canopy. |
| `EFFTB` | Daily mean temp [°C] | Light-use efficiency of single leaf [kg CO₂ ha⁻¹ hr⁻¹ / J m⁻² s⁻¹] | |
| `AMAXTB` | DVS [-] | Maximum leaf CO₂ assimilation rate [kg CO₂ ha⁻¹ hr⁻¹] | |
| `TMPFTB` | Average temp [°C] | Reduction factor for AMAX [-] | Temperature effect on photosynthesis. Range: 0–1. |
| `TMNFTB` | Min. temp [°C] | Reduction factor for gross assimilation [-] | Frost/damage effect. Range: 0–1. |
| `RFSETB` | DVS [-] | Reduction factor for senescence [-] | |

#### Conversion Efficiency Parameters

Biomass conversion efficiency from assimilates (CH₂O) to plant organs.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `CVL` | real | [kg/kg] | 0.60 – 0.80 | Conversion efficiency to **leaves**. |
| `CVO` | real | [kg/kg] | 0.60 – 0.80 | Conversion efficiency to **storage organs** (grain, tuber). |
| `CVR` | real | [kg/kg] | 0.60 – 0.80 | Conversion efficiency to **roots**. |
| `CVS` | real | [kg/kg] | 0.60 – 0.80 | Conversion efficiency to **stems**. |

#### Maintenance Respiration Parameters

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `Q10` | real | [-] | Relative increase in respiration rate per 10°C temperature increase. Typical: 2.0. |
| `RML` | real | [kg CH₂O/kg/d] | Relative maintenance respiration rate of **leaves**. |
| `RMO` | real | [kg CH₂O/kg/d] | Relative maintenance respiration rate of **storage organs**. |
| `RMR` | real | [kg CH₂O/kg/d] | Relative maintenance respiration rate of **roots**. |
| `RMS` | real | [kg CH₂O/kg/d] | Relative maintenance respiration rate of **stems**. |

#### Partitioning Tables

Define how total above-ground dry matter is distributed among organs as a function of development stage.

| Table | Description | Constraint |
|-------|-------------|------------|
| `FRTB` | Fraction of **total** dry matter partitioned to roots vs. DVS. | Independent; the remainder goes to above-ground. |
| `FLTB` | Fraction of **above-ground** dry matter partitioned to leaves vs. DVS. | FL + FS + FO = 1.0 at all DVS values. |
| `FSTB` | Fraction of **above-ground** dry matter partitioned to stems vs. DVS. | |
| `FOTB` | Fraction of **above-ground** dry matter partitioned to storage organs vs. DVS. | |

#### Death Rate Parameters

| Key/Table | Type | Units | Description |
|-----------|------|-------|-------------|
| `PERDL` | real | [d⁻¹] | Maximum relative death rate of leaves due to water stress. |
| `RDRRTB` | table | — | Relative death rate of roots as function of DVS [d⁻¹]. |
| `RDRSTB` | table | — | Relative death rate of stems as function of DVS [d⁻¹]. |

#### Water Use Parameters

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `CFET` | real | [-] | Correction factor for transpiration rate. Typically `1.0`. |
| `DEPNR` | real | [-] | Crop group number for soil water depletion (1–5). Higher values = greater drought tolerance. |
| `IAIRDU` | int | [-] | Air ducts in roots: `1` = present (wetland crops like rice), `0` = absent. |
| `IOX` | int | [-] | Oxygen stress switch: `0` = no oxygen stress, `1` = oxygen stress active. |

#### Rooting Parameters

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `RDI` | real | [cm] | 5 – 20 | Initial rooting depth at emergence. |
| `RRI` | real | [cm/d] | 0.5 – 2.0 | Maximum daily increase in rooting depth. |
| `RDMCR` | real | [cm] | 30 – 150 | Maximum rooting depth (crop-specific). Rice: ~80 cm; wheat: ~120 cm. |

#### Nutrient Parameters (NPK)

##### Concentration Limits

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `NMINSO` | real | [kg N/kg] | Minimum (residual) N concentration in storage organs. |
| `NMINVE` | real | [kg N/kg] | Minimum (residual) N concentration in vegetative organs. |
| `NMAXSO` | real | [kg N/kg] | Maximum N concentration in storage organs. |
| `NMAXVE` | real | [kg N/kg] | Maximum N concentration in leaves at DVS=0. |
| `PMINSO` | real | [kg P/kg] | Minimum P concentration in storage organs. |
| `PMINVE` | real | [kg P/kg] | Minimum P concentration in vegetative organs. |
| `PMAXSO` | real | [kg P/kg] | Maximum P concentration in storage organs. |
| `PMAXVE` | real | [kg P/kg] | Maximum P concentration in leaves at DVS=0. |
| `KMINSO` | real | [kg K/kg] | Minimum K concentration in storage organs. |
| `KMINVE` | real | [kg K/kg] | Minimum K concentration in vegetative organs. |
| `KMAXSO` | real | [kg K/kg] | Maximum K concentration in storage organs. |
| `KMAXVE` | real | [kg K/kg] | Maximum K concentration in leaves at DVS=0. |

##### Uptake and Translocation

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `RNUPTAKEMAX` | real | [kg N/ha/d] | Maximum daily nitrogen uptake rate. |
| `RPUPTAKEMAX` | real | [kg P/ha/d] | Maximum daily phosphorus uptake rate. |
| `RKUPTAKEMAX` | real | [kg K/ha/d] | Maximum daily potassium uptake rate. |
| `DVS_NPK_STOP` | real | [-] | DVS above which crop NPK uptake ceases. |
| `DVS_NPK_TRANSL` | real | [-] | DVS above which NPK translocation from vegetative to storage organs occurs. |

##### NPK Stress and Distribution

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `NLAI_NPK` | real | [-] | Coefficient for LAI reduction due to NPK stress. |
| `NSLA_NPK` | real | [-] | Coefficient for SLA reduction due to NPK stress. |
| `NPART` | real | [-] | Coefficient for leaf allocation effect due to N stress. |
| `NLUE_NPK` | real | [-] | Coefficient for CO₂ assimilation reduction due to NPK stress. |
| `NPK_TRANSLRT_FR` | real | [-] | Fraction of NPK translocated from roots. |
| `RDRLV_NPK` | real | [d⁻¹] | Maximum relative death rate of leaves due to NPK stress. |

##### Organ-Specific Concentration Fractions

| Key | Description |
|-----|-------------|
| `NMAXRT_FR` | Max N conc. in roots as fraction of max N in leaves. |
| `NMAXST_FR` | Max N conc. in stems as fraction of max N in leaves. |
| `NCRIT_FR` | Critical N conc. as fraction of max N conc. |
| `NRESIDRT` | Residual N fraction in roots [kg N/kg]. |
| `PMAXRT_FR`, `PMAXST_FR`, `PCRIT_FR`, `PRESIDRT` | Same for phosphorus. |
| `KMAXRT_FR`, `KMAXST_FR`, `KCRIT_FR`, `KRESIDRT` | Same for potassium. |

##### Translocation Time Coefficients

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `TCNT` | int | [d] | Time coefficient for N translocation. |
| `TCPT` | int | [d] | Time coefficient for P translocation. |
| `TCKT` | int | [d] | Time coefficient for K translocation. |

##### NPK Concentration Tables

| Table | x-variable | y-variable | Description |
|-------|-----------|-----------|-------------|
| `NMAXLV_TB` | DVS [-] | Max N conc. in leaves [kg N/kg] | |
| `PMAXLV_TB` | DVS [-] | Max P conc. in leaves [kg P/kg] | |
| `KMAXLV_TB` | DVS [-] | Max K conc. in leaves [kg K/kg] | |

##### Other

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `YZERO` | real | [kg/ha] | Maximum amount of vegetative organs at zero yield. |
| `NFIX` | real | [kg/kg] | Fraction of N uptake from biological fixation (e.g., legumes). |

#### Soil Physical Properties (root zone)

These parameters define the soil water retention properties used by WOFOST for water stress calculations.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `SMW` | real | [cm³/cm³] | 0.05 – 0.15 | Soil moisture at **wilting point** (ψ ≈ −1500 cm or −150 kPa). |
| `SMFCF` | real | [cm³/cm³] | 0.20 – 0.40 | Soil moisture at **field capacity** (ψ ≈ −100 to −340 cm). Southern clay: 0.35–0.40. |
| `SM0` | real | [cm³/cm³] | 0.35 – 0.50 | Soil moisture at **saturation** (≈ porosity Φ). |
| `CRAIRC` | real | [cm³/cm³] | 0.04 – 0.10 | **Critical air content** for aeration. Below this, oxygen stress limits root function. |

#### Example

```text
// --- Meta ---
CRPNAM : Gaoyou Rice (European Rice ecotype)

// --- Emergence ---
TBASEM : 8.0          // [°C]
TEFFMX : 22.0         // [°C]
TSUMEM : 100.0        // [°C·d]

// --- Phenology ---
IDSL   : 0            // Temperature-driven only
TSUM1  : 1490.0       // [°C·d]
TSUM2  : 610.0        // [°C·d]
DVSI   : 0.33         // Initial DVS
DVSEND : 2.00

DTSMTB : 0.00 0.00 10.00 0.00 25.00 15.00

// --- Initial ---
TDWI   : 130.0        // [kg/ha]
LAIEM  : 0.0007       // [ha/ha]
RGRLAI : 0.0070       // [ha/ha/d]

// --- Green Area ---
SPAN   : 50.0         // [d]
TBASE  : 0.0          // [°C]
SPA    : 0.000

SLATB  : 0.00 0.0017 0.10 0.0017 0.58 0.0019 0.86 0.00200 1.10 0.0015 1.40 0.0015 2.00 0.0010
SSATB  : 0.00 0.0000 2.00 0.0000

// --- Assimilation ---
KDIFTB : 0.00 0.600 2.00 0.600
EFFTB  : 0.0 0.450 40.0 0.450
AMAXTB : 0.00 40.00 1.00 40.00 1.30 40.00 2.00 40.00
TMPFTB : 0.0 0.010 10.0 0.000 25.0 1.000 35.0 1.000 42.0 0.000
TMNFTB : 0.0 0.00 3.0 1.00

// --- Conversion ---
CVL : 0.754
CVO : 0.684
CVR : 0.754
CVS : 0.754

// --- Respiration ---
Q10  : 2.0
RML  : 0.0300
RMO  : 0.0150
RMR  : 0.0150
RMS  : 0.0200
RFSETB : 0.00 0.300 2.00 0.400

// --- Partitioning ---
FRTB : 0.000 0.200 0.500 0.200 0.800 0.150 1.000 0.150 1.100 0.000 2.000 0.000
FLTB : 0.000 0.400 0.580 0.520 0.800 0.380 0.850 0.000 0.900 0.000 1.000 0.000 2.000 0.000
FSTB : 0.000 0.600 0.580 0.480 0.800 0.620 0.850 0.400 0.900 0.000 1.000 0.000 2.000 0.000
FOTB : 0.000 0.000 0.580 0.000 0.800 0.000 0.850 0.600 0.900 1.000 1.000 1.000 2.000 1.000

// --- Death ---
PERDL  : 0.030
RDRRTB : 0.00 0.000 1.50 0.000 1.51 0.020 2.00 0.020
RDRSTB : 0.00 0.000 1.50 0.000 1.51 0.020 2.00 0.020

// --- Water ---
CFET   : 1.00
DEPNR  : 3.5
IAIRDU : 1
IOX    : 0

// --- Rooting ---
RDI    : 10.0         // [cm]
RRI    : 1.2          // [cm/d]
RDMCR  : 40.0         // [cm]

// --- Nutrients ---
NMINSO : 0.0043
NMINVE : 0.0021
NMAXSO : 0.0144
NMAXVE : 0.0600
// ... (all nutrient parameters) ...
RNUPTAKEMAX : 7.2    // [kg N/ha/d]
RPUPTAKEMAX : 0.8    // [kg P/ha/d]
RKUPTAKEMAX : 7.4    // [kg K/ha/d]

// --- Soil Physics ---
SMW    : 0.085        // [cm^3/cm^3]
SMFCF  : 0.27
SM0    : 0.46
CRAIRC : 0.19
```

---

### 3.2 Meteorological Data — `cropmeteo.input`

**File path**: `<input_folder>/cropmeteo.input`
**Source code**: `src/cropsrc/MeteoState.h`, `src/cropsrc/MeteoInit.h`

Provides daily meteorological forcing data for the WOFOST crop model. The file is a whitespace-separated numeric table with at least nine columns; it is not YAML and it does not use the generic `time value` format from Section 4.

The file must also contain location metadata in comment or plain-text lines, including `LAT`, `LON`, and `ELEV`. Optional `ANGSTA` and `ANGSTB` metadata override the default Angstrom coefficients.

#### Format

Time series with one row per time step. Each row contains:

| Column | Type | Units | Description |
|--------|------|-------|-------------|
| 1 | real | [d or model time unit] | Time/index stored as `MeteoState::time`; the reader does not convert it. |
| 2 | int/real | [-] | Calendar year, e.g. `2004`. |
| 3 | int/real | [day of year] | DOY, normally 1–366. |
| 4 | real | [J/m²/d] | Daily global solar radiation (`IRRAD`). |
| 5 | real | [°C] | Daily minimum temperature (`TMIN`). |
| 6 | real | [°C] | Daily maximum temperature (`TMAX`). |
| 7 | real | [hPa] | Vapour pressure (`VAP`), not kPa. |
| 8 | real | [m/s] | Wind speed at 2 m (`WIND`). |
| 9 | real | [mm/d] | Daily precipitation (`RAIN`). |

The reader requires at least nine numeric columns. It stores computed open-water, bare-soil, and reference-crop ET as `E0`, `ES0`, and `ET0` in [cm/d] after converting the internal calculation result from mm/d.

---

### 3.3 Agricultural Management — `agro.input`

**File path**: `<input_folder>/agro.input`
**Format**: YAML

Defines crop variety selection and planting/harvest schedule.

#### Key Fields

| Key | Type | Description |
|-----|------|-------------|
| `crop` | string | Crop type name (e.g., `sugar-beet`, `rice`). |
| `variety` | string | Variety identifier matching a section in the crop YAML file. |
| `crop_start_type` | string | Start type: `sowing` or `emergence`. |
| `crop_start_date` | string | Crop start date (YYYY-MM-DD). |
| `crop_end_type` | string | End type: `harvest` or `maturity`. |
| `crop_end_date` | string | Crop end date (YYYY-MM-DD). |
| `max_duration` | int | Maximum crop growth duration [days]. |

#### Example

```yaml
crop: sugar-beet
variety: sugar-beet-601
crop_start_type: emergence
crop_start_date: 2004-01-02
crop_end_type: harvest
crop_end_date: 2006-07-20
max_duration: 300
```

---

### 3.4 WOFOST–SERGHEI Coupling Interface

**Source code**: `src/SourceSinkCrop.h`

The WOFOST crop growth model is two-way coupled with the SERGHEI groundwater model through the root zone.

#### Data Flow

```
┌─────────────────────────────┐
│       SERGHEI (GW)          │
│                             │
│  3D soil moisture field     │──── Root zone soil moisture ────┐
│                             │     (depth-weighted average     │
│                             │      over [0, RD])              │
└─────────────────────────────┘                                 │
                                                                ▼
┌──────────────────────────────────────────────────────────┐
│                    WOFOST Crop Model                      │
│                                                          │
│  Inputs:  soil moisture, meteorological data             │
│  Outputs: LAI, root depth (RD), biomass, DVS             │
│                                                          │
│  Internal processes:                                     │
│    - Phenology (development stage)                       │
│    - CO₂ assimilation & growth respiration               │
│    - Biomass partitioning (leaf, stem, root, storage)    │
│    - Transpiration & water stress                        │
│    - Nutrient (NPK) uptake & stress (Not currently available)                     │
└──────────────────────────────────────────────────────────┘
        │                                    │
        │ LAI (leaf area index)              │ RD (root depth)
        │                                    │
        ▼                                    ▼
┌──────────────────────────────────────────────────────────┐
│       SERGHEI (GW + SW)                                  │
│                                                          │
│  ET partitioning:                                        │
│    T_split = 1 - exp(-k_ext × LAI)    (transpiration)   │
│    E_split = 1 - T_split              (soil evaporation) │
│                                                          │
│  Root water extraction:                                  │
│    Depth range: [0, RD]                                  │
│    Rate: based on soil moisture and potential ET         │
└──────────────────────────────────────────────────────────┘
```

#### Coupling Variables

| Variable | Direction | Units | Description | Source Code |
|----------|-----------|-------|-------------|-------------|
| Root zone soil moisture (RZSM) | GW → WOFOST | [m³/m³] | Depth-weighted average volumetric water content over [0, RD]. | `compute_root_zone_moisture()` in `SourceSinkCrop.h` |
| LAI | WOFOST → SERGHEI | [m²/m²] | Leaf area index. Range: 0–7. Controls radiation interception. | `wofost_lai` in `SourceSinkCrop.h` |
| RD (root depth) | WOFOST → SERGHEI | [cm] | Current rooting depth. Range: RDI to RDMCR. Determines water extraction depth. | `wofost_rd` in `SourceSinkCrop.h` |
| T_split | Computed | [-] | Transpiration fraction of potential ET. | `compute_evap_partitioning()` in `SourceSinkCrop.h` |
| E_split | Computed | [-] | Soil evaporation fraction of potential ET. | `E_split = 1 - T_split` |

#### Water Stress Calculation

WOFOST uses the root zone soil moisture to compute a water stress factor:

```
RZSM = (Σ θ_i × Δz_i) / RD    for layers i within [0, RD]

Stress factor depends on:
  - SMW (wilting point)
  - SMFCF (field capacity)
  - SM0 (saturation)
  - DEPNR (crop group number)
  - IAIRDU (air ducts, for waterlogging tolerance)
```

- **Drought stress**: When RZSM approaches SMW, transpiration is reduced.
- **Oxygen stress**: When RZSM approaches SM0 (waterlogging), root function is impaired (if `IOX = 1`). Rice with `IAIRDU = 1` is tolerant.

---

## 4. Triggered Irrigation and Drainage Boundaries

### 4.1 Overview

**Source code**: `src/Parser.h` → `readExtBCFile()` (parsing and validation), `src/BC.h` (trigger evaluation and hydraulic application, `SWE_BC_TRIGGERED_IRRIGATION = 13`, `SWE_BC_TRIGGERED_DRAINAGE = 14`), `src/TimeIntegrator.h` (calls `ExtBC::apply()` each time step).

**Reference test cases**: `tests/t15-gaoyou-wofost-triggerIrrigation/` (depth-trigger irrigation), `tests/t16-gaoyou-wofost-triggerIrrigation/` (W1 shallow-wet-dry with a soil-moisture-trigger stage)

`bctype : 13` and `bctype : 14` implement rule-based, automatically controlled surface-water boundaries for paddy field management:

- **Triggered irrigation** (`bctype : 13`, `SWE_BC_TRIGGERED_IRRIGATION`): an inflow boundary. A management strategy table (per crop stage) defines target water depths; when the mean observed water depth falls below the stage lower limit, water is applied at the stage discharge until the upper limit is reached. Stages may alternatively (or additionally) trigger on **root-zone soil moisture** via the per-stage `trig_type`/`theta_trig` columns (see [4.3](#43-irrigation-strategy-file)), e.g. for shallow-wet-dry (浅湿露干) schedules.
- **Triggered drainage** (`bctype : 14`, `SWE_BC_TRIGGERED_DRAINAGE`): an outflow boundary. When the mean observed water depth exceeds the stage "drain-on" threshold, water is removed — either by critical (free) outflow or by a fixed pumping discharge — until the "drain-off" threshold is reached.

Both boundary types are evaluated **every hydrodynamic time step**: the stage is selected from the strategy table, the trigger condition is tested against the mean water depth over a user-defined observation region, and the boundary discharge is applied or closed accordingly. On/off hysteresis is built in (see [4.3](#43-irrigation-strategy-file) and [4.4](#44-drainage-strategy-file)).

Two stage-selection modes are supported:

- **DAT mode** (`use_DVS : 0`): the stage is selected by *Days After Transplanting/sowing*, `DAT = floor((etime - crop_start_time)/86400) + 1`, where `cropStartTime` is fixed at 0 in the current code, so `DAT` is effectively days since simulation start plus one.
- **DVS mode** (`use_DVS : 1`): available only in `CROP_GROWTH_MODEL` builds. The stage is selected by the **spatial mean WOFOST development stage (DVS)** averaged over the trigger observation cells. The mean DVS is recomputed every hydrodynamic step from the WOFOST state (initialised to `DVSI` from `cropparameter.input` at startup, advanced once per daily crop step). In non-crop builds `use_DVS : 1` cannot be served: the boundary prints a one-time warning ("DVS data are not connected …") at the first time step and irrigation/drainage stays disabled (not a fatal error). Note that stage selection starts at `DVSI`: with `DVSI : 0.33` the row covering 0.10–0.35 is matched from t = 0, so a regreening row (0.00–0.10) can never be reached — set `DVSI` accordingly.

### 4.2 Boundary Definition — `extbc.input`

**File path**: `<input_folder>/extbc.input`
**Source code**: `src/Parser.h` → `readExtBCFile()`

#### Boundary Condition Types

| Value | Constant Name | Class | Description |
|:---:|------|:---:|-------------|
| 13 | `SWE_BC_TRIGGERED_IRRIGATION` | Inflow | Rule-based irrigation inlet; discharge applied when depth trigger fires. |
| 14 | `SWE_BC_TRIGGERED_DRAINAGE` | Outflow | Rule-based drainage outlet; critical or fixed-Q outflow when depth trigger fires. |

#### Per-Boundary Keys

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `id` | string | **Yes** | Boundary identifier (must precede the other keys of its block). |
| `bctype` | int | **Yes** | `13` for triggered irrigation, `14` for triggered drainage. |
| `polygon` | string | **Yes** | Polygon file selecting the boundary cells where water enters/leaves the domain (same format as other hydraulic BCs; see [4.5](#45-trigger-polygons)). |
| `trigger_polygon` | string | No* | Observation region over which the mean trigger depth (and, in DVS mode, the mean crop DVS) is computed. Either a polygon file name or the keyword `ALL_DOMAIN`. If omitted, the boundary's own `polygon` is used and a warning is printed. |
| `direction` | real ×2 | **Yes** | Unit-normal vector of the boundary (normalized internally), e.g. `0.0 1.0`. Flow is applied along `+direction` for irrigation and along the outward normal for drainage. |
| `bcvals` | real ×1 | **Yes** | Exactly one value: the **global trigger mode** for irrigation (see [4.3](#43-irrigation-strategy-file)); per-stage `trig_type` columns in the strategy file override it row by row. For drainage boundaries it is parsed but unused; provide `1` for consistency. |
| `irrigation_strategy` | string | If bctype=13 | Filename of the irrigation strategy table (see [4.3](#43-irrigation-strategy-file)). Mandatory and fatal if missing for `bctype : 13`. |
| `drainage_strategy` | string | If bctype=14 | Filename of the drainage strategy table (see [4.4](#44-drainage-strategy-file)). Mandatory and fatal if missing for `bctype : 14`. |

#### Parsing and Validation Rules

The reader enforces the following (all violations are fatal errors unless noted):

- Triggered irrigation/drainage requires **exactly 1** value in `bcvals` (the trigger mode).
- `trigger_polygon` may be `ALL_DOMAIN` (case-insensitive) or a polygon file with **1 point** (nearest cell), **≥ 3 points** (polygon region), or 0 points meaning all cells; anything else is rejected.
- Strategy files must define `use_DVS` (0 or 1) and a positive `ndata`, followed by exactly `ndata` data rows; unknown keys inside a strategy file are rejected.
- Strategy rows must satisfy `max ≥ min` for DAT and DVS (and water-depth) ranges, and rows must be **sorted and non-overlapping** in both DAT and DVS.
- Irrigation rows with `irrigation_enabled : 1` must have `Qirr ≠ 0`; drainage rows with `drainage_enabled : 1` and `drain_mode : 2` must have `Qdrain ≠ 0`.
- Drainage rows require `drain_on_h ≥ drain_off_h`, and `drain_mode ∈ {1, 2}`.

#### Example — `tests/t15-gaoyou-wofost-triggerIrrigation/extbc.input`

```text
bccount : 2

id : trigger_drainage
bctype : 14
polygon : polygonoutlet.input
trigger_polygon : ALL_DOMAIN
direction : 0.0 1.0
bcvals : 1
drainage_strategy : test2.input

id : trigger_irrigation
bctype : 13
polygon : polygoninlet.input
trigger_polygon : ALL_DOMAIN
direction : 0.0 1.0
bcvals : 1
irrigation_strategy : test.input

//bcvals(0)：全局触发模式，`1` 为只用地表水深，`2` 为只用根系层土壤水（策略行需提供 theta_trig），`3` 为任一满足触发，`4` 为同时满足触发；策略行的 trig_type 列可逐阶段覆盖
```

> **Note**: Trigger modes 2–4 are the **root-zone soil-moisture triggers** described in [4.3](#43-irrigation-strategy-file). They require a build with both `CROP_GROWTH_MODEL` and `SERGHEI_RE_MODEL` enabled, and every stage row whose effective mode is 2–4 must provide `theta_trig` (a startup error is issued otherwise). In builds without those macros, modes 2–4 and any `theta_*` columns are rejected at startup.

### 4.3 Irrigation Strategy File

**Referenced from**: `extbc.input` key `irrigation_strategy` (file must be in the same directory as `extbc.input`)
**Source code**: `src/Parser.h` → `readExtBCFile()` (reading); `src/BC.h` → `selectIrrigationStage()`, `checkIrrigationTrigger()`, `applyTriggeredIrrigationStage()`

#### File Format

```text
use_DVS : <0|1>
ndata : <number_of_stage_rows>
// dat_min dat_max dvs_min dvs_max h_min h_max gw_depth_trigger Qirr irrigation_enabled drainage_enabled [trig_type] [theta_trig] [theta_stop]
<dat_min> <dat_max> <dvs_min> <dvs_max> <h_min> <h_max> <gw_depth_trigger> <Qirr> <irrigation_enabled> <drainage_enabled> [<trig_type>] [<theta_trig>] [<theta_stop>]
... (ndata rows)
```

`//` comment lines and blank lines are skipped. The comment header line documenting the column order is conventional but not required.

#### Keys and Columns

| Key/Column | Type | Units | Description |
|-----|------|-------|-------------|
| `use_DVS` | int | - | Stage-selection mode: `0` = by DAT, `1` = by mean crop DVS (crop builds only). |
| `ndata` | int | - | Number of strategy rows that follow. |
| `dat_min`, `dat_max` | real | [d] | DAT window of the stage (inclusive). Used when `use_DVS : 0`. |
| `dvs_min`, `dvs_max` | real | [-] | Mean-DVS window of the stage (inclusive). Used when `use_DVS : 1`. |
| `h_min` | real | [m] | Depth trigger: irrigation **starts** when the mean observed depth < `h_min` (while inactive). |
| `h_max` | real | [m] | Depth stop: irrigation **stops** when the mean observed depth ≥ `h_max` (while active). If `h_max < h_min` it is raised to `h_min` internally. |
| `gw_depth_trigger` | real | [m] | Groundwater-depth trigger threshold (reserved). Required by the format but **not yet connected**; enter a large dummy value such as `999.0`. |
| `Qirr` | real | [m³/s] | Irrigation discharge applied while the trigger is active. Must be non-zero if `irrigation_enabled : 1`. |
| `irrigation_enabled` | int (0/1) | - | Whether irrigation may run in this stage (`0` disables irrigation in the stage). |
| `drainage_enabled` | int (0/1) | - | Informational column recording whether the paired drainage strategy allows drainage in this stage. It is **not** used by the irrigation boundary itself; the drainage boundary is governed by its own strategy file. |
| `trig_type` | int (0–4) | No (default 0) | **Per-stage trigger type** appended after `drainage_enabled` (old 10-column files keep working): `0` = use the global `bcvals(0)` mode; `1` = surface-water depth only; `2` = root-zone soil moisture only; `3` = either condition opens; `4` = both conditions required. |
| `theta_trig` | real | If mode ≥ 2 | Soil-moisture trigger threshold as a **fraction of saturation**, ∈(0,1] (e.g. `0.70` = 70% of θ_sat). Irrigation opens when the trigger-cell mean relative saturation S = Σ(wc·dz)/Σ(θ_s·dz) falls below it. Mandatory whenever the effective mode is 2–4 (otherwise a startup error). |
| `theta_stop` | real | No (default −1) | Optional early-stop threshold (fraction of saturation, ∈[theta_trig,1]): irrigation stops once S ≥ `theta_stop` (e.g. rainfall has re-wetted the profile). `-1` disables early stop; `h_max` remains the hard upper bound in all cases. |

#### Trigger Mode (`bcvals(0)` and per-stage `trig_type`)

The effective mode of a stage is the row's `trig_type` if given, otherwise the global `bcvals(0)`:

| Value | Opening condition (while inactive) |
|:---:|---------|
| 1 | mean depth < `h_min` (surface-water depth trigger). |
| 2 | root-zone relative saturation S < `theta_trig` (soil-moisture trigger; the row must provide `theta_trig`). |
| 3 | either condition (depth **or** soil moisture). |
| 4 | both conditions (depth **and** soil moisture). |
| other | Rejected at runtime with a warning; falls back to the depth trigger. |

While active, **all modes share one continuation rule**: irrigation continues while `mean depth < h_max` **and** S < `theta_stop` (the soil-moisture clause only applies when `theta_stop > 0`). `h_max` is always a hard upper bound — a stage never keeps filling past it even if the soil is still dry.

S is recomputed **every hydrodynamic step** from the current groundwater moisture field `wc` over the trigger cells: per column it is depth-weighted over the WOFOST root zone [0, RD] and normalized by the layer saturations (θ_s from the van Genuchten soil table, `vg.input`). In the asynchronous coupling S lags the surface-water step by at most one groundwater sub-cycle. Modes 2–4 require a build with both `CROP_GROWTH_MODEL` and `SERGHEI_RE_MODEL` enabled; other builds reject such rows at startup.

#### Applied Hydraulics

While the trigger is active, the boundary first flattens the water surface across the boundary cross-section (`flattenWaterSurface`, as for `SWE_BC_Q_T`), then distributes the total discharge `Qirr` over the boundary cells (`distributeDischarge`):

- Wet boundary cells receive unit-width discharge weighted by their water depth, `hu = Q·(h/Σh)/ds · nx`, `hv = ... · ny`.
- If the whole boundary is dry, each cell is first raised to the critical depth `hc = (q²/g)^(1/3)` (at least `hmin`) and the unit discharge `q = Q/(N·ds)` is imposed, so momentum is not zeroed on startup.
- When the trigger is inactive (or the stage is disabled/absent), boundary `hu`, `hv` are set to zero.

#### Example — Irrigation Strategy (`test.input`)

```text
use_DVS : 0
ndata : 7
// dat_min dat_max dvs_min dvs_max h_min h_max gw_depth_trigger Qirr irrigation_enabled drainage_enabled
1 7 0.00 0.10 0.030 0.070 999.0 0.003 1 0
8 25 0.10 0.35 0.020 0.060 999.0 0.003 1 0
26 40 0.35 0.55 0.000 0.000 999.0 0.000 0 1
41 65 0.55 0.85 0.030 0.070 999.0 0.003 1 0
66 85 0.85 1.05 0.020 0.060 999.0 0.003 1 0
86 110 1.05 1.65 0.000 0.030 999.0 0.002 1 0
111 125 1.65 2.00 0.000 0.000 999.0 0.000 0 1
```

This example keeps a 2–7 cm ponded layer during vegetative growth (stages 1–5), reduces the target depth to 0–3 cm before heading (stage 6, Qirr lowered to 0.002 m³/s), and disables irrigation during the final drainage stage (days 111–125, `irrigation_enabled : 0`).

#### Example — W1 shallow-wet-dry (浅湿露干) full-season schedule (`tests/t16-gaoyou-wofost-triggerIrrigation/test.input`)

```text
use_DVS : 0
ndata : 7
// dat_min dat_max dvs_min dvs_max h_min h_max gw_depth_trigger Qirr irrigation_enabled drainage_enabled trig_type theta_trig theta_stop
1 7 0.00 0.10 0.005 0.020 999.0 0.003 1 0 1 -1 -1
8 25 0.10 0.35 0.000 0.010 999.0 0.003 1 0 2 0.70 -1
26 40 0.35 0.55 0.000 0.000 999.0 0.000 0 1
41 65 0.55 0.85 0.010 0.020 999.0 0.003 1 0 1 -1 -1
66 85 0.85 1.05 0.005 0.015 999.0 0.003 1 0 1 -1 -1
86 110 1.05 1.65 0.005 0.010 999.0 0.002 1 0 1 -1 -1
111 125 1.65 2.00 0.000 0.000 999.0 0.000 0 1
```

Stage by stage: regreening (days 1–7) keeps a 5–20 mm depth layer with the depth trigger; **early tillering (days 8–25) switches to the soil-moisture trigger** (`trig_type : 2`, `theta_trig : 0.70`) — irrigation opens when the root-zone relative saturation S drops below 70% of θ_sat and fills up to `h_max` = 10 mm (`h_min` is unused for opening in mode 2, but `h_max` remains the hard stop); days 26–40 are the mid-season drainage (晒田) stage with irrigation disabled; jointing–booting, heading–flowering and milk ripening return to shallow depth triggers (10→20, 5→15, 5→10 mm); days 111–125 are the final dry-off (落干) with irrigation disabled. The drainage half of the schedule (drain-to-dry in stages 3 and 7) lives in the paired drainage strategy file. The three optional columns may be omitted entirely on any row (the disabled rows above keep the old 10-column form), and `theta_stop` may replace `-1` with e.g. `0.95` on the moisture-trigger row to stop early once rainfall has re-wetted the profile. Note that the strategy parser does **not** strip inline `//` comments, so nothing may follow the last column of a row.

### 4.4 Drainage Strategy File

**Referenced from**: `extbc.input` key `drainage_strategy` (file must be in the same directory as `extbc.input`)
**Source code**: `src/Parser.h` → `readExtBCFile()` (reading); `src/BC.h` → `selectDrainageStage()`, `checkDrainageTrigger()`, `applyTriggeredDrainageStage()`

#### File Format

```text
use_DVS : <0|1>
ndata : <number_of_stage_rows>
// dat_min dat_max dvs_min dvs_max drain_on_h drain_off_h drain_mode Qdrain drainage_enabled
<dat_min> <dat_max> <dvs_min> <dvs_max> <drain_on_h> <drain_off_h> <drain_mode> <Qdrain> <drainage_enabled>
... (ndata rows)
```

#### Keys and Columns

| Key/Column | Type | Units | Description |
|-----|------|-------|-------------|
| `use_DVS` | int | - | Stage-selection mode: `0` = by DAT, `1` = by mean crop DVS (crop builds only). |
| `ndata` | int | - | Number of strategy rows that follow. |
| `dat_min`, `dat_max` | real | [d] | DAT window of the stage (inclusive). |
| `dvs_min`, `dvs_max` | real | [-] | Mean-DVS window of the stage (inclusive). |
| `drain_on_h` | real | [m] | Drainage **starts** when the mean observed depth > `drain_on_h` (while inactive). Must be ≥ `drain_off_h`. |
| `drain_off_h` | real | [m] | Drainage **stops** when the mean observed depth ≤ `drain_off_h` (while active). |
| `drain_mode` | int | - | Outflow mode: `1` = critical/free outflow (each boundary cell with `h ≥ hmin` discharges at Froude 1, `|u| = √(gh)`); `2` = fixed pumping discharge `Qdrain` distributed over wet boundary cells proportionally to their depth. |
| `Qdrain` | real | [m³/s] | Pumping discharge for `drain_mode : 2`. Unused for `drain_mode : 1` (enter any value, e.g. `0.005`). Must be non-zero if `drainage_enabled : 1` **and** `drain_mode : 2`. |
| `drainage_enabled` | int (0/1) | - | Whether drainage may run in this stage. |

When the trigger is inactive (or the stage is disabled/absent), boundary `hu`, `hv` are set to zero.

#### Example — Drainage Strategy (`test2.input`)

```text
use_DVS : 0
ndata : 7
// dat_min dat_max dvs_min dvs_max drain_on_h drain_off_h drain_mode Qdrain drainage_enabled
1 7 0.00 0.10 0.090 0.040 1 0.005 1
8 25 0.10 0.35 0.090 0.040 1 0.005 1
26 40 0.35 0.55 0.005 0.000 2 0.005 1
41 65 0.55 0.85 0.080 0.060 1 0.005 1
66 85 0.85 1.05 0.080 0.040 1 0.005 1
86 110 1.05 1.65 0.040 0.010 2 0.005 1
111 125 1.65 2.00 0.005 0.000 2 0.005 1
```

This example drains whenever the ponded depth exceeds 8–9 cm during the early stages (critical outflow), forces the field nearly dry during mid-season drainage and the final period (`drain_mode : 2`, fixed pumping `Qdrain = 0.005` m³/s down to ~0 depth), and relaxes the thresholds before heading.

### 4.5 Trigger Polygons

**Referenced from**: `extbc.input` key `trigger_polygon`
**Source code**: `src/Parser.h` (file reading); `src/BC.h` → `find_trigger_cells()`

The observation region over which the mean water depth (and mean DVS) is computed is selected by `trigger_polygon`:

| Option | Format | Meaning |
|--------|--------|---------|
| `ALL_DOMAIN` | keyword (case-insensitive) | All non-nodata cells of the domain are observation cells. |
| Polygon file, ≥ 3 points | `NPOINTS <n>` then n `x y` vertex lines | Cells whose centres lie inside the polygon. |
| Polygon file, 1 point | `NPOINTS 1` then one `x y` line | The single grid cell whose centre is nearest to the point. |

The polygon file format is the same `NPOINTS`-header format used by other boundary/source polygons, e.g.:

```text
NPOINTS 4
27.875 -0.125
27.875 0.625
28.625 0.625
28.625 -0.125
```

If no trigger cell is found at all, a warning is printed and the trigger stays permanently off. If `trigger_polygon` is omitted, the boundary's own `polygon` is reused as the observation region (warning printed once).

### 4.6 Runtime Behaviour and Output

- **Evaluation frequency.** Triggered boundaries are re-evaluated at every hydrodynamic time step (`ExtBC::apply` in `TimeIntegrator.h`), before halo exchange; discharges are integrated afterwards like any other boundary.
- **Discharge bookkeeping.** `bctype : 13` is classified as an *inflow* boundary and `bctype : 14` as an *outflow* boundary, so their fluxes appear in `domainTimeSeries.out` under `BoundaryInflow`/`BoundaryOutflow` and in the per-boundary `BoundaryFlow_<i>` / `BoundaryAccum_<i>` columns (index `i` = boundary order in `extbc.input`, 0-based; inflow is reported positive, outflow negative).
- **Solute transport interaction.** If surface reactive transport is enabled but `rtswbc.input` provides no per-species BC for a triggered boundary, gate-aware fallbacks apply (one-time warning each):
  - Triggered irrigation → `RTSW_BC_TRIGGERED_IRRIGATION` with concentration `0.0` while the irrigation is **active** (clear-water Dirichlet) and zero-gradient while **inactive** (boundary concentration follows the interior; no solute is forced when no water enters). The advective inflow term of the mass budget captures the delivered solute while active.
  - Triggered drainage → `RTSW_BC_TRIGGERED_DRAINAGE` (zero-gradient outflow at interior concentration; solute leaves only while drainage actually discharges).
  Specify `spec_X_bctype` entries in `rtswbc.input` (see [Section 1.2](#12-boundary-conditions--rtswbcinput)) to override the fallback. For fertilizer-loaded irrigation water use the gated type `spec_X_bctype : 5` with `spec_X_bcval` or a `spec_X_bcfile` timeseries, so the concentration BC follows the irrigation on/off cycle. Plain types 2/3/4 are still accepted (they apply unconditionally, independent of the trigger state).
- **Coupling with WOFOST.** In `CROP_GROWTH_MODEL` builds the crop-aware overload is used, enabling `use_DVS : 1` stage selection from the mean crop DVS over the trigger cells and providing the root depth (RD) for the soil-moisture triggers (modes 2–4). The DVS and RD are taken from the WOFOST state; no extra input file is required.

---

## 5. Timeseries File Format

Several input files use a common timeseries format for time-varying data. This includes boundary condition timeseries (`spec_X_bcfile`), source/sink timeseries (`spec_X_ssfile`), and forcing data.

#### Format

```text
ndata <N>
<time_1> <value_1>
<time_2> <value_2>
...
<time_N> <value_N>
```

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `ndata` | int | — | Number of data points. Must appear on the first line before any data. |
| `time` | real | [s] | Time in seconds since simulation start. Must be monotonically increasing. |
| `value` | real | varies | Data value. Units depend on context (see below). |

#### Context-Specific Units

| File Type | Value Units | Example |
|-----------|-------------|---------|
| Boundary condition timeseries (`spec_X_bcfile`) | [mg/L] | `rtirrigationNH4.input` — NH4 concentration in irrigation water |
| Source/sink timeseries (`spec_X_ssfile`) | [mg/m²/s] | `rtfertNH4-DON-NH4.input` — area mass flux of NH4 fertilizer |
| Evapotranspiration timeseries | [m/s] | `et.input` — potential ET rate |
| Rainfall timeseries | [m/s] | `rainfall.input` — rainfall rate |

#### Interpolation

- **Boundary-condition timeseries** (`spec_X_bcfile`, types 4/5): always linear interpolation.
- **sstype = 1** (linear): Linear interpolation between consecutive data points.
- **sstype = 2** (step): Piecewise constant (step function). The value holds until the next time point.
- A **step profile** can also be achieved under linear interpolation by repeating the same value at consecutive points.
- **End behaviour**: beyond the last data point the value is **held constant** (a nonzero last value stays in effect until the end of the simulation — end the series at the simulation length); before the first point the first segment is extended (possible backward extrapolation — start the series at t = 0).

#### Example — Fertilizer Application (`rtfertNH4-DON-NH4.input`)

```text
ndata 126
0.00000e+00 0.00000e+00
561600.00    1.31250e+00
838800.00    1.12500e+00
1080000.00   0.00000e+00
3240000.00   1.68750e+00
3499200.00   0.00000e+00
...
95040000.0   0.00000e+00
```

(Rows 2–6 correspond to: day 6.5 NH4 pulse starts, day 9.7 application, day 12.5 pulse ends, day 37.5 second application, day 40.5 ends; last row = end of simulation.)

> **Note 1**: The `ndata` keyword is parsed by reading one token before the count. The parser ignores the label and reads the number that follows.
> **Note 2**: Data files must be **two plain numeric columns** — the parser reads tokens from a raw stream, so **comment lines and inline comments are not allowed** (the comment rules of `rtswbc.input` and other key:value files do not apply to timeseries data files).

---

## 6. Output Files and Headers

Output files are written to the configured output directory. Text time series use whitespace-separated columns; the first line is the header and subsequent lines contain values in scientific notation. `spec0`, `spec1`, etc. are zero-based species indices, not the strings in `SPECIES_Name`.

### 6.1 Reactive-transport text outputs

#### `RTSubsurfaceTimeSeries.out`

One row is written at initialization and at each reactive-transport output time. The exact header depends on the run configuration:

| Header pattern | Unit | Meaning |
|---|---|---|
| `Time[s]` | s | Simulation time. |
| `specX_LiquidMass[mg]` | mg | Total dissolved mass of species X in groundwater cells. |
| `specX_SolidMass[mg]` | mg | Total adsorbed mass; present when species X has adsorption. |
| `specX_TotalMass[mg]` | mg | Dissolved plus solid mass. |
| `specX_StorageChange[mg/s]` | mg/s | Storage-change rate for the transport step. |
| `specX_MassExch[mg/s]` | mg/s | Mass exchange diagnostic used by the transport integrator. |
| `specX_ExchangeIntoGW[mg/s]` | mg/s | Surface-water to groundwater exchange contribution. |
| `specX_BCInflow[mg/s]`, `specX_BCOutflow[mg/s]` | mg/s | Total reactive-transport boundary inflow and outflow. |
| `specX_RTBC..._Inflow[mg/s]`, `..._Outflow[mg/s]` | mg/s | Inflow/outflow split for each configured groundwater RT boundary. The header encodes boundary index, location, direction, and type. |
| `specX_SSInflow[mg/s]`, `specX_SSOutflow[mg/s]` | mg/s | Source/sink mass entering or leaving the groundwater solute system. |
| `specX_Reaction[mg/s]` | mg/s | Net reaction contribution; present when the reaction module is active. |
| `specX_BalanceResidual[mg/s]` | mg/s | Absolute mass-balance residual. |
| `specX_BalanceResidualRelThroughput[-]` | - | Residual normalized by throughput. |

When `SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS=1`, the additional columns `BalanceRHS`, `BalanceResidualRel`, `TransportHalf1Residual`, `TransportHalf3Residual`, and `NegativeClippedMassRate` are inserted before the final residual columns. Their units are encoded in the header emitted by the program.

#### `RTSurfaceTimeSeries.out`

For each species, the output contains `SurfaceMass[mg]`, `StorageChange[mg/s]`, `BCInflow[mg/s]`, `BCOutflow[mg/s]`, `SSInflow[mg/s]`, `SSOutflow[mg/s]`, and the final balance residual columns. `MassExch[mg/s]` appears when surface–subsurface exchange is enabled, `Reaction[mg/s]` when reactions are enabled, and `RainInflow[mg/s]` when rainfall forcing is enabled. Optional mass-conservation metrics follow the same compile-time switch described above.

### 6.2 Hydrodynamic and coupling text outputs

| File | Header fields and meaning |
|---|---|
| `domainTimeSeries.out` | `Time` [s], `SurfaceVolume` [m3], optional boundary inflow/outflow and accumulated volumes [m3/s or m3], optional per-boundary `BoundaryFlow_i` and `BoundaryAccum_i`, optional `RainFlux`/`RainAccum`, `SurfaceEvapFlux`/`SurfaceEvapAccum`, infiltration `InfFlux`/`InfAccum`, and optional `SurfaceSolidVolume`. The actual columns depend on boundaries, rainfall, infiltration, and compile-time sediment settings. |
| `SubsurfaceTimeSeries.out` | `Time` [s], `SubSurfaceVolume [m3]`, `ExchangeRate [m3/s]`, `BoundaryInflow`, `BoundaryOutflow`, `Source/SinkInflow [m3/s]`, `Source/SinkOutflow [m3/s]`, `Source/SinkOutflow_RootTransp[m3/s]`, and `Source/SinkOutflow_SoilEvap[m3/s]`. |
| `RootZoneWaterContentTimeSeries.out` | Time followed by one `RootZoneWaterContent_Mean [-]` value for each configured groundwater source/sink object with root-zone coupling. |

### 6.3 NetCDF outputs and variable headers

| File | Dimensions | Main variables and units |
|---|---|---|
| `output_subsurface.nc` | `t,z,y,x` | `t` [s], coordinates `x,y,z` [m], `hd` [m], `wc` [-], `qx,qy,qz` [m3/s or model flux convention], and, when ET is enabled, `evapGW` and `transpGW`. |
| `output_transport.nc` | `t,n,z,y,x` | `t` [s], `n` zero-based species index, `x,y,z` [m], `c` [mg/L], `c_solid` [mg/kg]. Debug builds may add `c_difxz`, `c_difzx`, `c_difxx`, and `c_difzz`. |
| `output_transportSW.nc` | `t,n,y,x` | `t` [s], `n` zero-based species index, `x,y` [m], and `csw` [mg/L]. |
| `crop.nc` | `t,y,x` | `t` [s], `x,y` [m], `DVS` [-], `LAI` [m2/m2], `RD` [cm], `TAGP/TWSO/TWLV/TWST/TWRT` [kg/ha], `RFTRA` [-], and `RZSM` [cm3/cm3]. |

The transport NetCDF writers store the current state at each output index. Variable attributes are limited in the current source, so the units above are the model conventions and should be used when post-processing.

---

## 7. Surface Water Hydrodynamics (SWE)

The surface-water module integrates the 2-D shallow-water equations on a uniform, DEM-based grid. The horizontal discretisation (`nx` × `ny` cells, cell size, lower-left corner) is defined **exclusively** by `dem.input`; every other raster input (roughness, initial conditions, solute initial fields) must carry an identical header. Boundary inputs specific to the triggered irrigation/drainage types (`bctype` 13/14) are described in [Section 4](#4-triggered-irrigation-and-drainage-boundaries); solute inputs in [Section 1](#1-surface-water-solute-transport).

A minimum surface-hydrodynamics configuration requires `parameters.input`, `dem.input` and `sw.input`. `extbc.input`, `rainfall.input`, `infiltration.input` and `wind.input` are optional and may simply be omitted.

### 7.1 Simulation Parameters — `parameters.input`

**File path**: `<input_folder>/parameters.input`
**Source code**: `src/Parser.h` → `readParamsFile()`

Colon-separated `key : value` file. Lines starting with `//` and trailing `//` comments are allowed. **Unrecognised keys abort the run**, and all required keys must be present.

#### Parameters

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `simLength` | real | **Yes** | Total simulation length [s]. |
| `endTime` | real | No | Explicit stop time [s]; when given, the run stops at `endTime` (the final step is shortened accordingly). |
| `cfl` | real | **Yes** | Target Courant number for the adaptive hydrodynamic time step (e.g. 0.5). |
| `parNx` | int | **Yes** | Number of MPI ranks in x (must match the total number of ranks, `parNx × parNy`). |
| `parNy` | int | **Yes** | Number of MPI ranks in y. |
| `outFreq` | real | **Yes** | Interval [s] between field (NetCDF/VTK) outputs. |
| `obsFreq` | real | **Yes** | Interval [s] between time-series outputs. |
| `nScreen` | int | No | Console progress print interval in time steps (default 100). |
| `outFormat` | string | **Yes** | Output format: `NETCDF`, `VTK` or `BIN`. |
| `BCtype` | string | **Yes** | Default condition on the four domain edges (not covered by `extbc.input` polygons): `REFLECTIVE`, `TRANSMISSIVE` or `PERIODIC`. |
| `writeRain`, `writeLandUse`, `writeSoilMap`, `writeRoughness`, `writeInfParameters` | int | No | Set to 1 to include the corresponding field in the output files. |

#### Example — `tests/t9-gaoyou-Ncycle-WOFOST/parameters.input`

```text
simLength  : 86400  // How long to simulate (s)
cfl        : 0.5      // CFL value to use
outFreq    : 3600      // Frequency for output in seconds
obsFreq    : 3600     // Frequency for time series output
nScreen    : 60
parNx      : 2       // Number of parallel tasks in the x-direction
parNy      : 2        // Number of parallel tasks in the y-direction
outFormat  : NETCDF   // Output format: NETCDF or VTK
BCtype     : REFLECTIVE // outer boundaries: TRANSMISSIVE, REFLECTIVE or PERIODIC
```

### 7.2 Terrain — `dem.input`

**File path**: `<input_folder>/dem.input`
**Source code**: `src/Parser.h` → `readHeaderDEMFile()` / `readDEMFile()`; `src/rasterTools.h` → `readRasterField()`

ESRI-ASCII-style raster. Six header records (one `label value` pair per line — **only the numeric token is read, the labels are not validated**), followed by the elevation data:

| Record | Value |
|---|---|
| 1 | Number of columns `nx` |
| 2 | Number of rows `ny` |
| 3 | x-coordinate of the lower-left corner `xll` [m] |
| 4 | y-coordinate of the lower-left corner `yll` [m] |
| 5 | Cell size `dx` [m] (square, uniform cells) |
| 6 | NODATA value |

Ordering and validation rules:

- The **first data row corresponds to the northernmost row** (largest y); rows proceed southward, and within a row values run from low x to high x (standard ESRI ASCII order).
- Cells equal to the NODATA value are deactivated (excluded from the computation).
- Every other raster read through `readRasterField()` (`roughness.input`, `hini.input`, `uini.input`, `vini.input`, solute initial grids) repeats the same six header records; its header values must **exactly match** `dem.input` (nx, ny, xll, yll, dx), otherwise the run aborts.

#### Example — `tests/t9-gaoyou-Ncycle-WOFOST/dem.input`

```text
n_cols 61
n_rows 221
xllcorner 0.0
yllcorber 0.0
cellsize 0.5
nodata_value -9999
7.310 7.310 7.309 ...   (61 values per row, 221 rows)
```

(Typos in the labels — such as `yllcorber` above — are tolerated because parsing is positional.)

### 7.3 Initial State and Friction — `sw.input`

**File path**: `<input_folder>/sw.input`
**Source code**: `src/Parser.h` → `readSWFile()`

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `initialMode` | string | **Yes** | Surface initial condition: `dry` (dry domain), `h` (constant depth), `h+z` (constant free-surface elevation; h = max(wse − z, 0)) or `file` (read `hini.input`, `uini.input`, `vini.input` rasters). Lower-case. |
| `initialValue` | real | If mode = h or h+z | Constant value for mode `h` (water depth [m]) or `h+z` (free-surface elevation [m]). |
| `friction` | string | **Yes** | Friction model: `none`, `manning`, `darcyweisbach` or `chezy`. Must match the compile-time `SERGHEI_FRICTION_MODEL` (this build: Manning). |
| `roughness` | string | If friction ≠ none | Roughness coefficient value (uniform over the domain), or `file` to read `roughness.input` (raster with the same header as the DEM), or `landuse` to derive it from the land-use map. |
| `dryDepth` | real | **Yes** | Depth threshold [m] below which velocities are zeroed (dry-cell treatment), e.g. 0.001. |

#### Example — `tests/t9-gaoyou-Ncycle-WOFOST/sw.input`

```text
initialMode :  h+z      // dry, h, h+z, file
initialValue : 7.4     // depth or water surface elevation constant value

friction : manning       // none, manning
roughness : 0.15         // roughness coefficient value, or set to 'file' to read from file

dryDepth : 0.001
```

### 7.4 External Boundaries — `extbc.input`

**File path**: `<input_folder>/extbc.input`
**Source code**: `src/Parser.h` → `readExtBCFile()`; boundary application in `src/BC.h`

Defines external (domain-edge) surface boundaries as polygons. If the file is absent, the global `BCtype` from `parameters.input` applies to every edge. The file starts with `bccount` (number of boundary blocks) followed by one block per boundary; each block starts with `id`:

| Key | Type | Description |
|-----|------|-------------|
| `bccount` | int | Number of boundary blocks (must be first). |
| `id` | string | Boundary identifier (starts a new block). |
| `bctype` | int | Boundary type code (table below). |
| `direction` | real ×2 | Outward normal components `(normalx, normaly)`, e.g. `0.0 1.0`; the vector is normalised internally and fixes the orientation of prescribed fluxes. |
| `polygon` | string | 2-D polygon file selecting the boundary cells (`NPOINTS n` + `x y` rows). |
| `bcvals` | real(s) | Constant boundary value(s); meaning depends on `bctype` (table below). |
| `hydrograph` | string | Timeseries file for types 10/11/12 (format in [Section 5](#5-timeseries-file-format); plain numbers, no comments). |
| `trigger_polygon`, `irrigation_strategy`, `drainage_strategy` | string | Only for `bctype` 13/14 — see [Section 4](#4-triggered-irrigation-and-drainage-boundaries). |

Unrecognised keys abort the run, the number of blocks must equal `bccount`, and each polygon must select at least one cell.

#### Boundary condition types

| Code | Name | Value source | Description |
|------|------|--------------|-------------|
| 1–3 | `PERIODIC` / `REFLECTIVE` / `TRANSMISSIVE` | — | Global edge types; only valid for `BCtype` in `parameters.input`, not as `extbc` block types. |
| 5 | `CRITICAL` | — | Critical-flow outflow (Froude = 1 along the outward normal). |
| 6 | `H_CONST` | `bcvals : h` | Constant water depth [m]. |
| 7 | `Q_CONST` | `bcvals : <unused> qx qy` | Constant unit discharge [m²/s]: the **2nd and 3rd** values give the x/y components; the boundary depth is raised to the critical depth if the imposed discharge is supercritical. |
| 8 | `WSE_CONST` | `bcvals : wse` | Constant free-surface elevation [m]. |
| 9 | `FREE_OUTFLOW` | — | Zero-gradient free outflow. |
| 10 | `HZ_T_INLET` | `hydrograph` = WSE(t) [m] | Stage hydrograph at an inlet; boundary depth = max(wse(t) − z, 0). |
| 11 | `HZ_T_OUTLET` | `hydrograph` = WSE(t) [m] | Stage hydrograph at an outlet; reverse flow (e.g. tidal inflow) is permitted. |
| 12 | `Q_T` | `hydrograph` = Q(t) [m³/s] | Discharge hydrograph distributed over the boundary cells; positive values act along `direction`. |
| 13 | `TRIGGERED_IRRIGATION` | strategy file | Rule-based irrigation — [Section 4](#4-triggered-irrigation-and-drainage-boundaries). |
| 14 | `TRIGGERED_DRAINAGE` | strategy file | Rule-based drainage — [Section 4](#4-triggered-irrigation-and-drainage-boundaries). |

#### Example — `tests/t9-gaoyou-Ncycle-WOFOST/extbc.input` (irrigation inlet and outlet, both discharge hydrographs)

```text
bccount : 2
id : irrigation
bctype : 12
direction : 0.0 1.0
polygon : polygoninlet.input
hydrograph: inlet_data_rate-2-input.input
id : outlet
bctype : 12
direction : 0.0 -1.0
polygon : polygonoutlet.input
hydrograph: outlet_data_rate-2-input.input
```

### 7.5 Rainfall — `rainfall.input` (optional)

**Source code**: `src/Parser.h` → `readRainfallFile()`

If the file is absent the run is simply rain-free. Header:

```text
timeUnits s          // s or h
rainUnits mm/s       // mm/s or mm/h
TimeSeriesNp <np>    // number of time points
nx <nx>              // rain-field partitions in x
ny <ny>              // rain-field partitions in y
```

followed by `np` rows of `time` (converted to [s]) and `nx × ny` rainfall values (auto-converted to [m/s]). With `nx = ny = 1` the rainfall is uniform over the domain; for a distributed rain field the partitions are ordered x-fastest (partition `k = j*nx + i` of the DEM grid), each partition occupying one column of the table.

### 7.6 Infiltration — `infiltration.input` (optional)

**Source code**: `src/Parser.h` → `readInfiltrationFile()`

If the file is absent the surface is impervious (`INF_NONE`). **In SW–GW coupled runs (this build) keep `infmodel : none`** — infiltration is then driven by the exchange boundary of [Section 8.4](#84-boundary-conditions--gwbcinput), and an additional surface infiltration model would double-count the loss.

| Key | Type | Description |
|-----|------|-------------|
| `infmodel` | string | `none`, `constant`, `horton` or `greenampt`. |
| `spatialInput` | string | Spatial distribution of the parameters: `soilmap`, `landuse`, `netcdf` or `raster`. |
| `infiltrationclasses` | int | Number of parameter classes (label 0 is always impervious). |
| `rate` | real | Constant model: infiltration capacity per class [mm/s] (auto-converted to m/s). |
| `k`, `fc`, `f0` | real | Horton model: decay rate and initial/final capacities per class [mm/s]. |
| `ks`, `psi`, `dtheta` | real | Green-Ampt model: saturated conductivity [m/s], suction head [m], water-content deficit [-]. |

### 7.7 Wind — `wind.input` (optional)

**Source code**: `src/Parser.h` → `readWindFile()`

Three header records — number of data points, wind-stress coefficient `CwT`, and minimum depth `hwmin` for wind stress — followed by rows of `t vx vy` (wind velocity components [m/s]). If the file is absent, no wind forcing is applied.

---

## 8. Groundwater Hydrodynamics (Subsurface Flow)

The subsurface module integrates the 3-D mixed-form Richards equation on a terrain-following mesh that reuses the horizontal DEM grid (`nx` × `ny`) plus `ndepth` vertical layers hanging below the topography. It is compiled in with `-DSERGHEI_SWE_GW=ON` (macro `SERGHEI_RE_MODEL`; enabled in this build) and is then **always active**: `subsurface.input` is required. Its outputs (`SubsurfaceTimeSeries.out`, `output_subsurface.nc`) are listed in [Section 6](#6-output-files-and-headers).

| File | Role |
|---|---|
| `subsurface.input` | Vertical discretisation, soil count, initial mode, solver and coupling switches (**required**) |
| `vg.input` | van Genuchten–Mualem parameters per soil type |
| `soilID.input` | 3-D soil-type index field |
| `gwbc.input` | Subsurface boundary conditions (optional; uncovered faces default to no-flow) |
| `gwss.input` | Subsurface source/sink zones: ET, recharge, drainage, tile drains (optional) |
| `head.input` / `theta.input` / `wt.input` | Initial condition fields when `initialMode` = 2/3/4 |

### 8.1 Domain and Solver — `subsurface.input`

**File path**: `<input_folder>/subsurface.input`
**Source code**: `src/GwInit.h` → `readGwFile()`

Colon-separated `key : value` file with `//` comments. All required keys must be present or the run aborts; unknown keys are **silently ignored** (e.g. the legacy `BCtype` line that appears in older test cases has no effect).

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `ndepth` | int | **Yes** | Number of vertical cells (`nz`). |
| `height` | real | **Yes** | Domain thickness [m]; with `dz_multiplier : 1.0`, `dz = height / ndepth`. |
| `dz_multiplier` | real | **Yes** | Vertical stretching ratio; `1.0` = uniform `dz`; `> 1` = geometrically growing `dz` (the actual thickness then exceeds `height`). |
| `dz_base` | real | If multiplier > 1 | Finest (top) cell thickness [m] for the stretched grid. |
| `nSoilID` | int | **Yes** | Number of soil types = number of entries per key in `vg.input`. |
| `initialMode` | int | **Yes** | Initial condition: `1` saturated, `2` `head.input`, `3` `theta.input`, `4` `wt.input` — see [Section 8.6](#86-initial-conditions). |
| `dt_init` | real | **Yes** | Initial (and minimum) subsurface time step [s]. |
| `dt_max` | real | **Yes** | Maximum subsurface time step [s]. |
| `gw_scheme` | int | **Yes** | Nonlinear scheme: `1` = modified Picard (PCA), `2` = Picard iteration. |
| `cg_iter` | int | No | Maximum CG iterations (default 10⁶). |
| `cg_tol` | real | No | CG convergence tolerance (default 1e-8). |
| `aev` | real | **Yes** | Air-entry value [m] of the modified van Genuchten model, used to smooth conductivity/retention near saturation; `0.0` disables the smoothing. |
| `async` | int | **Yes** | SW–GW coupling mode: `1` = asynchronous (the groundwater module marches with its own adaptive `dt` between `dt_init` and `dt_max`, exchanging with the surface whenever it catches up; the surface step is capped at the groundwater step), `0` = synchronous (both modules use the same, smaller `dt` each step). |

#### Example — `tests/t9-gaoyou-Ncycle-WOFOST/subsurface.input`

```text
BCtype  : REFLECTIVE  // legacy, ignored
height : 0.4          // Thickness of the subsurface domain
dz_multiplier : 1.0   // Non-uniform dz
dz_base : 0.05
ndepth  : 8           // Number of cells in z-direction
nSoilID : 4           // Number of soil types
initialMode : 1       // saturated(1), file-h(2), file-theta(3), file-wt(4)
dt_init : 0.01        // initial (also minimum) dt for subsurface
dt_max  : 100.0       // maximum allowable dt
gw_scheme : 1         // numerical schemes, 1:PCA, 2:Picard
async : 1             // Asynchronous coupling
aev : 0.0
```

### 8.2 Soil Hydraulic Parameters — `vg.input`

**File path**: `<input_folder>/vg.input`
**Source code**: `src/GwInit.h` → `readVGParameters()`

One key per van Genuchten–Mualem parameter; the values for the `nSoilID` soil types are given as a **semicolon-separated list** (no spaces):

| Key | Units | Description |
|-----|-------|-------------|
| `Ks` | [m/s] | Saturated hydraulic conductivity. |
| `Phi` | [-] | Porosity (may equal `ThetaS`). |
| `ThetaS` | [-] | Saturated water content. |
| `ThetaR` | [-] | Residual water content. |
| `n` | [-] | van Genuchten `n` (> 1); `m = 1 − 1/n` is computed internally. |
| `alpha` | [1/m] | van Genuchten `alpha`. |

The soil index of each cell (from `soilID.input`) selects the column. **Always provide `vg.input`** — the fallback used when the file is missing relies on legacy scalar keys that are no longer parsed.

#### Example — `tests/t9-gaoyou-Ncycle-WOFOST/vg.input`

```text
alpha   : 0.56;0.68;1.69;1.93
n       : 1.65;1.61;1.33;1.31
Ks      : 4.86e-6;2.94e-6;5.61e-7;4.88e-7
Phi     : 0.51;0.47;0.43;0.39
ThetaR  : 0.048;0.046;0.036;0.035
ThetaS  : 0.51;0.47;0.43;0.39
```

### 8.3 Soil Map — `soilID.input`

**File path**: `<input_folder>/soilID.input`
**Source code**: `src/GwInit.h` → `readSoilID()`

First record `n_soil <N>` (must equal `nSoilID`), followed by `nx × ny × nz` integers (0-based soil indices referencing the `vg.input` columns). Ordering is **x fastest, then y, then z**: for each layer, DEM row-major from north to south; layers from top to bottom. When `nSoilID : 1` the data block (and the file itself) may be omitted — all cells then use soil type 0.

### 8.4 Boundary Conditions — `gwbc.input`

**File path**: `<input_folder>/gwbc.input`
**Source code**: `src/GwInit.h` → `readGwBCFile()`; boundary application in `src/GwBC.h`

Optional file; subsurface faces not covered by a boundary block are impermeable (no-flow). Structure as in `extbc.input`: `bccount` followed by blocks starting with `id`. **`direction` must appear before `bcvals`/`bcfile`** (it fixes the size of the boundary data array).

| Key | Type | Description |
|-----|------|-------------|
| `bccount` | int | Number of boundary blocks (first). |
| `id` | string | Boundary identifier (starts a block). |
| `direction` | int | Boundary face: `1` +x (east), `2` −x (west), `3` +y (north), `4` −y (south), `5` bottom (domain base), `6` top (ground surface). |
| `bctype` | int | Type code (table below). |
| `polygon` | string | 2-D polygon file; the polygon is projected onto the selected face. |
| `bcvals` | real | Constant value broadcast to every boundary face; meaning per type. |
| `timeseries` | string | Time-series file for types 5/6/7. |
| `bcfile` | string | Spatially distributed values for types 2/3/4. |

#### Boundary condition types

| Code | Name | Value | Description |
|------|------|-------|-------------|
| 1 | `NOFLOW` | — | Impermeable face. |
| 2 | `H_CONST` | `bcvals : h` [m] | Constant pressure head (Dirichlet). |
| 3 | `Q_CONST` | `bcvals : q` [m/s] | Constant Darcy flux (Neumann). |
| 4 | `WT_CONST` | `bcvals : wt` [m] | Constant water-table elevation; head = wt − z. |
| 5 | `H_T` | `timeseries` h(t) [m] | Pressure-head time series (linear interpolation). |
| 6 | `Q_T` | `timeseries` q(t) [m/s] | Flux time series. |
| 7 | `WT_T` | `timeseries` wt(t) [m] | Water-table time series; supports per-cell columns. |
| 8 | `SWE` | — | **Surface–groundwater exchange**; `direction` must be `6`. The ghost head equals the surface water head; three exchange regimes are distinguished (no exchange / infiltration / exfiltration capped at the saturated conductivity `Ks`). This is the standard top boundary of coupled runs. |
| 9 | `FD` | — | Free drainage (unit gravity gradient, q = −k); `direction` must be `5`. Standard bottom boundary. |

Timeseries files (types 5/6/7) start with `ndata <np> <ncell>`; each row is `t value` (`ncell = 1`) or `t v1 … v_ncell` (per-cell values, used by `WT_T`). `bcfile` files start with `ndata <nx> <ny> <nz>` followed by the values, with dimensions matching the face (`direction` 1/2: `ny × nz`; 3/4: `nx × nz`; 5/6: `nx × ny`).

#### Example — `tests/t9-gaoyou-Ncycle-WOFOST/gwbc.input` (top SWE exchange + bottom free drainage; lateral fixed heads commented out)

```text
bccount : 2
id : surface
direction : 6
bctype : 8
bcvals : 0.0
polygon : polygontop.input
id : bottom
direction : 5
bctype : 9
polygon : polygontop.input

//id : left
//direction : 2
//bctype : 2
//bcvals : 0
//polygon : polygonleft.input
//bcfile : headleft.input
```

### 8.5 Source/Sink Zones — `gwss.input`

**File path**: `<input_folder>/gwss.input`
**Source code**: `src/GwInit.h` → `readGwSSFile()`; application in `src/SourceSink.h`

Optional file (absent → no subsurface source/sink; rename the file to disable it, as done in `tests/t15`). Structure: `sscount` followed by blocks starting with `id`:

| Key | Type | Description |
|-----|------|-------------|
| `sscount` | int | Number of source/sink blocks (first). |
| `id` | string | Identifier (starts a block). |
| `sstype` | int | Type code (table below). |
| `polygon` | string | **3-D** polygon file: `NPOINTS n` + `x y z` rows; `z` is negative below the ground surface (e.g. `-0.3` = 0.3 m deep). 4 vertices = footprint + bottom depth (the whole column above that depth is affected); 8 vertices = closed box between two depths (used for tile drains). |
| `ssvals` | real | Constant value (type 2: fixed head [m]). |
| `timeseries` | string | Type 0: ET file (special format below); type 1: `ndata` + `t q` file. |
| `TileDrainageFile` | string | Type 3: drain parameter file (below). |
| `shapecorrection` | real | Type 2: conductance multiplier `Cpipe` [-]. |
| `direction` | int | Type 2: axis of the drainage specialisation (1–4; default z). |

#### Source/sink types

| Code | Type | Description |
|------|------|-------------|
| 0 | Evapotranspiration | Root-zone water uptake — see the ET file format below. |
| 1 | Flux | Volumetric source/sink rate [1/s] (change of water content per unit time; positive = recharge), from the `timeseries` file. |
| 2 | Head drainage | Internal drainage towards a fixed head `ssvals` with conductance factor `shapecorrection`. |
| 3 | Tile drainage | Drain sink `Γ·h` (active only where the pressure head is positive), parameters from `TileDrainageFile`. |

#### ET timeseries file (`timeseries` of `sstype : 0`)

```text
ndata <N>                     // number of data rows
f <0.5>                       // Beer's-law extinction coefficient for the E/T split
stress <h1> <h2> <h3> <h4>    // Feddes stress-reduction thresholds, pressure head [m]
rootx <xs> <ys> <zs>          // root-distribution parameters (parsed; see note)
rootp <px> <py> <pz>          // root-distribution parameters (parsed; see note)
<t> <ET0> <LAI> <RD>          // N rows: time [s], potential ET [m/s], LAI [-], rooting depth [m]
```

`E_pot = ET0·exp(−f·LAI)`, `T_pot = ET0 − E_pot`; transpiration is distributed over the root zone using the Feddes stress function α(h) between `h1…h4` and a normalised root-density shape. **In `CROP_GROWTH_MODEL` builds (this repository) LAI and RD are taken from WOFOST each step and only ET0 comes from this file** (hybrid coupling); the LAI/RD columns are then placeholders. The root-density shape currently uses fixed parameters (`zs = 0.2·RD`, `pz = 1`) — the `rootx`/`rootp` records are read for compatibility but do not shape the current distribution. The potential soil evaporation is also published to the surface module (evaporation of ponded water where present).

#### Tile-drain parameter file (`TileDrainageFile`, e.g. `drains.input`)

| Key | Units | Description |
|-----|-------|-------------|
| `diameter` | [m] | Effective drain diameter `d_e`. |
| `scale_factor` | [-] | Reduction factor for the grid-size effect (typically 1, 2 or 4). |

The dimensionless resistance correction `Cd` is computed internally from `d_e`, `scale_factor` and the grid size after the HYDRUS formula.

#### Example — `tests/t9-gaoyou-Ncycle-WOFOST/gwss.input`

```text
sscount : 1
id : root1
sstype : 0
ssvals : 0.0
polygon : polygonroot.input
timeseries : et.input

//id : TileDrains
//sstype : 3
//polygon : polygondrains.input
//TileDrainageFile : drains.input
```

### 8.6 Initial Conditions

**Source code**: `src/GwInit.h` → `setGwState()` / `readGwICFile()`

Selected by `initialMode` in `subsurface.input`:

| Mode | Source | Description |
|------|--------|-------------|
| `1` | — | Saturated: `wc = ThetaS`; the pressure head continues the surface water column (`h = h_surface + depth below the surface`). |
| `2` | `head.input` | 3-D pressure-head field [m]. |
| `3` | `theta.input` | 3-D water-content field [-]; heads derived via van Genuchten. |
| `4` | `wt.input` | 2-D water-table elevation [m]; `h = wt − z` per cell; NODATA cells are set saturated. |

File format (`head.input` / `theta.input` / `wt.input`): three header records (`nx`, `ny`, NODATA value — labels not validated), then the values with **x fastest, then y, then z** (2-D: `nx × ny`, north-to-south rows as in the DEM; 3-D: one `nx × ny` block per layer, top → bottom). The dimensions must match the DEM, and the file **must exist** for modes 2/3/4.

### 8.7 Surface–Subsurface Coupling (summary)

- With the `SWE` top boundary ([Section 8.4](#84-boundary-conditions--gwbcinput)) the exchange flux is computed from the surface and groundwater heads every groundwater step; exfiltration is capped at the saturated conductivity, and the exchanged water is applied to the surface depth directly. The exchange rate is reported as `ExchangeRate` in `SubsurfaceTimeSeries.out` ([Section 6.2](#62-hydrodynamic-and-coupling-text-outputs)).
- `async : 1` decouples the time steps (recommended; used by all reference cases); `async : 0` forces identical steps.
- Rainfall is always a surface forcing in this build; evaporation is split between ponded surface water and the top soil layer via the ET zone of [Section 8.5](#85-sourcesink-zones--gwssinput).

---

## Appendix A: Species Ordering Convention

The nitrogen cycle model uses a fixed species ordering. The species index (0-based) must be consistent across all input files.

| Index | Name | Chemical Formula | Role in Nitrogen Cycle |
|:---:|------|-----------------|----------------------|
| 0 | NH4 (Ammonium) | NH₄⁺ | Produced by mineralization; consumed by nitrification |
| 1 | NO3 (Nitrate) | NO₃⁻ | Produced by nitrification; consumed by denitrification |
| 2 | DO (Dissolved Oxygen) | O₂ (aq) | Consumed by nitrification and aerobic respiration; replenished by reaeration |
| 3 | DOC (Dissolved Organic Carbon) | C (org) | Electron donor for denitrification; consumed by aerobic respiration; released from sediment |
| 4 | DON (Dissolved Organic Nitrogen) | N (org) | Produced by DOC/DON release; consumed by mineralization to NH4 |

### Nitrogen Cycle Pathway

```
                    Mineralization
    DON ──────────────────────────→ NH4
                                      │
                        Nitrification │ (requires O2)
                                      ▼
    O2 ←── Reaeration             NO3
                                      │
                      Denitrification │ (consumes DOC, inhibited by O2)
                                      ▼
                                     N2 (lost)

    DOC + O2 ──→ CO2    (Aerobic Respiration)
    DOC ←── Sediment Release
    O2  ←── Atmospheric Reaeration
```

---

## Appendix B: Input File Quick-Reference Checklist

The following table lists all input files of the coupled model (hydrodynamics, solute transport, crop growth). Check that all required files are present in the input directory before running a simulation.

### Surface Water Hydrodynamics

| File | Required | Purpose |
|------|----------|---------|
| `parameters.input` | **Yes** | Simulation control: length, CFL, MPI layout, output frequency/format ([Section 7.1](#71-simulation-parameters--parametersinput)) |
| `dem.input` | **Yes** | Terrain raster; defines the horizontal grid ([Section 7.2](#72-terrain--deminput)) |
| `sw.input` | **Yes** | Initial state, friction, roughness ([Section 7.3](#73-initial-state-and-friction--swinput)) |
| `extbc.input` | No* | External surface boundaries, incl. `bctype : 13`/`14` triggered irrigation/drainage ([Sections 4](#4-triggered-irrigation-and-drainage-boundaries) and [7.4](#74-external-boundaries--extbcinput); *required if boundaries other than the global `BCtype` are used) |
| `roughness.input`, `hini.input`, `uini.input`, `vini.input` | If referenced in `sw.input` | Spatial roughness / initial-state rasters ([Section 7.3](#73-initial-state-and-friction--swinput)) |
| `rainfall.input` | No | Rainfall forcing ([Section 7.5](#75-rainfall--rainfallinput-optional)) |
| `infiltration.input` | No | Surface infiltration model — keep `none` in SW–GW coupled runs ([Section 7.6](#76-infiltration--infiltrationinput-optional)) |
| `wind.input` | No | Wind forcing ([Section 7.7](#77-wind--windinput-optional)) |

### Surface Water Solute Transport

| File | Required | Purpose |
|------|----------|---------|
| `rttransportsw.input` | **Yes** | Species count, transport properties, initial conditions, solver settings |
| `rtswbc.input` | No* | Solute boundary conditions (*required if external boundaries exist) |
| `rtswss.input` | No | Solute source/sink zones (fertilizer, etc.) |
| `reactionsw.input` | If ReactionModule=1 | Reaction kinetics, nitrogen cycle parameters |
| `rtinitialsw*.input` | If aq_mode=1 | Spatial initial concentration grids for each species |

### Groundwater Hydrodynamics

| File | Required | Purpose |
|------|----------|---------|
| `subsurface.input` | **Yes** | Discretisation, solver, coupling switches, initial mode ([Section 8.1](#81-domain-and-solver--subsurfaceinput); the module is always active in this build) |
| `vg.input` | **Yes** | van Genuchten–Mualem parameters per soil type ([Section 8.2](#82-soil-hydraulic-parameters--vginput)) |
| `soilID.input` | If nSoilID > 1 | 3-D soil-type index field ([Section 8.3](#83-soil-map--soilidinput)) |
| `gwbc.input` | No | Subsurface boundaries: SWE top exchange, FD bottom, heads/fluxes ([Section 8.4](#84-boundary-conditions--gwbcinput)) |
| `gwss.input` | No | ET / recharge / drainage / tile-drain zones ([Section 8.5](#85-sourcesink-zones--gwssinput)) |
| `head.input` / `theta.input` / `wt.input` | If initialMode = 2/3/4 | Initial-condition fields ([Section 8.6](#86-initial-conditions)) |
| `et.input` | If gwss sstype = 0 | ET0 series with Feddes/root parameters ([Section 8.5](#85-sourcesink-zones--gwssinput)) |
| `drains.input` | If gwss sstype = 3 | Tile-drain parameters ([Section 8.5](#85-sourcesink-zones--gwssinput)) |

### Groundwater Solute Transport

| File | Required | Purpose |
|------|----------|---------|
| `rttransportgw.input` | **Yes** | Species count, transport properties, initial conditions, solver settings |
| `rtgwbc.input` | No* | Solute boundary conditions |
| `rtgwss.input` | No | Solute source/sink zones (root uptake, etc.) |
| `reactiongw.input` | If ReactionModule=1 | Reaction kinetics, sorption, nitrogen cycle |
| `rtinitialgw*.input` | If aq_mode=1 | Spatial liquid-phase initial concentration |
| `rtinitialgwsolid*.input` | If solid_mode=1 | Spatial solid-phase (adsorbed) initial concentration |

### WOFOST Crop Growth

| File | Required | Purpose |
|------|----------|---------|
| `cropparameter.input` | **Yes** | Crop parameters, soil properties, nutrient parameters |
| `cropmeteo.input` | **Yes** | Daily meteorological forcing data |
| `agro.input` | **Yes** | Crop variety, planting/harvest dates |

### Supporting Files

| File | Required | Purpose |
|------|----------|---------|
| `extbc.input` | If surface boundaries used | Hydraulic external boundary conditions ([Section 7.4](#74-external-boundaries--extbcinput)); `bctype : 13`/`14` configure triggered irrigation/drainage (see [Section 4](#4-triggered-irrigation-and-drainage-boundaries)) |
| Irrigation/drainage strategy files (e.g., `test.input`, `test2.input`) | If bctype=13/14 | Stage tables for triggered irrigation/drainage (see [Section 4](#4-triggered-irrigation-and-drainage-boundaries)) |
| `polygon*.input` | As referenced | Polygon vertex coordinates for BC/SS spatial zones and trigger observation regions |
| Timeseries files (e.g., `rtfert*.input`) | As referenced | Time-varying boundary/source data |

### Main Output Files

| File | Produced when | Purpose |
|------|---------------|---------|
| `domainTimeSeries.out` | Surface hydrodynamics enabled | Domain-scale surface-water volume, boundary flow, rainfall, evaporation, and infiltration diagnostics. |
| `SubsurfaceTimeSeries.out` | Groundwater hydrodynamics enabled | Domain-scale groundwater volume, exchange, boundary flow, and source/sink diagnostics. |
| `RTSurfaceTimeSeries.out` | Surface solute transport enabled | Per-species surface solute mass and mass-balance diagnostics. |
| `RTSubsurfaceTimeSeries.out` | Groundwater solute transport enabled | Per-species liquid/solid/total mass, boundary/source/reaction terms, and mass-balance diagnostics. |
| `RootZoneWaterContentTimeSeries.out` | Root-zone coupling enabled | Root-zone water content for configured coupling objects. |
| `output_transportSW.nc` | Surface solute transport enabled | Time-dependent 2-D dissolved concentration fields `csw[t,n,y,x]`. |
| `output_transport.nc` | Groundwater solute transport enabled | Time-dependent 3-D dissolved and solid concentration fields `c[t,n,z,y,x]` and `c_solid[t,n,z,y,x]`. |
| `crop.nc` | WOFOST enabled | Time-dependent crop development, biomass, rooting depth, stress, and root-zone moisture fields. |

See [Section 6](#6-output-files-and-headers) for the complete text headers and NetCDF variable list. Optional debug and mass-conservation columns are controlled by compile-time flags and are therefore not present in every run.

---

*End of User Manual — Coupled Surface–Subsurface Hydrodynamics, Solute Transport & WOFOST Crop Growth*
