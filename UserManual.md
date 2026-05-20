# SERGHEI User Manual — Solute Transport & WOFOST Crop Growth Coupling

> **Scope**: This manual covers the **solute (reactive) transport** module (surface water and groundwater) and the **WOFOST 7.2 crop growth** module within SERGHEI. Hydrodynamics input files are not documented here.
>
> **Reference test case**: `tests/t9-gaoyou-Ncycle-WOFOST/` (paddy rice with nitrogen cycle)

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
- [4. Timeseries File Format](#4-timeseries-file-format)
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
| `aq_val` | real | If `aq_mode=0` | Uniform initial concentration [mg/L or kg/m³]. |
| `aq_file` | string | If `aq_mode=1` | Filename for spatially distributed initial concentration (ESRI ASCII Grid or simplified format; see [Section 2.5](#25-initial-condition-files)). |
| `diffusion_molecular` | real | **Yes** | Molecular diffusion coefficient [m²/s]. Typical range: 1.0×10⁻¹⁰ to 1.0×10⁻⁸. |
| `alpha_L` | real | **Yes** | Longitudinal dispersivity [m]. Typical range: 0.01 to 100. |
| `alpha_T` | real | **Yes** | Transverse dispersivity [m]. Typical range: 0.001 to 10. Often 1/10 of `alpha_L`. |
| `RainCon` | real | **Yes** | Solute concentration in rainfall for this species [mg/L or kg/m³]. Set to `0.0` if rainfall carries no solute. |

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
**Source code**: `src/RTInitSW.h` → `readRTBCFileSW()`, `src/BC.h` (lines 48–52)

Defines solute boundary conditions at external boundaries for each species. Boundary geometry (cell indices) is inherited from the hydrodynamic boundary conditions (`extbc.input`).

#### Boundary Condition Types

| Value | Constant Name | Description |
|:---:|------|-------------|
| 1 | `RTSW_BC_NONE` | No reactive transport BC applied at this boundary. |
| 2 | `RTSW_BC_CON_CONST` | **Constant concentration** (Dirichlet). The boundary cells are held at a fixed concentration specified by `spec_X_bcval`. |
| 3 | `RTSW_BC_ZEROGRAD` | **Zero-gradient** (Neumann, free outflow). The concentration at the boundary is extrapolated from the interior (zero concentration gradient normal to the boundary). Used for free-discharge outlets. |
| 4 | `RTSW_BC_CON_T` | **Time-varying concentration** (Dirichlet time series). The boundary concentration varies in time according to an external timeseries file specified by `spec_X_bcfile`. See [Section 4](#4-timeseries-file-format) for timeseries format. |

#### File Format

```text
bccount : <number_of_boundaries>

id : <boundary_name>            // Matches a boundary in extbc.input
direction : <nx> <ny>           // Normal vector (will be normalized to unit length)

// Per-species BC (X = species index, 0-based)
spec_X_bctype : <type>          // BC type (1–4, see table above)
spec_X_bcval  : <value>         // Constant concentration [mg/L] (only for type 2)
spec_X_bcfile : <filename>      // Timeseries file (only for type 4)
```

#### Cross-Validation Rules

- **Type 2** (constant): Use `spec_X_bcval`. Do **not** specify `spec_X_bcfile`.
- **Type 4** (time-varying): Use `spec_X_bcfile`. Do **not** specify `spec_X_bcval`.
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
| `Rate_Max_Nit` | real | [mg NH4/L/s] | 1×10⁻⁵ – 5×10⁻⁵ | Maximum nitrification rate. Corresponds to soil rate ~4–12 mg N/kg/d. Rice field measured range: 0.5–17 mg/kg/d (lower during flooding). |
| `K_Monod_NH4` | real | [mg NH4-N/L] | 0.5 – 2.0 | Half-saturation constant for NH4 in Monod kinetics. |
| `K_Monod_DO` | real | [mg O2/L] | 0.1 – 2.0 | Half-saturation constant for dissolved oxygen. Set very small (e.g., 0.01) to effectively ignore DO limitation during flooded conditions. |

#### Denitrification (NO3⁻ → N2)

Reduction of nitrate to nitrogen gas by heterotrophic bacteria under anaerobic conditions. This is the dominant nitrogen loss pathway in flooded rice paddies.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Denit` | real | [mg NO3/L/s] | 5×10⁻⁶ – 3×10⁻⁴ | Maximum denitrification rate. Corresponds to soil rate ~10–25 mg N/kg/d. Rice field measured range: 6–18 mg/kg/d (higher during flooding). |
| `K_Monod_NO3` | real | [mg NO3-N/L] | 0.2 – 1.0 | Half-saturation constant for NO3. |
| `K_Monod_DOC` | real | [mg DOC/L] | 0.5 – 10.0 | Half-saturation constant for dissolved organic carbon (electron donor). Set small to ignore DOC limitation. |
| `Ki_Inhib_DO` | real | [mg O2/L] | 0.1 – 2.0 | DO inhibition constant. Higher DO inhibits denitrification. Set small (e.g., 0.01) to allow denitrification even with low DO under flooded conditions. |

#### Aerobic Respiration (DOC + O2 → CO2)

Heterotrophic consumption of dissolved organic carbon with oxygen as electron acceptor.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Hetero` | real | [mg DOC/L/s] | 5×10⁻⁵ – 2×10⁻⁴ | Maximum aerobic respiration rate. Surface water heterotrophic respiration is typically more active than soil. |
| `K_Monod_DOC_Aerobic` | real | [mg DOC/L] | 0.5 – 10.0 | Half-saturation constant for DOC. |
| `K_Monod_DO_Aerobic` | real | [mg O2/L] | 0.1 – 2.0 | Half-saturation constant for DO. |

#### Mineralization (DON → NH4⁺)

Decomposition of dissolved organic nitrogen to ammonium.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Min` | real | [mg DON/L/s] | 1×10⁻⁷ – 1×10⁻⁵ | Maximum mineralization rate. |
| `K_Monod_DON` | real | [mg DON-N/L] | 0.1 – 5.0 | Half-saturation constant for DON. |
| `K_Monod_DO_Min` | real | [mg O2/L] | 0.1 – 2.0 | Half-saturation constant for DO. |

#### Surface-Water-Specific Source/Sink Parameters

These parameters govern atmospheric reaeration and sediment-water solute exchange processes unique to surface water.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Ka_sw` | real | [1/s] | 0 – 1×10⁻⁴ | Atmospheric reaeration rate constant. Set to `0.0` to disable reaeration for validation cases. Typical: 5×10⁻⁵. |
| `DO_sat` | real | [mg/L] | 8 – 12 | Saturated dissolved oxygen concentration. Temperature-dependent: ~8.1 mg/L at 25–28°C. |
| `DOC_eq_sw` | real | [mg C/L] | 2 – 20 | Equilibrium DOC concentration in surface water (sediment release endpoint). Fertile paddy fields: ~10 mg/L. |
| `DON_eq_sw` | real | [mg N/L] | 0.1 – 5.0 | Equilibrium DON concentration in surface water (sediment release endpoint). |
| `K_rel_sw` | real | [1/s] | 0 – 1×10⁻⁵ | Release rate constant for sediment-water solute exchange. Set to `0.0` to disable. |

#### Environmental Correction Factors

Temperature and pH correction for all biological reaction rates.

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Temp_Coeff_Theta` | real | [-] | 1.02 – 1.10 | Arrhenius temperature correction coefficient (θ). Universal for nitrogen transformations. |
| `Opt_Temp` | real | [°C] | 20 – 28 | Optimal (reference) temperature for biological reactions. Rice growing season: 25°C. |
| `Opt_pH` | real | [-] | 6.5 – 7.5 | Optimal pH for biological reactions. Rice paddies: 5.5–7.5. |

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
K_Monod_NH4   : 0.5           // [mg NH4-N/L]
K_Monod_DO    : 1.5           // [mg O2/L]

# === Denitrification ===
Rate_Max_Denit : 5.0e-6       // [mg NO3/L/s]
K_Monod_NO3    : 0.5          // [mg NO3-N/L]
K_Monod_DOC    : 5.0          // [mg DOC/L]
Ki_Inhib_DO    : 1.0          // [mg O2/L]

# === Aerobic Respiration ===
Rate_Max_Hetero       : 5.0e-5  // [mg DOC/L/s]
K_Monod_DOC_Aerobic   : 5.0     // [mg DOC/L]
K_Monod_DO_Aerobic    : 1.0     // [mg O2/L]

# === Mineralization ===
Rate_Max_Min   : 5.0e-6        // [mg DON/L/s]
K_Monod_DON    : 0.21015       // [mg DON-N/L]
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
| `solid_mode` | int | **Yes** | Solid-phase initial condition mode: `0` = constant value (`solid_val`), `1` = spatial file (`solid_file`). |
| `solid_val` | real | If `solid_mode=0` | Uniform initial solid-phase (adsorbed) concentration [mg/kg]. |
| `solid_file` | string | If `solid_mode=1` | Filename for spatially distributed solid-phase initial concentration. |
| `aq_phase` | int | No | Alias for `aq_mode` (accepted for backward compatibility). |

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
direction : 6                            // -z face (top of groundwater domain)
polygon : polygontop.input
spec_0_bctype : 8                        // NH4: surface-water exchange
spec_1_bctype : 8                        // NO3: surface-water exchange
spec_2_bctype : 8                        // DO: surface-water exchange
spec_3_bctype : 8                        // DOC: surface-water exchange
spec_4_bctype : 8                        // DON: surface-water exchange

// --- Bottom boundary: free drainage ---
id : rtbottom
direction : 5                            // +z face (bottom of domain)
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

#### Additional Global Parameters

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Bulk_Density` | real | [g/cm³] | 1.1 – 1.4 | Bulk density of the porous medium (ρb). Typical rice plow layer: 1.30 g/cm³. |

#### Sorption Parameters (per-species, semicolon-separated)

| Key | Type | Units | Description |
|-----|------|-------|-------------|
| `Sorption_Type` | int list | [-] | Sorption model for each species: `0` = None, `1` = Linear (uses `Kd`), `2` = Freundlich (uses `Kf`, `Nf`), `3` = Langmuir (uses `Alpha_D`, `Beta_D`). |
| `Kd` | real list | [L/kg or m³/kg] | Distribution coefficient for **Linear** sorption isotherm. Typical NH4 in paddy: 2.0–5.0. Set `0.0` for non-adsorbing species. |
| `Kf` | real list | [(mg/kg)/(mg/L)^Nf] | Freundlich constant. Set `0.0` if not using Freundlich. |
| `Nf` | real list | [-] | Freundlich nonlinearity exponent. Typical: 0.5–1.2. Set `1.0` if not used. |
| `Alpha_D` | real list | [L/mg] | Langmuir equilibrium constant. Set `0.0` if not using Langmuir. |
| `Beta_D` | real list | [mg/kg] | Langmuir sorption capacity. Set `0.0` if not using Langmuir. |
| `Beta` | real list | [1/s] | Non-equilibrium sorption parameter. Set `0.0` to disable non-equilibrium effects. |

#### Sorption Model Summary

| `Sorption_Type` | Model | Required Parameters | Isotherm Equation |
|:---:|-------|---------------------|-------------------|
| 0 | None (no sorption) | — | — |
| 1 | Linear | `Kd` | S = Kd × C |
| 2 | Freundlich | `Kf`, `Nf` | S = Kf × C^Nf |
| 3 | Langmuir | `Alpha_D`, `Beta_D` | S = (Alpha_D × Beta_D × C) / (1 + Alpha_D × C) |

where S = adsorbed concentration [mg/kg], C = dissolved concentration [mg/L].

#### Groundwater-Specific Source/Sink Parameters

| Key | Type | Units | Typical Range | Description |
|-----|------|-------|---------------|-------------|
| `Ka_DO` | real | [1/s] | 0 – 1×10⁻⁵ | Reaeration rate constant for groundwater. Typical: 1×10⁻⁵ (≈ 0.5/3600). |
| `DOC_eq` | real | [mg C/L] | 2 – 20 | DOC release equilibrium concentration from sediment. |
| `DON_eq` | real | [mg N/L] | 0.1 – 5.0 | DON release equilibrium concentration from sediment. |
| `K_rel_DOC` | real | [1/s] | 0 – 1×10⁻⁵ | DOC release rate constant. Set `0.0` to disable. |
| `K_rel_DON` | real | [1/s] | 0 – 1×10⁻⁵ | DON release rate constant. Set `0.0` to disable. |

#### Example

```text
# === Global ===
ReactModel_Type           : 1
Bulk_Density              : 1.30                // [g/cm^3]
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
K_Monod_NH4   : 0.5                    // [mg NH4-N/L]
K_Monod_DO    : 0.5                    // [mg O2/L]

# === Denitrification ===
Rate_Max_Denit : 5.0e-5                // [mg NO3/L/s]
K_Monod_NO3    : 0.5
K_Monod_DOC    : 2.0
Ki_Inhib_DO    : 0.2

# === Aerobic Respiration ===
Rate_Max_Hetero       : 1.0e-4
K_Monod_DOC_Aerobic   : 5.0
K_Monod_DO_Aerobic    : 0.5

# === Mineralization ===
Rate_Max_Min  : 3.0e-6
K_Monod_DON   : 0.21015
K_Monod_DO_Min: 0.48

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
sstype : <type>                  // Source/sink type
ssvals : <value>                 // Constant value (if applicable)
polygon : <polygon_file>         // Spatial extent polygon
timeseries : <timeseries_file>   // Time series data (if applicable)
```

#### Source/Sink Types (`sstype`)

| Value | Description |
|:---:|-------------|
| 0 | Root water uptake zone (coupled with WOFOST crop model). |
| 1 | Volumetric mass flux injection [mg/L/s]. |
| 2 | Concentration-specified point source coupled with water flow [mg/L]. |
| 3 | Area mass flux injection (fertilizer spreading) [mg/m²/s]. Automatically distributed over vertical grid thickness. |

#### Example — Root Zone for WOFOST Coupling

```text
sscount : 1
id : root1
sstype : 0                       // Root uptake zone
ssvals : 0.0
polygon : polygonroot.input      // Root zone polygon
timeseries : et.input            // Evapotranspiration time series
```

---

### 2.5 Initial Condition Files

**Source code**: `src/RTInitSW.h` → `readRtICFileSW()`, `src/RTInitGW.h` → `readRtICFileGW()`

Initial concentration fields can be specified as spatial grids. Both surface water and groundwater modules support two file formats.

#### File Naming Convention

| Domain | Phase | Pattern | Example |
|--------|-------|---------|---------|
| Surface water | Liquid (aqueous) | `rtinitialsw<Name>.input` | `rtinitialswNH4.input` |
| Groundwater | Liquid (aqueous) | `rtinitialgw<Name>.input` | `rtinitialgwNO3.input` |
| Groundwater | Solid (adsorbed) | `rtinitialgwsolid<Name>.input` | `rtinitialgwsolidNH4.input` |

#### Format 1: Simplified ASCII (3-line header)

```text
<ncols>
<nrows>
<nodata_value>
<data values...>     // Row-major, nrows × ncols values
```

#### Format 2: ESRI ASCII Grid (6-line header)

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
- Nodata values are passed through (not masked).
- Units: liquid phase [mg/L or kg/m³], solid phase [mg/kg].

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
| `CROP_NO` | int | Crop identification number. Example: `501` for rice. |

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
CROP_NO : 501

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

Provides daily meteorological forcing data for the WOFOST crop model.

#### Format

Time series with one row per time step. Each row contains:

| Column | Type | Units | Description |
|--------|------|-------|-------------|
| 1 | real | [s] | Time (seconds since simulation start) |
| 2 | real | [m/s] | Rainfall rate |
| 3 | real | [°C] | Minimum daily temperature |
| 4 | real | [°C] | Maximum daily temperature |
| 5 | real | [W/m²] | Shortwave radiation |
| 6 | real | [kPa] | Vapor pressure |
| 7 | real | [m/s] | Wind speed |
| 8+ | real | various | Additional ET/potential evaporation values |

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

## 4. Timeseries File Format

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
| Boundary condition timeseries (`spec_X_bcfile`) | [mg/L or kg/m³] | `rtirrigationNH4.input` — NH4 concentration in irrigation water |
| Source/sink timeseries (`spec_X_ssfile`) | [mg/m²/s] | `rtfertNH4-DON-NH4.input` — area mass flux of NH4 fertilizer |
| Evapotranspiration timeseries | [m/s] | `et.input` — potential ET rate |
| Rainfall timeseries | [m/s] | `rainfall.input` — rainfall rate |

#### Interpolation

- **sstype = 1** (linear): Linear interpolation between consecutive data points.
- **sstype = 2** (step): Piecewise constant (step function). The value holds until the next time point.

#### Example — Fertilizer Application (`rtfertNH4-DON-NH4.input`)

```text
ndata 126
0.00000e+00 0.00000e+00
561600.00    1.31250e+00         // Day 6.5: NH4 pulse starts
838800.00    1.12500e+00         // Day 9.7: NH4 application
1080000.00   0.00000e+00         // Day 12.5: pulse ends
3240000.00   1.68750e+00         // Day 37.5: second application
3499200.00   0.00000e+00         // Day 40.5: ends
...
95040000.0   0.00000e+00         // End of simulation
```

> **Note**: The `ndata` keyword is parsed by reading one token before the count. The parser ignores the label and reads the number that follows.

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

The following table lists all input files related to solute transport and crop growth. Check that all required files are present in the input directory before running a simulation.

### Surface Water Solute Transport

| File | Required | Purpose |
|------|----------|---------|
| `rttransportsw.input` | **Yes** | Species count, transport properties, initial conditions, solver settings |
| `rtswbc.input` | No* | Solute boundary conditions (*required if external boundaries exist) |
| `rtswss.input` | No | Solute source/sink zones (fertilizer, etc.) |
| `reactionsw.input` | If ReactionModule=1 | Reaction kinetics, nitrogen cycle parameters |
| `rtinitialsw*.input` | If aq_mode=1 | Spatial initial concentration grids for each species |

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
| `polygon*.input` | As referenced | Polygon vertex coordinates for BC/SS spatial zones |
| `drains.input` | If tile drains used | Tile drainage parameters |
| Timeseries files (e.g., `rtfert*.input`) | As referenced | Time-varying boundary/source data |

---

*End of User Manual — Solute Transport & WOFOST Crop Growth*
