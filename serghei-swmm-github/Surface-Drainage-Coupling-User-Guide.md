# SERGHEI-SWMM Surface–Drainage Coupling User Guide

This guide explains which input files are needed for different scenarios, what each file contains, and what the parameters mean. The drainage inputs are explained in detail; the surface and subsurface inputs follow the existing SERGHEI conventions and are only summarized. The subsurface module is disabled by default in the examples.

## 1. Required input files

To run a coupled surface–drainage case, the case directory needs the following files:

| File | Required | Purpose |
|------|----------|---------|
| `parameters.input` | yes | global time, CFL, output frequency, process decomposition |
| `dem.input` | yes | surface elevation raster |
| `sw.input` | yes | surface shallow-water initial conditions and roughness |
| `rainfall.input` | yes | rainfall time series |
| `Drainage-parameter.input` | yes (defaults used if missing) | drainage and 1D–2D exchange parameters |
| `Drainage.inp` | yes | pipe network topology, pumps, outfall boundaries |

The surface module also has optional files (`extbc.input`, `infiltration.input`, `wind.input`, `Ori_dem.input`, ...), which fall back to defaults when absent, exactly as in a surface-only case.

The differences between scenarios are almost entirely in which sections you write in `Drainage.inp`; see Section 3.

## 2. Input files and parameters

### 2.1 parameters.input (global, surface module)

`key : value` format; `//` starts a comment.

```text
simLength : 3600   // total simulation time (s)
cfl       : 0.4    // surface CFL number
outFreq   : 300    // surface NetCDF field output interval (s)
obsFreq   : 10     // time-series output interval (s); also the write interval of DrainageTimeSeries.out
parNx     : 1      // MPI processes in x
parNy     : 1      // MPI processes in y
outFormat : NETCDF // surface output format
BCtype    : REFLECTIVE  // 2D outer boundary type
```

### 2.2 Surface files (dem / sw / rainfall)

- `dem.input`: ESRI ASCII grid giving the ground elevation of every cell.
- `sw.input`: surface initial depth, roughness, dry/wet thresholds, same as a surface-only case.
- `rainfall.input`: rainfall time series (time unit s, intensity unit mm/h).

Readers familiar with the surface module already know these files, so they are not expanded here.

### 2.3 Drainage-parameter.input (drainage + 1D–2D exchange parameters)

`key : value` format; `//` starts a comment. If the file is missing, built-in defaults are used and the run is unaffected.

```text
inletLength : 0.5   // inlet length along the flow direction (m)
inletWidth  : 0.3   // inlet width (m)
junctionSize : 0.6  // manhole diameter (m), approximated as a square footprint in the exchange
Cw : 0.8            // weir discharge coefficient (surface water flowing freely into a node)
Co : 0.6            // orifice discharge coefficient (submerged/pressure exchange)
w : 1.0             // weir width weighting factor (default 1.0)
outputLinkDepth : 1 // also write LINK_*_depth[m] to the time-series output (0/1)
```

| Parameter | Unit | Meaning |
|-----------|------|---------|
| `inletLength` | m | length of an inlet along the flow direction |
| `inletWidth` | m | width of an inlet |
| `junctionSize` | m | manhole diameter, used to determine the footprint that exchanges with the surface |
| `Cw` | — | weir discharge coefficient, used for free surface inflow |
| `Co` | — | orifice discharge coefficient, used for submerged/pressure exchange |
| `w` | — | weir width weighting factor, default 1.0 |
| `outputLinkDepth` | — | when `1`, extra `LINK_*_depth[m]` columns are written |

### 2.4 Drainage.inp (network topology and boundaries)

Plain text; anything after `;` is a comment.

#### [JUNCTIONS] manholes / outfall nodes

```text
;; id  invert_elev  max_depth  init_depth  sur_depth  ponded_area  X  Y  type  sealed
0  7.0  2.0  0.0  0.0  300.0  50.0  0.0  junction  0
```

| Column | Unit | Meaning |
|--------|------|---------|
| `id` | — | node number, must be numbered 0,1,2,... in file reading order |
| `invert_elev` | m | absolute elevation of the node bottom (invert) |
| `max_depth` | m | full depth; rim elevation = invert + max depth |
| `init_depth` | m | 0 means the model derives it from the connected links |
| `sur_depth` | m | surcharge depth above max depth |
| `ponded_area` | m² | surface area available when flooding |
| `X` `Y` | m | node coordinates; they select the surface cell for exchange and must lie inside the domain |
| `type` | — | `junction` (manhole) or `outlet` (outfall node) |
| `sealed` | — | `1` disables surface inflow into this node, `0` enables it |

#### [INLETS] street inlets

```text
;; id  invert_elev  max_depth  init_depth  sur_depth  ponded_area  X  Y  connected_manhole
0  7.0  1.0  0.0  0.0  1.0  50.0  0.0  0
```

The inlet geometry comes from `inletLength` / `inletWidth` in `Drainage-parameter.input`; the last column is the manhole it drains into.

#### [STORAGE] storage unit (wet well)

```text
;; id  invert_elev  max_depth  init_depth  shape  p1  p2  p3  sur_depth  ponded_area  X  Y
0  5.0  3.0  0.0  FUNCTIONAL  25.0  0.0  0.0  0.0  0.0  150.0  0.0
```

| Column | Unit | Meaning |
|--------|------|---------|
| `shape` | — | `FUNCTIONAL` or `CYLINDRICAL` |
| `p1 p2 p3` | — | `FUNCTIONAL`: area `A(d)=p3+p1*d^p2`; `CYLINDRICAL`: `p1`, `p2` are the axes |

The example `p1=25, p2=0, p3=0` is a prismatic wet well with constant 25 m² area, 3 m deep, 75 m³ capacity.

#### [CONDUITS] pipes

```text
;; id  from  to  length  roughness  in_offset  out_offset  init_flow  max_flow  diameter
0  1  0  80.0  0.013  0.0  0.0  0.0  0.0  0.50
```

| Column | Unit | Meaning |
|--------|------|---------|
| `from` `to` | — | upstream/downstream node numbers; the `from` end invert should be higher than the `to` end |
| `length` | m | pipe length, given by input, not computed from coordinates |
| `roughness` | — | Manning roughness coefficient |
| `in_offset / out_offset` | m | elevation offset of each pipe end above the node invert |
| `init_flow` | m³/s | initial pipe flow; 0 derives it from node depths |
| `max_flow` | m³/s | flow limit |
| `diameter` | m | circular pipe diameter |

#### [PUMPS] pumps and [CURVES] pump curves

```text
[PUMPS]
;; id  from  to  curveId  init_status  startup_depth  shutoff_depth
0  0  2  0  0  0.4  0.15

[CURVES]
;; curveId  type  head_H(m)  flow_Q(m3/s)
0  pump3  0.0  0.50
0  pump3  1.5  0.40
0  pump3  3.0  0.30
0  pump3  5.0  0.000
```

[PUMPS] fields:

| Column | Unit | Meaning |
|--------|------|---------|
| `from` `to` | — | upstream/downstream node numbers; a pump is a link, not a node |
| `curveId` | — | an integer is a TYPE3 pump (references `[CURVES]`); `*` / `IDEAL` / `-1` is an ideal pump |
| `init_status` | — | `ON`/`1` on, `OFF`/`0` off |
| `startup_depth` | m | pump starts when the upstream node depth `> startup_depth`; 0 disables it |
| `shutoff_depth` | m | pump stops when the upstream node depth `< shutoff_depth`; 0 disables it |

[CURVES] fields:

| Column | Unit | Meaning |
|--------|------|---------|
| `curveId` | — | curve number, referenced by the pump `curveId` |
| `type` | — | always `pump3` / `PUMP3` |
| `head_H` | m | head = downstream absolute stage − upstream absolute stage |
| `flow_Q` | m³/s | flow at that head; linear interpolation between points |

Notes: a TYPE3 pump interpolates flow from the curve using the head; to stop pumping above some head, set the last point's flow to 0. An ideal pump ignores head, pumps the upstream-node inflow (plus overflow), and needs no `[CURVES]`.

#### [OUTFALLS] outfall boundary and [RIVERSTAGE] river-stage series

```text
[OUTFALLS]
;; nodeId  mode  param  gated
2  TIMESERIES  0  yes

[RIVERSTAGE]
;; seriesId  time(s)  absolute_stage(m)
0  0     6.5
0  1800  8.5
0  3600  6.5
```

[OUTFALLS] fields:

| Column | Meaning |
|--------|---------|
| `nodeId` | outfall node number (an `outlet` node in `[JUNCTIONS]`) |
| `mode` | `FREE` free outfall; `FIXED` fixed stage; `TIMESERIES` time-varying stage |
| `param` | unused for `FREE`; absolute stage for `FIXED`; `[RIVERSTAGE]` series number for `TIMESERIES` |
| `gated` | `yes`/`1` installs a flap gate that blocks backflow, `no`/`0` allows backflow |

[RIVERSTAGE] fields: `seriesId`, `time(s)`, `absolute_stage(m)`, linearly interpolated at the current time.

## 3. Which sections each scenario needs

Only write the relevant sections in `Drainage.inp`:

| Scenario | Drainage.inp sections needed |
|----------|------------------------------|
| Gravity drainage (no pump, no storage) | `[JUNCTIONS]`, `[CONDUITS]`, `[OUTFALLS]` |
| Lift station with a wet well | add `[STORAGE]`, `[PUMPS]`, and `[CURVES]` (for TYPE3) |
| Fixed downstream stage | `[OUTFALLS]` with `FIXED <stage>` |
| Time-varying river stage | `[OUTFALLS]` with `TIMESERIES <series>` plus `[RIVERSTAGE]` |
| Street-inlet inflow | add `[INLETS]` |
| Ideal pump | `[PUMPS]` with `curveId` set to `*`; no `[CURVES]` needed |

To put a wet well in front of a pump, model it yourself: define a `[STORAGE]` node, let the upstream pipes drain into it, and set the pump `from` to that node. The model does not attach a storage unit to a pump automatically.

## 4. Build and run

Enable the drainage module at CMake time (default OFF) together with the surface module:

```bash
cmake -S . -B build \
  -DSERGHEI_SWE_MODEL=ON \
  -DSERGHEI_SUBSURFACE_MODEL=OFF \
  -DSERGHEI_DRAINAGE_MODEL=ON
cmake --build build -j 4
```

Run:

```bash
build/bin/serghei input_dir/ output_dir/ n_threads
```

- The output directory must not already exist; the model creates it. Delete the old one before re-running.
- The screen prints `the swe dt and drainage dt are: ...` every step, a leftover debug print in the source; redirect it to `run.log`.
- The main results are in `output_dir/DrainageTimeSeries.out`; surface fields are written as NetCDF.

## 5. Examples

`examples/` contains four runnable cases covering "with/without pump" and "fixed/variable outfall stage" combinations. Each case directory has all input files (with commented `Drainage.inp`), a `schematic.png` structure diagram, and a `README.md`:

| Case | Scenario |
|------|----------|
| [case01_no_pump](examples/case01_no_pump/) | gravity drainage, no pump, FREE outfall |
| [case02_type3_fixed](examples/case02_type3_fixed/) | wet well + TYPE3 pump, FIXED outfall stage |
| [case03_type3_timeseries](examples/case03_type3_timeseries/) | wet well + TYPE3 pump, TIMESERIES river stage |
| [case04_ideal_timeseries](examples/case04_ideal_timeseries/) | IDEAL pump + TIMESERIES river stage + flap gate |

See [examples/README.md](examples/README.md) for details.
