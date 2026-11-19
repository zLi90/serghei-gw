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
junctionSize : 0.6  // manhole diameter (m); sets the exchange footprint, orifice area and weir perimeter
Cw : 0.8            // weir discharge coefficient (surface water flowing freely into a node)
Co : 0.6            // orifice discharge coefficient (submerged/pressure exchange)
outputLinkDepth : 1 // also write LINK_*_depth[m] to the time-series output (0/1)
writeNodeVolume : 0 // also write DrainageNodeVolume.out (per-node exchange accumulation, 1/0)
```

| Parameter | Unit | Meaning |
|-----------|------|---------|
| `inletLength` | m | length of an inlet along the flow direction |
| `inletWidth` | m | width of an inlet |
| `junctionSize` | m | manhole diameter; sets the exchange footprint, orifice area (junctionSize²) and weir perimeter (4×junctionSize) |
| `Cw` | — | weir discharge coefficient, used for free surface inflow |
| `Co` | — | orifice discharge coefficient, used for submerged/pressure exchange |
| `outputLinkDepth` | — | when `1`, extra `LINK_*_depth[m]` columns are written |
| `writeNodeVolume` | — | when `1`, also write `DrainageNodeVolume.out` (per-node exchange accumulation) |

The orifice area and weir perimeter are now computed per node type instead of using a single parameter: an inlet (INLET) uses `inletLength×inletWidth` and `2(inletLength+inletWidth)`, while a manhole (JUNCTION) uses `junctionSize²` and `4×junctionSize`.

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

## 5. Output files and results

After a run, the output directory (`output_dir/`) contains:

| File | Contents |
|------|----------|
| `DrainageTimeSeries.out` | drainage time series (main result, one row per `obsFreq` seconds) |
| `DrainageNodeVolume.out` | per-node exchange accumulation (requires `writeNodeVolume : 1`) |
| `domainTimeSeries.out` | surface-module time series |
| `output.nc` (or multiple NetCDF files written every `outFreq`) | surface 2D field output |
| `log.out` | run log |

`DrainageTimeSeries.out` is a whitespace-separated text file with the following columns in order:

| Column | Unit | Meaning |
|--------|------|---------|
| `Time` | s | simulation time |
| `two2oneAccum` | L³ | cumulative surface-to-drainage exchange volume over the whole run (sum over all exchanging nodes) |
| `one2twoAccum` | L³ | cumulative drainage-to-surface overflow volume over the whole run (sum over all nodes) |
| `OutfallInflow` | L³ | cumulative inflow summed over all outfall nodes (water arriving at the outfalls) |
| `OutfallBackflow` | L³ | cumulative backflow summed over all outfall nodes (high boundary flowing back into the network) |
| `OutfallDischarge` | L³/s | current discharge summed over all outfall nodes (instantaneous) |
| `OutfallBackflowDischarge` | L³/s | current backflow summed over all outfall nodes (instantaneous) |
| `JUNCTION_i_2D-1D` / `JUNCTION_i_1D-2D` | L³/s | current exchange between manhole i and the surface (into / out of the network) |
| `INLET_i_2D-1D` / `INLET_i_1D-2D` | L³/s | current exchange between inlet i and the surface |
| `OUTFALL_i` | L³/s | current inflow at outfall node i (single node, not summed) |
| `NODE_i_depth` | m | current depth at node i |
| `LINK_i_depth` | m | current depth in link i (requires `outputLinkDepth : 1`) |
| `PUMP_LINK_j_flow / head / setting / limited / offCurve` | see below | pump flow, head, switch, limited flag, off-curve flag |
| `PUMP_LINK_j_volume / runTime / startUps` | m³ / s / count | cumulative pumped volume, operating time, startups |
| `PUMP_LINK_j_hydraulicEnergy` | kWh | cumulative hydraulic energy (excluding pump/motor efficiency) |
| `PUMP_LINK_j_offCurveLow / offCurveHigh` | s | cumulative operating time below/above the curve |

Notes:

- `two2oneAccum`, `one2twoAccum`, `OutfallInflow`, and `OutfallBackflow` are **cumulative over the whole run**; columns whose names contain `Discharge` or `_2D-1D` / `_1D-2D` are **instantaneous flow rates** (L³/s).
- `OUTFALL_i` is the instantaneous inflow of a single outfall; `OutfallDischarge` is the instantaneous sum over all outfalls.

`DrainageNodeVolume.out` (written when `writeNodeVolume : 1`) records, for each output time, one line per node: `Time` (s), `Node` (node index), `two2oneAccum` (cumulative surface-to-network volume at this node, m³), and `one2twoAccum` (cumulative network-to-surface volume at this node, m³). These two columns are **per-node** accumulations, unlike the global totals in `DrainageTimeSeries.out`, and are useful for checking each node's water balance.

## 6. Examples

`examples/` contains four runnable cases covering "with/without pump" and "fixed/variable outfall stage" combinations. Each case directory has all input files (with commented `Drainage.inp`), a `schematic.png` structure diagram, and a `README.md`:

| Case | Scenario |
|------|----------|
| [case01_no_pump](examples/case01_no_pump/) | gravity drainage, no pump, FREE outfall |
| [case02_type3_fixed](examples/case02_type3_fixed/) | wet well + TYPE3 pump, FIXED outfall stage |
| [case03_type3_timeseries](examples/case03_type3_timeseries/) | wet well + TYPE3 pump, TIMESERIES river stage |
| [case04_ideal_timeseries](examples/case04_ideal_timeseries/) | IDEAL pump + TIMESERIES river stage + flap gate |

See [examples/README.md](examples/README.md) for details.
