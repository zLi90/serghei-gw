# Case 1: Gravity drainage (no pump)

## Scenario

Rainfall forms runoff on the 2D surface, enters the network through manhole J0, flows through two gravity pipes C0 and C1, and leaves at outfall O2 as a free outfall (FREE). There is no pump and no storage unit; this case introduces the basic surface → network → outfall coupling.

See `schematic.png`:

```text
rain -> surface (2D)
        | manhole J0 (sealed=0, exchanges with surface)
       C0 (gravity pipe, 0.6 m)
        | manhole J1 (sealed=1, 1D only)
       C1 (gravity pipe, 0.6 m)
        | outfall O2 (FREE)
```

## Input files

All files live in this directory; pass this directory to `serghei`.

- `parameters.input`: global time, CFL, output frequency, decomposition.
- `dem.input`: uniform 30×10 grid, 10 m cells, ground elevation 9.0 m.
- `sw.input`: surface shallow-water initial values.
- `rainfall.input`: constant 500 mm/h rain driving runoff.
- `Drainage-parameter.input`: drainage and 1D–2D exchange parameters.
- `Drainage.inp`: network topology, every field commented with `;`.

## Run

```bash
serghei this_dir/ this_dir/output/ 2
```

The `output/` directory must not already exist.

## Results

Open `output/DrainageTimeSeries.out`. Key columns:

- `JUNCTION_0_2D-1D[L^3/s]`: surface inflow into manhole J0.
- `NODE_*_depth[m]`, `LINK_*_depth[m]`: node/pipe depths.
- `OutfallInflow[L^3]`: cumulative volume discharged at the outfall (~490 m³ over 3600 s in this case).
- `OutfallDischarge[L^3/s]`: current outfall flow.

The network starts dry, fills as rain accumulates, and the pipes transition from dry to flowing, matching gravity drainage behavior.
