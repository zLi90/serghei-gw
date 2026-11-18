# Case 2: Wet well + TYPE3 pump + fixed outfall stage

## Scenario

Runoff enters manhole J1, flows by gravity through pipe C0 into a buried wet well (storage S0), and pump P0 lifts it to a fixed-stage (FIXED) outfall O2. The pump is a TYPE3 head–flow pump: its flow is interpolated from an H–Q curve using the instantaneous head, and it is switched on/off by the wet-well depth (`yOn=0.4`, `yOff=0.15`).

```text
rain -> surface (2D)
        | manhole J1 (sealed=0)
       C0 -> wet well S0 (25 m2, invert 5.0 m, depth 3 m)
              | pump P0 (TYPE3)
              | outfall O2 (FIXED 6.5 m)
```

See `schematic.png`.

## Key points

- The wet well is defined in `[STORAGE]`: `FUNCTIONAL 25.0 0.0 0.0` is a prismatic tank of constant 25 m² area and 75 m³ capacity.
- A pump is a link, `from` the wet well `to` the outfall; it is not a node, and the storage unit is not attached automatically.
- Startup/shutoff use only the upstream-node depth: `>0.4 m` starts, `<0.15 m` stops.
- Head = downstream absolute stage − upstream absolute stage; the downstream FIXED stage is constant 6.5 m.

## Run

```bash
serghei this_dir/ this_dir/output/ 2
```

## Results

In `output/DrainageTimeSeries.out`, the `PUMP_LINK_*` columns report:

- `flow`: actual pumped flow, ~0.49 m³/s at steady state here.
- `head`: ~0.19 m here.
- `setting`: current switch (1 = on).
- `limited`: whether flow was reduced by water availability (0 = no).
- `volume / runTime / startUps`: cumulative pumped volume, operating time, startups (~805 m³, ~1830 s, 24 times here).
- `hydraulicEnergy`: hydraulic energy (kWh), excluding pump/motor efficiency.

When the wet-well inflow exceeds the pump curve capacity, the wet-well level rises (~1.3 m here) and the pump runs continuously; early on, the level oscillates between 0.15–0.4 m and the pump cycles.
