# Case 4: IDEAL pump + time-varying river + flap gate

## Scenario

Manhole J0 receives surface runoff; the IDEAL pump P0 lifts it to manhole J1, and gravity pipe C0 carries it to outfall O2. The outfall uses a time-varying river stage (TIMESERIES) with a flap gate (`gated=yes`) that blocks backflow at high tide.

```text
rain -> surface (2D) -> manhole J0 -> IDEAL pump P0 -> manhole J1 -> pipe C0 -> outfall O2 (TIMESERIES, flap gate)
```

See `schematic.png`.

## Key points

- In `[PUMPS]`, setting `curveId` to `*` (or `IDEAL`, `-1`) selects an ideal pump; no `[CURVES]` is needed.
- An ideal pump ignores head: flow = current upstream-node inflow + overflow (minus competing outflow). It "pumps whatever arrives", not a prescribed constant flow.
- `gated=yes` installs a flap gate at the outfall that blocks backflow when the river is higher than the network.

## Run

```bash
serghei this_dir/ this_dir/output/ 2
```

## Results

- `PUMP_LINK_*_flow`: pump flow grows from 0 to ~0.66 m³/s as the upstream inflow grows.
- `PUMP_LINK_*_setting`: the pump turns on (0 → 1) once manhole J0 depth exceeds 0.8 m.
- `OutfallBackflowDischarge`: stays 0 throughout, showing the flap gate blocked high-tide backflow.

Change the outfall `yes` to `no` in `Drainage.inp` and re-run to observe the river backflowing through C0 into J1 at high tide.
