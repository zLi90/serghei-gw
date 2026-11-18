# Case 3: Wet well + TYPE3 pump + time-varying river stage (backwater)

## Scenario

Same as case 2, but the outfall is changed from a fixed stage to a time series (TIMESERIES): the river rises from 6.5 m to 8.5 m over 1200 s, then falls back to 6.5 m. The higher river raises the pump head, so the TYPE3 curve delivers less flow.

```text
rain -> surface (2D) -> manhole J1 -> C0 -> wet well S0 -> pump P0 -> outfall O2 (TIMESERIES river)
```

See `schematic.png`.

## Key points

- `[OUTFALLS] 2 TIMESERIES 0 no` references series 0 in `[RIVERSTAGE]`; `no` means no flap gate.
- The `[RIVERSTAGE]` series is linearly interpolated at the current time.
- The pump curve ends at `(2.5, 0)`: the higher the head, the lower the flow, reaching zero at 2.5 m.
- Flow drops sharply at high tide; raising the peak stage further would stop the pump entirely.

## Run

```bash
serghei this_dir/ this_dir/output/ 2
```

## Results

Look at `PUMP_LINK_*_flow` and `NODE_2_depth` (outfall depth):

- Low tide (river 6.5 m, outfall depth 0): pump flow ~0.48 m³/s.
- High tide (river 8.5 m, outfall depth 2.0 m): pump flow drops to ~0.05 m³/s.
- Flow recovers as the river falls. This is the backwater effect of a variable downstream boundary on the pump.
