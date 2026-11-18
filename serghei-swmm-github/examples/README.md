# Surface–Drainage Coupling Examples

This directory contains four self-contained demonstration cases. They all run the surface + drainage coupling with the subsurface module disabled, and they share the same surface grid and rainfall; only `Drainage.inp` differs. Each case directory contains every input file required to run it, a `schematic.png` structure diagram, and a `README.md`.

| Case | Scenario | Pump | Outfall boundary |
|------|----------|------|------------------|
| [case01_no_pump](case01_no_pump/) | gravity drainage only | none | FREE free outfall |
| [case02_type3_fixed](case02_type3_fixed/) | wet-well lift station | TYPE3 (head–flow curve) | FIXED stage |
| [case03_type3_timeseries](case03_type3_timeseries/) | wet-well station under backwater | TYPE3 | TIMESERIES river stage |
| [case04_ideal_timeseries](case04_ideal_timeseries/) | ideal pump + flap gate | IDEAL (pump whatever arrives) | TIMESERIES + flap gate |

## How to run

Using case 1 as an example (in WSL or Linux):

```bash
# After building, the serghei executable is at build/bin/serghei
SERGHEI_BIN=/path/to/serghei-swmm/build/bin/serghei
"$SERGHEI_BIN" /path/to/case01_no_pump/ /path/to/case01_no_pump/output/ 2
```

The three arguments are: input directory (trailing `/`), output directory (trailing `/`, must not already exist), and the number of OpenMP threads.

Run all four cases in one go:

```bash
cd examples
SERGHEI_BIN=/path/to/serghei-swmm/build/bin/serghei
for d in case01_no_pump case02_type3_fixed case03_type3_timeseries case04_ideal_timeseries; do
  rm -rf "$d/output"
  "$SERGHEI_BIN" "$d/" "$d/output/" 2 > "$d/run.log" 2>&1
done
```

> The screen prints `the swe dt and drainage dt are: ...` on every step (a leftover debug print); it is redirected to `run.log` above. A successful run ends with `[OK] SIMULATION FINISHED`.

## Files in each case

- `parameters.input`: global time, CFL, output frequency, MPI decomposition (`//` comments).
- `dem.input`: surface elevation raster (ESRI ASCII).
- `sw.input`: surface shallow-water initial conditions and roughness.
- `rainfall.input`: rainfall time series.
- `Drainage-parameter.input`: drainage and 1D–2D exchange parameters (`//` comments).
- `Drainage.inp`: network topology and boundaries (`;` comments explaining every section).
- `README.md`: case description.
- `schematic.png`: physical structure diagram.

For the full meaning of every input field, see the repository-root [Surface–Drainage Coupling User Guide](../Surface-Drainage-Coupling-User-Guide.md).

## Regenerate the schematics

```bash
cd examples
python3 make_schematics.py
```
