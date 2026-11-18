# SERGHEI-SWMM

SERGHEI-SWMM couples the SERGHEI two-dimensional shallow-water model with a one-dimensional drainage pipe network (dynamic-wave routing in the style of SWMM). Rainfall runs off the 2D surface and enters the network through street inlets and manholes; pipes, storage units (wet wells), pumps (TYPE3 / IDEAL), and free/fixed/time-series outfalls route the flow to a receiving water body.

## Documentation

- [Surface–Drainage Coupling User Guide](Surface-Drainage-Coupling-User-Guide.md): which input files are needed, what they contain, and what each drainage parameter means.
- [examples/README.md](examples/README.md): four runnable demonstration cases.

## Build

Drainage is disabled by default. To build the surface + drainage coupled model:

```bash
cmake -S . -B build \
  -DSERGHEI_SWE_MODEL=ON \
  -DSERGHEI_SUBSURFACE_MODEL=OFF \
  -DSERGHEI_DRAINAGE_MODEL=ON
cmake --build build -j 4
```

Then run an example:

```bash
build/bin/serghei examples/case02_type3_fixed/ out/ 2
```

Note: the output directory must not already exist.
