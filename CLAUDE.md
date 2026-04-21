# CLAUDE.md

FESOM2 Dwarf Tracer — a standalone advection-transport miniapp extracted from FESOM2.
Fortran 2008 with CMake build system. Builds a shared library (`libfesom_tracer`) and three driver executables.

## Build & verify

```bash
./configure.sh --compiler gnu --precision dp --clean --build
cd build_gnu_dp && ./run.sh 1 20 20 10 --periodic
```

Expected output: sum conserved at `0.3249000000E+05`, 0 NaN at every step.

Supported compilers: `gnu`, `intel`, `nvidia`. Precisions: `dp`, `sp`, `hp` (HP is NVIDIA-only).

## Directory layout

- `src/` — executable entry points (3 drivers)
- `lib/` — shared library Fortran sources
- `cmake/` — CMake modules
- `python/` — Python demo/test scripts
- `visualization/` — plotting scripts
- `docs/` — detailed docs (see [Further reading](#further-reading))

## Key files

| Purpose | File |
|---|---|
| Precision parameters (WP, MP) | `lib/oce_modules.F90` — see [PRECISION.md](docs/PRECISION.md) |
| Compiler flags | `CMakeLists.txt` → `apply_fesom_compile_flags()` |
| Advection entry point | `lib/oce_adv_tra_driver.F90` → `do_oce_adv_tra()` |
| Mesh type (T_MESH) | `lib/MOD_MESH.F90` |
| Tracer types | `lib/MOD_TRACER.F90` (T_TRACER_WORK=MP, T_TRACER_DATA=WP) |
| Analytic mesh generation | `lib/analytic_mesh.F90` |
| FP16 math wrappers | `lib/hp_math_intrinsics.F90` |
| C interface (optional) | `lib/tracer_c_interface.F90` — see [INTERFACES.md](docs/INTERFACES.md) |

## Coding conventions

- All reals use `real(kind=WP)` or `real(kind=MP)`, never hardcoded `real*8`.
- Literal constants use `_WP` suffix (e.g., `1.0e-3_WP`).
- Mesh/geometry arrays must be `real(kind=MP)` — `MP = max(WP, 4)`, at least single precision.
- Tracer data arrays use `real(kind=WP)` — the precision under test.
- Diagnostic sums wrap with `dble()` to avoid FP16 accumulation overflow.
- MPI datatype: use `MPI_WP`, not hardcoded `MPI_DOUBLE_PRECISION`.

## Precision pitfalls

- FP16 locals in mesh computation overflow (`b * r_earth` exceeds 65504) — always use MP for mesh locals.
- GNU and Intel do not support `real(kind=2)` — HP builds are NVIDIA-only.
- NVIDIA HP: intrinsic math functions lack FP16 overloads → `hp_math_intrinsics.F90` provides wrappers.
- Halo exchange is compiled out for HP (`#if !defined(USE_HALF_PRECISION)`).

## OpenMP (vertex-parallel)

Build with `--openmp` flag:
```bash
./configure.sh --compiler gnu --precision dp --openmp --clean --build
cd build_gnu_dp && OMP_NUM_THREADS=4 ./run.sh 1 200 200 10 --periodic
```

**Critical: MPI process binding.** OpenMPI defaults to `--bind-to core`, which pins all
OpenMP threads to a single core. The generated `run.sh` uses `--bind-to none` to fix this.
If running `mpirun` directly, always add `--bind-to none` (or `--bind-to socket`).

Edge-to-node scatter uses vertex-gather (each node gathers from incident edges).
See [VERTEX_PARALLEL.md](docs/VERTEX_PARALLEL.md) for details.

Key files:
- `lib/oce_node_edge_map.F90` — builds node-to-edge adjacency at init; also
  pre-packs static edge geometry into node-contiguous arrays for
  memory-bandwidth optimization
- Must call `build_node_edge_map()` and `build_node_edge_packed()` after mesh
  setup in all drivers (pre-packed arrays eliminate ~3 levels of indirect
  indexing in the vertex-gather hot loops)

## Compiler-specific notes

- **Intel**: uses `mpiifx`/`mpiicx`; source `setvars.sh` before building.
- **NVIDIA**: must use bundled MPI at `$NVHPC_ROOT/comm_libs/mpi/bin/mpifort`.
- Don't copy binaries between build dirs — rpath uses `$ORIGIN/../lib`.

## Further reading

- [ARCHITECTURE.md](docs/ARCHITECTURE.md) — directory layout, source file descriptions, build products
- [PRECISION.md](docs/PRECISION.md) — mixed precision design (WP/MP), compiler support, FP16 details
- [INTERFACES.md](docs/INTERFACES.md) — C-compatible API, Python wrapper, array layouts
- [VISUALIZATION.md](docs/VISUALIZATION.md) — plotting scripts setup, usage, multi-precision comparison
- [CHANGES.md](docs/CHANGES.md) — changelog

## Documentation maintenance

After implementing a feature that changes behavior, CLI options, file formats, or conventions, update the relevant doc in `docs/`. If no existing doc covers the new feature, ask the user whether a new doc should be created in `docs/` before proceeding. Verify documentation integrity periodically — if code and docs have diverged, fix the docs before moving on.
