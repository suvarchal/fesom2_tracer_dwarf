# Changelog

## 2026-04-21 — Pre-packed node-edge geometry (memory-bandwidth optimization)

### Motivation
Vertex-parallel OpenMP scaling saturated at ~2.5-3x on 500x500x3 regardless of
thread count. Benchmarks showed the kernel was memory-bandwidth-bound: the
vertex-gather loops chased `node_edge_idx → edge → mesh%edge_tri → nlevels`
through three levels of indirection per edge, hitting random cache lines.

### Changes
- `lib/oce_node_edge_map.F90`: added pre-packed arrays (`node_edge_dxdy`,
  `node_edge_elems`, `node_edge_enodes`, `node_edge_nl1/nu1/nl2/nu2`) and
  `build_node_edge_packed()` subroutine that copies static mesh geometry into
  node-contiguous layout at init time.
- `lib/oce_adv_tra_driver.F90`: converted three vertex-gather loops (FCT LO
  scatter, flux-to-tracer scatter, fused UPW1 kernel) to read pre-packed
  arrays instead of chasing `node_edge_idx → mesh%*(:, edge)`.
- `lib/oce_adv_tra_fct.F90`: converted FCT antidiffusive split loop.
- `src/fesom_analytic.F90`: call `build_node_edge_packed()` after
  `build_node_edge_map()`.

### Results (GNU gfortran `-O3`, 500x500x3, single MPI rank)

| Threads | Before (s) | After (s) | Absolute speedup |
|---------|-----------|-----------|------------------|
| 1       | 15.48     | 0.69      | 22.4x            |
| 4       |  6.42     | 0.30      | 21.3x            |
| 16      |  5.67     | 0.30      | 19.0x            |

The advection loop becomes ~20x faster in absolute terms. See
[PERFORMANCE.md](PERFORMANCE.md#pre-packed-node-edge-geometry-2026-04).

### Verification
- GNU DP: sum conserved at `0.3249000000E+05`, 0 NaN at every step, both
  with 1 and 4 threads.

---

## 2026-02-14 — Repository restructuring and dead code removal

### Directory renames
- `src/` renamed to `lib/` (shared library sources)
- `dwarf_ini/` renamed to `src/` (executable entry points)
- Python scripts moved from root to `python/`
- `tracer_c_interface.F90` and `tracer_init_from_mesh.F90` moved from root to `lib/`
- All documentation markdown files collected into `docs/`
- CMakeLists.txt updated: `src_home` variable renamed to `lib_home`, all paths adjusted

### Dead code removal
- Deleted `lib/MOD_ICE.F90` (849 lines) — full sea-ice module never used by tracer advection
- Deleted `diagnostics_stub.F90` (root and `lib/` copies) — provided `ldiag_DVD = .false.` making all DVD diagnostic branches dead
- Removed `use diagnostics` imports and both `ldiag_DVD` guarded blocks (~80 lines) from `lib/oce_adv_tra_driver.F90`
- Removed `use MOD_ICE`, `type(t_ice)` optional parameter, and ice read/write blocks from `lib/io_restart_derivedtype.F90` (interface + both implementations)
- Removed `USE MOD_ICE` and `type(t_ice)` declarations from `src/fesom.F90` and `lib/tracer_c_interface.F90`
- Removed unused `ice_rejected_salt` allocatable from `lib/oce_modules.F90`
- Removed `diagnostics_stub.F90` and `MOD_ICE.F90` from CMakeLists.txt build sources

### Verification
- GNU DP: builds clean, 10-step analytic run produces identical results (sum perfectly conserved, 0 NaN)

---

## 2026-01 — Mixed precision (MP) implementation

### Design
Introduced a second precision parameter `MP = max(WP, 4)` in `oce_modules.F90`. Mesh and tracer work arrays use `real(kind=MP)` (at least single precision) while tracer data arrays use `real(kind=WP)` (the precision under test). When `WP >= 4` (SP or DP builds), `MP = WP` and all MP changes are no-ops — full backward compatibility.

### Changes
- `lib/oce_modules.F90`: added `MP` parameter; physical constants `pi`, `r_earth`, `omega` changed to `real(kind=MP)`
- `lib/MOD_MESH.F90` (`T_MESH`): all coordinate, area, volume, and geometry arrays changed to `real(kind=MP)`
- `lib/MOD_TRACER.F90`: `T_TRACER_WORK` arrays (fluxes, FCT arrays, tendencies) changed to `real(kind=MP)`; `T_TRACER_DATA` remains `real(kind=WP)`
- `lib/oce_mesh.F90`: local variables in `mesh_areas`, `mesh_auxiliary_arrays`, `edge_center`, `elem_center` changed to MP to prevent FP16 overflow in `b * r_earth`
- `lib/gen_modules_rotate_grid.F90`: added `trim_cyclic_mp` subroutine for MP-kind cyclic trimming
- `lib/gen_halo_exchange.F90`: `exchange_nod`/`exchange_elem` calls on MP arrays guarded with `#if !defined(USE_HALF_PRECISION)` (no-op for npes=1)
- Diagnostic sums wrapped with `dble()` to avoid FP16 accumulation overflow

### Verification
- GNU DP/SP: sum perfectly conserved, 0 NaN
- NVIDIA DP/SP: sum perfectly conserved, 0 NaN
- NVIDIA HP: all values finite, 0 NaN, sum drifts slightly (expected at FP16 precision)

---

## 2026-01 — Half precision (FP16) support

### Changes
- `lib/hp_math_intrinsics.F90`: new module providing SP-promoted wrappers for 10 intrinsic math functions (`cos`, `sin`, `abs`, `tan`, `sqrt`, `asin`, `atan2`, `sign`, `exp`, `log`) that lack native FP16 overloads in nvfortran. Active only under `#ifdef USE_HALF_PRECISION`.
- CMakeLists.txt: added `USE_HALF_PRECISION` option with `FATAL_ERROR` for GNU and Intel (no `real(kind=2)` support); NVIDIA uses preprocessor define instead of `-r2` flag

### Compiler support
- GNU: DP, SP only (no `real(kind=2)`)
- Intel: DP, SP only (no `-r2` flag)
- NVIDIA: DP, SP, HP (`real(kind=2)` supported, intrinsic wrappers required)

---

## 2026-01 — Single and double precision support

### Changes
- `lib/oce_modules.F90`: introduced three-way preprocessor for WP: `USE_HALF_PRECISION` sets WP=2, `USE_SINGLE_PRECISION` sets WP=4, default WP=8
- `lib/gen_halo_exchange.F90`: replaced 58 hardcoded `real*8` declarations with `real(kind=WP)`
- `lib/oce_mesh.F90`, `lib/gen_modules_partitioning.F90`: replaced hardcoded MPI datatypes with precision-aware `MPI_WP`
- CMakeLists.txt: added `USE_SINGLE_PRECISION` option; configured precision-specific compiler flags for all compilers (Intel `-r4`/`-r8`, GNU `-fdefault-real-8`, NVHPC `-r4`/`-r8`, Cray `-s real32`/`-s real64`)
- `configure.sh`: added `--precision dp|sp|hp` flag; build directories named `build_<compiler>_<precision>`

---

## 2025-10-10 — Analytic mesh dwarf

### New files
- `src/fesom_analytic.F90`: standalone entry point that generates a rectangular triangular mesh in-memory with configurable grid size, vertical levels, and boundary conditions (closed-wall or doubly-periodic). No mesh files or restart data required.
- `lib/analytic_mesh.F90`: mesh generation module producing a structured triangular grid with topology, coordinates, areas, and connectivity arrays
- `lib/mesh_output.F90`: binary output of mesh arrays and per-step tracer scalar fields for Python visualization

### Features
- Command-line arguments: `nx`, `ny`, `nl`, `--periodic`, `--save-mesh`, `--save-scalars`, `--output-dir`
- Domain: 100 km x 100 km, depth 1000 m
- Initial conditions: sin(x)*cos(y) temperature field + Gaussian salinity blob
- Tracer sum is perfectly conserved at every step (validates advection correctness)

---

## 2025-10-10 — Mesh-init dwarf and C interface

### Mesh-init dwarf
- `src/fesom_mesh_init.F90`: entry point that reads pre-partitioned mesh files (not restart files) and sets custom tracer initial conditions before advection
- `lib/tracer_init_from_mesh.F90`: initialization module that calls `mesh_setup()`, allocates dynamics/tracer/work arrays, and initializes FCT

### C interface
- `lib/tracer_c_interface.F90`: C-compatible wrappers (`iso_c_binding`) exposing tracer library to Python/ctypes. Functions include `tracer_init_mpi`, `tracer_init`, `tracer_advect_step`, `tracer_get_stats`, `tracer_get_values`, `tracer_finalize`, `tracer_load_mesh_partition`, `tracer_load_mesh_from_files`, array-based initialization functions, and a complete workflow function
- Built conditionally via `-DENABLE_TRACER_C_INTERFACE=ON`

---

## 2025-10-09 — Initial shared library refactoring

### Build system
- Replaced `dwarf_linkfiles.sh` symlink workflow with direct CMake source references
- New top-level `CMakeLists.txt` producing shared library (`libfesom_tracer_Fortran.so`) and three executables
- `configure.sh` script for compiler/platform selection
- Compiler-specific flags for Intel, GNU, NVHPC, Cray
- OpenACC support for NVIDIA builds

### Stubs
- `diagnostics_stub.F90`: minimal module providing `ldiag_DVD = .false.` to avoid pulling in full `gen_modules_diag.F90` with 15+ dependencies

### Executable
- `fesom_tracer` (restart-based): reads binary restart files, runs 10 advection steps, writes updated restarts

### Structure
- Shared library contains all modules and advection routines
- Executables link against the library
- Standard `bin/` + `lib/` output layout with install targets
