# Vertex-Parallel OpenMP Strategy

## Overview

Replaces the lock-based edge-scatter OpenMP pattern with a vertex-gather
approach. Each vertex loops over its incident edges and gathers flux
contributions, writing only to its own RHS array element. This is
inherently race-free — no locks, atomics, or edge coloring needed.

## How it works

The original edge-scatter pattern:
```fortran
!$OMP PARALLEL DO
do edge = 1, myDim_edge2D
    ! compute flux on edge
    ! lock(node1); accumulate +flux into node1; unlock(node1)
    ! lock(node2); accumulate -flux into node2; unlock(node2)
end do
```

The vertex-gather replacement:
```fortran
!$OMP PARALLEL DO
do n = 1, myDim_nod2D
    do ie = 1, node_edge_num(n)
        edge = node_edge_idx(ie, n)
        sgn  = node_edge_sign(ie, n)   ! +1 or -1
        ! accumulate sgn * flux(edge) into node n
    end do
end do
```

Each thread writes only to its own vertex — no synchronization needed.

## Node-to-edge map

A CSR-like mapping built once during initialization (`build_node_edge_map`):

| Array              | Shape                        | Description                          |
|--------------------|------------------------------|--------------------------------------|
| `node_edge_num(n)` | `(total_nod)`                | Number of incident edges per node    |
| `node_edge_idx(ie,n)` | `(max_edges_per_node, total_nod)` | Edge indices for each node |
| `node_edge_sign(ie,n)` | `(max_edges_per_node, total_nod)` | +1 if node is 1st endpoint, -1 if 2nd |

For triangular meshes, `max_edges_per_node` is typically 6 (interior nodes).

## Converted scatter patterns

Three hot-path scatters in the advection code were converted:

1. **fct_LO scatter** (`oce_adv_tra_driver.F90`, `do_oce_adv_tra`):
   Edge flux accumulated into low-order solution at both endpoints.

2. **dttf_h scatter** (`oce_adv_tra_driver.F90`, `oce_tra_adv_flux2dtracer`):
   Horizontal flux contribution to tracer tendency at both endpoints.

3. **fct_plus/fct_minus scatter** (`oce_adv_tra_fct.F90`, `oce_tra_adv_fct`):
   Positive/negative antidiffusive flux accumulation for FCT limiting.

One init-time scatter was also converted:

4. **nn_pos/nboundary_lay** (`oce_muscl_adv.F90`, `muscl_adv_init`):
   Neighbor list and boundary layer detection.

OpenACC code paths are preserved unchanged (still use atomics).

## Build

```bash
./configure.sh --compiler gnu --precision dp --openmp --clean --build
cd build_gnu_dp && OMP_NUM_THREADS=4 ./run.sh 1 200 200 10 --periodic
```

## MPI process binding (critical)

OpenMPI defaults to `--bind-to core`, which pins all OpenMP threads to a
**single CPU core**. This completely kills thread scaling — all threads
compete for one core, showing 100% CPU in `top` instead of N×100%.

The generated `run.sh` uses `--bind-to none` to fix this. If calling
`mpirun` directly, always use one of:

```bash
mpirun --bind-to none -np 1 ./bin/fesom_tracer_analytic ...   # no binding
mpirun --bind-to socket -np 1 ./bin/fesom_tracer_analytic ... # bind to socket
```

You can verify binding with `--report-bindings`:
```bash
mpirun --report-bindings --bind-to none -np 1 ./bin/fesom_tracer_analytic ...
```

## Key properties

- **Bit-identical results** regardless of thread count (deterministic accumulation order)
- **No locks or atomics** — fully race-free by construction
- **No edge coloring** infrastructure needed
- **Non-OpenMP builds** work unchanged (gather loops run serially)

## Performance

### 200x200x3 mesh, 200 timesteps, GNU gfortran, `--bind-to none`

| Threads | advection_loop (s) | Speedup | adv_tra_hor (s) | Speedup |
|---------|-------------------|---------|-----------------|---------|
| 1       | 1.99              | 1.00x   | 0.67            | 1.00x   |
| 2       | 1.29              | 1.54x   | 0.39            | 1.73x   |
| 4       | 0.92              | 2.17x   | 0.23            | 2.95x   |
| 8       | 0.77              | 2.59x   | 0.16            | 4.26x   |

### 500x500x3 mesh, 200 timesteps, GNU gfortran, `--bind-to none`

| Threads | advection_loop (s) | Speedup | adv_tra_hor (s) | Speedup |
|---------|-------------------|---------|-----------------|---------|
| 1       | 12.86             | 1.00x   | 4.39            | 1.00x   |
| 2       | 9.23              | 1.39x   | 3.02            | 1.45x   |
| 4       | 6.35              | 2.03x   | 1.79            | 2.46x   |
| 8       | 5.88              | 2.19x   | 1.60            | 2.75x   |

`adv_tra_hor` (edge flux computation) scales near-linearly at smaller mesh
sizes. At 500x500, scaling flattens at 8 threads due to memory bandwidth
saturation. The overall `advection_loop` speedup is limited by serial
portions (halo exchange, unprofiled init sections).
