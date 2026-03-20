# OpenMP Scaling Performance

Benchmark comparing three OpenMP parallelization strategies for the
edge-to-node scatter operation in tracer advection.

## Setup

- **Mesh**: 500x500x3 (249,001 nodes, 2 vertical levels, doubly-periodic)
- **Timesteps**: 200 (400 tracer calls: 2 tracers x 200 steps)
- **Compiler**: GNU gfortran, `-O3 -fopenmp`, double precision
- **MPI binding**: `mpirun --bind-to none` (critical — see below)
- **Hardware**: Single MPI rank, varying OMP_NUM_THREADS

## Branches

| Branch | Strategy | Synchronization |
|--------|----------|-----------------|
| `main` | Edge loop with `omp_set_lock`/`omp_unset_lock` | Per-node locks |
| `check-openmp` | Edge loop grouped by edge coloring | None (colors avoid conflicts) |
| `vertex-parallel` | Node loop gathering from incident edges | None (each thread owns its node) |

## Results: advection_loop (total advection time)

| Threads | main (locks) | check-openmp (coloring) | vertex-parallel (gather) |
|---------|-------------|------------------------|--------------------------|
| 1       | 17.00s      | 12.35s                 | 12.86s                   |
| 2       | 11.91s      | 9.86s                  | 9.23s                    |
| 4       | 8.25s       | 7.37s                  | 6.35s                    |
| 8       | 7.17s       | 6.81s                  | 5.88s                    |

### Speedup (relative to own 1-thread baseline)

| Threads | main (locks) | check-openmp (coloring) | vertex-parallel (gather) |
|---------|-------------|------------------------|--------------------------|
| 1       | 1.00x       | 1.00x                  | 1.00x                    |
| 2       | 1.43x       | 1.25x                  | 1.39x                    |
| 4       | 2.06x       | 1.68x                  | 2.03x                    |
| 8       | 2.37x       | 1.81x                  | 2.19x                    |

## Results: adv_tra_hor (edge flux computation only)

This is the same code in all branches — only the scatter differs.

| Threads | main | check-openmp | vertex-parallel |
|---------|------|-------------|-----------------|
| 1       | 4.42s | 4.41s      | 4.39s           |
| 2       | 3.05s | 3.05s      | 3.02s           |
| 4       | 1.80s | 1.83s      | 1.79s           |
| 8       | 1.53s | 1.54s      | 1.60s           |

## Detailed profiling breakdown (vertex-parallel branch)

The advection loop has four distinct phases. The vertex-gather
(`flux2dtracer`) is the best-scaling component:

| Component | 1 thread | 4 threads | 8 threads | Speedup (8t) |
|-----------|----------|-----------|-----------|-------------|
| `adv_tra_hor` (edge flux computation) | 4.40s | 1.80s | 1.53s | **2.88x** |
| `flux2dtracer` (vertex-gather scatter) | 4.85s | 1.50s | 1.00s | **4.85x** |
| `adv_tra_ver` (vertical flux) | 0.64s | 0.37s | 0.35s | 1.82x |
| remainder (driver zero/accum loops) | 3.03s | 2.70s | 2.76s | 1.10x |

### Where the performance comes from

1. **Lock elimination** is the biggest single gain. At 1 thread, main takes
   17.0s vs 12.9s for vertex-parallel — a 4.1s penalty from `omp_set_lock`/
   `omp_unset_lock` even with zero contention.

2. **The vertex-gather scales 4.85x at 8 threads** — better than the edge
   flux computation (2.88x). This is because the gather is a pure streaming
   operation with perfect data locality: each thread writes only to its own
   nodes, no false sharing, no synchronization.

3. **The edge flux computation** (`adv_tra_hor`) is the same code in all
   branches and scales identically. It uses `!$OMP PARALLEL DO` over edges
   and is naturally race-free (each edge writes to its own flux slot).

4. **The driver loops** (zeroing, accumulation, tendency application) barely
   scale (1.10x) — they are simple streaming operations that are
   memory-bandwidth-limited.

### Potential optimization: fused vertex loop

Currently, advection has two separate phases:
1. **Edge loop** (`adv_tra_hor`): compute flux per edge → write to `flux(nz, edge)`
2. **Vertex loop** (`flux2dtracer`): gather from incident edges → accumulate to `rhs(nz, node)`

These could be fused into a single vertex loop that computes the flux
on-the-fly during gather, eliminating the intermediate flux array. Each
edge flux would be computed twice (once per endpoint), but the memory
savings from eliminating the flux array could more than compensate — especially
since the vertex-gather already scales 4.85x.

## Critical: MPI process binding

OpenMPI defaults to `--bind-to core`. With hybrid MPI+OpenMP, this pins
all threads to a single core, completely destroying scaling. This was the
original reason scaling appeared broken — not the parallelization strategy.

```bash
# BAD: all threads on 1 core (default)
mpirun -np 1 ./bin/fesom_tracer_analytic ...

# GOOD: threads can use all cores
mpirun --bind-to none -np 1 ./bin/fesom_tracer_analytic ...

# Verify binding:
mpirun --report-bindings --bind-to none -np 1 ./bin/fesom_tracer_analytic ...
```

![OpenMP Scaling](omp_scaling.png)
