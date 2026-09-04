# FESOM2 Dwarf Tracer: Mesh Data Structures Reference

## Overview

The FESOM2 tracer advection dwarf operates on an unstructured triangular mesh with
a structured vertical (z-level) discretization. The mesh is partitioned across MPI
ranks, with halo exchange for parallel communication.

This document describes the mesh data structures as consumed by the tracer advection
kernel, and serves as a reference for building an analytic mesh generator that removes
the dependency on pre-partitioned mesh files.

## Current Data Flow

For `npes=1`, `read_mesh()` derives identity node, element, and edge ownership
from the global mesh headers and initializes empty halo communication lists.
For `npes>1`, it reads the pre-partitioned metadata from `dist_N/`.

```
Pre-partitioned files (dist_N/)         Global mesh files
  rpart.out                               nod2d.out
  my_listXXXXX.out                        elem2d.out
  com_infoXXXXX.out                       aux3d.out / depth_zlev.out
                                          edges.out / edge_tri.out / edgenum.out
                                          elvls.out / nlvls.out
        |                                       |
        v                                       v
   read_mesh() -----> populate t_partit + t_mesh (raw)
        |
        v
   test_tri()              -- ensure CCW node ordering in triangles
   load_edges()            -- read edge/edge_tri, build elem_edges
   find_neighbors()        -- build elem_neighbors, nod_in_elem2D
   find_levels()           -- read nlevels, nlevels_nod2D, set ulevels
   find_levels_min_e2n()   -- nlevels_nod2D_min, ulevels_nod2D_max
   mesh_areas()            -- elem_area, area, areasvol, mesh_resolution
   mesh_auxiliary_arrays() -- edge_dxdy, edge_cross_dxdy, gradient_sca,
                              gradient_vec, metric_factor, elem_cos, coriolis
        |
        v
   tracer_init_mesh_and_arrays()  -- allocate hnode, helem, dynamics, tracers
        |
        v
   do_oce_adv_tra()               -- run advection kernel
```

## Key Data Structures

### `t_partit` -- Partition / MPI layout

| Field | Type | Description |
|---|---|---|
| `npes` | int | Total number of MPI ranks |
| `mype` | int | This rank's ID (0-based) |
| `MPI_COMM_FESOM` | int | MPI communicator |
| `part(npes+1)` | int(:) | Cumulative node count per rank |
| `myDim_nod2D` | int | Number of owned nodes |
| `eDim_nod2D` | int | Number of halo nodes |
| `myList_nod2D(:)` | int(:) | Global indices of local nodes (size: myDim+eDim) |
| `myDim_elem2D` | int | Number of owned elements |
| `eDim_elem2D` | int | Number of halo elements |
| `eXDim_elem2D` | int | Number of extended-halo elements |
| `myList_elem2D(:)` | int(:) | Global indices of local elements |
| `myDim_edge2D` | int | Number of owned edges |
| `eDim_edge2D` | int | Number of halo edges |
| `myList_edge2D(:)` | int(:) | Global indices of local edges |
| `com_nod2D` | com_struct | Halo exchange info for nodes |
| `com_elem2D` | com_struct | Halo exchange info for elements |
| `com_elem2D_full` | com_struct | Full-halo exchange info for elements |
| `s_mpitype_*, r_mpitype_*` | int arrays | MPI derived datatypes for exchange |

#### `com_struct` -- Communication pattern

| Field | Description |
|---|---|
| `rPEnum` | Number of ranks we receive from |
| `rPE(MAX_NEIGHBOR_PARTITIONS)` | Rank IDs we receive from |
| `rptr(MAX_NEIGHBOR_PARTITIONS+1)` | Pointer into rlist per rank |
| `rlist(:)` | Local indices of nodes/elems to receive |
| `sPEnum` | Number of ranks we send to |
| `sPE(MAX_NEIGHBOR_PARTITIONS)` | Rank IDs we send to |
| `sptr(MAX_NEIGHBOR_PARTITIONS+1)` | Pointer into slist per rank |
| `slist(:)` | Local indices of nodes/elems to send |

### `t_mesh` -- Mesh geometry

#### Scalars

| Field | Description |
|---|---|
| `nod2D` | Total (global) number of 2D nodes |
| `elem2D` | Total (global) number of 2D elements |
| `edge2D` | Total (global) number of 2D edges |
| `edge2D_in` | Number of internal (non-boundary) edges |
| `nl` | Number of vertical levels (layer interfaces) |
| `ocean_area` | Total ocean surface area (m^2) |

#### Horizontal topology (local indices)

| Field | Shape | Description |
|---|---|---|
| `coord_nod2D(2, nNod)` | real | Node coordinates in radians (rotated frame) |
| `geo_coord_nod2D(2, nNod)` | real | Node coordinates in geographic frame |
| `elem2D_nodes(3, nElem)` | int | 3 node indices per triangle (CCW ordering) |
| `edges(2, nEdge)` | int | 2 node indices per edge |
| `edge_tri(2, nEdge)` | int | 2 element indices per edge (0 = boundary) |
| `elem_edges(3, nElem)` | int | 3 edge indices per element |
| `elem_neighbors(3, nElem)` | int | 3 neighboring element indices (<=0 = boundary) |
| `nod_in_elem2D(maxAdj, nNod)` | int | Elements surrounding each node |
| `nod_in_elem2D_num(nNod)` | int | Count of elements per node |

Where `nNod = myDim_nod2D + eDim_nod2D`, `nElem = myDim_elem2D + eDim_elem2D + eXDim_elem2D`, `nEdge = myDim_edge2D + eDim_edge2D`.

#### Vertical structure

| Field | Shape | Description |
|---|---|---|
| `zbar(nl)` | real | Depth of layer interfaces (negative, surface=0) |
| `Z(nl-1)` | real | Mid-depth of each layer |
| `depth(nNod or nElem)` | real | Bottom depth at nodes or elements |
| `nlevels(nElem)` | int | Number of active levels at each element |
| `nlevels_nod2D(nNod)` | int | Number of active levels at each node |
| `nlevels_nod2D_min(nNod)` | int | Min nlevels among surrounding elements |
| `ulevels(nElem)` | int | Upper level index (1 = no cavity) |
| `ulevels_nod2D(nNod)` | int | Upper level index at nodes |
| `ulevels_nod2D_max(nNod)` | int | Max ulevels among surrounding elements |

#### Geometric quantities

| Field | Shape | Description |
|---|---|---|
| `elem_area(nElem)` | real | Triangle area (m^2) |
| `area(nl, nNod)` | real | Scalar cell upper/lower face area (m^2) |
| `areasvol(nl, nNod)` | real | Scalar cell mid-area for volume (m^2) |
| `area_inv(nl, nNod)` | real | 1/area |
| `areasvol_inv(nl, nNod)` | real | 1/areasvol |
| `mesh_resolution(nNod)` | real | Effective mesh resolution at nodes (m) |
| `elem_cos(nElem)` | real | cos(latitude) at element center |
| `metric_factor(nElem)` | real | tan(lat)/R_earth at element center |
| `edge_dxdy(2, nEdge)` | real | Edge vector (dx,dy) in radians |
| `edge_cross_dxdy(4, nEdge)` | real | Cross-distances from edge midpoint to element centers (m) |
| `gradient_sca(6, nElem)` | real | Scalar gradient reconstruction coefficients |
| `gradient_vec(6, nElem)` | real | Vector gradient reconstruction coefficients (least-squares) |
| `coriolis(nElem)` | real | Coriolis parameter at elements |
| `coriolis_node(nNod)` | real | Coriolis parameter at nodes |

#### ALE / layer thickness arrays (allocated in `tracer_init_mesh_and_arrays`)

| Field | Shape | Description |
|---|---|---|
| `hnode(nl-1, nNod)` | real | Layer thickness at nodes |
| `hnode_new(nl-1, nNod)` | real | Updated layer thickness |
| `helem(nl-1, nElem)` | real | Layer thickness at elements |
| `hbar(nNod)` | real | Sea surface elevation |
| `hbar_old(nNod)` | real | Previous elevation |
| `dhe(nElem)` | real | Depth increment on elements |
| `zbar_3d_n(nl, nNod)` | real | 3D layer boundary depths at nodes |
| `Z_3d_n(nl-1, nNod)` | real | 3D mid-layer depths at nodes |

### `t_dyn` -- Dynamics (velocity fields)

| Field | Shape | Description |
|---|---|---|
| `uv(2, nl-1, nElem)` | real | Horizontal velocity (u,v) at elements per layer |
| `w(nl, nNod)` | real | Vertical velocity at nodes per level |

### `t_tracer` -- Tracer fields

| Field | Description |
|---|---|
| `num_tracers` | Number of tracers (default 2: T, S) |
| `data(:)` | Array of `t_tracer_data` |
| `work` | `t_tracer_work` shared workspace |

#### `t_tracer_data` -- Per-tracer arrays

| Field | Shape | Description |
|---|---|---|
| `values(nl-1, nNod)` | real | Current tracer values |
| `valuesAB(nl-1, nNod)` | real | Adams-Bashforth interpolated values |
| `valuesold(nl-1, nNod, AB_order)` | real | Previous timestep values |
| `tra_adv_hor` | char(20) | Horizontal advection scheme (`UPW1`, `MUSCL`, `MFCT`) |
| `tra_adv_ver` | char(20) | Vertical advection scheme (`UPW1`, `QR4C`, `CDIFF`, `PPM`) |
| `tra_adv_lim` | char(20) | Limiter (`NONE`, `FCT`) |

#### `t_tracer_work` -- Shared work arrays

| Field | Shape | Description |
|---|---|---|
| `del_ttf(nl-1, nNod)` | real | Tracer tendency |
| `del_ttf_advhoriz(nl-1, nNod)` | real | Horizontal advection tendency |
| `del_ttf_advvert(nl-1, nNod)` | real | Vertical advection tendency |
| `fct_LO(nl-1, nNod)` | real | FCT low-order solution |
| `adv_flux_hor(nl-1, nEdge)` | real | Horizontal advective flux |
| `adv_flux_ver(nl, nNod)` | real | Vertical advective flux |
| `fct_ttf_max/min(nl-1, nNod)` | real | FCT bounds |
| `fct_plus/minus(nl-1, nNod)` | real | FCT limiters |

## Simplifications for Analytic Toy Mesh

The following simplifications reduce the mesh to a minimal but fully functional state:

| Simplification | Effect |
|---|---|
| `cartesian = .true.` | `elem_cos = 1`, `metric_factor = 0`, no grid rotation |
| `use_cavity = .false.` | `ulevels = 1` everywhere, no cavity arrays |
| `force_rotation = .false.` | Coordinates used as-is (no rotated/geographic transform) |
| Uniform depth | `nlevels = nl` and `nlevels_nod2D = nl` everywhere |
| `npes = 1` | No halo exchange needed: `eDim_* = 0`, `eXDim_elem2D = 0`, trivial `com_struct` |
| Flat bottom | No partial cells, `hnode = zbar(k) - zbar(k+1)` for all nodes |

## Mesh Files Reference (current)

| File | Content |
|---|---|
| `nod2d.out` | `nod2D` then per-node: `id lon lat boundary_flag` |
| `elem2d.out` | `elem2D` then per-elem: `n1 n2 n3` |
| `aux3d.out` | `nl`, `zbar(1:nl)`, then per-node bottom depth |
| `edges.out` | Per-edge: `n1 n2` |
| `edge_tri.out` | Per-edge: `e1 e2` |
| `edgenum.out` | `edge2D`, `edge2D_in` |
| `elvls.out` | Per-element: number of levels |
| `nlvls.out` | Per-node: number of levels |
| `dist_N/rpart.out` | `npes`, then node count per rank |
| `dist_N/my_listXXXXX.out` | Per-rank: local node/elem/edge lists |
| `dist_N/com_infoXXXXX.out` | Per-rank: halo communication patterns |
