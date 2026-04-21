!> @file oce_node_edge_map.F90
!! @brief Node-to-edge adjacency map for vertex-parallel OpenMP
!! @details Builds an inverse mapping from nodes to their incident edges,
!!          enabling race-free vertex-gather parallelization. Each vertex
!!          gathers flux contributions from its incident edges, writing
!!          only to its own RHS — no locks or atomics needed.
!!
!!          The pre-packed arrays (node_edge_dxdy, node_edge_elems, etc.)
!!          store static mesh geometry in node-contiguous layout, eliminating
!!          indirect indexing through node_edge_idx → edge → mesh arrays.

module oce_node_edge_map_module
  use o_PARAM, only: WP, MP
  implicit none
  save
  private
  public :: build_node_edge_map, build_node_edge_packed
  public :: node_edge_built, max_node_edges
  public :: node_edge_num, node_edge_idx, node_edge_sign
  public :: node_edge_dxdy, node_edge_elems, node_edge_enodes
  public :: node_edge_nl1, node_edge_nu1, node_edge_nl2, node_edge_nu2

  logical :: node_edge_built = .false.
  integer :: max_node_edges = 0

  !> Number of incident edges per node
  integer, allocatable :: node_edge_num(:)       ! (total_nod)
  !> Edge indices for each node
  integer, allocatable :: node_edge_idx(:,:)     ! (max_node_edges, total_nod)
  !> Sign: +1 if node is 1st endpoint of edge, -1 if 2nd endpoint
  integer, allocatable :: node_edge_sign(:,:)    ! (max_node_edges, total_nod)

  !> Pre-packed static geometry per (edge_slot, node) — eliminates indirect indexing
  real(kind=MP), allocatable :: node_edge_dxdy(:,:,:)    ! (4, max_node_edges, total_nod)
  integer,       allocatable :: node_edge_elems(:,:,:)   ! (2, max_node_edges, total_nod)
  integer,       allocatable :: node_edge_enodes(:,:,:)  ! (2, max_node_edges, total_nod)
  !> Pre-resolved per-element level bounds (eliminates edge_tri → nlevels double indirection)
  integer,       allocatable :: node_edge_nl1(:,:)       ! (max_node_edges, total_nod)
  integer,       allocatable :: node_edge_nu1(:,:)       ! (max_node_edges, total_nod)
  integer,       allocatable :: node_edge_nl2(:,:)       ! (max_node_edges, total_nod)
  integer,       allocatable :: node_edge_nu2(:,:)       ! (max_node_edges, total_nod)

contains

  subroutine build_node_edge_map(partit, mesh)
    use MOD_MESH
    use MOD_PARTIT
    implicit none
    type(t_mesh),   intent(in) :: mesh
    type(t_partit), intent(in) :: partit

    integer :: edge, n1, n2, total_nod
    integer, allocatable :: cnt(:)

    total_nod = partit%myDim_nod2D + partit%eDim_nod2D

    ! First pass: count edges per node
    allocate(cnt(total_nod))
    cnt = 0
    do edge = 1, partit%myDim_edge2D
      n1 = mesh%edges(1, edge)
      n2 = mesh%edges(2, edge)
      cnt(n1) = cnt(n1) + 1
      cnt(n2) = cnt(n2) + 1
    end do

    max_node_edges = maxval(cnt)

    ! Allocate mapping arrays
    if (allocated(node_edge_num))  deallocate(node_edge_num)
    if (allocated(node_edge_idx))  deallocate(node_edge_idx)
    if (allocated(node_edge_sign)) deallocate(node_edge_sign)

    allocate(node_edge_num(total_nod))
    allocate(node_edge_idx(max_node_edges, total_nod))
    allocate(node_edge_sign(max_node_edges, total_nod))

    node_edge_num  = 0
    node_edge_idx  = 0
    node_edge_sign = 0

    ! Second pass: fill the mapping
    do edge = 1, partit%myDim_edge2D
      n1 = mesh%edges(1, edge)
      n2 = mesh%edges(2, edge)

      node_edge_num(n1) = node_edge_num(n1) + 1
      node_edge_idx(node_edge_num(n1), n1)  = edge
      node_edge_sign(node_edge_num(n1), n1) = 1    ! 1st endpoint: +flux

      node_edge_num(n2) = node_edge_num(n2) + 1
      node_edge_idx(node_edge_num(n2), n2)  = edge
      node_edge_sign(node_edge_num(n2), n2) = -1   ! 2nd endpoint: -flux
    end do

    node_edge_built = .true.

    deallocate(cnt)

    if (partit%mype == 0) then
      write(*, '(A,I4,A,I8,A)') '  Node-edge map: max ', max_node_edges, &
        ' edges/node, ', total_nod, ' nodes'
    end if

  end subroutine build_node_edge_map

  !=========================================================================
  ! Pre-pack static mesh geometry into node-contiguous arrays.
  ! Must be called after build_node_edge_map.
  !=========================================================================
  subroutine build_node_edge_packed(partit, mesh)
    use MOD_MESH
    use MOD_PARTIT
    implicit none
    type(t_mesh),   intent(in) :: mesh
    type(t_partit), intent(in) :: partit

    integer :: n, e, edge, el1, el2, total_nod
    real(kind=8) :: mem_mb

    if (.not. node_edge_built) then
      write(*, '(A)') 'ERROR: build_node_edge_packed called before build_node_edge_map'
      stop 1
    end if

    total_nod = partit%myDim_nod2D + partit%eDim_nod2D

    ! Allocate pre-packed arrays
    if (allocated(node_edge_dxdy))   deallocate(node_edge_dxdy)
    if (allocated(node_edge_elems))  deallocate(node_edge_elems)
    if (allocated(node_edge_enodes)) deallocate(node_edge_enodes)
    if (allocated(node_edge_nl1))    deallocate(node_edge_nl1)
    if (allocated(node_edge_nu1))    deallocate(node_edge_nu1)
    if (allocated(node_edge_nl2))    deallocate(node_edge_nl2)
    if (allocated(node_edge_nu2))    deallocate(node_edge_nu2)

    allocate(node_edge_dxdy(4, max_node_edges, total_nod))
    allocate(node_edge_elems(2, max_node_edges, total_nod))
    allocate(node_edge_enodes(2, max_node_edges, total_nod))
    allocate(node_edge_nl1(max_node_edges, total_nod))
    allocate(node_edge_nu1(max_node_edges, total_nod))
    allocate(node_edge_nl2(max_node_edges, total_nod))
    allocate(node_edge_nu2(max_node_edges, total_nod))

    node_edge_dxdy   = 0.0_MP
    node_edge_elems  = 0
    node_edge_enodes = 0
    node_edge_nl1    = 0
    node_edge_nu1    = 0
    node_edge_nl2    = 0
    node_edge_nu2    = 0

    ! Pack geometry from edge-indexed arrays into node-contiguous layout
    do n = 1, total_nod
      do e = 1, node_edge_num(n)
        edge = node_edge_idx(e, n)

        ! Edge geometry
        node_edge_dxdy(:, e, n) = mesh%edge_cross_dxdy(:, edge)

        ! Adjacent elements
        node_edge_elems(:, e, n) = mesh%edge_tri(:, edge)

        ! Endpoint nodes
        node_edge_enodes(:, e, n) = mesh%edges(:, edge)

        ! Pre-resolve level bounds from element indices
        el1 = mesh%edge_tri(1, edge)
        el2 = mesh%edge_tri(2, edge)

        node_edge_nl1(e, n) = mesh%nlevels(el1) - 1
        node_edge_nu1(e, n) = mesh%ulevels(el1)

        if (el2 > 0) then
          node_edge_nl2(e, n) = mesh%nlevels(el2) - 1
          node_edge_nu2(e, n) = mesh%ulevels(el2)
        else
          node_edge_nl2(e, n) = 0
          node_edge_nu2(e, n) = 0
        end if
      end do
    end do

    if (partit%mype == 0) then
      mem_mb = dble(max_node_edges) * dble(total_nod) * &
               (4.0d0 * MP + 4.0d0 * 4 + 4.0d0 * 4) / (1024.0d0 * 1024.0d0)
      write(*, '(A,F8.1,A)') '  Node-edge packed arrays: ', mem_mb, ' MB'
    end if

  end subroutine build_node_edge_packed

end module oce_node_edge_map_module
