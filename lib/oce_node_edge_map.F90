!> @file oce_node_edge_map.F90
!! @brief Node-to-edge adjacency map for vertex-parallel OpenMP
!! @details Builds an inverse mapping from nodes to their incident edges,
!!          enabling race-free vertex-gather parallelization. Each vertex
!!          gathers flux contributions from its incident edges, writing
!!          only to its own RHS — no locks or atomics needed.

module oce_node_edge_map_module
  use o_PARAM, only: WP, MP
  implicit none
  save
  private
  public :: build_node_edge_map
  public :: node_edge_built, max_node_edges
  public :: node_edge_num, node_edge_idx, node_edge_sign

  logical :: node_edge_built = .false.
  integer :: max_node_edges = 0

  !> Number of incident edges per node
  integer, allocatable :: node_edge_num(:)       ! (total_nod)
  !> Edge indices for each node
  integer, allocatable :: node_edge_idx(:,:)     ! (max_node_edges, total_nod)
  !> Sign: +1 if node is 1st endpoint of edge, -1 if 2nd endpoint
  integer, allocatable :: node_edge_sign(:,:)    ! (max_node_edges, total_nod)

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

end module oce_node_edge_map_module
