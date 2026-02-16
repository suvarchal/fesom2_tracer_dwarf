!> @brief Edge coloring for lock-free OpenMP edge-to-node scatter
!> @details Partitions edges into color groups where no two same-color edges
!>          share a node. This allows lock-free parallel processing of each
!>          color group. Standard FEM assembly technique for unstructured meshes.
module oce_edge_coloring_module
    implicit none
    private
    public :: t_edge_coloring, compute_edge_coloring

    type t_edge_coloring
        integer :: num_colors = 0
        integer, allocatable :: edge_color(:)        ! (myDim_edge2D) color of each edge
        integer, allocatable :: color_edges(:,:)     ! (max_edges_per_color, num_colors) edge indices per color
        integer, allocatable :: color_edge_count(:)  ! (num_colors) number of edges in each color
    end type t_edge_coloring

contains

!> Compute edge coloring using greedy algorithm
subroutine compute_edge_coloring(coloring, myDim_edge2D, myDim_nod2D, eDim_nod2D, edges)
    implicit none
    type(t_edge_coloring), intent(inout) :: coloring
    integer, intent(in) :: myDim_edge2D
    integer, intent(in) :: myDim_nod2D
    integer, intent(in) :: eDim_nod2D
    integer, intent(in) :: edges(:,:)

    integer :: e, n1, n2, c, max_color, total_nodes
    integer, allocatable :: node_color_used(:,:)
    integer :: max_node_degree
    integer, allocatable :: node_degree(:)
    integer, allocatable :: temp_edges(:)
    integer :: i, j, nc, key_edge, key_node, cur_node

    if (myDim_edge2D == 0) return

    total_nodes = myDim_nod2D + eDim_nod2D

    ! Compute max node degree to size the color-used array
    allocate(node_degree(total_nodes))
    node_degree = 0
    do e = 1, myDim_edge2D
        n1 = edges(1, e)
        n2 = edges(2, e)
        if (n1 >= 1 .and. n1 <= total_nodes) node_degree(n1) = node_degree(n1) + 1
        if (n2 >= 1 .and. n2 <= total_nodes) node_degree(n2) = node_degree(n2) + 1
    end do
    max_node_degree = maxval(node_degree)
    deallocate(node_degree)

    ! Vizing's theorem: edge chromatic number <= max_degree + 1
    ! Double it for safety with greedy algorithm
    max_color = 2 * max_node_degree + 2

    ! node_color_used(c, node) = 1 if color c is used by an edge touching node
    allocate(node_color_used(max_color, total_nodes))
    node_color_used = 0

    ! Allocate edge_color array
    if (allocated(coloring%edge_color)) deallocate(coloring%edge_color)
    allocate(coloring%edge_color(myDim_edge2D))
    coloring%edge_color = 0

    ! Greedy coloring: assign each edge the smallest color not used by adjacent edges
    coloring%num_colors = 0
    do e = 1, myDim_edge2D
        n1 = edges(1, e)
        n2 = edges(2, e)

        ! Find smallest available color
        do c = 1, max_color
            if (node_color_used(c, n1) == 0 .and. node_color_used(c, n2) == 0) then
                coloring%edge_color(e) = c
                node_color_used(c, n1) = 1
                node_color_used(c, n2) = 1
                if (c > coloring%num_colors) coloring%num_colors = c
                exit
            end if
        end do
    end do

    deallocate(node_color_used)

    ! Build color_edge_count and color_edges arrays
    if (allocated(coloring%color_edge_count)) deallocate(coloring%color_edge_count)
    allocate(coloring%color_edge_count(coloring%num_colors))
    coloring%color_edge_count = 0

    ! Count edges per color
    do e = 1, myDim_edge2D
        c = coloring%edge_color(e)
        coloring%color_edge_count(c) = coloring%color_edge_count(c) + 1
    end do

    ! Allocate color_edges with max edges per color
    if (allocated(coloring%color_edges)) deallocate(coloring%color_edges)
    allocate(coloring%color_edges(maxval(coloring%color_edge_count), coloring%num_colors))
    coloring%color_edges = 0

    ! Fill color_edges
    allocate(temp_edges(coloring%num_colors))
    temp_edges = 0
    do e = 1, myDim_edge2D
        c = coloring%edge_color(e)
        temp_edges(c) = temp_edges(c) + 1
        coloring%color_edges(temp_edges(c), c) = e
    end do
    deallocate(temp_edges)

    ! Sort edges within each color by min node index for cache locality
    do c = 1, coloring%num_colors
        nc = coloring%color_edge_count(c)
        ! Insertion sort by min(node1, node2)
        do i = 2, nc
            key_edge = coloring%color_edges(i, c)
            key_node = min(edges(1, key_edge), edges(2, key_edge))
            j = i - 1
            do while (j >= 1)
                cur_node = min(edges(1, coloring%color_edges(j, c)), &
                               edges(2, coloring%color_edges(j, c)))
                if (cur_node <= key_node) exit
                coloring%color_edges(j + 1, c) = coloring%color_edges(j, c)
                j = j - 1
            end do
            coloring%color_edges(j + 1, c) = key_edge
        end do
    end do

    write(*, '(A,I4,A,I8,A)') '  Edge coloring: ', coloring%num_colors, &
        ' colors for ', myDim_edge2D, ' edges (sorted by node locality)'

end subroutine compute_edge_coloring

end module oce_edge_coloring_module
