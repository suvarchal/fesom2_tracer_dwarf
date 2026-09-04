!> @file atlas_fesom_mesh.F90
!! @brief Atlas-backed mesh setup wrapper for tracer dwarf
!! @details When ENABLE_ATLAS is enabled, this module initializes Atlas
!!          and generates a mesh-driven configuration path for the tracer
!!          dwarf. The resulting FESOM mesh arrays are populated through the
!!          existing in-memory analytic mesh generator. When Atlas is disabled,
!!          it falls back to the standard mesh_setup() file-based path.

module atlas_fesom_mesh_module
  use MOD_MESH
  use MOD_PARTIT
  use o_PARAM, only: WP, MP
  use oce_mesh_module, only: mesh_setup, test_tri, find_levels, find_levels_min_e2n, &
                              mesh_areas, mesh_auxiliary_arrays, read_vertical_grid
  use analytic_mesh_module, only: generate_analytic_mesh
  use g_config, only: force_rotation, MeshPath
  use g_rotate_grid, only: set_mesh_transform_matrix, g2r
  use par_support_interfaces, only: init_mpi_types, init_gatherLists
  use iso_fortran_env, only: output_unit

  use mpi
#ifdef ENABLE_ATLAS
  use atlas_module
  use, intrinsic :: iso_c_binding, only: c_int
#endif

  implicit none
  private

  public :: mesh_setup_with_atlas, compute_tracer_stats_atlas
  public :: atlas_fesom_enabled, atlas_fesom_active
  public :: atlas_halo_exchange_nodal
#ifdef ENABLE_ATLAS
  public :: atlas_mesh_to_fesom_mesh, set_atlas_stats_mesh, set_fesom_dist
  ! Module-level Atlas mesh (persists for stats computation)
  type(atlas_Mesh), save :: atlas_mesh_global
  type(atlas_functionspace_NodeColumns), save :: atlas_nodes_global
  logical, save :: atlas_mesh_active = .false.
  integer, save :: atlas_mesh_halo = 0
#endif

contains

  logical function atlas_fesom_enabled()
    character(len=32) :: value
    integer :: status
    logical, save :: checked = .false.
    logical, save :: enabled = .false.

#ifdef ENABLE_ATLAS
    if (.not. checked) then
      call get_environment_variable('ATLAS_FESOM', value, status=status)
      enabled = status == 0 .and. trim(value) == '1'
      checked = .true.
    end if
    atlas_fesom_enabled = enabled
#else
    atlas_fesom_enabled = .false.
#endif
  end function atlas_fesom_enabled

  logical function atlas_fesom_active()
#ifdef ENABLE_ATLAS
    atlas_fesom_active = atlas_mesh_active
#else
    atlas_fesom_active = .false.
#endif
  end function atlas_fesom_active

#ifdef ENABLE_ATLAS
  subroutine set_atlas_stats_mesh(mesh, halo)
    type(atlas_Mesh), intent(inout) :: mesh
    integer, intent(in) :: halo

    atlas_mesh_global = mesh
    atlas_mesh_halo = halo
    atlas_nodes_global = atlas_functionspace_NodeColumns(atlas_mesh_global, &
                                                         halo=atlas_mesh_halo)
    atlas_mesh_active = .true.
  end subroutine set_atlas_stats_mesh
#endif

  subroutine atlas_halo_exchange_nodal(nodal_data)
    real(kind=WP), intent(inout) :: nodal_data(:,:)
#ifdef ENABLE_ATLAS
    type(atlas_Field) :: field
    real(WP), pointer :: field_data(:,:)

    field = atlas_Field(nodal_data)
    call atlas_nodes_global%halo_exchange(field)

    call field%final()
#endif
  end subroutine atlas_halo_exchange_nodal

  !> @brief Mesh setup entry point with optional Atlas path
  !! @param[inout] partit Partition structure
  !! @param[inout] mesh Mesh structure
  subroutine mesh_setup_with_atlas(partit, mesh)
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh

#ifdef ENABLE_ATLAS
    type(atlas_Grid)             :: grid_obj
    type(atlas_MeshGenerator)   :: meshgen_obj
    type(atlas_Mesh)            :: mesh_obj
    type(atlas_GridDistribution) :: dist_obj
    type(atlas_Config)          :: meshgen_config

    integer :: io_stat
    character(len=256) :: atlas_grid_name
    integer :: ncells, nnodes, r
    integer :: nl_default
    real(kind=MP), allocatable :: zbar_default(:)
    type(atlas_mesh_Nodes) :: mesh_nodes_obj
    type(atlas_mesh_Cells) :: mesh_cells_obj
    logical :: have_dist

    if (.not. atlas_fesom_enabled()) then
      stop 'Atlas FESOM support not enabled'
    end if

    ! Use sensible defaults for fesom-pi grid
    nl_default = 10

    call read_vertical_grid(partit, trim(MeshPath)//'aux3d.out', nl_default, &
                            zbar_default, io_stat)
    if (io_stat /= 0) then
      nl_default = 10
      allocate(zbar_default(nl_default))
      zbar_default = [( -1000.0_MP * real(r-1, MP) / real(nl_default-1, MP), &
                        r=1,nl_default )]
    end if

    atlas_grid_name = 'fesom-pi'
    call get_environment_variable('ATLAS_GRID', atlas_grid_name, status=io_stat)
    if (io_stat /= 0 .or. len_trim(atlas_grid_name) == 0) then
      atlas_grid_name = 'fesom-pi'
    end if

    ! Create the configured Atlas grid.
    if (partit%mype == 0) then
      write(output_unit, '(3A)') '  --> Attempting to use Atlas grid "', &
        trim(atlas_grid_name), '"...'
    end if

    grid_obj = atlas_Grid(trim(atlas_grid_name))
    write(output_unit, *) '  --> Atlas grid found: ', grid_obj%name()
    call set_fesom_dist(partit, grid_obj, dist_obj, have_dist)

    ! Create mesh generator and generate mesh
    meshgen_config = atlas_Config()
    call meshgen_config%set('triangulate', .true.) ! no quad element support for FESOM mesh available
    if (have_dist) then
      mesh_obj = atlas_Mesh(grid_obj, dist_obj, meshgen_config)
      call dist_obj%FINAL()
    else
      mesh_obj = atlas_Mesh(grid_obj, meshgen_config)
    end if

    atlas_nodes_global = atlas_functionspace_NodeColumns(mesh_obj, halo=2)

    ! Get mesh dimensions for logging
    mesh_nodes_obj = mesh_obj%nodes()
    mesh_cells_obj = mesh_obj%cells()
    nnodes = int(mesh_nodes_obj%size())
    ncells = int(mesh_cells_obj%size())

    if (partit%mype == 0) then
      write(output_unit, '(3A)') '  --> Setting up mesh from Atlas grid "', &
        trim(atlas_grid_name), '"'
      write(output_unit, '(A,I0,A,I0)') '      Atlas mesh: ', ncells, ' cells, ', nnodes, ' nodes'
      write(output_unit, '(A)') '      Converting to FESOM format...'
    end if

    ! Convert atlas mesh to FESOM mesh
    call atlas_mesh_to_fesom_mesh(mesh_obj, zbar_default, partit, mesh)
    partit%myDim_elem2D_shrinked = mesh%elem2D
    deallocate(zbar_default)

    ! Store mesh for later use by compute_tracer_stats_atlas
    atlas_mesh_global = mesh_obj
    atlas_mesh_active = .true.
    call meshgen_config%FINAL()
    call grid_obj%FINAL()
    call mesh_obj%FINAL()
#else
    stop 'Atlas support not enabled'
#endif

  end subroutine mesh_setup_with_atlas

#ifdef ENABLE_ATLAS
  !> @brief Set Atlas distribution from the standard FESOM owner lists
  subroutine set_fesom_dist(partit, grid_obj, dist_obj, have_dist)
    type(t_partit), intent(in) :: partit
    type(atlas_Grid), intent(inout) :: grid_obj
    type(atlas_GridDistribution), intent(out) :: dist_obj
    logical, intent(out) :: have_dist

    character(len=32) :: use_fesom_dist_str
    character(len=256) :: dist_mesh_dir, rpart_file, owner_file
    character(len=10) :: npes_string
    character(len=5) :: rank_string
    integer :: funit, r, node_idx, npes_in_file, error_status, ierror, io_stat
    integer :: file_rank, owned_count, ghost_count
    integer, allocatable :: part_csr(:), part_per_node(:), local_node_list(:)

    have_dist = .false.
    call get_environment_variable('ATLAS_USE_FESOM_DIST', use_fesom_dist_str, &
                                  status=io_stat)
    if (io_stat /= 0) use_fesom_dist_str = '0'
    if (trim(use_fesom_dist_str) /= '1') return

    write(npes_string, "(I10)") partit%npes
    dist_mesh_dir = trim(MeshPath)//'dist_'//trim(adjustl(npes_string))//'/'
    rpart_file = trim(dist_mesh_dir)//'rpart.out'
    error_status = 0
    if (partit%mype == 0) then
      open(newunit=funit, file=trim(rpart_file), action='read', status='old', &
           iostat=io_stat)
      if (io_stat /= 0) then
        write(output_unit, '(3A)') &
          '  WARNING: ATLAS_USE_FESOM_DIST=1 but cannot open ', &
          trim(rpart_file), ' - using default Atlas distribution'
        error_status = 1
      else
        allocate(part_csr(partit%npes+1))
        read(funit, *) npes_in_file
        if (npes_in_file /= partit%npes) error_status = 1
        part_csr(1) = 1
        read(funit, *) part_csr(2:partit%npes+1)
        do r = 2, partit%npes+1
          part_csr(r) = part_csr(r-1) + part_csr(r)
        end do
        close(funit)
      end if
    end if

    if (partit%npes > 1) then
      call MPI_BCast(error_status, 1, MPI_INTEGER, 0, &
                     partit%MPI_COMM_FESOM, ierror)
      if (error_status == 0) then
        if (partit%mype /= 0) allocate(part_csr(partit%npes+1))
        call MPI_BCast(part_csr, partit%npes+1, MPI_INTEGER, 0, &
                       partit%MPI_COMM_FESOM, ierror)
      end if
    end if

    if (error_status == 0) then
      allocate(part_per_node(int(grid_obj%size())))
      part_per_node = -1
      if (partit%mype == 0) then
        do r = 0, partit%npes - 1
          write(rank_string, '(I5.5)') r
          owner_file = trim(dist_mesh_dir)//'my_list'//rank_string//'.out'
          open(newunit=funit, file=trim(owner_file), action='read', status='old', &
               iostat=io_stat)
          if (io_stat /= 0) error stop 'Cannot open FESOM partition owner list'
          read(funit, *) file_rank
          read(funit, *) owned_count
          read(funit, *) ghost_count
          if (file_rank /= r) error stop 'Unexpected rank in FESOM partition owner list'
          allocate(local_node_list(owned_count + ghost_count))
          read(funit, *) local_node_list
          close(funit)
          do node_idx = 1, owned_count
            part_per_node(local_node_list(node_idx)) = r
          end do
          deallocate(local_node_list)
        end do
        if (any(part_per_node < 0)) error stop 'Incomplete FESOM node ownership map'
      end if
      call MPI_BCast(part_per_node, int(grid_obj%size()), MPI_INTEGER, 0, &
                     partit%MPI_COMM_FESOM, ierror)
      deallocate(part_csr)
      dist_obj = atlas_GridDistribution(part_per_node)
      deallocate(part_per_node)
      have_dist = .true.
      if (partit%mype == 0) write(output_unit, '(3A,I0,A)') &
        '  --> ATLAS_USE_FESOM_DIST: using partition from ', &
        trim(rpart_file), ' (', partit%npes, ' partitions)'
    else
      if (allocated(part_csr)) deallocate(part_csr)
    end if
  end subroutine set_fesom_dist

  !> @brief Convert an Atlas mesh to a FESOM t_mesh
  !! @details Extracts node coordinates and triangle connectivity, asks Atlas
  !!          to build edge topology, maps it into the FESOM t_mesh arrays,
  !!          builds level and partition arrays, and calls the
  !!          standard mesh computation routines (mesh_areas,
  !!          mesh_auxiliary_arrays).
  !! @param[inout] atlas_msh  Atlas mesh object (RegularLonLat or similar)
  !! @param[in]    zbar       Vertical interface depths
  !! @param[inout] partit     FESOM partition structure (filled for npes=1)
  !! @param[inout] mesh3      Output FESOM t_mesh (must be uninitialised)
  subroutine atlas_mesh_to_fesom_mesh(atlas_msh, zbar, partit, mesh3)
    type(atlas_Mesh),   intent(inout) :: atlas_msh
    real(kind=MP),      intent(in)    :: zbar(:)
    type(t_partit),     intent(inout) :: partit
    type(t_mesh),       intent(inout) :: mesh3

    ! Atlas helper objects
    type(atlas_mesh_Nodes)             :: atlas_nodes
    type(atlas_mesh_Cells)             :: atlas_cells
    type(atlas_mesh_Edges)             :: atlas_edges
    type(atlas_MultiBlockConnectivity) :: atlas_edge_nodes
    type(atlas_MultiBlockConnectivity) :: atlas_edge_cells
    type(atlas_Field)                  :: node_global_index_field
    type(atlas_Field)                  :: cell_global_index_field
    type(atlas_Field)                  :: ghost_field
    type(atlas_Field)                  :: node_lonlat_field
    type(atlas_Field)                  :: cell_flags_field
    type(atlas_MultiBlockConnectivity) :: atlas_cell_nodes
    type(atlas_Connectivity)           :: atlas_node_edges
    integer(ATLAS_KIND_GIDX), pointer  :: node_global_index(:)
    integer(ATLAS_KIND_GIDX), pointer  :: cell_global_index(:)
    integer(c_int), pointer            :: node_ghost(:)
    integer(c_int), pointer            :: cell_flags(:)
    real(8), pointer                   :: node_lonlat(:,:)
    integer(ATLAS_KIND_IDX), pointer   :: cell_nodes(:,:), cell_node_cols(:)
    integer(ATLAS_KIND_IDX), pointer   :: node_edges(:,:), node_edge_cols(:)
    integer(ATLAS_KIND_IDX), pointer   :: edge_nodes(:,:), edge_node_cols(:)
    integer(ATLAS_KIND_IDX), pointer   :: edge_cells(:,:), edge_cell_cols(:)

    ! Mesh dimensions
    integer :: nnodes, ncells, atlas_ncells, nl, edge2D_local, edge2D_in_local
    integer :: global_nnodes, global_ncells, global_edges, global_internal_edges
    integer :: owned_edges, owned_internal_edges, owner_node, ierror
    integer :: atlas_edge_count
    integer, allocatable :: edge_old_to_new(:)
    logical, allocatable :: node_owned(:)

    integer :: i, j
    integer :: tmp_el
    logical :: patch_seen

    ! Edge classification
    integer :: n_internal
    integer :: e, n, k, q, e1, node_a, node_b, edge_id
    integer :: old_edge, new_edge, cell1, cell2
    integer :: elnodes(3), eledges(3)
    real(kind=WP) :: geographic_lon, geographic_lat, rotated_lon, rotated_lat

    ! Owned vs ghost node split
    integer :: n_owned, n_ghost

    real(MP), parameter :: deg2rad = acos(-1.0_MP) / 180.0_MP

    ! ------------------------------------------------------------------
    ! 1. Basic Atlas mesh dimensions
    ! ------------------------------------------------------------------
    atlas_nodes = atlas_msh%nodes()
    atlas_cells = atlas_msh%cells()
    nnodes = int(atlas_nodes%size())
    atlas_ncells = int(atlas_cells%size())
    ncells = atlas_ncells
    nl = size(zbar)

    cell_flags_field = atlas_cells%field('flags')
    call cell_flags_field%data(cell_flags)
    ncells = 0
    patch_seen = .false.
    do e = 1, atlas_ncells
      if (iand(cell_flags(e), 258_c_int) /= 0) then
        patch_seen = .true.
      else
        if (patch_seen) then
          error stop 'Atlas ghost/PATCH cells are not stored after owned cells'
        end if
        ncells = ncells + 1
      end if
    end do
    call cell_flags_field%final()

    node_global_index_field = atlas_nodes%global_index()
    cell_global_index_field = atlas_cells%global_index()
    call node_global_index_field%data(node_global_index)
    call cell_global_index_field%data(cell_global_index)

    ghost_field = atlas_nodes%ghost()
    n_owned = 0
    n_ghost = 0
    allocate(node_owned(nnodes))
    call ghost_field%data(node_ghost)
    do n=1, nnodes
      if (node_ghost(n) == 0_c_int) then
        node_owned(n) = .true.
        n_owned = n_owned + 1
        if (n_ghost > 0) then
          STOP 'Atlas mesh has owned nodes after ghost nodes; fesom_mesh requires all owned nodes first'
        end if
      else
        node_owned(n) = .false.
        n_ghost = n_ghost + 1
      end if
    end do
    call ghost_field%final()
    call MPI_Allreduce(n_owned, global_nnodes, 1, MPI_INTEGER, MPI_SUM, &
                       partit%MPI_COMM_FESOM, ierror)
    call MPI_Allreduce(ncells, global_ncells, 1, MPI_INTEGER, MPI_SUM, &
                       partit%MPI_COMM_FESOM, ierror)
    if (partit%mype == 0) then
      write(output_unit, '(A,I0,A,I0,A,I0)') &
        '  atlas_mesh_to_fesom_mesh: nnodes=', nnodes, ', ncells=', ncells, &
        ', Atlas cells=', atlas_ncells
    end if

    ! ------------------------------------------------------------------
    ! 2. Atlas node coordinates: lon/lat degrees -> radians
    ! ------------------------------------------------------------------
    mesh3%nod2D = global_nnodes
    allocate(mesh3%coord_nod2D(2, nnodes))
    call set_mesh_transform_matrix()
    node_lonlat_field = atlas_nodes%lonlat()
    call node_lonlat_field%data(node_lonlat)
    do n = 1, nnodes
      geographic_lon = real(node_lonlat(1, n) * deg2rad, WP)
      geographic_lat = real(node_lonlat(2, n) * deg2rad, WP)
      if (force_rotation) then
        call g2r(geographic_lon, geographic_lat, rotated_lon, rotated_lat)
        mesh3%coord_nod2D(:, n) = real([rotated_lon, rotated_lat], MP)
      else
        mesh3%coord_nod2D(:, n) = real([geographic_lon, geographic_lat], MP)
      end if
    end do
    call node_lonlat_field%final()

    ! ------------------------------------------------------------------
    ! 3. Atlas cell connectivity
    ! ------------------------------------------------------------------
    mesh3%elem2D = global_ncells
    allocate(mesh3%elem2D_nodes(3, ncells))
    atlas_cell_nodes = atlas_cells%node_connectivity()
    call atlas_cell_nodes%padded_data(cell_nodes, cell_node_cols)
    if (any(cell_node_cols /= 3)) then
      error stop 'Atlas grid conversion requires triangular cells'
    end if
    mesh3%elem2D_nodes = int(cell_nodes(1:3, 1:ncells))

    ! ------------------------------------------------------------------
    ! 4. Vertical structure from aux3d.out, with a flat bottom
    ! ------------------------------------------------------------------
    mesh3%nl = nl
    allocate(mesh3%zbar(nl))
    mesh3%zbar = zbar
    allocate(mesh3%Z(nl-1))
    mesh3%Z = 0.5_MP * (mesh3%zbar(1:nl-1) + mesh3%zbar(2:nl))
    allocate(mesh3%depth(nnodes))
    mesh3%depth = mesh3%zbar(nl)
    allocate(mesh3%elem_depth(ncells))
    mesh3%elem_depth = mesh3%zbar(nl)

    ! ------------------------------------------------------------------
    ! 5. Build edge topology with Atlas
    ! ------------------------------------------------------------------
    call atlas_build_edges(atlas_msh)
    call atlas_build_node_to_edge_connectivity(atlas_msh)
    atlas_edges = atlas_msh%edges()
    atlas_edge_nodes = atlas_edges%node_connectivity()
    atlas_edge_cells = atlas_edges%cell_connectivity()
    call atlas_edge_nodes%padded_data(edge_nodes, edge_node_cols)
    call atlas_edge_cells%padded_data(edge_cells, edge_cell_cols)
    atlas_node_edges = atlas_nodes%edge_connectivity()
    call atlas_node_edges%padded_data(node_edges, node_edge_cols)

    atlas_edge_count = int(atlas_edges%size())
    owned_edges = 0
    owned_internal_edges = 0
    do old_edge = 1, atlas_edge_count
      if (node_global_index(edge_nodes(1, old_edge)) <= &
          node_global_index(edge_nodes(2, old_edge))) then
        owner_node = int(edge_nodes(1, old_edge))
      else
        owner_node = int(edge_nodes(2, old_edge))
      end if
      if (node_owned(owner_node)) then
        owned_edges = owned_edges + 1
        cell1 = int(edge_cells(1, old_edge))
        cell2 = int(edge_cells(2, old_edge))
        if (cell1 > 0 .and. cell1 <= atlas_ncells .and. &
            cell2 > 0 .and. cell2 <= atlas_ncells) &
          owned_internal_edges = owned_internal_edges + 1
      end if
    end do
    call MPI_Allreduce(owned_edges, global_edges, 1, MPI_INTEGER, MPI_SUM, &
                       partit%MPI_COMM_FESOM, ierror)
    call MPI_Allreduce(owned_internal_edges, global_internal_edges, 1, &
                       MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierror)
    deallocate(node_owned)

    allocate(edge_old_to_new(atlas_edge_count))
    edge_old_to_new = 0
    edge2D_local = 0
    do old_edge = 1, atlas_edge_count
      cell1 = int(edge_cells(1, old_edge))
      cell2 = 0
      if (edge_cell_cols(old_edge) == 2) cell2 = int(edge_cells(2, old_edge))
      if (cell1 > ncells) cell1 = 0
      if (cell2 > ncells) cell2 = 0
      if (cell1 > 0 .and. cell2 > 0) then
        edge2D_local = edge2D_local + 1
        edge_old_to_new(old_edge) = edge2D_local
      end if
    end do
    edge2D_in_local = edge2D_local
    do old_edge = 1, atlas_edge_count
      cell1 = int(edge_cells(1, old_edge))
      cell2 = 0
      if (edge_cell_cols(old_edge) == 2) cell2 = int(edge_cells(2, old_edge))
      if (cell1 > ncells) cell1 = 0
      if (cell2 > ncells) cell2 = 0
      if ((cell1 > 0) .neqv. (cell2 > 0)) then
        edge2D_local = edge2D_local + 1
        edge_old_to_new(old_edge) = edge2D_local
      end if
    end do

    ! Preserve Atlas-local edge ordering and orientation.
    allocate(mesh3%edges(2, edge2D_local))
    allocate(mesh3%edge_tri(2, edge2D_local))
    allocate(mesh3%elem_edges(3, ncells))
    mesh3%elem_edges = 0

    do old_edge = 1, atlas_edge_count
      new_edge = edge_old_to_new(old_edge)
      if (new_edge == 0) cycle
      mesh3%edges(:, new_edge) = int(edge_nodes(:, old_edge))
      mesh3%edge_tri(1, new_edge) = int(edge_cells(1, old_edge))
      if (edge_cell_cols(old_edge) == 2) then
        mesh3%edge_tri(2, new_edge) = int(edge_cells(2, old_edge))
      else
        mesh3%edge_tri(2, new_edge) = 0
      end if
      if (mesh3%edge_tri(1, new_edge) > ncells) mesh3%edge_tri(1, new_edge) = 0
      if (mesh3%edge_tri(2, new_edge) > ncells) mesh3%edge_tri(2, new_edge) = 0
      if (mesh3%edge_tri(1, new_edge) == 0 .and. &
          mesh3%edge_tri(2, new_edge) > 0) then
        mesh3%edge_tri(1, new_edge) = mesh3%edge_tri(2, new_edge)
        mesh3%edge_tri(2, new_edge) = 0
      end if
    end do
    n_internal = edge2D_in_local

    mesh3%edge2D    = global_edges
    mesh3%edge2D_in = global_internal_edges

    ! ------------------------------------------------------------------
    ! 6. Build elem_edges: elem_edges(q,e) = edge opposite to node q
    ! ------------------------------------------------------------------
    ! Find the edge opposite each cell node through node-to-edge adjacency.
    do e = 1, ncells
      elnodes = mesh3%elem2D_nodes(:, e)
      do q = 1, 3
        node_a = elnodes(mod(q, 3) + 1)
        node_b = elnodes(mod(q + 1, 3) + 1)
        do k = 1, int(node_edge_cols(node_a))
          old_edge = int(node_edges(k, node_a))
          edge_id = edge_old_to_new(old_edge)
          if (edge_id == 0) cycle
          if ((mesh3%edges(1, edge_id) == node_a .and. &
               mesh3%edges(2, edge_id) == node_b) .or. &
              (mesh3%edges(1, edge_id) == node_b .and. &
               mesh3%edges(2, edge_id) == node_a)) then
            mesh3%elem_edges(q, e) = edge_id
            exit
          end if
        end do
        if (mesh3%elem_edges(q, e) == 0) then
          error stop 'Cannot find Atlas edge for triangular cell'
        end if
      end do
    end do
    deallocate(edge_old_to_new)

    ! ------------------------------------------------------------------
    ! 7. Partition arrays: Atlas stores owned nodes before ghost nodes.
    ! ------------------------------------------------------------------
    partit%myDim_nod2D  = n_owned;      partit%eDim_nod2D   = n_ghost
    partit%myDim_elem2D = ncells;       partit%eDim_elem2D  = 0
    partit%eXDim_elem2D = 0
    partit%myDim_edge2D = edge2D_local; partit%eDim_edge2D  = 0

    allocate(partit%myList_nod2D(n_owned))
    do n = 1, n_owned
      partit%myList_nod2D(n) = int(node_global_index(n))
    end do
    allocate(partit%myList_elem2D(ncells))
    do n = 1, ncells
      partit%myList_elem2D(n) = int(cell_global_index(n))
    end do
    allocate(partit%myList_edge2D(edge2D_local))
    do n = 1, edge2D_local
      partit%myList_edge2D(n) = n
    end do
    ! CSR partition vector from the actual owned-node count on each rank.
    allocate(partit%part(partit%npes+1))
    partit%part(1) = 1
    call MPI_Allgather(n_owned, 1, MPI_INTEGER, partit%part(2), 1, MPI_INTEGER, &
                       partit%MPI_COMM_FESOM, ierror)
    do n = 2, partit%npes + 1
      partit%part(n) = partit%part(n) + partit%part(n-1)
    end do
    if (partit%part(partit%npes+1) /= mesh3%nod2D + 1) then
      error stop 'Atlas node partition vector does not match global node count'
    end if

    ! Empty halo communication structures (npes=1: no halo exchange)
    partit%com_nod2D%rPEnum  = 0; partit%com_nod2D%sPEnum  = 0
    partit%com_elem2D%rPEnum = 0; partit%com_elem2D%sPEnum = 0
    partit%com_elem2D_full%rPEnum = 0; partit%com_elem2D_full%sPEnum = 0
    allocate(partit%com_nod2D%rlist(0))
    allocate(partit%com_nod2D%slist(0))
    allocate(partit%com_elem2D%rlist(0))
    allocate(partit%com_elem2D%slist(0))
    allocate(partit%com_elem2D_full%rlist(0))
    allocate(partit%com_elem2D_full%slist(0))
    allocate(partit%myInd_elem2D_shrinked(ncells))
    do n = 1, ncells; partit%myInd_elem2D_shrinked(n) = n; end do
    partit%myDim_elem2D_shrinked = ncells

    ! ------------------------------------------------------------------
    ! 8. Level arrays
    ! ------------------------------------------------------------------
    allocate(mesh3%bc_index_nod2D(nnodes));   mesh3%bc_index_nod2D = 0

    call find_levels(partit, mesh3)

    ! ------------------------------------------------------------------
    ! 9. elem_neighbors and nod_in_elem2D
    ! ------------------------------------------------------------------
    allocate(mesh3%elem_neighbors(3, ncells))
    mesh3%elem_neighbors = 0
    do e = 1, ncells
      eledges = mesh3%elem_edges(:, e)
      do k = 1, 3
        e1 = mesh3%edge_tri(1, eledges(k))
        if (e1 == e) e1 = mesh3%edge_tri(2, eledges(k))
        mesh3%elem_neighbors(k, e) = e1
      end do
    end do

    allocate(mesh3%nod_in_elem2D_num(nnodes))
    mesh3%nod_in_elem2D_num = 0
    do e = 1, ncells
      do k = 1, 3
        n = mesh3%elem2D_nodes(k, e)
        mesh3%nod_in_elem2D_num(n) = mesh3%nod_in_elem2D_num(n) + 1
      end do
    end do
    allocate(mesh3%nod_in_elem2D(maxval(mesh3%nod_in_elem2D_num), nnodes))
    mesh3%nod_in_elem2D = 0
    mesh3%nod_in_elem2D_num = 0
    do e = 1, ncells
      do k = 1, 3
        n = mesh3%elem2D_nodes(k, e)
        mesh3%nod_in_elem2D_num(n) = mesh3%nod_in_elem2D_num(n) + 1
        mesh3%nod_in_elem2D(mesh3%nod_in_elem2D_num(n), n) = e
      end do
    end do

    do n = 1, nnodes
      do i = 2, mesh3%nod_in_elem2D_num(n)
        tmp_el = mesh3%nod_in_elem2D(i, n)
        j = i - 1
        do while (j >= 1 .and. &
                  cell_global_index(mesh3%nod_in_elem2D(j, n)) > cell_global_index(tmp_el))
          mesh3%nod_in_elem2D(j+1, n) = mesh3%nod_in_elem2D(j, n)
          j = j - 1
        end do
        mesh3%nod_in_elem2D(j+1, n) = tmp_el
      end do
    end do

    call node_global_index_field%final()
    call cell_global_index_field%final()

    ! ------------------------------------------------------------------
    ! 10. Mesh computation routines
    ! ------------------------------------------------------------------
    call init_mpi_types(partit, mesh3)
    call init_gatherLists(partit)
    call test_tri(partit, mesh3)
    call find_levels_min_e2n(partit, mesh3)
    call mesh_areas(partit, mesh3)
    call mesh_auxiliary_arrays(partit, mesh3)

    if (partit%mype == 0) then
      write(output_unit, '(A,I0,A,I0,A,I0)') &
        '  atlas_mesh_to_fesom_mesh complete: ', mesh3%nod2D, ' nodes, ', &
        mesh3%elem2D, ' cells, ', mesh3%edge2D, ' edges'
    end if

  end subroutine atlas_mesh_to_fesom_mesh
#endif

  !> @brief Compute tracer min/max/sum using Atlas NodeColumns
  !! @details Atlas mesh overlap nodes are marked as ghosts during conversion,
  !!          so NodeColumns reductions include each global node exactly once.
  !! @param[in] tracer_data 2D tracer array (nz, npoints)
  !! @param[in] n_owned Number of owned nodes on this rank
  !! @param[in] partit Partition structure (for MPI communicator)
  !! @param[out] tmin Minimum value (real(8))
  !! @param[out] tmax Maximum value (real(8))
  !! @param[out] tsum Sum of all values (real(8))
  subroutine compute_tracer_stats_atlas(tracer_data, tmin, tmax, tsum)
    real(kind=WP), intent(in) :: tracer_data(:,:)
    real(8), intent(out) :: tmin, tmax, tsum

#ifdef ENABLE_ATLAS
    type(atlas_Field) :: tracer_field
    real(WP), pointer :: atlas_tracer(:,:)
    real(WP) :: tmin_wp, tmax_wp, tsum_wp

    if (.not. atlas_fesom_active()) then
      stop 'compute_tracer_stats_atlas: Atlas FESOM not active'
    end if

    tracer_field = atlas_nodes_global%create_field(name='tracer', &
                             kind=atlas_real(WP), &
                             levels=size(tracer_data, 1))
    call tracer_field%data(atlas_tracer)
    if (size(atlas_tracer, 1) /= size(tracer_data, 1) .or. &
        size(atlas_tracer, 2) /= size(tracer_data, 2)) then
      error stop 'compute_tracer_stats_atlas: tracer shape does not match Atlas mesh'
    end if
    atlas_tracer = tracer_data

    call atlas_nodes_global%minimum(tracer_field, tmin_wp)
    call atlas_nodes_global%maximum(tracer_field, tmax_wp)
    call atlas_nodes_global%order_independent_sum(tracer_field, tsum_wp)
    tmin = dble(tmin_wp)
    tmax = dble(tmax_wp)
    tsum = dble(tsum_wp)

    call tracer_field%final()
#endif
  end subroutine compute_tracer_stats_atlas

end module atlas_fesom_mesh_module
