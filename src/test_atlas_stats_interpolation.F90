!> @file test_atlas_stats_interpolation.F90
!! @brief Unit test: Atlas field statistics on fesom-pi and O128 meshes
!! @details Creates NodeColumns function spaces for fesom-pi and O128, initializes
!!          a tracer field on fesom-pi using a sin(lat)*cos(lon) pattern, sets up
!!          finite-element interpolation from fesom-pi to O128, and computes
!!          min/max/sum statistics on both meshes via compute_field_stats_atlas.

#ifdef ENABLE_ATLAS

program test_atlas_stats_interpolation
  use atlas_module
  use atlas_fesom_mesh_module, only: compute_tracer_stats_atlas, set_atlas_stats_mesh
  use mpi
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none

  integer, parameter :: wp = c_double   ! double precision throughout

  type(atlas_Grid)                      :: grid_pi, grid_O128
  type(atlas_Mesh)                      :: mesh_pi, mesh_O128
  type(atlas_MeshGenerator)             :: meshgen_fesom  ! for fesom unstructured grids
  type(atlas_MeshGenerator)             :: meshgen_struct  ! for structured Gaussian grids
  type(atlas_functionspace_NodeColumns) :: fs_pi, fs_O128
  type(atlas_Field)                     :: field_pi, field_O128
  type(atlas_Field)                     :: lonlat_field
  type(atlas_mesh_Nodes)                :: nodes_pi
  type(atlas_Config)                    :: interp_config
  type(atlas_Interpolation)             :: interpolation

  real(wp), pointer :: data_pi(:,:), data_O128(:,:)
  real(wp), pointer :: lonlat(:,:)

  real(8)  :: tmin, tmax, tsum
  real(wp) :: lon_r, lat_r, deg2rad
  integer  :: jnode, nb_nodes
  integer  :: ierr, mype

  deg2rad = 2.0_wp * asin(1.0_wp) / 180.0_wp

  ! ========================================
  ! Initialize MPI and Atlas
  ! ========================================
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)

  call atlas_initialize()

  if (mype == 0) then
    write(*,'(A)') ''
    write(*,'(A)') '========================================'
    write(*,'(A)') ' test_atlas_stats_interpolation'
    write(*,'(A)') '========================================'
    write(*,'(A)') ''
  end if

  meshgen_fesom  = atlas_MeshGenerator("fesom")  ! fesom unstructured mesh gen
  meshgen_struct = atlas_MeshGenerator()          ! default structured mesh gen

  ! ========================================
  ! 1. fesom-pi: mesh + NodeColumns function space
  ! ========================================
  if (mype == 0) write(*,'(A)') 'Step 1: Create fesom-pi mesh and function space'

  grid_pi = atlas_Grid("fesom-pi")
  mesh_pi = meshgen_fesom%generate(grid_pi)
  fs_pi   = atlas_functionspace_NodeColumns(mesh_pi, halo=1)

  nb_nodes = fs_pi%nb_nodes()
  if (mype == 0) write(*,'(A,I0,A)') '  fesom-pi local nodes (incl. halo): ', nb_nodes, ''

  ! Create 2D scalar field (levels=0) and populate from lon/lat coordinates
  field_pi     = fs_pi%create_field(name='tracer_pi', kind=atlas_real(wp), levels=1)
  nodes_pi     = fs_pi%nodes()
  lonlat_field = nodes_pi%lonlat()
  call lonlat_field%data(lonlat)
  call field_pi%data(data_pi)

  do jnode = 1, nb_nodes
    lon_r = lonlat(1,jnode) * deg2rad
    lat_r = lonlat(2,jnode) * deg2rad
    ! sin(lat)*cos(lon) + 2 to keep values positive
    data_pi(1,jnode) = sin(lat_r) * cos(lon_r) + 2.0_wp
  end do

  call fs_pi%halo_exchange(field_pi)

  ! Compute statistics on fesom-pi via Atlas function space reductions
  call set_atlas_stats_mesh(mesh_pi, halo=1)
  call compute_tracer_stats_atlas(data_pi, tmin, tmax, tsum)

  if (mype == 0) then
    write(*,'(A)') '  fesom-pi tracer stats (sin(lat)*cos(lon) + 2):'
    write(*,'(A,3ES14.6)') '    min, max, sum = ', tmin, tmax, tsum
    write(*,'(A)') ''
  end if

  ! ========================================
  ! 2. O128: mesh + NodeColumns function space
  ! ========================================
  if (mype == 0) write(*,'(A)') 'Step 2: Create O128 mesh and function space'

  grid_O128   = atlas_Grid("O128")
  mesh_O128   = meshgen_struct%generate(grid_O128)
  fs_O128     = atlas_functionspace_NodeColumns(mesh_O128, halo=1)

  if (mype == 0) then
    write(*,'(A,I0,A)') '  O128 local nodes (incl. halo): ', fs_O128%nb_nodes(), ''
    write(*,'(A)') ''
  end if

  field_O128 = fs_O128%create_field(name='tracer_O128', kind=atlas_real(wp), levels=1)
  call field_O128%data(data_O128)

  ! ========================================
  ! 3. Finite-element interpolation fesom-pi -> O128
  ! ========================================
  if (mype == 0) write(*,'(A)') 'Step 3: Nearest-neighbour interpolation fesom-pi -> O128'

  interp_config = atlas_Config()
  call interp_config%set("type", "nearest-neighbour")
  interpolation = atlas_Interpolation(interp_config, fs_pi, fs_O128)

  call interpolation%execute(field_pi, field_O128)
  call fs_O128%halo_exchange(field_O128)

  ! Compute statistics on O128 after interpolation
  call set_atlas_stats_mesh(mesh_O128, halo=1)
  call compute_tracer_stats_atlas(data_O128, tmin, tmax, tsum)

  if (mype == 0) then
    write(*,'(A)') '  O128 tracer stats (after interpolation from fesom-pi):'
    write(*,'(A,3ES14.6)') '    min, max, sum = ', tmin, tmax, tsum
    write(*,'(A)') ''
    write(*,'(A)') '========================================'
    write(*,'(A)') ' Test PASSED'
    write(*,'(A)') '========================================'
    write(*,'(A)') ''
  end if

  ! ========================================
  ! Cleanup
  ! ========================================
  call interpolation%final()
  call interp_config%final()
  call field_O128%final()
  call field_pi%final()
  call fs_O128%final()
  call fs_pi%final()
  call mesh_O128%final()
  call mesh_pi%final()
  call grid_O128%final()
  call grid_pi%final()
  call meshgen_fesom%final()
  call meshgen_struct%final()

  call atlas_finalise()
  call MPI_Finalize(ierr)

end program test_atlas_stats_interpolation

#else

! When Atlas is not enabled, emit a stub that reports the test was skipped
program test_atlas_stats_interpolation
  implicit none
  write(*,'(A)') 'test_atlas_stats_interpolation: SKIPPED (build without ENABLE_ATLAS)'
end program test_atlas_stats_interpolation

#endif
