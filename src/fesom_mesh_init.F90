!> @file fesom_mesh_init.F90
!! @brief Tracer dwarf main program using mesh file initialization
!! @details This version initializes from mesh files (not restart files)
!!          allowing custom tracer values to be set before advection.

program tracer_dwarf_mesh_init
  use mpi
  use MOD_MESH
  use MOD_PARTIT
  use MOD_DYN
  use MOD_TRACER
  use MOD_PARSUP
  use oce_adv_tra_driver_module
  use tracer_init_from_mesh_module
  use atlas_fesom_mesh_module, only: atlas_fesom_enabled, compute_tracer_stats_atlas
  use g_config
  use o_PARAM
#ifdef ENABLE_ATLAS
  use atlas_module, only: atlas_initialize, atlas_finalize
#endif
  
  implicit none
  
  type(t_partit) :: partit
  type(t_mesh)   :: mesh
  type(t_dyn)    :: dynamics
  type(t_tracer) :: tracers
  
  integer :: n, nz, istep, nsteps, ivar,ierr, elem, kz, nzmin, nzmax
  integer :: elem_nodes(3)
  real(kind=WP) :: dt_local, x_coord, y_coord, z_coord
  real(8), dimension(2) :: trmin_loc, trmax_loc, trsum_loc, trmin, trmax, trsum
  
  ! ========================================
  ! Initialize MPI first
  ! ========================================
  call MPI_Init(ierr)
#ifdef ENABLE_ATLAS
  if (atlas_fesom_enabled()) call atlas_initialize()
#endif
  partit%MPI_COMM_FESOM = MPI_COMM_WORLD
  
  ! Then call par_init to set up partition info
  call par_init(partit)
  
  if (partit%mype == 0) then
    write(*, '(A)') ''
    write(*, '(A)') '================================================'
    write(*, '(A)') 'FESOM2 Tracer Advection Dwarf'
    write(*, '(A)') 'Mesh File Initialization Version'
    write(*, '(A)') '================================================'
    write(*, '(A,I6)') 'Running on ', partit%npes, ' processes'
    write(*, '(A)') '================================================'
    write(*, '(A)') ''
    write(*, '(A,I2,A)') 'Working Precision (WP) = ', WP, ' bytes'
    if (WP == 2) then
      write(*, '(A)') '*** HALF PRECISION MODE (FP16) ***'
      write(*, '(A)') '*** EXPERIMENTAL - GPU ONLY ***'
    else if (WP == 4) then
      write(*, '(A)') '*** SINGLE PRECISION MODE ***'
    else if (WP == 8) then
      write(*, '(A)') '*** DOUBLE PRECISION MODE ***'
    end if
    write(*, '(A)') ''
  end if
  
  ! ========================================
  ! Set basic configuration
  ! ========================================
  if (partit%mype == 0) write(*, '(A)') 'Setting configuration...'
  
  ! Set mesh path - adjust this to point to your mesh
  MeshPath = './tests/data/pi/'
  
  ! Set basic parameters
  use_ice = .false.
  use_sw_pene = .false.
  use_cavity = .false.
  flag_debug = .false.
  
  ! ========================================
  ! Initialize mesh and arrays from mesh files
  ! ========================================
  call tracer_init_mesh_and_arrays(partit, mesh, dynamics, tracers)

  ! Use a smooth, spatially varying horizontal velocity field.
  ! Need to change zero-velocity field to a spatially varying field as advection would not be properly tested in using the halo exchange
  nz = mesh%nl - 1
  do elem = 1, partit%myDim_elem2D
    elem_nodes = mesh%elem2D_nodes(:, elem)
    x_coord = sum(real(mesh%coord_nod2D(1, elem_nodes), WP)) / 3.0_WP
    y_coord = sum(real(mesh%coord_nod2D(2, elem_nodes), WP)) / 3.0_WP
    do kz = 1, nz
      z_coord = real(kz-1, WP) / real(nz-1, WP)
      dynamics%uv(1, kz, elem) = 0.1_WP * &
        (1.0_WP + 0.25_WP*sin(x_coord)*cos(y_coord)) * (1.0_WP - 0.2_WP*z_coord)
      dynamics%uv(2, kz, elem) = 0.05_WP * &
        cos(x_coord)*sin(y_coord) * (1.0_WP - 0.2_WP*z_coord)
    end do
  end do
  
  ! ========================================
  ! Set custom tracer values
  ! ========================================
  if (partit%mype == 0) then
    write(*, '(A)') '================================================'
    write(*, '(A)') 'Setting custom tracer values'
    write(*, '(A)') '================================================'
  end if
  
  ! Set tracers to affine functions of horizontal position and depth.
  do n = 1, partit%myDim_nod2D
    x_coord = real(mesh%coord_nod2D(1, n), WP)
    y_coord = real(mesh%coord_nod2D(2, n), WP)
    do istep = 1, nz
      z_coord = real(istep-1, WP) / real(nz-1, WP)
      tracers%data(1)%values(istep, n) = 15.0_WP + 2.0_WP*x_coord &
                                      - 1.0_WP*y_coord - 10.0_WP*z_coord
      tracers%data(2)%values(istep, n) = 35.0_WP - 0.5_WP*x_coord &
                                      + 0.25_WP*y_coord + 1.0_WP*z_coord
    end do
  end do
  
  ! Initialize AB and old values
  do n = 1, tracers%num_tracers
    tracers%data(n)%valuesAB = tracers%data(n)%values
    do istep = 1, tracers%data(n)%AB_order
      tracers%data(n)%valuesold(:, :, istep) = tracers%data(n)%values
    end do
  end do
  
  if (partit%mype == 0) then
    write(*, '(A)') '  Temperature: linear in x, y, and depth'
    write(*, '(A)') '  Salinity: linear in x, y, and depth'
    write(*, '(A)') ''
  end if
  
  ! ========================================
  ! Print initial statistics
  ! ========================================
  if (partit%mype == 0) then
    write(*, '(A)') '================================================'
    write(*, '(A)') 'Initial tracer statistics'
    write(*, '(A)') '================================================'
  end if
  
  ! Compute tracer statistics using Atlas (if available) or MPI_Allreduce
  if (atlas_fesom_enabled()) then
    call compute_tracer_stats_atlas(tracers%data(1)%values, trmin(1), trmax(1), trsum(1))
    call compute_tracer_stats_atlas(tracers%data(2)%values, trmin(2), trmax(2), trsum(2))
  else
    do ivar = 1, 2
      trmin_loc(ivar) = dble(minval(tracers%data(ivar)%values(:, 1:partit%myDim_nod2D)))
      trmax_loc(ivar) = dble(maxval(tracers%data(ivar)%values(:, 1:partit%myDim_nod2D)))
      trsum_loc(ivar) = sum(dble(tracers%data(ivar)%values(:, 1:partit%myDim_nod2D)))
    end do
    call MPI_Allreduce(trmin_loc, trmin, 2, MPI_DOUBLE_PRECISION, MPI_MIN, &
                      partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(trmax_loc, trmax, 2, MPI_DOUBLE_PRECISION, MPI_MAX, &
                      partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(trsum_loc, trsum, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      partit%MPI_COMM_FESOM, ierr)
  end if


  if (partit%mype == 0) then
    write(*, '(A,3E18.10)') '  Temperature: min, max, sum = ', trmin(1), trmax(1), trsum(1)
    write(*, '(A,3E18.10)') '  Salinity:    min, max, sum = ', trmin(2), trmax(2), trsum(2)
    write(*, '(A)') ''
  end if
  
  ! ========================================
  ! Run advection test
  ! ========================================
  if (partit%mype == 0) then
    write(*, '(A)') '================================================'
    write(*, '(A)') 'Running advection test'
    write(*, '(A)') '================================================'
  end if
  
  nsteps = 10
  dt_local = 360.0_WP
  
  if (partit%mype == 0) then
    write(*, '(A,I6,A,E12.4)') '  Steps: ', nsteps, ', dt = ', dt_local
    write(*, '(A)') ''
  end if
  
  do istep = 1, nsteps
    ! Advect each tracer
    do n = 1, tracers%num_tracers
      tracers%work%del_ttf = 0.0_MP
      tracers%work%del_ttf_advhoriz = 0.0_MP
      tracers%work%del_ttf_advvert = 0.0_MP
      tracers%data(n)%valuesAB = tracers%data(n)%values

      call do_oce_adv_tra(dt_local, dynamics%uv, dynamics%w, dynamics%w, dynamics%w, n, &
                          dynamics, tracers, partit, mesh)

      tracers%work%del_ttf = tracers%work%del_ttf_advhoriz &
                            + tracers%work%del_ttf_advvert
      do kz = 1, partit%myDim_nod2D
        nzmin = mesh%ulevels_nod2D(kz)
        nzmax = mesh%nlevels_nod2D(kz) - 1
        tracers%data(n)%values(nzmin:nzmax, kz) = &
          tracers%data(n)%values(nzmin:nzmax, kz) &
          + tracers%work%del_ttf(nzmin:nzmax, kz) &
          / mesh%hnode_new(nzmin:nzmax, kz)
      end do
    end do
    
    ! Print statistics every step (computed globally across all ranks)
    if (atlas_fesom_enabled()) then
      call compute_tracer_stats_atlas(tracers%data(1)%values, trmin(1), trmax(1), trsum(1))
    else
      ! Compute local min/max/sum on this rank
      trmin_loc(1) = dble(minval(tracers%data(1)%values(:, 1:partit%myDim_nod2D)))
      trmax_loc(1) = dble(maxval(tracers%data(1)%values(:, 1:partit%myDim_nod2D)))
      trsum_loc(1) = sum(dble(tracers%data(1)%values(:, 1:partit%myDim_nod2D)))

      ! Reduce across all ranks
      call MPI_Allreduce(trmin_loc(1), trmin(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                        partit%MPI_COMM_FESOM, ierr)
      call MPI_Allreduce(trmax_loc(1), trmax(1), 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                        partit%MPI_COMM_FESOM, ierr)
      call MPI_Allreduce(trsum_loc(1), trsum(1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        partit%MPI_COMM_FESOM, ierr)
    end if
    if (partit%mype == 0) then
      write(*, '(A,I4,A,3E18.10)') '  Step ', istep, ': T min, max, sum = ', trmin(1), trmax(1), trsum(1)
    end if
  end do
  
  if (partit%mype == 0) then
    write(*, '(A)') ''
    write(*, '(A)') '================================================'
    write(*, '(A)') 'Advection test complete!'
    write(*, '(A)') '================================================'
    write(*, '(A)') ''
  end if
  
  ! ========================================
  ! Finalize MPI
  ! ========================================
#ifdef ENABLE_ATLAS
  if (atlas_fesom_enabled()) call atlas_finalize()
#endif
  call par_ex(partit%MPI_COMM_FESOM, partit%mype)
  
  if (partit%mype == 0) then
    write(*, '(A)') 'Program finished successfully'
  end if

end program tracer_dwarf_mesh_init
