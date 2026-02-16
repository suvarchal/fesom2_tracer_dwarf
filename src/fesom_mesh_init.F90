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
  use oce_edge_coloring_module, only: compute_edge_coloring
  use g_config
  use o_PARAM
  
  implicit none
  
  type(t_partit) :: partit
  type(t_mesh)   :: mesh
  type(t_dyn)    :: dynamics
  type(t_tracer) :: tracers
  
  integer :: n, nz, istep, nsteps, ierr
  real(kind=WP) :: dt_local
  real(kind=WP) :: tmin, tmax, tsum, smin, smax, ssum
  
  ! ========================================
  ! Initialize MPI first
  ! ========================================
  call MPI_Init(ierr)
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

  ! Compute edge coloring for lock-free OpenMP edge-to-node scatter
  if (partit%mype == 0) write(*, '(A)') '  --> Computing edge coloring'
  call compute_edge_coloring(partit%edge_coloring, partit%myDim_edge2D, &
                             partit%myDim_nod2D, partit%eDim_nod2D, mesh%edges)

  ! ========================================
  ! Set custom tracer values
  ! ========================================
  if (partit%mype == 0) then
    write(*, '(A)') '================================================'
    write(*, '(A)') 'Setting custom tracer values'
    write(*, '(A)') '================================================'
  end if
  
  nz = mesh%nl - 1
  
  ! Set temperature: warm at surface, cold at depth
  do n = 1, partit%myDim_nod2D
    do istep = 1, nz
      ! Simple vertical gradient: 20°C at surface, 5°C at bottom
      tracers%data(1)%values(istep, n) = 20.0_WP - 15.0_WP * real(istep-1, WP) / real(nz-1, WP)
    end do
  end do
  
  ! Set salinity: constant 35 PSU
  tracers%data(2)%values(:, 1:partit%myDim_nod2D) = 35.0_WP
  
  ! Initialize AB and old values
  do n = 1, tracers%num_tracers
    tracers%data(n)%valuesAB = tracers%data(n)%values
    do istep = 1, tracers%data(n)%AB_order
      tracers%data(n)%valuesold(:, :, istep) = tracers%data(n)%values
    end do
  end do
  
  if (partit%mype == 0) then
    write(*, '(A)') '  Temperature: 20°C (surface) -> 5°C (bottom)'
    write(*, '(A)') '  Salinity: 35 PSU (constant)'
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
  
  tmin = minval(tracers%data(1)%values(:, 1:partit%myDim_nod2D))
  tmax = maxval(tracers%data(1)%values(:, 1:partit%myDim_nod2D))
  tsum = sum(tracers%data(1)%values(:, 1:partit%myDim_nod2D))
  
  smin = minval(tracers%data(2)%values(:, 1:partit%myDim_nod2D))
  smax = maxval(tracers%data(2)%values(:, 1:partit%myDim_nod2D))
  ssum = sum(tracers%data(2)%values(:, 1:partit%myDim_nod2D))
  
  if (partit%mype == 0) then
    write(*, '(A,3E14.6)') '  Temperature: min, max, sum = ', tmin, tmax, tsum
    write(*, '(A,3E14.6)') '  Salinity:    min, max, sum = ', smin, smax, ssum
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
  dt_local = 1.0e-3_WP
  
  if (partit%mype == 0) then
    write(*, '(A,I6,A,E12.4)') '  Steps: ', nsteps, ', dt = ', dt_local
    write(*, '(A)') ''
  end if
  
  do istep = 1, nsteps
    ! Advect each tracer
    do n = 1, tracers%num_tracers
      call do_oce_adv_tra(dt_local, dynamics%uv, dynamics%w, dynamics%w, dynamics%w, n, &
                          dynamics, tracers, partit, mesh)
    end do
    
    ! Print statistics every step
    if (partit%mype == 0) then
      tmin = minval(tracers%data(1)%values(:, 1:partit%myDim_nod2D))
      tmax = maxval(tracers%data(1)%values(:, 1:partit%myDim_nod2D))
      tsum = sum(tracers%data(1)%values(:, 1:partit%myDim_nod2D))
      
      write(*, '(A,I4,A,3E14.6)') '  Step ', istep, ': T min, max, sum = ', tmin, tmax, tsum
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
  call par_ex(partit%MPI_COMM_FESOM, partit%mype)
  
  if (partit%mype == 0) then
    write(*, '(A)') 'Program finished successfully'
  end if

end program tracer_dwarf_mesh_init
