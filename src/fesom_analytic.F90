!> @file fesom_analytic.F90
!! @brief Tracer dwarf main program using analytic mesh generation
!! @details This version generates a rectangular triangular mesh entirely
!!          in-memory (no file I/O), enabling single-process (npes=1) testing
!!          with arbitrary mesh sizes.
!!
!! Usage:
!!   mpirun -np 1 fesom_tracer_analytic [nx] [ny] [nl] [--save-mesh] [--save-scalars] [--periodic]
!!
!! Arguments (all optional, positional):
!!   nx  - Number of grid points in x direction (default: 20)
!!   ny  - Number of grid points in y direction (default: 20)
!!   nl  - Number of vertical levels (default: 10, must be >= 3)
!!
!! Options:
!!   --save-mesh     Write mesh arrays as raw binary files for Python
!!                   visualization (see plot_mesh.py)
!!   --save-scalars  Write tracer fields (temperature, salinity) at every
!!                   advection step (see plot_scalars.py)
!!   --periodic      Use doubly-periodic boundary conditions (default: closed-wall)
!!   --output-dir DIR  Output directory for mesh/scalar files (default: mesh_output)
!!
!! Examples:
!!   mpirun -np 1 fesom_tracer_analytic              # 400 nodes (20x20)
!!   mpirun -np 1 fesom_tracer_analytic 40 40        # 1600 nodes (40x40)
!!   mpirun -np 1 fesom_tracer_analytic 40 40 20     # 1600 nodes, 20 levels
!!   mpirun -np 1 fesom_tracer_analytic 20 20 --save-mesh  # save mesh output
!!   mpirun -np 1 fesom_tracer_analytic 20 20 --save-scalars  # save tracer output
!!   mpirun -np 1 fesom_tracer_analytic 20 20 --periodic  # periodic boundaries
!!   mpirun -np 1 fesom_tracer_analytic 20 20 --save-scalars --output-dir output_sp
!!
!! Domain size (Lx=Ly=100km) and depth (1000m) are fixed.

program tracer_dwarf_analytic
  use mpi
  use MOD_MESH
  use MOD_PARTIT
  use MOD_DYN
  use MOD_TRACER
  use MOD_PARSUP
  use oce_adv_tra_driver_module
  use oce_adv_tra_fct_module
  use analytic_mesh_module
  use mesh_output_module
  use oce_node_edge_map_module
  use g_config
  use o_PARAM
  use fesom_profiler
#ifdef USE_HALF_PRECISION
  use hp_math_intrinsics
#endif

  implicit none

  type(t_partit) :: partit
  type(t_mesh)   :: mesh
  type(t_dyn)    :: dynamics
  type(t_tracer) :: tracers

  integer :: n, nz, istep, nsteps, ierr
  real(kind=WP) :: dt_local
  real(8) :: tmin, tmax, tsum, smin, smax, ssum
  integer :: nan_count
  real(kind=MP) :: xn, yn, xc, yc, r
  real(kind=WP) :: zn
  integer :: nzmin, nzmax, k, kz

  ! Mesh parameters with defaults
  integer :: nx, ny, nl
  real(kind=MP), parameter :: Lx = 100000.0_MP   ! 100 km
  real(kind=MP), parameter :: Ly = 100000.0_MP   ! 100 km
  real(kind=MP), parameter :: max_depth = 1000.0_MP  ! 1000 m

  ! Command-line parsing
  integer :: nargs, iarg, io_stat, pos_count
  character(len=256) :: arg_str
  character(len=256) :: output_dir
  logical :: save_mesh
  logical :: save_scalars
  logical :: periodic
  logical :: next_is_output_dir
  integer :: scalar_units(2)

  ! Defaults
  nx = 20
  ny = 20
  nl = 10
  save_mesh = .false.
  save_scalars = .false.
  periodic = .false.
  output_dir = 'mesh_output'
  next_is_output_dir = .false.

  ! Parse command-line arguments: [nx] [ny] [nl] [--save-mesh] [--save-scalars]
  ! Flags can appear in any position
  nargs = command_argument_count()
  pos_count = 0
  do iarg = 1, nargs
    call get_command_argument(iarg, arg_str)
    if (next_is_output_dir) then
      output_dir = trim(arg_str)
      next_is_output_dir = .false.
      cycle
    end if
    if (trim(arg_str) == '--save-mesh') then
      save_mesh = .true.
      cycle
    end if
    if (trim(arg_str) == '--save-scalars') then
      save_scalars = .true.
      cycle
    end if
    if (trim(arg_str) == '--periodic') then
      periodic = .true.
      cycle
    end if
    if (trim(arg_str) == '--output-dir') then
      next_is_output_dir = .true.
      cycle
    end if
    pos_count = pos_count + 1
    select case (pos_count)
    case (1)
      read(arg_str, *, iostat=io_stat) nx
      if (io_stat /= 0) then
        write(*, '(A,A)') 'Error: invalid nx argument: ', trim(arg_str)
        stop 1
      end if
    case (2)
      read(arg_str, *, iostat=io_stat) ny
      if (io_stat /= 0) then
        write(*, '(A,A)') 'Error: invalid ny argument: ', trim(arg_str)
        stop 1
      end if
    case (3)
      read(arg_str, *, iostat=io_stat) nl
      if (io_stat /= 0) then
        write(*, '(A,A)') 'Error: invalid nl argument: ', trim(arg_str)
        stop 1
      end if
    case default
      write(*, '(A,A)') 'Error: unexpected argument: ', trim(arg_str)
      stop 1
    end select
  end do

  ! Validate parameters
  if (nx < 2 .or. ny < 2) then
    write(*, '(A)') 'Error: nx and ny must be >= 2'
    stop 1
  end if
  if (nl < 3) then
    write(*, '(A)') 'Error: nl must be >= 3'
    stop 1
  end if

  ! ========================================
  ! Initialize MPI
  ! ========================================
  call MPI_Init(ierr)
  partit%MPI_COMM_FESOM = MPI_COMM_WORLD
  call par_init(partit)

  ! Initialize profiler
  call fesom_profiler_init()

  if (partit%mype == 0) then
    write(*, '(A)') ''
    write(*, '(A)') '================================================'
    write(*, '(A)') 'FESOM2 Tracer Advection Dwarf'
    write(*, '(A)') 'Analytic Mesh Version (no file I/O)'
    write(*, '(A)') '================================================'
    write(*, '(A,I6)') 'Running on ', partit%npes, ' processes'
    write(*, '(A)') '================================================'
    write(*, '(A)') ''
    write(*, '(A,I2,A)') 'Working Precision (WP) = ', WP, ' bytes'
    if (WP == 4) then
      write(*, '(A)') '*** SINGLE PRECISION MODE ***'
    else if (WP == 8) then
      write(*, '(A)') '*** DOUBLE PRECISION MODE ***'
    end if
    write(*, '(A)') ''
  end if

  ! ========================================
  ! Generate analytic mesh
  ! ========================================
  call fesom_profiler_start("mesh_init")
  call generate_analytic_mesh(nx, ny, nl, Lx, Ly, max_depth, partit, mesh, periodic)

  ! Build node-to-edge map for vertex-parallel OpenMP
  call build_node_edge_map(partit, mesh)
  call build_node_edge_packed(partit, mesh)

  ! ========================================
  ! Optional: save mesh for Python visualization
  ! ========================================
  if (save_mesh .and. partit%mype == 0) then
    write(*, '(A,A)') '  --> Writing mesh to ', trim(output_dir)
    call write_mesh_for_python(mesh, trim(output_dir))
  end if

  call fesom_profiler_end("mesh_init")

  ! ========================================
  ! Initialize ALE arrays
  ! ========================================
  call fesom_profiler_start("setup")
  if (partit%mype == 0) write(*, '(A)') '  --> Initializing ALE arrays'

  nz = mesh%nl - 1

  allocate(mesh%hnode(nz, partit%myDim_nod2D+partit%eDim_nod2D))
  allocate(mesh%hnode_new(nz, partit%myDim_nod2D+partit%eDim_nod2D))
  allocate(mesh%helem(nz, partit%myDim_elem2D+partit%eDim_elem2D))
  allocate(mesh%hbar(partit%myDim_nod2D+partit%eDim_nod2D))
  allocate(mesh%hbar_old(partit%myDim_nod2D+partit%eDim_nod2D))
  allocate(mesh%dhe(partit%myDim_elem2D))
  allocate(mesh%zbar_3d_n(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D))
  allocate(mesh%Z_3d_n(nz, partit%myDim_nod2D+partit%eDim_nod2D))

  mesh%hbar = 0.0_MP
  mesh%hbar_old = 0.0_MP
  mesh%dhe = 0.0_MP
  mesh%hnode = 0.0_MP
  mesh%hnode_new = 0.0_MP
  mesh%helem = 0.0_MP
  mesh%zbar_3d_n = 0.0_MP
  mesh%Z_3d_n = 0.0_MP

  ! Initialize layer thicknesses from zbar
  do n = 1, partit%myDim_nod2D+partit%eDim_nod2D
    do istep = 1, nz
      mesh%hnode(istep, n) = mesh%zbar(istep) - mesh%zbar(istep+1)
      mesh%hnode_new(istep, n) = mesh%hnode(istep, n)
    end do
  end do

  do n = 1, partit%myDim_elem2D+partit%eDim_elem2D
    do istep = 1, nz
      mesh%helem(istep, n) = mesh%zbar(istep) - mesh%zbar(istep+1)
    end do
  end do

  ! Initialize 3D depth arrays
  do n = 1, partit%myDim_nod2D+partit%eDim_nod2D
    mesh%zbar_3d_n(:, n) = mesh%zbar(:)
    mesh%Z_3d_n(:, n) = mesh%Z(:)
  end do

  ! ========================================
  ! Allocate dynamics arrays
  ! ========================================
  if (partit%mype == 0) write(*, '(A)') '  --> Allocating dynamics arrays'

  allocate(dynamics%uv(2, nz, partit%myDim_elem2D+partit%eDim_elem2D))
  allocate(dynamics%w(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D))
  ! Uniform horizontal current: 0.1 m/s in x-direction
  dynamics%uv(1, :, :) = 0.1_WP
  dynamics%uv(2, :, :) = 0.0_WP
  dynamics%w  = 0.0_WP

  ! ========================================
  ! Allocate tracer arrays
  ! ========================================
  if (partit%mype == 0) write(*, '(A)') '  --> Allocating tracer arrays'

  tracers%num_tracers = 2
  allocate(tracers%data(tracers%num_tracers))

  do n = 1, tracers%num_tracers
    allocate(tracers%data(n)%values(nz, partit%myDim_nod2D+partit%eDim_nod2D))
    allocate(tracers%data(n)%valuesAB(nz, partit%myDim_nod2D+partit%eDim_nod2D))
    allocate(tracers%data(n)%valuesold(nz, partit%myDim_nod2D+partit%eDim_nod2D, tracers%data(n)%AB_order))

    tracers%data(n)%values    = 0.0_WP
    tracers%data(n)%valuesAB  = 0.0_WP
    tracers%data(n)%valuesold = 0.0_WP

    tracers%data(n)%tra_adv_hor = 'UPW1'
    tracers%data(n)%tra_adv_ver = 'UPW1'
    tracers%data(n)%tra_adv_lim = 'NONE'
  end do

  ! ========================================
  ! Allocate work arrays and init FCT
  ! ========================================
  if (partit%mype == 0) write(*, '(A)') '  --> Allocating work arrays'

  allocate(tracers%work%del_ttf(nz, partit%myDim_nod2D+partit%eDim_nod2D))
  allocate(tracers%work%del_ttf_advhoriz(nz, partit%myDim_nod2D+partit%eDim_nod2D))
  allocate(tracers%work%del_ttf_advvert(nz, partit%myDim_nod2D+partit%eDim_nod2D))

  tracers%work%del_ttf = 0.0_MP
  tracers%work%del_ttf_advhoriz = 0.0_MP
  tracers%work%del_ttf_advvert = 0.0_MP

  if (partit%mype == 0) write(*, '(A)') '  --> Initializing FCT arrays'
  call oce_adv_tra_fct_init(tracers%work, partit, mesh)

  ! ========================================
  ! Set custom tracer values
  ! ========================================
  if (partit%mype == 0) then
    write(*, '(A)') ''
    write(*, '(A)') '================================================'
    write(*, '(A)') 'Setting custom tracer values'
    write(*, '(A)') '================================================'
  end if

  xc = 0.5_MP * Lx   ! domain center x (meters)
  yc = 0.5_MP * Ly   ! domain center y (meters)

  ! Temperature: sinusoidal pattern — warm ridge + vertical stratification
  do n = 1, partit%myDim_nod2D
    xn = mesh%coord_nod2D(1, n) * r_earth   ! x in meters
    yn = mesh%coord_nod2D(2, n) * r_earth   ! y in meters
    do istep = 1, nz
      zn = real(istep-1, WP) / real(nz-1, WP)   ! 0 at surface, 1 at bottom
      tracers%data(1)%values(istep, n) = 15.0_WP &
        + 5.0_WP * sin(2.0_WP * pi * xn / Lx) * cos(2.0_WP * pi * yn / Ly) &
        - 10.0_WP * zn
    end do
  end do

  ! Salinity: Gaussian blob centered in the domain + background
  do n = 1, partit%myDim_nod2D
    xn = mesh%coord_nod2D(1, n) * r_earth
    yn = mesh%coord_nod2D(2, n) * r_earth
    r = sqrt((xn - xc)**2 + (yn - yc)**2)
    do istep = 1, nz
      zn = real(istep-1, WP) / real(nz-1, WP)
      tracers%data(2)%values(istep, n) = 34.0_WP &
        + 2.0_WP * exp(-r**2 / (0.1_WP * Lx)**2) * (1.0_WP - zn)
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
    write(*, '(A)') '  Temperature: sin(x)*cos(y) warm ridge + vertical stratification'
    write(*, '(A)') '  Salinity: Gaussian blob at domain center + depth decay'
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
  tsum = sum(dble(tracers%data(1)%values(:, 1:partit%myDim_nod2D)))

  smin = minval(tracers%data(2)%values(:, 1:partit%myDim_nod2D))
  smax = maxval(tracers%data(2)%values(:, 1:partit%myDim_nod2D))
  ssum = sum(dble(tracers%data(2)%values(:, 1:partit%myDim_nod2D)))

  nan_count = count(.not.(tracers%data(1)%values(:, 1:partit%myDim_nod2D) &
                   == tracers%data(1)%values(:, 1:partit%myDim_nod2D)))

  if (partit%mype == 0) then
    write(*, '(A,2E14.6,E18.10)') '  Temperature: min, max, sum = ', tmin, tmax, tsum
    write(*, '(A,2E14.6,E18.10)') '  Salinity:    min, max, sum = ', smin, smax, ssum
    write(*, '(A,I8)') '  Temperature NaN count = ', nan_count
    write(*, '(A)') ''
  end if

  call fesom_profiler_end("setup")

  ! ========================================
  ! Run advection test
  ! ========================================
  if (partit%mype == 0) then
    write(*, '(A)') '================================================'
    write(*, '(A)') 'Running advection test'
    write(*, '(A)') '================================================'
  end if

  nsteps = 10
  dt_local = 1000.0_WP  ! CFL ~ dt*u/dx = 1000*1/5263 ~ 0.19

  if (partit%mype == 0) then
    write(*, '(A,I6,A,E12.4)') '  Steps: ', nsteps, ', dt = ', dt_local
    write(*, '(A)') ''
  end if

  ! Open scalar output files and write initial condition (step 0)
  if (save_scalars .and. partit%mype == 0) then
    write(*, '(A,A)') '  --> Writing scalar output to ', trim(output_dir)
    call open_scalar_output(mesh, partit%myDim_nod2D, nz, nsteps, &
                            tracers%num_tracers, trim(output_dir), scalar_units)
    call write_scalar_step(scalar_units, tracers%num_tracers, tracers, nz, &
                           partit%myDim_nod2D)
  end if

  call fesom_profiler_start("advection_loop")
  call fesom_profiler_set_timesteps(nsteps)
  do istep = 1, nsteps
    do n = 1, tracers%num_tracers
      ! Zero advection tendencies before each tracer
!$OMP PARALLEL DO COLLAPSE(2)
      do k = 1, partit%myDim_nod2D + partit%eDim_nod2D
        do nzmin = 1, nz
          tracers%work%del_ttf_advhoriz(nzmin, k) = 0.0_MP
          tracers%work%del_ttf_advvert(nzmin, k)  = 0.0_MP
          tracers%work%del_ttf(nzmin, k)          = 0.0_MP
        end do
      end do
!$OMP END PARALLEL DO

      call do_oce_adv_tra(dt_local, dynamics%uv, dynamics%w, dynamics%w, dynamics%w, n, &
                          dynamics, tracers, partit, mesh)

      ! Accumulate horizontal + vertical tendencies into del_ttf
!$OMP PARALLEL DO PRIVATE(nzmax, nzmin, kz)
      do k = 1, partit%myDim_nod2D + partit%eDim_nod2D
        nzmax = mesh%nlevels_nod2D(k) - 1
        nzmin = mesh%ulevels_nod2D(k)
        do kz = nzmin, nzmax
          tracers%work%del_ttf(kz, k) = tracers%work%del_ttf(kz, k) &
            + tracers%work%del_ttf_advhoriz(kz, k) &
            + tracers%work%del_ttf_advvert(kz, k)
        end do
      end do
!$OMP END PARALLEL DO

      ! Apply tendency to tracer values
!$OMP PARALLEL DO PRIVATE(nzmax, nzmin, kz)
      do k = 1, partit%myDim_nod2D
        nzmax = mesh%nlevels_nod2D(k) - 1
        nzmin = mesh%ulevels_nod2D(k)
        do kz = nzmin, nzmax
          tracers%data(n)%values(kz, k) = tracers%data(n)%values(kz, k) &
            + tracers%work%del_ttf(kz, k) / mesh%hnode_new(kz, k)
        end do
      end do
!$OMP END PARALLEL DO
    end do

    if (partit%mype == 0) then
      tmin = minval(tracers%data(1)%values(:, 1:partit%myDim_nod2D))
      tmax = maxval(tracers%data(1)%values(:, 1:partit%myDim_nod2D))
      tsum = sum(dble(tracers%data(1)%values(:, 1:partit%myDim_nod2D)))
      nan_count = count(.not.(tracers%data(1)%values(:, 1:partit%myDim_nod2D) &
                       == tracers%data(1)%values(:, 1:partit%myDim_nod2D)))

      write(*, '(A,I4,A,2E14.6,E18.10,A,I6)') '  Step ', istep, &
        ': T min, max, sum = ', tmin, tmax, tsum, '  NaN:', nan_count
    end if

    ! Write scalar frame for this step
    if (save_scalars .and. partit%mype == 0) then
      call write_scalar_step(scalar_units, tracers%num_tracers, tracers, nz, &
                             partit%myDim_nod2D)
    end if
  end do

  ! Close scalar output files
  if (save_scalars .and. partit%mype == 0) then
    call close_scalar_output(scalar_units, tracers%num_tracers)
  end if
  call fesom_profiler_end("advection_loop")

  if (partit%mype == 0) then
    write(*, '(A)') ''
    write(*, '(A)') '================================================'
    write(*, '(A)') 'Advection test complete!'
    write(*, '(A)') '================================================'
    write(*, '(A)') ''
  end if

  ! ========================================
  ! Finalize profiler and MPI
  ! ========================================
  call fesom_profiler_finalize(partit%MPI_COMM_FESOM, partit%mype)
  call par_ex(partit%MPI_COMM_FESOM, partit%mype)

  if (partit%mype == 0) then
    write(*, '(A)') 'Program finished successfully'
  end if

end program tracer_dwarf_analytic
