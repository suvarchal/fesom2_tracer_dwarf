!> @file tracer_c_interface.F90
!> @brief C-compatible interface for tracer advection library
!> @details
!> This module provides C-compatible wrapper functions using iso_c_binding
!> to enable Python and other languages to use the FESOM tracer advection library.
!>
!> The module wraps the main workflow from fesom.F90 into callable functions
!> that can be accessed via ctypes or similar foreign function interfaces.

module tracer_c_interface
  use iso_c_binding
  use iso_fortran_env, only: error_unit
  use o_PARAM, only: MAX_PATH, WP
  use MOD_MESH, only: t_mesh
  use MOD_PARTIT, only: t_partit
  use MOD_TRACER, only: t_tracer
  use MOD_DYN, only: t_dyn
  use mod_parsup, only: par_ex
  use par_support_interfaces, only: par_init, init_mpi_types, init_gatherLists
  use restart_derivedtype_module, only: read_all_bin_restarts, write_all_bin_restarts
  use fortran_utils, only: int_to_txt
  use oce_adv_tra_driver_module, only: do_oce_adv_tra
  use g_comm_auto, only: exchange_nod

  implicit none

  ! Module-level storage for derived types
  type(t_mesh),   target, save :: global_mesh
  type(t_partit), target, save :: global_partit
  type(t_tracer), target, save :: global_tracers
  type(t_dyn),    target, save :: global_dyn

  ! Status flags
  logical, save :: tracer_initialized = .false. !.true. ! .false.
  logical, save :: mpi_initialized_by_us = .false.

  ! Paths
  character(len=MAX_PATH), save :: restart_path = ''
  character(len=MAX_PATH), save :: output_path = ''

contains

  !> @brief Initialize MPI (if not already initialized)
  !> @param[in] mpi_comm MPI communicator (optional, use MPI_COMM_WORLD if 0)
  !> @return 0 on success, non-zero on error
  function tracer_init_mpi(mpi_comm) result(status) bind(C, name="tracer_init_mpi")
    use mpi
    integer(c_int), intent(in), value :: mpi_comm
    integer(c_int) :: status
    logical :: mpi_is_initialized
    integer :: ierr

    status = 0

    ! Check if MPI is already initialized
    call MPI_Initialized(mpi_is_initialized, ierr)
    if (ierr /= 0) then
      write(unit=error_unit, fmt='(A)') '### error: MPI_Initialized failed'
      status = ierr
      return
    end if

    if (.not. mpi_is_initialized) then
      ! MPI not initialized, we need to initialize it
      call MPI_Init(ierr)
      if (ierr /= 0) then
        write(unit=error_unit, fmt='(A)') '### error: MPI_Init failed'
        status = ierr
        return
      end if
      mpi_initialized_by_us = .true.
    end if

    ! Set up the partit structure with MPI info
    if (mpi_comm == 0) then
      global_partit%MPI_COMM_FESOM = MPI_COMM_WORLD
    else
      global_partit%MPI_COMM_FESOM = mpi_comm
    end if

    ! Initialize partitioning parameters
    call par_init(global_partit)

  end function tracer_init_mpi


  !> @brief Finalize MPI (if we initialized it)
  !> @return 0 on success, non-zero on error
  function tracer_finalize_mpi() result(status) bind(C, name="tracer_finalize_mpi")
    use mpi
    integer(c_int) :: status
    integer :: ierr

    status = 0

    if (mpi_initialized_by_us) then
      call MPI_Finalize(ierr)
      if (ierr /= 0) then
        write(unit=error_unit, fmt='(A)') '### error: MPI_Finalize failed'
        status = ierr
        return
      end if
      mpi_initialized_by_us = .false.
    end if

  end function tracer_finalize_mpi


  !> @brief Initialize tracer advection from binary restart files
  !> @param[in] restart_dir Path to restart directory (C string)
  !> @param[in] npes Number of MPI processes (must match restart files)
  !> @return 0 on success, non-zero on error
  function tracer_init(restart_dir) result(status) bind(C, name="tracer_init")
    character(kind=c_char), dimension(*), intent(in) :: restart_dir
    integer(c_int) :: status

    character(len=MAX_PATH) :: npepath
    logical :: dir_exist
    integer :: i

    status = 0

    ! Convert C string to Fortran string
    i = 1
    do while (restart_dir(i) /= c_null_char .and. i <= MAX_PATH)
      restart_path(i:i) = restart_dir(i)
      i = i + 1
    end do
    restart_path(i:) = ' '

    ! Construct path to restart files
    npepath = trim(restart_path) // "/np" // int_to_txt(global_partit%npes)

    ! Check if restart directory exists
    inquire(file=trim(npepath), exist=dir_exist)
    if (.not. dir_exist) then
      write(unit=error_unit, fmt='(3A)') &
        '### error: could not find restart directory: ', trim(npepath)
      status = -1
      return
    end if

    ! Read binary restart files
    call read_all_bin_restarts(npepath, partit=global_partit, mesh=global_mesh, &
                               dynamics=global_dyn, tracers=global_tracers)

    ! Initialize MPI types
    call init_mpi_types(global_partit, global_mesh)
    call init_gatherLists(global_partit)

    tracer_initialized = .true.

    if (global_partit%mype == 0) then
      write(*, '(A,A)') 'Tracer initialized from: ', trim(npepath)
    end if

  end function tracer_init


  !> @brief Perform one tracer advection time step
  !> @param[in] dt Time step size
  !> @return 0 on success, non-zero on error
  function tracer_advect_step(dt) result(status) bind(C, name="tracer_advect_step")
    real(c_double), intent(in), value :: dt
    integer(c_int) :: status
    integer :: n, nz, nzmax, nzmin

    if (.not. tracer_initialized) then
      write(unit=error_unit, fmt='(A)') &
        '### error: tracer not initialized. Call tracer_init first.'
      status = -1
      return
    end if

    status = 0

    ! Zero out advection flux arrays
    !$OMP PARALLEL DO PRIVATE(n, nz)
    do n = 1, global_partit%myDim_nod2D + global_partit%eDim_nod2D
      do nz = 1, global_mesh%nl - 1
        global_tracers%work%del_ttf_advhoriz(nz, n) = 0.0_WP
        global_tracers%work%del_ttf_advvert(nz, n)  = 0.0_WP
      end do
    end do
    !$OMP END PARALLEL DO

    ! Perform tracer advection
    call do_oce_adv_tra(dt, global_dyn%uv, global_dyn%w, global_dyn%w_i, &
                        global_dyn%w_e, 1, global_dyn, global_tracers, &
                        global_partit, global_mesh)

    ! Update tracer values with total flux
    !$OMP PARALLEL DO PRIVATE(n, nz, nzmax, nzmin)
    do n = 1, global_partit%myDim_nod2D + global_partit%eDim_nod2D
      nzmax = global_mesh%nlevels_nod2D(n) - 1
      nzmin = global_mesh%ulevels_nod2D(n)
      do nz = nzmin, nzmax
        global_tracers%work%del_ttf(nz, n) = global_tracers%work%del_ttf(nz, n) + &
                                             global_tracers%work%del_ttf_advhoriz(nz, n) + &
                                             global_tracers%work%del_ttf_advvert(nz, n)
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(n, nz, nzmax, nzmin)
    do n = 1, global_partit%myDim_nod2D
      nzmax = global_mesh%nlevels_nod2D(n) - 1
      nzmin = global_mesh%ulevels_nod2D(n)
      do nz = nzmin, nzmax
        global_tracers%data(1)%values(nz, n) = global_tracers%data(1)%values(nz, n) + &
                                                global_tracers%work%del_ttf(nz, n) / &
                                                global_mesh%hnode_new(nz, n)
      end do
    end do
    !$OMP END PARALLEL DO

    ! Exchange halo values
    call exchange_nod(global_tracers%data(1)%values(:,:), global_partit, luse_g2g=.true.)

  end function tracer_advect_step


  !> @brief Get tracer statistics (min, max, sum)
  !> @param[out] tmin Minimum tracer value
  !> @param[out] tmax Maximum tracer value
  !> @param[out] tsum Sum of tracer values
  !> @return 0 on success, non-zero on error
  function tracer_get_stats(tmin, tmax, tsum) result(status) &
       bind(C, name="tracer_get_stats")
    real(c_double), intent(out) :: tmin, tmax, tsum
    integer(c_int) :: status

    if (.not. tracer_initialized) then
      write(unit=error_unit, fmt='(A)') &
        '### error: tracer not initialized.'
      status = -1
      tmin = 0.0_c_double
      tmax = 0.0_c_double
      tsum = 0.0_c_double
      return
    end if

    tmin = real(minval(global_tracers%data(1)%values), c_double)
    tmax = real(maxval(global_tracers%data(1)%values), c_double)
    tsum = real(sum(global_tracers%data(1)%values), c_double)
    status = 0

  end function tracer_get_stats


  !> @brief Get tracer array size
  !> @param[out] nz Number of vertical levels
  !> @param[out] nn Number of nodes
  !> @return 0 on success, non-zero on error
  function tracer_get_size(nz, nn) result(status) bind(C, name="tracer_get_size")
    integer(c_int), intent(out) :: nz, nn
    integer(c_int) :: status

    if (.not. tracer_initialized) then
      write(unit=error_unit, fmt='(A)') &
        '### error: tracer not initialized.'
      status = -1
      nz = 0
      nn = 0
      return
    end if

    nz = int(global_mesh%nl - 1, c_int)
    nn = int(global_partit%myDim_nod2D, c_int)
    status = 0

  end function tracer_get_size


  !> @brief Get tracer values (column-major for Fortran/Python compatibility)
  !> @param[out] values Tracer values array (nz x nn)
  !> @param[in] nz Number of vertical levels
  !> @param[in] nn Number of nodes
  !> @return 0 on success, non-zero on error
  function tracer_get_values(values, nz, nn) result(status) bind(C, name="tracer_get_values")
    integer(c_int), intent(in), value :: nz, nn
    real(c_double), dimension(nz, nn), intent(out) :: values
    integer(c_int) :: status
    integer :: i, j

    if (.not. tracer_initialized) then
      write(unit=error_unit, fmt='(A)') &
        '### error: tracer not initialized.'
      status = -1
      return
    end if

    ! Copy tracer values
    do j = 1, nn
      do i = 1, nz
        values(i, j) = real(global_tracers%data(1)%values(i, j), c_double)
      end do
    end do

    status = 0

  end function tracer_get_values


  !> @brief Write restart files and cleanup
  !> @param[in] output_dir Path to output directory (C string)
  !> @return 0 on success, non-zero on error
  function tracer_finalize(output_dir) result(status) bind(C, name="tracer_finalize")
    character(kind=c_char), dimension(*), intent(in) :: output_dir
    integer(c_int) :: status

    character(len=MAX_PATH) :: npepath, meta
    integer, dimension(3) :: time = (/0, 0, 0/)
    integer :: i

    if (.not. tracer_initialized) then
      write(unit=error_unit, fmt='(A)') &
        '### error: tracer not initialized.'
      status = -1
      return
    end if

    status = 0

    ! Convert C string to Fortran string
    i = 1
    do while (output_dir(i) /= c_null_char .and. i <= MAX_PATH)
      output_path(i:i) = output_dir(i)
      i = i + 1
    end do
    output_path(i:) = ' '

    ! Construct output paths
    npepath = trim(output_path) // "/np" // int_to_txt(global_partit%npes)
    meta = trim(output_path) // "/meta.time"

    ! Write binary restart files
    call write_all_bin_restarts(time, npepath, meta, partit=global_partit, &
                                mesh=global_mesh, dynamics=global_dyn, &
                                tracers=global_tracers)

    if (global_partit%mype == 0) then
      write(*, '(A,A)') 'Restart files written to: ', trim(npepath)
    end if

    tracer_initialized = .false.

  end function tracer_finalize


  !> @brief Complete workflow: init, run N steps, finalize
  !> @param[in] restart_dir Path to restart directory (C string)
  !> @param[in] nsteps Number of time steps to run
  !> @param[in] dt Time step size
  !> @return 0 on success, non-zero on error
  function tracer_run_workflow(restart_dir, nsteps, dt) result(status) &
       bind(C, name="tracer_run_workflow")
    character(kind=c_char), dimension(*), intent(in) :: restart_dir
    integer(c_int), intent(in), value :: nsteps
    real(c_double), intent(in), value :: dt
    integer(c_int) :: status
    integer :: step
    real(c_double) :: tmin, tmax, tsum

    ! Initialize
    status = tracer_init(restart_dir)
    if (status /= 0) return

    ! Run time steps
    do step = 1, nsteps
      status = tracer_advect_step(dt)
      if (status /= 0) return

      ! Get statistics
      status = tracer_get_stats(tmin, tmax, tsum)
      if (status /= 0) return

      if (global_partit%mype == 0) then
        write(*, '(A,I3,A,3ES15.6)') 'Step ', step, ': min/max/sum = ', tmin, tmax, tsum
      end if
    end do

    ! Finalize (write to same directory)
    status = tracer_finalize(restart_dir)

  end function tracer_run_workflow


  !> @brief Reset/cleanup module state
  !> @return 0 on success, non-zero on error
  function tracer_cleanup() result(status) bind(C, name="tracer_cleanup")
    integer(c_int) :: status

    tracer_initialized = .false.
    restart_path = ''
    output_path = ''
    status = 0

    ! Note: Actual deallocation of derived type structures should be done carefully
    ! as they may contain many allocated arrays

  end function tracer_cleanup


  !============================================================================
  ! GENERIC INITIALIZATION FUNCTIONS (for Python array initialization)
  !============================================================================

  !> @brief Set mesh dimensions
  !> @param[in] nl Number of vertical levels
  !> @param[in] nod2D Number of 2D nodes
  !> @param[in] edge2D Number of edges
  !> @param[in] elem2D Number of elements
  !> @return 0 on success, non-zero on error
  function tracer_set_mesh_dims(nl, nod2D, edge2D, elem2D) result(status) &
       bind(C, name="tracer_set_mesh_dims")
    integer(c_int), intent(in), value :: nl, nod2D, edge2D, elem2D
    integer(c_int) :: status

    status = 0

    ! Set mesh dimensions
    global_mesh%nl = nl
    global_mesh%nod2D = nod2D
    global_mesh%edge2D = edge2D
    global_mesh%elem2D = elem2D

    if (global_partit%mype == 0) then
      write(*, '(A,4I8)') 'Set mesh dimensions: nl, nod2D, edge2D, elem2D = ', &
                          nl, nod2D, edge2D, elem2D
    end if

  end function tracer_set_mesh_dims


  !> @brief Set vertical level arrays
  !> @param[in] nlevels Bottom level index for each node (nod2D)
  !> @param[in] ulevels Top level index for each node (nod2D)
  !> @param[in] nod2D Number of 2D nodes
  !> @return 0 on success, non-zero on error
  function tracer_set_levels(nlevels, ulevels, nod2D) result(status) &
       bind(C, name="tracer_set_levels")
    integer(c_int), intent(in), value :: nod2D
    integer(c_int), dimension(nod2D), intent(in) :: nlevels, ulevels
    integer(c_int) :: status
    integer :: i

    status = 0

    ! Allocate level arrays if not already allocated
    if (.not. allocated(global_mesh%nlevels_nod2D)) then
      allocate(global_mesh%nlevels_nod2D(nod2D))
    end if
    if (.not. allocated(global_mesh%ulevels_nod2D)) then
      allocate(global_mesh%ulevels_nod2D(nod2D))
    end if

    ! Copy level arrays
    do i = 1, nod2D
      global_mesh%nlevels_nod2D(i) = nlevels(i)
      global_mesh%ulevels_nod2D(i) = ulevels(i)
    end do

    if (global_partit%mype == 0) then
      write(*, '(A,I8,A,2I6)') 'Set level arrays for ', nod2D, ' nodes, range: ', &
                               minval(nlevels), maxval(nlevels)
    end if

  end function tracer_set_levels


  !> @brief Set layer thickness array
  !> @param[in] hnode Layer thickness at nodes (nl-1, nod2D)
  !> @param[in] nl Number of vertical levels
  !> @param[in] nod2D Number of 2D nodes
  !> @return 0 on success, non-zero on error
  function tracer_set_thickness(hnode, nl, nod2D) result(status) &
       bind(C, name="tracer_set_thickness")
    integer(c_int), intent(in), value :: nl, nod2D
    real(c_double), dimension(nl-1, nod2D), intent(in) :: hnode
    integer(c_int) :: status
    integer :: i, j

    status = 0

    ! Allocate thickness array if not already allocated
    if (.not. allocated(global_mesh%hnode_new)) then
      allocate(global_mesh%hnode_new(nl-1, nod2D))
    end if

    ! Copy thickness array (convert from c_double to WP)
    do j = 1, nod2D
      do i = 1, nl-1
        global_mesh%hnode_new(i, j) = real(hnode(i, j), WP)
      end do
    end do

    if (global_partit%mype == 0) then
      write(*, '(A,2I8,A,2ES12.4)') 'Set thickness array: ', nl-1, nod2D, &
                                    ', range: ', minval(hnode), maxval(hnode)
    end if

  end function tracer_set_thickness


  !> @brief Set partition information
  !> @param[in] myDim_nod2D Number of owned nodes
  !> @param[in] eDim_nod2D Number of halo nodes
  !> @return 0 on success, non-zero on error
  function tracer_set_partition(myDim_nod2D, eDim_nod2D) result(status) &
       bind(C, name="tracer_set_partition")
    integer(c_int), intent(in), value :: myDim_nod2D, eDim_nod2D
    integer(c_int) :: status

    status = 0

    ! Set partition dimensions
    global_partit%myDim_nod2D = myDim_nod2D
    global_partit%eDim_nod2D = eDim_nod2D

    ! Initialize MPI types (needed for halo exchange)
    call init_mpi_types(global_partit, global_mesh)
    call init_gatherLists(global_partit)

    if (global_partit%mype == 0) then
      write(*, '(A,2I8)') 'Set partition: myDim_nod2D, eDim_nod2D = ', &
                          myDim_nod2D, eDim_nod2D
    end if

  end function tracer_set_partition


  !> @brief Allocate tracer arrays
  !> @param[in] num_tracers Number of tracers (e.g., 1 for T only, 2 for T+S)
  !> @param[in] nl Number of vertical levels
  !> @param[in] nod2D Number of 2D nodes
  !> @param[in] AB_order Adams-Bashforth order (default 2)
  !> @return 0 on success, non-zero on error
  function tracer_allocate_tracers(num_tracers, nl, nod2D, AB_order) result(status) &
       bind(C, name="tracer_allocate_tracers")
    integer(c_int), intent(in), value :: num_tracers, nl, nod2D, AB_order
    integer(c_int) :: status
    integer :: n, node_size

    status = 0
    node_size = global_partit%myDim_nod2D + global_partit%eDim_nod2D

    ! Set number of tracers
    global_tracers%num_tracers = num_tracers

    ! Allocate tracer data array
    if (allocated(global_tracers%data)) then
      deallocate(global_tracers%data)
    end if
    allocate(global_tracers%data(num_tracers))

    ! Allocate arrays for each tracer
    do n = 1, num_tracers
      allocate(global_tracers%data(n)%values(nl-1, node_size))
      allocate(global_tracers%data(n)%valuesAB(nl-1, node_size))
      allocate(global_tracers%data(n)%valuesold(AB_order-1, nl-1, node_size))
      
      ! Initialize to zero
      global_tracers%data(n)%values = 0.0_WP
      global_tracers%data(n)%valuesAB = 0.0_WP
      global_tracers%data(n)%valuesold = 0.0_WP
      
      ! Set default parameters (simple upwind scheme for testing)
      global_tracers%data(n)%AB_order = AB_order
      global_tracers%data(n)%ID = n
      global_tracers%data(n)%tra_adv_hor = 'UPW1'
      global_tracers%data(n)%tra_adv_ver = 'UPW1'
      global_tracers%data(n)%tra_adv_lim = 'NONE'
      global_tracers%data(n)%tra_adv_ph = 1.0_WP
      global_tracers%data(n)%tra_adv_pv = 1.0_WP
      global_tracers%data(n)%smooth_bh_tra = .false.
      global_tracers%data(n)%i_vert_diff = .false.
    end do

    ! Allocate work arrays
    if (allocated(global_tracers%work%del_ttf)) then
      deallocate(global_tracers%work%del_ttf)
    end if
    if (allocated(global_tracers%work%del_ttf_advhoriz)) then
      deallocate(global_tracers%work%del_ttf_advhoriz)
    end if
    if (allocated(global_tracers%work%del_ttf_advvert)) then
      deallocate(global_tracers%work%del_ttf_advvert)
    end if
    
    allocate(global_tracers%work%del_ttf(nl-1, node_size))
    allocate(global_tracers%work%del_ttf_advhoriz(nl-1, node_size))
    allocate(global_tracers%work%del_ttf_advvert(nl-1, node_size))
    
    global_tracers%work%del_ttf = 0.0_WP
    global_tracers%work%del_ttf_advhoriz = 0.0_WP
    global_tracers%work%del_ttf_advvert = 0.0_WP

    if (global_partit%mype == 0) then
      write(*, '(A,I3,A,I3,A,2I8)') 'Allocated ', num_tracers, ' tracers with AB_order=', &
                                    AB_order, ', dims: ', nl-1, node_size
    end if

  end function tracer_allocate_tracers


  !> @brief Set tracer values from array
  !> @param[in] tracer_id Tracer ID (1-based, e.g., 1=temperature, 2=salinity)
  !> @param[in] values Tracer values (nl-1, nod2D)
  !> @param[in] nl Number of vertical levels
  !> @param[in] nod2D Number of 2D nodes
  !> @return 0 on success, non-zero on error
  function tracer_set_values(tracer_id, values, nl, nod2D) result(status) &
       bind(C, name="tracer_set_values")
    integer(c_int), intent(in), value :: tracer_id, nl, nod2D
    real(c_double), dimension(nl-1, nod2D), intent(in) :: values
    integer(c_int) :: status
    integer :: i, j

    status = 0

    ! Validate tracer_id
    if (tracer_id < 1 .or. tracer_id > global_tracers%num_tracers) then
      write(unit=error_unit, fmt='(A,I3,A,I3)') &
        '### error: invalid tracer_id ', tracer_id, ', num_tracers=', global_tracers%num_tracers
      status = -1
      return
    end if

    ! Copy values (convert from c_double to WP)
    do j = 1, nod2D
      do i = 1, nl-1
        global_tracers%data(tracer_id)%values(i, j) = real(values(i, j), WP)
        global_tracers%data(tracer_id)%valuesAB(i, j) = real(values(i, j), WP)
        global_tracers%data(tracer_id)%valuesold(1, i, j) = real(values(i, j), WP)
      end do
    end do

    if (global_partit%mype == 0) then
      write(*, '(A,I3,A,2ES12.4)') 'Set tracer ', tracer_id, ' values, range: ', &
                                   minval(values), maxval(values)
    end if

  end function tracer_set_values


  !> @brief Set velocity fields to zero (for testing)
  !> @param[in] nl Number of vertical levels
  !> @param[in] nod2D Number of 2D nodes
  !> @param[in] edge2D Number of edges
  !> @return 0 on success, non-zero on error
  function tracer_set_velocity_zero(nl, nod2D, edge2D) result(status) &
       bind(C, name="tracer_set_velocity_zero")
    integer(c_int), intent(in), value :: nl, nod2D, edge2D
    integer(c_int) :: status

    status = 0

    ! Allocate velocity arrays if not already allocated
    if (.not. allocated(global_dyn%uv)) then
      allocate(global_dyn%uv(nl-1, edge2D, 2))
    end if
    if (.not. allocated(global_dyn%w)) then
      allocate(global_dyn%w(nl, nod2D))
    end if
    if (.not. allocated(global_dyn%w_e)) then
      allocate(global_dyn%w_e(nl, edge2D))
    end if
    if (.not. allocated(global_dyn%w_i)) then
      allocate(global_dyn%w_i(nl, nod2D))
    end if

    ! Set all velocities to zero
    global_dyn%uv = 0.0_WP
    global_dyn%w = 0.0_WP
    global_dyn%w_e = 0.0_WP
    global_dyn%w_i = 0.0_WP

    if (global_partit%mype == 0) then
      write(*, '(A,3I8)') 'Set velocity fields to zero: nl, nod2D, edge2D = ', &
                          nl, nod2D, edge2D
    end if

  end function tracer_set_velocity_zero


  !> @brief Set velocity fields from arrays
  !> @param[in] uv Horizontal velocity (nl-1, edge2D, 2)
  !> @param[in] w Vertical velocity (nl, nod2D)
  !> @param[in] nl Number of vertical levels
  !> @param[in] nod2D Number of 2D nodes
  !> @param[in] edge2D Number of edges
  !> @return 0 on success, non-zero on error
  function tracer_set_velocity(uv, w, nl, nod2D, edge2D) result(status) &
       bind(C, name="tracer_set_velocity")
    integer(c_int), intent(in), value :: nl, nod2D, edge2D
    real(c_double), dimension(nl-1, edge2D, 2), intent(in) :: uv
    real(c_double), dimension(nl, nod2D), intent(in) :: w
    integer(c_int) :: status
    integer :: i, j, k

    status = 0

    ! Allocate velocity arrays if not already allocated
    if (.not. allocated(global_dyn%uv)) then
      allocate(global_dyn%uv(nl-1, edge2D, 2))
    end if
    if (.not. allocated(global_dyn%w)) then
      allocate(global_dyn%w(nl, nod2D))
    end if
    if (.not. allocated(global_dyn%w_e)) then
      allocate(global_dyn%w_e(nl, edge2D))
    end if
    if (.not. allocated(global_dyn%w_i)) then
      allocate(global_dyn%w_i(nl, nod2D))
    end if

    ! Copy horizontal velocity (convert from c_double to WP)
    do k = 1, 2
      do j = 1, edge2D
        do i = 1, nl-1
          global_dyn%uv(i, j, k) = real(uv(i, j, k), WP)
        end do
      end do
    end do

    ! Copy vertical velocity
    do j = 1, nod2D
      do i = 1, nl
        global_dyn%w(i, j) = real(w(i, j), WP)
      end do
    end do

    ! Set w_e and w_i to zero (can be derived later if needed)
    global_dyn%w_e = 0.0_WP
    global_dyn%w_i = 0.0_WP

    if (global_partit%mype == 0) then
      write(*, '(A,2ES12.4,A,2ES12.4)') 'Set velocity: uv range=', &
                                        minval(uv), maxval(uv), ', w range=', &
                                        minval(w), maxval(w)
    end if

  end function tracer_set_velocity


  !> @brief Mark initialization as complete
  !> @return 0 on success, non-zero on error
  function tracer_init_complete() result(status) bind(C, name="tracer_init_complete")
    integer(c_int) :: status

    status = 0
    tracer_initialized = .true.

    if (global_partit%mype == 0) then
      write(*, '(A)') 'Tracer initialization from arrays complete!'
    end if

  end function tracer_init_complete


  !> @brief Load mesh from partition files (alternative to array initialization)
  !> @param[in] restart_dir Path to restart directory containing mesh partition
  !> @return 0 on success, non-zero on error
  !> @details This function loads a complete mesh from existing partition files,
  !>          providing all connectivity and communication arrays needed for advection.
  !>          Use this instead of tracer_init() when you want to set tracer values
  !>          from arrays but need a real mesh structure.
  function tracer_load_mesh_partition(restart_dir) result(status) &
       bind(C, name="tracer_load_mesh_partition")
    character(kind=c_char), dimension(*), intent(in) :: restart_dir
    integer(c_int) :: status

    character(len=MAX_PATH) :: npepath
    logical :: dir_exist
    integer :: i

    status = 0

    ! Convert C string to Fortran string
    i = 1
    do while (restart_dir(i) /= c_null_char .and. i <= MAX_PATH)
      restart_path(i:i) = restart_dir(i)
      i = i + 1
    end do
    restart_path(i:) = ' '

    ! Construct path to restart files
    npepath = trim(restart_path) // "/np" // int_to_txt(global_partit%npes)

    ! Check if restart directory exists
    inquire(file=trim(npepath), exist=dir_exist)
    if (.not. dir_exist) then
      write(unit=error_unit, fmt='(3A)') &
        '### error: could not find mesh partition directory: ', trim(npepath)
      status = -1
      return
    end if

    ! Read mesh and partition info from binary restart files
    ! This loads all connectivity arrays needed for advection
    call read_all_bin_restarts(npepath, partit=global_partit, mesh=global_mesh, &
                               dynamics=global_dyn)

    ! Initialize MPI types for communication
    call init_mpi_types(global_partit, global_mesh)
    call init_gatherLists(global_partit)

    ! Mark as initialized so other functions can work
    tracer_initialized = .true.

    if (global_partit%mype == 0) then
      write(*, '(A,A)') 'Loaded mesh partition from: ', trim(npepath)
      write(*, '(A,I8,A,I8,A,I8)') '  Mesh: nl=', global_mesh%nl, &
                                   ', nod2D=', global_mesh%nod2D, &
                                   ', edge2D=', global_mesh%edge2D
      write(*, '(A,I8,A,I8)') '  Partition: myDim=', global_partit%myDim_nod2D, &
                              ', eDim=', global_partit%eDim_nod2D
    end if

  end function tracer_load_mesh_partition


  !> @brief Load mesh from mesh partition files (not restart files)
  !> @param[in] mesh_path Path to mesh directory
  !> @return 0 on success, non-zero on error
  !> @details This function loads mesh from partition files using mesh_setup(),
  !>          like fesom_mesh_init.F90 does. This reads from mesh files, not restart files.
  !>          The MeshPath from g_config is used.
  function tracer_load_mesh_from_files(mesh_path) result(status) &
       bind(C, name="tracer_load_mesh_from_files")
    use g_config, only: MeshPath
    use atlas_fesom_mesh_module, only: mesh_setup_with_atlas
    character(kind=c_char), dimension(*), intent(in) :: mesh_path
    integer(c_int) :: status
    
    integer :: i
    
    status = 0
    
    ! Convert C string to Fortran string and set MeshPath
    i = 1
    do while (mesh_path(i) /= c_null_char .and. i <= MAX_PATH)
      MeshPath(i:i) = mesh_path(i)
      i = i + 1
    end do
    MeshPath(i:) = ' '
    
    if (global_partit%mype == 0) then
      write(*, '(A,A)') 'Loading mesh from partition files: ', trim(MeshPath)
    end if
    
    ! Load mesh using Atlas-backed setup when enabled, otherwise fallback to mesh_setup
    if (atlas_fesom_enabled()) then
        call mesh_setup_with_atlas(global_partit, global_mesh)
    else
        call mesh_setup(global_partit, global_mesh)
    end if
    
    ! Mark as initialized
    tracer_initialized = .true.
    
    if (global_partit%mype == 0) then
      write(*, '(A)') 'Mesh loaded successfully from partition files'
      write(*, '(A,I8,A,I8,A,I8)') '  Mesh: nl=', global_mesh%nl, &
                                   ', nod2D=', global_mesh%nod2D, &
                                   ', edge2D=', global_mesh%edge2D
      write(*, '(A,I8,A,I8)') '  Partition: myDim=', global_partit%myDim_nod2D, &
                              ', eDim=', global_partit%eDim_nod2D
    end if
    
  end function tracer_load_mesh_from_files


  !> @brief Initialize ALE arrays (layer thickness) needed for advection
  !> @return 0 on success, non-zero on error
  !> @details This initializes layer thickness arrays from mesh%zbar.
  !>          Must be called after mesh is loaded and before advection.
  function tracer_init_ale_arrays() result(status) &
       bind(C, name="tracer_init_ale_arrays")
    integer(c_int) :: status
    integer :: n, k
    
    status = 0
    
    if (.not. tracer_initialized) then
      write(unit=error_unit, fmt='(A)') &
        '### error: tracer_init_ale_arrays called before initialization'
      status = -1
      return
    end if
    
    ! Allocate layer thickness arrays at nodes
    if (.not. allocated(global_mesh%hnode)) then
      allocate(global_mesh%hnode(global_mesh%nl-1, &
               global_partit%myDim_nod2D+global_partit%eDim_nod2D))
      allocate(global_mesh%hnode_new(global_mesh%nl-1, &
               global_partit%myDim_nod2D+global_partit%eDim_nod2D))
    end if
    
    ! Allocate layer thickness at elements
    if (.not. allocated(global_mesh%helem)) then
      allocate(global_mesh%helem(global_mesh%nl-1, &
               global_partit%myDim_elem2D+global_partit%eDim_elem2D))
    end if
    
    ! Initialize layer thicknesses from zbar
    do n = 1, global_partit%myDim_nod2D+global_partit%eDim_nod2D
      do k = 1, global_mesh%nl-1
        global_mesh%hnode(k, n) = global_mesh%zbar(k) - global_mesh%zbar(k+1)
        global_mesh%hnode_new(k, n) = global_mesh%hnode(k, n)
      end do
    end do
    
    do n = 1, global_partit%myDim_elem2D+global_partit%eDim_elem2D
      do k = 1, global_mesh%nl-1
        global_mesh%helem(k, n) = global_mesh%zbar(k) - global_mesh%zbar(k+1)
      end do
    end do
    
    if (global_partit%mype == 0) then
      write(*, '(A)') 'ALE arrays initialized for advection'
    end if
    
  end function tracer_init_ale_arrays


  !> @brief Complete initialization for advection (dynamics + work arrays + FCT)
  !> @return 0 on success, non-zero on error
  !> @details This performs remaining initialization needed for advection:
  !>          - Allocates dynamics arrays (uv, w, w_e, w_i)
  !>          - Allocates tracer work arrays
  !>          - Initializes FCT arrays
  !>          Call this after loading mesh, allocating tracers, and setting velocities.
  function tracer_init_for_advection() result(status) &
       bind(C, name="tracer_init_for_advection")
    use oce_adv_tra_fct_module, only: oce_adv_tra_fct_init
    integer(c_int) :: status
    integer :: n
    
    status = 0
    
    if (.not. tracer_initialized) then
      write(unit=error_unit, fmt='(A)') &
        '### error: tracer_init_for_advection called before mesh loaded'
      status = -1
      return
    end if
    
    if (global_partit%mype == 0) then
      write(*, '(A)') 'Initializing remaining arrays for advection...'
    end if
    
    ! Allocate dynamics arrays if not already done
    if (.not. allocated(global_dyn%uv)) then
      allocate(global_dyn%uv(2, global_mesh%nl-1, &
               global_partit%myDim_elem2D+global_partit%eDim_elem2D))
      global_dyn%uv = 0.0_WP
    end if
    
    if (.not. allocated(global_dyn%w)) then
      allocate(global_dyn%w(global_mesh%nl, &
               global_partit%myDim_nod2D+global_partit%eDim_nod2D))
      global_dyn%w = 0.0_WP
    end if
    
    ! Allocate w_e and w_i (vertical velocities at edges/interfaces)
    if (.not. allocated(global_dyn%w_e)) then
      allocate(global_dyn%w_e(global_mesh%nl, &
               global_partit%myDim_edge2D+global_partit%eDim_edge2D))
      global_dyn%w_e = 0.0_WP
    end if
    
    if (.not. allocated(global_dyn%w_i)) then
      allocate(global_dyn%w_i(global_mesh%nl, &
               global_partit%myDim_nod2D+global_partit%eDim_nod2D))
      global_dyn%w_i = 0.0_WP
    end if
    
    ! Allocate work arrays if not already done
    if (.not. allocated(global_tracers%work%del_ttf)) then
      allocate(global_tracers%work%del_ttf(global_mesh%nl-1, &
               global_partit%myDim_nod2D+global_partit%eDim_nod2D))
      allocate(global_tracers%work%del_ttf_advhoriz(global_mesh%nl-1, &
               global_partit%myDim_nod2D+global_partit%eDim_nod2D))
      allocate(global_tracers%work%del_ttf_advvert(global_mesh%nl-1, &
               global_partit%myDim_nod2D+global_partit%eDim_nod2D))
      
      global_tracers%work%del_ttf = 0.0_WP
      global_tracers%work%del_ttf_advhoriz = 0.0_WP
      global_tracers%work%del_ttf_advvert = 0.0_WP
    end if
    
    ! Initialize FCT arrays
    call oce_adv_tra_fct_init(global_tracers%work, global_partit, global_mesh)
    
    if (global_partit%mype == 0) then
      write(*, '(A)') 'Full initialization complete - ready for advection!'
    end if
    
  end function tracer_init_for_advection

end module tracer_c_interface
