!> @file tracer_init_from_mesh.F90
!! @brief Initialize tracer dwarf from mesh files (not restart files)
!! @details This module provides initialization routines that read mesh from
!!          mesh files (like main FESOM in initial mode) rather than from restart files.
!!          This is cleaner for testing with custom tracer values.

module tracer_init_from_mesh_module
  use MOD_MESH
  use MOD_PARTIT
  use MOD_DYN
  use MOD_TRACER
  use MOD_PARSUP
  use oce_mesh_module, only: mesh_setup
  use atlas_fesom_mesh_module, only: mesh_setup_with_atlas, atlas_fesom_enabled
  use oce_adv_tra_fct_module
  use oce_muscl_adv_module
  use o_PARAM
  
  implicit none
  private
  
  public :: tracer_init_mesh_and_arrays
  
contains

  !> @brief Initialize mesh from mesh files and set up tracer arrays
  !! @details This follows the pattern of fesom_init but simplified for tracer advection only
  !! @param[inout] partit Partition structure
  !! @param[inout] mesh Mesh structure  
  !! @param[inout] dynamics Dynamics structure
  !! @param[inout] tracers Tracer structure
  subroutine tracer_init_mesh_and_arrays(partit, mesh, dynamics, tracers)
    use iso_fortran_env, only: output_unit
    
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh
    type(t_dyn),    intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    
    integer :: n, nz, istep
    
    if (partit%mype == 0) then
      write(output_unit, '(A)') ''
      write(output_unit, '(A)') '=========================================='
      write(output_unit, '(A)') 'Tracer Dwarf: Initializing from mesh files'
      write(output_unit, '(A)') '=========================================='
    end if
    
    ! ========================================
    ! 1. Read configuration from namelists
    ! ========================================
    if (partit%mype == 0) write(output_unit, '(A)') '  --> Reading configuration'
    ! Note: read_config is called by the main program before this routine
    
    ! ========================================
    ! 2. Set up mesh from mesh files
    ! ========================================
    if (partit%mype == 0) write(output_unit, '(A)') '  --> Setting up mesh from files'
    if (atlas_fesom_enabled()) then
      call mesh_setup_with_atlas(partit, mesh)
    else
      call mesh_setup(partit, mesh)
    end if
    
    if (partit%mype == 0) then
      write(output_unit, '(A,I8)') '      Mesh nl     = ', mesh%nl
      write(output_unit, '(A,I8)') '      Mesh nod2D  = ', mesh%nod2D
      write(output_unit, '(A,I8)') '      Mesh edge2D = ', mesh%edge2D
      write(output_unit, '(A,I8)') '      Partit myDim_nod2D = ', partit%myDim_nod2D
      write(output_unit, '(A,I8)') '      Partit eDim_nod2D  = ', partit%eDim_nod2D
    end if
    
    ! ========================================
    ! 3. Initialize ALE arrays (simplified version of init_ale)
    ! ========================================
    if (partit%mype == 0) write(output_unit, '(A)') '  --> Initializing ALE arrays'
    
    nz = mesh%nl - 1
    
    ! Allocate layer thickness arrays at nodes
    allocate(mesh%hnode(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D))
    allocate(mesh%hnode_new(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D))
    
    ! Allocate layer thickness at elements
    allocate(mesh%helem(mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D))
    
    ! Allocate elevation arrays
    allocate(mesh%hbar(partit%myDim_nod2D+partit%eDim_nod2D))
    allocate(mesh%hbar_old(partit%myDim_nod2D+partit%eDim_nod2D))
    
    ! Allocate dhe (depth increment on elements)
    allocate(mesh%dhe(partit%myDim_elem2D))
    
    ! Allocate 3D depth arrays
    allocate(mesh%zbar_3d_n(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D))
    allocate(mesh%Z_3d_n(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D))
    
    ! Initialize to zero
    mesh%hbar = 0.0_WP
    mesh%hbar_old = 0.0_WP
    mesh%dhe = 0.0_WP
    mesh%hnode = 0.0_WP
    mesh%hnode_new = 0.0_WP
    mesh%helem = 0.0_WP
    mesh%zbar_3d_n = 0.0_WP
    mesh%Z_3d_n = 0.0_WP
    
    ! Initialize layer thicknesses from zbar (constant in time for this test)
    ! For nodes
    do n = 1, partit%myDim_nod2D+partit%eDim_nod2D
      do istep = 1, mesh%nl-1
        mesh%hnode(istep, n) = mesh%zbar(istep) - mesh%zbar(istep+1)
        mesh%hnode_new(istep, n) = mesh%hnode(istep, n)
      end do
    end do
    
    ! For elements
    do n = 1, partit%myDim_elem2D+partit%eDim_elem2D
      do istep = 1, mesh%nl-1
        mesh%helem(istep, n) = mesh%zbar(istep) - mesh%zbar(istep+1)
      end do
    end do
    
    ! Initialize 3D depth arrays (simplified - no partial cells)
    do n = 1, partit%myDim_nod2D+partit%eDim_nod2D
      mesh%zbar_3d_n(:, n) = mesh%zbar(:)
      mesh%Z_3d_n(:, n) = mesh%Z(:)
    end do
    
    if (partit%mype == 0) write(output_unit, '(A)') '      ALE arrays initialized'
    
    ! ========================================
    ! 4. Allocate dynamics arrays
    ! ========================================
    if (partit%mype == 0) write(output_unit, '(A)') '  --> Allocating dynamics arrays'
    
    ! Allocate velocity arrays
    ! uv is on elements, w is on nodes
    allocate(dynamics%uv(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D))
    allocate(dynamics%w(mesh%nl, partit%myDim_nod2D+partit%eDim_nod2D))
    
    ! Initialize to zero
    dynamics%uv = 0.0_WP
    dynamics%w  = 0.0_WP
    
    if (partit%mype == 0) then
      write(output_unit, '(A)') '      Velocity arrays allocated and zeroed'
      write(output_unit, '(A,I8)') '        uv: myDim_elem2D+eDim_elem2D = ', &
                                    partit%myDim_elem2D+partit%eDim_elem2D
      write(output_unit, '(A,I8)') '        w:  myDim_nod2D+eDim_nod2D   = ', &
                                    partit%myDim_nod2D+partit%eDim_nod2D
    end if
    
    ! ========================================
    ! 5. Allocate tracer arrays
    ! ========================================
    if (partit%mype == 0) write(output_unit, '(A)') '  --> Allocating tracer arrays'
    
    ! Set number of tracers (temperature and salinity)
    tracers%num_tracers = 2
    
    ! Allocate tracer data
    allocate(tracers%data(tracers%num_tracers))
    
    ! Allocate arrays for each tracer
    do n = 1, tracers%num_tracers
      allocate(tracers%data(n)%values(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D))
      allocate(tracers%data(n)%valuesAB(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D))
      allocate(tracers%data(n)%valuesold(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, tracers%data(n)%AB_order))
      
      ! Initialize to zero
      tracers%data(n)%values    = 0.0_WP
      tracers%data(n)%valuesAB  = 0.0_WP
      tracers%data(n)%valuesold = 0.0_WP
      
      ! Set default advection schemes
      tracers%data(n)%tra_adv_hor = 'UPW1'
      tracers%data(n)%tra_adv_ver = 'UPW1'  ! First-order upwind for vertical
      tracers%data(n)%tra_adv_lim = 'NONE'  ! No FCT limiter
    end do
    
    if (partit%mype == 0) then
      write(output_unit, '(A,I2,A)') '      Allocated ', tracers%num_tracers, ' tracers'
    end if
    
    ! ========================================
    ! 6. Allocate work arrays for advection
    ! ========================================
    if (partit%mype == 0) write(output_unit, '(A)') '  --> Allocating work arrays'
    
    ! Work arrays are in tracers%work, not in individual tracer data
    allocate(tracers%work%del_ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D))
    allocate(tracers%work%del_ttf_advhoriz(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D))
    allocate(tracers%work%del_ttf_advvert(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D))
    
    tracers%work%del_ttf = 0.0_WP
    tracers%work%del_ttf_advhoriz = 0.0_WP
    tracers%work%del_ttf_advvert = 0.0_WP
    
    if (partit%mype == 0) write(output_unit, '(A)') '      Work arrays allocated'
    
    ! ========================================
    ! 7. Initialize FCT and MUSCL arrays
    ! ========================================
    if (partit%mype == 0) write(output_unit, '(A)') '  --> Initializing FCT arrays'
    call oce_adv_tra_fct_init(tracers%work, partit, mesh)
    
    ! NOTE: Skipping muscl_adv_init for now as it requires SSH stiffness matrix
    ! which needs init_stiff_mat_ale. UPW1 scheme doesn't need MUSCL arrays.
    ! if (partit%mype == 0) write(output_unit, '(A)') '  --> Initializing MUSCL arrays'
    ! call muscl_adv_init(tracers%work, partit, mesh)
    
    ! ========================================
    ! Done
    ! ========================================
    if (partit%mype == 0) then
      write(output_unit, '(A)') ''
      write(output_unit, '(A)') '=========================================='
      write(output_unit, '(A)') 'Initialization complete!'
      write(output_unit, '(A)') '=========================================='
      write(output_unit, '(A)') ''
    end if
    
  end subroutine tracer_init_mesh_and_arrays

end module tracer_init_from_mesh_module
