module oce_adv_tra_driver_module
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use g_comm_auto
    use oce_adv_tra_hor_interfaces
    use oce_adv_tra_ver_interfaces
    use oce_adv_tra_fct_module, only: oce_tra_adv_fct
    use oce_node_edge_map_module
    use fesom_profiler

    implicit none

    private
    public :: do_oce_adv_tra, oce_tra_adv_flux2dtracer

contains

!
!
!===============================================================================
subroutine do_oce_adv_tra(dt, vel, w, wi, we, tr_num, dynamics, tracers, partit, mesh)
    use MOD_MESH
    use MOD_TRACER
    use O_PARAM, only: WP, MP
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    use g_comm_auto
    use oce_adv_tra_hor_interfaces
    use oce_adv_tra_ver_interfaces
    use oce_adv_tra_fct_module, only: oce_tra_adv_fct
    ! oce_tra_adv_flux2dtracer is now in the same module
    implicit none
    real(kind=WP),  intent(in),    target :: dt
    integer,        intent(in)            :: tr_num
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(inout), target :: tracers
    type(t_dyn)   , intent(inout), target :: dynamics
    real(kind=WP),  intent(in)            :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
    real(kind=WP),  intent(in), target    :: W(mesh%nl,    partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP),  intent(in), target    :: WI(mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP),  intent(in), target    :: WE(mesh%nl,   partit%myDim_nod2D+partit%eDim_nod2D)

    real(kind=WP),  pointer, dimension (:,:)   :: pwvel
    real(kind=WP),  pointer, dimension (:,:)   :: ttf, ttfAB
    real(kind=MP),  pointer, dimension (:,:)   :: fct_LO
    real(kind=MP),  pointer, dimension (:,:)   :: adv_flux_hor, adv_flux_ver, dttf_h, dttf_v
    real(kind=MP),  pointer, dimension (:,:)   :: fct_ttf_min, fct_ttf_max
    real(kind=MP),  pointer, dimension (:,:)   :: fct_plus, fct_minus

    integer,        pointer, dimension (:)     :: nboundary_lay
    real(kind=MP),  pointer, dimension (:,:,:) :: edge_up_dn_grad

    integer       :: el(2), enodes(2), nz, n, e
    integer       :: nl12, nu12, nl1, nl2, nu1, nu2
    real(kind=WP) :: cLO, cHO
    real(kind=MP) :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP) :: qc, qu, qd
    real(kind=WP) :: tvert(mesh%nl), tvert_e(mesh%nl), b, c, d, da, db, dg, Tupw1
    real(kind=MP) :: a, vflux
    real(kind=WP) :: Tmean, Tmean1, Tmean2, num_ord
    real(kind=WP) :: opth, optv
    logical       :: do_zero_flux

!#include "associate_part_def.h"
!#include "associate_mesh_def.h"
!#include "associate_part_ass.h"
!#include "associate_mesh_ass.h"
    ttf             => tracers%data(tr_num)%values
    ttfAB           => tracers%data(tr_num)%valuesAB
    opth            =  tracers%data(tr_num)%tra_adv_ph
    optv            =  tracers%data(tr_num)%tra_adv_pv
    fct_LO          => tracers%work%fct_LO
    adv_flux_ver    => tracers%work%adv_flux_ver
    adv_flux_hor    => tracers%work%adv_flux_hor
    edge_up_dn_grad => tracers%work%edge_up_dn_grad
    nboundary_lay   => tracers%work%nboundary_lay
    fct_ttf_min     => tracers%work%fct_ttf_min
    fct_ttf_max     => tracers%work%fct_ttf_max
    fct_plus        => tracers%work%fct_plus
    fct_minus       => tracers%work%fct_minus
    dttf_h          => tracers%work%del_ttf_advhoriz
    dttf_v          => tracers%work%del_ttf_advvert
    
    !___________________________________________________________________________
    ! compute FCT horzontal and vertical low order solution as well as lw order
    ! part of antidiffusive flux
    if (trim(tracers%data(tr_num)%tra_adv_lim)=='FCT') then
        ! compute the low order upwind horizontal flux
        ! o_init_zero=.true.  : zero the horizontal flux before computation
        ! o_init_zero=.false. : input flux will be substracted
        call adv_tra_hor_upw1(vel, ttf, partit, mesh, adv_flux_hor, o_init_zero=.true.)
        ! update the LO solution for horizontal contribution
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
        do n=1, partit%myDim_nod2D+partit%eDim_nod2D
           do nz=1, mesh%nl - 1
              fct_LO(nz,n) = 0.0_MP
           end do
        end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
        !$ACC END PARALLEL LOOP
#endif

#ifndef ENABLE_OPENACC
        ! Vertex-gather: race-free OpenMP parallelization over nodes
        ! Uses pre-packed level bounds; node_edge_idx still needed for dynamic flux array
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, e, nl1, nu1, nl2, nu2, nu12, nl12, nz)
        do n=1, partit%myDim_nod2D
            do e=1, node_edge_num(n)
                nl1=node_edge_nl1(e, n)
                nu1=node_edge_nu1(e, n)
                nl2=node_edge_nl2(e, n)
                nu2=node_edge_nu2(e, n)
                nl12 = max(nl1,nl2)
                nu12 = nu1
                if (nu2>0) nu12 = min(nu1,nu2)
                do nz=nu12, nl12
                    fct_LO(nz, n)=fct_LO(nz, n) + node_edge_sign(e, n) * adv_flux_hor(nz, node_edge_idx(e, n))
                end do
            end do
        end do
!$OMP END PARALLEL DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC PARALLEL LOOP GANG PRIVATE(enodes, el) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#else
        !$ACC UPDATE SELF(fct_lo, adv_flux_hor)
#endif
        do e=1, partit%myDim_edge2D
            enodes=mesh%edges(:,e)
            el=mesh%edge_tri(:,e)
            nl1=mesh%nlevels(el(1))-1
            nu1=mesh%ulevels(el(1))
            nl2=0
            nu2=0
            if(el(2)>0) then
                nl2=mesh%nlevels(el(2))-1
                nu2=mesh%ulevels(el(2))
            end if
            nl12 = max(nl1,nl2)
            nu12 = nu1
            if (nu2>0) nu12 = min(nu1,nu2)
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC LOOP VECTOR
#endif
            do nz=nu12, nl12
#if !defined(DISABLE_OPENACC_ATOMICS)
               !$ACC ATOMIC UPDATE
#endif
               fct_LO(nz, enodes(1))=fct_LO(nz, enodes(1))+adv_flux_hor(nz, e)
#if !defined(DISABLE_OPENACC_ATOMICS)
               !$ACC ATOMIC UPDATE
#endif
               fct_LO(nz, enodes(2))=fct_LO(nz, enodes(2))-adv_flux_hor(nz, e)
            end do
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC END LOOP
#endif
        end do
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC END PARALLEL LOOP
#else
        !$ACC UPDATE DEVICE(fct_lo)
#endif
#endif

        ! compute the low order upwind vertical flux (explicit part only)
        ! zero the input/output flux before computation
        call adv_tra_ver_upw1(we, ttf, partit, mesh, adv_flux_ver, o_init_zero=.true.)
        ! update the LO solution for vertical contribution

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nu1, nl1, nz)
#else
        !$ACC PARALLEL LOOP GANG PRESENT(fct_LO) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
        do n=1, partit%myDim_nod2D
            nu1 = mesh%ulevels_nod2D(n)
            nl1 = mesh%nlevels_nod2D(n)
            !!PS do  nz=1, nlevels_nod2D(n)-1
            !$ACC LOOP VECTOR
            do  nz= nu1, nl1-1
                fct_LO(nz,n)=(ttf(nz,n)*mesh%hnode(nz,n)+(fct_LO(nz,n)+(adv_flux_ver(nz, n)-adv_flux_ver(nz+1, n)))*dt/mesh%areasvol(nz,n))/mesh%hnode_new(nz,n)
            end do
            !$ACC END LOOP
        end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
        !$ACC END PARALLEL LOOP
#endif


        !_______________________________________________________________________

        if (dynamics%use_wsplit) then !wvel/=wvel_e
            ! update for implicit contribution (w_split option)
!when adv_tra_vert_impl is ported to ACC the UPDATEs below wont be needed!
!$ACC UPDATE HOST(fct_LO)
            call adv_tra_vert_impl(dt, wi, fct_LO, partit, mesh)
!$ACC UPDATE DEVICE(fct_LO)
            ! compute the low order upwind vertical flux (full vertical velocity)
            ! zero the input/output flux before computation
            ! --> compute here low order part of vertical anti diffusive fluxes,
            !     has to be done on the full vertical velocity w
            call adv_tra_ver_upw1(w, ttf, partit, mesh, adv_flux_ver, o_init_zero=.true.)
        end if
#if !defined(USE_HALF_PRECISION)
        call exchange_nod(fct_LO, partit, luse_g2g = .true.)
#endif
!$OMP BARRIER
    end if !--> if (trim(tracers%data(tr_num)%tra_adv_lim)=='FCT') then
    
    do_zero_flux=.true.
    if (trim(tracers%data(tr_num)%tra_adv_lim)=='FCT') do_zero_flux=.false.

    !___________________________________________________________________________
    ! Fused path: for non-FCT UPW1+UPW1, compute flux and scatter in one vertex loop
    if (trim(tracers%data(tr_num)%tra_adv_lim)/='FCT' .and. &
        trim(tracers%data(tr_num)%tra_adv_hor)=='UPW1' .and. &
        trim(tracers%data(tr_num)%tra_adv_ver)=='UPW1') then

        pwvel=>we
        ! Compute vertical fluxes (node-indexed, no scatter needed)
        call adv_tra_ver_upw1(pwvel, ttfAB, partit, mesh, adv_flux_ver, o_init_zero=.true.)

        ! Fused: horizontal flux computation + horizontal scatter + vertical scatter
        ! in one vertex loop — no intermediate flux array needed
        call fesom_profiler_start("fused_vertex")
        call adv_tra_upw1_vertex_fused(vel, ttfAB, dt, dttf_h, dttf_v, adv_flux_ver, partit, mesh)
        call fesom_profiler_end("fused_vertex")

    else
    !___________________________________________________________________________
    ! Original path: separate flux computation + scatter (for FCT, MUSCL, etc.)
    call fesom_profiler_start("adv_tra_hor")
    SELECT CASE(trim(tracers%data(tr_num)%tra_adv_hor))
        CASE('MUSCL')
            call adv_tra_hor_muscl(vel, ttfAB, partit, mesh, opth,  adv_flux_hor, edge_up_dn_grad, nboundary_lay, o_init_zero=do_zero_flux)
        CASE('MFCT')
             call adv_tra_hor_mfct(vel, ttfAB, partit, mesh, opth,  adv_flux_hor, edge_up_dn_grad,                o_init_zero=do_zero_flux)
        CASE('UPW1')
             call adv_tra_hor_upw1(vel, ttfAB, partit, mesh,        adv_flux_hor,                                 o_init_zero=do_zero_flux)
        CASE DEFAULT !unknown
            if (partit%mype==0) write(*,*) 'Unknown horizontal advection type ',  trim(tracers%data(tr_num)%tra_adv_hor), '! Check your namelists!'
            call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
    END SELECT
    call fesom_profiler_end("adv_tra_hor")
    if (trim(tracers%data(tr_num)%tra_adv_lim)=='FCT') then
       pwvel=>w
    else
       pwvel=>we
    end if

    call fesom_profiler_start("adv_tra_ver")
    SELECT CASE(trim(tracers%data(tr_num)%tra_adv_ver))
        CASE('QR4C')
            call adv_tra_ver_qr4c (   pwvel, ttfAB, partit, mesh, optv, adv_flux_ver, o_init_zero=do_zero_flux)
        CASE('CDIFF')
            call adv_tra_ver_cdiff(   pwvel, ttfAB, partit, mesh,       adv_flux_ver, o_init_zero=do_zero_flux)
        CASE('PPM')
            call adv_tra_vert_ppm(dt, pwvel, ttfAB, partit, mesh,       adv_flux_ver, o_init_zero=do_zero_flux)
        CASE('UPW1')
            call adv_tra_ver_upw1 (   pwvel, ttfAB, partit, mesh,       adv_flux_ver, o_init_zero=do_zero_flux)
        CASE DEFAULT !unknown
            if (partit%mype==0) write(*,*) 'Unknown vertical advection type ',  trim(tracers%data(tr_num)%tra_adv_ver), '! Check your namelists!'
            call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
    END SELECT
    call fesom_profiler_end("adv_tra_ver")

    call fesom_profiler_start("flux2dtracer")
    if (trim(tracers%data(tr_num)%tra_adv_lim)=='FCT') then
       call oce_tra_adv_fct(dt, ttf, fct_LO, adv_flux_hor, adv_flux_ver, fct_ttf_min, fct_ttf_max, fct_plus, fct_minus, edge_up_dn_grad, partit, mesh)
       call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, adv_flux_hor, adv_flux_ver, partit, mesh, use_lo=.TRUE., ttf=ttf, lo=fct_LO)
    else
       call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, adv_flux_hor, adv_flux_ver, partit, mesh)
    end if
    call fesom_profiler_end("flux2dtracer")
    end if ! fused vs original path

end subroutine do_oce_adv_tra
!
!
!===============================================================================
subroutine oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, flux_h, flux_v, partit, mesh, use_lo, ttf, lo)
    use MOD_MESH
    use o_ARRAYS
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto
    implicit none
    real(kind=WP), intent(in),    target :: dt
    type(t_partit),intent(inout), target :: partit
    type(t_mesh),  intent(in),    target :: mesh
    real(kind=MP), intent(inout)      :: dttf_h(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=MP), intent(inout)      :: dttf_v(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=MP), intent(inout)      :: flux_h(mesh%nl-1, partit%myDim_edge2D)
    real(kind=MP), intent(inout)      :: flux_v(mesh%nl,   partit%myDim_nod2D)
    logical,       optional           :: use_lo
    real(kind=MP), optional           :: lo (mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), optional           :: ttf(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    integer                           :: n, nz, k, elem, enodes(3), num, el(2), nu12, nl12, nu1, nu2, nl1, nl2, edge
!#include "associate_part_def.h"
!#include "associate_mesh_def.h"
!#include "associate_part_ass.h"
!#include "associate_mesh_ass.h"
    !___________________________________________________________________________
    ! c. Update the solution
    ! Vertical
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, nz, el, nu12, nl12, nu1, nu2, nl1, nl2, edge)
#endif
    if (present(use_lo)) then
       if (use_lo) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
          do n=1, partit%myDim_nod2D
             nu1 = mesh%ulevels_nod2D(n)
             nl1 = mesh%nlevels_nod2D(n)
             !!PS do nz=1,nlevels_nod2D(n)-1
             !$ACC LOOP VECTOR
             do nz=nu1, nl1-1
                dttf_v(nz,n)=dttf_v(nz,n)-ttf(nz,n)*mesh%hnode(nz,n)+LO(nz,n)*mesh%hnode_new(nz,n)
             end do
             !$ACC END LOOP
          end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
         !$ACC END PARALLEL LOOP
#endif
       end if
    end if
#ifndef ENABLE_OPENACC
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
    do n=1, partit%myDim_nod2D
        nu1 = mesh%ulevels_nod2D(n)
        nl1 = mesh%nlevels_nod2D(n)
        !$ACC LOOP VECTOR
        do nz=nu1,nl1-1
            dttf_v(nz,n)=dttf_v(nz,n) + (flux_v(nz,n)-flux_v(nz+1,n))*dt/mesh%areasvol(nz,n)
        end do
        !$ACC END LOOP
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
    !$ACC END PARALLEL LOOP
#endif
    ! Horizontal
#ifndef ENABLE_OPENACC
    ! Vertex-gather: race-free OpenMP parallelization over nodes
    ! Uses pre-packed level bounds; node_edge_idx still needed for dynamic flux array
!$OMP DO
    do n=1, partit%myDim_nod2D
        do edge=1, node_edge_num(n)
            nl1=node_edge_nl1(edge, n)
            nu1=node_edge_nu1(edge, n)
            nl2=node_edge_nl2(edge, n)
            nu2=node_edge_nu2(edge, n)
            nl12 = max(nl1,nl2)
            nu12 = nu1
            if (nu2>0) nu12 = min(nu1,nu2)
            do nz=nu12, nl12
                dttf_h(nz,n)=dttf_h(nz,n) + node_edge_sign(edge, n) * flux_h(nz, node_edge_idx(edge, n))*dt/mesh%areasvol(nz,n)
            end do
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC PARALLEL LOOP GANG PRIVATE(enodes, el) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#else
    !$ACC UPDATE SELF(dttf_h, flux_h)
#endif
    do edge=1, partit%myDim_edge2D
        enodes(1:2)=mesh%edges(:,edge)
        el=mesh%edge_tri(:,edge)
        nl1=mesh%nlevels(el(1))-1
        nu1=mesh%ulevels(el(1))
        nl2=0
        nu2=0
        if(el(2)>0) then
            nl2=mesh%nlevels(el(2))-1
            nu2=mesh%ulevels(el(2))
        end if
        nl12 = max(nl1,nl2)
        nu12 = nu1
        if (nu2>0) nu12 = min(nu1,nu2)
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC LOOP VECTOR
#endif
        do nz=nu12, nl12
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC UPDATE
#endif
            dttf_h(nz,enodes(1))=dttf_h(nz,enodes(1))+flux_h(nz,edge)*dt/mesh%areasvol(nz,enodes(1))
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC UPDATE
#endif
            dttf_h(nz,enodes(2))=dttf_h(nz,enodes(2))-flux_h(nz,edge)*dt/mesh%areasvol(nz,enodes(2))
        end do
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC END LOOP
#endif
    end do
#if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC END PARALLEL LOOP
#else
    !$ACC UPDATE DEVICE(dttf_h)
#endif
#endif

end subroutine oce_tra_adv_flux2dtracer

!
!
!===============================================================================
! Fused UPW1: compute horizontal flux + scatter horizontal + scatter vertical
! in one vertex loop. Eliminates the intermediate flux(nz,edge) array.
! Each edge flux is computed twice (once per endpoint) but this is offset by
! eliminating the flux array write/read and having one parallel region.
!===============================================================================
subroutine adv_tra_upw1_vertex_fused(vel, ttf, dt, dttf_h, dttf_v, flux_v, partit, mesh)
    use MOD_MESH
    use O_PARAM, only: WP, MP, r_earth
#ifdef USE_HALF_PRECISION
    use hp_math_intrinsics
#endif
    USE MOD_PARTIT
    USE MOD_PARSUP
    implicit none
    type(t_partit),intent(in),    target :: partit
    type(t_mesh),  intent(in),    target :: mesh
    real(kind=WP), intent(in)            :: dt
    real(kind=WP), intent(in)            :: ttf(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)            :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
    real(kind=MP), intent(inout)         :: dttf_h(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=MP), intent(inout)         :: dttf_v(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=MP), intent(in)            :: flux_v(mesh%nl,   partit%myDim_nod2D)

    integer       :: n, nz, ie, enodes(2), el(2), sgn
    integer       :: nl1, nl2, nu1, nu2, nl12, nu12
    real(kind=MP) :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=MP) :: vflux, flux_val

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, ie, enodes, el, sgn, &
!$OMP              nl1, nl2, nu1, nu2, nl12, nu12, deltaX1, deltaY1, deltaX2, deltaY2, vflux, flux_val)
    do n = 1, partit%myDim_nod2D
        ! --- Vertical scatter (from pre-computed vertical flux) ---
        nu1 = mesh%ulevels_nod2D(n)
        nl1 = mesh%nlevels_nod2D(n)
        do nz = nu1, nl1-1
            dttf_v(nz,n) = dttf_v(nz,n) + (flux_v(nz,n) - flux_v(nz+1,n)) * dt / mesh%areasvol(nz,n)
        end do

        ! --- Fused horizontal: compute edge flux on-the-fly + accumulate ---
        ! Uses pre-packed node-contiguous arrays (no indirect indexing)
        do ie = 1, node_edge_num(n)
            sgn    = node_edge_sign(ie, n)
            enodes = node_edge_enodes(:, ie, n)
            el     = node_edge_elems(:, ie, n)

            ! Geometry from pre-packed arrays
            deltaX1 = node_edge_dxdy(1, ie, n)
            deltaY1 = node_edge_dxdy(2, ie, n)

            nl1 = node_edge_nl1(ie, n)
            nu1 = node_edge_nu1(ie, n)

            ! Element 2 (boundary edges have el(2)==0, nl2/nu2==0)
            nl2 = node_edge_nl2(ie, n)
            nu2 = node_edge_nu2(ie, n)
            deltaX2 = node_edge_dxdy(3, ie, n)
            deltaY2 = node_edge_dxdy(4, ie, n)

            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)

            ! (A) Only el(1) contributes (cavity surface levels)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                flux_val = -0.5_WP * (ttf(nz,enodes(1))*(vflux+abs(vflux)) + ttf(nz,enodes(2))*(vflux-abs(vflux)))
                dttf_h(nz, n) = dttf_h(nz, n) + sgn * flux_val * dt / mesh%areasvol(nz, n)
            end do

            ! (B) Only el(2) contributes (cavity surface levels)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    flux_val = -0.5_WP * (ttf(nz,enodes(1))*(vflux+abs(vflux)) + ttf(nz,enodes(2))*(vflux-abs(vflux)))
                    dttf_h(nz, n) = dttf_h(nz, n) + sgn * flux_val * dt / mesh%areasvol(nz, n)
                end do
            end if

            ! (C) Both elements contribute (main case)
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                flux_val = -0.5_WP * (ttf(nz,enodes(1))*(vflux+abs(vflux)) + ttf(nz,enodes(2))*(vflux-abs(vflux)))
                dttf_h(nz, n) = dttf_h(nz, n) + sgn * flux_val * dt / mesh%areasvol(nz, n)
            end do

            ! (D) Remaining levels, only el(1)
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                flux_val = -0.5_WP * (ttf(nz,enodes(1))*(vflux+abs(vflux)) + ttf(nz,enodes(2))*(vflux-abs(vflux)))
                dttf_h(nz, n) = dttf_h(nz, n) + sgn * flux_val * dt / mesh%areasvol(nz, n)
            end do

            ! (E) Remaining levels, only el(2)
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                flux_val = -0.5_WP * (ttf(nz,enodes(1))*(vflux+abs(vflux)) + ttf(nz,enodes(2))*(vflux-abs(vflux)))
                dttf_h(nz, n) = dttf_h(nz, n) + sgn * flux_val * dt / mesh%areasvol(nz, n)
            end do
        end do
    end do
!$OMP END PARALLEL DO
end subroutine adv_tra_upw1_vertex_fused

end module oce_adv_tra_driver_module
