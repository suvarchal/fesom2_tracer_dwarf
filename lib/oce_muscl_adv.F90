module oce_muscl_adv_module
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_TRACER
    use o_ARRAYS
    use o_PARAM
    use g_comm_auto
    use g_config
    use oce_node_edge_map_module
#ifdef USE_HALF_PRECISION
    use hp_math_intrinsics
#endif

    implicit none
    
    private
    public :: muscl_adv_init, find_up_downwind_triangles, fill_up_dn_grad

contains

! A set of routines to implement MUSCL-type of advection
! For description, see Abalakin, I., Dervieux, A., Kozubskaya, T., 2002. A
! vertex-centered high-order MUSCL scheme applying to linearized Euler acoustics.
! INRIA, Rapport de recherche 4459.
!
!  The advection routine is solve_tracer_muscl
!  muscl_adv_init initializes several arrays needed for the algorithm
!  
!  The algorithm works with the concept of upwind and downwind triangles to a given
!  edge. This introduces additional halo communications.
!  sergey.danilov@awi.de 2012
!
!Contains:
!	muscl_adv_init
!	find_up_downwind_triangles
!	fill_up_dn_grad
!	adv_tracer_muscl
subroutine muscl_adv_init(twork, partit, mesh)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_TRACER
    use o_ARRAYS
    use o_PARAM
    use g_comm_auto
    use g_config
    ! find_up_downwind_triangles is now in the same module
    IMPLICIT NONE
    integer     :: n, k, n1, n2

    type(t_mesh),        intent(inout), target :: mesh
    type(t_partit),      intent(inout), target :: partit
    type(t_tracer_work), intent(inout), target :: twork

!#include "associate_part_def.h"
!#include "associate_mesh_def.h"
!#include "associate_part_ass.h"
!#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    ! find upwind and downwind triangle for each local edge 
    call find_up_downwind_triangles(twork, partit, mesh)
    
    !___________________________________________________________________________
    mesh%nn_size=0
    k=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n)
!$OMP DO REDUCTION(max: k)
    do n=1, partit%myDim_nod2D
        ! get number of  neighbouring nodes from sparse stiffness matrix
        ! stiffnes matrix filled up in subroutine init_stiff_mat_ale
        ! --> mesh%SSH_stiff%rowptr... compressed row index of sparse matrix
        ! --> mesh%SSH_stiff%values... array with values at row column location, has legth nod2d+1
        ! --> mesh%SSH_stiff%rowptr(n) ... gives index location in mesh%SSH_stiff%values where the 
        !                             next value switches to a new row
        ! --> mesh%SSH_stiff%rowptr(n+1)-mesh%SSH_stiff%rowptr(n) gives maximum number of 
        !     neighbouring nodes within a single row of the sparse matrix
        k=max(k, mesh%SSH_stiff%rowptr(n+1)-mesh%SSH_stiff%rowptr(n))
    end do
!$OMP END DO    
!$OMP END PARALLEL
    mesh%nn_size=k
    !___________________________________________________________________________
    allocate(mesh%nn_num(partit%myDim_nod2D), mesh%nn_pos(mesh%nn_size,partit%myDim_nod2D))
    ! These are the same arrays that we also use in quadratic reconstruction
    !MOVE IT TO SOMEWHERE ELSE
!$OMP PARALLEL DO
    do n=1, partit%myDim_nod2D
        ! number of neigbouring nodes to node n
        mesh%nn_num(n)=1
        ! local position of neigbouring nodes
        mesh%nn_pos(1,n)=n
    end do   
!$OMP END PARALLEL DO
    !___________________________________________________________________________
    allocate(twork%nboundary_lay(partit%myDim_nod2D+partit%eDim_nod2D)) !node n becomes a boundary node after layer twork%nboundary_lay(n)
    twork%nboundary_lay=mesh%nl-1
    ! Vertex-gather: race-free OpenMP parallelization over nodes
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, k, n1, n2)
!$OMP DO
    do n=1, partit%myDim_nod2D
        ! Build nn_pos by gathering neighbor nodes from incident edges
        do k=1, node_edge_num(n)
            ! Determine the other endpoint of this edge
            if (node_edge_sign(k, n) == 1) then
                n2 = mesh%edges(2, node_edge_idx(k, n))
            else
                n2 = mesh%edges(1, node_edge_idx(k, n))
            end if
            mesh%nn_pos(mesh%nn_num(n)+1, n) = n2
            mesh%nn_num(n) = mesh%nn_num(n) + 1
        end do

        ! Update nboundary_lay from incident edges
        do k=1, node_edge_num(n)
            n1 = node_edge_idx(k, n)  ! edge index
            if (any(mesh%edge_tri(:, n1) <= 0)) then
                twork%nboundary_lay(n) = 0
            else
                twork%nboundary_lay(n) = min(twork%nboundary_lay(n), &
                    minval(mesh%nlevels(mesh%edge_tri(:, n1))) - 1)
            end if
        end do
    end do
!$OMP END DO
    ! Also update nboundary_lay for extended nodes
!$OMP DO
    do n=partit%myDim_nod2D+1, partit%myDim_nod2D+partit%eDim_nod2D
        do k=1, node_edge_num(n)
            n1 = node_edge_idx(k, n)  ! edge index
            if (any(mesh%edge_tri(:, n1) <= 0)) then
                twork%nboundary_lay(n) = 0
            else
                twork%nboundary_lay(n) = min(twork%nboundary_lay(n), &
                    minval(mesh%nlevels(mesh%edge_tri(:, n1))) - 1)
            end if
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
end SUBROUTINE muscl_adv_init
!
!
!_______________________________________________________________________________
SUBROUTINE find_up_downwind_triangles(twork, partit, mesh)
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
USE MOD_TRACER
USE o_ARRAYS
USE o_PARAM
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
integer                    :: n, k, ednodes(2), elem, el
real(kind=MP)              :: x(2),b(2), c(2), cr, bx, by, xx, xy, ab, ax
real(kind=MP), allocatable :: coord_elem(:, :,:)
real(kind=WP), allocatable :: temp(:)
integer, allocatable       :: temp_i(:), e_nodes(:,:)

type(t_mesh),        intent(in)   , target :: mesh
type(t_partit),      intent(inout), target :: partit
type(t_tracer_work), intent(inout), target :: twork
!#include "associate_part_def.h"
!#include "associate_mesh_def.h"
!#include "associate_part_ass.h"
!#include "associate_mesh_ass.h"

allocate(twork%edge_up_dn_tri(2,partit%myDim_edge2D))
allocate(twork%edge_up_dn_grad(4,mesh%nl-1,partit%myDim_edge2D))
twork%edge_up_dn_tri=0
! =====
! In order that this procedure works, we need to know nodes and their coordinates 
! on the extended set of elements (not only my, but myDim+eDim+eXDim) 
! =====
allocate(coord_elem(2, 3, partit%myDim_elem2D+partit%eDim_elem2D+partit%eXDim_elem2D))
allocate(temp(partit%myDim_elem2D+partit%eDim_elem2D+partit%eXDim_elem2D))
   DO n=1,3
        DO k=1,2
!$OMP PARALLEL
!$OMP DO
           DO el=1,partit%myDim_elem2D
              temp(el)=mesh%coord_nod2D(k,mesh%elem2D_nodes(n,el))
           END DO
!$OMP END DO
!$OMP MASTER
	   call exchange_elem(temp, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
           DO el=1, partit%myDim_elem2D+partit%eDim_elem2D+partit%eXDim_elem2D
              coord_elem(k,n,el)=temp(el)
           END DO
!$OMP END DO
!$OMP END PARALLEL
	END DO
   END DO
deallocate(temp)

allocate(e_nodes(3,partit%myDim_elem2D+partit%eDim_elem2D+partit%eXDim_elem2D))
allocate(temp_i(partit%myDim_elem2D+partit%eDim_elem2D+partit%eXDim_elem2D))
    DO n=1,3
!$OMP PARALLEL
!$OMP DO
       do el=1,partit%myDim_elem2D
          temp_i(el)=partit%myList_nod2D(mesh%elem2D_nodes(n,el))
       end do
!$OMP END DO
!$OMP MASTER
       call exchange_elem(temp_i, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
       DO el=1, partit%myDim_elem2D+partit%eDim_elem2D+partit%eXDim_elem2D
          e_nodes(n, el)=temp_i(el)
       END DO
!$OMP END DO
!$OMP END PARALLEL
    END DO   
deallocate(temp_i)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, k, ednodes, elem, el, x,b, c, cr, bx, by, xx, xy, ab, ax)
!$OMP DO
DO n=1, partit%myDim_edge2D
   ednodes=mesh%edges(:,n) 
   x=mesh%coord_nod2D(:,ednodes(2))-mesh%coord_nod2D(:,ednodes(1))
      	 if(x(1)>cyclic_length/2._WP) x(1)=x(1)-cyclic_length
         if(x(1)<-cyclic_length/2._WP) x(1)=x(1)+cyclic_length
	
   ! Find upwind (in the sense of x) triangle, i. e. 
   ! find which triangle contains -x:
   x=-x
   DO k=1,mesh%nod_in_elem2D_num(ednodes(1))
      elem=mesh%nod_in_elem2D(k,ednodes(1))
      if(e_nodes(1,elem)==partit%myList_nod2D(ednodes(1))) then
	 b=coord_elem(:,2,elem)-coord_elem(:,1,elem)
	 c=coord_elem(:,3,elem)-coord_elem(:,1,elem)
      elseif(e_nodes(2,elem)==partit%myList_nod2D(ednodes(1))) then
	 b=coord_elem(:,1,elem)-coord_elem(:,2,elem)
	 c=coord_elem(:,3,elem)-coord_elem(:,2,elem)
      else	 
	 b=coord_elem(:,1,elem)-coord_elem(:,3,elem)
	 c=coord_elem(:,2,elem)-coord_elem(:,3,elem)
      end if
      	 if(b(1)>cyclic_length/2._WP) b(1)=b(1)-cyclic_length
         if(b(1)<-cyclic_length/2._WP) b(1)=b(1)+cyclic_length
	 if(c(1)>cyclic_length/2._WP) c(1)=c(1)-cyclic_length
         if(c(1)<-cyclic_length/2._WP) c(1)=c(1)+cyclic_length
      ! the vector x has to be between b and c
      ! Decompose b and x into parts along c and along (-cy,cx), i.e.
      ! 90 degree counterclockwise
      cr=sum(c*c)
      bx=sum(b*c)/cr
      by=(-b(1)*c(2)+b(2)*c(1))/cr
      xx=sum(x*c)/cr
      xy=(-x(1)*c(2)+x(2)*c(1))/cr
      ab=atan2(by,bx)
      ax=atan2(xy,xx)
      ! Since b and c are the sides of triangle, |ab|<pi, and atan2 should 
      ! be what is needed
      if((ab>0.0_WP).and.(ax>0.0_WP).and.(ax<ab)) then
      twork%edge_up_dn_tri(1,n)=elem
      cycle
      endif
      if((ab<0.0_WP).and.(ax<0.0_WP).and.(ax>ab)) then
      twork%edge_up_dn_tri(1,n)=elem
      cycle
      endif
      if((ab==ax).or.(ax==0.0_WP)) then
      twork%edge_up_dn_tri(1,n)=elem
      cycle
      endif
   END DO
   ! Find downwind element
   x=-x
   DO k=1,mesh%nod_in_elem2D_num(ednodes(2))
      elem=mesh%nod_in_elem2D(k,ednodes(2))
      if(e_nodes(1,elem)==partit%myList_nod2D(ednodes(2))) then
      	 b=coord_elem(:,2,elem)-coord_elem(:,1,elem)
	 c=coord_elem(:,3,elem)-coord_elem(:,1,elem)
      elseif(e_nodes(2, elem)==partit%myList_nod2D(ednodes(2))) then
	 b=coord_elem(:,1,elem)-coord_elem(:,2,elem)
	 c=coord_elem(:,3,elem)-coord_elem(:,2,elem)
      else	 
	 b=coord_elem(:,1,elem)-coord_elem(:,3,elem)
	 c=coord_elem(:,2,elem)-coord_elem(:,3,elem)
      end if
      	 if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
         if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
	 if(c(1)>cyclic_length/2.) c(1)=c(1)-cyclic_length
         if(c(1)<-cyclic_length/2.) c(1)=c(1)+cyclic_length
      ! the vector x has to be between b and c
      ! Decompose b and x into parts along c and along (-cy,cx), i.e.
      ! 90 degree counterclockwise
      cr=sum(c*c)
      bx=sum(b*c)/cr
      by=(-b(1)*c(2)+b(2)*c(1))/cr
      xx=sum(x*c)/cr
      xy=(-x(1)*c(2)+x(2)*c(1))/cr
      ab=atan2(by,bx)
      ax=atan2(xy,xx)
      ! Since b and c are the sides of triangle, |ab|<pi, and atan2 should 
      ! be what is needed
      if((ab>0.0_WP).and.(ax>0.0_WP).and.(ax<ab)) then
      twork%edge_up_dn_tri(2,n)=elem
      cycle
      endif
      if((ab<0.0_WP).and.(ax<0.0_WP).and.(ax>ab)) then
      twork%edge_up_dn_tri(2,n)=elem
      cycle
      endif
      if((ab==ax).or.(ax==0.0)) then
      twork%edge_up_dn_tri(2,n)=elem
      cycle
      endif
   END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! For edges touching the boundary --- up or downwind elements may be absent.  
! We return to the standard Miura at nodes that
! belong to such edges. Same at the depth.
! Count the number of 'good' edges:
!k=0
!DO n=1, partit%myDim_edge2D
!   if((twork%edge_up_dn_tri(1,n).ne.0).and.(twork%edge_up_dn_tri(2,n).ne.0)) k=k+1
!END DO

!$OMP PARALLEL DO
DO n=1, partit%myDim_edge2D
   twork%edge_up_dn_grad(:, :, n)=0.0_MP
END DO
!$OMP END PARALLEL DO
deallocate(e_nodes, coord_elem)
end SUBROUTINE find_up_downwind_triangles
!
!
!_______________________________________________________________________________
SUBROUTINE fill_up_dn_grad(twork, partit, mesh)
! ttx, tty  elemental gradient of tracer 
USE o_PARAM
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
USE MOD_TRACER
USE o_ARRAYS
IMPLICIT NONE
integer                  :: edge, n, nz, elem, k, ednodes(2), nzmin, nzmax
real(kind=MP)            :: tvol, tx, ty
type(t_mesh),        intent(in),    target :: mesh
type(t_partit),      intent(inout), target :: partit
type(t_tracer_work), intent(inout), target :: twork
!#include "associate_part_def.h"
!#include "associate_mesh_def.h"
!#include "associate_part_ass.h"
!#include "associate_mesh_ass.h"
	!___________________________________________________________________________
	! loop over edge segments
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, n, nz, elem, k, ednodes, nzmin, nzmax, tvol, tx, ty)
!$OMP DO
	DO edge=1, partit%myDim_edge2D
		ednodes=mesh%edges(:,edge)
		!_______________________________________________________________________
		! case when edge has upwind and downwind triangle on the surface
		if((twork%edge_up_dn_tri(1,edge).ne.0.0_WP).and.(twork%edge_up_dn_tri(2,edge).ne.0.0_WP)) then
			nzmin = maxval(mesh%ulevels_nod2D_max(ednodes))
			nzmax = minval(mesh%nlevels_nod2D_min(ednodes))
			
			!___________________________________________________________________
			! loop over not shared depth levels of edge node 1 (ednodes(1))
			DO nz=mesh%ulevels_nod2D(ednodes(1)), nzmin-1
				tvol=0.0_MP
				tx=0.0_MP
				ty=0.0_MP
				! loop over number triangles that share the nodeedge points ednodes(1)
				! --> calculate mean gradient at ednodes(1) over the sorounding 
				!     triangle gradients
				DO k=1, mesh%nod_in_elem2D_num(ednodes(1))
					elem=mesh%nod_in_elem2D(k,ednodes(1))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(mesh%nlevels(elem)-1<nz .or. nz<mesh%ulevels(elem)) cycle
					tvol=tvol+mesh%elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*mesh%elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*mesh%elem_area(elem)
				END DO
				twork%edge_up_dn_grad(1,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(3,nz,edge)=ty/tvol
			END DO
			
			!___________________________________________________________________
			! loop over not shared depth levels of edge node 2 (ednodes(2))
			DO nz=mesh%ulevels_nod2D(ednodes(2)),nzmin-1
				tvol=0.0_MP
				tx=0.0_MP
				ty=0.0_MP
				! loop over number triangles that share the nodeedge points ednodes(2)
				! --> calculate mean gradient at ednodes(2) over the sorounding 
				!     triangle gradients
				DO k=1, mesh%nod_in_elem2D_num(ednodes(2))
					elem=mesh%nod_in_elem2D(k,ednodes(2))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(mesh%nlevels(elem)-1<nz .or. nz<mesh%ulevels(elem)) cycle
					tvol=tvol+mesh%elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*mesh%elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*mesh%elem_area(elem)
				END DO
				twork%edge_up_dn_grad(2,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(4,nz,edge)=ty/tvol
			END DO
			
			!___________________________________________________________________
			! loop over shared depth levels
			!!PS DO nz=1, minval(nlevels_nod2D_min(ednodes))-1
			DO nz=nzmin, nzmax-1
				! tracer gradx for upwind and downwind tri
				twork%edge_up_dn_grad(1:2,nz,edge)=tr_xy(1,nz,twork%edge_up_dn_tri(:,edge))
				! tracer grady for upwind and downwind tri
				twork%edge_up_dn_grad(3:4,nz,edge)=tr_xy(2,nz,twork%edge_up_dn_tri(:,edge))
			END DO
			
			!___________________________________________________________________
			! loop over not shared depth levels of edge node 1 (ednodes(1))
			!!PS DO nz=minval(nlevels_nod2D_min(ednodes)),nlevels_nod2D(ednodes(1))-1
			DO nz=nzmax, mesh%nlevels_nod2D(ednodes(1))-1
				tvol=0.0_MP
				tx=0.0_MP
				ty=0.0_MP
				! loop over number triangles that share the nodeedge points ednodes(1)
				! --> calculate mean gradient at ednodes(1) over the sorounding 
				!     triangle gradients
				DO k=1, mesh%nod_in_elem2D_num(ednodes(1))
					elem=mesh%nod_in_elem2D(k,ednodes(1))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(mesh%nlevels(elem)-1<nz .or. nz<mesh%ulevels(elem)) cycle
					tvol=tvol+mesh%elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*mesh%elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*mesh%elem_area(elem)
				END DO
				twork%edge_up_dn_grad(1,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(3,nz,edge)=ty/tvol
			END DO
			!___________________________________________________________________
			! loop over not shared depth levels of edge node 2 (ednodes(2))
			!!PS DO nz=minval(nlevels_nod2D_min(ednodes)),nlevels_nod2D(ednodes(2))-1
			DO nz=nzmax, mesh%nlevels_nod2D(ednodes(2))-1
				tvol=0.0_MP
				tx=0.0_MP
				ty=0.0_MP
				! loop over number triangles that share the nodeedge points ednodes(2)
				! --> calculate mean gradient at ednodes(2) over the sorounding 
				!     triangle gradients
				DO k=1, mesh%nod_in_elem2D_num(ednodes(2))
					elem=mesh%nod_in_elem2D(k,ednodes(2))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(mesh%nlevels(elem)-1<nz .or. nz<mesh%ulevels(elem)) cycle
					tvol=tvol+mesh%elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*mesh%elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*mesh%elem_area(elem)
				END DO
				twork%edge_up_dn_grad(2,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(4,nz,edge)=ty/tvol
			END DO
		!_______________________________________________________________________
		! case when edge either upwind or downwind triangle on the surface
		! --> surface boundary edge
		else
			! Only linear reconstruction part
			nzmin = mesh%ulevels_nod2D(ednodes(1))
			nzmax = mesh%nlevels_nod2D(ednodes(1))
			!!PS DO nz=1,nlevels_nod2D(ednodes(1))-1
			DO nz=nzmin,nzmax-1
				tvol=0.0_MP
				tx=0.0_MP
				ty=0.0_MP
				DO k=1, mesh%nod_in_elem2D_num(ednodes(1))
					elem=mesh%nod_in_elem2D(k,ednodes(1))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(mesh%nlevels(elem)-1 < nz .or. nz<mesh%ulevels(elem) ) cycle
					tvol=tvol+mesh%elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*mesh%elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*mesh%elem_area(elem)
				END DO
				twork%edge_up_dn_grad(1,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(3,nz,edge)=ty/tvol
			END DO
			nzmin = mesh%ulevels_nod2D(ednodes(2))
			nzmax = mesh%nlevels_nod2D(ednodes(2))
			!!PS DO nz=1,nlevels_nod2D(ednodes(2))-1
			DO nz=nzmin,nzmax-1
				tvol=0.0_MP
				tx=0.0_MP
				ty=0.0_MP
				DO k=1, mesh%nod_in_elem2D_num(ednodes(2))
					elem=mesh%nod_in_elem2D(k,ednodes(2))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(mesh%nlevels(elem)-1 < nz .or. nz<mesh%ulevels(elem) ) cycle
					tvol=tvol+mesh%elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*mesh%elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*mesh%elem_area(elem)
				END DO
				twork%edge_up_dn_grad(2,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(4,nz,edge)=ty/tvol
			END DO
		end if  
	END DO
!$OMP END DO
!$OMP END PARALLEL
END SUBROUTINE fill_up_dn_grad

end module oce_muscl_adv_module
