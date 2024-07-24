! Copyright (c) 2023-2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements an inverse distance weighting based source term
module idw_source_term
  use num_types, only : rp, dp
  use field_list, only : field_list_t
  use json_module, only : json_file, json_value, json_core
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use field, only : field_t
  use utils, only : neko_error
  use tri_mesh, only : tri_mesh_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use file, only : file_t
  use comm, only : MPI_REAL_PRECISION, NEKO_COMM
  use intersection_detector, only : intersect_detector_t
  use mesh, only : mesh_t
  use stack, only : stack_i4_t, stack_pt_t   
  use global_interpolation, only : global_interpolation_t
  use point, only : point_t
  use math, only : NEKO_EPS         
  use gather_scatter
  use mpi_f08
  use logger
  implicit none
  private

  !> Inverse distance weighting source term.
  type, public, extends(source_term_t) :: idw_source_term_t
     !> Smallest distance between between points and dofs     
     real(kind=dp), private :: ds_min
     type(intersect_detector_t), private :: intersect
     type(global_interpolation_t), private :: global_interp
     type(stack_pt_t), private ::  lagrangian_points
     type(stack_i4_t), private, allocatable :: lag_el(:)
     real(kind=dp), private, allocatable :: xyz(:,:)
     real(kind=rp), private, allocatable :: fu_ib(:)
     real(kind=rp), private, allocatable :: fv_ib(:)
     real(kind=rp), private, allocatable :: fw_ib(:)
     real(kind=rp), private :: pwr_param
     type(field_t), private :: w
     type(gs_t), private :: gs
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => idw_source_term_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => idw_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => idw_source_term_compute
     !> Initialise lagrangian from a boundary mesh
     procedure, private, pass(this) :: init_boundary_mesh => idw_init_boundary_mesh
  end type idw_source_term_t

contains

  subroutine idw_source_term_init_from_json(this, json, fields, coef)
    class(idw_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: start_time, end_time
    type(json_value), pointer :: json_object_list
    type(json_core) :: core
    type(json_file) :: object_settings
    character(len=:), allocatable :: object_type
    integer :: n_regions, i, j, k, e
    character(len=LOG_SIZE) :: log_buf
    real(kind=dp) :: dx_min, dy_min, dz_min, aabb_padding
    type(point_t), pointer :: pt(:)
    type(stack_i4_t) :: overlaps

    ! Mandatory fields for the general source term
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    call neko_log%section('Inverse distance weighting')

    call json_get_or_default(json, "power_parameter", this%pwr_param, 0.5_rp)
    write(log_buf, '(A,f5.2)')  'IDW Power  : ', this%pwr_param
    call neko_log%message(log_buf)
    
    ! Naive apporach to find the smallest distance between two dofs in the mesh
    dx_min = huge(0.0_rp)
    dy_min = huge(0.0_rp)
    dz_min = huge(0.0_rp)
    associate (x => coef%dof%x, y => coef%dof%y, z => coef%dof%z, lx => coef%Xh%lx)
      do e = 1, coef%msh%nelv
         dx_min = min(dx_min, &
              (minval(abs(cshift(x(:,:,:,e), dim = 1, shift = 1) - &
              cshift(x(:,:,:,e), dim = 1, shift = -1))) + &
              minval(abs(cshift(x(:,:,:,e), dim = 1, shift = 1) - &
              cshift(x(:,:,:,e), dim = 2, shift = -1))) + &
              minval(abs(cshift(x(:,:,:,e), dim = 3, shift = 1) - &
              cshift(x(:,:,:,e), dim = 3, shift = -1))))    / 3.0_dp)

         dy_min = min(dy_min, &
              (minval(abs(cshift(y(:,:,:,e), dim = 1, shift = 1) - &
              cshift(y(:,:,:,e), dim = 1, shift = -1))) + &
              minval(abs(cshift(y(:,:,:,e), dim = 1, shift = 1) - &
              cshift(y(:,:,:,e), dim = 2, shift = -1))) + &
              minval(abs(cshift(y(:,:,:,e), dim = 3, shift = 1) - &
              cshift(y(:,:,:,e), dim = 3, shift = -1))))    / 3.0_dp)

         dz_min = min(dz_min, &
              (minval(abs(cshift(z(:,:,:,e), dim = 1, shift = 1) - &
              cshift(z(:,:,:,e), dim = 1, shift = -1))) + &
              minval(abs(cshift(z(:,:,:,e), dim = 1, shift = 1) - &
              cshift(z(:,:,:,e), dim = 2, shift = -1))) + &
              minval(abs(cshift(z(:,:,:,e), dim = 3, shift = 1) - &
              cshift(z(:,:,:,e), dim = 3, shift = -1))))    / 3.0_dp)
        
      end do
    end associate

    this%ds_min = (dx_min + dy_min + dz_min) / 3.0_dp
    call MPI_Allreduce(MPI_IN_PLACE, this%ds_min, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM)
    write(log_buf, '(A,ES13.6)') 'Minimum ds :',  this%ds_min
    call neko_log%message(log_buf)

    aabb_padding = 4 * this%ds_min
    
    call this%intersect%init(coef%msh, aabb_padding)
    call this%lagrangian_points%init()
    
    call json%get('objects', json_object_list)
    call json%info('objects', n_children=n_regions)
    call json%get_core(core)

    if (n_regions .lt. 10) then
       write(log_buf, '(A, I1)') 'Objects    : ', n_regions
    else if (n_regions .ge. 100) then
       write(log_buf, '(A, I2)') 'Objects    : ', n_regions
    else
       write(log_buf, '(A, I3)') 'Objects    : ', n_regions
    end if
    call neko_log%message(log_buf)

    call neko_log%begin()    
    do i = 1, n_regions
       call neko_log%begin()
       call json_extract_item(core, json_object_list, i , object_settings)
       call json_get_or_default(object_settings, 'type', object_type, 'none')

       call neko_log%message('Type       : '// trim(object_type))
       select case (object_type)
       case ('boundary_mesh')
          call this%init_boundary_mesh(object_settings)
       case ('none')
          call neko_error('IDW source term objects require a region type')
       case default
          call neko_error('IDW source term unkown region type')          
       end select
       call neko_log%end()
    end do
    call neko_log%end()


    ! Report total number of lagrangian points generated, this might differ
    ! from the STL sources due to refinement (not implemented yet!)   
    if (this%lagrangian_points%size() .lt. 1e1) then
           call neko_log%message('Type       : '// trim(object_type))
       write(log_buf, '(A, I1)') 'Tot lagpts : ', this%lagrangian_points%size()
    else if (this%lagrangian_points%size() .lt. 1e2) then
       write(log_buf, '(A, I2)') 'Tot lagpts : ', this%lagrangian_points%size()
    else if (this%lagrangian_points%size() .lt. 1e3) then
       write(log_buf, '(A, I3)') 'Tot lagpts : ', this%lagrangian_points%size()
    else if (this%lagrangian_points%size() .lt. 1e4) then
       write(log_buf, '(A, I4)') 'Tot lagpts : ', this%lagrangian_points%size()
    else if (this%lagrangian_points%size() .lt. 1e5) then
       write(log_buf, '(A, I5)') 'Tot lagpts : ', this%lagrangian_points%size()
    else if (this%lagrangian_points%size() .ge. 1e6) then
       write(log_buf, '(A, I6)') 'Tot lagpts : ', this%lagrangian_points%size()
    else if (this%lagrangian_points%size() .ge. 1e7) then
       write(log_buf, '(A, I7)') 'Tot lagpts : ', this%lagrangian_points%size()
    else if (this%lagrangian_points%size() .ge. 1e8) then
       write(log_buf, '(A, I8)') 'Tot lagpts : ', this%lagrangian_points%size()
    end if
    call neko_log%message(log_buf)
    
    allocate(this%xyz(3, this%lagrangian_points%size()))
    allocate(this%fu_ib(this%lagrangian_points%size()))
    allocate(this%fv_ib(this%lagrangian_points%size()))
    allocate(this%fw_ib(this%lagrangian_points%size()))

    pt => this%lagrangian_points%array()
    do i = 1, this%lagrangian_points%size()
       this%xyz(1, i) = pt(i)%x(1)
       this%xyz(2, i) = pt(i)%x(2)
       this%xyz(3, i) = pt(i)%x(3)
    end do
    
    call this%global_interp%init(coef%dof)

    call this%global_interp%find_points_xyz(this%xyz, &
         this%lagrangian_points%size())

    ! Construct list of overlapping elements for each lagrangian particle    
    allocate(this%lag_el(this%lagrangian_points%size()))
    call overlaps%init()
    pt => this%lagrangian_points%array()
    do i = 1, this%lagrangian_points%size()
       call this%lag_el(i)%init()
       call this%intersect%overlap(pt(i), overlaps)
       do while(.not. overlaps%is_empty())
          e = overlaps%pop()
          call this%lag_el(i)%push(e)
       end do
    end do
    call overlaps%free()

    ! Construct weight field
    call this%w%init(coef%dof, "ib_weight")

    call this%gs%init(coef%dof)

    call idw_compute_weight(this%w, this%lagrangian_points, this%lag_el, &
         coef%dof%x, coef%dof%y, coef%dof%z, aabb_padding, &
         this%pwr_param, this%gs, coef%Xh%lx,coef%msh%nelv)

  end subroutine idw_source_term_init_from_json

  subroutine idw_source_term_free(this)
    class(idw_source_term_t), intent(inout) :: this
    integer :: i
    
    call this%free_base()

    call this%lagrangian_points%free()

    if (allocated(this%xyz)) then
       deallocate(this%xyz)
    end if

    if (allocated(this%fu_ib)) then
       deallocate(this%fu_ib)
    end if

    if (allocated(this%fv_ib)) then
       deallocate(this%fv_ib)
    end if

    if (allocated(this%fw_ib)) then
       deallocate(this%fw_ib)
    end if

    if (allocated(this%lag_el)) then
       do i = 1, size(this%lag_el)
          call this%lag_el(i)%free()
       end do
       deallocate(this%lag_el)
    end if

    call this%gs%free()

  end subroutine idw_source_term_free

  subroutine idw_source_term_compute(this, t, tstep)
    class(idw_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u, v, w, fu, fv, fw
    integer :: n

    n = this%fields%item_size(1)

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)
    
    associate(global_interp => this%global_interp, &
         fu_ib => this%fu_ib, fv_ib => this%fv_ib, fw_ib => this%fw_ib)

      call global_interp%evaluate(fu_ib, u%x)
      call global_interp%evaluate(fv_ib, v%x)
      call global_interp%evaluate(fw_ib, w%x)

    end associate
    
  end subroutine idw_source_term_compute

  subroutine idw_init_boundary_mesh(this, json)
    class(idw_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(file_t) :: mesh_file
    type(tri_mesh_t) :: boundary_mesh
    character(len=:), allocatable :: mesh_file_name
    character(len=LOG_SIZE) :: log_buf
    integer :: i, el_idx
    type(stack_i4_t) :: overlaps

    call json_get(json, 'name', mesh_file_name)
    mesh_file = file_t(mesh_file_name)
    call neko_log%message('Filename   : '// trim(mesh_file_name))
    call mesh_file%read(boundary_mesh)
    if (boundary_mesh%mpts .lt. 1e1) then
       write(log_buf, '(A, I1)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e2) then
       write(log_buf, '(A, I2)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e3) then
       write(log_buf, '(A, I3)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e4) then
       write(log_buf, '(A, I4)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e5) then
       write(log_buf, '(A, I5)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .ge. 1e6) then
       write(log_buf, '(A, I6)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .ge. 1e7) then
       write(log_buf, '(A, I7)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .ge. 1e8) then
       write(log_buf, '(A, I8)') ' `-Points  : ', boundary_mesh%mpts
    end if
    call neko_log%message(log_buf)

    if (boundary_mesh%nelv .eq. 0) then
       call neko_error('No elements in the boundary mesh')
    end if

    call overlaps%init()
    
    do i = 1, boundary_mesh%mpts
       call this%intersect%overlap(boundary_mesh%points(i), overlaps)

       if (.not. overlaps%is_empty()) then
          call this%lagrangian_points%push(boundary_mesh%points(i))
       end if

       do while (.not. overlaps%is_empty())
          el_idx = overlaps%pop()
       end do
       call overlaps%clear()

    end do
    
    call boundary_mesh%free()
    
  end subroutine idw_init_boundary_mesh

  !> Compute IB weight field
  subroutine idw_compute_weight(w, lag_pts, lag_el, x, y, z, rmax, p, gs, lx, ne)
    type(field_t), intent(inout) :: w
    type(stack_pt_t), intent(inout) :: lag_pts
    type(stack_i4_t), intent(inout) :: lag_el(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne) :: x, y, z
    real(kind=rp), intent(inout) :: rmax, p
    type(gs_t), intent(inout) :: gs
    integer, pointer :: el(:)
    type(point_t), pointer :: pt(:)
    integer :: i, j, k, l, e, ee
    real(kind=rp) :: r

    w%x = 0.0_rp

    pt => lag_pts%array()
    do i = 1, lag_pts%size()
       el => lag_el(i)%array()
       do ee = 1, lag_el(i)%size()
          e = el(ee)
          do l = 1, lx 
             do k = 1, lx
                do j = 1, lx
                   r = sqrt((x(j,k,l,e) - pt(i)%x(1))**2 &
                        + (y(j,k,l,e) - pt(i)%x(2))**2 &
                        + (z(j,k,l,e) - pt(i)%x(3))**2)
                   
                   w%x(j, k, l, e) = w%x(j, k, l, e) &
                        + inv_dist_weight(r, rmax, p)
                end do
             end do
          end do
       end do
    end do

    call gs%op(w, GS_OP_ADD)
    
  end subroutine idw_compute_weight
  
  !> Inverse distance weighting coefficient
  !! @param r Radial distance to Lagrangian point.
  !! @param rmax Radial distance for sphere of influence.
  pure function inv_dist_weight(r, rmax, p) result(idw)    
    real(kind=rp), intent(in) :: r
    real(kind=rp), intent(in) :: rmax
    real(kind=rp), intent(in) :: p
    real(kind=rp) :: idw
    
    if(r .ge. rmax) then
       idw = 0.0_rp
    else
       idw = ((rmax-r)/(rmax * r + NEKO_EPS))**p
    end if
    
  end function inv_dist_weight

end module idw_source_term
