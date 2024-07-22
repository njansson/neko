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
  use mpi_f08
  use logger
  ! use point, only : point_t
  ! use global_interpolation, only : global_interpolation_t
  
  ! use uset, only : uset_i4_t
  ! use gather_scatter
  ! use file
  ! use stack
  ! use math
  implicit none
  private

  !> Inverse distance weighting source term.
  type, public, extends(source_term_t) :: idw_source_term_t
!     type(tri_mesh_t) :: model
!     type(global_interpolation_t) :: global_interp
!     real(kind=rp), allocatable :: xyz(:,:)
!     real(kind=rp), allocatable :: F_ib(:,:)
!     type(stack_i4_t), allocatable :: neigh_el(:)
!     type(field_t) :: w
!     real(kind=rp), allocatable :: rmax(:)
!     real(kind=rp) :: pwr_param
!     type(field_list_t) :: sampled_fields
!     logical :: stationary = .true.
!     logical :: w_computed = .false.
     !     type(gs_t) :: gs
     ! Smallest distance between between points and dofs     
     real(kind=dp) :: ds_min
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => idw_source_term_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => idw_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => idw_source_term_compute
     
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
    real(kind=dp) :: dx_min, dy_min, dz_min
    type(field_t), pointer :: u

    ! Mandatory fields for the general source term
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    call neko_log%section('Inverse distance weighting')
    
    ! Naive apporach to find the smallest distance between two dofs in the mesh
    u => neko_field_registry%get_field('u')
    dx_min = huge(0.0_rp)
    dy_min = huge(0.0_rp)
    dz_min = huge(0.0_rp)
    associate (x => u%dof%x, y=>u%dof%y, z=>u%dof%z, lx => u%Xh%lx)
      do e = 1, u%msh%nelv
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

  end subroutine idw_source_term_init_from_json

  subroutine idw_source_term_free(this)
    class(idw_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine idw_source_term_free

  subroutine idw_source_term_compute(this, t, tstep)
    class(idw_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u, v, w
    
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

  end subroutine idw_source_term_compute

  subroutine idw_init_boundary_mesh(this, json)
    class(idw_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(file_t) :: mesh_file
    type(tri_mesh_t) :: boundary_mesh
    character(len=:), allocatable :: mesh_file_name
    character(len=LOG_SIZE) :: log_buf
    
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
  
    
  end subroutine idw_init_boundary_mesh

end module idw_source_term
