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
     real(kind=dp) :: ds_min
     type(intersect_detector_t) :: intersect
     type(global_interpolation_t) :: global_interp
     type(point_t), allocatable :: lag_pts(:)
     type(stack_i4_t), allocatable :: lag_el(:)
     real(kind=dp), allocatable :: xyz(:,:)
     real(kind=rp), allocatable :: fu_ib(:)
     real(kind=rp), allocatable :: fv_ib(:)
     real(kind=rp), allocatable :: fw_ib(:)
     real(kind=rp)  :: pwr_param
     real(kind=rp) :: rmax
     type(field_t)  :: w
     type(field_t)  :: ds 
     type(gs_t)  :: gs
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
    real(kind=dp) :: aabb_padding,dx_max, dy_max, dz_max, ds_max, ds_min
    type(stack_i4_t) :: overlaps
    type(stack_pt_t) ::  lagrangian_points
         

    type(file_t) :: wf, dsf

    ! Mandatory fields for the general source term
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    call neko_log%section('Inverse distance weighting')

    call json_get_or_default(json, "rmax", this%rmax, 1.0_rp)
    write(log_buf, '(A,f5.2)')  'Rmax       : ', this%rmax
    call neko_log%message(log_buf)
    call json_get_or_default(json, "power_parameter", this%pwr_param, 0.5_rp)
    write(log_buf, '(A,f5.2)')  'IDW Power  : ', this%pwr_param
    call neko_log%message(log_buf)
    
    ! Naive apporach to find the smallest distance between two dofs in the mesh

    ds_min = huge(0.0_rp)
    ds_max = -huge(0.0_rp)

    call this%ds%init(coef%dof)

    associate (x => coef%dof%x, y => coef%dof%y, z => coef%dof%z, &
         lx => coef%Xh%lx, ds => this%ds%x)

      do e = 1, coef%msh%nelv
         do k = 2, lx-1
            do j = 2, lx-1
               do i = 2, lx-1
                  dx_max = max(abs(x(i,j,k,e) - x(i+1,j,k,e)), &
                       abs(x(i-1,j,k,e) - x(i,j,k,e)), &
                       abs(x(i,j,k,e) - x(i,j+1,k,e)), &
                       abs(x(i,j-1,k,e) - x(i,j,k,e)), &
                       abs(x(i,j,k,e) - x(i,j,k+1,e)), &
                       abs(x(i,j,k-1,e) - x(i,j,k,e)), &
                       abs(x(i-1, j-1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j-1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j+1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j+1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j-1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j-1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j+1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j+1, k-1, e) - x(i,j,k,e)))

                  dy_max = max(abs(y(i,j,k,e) - y(i+1,j,k,e)), &
                       abs(y(i-1,j,k,e) - y(i,j,k,e)), &
                       abs(y(i,j,k,e) - y(i,j+1,k,e)), &
                       abs(y(i,j-1,k,e) - y(i,j,k,e)), &
                       abs(y(i,j,k,e) - y(i,j,k+1,e)), &
                       abs(y(i,j,k-1,e) - y(i,j,k,e)), &
                       abs(y(i-1, j-1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j-1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j+1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j+1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j-1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j-1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j+1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j+1, k-1, e) - y(i,j,k,e)))


                  dz_max = max(abs(z(i,j,k,e) - z(i+1,j,k,e)), &
                       abs(z(i-1,j,k,e) - z(i,j,k,e)), &
                       abs(z(i,j,k,e) - z(i,j+1,k,e)), &
                       abs(z(i,j-1,k,e) - z(i,j,k,e)), &
                       abs(z(i,j,k,e) - z(i,j,k+1,e)), &
                       abs(z(i,j,k-1,e) - z(i,j,k,e)), &
                       abs(z(i-1, j-1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j-1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j+1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j+1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j-1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j-1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j+1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j+1, k-1, e) - z(i,j,k,e)))
                  ds(i,j,k,e) = (dx_max + dy_max + dz_max) / 3.0_rp
               end do
            end do
         end do

         
         ds(1,:,:,e) = ds(2,:,:,e)
         ds(lx,:,:,e) = ds(lx-1,:,:,e)

         ds(:,1,:,e) = ds(:,2,:,e)
         ds(:,lx,:,e) = ds(:,lx-1,:,e)

         ds(:,:,1,e) = ds(:,:,2,e)
         ds(:,:,lx,e) = ds(:,:,lx-1,e)


         ds_max = max(ds_max, maxval(ds(:,:,:,e)))
         ds_min = min(ds_min, minval(ds(:,:,:,e)))
                 
      end do
    end associate

    this%ds_min = ds_min
    
    call MPI_Allreduce(MPI_IN_PLACE, this%ds_min, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM)
    write(log_buf, '(A,ES13.6)') 'Minimum ds :',  this%ds_min
    call neko_log%message(log_buf)

    aabb_padding = 1 * ds_max
    
    call this%intersect%init(coef%msh, aabb_padding)
    call lagrangian_points%init()
    
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
          call this%init_boundary_mesh(lagrangian_points, object_settings)
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
    if (lagrangian_points%size() .lt. 1e1) then
           call neko_log%message('Type       : '// trim(object_type))
       write(log_buf, '(A, I1)') 'Tot lagpts : ', lagrangian_points%size()
    else if (lagrangian_points%size() .lt. 1e2) then
       write(log_buf, '(A, I2)') 'Tot lagpts : ', lagrangian_points%size()
    else if (lagrangian_points%size() .lt. 1e3) then
       write(log_buf, '(A, I3)') 'Tot lagpts : ', lagrangian_points%size()
    else if (lagrangian_points%size() .lt. 1e4) then
       write(log_buf, '(A, I4)') 'Tot lagpts : ', lagrangian_points%size()
    else if (lagrangian_points%size() .lt. 1e5) then
       write(log_buf, '(A, I5)') 'Tot lagpts : ', lagrangian_points%size()
    else if (lagrangian_points%size() .ge. 1e6) then
       write(log_buf, '(A, I6)') 'Tot lagpts : ', lagrangian_points%size()
    else if (lagrangian_points%size() .ge. 1e7) then
       write(log_buf, '(A, I7)') 'Tot lagpts : ', lagrangian_points%size()
    else if (lagrangian_points%size() .ge. 1e8) then
       write(log_buf, '(A, I8)') 'Tot lagpts : ', lagrangian_points%size()
    end if
    call neko_log%message(log_buf)
    
    allocate(this%xyz(3, lagrangian_points%size()))
    allocate(this%fu_ib(lagrangian_points%size()))
    allocate(this%fv_ib(lagrangian_points%size()))
    allocate(this%fw_ib(lagrangian_points%size()))
    allocate(this%lag_pts(lagrangian_points%size()))
    
    select type(pt => lagrangian_points%data)
    type is (point_t)
       do i = 1, lagrangian_points%size()
          this%xyz(1, i) = pt(i)%x(1)
          this%xyz(2, i) = pt(i)%x(2)
          this%xyz(3, i) = pt(i)%x(3)
          this%lag_pts(i) = pt(i)
       end do
    end select

    
    call this%global_interp%init(coef%dof)

    call this%global_interp%find_points_xyz(this%xyz, &
         lagrangian_points%size())

    ! Construct list of overlapping elements for each lagrangian particle    
    allocate(this%lag_el(lagrangian_points%size()))
    call overlaps%init()
    do i = 1, size(this%lag_pts)
       call this%lag_el(i)%init()
       call this%intersect%overlap(this%lag_pts(i), overlaps)
       do while(.not. overlaps%is_empty())
          e = overlaps%pop()
          call this%lag_el(i)%push(e)
       end do
    end do
    call overlaps%free()

    ! Construct weight field
    call this%w%init(coef%dof, "ib_weight")

    call this%gs%init(coef%dof)

    this%ds%x = this%ds%x * coef%mult    
    call this%gs%op(this%ds, GS_OP_ADD)


    call idw_compute_weight(this%w, this%lag_pts, this%lag_el, &
         coef%dof%x, coef%dof%y, coef%dof%z, this%ds%x, this%rmax, &
         this%pwr_param, this%gs, coef%Xh%lx,coef%msh%nelv)

    this%w%x = this%w%x * coef%mult

    call this%gs%op(this%w, GS_OP_ADD)
    
    wf = file_t("w.fld")
    call wf%write(this%w)

    dsf = file_t("ds.fld")
    call dsf%write(this%ds)
    
    call lagrangian_points%free()
    
  end subroutine idw_source_term_init_from_json

  subroutine idw_source_term_free(this)
    class(idw_source_term_t), intent(inout) :: this
    integer :: i
    
    call this%free_base()



    call this%ds%free()

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

    if (allocated(this%lag_pts)) then
       deallocate(this%lag_pts)
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
    integer :: i, j, k, l, e, ee, n
    real(kind=rp) :: r, idw, dt
    n = this%fields%item_size(1)

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)

    !> @todo Change this once we have variable time-stepping
    dt = t / tstep

    associate(global_interp => this%global_interp, &
         fu_ib => this%fu_ib, fv_ib => this%fv_ib, fw_ib => this%fw_ib, &
         lag_pts => this%lag_pts, &
         ds => this%ds%x, c => this%w%x, x => this%w%dof%x, &
         y => this%w%dof%y, z => this%w%dof%z, lx => this%w%Xh%lx)


      fu_ib = 0.0_rp
      fv_ib = 0.0_rp
      fw_ib = 0.0_rp
      
      call global_interp%evaluate(fu_ib, u%x)
      call global_interp%evaluate(fv_ib, v%x)
      call global_interp%evaluate(fw_ib, w%x)

      do i = 1, size(this%lag_pts)
         select type (el => this%lag_el(i)%data)
         type is (integer)
            do ee = 1, this%lag_el(i)%size()
               e = el(ee)
               do l = 1, lx 
                  do k = 1, lx
                     do j = 1, lx
                        if (abs(this%w%x(j,k,l,e)) .gt. 1e-4_rp) then
                           r = sqrt((x(j,k,l,e) - lag_pts(i)%x(1))**2 &
                                + (y(j,k,l,e) - lag_pts(i)%x(2))**2 &
                                + (z(j,k,l,e) - lag_pts(i)%x(3))**2)
                           r = r / ds(j,k,l,e)
                           idw = inv_dist_weight(r, this%rmax, this%pwr_param)

                           fu%x(j,k,l,e) = fu%x(j,k,l,e) &
                                + (-fu_ib(i) * idw) / (this%w%x(j,k,l,e) * dt)
                           
                           fv%x(j,k,l,e) = fv%x(j,k,l,e) &
                                + (-fv_ib(i) * idw) / (this%w%x(j,k,l,e) * dt)
                           
                           fw%x(j,k,l,e) = fw%x(j,k,l,e) &
                                + (-fw_ib(i) * idw) / (this%w%x(j,k,l,e) * dt)
                           

                        end if
                     end do
                  end do
               end do
            end do
         end select
      end do
      

      fu%x = fu%x * this%coef%mult
      fv%x = fv%x * this%coef%mult
      fw%x = fw%x * this%coef%mult

      call this%gs%op(fu, GS_OP_ADD)
      call this%gs%op(fv, GS_OP_ADD)
      call this%gs%op(fw, GS_OP_ADD)

    end associate
    
  end subroutine idw_source_term_compute

  subroutine idw_init_boundary_mesh(this, lag_pts, json)
    class(idw_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(stack_pt_t), intent(inout) :: lag_pts
    type(file_t) :: mesh_file
    type(tri_mesh_t) :: boundary_mesh
    character(len=:), allocatable :: mesh_file_name
    character(len=LOG_SIZE) :: log_buf
    integer :: i, j, el_idx
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

!       if (overlaps%size() .gt. 0) then
          call lag_pts%push(boundary_mesh%points(i))

          do while (.not. overlaps%is_empty())
             el_idx = overlaps%pop()
          end do
 !      end if
       call overlaps%clear()

    end do
    
    call boundary_mesh%free()
    
  end subroutine idw_init_boundary_mesh

  !> Compute IB weight field
  subroutine idw_compute_weight(w, lag_pts, lag_el, x, y, z, ds, rmax, p, gs, lx, ne)
    type(field_t), intent(inout) :: w
    type(point_t), allocatable, intent(inout) :: lag_pts(:)
    type(stack_i4_t), allocatable, intent(inout) :: lag_el(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne) :: x, y, z, ds
    real(kind=rp), intent(inout) :: p
    real(kind=rp), intent(inout) :: rmax
    type(gs_t), intent(inout) :: gs
    integer :: i, j, k, l, e, ee
    real(kind=rp) :: r

    w%x = 0.0_rp

    do i = 1, size(lag_pts)
       select type(el => lag_el(i)%data)
       type is (integer)
          do ee = 1, lag_el(i)%size()
             e = el(ee)
             do l = 1, lx 
                do k = 1, lx
                   do j = 1, lx
                      r = sqrt((x(j,k,l,e) - lag_pts(i)%x(1))**2 &
                           + (y(j,k,l,e) - lag_pts(i)%x(2))**2 &
                           + (z(j,k,l,e) - lag_pts(i)%x(3))**2)
                      r = r / ds(j,k,l,e)
                      w%x(j, k, l, e) = w%x(j, k, l, e) &
                           + inv_dist_weight(r, rmax, p)
                   end do
                end do
             end do
          end do
       end select
    end do
    
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
