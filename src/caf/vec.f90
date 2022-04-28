! Copyright (c) 2022, The Neko Authors
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
!> Defines a coarray vector
module vec
  use num_types
  use utils
  use tuple
  use stack
  implicit none
  private
  
  type :: vec_data_t
     real(kind=dp), allocatable :: x(:) !< Local part of global vector
  end type vec_data_t

  type :: vec_oimg_t
     type(tuple_i4r8_t), allocatable :: data(:)
     integer :: top_
     integer :: size_
  end type vec_oimg_t

  type, public :: vec_t
     type(vec_data_t), allocatable :: X[:]  !< Coarray of local vectors
     integer :: n                           !< Size of local vector
     integer :: m                           !< Size of global vector
     integer :: range(2)                    !< Ownership range
     integer, allocatable :: glb_range(:)   !< Start offset for each image
     type(vec_oimg_t), allocatable :: oimg(:)[:] !< Off image buffer
     type(stack_i4_t), private :: neigh_img   !< Neigh. images
     type(tuple_i4r8_t), private, allocatable :: tmp_buf(:)
     contains
     procedure, pass(this) :: init => vec_init
     procedure, pass(this) :: free => vec_free
     procedure, pass(this) :: zero => vec_zero
     procedure, pass(this) :: finalize => vec_finalize
     procedure, pass(this) :: set_scalar => vec_set_scalar
     procedure, pass(this) :: set_block => vec_set_block
     generic :: set => set_scalar, set_block
     procedure, pass(this) :: add_scalar => vec_add_scalar
     procedure, pass(this) :: add_block => vec_add_block
     generic :: add => add_scalar, add_block
  end type vec_t
  
contains

  !> Initialise a vector with @a n local entries
  subroutine vec_init(this, n)
    class(vec_t), intent(inout) :: this
    integer, intent(in) :: n
    integer, allocatable :: tmp[:]
    integer :: acc, i 
    
    call vec_free(this)

    this%n = n
    this%m = n
    call co_sum(this%m)

    allocate(tmp[*])
    allocate(this%X[*])
    sync all
    allocate(this%X%x(n))

    tmp[this_image()] = this%n
    sync all
    
    acc = 0
    do i = 1, this_image()
       acc = acc + tmp[i]
    end do
    sync all
    
    this%range(1) = acc - this%n + 1
    this%range(2) = this%range(1) + this%n
    
    tmp[this_image()] = this%range(1)
    sync memory
    
    allocate(this%glb_range(num_images() + 1))
    
    do i = 1, num_images()
       this%glb_range(i) = tmp[i]
    end do
    this%glb_range(num_images() + 1) = this%m

    allocate(this%oimg(num_images())[*])
    sync all
    
    do i = 1, num_images()
       allocate(this%oimg(i)%data(32))
       this%oimg(i)%size_ = 32
       this%oimg(i)%top_ = 0
    end do
    sync all
    deallocate(tmp)
    
  end subroutine vec_init

  !> Destroy a vector
  subroutine vec_free(this)
    class(vec_t), intent(inout) :: this
    integer :: i
    
    this%n = 0
    this%m = 0

    if (allocated(this%X)) then
       if (allocated(this%X%x)) then
          deallocate(this%X%x)
       end if
       deallocate(this%X)
    end if

    if (allocated(this%glb_range)) then
       deallocate(this%glb_range)
    end if

    if (allocated(this%oimg)) then
       do i = 1, num_images()
          if (allocated(this%oimg(i)%data)) then
             deallocate(this%oimg(i)%data)
          end if
       end do
       deallocate(this%oimg)
    end if

    call this%neigh_img%free()

    if (allocated(this%tmp_buf)) then
       deallocate(this%tmp_buf)
    end if
    
  end subroutine vec_free

  !> Zero a vector
  subroutine vec_zero(this)
    class(vec_t), intent(inout) :: this

    this%X%x = 0.0_rp
    sync memory
    
  end subroutine vec_zero

  subroutine vec_finalize(this)
    class(vec_t), intent(inout) :: this
    type(tuple_i4r8_t), allocatable :: tmp(:)
    integer :: i, j, src, osize, offset

    sync all
    
    if (.not. allocated(this%tmp_buf)) then
       call this%neigh_img%init()
       do i = 1, num_images()
          src = mod((this_image() + i), num_images()) + 1

          if (src .eq. this_image()) then
             cycle
          end if

          ! obtain remote size
          osize = this%oimg(this_image())[src]%top_

          if (osize .gt. 0) then
               
             if (.not. allocated(this%tmp_buf)) then
                allocate(this%tmp_buf(osize))
             else if (osize .gt. size(this%tmp_buf)) then
                allocate(tmp(osize))
                tmp(1:size(this%tmp_buf)) = this%tmp_buf
                call move_alloc(tmp, this%tmp_buf)
             end if

             ! Only neighbours with non-zero data (assume static pattern)
             call this%neigh_img%push(src)

             ! Copy data from remote image
             this%tmp_buf(1:osize) = this%oimg(this_image())[src]%data(1:osize)
             do j = 1, osize
                offset = this%tmp_buf(j)%x - this%range(1)
                this%X%x(offset) = this%X%x(offset) + this%tmp_buf(j)%y
             end do

          end if
       end do
    else
       select type(neighp=>this%neigh_img%data)
       type is (integer)
          do i = 1, this%neigh_img%top_
             src = neighp(i)
             osize = this%oimg(this_image())[src]%top_
             if (osize .gt. 0) then
                this%tmp_buf(1:osize) = &
                     this%oimg(this_image())[src]%data(1:osize)
                do j = 1, osize
                   offset = this%tmp_buf(j)%x - this%range(1)
                   this%X%x(offset) = this%X%x(offset) + this%tmp_buf(j)%y
                end do
             end if
          end do
       end select
    end if

    sync all

    ! Reset off-image buffer
    do i = 0, num_images() - 1
       this%oimg(i)%top_ = 0
    end do
    
  end subroutine vec_finalize

  subroutine vec_set_scalar(this, s, j)
    class(vec_t), intent(inout) :: this
    real(kind=rp), intent(in) :: s
    integer, intent(in) :: j
    integer :: owner
    type(tuple_i4r8_t) :: vec_tuple    
    
    if (j .ge. this%range(1) .and. j .lt. this%range(2)) then
       this%X%x(j - this%range(1)) = s
    else
       vec_tuple%x = j
       vec_tuple%y = s
       owner = bsearch(this%glb_range, j, 1, num_images())
       call vec_oimg_push(this%oimg(owner), vec_tuple)
    end if
  end subroutine vec_set_scalar

  subroutine vec_set_block(this, b, j, n)
    class(vec_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: b(n)
    integer, intent(in) :: j(n)
    integer :: i, owner
    type(tuple_i4r8_t) :: vec_tuple

    do i = 1, n
       if (j(i) .ge. this%range(1) .and. j(i) .lt. this%range(2)) then
          this%X%x(j(i) - this%range(1)) = b(i)
       else
          vec_tuple%x = j(i)
          vec_tuple%y = b(i)
          owner = bsearch(this%glb_range, j(i), 1, num_images())
          call vec_oimg_push(this%oimg(owner), vec_tuple)
       end if
    end do
    
  end subroutine vec_set_block

  subroutine vec_add_scalar(this, s, j)
    class(vec_t), intent(inout) :: this
    real(kind=rp), intent(in) :: s
    integer, intent(in) :: j
    integer :: offset, owner
    type(tuple_i4r8_t) :: vec_tuple

    if (j .ge. this%range(1) .and. j .lt. this%range(2)) then
       offset = j - this%range(1) 
       this%X%x(offset) = this%X%x(offset) + s
    else
       vec_tuple%x = j
       vec_tuple%y = s
       owner = bsearch(this%glb_range, j, 1, num_images())
       call vec_oimg_push(this%oimg(owner), vec_tuple)
    end if
  end subroutine vec_add_scalar
  
  subroutine vec_add_block(this, b, j, n)
    class(vec_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: b(n)
    integer, intent(in) :: j(n)
    integer :: i, offset, owner
    type(tuple_i4r8_t) :: vec_tuple
    
    do i = 1, n
       if (j(i) .ge. this%range(1) .and. j(i) .lt. this%range(2)) then
          offset = j(i) - this%range(1)
          this%X%x(offset) = this%X%x(offset) + b(i)
       else
          vec_tuple%x = j(i)
          vec_tuple%y = b(i)
          owner = bsearch(this%glb_range, j(i), 1, num_images())
          call vec_oimg_push(this%oimg(owner), vec_tuple)
       end if
    end do
    
  end subroutine vec_add_block

  subroutine vec_oimg_push(oimg, vec_tuple)
    type(vec_oimg_t), intent(inout) ::oimg
    type(tuple_i4r8_t), intent(inout) :: vec_tuple
    type(tuple_i4r8_t), allocatable :: tmp(:)
    integer :: i

    if (oimg%top_ .eq. oimg%size_) then
       oimg%size_ = ishft(oimg%size_, 1)
       allocate(tmp(oimg%size_))
       tmp(1:oimg%top_) = oimg%data
       call move_alloc(tmp, oimg%data)
    end if

    oimg%top_ = oimg%top_ + 1
    oimg%data(oimg%top_) = vec_tuple
    
  end subroutine vec_oimg_push
  
end module vec
