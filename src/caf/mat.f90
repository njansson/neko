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
!> Defines a coarray sparse matrix
module mat
  use num_types
  use utils
  use tuple
  use stack
  implicit none
  private

  type :: mat_csr_t
     integer, allocatable :: col(:)
     integer, allocatable :: rpt(:)
     real(kind=dp), allocatable :: val(:)
     integer :: nz
  end type mat_csr_t

  type :: mat_oimg_t
     type(tuple_2i4r8_t), allocatable :: data(:)
     integer :: top_
     integer :: size_ 
  end type mat_oimg_t

  type, public :: mat_t
     type(mat_csr_t), allocatable :: D[:] !< Diagonal part 
     type(mat_csr_t), allocatable :: O[:] !< Off-diagonal part
     integer :: m                         !< Local number of rows
     integer :: n                         !< Local number of columns
     integer :: m_glb                     !< Global number of rows     
     integer :: n_glb                     !< Global number of columns
     integer :: range(2)                  !< Ownership range
     integer, allocatable :: glb_range(:) !< Start offset for each image
     type(mat_oimg_t), allocatable :: oimg(:)[:] !< Off image buffer
     type(stack_i4_t), private :: neigh_img      !< Neigh. images
     type(stack_i4r8t2_t), private, allocatable :: rs(:) !< Row-stacks
     type(tuple_2i4r8_t), private, allocatable :: tmp_buf(:)
     logical, private :: assembled = .false.
   contains
     procedure, pass(A) :: init => mat_init
     procedure, pass(A) :: free => mat_free
     procedure, pass(A) :: finalize => mat_finalize
     procedure, pass(A) :: add_scalar => mat_add_scalar
     procedure, pass(A) :: add_block => mat_add_block
     generic :: add => add_scalar, add_block
  end type mat_t

contains

  !> Initialise a matrix with @a m local rows and @a n local columns
  subroutine mat_init(A, m, n)
    class(mat_t), intent(inout) :: A
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, allocatable :: tmp[:]
    integer :: acc, i
    
    call mat_free(A)

    A%m = m
    A%m_glb = m
    call co_sum(A%m_glb)


    A%n = n
    A%n_glb = n
    call co_sum(A%n_glb)

    allocate(tmp[*])
    allocate(A%D[*])
    allocate(A%O[*])
    sync all

    tmp[this_image()] = A%m
    sync all

    acc = 0
    do i = 1, this_image()
       acc = acc + tmp[i]
    end do
    sync all

    A%range(1) = acc - A%m + 1
    A%range(2) = A%range(1) + A%m

    tmp[this_image()] = A%range(1)
    sync memory

    allocate(A%glb_range(num_images() + 1))

    do i = 1, num_images()
       A%glb_range(i) = tmp[i]
    end do
    A%glb_range(num_images() + 1) = A%m_glb

    allocate(A%oimg(num_images())[*])
    sync all

    do i = 1, num_images()
       allocate(A%oimg(i)%data(32))
       A%oimg(i)%size_ = 32
       A%oimg(i)%top_ = 0
    end do
    
    sync all
    deallocate(tmp)

    allocate(A%rs(m))

    do i = 1, m
       call A%rs(i)%init()
    end do
    
  end subroutine mat_init

  !> Destroy a matrix
  subroutine mat_free(A)
    class(mat_t), intent(inout) :: A
    integer :: i
    
    A%m = 0
    A%n = 0
    
    A%m_glb = 0
    A%n_glb = 0

    if (allocated(A%D)) then       
       call mat_csr_free(A%D)
    end if

    if (allocated(A%O)) then
       call mat_csr_free(A%O)
    end if

    if (allocated(A%glb_range)) then
       deallocate(A%glb_range)
    end if

    if (allocated(A%oimg)) then
       do i = 1, num_images()
          if (allocated(A%oimg(i)%data)) then
             deallocate(A%oimg(i)%data)
          end if
       end do
       deallocate(A%oimg)
    end if

    call A%neigh_img%free()
    
    if (allocated(A%rs)) then
       deallocate(A%rs)
    end if

    if (allocated(A%tmp_buf)) then
       deallocate(A%tmp_buf)
    end if
    
  end subroutine mat_free

  !> Destroy a CSR structure
  subroutine mat_csr_free(csr)
    type(mat_csr_t), intent(inout) :: csr

    if (allocated(csr%col)) then
       deallocate(csr%col)
    end if

    if (allocated(csr%rpt)) then
       deallocate(csr%rpt)
    end if

    if (allocated(csr%val)) then
       deallocate(csr%val)
    end if
    
  end subroutine mat_csr_free

  subroutine mat_finalize(A)
    class(mat_t), intent(inout) :: A
    type(tuple_2i4r8_t), allocatable :: tmp(:)
    integer :: i, j, k, src, osize
    type(tuple_i4r8_t) :: col_data
    integer :: o_nz, d_nz, o_v, d_v

    sync all

    if (.not. allocated(A%tmp_buf)) then
       ! Dummy allocation to mark initial finalisation
       allocate(A%tmp_buf(16))
       call A%neigh_img%init()
       do i = 1, num_images()
          src = mod((this_image() + i), num_images()) + 1

          if (src .eq. this_image()) then
             cycle
          end if

          ! obtain remote size
          osize = A%oimg(this_image())[src]%top_

          if (osize .gt. 0) then

             if (.not. allocated(A%tmp_buf)) then
                allocate(A%tmp_buf(osize))
             else if (osize .gt. size(A%tmp_buf)) then
                allocate(tmp(osize))
                tmp(1:size(A%tmp_buf)) = A%tmp_buf
                call move_alloc(tmp, A%tmp_buf)
             end if
             
             ! Only neighbours with non-zero data (assume static pattern)
             call A%neigh_img%push(src)

             ! Copy data from remote image
             A%tmp_buf(1:osize) = A%oimg(this_image())[src]%data(1:osize)
             do j = 1, osize
                call A%add(A%tmp_buf(j)%x, A%tmp_buf(j)%y, A%tmp_buf(j)%z)
             end do
          end if
       end do

       o_nz = 0
       d_nz = 0
       
       ! Sort rows using insertion sort and count non-zeros       
       do i = 1, A%m          
          select type (ep => A%rs(i)%data)
          type is (tuple_i4r8_t)
             do j = 2, A%rs(i)%top_
                col_data = ep(j)
                k = j - 1
                do while(k .ge. 1 .and. ep(k)%x .gt. col_data%x)
                   ep(k+1) = ep(k)
                   k = k - 1
                end do
                ep(k+1) = col_data
             end do

             do j = 1, A%rs(i)%top_
                if (ep(j)%x .lt. A%range(1) .or. ep(j)%x .ge. A%range(2)) then
                   o_nz = o_nz + 1
                else
                   d_nz = d_nz + 1
                end if
             end do
          end select
       end do

       ! Construct diagonal and off-diagonal CSR blocks
       allocate(A%O%val(o_nz), A%O%col(o_nz), A%O%rpt(A%m + 1))
       allocate(A%D%val(d_nz), A%D%col(d_nz), A%D%rpt(A%m + 1))

       ! counters
       o_v = 1
       d_v = 1
              
       do i = 1, A%m

          A%O%rpt(i) = o_v
          A%D%rpt(i) = d_v
                       
          select type (ep => A%rs(i)%data)
          type is (tuple_i4r8_t)
             do j = 1, A%rs(i)%top_
                if (ep(j)%x .lt. A%range(1) .or. ep(j)%x .ge. A%range(2)) then
                   A%O%col(o_v) = ep(j)%x
                   A%O%val(o_v) = ep(j)%y
                   o_v = o_v + 1                            
                else
                   A%D%col(d_v) = ep(j)%x
                   A%D%val(d_v) = ep(j)%y
                   d_v = d_v + 1                                      
                end if
             end do
          end select
       end do

       A%O%rpt(A%m + 1) = o_nz
       A%D%rpt(A%m + 1) = d_nz

       A%O%nz = o_nz
       A%D%nz = d_nz

    else
       select type(neighp => A%neigh_img%data)
       type is(integer)
          do i = 1, A%neigh_img%top_
             src = neighp(i)
             osize = A%oimg(this_image())[src]%top_
             if (osize .gt. 0) then
                A%tmp_buf(1:osize) = A%oimg(this_image())[src]%data(1:osize)
                do j = 1, osize
                   call A%add(A%tmp_buf(j)%x, A%tmp_buf(j)%y, A%tmp_buf(j)%z)
                end do
             end if
          end do
       end select
    end if

    sync all

    ! Reset off-image buffer
    do i = 1, num_images()
       A%oimg(i)%top_ = 0
    end do

    A%assembled = .true.
  end subroutine mat_finalize

  subroutine mat_zero(A)
    class(mat_t), intent(inout) :: A

    if (allocated(A%O%val)) then
       A%O%val = 0d0
    end if

    if (allocated(A%D%val)) then
       A%D%val = 0d0
    end if

  end subroutine mat_zero

  subroutine mat_add_scalar(A, r, c, v)
    class(mat_t), intent(inout) :: A
    integer, intent(in) :: r
    integer, intent(in) :: c
    real(kind=dp) :: v
    type(tuple_i4r8_t) :: mat_tuple
    type(tuple_2i4r8_t) :: oimg_mat_tuple
    integer :: i, lr, owner

    if (r .ge. A%range(1) .and. r .lt. A%range(2)) then

       lr = (r - A%range(1)) + 1 ! Local row

       if (A%assembled) then
          if (c .lt. A%range(1) .or. c .ge. A%range(2)) then
             associate(rpt => A%O%rpt, col =>A%O%col, val => A%O%val)
               do i = rpt(lr), rpt(lr+1)
                  if (col(i) .eq. c) then
                     val(i) = val(i) + v
                     return
                  end if
               end do
             end associate
          else
             associate(rpt => A%D%rpt, col =>A%D%col, val => A%D%val)
               do i = rpt(lr), rpt(lr+1)
                  if (col(i) .eq. c) then
                     val(i) = val(i) + v
                     return
                  end if
               end do
             end associate
          end if
       else       
          select type(ep => A%rs(lr)%data)
          type is (tuple_i4r8_t)
             do i = 1, A%rs(lr)%top_
                if (ep(i)%x .eq. c) then
                   ep(i)%y = ep(i)%y + v
                   return
                end if
             end do
          class default
             call neko_error('Invalid type (mat_add_scalar)')
          end select
       
          mat_tuple%x = c
          mat_tuple%y = v
          call A%rs(lr)%push(mat_tuple)
       end if
    else
       oimg_mat_tuple%x = r
       oimg_mat_tuple%y = c
       oimg_mat_tuple%z = v
       owner = bsearch(A%glb_range, r, 1, num_images())
       call mat_oimg_push(A%oimg(owner), oimg_mat_tuple)
    end if
    
  end subroutine mat_add_scalar

  subroutine mat_add_block(A, m, r, n, c, b)
    class(mat_t), intent(inout) :: A
    integer, intent(in) :: m
    integer, intent(in) :: r(m)
    integer, intent(in) :: n
    integer, intent(in) :: c(n)
    real(kind=dp), intent(in) :: b(m*n)
    type(tuple_i4r8_t) :: mat_tuple
    type(tuple_2i4r8_t) :: oimg_mat_tuple
    integer :: i, j, k, l, lr, owner

    l = 0
    do i = 1, m
       do j = 1, n
          l = l + 1
          if (r(i) .ge. A%range(1) .and. r(i) .lt. A%range(2)) then

             lr= (r(i) - A%range(1)) + 1 ! local row
             if (A%assembled) then
                if (c(j) .lt. A%range(1) .or. c(j) .ge. A%range(2)) then
                   associate(rpt => A%O%rpt, col =>A%O%col, val => A%O%val)
                     do k = rpt(lr), rpt(lr+1)
                        if (col(k) .eq. c(j)) then
                           val(k) = val(k) + b(l)
                           return
                        end if
                     end do
                   end associate
                else
                   associate(rpt => A%D%rpt, col =>A%D%col, val => A%D%val)
                     do k = rpt(lr), rpt(lr+1)
                        if (col(k) .eq. c(j)) then
                           val(k) = val(k) + b(l)
                           return
                        end if
                     end do
                   end associate
                end if
             else
                select type(ep => A%rs(lr)%data)
                type is (tuple_i4r8_t)
                   do k = 1, A%rs(lr)%top_
                      if (ep(k)%x .eq. c(j)) then
                         ep(k)%y = ep(k)%y + b(l)
                         goto 42
                      end if
                   end do
                class default
                   call neko_error('Invalid type (mat_add_block)')
                end select
             
                mat_tuple%x = c(j)
                mat_tuple%y = b(l)
                call A%rs(lr)%push(mat_tuple)
             end if
          else
             oimg_mat_tuple%x = r(i)
             oimg_mat_tuple%y = c(j)
             oimg_mat_tuple%z = b(l)
             owner = bsearch(A%glb_range, r(i), 1, num_images())
             call mat_oimg_push(A%oimg(owner), oimg_mat_tuple)
42        end if

       end do
    end do
    
  end subroutine mat_add_block
  
  subroutine mat_oimg_push(oimg, oimg_mat_tuple)
    type(mat_oimg_t), intent(inout) :: oimg
    type(tuple_2i4r8_t), intent(inout) :: oimg_mat_tuple
    type(tuple_2i4r8_t), allocatable :: tmp(:)

    if (oimg%top_ .eq. oimg%size_) then
       oimg%size_ = ishft(oimg%size_, 1)
       allocate(tmp(oimg%size_))
       tmp(1:oimg%top_) = oimg%data
       call move_alloc(tmp, oimg%data)
    end if

    oimg%top_ = oimg%top_ + 1
    oimg%data(oimg%top_) = oimg_mat_tuple
    
  end subroutine mat_oimg_push
  
end module mat
