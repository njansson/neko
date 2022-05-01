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
     integer, allocatable :: col
     integer, allocatable :: rpt
     real(kind=dp), allocatable :: val
  end type mat_csr_t

  type, public :: mat_t
     type(mat_csr_t), allocatable :: D[:] !< Diagonal part 
     type(mat_csr_t), allocatable :: O[:] !< Off-diagonal part
     integer :: m                         !< Local number of rows
     integer :: n                         !< Local number of columns
     integer :: m_glb                     !< Global number of rows     
     integer :: n_glb                     !< Global number of columns
     integer :: range(2)                  !< Ownership range
     integer, allocatable :: glb_range(:) !< Start offset for each image
   contains
     procedure, pass(this) :: init => mat_init
     procedure, pass(this) :: free => mat_free
  end type mat_t

contains

  !> Initialise a matrix with @a m local rows and @a n local columns
  subroutine mat_init(this, m, n)
    class(mat_t), intent(inout) :: this
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, allocatable :: tmp[:]
    integer :: acc, i
    
    call mat_free(this)

    this%m = m
    this%m_glb = m
    call co_sum(this%m_glb)


    this%n = n
    this%n_glb = n
    call co_sum(this%n_glb)

    allocate(tmp[*])
    allocate(this%D[*])
    allocate(this%O[*])
    sync all

    tmp[this_image()] = this%m
    sync all

    acc = 0
    do i = 1, this_image()
       acc = acc + tmp[i]
    end do
    sync all

    this%range(1) = acc - this%m + 1
    this%range(2) = this%range(1) + this%m

    tmp[this_image()] = this%range(1)
    sync memory

    allocate(this%glb_range(num_images() + 1))

    do i = 1, num_images()
       this%glb_range(i) = tmp[i]
    end do
    this%glb_range(num_images() + 1) = this%m_glb


    sync all
    deallocate(tmp)
    
  end subroutine mat_init

  subroutine mat_free(this)
    class(mat_t), intent(inout) :: this

    this%m = 0
    this%n = 0
    
    this%m_glb = 0
    this%n_glb = 0

    if (allocated(this%D)) then       
       call mat_csr_free(this%D)
    end if

    if (allocated(this%O)) then
       call mat_csr_free(this%O)
    end if

    if (allocated(this%glb_range)) then
       deallocate(this%glb_range)
    end if
    
  end subroutine mat_free

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
  
end module mat
