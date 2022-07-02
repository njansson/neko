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
!> Defines sparse matrix vector multiplication
module spmv
  use num_types
  use mat
  use vec
  implicit none

  
contains

  !> SPMV \f$ y = Ax \f$
  subroutine caf_spmv(A, x, y)
    type(mat_t), intent(in) :: A
    type(vec_t), intent(in) :: x
    type(vec_t), intent(inout) :: y
    real(kind=dp) :: tmp
    integer :: r, i

    ! Diagonal CSR block
    if (A%D%nz .gt. 0) then
       associate(rpt => A%D%rpt, col => A%D%col, val => A%D%val, &
            xp => x%X%x, yp => y%X%x)
         do r = 1, A%m
            tmp = 0d0
            do i = rpt(r), rpt(r+1) - 1
               tmp  = tmp + val(i) * xp(col(i))
            end do
            yp(r) = tmp
         end do
       end associate
    end if

    ! Off-Diagonal CSR block    
    if (A%O%nz .gt. 0) then
       associate(rpt => A%O%rpt, col => A%O%col, val => A%O%val, &
            xp => x%X%x, yp => y%X%x)
         do r = 1, A%m
            tmp = 0d0
            do i = rpt(r), rpt(r+1) - 1
               tmp  = tmp + val(i) * xp(col(i))
            end do
            yp(r) = yp(r) + tmp
         end do
       end associate
    end if
    
  end subroutine caf_spmv
  
end module spmv
