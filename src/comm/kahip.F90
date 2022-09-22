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
!> Interface to KaHIP
module kahip
  use comm
  use utils
  use mesh, only : mesh_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  interface
     subroutine parhip_partition_kway &
          (vtxdist, xadj, adjncy, vwgt, adjwgt, nparts, imbalance, &
          suppress_output, seed, mode, edgecut, part) &
          bind(c, name='ParHIPPartitionKWay_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       ! idxtype variables
       type(c_ptr), value :: vtxdist, xadj, adjncy, vwgt, adjwgt, nparts,&
            edgecut, part
       ! double variables
       type(c_ptr), value :: imbalance
       ! int (and bool) variables
       integer(c_int) :: suppress_output, seed, mode
     end subroutine parhip_partition_kway
  end interface

contains

#ifdef HAVE_KAHIP
#else
#endif
  
end module kahip
