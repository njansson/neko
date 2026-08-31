! Copyright (c) 2020-2025, The Neko Authors
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
!> Defines GPU-aware OpenSHMEM gather-scatter communication
module gs_device_cray_shmem
  use num_types, only : rp, c_rp, i8
  use gs_comm, only : gs_comm_t, GS_VEC_NC
  use stack, only : stack_i4_t
  use htable, only : htable_i4_t
  use comm, only : NEKO_COMM, pe_rank, pe_size, global_pe_size
  use mpi_f08, only : MPI_Allreduce, MPI_Alltoall, MPI_INTEGER, MPI_MAX
  use device, only : device_alloc, device_free, device_memcpy, device_sync, &
       device_get_ptr, device_stream_wait_event, HOST_TO_DEVICE, &
       DEVICE_TO_HOST
  use utils, only : neko_error
#ifdef HAVE_OPENSHMEM
  use shmem, only : shmem_calloc, shmem_free, shmem_barrier_all, &
       shmem_putmem_signal_nbi, shmem_signal_wait_until, &
       shmemx_space_create, shmemx_space_malloc, &
       shmemx_space_free, shmemx_query_gpu_awareness, shmem_space_config_t, &
       SHMEM_MTYPE_GPU, SHMEM_TEAM_WORLD, SHMEM_SIGNAL_SET, SHMEM_CMP_GE
#endif
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_loc, &
       c_f_pointer, c_associated, c_sizeof, c_size_t, c_int, c_int32_t, &
       c_int64_t, c_long
  implicit none
  private

  !> Whether this build can use the backend at all: a native OpenSHMEM
  !! library (--with-openshmem) plus a CUDA or HIP device backend for the
  !! pack/unpack kernels. Whether the *library* is GPU-aware is a runtime
  !! property and is checked in init, see shmemx_query_gpu_awareness.
#if defined(HAVE_OPENSHMEM) && (defined(HAVE_HIP) || defined(HAVE_CUDA))
  logical, parameter, public :: GS_DEVICE_CRAY_SHMEM_AVAIL = .true.
#else
  logical, parameter, public :: GS_DEVICE_CRAY_SHMEM_AVAIL = .false.
#endif

  !> Round-robin depth of the two buffers. Cycling them is what lets the
  !! back-pressure handshake be implicit, costing no messages at all:
  !! with peer sets symmetric (checked in init) every peer both sends to
  !! and receives from us each round, and the data signal already carries
  !! everything an explicit ack would.
  !!
  !! Two recv slabs. Peer p writes our slab r mod 2 in round r, and last
  !! wrote it in round r-2, which we unpacked in nbwait(r-2). p's round r
  !! put follows p's nbwait(r-1), where p waited for our round r-1 data;
  !! we put that only after nbsend(r-1)'s device_sync, which completes
  !! the nbwait(r-2) unpacks. So p cannot overwrite a slab we are still
  !! reducing.
  !!
  !! Three send slabs. Repacking slab r mod 3 in round r needs our round
  !! r-3 put drained. We received p's round r-2 data in nbwait(r-2),
  !! which p sent after its own nbsend(r-2) sync, which completes p's
  !! nbwait(r-3) unpack of our round r-3 slab. Two slabs would rest on
  !! p's round r-1 data, which we have not yet waited for when we pack.
  integer, parameter :: GS_CRAY_SHMEM_SEND_SLABS = 3
  integer, parameter :: GS_CRAY_SHMEM_RECV_SLABS = 2

  !> Headroom factor applied to the first instance's need when sizing the
  !! process-wide space, so the gs instances created later (the scalar,
  !! the multigrid levels, all smaller than the fine fluid gs that comes
  !! first) fit as well without a second space.
  integer(c_size_t), parameter :: GS_CRAY_SHMEM_HEADROOM = 4

  !> Floor and default unit for the space size
  integer(c_size_t), parameter :: GS_CRAY_SHMEM_MIB = 1048576
  integer(c_size_t), parameter :: GS_CRAY_SHMEM_MIN = 32*GS_CRAY_SHMEM_MIB

  !> Whether to stage the packed halo through host memory before putting
  !! it, set from NEKO_GS_CRAY_SHMEM_HOST_SEND.
  !!
  !! The send buffer is only ever a local source, never a remote target,
  !! and the library accepts host memory there, so this trades the
  !! library's per-message copy out of GPU memory into the NIC command
  !! queue (the INJECT path, which only exists for GPU-attached sources)
  !! for a single contiguous device-to-host copy of the packed slab. The
  !! receive side stays on the GPU either way.
  !!
  !! Measured on Dardel (MI250X, Slingshot 11): no difference at all, so
  !! the sender's GPU involvement is not what makes this backend slower
  !! than device MPI there -- the cost tracks the destination being GPU
  !! memory, not the source. Kept because the HPE CUG 2025 paper measures
  !! that per-message copy as significant on PCIe-attached discrete GPUs,
  !! where this may yet pay for itself. Off by default.
  logical :: gs_cray_shmem_host_send = .false.
  logical :: gs_cray_shmem_host_send_read = .false.

  !> The process-wide GPU memory space every instance allocates from.
  !!
  !! HPE Slingshot SHMEM 12 permits a single space per job ("LIBSMA ERROR:
  !! cannot create more than 1 new space"), so this is a device symmetric
  !! heap shared by every gs instance rather than one space each, and it
  !! is sized up front for all of them (see gs_cray_shmem_space_init).
  !!
  !! It is deliberately never destroyed. The documented cap is on spaces,
  !! and whether destroying one returns the quota is unspecified, so a
  !! teardown followed by a new gs would risk an unrecoverable abort. The
  !! buffers inside it are freed normally, and the space itself goes away
  !! with the process.
  type(c_ptr) :: gs_cray_shmem_space = C_NULL_PTR

  !> Symmetric device buffer for one direction of communication
  type, private :: gs_device_cray_shmem_buf_t
     !> Number of dofs per neighbor
     integer, allocatable :: ndofs(:)
     !> Local offset in this buffer for each neighbor
     integer, allocatable :: offset(:)
     !> For send_buf: offset in the remote PE's recv buffer where our data
     !! should land. Unused for recv_buf.
     integer, allocatable :: remote_offset(:)
     !> Total number of dofs in this buffer on this PE
     integer :: total = 0
     !> Maximum total across all PEs; also the scalar slab stride, the
     !! vector slab stride being GS_VEC_NC times it
     integer :: max_total = 0
     !> Number of round-robin slabs in the buffer
     integer :: nslabs = 1
     !> Symmetric buffer in the GPU memory space. A device address: it is
     !! only ever passed to kernels and to SHMEM, never dereferenced here.
     type(c_ptr) :: buf_d = C_NULL_PTR
     !> Dof mapping for the pack/unpack kernels (ordinary device memory)
     type(c_ptr) :: dof_d = C_NULL_PTR
   contains
     procedure, pass(this) :: init => gs_device_cray_shmem_buf_init
     procedure, pass(this) :: alloc => gs_device_cray_shmem_buf_alloc
     procedure, pass(this) :: free => gs_device_cray_shmem_buf_free
  end type gs_device_cray_shmem_buf_t

  !> Gather-scatter communication using GPU-aware OpenSHMEM (HPE Slingshot
  !! SHMEM 12 and later).
  !!
  !! The protocol is gs_shmem's, with its explicit acks removed. Each PE
  !! holds symmetric send and recv buffers sized to the global maximum dof
  !! count, cycled over several slabs, and one symmetric uint64 array of
  !! pe_size slots: data_signals[r] is set to `iter` by PE r when it has
  !! put data into our recv buffer, via shmem_putmem_signal_nbi.
  !!
  !! There is no ack message. gs_shmem needs one so a fast sender cannot
  !! overwrite a buffer the receiver has not read, but with rotating slabs
  !! and symmetric peer sets the data signal already proves it: a peer's
  !! round r data can only have been sent after that peer's round r-1
  !! unpack completed. The slab-count parameters carry the full argument.
  !! Per round this leaves exactly one message per peer, as MPI has.
  !!
  !! What differs from the host backend is where the data lives and who
  !! moves it. Both buffers are allocated from a GPU-attached memory space
  !! (shmemx_space_create with SHMEM_MTYPE_GPU), so the puts move data
  !! device-to-device without staging through the host, while the packing
  !! and the reduction are done by the same kernels the device MPI backend
  !! uses. The signal words stay in the default (CPU) symmetric heap: every
  !! wait is a host-side poll, and polling device memory from the host would
  !! be far slower.
  !!
  !! Cray SHMEM is GPU-aware but not GPU-initiated -- there is no
  !! stream-ordered put and no call that can be issued from inside a kernel
  !! -- so the host has to synchronize the stream between the pack and the
  !! first put, and again between the unpacks and the acks. That is the
  !! structural cost of this backend against gs_device_shmem (NVSHMEM),
  !! which issues its puts from the push kernel.
  type, public, extends(gs_comm_t) :: gs_device_cray_shmem_t
     type(gs_device_cray_shmem_buf_t) :: send_buf
     type(gs_device_cray_shmem_buf_t) :: recv_buf
     !> Host mirror of the send buffer, same slab layout; allocated only
     !! when gs_cray_shmem_host_send is on
     real(kind=rp), allocatable :: host_send(:)
     !> Symmetric data-arrival signals, pe_size slots, in the CPU heap
     type(c_ptr) :: data_signals_ptr = C_NULL_PTR
     !> Monotonically increasing round counter; the sender writes it into
     !! the receiver's data signal slot with every put. All waits use
     !! CMP_GE, so peer counts and orderings may differ between ranks.
     !!
     !! The counter must stay in lockstep across ranks, since a sender
     !! picks the receiver's slab from its own value: every rank runs the
     !! same sequence of gs ops on a given instance (SPMD), which is the
     !! same assumption gs_device_shmem makes for its parity slabs.
     integer(kind=i8) :: iter = 0
   contains
     procedure, pass(this) :: init => gs_device_cray_shmem_init
     procedure, pass(this) :: free => gs_device_cray_shmem_free
     procedure, pass(this) :: nbsend => gs_device_cray_shmem_nbsend
     procedure, pass(this) :: nbrecv => gs_device_cray_shmem_nbrecv
     procedure, pass(this) :: nbwait => gs_device_cray_shmem_nbwait
     procedure, pass(this) :: nbsend_vec => gs_device_cray_shmem_nbsend_vec
     procedure, pass(this) :: nbrecv_vec => gs_device_cray_shmem_nbrecv_vec
     procedure, pass(this) :: nbwait_vec => gs_device_cray_shmem_nbwait_vec
  end type gs_device_cray_shmem_t

#ifdef HAVE_HIP
  interface
     subroutine hip_gs_pack(u_d, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'hip_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, offset
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_pack
  end interface

  interface
     subroutine hip_gs_unpack(u_d, op, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'hip_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_unpack
  end interface

  interface
     subroutine hip_gs_pack_vec(u_d, buf_d, dof_d, offset, n, nc, ns, stream) &
          bind(c, name = 'hip_gs_pack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_pack_vec
  end interface

  interface
     subroutine hip_gs_unpack_vec(u_d, op, buf_d, dof_d, offset, n, nc, ns, &
          stream) bind(c, name = 'hip_gs_unpack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_unpack_vec
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_gs_pack(u_d, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'cuda_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, offset
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_pack
  end interface

  interface
     subroutine cuda_gs_unpack(u_d, op, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'cuda_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_unpack
  end interface

  interface
     subroutine cuda_gs_pack_vec(u_d, buf_d, dof_d, offset, n, nc, ns, stream) &
          bind(c, name = 'cuda_gs_pack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_pack_vec
  end interface

  interface
     subroutine cuda_gs_unpack_vec(u_d, op, buf_d, dof_d, offset, n, nc, ns, &
          stream) bind(c, name = 'cuda_gs_unpack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_unpack_vec
  end interface
#endif

contains

  !> Create the process-wide GPU memory space, once. Sized to a multiple
  !! of what the first gs instance needs, which is the fine fluid gs and
  !! so the largest, with a floor; override with NEKO_GS_CRAY_SHMEM_SIZE,
  !! in MiB. Collective: every PE reaches this with the same size, since
  !! the buffers are sized by a global maximum.
  !! @param bytes_needed what the calling instance is about to allocate
  subroutine gs_cray_shmem_space_init(bytes_needed)
    integer(c_size_t), intent(in) :: bytes_needed
    character(len=255) :: env_size
    integer :: env_len, env_mib, ierr
    integer(c_size_t) :: sz
#ifdef HAVE_OPENSHMEM
    type(shmem_space_config_t) :: config

    if (.not. gs_cray_shmem_host_send_read) then
       call get_environment_variable("NEKO_GS_CRAY_SHMEM_HOST_SEND", &
            env_size, env_len)
       gs_cray_shmem_host_send = (env_len .gt. 0)
       if (env_len .gt. 0) then
          gs_cray_shmem_host_send = (env_size(1:env_len) .ne. "0")
       end if
       gs_cray_shmem_host_send_read = .true.
    end if

    if (c_associated(gs_cray_shmem_space)) return

    sz = max(GS_CRAY_SHMEM_HEADROOM * bytes_needed, GS_CRAY_SHMEM_MIN)

    call get_environment_variable("NEKO_GS_CRAY_SHMEM_SIZE", env_size, env_len)
    if (env_len .gt. 0) then
       read(env_size(1:env_len), *, iostat = ierr) env_mib
       if (ierr .ne. 0 .or. env_mib .le. 0) then
          call neko_error('NEKO_GS_CRAY_SHMEM_SIZE must be a size in MiB')
       end if
       sz = int(env_mib, c_size_t) * GS_CRAY_SHMEM_MIB
    end if

    config%size = sz
    config%mtype = SHMEM_MTYPE_GPU
    if (shmemx_space_create(SHMEM_TEAM_WORLD, 0_c_long, config, &
         gs_cray_shmem_space) .ne. 0) then
       call neko_error('shmemx_space_create failed for the GPU memory space')
    end if
#endif
  end subroutine gs_cray_shmem_space_init

  !> Per-neighbor bookkeeping and the dof map for the pack/unpack kernels.
  !! Does not touch symmetric memory; the buffer itself is allocated in
  !! alloc() once the size of the space is known for both directions.
  !! @param pe_order ranks in send_pe / recv_pe order (1-based)
  !! @param dof_stack indexed by rank, lower bound = 0 (per init_dofs)
  !! @param mark_dupes mark dofs appearing more than once, for the unpack
  subroutine gs_device_cray_shmem_buf_init(this, pe_order, dof_stack, &
       mark_dupes)
    class(gs_device_cray_shmem_buf_t), intent(inout) :: this
    integer, allocatable, intent(inout) :: pe_order(:)
    type(stack_i4_t), allocatable, intent(inout) :: dof_stack(:)
    logical, intent(in) :: mark_dupes
    integer, allocatable :: dofs(:)
    integer :: i, j, k, n, dupe, marked
    integer(c_size_t) :: sz
    type(htable_i4_t) :: doftable
    integer(c_int32_t) :: i4_dummy

    n = size(pe_order)

    allocate(this%ndofs(n))
    allocate(this%offset(n))
    allocate(this%remote_offset(n))

    do i = 1, n
       this%remote_offset(i) = -1
    end do

    this%total = 0
    do i = 1, n
       this%ndofs(i) = dof_stack(pe_order(i))%size()
       this%offset(i) = this%total
       this%total = this%total + this%ndofs(i)
    end do

    sz = c_sizeof(i4_dummy) * int(max(this%total, 1), c_size_t)
    call device_alloc(this%dof_d, sz)

    if (mark_dupes) call doftable%init(2*max(this%total, 1))
    allocate(dofs(max(this%total, 1)))

    ! Copy from dof_stack into dofs, optionally marking duplicates
    marked = 0
    do i = 1, n
       ! %array() breaks on cray
       select type (arr => dof_stack(pe_order(i))%data)
       type is (integer)
          do j = 1, this%ndofs(i)
             k = this%offset(i) + j
             if (mark_dupes) then
                if (doftable%get(arr(j), dupe) .eq. 0) then
                   if (dofs(dupe) .gt. 0) then
                      dofs(dupe) = -dofs(dupe)
                      marked = marked + 1
                   end if
                   dofs(k) = -arr(j)
                   marked = marked + 1
                else
                   call doftable%set(arr(j), k)
                   dofs(k) = arr(j)
                end if
             else
                dofs(k) = arr(j)
             end if
          end do
       end select
    end do

    if (this%total .gt. 0) then
       call device_memcpy(dofs, this%dof_d, this%total, HOST_TO_DEVICE, &
            sync = .true.)
    end if

    deallocate(dofs)
    call doftable%free()

  end subroutine gs_device_cray_shmem_buf_init

  !> Allocate the symmetric buffer from the GPU memory space @a space.
  !! Sized for up to GS_VEC_NC components so the fused vector path can
  !! reuse it; the scalar path uses the first max_total of each slab.
  !! @param max_total global maximum dof count for this direction
  !! @param nslabs number of round-robin slabs
  subroutine gs_device_cray_shmem_buf_alloc(this, space, max_total, nslabs)
    class(gs_device_cray_shmem_buf_t), intent(inout) :: this
    type(c_ptr), intent(in) :: space
    integer, intent(in) :: max_total, nslabs
    integer(c_size_t) :: sz
    real(c_rp) :: rp_dummy
#ifndef HAVE_OPENSHMEM
    call neko_error('Neko was not built with OpenSHMEM support')
#else

    this%max_total = max_total
    this%nslabs = nslabs

    sz = c_sizeof(rp_dummy) * &
         int(max(nslabs*GS_VEC_NC*max_total, 1), c_size_t)
    this%buf_d = shmemx_space_malloc(space, sz)
    if (.not. c_associated(this%buf_d)) then
       call neko_error('shmemx_space_malloc failed, raise ' // &
            'NEKO_GS_CRAY_SHMEM_SIZE (MiB)')
    end if
#endif
  end subroutine gs_device_cray_shmem_buf_alloc

  !> Release the symmetric buffer and the bookkeeping.
  !! @param space the space @a buf_d was allocated from
  subroutine gs_device_cray_shmem_buf_free(this, space)
    class(gs_device_cray_shmem_buf_t), intent(inout) :: this
    type(c_ptr), intent(in) :: space

    if (allocated(this%ndofs)) deallocate(this%ndofs)
    if (allocated(this%offset)) deallocate(this%offset)
    if (allocated(this%remote_offset)) deallocate(this%remote_offset)

#ifdef HAVE_OPENSHMEM
    if (c_associated(this%buf_d) .and. c_associated(space)) then
       call shmemx_space_free(space, this%buf_d)
    end if
#endif
    this%buf_d = C_NULL_PTR

    if (c_associated(this%dof_d)) call device_free(this%dof_d)
    this%dof_d = C_NULL_PTR

    this%total = 0
    this%max_total = 0

  end subroutine gs_device_cray_shmem_buf_free

  !> Initialise GPU-aware OpenSHMEM based communication method
  subroutine gs_device_cray_shmem_init(this, send_pe, recv_pe)
    class(gs_device_cray_shmem_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer :: i, ierr
    integer :: totals(2), max_totals(2)
    integer(c_size_t) :: i64_size, bytes_send, bytes_recv
    integer(c_int64_t) :: i64_dummy
    integer, allocatable :: local_offsets(:), remote_offsets(:)
    real(c_rp) :: rp_dummy
#if !defined(HAVE_OPENSHMEM) || (!defined(HAVE_HIP) && !defined(HAVE_CUDA))
    call neko_error('gs_device_cray_shmem: no OpenSHMEM or no device backend')
#else

    ! SHMEM PEs are the whole job, so a gs on a subset of the ranks cannot
    ! address its peers by PE number (same restriction as gs_shmem).
    if (pe_size .ne. global_pe_size) then
       call neko_error('gs_device_cray_shmem requires all ranks in the job')
    end if

    ! GPU-awareness is a runtime property of the library: without it the
    ! puts below would be handed device pointers the library will not
    ! recognise. Requires SHMEM_GPU_SUPPORT_ENABLED=1 in the environment.
    if (shmemx_query_gpu_awareness(SHMEM_TEAM_WORLD) .ne. 1) then
       call neko_error('gs_device_cray_shmem: library is not GPU-aware, ' // &
            'set SHMEM_GPU_SUPPORT_ENABLED=1')
    end if

    call this%init_order(send_pe, recv_pe)

    ! The implicit handshake (see the slab-count parameters) needs every
    ! peer to appear in both directions, so that receiving a peer's data
    ! is evidence about its progress. The gs schedule registers every
    ! sharing peer both ways; fail loudly if that ever changes. The order
    ! may differ -- nothing here indexes the two lists together.
    if (size(this%send_pe) .ne. size(this%recv_pe)) then
       call neko_error('gs_device_cray_shmem requires symmetric peer sets')
    end if
    do i = 1, size(this%send_pe)
       if (.not. any(this%recv_pe .eq. this%send_pe(i))) then
          call neko_error('gs_device_cray_shmem requires symmetric peer sets')
       end if
    end do

    call this%send_buf%init(this%send_pe, this%send_dof, .false.)
    call this%recv_buf%init(this%recv_pe, this%recv_dof, .true.)

    ! Symmetric memory must have the same size on every PE, and the space
    ! has to be created before anything can be allocated from it, so both
    ! directions are reduced in one go here rather than inside buf%init.
    totals(1) = this%send_buf%total
    totals(2) = this%recv_buf%total
    call MPI_Allreduce(totals, max_totals, 2, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)

    ! Both direction buffers come out of the one process-wide space. Only
    ! the recv buffer is ever a remote target -- the send buffer is purely
    ! a local source, which the library also accepts as ordinary device
    ! memory -- but keeping both in the space means both are registered
    ! with the NIC up front.
    bytes_send = c_sizeof(rp_dummy) * &
         int(max(GS_CRAY_SHMEM_SEND_SLABS*GS_VEC_NC*max_totals(1), 1), &
         c_size_t)
    bytes_recv = c_sizeof(rp_dummy) * &
         int(max(GS_CRAY_SHMEM_RECV_SLABS*GS_VEC_NC*max_totals(2), 1), &
         c_size_t)
    call gs_cray_shmem_space_init(bytes_send + bytes_recv)

    call this%send_buf%alloc(gs_cray_shmem_space, max_totals(1), &
         GS_CRAY_SHMEM_SEND_SLABS)
    call this%recv_buf%alloc(gs_cray_shmem_space, max_totals(2), &
         GS_CRAY_SHMEM_RECV_SLABS)

    if (gs_cray_shmem_host_send) then
       allocate(this%host_send(max(GS_CRAY_SHMEM_SEND_SLABS*GS_VEC_NC * &
            max_totals(1), 1)))
    end if

    ! Per-rank symmetric signal arrays, in the default (CPU) heap: every
    ! wait on them is a host-side poll. Size pe_size on every PE, so no
    ! slot handshake is needed -- each PE writes at the remote PE's rank.
    i64_size = c_sizeof(i64_dummy)
    this%data_signals_ptr = shmem_calloc(int(pe_size, c_size_t), i64_size)
    if (.not. c_associated(this%data_signals_ptr)) then
       call neko_error('shmem_calloc failed for gs data signals')
    end if

    ! Exchange, for every send peer, the offset in the receiver's recv
    ! buffer where our slab must land. A single Alltoall keeps this
    ! deadlock-free for arbitrary, non-uniform peer sets.
    allocate(local_offsets(0:pe_size - 1))
    allocate(remote_offsets(0:pe_size - 1))
    local_offsets = -1
    do i = 1, size(this%recv_pe)
       local_offsets(this%recv_pe(i)) = this%recv_buf%offset(i)
    end do
    call MPI_Alltoall(local_offsets, 1, MPI_INTEGER, &
         remote_offsets, 1, MPI_INTEGER, NEKO_COMM, ierr)
    do i = 1, size(this%send_pe)
       this%send_buf%remote_offset(i) = remote_offsets(this%send_pe(i))
    end do
    deallocate(local_offsets)
    deallocate(remote_offsets)

    this%iter = 0
    this%vec_supported = .true.

    ! No one-sided communication before every PE has its buffers.
    call shmem_barrier_all()
#endif
  end subroutine gs_device_cray_shmem_init

  !> Deallocate GPU-aware OpenSHMEM based communication method
  subroutine gs_device_cray_shmem_free(this)
    class(gs_device_cray_shmem_t), intent(inout) :: this

#ifdef HAVE_OPENSHMEM
    ! The frees below are collective; synchronize first so no in-flight put
    ! targets a buffer we are about to release.
    call shmem_barrier_all()

    if (c_associated(this%data_signals_ptr)) then
       call shmem_free(this%data_signals_ptr)
    end if
    this%data_signals_ptr = C_NULL_PTR
    this%iter = 0
#endif

    ! The buffers are released, the space they came from is not; see the
    ! gs_cray_shmem_space declaration.
    call this%send_buf%free(gs_cray_shmem_space)
    call this%recv_buf%free(gs_cray_shmem_space)

    if (allocated(this%host_send)) deallocate(this%host_send)

    call this%free_order()
    call this%free_dofs()

  end subroutine gs_device_cray_shmem_free

  !> Pack the shared dofs into this round's send slab and put them into
  !! each neighbor's recv slab with a signal.
  !!
  !!
  !! The device_sync here is load-bearing twice over: it completes the
  !! pack before the first host-issued put, and it completes the previous
  !! round's unpacks (the caller's stream was joined with every unpack
  !! stream in nbwait), which is what makes the put that follows evidence
  !! to our peers that we are done with their last slab.
  subroutine gs_device_cray_shmem_nbsend(this, u, n, tag, deps, strm)
    class(gs_device_cray_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, dst, sbase, rbase
    integer(c_size_t) :: nbytes
    type(c_ptr) :: u_d, src_p, dst_p
    real(kind=rp), pointer :: send_data(:), recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:)
    real(c_rp) :: rp_dummy
#ifdef HAVE_OPENSHMEM

    u_d = device_get_ptr(u)

    ! Each gs op gets a fresh signal value so receivers can tell one
    ! round's put from the next.
    this%iter = this%iter + 1

    call gs_device_cray_shmem_bufs(this, send_data, recv_data, data_signals)
    sbase = gs_cray_shmem_slab(this%send_buf, this%iter)
    rbase = gs_cray_shmem_slab(this%recv_buf, this%iter)

    ! Bulk-pack every peer's slab in one kernel, ordered after the gather.
    ! Doing it here, before any per-peer work, keeps every read of u ahead
    ! of the unpack writes to u in nbwait. The slab is free: it was last
    ! put in round iter-3, and every peer's round iter-2 data, which we
    ! waited for in nbwait, implies that put drained.
    call device_stream_wait_event(strm, deps, 0)
    call gs_device_cray_shmem_pack(u_d, c_loc(send_data(sbase + 1)), &
         this%send_buf%dof_d, 0, this%send_buf%total, strm)

    ! One contiguous copy of the packed slab, in place of the library's
    ! per-message copy out of GPU memory; see gs_cray_shmem_host_send.
    if (gs_cray_shmem_host_send .and. this%send_buf%total .gt. 0) then
       src_p = c_loc(send_data(sbase + 1))
       dst_p = c_loc(this%host_send(sbase + 1))
       call device_memcpy(dst_p, src_p, &
            int(this%send_buf%total, c_size_t) * c_sizeof(rp_dummy), &
            DEVICE_TO_HOST, sync = .false., strm = strm)
    end if

    ! Cray SHMEM has no stream-ordered put, so the pack must be finished
    ! before the first one; see the routine comment for what else this
    ! sync carries.
    call device_sync(strm)

    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)

       if (gs_cray_shmem_host_send) then
          src_p = c_loc(this%host_send(sbase + this%send_buf%offset(i) + 1))
       else
          src_p = c_loc(send_data(sbase + this%send_buf%offset(i) + 1))
       end if

       nbytes = int(this%send_buf%ndofs(i), c_size_t) * c_sizeof(rp_dummy)
       call shmem_putmem_signal_nbi( &
            c_loc(recv_data(rbase + this%send_buf%remote_offset(i) + 1)), &
            src_p, &
            nbytes, &
            c_loc(data_signals(pe_rank + 1)), &
            this%iter, SHMEM_SIGNAL_SET, dst)
    end do
#endif
  end subroutine gs_device_cray_shmem_nbsend

  !> No-op: receives are completed by the remote put-with-signal.
  subroutine gs_device_cray_shmem_nbrecv(this, tag)
    class(gs_device_cray_shmem_t), intent(inout) :: this
    integer, intent(in) :: tag

  end subroutine gs_device_cray_shmem_nbrecv

  !> Wait per-neighbor for the signal that its data has landed and reduce
  !! the slab into u on that neighbor's stream. The unpacks are only
  !! joined into the caller's stream, never waited for on the host: the
  !! next nbsend's device_sync is what completes them.
  subroutine gs_device_cray_shmem_nbwait(this, u, n, op, strm)
    class(gs_device_cray_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op, i, src, rbase
    integer(c_int64_t) :: dummy
    type(c_ptr) :: u_d
    real(kind=rp), pointer :: send_data(:), recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:)
#ifdef HAVE_OPENSHMEM

    u_d = device_get_ptr(u)

    call gs_device_cray_shmem_bufs(this, send_data, recv_data, data_signals)
    rbase = gs_cray_shmem_slab(this%recv_buf, this%iter)

    ! Wait for every peer's slab, then reduce the whole recv buffer in one
    ! kernel. Per-peer unpacks on per-peer streams, launched as each slab
    ! lands, overlap the reduction with the waiting but cost a launch, an
    ! event record and a join per peer; at Neko's neighbour counts that
    ! overhead dominates the overlap, which is why the device MPI backend
    ! reduces in bulk too. A dof shared with several peers is handled by
    ! the duplicate marking in the recv dof map.
    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       dummy = shmem_signal_wait_until(c_loc(data_signals(src + 1)), &
            SHMEM_CMP_GE, this%iter)
    end do

    call gs_device_cray_shmem_unpack(u_d, op, c_loc(recv_data(rbase + 1)), &
         this%recv_buf%dof_d, 0, this%recv_buf%total, strm)
#endif
  end subroutine gs_device_cray_shmem_nbwait

  !> Fused nc-component send. @a u is the compact shared device buffer
  !! (component-outer, per-component stride n). The packed layout is
  !! component-minor per dof, so a peer's slab stays contiguous and the
  !! offsets within a slab are simply scaled by nc. See the scalar
  !! nbsend for the slab and ack reasoning.
  subroutine gs_device_cray_shmem_nbsend_vec(this, u, n, nc, tag, deps, strm)
    class(gs_device_cray_shmem_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, dst, sbase, rbase
    integer(c_size_t) :: nbytes
    type(c_ptr) :: u_d, src_p, dst_p
    real(kind=rp), pointer :: send_data(:), recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:)
    real(c_rp) :: rp_dummy
#ifdef HAVE_OPENSHMEM

    u_d = device_get_ptr(u)

    this%iter = this%iter + 1

    call gs_device_cray_shmem_bufs(this, send_data, recv_data, data_signals)
    sbase = gs_cray_shmem_slab(this%send_buf, this%iter)
    rbase = gs_cray_shmem_slab(this%recv_buf, this%iter)

    call device_stream_wait_event(strm, deps, 0)
    call gs_device_cray_shmem_pack_vec(u_d, c_loc(send_data(sbase + 1)), &
         this%send_buf%dof_d, 0, this%send_buf%total, nc, n, strm)

    if (gs_cray_shmem_host_send .and. this%send_buf%total .gt. 0) then
       src_p = c_loc(send_data(sbase + 1))
       dst_p = c_loc(this%host_send(sbase + 1))
       call device_memcpy(dst_p, src_p, &
            int(nc * this%send_buf%total, c_size_t) * c_sizeof(rp_dummy), &
            DEVICE_TO_HOST, sync = .false., strm = strm)
    end if

    call device_sync(strm)

    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)

       if (gs_cray_shmem_host_send) then
          src_p = c_loc(this%host_send(sbase + &
               nc*this%send_buf%offset(i) + 1))
       else
          src_p = c_loc(send_data(sbase + nc*this%send_buf%offset(i) + 1))
       end if

       nbytes = int(nc * this%send_buf%ndofs(i), c_size_t) * &
            c_sizeof(rp_dummy)
       call shmem_putmem_signal_nbi( &
            c_loc(recv_data(rbase + nc*this%send_buf%remote_offset(i) + 1)), &
            src_p, &
            nbytes, &
            c_loc(data_signals(pe_rank + 1)), &
            this%iter, SHMEM_SIGNAL_SET, dst)
    end do
#endif
  end subroutine gs_device_cray_shmem_nbsend_vec

  !> No-op: receives are completed by the remote put-with-signal.
  subroutine gs_device_cray_shmem_nbrecv_vec(this, tag, nc)
    class(gs_device_cray_shmem_t), intent(inout) :: this
    integer, intent(in) :: tag, nc

  end subroutine gs_device_cray_shmem_nbrecv_vec

  !> Fused nc-component wait and reduction, see the scalar nbwait.
  subroutine gs_device_cray_shmem_nbwait_vec(this, u, n, nc, op, strm)
    class(gs_device_cray_shmem_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op, i, src, rbase
    integer(c_int64_t) :: dummy
    type(c_ptr) :: u_d
    real(kind=rp), pointer :: send_data(:), recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:)
#ifdef HAVE_OPENSHMEM

    u_d = device_get_ptr(u)

    call gs_device_cray_shmem_bufs(this, send_data, recv_data, data_signals)
    rbase = gs_cray_shmem_slab(this%recv_buf, this%iter)

    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       dummy = shmem_signal_wait_until(c_loc(data_signals(src + 1)), &
            SHMEM_CMP_GE, this%iter)
    end do

    call gs_device_cray_shmem_unpack_vec(u_d, op, &
         c_loc(recv_data(rbase + 1)), this%recv_buf%dof_d, 0, &
         this%recv_buf%total, nc, n, strm)
#endif
  end subroutine gs_device_cray_shmem_nbwait_vec

  !> First element of the slab round @a iter uses in @a buf, in elements
  !! of the vector-sized layout. Rounds cycle through the slabs, so a
  !! slab is only rewritten every nslabs rounds.
  pure function gs_cray_shmem_slab(buf, iter) result(base)
    type(gs_device_cray_shmem_buf_t), intent(in) :: buf
    integer(kind=i8), intent(in) :: iter
    integer :: base

    base = int(mod(iter, int(buf%nslabs, i8))) * GS_VEC_NC * buf%max_total

  end function gs_cray_shmem_slab

#ifdef HAVE_OPENSHMEM
  !> Bind Fortran views to the symmetric buffers and signal arrays. The
  !! buffer views are device addresses: they exist so the puts below can
  !! name an element offset with c_loc, and are never dereferenced here.
  subroutine gs_device_cray_shmem_bufs(this, send_data, recv_data, &
       data_signals)
    class(gs_device_cray_shmem_t), intent(inout) :: this
    real(kind=rp), pointer, intent(out) :: send_data(:), recv_data(:)
    integer(c_int64_t), pointer, intent(out) :: data_signals(:)

    call c_f_pointer(this%send_buf%buf_d, send_data, &
         [max(this%send_buf%nslabs*GS_VEC_NC*this%send_buf%max_total, 1)])
    call c_f_pointer(this%recv_buf%buf_d, recv_data, &
         [max(this%recv_buf%nslabs*GS_VEC_NC*this%recv_buf%max_total, 1)])
    call c_f_pointer(this%data_signals_ptr, data_signals, [pe_size])

  end subroutine gs_device_cray_shmem_bufs
#endif

  !> Pack @a n dofs from @a u_d into @a buf_d at @a offset
  subroutine gs_device_cray_shmem_pack(u_d, buf_d, dof_d, offset, n, strm)
    type(c_ptr), intent(in) :: u_d, buf_d, dof_d, strm
    integer, intent(in) :: offset, n

#ifdef HAVE_HIP
    call hip_gs_pack(u_d, buf_d, dof_d, offset, n, strm)
#elif HAVE_CUDA
    call cuda_gs_pack(u_d, buf_d, dof_d, offset, n, strm)
#else
    call neko_error('gs_device_cray_shmem: no device backend')
#endif

  end subroutine gs_device_cray_shmem_pack

  !> Reduce @a n dofs from @a buf_d at @a offset into @a u_d with @a op
  subroutine gs_device_cray_shmem_unpack(u_d, op, buf_d, dof_d, offset, n, &
       strm)
    type(c_ptr), intent(in) :: u_d, buf_d, dof_d, strm
    integer, intent(in) :: op, offset, n

#ifdef HAVE_HIP
    call hip_gs_unpack(u_d, op, buf_d, dof_d, offset, n, strm)
#elif HAVE_CUDA
    call cuda_gs_unpack(u_d, op, buf_d, dof_d, offset, n, strm)
#else
    call neko_error('gs_device_cray_shmem: no device backend')
#endif

  end subroutine gs_device_cray_shmem_unpack

  !> Fused nc-component pack, see gs_device_cray_shmem_pack
  subroutine gs_device_cray_shmem_pack_vec(u_d, buf_d, dof_d, offset, n, nc, &
       ns, strm)
    type(c_ptr), intent(in) :: u_d, buf_d, dof_d, strm
    integer, intent(in) :: offset, n, nc, ns

#ifdef HAVE_HIP
    call hip_gs_pack_vec(u_d, buf_d, dof_d, offset, n, nc, ns, strm)
#elif HAVE_CUDA
    call cuda_gs_pack_vec(u_d, buf_d, dof_d, offset, n, nc, ns, strm)
#else
    call neko_error('gs_device_cray_shmem: no device backend')
#endif

  end subroutine gs_device_cray_shmem_pack_vec

  !> Fused nc-component unpack, see gs_device_cray_shmem_unpack
  subroutine gs_device_cray_shmem_unpack_vec(u_d, op, buf_d, dof_d, offset, &
       n, nc, ns, strm)
    type(c_ptr), intent(in) :: u_d, buf_d, dof_d, strm
    integer, intent(in) :: op, offset, n, nc, ns

#ifdef HAVE_HIP
    call hip_gs_unpack_vec(u_d, op, buf_d, dof_d, offset, n, nc, ns, strm)
#elif HAVE_CUDA
    call cuda_gs_unpack_vec(u_d, op, buf_d, dof_d, offset, n, nc, ns, strm)
#else
    call neko_error('gs_device_cray_shmem: no device backend')
#endif

  end subroutine gs_device_cray_shmem_unpack_vec

end module gs_device_cray_shmem
