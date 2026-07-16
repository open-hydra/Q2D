!>@brief MPI communication of ghost cells across inter-block boundaries (2D).
!> Builds a communication schedule from the BC connectivity array, then
!> provides non-blocking exchange routines for the P and R arrays.
!> Messages are aggregated per remote rank to minimize latency.
!> In serial mode (USE_MPI not defined) all routines are no-ops.
!>
!> 2D data structures:
!> - `ghost_entry_type` with (i, j) coordinates only.
!> - `face_group_type` with guide offsets (di, dj) only.
!> - `P/R` of rank 3 (nprim, i, j).
!> - `Pg(nprim, 4)`: 2 ghost columns + 2 stencil columns.
!> - `PG_NCOLS = 2` (columns 3..4 of Pg).
!> - P and R exchange based on persistent MPI requests.
module Q2D_Mod_GhostExchange
  use iso_fortran_env, only: R8 => real64, I4 => int32
  use Q2D_Mod_MPI

  implicit none
  private

  !> Single ghost-cell exchange entry (2D coordinates i, j).
  type :: ghost_entry_type
    integer :: bc_idx          !< Index into domain%bc(:)
    integer :: local_block     !< Block on this rank
    integer :: remote_block    !< Block on the remote rank
    integer :: remote_rank     !< MPI rank of the remote block
    integer :: src_i, src_j, src_f  !< Source cell coordinates and face
    integer :: dst_i, dst_j, dst_f  !< Destination cell coordinates and face
  end type ghost_entry_type

  !> Per-rank aggregation data
  type :: rank_group_type
    integer :: rank = -1       !< Remote MPI rank
    integer :: offset = 0      !< Starting index in the entry list (1-based)
    integer :: count  = 0      !< Number of entries for this rank
  end type rank_group_type

  !> Face group: contiguous sequence of entries sharing (block, face).
  !> Precomputes the guide offsets (di, dj).
  type :: face_group_type
    integer :: block            !< Local block index
    integer :: face             !< Face number (1-4)
    integer :: entry_start      !< Index of the first entry in the send/recv list
    integer :: entry_count      !< Number of entries
    integer :: di, dj           !< Precomputed guide(face,1..2)
  end type face_group_type

  !> Communication schedule: entries grouped by remote rank for aggregated messages.
  type :: ghost_schedule_type
    integer :: n_send = 0, n_recv = 0
    type(ghost_entry_type), allocatable :: send_list(:)
    type(ghost_entry_type), allocatable :: recv_list(:)
    logical :: built = .false.

    ! Per-rank aggregation of sends
    integer :: n_send_ranks = 0
    type(rank_group_type), allocatable :: send_groups(:)

    ! Per-rank aggregation of receives
    integer :: n_recv_ranks = 0
    type(rank_group_type), allocatable :: recv_groups(:)

    ! Face groups for aggregated pack/unpack (reduce per-cell overhead)
    integer :: n_send_faces = 0
    type(face_group_type), allocatable :: send_faces(:)
    integer :: n_recv_faces = 0
    type(face_group_type), allocatable :: recv_faces(:)

    ! Preallocated MPI buffers (one contiguous buffer per rank group)
    real(R8), allocatable :: P_send_buf(:)   !< (nprim*gc * n_send)
    real(R8), allocatable :: P_recv_buf(:)   !< (nprim*gc * n_recv)
    real(R8), allocatable :: Pg_send_buf(:)  !< (nprim*PG_NCOLS * n_send) — 2 columns
    real(R8), allocatable :: Pg_recv_buf(:)  !< (nprim*PG_NCOLS * n_recv)
    integer, allocatable  :: send_req(:)     !< (n_send_ranks)
    integer, allocatable  :: recv_req(:)     !< (n_recv_ranks)
    integer, allocatable  :: mpi_stat(:,:)   !< (MPI_STATUS_SIZE, max(n_send_ranks, n_recv_ranks))
    ! Persistent MPI requests for the exchanges (initialized once, started at every step)
    integer, allocatable  :: P_send_req_pers(:)  !< (n_send_ranks)
    integer, allocatable  :: P_recv_req_pers(:)  !< (n_recv_ranks)
    integer, allocatable  :: R_send_req_pers(:)  !< (n_send_ranks)
    integer, allocatable  :: R_recv_req_pers(:)  !< (n_recv_ranks)
    logical :: persistent_P_init = .false.
    logical :: persistent_R_init = .false.
  end type ghost_schedule_type

  type(ghost_schedule_type), public :: ghost_sched

  public :: build_ghost_schedule
  public :: build_local_bc_index
  public :: cleanup_ghost_schedule
  public :: exchange_ghost_P_post_recv, exchange_ghost_P_pack
  public :: exchange_ghost_P_post_send, exchange_ghost_P_wait_unpack
  public :: exchange_ghost_P_wait_send
  public :: exchange_ghost_R
  public :: exchange_ghost_Pg
  public :: gather_P_to_root
  public :: gather_diagnostic_to_root
  public :: scatter_P_from_root
  public :: mpi_io_barrier
  public :: Ghost_Interrank

contains


  !> Builds the communication schedule by scanning BC entries of type 101
  !> (block-to-block connection) and 201 (periodic). Entries are sorted by
  !> remote rank for aggregated MPI sends. BC 103 (multi-solver coupling)
  !> does not take part in the inter-block ghost exchange: it is handled
  !> separately as a multi-solver coupling stub in the BC dispatch.
  !> To be called after partition_blocks and the domain construction.
  subroutine build_ghost_schedule(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: gc, nprim

    implicit none
    type(Q2D_domain_type), intent(in) :: domain
    ! Local variables
    integer :: i, ns, nr, bm, bs

    if (mpi_size_ <= 1) then
      ghost_sched%built = .true.
      return
    end if

#ifdef USE_MPI
    ! Count send and receive entries
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 101 .and. domain%bc(i)%type /= 201) cycle  ! inter-rank BC 201 (periodic) is a plain copy like 101
      bm = domain%bc(i)%b   ! destination block
      bs = domain%bc(i)%bs  ! source block
      if (.not. allocated(block_owner)) cycle
      ! If the rank owns the destination but not the source => it must RECEIVE
      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
      end if
      ! If the rank owns the source but not the destination => it must SEND
      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
      end if
    end do

    ghost_sched%n_send = ns
    ghost_sched%n_recv = nr
    if (allocated(ghost_sched%send_list)) deallocate(ghost_sched%send_list)
    if (allocated(ghost_sched%recv_list)) deallocate(ghost_sched%recv_list)
    allocate(ghost_sched%send_list(ns))
    allocate(ghost_sched%recv_list(nr))

    ! Fill the send and receive lists
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 101 .and. domain%bc(i)%type /= 201) cycle  ! inter-rank BC 201 (periodic) is a plain copy like 101
      bm = domain%bc(i)%b
      bs = domain%bc(i)%bs

      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
        ghost_sched%recv_list(nr)%bc_idx = i
        ghost_sched%recv_list(nr)%local_block  = bm
        ghost_sched%recv_list(nr)%remote_block = bs
        ghost_sched%recv_list(nr)%remote_rank  = block_owner(bs)
        ghost_sched%recv_list(nr)%src_i = domain%bc(i)%is
        ghost_sched%recv_list(nr)%src_j = domain%bc(i)%js
        ghost_sched%recv_list(nr)%src_f = domain%bc(i)%fs
        ghost_sched%recv_list(nr)%dst_i = domain%bc(i)%i
        ghost_sched%recv_list(nr)%dst_j = domain%bc(i)%j
        ghost_sched%recv_list(nr)%dst_f = domain%bc(i)%f
      end if

      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
        ghost_sched%send_list(ns)%bc_idx = i
        ghost_sched%send_list(ns)%local_block  = bs
        ghost_sched%send_list(ns)%remote_block = bm
        ghost_sched%send_list(ns)%remote_rank  = block_owner(bm)
        ghost_sched%send_list(ns)%src_i = domain%bc(i)%is
        ghost_sched%send_list(ns)%src_j = domain%bc(i)%js
        ghost_sched%send_list(ns)%src_f = domain%bc(i)%fs
        ghost_sched%send_list(ns)%dst_i = domain%bc(i)%i
        ghost_sched%send_list(ns)%dst_j = domain%bc(i)%j
        ghost_sched%send_list(ns)%dst_f = domain%bc(i)%f
      end if
    end do

    ! Sort the entries by remote rank and build the rank groups
    call sort_entries_by_rank(ghost_sched%send_list, ns)
    call sort_entries_by_rank(ghost_sched%recv_list, nr)
    call build_rank_groups(ghost_sched%send_list, ns, ghost_sched%send_groups, ghost_sched%n_send_ranks)
    call build_rank_groups(ghost_sched%recv_list, nr, ghost_sched%recv_groups, ghost_sched%n_recv_ranks)

    ! Build the face groups for aggregated pack/unpack
    call build_face_groups_send(ghost_sched%send_list, ns, &
                                ghost_sched%send_faces, ghost_sched%n_send_faces)
    call build_face_groups_recv(ghost_sched%recv_list, nr, &
                                ghost_sched%recv_faces, ghost_sched%n_recv_faces)

    ! Preallocate the MPI buffers (reused at every RK step). PG_NCOLS=2 in 2D.
    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns > 0 .or. nr > 0) then
      block
        use mpi, only: MPI_STATUS_SIZE
        integer :: max_ranks
        integer, parameter :: PG_NCOLS = 2   ! columns 3..4 of Pg
        max_ranks = max(ghost_sched%n_send_ranks, ghost_sched%n_recv_ranks, 1)
        allocate(ghost_sched%P_send_buf(nprim*gc * max(ns,1)))
        allocate(ghost_sched%P_recv_buf(nprim*gc * max(nr,1)))
        allocate(ghost_sched%Pg_send_buf(nprim*PG_NCOLS * max(ns,1)))
        allocate(ghost_sched%Pg_recv_buf(nprim*PG_NCOLS * max(nr,1)))
        allocate(ghost_sched%send_req(max(ghost_sched%n_send_ranks,1)))
        allocate(ghost_sched%recv_req(max(ghost_sched%n_recv_ranks,1)))
        allocate(ghost_sched%mpi_stat(MPI_STATUS_SIZE, max_ranks))
      end block
    end if

    ! Initialize the persistent MPI requests for the P and R exchanges
    call init_persistent_P_requests()
    call init_persistent_R_requests()

    ! Diagnostics: n_send_ranks/n_recv_ranks are PER-RANK values (number of
    ! distinct neighbours of THIS rank). Printing only the root value would be
    ! misleading: it would read as "only N ranks communicate" whereas every
    ! rank communicates and root simply has N neighbours. Hence reduce over
    ! ALL ranks and report min/max/average + how many ranks have zero
    ! neighbours (the actual red flag).
    block
      use mpi, only: MPI_INTEGER, MPI_MIN, MPI_MAX, MPI_SUM, MPI_COMM_WORLD
      integer :: nloc, nmin, nmax, nsum, nzero, nzloc, ierr
      nloc  = ghost_sched%n_send_ranks
      nzloc = merge(1, 0, (ghost_sched%n_send_ranks == 0 .and. &
                           ghost_sched%n_recv_ranks == 0))
      call MPI_REDUCE(nloc,  nmin,  1, MPI_INTEGER, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(nloc,  nmax,  1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(nloc,  nsum,  1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(nzloc, nzero, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (mpi_is_root) then
        write(*,'(A,I0,A,I0,A,F5.1,A)') ' Ghost schedule: neighbours/rank min=', &
          nmin, ' max=', nmax, ' avg=', real(nsum)/real(mpi_size_), &
          ' (per-rank; every rank exchanges over the whole shared boundary)'
        if (nzero > 0) then
          write(*,'(A,I0,A,I0,A)') ' [WARNING] ', nzero, ' of ', mpi_size_, &
            ' ranks have ZERO ghost neighbours -- check partition/connectivity'
        end if
      end if
    end block
#endif

    ghost_sched%built = .true.
  end subroutine build_ghost_schedule


  !> Builds the prefiltered list of BC indices with Bm local to this rank.
  !> Avoids the O(nbound) scan in hot paths (Fill_Ghost_Cell, etc.).
  subroutine build_local_bc_index(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: i, n, ns

    ! Count local BC entries (local Bm) and local source entries (local Bs, type 101/201)
    n = 0; ns = 0
    do i = 1, domain%nbound
      if (is_local_block(domain%bc(i)%b)) n = n + 1
      if (domain%bc(i)%type == 101 .or. domain%bc(i)%type == 201) then
        if (is_local_block(domain%bc(i)%bs)) ns = ns + 1
      end if
    end do

    ! Build local_bc_idx (local Bm)
    domain%n_local_bc = n
    if (allocated(domain%local_bc_idx)) deallocate(domain%local_bc_idx)
    allocate(domain%local_bc_idx(n))

    ! Build local_bs_idx (local Bs, type 101/201 only — for the Pg computation)
    domain%n_local_bs = ns
    if (allocated(domain%local_bs_idx)) deallocate(domain%local_bs_idx)
    allocate(domain%local_bs_idx(ns))

    ! Fill both index arrays in a single pass
    n = 0; ns = 0
    do i = 1, domain%nbound
      if (is_local_block(domain%bc(i)%b)) then
        n = n + 1
        domain%local_bc_idx(n) = i
      end if
      if (domain%bc(i)%type == 101 .or. domain%bc(i)%type == 201) then
        if (is_local_block(domain%bc(i)%bs)) then
          ns = ns + 1
          domain%local_bs_idx(ns) = i
        end if
      end if
    end do
  end subroutine build_local_bc_index


  !> Frees the persistent MPI requests and deallocates the schedule arrays.
  subroutine cleanup_ghost_schedule()
    implicit none
#ifdef USE_MPI
    call cleanup_persistent_requests()
#endif
  end subroutine cleanup_ghost_schedule


#ifdef USE_MPI
  !> Initializes the persistent MPI requests for the P field exchange.
  !> Called once from build_ghost_schedule. Uses MPI_SEND_INIT/RECV_INIT.
  subroutine init_persistent_P_requests()
    use Q2D_Global_m, only: gc, nprim
    use mpi

    implicit none
    integer :: r, ierr, tag, buf_pos, entry_size

    entry_size = nprim * gc

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(ghost_sched%P_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%P_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

    ! Initialize the persistent receives
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank
      call MPI_RECV_INIT(ghost_sched%P_recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%P_recv_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Initialize the persistent sends
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_
      call MPI_SEND_INIT(ghost_sched%P_send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%P_send_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ghost_sched%persistent_P_init = .true.
  end subroutine init_persistent_P_requests


  !> Initializes the persistent MPI requests for the R field exchange.
  !> Reuses P_send_buf/P_recv_buf (same size). Uses distinct tags (offset 30000).
  subroutine init_persistent_R_requests()
    use Q2D_Global_m, only: gc, nprim
    use mpi

    implicit none
    integer :: r, ierr, tag, buf_pos, entry_size
    integer, parameter :: TAG_OFFSET = 30000

    entry_size = nprim * gc

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(ghost_sched%R_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%R_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

    ! Initialize the persistent receives for R
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + TAG_OFFSET
      call MPI_RECV_INIT(ghost_sched%P_recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%R_recv_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Initialize the persistent sends for R
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + TAG_OFFSET
      call MPI_SEND_INIT(ghost_sched%P_send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%R_send_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ghost_sched%persistent_R_init = .true.
  end subroutine init_persistent_R_requests


  !> Frees all persistent MPI requests.
  subroutine cleanup_persistent_requests()
    use mpi
    implicit none
    integer :: ierr, r

    if (ghost_sched%persistent_P_init) then
      do r = 1, ghost_sched%n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%P_send_req_pers(r), ierr)
      end do
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%P_recv_req_pers(r), ierr)
      end do
      ghost_sched%persistent_P_init = .false.
    end if

    if (ghost_sched%persistent_R_init) then
      do r = 1, ghost_sched%n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%R_send_req_pers(r), ierr)
      end do
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%R_recv_req_pers(r), ierr)
      end do
      ghost_sched%persistent_R_init = .false.
    end if
  end subroutine cleanup_persistent_requests


  !> Sorts the ghost entries by remote_rank with insertion sort (the lists are small).
  subroutine sort_entries_by_rank(list, n)
    implicit none
    type(ghost_entry_type), intent(inout) :: list(:)
    integer, intent(in) :: n
    ! Local variables
    integer :: i, j
    type(ghost_entry_type) :: tmp

    do i = 2, n
      tmp = list(i)
      j = i - 1
      do while (j >= 1 .and. list(j)%remote_rank > tmp%remote_rank)
        list(j+1) = list(j)
        j = j - 1
      end do
      list(j+1) = tmp
    end do
  end subroutine sort_entries_by_rank


  !> Builds the rank groups from a sorted entry list.
  subroutine build_rank_groups(list, n, groups, n_groups)
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(rank_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    ! Local variables
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count the distinct ranks
    ng = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%rank   = list(1)%remote_rank
    groups(1)%offset = 1
    groups(1)%count  = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) then
        ng = ng + 1
        groups(ng)%rank   = list(i)%remote_rank
        groups(ng)%offset = i
        groups(ng)%count  = 1
      else
        groups(ng)%count = groups(ng)%count + 1
      end if
    end do
  end subroutine build_rank_groups


  !> Builds the face groups from the send entry list.
  !> Groups consecutive entries sharing (local_block, src_f).
  !> Guide offsets di, dj taken from guide(:, 1..2).
  subroutine build_face_groups_send(list, n, groups, n_groups)
    use Q2D_Parameters_m, only: guide
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(face_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count the face groups
    ng = 1
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%src_f /= list(i-1)%src_f) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%block = list(1)%local_block
    groups(1)%face  = list(1)%src_f
    groups(1)%entry_start = 1
    groups(1)%entry_count = 1
    groups(1)%di = guide(list(1)%src_f, 1)
    groups(1)%dj = guide(list(1)%src_f, 2)
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%src_f /= list(i-1)%src_f) then
        ng = ng + 1
        groups(ng)%block = list(i)%local_block
        groups(ng)%face  = list(i)%src_f
        groups(ng)%entry_start = i
        groups(ng)%entry_count = 1
        groups(ng)%di = guide(list(i)%src_f, 1)
        groups(ng)%dj = guide(list(i)%src_f, 2)
      else
        groups(ng)%entry_count = groups(ng)%entry_count + 1
      end if
    end do
  end subroutine build_face_groups_send


  !> Builds the face groups from the receive entry list.
  !> Groups consecutive entries sharing (local_block, dst_f).
  subroutine build_face_groups_recv(list, n, groups, n_groups)
    use Q2D_Parameters_m, only: guide
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(face_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count the face groups
    ng = 1
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%dst_f /= list(i-1)%dst_f) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%block = list(1)%local_block
    groups(1)%face  = list(1)%dst_f
    groups(1)%entry_start = 1
    groups(1)%entry_count = 1
    groups(1)%di = guide(list(1)%dst_f, 1)
    groups(1)%dj = guide(list(1)%dst_f, 2)
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%dst_f /= list(i-1)%dst_f) then
        ng = ng + 1
        groups(ng)%block = list(i)%local_block
        groups(ng)%face  = list(i)%dst_f
        groups(ng)%entry_start = i
        groups(ng)%entry_count = 1
        groups(ng)%di = guide(list(i)%dst_f, 1)
        groups(ng)%dj = guide(list(i)%dst_f, 2)
      else
        groups(ng)%entry_count = groups(ng)%entry_count + 1
      end if
    end do
  end subroutine build_face_groups_recv
#endif


  !> Starts the persistent P receives. To be called from !$omp single.
  subroutine exchange_ghost_P_post_recv(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_post_recv()
#endif
  end subroutine exchange_ghost_P_post_recv


  !> Packs the P send buffer. To be called from !$omp do over the face groups.
  !> n_send_faces is accessible via ghost_sched%n_send_faces.
  subroutine exchange_ghost_P_pack(domain, fg_start, fg_end)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    implicit none
    type(Q2D_domain_type), intent(in) :: domain
    integer, intent(in) :: fg_start, fg_end
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_pack_faces(domain, fg_start, fg_end)
#endif
  end subroutine exchange_ghost_P_pack


  !> Starts the persistent P sends. To be called from !$omp single once the pack is complete.
  subroutine exchange_ghost_P_post_send(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_post_send()
#endif
  end subroutine exchange_ghost_P_post_send


  !> Waits for the P receives and unpacks. Wait from !$omp single, unpack from !$omp do.
  subroutine exchange_ghost_P_wait_unpack(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_wait_recv()
    call exchange_P_unpack_faces(domain, 1, ghost_sched%n_recv_faces)
#endif
  end subroutine exchange_ghost_P_wait_unpack


  !> Waits for the P sends to complete. To be called from !$omp single.
  subroutine exchange_ghost_P_wait_send(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_wait_send()
#endif
  end subroutine exchange_ghost_P_wait_send


  !> Exchanges the R (residual) ghost-cell data between MPI ranks.
  subroutine exchange_ghost_R(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_R_field_begin(domain)
    call exchange_R_field_end(domain)
#endif
  end subroutine exchange_ghost_R


  !> Gathers on rank 0 the interior P data of all blocks for I/O.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL to overlap the communications.
  !> On root, temporarily allocates P for the remote blocks if not already allocated.
  !> P of rank 3 (nprim, n1, n2).
  subroutine gather_P_to_root(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nprim, gc

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: b, n1, n2

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute the per-block sizes and the total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        blk_ncells(b) = nprim * n1 * n2
        blk_offset(b) = total_buf_size
        if (mpi_is_root) then
          if (.not. is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        else
          if (is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        end if
      end do
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! Recompute the offsets contiguously
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do

      nreq = 0
      if (mpi_is_root) then
        ! Root: start all receives for the remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack the local blocks and start all sends
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%P(:, 1:n1, 1:n2), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all requests to complete
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack the received data (allocate P if it had been deallocated)
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            if (.not. allocated(domain%blk(b)%P)) &
              allocate(domain%blk(b)%P(nprim, 1-gc:n1+gc, 1-gc:n2+gc))
            domain%blk(b)%P(:, 1:n1, 1:n2) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [nprim, n1, n2])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine gather_P_to_root


  !> Gathers on root the diagnostic arrays (R, dtlocal, beta) from all ranks.
  !> On root, the remote blocks may have these arrays deallocated by
  !> deallocate_remote_computation_data: here they are reallocated to receive
  !> the data. All ranks must call this routine (collective).
  !> R of rank 3 (nprim, n1, n2); dtlocal and beta of rank 2 (n1, n2).
  subroutine gather_diagnostic_to_root(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nprim

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: b, n1, n2

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, offset, total_buf_size
      integer :: ncells_R, ncells_dt
      integer, allocatable :: reqs(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)
      ! Per block: nprim*n1*n2 for R, n1*n2 for dtlocal, n1*n2 for beta
      integer, allocatable :: blk_offset(:), blk_size(:)

      allocate(blk_offset(domain%nb), blk_size(domain%nb))

      ! Compute the buffer sizes
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        blk_size(b) = nprim * n1 * n2 + 2 * n1 * n2  ! R + dtlocal + beta
      end do

      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_size(b)
        else
          if (is_local_block(b)) offset = offset + blk_size(b)
        end if
      end do
      total_buf_size = offset
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! On root: ensure R, dtlocal, beta are allocated for the remote blocks
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            if (.not. allocated(domain%blk(b)%R)) &
              allocate(domain%blk(b)%R(nprim, n1, n2))
            if (.not. allocated(domain%blk(b)%dtlocal)) &
              allocate(domain%blk(b)%dtlocal(n1, n2))
            if (.not. allocated(domain%blk(b)%beta)) &
              allocate(domain%blk(b)%beta(n1, n2))
          end if
        end do
      end if

      nreq = 0
      if (mpi_is_root) then
        ! Root: start the receives for the remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_size(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), 1000+b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack the local blocks and send
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            ncells_R  = nprim * n1 * n2
            ncells_dt = n1 * n2
            offset = blk_offset(b)
            buf(offset+1 : offset+ncells_R) = &
              reshape(domain%blk(b)%R(:, 1:n1, 1:n2), [ncells_R])
            offset = offset + ncells_R
            buf(offset+1 : offset+ncells_dt) = &
              reshape(domain%blk(b)%dtlocal(1:n1, 1:n2), [ncells_dt])
            offset = offset + ncells_dt
            buf(offset+1 : offset+ncells_dt) = &
              reshape(domain%blk(b)%beta(1:n1, 1:n2), [ncells_dt])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_size(b), MPI_DOUBLE_PRECISION, &
                           0, 1000+b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all requests to complete
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack the received data
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            ncells_R  = nprim * n1 * n2
            ncells_dt = n1 * n2
            offset = blk_offset(b)
            domain%blk(b)%R(:, 1:n1, 1:n2) = &
              reshape(buf(offset+1 : offset+ncells_R), [nprim, n1, n2])
            offset = offset + ncells_R
            domain%blk(b)%dtlocal(1:n1, 1:n2) = &
              reshape(buf(offset+1 : offset+ncells_dt), [n1, n2])
            offset = offset + ncells_dt
            domain%blk(b)%beta(1:n1, 1:n2) = &
              reshape(buf(offset+1 : offset+ncells_dt), [n1, n2])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_size)
    end block
#endif
  end subroutine gather_diagnostic_to_root


  !> Scatters the block P data from rank 0 to the owner ranks after reading.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL to overlap the communications.
  !> P of rank 3 (nprim, n1, n2).
  subroutine scatter_P_from_root(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nprim

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: b, n1, n2

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute the per-block sizes and the total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        blk_ncells(b) = nprim * n1 * n2
      end do

      ! Compute the contiguous buffer offsets
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do
      total_buf_size = offset
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      nreq = 0
      if (mpi_is_root) then
        ! Root: pack the remote blocks and start all sends
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%P(:, 1:n1, 1:n2), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: start all receives for the local blocks
        do b = 1, domain%nb
          if (is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all requests to complete
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Non-root: unpack the received data
      if (.not. mpi_is_root) then
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            domain%blk(b)%P(:, 1:n1, 1:n2) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [nprim, n1, n2])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine scatter_P_from_root


  !> Exchanges the Pg(:,3:4) stencil data for inter-rank viscous connection BCs.
  !> Messages are aggregated per remote rank.
  !> PG_NCOLS=2: columns 3..4 of Pg.
  subroutine exchange_ghost_Pg(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_Pg_field(domain)
#endif
  end subroutine exchange_ghost_Pg


#ifdef USE_MPI
  !> Internal: exchanges the Pg(:,3:4) stencil with aggregated messages.
  subroutine exchange_Pg_field(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nprim
    use mpi

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: ns, nr, i, r, ierr, tag, bc_idx, buf_pos
    integer :: entry_size
    integer, parameter :: PG_NCOLS = 2   ! columns 3..4 of Pg

    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns == 0 .and. nr == 0) return

    entry_size = nprim * PG_NCOLS

    ! Start one receive per remote rank (aggregated)
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + 10000  ! unique per-rank tag for Pg
      call MPI_IRECV(ghost_sched%Pg_recv_buf(buf_pos), &
                     ghost_sched%recv_groups(r)%count * entry_size, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%recv_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%recv_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Pack all entries of each rank into a contiguous buffer, then one send per rank
    do r = 1, ghost_sched%n_send_ranks
      do i = ghost_sched%send_groups(r)%offset, &
             ghost_sched%send_groups(r)%offset + ghost_sched%send_groups(r)%count - 1
        bc_idx = ghost_sched%send_list(i)%bc_idx
        buf_pos = (i - 1) * entry_size
        ghost_sched%Pg_send_buf(buf_pos+1         : buf_pos+  nprim) = domain%bc(bc_idx)%Pg(:,3)
        ghost_sched%Pg_send_buf(buf_pos+  nprim+1 : buf_pos+2*nprim) = domain%bc(bc_idx)%Pg(:,4)
      end do
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + 10000  ! sender rank as tag for Pg
      call MPI_ISEND(ghost_sched%Pg_send_buf(buf_pos), &
                     ghost_sched%send_groups(r)%count * entry_size, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%send_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%send_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Wait for all receives
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%recv_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if

    ! Unpack into bc%Pg(:,3:4)
    do i = 1, nr
      bc_idx = ghost_sched%recv_list(i)%bc_idx
      buf_pos = (i - 1) * entry_size
      domain%bc(bc_idx)%Pg(:,3) = ghost_sched%Pg_recv_buf(buf_pos+1         : buf_pos+  nprim)
      domain%bc(bc_idx)%Pg(:,4) = ghost_sched%Pg_recv_buf(buf_pos+  nprim+1 : buf_pos+2*nprim)
    end do

    ! Wait for all sends to complete
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%send_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_Pg_field


  !> Starts the persistent receives for the P field. To be called from !$omp single.
  subroutine exchange_P_post_recv()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%P_recv_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_post_recv


  !> Packs the P send buffer. Safe to call from !$omp do over the
  !> face groups. Source: P(:, Is+di*(g-1), Js+dj*(g-1)).
  subroutine exchange_P_pack_faces(domain, fg_start, fg_end)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: gc, nprim
    implicit none
    type(Q2D_domain_type), intent(in) :: domain
    integer, intent(in) :: fg_start, fg_end
    ! Local variables
    integer :: fg, i, g, b, buf_pos
    integer :: Is, Js, di, dj
    integer :: entry_size

    entry_size = nprim * gc

    do fg = fg_start, fg_end
      b  = ghost_sched%send_faces(fg)%block
      di = ghost_sched%send_faces(fg)%di
      dj = ghost_sched%send_faces(fg)%dj
      do i = ghost_sched%send_faces(fg)%entry_start, &
             ghost_sched%send_faces(fg)%entry_start + ghost_sched%send_faces(fg)%entry_count - 1
        Is = ghost_sched%send_list(i)%src_i
        Js = ghost_sched%send_list(i)%src_j
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          ghost_sched%P_send_buf(buf_pos + (g-1)*nprim+1 : buf_pos + g*nprim) = &
            domain%blk(b)%P(:, Is+di*(g-1), Js+dj*(g-1))
        end do
      end do
    end do
  end subroutine exchange_P_pack_faces


  !> Starts the persistent sends for the P field. To be called from !$omp single after the pack.
  subroutine exchange_P_post_send()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%P_send_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_post_send


  !> Waits for the P receives. To be called from !$omp single.
  subroutine exchange_P_wait_recv()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%P_recv_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_wait_recv


  !> Unpacks the P receive buffer. Safe to call from !$omp do
  !> over the face groups. Destination: P(:, Id-di*g, Jd-dj*g).
  subroutine exchange_P_unpack_faces(domain, fg_start, fg_end)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: gc, nprim
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer, intent(in) :: fg_start, fg_end
    ! Local variables
    integer :: fg, i, g, b, buf_pos
    integer :: Id, Jd, di, dj
    integer :: entry_size

    entry_size = nprim * gc

    do fg = fg_start, fg_end
      b  = ghost_sched%recv_faces(fg)%block
      di = ghost_sched%recv_faces(fg)%di
      dj = ghost_sched%recv_faces(fg)%dj
      do i = ghost_sched%recv_faces(fg)%entry_start, &
             ghost_sched%recv_faces(fg)%entry_start + ghost_sched%recv_faces(fg)%entry_count - 1
        Id = ghost_sched%recv_list(i)%dst_i
        Jd = ghost_sched%recv_list(i)%dst_j
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          domain%blk(b)%P(:, Id-di*g, Jd-dj*g) = &
            ghost_sched%P_recv_buf(buf_pos + (g-1)*nprim+1 : buf_pos + g*nprim)
        end do
      end do
    end do
  end subroutine exchange_P_unpack_faces


  !> Waits for the P sends to complete. To be called from !$omp single.
  subroutine exchange_P_wait_send()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%P_send_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_wait_send


  !> Internal: starts the non-blocking R field exchange with persistent requests.
  subroutine exchange_R_field_begin(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: gc, nprim
    use mpi

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: ns, nr, i, g, b, fg, ierr, buf_pos
    integer :: Is, Js, di, dj
    integer :: entry_size

    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns == 0 .and. nr == 0) return

    entry_size = nprim * gc

    ! Start the persistent receives for R
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%R_recv_req_pers, ierr)
      call check_mpi_error(ierr)
    end if

    ! Pack using the face groups (reads R instead of P)
    do fg = 1, ghost_sched%n_send_faces
      b  = ghost_sched%send_faces(fg)%block
      di = ghost_sched%send_faces(fg)%di
      dj = ghost_sched%send_faces(fg)%dj
      do i = ghost_sched%send_faces(fg)%entry_start, &
             ghost_sched%send_faces(fg)%entry_start + ghost_sched%send_faces(fg)%entry_count - 1
        Is = ghost_sched%send_list(i)%src_i
        Js = ghost_sched%send_list(i)%src_j
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          ghost_sched%P_send_buf(buf_pos + (g-1)*nprim+1 : buf_pos + g*nprim) = &
            domain%blk(b)%R(:, Is+di*(g-1), Js+dj*(g-1))
        end do
      end do
    end do

    ! Start the persistent sends for R
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%R_send_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_R_field_begin


  !> Internal: completes the non-blocking exchange and unpacks R into the ghost cells.
  subroutine exchange_R_field_end(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: gc, nprim
    use mpi

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: nr, i, g, b, fg, ierr, buf_pos
    integer :: Id, Jd, di, dj
    integer :: entry_size

    nr = ghost_sched%n_recv
    if (ghost_sched%n_send == 0 .and. nr == 0) return

    entry_size = nprim * gc

    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%R_recv_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if

    ! Unpack using the face groups
    do fg = 1, ghost_sched%n_recv_faces
      b  = ghost_sched%recv_faces(fg)%block
      di = ghost_sched%recv_faces(fg)%di
      dj = ghost_sched%recv_faces(fg)%dj
      do i = ghost_sched%recv_faces(fg)%entry_start, &
             ghost_sched%recv_faces(fg)%entry_start + ghost_sched%recv_faces(fg)%entry_count - 1
        Id = ghost_sched%recv_list(i)%dst_i
        Jd = ghost_sched%recv_list(i)%dst_j
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          domain%blk(b)%R(:, Id-di*g, Jd-dj*g) = &
            ghost_sched%P_recv_buf(buf_pos + (g-1)*nprim+1 : buf_pos + g*nprim)
        end do
      end do
    end do

    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%R_send_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_R_field_end
#endif


  !> MPI barrier to synchronize the ranks after I/O operations.
  !> Prevents non-root ranks from advancing to the next MPI collective
  !> while root is still writing to file.
  subroutine mpi_io_barrier()
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    end block
#endif
  end subroutine mpi_io_barrier


  !> Fills the inter-rank ghost stencil Pg(:, 1:gc) for the source face Fm.
  subroutine Ghost_Interrank( Im, Jm, Fm, blkm, Pg )
    use Q2D_Global_m,         only: gc, nprim
    use Q2D_Parameters_m,     only: guide
    use Q2D_Advanced_Types_m, only: Q2D_block_type
    implicit none
    integer, intent(in)                  :: Im, Jm, Fm
    type(Q2D_block_type), intent(inout)  :: blkm
    real(R8), intent(inout)              :: Pg(nprim, 4)
    integer                              :: g, Ig, Jg

    do g = 1, gc
        Ig = Im - guide(Fm,1)*g
        Jg = Jm - guide(Fm,2)*g
        Pg(:,g) = blkm % P(:,Ig,Jg)
    enddo
  end subroutine Ghost_Interrank


end module Q2D_Mod_GhostExchange
