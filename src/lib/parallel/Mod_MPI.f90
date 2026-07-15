!>@brief Modulo di infrastruttura MPI per Q2D.
!> Fornisce la gestione dell'ambiente MPI, il partizionamento dei blocchi e
!> i wrapper delle riduzioni globali. Compilato senza USE_MPI, tutte le
!> routine sono no-op o operazioni identità (fallback seriale).
!>
!> L'infrastruttura MPI è indipendente dalla dimensionalità del solutore.
!> Le procedure pubbliche sono sotto `#ifdef USE_MPI` con fallback identità
!> nel ramo #else.
module Q2D_Mod_MPI
#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env, only: R8 => real64, I4 => int32

  implicit none
  private

  ! --- Stato pubblico ---
  integer, public :: mpi_rank_ = 0         !< Rank di questo processo (0 in seriale)
  integer, public :: mpi_size_ = 1         !< Numero totale di processi MPI
  logical, public :: mpi_is_root = .true.  !< Vero su rank 0

  ! --- Mappa blocco → rank ---
  integer, allocatable, public :: block_owner(:)     !< block_owner(b) = rank proprietario del blocco b
  integer, allocatable, public :: local_block_ids(:) !< ID globali dei blocchi posseduti da questo rank
  integer, public :: n_local_blocks = 0

  ! --- Procedure pubbliche ---
  public :: mpi_init_env, mpi_finalize_env
  public :: is_local_block
  public :: partition_blocks
  public :: mpi_allreduce_sum_r8, mpi_allreduce_min_r8
  public :: mpi_allreduce_sum_r8_array
  public :: mpi_reduce_sum_r8, mpi_reduce_sum_r8_array
  public :: mpi_allreduce_norm2
  public :: mpi_bcast_logical, mpi_bcast_integer
  public :: mpi_abort_all
  public :: check_mpi_error

contains


  !> Inizializza l'ambiente MPI con supporto ai thread. Da chiamare all'inizio del main.
  !> Usa MPI_THREAD_FUNNELED: solo il thread master effettua chiamate MPI.
  subroutine mpi_init_env()
#ifdef USE_MPI
    integer :: ierr, provided
    logical :: already

    call MPI_Initialized(already, ierr)
    if (.not. already) then
      call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierr)
      call check_mpi_error(ierr)
      if (provided < MPI_THREAD_FUNNELED) then
        write(*,'(A)') ' ERROR: MPI does not support the required threading level (MPI_THREAD_FUNNELED)'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end if
    end if
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank_, ierr)
    call check_mpi_error(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size_, ierr)
    call check_mpi_error(ierr)
    mpi_is_root = (mpi_rank_ == 0)
#else
    mpi_rank_ = 0
    mpi_size_ = 1
    mpi_is_root = .true.
#endif
  end subroutine mpi_init_env


  !> Finalizza l'ambiente MPI. Da chiamare alla fine del main.
  subroutine mpi_finalize_env()
#ifdef USE_MPI
    integer :: ierr
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
#endif
  end subroutine mpi_finalize_env


  !> Verifica se il blocco b è posseduto da questo rank.
  logical function is_local_block(b)
    integer, intent(in) :: b

    if (.not. allocated(block_owner)) then
      is_local_block = .true.  ! modalità seriale o partizionamento non ancora eseguito
      return
    end if
    is_local_block = (block_owner(b) == mpi_rank_)
  end function is_local_block


  !> Partiziona nb blocchi tra i rank MPI con bilanciamento di carico greedy.
  !> blk_ncells(b) = numero totale di celle del blocco b.
  subroutine partition_blocks(nb, blk_ncells)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    ! Variabili locali
    integer :: b, r, nloc
    integer, allocatable :: rank_load(:)

    if (allocated(block_owner)) deallocate(block_owner)
    if (allocated(rank_load)) deallocate(rank_load)
    allocate(block_owner(nb))
    allocate(rank_load(0:mpi_size_-1))
    rank_load = 0

    ! Greedy: assegna ogni blocco al rank con il carico minore
    r = 0
    do b = 1, nb
      r = minloc(rank_load, dim=1) - 1   ! rank con carico minimo (indicizzato da 0)
      block_owner(b) = r
      rank_load(r) = rank_load(r) + blk_ncells(b)
    end do

    ! Costruisce la lista dei blocchi locali
    nloc = count(block_owner == mpi_rank_)
    n_local_blocks = nloc
    if (allocated(local_block_ids)) deallocate(local_block_ids)
    allocate(local_block_ids(nloc))
    nloc = 0
    do b = 1, nb
      if (block_owner(b) == mpi_rank_) then
        nloc = nloc + 1
        local_block_ids(nloc) = b
      end if
    end do

    deallocate(rank_load)
  end subroutine partition_blocks


  !> MPI_ALLREDUCE con MPI_SUM per uno scalare real(R8).
  subroutine mpi_allreduce_sum_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_sum_r8


  !> MPI_ALLREDUCE con MPI_MIN per uno scalare real(R8).
  subroutine mpi_allreduce_min_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_min_r8


  !> MPI_ALLREDUCE con MPI_SUM per un array di real(R8).
  !> Opera in-place: arr è sovrascritto con il risultato globale.
  subroutine mpi_allreduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(MPI_IN_PLACE, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_allreduce_sum_r8_array


  !> Norma L2 distribuita: somma locale dei quadrati + MPI_ALLREDUCE SUM + sqrt.
  function mpi_allreduce_norm2(x, n) result(global_norm)
    integer, intent(in) :: n
    real(R8), intent(in) :: x(n)
    real(R8) :: global_norm
    real(R8) :: local_ss

    local_ss = dot_product(x, x)
#ifdef USE_MPI
    call mpi_allreduce_sum_r8(local_ss, global_norm)
    global_norm = sqrt(global_norm)
#else
    global_norm = sqrt(local_ss)
#endif
  end function mpi_allreduce_norm2


  !> MPI_REDUCE con MPI_SUM per uno scalare real(R8) verso il root (rank 0).
  subroutine mpi_reduce_sum_r8(val)
    real(R8), intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    if (mpi_rank_ == 0) then
      call MPI_REDUCE(MPI_IN_PLACE, val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    else
      call MPI_REDUCE(val, val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8


  !> MPI_REDUCE con MPI_SUM per un array di real(R8) verso il root (rank 0).
  !> Il risultato è valido solo su root. Più efficiente di ALLREDUCE quando
  !> solo root necessita del risultato (es. per l'output dei residui).
  subroutine mpi_reduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    if (mpi_rank_ == 0) then
      call MPI_REDUCE(MPI_IN_PLACE, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    else
      call MPI_REDUCE(arr, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8_array


  !> Trasmette un valore logico da root (rank 0) a tutti i rank.
  subroutine mpi_bcast_logical(val)
    logical, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_logical


  !> Trasmette un valore intero da root (rank 0) a tutti i rank.
  subroutine mpi_bcast_integer(val)
    integer, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_integer


  !> Controlla il codice di ritorno MPI. In caso di errore abortisce tutti i rank.
  subroutine check_mpi_error(ierr)
    integer, intent(in) :: ierr
#ifdef USE_MPI
    integer :: abort_ierr
    if (ierr /= 0) then
      write(*,'(A,I4,A,I6)') ' [RANK ', mpi_rank_, '] MPI error code: ', ierr
      call MPI_ABORT(MPI_COMM_WORLD, 1, abort_ierr)
    end if
#endif
  end subroutine check_mpi_error


  !> Abort consapevole di MPI: stampa il messaggio, poi abortisce tutti i rank.
  subroutine mpi_abort_all(message)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer :: ierr
    write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] ABORT: ', trim(message)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#else
    write(*,'(A,A)') ' ABORT: ', trim(message)
    error stop
#endif
  end subroutine mpi_abort_all


end module Q2D_Mod_MPI
