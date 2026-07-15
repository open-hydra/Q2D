module Q2D_IO_Probes
  use iso_fortran_env,      only: I4 => int32, R8 => real64
  use Q2D_Parameters_m,     only: clen, llen
  use Q2D_Base_Types_m,     only: Q2D_vector_2D_type
  use Q2D_Config_Types_m,   only: obj_io_probes
  use Q2D_Advanced_Types_m, only: Q2D_domain_type

  implicit none
  private
  public :: Setup_Probes
  public :: Write_Probes_Data
  public :: nprobes

  type :: real_ptr
    real(R8), pointer :: p => null()
  end type real_ptr

  type :: obj_probe_t
    type(Q2D_vector_2D_type)         :: location          ! coordinate (x, y)
    integer, dimension(3)            :: ilocation         ! indici (b, i, j)
    integer                          :: nvar
    character(len=llen), allocatable :: names(:)
    type(real_ptr),      allocatable :: variables(:)
    real(R8),            pointer     :: P(:)
    real(R8)                         :: T, M
    real(R8)                         :: dtime
    integer                          :: ntime
    integer                          :: diter
    integer                          :: unit
  contains
    private
    procedure, pass(self) :: Place
    procedure, pass(self) :: Update
  end type obj_probe_t

  integer :: nprobes = 0
  type(obj_probe_t), allocatable, target :: probe(:)

contains

  subroutine Setup_Probes ( domain, newrun )
    use Q2D_Mod_MPI, only: is_local_block
    use strings,     only: parse
    implicit none
    type(Q2D_domain_type), intent(in) :: domain
    logical,               intent(in) :: newrun
    ! Variabili locali
    integer             :: i, ii, error
    character(len=clen) :: string(8)

    if (.not. allocated(obj_io_probes)) then
      nprobes = 0
      return
    end if

    nprobes = size(obj_io_probes)
    if (nprobes == 0) return

    if (allocated(probe)) deallocate(probe)
    allocate(probe(1:nprobes))

    do i = 1, nprobes

      obj_io_probes(i)%file = 'OUTPUT/'//trim(obj_io_probes(i)%file)//'.txt'

      ! Nomi delle variabili
      string = ''
      probe(i)%nvar = 0
      call parse(obj_io_probes(i)%varnames, ' ', string)
      do ii = 1, size(string)
        if (string(ii) /= '') probe(i)%nvar = probe(i)%nvar + 1
      end do
      allocate(probe(i)%names(1:probe(i)%nvar))
      do ii = 1, probe(i)%nvar
        probe(i)%names(ii) = trim(string(ii))
      end do

      ! Frequenze di scrittura
      probe(i)%dtime = obj_io_probes(i)%dtime
      probe(i)%diter = obj_io_probes(i)%diter
      probe(i)%ntime = 1

      ! Posizione — 3 componenti per ilocation, 2 per location
      probe(i)%location%c = obj_io_probes(i)%loc
      probe(i)%ilocation  = obj_io_probes(i)%iloc

      if (sum(obj_io_probes(i)%iloc) == 0) call probe(i)%Place(domain)

      if (.not. is_local_block(probe(i)%ilocation(1))) cycle

      call Assign_Variables(probe(i), domain)

      if (.not. newrun) then
        open(newunit=probe(i)%unit, file=trim(obj_io_probes(i)%file), status='OLD', iostat=error)
        if (error /= 0) then
          obj_io_probes(i)%error_message = &
            "[ERROR] You restarted from an old solution but the probes files were not found."
          return
        end if
        error = 0
        do while (error == 0)
          read(probe(i)%unit, *, iostat=error)
        end do
        backspace(probe(i)%unit)
      else
        open(newunit=probe(i)%unit, file=trim(obj_io_probes(i)%file), status='REPLACE', iostat=error)
      end if

    end do

  end subroutine Setup_Probes


  !> Determina gli indici di cella della sonda a partire dalle coordinate (x, y).
  subroutine Place ( self, domain )
    implicit none
    class(obj_probe_t)                 :: self
    type(Q2D_domain_type), intent(in)  :: domain
    ! Variabili locali
    integer  :: i, j, b
    real(R8) :: d0, d

    d0 = huge(1.0_R8)
    do b = 1, domain%nb
      do j = 0, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          d = norm2(self%location%c - domain%blk(b)%node(i, j)%c)
          if (d < d0) then
            d0 = d
            self%ilocation = [b, i, j]
          end if
        end do
      end do
    end do
  end subroutine Place


  subroutine Assign_Variables ( probe_, domain )
    use Q2D_Global_m, only: nsc, nu, nv, np
    use IR_Precision, only: str
    implicit none
    type(obj_probe_t),     intent(inout), target :: probe_
    type(Q2D_domain_type), intent(in),    target :: domain
    ! Variabili locali
    integer :: v, s, b, i, j

    b = probe_%ilocation(1)
    i = probe_%ilocation(2)
    j = probe_%ilocation(3)

    probe_%P => domain%blk(b)%P(:, i, j)
    allocate(probe_%variables(1:probe_%nvar))

    do v = 1, probe_%nvar
      do s = 1, nsc
        if (probe_%names(v) == 'rho('//trim(str(.true.,s))//')') &
          probe_%variables(v)%p => probe_%P(s)
      end do
      ! Nessuna componente "w".
      if (probe_%names(v) == 'u') probe_%variables(v)%p => probe_%P(nu)
      if (probe_%names(v) == 'v') probe_%variables(v)%p => probe_%P(nv)
      if (probe_%names(v) == 'p') probe_%variables(v)%p => probe_%P(np)
      if (probe_%names(v) == 'T') probe_%variables(v)%p => probe_%T
      if (probe_%names(v) == 'M') probe_%variables(v)%p => probe_%M
    end do

  end subroutine Assign_Variables


  subroutine Update ( self )
    use Q2D_Global_m,            only: nsc, nu, nv, np
    use FLINT_Lib_Thermodynamic, only: f_Rtot, f_ss, EOS
    implicit none
    class(obj_probe_t), intent(inout) :: self
    ! Variabili locali
    real(R8) :: Rtot, speed_of_sound

    Rtot = f_Rtot(self%P(1:nsc))
    speed_of_sound = f_ss(self%P(1:nsc), self%P(np), sum(self%P(1:nsc)), Rtot)
    self%T = EOS(p=self%P(np), rho=sum(self%P(1:nsc)), R=Rtot)
    ! Mach dalle sole componenti di velocità (u, v); nessuna componente w.
    self%M = norm2(self%P(nu:nv)) / speed_of_sound

  end subroutine Update


  subroutine Write_Probes_Data ( iter, time )
    use Q2D_Mod_MPI, only: is_local_block
    implicit none
    integer,  intent(in) :: iter
    real(R8), intent(in) :: time
    ! Variabili locali
    integer :: i, v

    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      if (mod(iter, probe(i)%diter) == 0) then
        call probe(i)%update
        write(probe(i)%unit, *) iter, (probe(i)%variables(v)%p, v=1, probe(i)%nvar)
      else if (time >= probe(i)%dtime * probe(i)%ntime) then
        call probe(i)%update
        write(probe(i)%unit, *) time, (probe(i)%variables(v)%p, v=1, probe(i)%nvar)
        probe(i)%ntime = probe(i)%ntime + 1
      end if
    end do

  end subroutine Write_Probes_Data

end module Q2D_IO_Probes
