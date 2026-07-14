module Q2D_IO_LineProbes
  ! Q2D-extra (peculiarità #5): line-probes (sampling along grid lines i/j
  ! fixed-index). MOSE has no line-probes. Registry-driven like the regular probes:
  ! the per-line-probe config is registered (Q2D_Read_IO::Register_LineProbes) and
  ! loaded (Backend_INI::Load_Ini) into obj_io_lineprobes (Config_Types); the
  ! section-level restart goes to obj_io%lineprobe_restart (Scan_Ini). Setup_LineProbes
  ! reads the loaded values — no fini%get bypass (resolves TODO C-11.d).
  use iso_fortran_env,      only: I4 => int32, R8 => real64
  use Q2D_Parameters_m,     only: clen, llen, hlen
  use Q2D_Advanced_Types_m, only: Q2D_domain_type

  implicit none
  private
  public :: Setup_LineProbes
  public :: Write_LineProbes_Data
  public :: Close_LineProbes
  public :: nlineprobes

  type :: real_array_ptr
    real(R8), pointer :: p(:) => null()
  end type real_array_ptr

  type :: obj_lineprobe_t
    character(len=llen)                :: name
    integer, dimension(2)              :: ilocation       ! (block, direction)  direction: 1 = i fixed, 2 = j fixed
    integer                            :: fixed_index     ! the constant i or j index defining the line
    integer                            :: npoints
    integer                            :: nvar
    character(len=clen), allocatable   :: names(:)
    type(real_array_ptr), allocatable  :: variables(:)
    real(R8), pointer                  :: P(:,:)          ! Pointer to (nprim, npoints) slice
    real(R8), allocatable              :: T(:), M(:)
    real(R8)                           :: dtime
    integer                            :: ntime
    integer                            :: diter
    real(R8), allocatable              :: nodes_lower(:,:) ! degenerate-strip node coords j=0 row (2, 0:np)
    real(R8), allocatable              :: nodes_upper(:,:) ! degenerate-strip node coords j=1 row (2, 0:np)
    integer                            :: unit
    logical                            :: header_written
  contains
    private
    procedure, pass(self) :: Update => Update_LineProbe
  end type obj_lineprobe_t

  integer :: nlineprobes = 0
  type(obj_lineprobe_t), allocatable, target :: lineprobe(:)

contains

  subroutine Setup_LineProbes ( domain )
    use strings,            only: parse
    use Q2D_Mod_MPI,        only: is_local_block
    use Q2D_Config_Types_m, only: obj_io, obj_io_lineprobes
    implicit none
    type(Q2D_domain_type), intent(in), target :: domain
    ! Local
    integer             :: i, j, b, dir, idx, ios
    character(len=llen) :: name, filepath
    character(len=clen) :: string(64)
    logical             :: restart, file_exists

    ! Registry-driven (peculiarità #5): the per-line-probe config is loaded into
    ! obj_io_lineprobes by Register_LineProbes + Load_Ini, and the section-level
    ! restart/newrun into obj_io%lineprobe_restart by Scan_Ini — no fini%get here.
    ! Absent input.ini -> obj_io_lineprobes never allocated (Read_Inifile bootstrap
    ! guard returns early), so the allocated-check reproduces the legacy no-op.
    if (.not. allocated(obj_io_lineprobes)) then
      nlineprobes = 0
      return
    end if
    nlineprobes = size(obj_io_lineprobes)
    if (nlineprobes == 0) return

    restart = obj_io%lineprobe_restart

    if (allocated(lineprobe)) deallocate(lineprobe)
    allocate(lineprobe(nlineprobes))

    do i = 1, nlineprobes
      name = obj_io_lineprobes(i)%name
      lineprobe(i)%name = name
      lineprobe(i)%unit = -1
      lineprobe(i)%header_written = .false.

      ! Variable list — parse the registered string (mirror IO_Probes::Setup_Probes)
      string = ''
      lineprobe(i)%nvar = 0
      call parse(obj_io_lineprobes(i)%varnames, ' ', string)
      ! The registry default 'none' is the no-variables sentinel: a coordinate-only
      ! probe (legal in legacy, where count_values on an absent key returned 0).
      ! Treat it as nvar=0 so Assign_Variables/write loops are zero-trip; a genuine
      ! misspelling ('rh0', or 'none' mixed with real vars) keeps nvar>=1 and is
      ! still caught by the unrecognized-variable error stop in Assign_Variables.
      if (.not. (size(string) >= 1 .and. trim(adjustl(obj_io_lineprobes(i)%varnames)) == 'none')) then
        do j = 1, size(string)
          if (string(j) /= '') lineprobe(i)%nvar = lineprobe(i)%nvar + 1
        end do
      end if
      allocate(lineprobe(i)%names(1:lineprobe(i)%nvar))
      do j = 1, lineprobe(i)%nvar
        lineprobe(i)%names(j) = trim(string(j))
      end do

      ! Frequency (registered; defaults 1e30 / 1e9 ≈ legacy huge = "write only if set")
      lineprobe(i)%dtime = obj_io_lineprobes(i)%dtime
      lineprobe(i)%diter = obj_io_lineprobes(i)%diter
      lineprobe(i)%ntime = 1

      ! Location (mandatory) — the '0' sentinel default reproduces the legacy
      ! `error stop` when index-position / fixed-index are omitted.
      lineprobe(i)%ilocation   = obj_io_lineprobes(i)%iloc
      lineprobe(i)%fixed_index = obj_io_lineprobes(i)%fixed_index
      if (lineprobe(i)%ilocation(1) < 1 .or. lineprobe(i)%ilocation(2) < 1) &
        error stop "[Q2D_IO_LineProbes] Missing index-position for lineprobe"
      if (lineprobe(i)%fixed_index < 1) &
        error stop "[Q2D_IO_LineProbes] Missing fixed-index for lineprobe"

      b   = lineprobe(i)%ilocation(1)
      dir = lineprobe(i)%ilocation(2)
      idx = lineprobe(i)%fixed_index

      ! MPI: skip lineprobes whose target block is remote (consistent with the
      ! MOSE IO_Probes::Setup_Probes pattern). domain%blk(b)%P is deallocated for
      ! non-local blocks by deallocate_remote_computation_data.
      if (.not. is_local_block(b)) cycle

      if (dir == 1) then
        ! i fixed, extract along j
        lineprobe(i)%npoints = domain%blk(b)%dim(2)
        lineprobe(i)%P => domain%blk(b)%P(:, idx, 1:domain%blk(b)%dim(2))
      else
        ! j fixed, extract along i
        lineprobe(i)%npoints = domain%blk(b)%dim(1)
        lineprobe(i)%P => domain%blk(b)%P(:, 1:domain%blk(b)%dim(1), idx)
      end if

      allocate(lineprobe(i)%T(lineprobe(i)%npoints))
      allocate(lineprobe(i)%M(lineprobe(i)%npoints))

      call Assign_Variables(lineprobe(i))
      call Setup_Nodes(lineprobe(i), domain)

      if (restart) then
        filepath = 'OUTPUT/'//trim(name)//'.tec'
        inquire(file=trim(filepath), exist=file_exists)
        if (file_exists) then
          open(newunit=lineprobe(i)%unit, file=trim(filepath), status='OLD', &
               form='formatted', position='append', iostat=ios)
          if (ios == 0) lineprobe(i)%header_written = .true.
        end if
      end if

    end do

  end subroutine Setup_LineProbes


  subroutine Setup_Nodes ( lp, domain )
    implicit none
    type(obj_lineprobe_t),  intent(inout) :: lp
    type(Q2D_domain_type),  intent(in)    :: domain
    ! Local
    integer :: b, dir, idx, ip, np_

    b   = lp%ilocation(1)
    dir = lp%ilocation(2)
    idx = lp%fixed_index
    np_ = lp%npoints

    allocate(lp%nodes_lower(1:2, 0:np_))
    allocate(lp%nodes_upper(1:2, 0:np_))

    if (dir == 1) then
      ! i fixed, extract along j
      do ip = 0, np_
        lp%nodes_lower(1:2, ip) = domain%blk(b)%node(idx-1, ip)%c(1:2)
        lp%nodes_upper(1:2, ip) = domain%blk(b)%node(idx,   ip)%c(1:2)
      end do
    else
      ! j fixed, extract along i
      do ip = 0, np_
        lp%nodes_lower(1:2, ip) = domain%blk(b)%node(ip, idx-1)%c(1:2)
        lp%nodes_upper(1:2, ip) = domain%blk(b)%node(ip, idx  )%c(1:2)
      end do
    end if

  end subroutine Setup_Nodes


  subroutine Assign_Variables ( lp )
    use Q2D_Global_m, only: nsc, nu, nv, np
    use IR_Precision, only: str
    implicit none
    type(obj_lineprobe_t), intent(inout), target :: lp
    ! Local
    integer :: v, s

    allocate(lp%variables(1:lp%nvar))

    do v = 1, lp%nvar
      do s = 1, nsc
        if (lp%names(v) == 'rho('//trim(str(.true.,s))//')') lp%variables(v)%p => lp%P(s, :)
      end do
      ! No "w" component.
      if (lp%names(v) == 'u') lp%variables(v)%p => lp%P(nu, :)
      if (lp%names(v) == 'v') lp%variables(v)%p => lp%P(nv, :)
      if (lp%names(v) == 'p') lp%variables(v)%p => lp%P(np, :)
      if (lp%names(v) == 'T') lp%variables(v)%p => lp%T
      if (lp%names(v) == 'M') lp%variables(v)%p => lp%M
      if (.not. associated(lp%variables(v)%p)) then
        write(*,*) '[Q2D_IO_LineProbes] Unrecognized variable: ', trim(lp%names(v))
        error stop "[Q2D_IO_LineProbes] Unknown lineprobe variable"
      end if
    end do

  end subroutine Assign_Variables


  subroutine Update_LineProbe ( self )
    use Q2D_Global_m,            only: nsc, nu, nv, np
    use FLINT_Lib_Thermodynamic, only: f_Rtot, f_ss, EOS
    implicit none
    class(obj_lineprobe_t), intent(inout) :: self
    ! Local
    integer  :: ip
    real(R8) :: Rtot, speed_of_sound, rho_tot

    do ip = 1, self%npoints
      rho_tot = sum(self%P(1:nsc, ip))
      Rtot = f_Rtot(self%P(1:nsc, ip))
      speed_of_sound = f_ss(self%P(1:nsc, ip), self%P(np, ip), rho_tot, Rtot)
      self%T(ip) = EOS(p=self%P(np, ip), rho=rho_tot, R=Rtot)
      ! Mach 2D (u, v).
      self%M(ip) = norm2(self%P(nu:nv, ip)) / speed_of_sound
    end do
  end subroutine Update_LineProbe


  subroutine Write_LineProbes_Data ( iter, time )
    use Q2D_Mod_MPI, only: is_local_block
    implicit none
    integer,  intent(in) :: iter
    real(R8), intent(in) :: time
    ! Local
    integer :: lp
    logical :: should_write

    if (nlineprobes == 0) return

    do lp = 1, nlineprobes
      ! MPI: skip lineprobes whose target block is remote (mirror the MOSE
      ! IO_Probes::Write_Probes_Data pattern).
      if (.not. is_local_block(lineprobe(lp)%ilocation(1))) cycle
      should_write = .false.
      if (mod(iter, lineprobe(lp)%diter) == 0) then
        should_write = .true.
      else if (time >= lineprobe(lp)%dtime * lineprobe(lp)%ntime) then
        should_write = .true.
        lineprobe(lp)%ntime = lineprobe(lp)%ntime + 1
      end if
      if (.not. should_write) cycle
      call Write_LineProbe_Tecplot(lineprobe(lp), time)
    end do
  end subroutine Write_LineProbes_Data


  subroutine Write_LineProbe_Tecplot ( lp, time )
    use IR_Precision, only: str, FR_P
    implicit none
    type(obj_lineprobe_t), intent(inout) :: lp
    real(R8), intent(in) :: time
    ! Local
    integer             :: i, v, Ni, NvarTot
    character(len=llen) :: filepath
    character(len=500)  :: varloc_str, zoneheader
    character(len=hlen) :: varheader

    call lp%Update()

    Ni = lp%npoints
    NvarTot = 2 + lp%nvar  ! x, y + variables

    if (.not. lp%header_written) then
      filepath = 'OUTPUT/'//trim(lp%name)//'.tec'
      open(newunit=lp%unit, file=trim(filepath), status='REPLACE', form='formatted')
      varheader = ' VARIABLES ="x" "y"'
      do v = 1, lp%nvar
        varheader = trim(varheader)//' "'//trim(lp%names(v))//'"'
      end do
      write(lp%unit, '(A)') trim(varheader)
      lp%header_written = .true.
    end if

    if (lp%nvar == 1) then
      varloc_str = ', VARLOCATION=([1-2]=NODAL,[3]=CELLCENTERED)'
    else
      varloc_str = ', VARLOCATION=([1-2]=NODAL,[3-'//trim(str(.true.,NvarTot))//']=CELLCENTERED)'
    end if

    zoneheader = ' ZONE  T = '//trim(lp%name)// &
                 ', I='//trim(str(no_sign=.true., n=Ni+1))// &
                 ', J=2, K=1'// &
                 ', DATAPACKING=BLOCK'//adjustl(trim(varloc_str))
    zoneheader = trim(zoneheader)//', SOLUTIONTIME='//trim(str(no_sign=.true., n=time))
    write(lp%unit, '(A)') trim(zoneheader)

    ! x coordinates: lower row (j=0) then upper row (j=1) — Tecplot BLOCK packing.
    write(lp%unit, FR_P) (lp%nodes_lower(1, i), i=0, Ni), (lp%nodes_upper(1, i), i=0, Ni)
    write(lp%unit, FR_P) (lp%nodes_lower(2, i), i=0, Ni), (lp%nodes_upper(2, i), i=0, Ni)

    do v = 1, lp%nvar
      write(lp%unit, FR_P) (lp%variables(v)%p(i), i=1, Ni)
    end do

    flush(lp%unit)

  end subroutine Write_LineProbe_Tecplot


  subroutine Close_LineProbes ()
    implicit none
    integer :: i
    logical :: is_open

    if (.not. allocated(lineprobe)) return

    do i = 1, nlineprobes
      if (lineprobe(i)%unit > 0) then
        inquire(unit=lineprobe(i)%unit, opened=is_open)
        if (is_open) close(lineprobe(i)%unit)
      end if
      if (allocated(lineprobe(i)%nodes_lower)) deallocate(lineprobe(i)%nodes_lower)
      if (allocated(lineprobe(i)%nodes_upper)) deallocate(lineprobe(i)%nodes_upper)
    end do

  end subroutine Close_LineProbes

end module Q2D_IO_LineProbes
