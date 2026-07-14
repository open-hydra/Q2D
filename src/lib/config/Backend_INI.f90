module Q2D_Backend_INI
  implicit none
  private
  public :: Load_Ini, Scan_Ini, Scan_NoChem

contains

  subroutine Load_Ini(fini)
    use Q2D_Registry_m, only: reg
    use Finer, only: file_ini
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i, error

    do i = 1, reg%size
      error = 1

      ! Each registry entry has exactly one associated typed pointer.
      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
      else if (associated(reg%params(i)%value%r)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%r, error=error)
      else if (associated(reg%params(i)%value%l)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%l, error=error)
      else if (associated(reg%params(i)%value%s)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%s, error=error)
      else if (associated(reg%params(i)%value%iarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%iarr, error=error)
      else if (associated(reg%params(i)%value%rarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%rarr, error=error)
      end if

      if (error == 0) reg%params(i)%is_set = .true.
    end do

  end subroutine


  subroutine Scan_Ini(fini, nprobes, probes_name, nlineprobes, lineprobes_names, nmgl)
    use Finer,                only: file_ini
    use Q2D_Config_Types_m,   only: obj_rot, obj_sim_param, obj_io
    use Q2D_Parameters_m,     only: codename, llen
    use IR_Precision,         only: str
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(out) :: nprobes, nlineprobes, nmgl
    character(len=16),   allocatable, intent(out) :: probes_name(:)
    character(len=llen), allocatable, intent(out) :: lineprobes_names(:)
    ! Local
    character(len=llen) :: wholestring
    integer :: error, i

    ! Read the MGL — the HYDRA_MG branch mirrors MOSE. ATLAS may set
    ! `obj_sim_param%HYDRA_MG = .true.` to read the `[HYDRA-Multigrid]` section
    ! instead of `[<codename>-Multigrid]`. Default `.false.` preserves the Q2D
    ! standalone-build INI parsing.
    nmgl = 1
    if (.not. obj_sim_param%HYDRA_MG) then
      call fini%get(section_name=trim(codename)//'-Multigrid', option_name='levels', val=nmgl, error=error)
    else
      call fini%get(section_name='HYDRA-Multigrid', option_name='levels', val=nmgl, error=error)
    end if

    ! Count the probes
    nprobes = 0
    do
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,nprobes+1)), val=wholestring, error=error)
      if (error/=0) exit
      nprobes = nprobes+1
    enddo

    allocate(probes_name(nprobes))
    do i = 1, nprobes
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,i)), val=probes_name(i), error=error)
    enddo

    ! Count the line-probes (Q2D-extra, peculiarità #5; MOSE has none) — mirror the
    ! probe count. The per-line-probe config is registry-driven (Register_LineProbes);
    ! only the count + section names + the section-level restart are scanned here.
    nlineprobes = 0
    do
      call fini%get(section_name=trim(codename)//'-LineProbes', option_name='lineprobe'//trim(str(.true.,nlineprobes+1)), val=wholestring, error=error)
      if (error/=0) exit
      nlineprobes = nlineprobes+1
    enddo

    allocate(lineprobes_names(nlineprobes))
    do i = 1, nlineprobes
      call fini%get(section_name=trim(codename)//'-LineProbes', option_name='lineprobe'//trim(str(.true.,i)), val=lineprobes_names(i), error=error)
    enddo

    ! Line-probe restart/newrun (section-level): the `restart` key takes priority,
    ! else `newrun` inverted, default .false. — verbatim-preserved from the legacy
    ! Q2D_IO_LineProbes::Setup_LineProbes parsing (peculiarità #5).
    obj_io%lineprobe_restart = .false.
    call fini%get(section_name=trim(codename)//'-LineProbes', option_name='restart', val=obj_io%lineprobe_restart, error=error)
    if (error /= 0) then
      call fini%get(section_name=trim(codename)//'-LineProbes', option_name='newrun', val=obj_io%lineprobe_restart, error=error)
      if (error == 0) obj_io%lineprobe_restart = .not. obj_io%lineprobe_restart
    end if

    ! Count the stationary faces for the rotating frame
    obj_rot%n_stationary = 0
    i = 0
    do
      call fini%get(section_name=trim(codename)//'-Rotating-Frame',option_name='stationary-face-'//trim(str(.true.,i+1)), val=wholestring, error=error)
      if (error /= 0) exit
      obj_rot%n_stationary = obj_rot%n_stationary + 1
    enddo
    i = i + 1
    if (obj_rot%n_stationary > 0) allocate(obj_rot%stationary_face_str(obj_rot%n_stationary))

    ! No-chemistry block list: unbounded read (no length/count/digit cap), mirror of
    ! legacy Q2D_IO_chemistry::read_ini. The standalone path reaches it through here.
    call Scan_NoChem(fini)

  end subroutine Scan_Ini


  ! Reads the unbounded [<codename>-Chemistry] exclude-blocks integer list from the
  ! FiNeR object. The block COUNT is unlimited (count_values sizes the token array;
  ! no fixed cap, unlike a Registry string param capped at llen). Tokens are read as
  ! a character array (FiNeR's array get supports character(*), not type(string)) and
  ! each non-empty token is validated with StringiFor's is_integer before conversion:
  ! FiNeR's integer get silently returns garbage on a non-integer token (to_number
  ! skips the read, errd stays 0), so this guard is stricter than the legacy
  ! fini%get(integer-array) path it mirrors. count_values is an upper bound (it counts
  ! raw separators, so multiple blanks leave empty trailing slots); empty tokens are
  ! benign and skipped, while a non-empty non-integer token rejects the whole list.
  ! Per-token width clen (16) covers any valid integer block id; the number of blocks
  ! is unbounded. Called from BOTH INI read paths (standalone + runtime/coupled).
  subroutine Scan_NoChem(fini)
    use Finer,              only: file_ini
    use StringiFor,         only: string
    use Q2D_Config_Types_m, only: obj_chemistry
    use Q2D_Parameters_m,   only: codename, clen
    use iso_fortran_env,    only: I4 => int32
    implicit none
    type(file_ini), intent(in) :: fini
    character(len=clen), allocatable :: tokens(:)
    integer(I4), allocatable :: nochem_tmp(:)
    type(string) :: tok
    integer(I4) :: nochem_count, v
    integer :: error
    logical :: reject_list

    if (allocated(obj_chemistry%no_chem_list)) deallocate(obj_chemistry%no_chem_list)

    nochem_count = 0
    if (fini%index(section_name=trim(codename)//'-Chemistry', option_name='exclude-blocks') > 0) &
      nochem_count = fini%count_values(section_name=trim(codename)//'-Chemistry', option_name='exclude-blocks')

    if (nochem_count > 0) then
      ! Pre-init empty: count_values over-counts when blanks repeat, and FiNeR's get
      ! writes only the real tokens (current finer_option_t get_a_option behaviour:
      ! it loops to size(split tokens), not size(val) — not a documented API contract),
      ! leaving the surplus upper-bound slots at the pre-init ''. Empty slots are then
      ! the multiple-blank artefact (skipped), distinct from a genuinely non-integer
      ! token (rejected). If FiNeR ever filled the surplus slots, the garbage would be
      ! non-empty and would trip is_integer below -> safe rejection, never silent.
      allocate(tokens(nochem_count))
      tokens = ''
      call fini%get(section_name=trim(codename)//'-Chemistry', option_name='exclude-blocks', &
                    val=tokens, error=error)
      allocate(nochem_tmp(nochem_count), source=0_I4)
      reject_list = (error /= 0)
      do v = 1, nochem_count
        if (len_trim(tokens(v)) == 0) cycle   ! empty upper-bound slot (collapsed blanks)
        tok = trim(adjustl(tokens(v)))
        if (tok%is_integer()) then
          nochem_tmp(v) = tok%to_number(kind=1_I4)
        else
          reject_list = .true.                ! non-empty non-integer token -> reject list
        end if
      end do
      if (reject_list) then
        obj_chemistry%warning_message = '[WARNING] Non-integer token in exclude-blocks; list ignored.'
        allocate(obj_chemistry%no_chem_list(0))
      else
        obj_chemistry%no_chem_list = pack(nochem_tmp, nochem_tmp >= 1)
      end if
    else
      allocate(obj_chemistry%no_chem_list(0))
    end if

  end subroutine Scan_NoChem

end module Q2D_Backend_INI
