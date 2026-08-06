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
    character(len=llen), allocatable, intent(out) :: probes_name(:)
    character(len=llen), allocatable, intent(out) :: lineprobes_names(:)
    ! Local variables
    character(len=llen) :: wholestring
    integer :: error, i

    ! Read the number of multigrid levels. ATLAS may set
    ! `obj_sim_param%HYDRA_MG = .true.` to read the `[HYDRA-Multigrid]` section
    ! instead of `[<codename>-Multigrid]`. The default `.false.` corresponds to
    ! the INI parsing of the standalone Q2D run.
    nmgl = 1
    if (.not. obj_sim_param%HYDRA_MG) then
      call fini%get(section_name=trim(codename)//'-Multigrid', option_name='levels', val=nmgl, error=error)
    else
      call fini%get(section_name='HYDRA-Multigrid', option_name='levels', val=nmgl, error=error)
    end if

    ! Probe count
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

    ! Line-probe count, with the same logic as the probe count.
    ! The configuration of the individual line-probes goes through the registry
    ! (Register_LineProbes); here only the count, the section names and the
    ! section-level restart are read.
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

    ! Line-probe restart/newrun (section level): the `restart` key takes
    ! priority; in its absence, negated `newrun` is used; default .false.
    obj_io%lineprobe_restart = .false.
    call fini%get(section_name=trim(codename)//'-LineProbes', option_name='restart', val=obj_io%lineprobe_restart, error=error)
    if (error /= 0) then
      call fini%get(section_name=trim(codename)//'-LineProbes', option_name='newrun', val=obj_io%lineprobe_restart, error=error)
      if (error == 0) obj_io%lineprobe_restart = .not. obj_io%lineprobe_restart
    end if

    ! Count of the stationary faces for the rotating frame
    obj_rot%n_stationary = 0
    i = 0
    do
      call fini%get(section_name=trim(codename)//'-Rotating-Frame',option_name='stationary-face-'//trim(str(.true.,i+1)), val=wholestring, error=error)
      if (error /= 0) exit
      obj_rot%n_stationary = obj_rot%n_stationary + 1
    enddo
    i = i + 1
    if (obj_rot%n_stationary > 0) allocate(obj_rot%stationary_face_str(obj_rot%n_stationary))

    ! List of blocks excluded from chemistry: read without constraints on length,
    ! count or number of digits. The standalone path reaches it from here.
    call Scan_NoChem(fini)

  end subroutine Scan_Ini


  ! Reads from the FiNeR object the exclude-blocks integer list of the
  ! [<codename>-Chemistry] section, with no limit on the number of blocks (count_values
  ! sizes the token array; no fixed cap, unlike a Registry string parameter truncated
  ! to llen). Tokens are read as a character array (FiNeR's array get supports
  ! character(*), not type(string)) and each non-empty token is validated with
  ! StringiFor's is_integer before conversion: FiNeR's integer get silently returns
  ! spurious values on a non-integer token (to_number skips the read, errd stays 0),
  ! so this guard prevents unreported wrong reads.
  ! count_values is an upper bound (it counts raw separators, so multiple blanks
  ! leave empty trailing slots); empty tokens are harmless and are skipped, while a
  ! non-empty non-integer token makes the whole list be rejected. The per-token width
  ! clen (16) covers any valid integer block id; the number of blocks is unlimited.
  ! Called by BOTH INI reading paths (standalone + runtime/coupled).
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
      ! Pre-initialization to the empty string: count_values overestimates in the
      ! presence of repeated blanks and FiNeR's get writes only the real tokens
      ! (current behavior of get_a_option in finer_option_t: it loops over
      ! size(separated tokens), not over size(val) — not a documented API contract),
      ! leaving the excess slots at the pre-initialized value ''. Empty slots are
      ! therefore the artifact of multiple blanks (skipped), distinct from a truly
      ! non-integer token (rejected). If FiNeR filled the excess slots, the spurious
      ! content would be non-empty and would trigger is_integer below -> safe
      ! rejection, never silent.
      allocate(tokens(nochem_count))
      tokens = ''
      call fini%get(section_name=trim(codename)//'-Chemistry', option_name='exclude-blocks', &
                    val=tokens, error=error)
      allocate(nochem_tmp(nochem_count), source=0_I4)
      reject_list = (error /= 0)
      do v = 1, nochem_count
        if (len_trim(tokens(v)) == 0) cycle   ! empty slot from the upper bound (collapsed blanks)
        tok = trim(adjustl(tokens(v)))
        if (tok%is_integer()) then
          nochem_tmp(v) = tok%to_number(kind=1_I4)
        else
          reject_list = .true.                ! non-empty non-integer token -> list rejected
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
