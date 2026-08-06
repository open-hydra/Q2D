module Q2D_Read_Ini

  implicit none
  private
  public :: Read_Inifile, Read_Inifile_Runtime

contains

  subroutine Read_Inifile()
    use Q2D_Global_m,       only: nsc
    use Finer,              only: file_ini
    use Q2D_Read_Sim_Param, only: Register_Sim_Param
    use Q2D_Read_IO,        only: Register_IO_Fields, Register_Probes, Register_LineProbes
    use Q2D_Read_Numerics,  only: Register_Numerics
    use Q2D_Read_Physics,   only: Register_Physics
    use Q2D_Read_Sources,   only: Register_Sources
    use Q2D_Backend_INI,    only: Load_Ini, Scan_Ini
    use Q2D_Config_Types_m, only: obj_chemistry, obj_sim_param
    use Q2D_Parameters_m,   only: llen
    use Q2D_Mod_MPI,        only: mpi_is_root
    implicit none
    ! Local variables
    type(file_ini) :: fini
    integer :: nprobes, nlineprobes, nmgl
    character(len=llen), allocatable :: probes_name(:)
    character(len=llen), allocatable :: lineprobes_names(:)
    logical :: have_ini

    ! Bootstrap guard: without input.ini the setup is a no-op consistent with the
    ! main loop (obj_sim_param%TODO=3 → do-while bypass).
    inquire(file='input.ini', exist=have_ini)
    if (.not. have_ini) then
      if (mpi_is_root) write(*,'(A)') &
        ' [Q2D] input.ini not found — driver wired, body is no-op until INPUT/ valid.'
      obj_sim_param%TODO = 3
      return
    end if

    allocate(obj_chemistry%RT(1:nsc+1), obj_chemistry%AT(1:nsc+1))

    ! Load input.ini
    call fini%load(filename='input.ini')

    ! Scan input.ini for the a-priori-unknown number of probes and of multigrid
    ! levels (plus the stationary faces of obj_rot)
    call Scan_Ini(fini, nprobes, probes_name, nlineprobes, lineprobes_names, nmgl)

    ! Build the registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(nprobes, probes_name)
    call Register_LineProbes(nlineprobes, lineprobes_names)
    call Register_Numerics(nmgl)
    call Register_Physics()
    ! Register the geometric sources, after the standard singletons so that
    ! Register_Sources finds the base types already active.
    call Register_Sources()

    ! Registry built: load the values from the ini file
    call Load_Ini(fini)

  end subroutine Read_Inifile


  subroutine Read_Inifile_Runtime()
    use Finer,           only: file_ini
    use Q2D_Backend_INI, only: Load_Ini, Scan_NoChem
    use Q2D_Registry_m,  only: Validate_Registry
    implicit none
    ! Local variables
    type(file_ini) :: fini
    character(len=1024) :: out

    ! Load input.ini
    call fini%load(filename='input.ini')

    ! Registry already built: load the values from the ini file
    call Load_Ini(fini)

    ! List of blocks excluded from chemistry: read without limits (it is not a
    ! Registry parameter), so that the runtime/coupled path also acquires it
    ! (the standalone path reads it via Scan_Ini).
    call Scan_NoChem(fini)

    ! Registry validation
    out = Validate_Registry()

  end subroutine Read_Inifile_Runtime


end module Q2D_Read_Ini
