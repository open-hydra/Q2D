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
    ! Local
    type(file_ini) :: fini
    integer :: nprobes, nlineprobes, nmgl
    character(len=16),   allocatable :: probes_name(:)
    character(len=llen), allocatable :: lineprobes_names(:)
    logical :: have_ini

    ! Bootstrap guard: without input.ini the setup is a no-op consistent with the
    ! main loop (obj_sim_param%TODO=3 → bypass do-while).
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

    ! Scan input.ini for unknown number of probes and multigrid levels (also obj_rot stationary faces)
    call Scan_Ini(fini, nprobes, probes_name, nlineprobes, lineprobes_names, nmgl)

    ! Build registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(nprobes, probes_name)
    call Register_LineProbes(nlineprobes, lineprobes_names)   ! Q2D-extra (peculiarità #5)
    call Register_Numerics(nmgl)
    call Register_Physics()
    ! Q2D-specific: register the geometric sources. After the standard singletons
    ! so that Register_Sources sees the base types already active.
    call Register_Sources()

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

  end subroutine Read_Inifile


  subroutine Read_Inifile_Runtime()
    use Finer,           only: file_ini
    use Q2D_Backend_INI, only: Load_Ini, Scan_NoChem
    use Q2D_Registry_m,  only: Validate_Registry
    implicit none
    ! Local
    type(file_ini) :: fini
    character(len=1024) :: out

    ! Load input.ini
    call fini%load(filename='input.ini')

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

    ! No-chemistry block list: unbounded read (not a Registry param), so the
    ! coupled/runtime path picks it up too (the standalone path reads it via Scan_Ini).
    call Scan_NoChem(fini)

    ! Validate registry
    out = Validate_Registry()

  end subroutine Read_Inifile_Runtime


end module Q2D_Read_Ini
