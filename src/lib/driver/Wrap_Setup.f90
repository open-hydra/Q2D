module Q2D_Wrap_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Q2D_setup

contains

  subroutine Q2D_setup(simulation)
    use Q2D_Advanced_Types_m,     only: Q2D_simulation_type
    use Q2D_Config_Types_m,       only: obj_sim_param, obj_io, obj_time_scheme, &
                                        obj_multigrid, obj_source
    use Q2D_Global_m,             only: Q2D_phase_prefix, nres
    use Q2D_Read_Ini,             only: Read_Inifile
    use Q2D_Load_ThermoTransport, only: Load_ThermoTransport
    use Q2D_Load_Chemistry,       only: Load_Chemistry
    use Q2D_Assign_Setup,         only: Assign_Setup
    use Q2D_IO_Solution,          only: Read_IC, Setup_Output_Solution
    use Q2D_Mod_Allocate_Data,    only: Setup_Data_Structure, &
                                        deallocate_remote_computation_data
    use Q2D_Mod_Multigrid,        only: Setup_Multigrid, Restriction
    use Q2D_IO_BC,                only: Setup_BC
    use Q2D_Mod_Metrics,          only: Setup_Metrics
    use Q2D_IO_Probes,            only: Setup_Probes
    use Q2D_IO_LineProbes,        only: Setup_LineProbes
    use Q2D_IO_Wall,              only: Initialize_Wall_File
    use Q2D_Lib_Ghost,            only: Fill_Ghost_Cell
    use Q2D_Mod_MPI,              only: mpi_is_root, partition_blocks
    use Q2D_Mod_GhostExchange,    only: build_ghost_schedule, build_local_bc_index
    use Q2D_Mod_Diagnostic,       only: Print_Header, Check_Tables, Check_Input, &
                                        Check_IC, Check_BC, Print_Warnings, &
                                        Stop_Simulation, Print_Shell_Info
    use Q2D_Mod_Sources,          only: Initialize_Sources  ! sorgenti geometriche (inizializzazione a setup)
    use Q2D_Lib_Source_Areavariation, only: Read_Block_Q2D_Area
    implicit none
    type(Q2D_simulation_type), intent(inout) :: simulation
    ! Locali
    integer :: m, ios

    ! ---- Struttura del setup -------------------------------------------------------
    !
    ! Il guard sul file di bootstrap risiede in Q2D_Read_Inifile.
    ! Collocazione delle inizializzazioni delle sorgenti geometriche e delle probe:
    !   - Read_Block_Q2D_Area (area-variation) → lettura UNICA su domain(1) più
    !     sotto (dopo Setup_Data_Structure, prima di Setup_Multigrid); la
    !     restrizione fine→coarse avviene via Coarse_Cross_Section in Setup_Multigrid.
    !   - Precompute_Area_Geometry (area-variation) → Q2D_Lib_Source_Areavariation,
    !     invocata dall'hook in Setup_Metrics.
    !   - Precompute_Axisymmetric_Geometry (assialsimmetria) →
    !     Q2D_Lib_Source_Axisymmetry, invocata dall'hook in Setup_Metrics.
    !   - Setup_LineProbes (line-probes) → Q2D_IO_LineProbes::Setup_LineProbes.

    if (mpi_is_root) call Print_Header()

    call Load_ThermoTransport()
    call Load_Chemistry()
    if (mpi_is_root) call Check_Tables()

    call Read_Inifile()
    if (obj_sim_param%TODO == 3) return   ! guard di bootstrap scattato in Read_Inifile (input.ini assente)
    if (mpi_is_root) call Check_Input()

    ! Aggancia il dispatcher delle sorgenti geometriche (assialsimmetria /
    ! area-variation) dopo il parse dell'INI e prima delle allocazioni.
    call Initialize_Sources()

    allocate(simulation%domain (obj_multigrid%MGL))
    allocate(simulation%IOfield(obj_multigrid%MGL))
    simulation%domain(1)%itermax = obj_sim_param%iter_threshold

    call Assign_Setup()

    call Read_IC(simulation%IOfield(1))
    if (mpi_is_root) call Check_IC()

    call Setup_Data_Structure(simulation%domain(1), simulation%IOfield(1))
    allocate(obj_sim_param%residuotot(nres))

    ! area_variation: Cross_Section letta UNA sola volta sulla mesh fine
    ! (domain(1)); la restrizione fine→coarse avviene in Setup_Multigrid via
    ! Coarse_Cross_Section (sotto-campionamento R1).
    if (obj_source%area_variation) call Read_Block_Q2D_Area(simulation%domain(1))

    if (obj_multigrid%MGL > 1) call Setup_Multigrid(simulation)

    call Setup_BC(simulation%domain)
    if (mpi_is_root) call Check_BC()

    do m = 1, obj_multigrid%MGL
      call partition_blocks(simulation%domain(m)%nb, &
        [(product(simulation%domain(m)%blk(ios)%dim(1:2)), ios=1, simulation%domain(m)%nb)])
      call build_ghost_schedule(simulation%domain(m))
      call build_local_bc_index(simulation%domain(m))
    end do

    do m = 1, obj_multigrid%MGL
      call Setup_Metrics(simulation%domain(m))
    end do

    do m = 1, obj_multigrid%MGL
      call deallocate_remote_computation_data(simulation%domain(m))
    end do

    call Setup_Probes(simulation%domain(1), obj_sim_param%newrun)
    call Setup_LineProbes(simulation%domain(1))
    call Setup_Output_Solution(simulation%IOfield)

    if (mpi_is_root .and. obj_io%write_wall) &
      call Initialize_Wall_File(simulation%domain(1), simulation%IOfield(1)%tec%extension)

    if (mpi_is_root) call Print_Shell_Info(simulation)

    if (mpi_is_root) then
      open(newunit=obj_io%unitRES, &
           file='OUTPUT/'//trim(Q2D_phase_prefix)//'residual-history.dat', &
           status='unknown', form='formatted')
      if (.not. obj_sim_param%newrun) then
        ios = 0
        do while (ios == 0)
          read(obj_io%unitRES, *, iostat=ios)
        end do
        backspace(obj_io%unitRES)
      end if
    end if

    write(obj_io%unitRES_format,'(A11,I0,A7)') '(I8,E20.10,', nres, 'E20.10)'

    select case (trim(obj_time_scheme%solver_type))
      case ('euler', 'RK2', 'RK3')
        obj_multigrid%MG_level = obj_multigrid%MGL
        do m = 1, obj_multigrid%MGL
          simulation%domain(m)%iter = 0
          if (obj_sim_param%newrun) then
            if (obj_time_scheme%time_accurate) then
              simulation%domain(m)%time = 0.0_R8
            else
              simulation%domain(m)%time = -1.0_R8
            end if
            obj_sim_param%iter_general = 0
          else
            if (obj_time_scheme%time_accurate) then
              simulation%domain(m)%time = simulation%IOfield(m)%solutiontime
              obj_sim_param%iter_general = 0
            else
              obj_sim_param%iter_general = int(simulation%IOfield(m)%solutiontime)
              simulation%domain(m)%time = -1.0_R8
            end if
          end if
          obj_sim_param%time_from_call = simulation%domain(m)%time
          obj_sim_param%iter_from_call = 0
        end do
        call Fill_Ghost_Cell(simulation%domain(1))
        if (obj_multigrid%MGL > 1) then
          do m = 2, obj_multigrid%MGL
            call Restriction(Fine=simulation%domain(m-1), Coarse=simulation%domain(m))
            call Fill_Ghost_Cell(simulation%domain(m))
          end do
        end if
    end select

    if (mpi_is_root) call Print_Warnings()
    if (mpi_is_root) call Stop_Simulation()

    call Cpu_Time(obj_sim_param%cputime(1))

  end subroutine Q2D_setup

end module Q2D_Wrap_Setup
