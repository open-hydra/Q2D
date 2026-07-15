module Q2D_Wrap_Postprocess
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Q2D_postprocess

  integer,        private :: id_stampa = 0
  character(81),  private :: RANS_shell_format     = "('Q2D | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(69),  private :: URANS_shell_format    = "('Q2D | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"
  character(96),  private :: RANS_shell_format_MG  = "('Q2D Grid Level', i2, ' | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(84),  private :: URANS_shell_format_MG = "('Q2D Grid Level', i2, ' | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"

contains

  subroutine Q2D_postprocess ( simulation )
    use Q2D_Advanced_Types_m, only: Q2D_simulation_type
    use Q2D_Config_Types_m,   only: obj_sim_param, obj_io, obj_time_scheme, obj_multigrid, obj_rot
    use Q2D_Global_m,         only: Q2D_phase_prefix
    use Q2D_Parameters_m,     only: llen
    use Q2D_IO_Wall,          only: IOwall, Write_Wall_Solution
    use Q2D_IO_Solution,      only: Write_solution
    use Q2D_Mod_Diagnostic,   only: Write_Diagnostic
    use Q2D_IO_Probes,        only: Write_Probes_Data
    use Q2D_IO_LineProbes,    only: Write_LineProbes_Data, Close_LineProbes
    use Q2D_Mod_Multigrid,    only: Prolongation
    use Q2D_Read_Ini,         only: Read_Inifile_Runtime
    use Q2D_Mod_MPI,          only: mpi_is_root
    use Q2D_Lib_RotatingFrame, only: RF_Torque_Walls
    use IR_Precision,         only: str
    implicit none
    type(Q2D_simulation_type), intent(inout) :: simulation
    ! Locali
    character(len=llen) :: solfile, dgsfile, wallfile, mgsol, mgwall
    real(R8) :: sim_time
    integer  :: m, level

    level = obj_multigrid%MG_level

    ! Aggiorna il solutiontime di IOfield/IOwall
    if (obj_time_scheme%time_accurate) then
      simulation%IOfield(1)%solutiontime = simulation%domain(1)%time
      IOwall%solutiontime = simulation%domain(1)%time
    else
      simulation%IOfield(1)%solutiontime = -real(obj_sim_param%iter_general, R8)
      IOwall%solutiontime = -real(obj_sim_param%iter_general, R8)
    end if

    if (level == 1) then
      solfile  = '/'//trim(Q2D_phase_prefix)//'field'
      wallfile = '/'//trim(Q2D_phase_prefix)//'wall'
      dgsfile  = '/'//trim(Q2D_phase_prefix)//'diagnostic'
    else
      solfile  = '/'//trim(Q2D_phase_prefix)//'field-level'//trim(str(.true.,level))
      ! wallfile/dgsfile vanno definiti a OGNI livello: il ramo TODO==3 più sotto
      ! li scrive incondizionatamente e nomi non inizializzati farebbero fallire
      ! i livelli multigrid coarse con "Cannot open file 'OUTPUT/'".
      wallfile = '/'//trim(Q2D_phase_prefix)//'wall-level'//trim(str(.true.,level))
      dgsfile  = '/'//trim(Q2D_phase_prefix)//'diagnostic-level'//trim(str(.true.,level))
      mgsol    = '/'//trim(Q2D_phase_prefix)//'field-prolongated'
      mgwall   = '/'//trim(Q2D_phase_prefix)//'wall-prolongated'
    end if

    ! ===========================================================================
    ! FINE SIMULAZIONE (TODO == 3)
    ! ===========================================================================
    if (obj_sim_param%TODO == 3) then

      if (mpi_is_root) then
        if (mod(simulation%domain(1)%iter, obj_io%res_diter) == 0) then
          write(obj_io%unitRES, trim(obj_io%unitRES_format)) &
            obj_sim_param%iter_general, simulation%domain(1)%time, obj_sim_param%residuotot
        end if
      end if

      ! Coppia di pressione sulle pareti rotanti.
      if (obj_rot%enabled .and. level == 1 .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0) &
        call RF_Torque_Walls( simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time )

      if (mpi_is_root) then
        if (obj_time_scheme%time_accurate) then
          if (obj_multigrid%MGL > 1) then
            write(*, URANS_shell_format_MG) level, obj_sim_param%iter_from_call, &
              simulation%domain(1)%time, simulation%domain(1)%dtglobal
          else
            write(*, URANS_shell_format) obj_sim_param%iter_from_call, &
              simulation%domain(1)%time, simulation%domain(1)%dtglobal
          end if
        else
          if (obj_multigrid%MGL > 1) then
            write(*, RANS_shell_format_MG) level, obj_sim_param%iter_from_call, &
              obj_sim_param%iter_general, obj_sim_param%residuotot(1)
          else
            write(*, RANS_shell_format) obj_sim_param%iter_from_call, &
              obj_sim_param%iter_general, obj_sim_param%residuotot(1)
          end if
        end if

        call cpu_time(obj_sim_param%cputime(2))
        sim_time = (obj_sim_param%cputime(2) - obj_sim_param%cputime(1)) / max(obj_sim_param%nthreads, 1)
        write(*,*)
        write(*,*) '  Time of operation was', sim_time/60.0_R8, 'min'
      end if

      if (associated(Write_solution)) then
        call Write_solution(simulation%domain(1), simulation%IOfield(1), solfile)
      end if
      if (obj_io%write_wall) call Write_Wall_Solution(simulation%domain(1), wallfile)
      if (.not. obj_time_scheme%time_accurate) &
        call Write_Diagnostic(simulation%domain(1), simulation%IOfield(1), dgsfile)

      ! Chiusura ordinata dei file delle line-probes.
      call Close_LineProbes()


    ! ===========================================================================
    ! SOLUZIONE INTERMEDIA (TODO == 2)
    ! ===========================================================================
    else if (obj_sim_param%TODO == 2) then

      if (level == 1) then
        id_stampa = id_stampa + 1
        if (.not. obj_io%sol_overwrite) solfile  = trim(solfile) //trim(str(.true.,id_stampa))
        if (.not. obj_io%sol_overwrite) wallfile = trim(wallfile)//trim(str(.true.,id_stampa))
        if (.not. obj_io%sol_overwrite) dgsfile  = trim(dgsfile) //trim(str(.true.,id_stampa))
      end if

      if (mpi_is_root) then
        if (mod(simulation%domain(level)%iter, obj_io%res_diter) == 0) then
          write(obj_io%unitRES, trim(obj_io%unitRES_format)) &
            obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
        end if
      end if

      if (level == 1) then
        call Write_Probes_Data(obj_sim_param%iter_general, simulation%domain(1)%time)
        call Write_LineProbes_Data(obj_sim_param%iter_general, simulation%domain(1)%time)
      end if

      ! Coppia di pressione sulle pareti rotanti.
      if (level == 1 .and. obj_rot%enabled .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0) &
        call RF_Torque_Walls( simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time )

      if (mpi_is_root) then
        if (mod(simulation%domain(level)%iter, obj_io%shell_diter) == 0) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1) then
              write(*, URANS_shell_format_MG) level, obj_sim_param%iter_from_call, &
                simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write(*, URANS_shell_format) obj_sim_param%iter_from_call, &
                simulation%domain(1)%time, simulation%domain(1)%dtglobal
            end if
          else
            if (obj_multigrid%MGL > 1) then
              write(*, RANS_shell_format_MG) level, obj_sim_param%iter_from_call, &
                obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            else
              write(*, RANS_shell_format) obj_sim_param%iter_from_call, &
                obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            end if
          end if
        end if
      end if

      ! Output della soluzione: a cadenza di iterazioni oppure di tempo fisico.
      if (mod(simulation%domain(level)%iter, obj_io%sol_diter) == 0) then
        if (mpi_is_root) write(*,*) ' ... writing iter-based solution'
        if (associated(Write_solution)) &
          call Write_solution(simulation%domain(level), simulation%IOfield(level), solfile)
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation(Fine=simulation%domain(m-1), Coarse=simulation%domain(m))
          end do
          if (associated(Write_solution)) &
            call Write_solution(simulation%domain(1), simulation%IOfield(1), mgsol)
          if (obj_io%write_wall) call Write_Wall_Solution(simulation%domain(1), mgwall)
        end if
        if (level == 1) then
          if (obj_io%write_wall) call Write_Wall_Solution(simulation%domain(level), wallfile)
          if (.not. obj_time_scheme%time_accurate) &
            call Write_Diagnostic(simulation%domain(level), simulation%IOfield(level), dgsfile)
        end if
      else if (simulation%domain(level)%time >= obj_sim_param%time_from_call + obj_io%sol_dtime) then
        obj_sim_param%time_from_call = simulation%domain(level)%time
        if (mpi_is_root) write(*,*) ' ... writing time-based solution'
        if (associated(Write_solution)) &
          call Write_solution(simulation%domain(level), simulation%IOfield(level), solfile)
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation(Fine=simulation%domain(m-1), Coarse=simulation%domain(m))
          end do
          if (associated(Write_solution)) &
            call Write_solution(simulation%domain(1), simulation%IOfield(1), mgsol)
          if (obj_io%write_wall) call Write_Wall_Solution(simulation%domain(1), mgwall)
        end if
        if (level == 1) then
          if (obj_io%write_wall) call Write_Wall_Solution(simulation%domain(level), wallfile)
          if (.not. obj_time_scheme%time_accurate) &
            call Write_Diagnostic(simulation%domain(level), simulation%IOfield(level), dgsfile)
        end if
      end if

      if (obj_multigrid%change_MG) then
        if (mpi_is_root) write(*,*) ' ... writing solution of grid-level ', level
        if (associated(Write_solution)) &
          call Write_solution(simulation%domain(level), simulation%IOfield(level), solfile)
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation(Fine=simulation%domain(m-1), Coarse=simulation%domain(m))
          end do
          if (associated(Write_solution)) &
            call Write_solution(simulation%domain(1), simulation%IOfield(1), mgsol)
          if (obj_io%write_wall) call Write_Wall_Solution(simulation%domain(1), mgwall)
        end if
        obj_multigrid%MG_level = obj_multigrid%MG_level - 1
        obj_sim_param%iter_from_call = 0
        obj_sim_param%time_from_call = 0.0_R8
      end if

      ! Scrittura per l'accoppiamento HYDRA. L'orchestratore top-level ATLAS può
      ! impostare `obj_sim_param%HYDRA_postprocess = .true.`; il default `.false.`
      ! preserva il comportamento del build standalone di Q2D (blocco no-op).
      if ( obj_sim_param%HYDRA_postprocess ) then
        if (mpi_is_root) write(*,*) ' ... writing HYDRA-coupling solution'
        if (associated(Write_solution)) &
          call Write_solution( simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          end do
          if (associated(Write_solution)) &
            call Write_solution( simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io%write_wall) call Write_Wall_Solution( simulation%domain(1), mgwall )
        end if
      end if


    ! ===========================================================================
    ! OPERAZIONI AUSILIARIE A OGNI ITERAZIONE (TODO == 1)
    ! ===========================================================================
    else if (obj_sim_param%TODO == 1) then

      if (mpi_is_root) then
        if (mod(simulation%domain(level)%iter, obj_io%res_diter) == 0) then
          write(obj_io%unitRES, trim(obj_io%unitRES_format)) &
            obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
        end if
      end if

      if (level == 1) then
        call Write_Probes_Data(obj_sim_param%iter_general, simulation%domain(1)%time)
        call Write_LineProbes_Data(obj_sim_param%iter_general, simulation%domain(1)%time)
      end if

      ! Coppia di pressione sulle pareti rotanti.
      if (level == 1 .and. obj_rot%enabled .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0) &
        call RF_Torque_Walls( simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time )

      if (mpi_is_root) then
        if (mod(simulation%domain(level)%iter, obj_io%shell_diter) == 0) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1) then
              write(*, URANS_shell_format_MG) level, obj_sim_param%iter_from_call, &
                simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write(*, URANS_shell_format) obj_sim_param%iter_from_call, &
                simulation%domain(1)%time, simulation%domain(1)%dtglobal
            end if
          else
            if (obj_multigrid%MGL > 1) then
              write(*, RANS_shell_format_MG) level, obj_sim_param%iter_from_call, &
                obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            else
              write(*, RANS_shell_format) obj_sim_param%iter_from_call, &
                obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            end if
          end if
        end if
      end if

    end if

    ! Rilettura periodica del file di input a runtime.
    if (mod(simulation%domain(level)%iter, obj_io%ini_diter) == 0) then
      call Read_Inifile_Runtime()
    end if

  end subroutine Q2D_postprocess

end module Q2D_Wrap_Postprocess
