program Q2D_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use Q2D_Procedures_m,      only: Q2D_type
  use Q2D_Advanced_Types_m,  only: Q2D_simulation_type
  use Q2D_Config_Types_m,    only: obj_sim_param
  use Q2D_Mod_MPI,           only: mpi_init_env, mpi_finalize_env, mpi_is_root, mpi_size_
#ifdef USE_MPI
  use Q2D_Mod_GhostExchange, only: cleanup_ghost_schedule
#endif

  implicit none
  type(Q2D_type)            :: Q2D
  type(Q2D_simulation_type) :: simulation

  ! Initialize the MPI environment (no-op if USE_MPI is not defined).
  call mpi_init_env()

#if defined (_OPENMP)
  !$omp parallel
  obj_sim_param%nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  if (mpi_is_root) then
    write(*,'(A)')    ' Parallel execution'
    write(*,'(A)')    ' OpenMP:'
    write(*,'(A,I4)') ' -  Number of threads --> ', obj_sim_param%nthreads
  end if
#else
  if (mpi_is_root) write(*,'(A)') ' Serial execution'
  obj_sim_param%nthreads = 1
#endif

#ifdef USE_MPI
  if (mpi_is_root) then
    write(*,'(A)')    ' MPI:'
    write(*,'(A,I4)') ' -  Number of ranks   --> ', mpi_size_
  end if
#endif

  ! Solve with Q2D
  call Q2D%setup(simulation)

  ! Early-exit guard: if setup has set TODO=3 the main loop is bypassed.
  if (obj_sim_param%TODO < 3) then
    obj_sim_param%TODO = 1
    do while (obj_sim_param%TODO <= 2)
      call Q2D%solve(simulation, Dummy_Function)
      if (obj_sim_param%TODO <= 2) call Q2D%postprocess(simulation)
    end do

    call Q2D%postprocess(simulation)
  end if

  ! Release the persistent MPI requests before finalization
#ifdef USE_MPI
  call cleanup_ghost_schedule()
#endif

  ! Finalize the MPI environment
  call mpi_finalize_env()

contains

  subroutine Dummy_Function
    ! Empty subroutine to be passed as argument to Q2D%solve.
    ! It can be used for user-defined operations during the solution process.
    ! It is employed in the HYDRA coupling procedures.
  end subroutine Dummy_Function

end program Q2D_program
