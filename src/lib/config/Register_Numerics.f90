module Q2D_Read_Numerics
  use iso_fortran_env,    only: I4 => int32, R8 => real64
  use Q2D_Config_Types_m, only: obj_time_scheme, obj_irs, obj_space_scheme, &
                                 obj_riemann, obj_multigrid, obj_sim_param
  use Q2D_Registry_m,     only: reg

  implicit none
  private
  public :: Register_Numerics

contains

  subroutine Register_Numerics(nmgl)
    use Q2D_Parameters_m, only: llen, codename
    use IR_Precision,     only: str
    implicit none
    ! Variabili locali
    integer, intent(in) :: nmgl
    character(len=llen) :: section

    section = trim(codename)//'-Numerics'

    !! ------------------------------------------------------
    !! Schema temporale -------------------------------------
    !! ------------------------------------------------------
    obj_time_scheme%warning_message = 'none'
    obj_time_scheme%error_message   = 'none'
    obj_time_scheme%description     = 'none'

    ! Tipo di solutore
    call reg%add( trim(section), 'time-scheme', obj_time_scheme%solver_type, 'euler', 'Time integration solver', 'euler, RK2, RK3', .true. )

    ! Coefficienti di stabilità e opzioni correlate
    call reg%add( trim(section), 'cfl', obj_time_scheme%cfl, '0.5', 'CFL number', '> 0', .true. )
    call reg%add( trim(section), 'vnn', obj_time_scheme%vnn, '0.3', 'VNN parameter', '> 0', .false. )
    call reg%add( trim(section), 'cfl-rise-threshold', obj_time_scheme%rampa_cfl_iter, '0', 'CFL rise threshold', '>= 0', .false. )

    ! Interruttore time-accurate
    call reg%add( trim(section), 'time-accurate', obj_time_scheme%time_accurate, '.false.', 'Time accurate switch', 'logical', .true. )

    ! Variabili di integrazione del nuovo stato
    call reg%add( trim(section), 'integration-variables', obj_time_scheme%integration_variables, 'cons', 'Integration variables (cons/prim)', 'cons ,  prim', .false. )

    ! Strategia di parallelismo OMP ('do' | 'task')
    call reg%add( trim(section), 'OMP-parallel', obj_time_scheme%omp_strategy, 'do', 'OMP parallelism strategy', 'do, task', .false. )


    ! Implicit residual smoothing --------------------------
    obj_irs%description     = 'none'
    obj_irs%warning_message = 'none'
    obj_irs%error_message   = 'none'
    call reg%add( trim(section), 'irs', obj_irs%enabled, '.false.', 'Implicit Residual Smoothing', 'logical', .false. )
    call reg%add( trim(section), 'irs-beta', obj_irs%beta, '0.0', 'IRS beta parameter', '>= 0', .false. )

    !! ------------------------------------------------------
    !! ------------------------------------------------------


    !! ------------------------------------------------------
    !! Schema spaziale --------------------------------------
    !! ------------------------------------------------------

    ! Discretizzazione spaziale ----------------------------
    obj_space_scheme%description     = 'none'
    obj_space_scheme%warning_message = 'none'
    obj_space_scheme%error_message   = 'none'
    call reg%add( trim(section), 'space-reconstruction', obj_space_scheme%space_reconstruction, '', 'Space reconstruction method', 'MUSCL-SD, MUSCL, first-order', .true. )
    call reg%add( trim(section), 'flux-limiter', obj_space_scheme%flux_limiter, '', 'Flux limiter for space reconstruction', 'vanalbada, minmod, superbee, vanleer, mc, LIMO3', .false. )


    ! Solutore di Riemann ----------------------------------
    obj_riemann%description     = 'none'
    obj_riemann%warning_message = 'none'
    obj_riemann%error_message   = 'none'
    call reg%add( trim(section), 'riemann-solver', obj_riemann%description, 'HLLC', 'Riemann solver', 'SLAU, SLAU2, HLLC+, HLLE++, HLLE, HLLEM, HLLC, AUSM+, AUSM+-up, AUSM+-up2, exact', .false. )
    call reg%add( trim(section), 'riemann-options-Minf', obj_riemann%Minf, '0.0', 'Mach infinity for AUSM+-up', '>= 0', .false. )


    ! Livelli multigrid -------------------------------------
    call Register_Multigrid_Levels(nmgl)

    !! ------------------------------------------------------
    !! ------------------------------------------------------

  end subroutine Register_Numerics


  subroutine Register_Multigrid_Levels(nmgl)
    use Q2D_Parameters_m,    only: llen, codename
    use Q2D_Config_Types_m,  only: obj_sim_param
    use IR_Precision,        only: str
    implicit none
    integer, intent(in) :: nmgl
    integer :: m
    character(len=llen) :: option

    ! Allocazione dell'array delle soglie di iterazione per ciascun livello
    obj_multigrid%MGL = max(nmgl, 1)
    allocate( obj_multigrid%iter_threshold(obj_multigrid%MGL) )
    obj_multigrid%iter_threshold = 0

    ! Registrazione dei parametri di iterazione per livello (level1-iter,
    ! level2-iter, ...). ATLAS può impostare `obj_sim_param%HYDRA_MG = .true.`
    ! per registrare la sezione `[HYDRA-Multigrid]` al posto di
    ! `[<codename>-Multigrid]`. Il default `.false.` corrisponde al layout del
    ! registry dell'esecuzione standalone di Q2D.
    if (.not. obj_sim_param%HYDRA_MG) then
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( trim(codename)//'-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .false. )
      enddo
    else
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( 'HYDRA-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .false. )
      enddo
    end if

  end subroutine Register_Multigrid_Levels

end module Q2D_Read_Numerics
