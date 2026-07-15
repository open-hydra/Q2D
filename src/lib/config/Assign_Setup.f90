module Q2D_Assign_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use IR_Precision,    only: str

  implicit none
  private
  public :: Assign_Setup

contains

  subroutine Assign_Setup()
    use Q2D_Config_Types_m, only: obj_sim_param, obj_rans, obj_riemann, obj_space_scheme, &
                                   obj_shock_detector, obj_time_scheme, obj_irs, obj_transport
    use Q2D_Global_m,       only: model
    use Q2D_Lib_Chemistry,  only: Setup_Chemistry
    use Q2D_Lib_Strang,     only: Setup_Strang_Splitting
    use Q2D_Mod_RANS,       only: Setup_RANS_Model
    use Q2D_Mod_Space,      only: Setup_Space_Scheme
    use Q2D_Mod_Fluxes,     only: Assign_Fluxes
    use Q2D_Mod_Riemann,    only: Assign_Riemann_Solver
    use Q2D_Mod_Newstate,   only: Assign_Integration_Variables
    use Q2D_IO_Solution,    only: Setup_Input_Solution
    use Q2D_Mod_Soot,       only: Setup_Soot
    use Q2D_Lib_RotatingFrame, only: Setup_RotatingFrame
    implicit none

    ! Set the simulation type
    if (obj_sim_param%simulation_type=='euler') then
      model = 0
    elseif (obj_sim_param%simulation_type=='navier-stokes') then
      if (obj_rans%model == 'none' .or. obj_rans%model == '') then
        model = 1
        obj_rans%model = 'none'
      else
        model = 2
      end if
    endif

    ! Set up the input solution
    call Setup_Input_Solution()

    ! Space
    call Setup_Space_Scheme()
    call Assign_Riemann_Solver()
    ! OMP dispatcher: assigns the procedure pointer `Fluxes` to Fluxes_omp_do or
    ! Fluxes_omp_task according to obj_time_scheme%omp_strategy ('do' default;
    ! 'task' uses block-level parallelism + intra-block taskloop).
    call Assign_Fluxes( obj_time_scheme%omp_strategy )
    if ( obj_riemann%SD .or. obj_space_scheme%SD) then
      obj_shock_detector%SD = .true.
    else
      obj_shock_detector%SD = .false.
    end if

    ! Time
    if (obj_time_scheme%solver_type /= 'euler') then
      read(obj_time_scheme%solver_type(3:3), *) obj_time_scheme%n_rk
    else
      obj_time_scheme%n_rk = 1
    end if
    if (obj_irs%beta>0d0) obj_irs%enabled = .true.
    call Assign_Integration_Variables()

    ! Chemistry assignment.
    ! Setup_Chemistry must precede Setup_Strang_Splitting because it sets
    ! obj_chemistry%use_strang from obj_chemistry%model (cf. Lib_Strang).
    call Setup_Chemistry()
    call Setup_Strang_Splitting()

    ! Soot model assignment.
    call Setup_Soot()

    ! RANS model assignment.
    ! Called after Setup_Chemistry/Setup_Strang_Splitting to keep a consistent
    ! setup ordering.
    call Setup_RANS_Model()

    ! Rotating-frame assignment. Order: after Setup_RANS_Model.
    ! Setup_RotatingFrame applies a double check (axisymmetric=.false. + axis=[0,0,1]);
    ! the rotating frame in 2D axisymmetric mode is not supported.
    call Setup_RotatingFrame()

    !! Descriptions, warnings and errors

    ! Simulation type
    if (obj_sim_param%simulation_type == 'euler') then
      obj_sim_param%description = 'Euler'
    else if (obj_sim_param%simulation_type == 'navier-stokes') then
      obj_sim_param%description = 'Navier-Stokes'
    end if
    ! Time scheme
    if (obj_time_scheme%solver_type == 'euler') then
      obj_time_scheme%description = 'Explicit Euler'
    else if (obj_time_scheme%solver_type == 'RK2') then
      obj_time_scheme%description = 'Second-order Runge-Kutta'
    else if (obj_time_scheme%solver_type == 'RK3') then
      obj_time_scheme%description = 'Third-order Runge-Kutta'
    end if
    if (obj_time_scheme%time_accurate) then
      obj_time_scheme%description = trim(obj_time_scheme%description)//' with time-accurate switch enabled'
    end if
    if (obj_irs%enabled) then
      obj_irs%description = 'Beta set to '//trim(str(.true.,real(obj_irs%beta)))
    end if
    ! Space scheme
    ! ... written in Mod_Space ...
    ! Riemann solver
    if (trim(obj_riemann%description) == 'AUSM+-up' .or. trim(obj_riemann%description) == 'AUSM+-up2') then
      if (obj_riemann%Minf == 0.0d0) then
        obj_riemann%error_message = '[ERROR] AUSM+-up solver selected. Minf must be defined in input.'
      end if
    endif
    ! Transport properties
    if (model>0 .and. obj_transport%description=='Unavailable') &
    write(*,'(A)') '[ERROR] Transport properties are unavailable for the selected phase: cannot run Navier-Stokes simulation'

  end subroutine Assign_Setup

end module Q2D_Assign_Setup
