module Q2D_Read_Physics
  use iso_fortran_env,    only: I4 => int32, R8 => real64
  use Q2D_Parameters_m,   only: llen, codename
  use Q2D_Config_Types_m, only: obj_chemistry, obj_soot, obj_rans, obj_rot, obj_sim_param
  use Q2D_Registry_m,     only: reg

  implicit none
  private
  public :: Register_Physics

contains

  subroutine Register_Physics()
    use Q2D_Global_m,  only: nsc
    use IR_Precision,  only: str
    implicit none
    integer :: m  ! Stationary-face counter (NOT a spatial loop index — no `k` as a 3D index).
    character(len=:), allocatable :: section

    !! ------------------------------------------------------
    !! Physical models --------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Physics'
    call reg%add(section, 'equations', obj_sim_param%simulation_type, '', 'Gas dynamics equations', 'euler , navier-stokes', .true.)
    ! Active Q2D RANS models: SA, SST, Wilcox2006, SSGLRR (SSG/LRR-omega RSM,
    ! 4-stress + omega; -SD selects the shock-detector diffusive flux variant).
    call reg%add(section, 'turbulence', obj_rans%model, 'none', 'RANS turbulence model', &
                 'SA, SA-R, SA-RC, SAcomp, SST, Wilcox2006, SSGLRR, SSGLRR-SD, none', .false.)
    call reg%add(section, 'chemistry', obj_chemistry%model, 'frozen', 'Chemistry model', 'frozen, finite-rate, equilibrium', .false.)
    call reg%add(section, 'soot-generation', obj_soot%model, 'none', 'Soot generation model', 'LL91, LIN, none', .false.)
    call reg%add(section, 'rotational-frame', obj_rot%model, 'none', 'Rotational frame model', 'rigid-body, none', .false.)

    !! ------------------------------------------------------
    !! Chemistry --------------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Chemistry'
    ! Exclude blocks from chemistry: read unbounded (no length cap) by Scan_NoChem
    ! (Q2D_Backend_INI), NOT via a Registry string param which would truncate at llen.
    ! ODE solver selection
    call reg%add(section, 'ode-solver', obj_chemistry%ode_name, 'H-radau5', 'ODE solver for chemistry', 'H-radau5, sdirk4b, ros4', .false.)
    ! ODE solver parameters
    call reg%add(section, 'ode-max-steps', obj_chemistry%max_ode_steps, '100000', 'Maximum ODE integration steps', '> 0', .false.)
    call reg%add(section, 'ode-relative-tol-species', obj_chemistry%RT(1:nsc), '1e-5', 'ODE relative tolerance for species', '> 0', .false.)
    call reg%add(section, 'ode-relative-tol-temperature', obj_chemistry%RT(nsc+1), '1e-5', 'ODE relative tolerance for temperature', '> 0', .false.)
    call reg%add(section, 'ode-absolute-tol-species', obj_chemistry%AT(1:nsc), '1e-5', 'ODE absolute tolerance for species', '> 0', .false.)
    call reg%add(section, 'ode-absolute-tol-temperature', obj_chemistry%AT(nsc+1), '1e-5', 'ODE absolute tolerance for temperature', '> 0', .false.)


    !! ------------------------------------------------------
    !! Turbulence -------------------------------------------
    !! ------------------------------------------------------
    ! Sezione [Q2D-Turbulence]. Q2D legacy IO_RANS.f90 leggeva da [Q2D-RANS]
    ! con auto-discovery: contenuto assorbito qui.
    section = trim(codename)//'-Turbulence'
    call reg%add(section, 'Prt', obj_rans%Prt, '0.85', 'Turbulent Prandtl number', '> 0', .false.)
    call reg%add(section, 'Sct', obj_rans%Sct, '0.90', 'Turbulent Schmidt number', '> 0', .false.)
    call reg%add(section, 'Sc', obj_rans%Sc, '0.7', 'Schmidt number', '> 0', .false.)
    ! k-coupling Q2D-extra registry: the MOSE peer leaves the field commented out
    ! (`k_energy_coupling` hardcoded `.false.`). Q2D makes the option
    ! user-configurable via INI. The feature is functionally identical in both
    ! (cf. Lib_SST).
    call reg%add(section, 'k-coupling', obj_rans%k_energy_coupling, '.false.', 'Turbulent kinetic energy coupling', 'logical', .false.)


    !! ------------------------------------------------------
    !! Rotating Frame ---------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Rotating-Frame'
    call reg%add(section, 'omega',  obj_rot%omega,      '0.0',         'Angular speed [rad/s]', '>= 0', .false.)
    ! axis(3)/origin(3) kept as placeholders (in 2D planar: axis perpendicular to the plane).
    call reg%add(section, 'axis',   obj_rot%axis_str,   '0.0 0.0 1.0', 'Rotation axis direction (3 components)', '', .false.)
    call reg%add(section, 'origin', obj_rot%origin_str, '0.0 0.0 0.0', 'Point on the rotation axis [m]', '', .false.)
    ! Stationary face entries (count scanned in Scan_Ini, strings allocated there)
    do m = 1, obj_rot%n_stationary
      call reg%add(section, 'stationary-face-'//trim(str(.true.,m)), obj_rot%stationary_face_str(m),'', 'Stationary face: block_index face_direction', '', .false.)
    end do

    !! ------------------------------------------------------
    !! GSI --------------------------------------------------
    !! ------------------------------------------------------
    !! TODO: Refactor GSI material properties to registry
    !! For now: HDPE, PP, Paraffin, HTPB, Carbon properties at module-level defaults
    !! in Q2D_Lib_GSI (including the Q2D-extra GSI_Paraffin).
    !! These will require dedicated sections in config types — TODO.
    !! ======================================================

  end subroutine Register_Physics

end module Q2D_Read_Physics
