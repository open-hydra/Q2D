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
    integer :: m  ! Contatore delle facce stazionarie (non è un indice di ciclo spaziale).
    character(len=:), allocatable :: section

    !! ------------------------------------------------------
    !! Modelli fisici ---------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Physics'
    call reg%add(section, 'equations', obj_sim_param%simulation_type, '', 'Gas dynamics equations', 'euler , navier-stokes', .true.)
    ! Modelli RANS attivi: SA, SST, Wilcox2006, SSGLRR (RSM SSG/LRR-omega,
    ! 4 stress + omega; il suffisso -SD seleziona la variante del flusso
    ! diffusivo con shock detector).
    call reg%add(section, 'turbulence', obj_rans%model, 'none', 'RANS turbulence model', &
                 'SA, SA-R, SA-RC, SAcomp, SST, Wilcox2006, SSGLRR, SSGLRR-SD, none', .false.)
    call reg%add(section, 'chemistry', obj_chemistry%model, 'frozen', 'Chemistry model', 'frozen, finite-rate, equilibrium', .false.)
    call reg%add(section, 'soot-generation', obj_soot%model, 'none', 'Soot generation model', 'LL91, LIN, none', .false.)
    call reg%add(section, 'rotational-frame', obj_rot%model, 'none', 'Rotational frame model', 'rigid-body, none', .false.)

    !! ------------------------------------------------------
    !! Chimica ----------------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Chemistry'
    ! Blocchi esclusi dalla chimica: lettura senza limite di lunghezza in Scan_NoChem
    ! (Q2D_Backend_INI), NON tramite un parametro stringa del Registry che troncherebbe a llen.
    ! Selezione del solutore ODE
    call reg%add(section, 'ode-solver', obj_chemistry%ode_name, 'H-radau5', 'ODE solver for chemistry', 'H-radau5, sdirk4b, ros4', .false.)
    ! Parametri del solutore ODE
    call reg%add(section, 'ode-max-steps', obj_chemistry%max_ode_steps, '100000', 'Maximum ODE integration steps', '> 0', .false.)
    call reg%add(section, 'ode-relative-tol-species', obj_chemistry%RT(1:nsc), '1e-5', 'ODE relative tolerance for species', '> 0', .false.)
    call reg%add(section, 'ode-relative-tol-temperature', obj_chemistry%RT(nsc+1), '1e-5', 'ODE relative tolerance for temperature', '> 0', .false.)
    call reg%add(section, 'ode-absolute-tol-species', obj_chemistry%AT(1:nsc), '1e-5', 'ODE absolute tolerance for species', '> 0', .false.)
    call reg%add(section, 'ode-absolute-tol-temperature', obj_chemistry%AT(nsc+1), '1e-5', 'ODE absolute tolerance for temperature', '> 0', .false.)


    !! ------------------------------------------------------
    !! Turbolenza -------------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Turbulence'
    call reg%add(section, 'Prt', obj_rans%Prt, '0.85', 'Turbulent Prandtl number', '> 0', .false.)
    call reg%add(section, 'Sct', obj_rans%Sct, '0.90', 'Turbulent Schmidt number', '> 0', .false.)
    call reg%add(section, 'Sc', obj_rans%Sc, '0.7', 'Schmidt number', '> 0', .false.)
    ! Accoppiamento dell'energia cinetica turbolenta k nell'equazione
    ! dell'energia, configurabile dall'utente via INI (cf. Lib_SST).
    call reg%add(section, 'k-coupling', obj_rans%k_energy_coupling, '.false.', 'Turbulent kinetic energy coupling', 'logical', .false.)


    !! ------------------------------------------------------
    !! Sistema di riferimento rotante -----------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-Rotating-Frame'
    call reg%add(section, 'omega',  obj_rot%omega,      '0.0',         'Angular speed [rad/s]', '>= 0', .false.)
    ! axis(3)/origin(3) mantenuti come placeholder (nel 2D piano l'asse è perpendicolare al piano).
    call reg%add(section, 'axis',   obj_rot%axis_str,   '0.0 0.0 1.0', 'Rotation axis direction (3 components)', '', .false.)
    call reg%add(section, 'origin', obj_rot%origin_str, '0.0 0.0 0.0', 'Point on the rotation axis [m]', '', .false.)
    ! Voci delle facce stazionarie (conteggio letto in Scan_Ini, stringhe allocate lì)
    do m = 1, obj_rot%n_stationary
      call reg%add(section, 'stationary-face-'//trim(str(.true.,m)), obj_rot%stationary_face_str(m),'', 'Stationary face: block_index face_direction', '', .false.)
    end do

    !! ------------------------------------------------------
    !! GSI --------------------------------------------------
    !! ------------------------------------------------------
    !! TODO: portare le proprietà dei materiali GSI nel registry.
    !! Per ora: proprietà di HDPE, PP, paraffina, HTPB e carbonio come default a
    !! livello di modulo in Q2D_Lib_GSI (inclusa GSI_Paraffin).
    !! Richiederanno sezioni dedicate nei config types — TODO.
    !! ======================================================

  end subroutine Register_Physics

end module Q2D_Read_Physics
