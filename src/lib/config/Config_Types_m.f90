module Q2D_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Parameters_m, only: hlen, llen, clen

  implicit none
  private

  !! ------------------------------------------------------
  !! Simulation parameters --------------------------------
  !! ------------------------------------------------------
  type :: simulation_parameters_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    logical             :: newrun                 ! Restart flag
    real(R8)            :: res_threshold          ! Residual threshold for stopping
    real(R8)            :: time_threshold         ! Maximum simulation time
    integer             :: iter_threshold         ! Maximum number of iterations
    character(len=clen) :: simulation_type        ! Simulation type (euler, navier-stokes)
    ! Service variables
    integer             :: iter_general           ! Number of iterations, including all MG levels
    integer             :: iter_from_call
    real(R8)            :: time_from_call
    real(R8), allocatable :: residuotot(:)
    integer             :: nthreads               ! Number of threads for the simulation
    real(R8)            :: cputime(2)             ! Simulation duration
    integer             :: TODO                   ! Selects solve and/or postprocess (status flag, Wrap_Setup)
    ! HYDRA coupling flags: Q2D can run either standalone or coupled via HYDRA.
    ! The `.false.` defaults correspond to the standalone behaviour; the
    ! high-level ATLAS orchestrator may set them to `.true.` for HYDRA coupling.
    logical             :: HYDRA_time_accurate = .false. ! dtglobal imposed by ATLAS via Set_Global_dt
    logical             :: HYDRA_postprocess   = .false. ! solution writing driven by ATLAS
    logical             :: HYDRA_MG            = .false. ! HYDRA-Multigrid INI section driven by ATLAS
  end type simulation_parameters_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Input-output -----------------------------------------
  !! ------------------------------------------------------
  type :: io_t
    character(len=llen)  :: warning_message
    character(len=llen)  :: error_message
    character(len=llen)  :: description
    ! USER INPUT
    integer              :: sol_diter, res_diter   ! Iteration interval to save the solution / the residuals
    real(R8)             :: sol_dtime              ! Time interval to save the solution
    logical              :: sol_overwrite          ! Switch to overwrite the solution
    character(len=llen)  :: sol_format, ini_format ! Solution format (tecplot, vtk) and encoding (ascii, raw)
    character(len=llen)  :: sol_variables, wall_variables
    integer              :: shell_diter            ! Shell update
    integer              :: ini_diter              ! input.ini update
    logical              :: lineprobe_restart      ! Line-probe restart/append flag (section [<codename>-LineProbes])
    ! Service variables
    character(len=llen)  :: nameinit               ! Name of the initial file
    logical              :: write_thermo, write_transport, write_composition
    character(len=hlen)  :: Ovarnames, ORANSname
    integer              :: Onvar
    real(R8)             :: IOtime
    logical              :: write_wall, write_wall_mechanical, write_wall_thermal, write_wall_mass
    character(len=hlen)  :: Owallnames
    integer              :: Onwall
    integer              :: unitRES                ! Residual file unit
    character(len=llen)  :: unitRES_format         ! Format for residual_history_file
  end type io_t

  type :: io_probes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: file
    character(len=hlen) :: varnames
    integer             :: iloc(3)                 ! 2D indices (b, i, j)
    real(R8)            :: loc(2)                  ! 2D coordinates (x, y)
    ! Service variables
    ! ...
  end type io_probes_t


  ! Line-probe configuration: sampling along a line at fixed i or j index.
  ! Structure analogous to io_probes_t, with the line-probe specifics
  ! (name + 2-component iloc [block, direction] + scalar fixed_index in place of loc).
  type :: io_lineprobes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: name                    ! line-probe section name (= base of the output file name)
    character(len=hlen) :: varnames                ! whitespace-separated list of variables
    integer             :: iloc(2)                 ! [block, direction] (direction=1 i fixed, 2 j fixed)
    integer             :: fixed_index             ! sampling index on the fixed axis of the line
  end type io_lineprobes_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! NUMERICAL SCHEME !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Time scheme ------------------------------------------
  !! ------------------------------------------------------
  type :: time_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    character(len=llen) :: solver_type             ! Solver type (explicit/implicit, euler, RK2, RK3)
    real(R8)            :: cfl                     ! Convective stability parameter
    real(R8)            :: vnn                     ! Diffusive stability parameter
    integer             :: rampa_cfl_iter          ! Number of CFL ramp iterations
    logical             :: time_accurate           ! Time-accurate switch
    integer             :: n_RK
    character(len=llen) :: integration_variables   ! Integration variables (cons/prim)
    character(len=llen) :: omp_strategy            ! OMP parallelism strategy ('do' | 'task')
    ! Service variables
    ! ...
  end type time_scheme_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Implicit residual smoothing --------------------------
  !! ------------------------------------------------------
  type irs_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    real(R8)  :: beta
    ! Service variables
    logical   :: enabled
  end type irs_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Spatial discretization -------------------------------
  !! ------------------------------------------------------
  type :: space_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    character(len=llen) :: space_reconstruction    ! Spatial reconstruction method
    character(len=llen) :: flux_limiter            ! Flux limiter for the spatial reconstruction
    ! Service variables
    logical :: SD
  end type space_scheme_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Riemann solver ---------------------------------------
  !! ------------------------------------------------------
  type :: riemann_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    real(R8) :: Minf
    ! Service variables
    logical  :: SD
  end type riemann_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Shock detector ---------------------------------------
  !! ------------------------------------------------------
  type :: shock_detector_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    ! ...
    ! Service variables
    logical :: SD
  end type shock_detector_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Multigrid --------------------------------------------
  !! ------------------------------------------------------
  type :: multigrid_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    integer              :: MGL                    ! Number of multigrid levels
    integer, allocatable :: iter_threshold(:)      ! Number of iterations for each level
    ! Service variables
    integer              :: MG_level               ! MG level currently being solved
    logical              :: change_MG              ! Logical variable to change MG level
  end type multigrid_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! BC ---------------------------------------------------
  !! ------------------------------------------------------
  type :: io_bc_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    ! ...
    ! Service variables
    logical, allocatable :: viscous_flag(:,:)      ! 2D shape (nb, 4 faces)
    logical, allocatable :: coupling_flag(:,:)     ! 2D shape (nb, 4 faces)
  end type io_bc_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! PHYSICS !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Chemistry --------------------------------------------
  !! ------------------------------------------------------
  type :: chemistry_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! USER INPUT
    character(len=llen)   :: model
    character(len=llen)   :: ode_name
    integer               :: max_ode_steps
    integer               :: iopt(3)
    real(R8), allocatable :: RT(:), AT(:)
    integer(I4), allocatable :: no_chem_list(:)
    ! Service variables
    character(len=llen)   :: mechanism_name
    logical               :: use_strang
    integer               :: imodel               !> 0 without chemistry; 1 finite-rate; 2 equilibrium
    integer               :: N_strang             !> 1 without strang;  2 with strang
    real(R8)              :: strangcoeff          !> 1 without strang; 0.5 with strang
  end type chemistry_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Soot -------------------------------------------------
  !! ------------------------------------------------------
  type :: soot_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! USER INPUT
    character(len=llen)   :: model
    logical               :: enabled
    logical               :: LL91
    logical               :: LIN
    ! Service variables
    real(R8)              :: rho_soot
  end type soot_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Turbulence closure -----------------------------------
  !! ------------------------------------------------------
  type :: rans_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! USER INPUT
    real(R8) :: Sc      ! Laminar Schmidt number
    real(R8) :: Sct     ! Turbulent Schmidt number
    real(R8) :: Prt     ! Turbulent Prandtl number
    character(len=llen) :: model
    logical :: SAcomp, SpalartShur, SAR
    logical :: QCR2000, blowing_corr, k_energy_coupling
    ! Service variables
    logical :: RSM, SD
  end type rans_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Rotating reference frame -----------------------------
  !! ------------------------------------------------------
  ! axis(3)/origin(3)/Omega_vec(3) kept as placeholders (in planar 2D the
  ! rotation axis is perpendicular to the plane).
  type :: rot_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! USER INPUT
    character(len=llen)   :: model
    real(R8)              :: omega                   !> Angular speed [rad/s]
    real(R8)              :: axis(3)                 !> Rotation axis unit vector (in planar 2D it is (0,0,1))
    real(R8)              :: origin(3)               !> A point on the axis [m]
    real(R8)              :: Omega_vec(3) = 0.0_R8   !> omega * axis  [rad/s]
    character(len=128)    :: axis_str
    character(len=128)    :: origin_str
    ! Service variables
    logical               :: enabled
    integer               :: n_stationary = 0
    character(len=32), allocatable :: stationary_face_str(:)
    integer, allocatable           :: stationary_face(:,:)
  end type rot_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! SOURCE TERMS !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Geometric source terms:
  !! `axisymmetric` (1/r terms), `area_variation` (G = ∇h, V·G).
  !! ------------------------------------------------------
  type :: source_t
    character(len=llen)              :: warning_message
    character(len=llen)              :: error_message
    character(len=llen)              :: description
    ! USER INPUT
    character(len=llen)              :: active_str    ! Raw INI string from [Q2D-Sources] active=
    character(len=llen), allocatable :: models(:)     ! Source model names after parsing
    ! Service variables
    integer                          :: n_models = 0
    logical                          :: axisymmetric = .false.
    logical                          :: area_variation = .false.
  end type source_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! OTHER TYPES !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Basic placeholder (stub for thermo/transport) --------
  !! ------------------------------------------------------
  type :: basic_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER INPUT
    ! ...
    ! Service variables
    ! ...
  end type basic_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! SINGLETON OBJECTS !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(simulation_parameters_t), save, public :: obj_sim_param
  type(io_t),                    save, public :: obj_io
  type(io_probes_t), allocatable,      public :: obj_io_probes(:)
  type(io_lineprobes_t), allocatable,  public :: obj_io_lineprobes(:)   ! line-probe configurations
  type(io_bc_t),                 save, public :: obj_io_bc
  type(time_scheme_t),           save, public :: obj_time_scheme
  type(irs_t),                   save, public :: obj_irs
  type(space_scheme_t),          save, public :: obj_space_scheme
  type(riemann_t),               save, public :: obj_riemann
  type(shock_detector_t),        save, public :: obj_shock_detector
  type(multigrid_t),             save, public :: obj_multigrid
  type(chemistry_t),             save, public :: obj_chemistry
  type(soot_t),                  save, public :: obj_soot
  type(rans_t),                  save, public :: obj_rans
  type(rot_t),                   save, public :: obj_rot
  type(basic_t),                 save, public :: obj_thermo
  type(basic_t),                 save, public :: obj_transport
  type(source_t),                save, public :: obj_source
  !! ------------------------------------------------------
  !! ------------------------------------------------------

end module Q2D_Config_Types_m
