module Q2D_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Parameters_m, only: hlen, llen, clen

  implicit none
  private

  !! ------------------------------------------------------
  !! Simulation Parameters --------------------------------
  !! ------------------------------------------------------
  type :: simulation_parameters_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    logical             :: newrun                 ! Restart flag
    real(R8)            :: res_threshold          ! Residual threshold for stopping
    real(R8)            :: time_threshold         ! Maximum simulation time
    integer             :: iter_threshold         ! Maximum iteration count
    character(len=clen) :: simulation_type        ! Type of simulation (euler, navier-stokes)
    ! Useful variables
    integer             :: iter_general           ! Number of iterations - including all MG levels
    integer             :: iter_from_call
    real(R8)            :: time_from_call
    real(R8), allocatable :: residuotot(:)
    integer             :: nthreads               ! Number of threads for simulation
    real(R8)            :: cputime(2)             ! Simulation time duration
    integer             :: TODO                   ! Decide solve and/or postprocess (state-machine flag, Wrap_Setup)
    ! HYDRA-coupling flags. These are COUPLING flags (NOT MASTER legacy): Q2D can
    ! be both standalone-build capable AND HYDRA-coupling capable. Defaults
    ! `.false.` preserve the Q2D standalone-build behavior; the ATLAS top-level
    ! orchestrator may set them `.true.` for HYDRA coupling.
    logical             :: HYDRA_time_accurate = .false. ! ATLAS-driven dtglobal via Set_Global_dt
    logical             :: HYDRA_postprocess   = .false. ! ATLAS-driven solution write
    logical             :: HYDRA_MG            = .false. ! ATLAS-driven HYDRA-Multigrid INI section
  end type simulation_parameters_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Input-Output -----------------------------------------
  !! ------------------------------------------------------
  type :: io_t
    character(len=llen)  :: warning_message
    character(len=llen)  :: error_message
    character(len=llen)  :: description
    ! USER-DEFINED INPUTS
    integer              :: sol_diter, res_diter   ! Iteration interval to save the solution / residuals
    real(R8)             :: sol_dtime              ! Time interval to save the solution
    logical              :: sol_overwrite          ! Switch to overwrite the solution
    character(len=llen)  :: sol_format, ini_format ! Solution format (tecplot, vtk) and formatting (ascii, raw)
    character(len=llen)  :: sol_variables, wall_variables
    integer              :: shell_diter            ! Shell update
    integer              :: ini_diter              ! input.ini update
    logical              :: lineprobe_restart      ! Q2D-extra (peculiarità #5): line-probe restart/append flag, read once at [<codename>-LineProbes] level
    ! Useful variables
    character(len=llen)  :: nameinit               ! Initial file name
    logical              :: write_thermo, write_transport, write_composition
    character(len=hlen)  :: Ovarnames, ORANSname
    integer              :: Onvar
    real(R8)             :: IOtime
    logical              :: write_wall, write_wall_mechanical, write_wall_thermal, write_wall_mass
    character(len=hlen)  :: Owallnames
    integer              :: Onwall
    integer              :: unitRES                ! Residuals file unit
    character(len=llen)  :: unitRES_format         ! Format for residual_history_file
  end type io_t

  type :: io_probes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: file
    character(len=hlen) :: varnames
    integer             :: iloc(3)                 ! 2D index (b, i, j) (vs the 3D peer's iloc(4) = (b,i,j,k))
    real(R8)            :: loc(2)                  ! 2D coordinates (x, y) (vs the 3D peer's loc(3))
    ! Useful variables
    ! ...
  end type io_probes_t


  ! Q2D-extra (peculiarità #5): line-probe config — sampling along an i/j fixed-index
  ! line. MOSE has no line-probes. Mirror of io_probes_t with line-probe specifics
  ! (name + 2-component iloc=[block,direction] + scalar fixed-index instead of loc).
  type :: io_lineprobes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: name                    ! line-probe section name (= output file base)
    character(len=hlen) :: varnames                ! whitespace-separated variable list
    integer             :: iloc(2)                 ! [block, direction] (direction=1 i-fixed, 2 j-fixed)
    integer             :: fixed_index             ! fixed-axis sample index along the line
  end type io_lineprobes_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! NUMERICAL SCHEME !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Time Scheme ------------------------------------------
  !! ------------------------------------------------------
  type :: time_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    character(len=llen) :: solver_type             ! Solver type (explicit/implicit, euler, RK2, RK3)
    real(R8)            :: cfl                     ! Convective stability parameter
    real(R8)            :: vnn                     ! Diffusive stability parameter
    integer             :: rampa_cfl_iter          ! Number of iterations for CFL ramp
    logical             :: time_accurate           ! Time-accurate switch
    integer             :: n_RK
    character(len=llen) :: integration_variables   ! Integration variables (cons/prim)
    character(len=llen) :: omp_strategy            ! OMP parallelism strategy ('do' | 'task'), Q2D-extra
    ! Useful variables
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
    ! USER-DEFINED INPUTS
    real(R8)  :: beta
    ! Useful variables
    logical   :: enabled
  end type irs_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Space Discretization ---------------------------------
  !! ------------------------------------------------------
  type :: space_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    character(len=llen) :: space_reconstruction    ! Space reconstruction method
    character(len=llen) :: flux_limiter            ! Flux limiter for space reconstruction
    ! Useful variables
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
    ! USER-DEFINED INPUTS
    real(R8) :: Minf
    ! Useful variables
    logical  :: SD
  end type riemann_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Shock Detector ---------------------------------------
  !! ------------------------------------------------------
  type :: shock_detector_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    ! ...
    ! Useful variables
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
    ! USER-DEFINED INPUTS
    integer              :: MGL                    ! Number of multigrid levels
    integer, allocatable :: iter_threshold(:)      ! Number of iterations for each level
    ! Useful variables
    integer              :: MG_level               ! Current MG-level being solved
    logical              :: change_MG              ! Logical variable to change MG-level
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
    ! USER-DEFINED INPUTS
    ! ...
    ! Useful variables
    logical, allocatable :: viscous_flag(:,:)      ! 2D shape (nb, 4 faces) (vs the 3D peer's 6 faces)
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
    ! USER-DEFINED INPUTS
    character(len=llen)   :: model
    character(len=llen)   :: ode_name
    integer               :: max_ode_steps
    integer               :: iopt(3)
    real(R8), allocatable :: RT(:), AT(:)
    integer(I4), allocatable :: no_chem_list(:)
    ! Useful variables
    character(len=llen)   :: mechanism_name
    logical               :: use_strang
    integer               :: imodel               !> 0 no chemistry; 1 finite-rate; 2 equilibrium
    integer               :: N_strang             !> 1 w/o strang;   2 w strang
    real(R8)              :: strangcoeff          !> 1 w/o strang; 0.5 w strang
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
    ! USER-DEFINED INPUTS
    character(len=llen)   :: model
    logical               :: enabled
    logical               :: LL91
    logical               :: LIN
    ! Useful variables
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
    ! USER-DEFINED INPUTS
    real(R8) :: Sc      ! Schmidt laminare
    real(R8) :: Sct     ! Schmidt turbolento
    real(R8) :: Prt     ! Prandtl turbolento
    character(len=llen) :: model
    logical :: SAcomp, SpalartShur, SAR
    logical :: QCR2000, blowing_corr, k_energy_coupling
    ! Useful variables
    logical :: RSM, SD
  end type rans_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Rotating-frame ---------------------------------------
  !! ------------------------------------------------------
  ! axis(3)/origin(3)/Omega_vec(3) kept as placeholders (in 2D planar the
  ! rotation axis is perpendicular to the plane).
  type :: rot_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! USER-DEFINED INPUTS
    character(len=llen)   :: model
    real(R8)              :: omega                   !> Angular speed [rad/s]
    real(R8)              :: axis(3)                 !> Unit rotation axis (in 2D planar it is (0,0,1))
    real(R8)              :: origin(3)               !> A point on the axis [m]
    real(R8)              :: Omega_vec(3) = 0.0_R8   !> omega * axis  [rad/s]
    character(len=128)    :: axis_str
    character(len=128)    :: origin_str
    ! Useful variables
    logical               :: enabled
    integer               :: n_stationary = 0
    character(len=32), allocatable :: stationary_face_str(:)
    integer, allocatable           :: stationary_face(:,:)
  end type rot_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! Q2D EXTRA !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Source terms (Q2D-specific)
  !! Geometric sources: `axisymmetric` (1/r), `area_variation` (G = ∇h, V·G).
  !! ------------------------------------------------------
  type :: source_t
    character(len=llen)              :: warning_message
    character(len=llen)              :: error_message
    character(len=llen)              :: description
    ! USER-DEFINED INPUTS
    character(len=llen)              :: active_str    ! Raw INI string from [Q2D-Sources] active=
    character(len=llen), allocatable :: models(:)     ! Parsed source-model names
    ! Useful variables
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
  !! Basic placeholder (thermo/transport stubs) ----
  !! ------------------------------------------------------
  type :: basic_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    ! ...
    ! Useful variables
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
  type(io_lineprobes_t), allocatable,  public :: obj_io_lineprobes(:)   ! Q2D-extra (peculiarità #5)
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
