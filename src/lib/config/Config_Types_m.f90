module Q2D_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Parameters_m, only: hlen, llen, clen

  implicit none
  private

  !! ------------------------------------------------------
  !! Parametri di simulazione -----------------------------
  !! ------------------------------------------------------
  type :: simulation_parameters_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! INPUT UTENTE
    logical             :: newrun                 ! Flag di restart
    real(R8)            :: res_threshold          ! Soglia sui residui per l'arresto
    real(R8)            :: time_threshold         ! Tempo massimo di simulazione
    integer             :: iter_threshold         ! Numero massimo di iterazioni
    character(len=clen) :: simulation_type        ! Tipo di simulazione (euler, navier-stokes)
    ! Variabili di servizio
    integer             :: iter_general           ! Numero di iterazioni, inclusi tutti i livelli MG
    integer             :: iter_from_call
    real(R8)            :: time_from_call
    real(R8), allocatable :: residuotot(:)
    integer             :: nthreads               ! Numero di thread per la simulazione
    real(R8)            :: cputime(2)             ! Durata della simulazione
    integer             :: TODO                   ! Seleziona solve e/o postprocess (flag di stato, Wrap_Setup)
    ! Flag di accoppiamento HYDRA: Q2D può funzionare sia standalone sia
    ! accoppiato via HYDRA. I default `.false.` corrispondono al comportamento
    ! standalone; l'orchestratore di alto livello ATLAS può impostarli `.true.`
    ! per l'accoppiamento HYDRA.
    logical             :: HYDRA_time_accurate = .false. ! dtglobal imposto da ATLAS via Set_Global_dt
    logical             :: HYDRA_postprocess   = .false. ! scrittura della soluzione pilotata da ATLAS
    logical             :: HYDRA_MG            = .false. ! sezione INI HYDRA-Multigrid pilotata da ATLAS
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
    ! INPUT UTENTE
    integer              :: sol_diter, res_diter   ! Intervallo di iterazioni per salvare la soluzione / i residui
    real(R8)             :: sol_dtime              ! Intervallo di tempo per salvare la soluzione
    logical              :: sol_overwrite          ! Interruttore per sovrascrivere la soluzione
    character(len=llen)  :: sol_format, ini_format ! Formato della soluzione (tecplot, vtk) e codifica (ascii, raw)
    character(len=llen)  :: sol_variables, wall_variables
    integer              :: shell_diter            ! Aggiornamento della shell
    integer              :: ini_diter              ! Aggiornamento di input.ini
    logical              :: lineprobe_restart      ! Flag restart/append delle line-probe (sezione [<codename>-LineProbes])
    ! Variabili di servizio
    character(len=llen)  :: nameinit               ! Nome del file iniziale
    logical              :: write_thermo, write_transport, write_composition
    character(len=hlen)  :: Ovarnames, ORANSname
    integer              :: Onvar
    real(R8)             :: IOtime
    logical              :: write_wall, write_wall_mechanical, write_wall_thermal, write_wall_mass
    character(len=hlen)  :: Owallnames
    integer              :: Onwall
    integer              :: unitRES                ! Unità del file dei residui
    character(len=llen)  :: unitRES_format         ! Formato per residual_history_file
  end type io_t

  type :: io_probes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! INPUT UTENTE
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: file
    character(len=hlen) :: varnames
    integer             :: iloc(3)                 ! Indici 2D (b, i, j)
    real(R8)            :: loc(2)                  ! Coordinate 2D (x, y)
    ! Variabili di servizio
    ! ...
  end type io_probes_t


  ! Configurazione line-probe: campionamento lungo una linea a indice i o j fissato.
  ! Struttura analoga a io_probes_t, con le specificità delle line-probe
  ! (name + iloc a 2 componenti [blocco, direzione] + fixed_index scalare al posto di loc).
  type :: io_lineprobes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! INPUT UTENTE
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: name                    ! nome della sezione line-probe (= base del nome file di output)
    character(len=hlen) :: varnames                ! lista di variabili separate da spazi
    integer             :: iloc(2)                 ! [blocco, direzione] (direzione=1 i fissato, 2 j fissato)
    integer             :: fixed_index             ! indice di campionamento sull'asse fissato della linea
  end type io_lineprobes_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! SCHEMA NUMERICO !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Schema temporale -------------------------------------
  !! ------------------------------------------------------
  type :: time_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! INPUT UTENTE
    character(len=llen) :: solver_type             ! Tipo di solutore (esplicito/implicito, euler, RK2, RK3)
    real(R8)            :: cfl                     ! Parametro di stabilità convettivo
    real(R8)            :: vnn                     ! Parametro di stabilità diffusivo
    integer             :: rampa_cfl_iter          ! Numero di iterazioni della rampa CFL
    logical             :: time_accurate           ! Interruttore time-accurate
    integer             :: n_RK
    character(len=llen) :: integration_variables   ! Variabili di integrazione (cons/prim)
    character(len=llen) :: omp_strategy            ! Strategia di parallelismo OMP ('do' | 'task')
    ! Variabili di servizio
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
    ! INPUT UTENTE
    real(R8)  :: beta
    ! Variabili di servizio
    logical   :: enabled
  end type irs_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Discretizzazione spaziale ----------------------------
  !! ------------------------------------------------------
  type :: space_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! INPUT UTENTE
    character(len=llen) :: space_reconstruction    ! Metodo di ricostruzione spaziale
    character(len=llen) :: flux_limiter            ! Flux limiter per la ricostruzione spaziale
    ! Variabili di servizio
    logical :: SD
  end type space_scheme_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Solutore di Riemann ----------------------------------
  !! ------------------------------------------------------
  type :: riemann_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! INPUT UTENTE
    real(R8) :: Minf
    ! Variabili di servizio
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
    ! INPUT UTENTE
    ! ...
    ! Variabili di servizio
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
    ! INPUT UTENTE
    integer              :: MGL                    ! Numero di livelli multigrid
    integer, allocatable :: iter_threshold(:)      ! Numero di iterazioni per ciascun livello
    ! Variabili di servizio
    integer              :: MG_level               ! Livello MG attualmente risolto
    logical              :: change_MG              ! Variabile logica per cambiare livello MG
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
    ! INPUT UTENTE
    ! ...
    ! Variabili di servizio
    logical, allocatable :: viscous_flag(:,:)      ! Shape 2D (nb, 4 facce)
    logical, allocatable :: coupling_flag(:,:)     ! Shape 2D (nb, 4 facce)
  end type io_bc_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FISICA !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Chimica ----------------------------------------------
  !! ------------------------------------------------------
  type :: chemistry_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! INPUT UTENTE
    character(len=llen)   :: model
    character(len=llen)   :: ode_name
    integer               :: max_ode_steps
    integer               :: iopt(3)
    real(R8), allocatable :: RT(:), AT(:)
    integer(I4), allocatable :: no_chem_list(:)
    ! Variabili di servizio
    character(len=llen)   :: mechanism_name
    logical               :: use_strang
    integer               :: imodel               !> 0 senza chimica; 1 finite-rate; 2 equilibrio
    integer               :: N_strang             !> 1 senza strang;  2 con strang
    real(R8)              :: strangcoeff          !> 1 senza strang; 0.5 con strang
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
    ! INPUT UTENTE
    character(len=llen)   :: model
    logical               :: enabled
    logical               :: LL91
    logical               :: LIN
    ! Variabili di servizio
    real(R8)              :: rho_soot
  end type soot_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Chiusura di turbolenza -------------------------------
  !! ------------------------------------------------------
  type :: rans_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! INPUT UTENTE
    real(R8) :: Sc      ! Numero di Schmidt laminare
    real(R8) :: Sct     ! Numero di Schmidt turbolento
    real(R8) :: Prt     ! Numero di Prandtl turbolento
    character(len=llen) :: model
    logical :: SAcomp, SpalartShur, SAR
    logical :: QCR2000, blowing_corr, k_energy_coupling
    ! Variabili di servizio
    logical :: RSM, SD
  end type rans_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Sistema di riferimento rotante -----------------------
  !! ------------------------------------------------------
  ! axis(3)/origin(3)/Omega_vec(3) mantenuti come placeholder (nel 2D piano
  ! l'asse di rotazione è perpendicolare al piano).
  type :: rot_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! INPUT UTENTE
    character(len=llen)   :: model
    real(R8)              :: omega                   !> Velocità angolare [rad/s]
    real(R8)              :: axis(3)                 !> Versore dell'asse di rotazione (nel 2D piano è (0,0,1))
    real(R8)              :: origin(3)               !> Un punto sull'asse [m]
    real(R8)              :: Omega_vec(3) = 0.0_R8   !> omega * axis  [rad/s]
    character(len=128)    :: axis_str
    character(len=128)    :: origin_str
    ! Variabili di servizio
    logical               :: enabled
    integer               :: n_stationary = 0
    character(len=32), allocatable :: stationary_face_str(:)
    integer, allocatable           :: stationary_face(:,:)
  end type rot_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!! TERMINI SORGENTE !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Termini sorgente geometrici:
  !! `axisymmetric` (termini 1/r), `area_variation` (G = ∇h, V·G).
  !! ------------------------------------------------------
  type :: source_t
    character(len=llen)              :: warning_message
    character(len=llen)              :: error_message
    character(len=llen)              :: description
    ! INPUT UTENTE
    character(len=llen)              :: active_str    ! Stringa INI grezza da [Q2D-Sources] active=
    character(len=llen), allocatable :: models(:)     ! Nomi dei modelli sorgente dopo il parsing
    ! Variabili di servizio
    integer                          :: n_models = 0
    logical                          :: axisymmetric = .false.
    logical                          :: area_variation = .false.
  end type source_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! ALTRI TIPI !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Placeholder di base (stub per thermo/transport) ------
  !! ------------------------------------------------------
  type :: basic_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! INPUT UTENTE
    ! ...
    ! Variabili di servizio
    ! ...
  end type basic_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! OGGETTI SINGLETON !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(simulation_parameters_t), save, public :: obj_sim_param
  type(io_t),                    save, public :: obj_io
  type(io_probes_t), allocatable,      public :: obj_io_probes(:)
  type(io_lineprobes_t), allocatable,  public :: obj_io_lineprobes(:)   ! configurazioni line-probe
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
