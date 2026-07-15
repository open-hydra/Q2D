module Q2D_Advanced_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Lib_ORION_data,    only: ORION_data
  use Q2D_Base_Types_m,  only: Q2D_tensor_2D_type, Q2D_vector_2D_type, Q2D_d_metrics_type
  use Q2D_Series_Data_m, only: time_series_type

  implicit none
  private

  public :: block_type
  public :: bc_type
  public :: Q2D_block_type
  public :: Q2D_bc_type
  public :: Q2D_domain_type
  public :: Q2D_simulation_type

  type :: block_type
    integer                                :: dim(2)             ! Numero di celle in i-j (ghost escluse)
    real(R8), allocatable                  :: sur(:,:)           ! Superficie di cella (misura 2D, analoga al volume in 3D)
    type(Q2D_vector_2D_type), allocatable  :: node(:,:)          ! Nodi della griglia (ghost incluse)
    type(Q2D_tensor_2D_type), allocatable  :: M(:,:)             ! Tensore metrico di trasformazione
    type(Q2D_vector_2D_type), allocatable  :: dl(:,:)            ! Lunghezza media di cella (nelle direzioni i/j)
    type(Q2D_d_metrics_type)               :: dir(2)             ! Oggetto direzione: facce i, facce j (4 facce per cella in 2D)
    real(R8), allocatable                  :: yn(:,:)            ! Distanza dalla parete più vicina
    real(R8), allocatable                  :: Cross_Section(:,:) ! Sezione trasversale alle facce i (sorgente di variazione d'area)
  end type block_type

  type :: bc_type
    integer                  :: i, j, b, f                                ! Coordinate ij, blocco e faccia (4 facce in 2D)
    integer                  :: type                                      ! Tipo di BC
    integer                  :: bs, is, js, ks, fs, d11, d12, d21, d22    ! Specifiche BC 1 (connessione)
    type(Q2D_tensor_2D_type) :: Mg(2)                                     ! Tensore metrico delle celle ghost
    type(Q2D_vector_2D_type) :: dlg(2)                                    ! Lunghezza media delle celle ghost
    real(R8)                 :: surg(2)                                   ! Superficie delle celle ghost
    real(R8), allocatable    :: Pg(:,:)                                   ! Stencil delle variabili primitive nelle celle ghost
    integer                  :: ni(2)                                     ! BC chimera
    integer, allocatable     :: donorID(:,:)                              ! BC chimera
    real(R8), allocatable    :: volume_fraction(:)                        ! BC chimera
    real(R8), allocatable    :: ext_flux(:)                               ! Accoppiamento multi-solutore (BC tipo 999)
  end type bc_type

  ! Estensione a gas perfetto di un blocco.
  type, extends(block_type) :: Q2D_block_type
    real(R8), dimension(:, :, :), allocatable :: P, PO                    ! Variabili primitive ai tempi n e n-1: (i, j, nprim)
    real(R8), dimension(:, :, :), allocatable :: R                        ! Residui: (i, j, nres)
    real(R8), dimension(:, :, :), allocatable :: RS1, RS2                 ! Residui dello smoothing implicito (memoria temporanea)
    real(R8), dimension(:, :),    allocatable :: dtlocal                  ! Passo temporale locale
    type(Q2D_tensor_2D_type),     allocatable :: vel_gradient(:, :)       ! Gradiente di velocità
    real(R8),                     allocatable :: rc_term1(:, :)           ! Termine correzione rotazione/curvatura di Spalart-Shur
    real(R8),                     allocatable :: rc_term2(:, :)           ! Termine correzione rotazione/curvatura di Spalart-Shur
    real(R8),                     allocatable :: axis_r_ave(:, :)         ! Raggio medio precalcolato (sorgenti assialsimmetriche)
    real(R8),                     allocatable :: axis_inv_r(:, :)         ! 1/r precalcolato (sorgenti assialsimmetriche)
    real(R8),                     allocatable :: area_geom_term_x(:, :)   ! G_x = -(1/h)*dh/dx per le sorgenti di variazione d'area
    real(R8),                     allocatable :: area_geom_term_y(:, :)   ! G_y = -(1/h)*dh/dy per le sorgenti di variazione d'area
    real(R8),                     allocatable :: beta(:,:)                ! Flag del rilevatore d'urto (campo rank-2 in 2D)
    logical                                   :: no_chem = .false.        ! Disabilita la chimica su un singolo blocco
    logical                                   :: area_variation_active = .false. ! Attiva il termine sorgente di variazione d'area
  end type Q2D_block_type

  type, extends(bc_type) :: Q2D_bc_type
    real(R8)              :: qw, Tw, Taw, hg, qrad             ! Specifiche BC di parete viscosa
    real(R8)              :: k_rough, eps_wall                 ! Rugosità di parete ks + emissività (codici BC 301/302/303/304)
    real(R8)              :: T0, p0, alpha, beta, mach, pamb   ! Specifiche BC 4 (inflow/outflow)
    real(R8)              :: mdot                              ! Specifica di inflow a portata imposta (BC 403/404)
    real(R8)              :: rel_fac                           ! Fattore di rilassamento BC 4 (inflow/outflow)
    real(R8), allocatable :: ci(:)                             ! Frazioni massiche delle specie, BC 4 (inflow/outflow)
    real(R8)              :: SF                                ! Segnaposto inutilizzato per specifiche BC
    real(R8)              :: aCoeff, n, pRef, Taf, haf         ! Legge di combustione di Vieille del grano SRM (BC 502)
    real(R8)              :: rhoGrain = 0.0_R8                 ! Densità del grano SRM [kg/m3]
    real(R8)              :: SF_geo   = 1.0_R8                 ! Fattore di scala geometrico SRM
    type(time_series_type) :: p0time                           ! Pressione totale BC dipendente dal tempo
  end type Q2D_bc_type

  type :: Q2D_domain_type
    real(R8)                                          :: time             ! Tempo di soluzione
    real(R8)                                          :: dtglobal         ! dt globale per simulazioni time-accurate
    integer                                           :: iter, itermax    ! Numero di iterazione corrente e massimo
    integer                                           :: nb, nbound       ! Numero di blocchi, numero di facce di contorno
    integer, dimension(:, :), allocatable             :: n_bf             ! Numero di elementi bc per faccia per blocco
    type(Q2D_block_type), dimension(:), allocatable   :: blk              ! Array allocabile dei blocchi
    type(Q2D_bc_type),    dimension(:), allocatable   :: bc               ! Array allocabile delle bc
    integer                                           :: n_local_bc = 0   ! Numero di voci BC locali
                                                                          ! (usato da Ghost_Residual_Extrapolation in Lib_IRS;
                                                                          ! prerequisito condiviso con MPI)
    integer                                           :: n_local_bs = 0   ! Numero di voci BC sorgente locali
    integer, dimension(:),    allocatable             :: local_bc_idx     ! Indici delle BC locali, costruiti da
                                                                          ! Mod_GhostExchange::build_local_bc_index
    integer, dimension(:),    allocatable             :: local_bs_idx     ! Indici delle BC sorgente locali
  end type Q2D_domain_type

  type :: Q2D_simulation_type
    type(Q2D_domain_type), dimension(:), allocatable :: domain
    type(ORION_data),      dimension(:), allocatable :: IOfield
  end type Q2D_simulation_type

end module Q2D_Advanced_Types_m
