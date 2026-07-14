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
    integer                                :: dim(2)             ! Number of cells in i-j (ghost not included)
    real(R8), allocatable                  :: sur(:,:)           ! Cell surface (sur replaces the 3D vol)
    type(Q2D_vector_2D_type), allocatable  :: node(:,:)          ! Mesh grid points (including ghost)
    type(Q2D_tensor_2D_type), allocatable  :: M(:,:)             ! Metric transformation tensor
    type(Q2D_vector_2D_type), allocatable  :: dl(:,:)            ! Average cell length (in i/j direction)
    type(Q2D_d_metrics_type)               :: dir(2)             ! Direction object: i-faces, j-faces (4 faces 2D vs the 3D peer's 6)
    real(R8), allocatable                  :: yn(:,:)            ! Nearest wall distance
    real(R8), allocatable                  :: Cross_Section(:,:) ! Cross section at i faces (area-variation source)
  end type block_type

  type :: bc_type
    integer                  :: i, j, b, f                                ! ij coordinates, block and face (4 faces 2D)
    integer                  :: type                                      ! BC type
    integer                  :: bs, is, js, ks, fs, d11, d12, d21, d22    ! BC 1 (connection) specifications
    type(Q2D_tensor_2D_type) :: Mg(2)                                     ! Ghost cell metric tensor
    type(Q2D_vector_2D_type) :: dlg(2)                                    ! Ghost cell average cell length
    real(R8)                 :: surg(2)                                   ! Ghost cell surface
    real(R8), allocatable    :: Pg(:,:)                                   ! Ghost cell primitive stencil
    integer                  :: ni(2)                                     ! BC chimera
    integer, allocatable     :: donorID(:,:)                              ! BC chimera
    real(R8), allocatable    :: volume_fraction(:)                        ! BC chimera
    real(R8), allocatable    :: ext_flux(:)                               ! Multi-solver coupling (Q2D 999 == MOSE 103)
  end type bc_type

  ! Perfect gas extension for a block.
  type, extends(block_type) :: Q2D_block_type
    real(R8), dimension(:, :, :), allocatable :: P, PO                    ! Primitive variables at time n and n-1: (i, j, nprim)
    real(R8), dimension(:, :, :), allocatable :: R                        ! Residuals: (i, j, nres)
    real(R8), dimension(:, :, :), allocatable :: RS1, RS2                 ! Implicit smoothing residuals (temporary storage)
    real(R8), dimension(:, :),    allocatable :: dtlocal                  ! Local time step
    type(Q2D_tensor_2D_type),     allocatable :: vel_gradient(:, :)       ! Gradient of velocity
    real(R8),                     allocatable :: rc_term1(:, :)           ! Spalart-Shur rotation/curvature correction term
    real(R8),                     allocatable :: rc_term2(:, :)           ! Spalart-Shur rotation/curvature correction term
    real(R8),                     allocatable :: axis_r_ave(:, :)         ! Precomputed average radius for axisymmetric source terms
    real(R8),                     allocatable :: axis_inv_r(:, :)         ! Precomputed 1/r for axisymmetric source terms
    real(R8),                     allocatable :: area_geom_term_x(:, :)   ! G_x = -(1/h)*dh/dx for area-variation sources
    real(R8),                     allocatable :: area_geom_term_y(:, :)   ! G_y = -(1/h)*dh/dy for area-variation sources
    real(R8),                     allocatable :: beta(:,:)                ! Shock detector flag (rank-3 reduced to rank-2 for 2D)
    logical                                   :: no_chem = .false.        ! Disable chemistry on a single block
    logical                                   :: area_variation_active = .false. ! Activate area-variation source term
  end type Q2D_block_type

  type, extends(bc_type) :: Q2D_bc_type
    real(R8)              :: qw, Tw, Taw, hg, qrad             ! BC viscous wall specifications
    real(R8)              :: k_rough, eps_wall                 ! Wall roughness ks + emissivity (BC codes 301/302/303/304)
    real(R8)              :: T0, p0, alpha, beta, mach, pamb   ! BC 4 (inflow/outflow) specifications
    real(R8)              :: mdot                              ! BC 403/404 mass-flux inflow specification
    real(R8)              :: rel_fac                           ! BC 4 (inflow/outflow) relaxation factor
    real(R8), allocatable :: ci(:)                             ! BC 4 (inflow/outflow) species mass fractions
    real(R8)              :: SF                                ! BC ? specifications (unused placeholder)
    real(R8)              :: aCoeff, n, pRef, Taf, haf         ! BC 502 SRM grain Vieille burn-rate
    real(R8)              :: rhoGrain = 0.0_R8                 ! SRM grain density [kg/m3]
    real(R8)              :: SF_geo   = 1.0_R8                 ! SRM geometric scale factor
    type(time_series_type) :: p0time                           ! BC time-dependent total pressure
  end type Q2D_bc_type

  type :: Q2D_domain_type
    real(R8)                                          :: time             ! Solution time
    real(R8)                                          :: dtglobal         ! Global dt for time-accurate simulation
    integer                                           :: iter, itermax    ! Iteration number
    integer                                           :: nb, nbound       ! Number of blocks, number of boundary faces
    integer, dimension(:, :), allocatable             :: n_bf             ! Number of bc elements per face per block
    type(Q2D_block_type), dimension(:), allocatable   :: blk              ! Allocatable block array
    type(Q2D_bc_type),    dimension(:), allocatable   :: bc               ! Allocatable bc array
    integer                                           :: n_local_bc = 0   ! Number of local BC entries (Lib_IRS Ghost_Residual_Extrapolation; shared prerequisite with MPI)
    integer                                           :: n_local_bs = 0   ! Number of local BC-source entries
    integer, dimension(:),    allocatable             :: local_bc_idx     ! Local BC indices (Mod_GhostExchange::build_local_bc_index)
    integer, dimension(:),    allocatable             :: local_bs_idx     ! Local BC-source indices
  end type Q2D_domain_type

  type :: Q2D_simulation_type
    type(Q2D_domain_type), dimension(:), allocatable :: domain
    type(ORION_data),      dimension(:), allocatable :: IOfield
  end type Q2D_simulation_type

end module Q2D_Advanced_Types_m
