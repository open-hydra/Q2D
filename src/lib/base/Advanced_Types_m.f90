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
    integer                                :: dim(2)             ! Number of cells in i-j (ghosts excluded)
    real(R8), allocatable                  :: sur(:,:)           ! Cell surface (2D measure, analogous to the volume in 3D)
    type(Q2D_vector_2D_type), allocatable  :: node(:,:)          ! Grid nodes (ghosts included)
    type(Q2D_tensor_2D_type), allocatable  :: M(:,:)             ! Transformation metric tensor
    type(Q2D_vector_2D_type), allocatable  :: dl(:,:)            ! Mean cell length (in the i/j directions)
    type(Q2D_d_metrics_type)               :: dir(2)             ! Direction object: i faces, j faces (4 faces per cell in 2D)
    real(R8), allocatable                  :: yn(:,:)            ! Distance to the nearest wall
    real(R8), allocatable                  :: Cross_Section(:,:) ! Cross-section at the i faces (area-variation source)
  end type block_type

  type :: bc_type
    integer                  :: i, j, b, f                                ! ij coordinates, block and face (4 faces in 2D)
    integer                  :: type                                      ! BC type
    integer                  :: bs, is, js, ks, fs, d11, d12, d21, d22    ! BC 1 specifications (connection)
    type(Q2D_tensor_2D_type) :: Mg(2)                                     ! Metric tensor of the ghost cells
    type(Q2D_vector_2D_type) :: dlg(2)                                    ! Mean length of the ghost cells
    real(R8)                 :: surg(2)                                   ! Surface of the ghost cells
    real(R8), allocatable    :: Pg(:,:)                                   ! Stencil of the primitive variables in the ghost cells
    integer                  :: ni(2)                                     ! Chimera BC
    integer, allocatable     :: donorID(:,:)                              ! Chimera BC
    real(R8), allocatable    :: volume_fraction(:)                        ! Chimera BC
    real(R8), allocatable    :: ext_flux(:)                               ! Multi-solver coupling (BC type 999)
  end type bc_type

  ! Perfect-gas extension of a block.
  type, extends(block_type) :: Q2D_block_type
    real(R8), dimension(:, :, :), allocatable :: P, PO                    ! Primitive variables at times n and n-1: (i, j, nprim)
    real(R8), dimension(:, :, :), allocatable :: R                        ! Residuals: (i, j, nres)
    real(R8), dimension(:, :, :), allocatable :: RS1, RS2                 ! Implicit-smoothing residuals (temporary storage)
    real(R8), dimension(:, :),    allocatable :: dtlocal                  ! Local time step
    type(Q2D_tensor_2D_type),     allocatable :: vel_gradient(:, :)       ! Velocity gradient
    real(R8),                     allocatable :: rc_term1(:, :)           ! Spalart-Shur rotation/curvature correction term
    real(R8),                     allocatable :: rc_term2(:, :)           ! Spalart-Shur rotation/curvature correction term
    real(R8),                     allocatable :: axis_r_ave(:, :)         ! Precomputed mean radius (axisymmetric sources)
    real(R8),                     allocatable :: axis_inv_r(:, :)         ! Precomputed 1/r (axisymmetric sources)
    real(R8),                     allocatable :: area_geom_term_x(:, :)   ! G_x = -(1/h)*dh/dx for the area-variation sources
    real(R8),                     allocatable :: area_geom_term_y(:, :)   ! G_y = -(1/h)*dh/dy for the area-variation sources
    real(R8),                     allocatable :: beta(:,:)                ! Shock-detector flag (rank-2 field in 2D)
    logical                                   :: no_chem = .false.        ! Disable chemistry on a single block
    logical                                   :: area_variation_active = .false. ! Enable the area-variation source term
  end type Q2D_block_type

  type, extends(bc_type) :: Q2D_bc_type
    real(R8)              :: qw, Tw, Taw, hg, qrad             ! Viscous-wall BC specifications
    real(R8)              :: k_rough, eps_wall                 ! Wall roughness ks + emissivity (BC codes 301/302/303/304)
    real(R8)              :: T0, p0, alpha, beta, mach, pamb   ! BC 4 specifications (inflow/outflow)
    real(R8)              :: mdot                              ! Imposed-mass-flow inflow specification (BC 403/404)
    real(R8)              :: rel_fac                           ! BC 4 relaxation factor (inflow/outflow)
    real(R8), allocatable :: ci(:)                             ! Species mass fractions, BC 4 (inflow/outflow)
    real(R8)              :: SF                                ! Unused placeholder for BC specifications
    real(R8)              :: aCoeff, n, pRef, Taf              ! SRM grain (BC 502): mdot_ref, n, p_ref, T_af (total)
    logical               :: srm_inflow = .false.              ! SRM grain of mode 2: applied as 404 (IO_BC)
    type(time_series_type) :: p0time                           ! Time-dependent BC total pressure
  end type Q2D_bc_type

  type :: Q2D_domain_type
    real(R8)                                          :: time             ! Solution time
    real(R8)                                          :: dtglobal         ! Global dt for time-accurate simulations
    integer                                           :: iter, itermax    ! Current and maximum iteration number
    integer                                           :: nb, nbound       ! Number of blocks, number of boundary faces
    integer, dimension(:, :), allocatable             :: n_bf             ! Number of bc elements per face per block
    type(Q2D_block_type), dimension(:), allocatable   :: blk              ! Allocatable array of blocks
    type(Q2D_bc_type),    dimension(:), allocatable   :: bc               ! Allocatable array of bcs
    integer                                           :: n_local_bc = 0   ! Number of local BC entries
                                                                          ! (used by Ghost_Residual_Extrapolation in Lib_IRS;
                                                                          ! shared prerequisite with MPI)
    integer                                           :: n_local_bs = 0   ! Number of local source BC entries
    integer, dimension(:),    allocatable             :: local_bc_idx     ! Indices of the local BCs, built by
                                                                          ! Mod_GhostExchange::build_local_bc_index
    integer, dimension(:),    allocatable             :: local_bs_idx     ! Indices of the local source BCs
  end type Q2D_domain_type

  type :: Q2D_simulation_type
    type(Q2D_domain_type), dimension(:), allocatable :: domain
    type(ORION_data),      dimension(:), allocatable :: IOfield
  end type Q2D_simulation_type

end module Q2D_Advanced_Types_m
