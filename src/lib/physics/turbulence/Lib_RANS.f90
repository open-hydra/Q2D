!>@brief: RANS abstract interfaces and procedure pointers (Q2D 2D-reduced superset).
!
! 7 procedure pointers + 7 abstract interfaces. Q2D superset:
!   +1 procedure pointer RANS_Diffusive_SourceQ2D (area-variation diffusive)
!   +1 abstract interface RANS_diffusive_SourceQ2D_if (idem)
!   Eddy_Viscosity_if 9-args (vs the 3D peer's 6-args): extra vr/r_ave + V_dot_G (optional)
!   RANS_diffusive_flux_if 12-args (vs the 3D peer's 9-args): extra int_l (2D analogue of
!   area) + vr/r_ave + mut optional
!
! `RSM_Symmetry` (2D restriction of the MOSE peer routine) was stripped as dead code
! pre-SSGLRR and reintroduced with the SSGLRR port (Sess 26): consumed by
! Lib_BC_Fluxes_Symmetry::BC_Symmetry_Visc `if(RSM)` block (theory §6.7).
module Q2D_Lib_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: Eddy_Viscosity, RANS_Diffusive_Flux, RANS_Diffusive_SourceQ2D
  public :: RANS_Source_Terms, RANS_Set_Wall_Values, RANS_Set_Blowing_Wall
  public :: RANS_Extrapolate_Wall, RANS_Enforce_Realizability
  public :: RSM_Symmetry

  ! Procedure pointers (8 — Q2D superset: 7 + 1 Q2D-extra)
  procedure(Eddy_Viscosity_if),             pointer :: Eddy_Viscosity
  procedure(RANS_Diffusive_Flux_if),        pointer :: RANS_Diffusive_Flux
  procedure(RANS_Diffusive_SourceQ2D_if),   pointer :: RANS_Diffusive_SourceQ2D => null()
  procedure(RANS_Source_Terms_if),          pointer :: RANS_Source_Terms
  procedure(RANS_Set_Wall_Values_if),       pointer :: RANS_Set_Wall_Values
  procedure(RANS_Set_Blowing_Wall_if),      pointer :: RANS_Set_Blowing_Wall
  procedure(RANS_Extrapolate_Wall_if),      pointer :: RANS_Extrapolate_Wall
  procedure(RANS_Enforce_Realizability_if), pointer :: RANS_Enforce_Realizability

  !> Eddy_Viscosity_if 9-args (Q2D extras).
  !> V_dot_G optional: present -> propagated to Strain_Tensor for S_33 = -V·G correction
  !> (modifies |S|^2 = |S|^2_2D + 2(V·G)^2 in Boussinesq); absent -> Definition B.
  abstract interface
    subroutine Eddy_Viscosity_if ( mut, rans_variables, mul, rho, vel_gradient, vr, r_ave, walldist, V_dot_G )
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(in), dimension(nRANS) :: rans_variables
      real(R8), intent(in)                   :: mul
      real(R8), intent(in)                   :: rho
      real(R8), intent(in), dimension(2,2)   :: vel_gradient
      real(R8), intent(in)                   :: walldist
      real(R8), intent(in)                   :: vr, r_ave
      real(R8), intent(out)                  :: mut
      real(R8), intent(in), optional         :: V_dot_G
    end subroutine Eddy_Viscosity_if
  end interface

  !> RANS_Diffusive_Flux_if 12-args (Q2D extra).
  !> mut optional: if present, eddy viscosity from caller is reused (skip recomputation
  !> in SST, ensures consistency with momentum stress mu_t under area variation).
  !> int_l is the 2D analogue of the 3D face area.
  abstract interface
    subroutine RANS_Diffusive_Flux_if(flux, rans_variables, mul, vel_gradient, rans_gradient, rho, &
                                      int_l, normal, vr, r_ave, dist, mut)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(in), dimension(nRANS)   :: rans_variables
      real(R8), intent(in)                     :: mul
      real(R8), intent(in), dimension(2,2)     :: vel_gradient
      real(R8), intent(in), dimension(nRANS,2) :: rans_gradient
      real(R8), intent(in)                     :: rho
      real(R8), intent(in)                     :: int_l, normal(2), dist
      real(R8), intent(in)                     :: vr, r_ave
      real(R8), intent(out), dimension(nRANS)  :: flux
      real(R8), intent(in), optional           :: mut
    end subroutine RANS_Diffusive_Flux_if
  end interface

  !> RANS_Diffusive_SourceQ2D_if (Q2D-extra, area-variation).
  !> Returns source(nrans,2); the caller projects onto the geometry (G or 1/r).
  abstract interface
    subroutine RANS_Diffusive_SourceQ2D_if(source, rans_variables, mul, mut, vel_gradient, rans_gradient, rho, dist)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(in), dimension(nRANS)    :: rans_variables
      real(R8), intent(in)                      :: mul
      real(R8), intent(in)                      :: mut
      real(R8), intent(in), dimension(2,2)      :: vel_gradient
      real(R8), intent(in), dimension(nRANS,2)  :: rans_gradient
      real(R8), intent(in)                      :: rho
      real(R8), intent(in)                      :: dist
      real(R8), intent(out), dimension(nRANS,2) :: source
    end subroutine RANS_Diffusive_SourceQ2D_if
  end interface

  abstract interface
    subroutine RANS_Source_Terms_if(domain)
      use Q2D_Advanced_Types_m, only: Q2D_domain_type
      implicit none
      type(Q2D_domain_type), intent(inout) :: domain
    end subroutine RANS_Source_Terms_if
  end interface

  abstract interface
    subroutine RANS_Set_Wall_Values_if(mil, rans_variables, dist)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(out), dimension(nRANS)  :: rans_variables
      real(R8), intent(in)                     :: mil
      real(R8), intent(in)                     :: dist
    end subroutine RANS_Set_Wall_Values_if
  end interface

  !> tau(2) (2D-reduced); the MOSE peer has tau(3).
  abstract interface
    subroutine RANS_Set_Blowing_Wall_if(rho, mil, tau, rans_variables, mdot, dist)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(out), dimension(nRANS)  :: rans_variables
      real(R8), intent(in)                     :: rho, mil, mdot, dist
      real(R8), intent(in), dimension(2)       :: tau
    end subroutine RANS_Set_Blowing_Wall_if
  end interface

  abstract interface
    subroutine RANS_Extrapolate_Wall_if ( prim, wall_rans, rho, wall_rho, ghost_rans )
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(in),  dimension(nRANS) :: prim, wall_rans
      real(R8), intent(in) :: rho, wall_rho
      real(R8), intent(out), dimension(nRANS) :: ghost_rans
    end subroutine RANS_Extrapolate_Wall_if
  end interface

  abstract interface
    subroutine RANS_Enforce_Realizability_if ( var )
      use iso_fortran_env, only: R8 => real64
      implicit none
      real(R8), intent(inout),  dimension(:) :: var
    end subroutine RANS_Enforce_Realizability_if
  end interface

contains

  ! Reynolds-stress values at a symmetry interface and in the ghost cell
  ! (theory §6.7). Q2D 2D reduction (4-stress + omega; R13=R23 dropped, 7->5):
  ! the symmetry plane is an x- or y-plane (no z-plane in 2D). The normal stresses
  ! R11,R22,R33 and the scalar omega are EVEN (zero-gradient ghost); the in-plane
  ! shear R12 is ODD at any in-plane symmetry plane (one index along the plane
  ! normal), so its face value is zero and its ghost sign-flips. This is the 2D
  ! restriction of the MOSE peer Lib_RANS::RSM_Symmetry (x-plane: R12,R13 odd;
  ! y-plane: R12,R23 odd) with the out-of-plane cross-stresses absent.
  ! var = (rho*R11, rho*R22, rho*R33, rho*R12, rho*omega); n is the face normal.
  subroutine RSM_Symmetry ( var, var_w, var_g, n )
    use Q2D_Global_m, only: nRANS
    implicit none
    real(R8), intent(in),  dimension(nRANS) :: var
    real(R8), intent(out), dimension(nRANS) :: var_w, var_g
    real(R8), intent(in),  dimension(2)     :: n

    ! Normal stresses + omega: zero-gradient (even)
    var_w(1:3) = var(1:3)
    var_w(5)   = var(5)
    var_g(1:3) = var(1:3)
    var_g(5)   = var(5)

    ! In-plane shear R12: odd at any in-plane symmetry plane
    var_w(4) = 0d0
    var_g(4) = - var(4)

  end subroutine RSM_Symmetry

end module Q2D_Lib_RANS
