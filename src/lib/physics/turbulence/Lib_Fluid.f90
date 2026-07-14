!>@brief: Fluid mechanics relations (Q2D 2D-reduced + Definition A).
!
! Q2D superset features:
!   - 4×Stress_Vector + 4×Stress_Tensor + 2×Strain_Tensor realizations (Q2D
!     superset for the axisymmetric and Definition A regimes);
!   - Stress_t(3,3) and Strain E(3,3) even in 2D — a legitimate exception to the
!     No-w principle (S_33 = -V·G is dilatational, NOT an out-of-plane velocity
!     component);
!   - V_dot_G optional: present -> Definition A (S_33 = -V·G); absent ->
!     Definition B (default 0, backward compat).
!
! Theory: S_33 = -⟨∂w/∂z⟩ = -(V·G), the BC-consistent depth-averaged out-of-plane
! strain (Definition A); see the q2D area-variation derivation.
module Q2D_Lib_Fluid
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  real(R8), private, parameter :: negativetwothird = -2.d0/3.d0

  ! Procedure pointers (3 — vs the MOSE peer's single Stress_Vector)
  procedure(Stress_Vector_if), pointer, public :: Stress_Vector
  procedure(Stress_Tensor_if), pointer, public :: Stress_Tensor
  procedure(Strain_Tensor_if), pointer, public :: Strain_Tensor

  public :: Stress_Vector_Std, Stress_Vector_Axi
  public :: Stress_Vector_2eq_Std, Stress_Vector_2eq_Axi
  public :: Stress_Vector_RSM_Std, Stress_Vector_RSM_Axi
  public :: Stress_Tensor_Std, Stress_Tensor_Axi
  public :: Stress_Tensor_2eq_Std, Stress_Tensor_2eq_Axi
  public :: Stress_Tensor_RSM_Std, Stress_Tensor_RSM_Axi
  public :: Strain_Tensor_Std, Strain_Tensor_Axi
  public :: Vorticity_Vector, Vorticity_Tensor   ! used by Lib_Spalart/Lib_SST (Q2D legacy esponeva come accessibili)
  public :: V_dot_G_face_avg

  ! Abstract interface Stress_Vector:
  ! V_dot_G optional: present -> Definition A (rigorous z-reduction with S_33 = -V·G);
  ! absent -> Definition B (default 0, backward compat with legacy call sites).
  abstract interface
    pure function Stress_Vector_if ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
      use iso_fortran_env, only: R8 => real64
      implicit none
      real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
      real(R8), intent(in), optional :: V_dot_G
      real(R8) :: Stress(2)
    end function Stress_Vector_if
  end interface

  abstract interface
    pure function Stress_Tensor_if ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
      use iso_fortran_env, only: R8 => real64
      implicit none
      real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
      real(R8), intent(in), optional :: V_dot_G
      real(R8) :: Stress_t(3,3)
    end function Stress_Tensor_if
  end interface

  abstract interface
    pure function Strain_Tensor_if ( Gradvel, vr, r_ave, V_dot_G ) result ( E )
      use iso_fortran_env, only: R8 => real64
      implicit none
      real(R8), intent(in) :: Gradvel(2,2), r_ave, vr
      real(R8), intent(in), optional :: V_dot_G
      real(R8) :: E(3,3)
    end function Strain_Tensor_if
  end interface

contains

  pure function Stress_Vector_Std ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! Q2D Definition A: V_dot_G correction enables (V·G) terms in Divervel and tau_zz.
    ! Default V_dot_G = 0 preserves Definition B (legacy backward compat).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy
    real(R8) :: mi, lam, Divervel, lamDivervel, VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1)
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2)
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_Std


  pure function Stress_Vector_Axi ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! V_dot_G optional argument for interface uniformity (NOT used in Axi:
    ! axisymmetric flows have no virtual top/bot wall in Def A sense).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy
    real(R8) :: mi, lam, Divervel, lamDivervel

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    if (r_ave .le. 1.d-6) then
      Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(2,2) ! de l'Hopital vr/r --> dv/dr = Gradvel(2,2)
    else
      Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave     ! cylindrical coordinate divergence
    endif
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1)
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2)
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_Axi


  pure function Stress_Tensor_Std ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! Q2D Definition A (when V_dot_G present): tau_zz = -(4/3)*mi*V_dot_G - (2/3)*mi*div_2D
    ! and Divervel = div_2D - V_dot_G to capture the z-divergence contribution.
    ! V_dot_G absent (default 0) recovers Definition B (legacy: tau_zz = lam*div_2D).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy, Tau_zz
    real(R8) :: mi, lam, Divervel, lamDivervel, VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1)
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2)
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )
    Tau_zz = lamDivervel + 2d0 * mi * (-VG)

    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz

    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy

    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_Std


  pure function Stress_Tensor_Axi ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! V_dot_G optional argument for interface uniformity (NOT used in Axi).
    ! Tensor evaluated cell-center (r=dr/2 > 0 by primal/dual mesh construction)
    ! -> no de l'Hopital fallback needed.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy, Tau_zz
    real(R8) :: mi, lam, Divervel, lamDivervel

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave ! cylindrical coordinate divergence
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1)
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2)
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )
    Tau_zz = 2d0 * mi * vr / r_ave + lamDivervel

    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz

    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy

    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_Axi


  pure function Stress_Vector_2eq_Std ( Gradvel, Normal, mi_l, mi_t, rhok, r_ave, vr, V_dot_G ) result ( Stress )
    ! Q2D Definition A (V_dot_G present): Divervel = div_2D - V_dot_G to capture
    ! z-divergence at top/bot walls. Boussinesq-2eq adds rhokdelta = -(2/3)ρk to
    ! diagonal. V_dot_G absent (default 0) recovers Definition B (legacy).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, rhok(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy
    real(R8) :: mi, lam, Divervel, lamDivervel, rhokdelta, VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG
    lamDivervel = lam*Divervel
    rhokdelta = negativetwothird * rhok(1)

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1) + rhokdelta
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2) + rhokdelta
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_2eq_Std


  pure function Stress_Vector_2eq_Axi ( Gradvel, Normal, mi_l, mi_t, rhok, r_ave, vr, V_dot_G ) result ( Stress )
    ! V_dot_G optional argument for interface uniformity (NOT used in Axi).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, rhok(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy
    real(R8) :: mi, lam, Divervel, lamDivervel, rhokdelta

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    if (r_ave .le. 1.d-6) then
      Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(2,2) ! de l'Hopital
    else
      Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave
    endif
    lamDivervel = lam*Divervel
    rhokdelta = negativetwothird * rhok(1)

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1) + rhokdelta
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2) + rhokdelta
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_2eq_Axi


  pure function Stress_Tensor_2eq_Std ( Gradvel, mi_l, mi_t, rhok, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! Q2D Definition A: Divervel = div_2D - V_dot_G,
    ! Tau_zz = lam*Divervel - 2*mi*VG + rhokdelta = -(4mi/3)*VG -(2mi/3)*div_2D -(2/3)rhok.
    ! Boussinesq-2eq adds rhokdelta = -(2/3)rho*k to diagonal.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, rhok(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy, Tau_zz
    real(R8) :: mi, lam, Divervel, lamDivervel, rhokdelta, VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG
    lamDivervel = lam*Divervel
    rhokdelta = negativetwothird * rhok(1)

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1) + rhokdelta
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2) + rhokdelta
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )
    Tau_zz = lamDivervel + 2d0 * mi * (-VG) + rhokdelta

    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz

    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy
    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_2eq_Std


  pure function Stress_Tensor_2eq_Axi ( Gradvel, mi_l, mi_t, rhok, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! V_dot_G optional argument for interface uniformity (NOT used in Axi).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, rhok(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy, Tau_zz
    real(R8) :: mi, lam, Divervel, lamDivervel, rhokdelta

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave
    lamDivervel = lam*Divervel
    rhokdelta = negativetwothird * rhok(1)

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1) + rhokdelta
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2) + rhokdelta
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )
    Tau_zz = 2d0 * mi * vr / r_ave + lamDivervel + rhokdelta

    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz

    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy
    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_2eq_Axi


  pure function Stress_Vector_RSM_Std ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! SSG/LRR-omega RSM momentum coupling (no Boussinesq, theory §6.8): the total
    ! face stress is the LAMINAR Newtonian tensor (mi_l only, mi_t ignored) MINUS
    ! the transported Reynolds stress rho*R_ij. var = Prim(nt:) so
    ! var(1)=rho*R11, var(2)=rho*R22, var(4)=rho*R12 (in-plane block; R13=R23=0).
    ! Definition A: Divervel = div_2D - V·G (z-divergence at virtual walls).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy
    real(R8) :: lam, Divervel, lamDivervel, VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    lam = negativetwothird * mi_l
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) - var(1)
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) - var(2)
    Tau_xy = mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) - var(4)

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_RSM_Std


  pure function Stress_Vector_RSM_Axi ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! Axisymmetric RSM momentum coupling: laminar (mi_l) cylindrical divergence
    ! (de l'Hopital at the axis) MINUS rho*R_ij. V_dot_G accepted for interface
    ! uniformity, NOT used in Axi.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy
    real(R8) :: lam, Divervel, lamDivervel

    lam = negativetwothird * mi_l
    if (r_ave .le. 1.d-6) then
      Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(2,2) ! de l'Hopital vr/r --> dv/dr
    else
      Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave     ! cylindrical divergence
    endif
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) - var(1)
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) - var(2)
    Tau_xy = mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) - var(4)

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_RSM_Axi


  pure function Stress_Tensor_RSM_Std ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! 3x3 total stress for the RSM geometric momentum source (Source_Area_Visc_Blk
    ! G·(tau_33 - tau_ii) and Source_Axi_Visc_Blk (tau_rr - tau_theta)/r): laminar
    ! (mi_l) tensor MINUS rho*R_ij. var(3)=rho*R33 enters tau_33. Definition A:
    ! tau_33^lam = lam*Divervel - 2*mi_l*VG, Divervel = div_2D - V·G.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy, Tau_zz
    real(R8) :: lam, Divervel, lamDivervel, VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    lam = negativetwothird * mi_l
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) - var(1)
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) - var(2)
    Tau_xy = mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) - var(4)
    Tau_zz = lamDivervel + 2d0 * mi_l * (-VG) - var(3)

    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz
    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy
    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_RSM_Std


  pure function Stress_Tensor_RSM_Axi ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! Axisymmetric 3x3 total stress: laminar (mi_l) cylindrical tensor MINUS
    ! rho*R_ij; tau_33^lam = 2*mi_l*vr/r + lam*Divervel (theta-theta). V_dot_G
    ! accepted for interface uniformity, NOT used in Axi.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    real(R8) :: Tau_xx, Tau_xy, Tau_yy, Tau_zz
    real(R8) :: lam, Divervel, lamDivervel

    lam = negativetwothird * mi_l
    Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave ! cylindrical divergence
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) - var(1)
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) - var(2)
    Tau_xy = mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) - var(4)
    Tau_zz = 2d0 * mi_l * vr / r_ave + lamDivervel - var(3)

    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz
    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy
    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_RSM_Axi


  pure function Vorticity_Vector ( Gradvel ) result ( omega )
    ! 2D-reduced: omega(1)=omega(2)=0; only omega(3) = G(2,1)-G(1,2).
    implicit none
    real(R8), intent(in), dimension(2,2) :: Gradvel
    real(R8), dimension(3) :: omega

    omega(1) = 0.d0
    omega(2) = 0.d0
    omega(3) = Gradvel(2,1) - Gradvel(1,2)
  end function Vorticity_Vector


  pure function Vorticity_Tensor ( omega ) result ( W )
    ! Convention: W_ij = (du_i/dx_j - du_j/dx_i)/2. Sign-flipped Hodge skew-symmetric;
    ! W is used only via |W|^2 = 2*sum(W_ij^2) in RANS models, sign convention immaterial.
    implicit none
    real(R8), intent(in), dimension(3) :: omega
    real(R8), dimension(3,3) :: W

    W(1,1) = 0.d0 ; W(2,2) = 0.d0 ; W(3,3) = 0.d0

    W(1,2) = -0.5d0*omega(3)
    W(2,1) =  0.5d0*omega(3)

    W(1,3) =  0.5d0*omega(2)
    W(3,1) = -0.5d0*omega(2)

    W(2,3) =  0.5d0*omega(1)
    W(3,2) = -0.5d0*omega(1)
  end function Vorticity_Tensor


  pure function Strain_Tensor_Std ( Gradvel, vr, r_ave, V_dot_G ) result ( E )
    ! Q2D Definition A (Model A, BC-consistent depth-averaged): S_33 = -V·G
    ! represents <∂w/∂z> from depth-averaging the 3D continuity equation under
    ! (H4) impermeability + slip at virtual top/bot walls. The formal identity
    ! <∂w/∂z> = -(V·G) follows from Leibniz + impermeability BC, without
    ! specifying any z-profile of w. When V_dot_G absent (default 0), reduces to
    ! Definition B (free-slip plug, no z-divergence correction).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8), dimension(3,3) :: E
    real(R8) :: VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    E(1,1) = Gradvel(1,1)
    E(2,2) = Gradvel(2,2)
    E(3,3) = -VG

    E(1,2) = 0.5d0*(Gradvel(1,2) + Gradvel(2,1))
    E(2,1) = E(1,2)

    E(3,1) = 0.d0
    E(1,3) = E(3,1)
    E(3,2) = 0.d0
    E(2,3) = E(3,2)
  end function Strain_Tensor_Std


  pure function Strain_Tensor_Axi ( Gradvel, vr, r_ave, V_dot_G ) result ( E )
    ! Axisymmetric variant: E(3,3) = vr/r_ave (theta-theta strain, with de l'Hopital
    ! fallback at axis). V_dot_G optional accepted for interface uniformity but NOT
    ! used in Axi.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8), dimension(3,3) :: E

    E(1,1) = Gradvel(1,1)
    E(2,2) = Gradvel(2,2)
    if (r_ave .le. 1.d-6) then
      E(3,3) = Gradvel(2,2) ! de l'Hopital
    else
      E(3,3) = vr/r_ave
    endif

    E(1,2) = 0.5d0*(Gradvel(1,2) + Gradvel(2,1))
    E(2,1) = E(1,2)

    E(3,1) = 0.d0
    E(1,3) = E(3,1)
    E(3,2) = 0.d0
    E(2,3) = E(3,2)
  end function Strain_Tensor_Axi


  pure function V_dot_G_face_avg ( V_loc, V_neigh, G_loc, G_neigh ) result ( VG )
    ! Canonical convention for V·G at the face (Q2D Def A).
    ! Two-average approximation: ⟨V⟩_face · ⟨G⟩_face, with ⟨V⟩=½(V_loc+V_neigh) and idem G.
    ! On a slip-symmetry BC (V_n_face=0 exactly), a cell-centered average would introduce
    ! a spurious O(1) term V_n_loc·G_n_loc; this convention gives exactly
    ! V_tan·G_tan (rigorous independently of h at the boundary).
    implicit none
    real(R8), intent(in) :: V_loc(2), V_neigh(2), G_loc(2), G_neigh(2)
    real(R8) :: VG
    real(R8) :: V_avg(2), G_avg(2)

    V_avg = 0.5d0 * ( V_loc + V_neigh )
    G_avg = 0.5d0 * ( G_loc + G_neigh )
    VG    = V_avg(1)*G_avg(1) + V_avg(2)*G_avg(2)
  end function V_dot_G_face_avg

end module Q2D_Lib_Fluid
