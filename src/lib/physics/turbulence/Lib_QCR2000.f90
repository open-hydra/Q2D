!>@brief Quadratic Constitutive Relations (Spalart 2000) — Q2D 2D-reduced.
!
! Reference: Spalart P.R. (2000), "Strategies for turbulence modelling and simulations",
!            Int. J. Heat Fluid Flow 21(3):252-263.
!
! 2D applicability caveat: QCR is motivated for 3D corner separation (junction
! flows wing-body, fan blade hub interactions). In Q2D (Tau_xz=Tau_yz=0 a priori
! because out-of-plane velocity gradients are absent), QCR acts only on the
! in-plane anisotropy (Tau_xx, Tau_yy, Tau_xy). Reduced physical payoff vs the
! target 3D use case. For 2D RANS engineering use SST/SA as the primary models;
! QCR is recommended as an optional correction for cases with significant
! anisotropy.
!
! Q2D provides 2 realizations (Stress_Vector_QCR2000_Std/Axi) following the
! Lib_Fluid Stress_Vector_Std/Axi axi-geometry switch pattern. The abstract
! Q2D Stress_Vector_if interface requires 8-args (Gradvel(2,2), Stress(2)),
! so the single 5-arg form (Gradvel(3,3), Normal(3), Stress(3)) cannot be
! assigned to the procedure pointer — hence the Std/Axi split.
!
! Dispatch: QCR2000 binds only in the SA branch; for SST/Wilcox the
! Stress_Vector_2eq path is used when k_energy_coupling is active, so QCR is
! mutually exclusive with it.
module Q2D_Lib_QCR2000
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  real(R8), parameter, private :: C_cr1             = 0.3d0
  real(R8), parameter, private :: negativetwothird  = -2.d0/3.d0

  public :: Stress_Vector_QCR2000_Std, Stress_Vector_QCR2000_Axi
  public :: Stress_Tensor_QCR2000_Std, Stress_Tensor_QCR2000_Axi

contains

  pure function Stress_Vector_QCR2000_Std ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! Signature 8-args matching Stress_Vector_if (Lib_Fluid).
    ! Q2D Def A: Divervel = G(1,1) + G(2,2) - VG (baked in via VG).
    ! 2D planar: Tau_xz_turb = Tau_yz_turb = 0 a priori
    ! (no out-of-plane velocity gradients) -> QCR correction on Tau_xx, Tau_yy, Tau_xy
    ! only, and Tau_zz_QCR = Tau_zz (no out-of-plane mixing).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t                ! turbulent (Boussinesq deviatoric)
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR          ! QCR-corrected turbulent
    real(R8) :: Tau_xx, Tau_yy, Tau_xy                      ! laminar + turbulent total
    real(R8) :: lam_l, Divervel, lamDivervel, VG
    real(R8) :: omega_z, Den, O_off

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    ! Q2D Def A divergence
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG

    ! Turbulent stress (Boussinesq deviatoric, with mi_t). 2D planar.
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! Vorticity z-component (omega(1)=omega(2)=0)
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den = sqrt(sum(Gradvel_full**2)) includes the (3,3)=-VG component.
    ! Avoid division by zero in laminar regions.
    Den = sqrt ( sum ( Gradvel**2 ) + VG**2 ) + 1d-40

    ! O_ij(1,2) = 2*W_ij(1,2)/Den = -omega_z/Den; O_ij(2,1) = +omega_z/Den.
    ! Off-diagonal magnitude in z-component sense (only nonzero entries in 2D).
    O_off = omega_z / Den   ! = -O_ij(1,2) = +O_ij(2,1)

    ! QCR correction (reduces 6 components to in-plane 3):
    ! T_ij_QCR = T_ij - C_cr1*[O_ik*T_jk + O_jk*T_ik], with sums on k = {1,2,3}.
    ! In 2D planar: O_ik = 0 except (1,2),(2,1); T_ik = 0 for k=3 turbulent stress.
    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1   ! = Tau_xx_t + 2*O_off*Tau_xy_t*C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1
                                                                ! = Tau_xy_t + O_off*(Tau_yy_t - Tau_xx_t)*C_cr1

    ! Laminar (linear Newtonian, no QCR correction). Mirror Stress_Vector_Std.
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_QCR2000_Std


  pure function Stress_Vector_QCR2000_Axi ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! Signature 8-args matching Stress_Vector_if. Cylindrical divergence.
    ! V_dot_G optional argument for interface uniformity (NOT used in Axi:
    ! axisymmetric flows have no virtual top/bot wall in Def A sense).
    ! 2D-axi No-w (no swirl): Tau_xz_turb = Tau_yz_turb = 0, and
    ! Tau_zz_turb (theta-theta) does NOT enter in-plane Stress(1:2).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR
    real(R8) :: Tau_xx, Tau_yy, Tau_xy
    real(R8) :: lam_l, Divervel, lamDivervel, G_zz
    real(R8) :: omega_z, Den, O_off

    ! Cylindrical coordinate divergence (de l'Hopital fallback at axis r_ave -> 0)
    if (r_ave .le. 1.d-6) then
      Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(2,2)  ! de l'Hopital vr/r --> dv/dr
      G_zz     = Gradvel(2,2)
    else
      Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave
      G_zz     = vr/r_ave
    end if

    ! Turbulent stress
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! Vorticity z-component
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den includes the (3,3) = vr/r_ave component (axi)
    Den = sqrt ( sum ( Gradvel**2 ) + G_zz**2 ) + 1d-40
    O_off = omega_z / Den

    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1

    ! Laminar
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_QCR2000_Axi


  pure function Stress_Tensor_QCR2000_Std ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! 3x3 stress tensor with QCR2000 correction (Spalart 2000), Std variant.
    ! Verified by the term_by_term_sympy audit (multi-method).
    ! Mirror Stress_Vector_QCR2000_Std pattern with:
    !   - 7-args signature matching Stress_Tensor_if (no Normal, output 3x3 instead of vector(2)·n).
    !   - Q2D Def A: Divervel = G(1,1) + G(2,2) - VG (baked via VG).
    !   - QCR correction on in-plane (xx, yy, xy); Tau_QCR_zz = Tau_zz UNCHANGED
    !     (No-w → O_θ·=0 → no QCR on (3,3) component).
    !   - Tau_xz, Tau_yz = 0 (No-w).
    ! Sanity check: O_off=0 → Tau_QCR_xx/yy/xy = Tau_xx/yy/xy non-QCR ✓.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t          ! turbulent (Boussinesq deviatoric, mi_t)
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR    ! QCR-corrected turbulent
    real(R8) :: Tau_xx, Tau_yy, Tau_xy, Tau_zz        ! full (laminar + turbulent_QCR)
    real(R8) :: lam_l, Divervel, lamDivervel, VG
    real(R8) :: omega_z, Den, O_off
    real(R8) :: mi_combined

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    ! Q2D Def A divergence (baked via VG)
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG

    ! Turbulent stress (Boussinesq deviatoric, with mi_t). 2D planar.
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! Vorticity z-component (omega(1)=omega(2)=0)
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den includes the (3,3) = -VG component (area-var Definition A).
    Den = sqrt ( sum ( Gradvel**2 ) + VG**2 ) + 1d-40
    O_off = omega_z / Den   ! = -O_ij(1,2) = +O_ij(2,1)

    ! QCR correction (reduces 6 components to in-plane 3):
    ! T_ij_QCR = T_ij - C_cr1*[O_ik*T_jk + O_jk*T_ik], with sums on k.
    ! O_θ·=0 (No-w) → no QCR correction on (3,3) component.
    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1

    ! Laminar (linear Newtonian, no QCR correction).
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    ! In-plane: laminar + QCR-corrected turbulent
    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    ! Tau_zz: V·G baked + NO QCR correction (No-w → O_θ·=0).
    ! Mirror Stress_Tensor_Std: Tau_zz = lamDivervel_combined + 2*mi*(-VG).
    mi_combined = mi_l + mi_t
    Tau_zz = (negativetwothird * mi_combined) * Divervel + 2d0 * mi_combined * (-VG)

    ! Assemble 3x3 tensor (τ_xz = τ_yz = 0)
    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz

    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy

    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_QCR2000_Std


  pure function Stress_Tensor_QCR2000_Axi ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! 3x3 stress tensor with QCR2000 correction, Axi variant.
    ! Mirror Stress_Tensor_Axi (Lib_Fluid) + QCR on in-plane (xx, yy, xy);
    ! Tau_QCR_θθ = Tau_θθ UNCHANGED (No-w → O_θ·=0).
    ! V_dot_G optional argument for interface uniformity (NOT used in Axi).
    ! Tensor evaluated cell-center (r_ave > 0 by mesh construction) → no de l'Hopital fallback
    ! (same convention as Stress_Tensor_Axi).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR
    real(R8) :: Tau_xx, Tau_yy, Tau_xy, Tau_zz
    real(R8) :: lam_l, Divervel, lamDivervel, G_zz
    real(R8) :: omega_z, Den, O_off
    real(R8) :: mi_combined

    ! Cylindrical coordinate divergence (V_r/r in (3,3) component).
    Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave
    G_zz     = vr/r_ave

    ! Turbulent stress (Boussinesq deviatoric, with mi_t).
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! Vorticity z-component
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den includes the (3,3) = vr/r_ave component (cylindrical Christoffel).
    Den = sqrt ( sum ( Gradvel**2 ) + G_zz**2 ) + 1d-40
    O_off = omega_z / Den

    ! QCR correction (in-plane only, O_θ·=0)
    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1

    ! Laminar
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    ! Tau_zz (= τ_θθ in cylindrical reading): vr/r baked + NO QCR (No-w).
    ! Mirror Stress_Tensor_Axi: Tau_zz = 2*mi*vr/r_ave + lamDivervel_combined.
    mi_combined = mi_l + mi_t
    Tau_zz = (negativetwothird * mi_combined) * Divervel + 2d0 * mi_combined * G_zz

    ! Assemble 3x3 tensor (τ_zr/rz off-diagonals zero)
    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz

    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy

    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_QCR2000_Axi


end module Q2D_Lib_QCR2000
