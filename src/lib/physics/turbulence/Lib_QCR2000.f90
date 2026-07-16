!>@brief Quadratic Constitutive Relations (Spalart 2000) — 2D reduction for Q2D.
!
! Reference: Spalart P.R. (2000), "Strategies for turbulence modelling and simulations",
!            Int. J. Heat Fluid Flow 21(3):252-263.
!
! Applicability limit in 2D: QCR was devised for 3D corner separations
! (wing-body junctions, hub interactions of blade rows). In Q2D
! (Tau_xz=Tau_yz=0 a priori because the out-of-plane velocity gradients are
! absent) QCR acts only on the in-plane anisotropy (Tau_xx, Tau_yy,
! Tau_xy): the physical benefit is therefore reduced with respect to the
! original 3D use case. For engineering 2D RANS the primary models remain SST/SA;
! QCR is recommended as an optional correction for cases with significant
! anisotropy.
!
! Q2D provides 2 realizations (Stress_Vector_QCR2000_Std/Axi) following the
! Std/Axi geometric switching scheme of Stress_Vector_Std/Axi in
! Lib_Fluid. The Q2D abstract interface Stress_Vector_if requires 8
! arguments (Gradvel(2,2), Stress(2)), so a possible single 5-argument
! form (Gradvel(3,3), Normal(3), Stress(3)) is not assignable to the
! procedure pointer — hence the Std/Axi separation.
!
! Dispatch: QCR2000 hooks in only on the SA branch; for SST/Wilcox, when
! k_energy_coupling is active, the Stress_Vector_2eq path is used, so
! QCR is mutually exclusive with it.
module Q2D_Lib_QCR2000
  use iso_fortran_env,    only: I4 => int32, R8 => real64
  use Q2D_Config_Types_m, only: obj_rot   ! accessed read-only in the pure functions below

  implicit none
  private

  real(R8), parameter, private :: C_cr1             = 0.3d0
  real(R8), parameter, private :: negativetwothird  = -2.d0/3.d0

  public :: Stress_Vector_QCR2000_Std, Stress_Vector_QCR2000_Axi
  public :: Stress_Tensor_QCR2000_Std, Stress_Tensor_QCR2000_Axi

contains

  pure function Stress_Vector_QCR2000_Std ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! 8-argument signature conforming to Stress_Vector_if (Lib_Fluid).
    ! Q2D Definition A: Divervel = G(1,1) + G(2,2) - VG (baked in via VG).
    ! 2D plane: Tau_xz_turb = Tau_yz_turb = 0 a priori
    ! (no out-of-plane velocity gradient) -> QCR correction only on
    ! Tau_xx, Tau_yy, Tau_xy, and Tau_zz_QCR = Tau_zz (no out-of-plane mixing).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t                ! turbulent (Boussinesq deviatoric)
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR          ! QCR-corrected turbulent
    real(R8) :: Tau_xx, Tau_yy, Tau_xy                      ! total laminar + turbulent
    real(R8) :: lam_l, Divervel, lamDivervel, VG
    real(R8) :: omega_z, Den, O_off

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    ! Divergence according to Q2D Definition A
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG

    ! Turbulent stress (Boussinesq deviatoric, with mi_t). 2D plane.
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! z-component of vorticity (omega(1)=omega(2)=0)
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den = sqrt(sum(Gradvel_full**2)) includes the (3,3)=-VG component.
    ! Avoids division by zero in laminar regions.
    Den = sqrt ( sum ( Gradvel**2 ) + VG**2 ) + 1d-40
    ! ── Objectivity in a rotating frame (verified on an analytic vortex,
    ! fixed axes vs rotating axes): O_ik = 2 W_ik/||grad u|| must be built with
    ! quantities in ABSOLUTE axes — omega_z shifted by 2*Oz (like every other
    ! vorticity consumer) AND Den from the absolute gradient = relative +
    ! entrainment gradient d(Omega x r)_a/dx_b = [[0,-Oz],[Oz,0]] in the
    ! plane (the (3,3) component is frame-invariant). Shifting omega alone
    ! is NOT enough: a Den from the relative gradient leaves tau_QCR
    ! frame-dependent and over-corrects; with both corrections
    ! tau_QCR is exactly objective (machine zero).
    if (obj_rot%enabled) then
      omega_z = omega_z + 2.0_R8*obj_rot%Omega_vec(3)
      Den = sqrt ( Gradvel(1,1)**2 + (Gradvel(1,2)-obj_rot%Omega_vec(3))**2 &
                 + (Gradvel(2,1)+obj_rot%Omega_vec(3))**2 + Gradvel(2,2)**2 &
                 + VG**2 ) + 1d-40
    end if

    ! O_ij(1,2) = 2*W_ij(1,2)/Den = -omega_z/Den; O_ij(2,1) = +omega_z/Den.
    ! Off-diagonal magnitude in the z-component sense (the only nonzero entries in 2D).
    O_off = omega_z / Den   ! = -O_ij(1,2) = +O_ij(2,1)

    ! QCR correction (the 6 components reduce to the 3 in-plane ones):
    ! T_ij_QCR = T_ij - C_cr1*[O_ik*T_jk + O_jk*T_ik], sums over k = {1,2,3}.
    ! In the 2D plane: O_ik = 0 except (1,2),(2,1); T_ik = 0 for k=3 in the turbulent stress.
    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1   ! = Tau_xx_t + 2*O_off*Tau_xy_t*C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1
                                                                ! = Tau_xy_t + O_off*(Tau_yy_t - Tau_xx_t)*C_cr1

    ! Laminar (linear Newtonian, without QCR correction). Mirrors Stress_Vector_Std.
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_QCR2000_Std


  pure function Stress_Vector_QCR2000_Axi ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! 8-argument signature conforming to Stress_Vector_if. Cylindrical divergence.
    ! Optional argument V_dot_G for interface uniformity (NOT used in
    ! Axi: axisymmetric flows have no top/bot virtual wall in the
    ! sense of Definition A).
    ! 2D-axi No-w (without swirl): Tau_xz_turb = Tau_yz_turb = 0, and
    ! Tau_zz_turb (theta-theta) does NOT enter Stress(1:2) in the plane.
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

    ! Divergence in cylindrical coordinates (de l'Hopital fallback at the axis r_ave -> 0)
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

    ! z-component of vorticity
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den includes the (3,3) = vr/r_ave component (axi)
    Den = sqrt ( sum ( Gradvel**2 ) + G_zz**2 ) + 1d-40
    ! ── Objectivity in a rotating frame — see the comment of the Std variant
    ! (objectivity verified on an analytic vortex): absolute omega AND
    ! absolute gradient in Den; the cylindrical (3,3) = vr/r component is
    ! frame-invariant. (The 2D-axi + rotating-frame combination is currently
    ! blocked by the compatibility check in Setup_RotatingFrame; this
    ! branch keeps the Axi realization frame-consistent
    ! in view of the 2D-axi extension with swirl.)
    if (obj_rot%enabled) then
      omega_z = omega_z + 2.0_R8*obj_rot%Omega_vec(3)
      Den = sqrt ( Gradvel(1,1)**2 + (Gradvel(1,2)-obj_rot%Omega_vec(3))**2 &
                 + (Gradvel(2,1)+obj_rot%Omega_vec(3))**2 + Gradvel(2,2)**2 &
                 + G_zz**2 ) + 1d-40
    end if
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
    ! Mirrors the scheme of Stress_Vector_QCR2000_Std with:
    !   - 7-argument signature conforming to Stress_Tensor_if (no Normal, 3x3 output
    !     instead of vector(2)·n).
    !   - Q2D Definition A: Divervel = G(1,1) + G(2,2) - VG (baked in via VG).
    !   - QCR correction on the in-plane components (xx, yy, xy); Tau_QCR_zz = Tau_zz
    !     UNCHANGED (No-w → O_θ·=0 → no QCR on the (3,3) component).
    !   - Tau_xz, Tau_yz = 0 (No-w).
    ! Consistency check: O_off=0 → Tau_QCR_xx/yy/xy = Tau_xx/yy/xy non-QCR ✓.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t          ! turbulent (Boussinesq deviatoric, mi_t)
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR    ! QCR-corrected turbulent
    real(R8) :: Tau_xx, Tau_yy, Tau_xy, Tau_zz        ! complete (laminar + turbulent_QCR)
    real(R8) :: lam_l, Divervel, lamDivervel, VG
    real(R8) :: omega_z, Den, O_off
    real(R8) :: mi_combined

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    ! Divergence according to Q2D Definition A (baked in via VG)
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG

    ! Turbulent stress (Boussinesq deviatoric, with mi_t). 2D plane.
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! z-component of vorticity (omega(1)=omega(2)=0)
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den includes the (3,3) = -VG component (Definition A with area variation).
    Den = sqrt ( sum ( Gradvel**2 ) + VG**2 ) + 1d-40
    ! ── Objectivity in a rotating frame (verified on an analytic vortex,
    ! fixed axes vs rotating axes): O_ik = 2 W_ik/||grad u|| must be built with
    ! quantities in ABSOLUTE axes — omega_z shifted by 2*Oz (like every other
    ! vorticity consumer) AND Den from the absolute gradient = relative +
    ! entrainment gradient d(Omega x r)_a/dx_b = [[0,-Oz],[Oz,0]] in the
    ! plane (the (3,3) component is frame-invariant). Shifting omega alone
    ! is NOT enough: a Den from the relative gradient leaves tau_QCR
    ! frame-dependent and over-corrects; with both corrections
    ! tau_QCR is exactly objective (machine zero).
    if (obj_rot%enabled) then
      omega_z = omega_z + 2.0_R8*obj_rot%Omega_vec(3)
      Den = sqrt ( Gradvel(1,1)**2 + (Gradvel(1,2)-obj_rot%Omega_vec(3))**2 &
                 + (Gradvel(2,1)+obj_rot%Omega_vec(3))**2 + Gradvel(2,2)**2 &
                 + VG**2 ) + 1d-40
    end if
    O_off = omega_z / Den   ! = -O_ij(1,2) = +O_ij(2,1)

    ! QCR correction (the 6 components reduce to the 3 in-plane ones):
    ! T_ij_QCR = T_ij - C_cr1*[O_ik*T_jk + O_jk*T_ik], sums over k.
    ! O_θ·=0 (No-w) → no QCR correction on the (3,3) component.
    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1

    ! Laminar (linear Newtonian, without QCR correction).
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    ! In the plane: laminar + QCR-corrected turbulent
    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    ! Tau_zz: V·G baked in + NO QCR correction (No-w → O_θ·=0).
    ! Mirrors Stress_Tensor_Std: Tau_zz = lamDivervel_combined + 2*mi*(-VG).
    mi_combined = mi_l + mi_t
    Tau_zz = (negativetwothird * mi_combined) * Divervel + 2d0 * mi_combined * (-VG)

    ! Assembly of the 3x3 tensor (τ_xz = τ_yz = 0)
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
    ! Mirrors Stress_Tensor_Axi (Lib_Fluid) + QCR on the in-plane components
    ! (xx, yy, xy); Tau_QCR_θθ = Tau_θθ UNCHANGED (No-w → O_θ·=0).
    ! Optional argument V_dot_G for interface uniformity (NOT used in Axi).
    ! Tensor evaluated at cell center (r_ave > 0 by mesh construction) →
    ! no de l'Hopital fallback (same convention as Stress_Tensor_Axi).
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

    ! Divergence in cylindrical coordinates (V_r/r in the (3,3) component).
    Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave
    G_zz     = vr/r_ave

    ! Turbulent stress (Boussinesq deviatoric, with mi_t).
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! z-component of vorticity
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den includes the (3,3) = vr/r_ave component (cylindrical Christoffel term).
    Den = sqrt ( sum ( Gradvel**2 ) + G_zz**2 ) + 1d-40
    ! ── Objectivity in a rotating frame — see the comment of the Std variant
    ! (objectivity verified on an analytic vortex): absolute omega AND
    ! absolute gradient in Den; the cylindrical (3,3) = vr/r component is
    ! frame-invariant. (The 2D-axi + rotating-frame combination is currently
    ! blocked by the compatibility check in Setup_RotatingFrame; this
    ! branch keeps the Axi realization frame-consistent
    ! in view of the 2D-axi extension with swirl.)
    if (obj_rot%enabled) then
      omega_z = omega_z + 2.0_R8*obj_rot%Omega_vec(3)
      Den = sqrt ( Gradvel(1,1)**2 + (Gradvel(1,2)-obj_rot%Omega_vec(3))**2 &
                 + (Gradvel(2,1)+obj_rot%Omega_vec(3))**2 + Gradvel(2,2)**2 &
                 + G_zz**2 ) + 1d-40
    end if
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

    ! Tau_zz (= τ_θθ in cylindrical reading): vr/r baked in + NO QCR (No-w).
    ! Mirrors Stress_Tensor_Axi: Tau_zz = 2*mi*vr/r_ave + lamDivervel_combined.
    mi_combined = mi_l + mi_t
    Tau_zz = (negativetwothird * mi_combined) * Divervel + 2d0 * mi_combined * G_zz

    ! Assembly of the 3x3 tensor (off-diagonal τ_zr/rz null)
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
