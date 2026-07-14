!>@brief: Kappa-Omega SST (2003) Menter. https://turbmodels.larc.nasa.gov/sst.html
!
! 12 procedures plus 1 Q2D-extra procedure
! (SST_Source_Diffusive_Q2D). 2D loops, 9/12-arg signatures with the
! extra vr/r_ave (axisymmetric) and V_dot_G (Definition A) arguments.
!
! When obj_rans%SpalartShur is set, SST_Source_Terms invokes
! Compute_Velocity_Gradient + Compute_RC_Terms and SST_Blk computes
! fr1 = max(min(frot,1.25),0) with frot=(1+cr1)*(2*rstar/(1+rstar))*
! (1-cr3*atan(cr2*rtilde))-cr1 and the SST-RC denominator D=sqrt(max(S²,0.09*omega²))
! obj_rot%enabled wires the absolute-vorticity correction in the
! SST_Blk body only (NOT SST_Eddy_Viscosity).
module Q2D_Lib_SST
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  real(R8), parameter, private :: &
    sigma_k1 = 0.85d0, &
    sigma_k2 = 1.0d0, &
    sigma_w1 = 0.5d0, &
    sigma_w2 = 0.856d0, &
    beta_1 = 0.075d0, &
    beta_2 = 0.0828d0, &
    beta_star = 0.09d0, &
    karman = 0.41d0, &
    alpha_1 = 0.31d0, &
    gamma_1 = 5d0/9d0, &
    gamma_2 = 0.44d0, &
    cr1 = 1.d0, &
    cr2 = 2.d0, &
    cr3 = 1.d0

  public :: SST_Source_Terms
  public :: SST_Set_Wall_Values, SST_Blowing_Correction, SST_Blowing_noCorrection
  public :: SST_Extrapolate_Wall, SST_Enforce_Realizability
  public :: SST_Eddy_Viscosity, SST_RANS_Diffusive_Flux
  public :: SST_Source_Diffusive_Q2D

contains

  subroutine SST_Source_Terms ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Lib_SpalartShur,  only: Compute_Velocity_Gradient, Compute_RC_Terms
    use Q2D_Config_Types_m,   only: obj_rans
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: b

    ! Spalart-Shur rotation/curvature precomputation.
    ! Compute_Velocity_Gradient + Compute_RC_Terms populate vel_gradient +
    ! rc_term1 + rc_term2 BEFORE the block loop. Conditional allocation in
    ! Mod_Allocate_Data guarantees the arrays exist iff SpalartShur.
    if ( obj_rans%SpalartShur ) then
      call Compute_Velocity_Gradient ( domain )
      call Compute_RC_Terms ( domain )
    end if

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks

      call SST_Blk ( domain % blk(b) % P,                 &
                     domain % blk(b) % r,                 &
                     domain % blk(b) % M,                 &
                     domain % blk(b) % sur,               &
                     domain % blk(b) % yn,                &
                     domain % blk(b) % vel_gradient,      &
                     domain % blk(b) % rc_term1,          &
                     domain % blk(b) % rc_term2,          &
                     domain % blk(b) % node,              &
                     domain % blk(b) % dim,               &
                     domain % blk(b) % area_geom_term_x,  &
                     domain % blk(b) % area_geom_term_y   )

    end do

  end subroutine SST_Source_Terms


  subroutine SST_Blk ( Prim, Res, M, sur, WDist, gradv, rc1, rc2, r, n, geom_x, geom_y )
    ! Signature 12-args. sur is the 2D analogue of the 3D cell volume.
    ! Body 2D loop; Q2D Def A V_dot_G at cell centre for the
    ! P_k = mu_t·[|S|^2_2D + 2(V·G)^2] correction (v7 §N16).
    ! k_energy_coupling read from obj_rans (vs the MOSE peer's arg).
    use Q2D_Base_Types_m,        only: Q2D_tensor_2D_type, Q2D_vector_2D_type
    use Q2D_Global_m,            only: nprim, nsc, nu, nv, np, nt, gc
    use Q2D_Config_Types_m,      only: obj_rans, obj_rot
    use Q2D_Lib_Fluid,           only: Vorticity_Vector, Vorticity_Tensor, Strain_Tensor
    use FLINT_Lib_Thermodynamic, only: Co_rotot_Rtot, f_laminarViscosity

    implicit none
    integer, intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: sur
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: WDist
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: M
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: gradv
    real(R8), dimension(n(1), n(2)), intent(in) :: rc1, rc2
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)), intent(in) :: r
    real(R8), dimension(1:n(1), 1:n(2)), intent(in) :: geom_x, geom_y
    integer :: i, j, ii, jj
    real(R8) :: rho, Rgas, mil, kap, ome, dist, Gradvel(2,2), Divel, Sij(3,3), S, Wij(3,3), diag, mi_t, F(2)
    real(R8) :: vort(3), O, rstar, D, rtilde, frot, fr1
    real(R8) :: Tij(2,2), Prod(2), Grad(2,2), dkDotdw, Diff, beta, gamma, Diss(2), Source(2)
    real(R8) :: r_ave, V_dot_G

    !$omp do collapse (2) private ( rho, Rgas, mil, kap, ome, dist, Gradvel, Divel, Sij, S, Wij, diag, mi_t, F ), &
    !$omp private ( Tij, Prod, Grad, dkDotdw, Diff, beta, gamma, Diss, Source, r_ave, V_dot_G, i, j, ii, jj ), &
    !$omp private ( vort, O, rstar, D, rtilde, frot, fr1 )

    do j = 1, n(2)
    do i = 1, n(1)

      r_ave = 0.25d0 * (r(i  ,j  )%c(2) + r(i-1  ,j  )%c(2) + r(i  ,j-1  )%c(2) + r(i-1  ,j-1  )%c(2))
      call Co_rotot_Rtot ( Prim(1:nsc,i,j), rho, Rgas )
      mil = f_laminarViscosity ( Prim(1:nsc,i,j), Prim(np,i,j), rho, Rgas )
      kap = Prim(nt,i,j) / rho
      ome = Prim(nt+1,i,j) / rho
      dist = WDist(i,j)

      Gradvel(:,1) = ( Prim(nu:nv,i+1,j) - Prim(nu:nv,i-1,j) ) * 0.5d0
      Gradvel(:,2) = ( Prim(nu:nv,i,j+1) - Prim(nu:nv,i,j-1) ) * 0.5d0
      Gradvel = Matmul ( Gradvel, M(i,j) % c )
      Divel = Gradvel(1,1) + Gradvel(2,2)

      ! Q2D Def A: V·G at cell centre for P_k = mu_t·[|S|^2_2D + 2(V·G)^2] (v7 §N16).
      ! Backward-compat: with area_geom_term=0 (universal zero-fill) -> V_dot_G=0 -> classic SST.
      V_dot_G = Prim(nu,i,j)*geom_x(i,j) + Prim(nv,i,j)*geom_y(i,j)

      Grad(:,1) = ( Prim(nt:nt+1,i+1,j) / sum( Prim(1:nsc,i+1,j) ) - &
                    Prim(nt:nt+1,i-1,j) / sum( Prim(1:nsc,i-1,j) ) ) * 0.5d0
      Grad(:,2) = ( Prim(nt:nt+1,i,j+1) / sum( Prim(1:nsc,i,j+1) ) - &
                    Prim(nt:nt+1,i,j-1) / sum( Prim(1:nsc,i,j-1) ) ) * 0.5d0
      Grad = Matmul ( Grad, M(i,j) % c )

      F = compute_F ( rho=rho, w=ome, k=kap, mil=mil, gradk=Grad(1,:), gradw=Grad(2,:), d=dist )
      gamma = gamma_1*F(1) + gamma_2*(1d0-F(1))
      beta = beta_1*F(1) + beta_2*(1d0-F(1))

      ! Q2D Def A: Strain_Tensor with V_dot_G -> S_33 = -V·G
      Sij = Strain_Tensor ( Gradvel, Prim(nv,i,j), r_ave, V_dot_G=V_dot_G )
      S = Sqrt ( 2d0 * sum ( Sij**2 ) )

      vort = Vorticity_Vector ( Gradvel )
      ! Absolute vorticity omega_abs =
      ! omega_rel + 2*Omega_frame, needed for O in the SST-RC denominator
      ! and for Wij passed to the Spalart-Shur fr1. Q2D Vorticity_Vector
      ! returns vector(3) (omega(1)=omega(2)=0, omega(3)=Gradvel(2,1)
      ! -Gradvel(1,2)); Omega_vec axis [0,0,1] strict.
      if (obj_rot%enabled) vort = vort + 2.0_R8 * obj_rot%Omega_vec
      Wij = Vorticity_Tensor ( vort )
      O = Sqrt ( 2d0 * sum ( Wij**2 ) )
      mi_t = rho*kap*alpha_1 / Max ( alpha_1*ome, S*F(2) )

      Prod(1) = mi_t * S**2 ! SST-2003m

      if ( obj_rans%SpalartShur ) then
        ! SST-RC (Spalart & Shur 1997): rotation/curvature correction.
        ! fr1 multiplies production; rstar=S/Omega from
        ! Compute_RC_Terms (cylindrical+area-var formal handling); rc_term2 is the
        ! numerator of rtilde (SST-specific denominator applied here: O*D^3 with D=max(S²,0.09*ω²)^½).
        rstar  = rc1(i,j)
        D      = sqrt( max( S**2, 0.09d0*ome**2) )
        rtilde = rc2(i,j) / ( O*D**3 )
        frot = (1.d0 + cr1)*(2.d0*rstar/(1.d0 + rstar))*(1.d0 - cr3*atan(cr2*rtilde)) - cr1
        fr1  = max( min(frot, 1.25d0), 0d0 )    ! limit fr1 in [0, 1.25]
        Prod(1) = fr1*Prod(1)
      end if

      Prod(1) = Min( Prod(1), 10d0*beta_star*rho*ome*kap )
      Prod(2) = ( gamma*rho / mi_t ) * Prod(1)

      dkDotdw = dot_product ( Grad(1,:), Grad(2,:) )
      Diff = 2d0*(1d0-F(1))*rho*sigma_w2/ome*dkDotdw

      Diss(1) = beta_star*rho*kap*ome
      Diss(2) = beta*rho*ome**2

      Source(1) = Prod(1) - Diss(1)
      Source(2) = Prod(2) - Diss(2) + Diff

      Res(nt:nt+1,i,j) = Res(nt:nt+1,i,j) - Source * sur(i,j)
      if (obj_rans%k_energy_coupling) Res(np,i,j) = Res(np,i,j) + Source(1) * sur(i,j)

    enddo ; enddo

  end subroutine SST_Blk


  subroutine SST_Eddy_Viscosity ( mi_t, rkw, mi_l, rho, Gradvel, vr, r_ave, dist, V_dot_G )
    ! Signature 9-args Q2D-form. Body uses Strain_Tensor with V_dot_G ->
    ! S_33 = -V·G modifies |S| and hence mi_t (Def A baked into the body).
    use Q2D_Global_m,  only: nRANS
    use Q2D_Lib_Fluid, only: Vorticity_Vector, Vorticity_Tensor, Strain_Tensor

    implicit none
    real(R8), intent(in), dimension(nRANS) :: rkw
    real(R8), intent(in)                   :: mi_l
    real(R8), intent(in)                   :: rho
    real(R8), intent(in)                   :: dist
    real(R8), intent(in), dimension(2,2)   :: Gradvel
    real(R8), intent(in)                   :: vr, r_ave
    real(R8), intent(out)                  :: mi_t
    real(R8), intent(in), optional         :: V_dot_G
    real(R8) :: kap, ome, Sij(3,3), S, F2, vort(3), Wij(3,3), O, VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    kap = rkw(1) / rho
    ome = rkw(2) / rho
    F2 = compute_F2 ( rho=rho, w=ome, k=kap, mil=mi_l, d=dist )
    Sij = Strain_Tensor ( Gradvel, vr, r_ave, V_dot_G=VG )
    S = Sqrt ( 2d0 * sum ( Sij**2 ) )
    vort = Vorticity_Vector ( Gradvel )
    Wij = Vorticity_Tensor ( vort )
    O = Sqrt ( 2d0 * sum ( Wij**2 ) )
    mi_t = rho*kap*alpha_1 / Max ( alpha_1*ome, S*F2 )

  end subroutine SST_Eddy_Viscosity


  subroutine SST_RANS_Diffusive_Flux ( Flux, rkw, mi_l, Gradvel, Gradkw, rho, int_l, Normal, vr, r_ave, Dist, mut )
    ! Signature 12-args. int_l is the 2D analogue of the 3D face area.
    ! Branch `if (present(mut))` Q2D-extra: reuse the caller's mu_t for consistency
    ! with the momentum-stress mu_t under area variation (V_dot_G modifies |S| and mu_t).
    ! The caller omits mut only for BC_Wall_Heat/Wall_Temperature where V=0 -> S_33=0 exactly.
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS)   :: rkw
    real(R8), intent(in)                     :: mi_l
    real(R8), intent(in), dimension(2,2)     :: Gradvel
    real(R8), intent(in), dimension(nRANS,2) :: Gradkw
    real(R8), intent(in)                     :: rho
    real(R8), intent(in)                     :: int_l, Normal(2), Dist
    real(R8), intent(in)                     :: vr, r_ave
    real(R8), intent(out), dimension(nRANS)  :: Flux
    real(R8), intent(in), optional           :: mut
    real(R8) :: kap, ome, F1, sigma_k, sigma_w, mi_t

    kap = rkw(1)/rho
    ome = rkw(2)/rho
    F1 = compute_F1 (rho=rho, w=ome, k=kap, mil=mi_l, &
                     gradk=Gradkw(1,:), gradw=Gradkw(2,:), d=Dist)
    sigma_k = sigma_k1*F1 + sigma_k2*(1d0-F1)
    sigma_w = sigma_w1*F1 + sigma_w2*(1d0-F1)
    if (present(mut)) then
      mi_t = mut
    else
      call SST_Eddy_Viscosity( mi_t=mi_t, rkw=rkw, mi_l=mi_l, rho=rho, &
                               Gradvel=Gradvel, vr=vr, r_ave=r_ave, dist=Dist )
    end if
    Flux(1) = ( mi_l + sigma_k*mi_t ) * dot_product ( Gradkw(1,:), Normal )
    Flux(2) = ( mi_l + sigma_w*mi_t ) * dot_product ( Gradkw(2,:), Normal )
    Flux = Flux * int_l

  end subroutine SST_RANS_Diffusive_Flux


  subroutine SST_Source_Diffusive_Q2D (source, rkw, mi_l, mut, gradvel, gradkw, rho, dist)
    ! Q2D-extra: SST diffusion projected onto G in Source_Area_Visc_Blk.
    ! gradvel accepted for interface conformance.
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS)    :: rkw
    real(R8), intent(in)                      :: mi_l
    real(R8), intent(in)                      :: mut
    real(R8), intent(in), dimension(2,2)      :: gradvel
    real(R8), intent(in), dimension(nRANS,2)  :: gradkw
    real(R8), intent(in)                      :: rho
    real(R8), intent(in)                      :: dist
    real(R8), intent(out), dimension(nRANS,2) :: source
    real(R8) :: kap, ome, F1, sigma_k, sigma_w

    kap = rkw(1)/rho
    ome = rkw(2)/rho
    F1 = compute_F1(rho=rho, w=ome, k=kap, mil=mi_l, &
                    gradk=gradkw(1,:), gradw=gradkw(2,:), d=dist)
    sigma_k = sigma_k1*F1 + sigma_k2*(1d0-F1)
    sigma_w = sigma_w1*F1 + sigma_w2*(1d0-F1)
    source(1,1) = ( mi_l + sigma_k*mut ) * gradkw(1,1)
    source(1,2) = ( mi_l + sigma_k*mut ) * gradkw(1,2)
    source(2,1) = ( mi_l + sigma_w*mut ) * gradkw(2,1)
    source(2,2) = ( mi_l + sigma_w*mut ) * gradkw(2,2)

  end subroutine SST_Source_Diffusive_Q2D


  subroutine SST_Set_Wall_Values ( mi_l, rkw, dist )
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: rkw
    real(R8), intent(in)                     :: mi_l
    real(R8), intent(in)                     :: dist

    rkw(1) = 0d0
    rkw(2) = 8d2 * mi_l / ( dist**2 ) ! Menter kw-SST smooth surface BC

  end subroutine SST_Set_Wall_Values


  subroutine SST_Blowing_Correction ( rho, mil, tau, rkw, mdot, dist )
    ! tau(2) (2D-reduced); the MOSE peer has tau(3).
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: rkw
    real(R8), intent(in)                     :: rho, mil, mdot, dist
    real(R8), dimension(2), intent(in)       :: tau
    real(R8) :: tauw, uT, v_w, v_plus, SB

    rkw(1) = 0d0
    tauw = sqrt( sum(tau**2) )
    uT = sqrt ( tauw / rho )
    v_w = abs(mdot)/rho
    v_plus = v_w / uT
    SB = 25d0 / ( v_plus*( 1d0 + 5d0 * v_plus ) )
    rkw(2) = (rho*uT)**2 * SB / mil ! Wilcox BC for mass injection

  end subroutine SST_Blowing_Correction


  subroutine SST_Blowing_noCorrection ( rho, mil, tau, rkw, mdot, dist )
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: rkw
    real(R8), intent(in)                     :: rho, mil, mdot, dist
    real(R8), dimension(2), intent(in)       :: tau

    rkw(1) = 0d0
    rkw(2) = 8d2 * mil / ( dist**2 )

  end subroutine SST_Blowing_noCorrection


  subroutine SST_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in),  dimension(nRANS) :: prim, wall_prim
    real(R8), intent(in) :: rho, wall_rho
    real(R8), intent(out), dimension(nRANS) :: ghost_prim

    ghost_prim(1) = - prim(1) ! enforcing k_wall == 0
    ghost_prim(2) = rho * 2d0 * wall_prim(2)/wall_rho - prim(2)

  end subroutine SST_Extrapolate_Wall


  subroutine SST_Enforce_Realizability ( rhokrhow )
    implicit none
    real(R8), intent(inout),  dimension(:) :: rhokrhow
    if ( rhokrhow(1) < 1d-10 ) rhokrhow(1) = 1d-10
    if ( rhokrhow(2) < 1d-10 ) rhokrhow(2) = 1d-10
  end subroutine SST_Enforce_Realizability


  function compute_F ( rho, w, k, mil, Gradw, Gradk, d ) result ( F )
    ! Gradw(2), Gradk(2) (2D); the MOSE peer has (3).
    implicit none
    real(R8), intent(in) :: rho, w, k, mil, Gradw(2), Gradk(2), d
    real(R8) :: CD, arg1, arg2, F(2)

    CD = Max ( 2d0*rho*sigma_w2/w*dot_product( Gradw, Gradk ), 1d-10 )
    arg1 = Min ( Max ( Sqrt(k)/(beta_star*w*d), &
                       5d2*mil/(rho*w*d**2) ), &
                 4d0*sigma_w2*rho*k/(CD*d**2) )
    arg2 = Max ( 2d0*Sqrt(k)/(beta_star*w*d), 5d2*mil/(rho*w*d**2) )
    F(1) = Tanh ( arg1**4 )
    F(2) = Tanh ( arg2**2 )
  end function compute_F


  function compute_F1 ( rho, w, k, mil, Gradw, Gradk, d ) result ( F1 )
    implicit none
    real(R8), intent(in) :: rho, w, k, mil, Gradw(2), Gradk(2), d
    real(R8) :: CD, arg1, F1

    CD = Max ( 2d0*rho*sigma_w2/w*dot_product( Gradw, Gradk ), 1d-10 )
    arg1 = Min ( Max ( Sqrt(k)/(beta_star*w*d), &
                       5d2*mil/(rho*w*d**2) ), &
                 4d0*sigma_w2*rho*k/(CD*d**2) )
    F1 = Tanh ( arg1**4 )
  end function compute_F1


  function compute_F2 ( rho, w, k, mil, d ) result ( F2 )
    implicit none
    real(R8), intent(in) :: rho, w, k, mil, d
    real(R8) :: arg2, F2

    arg2 = Max ( 2d0*Sqrt(k)/(beta_star*w*d), 5d2*mil/(rho*w*d**2) )
    F2 = Tanh ( arg2**2 )

  end function compute_F2

end module Q2D_Lib_SST
