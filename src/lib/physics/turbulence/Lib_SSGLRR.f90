!>@brief: SSG/LRR-omega Reynolds-stress transport (Q2D 2D reduction, 4 stresses + omega).
!         see https://turbmodels.larc.nasa.gov/rsm-ssglrr.html.
!
! Model identity (verified against the primary source Eisfeld & Rumsey,
! "Length-Scale Correction for Reynolds-Stress Modeling", AIAA J. 2020, Appendix A
! Eqs 17-37 + Table 2): the algebra below is the *baseline* SSG/LRR-omega
! ("w2012") — the omega destruction is beta*rho*omega**2 (Eq 7/29 with the
! length-scale correction F^(LSC) == 0). The w2019 variant adds to the omega
! destruction only the LSC term (Eq 29), active solely near stagnation/
! reattachment (chi != 0 for L_t/l_log > 1.466) — absent here; on attached flows
! (e.g. flat plate) the LSC is inactive, so there baseline == w2019.
! Coefficient convention: the anisotropy is a_ij = R_ij/k - 2/3 delta = 2 b_ij
! (with b_ij the Table 2 normalization, R_ij/2k - 1/3 delta), so the constants
! below are Table 2 rescaled (C1,C1*,C3*,C4,C5 halved; C2 divided by four;
! C3,D unchanged) — algebraically identical to the paper.
!
! Per-regime Q2D reduction, formally derived in:
!   docs/theory/q2d-rans-closures.md §6 + Appendix A.
!
! The 3D model transports 7 variables (rho*R11,R22,R33,R12,R13,R23, rho*omega).
! Every Q2D geometric regime eliminates the two cross stresses R13,R23, which
! vanish by a mirror symmetry of the regime (reflection in theta in 2D-axi;
! parity in z in area variation; partial_z==0 in 2D-plane) and are never
! regenerated (7->5 reduction proven closed, Appendix A.6).
! Hence nRANS=5:
!   var(1)=rho*R11  var(2)=rho*R22  var(3)=rho*R33  var(4)=rho*R12  var(5)=rho*omega.
! R33 is RETAINED in every regime (normal stress fed by the pressure-strain even
! where its direct production P33 vanishes — "four stresses, not three").
!
! The three regimes share ONE SINGLE structure (theory §6.1):
!   - the SAME 3x3 SSG+LRR algebra for pressure-strain / dissipation /
!     realizability (full traces over 3x3);
!   - the SAME exact production P_ij = -rho(R_ik d_k U_j + R_jk d_k U_i)
!     fed by a BLOCK-DIAGONAL mean-velocity gradient whose only out-of-plane
!     component is the regime dilatation
!     (grad U)_33 = S_33 in {0, u_r/r, -V·G}, exactly the value already provided
!     by Strain_Tensor_Std/_Axi (no 2x2->3x3 widening). The out-of-plane
!     curvature/dilatation production is therefore P33 = -2 rho R33 S_33.
!   - geometric transport/diffusion sources IDENTICAL IN FORM to the scalar
!     mean-flow RANS ones: the convective -rho*u_r*R_ij/r (2D-axi) /
!     +rho*R_ij*(V·G) (area variation) are applied generically to the
!     post-pressure indices by Source_Axi_Eul_Blk / Source_Area_Eul_Blk; the
!     diffusive +F_r/r (2D-axi) / -F·G (area variation) is routed through
!     RANS_Diffusive_SourceQ2D (SSGLRR_Source_Diffusive_Q2D) by
!     Source_Axi/Area_Visc_Blk.
!
! Momentum coupling (no Boussinesq, theory §6.8): total stress =
! laminar(mu_l) minus rho*R_ij, via Stress_Vector_RSM /
! Stress_Tensor_RSM (Q2D_Lib_Fluid).
!
! Interface conformity: the bodies below follow the abstract 2D interfaces
! of Q2D's Lib_RANS (9/12 arguments, with the extra vr/r_ave + V_dot_G).
module Q2D_Lib_SSGLRR
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  real(R8), dimension(3,3), parameter :: &
  deltaij = reshape((/ 1d0, 0d0, 0d0, &
                       0d0, 1d0, 0d0, &
                       0d0, 0d0, 1d0  /), shape(deltaij), order=(/2,1/) )

  real(R8), parameter, private :: &
  C_mu = 0.09d0, &
  ! omega-branch coefficients (near-wall LRR-omega)
  alphaw_w = 0.5556d0, &
  betaw_w = 0.075d0, &
  sigmaw_w = 0.5d0, &
  sigmad_w = 0d0, &
  C1_w = 1.8d0, &
  C1star_w = 0d0, &
  C2_w = 0d0, &
  C3_w = 0.8d0, &
  C3star_w = 0d0, &
  C2_LRR = 0.52d0, &
  C4_w = 0.5d0*(18d0*C2_LRR + 12d0)/11d0, &
  C5_w = 0.5d0*(-14d0*C2_LRR+20d0)/11d0, &
  D_w = 0.75d0*C_mu, &
  ! epsilon-branch coefficients (free-stream SSG-epsilon)
  alphaw_e = 0.44d0, &
  betaw_e = 0.0828d0, &
  sigmaw_e = 0.856d0, &
  sigmad_e = 1.712d0, &
  C1_e = 1.7d0, &
  C1star_e = 0.9d0, &
  C2_e = 1.05d0, &
  C3_e = 0.8d0, &
  C3star_e = 0.65d0, &
  C4_e = 0.625d0, &
  C5_e = 0.2d0, &
  D_e = 0.22d0

  public :: SSGLRR_Source_Terms
  public :: SSGLRR_Eddy_Viscosity
  public :: SSGLRR_RANS_Diffusive_Flux, SSGLRR_SD_RANS_Diffusive_Flux
  public :: SSGLRR_Source_Diffusive_Q2D
  public :: SSGLRR_Set_Wall_Values, SSGLRR_Extrapolate_Wall
  public :: SSGLRR_Enforce_Realizability

contains

  subroutine SSGLRR_Source_Terms ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: b

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks

      call SSGLRR_Blk ( domain % blk(b) % P,                 &
                        domain % blk(b) % r,                 &
                        domain % blk(b) % M,                 &
                        domain % blk(b) % sur,               &
                        domain % blk(b) % yn,                &
                        domain % blk(b) % node,              &
                        domain % blk(b) % dim,               &
                        domain % blk(b) % area_geom_term_x,  &
                        domain % blk(b) % area_geom_term_y   )

    end do

  end subroutine SSGLRR_Source_Terms


  subroutine SSGLRR_Blk ( Prim, Res, M, sur, WDist, r, n, geom_x, geom_y )
    ! 2D block loop. sur is the 2D analogue of the 3D cell volume.
    ! The block-diagonal velocity gradient is built from the stored 2x2
    ! gradient plus the regime out-of-plane strain S_33 supplied by
    ! Strain_Tensor_Std/_Axi (0 / u_r/r with de l'Hopital / -V·G Definition A).
    ! Production, pressure-strain, dissipation and the omega source are
    ! the full 3x3 model algebra evaluated on the 5 retained variables;
    ! the cross stresses R13,R23 are identically zero and never regenerated
    ! (7->5 reduction, Appendix A.6), so only 5 residuals are assembled.
    ! The convective/diffusive geometric sources are added separately by the
    ! Source_Axi/Area_{Eul,Visc}_Blk wiring (theory §6.1).
    use Q2D_Base_Types_m,        only: Q2D_tensor_2D_type, Q2D_vector_2D_type
    use Q2D_Global_m,            only: nprim, nsc, nu, nv, np, nt, gc
    use Q2D_Config_Types_m,      only: obj_rot
    use Q2D_Lib_Fluid,           only: Vorticity_Vector, Vorticity_Tensor, Strain_Tensor
    use FLINT_Lib_Thermodynamic, only: Co_rotot_Rtot, f_laminarViscosity

    implicit none
    integer, intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: sur
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: WDist
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: M
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)), intent(in) :: r
    real(R8), dimension(1:n(1), 1:n(2)), intent(in) :: geom_x, geom_y
    ! Local
    integer  :: i, j, ii, jj, kk
    real(R8) :: rho, Rgas, mil, rhoRij(3,3), omega, kappa, rhoe, dist
    real(R8) :: Gradvel(2,2), Gradvel3(3,3), Grad(5,2), gradk(2), Divel, Sij(3,3), Sij_star(3,3), Wij(3,3)
    real(R8) :: omega_w(3), r_ave, V_dot_G
    real(R8) :: F1, F2, alphaw, betaw, sigmad, C1, C2, C3, C4, C5, C1star, C3star
    real(R8), dimension(3,3) :: rhoP, rhoeij, aij, C1_term, C2_term, C3_term, C4_term, C5_term, aikakj, rhoPI
    real(R8) :: rhoPkk, aklakl, aklSkl, Diff, Source(5)

    !$omp do collapse (2) private (rho, Rgas, mil, rhoRij, omega, kappa, rhoe, dist) &
    !$omp private (Gradvel, Gradvel3, Grad, gradk, Divel, Sij, Sij_star, Wij, omega_w, r_ave, V_dot_G) &
    !$omp private (F1, F2, alphaw, betaw, sigmad, C1, C2, C3, C4, C5, C1star, C3star) &
    !$omp private (rhoP, rhoeij, aij, C1_term, C2_term, C3_term, C4_term, C5_term, aikakj, rhoPI) &
    !$omp private (rhoPkk, aklakl, aklSkl, Diff, Source, i, j, ii, jj, kk)

    do j = 1, n(2)
    do i = 1, n(1)

      r_ave = 0.25d0 * ( r(i,j)%c(2) + r(i-1,j)%c(2) + r(i,j-1)%c(2) + r(i-1,j-1)%c(2) )

      ! Local variables (var(1:4)=rho*Rij, var(5)=rho*omega; R13=R23 eliminated)
      call Co_rotot_Rtot ( Prim(1:nsc,i,j), rho, Rgas )
      mil = f_laminarViscosity ( Prim(1:nsc,i,j), Prim(np,i,j), rho, Rgas )
      rhoRij(1,1) = Prim(nt  ,i,j)
      rhoRij(2,2) = Prim(nt+1,i,j)
      rhoRij(3,3) = Prim(nt+2,i,j)
      rhoRij(1,2) = Prim(nt+3,i,j)
      rhoRij(2,1) = rhoRij(1,2)
      rhoRij(1,3) = 0d0 ; rhoRij(3,1) = 0d0   ! R13 eliminated (regime mirror symmetry)
      rhoRij(2,3) = 0d0 ; rhoRij(3,2) = 0d0   ! R23 eliminated
      omega = Prim(nt+4,i,j) / rho
      kappa = 0.5d0 * ( rhoRij(1,1) + rhoRij(2,2) + rhoRij(3,3) ) / rho
      rhoe = rho * C_mu * kappa * omega
      dist = WDist(i,j)

      ! In-plane velocity gradient d u_a / d x_b (a,b in {1,2})
      Gradvel(:,1) = ( Prim(nu:nv,i+1,j) - Prim(nu:nv,i-1,j) ) * 0.5d0
      Gradvel(:,2) = ( Prim(nu:nv,i,j+1) - Prim(nu:nv,i,j-1) ) * 0.5d0
      Gradvel = Matmul ( Gradvel, M(i,j) % c )

      ! Q2D Definition A: V·G at cell center -> Strain_Tensor S_33 = -V·G (area variation).
      V_dot_G = Prim(nu,i,j)*geom_x(i,j) + Prim(nv,i,j)*geom_y(i,j)

      ! Gradient of (R11,R22,R33,R12,omega) and of k = 0.5(R11+R22+R33).
      Grad(:,1) = ( Prim(nt:nt+4,i+1,j) / sum( Prim(1:nsc,i+1,j) ) - &
                    Prim(nt:nt+4,i-1,j) / sum( Prim(1:nsc,i-1,j) ) ) * 0.5d0
      Grad(:,2) = ( Prim(nt:nt+4,i,j+1) / sum( Prim(1:nsc,i,j+1) ) - &
                    Prim(nt:nt+4,i,j-1) / sum( Prim(1:nsc,i,j-1) ) ) * 0.5d0
      Grad = Matmul ( Grad, M(i,j) % c )
      gradk = 0.5d0 * ( Grad(1,:) + Grad(2,:) + Grad(3,:) )

      ! Block-diagonal mean velocity gradient (Appendix A.1):
      ! (grad U)_33 = S_33 from Strain_Tensor (procedure pointer _Std/_Axi -> 0 / vr/r / -V·G).
      Sij = Strain_Tensor ( Gradvel, Prim(nv,i,j), r_ave, V_dot_G=V_dot_G )
      Gradvel3 = 0d0
      Gradvel3(1:2,1:2) = Gradvel
      Gradvel3(3,3) = Sij(3,3)
      Divel = Gradvel3(1,1) + Gradvel3(2,2) + Gradvel3(3,3)
      Sij_star = Sij - 1d0/3d0 * Divel * deltaij

      ! In-plane vorticity (only omega_3 = dv/dx - du/dy survives), absolute axes.
      omega_w = Vorticity_Vector ( Gradvel )
      if (obj_rot%enabled) omega_w = omega_w + 2.0_R8 * obj_rot%Omega_vec
      Wij = Vorticity_Tensor ( omega_w )

      ! Coefficient blending (F1 -> 1 at the wall, -> 0 in the free stream)
      F1 = compute_F1 ( rho=rho, w=omega, k=kappa, mil=mil, gradw=Grad(5,:), gradk=gradk, d=dist )
      F2 = 1d0 - F1
      alphaw = F1*alphaw_w + F2*alphaw_e
      betaw  = F1*betaw_w  + F2*betaw_e
      sigmad = F1*sigmad_w + F2*sigmad_e
      C1     = F1*C1_w     + F2*C1_e
      C2     = F1*C2_w     + F2*C2_e
      C3     = F1*C3_w     + F2*C3_e
      C4     = F1*C4_w     + F2*C4_e
      C5     = F1*C5_w     + F2*C5_e
      C1star = F1*C1star_w + F2*C1star_e
      C3star = F1*C3star_w + F2*C3star_e

      ! Production rho*P_ij = -rho(R_ik d_k U_j + R_jk d_k U_i), block-diagonal gradient.
      rhoP = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        rhoP(ii,jj) = rhoP(ii,jj) + rhoRij(ii,kk)*Gradvel3(jj,kk)
      enddo
      enddo
      enddo
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        rhoP(ii,jj) = rhoP(ii,jj) + rhoRij(jj,kk)*Gradvel3(ii,kk)
      enddo
      enddo
      enddo
      rhoP = -rhoP

      ! Dissipation (isotropic)
      rhoeij = 2d0/3d0 * rhoe * deltaij

      ! Pressure-strain correlation (SSG slow part + LRR rapid part). Full 3x3 traces.
      aij = rhoRij/(rho*kappa) - 2d0/3d0 * deltaij
      aklakl = Sum (aij*aij)

      ! C1 term
      rhoPkk = rhoP(1,1) + rhoP(2,2) + rhoP(3,3)
      C1_term = -( C1*rhoe + 0.5d0*C1star*rhoPkk )*aij

      ! C2 term
      aikakj = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        aikakj(ii,jj) = aikakj(ii,jj) + aij(ii,kk)*aij(kk,jj)
      enddo
      enddo
      enddo
      C2_term = C2*rhoe * ( aikakj - 1d0/3d0*aklakl*deltaij )

      ! C3 term
      C3_term = ( C3 - C3star * Sqrt( aklakl ) ) * rho*kappa*Sij_star

      ! C4 term
      C4_term = 0d0
      aklSkl = Sum ( aij*Sij )
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C4_term(ii,jj) = C4_term(ii,jj) + aij(ii,kk)*Sij(jj,kk)
      enddo
      enddo
      enddo
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C4_term(ii,jj) = C4_term(ii,jj) + aij(jj,kk)*Sij(ii,kk)
      enddo
      enddo
      enddo
      C4_term = C4_term - 2d0/3d0*aklSkl*deltaij
      C4_term = C4*rho*kappa * C4_term

      ! C5 term
      C5_term = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C5_term(ii,jj) = C5_term(ii,jj) + aij(ii,kk)*Wij(jj,kk)
      enddo
      enddo
      enddo
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C5_term(ii,jj) = C5_term(ii,jj) + aij(jj,kk)*Wij(ii,kk)
      enddo
      enddo
      enddo
      C5_term = C5*rho*kappa * C5_term

      rhoPI = C1_term + C2_term + C3_term + C4_term + C5_term

      ! Cross-diffusion of omega
      Diff = sigmad * rho/omega * Max( Dot_Product ( gradk, Grad(5,:) ), 0d0 )

      ! Source terms (5 retained: R11,R22,R33,R12,omega). Source(13)=Source(23)=0
      ! identically (7->5 reduction, Appendix A.6), so they are not assembled.
      Source(1) = rhoP(1,1) + rhoPI(1,1) - rhoeij(1,1)
      Source(2) = rhoP(2,2) + rhoPI(2,2) - rhoeij(2,2)
      Source(3) = rhoP(3,3) + rhoPI(3,3) - rhoeij(3,3)
      Source(4) = rhoP(1,2) + rhoPI(1,2) - rhoeij(1,2)
      Source(5) = alphaw*omega/kappa*0.5d0*rhoPkk - betaw*rho*omega**2 + Diff

      ! ── Rotating frame: exact Coriolis (rotational) production ─────────────
      ! In a frame rotating at constant Omega = (0,0,Oz) the second-moment
      ! transport equation contains the EXACT term (no closure)
      !   G_ij = -2 Omega_k ( eps_ikl R_lj + eps_jkl R_li ),
      ! obtained by taking the moments of the fluctuating-momentum equation,
      ! whose Coriolis acceleration is -2 eps_ikl Omega_k u'_l.
      ! It is traceless (G_kk = 0: the k content and the omega row are not
      ! touched) and is NOT absorbable into the production built from the
      ! relative gradients: P(grad_abs) - P(grad_rel) = G/2 exactly, so
      ! G must be added explicitly. It coexists without double counting with the
      ! W_abs shift in the rapid pressure-strain C5 above (non-inertial
      ! form of Speziale-Sarkar-Gatski 1991).
      ! With axis [0,0,1] (constraint imposed by Setup_RotatingFrame) and R13 = R23 = 0:
      !   G11 = +4 Oz R12,  G22 = -4 Oz R12,  G33 = 0,  G12 = 2 Oz (R22 - R11),
      !   G13 = G23 = 0  -> the 7->5 reduction stays closed also under rotation.
      ! Frame invariance: with G the solution in rotating axes coincides with
      ! the rotated inertial solution; without G an O(1) error remains,
      ! independent of the time step. Solver-level regression:
      ! the V&V case rotating-channel-rsm.
      if (obj_rot%enabled) then
        Source(1) = Source(1) + 4d0*obj_rot%Omega_vec(3)*rhoRij(1,2)
        Source(2) = Source(2) - 4d0*obj_rot%Omega_vec(3)*rhoRij(1,2)
        Source(4) = Source(4) + 2d0*obj_rot%Omega_vec(3)*(rhoRij(2,2) - rhoRij(1,1))
      end if

      Res(nt:nt+4,i,j) = Res(nt:nt+4,i,j) - Source * sur(i,j)

    enddo ; enddo

  end subroutine SSGLRR_Blk


  subroutine SSGLRR_Eddy_Viscosity ( mi_t, var, mi_l, rho, Gradvel, vr, r_ave, dist, V_dot_G )
    ! mu_t = rho*k/omega — the turbulent diffusivity for the scalar (species/heat)
    ! transport; the Reynolds stress enters momentum directly (no Boussinesq).
    ! Signature 9-arg Q2D-form (vr/r_ave/V_dot_G accepted for interface uniformity,
    ! not used: mu_t is a pure algebraic function of k and omega, Appendix A.2).
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS) :: var       ! : [rho*Rij rho*omega]
    real(R8), intent(in)                   :: mi_l
    real(R8), intent(in)                   :: rho
    real(R8), intent(in)                   :: dist
    real(R8), intent(in), dimension(2,2)   :: Gradvel
    real(R8), intent(in)                   :: vr, r_ave
    real(R8), intent(out)                  :: mi_t
    real(R8), intent(in), optional         :: V_dot_G
    real(R8) :: rhok, w

    w    = var(5) / rho
    rhok = 0.5d0 * Sum( var(1:3) )
    mi_t = rhok / w

  end subroutine SSGLRR_Eddy_Viscosity


  subroutine SSGLRR_RANS_Diffusive_Flux ( Flux, var, mi_l, Gradvel, Grad, rho, int_l, Normal, vr, r_ave, Dist, mut )
    ! Daly-Harlow generalized-gradient (GGDH) diffusive flux of each stress
    ! component: F^(ij)_m = mu_l d_m R_ij + D_* rho*R_ml d_l R_ij (the stress
    ! tensor itself is the anisotropic diffusivity), theory §6.4. In 2D with
    ! R13=R23=0 only the in-plane R_ml block (var1=rhoR11, var2=rhoR22,
    ! var4=rhoR12) drives the GGDH; omega diffuses with (mu_l + sigmaw*rho*k/omega).
    ! Signature 12-arg Q2D-form; int_l is the 2D analogue of the face area.
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS)   :: var       ! : [rho*Rij rho*omega]
    real(R8), intent(in)                     :: mi_l
    real(R8), intent(in), dimension(2,2)     :: Gradvel
    real(R8), intent(in), dimension(nRANS,2) :: Grad      ! : gradient of [ Rij omega ]
    real(R8), intent(in)                     :: rho
    real(R8), intent(in)                     :: int_l, Normal(2), Dist
    real(R8), intent(in)                     :: vr, r_ave
    real(R8), intent(out), dimension(nRANS)  :: Flux
    real(R8), intent(in), optional           :: mut
    ! Local
    integer  :: i
    real(R8) :: omega, kappa, gradk(2), F1, F2, D_, sigmaw, Diffx(4), Diffy(4)

    omega = var(5)/rho
    kappa = 0.5d0 * Sum ( var(1:3) ) / rho
    gradk = 0.5d0 * ( Grad(1,:) + Grad(2,:) + Grad(3,:) )

    F1 = compute_F1 ( rho=rho, w=omega, k=kappa, mil=mi_l, gradw=Grad(5,:), gradk=gradk, d=Dist )
    F2 = 1d0 - F1
    D_ = ( D_w * F1 + D_e * F2 ) / ( C_mu*omega )
    sigmaw = F1 * sigmaw_w + F2 * sigmaw_e

    Diffx = mi_l*Grad(1:4,1) + D_*( var(1)*Grad(1:4,1) + var(4)*Grad(1:4,2) )
    Diffy = mi_l*Grad(1:4,2) + D_*( var(4)*Grad(1:4,1) + var(2)*Grad(1:4,2) )
    do i = 1, 4
      Flux(i) = Diffx(i)*Normal(1) + Diffy(i)*Normal(2)
    end do
    Flux(5) = ( mi_l + sigmaw*rho*kappa/omega ) * Dot_Product ( Grad(5,:), Normal )
    Flux = Flux * int_l

  end subroutine SSGLRR_RANS_Diffusive_Flux


  subroutine SSGLRR_SD_RANS_Diffusive_Flux ( Flux, var, mi_l, Gradvel, Grad, rho, int_l, Normal, vr, r_ave, Dist, mut )
    ! Shock-detector "-SD" variant: the anisotropic GGDH tensor is replaced by the
    ! isotropic (mu_l + D_/C_mu*mu_t) gradient diffusion (mu_t = rho*k/omega), more
    ! robust across shocks (theory §6.4). Signature 12-arg Q2D-form.
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS)   :: var
    real(R8), intent(in)                     :: mi_l
    real(R8), intent(in), dimension(2,2)     :: Gradvel
    real(R8), intent(in), dimension(nRANS,2) :: Grad
    real(R8), intent(in)                     :: rho
    real(R8), intent(in)                     :: int_l, Normal(2), Dist
    real(R8), intent(in)                     :: vr, r_ave
    real(R8), intent(out), dimension(nRANS)  :: Flux
    real(R8), intent(in), optional           :: mut
    ! Local
    integer  :: i
    real(R8) :: omega, kappa, gradk(2), F1, F2, D_, sigmaw, mi_t

    omega = var(5)/rho
    kappa = 0.5d0 * Sum ( var(1:3) ) / rho
    gradk = 0.5d0 * ( Grad(1,:) + Grad(2,:) + Grad(3,:) )

    F1 = compute_F1 ( rho=rho, w=omega, k=kappa, mil=mi_l, gradw=Grad(5,:), gradk=gradk, d=Dist )
    F2 = 1d0 - F1
    D_ = 0.5d0*C_mu*F1 + 2d0/3d0*0.22d0*F2
    sigmaw = F1 * sigmaw_w + F2 * sigmaw_e
    mi_t = rho*kappa/omega

    do i = 1, 4
      Flux(i) = (mi_l + D_/C_mu*mi_t) * Dot_Product ( Grad(i,:), Normal )
    end do
    Flux(5) = ( mi_l + sigmaw*mi_t ) * Dot_Product ( Grad(5,:), Normal )
    Flux = Flux * int_l

  end subroutine SSGLRR_SD_RANS_Diffusive_Flux


  subroutine SSGLRR_Source_Diffusive_Q2D ( source, var, mi_l, mut, Gradvel, Grad, rho, dist )
    ! GGDH diffusive flux components projected onto the geometry by
    ! Source_Axi_Visc_Blk (+F_r/r) and Source_Area_Visc_Blk (-F·G). Returns
    ! source(i,d) = F^(i)_d (flux of variable i in direction d, no normal/int_l).
    ! Mirrors SSGLRR_RANS_Diffusive_Flux (std GGDH); the "-SD" geometric variant
    ! is the same projection with the isotropic tensor and is only reachable in
    ! 2D-axi/area-variation (the 2D-plane V&V has no geometric source).
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS)    :: var
    real(R8), intent(in)                      :: mi_l
    real(R8), intent(in)                      :: mut
    real(R8), intent(in), dimension(2,2)      :: Gradvel
    real(R8), intent(in), dimension(nRANS,2)  :: Grad
    real(R8), intent(in)                      :: rho
    real(R8), intent(in)                      :: dist
    real(R8), intent(out), dimension(nRANS,2) :: source
    ! Local
    real(R8) :: omega, kappa, gradk(2), F1, F2, D_, sigmaw, Dvisc

    omega = var(5)/rho
    kappa = 0.5d0 * Sum ( var(1:3) ) / rho
    gradk = 0.5d0 * ( Grad(1,:) + Grad(2,:) + Grad(3,:) )

    F1 = compute_F1 ( rho=rho, w=omega, k=kappa, mil=mi_l, gradw=Grad(5,:), gradk=gradk, d=dist )
    F2 = 1d0 - F1
    D_ = ( D_w * F1 + D_e * F2 ) / ( C_mu*omega )
    sigmaw = F1 * sigmaw_w + F2 * sigmaw_e

    source(1:4,1) = mi_l*Grad(1:4,1) + D_*( var(1)*Grad(1:4,1) + var(4)*Grad(1:4,2) )
    source(1:4,2) = mi_l*Grad(1:4,2) + D_*( var(4)*Grad(1:4,1) + var(2)*Grad(1:4,2) )
    Dvisc = mi_l + sigmaw*rho*kappa/omega
    source(5,1) = Dvisc * Grad(5,1)
    source(5,2) = Dvisc * Grad(5,2)

  end subroutine SSGLRR_Source_Diffusive_Q2D


  subroutine SSGLRR_Set_Wall_Values ( mi_l, var, dist )
    ! Solid wall: no-slip kills all fluctuations (R_ij=0); Menter smooth-wall omega.
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: var
    real(R8), intent(in)                     :: mi_l
    real(R8), intent(in)                     :: dist

    var(1:4) = 0d0                       ! rho*Rij = 0 at the wall
    var(5)   = 8d2 * mi_l / ( dist**2 )  ! Menter kw-SST smooth-surface omega

  end subroutine SSGLRR_Set_Wall_Values


  subroutine SSGLRR_Extrapolate_Wall ( var, wall_var, rho, wall_rho, ghost_var )
    ! Ghost enforces rho*Rij_wall == 0 (sign-flip) and extrapolates omega.
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in),  dimension(nRANS) :: var, wall_var
    real(R8), intent(in) :: rho, wall_rho
    real(R8), intent(out), dimension(nRANS) :: ghost_var

    ghost_var(1:4) = - var(1:4)                                    ! face value rho*Rij == 0
    ghost_var(5)   = rho * 2d0 * wall_var(5)/wall_rho - var(5)     ! extrapolate omega_wall

  end subroutine SSGLRR_Extrapolate_Wall


  subroutine SSGLRR_Enforce_Realizability ( rhoRij )
    ! Applied after each update (regime-invariant, theory §6.3): floor the normal
    ! stresses and Schwarz-limit the retained shear. R13,R23 are absent (7->5).
    implicit none
    real(R8), intent(inout),  dimension(:) :: rhoRij

    ! 1st condition: u'u', v'v', w'w' > 0
    if ( rhoRij(1) < 1d-16 ) rhoRij(1) = 1d-16
    if ( rhoRij(2) < 1d-16 ) rhoRij(2) = 1d-16
    if ( rhoRij(3) < 1d-16 ) rhoRij(3) = 1d-16

    ! 2nd condition: |u'v'| <= sqrt( u'u' * v'v' )
    if ( Abs( rhoRij(4) ) > Sqrt( rhoRij(1)*rhoRij(2) ) ) then
      rhoRij(4) = Sign( 1d0, rhoRij(4) ) * Sqrt( rhoRij(1)*rhoRij(2) )
    end if

  end subroutine SSGLRR_Enforce_Realizability


  function compute_F1 ( rho, w, k, mil, Gradw, Gradk, d ) result ( F1 )
    ! SST blending argument (theory §6.5): F1 = tanh(zeta^4), 2D gradients.
    implicit none
    real(R8), intent(in) :: rho, w, k, mil, Gradw(2), Gradk(2), d
    real(R8) :: CD, zeta, F1

    CD = sigmad_e * rho / w * Max ( Dot_Product( Gradw, Gradk ), 1d-20 )
    zeta = Min ( Max( Sqrt(k)/(C_mu*w*d), 5d2*mil/(rho*w*d**2) ), &
                 4d0*sigmaw_e*rho*k/(CD*d**2) )
    F1 = Tanh ( zeta**4 )

  end function compute_F1

end module Q2D_Lib_SSGLRR
