!>@brief: Modello di turbolenza Kappa-Omega Wilcox 2006, variante TMR "Wilcox2006-klim"
!         (Wilcox 2006 baseline + limitatore di produzione dell'equazione di k
!         min(P, 20*beta_star*rho*omega*k); l'equazione di omega usa la produzione
!         ESATTA non limitata gamma*omega/k*P, come prescrive la variante klim —
!         per questo Prod(2) è calcolato PRIMA che Prod(1) venga limitato più sotto).
!         https://turbmodels.larc.nasa.gov/wilcox.html (now tmbwg.github.io/turbmodels)
!
! 6 procedure più 1 procedura extra di Q2D
! (Wilcox2006_Source_Diffusive_Q2D).
!
! Nota fisica: la turbolenza è intrinsecamente 3D (Pope 2000, Turbulent Flows
!       ch. 6: la cascata di Kolmogorov richiede il vortex stretching
!       omega·grad(V), esclusivamente 3D). Il modello RANS 2D Wilcox 2006 è uno
!       strumento ingegneristico standard (modello primario della NASA TMR), ma
!       con accuratezza degradata rispetto a DNS/LES 3D. La RANS 2D di Q2D non
!       altera questo limite intrinseco di modellazione.
!
! Caratteristiche della realizzazione 2D:
!   - firma con Gradvel(2,2) (nessun gradiente di velocità fuori piano);
!   - Strain_Tensor restituisce Eij(3,3) con E(3,3) = -V·G (Std) o vr/r_ave (Axi),
!     incorporato tramite il dispatch a procedure pointer di Lib_Fluid;
!   - Vorticity_Vector restituisce vector(3) con omega(1)=omega(2)=0;
!   - somma tripla sugli INDICI TENSORIALI (1..3) mantenuta nella correzione
!     f_beta (Wij(ii,jj)*Wij(jj,kk)*Sij_hat(kk,ii)) — NON spaziale;
!   - int_l è l'analogo 2D dell'area di faccia 3D, sur del volume 3D;
!   - il flag obj_rans%blowing_corr è letto dal parsing di Setup_RANS_Model;
!   - correzione a vorticità assoluta attivata da obj_rot%enabled.
!
! Riferimenti di letteratura primaria:
!   - Wilcox D.C. (2006) "Turbulence Modeling for CFD" 3rd ed., DCW Industries.
!   - NASA TMR 2D NACA0012 Wilcox 2006 https://turbmodels.larc.nasa.gov/wilcox.html
module Q2D_Lib_Wilcox2006
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  real(R8), parameter, private :: &
    sigma_k   = 0.6d0, &
    sigma_w   = 0.5d0, &
    beta_star = 0.09d0, &
    beta_0    = 0.0708d0, &
    gamma     = 13d0/25d0, &
    Clim      = 7d0/8d0

  public :: Wilcox2006_Source_Terms
  public :: Wilcox2006_Set_Wall_Values
  public :: Wilcox2006_Extrapolate_Wall, Wilcox2006_Enforce_Realizability
  public :: Wilcox2006_Eddy_Viscosity, Wilcox2006_RANS_Diffusive_Flux
  public :: Wilcox2006_Source_Diffusive_Q2D

contains

  subroutine Wilcox2006_Source_Terms ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: b

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: salta i blocchi non locali

      call Wilcox_Blk ( domain % blk(b) % P,                 &
                        domain % blk(b) % r,                 &
                        domain % blk(b) % M,                 &
                        domain % blk(b) % sur,               &
                        domain % blk(b) % node,              &
                        domain % blk(b) % dim,               &
                        domain % blk(b) % area_geom_term_x,  &
                        domain % blk(b) % area_geom_term_y   )

    end do

  end subroutine Wilcox2006_Source_Terms


  subroutine Wilcox_Blk ( Prim, Res, M, sur, r, n, geom_x, geom_y )
    ! Firma a 8 argomenti. Corpo con loop 2D (do j; do i — nessun do k).
    ! sur è l'analogo 2D del volume di cella 3D (metrica di faccia vs volume di cella).
    ! k_energy_coupling letto da obj_rans.
    use Q2D_Base_Types_m,        only: Q2D_tensor_2D_type, Q2D_vector_2D_type
    use Q2D_Global_m,            only: nprim, nsc, nu, nv, np, nt, gc
    use Q2D_Config_Types_m,      only: obj_rans, obj_rot
    use Q2D_Lib_Fluid,           only: Vorticity_Vector, Vorticity_Tensor, Strain_Tensor
    use FLINT_Lib_Thermodynamic, only: Co_rotot_Rtot, f_laminarViscosity

    implicit none
    integer, intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc),        intent(in)    :: sur
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: M
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)),    intent(in)    :: r
    real(R8), dimension(1:n(1), 1:n(2)),                    intent(in)    :: geom_x, geom_y
    integer :: i, j, ii, jj, kk
    real(R8) :: rho, Rgas, mil, kap, ome, Gradvel(2,2), Divel, Sij(3,3), Sij_bar(3,3), diag, ome_hat
    real(R8) :: mi_t, Tij(3,3), Prod(2)
    real(R8) :: Gradkw(2,2), dkDotdw, sigma_d, Diff
    real(R8) :: Sij_hat(3,3), Wij(3,3), Num, X_w, f_beta, beta, Diss(2), Source(2)
    real(R8) :: omega_w(3), r_ave, V_dot_G
    logical  :: k_energy_coupling

    k_energy_coupling = obj_rans%k_energy_coupling

    !$omp do collapse (2) private ( rho, Rgas, mil, kap, ome, Gradvel, Divel, Sij, Sij_bar, diag, ome_hat ), &
    !$omp private ( mi_t, Tij, Prod, Gradkw, dkDotdw, sigma_d, Diff, Sij_hat, Wij, Num, X_w, f_beta, beta ), &
    !$omp private ( Diss, Source, omega_w, r_ave, V_dot_G, i, j, ii, jj, kk )

    do j = 1, n(2)
    do i = 1, n(1)

      r_ave = 0.25d0 * (r(i  ,j  )%c(2) + r(i-1  ,j  )%c(2) + r(i  ,j-1  )%c(2) + r(i-1  ,j-1  )%c(2))
      call Co_rotot_Rtot ( Prim(1:nsc,i,j), rho, Rgas )
      mil = f_laminarViscosity ( Prim(1:nsc,i,j), Prim(np,i,j), rho, Rgas )
      kap = Prim(nt,i,j)   / rho
      ome = Prim(nt+1,i,j) / rho

      ! Gradiente di velocità (solo nu:nv; nessun d/dz fuori piano)
      Gradvel(:,1) = ( Prim(nu:nv,i+1,j) - Prim(nu:nv,i-1,j) ) * 0.5d0
      Gradvel(:,2) = ( Prim(nu:nv,i,j+1) - Prim(nu:nv,i,j-1) ) * 0.5d0
      Gradvel = Matmul ( Gradvel, M(i,j) % c )

      ! Definizione A di Q2D: V·G al centro cella (S_33 = -V·G in Lib_Fluid::Strain_Tensor_Std).
      ! Con area_geom_term=0 (riempimento nullo universale) -> V_dot_G=0 -> Definizione B.
      V_dot_G = Prim(nu,i,j)*geom_x(i,j) + Prim(nv,i,j)*geom_y(i,j)

      ! Tensore di deformazione 3×3 (E(3,3) = -V·G Std o vr/r_ave Axi)
      Sij = Strain_Tensor ( Gradvel, Prim(nv,i,j), r_ave, V_dot_G=V_dot_G )

      ! Divergenza universale dalla traccia di Sij (Std o Axi via dispatch)
      Divel = Sij(1,1) + Sij(2,2) + Sij(3,3)

      ! Deformazione deviatorica Sij_bar
      Sij_bar = Sij
      diag = 1d0/3d0 * Divel
      Sij_bar(1,1) = Sij_bar(1,1) - diag
      Sij_bar(2,2) = Sij_bar(2,2) - diag
      Sij_bar(3,3) = Sij_bar(3,3) - diag

      ! Eddy viscosity with Clim limiter on omega
      ome_hat = Max ( ome, Clim * Sqrt( 2d0 * Sum( Sij_bar**2 ) / beta_star ) )
      mi_t = rho * kap / ome_hat

      ! Turbulent stress tensor (Boussinesq + 2/3 rho*k diagonal)
      Tij = 2d0 * mi_t * Sij
      diag = 2d0/3d0 * mi_t * Divel + 2d0/3d0 * rho * kap
      Tij(1,1) = Tij(1,1) - diag
      Tij(2,2) = Tij(2,2) - diag
      Tij(3,3) = Tij(3,3) - diag

      ! Production = Tij : Gradvel = Tij : Sij (Tij symmetric, antisym part null)
      Prod(1) = sum ( Tij * Sij )
      Prod(2) = ( gamma * ome / kap ) * Prod(1)
      Prod(1) = Min( Prod(1), 20d0*beta_star*rho*ome*kap ) ! k-eqn production limiter

      ! Gradient of kappa and omega (2 spatial dims)
      Gradkw(:,1) = ( Prim(nt:nt+1,i+1,j) / sum( Prim(1:nsc,i+1,j) ) - &
                      Prim(nt:nt+1,i-1,j) / sum( Prim(1:nsc,i-1,j) ) ) * 0.5d0
      Gradkw(:,2) = ( Prim(nt:nt+1,i,j+1) / sum( Prim(1:nsc,i,j+1) ) - &
                      Prim(nt:nt+1,i,j-1) / sum( Prim(1:nsc,i,j-1) ) ) * 0.5d0
      Gradkw = Matmul ( Gradkw, M(i,j) % c )

      ! Cross-Diffusion term
      dkDotdw = Dot_Product ( Gradkw(1,:), Gradkw(2,:) )
      sigma_d = 1d0/8d0
      if ( dkDotdw <= 0d0 ) sigma_d = 0d0
      Diff = sigma_d * rho / ome * dkDotdw

      ! Cross-strain Sij_hat for f_beta dissipation correction
      Sij_hat = Sij
      Sij_hat(1,1) = Sij(1,1) - 0.5d0 * Divel
      Sij_hat(2,2) = Sij(2,2) - 0.5d0 * Divel
      Sij_hat(3,3) = Sij(3,3) - 0.5d0 * Divel

      omega_w = Vorticity_Vector ( Gradvel )    ! 2D-reduced (only omega(3) ≠ 0)
      ! Absolute vorticity for f_beta:
      ! omega_abs = omega_rel + 2*Omega_frame. Q2D Vorticity_Vector returns vector(3)
      ! with omega(1)=omega(2)=0; Omega_vec axis [0,0,1] strict.
      if (obj_rot%enabled) omega_w = omega_w + 2.0_R8 * obj_rot%Omega_vec
      Wij = Vorticity_Tensor ( omega_w )

      ! Triple sum over TENSOR INDICES (1..3): NON-spatial, preserve.
      Num = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        Num = Num + Wij(ii,jj) * Wij(jj,kk) * Sij_hat(kk,ii)
      end do
      end do
      end do

      X_w = Abs ( Num / ( beta_star * ome )**3 )
      f_beta = ( 1d0 + 85d0 * X_w ) / ( 1d0 + 100d0 * X_w )
      beta = beta_0 * f_beta

      ! Dissipation
      Diss(1) = beta_star * rho * kap * ome
      Diss(2) = beta * rho * ome**2

      ! Source (sign convention: Res integrates RHS; -= Source · sur)
      Source(1) = Prod(1) - Diss(1)
      Source(2) = Prod(2) - Diss(2) + Diff

      Res(nt:nt+1,i,j) = Res(nt:nt+1,i,j) - Source * sur(i,j)
      if (k_energy_coupling) Res(np,i,j) = Res(np,i,j) + Source(1) * sur(i,j)

    end do ; end do

  end subroutine Wilcox_Blk


  subroutine Wilcox2006_Eddy_Viscosity ( mi_t, rkw, mi_l, rho, Gradvel, vr, r_ave, dist, V_dot_G )
    ! Signature 9-args Q2D-form. Body uses Strain_Tensor with V_dot_G ->
    ! S_33 = -V·G modifies |Sij_bar| and hence mi_t (Def A baked into the body).
    use Q2D_Global_m,  only: nRANS
    use Q2D_Lib_Fluid, only: Strain_Tensor

    implicit none
    real(R8), intent(in), dimension(nRANS) :: rkw
    real(R8), intent(in)                   :: mi_l
    real(R8), intent(in)                   :: rho
    real(R8), intent(in)                   :: dist
    real(R8), intent(in), dimension(2,2)   :: Gradvel
    real(R8), intent(in)                   :: vr, r_ave
    real(R8), intent(out)                  :: mi_t
    real(R8), intent(in), optional         :: V_dot_G
    real(R8) :: ome, Sij(3,3), Divel, diag, ome_hat, VG

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    ome = rkw(2) / rho
    Sij = Strain_Tensor ( Gradvel, vr, r_ave, V_dot_G=VG )
    Divel = Sij(1,1) + Sij(2,2) + Sij(3,3)
    diag = 1d0/3d0 * Divel
    Sij(1,1) = Sij(1,1) - diag
    Sij(2,2) = Sij(2,2) - diag
    Sij(3,3) = Sij(3,3) - diag
    ome_hat = Max ( ome, Clim * Sqrt( 2d0 * Sum( Sij**2 ) / beta_star ) )
    mi_t = rkw(1) / ome_hat

  end subroutine Wilcox2006_Eddy_Viscosity


  subroutine Wilcox2006_RANS_Diffusive_Flux ( Flux, rkw, mi_l, Gradvel, Gradkw, rho, int_l, Normal, vr, r_ave, Dist, mut )
    ! Signature 12-args. int_l is the 2D analogue of the 3D face area.
    ! Wilcox 2006 simpler than SST: constant sigmas, no F1 blending.
    ! Diffusive viscosity = mi_l + sigma * rho*k/omega (Wilcox 2006 cap. 4 baseline).
    ! mut/Gradvel/vr/r_ave accepted for interface conformance (Wilcox uses rho*k/omega).
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
    real(R8) :: rhok, ome

    rhok = rkw(1)
    ome  = rkw(2) / rho
    Flux(1) = ( mi_l + sigma_k * rhok / ome ) * Dot_Product ( Gradkw(1,:), Normal )
    Flux(2) = ( mi_l + sigma_w * rhok / ome ) * Dot_Product ( Gradkw(2,:), Normal )
    Flux = Flux * int_l

  end subroutine Wilcox2006_RANS_Diffusive_Flux


  subroutine Wilcox2006_Source_Diffusive_Q2D ( source, rkw, mi_l, mut, gradvel, gradkw, rho, dist )
    ! Q2D-extra: Wilcox diffusion projected onto G in Source_Area_Visc_Blk.
    ! Mirror SST_Source_Diffusive_Q2D pattern; constant sigmas (no F1 blending).
    ! mut/gradvel accepted for interface conformance (Wilcox uses rho*k/omega).
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
    real(R8) :: rhok, ome

    rhok = rkw(1)
    ome  = rkw(2) / rho
    source(1,1) = ( mi_l + sigma_k * rhok / ome ) * gradkw(1,1)
    source(1,2) = ( mi_l + sigma_k * rhok / ome ) * gradkw(1,2)
    source(2,1) = ( mi_l + sigma_w * rhok / ome ) * gradkw(2,1)
    source(2,2) = ( mi_l + sigma_w * rhok / ome ) * gradkw(2,2)

  end subroutine Wilcox2006_Source_Diffusive_Q2D


  subroutine Wilcox2006_Set_Wall_Values ( mi_l, rkw, dist )
    ! Smooth-wall BC for k=0, omega=8e2*mi_l/dist^2 (Menter approx).
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: rkw
    real(R8), intent(in)                     :: mi_l
    real(R8), intent(in)                     :: dist

    rkw(1) = 0d0                       ! Solid surface condition on k
    rkw(2) = 8d2 * mi_l / ( dist**2 )  ! Smooth surface BC for omega (Menter kw-SST)

  end subroutine Wilcox2006_Set_Wall_Values


  subroutine Wilcox2006_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in),  dimension(nRANS) :: prim, wall_prim
    real(R8), intent(in) :: rho, wall_rho
    real(R8), intent(out), dimension(nRANS) :: ghost_prim

    ghost_prim(1) = - prim(1)                                      ! enforcing k_wall == 0
    ghost_prim(2) = rho * 2d0 * wall_prim(2)/wall_rho - prim(2)    ! extrapolating from omega_wall

  end subroutine Wilcox2006_Extrapolate_Wall


  subroutine Wilcox2006_Enforce_Realizability ( rhokrhow )
    implicit none
    real(R8), intent(inout), dimension(:) :: rhokrhow

    if ( rhokrhow(1) < 1d-10 ) rhokrhow(1) = 1d-10
    if ( rhokrhow(2) < 1d-10 ) rhokrhow(2) = 1d-10

  end subroutine Wilcox2006_Enforce_Realizability

end module Q2D_Lib_Wilcox2006
