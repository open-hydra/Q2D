!>@brief: Spalart-Allmaras model (without the ft2 term). https://turbmodels.larc.nasa.gov/spalart.html
!
! 12 procedures plus Spalart_Source_Diffusive_Q2D (area-variation diffusive
! source). 2D loop, 9/12-argument signatures with the additional arguments
! vr/r_ave (axisymmetric) and V_dot_G (Definition A).
!
! With obj_rans%SpalartShur active, Spalart_Source_Terms invokes
! Compute_Velocity_Gradient + Compute_RC_Terms and SA_Source_Blk computes
! fr1 = (1+cr1)*(2*rstar/(1+rstar))*(1-cr3*atan(cr2*rtilde)) - cr1, so that
! Production = fr1*cb1*Stilde*nit. The SAR and SAcomp variants are read from the
! obj_rans flags (set in Setup_RANS_Model); obj_rot%enabled activates the
! absolute-vorticity correction.
module Q2D_Lib_Spalart
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  real(R8), parameter, private :: &
    cb1 = 0.1355d0, &
    cb2 = 0.622d0, &
    cw2 = 0.3d0, &
    cw3 = 2.0d0, &
    cv1 = 7.1d0, &
    cK = 0.41d0, &
    ! ct1..ct4 intentionally absent: this is the TMR SA-noft2 variant
    ! (no trip, f_t2 == 0); the trip-function constants would be
    ! dead code here.
    sigma = 2d0/3d0, &
    c2 = 0.7d0, &
    c3 = 0.9d0, &
    cr1 = 1.d0, &
    cr2 = 12.d0, &
    cr3 = 1.d0, &
    Crot = 2.d0, &
    cw1 = cb1/cK**2 + (1d0+cb2)/sigma, &
    ! (the QCR constant C_cr1 is defined in Q2D_Lib_QCR2000, not here.)
    tolerance = 1.d-2

  public :: Spalart_Source_Terms
  public :: Spalart_Set_Wall_Values, Spalart_Set_Blowing_Wall
  public :: Spalart_Extrapolate_Wall, Spalart_Enforce_Realizability
  public :: Spalart_Eddy_Viscosity, Spalart_RANS_Diffusive_Flux
  public :: Spalart_Source_Diffusive_Q2D

contains

  subroutine Spalart_Source_Terms ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Lib_Fluid,        only: Vorticity_Vector, Vorticity_Tensor
    use Q2D_Lib_SpalartShur,  only: Compute_Velocity_Gradient, Compute_RC_Terms
    use Q2D_Config_Types_m,   only: obj_rans
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: b

    ! Spalart-Shur rotation/curvature precomputation.
    ! Compute_Velocity_Gradient + Compute_RC_Terms populate vel_gradient +
    ! rc_term1 + rc_term2 BEFORE the block loop. The conditional allocation in
    ! Mod_Allocate_Data ensures the arrays exist if and only if SpalartShur.
    if ( obj_rans%SpalartShur ) then
      call Compute_Velocity_Gradient ( domain )
      call Compute_RC_Terms ( domain )
    end if

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks

      call SA_Source_Blk ( domain % blk(b) % P,                 &
                           domain % blk(b) % r,                 &
                           domain % blk(b) % M,                 &
                           domain % blk(b) % yn,                &
                           domain % blk(b) % sur,               &
                           domain % blk(b) % vel_gradient,      &
                           domain % blk(b) % rc_term1,          &
                           domain % blk(b) % rc_term2,          &
                           domain % blk(b) % node,              &
                           domain % blk(b) % dim,               &
                           domain % blk(b) % area_geom_term_x,  &
                           domain % blk(b) % area_geom_term_y   )

    end do

  end subroutine Spalart_Source_Terms


  subroutine SA_Source_Blk ( Prim, Res, M, yn, sur, gradv, rc1, rc2, r, n, geom_x, geom_y )
    ! 12-argument signature (sur, r, geom_x, geom_y for the area variation);
    ! body with a 2D loop (do j; do i — no do k).
    ! sur is the 2D analogue of the 3D cell volume (cell surface as a volume
    ! metric).
    use Q2D_Base_Types_m,        only: Q2D_tensor_2D_type, Q2D_vector_2D_type
    use Q2D_Global_m,            only: nprim, nsc, nu, nv, np, nt, gc
    use Q2D_Config_Types_m,      only: obj_rans, obj_rot
    use Q2D_Lib_Fluid,           only: Vorticity_Vector, Vorticity_Tensor, Strain_Tensor
    use FLINT_Lib_Thermodynamic, only: Co_rotot_Rtot, f_laminarViscosity

    implicit none
    integer, intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: yn
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: sur
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: M
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: gradv
    real(R8), dimension(n(1), n(2)), intent(in) :: rc1, rc2
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)), intent(in) :: r
    ! Q2D area-variation geometry (Definition A): cell-centered,
    ! universally zeroed when the area variation is off.
    ! Relevant for the strain-based branches (SAR: |S| in the Production; SpalartShur:
    ! D2 denominator); for the standard SA-noft2 (vorticity-based |omega|) the
    ! 2(V·G)^2 correction does not apply -> Eij (and hence V_dot_G) does not enter
    ! the Production.
    real(R8), dimension(1:n(1), 1:n(2)), intent(in) :: geom_x, geom_y
    integer :: i, j
    real(R8) :: rho, Rgas, nit, mil, Gradvel(2,2), Gradnit(2), omega(3), Om, abs_Eij, Wij(3,3)
    real(R8) :: Eij(3,3), chi, Stilde, Sbar, fw, Gradnit2, r_Fun, g_Fun, fr1, rstar, D2, rtilde
    real(R8) :: Production, Diffusion, Destruction, Source
    real(R8) :: r_ave, V_dot_G

    !$omp do collapse (2) private ( rho, Rgas, nit, mil, Gradvel, Gradnit, omega, Om, abs_Eij, Wij ), &
    !$omp private ( Eij, chi, Stilde, Sbar, fw, Gradnit2, r_Fun, g_Fun, fr1, rstar, D2, rtilde ), &
    !$omp private ( Production, Diffusion, Destruction, Source, r_ave, V_dot_G, i, j )

    do j = 1, n(2)
    do i = 1, n(1)

      r_ave = 0.25d0 * (r(i  ,j  )%c(2) + r(i-1  ,j  )%c(2) + r(i  ,j-1  )%c(2) + r(i-1  ,j-1  )%c(2))
      call Co_rotot_Rtot ( Prim(1:nsc,i,j), rho, Rgas )
      mil = f_laminarViscosity ( Prim(1:nsc,i,j), Prim(np,i,j), rho, Rgas )
      nit = Prim (nt,i,j) / rho  ! nu-tilde
      chi = Prim (nt,i,j) / mil  ! Spalart chi variable

      Gradvel(:,1) = ( Prim(nu:nv,i+1,j) - Prim(nu:nv,i-1,j) ) * 0.5d0
      Gradvel(:,2) = ( Prim(nu:nv,i,j+1) - Prim(nu:nv,i,j-1) ) * 0.5d0
      Gradvel = Matmul ( Gradvel, M(i,j) % c )

      omega = Vorticity_Vector ( Gradvel )    ! 2D reduction (only omega(3) ≠ 0)
      ! Absolute vorticity for SA:
      ! omega_abs = omega_rel + 2*Omega_frame, physically correct for Om,
      ! Stilde and D2 (SA-RC denominator). Vorticity_Vector returns
      ! vector(3) (omega(1)=omega(2)=0); Omega_vec axis strictly [0,0,1].
      if (obj_rot%enabled) omega = omega + 2.0_R8 * obj_rot%Omega_vec
      Om    = Sqrt ( sum ( omega**2 ) )

      Sbar = nit / ( cK**2 * yn(i,j)**2 ) * fv2 (chi)
      if ( Sbar >= (-c2*Om) ) then
        Stilde = Om + Sbar
      else
        Stilde = Om + Om * (c2**2*Om + c3*Sbar) / ((c3 - 2*c2)*Om - Sbar)
      endif
      if ( Stilde == 0d0 ) Stilde = 1d-40

      ! Strain-rate tensor (Q2D Def A: V_dot_G correction; Eij(3,3)=-V·G).
      ! Always computed: SAR uses abs_Eij in the Production, SpalartShur uses abs_Eij in the
      ! D2 denominator for rtilde. Standard SA does not use it (negligible overhead).
      V_dot_G = Prim(nu,i,j)*geom_x(i,j) + Prim(nv,i,j)*geom_y(i,j)
      Eij     = Strain_Tensor ( Gradvel, Prim(nv,i,j), r_ave, V_dot_G=V_dot_G )
      abs_Eij = sqrt ( 2*sum (Eij**2) )

      if ( obj_rans%SpalartShur ) then
        ! SA-RC (Spalart & Shur 1997): rotation/curvature correction.
        ! fr1 multiplies the production; rc_term1=rstar=S/Omega
        ! from Compute_RC_Terms (formal cylindrical+area-var treatment); rc_term2 is the
        ! numerator of rtilde (SA-specific denominator applied here: D2^2).
        Wij = Vorticity_Tensor ( omega )    ! ABSOLUTE omega (set above)
        rstar  = rc1(i,j)
        D2     = 0.5d0*(abs_Eij**2 + Om**2)        ! SA-RC denominator (Spalart-Shur 1997 §3.2)
        rtilde = rc2(i,j) / D2**2
        fr1 = (1.d0 + cr1)*(2.d0*rstar/(1.d0 + rstar))*(1.d0 - cr3*datan(cr2*rtilde)) - cr1
        Production = fr1 * cb1 * Stilde * nit

      else if ( obj_rans%SAR ) then
        ! SA-R (Spalart 2000, strain-based): adds |Sij| on top of |omega|.
        Production = cb1 * ( Stilde + Crot * Min(0d0, abs_Eij - Om) ) * nit

      else
        ! Standard SA (only |omega|, vorticity-based)
        Production = cb1 * Stilde * nit
      end if

      Gradnit(1) = ( Prim(nt,i+1,j) / sum( Prim(1:nsc,i+1,j) ) - &
                     Prim(nt,i-1,j) / sum( Prim(1:nsc,i-1,j) ) ) * 0.5d0
      Gradnit(2) = ( Prim(nt,i,j+1) / sum( Prim(1:nsc,i,j+1) ) - &
                     Prim(nt,i,j-1) / sum( Prim(1:nsc,i,j-1) ) ) * 0.5d0
      Gradnit = Matmul ( Gradnit, M(i,j) % c )

      Gradnit2 = Sum ( Gradnit**2 )
      Diffusion = cb2 / sigma * Gradnit2

      r_Fun = nit / ( Stilde * cK**2 * yn(i,j)**2 )
      r_Fun = Min ( r_Fun, 10d0 )
      g_Fun = r_Fun + cw2 * ( r_Fun**6 - r_Fun )
      fw = g_Fun * ( ( 1d0 + cw3**6 ) / ( g_Fun**6 + cw3**6 ) )**(1.d0/6.d0)
      Destruction = cw1 * fw * ( nit / yn(i,j) )**2

      if ( obj_rans%SAcomp ) call Compressibility_Correction &
        ( nit, rho, chi, Prim(1:nsc,i,j), Prim(np,i,j), Rgas, Om, Production )

      Source = rho * ( Production + Diffusion - Destruction ) * sur(i,j)
      Res(nt,i,j) = Res(nt,i,j) - Source

    enddo ; enddo

  end subroutine SA_Source_Blk


  function func (yy,rr) result (result)
    implicit none
    real(R8), intent(in) :: yy, rr
    real(R8) :: result, ff
    ff = 0.44d0 / ( 1d0 + 14d0 * yy**5d0 ) + 0.56d0
    result = yy * yy * ff - rr
  end function func


  pure function fv1 (chi) result (result)
    implicit none
    real(R8), intent(in) :: chi
    real(R8) :: result
    result = ( chi**3 ) / ( chi**3 + cv1**3 )
  end function fv1


  function fv2 (chi) result (result)
    implicit none
    real(R8) :: result
    real(R8), intent(in) :: chi
    result = 1d0 - chi / ( 1d0 + chi * fv1(chi) )
  end function fv2


  ! Paciorri-Sabetta compressibility correction https://doi.org/10.2514/2.3967
  subroutine Compressibility_Correction ( nit, rho, chi, rhoi, p, Rgas, om, prod )
    use Q2D_Global_m,            only: nsc
    use FLINT_Lib_Thermodynamic, only: f_ss

    implicit none
    real(R8), intent(in) :: nit, rho, chi, rhoi(nsc), p, Rgas, om
    real(R8), intent(inout) :: prod
    real(R8) :: sound, c_mach, staux, err, deltamach
    real(R8) :: aux_num, aux_den, c_plus, c_minus, res, f_cor1, f_cor2, f_cor

    sound = f_ss ( rhoi, p, rho, Rgas )
    c_mach = 1.d0
    staux = 25.d0*(om*nit*fv1(chi))/sound
    err = 1.d0
    do while (err>tolerance)
      deltamach=.1d0*c_mach
      aux_num=func(c_mach,staux)*2.*deltamach
      c_plus=c_mach+deltamach
      c_minus=c_mach-deltamach
      aux_den=func(c_plus,staux)-func(c_minus,staux)
      res=aux_num/aux_den
      err=dabs(res)
      c_mach=c_mach-res
    end do
    f_cor1=0.6d0*(1.d0/(1.d0+9.d0*(c_mach)**6.))+0.4d0
    f_cor2=0.44d0*(1.d0/(1.d0+14.d0*(c_mach)**5.))+0.56d0
    f_cor=f_cor1*f_cor2
    prod=prod*f_cor

  end subroutine Compressibility_Correction


  subroutine Spalart_Set_Wall_Values ( mi_l, mitilde, dist )
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: mitilde
    real(R8), intent(in)                     :: mi_l
    real(R8), intent(in)                     :: dist

    mitilde = 0d0 ! smooth solid wall condition

  end subroutine Spalart_Set_Wall_Values


  subroutine Spalart_Set_Blowing_Wall ( rho, mil, tau, mitilde, mdot, dist )
    ! tau(2): 2D-reduced wall stress vector.
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(out), dimension(nRANS)  :: mitilde
    real(R8), intent(in)                     :: rho, mil, mdot, dist
    real(R8), dimension(2), intent(in)       :: tau

    mitilde = 0d0 ! smooth solid wall condition

  end subroutine Spalart_Set_Blowing_Wall


  subroutine Spalart_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in),  dimension(nRANS) :: prim, wall_prim
    real(R8), intent(in) :: rho, wall_rho
    real(R8), intent(out), dimension(nRANS) :: ghost_prim

    ghost_prim(1) = - prim(1) ! enforces mit_wall == 0

  end subroutine Spalart_Extrapolate_Wall


  subroutine Spalart_Enforce_Realizability ( mitilde )

    implicit none
    real(R8), intent(inout),  dimension(:) :: mitilde
    if ( mitilde(1) < 1d-10 ) mitilde(1) = 1d-10

  end subroutine Spalart_Enforce_Realizability


  subroutine Spalart_Eddy_Viscosity ( mut, mitilde, milam, rho, Gradvel, vr, r_ave, dist, V_dot_G )
    ! 9-argument signature (vr, r_ave, optional V_dot_G accepted for interface
    ! conformity; standard SA-noft2 is chi-based and does not use them in the body).
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS) :: mitilde
    real(R8), intent(in)                   :: milam
    real(R8), intent(in)                   :: rho
    real(R8), intent(in)                   :: dist
    real(R8), intent(in), dimension(2,2)   :: Gradvel
    real(R8), intent(in)                   :: vr, r_ave
    real(R8), intent(out)                  :: mut
    real(R8), intent(in), optional         :: V_dot_G
    real(R8) :: chi

    chi = mitilde(1)/milam
    mut = fv1( chi )*mitilde(1)

  end subroutine Spalart_Eddy_Viscosity


  subroutine Spalart_RANS_Diffusive_Flux (flux, mitilde, milam, gradvel, gradnit, rho, int_l, normal, vr, r_ave, dist, mut)
    ! 12-argument signature. int_l is the 2D analogue of the 3D face area.
    ! gradvel/vr/r_ave/dist/mut accepted for interface conformity (the SA flux
    ! uses mitilde directly).
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS)   :: mitilde
    real(R8), intent(in)                     :: milam
    real(R8), intent(in), dimension(2,2)     :: gradvel
    real(R8), intent(in), dimension(nRANS,2) :: gradnit
    real(R8), intent(in)                     :: rho
    real(R8), intent(in)                     :: int_l, normal(2), dist
    real(R8), intent(in)                     :: vr, r_ave
    real(R8), intent(out), dimension(nRANS)  :: flux
    real(R8), intent(in), optional           :: mut

    flux(1) = int_l*(milam + mitilde(1))*dot_product(gradnit(1,:), normal)/sigma

  end subroutine Spalart_RANS_Diffusive_Flux


  subroutine Spalart_Source_Diffusive_Q2D (source, mitilde, milam, mut, gradvel, gradnit, rho, dist)
    ! SA diffusion projected onto G in Source_Area_Visc_Blk (area-variation source).
    ! mut/gradvel accepted for interface conformity (SA does not use them).
    use Q2D_Global_m, only: nRANS

    implicit none
    real(R8), intent(in), dimension(nRANS)    :: mitilde
    real(R8), intent(in)                      :: milam
    real(R8), intent(in)                      :: mut
    real(R8), intent(in), dimension(2,2)      :: gradvel
    real(R8), intent(in), dimension(nRANS,2)  :: gradnit
    real(R8), intent(in)                      :: rho
    real(R8), intent(in)                      :: dist
    real(R8), intent(out), dimension(nRANS,2) :: source

    source(1,1) = (milam + mitilde(1)) * gradnit(1,1) / sigma
    source(1,2) = (milam + mitilde(1)) * gradnit(1,2) / sigma

  end subroutine Spalart_Source_Diffusive_Q2D

end module Q2D_Lib_Spalart
