module Q2D_Lib_BC_Fluxes_SRM
  ! BC 502: SRM grain combustion, as MOSE applies it (MOSE_Lib_BC_Fluxes_GSI_SRM). Burn law
  ! `mdot = aCoeff * (P/pRef)^n` (Saint-Robert/Vieille, aCoeff the mass flux of the propellant gas at pRef),
  ! the flame temperature Taf of the propellant as the total temperature of the products and their mass
  ! fractions BC_ci.
  !
  ! The face carries the complete flux of the grain: the pressure of the boundary cell on the face (the
  ! boundary layer assumption p(wall) = p(cell)), the momentum and the enthalpy convected in with the blown
  ! mass, and the viscous terms. No reflective (symmetry) flux is applied alongside it while the grain burns:
  ! it would impose u.n = 0 at a permeable face and return the impermeable-wall pressure, below the pressure
  ! of the cell since the gas moves away from the face (by about rho*c*u_blow in the exact reflective Riemann
  ! problem; SLAU, 6 % of that at M = 0.03). A surface that does not burn (mdot = 0: aCoeff = 0) leaves the
  ! residual alone and reports it (burning): the caller closes the face as an impermeable wall at Taf.
  !
  ! The products leave the wall at the static temperature T_w of h(T_w) + u_b^2/2 = h(Taf), u_b = mdot/rho_wall
  ! (SRM_Static_Temperature), so that the gas gets back the total enthalpy h(Taf) of the tables per unit mass,
  ! that of the propellant; the wall state takes T_w, the composition BC_ci and its gas constant.
  !
  ! 2D formulation:
  !   * ghost-cell indices Ig/Jg; Normal(2), M(2,2), Gradient(nprim,2),
  !     Stress(2), Blowing(2);
  !   * int_l: integrated face length;
  !   * rank-3 residual.
  !
  ! V_dot_G: Stress_Vector is invoked without the optional argument V_dot_G
  ! (default 0), consistent with the absence of local area variation at the wall.
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index, Compute_Wall_Properties
  use Q2D_Lib_Fluid,              only: Stress_Vector
  use Q2D_Lib_RANS,               only: RANS_Set_Blowing_Wall, RANS_Diffusive_Flux, RANS_Extrapolate_Wall
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, co_k_mi_lam_Wilke, f_tabT, h_tab, cp_tab, Ri_tab

  implicit none
  private

  public :: BC_SRM, SRM_Static_Temperature

contains


  subroutine BC_SRM ( Im, Jm, Fm, Blk, aCoeff, n, pRef, Taf, BC_ci, burning, Ovar )
    implicit none
    integer,  intent(in) :: Im, Jm, Fm
    real(R8), intent(in) :: aCoeff, n, pRef, Taf, BC_ci(nsc)
    type(Q2D_block_type), intent(inout) :: Blk
    logical,  intent(out), optional :: burning     ! .false.: the surface does not burn, the caller closes the face
    real(R8), optional, dimension(7), intent(inout) :: Ovar
    ! Locals
    integer  :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Ig, Jg, s
    real(R8) :: Normal(2), int_l, Dist, dl, M(2,2), Prim(nprim), rho, Rgas, Gradient(nprim,2)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, mdot, T, T_w, h_w
    real(R8) :: Stress(2), Blowing(2), Flux(nprim)
    real(R8) :: q_conv
    real(R8) :: r_ave
    logical  :: burns

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    ! Face metrics
    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l
    Dist   = Blk % yn(Im,Jm)
    M      = Blk % M (Im,Jm) % c

    ! Boundary-cell variables
    Prim = Blk % P(:,Im,Jm)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )

    ! Saint-Robert's law, signed by the face orientation: (-modfm2) is +1 on a low face and -1 on a high one, so
    ! the mass is injected into the domain whichever side of the block the grain sits on.
    mdot = aCoeff * (Prim(np)/pRef)**n * (-modfm2)

    ! A surface that does not burn (aCoeff = 0) injects nothing, and a failed state must not be turned into a
    ! flux: leave the residual alone, the caller closes the face (an impermeable wall at Taf).
    burns = ieee_is_finite(mdot) .and. ieee_is_finite(Taf) .and. Taf > 0d0 .and. mdot * (-modfm2) > 0d0
    if (present(burning)) burning = burns
    if (.not. burns) return

    ! The wall state: the products at the static temperature T_w of total Taf, with their composition and its
    ! gas constant, at the pressure of the cell
    T_w      = SRM_Static_Temperature ( Taf, abs(mdot), Prim(np), BC_ci )
    rho_wall = Prim(np) / ( sum ( BC_ci * Ri_tab(1:nsc) ) * T_w )
    Blowing  = (mdot/rho_wall) * Normal ! Blowing velocity

    ! r_ave for Stress_Vector (2D-axisymmetric terms / area variation)
    select case(Fm)
      case(1)
        r_ave = 0.5d0 * ( blk % node(im-1,jm)%c(2) + blk % node(im-1,jm-1)%c(2) )
      case(2)
        r_ave = 0.5d0 * ( blk % node(im,jm)%c(2)   + blk % node(im,jm-1)%c(2)   )
      case(3)
        r_ave = 0.5d0 * ( blk % node(im,jm-1)%c(2) + blk % node(im-1,jm-1)%c(2) )
      case(4)
        r_ave = 0.5d0 * ( blk % node(im,jm)%c(2)   + blk % node(im-1,jm)%c(2)   )
    end select

    q_conv   = 0.d0
    Stress   = 0.d0
    Gradient = 0d0
    Prim_wall(1:nsc) = BC_ci * rho_wall ! Partial densities at the wall
    mil = 0d0
    if (model>0) then
      T = Prim(np) / ( rho * Rgas )
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_w, mil, kl )
      Gradient(np,Dir) = ( T - T_w ) * modfm3
      Gradient(np,:)   = matmul ( Gradient(np,:), M )
      q_conv           = kl * dot_product( Gradient(np,:), Normal )
      Gradient(nu:nv,Dir) = ( Prim(nu:nv) - Blowing ) * modfm3    ! Velocity gradient due to blowing
      Gradient(nu:nv,:)   = matmul ( Gradient(nu:nv,:), M )       ! Velocity-gradient transformation
      Stress = Stress_Vector ( Gradient(nu:nv,:), Normal, mil, 0d0, Prim(nt:), r_ave, 0.d0 )
    endif

    ! Enthalpy of the injected products at T_w, from the thermodynamic tables
    h_w = 0d0
    do s = 1, nsc
      h_w = h_w + BC_ci(s) * f_tabT ( T_w, s, h_tab )
    end do

    ! Fluxes: the complete flux of the face, the pressure of the cell on it included; the energy
    ! mdot (h(T_w) + u_b^2/2) is mdot h(Taf), the total enthalpy of the products
    Flux        = 0d0
    Flux(1:nsc) = - mdot * int_l * BC_ci
    Flux(nu:nv) = ( Stress - mdot * Blowing - Prim(np) * Normal ) * int_l
    Flux(np)    = int_l * ( q_conv - mdot * ( h_w + 0.5d0 * sum( Blowing**2 ) ) &
                             + dot_product( Stress, Blowing ) )

    if (model==2) then
      dl = blk % yn(im,jm)
      call RANS_Set_Blowing_Wall ( rho=rho_wall, mil=mil,                 &
                                    tau=Stress, rans_variables=Prim_Wall(nt:nprim), &
                                    mdot=mdot, dist=dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
      Gradient(nt:nprim,:)   = matmul ( Gradient(nt:nprim,:), M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim),                  &
                                 rans_variables=Prim_Wall(nt:nprim),   &
                                 vel_gradient=Gradient(nu:nv,:),       &
                                 rans_gradient=Gradient(nt:nprim,:),   &
                                 mul=mil, rho=rho_wall,                &
                                 int_l=int_l, normal=Normal,           &
                                 vr=0.d0, r_ave=r_ave, dist=1d-20 )
      ! Ghost-cell extrapolation of the RANS variables
      Ig = Im - guide(Fm,1)
      Jg = Jm - guide(Fm,2)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm) = Blk % r (:,Im,Jm) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties ( stress=Stress, pw=prim(np), qw=q_conv,  &
                                                       mdot=mdot, y=blk%dl(im,jm)%c(dir)*0.5d0, &
                                                       Tw=T_w, rhow=rho_wall, mu=mil,          &
                                                       exit_array=Ovar )

  end subroutine BC_SRM


  pure real(R8) function SRM_Static_Temperature ( T0, g, p, y ) result ( T )
    ! Static temperature of the products of mass fractions `y` that cross a face with the mass flux `g` >= 0
    ! (kg/(m^2 s)) at the pressure `p` and have the total temperature `T0` (the Taf of a grain, 502): the root of
    ! h(T) + u(T)^2/2 = h(T0), u = g R T / p, with h and R of the tables of the gas (f_tabT, Ri_tab), so that the gas
    ! carries the total enthalpy h(T0). The left side grows with T and is not below h(T0) at T0: Newton from T0 comes
    ! down to the root, in two or three steps; cp of the tables as the slope. g = 0 gives T0. As in MOSE.
    implicit none
    real(R8), intent(in) :: T0, g, p, y(nsc)
    ! Local
    integer  :: it, s
    real(R8) :: R, a2, h0, F, dF, dT

    T = T0
    if ( g == 0d0 ) return
    R  = sum ( y * Ri_tab(1:nsc) )
    a2 = ( g * R / p )**2                       ! u^2 = a2 T^2
    h0 = 0d0
    do s = 1, nsc
      h0 = h0 + y(s) * f_tabT ( T0, s, h_tab )
    end do
    do it = 1, 50
      F  = 0.5d0 * a2 * T * T - h0
      dF = a2 * T
      do s = 1, nsc
        F  = F  + y(s) * f_tabT ( T, s, h_tab )
        dF = dF + y(s) * f_tabT ( T, s, cp_tab )
      end do
      dT = F / dF
      T  = T - dT
      if ( abs(dT) <= 1d-12 * T0 ) exit
    end do

  end function SRM_Static_Temperature


end module Q2D_Lib_BC_Fluxes_SRM
