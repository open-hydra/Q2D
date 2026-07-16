module Q2D_Lib_Riemann_LF
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Riemann_PLLF, riemann_LLF

contains

  ! Local Lax-Friedrichs (alias Rusanov) - 2D
  subroutine riemann_LLF(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: E0
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: uln, urn, A
    real(R8) :: Frr, Fur, Fvr, Fer, Frl, Ful, Fvl, Fel
    real(R8) :: dummy_beta

    dummy_beta = beta

    uln = ul*nx + vl*ny
    urn = ur*nx + vr*ny

    A = MAX(ABS(uln-al),ABS(urn-ar),ABS(uln+al),ABS(urn+ar))

    call fluxes(pr,dr,ur,vr,nx,ny,Frr,Fur,Fvr,Fer)
    call fluxes(pl,dl,ul,vl,nx,ny,Frl,Ful,Fvl,Fel)

    F_r = 0.5d0*(Frl + Frr - A*(drtot - dltot))
    F_u = 0.5d0*(Fur + Ful - A*(drtot*ur - dltot*ul))
    F_v = 0.5d0*(Fvr + Fvl - A*(drtot*vr - dltot*vl))
    F_E = 0.5d0*(Fel + Fer - A*(drtot*E0(pr,dr,sqrt(ur*ur+vr*vr)) &
                              - dltot*E0(pl,dl,sqrt(ul*ul+vl*vl))))

  end subroutine riemann_LLF


  ! Preconditioned local Lax-Friedrichs - 2D
  subroutine Riemann_PLLF(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: f_Rtot, f_cp, E0
    implicit none
    real(R8), intent(in)  :: dl(nsc), ul, vl, pl, al
    real(R8), intent(in)  :: dr(nsc), ur, vr, pr, ar
    real(R8), intent(in)  :: dltot, drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: Rgasl, Rgasr
    real(R8) :: uln, urn, Tl, Tr, rho_Tl, rho_Tr, cpl, cpr, U_rl, U_rr
    real(R8) :: thetal, thetar, alphal, alphar, betal, betar, ul_mod, ur_mod
    real(R8) :: al_mod, ar_mod, A, Frr, Fur, Fvr, Fer, Frl, Ful, Fvl, Fel
    real(R8) :: dummy_beta

    dummy_beta = beta

    Rgasl = f_Rtot(dl)
    Rgasr = f_Rtot(dr)

    uln = ul*nx + vl*ny
    urn = ur*nx + vr*ny

    Tl = pl / ( dltot * Rgasl )
    Tr = pr / ( drtot * Rgasr )
    rho_Tl = - dltot / Tl
    rho_Tr = - drtot / Tr
    cpl = f_cp ( dl, Tl, dltot )
    cpr = f_cp ( dr, Tr, drtot )

    U_rl = max( 1d-5*al, min( abs(uln) , al ) )
    U_rr = max( 1d-5*ar, min( abs(urn) , ar ) )
    thetal = 1d0/U_rl**2 - rho_Tl/( dltot * cpl )
    thetar = 1d0/U_rr**2 - rho_Tr/( drtot * cpr )
    betal = 1d0/al**2
    betar = 1d0/ar**2
    alphal = 0.5d0*( 1d0 - betal*U_rl**2 )
    alphar = 0.5d0*( 1d0 - betar*U_rr**2 )
    ul_mod = abs( uln*( 1d0 - alphal ) )
    ur_mod = abs( urn*( 1d0 - alphar ) )
    al_mod = sqrt( alphal**2 * uln**2 + U_rl**2 )
    ar_mod = sqrt( alphar**2 * urn**2 + U_rr**2 )

    A = max( ul_mod + al_mod, ur_mod + ar_mod )

    call fluxes(pr,dr,ur,vr,nx,ny,Frr,Fur,Fvr,Fer)
    call fluxes(pl,dl,ul,vl,nx,ny,Frl,Ful,Fvl,Fel)

    F_r = 0.5d0*(Frl + Frr - A*(drtot - dltot))
    F_u = 0.5d0*(Fur + Ful - A*(drtot*ur - dltot*ul))
    F_v = 0.5d0*(Fvr + Fvl - A*(drtot*vr - dltot*vl))
    F_E = 0.5d0*(Fel + Fer - A*(drtot*E0(pr,dr,sqrt(ur*ur+vr*vr)) &
                              - dltot*E0(pl,dl,sqrt(ul*ul+vl*vl))))

  end subroutine Riemann_PLLF


  pure subroutine fluxes(p,r,u,v,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)  :: p
    real(R8), intent(in)  :: r(nsc)
    real(R8), intent(in)  :: u, v
    real(R8), intent(in)  :: nx, ny
    real(R8), intent(out) :: F_r
    real(R8), intent(out) :: F_u, F_v
    real(R8), intent(out) :: F_E

    F_r = sum(r)*(u*nx + v*ny)
    F_u = F_r*u + p*nx
    F_v = F_r*v + p*ny
    F_E = F_r*H0(p,r,sqrt(u*u + v*v))

  end subroutine fluxes

end module Q2D_Lib_Riemann_LF
