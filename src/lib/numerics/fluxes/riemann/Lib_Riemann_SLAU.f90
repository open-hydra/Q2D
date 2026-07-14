module Q2D_Lib_Riemann_SLAU
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_SLAU, riemann_SLAU2

contains

  subroutine riemann_SLAU_base(dltot,drtot,ul,vl,ur,vr,pl,pr,al,ar,nx,ny, &
                                F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF)
    implicit none
    real(R8), intent(in)  :: dltot,drtot,ul,vl,ur,vr,pl,pr,al,ar,nx,ny
    real(R8), intent(out) :: F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF
    real(R8) :: unl, unr, mL, mR, Mach_tilde, f_rho
    real(R8) :: Vn_mag, Vn_magL, Vn_magR

    unl = ul*nx + vl*ny
    unr = ur*nx + vr*ny
    sq_veli = ul**2 + vl**2
    sq_velj = ur**2 + vr**2

    aF = 0.5d0*(al + ar)

    mL = unl / aF
    mR = unr / aF

    Mach_tilde = min(1.d0, sqrt(0.5d0*(sq_veli + sq_velj)) / aF)
    Chi = (1.d0 - Mach_tilde)**2

    Vn_mag  = (dltot*abs(unl) + drtot*abs(unr)) / (dltot + drtot)
    f_rho   = -max(min(mL,0.d0), -1.d0) * min(max(mR,0.d0), 1.d0)
    Vn_magL = (1.d0 - f_rho)*Vn_mag + f_rho*abs(unl)
    Vn_magR = (1.d0 - f_rho)*Vn_mag + f_rho*abs(unr)

    F_r = 0.5d0 * ( dltot*(unl + Vn_magL) + drtot*(unr - Vn_magR) - (Chi/aF)*(pr - pl) )

    if (abs(mL) < 1.d0) then
      BetaL = 0.25d0*(2.d0 - mL)*(mL + 1.d0)**2
    else if (mL >= 0.d0) then
      BetaL = 1.d0
    else
      BetaL = 0.d0
    endif

    if (abs(mR) < 1.d0) then
      BetaR = 0.25d0*(2.d0 + mR)*(mR - 1.d0)**2
    else if (mR >= 0.d0) then
      BetaR = 0.d0
    else
      BetaR = 1.d0
    endif

  end subroutine riemann_SLAU_base


  subroutine riemann_SLAU(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: BetaL, BetaR, Chi, sq_veli, sq_velj, aF, p_half, H_L, H_R
    real(R8) :: dummy_beta

    dummy_beta = beta

    call riemann_SLAU_base(dltot,drtot,ul,vl,ur,vr,pl,pr,al,ar,nx,ny, &
                           F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF)

    p_half = 0.5d0*(pl + pr) + 0.5d0*(BetaL - BetaR)*(pl - pr) &
           + (1.d0 - Chi)*(BetaL + BetaR - 1.d0)*0.5d0*(pl + pr)

    H_L = H0(pl, dl, sqrt(sq_veli))
    H_R = H0(pr, dr, sqrt(sq_velj))

    if (F_r >= 0.d0) then
      F_u = F_r*ul + nx*p_half
      F_v = F_r*vl + ny*p_half
      F_E = F_r*H_L
    else
      F_u = F_r*ur + nx*p_half
      F_v = F_r*vr + ny*p_half
      F_E = F_r*H_R
    endif

  end subroutine riemann_SLAU


  subroutine riemann_SLAU2(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: BetaL, BetaR, Chi, sq_veli, sq_velj, aF, p_half, H_L, H_R
    real(R8) :: dummy_beta

    dummy_beta = beta

    call riemann_SLAU_base(dltot,drtot,ul,vl,ur,vr,pl,pr,al,ar,nx,ny, &
                           F_r,BetaL,BetaR,Chi,sq_veli,sq_velj,aF)

    p_half = 0.5d0*(pl + pr) + 0.5d0*(BetaL - BetaR)*(pl - pr) &
           + sqrt(0.5d0*(sq_veli + sq_velj))*(BetaL + BetaR - 1.d0)*aF*0.5d0*(dltot + drtot)

    H_L = H0(pl, dl, sqrt(sq_veli))
    H_R = H0(pr, dr, sqrt(sq_velj))

    if (F_r >= 0.d0) then
      F_u = F_r*ul + nx*p_half
      F_v = F_r*vl + ny*p_half
      F_E = F_r*H_L
    else
      F_u = F_r*ur + nx*p_half
      F_v = F_r*vr + ny*p_half
      F_E = F_r*H_R
    endif

  end subroutine riemann_SLAU2

end module Q2D_Lib_Riemann_SLAU
