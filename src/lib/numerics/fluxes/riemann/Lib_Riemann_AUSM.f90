! Riemann solvers of the AUSM family (2D formulation)
! References:
!   Liou, M.-S., and Steffen, C. J. Jr. (1993), JCP, 107, 23-39.
!   Liou, M.-S. (1996), JCP, 129, 364-382.
!   Liou, M.-S. (2006), JCP, 214, 137-170.
!   Kitamura, K., Shima, E. (AUSM+-up2 hypersonic).
module Q2D_Lib_Riemann_AUSM
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Config_Types_m, only: obj_riemann

  implicit none
  private
  public :: riemann_Hanel, riemann_AUSMp
  public :: riemann_AUSMp_up, riemann_AUSMp_up2

contains

  subroutine riemann_Hanel(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: unl, unr, H_L, H_R, uLp, uRm, pLp, pRm, phalf
    real(R8) :: dummy_beta

    dummy_beta = beta

    unl = ul*nx + vl*ny
    unr = ur*nx + vr*ny

    H_L = H0(pl, dl, sqrt(ul*ul + vl*vl))
    H_R = H0(pr, dr, sqrt(ur*ur + vr*vr))

    if (abs(unl) <= al) then
      uLp = (unl + al)**2/(4.d0*al)
    else
      uLp = 0.5d0*(unl + abs(unl))
    endif

    if (abs(unr) <= ar) then
      uRm = (-(unr - ar)**2)/(4.d0*ar)
    else
      uRm = 0.5d0*(unr - abs(unr))
    endif

    if (abs(unl) <= al) then
      pLp = pl*uLp*(2d0 - unl/al)/al
    else
      pLp = pl*uLp/unl
    endif

    if (abs(unr) <= ar) then
      pRm = pr*uRm*(-2d0 - unr/ar)/ar
    else
      pRm = pr*uRm/unr
    endif

    phalf = pLp + pRm

    F_r = uLp*dltot + uRm*drtot

    F_u = uLp*(dltot*ul) + uRm*(drtot*ur) + phalf*nx
    F_v = uLp*(dltot*vl) + uRm*(drtot*vr) + phalf*ny

    F_E = uLp*(dltot*H_L) + uRm*(drtot*H_R)

  end subroutine riemann_Hanel


  subroutine riemann_AUSMp(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8), parameter :: betap = 0.125d0, alpha = 0.1875d0
    real(R8) :: unl, unr, ML, MR
    real(R8) :: MplusL, MminusR, Mhalf
    real(R8) :: pplusL, pminusR, phalf
    real(R8) :: H_L, H_R, a_half, mdotm, mdotp
    real(R8) :: dummy_beta

    dummy_beta = beta

    unl = ul*nx + vl*ny
    unr = ur*nx + vr*ny

    H_L = H0(pl, dl, sqrt(ul**2 + vl**2))
    H_R = H0(pr, dr, sqrt(ur**2 + vr**2))

    a_half = sqrt( al * ar )

    ML = unl / a_half
    MR = unr / a_half

    if (abs(ML) < 1.d0) then
      MplusL = 0.25d0*(ML + 1.d0)**2 + betap*(ML**2 - 1.d0)**2
      pplusL = 0.25d0*(ML + 1.d0)**2 * (2.d0 - ML) + alpha*ML*(ML**2 - 1.d0)**2
    else
      MplusL = 0.5d0*(ML + abs(ML))
      pplusL = 0.5d0*(1.d0 + sign(1.d0, ML))
    endif

    if (abs(MR) < 1.d0) then
      MminusR = -0.25d0*(MR - 1.d0)**2 - betap*(MR**2 - 1.d0)**2
      pminusR = 0.25d0*(MR - 1.d0)**2 * (2.d0 + MR) - alpha*MR*(MR**2 - 1.d0)**2
    else
      MminusR = 0.5d0*(MR - abs(MR))
      pminusR = 0.5d0*(1.d0 - sign(1.d0, MR))
    endif

    Mhalf = MplusL + MminusR
    phalf = pplusL*pl + pminusR*pr

    mdotp = dltot * a_half * 0.5d0 * (Mhalf + abs(Mhalf))
    mdotm = drtot * a_half * 0.5d0 * (Mhalf - abs(Mhalf))

    F_r = mdotp + mdotm
    F_u = mdotp * ul + mdotm * ur + nx * phalf
    F_v = mdotp * vl + mdotm * vr + ny * phalf
    F_E = mdotp * H_L + mdotm * H_R

  end subroutine riemann_AUSMp


  subroutine riemann_AUSMp_up(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0, f_gamma, f_Rtot
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8), parameter :: Kp = 0.25d0, Ku = 0.75d0, sigma = 1.0d0, betap = 0.125d0
    real(R8) :: unl, unr, H_L, H_R
    real(R8) :: gamma_L, gamma_R, Rgasl, Rgasr
    real(R8) :: astarL, astarR, astarL2, astarR2
    real(R8) :: atildeL, atildeR, aF
    real(R8) :: mL, mR, MFsq, Mrefsq, fa, alpha
    real(R8) :: mLP, mRM, betaLP, betaRM
    real(R8) :: rhoF, Mp, Pu, mF, mdot
    real(R8) :: upwind_u, upwind_v, upwind_H
    real(R8) :: phalf, p1, p2
    real(R8) :: dummy_beta

    dummy_beta = beta

    unl = ul*nx + vl*ny
    unr = ur*nx + vr*ny

    H_L = H0(pl, dl, sqrt(ul**2 + vl**2))
    H_R = H0(pr, dr, sqrt(ur**2 + vr**2))

    Rgasl = f_Rtot(dl)
    Rgasr = f_Rtot(dr)
    gamma_L = f_gamma(dl, pl, dltot, Rgasl)
    gamma_R = f_gamma(dr, pr, drtot, Rgasr)
    astarL2 = 2d0*(gamma_L - 1d0)/(gamma_L + 1.d0)*H_L
    astarR2 = 2d0*(gamma_R - 1d0)/(gamma_R + 1.d0)*H_R
    astarL = sqrt(astarL2)
    astarR = sqrt(astarR2)

    atildeL = astarL2 / max(astarL,  unl)
    atildeR = astarR2 / max(astarR, -unr)

    aF = min(atildeL, atildeR)

    mL = unl / aF
    mR = unr / aF

    MFsq = 0.5d0*(mL*mL + mR*mR)
    Mrefsq = min(1.0d0, max(MFsq, obj_riemann%Minf*obj_riemann%Minf))
    fa = 2.d0*sqrt(Mrefsq) - Mrefsq
    alpha = 0.1875d0 * (-4.d0 + 5.d0 * fa * fa)

    if (abs(mL) <= 1.d0) then
      p1 = 0.25d0 * (mL + 1.d0)**2
      p2 = (mL*mL - 1.d0)**2
      mLP    = p1 + betap * p2
      betaLP = p1 * (2.d0 - mL) + alpha * mL * p2
    else
      mLP    = 0.5d0 * (mL + abs(mL))
      betaLP = mLP/mL
    endif

    if (abs(mR) <= 1.d0) then
      p1 = 0.25d0 * (mR - 1.d0)**2
      p2 = (mR*mR - 1.d0)**2
      mRM    = -p1 - betap * p2
      betaRM = p1 * (2.d0 + mR) - alpha * mR * p2
    else
      mRM    = 0.5d0 * (mR - abs(mR))
      betaRM = mRM/mR
    endif

    rhoF = 0.5d0 * (dltot + drtot)
    Mp = - (Kp / fa) * max(0.d0, 1.d0 - sigma * MFsq) * (pr - pl) / (rhoF * aF * aF)
    Pu = - Ku * fa * betaLP * betaRM * 2.d0 * rhoF * aF * (unr - unl)

    mF = mLP + mRM + Mp
    mdot = aF * ( max(mF, 0.d0) * dltot + min(mF, 0.d0) * drtot )

    phalf = betaLP * pl + betaRM * pr + Pu

    if (mF >= 0.d0) then
      upwind_u = ul
      upwind_v = vl
      upwind_H = H_L
    else
      upwind_u = ur
      upwind_v = vr
      upwind_H = H_R
    endif

    F_r = mdot
    F_u = mdot * upwind_u + nx * phalf
    F_v = mdot * upwind_v + ny * phalf
    F_E = mdot * upwind_H

  end subroutine riemann_AUSMp_up


  subroutine riemann_AUSMp_up2(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: H0, f_gamma, f_Rtot
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al,dltot
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8), parameter :: Kp = 0.25d0, sigma = 1.0d0, betap = 0.125d0
    real(R8) :: unl, unr, sq_vel
    real(R8) :: H_L, H_R
    real(R8) :: gamma_L, gamma_R, Rgasl, Rgasr
    real(R8) :: astarL2, astarR2, astarL, astarR
    real(R8) :: ahatL, ahatR, aF
    real(R8) :: mL, mR, MFsq, param1, Mrefsq, fa, alpha
    real(R8) :: mLP, mRM, pLP, pRM
    real(R8) :: rhoF, Mp, mF, mdot, phalf
    real(R8) :: upwind_u, upwind_v, upwind_H
    real(R8) :: p1, p2
    real(R8) :: dummy_beta

    dummy_beta = beta

    unl = ul*nx + vl*ny
    unr = ur*nx + vr*ny
    sq_vel = 0.5d0*((ul**2 + vl**2) + (ur**2 + vr**2))

    H_L = H0(pl, dl, sqrt(ul**2 + vl**2))
    H_R = H0(pr, dr, sqrt(ur**2 + vr**2))

    Rgasl = f_Rtot(dl)
    Rgasr = f_Rtot(dr)
    gamma_L = f_gamma(dl, pl, dltot, Rgasl)
    gamma_R = f_gamma(dr, pr, drtot, Rgasr)
    astarL2 = 2d0*(gamma_L - 1d0)/(gamma_L + 1.d0)*H_L
    astarR2 = 2d0*(gamma_R - 1d0)/(gamma_R + 1.d0)*H_R
    astarL = sqrt(astarL2)
    astarR = sqrt(astarR2)

    ahatL = astarL2 / max(astarL, unl)
    ahatR = astarR2 / max(astarR, -unr)
    aF = min(ahatL, ahatR)

    mL = unl / aF
    mR = unr / aF
    MFsq = 0.5d0*(mL*mL + mR*mR)
    param1 = max(MFsq, obj_riemann%Minf*obj_riemann%Minf)
    Mrefsq = min(1.0d0, param1)
    fa = 2.d0*sqrt(Mrefsq) - Mrefsq
    alpha = 3.d0/16.d0*(-4.d0 + 5.d0*fa*fa)

    if (abs(mL) <= 1.d0) then
      p1 = 0.25d0*(mL+1.d0)**2
      p2 = (mL*mL-1.d0)**2
      mLP = p1 + betap*p2
      pLP = p1*(2.d0 - mL) + alpha*mL*p2
    else
      mLP = 0.5d0*(mL + abs(mL))
      pLP = mLP / mL
    endif

    if (abs(mR) <= 1.d0) then
      p1 = 0.25d0*(mR-1.d0)**2
      p2 = (mR*mR-1.d0)**2
      mRM = -p1 - betap*p2
      pRM =  p1*(2.d0 + mR) - alpha*mR*p2
    else
      mRM = 0.5d0*(mR - abs(mR))
      pRM = mRM / mR
    endif

    rhoF = 0.5d0*(dltot + drtot)
    Mp = -(Kp/fa)*max(0.d0, 1.d0 - sigma*MFsq)*(pr - pl)/(rhoF*aF*aF)
    mF = mLP + mRM + Mp
    mdot = aF*( max(mF,0.d0)*dltot + min(mF,0.d0)*drtot )

    phalf = 0.5d0*(pr + pl) + 0.5d0*(pLP - pRM)*(pl - pr) + sqrt(sq_vel)*(pLP + pRM - 1.d0)*rhoF*aF

    if (mF >= 0.d0) then
      upwind_u = ul
      upwind_v = vl
      upwind_H = H_L
    else
      upwind_u = ur
      upwind_v = vr
      upwind_H = H_R
    endif

    F_r = mdot
    F_u = mdot * upwind_u + nx * phalf
    F_v = mdot * upwind_v + ny * phalf
    F_E = mdot * upwind_H

  end subroutine riemann_AUSMp_up2

end module Q2D_Lib_Riemann_AUSM
