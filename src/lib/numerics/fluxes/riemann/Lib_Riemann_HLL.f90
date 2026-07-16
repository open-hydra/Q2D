!>@brief: Riemann solvers of the HLL family (2D formulation).
!
! 7 public solvers + 3 helper routines; in 2D the face tangent is
! unique: `t = (-ny,nx)`.
!
! Public solvers:
!   - riemann_HLLE   (Einfeldt-Roe HLL flux)
!   - riemann_HLLEM  (HLLEM = HLL + antidiffusion on contact/shear, Einfeldt 1991)
!   - riemann_HLLC   (Batten intermediate star wave, HLLC formulation)
!   - riemann_HLLEpp (HLLE++ blend with Roe in the Tramel/OVERFLOW style, Tramel 2009)
!   - riemann_HLLEMSD (HLLEM+HLLE blend driven by the shock detector)
!   - riemann_HLLCSD  (HLLC+HLLE blend driven by the shock detector)
!   - riemann_HLLCHLLE (rotated HLLC/HLLE for faces not aligned with the discontinuity)
!
! Reference:
!   Tramel, R., Nichols, R., & Buning, P. (2009). Addition of improved
!   shock-capturing schemes to OVERFLOW 2.1. AIAA-2009-3988.
!
! HLLEMSD/HLLCSD are the main consumers of the `beta` computed by
! Lib_Shock_Detector: without them the shock detector remains unused.
!
module Q2D_Lib_Riemann_HLL
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_HLLE, riemann_HLLEM, riemann_HLLC
  public :: riemann_HLLEpp
  public :: riemann_HLLEMSD, riemann_HLLCSD
  public :: riemann_HLLCHLLE

contains

  subroutine riemann_HLLC(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    implicit none
    real(R8), intent(in) :: dl(nsc), dr(nsc)
    real(R8), intent(in) :: ul, vl, ur, vr
    real(R8), intent(in) :: pl, pr, al, ar
    real(R8), intent(in) :: dltot, drtot
    real(R8), intent(in) :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: rho_roe, u_roe, v_roe, a_roe, un_roe, h0_roe
    real(R8) :: h0l, h0r, unl, unr
    real(R8) :: S1, S4, Sstar, pstar, inv_denom
    real(R8) :: U1S, U2S, U3S, U5S
    real(R8) :: e0l, e0r


    unl = ul*nx + vl*ny
    unr = ur*nx + vr*ny

    call roe_averages(nsc, dl, dr, dltot, drtot, ul, vl, ur, vr, pl, pr, &
                      h0l, h0r, rho_roe, u_roe, v_roe, a_roe, h0_roe)

    un_roe = u_roe*nx + v_roe*ny

    e0l = h0l - pl/dltot
    e0r = h0r - pr/drtot

    S1 = min(0.d0, unl - al, un_roe - a_roe)
    S4 = max(0.d0, unr + ar, un_roe + a_roe)

    inv_denom = 1.d0 / ( dltot*(S1-unl) - drtot*(S4-unr) )

    Sstar = ( pr - pl + dltot*unl*(S1-unl) - drtot*unr*(S4-unr) ) * inv_denom
    pstar = dltot*(unl-S1)*(unl-Sstar) + pl

    if (S1 >= 0.d0) then

      F_r = dltot * unl
      F_u = F_r*ul + pl*nx
      F_v = F_r*vl + pl*ny
      F_E = F_r*h0l

    elseif (Sstar >= 0.d0) then

      F_r = dltot * unl
      F_u = F_r*ul + pl*nx
      F_v = F_r*vl + pl*ny
      F_E = F_r*h0l

      inv_denom = 1.d0 / (S1 - Sstar)

      U1S = dltot*(S1-unl)*inv_denom
      U2S = ((S1-unl)*dltot*ul + (pstar-pl)*nx)*inv_denom
      U3S = ((S1-unl)*dltot*vl + (pstar-pl)*ny)*inv_denom
      U5S = ((S1-unl)*dltot*e0l - pl*unl + pstar*Sstar)*inv_denom

      F_r = F_r + S1*(U1S - dltot)
      F_u = F_u + S1*(U2S - dltot*ul)
      F_v = F_v + S1*(U3S - dltot*vl)
      F_E = F_E + S1*(U5S - dltot*e0l)

    elseif (S4 > 0.d0) then

      F_r = drtot * unr
      F_u = F_r*ur + pr*nx
      F_v = F_r*vr + pr*ny
      F_E = F_r*h0r

      inv_denom = 1.d0 / (S4 - Sstar)

      U1S = drtot*(S4-unr)*inv_denom
      U2S = ((S4-unr)*drtot*ur + (pstar-pr)*nx)*inv_denom
      U3S = ((S4-unr)*drtot*vr + (pstar-pr)*ny)*inv_denom
      U5S = ((S4-unr)*drtot*e0r - pr*unr + pstar*Sstar)*inv_denom

      F_r = F_r + S4*(U1S - drtot)
      F_u = F_u + S4*(U2S - drtot*ur)
      F_v = F_v + S4*(U3S - drtot*vr)
      F_E = F_E + S4*(U5S - drtot*e0r)

    else

      F_r = drtot * unr
      F_u = F_r*ur + pr*nx
      F_v = F_r*vr + pr*ny
      F_E = F_r*h0r

    endif

  end subroutine riemann_HLLC


  subroutine riemann_HLLE(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    implicit none
    real(R8), intent(in) :: dl(nsc), dr(nsc)
    real(R8), intent(in) :: ul, vl, ur, vr
    real(R8), intent(in) :: pl, pr, al, ar
    real(R8), intent(in) :: dltot, drtot
    real(R8), intent(in) :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: rho_roe, u_roe, v_roe, a_roe, un_roe, h0_roe
    real(R8) :: h0l, h0r, unl, unr
    real(R8) :: S1, S4, inv_denom
    real(R8) :: FrL, FrR


    unl = ul*nx + vl*ny
    unr = ur*nx + vr*ny

    call roe_averages(nsc, dl, dr, dltot, drtot, ul, vl, ur, vr, pl, pr, &
                      h0l, h0r, rho_roe, u_roe, v_roe, a_roe, h0_roe)

    un_roe = u_roe*nx + v_roe*ny

    S1 = min(0.d0, unl - al, un_roe - a_roe)
    S4 = max(0.d0, unr + ar, un_roe + a_roe)

    inv_denom = 1.d0 / (S4 - S1)

    FrL = dltot*unl
    FrR = drtot*unr

    F_r = ( S4*FrL - S1*FrR + S1*S4*(drtot-dltot) ) * inv_denom

    F_u = ( S4*(FrL*ul + pl*nx) - S1*(FrR*ur + pr*nx) &
          + S1*S4*(drtot*ur - dltot*ul) ) * inv_denom

    F_v = ( S4*(FrL*vl + pl*ny) - S1*(FrR*vr + pr*ny) &
          + S1*S4*(drtot*vr - dltot*vl) ) * inv_denom

    F_E = ( S4*(FrL*h0l) - S1*(FrR*h0r) &
          + S1*S4*( drtot*(h0r - pr/drtot) &
                  - dltot*(h0l - pl/dltot) ) ) * inv_denom

  end subroutine riemann_HLLE


  subroutine riemann_HLLCHLLE(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: nx1, ny1, nx2, ny2
    real(R8) :: alfa1, alfa2
    real(R8) :: F_rHLLE, F_uHLLE, F_vHLLE, F_EHLLE
    real(R8) :: F_rHLLC, F_uHLLC, F_vHLLC, F_EHLLC
    real(R8) :: abs_dv, dum

    abs_dv = sqrt((ur-ul)**2 + (vr-vl)**2)
    dum = 0.5d0*( sqrt(ul**2+vl**2) + sqrt(ur**2+vr**2) )

    if (abs_dv/(dum + 1.0d-10) > 1.0d-12) then

      nx1 = (ur-ul) / abs_dv
      ny1 = (vr-vl) / abs_dv

      alfa1 = nx*nx1 + ny*ny1

      nx1 = sign(1.d0, alfa1)*nx1
      ny1 = sign(1.d0, alfa1)*ny1
      alfa1 = sign(1.d0, alfa1)*alfa1

      nx2 = -ny1
      ny2 =  nx1
      alfa2 = nx*nx2 + ny*ny2

      nx2 = sign(1.d0, alfa2)*nx2
      ny2 = sign(1.d0, alfa2)*ny2
      alfa2 = sign(1.d0, alfa2)*alfa2

      call riemann_HLLE(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx1,ny1,F_rHLLE,F_uHLLE,F_vHLLE,F_EHLLE)
      call riemann_HLLC(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx2,ny2,F_rHLLC,F_uHLLC,F_vHLLC,F_EHLLC)

      F_r = F_rHLLE*alfa1 + F_rHLLC*alfa2
      F_u = F_uHLLE*alfa1 + F_uHLLC*alfa2
      F_v = F_vHLLE*alfa1 + F_vHLLC*alfa2
      F_E = F_EHLLE*alfa1 + F_EHLLC*alfa2

    else

      call riemann_HLLC(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_rHLLC,F_uHLLC,F_vHLLC,F_EHLLC)
      F_r = F_rHLLC
      F_u = F_uHLLC
      F_v = F_vHLLC
      F_E = F_EHLLC

    endif

  end subroutine riemann_HLLCHLLE


  pure subroutine roe_averages(nsc, dl, dr, dltot, drtot, ul, vl, ur, vr, pl, pr, &
                                h0l, h0r, rho_roe, u_roe, v_roe, a_roe, h0_roe)
    use FLINT_Lib_Thermodynamic, only: f_Rtot, h_tab, cp_tab, Wm_tab, Runiv
    implicit none
    integer, intent(in)   :: nsc
    real(R8), intent(in)  :: dl(nsc), dr(nsc)
    real(R8), intent(in)  :: dltot, drtot
    real(R8), intent(in)  :: ul, vl, ur, vr
    real(R8), intent(in)  :: pl, pr
    real(R8), intent(out) :: h0l, h0r
    real(R8), intent(out) :: rho_roe, u_roe, v_roe, a_roe, h0_roe
    integer  :: s, Til, Tir
    real(R8) :: Rl, Rr, Tl, Tr, dTl, dTr
    real(R8) :: inv_dltot, inv_drtot
    real(R8) :: srL, srR, inv_sr
    real(R8) :: cv_roe, sum_ei
    real(R8) :: R_roe, T_roe
    real(R8) :: gam_roe, vel2
    real(R8) :: hl, hr, el, er, invW, d_roe

    Rl = f_Rtot(dl)
    Rr = f_Rtot(dr)

    inv_dltot = 1.d0 / dltot
    inv_drtot = 1.d0 / drtot

    Tl = pl * inv_dltot / Rl
    Tr = pr * inv_drtot / Rr

    Til = int(Tl)
    Tir = int(Tr)
    dTl = Tl - Til
    dTr = Tr - Tir

    srL = sqrt(dltot)
    srR = sqrt(drtot)
    inv_sr = 1.d0 / (srL + srR)

    rho_roe = srL * srR

    u_roe = (srR*ur + srL*ul) * inv_sr
    v_roe = (srR*vr + srL*vl) * inv_sr

    R_roe = (srR*Rr + srL*Rl) * inv_sr
    T_roe = (srR*Tr + srL*Tl) * inv_sr

    h0l = 0.d0
    h0r = 0.d0
    cv_roe = 0.d0
    sum_ei = 0.d0

    do s = 1, nsc

      invW = Runiv / Wm_tab(s)

      hl = h_tab(Til, s) + (h_tab(Til+1, s) - h_tab(Til, s)) * dTl
      hr = h_tab(Tir, s) + (h_tab(Tir+1, s) - h_tab(Tir, s)) * dTr

      el = hl - invW * Tl
      er = hr - invW * Tr

      h0l = h0l + hl * dl(s) * inv_dltot
      h0r = h0r + hr * dr(s) * inv_drtot

      d_roe = (srR*dr(s)*inv_drtot + srL*dl(s)*inv_dltot) * inv_sr

      cv_roe = cv_roe + d_roe * ( &
          0.5d0 * ( &
            cp_tab(Til, s) + (cp_tab(Til+1, s)-cp_tab(Til, s))*dTl + &
            cp_tab(Tir, s) + (cp_tab(Tir+1, s)-cp_tab(Tir, s))*dTr ) &
          - invW )

      sum_ei = sum_ei + d_roe * (srR*er + srL*el) * inv_sr

    end do

    h0l = h0l + 0.5d0*(ul*ul + vl*vl)
    h0r = h0r + 0.5d0*(ur*ur + vr*vr)

    h0_roe = (srR*h0r + srL*h0l) * inv_sr

    gam_roe = 1.d0 + R_roe / cv_roe

    vel2 = u_roe*u_roe + v_roe*v_roe

    a_roe = sqrt( (gam_roe - 1.d0) * &
          ( h0_roe - 0.5d0*vel2 + cv_roe*T_roe - sum_ei ) )

  end subroutine roe_averages


  pure function hartenHyman(v,e) result(l)
    ! Harten-Hyman entropy fix used by riemann_HLLEpp.
    implicit none
    real(R8), intent(in) :: v, e
    real(R8) :: l
    real(R8) :: al

    al = abs(v)
    l = al
    if (al<=e) l = (v*v + e*e)/(2d0*e)

  end function hartenHyman


  pure subroutine fluxes_helper ( p, r, u, v, nx, ny, F_r, F_u, F_v, F_E )
    ! Inviscid flux in conservative form from the primitives (2D).
    use Q2D_Global_m,            only: nsc
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
    F_E = F_r*H0(p, r, sqrt(u*u + v*v))

  end subroutine fluxes_helper


  subroutine riemann_HLLEM ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                             beta, nx, ny, F_r, F_u, F_v, F_E )
    ! HLLEM (Einfeldt et al. 1991): HLL + antidiffusion on contact and shear.
    ! 2D formulation: unique tangent `t = (-ny, nx)`, a single shear wave.
    ! `beta` is ignored (HLLEM does not use the shock detector; the argument exists
    ! only for uniformity with the common abstract interface).
    use Q2D_Global_m, only: nsc
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al
    real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    ! Local variables
    real(R8) :: h0r, h0l
    real(R8) :: S1, S4, C, delta, inv_denom
    real(R8) :: Frl, Ful, Fvl, Fel
    real(R8) :: Frr, Fur, Fvr, Fer
    real(R8) :: rho_ROE, u_ROE, v_ROE, a_ROE, un_ROE, h0_roe
    real(R8) :: drho, dp, du, dv
    real(R8) :: tx, ty
    real(R8) :: alpha_c, alpha_s
    real(R8) :: e0l_cons, e0r_cons

    ! Roe averages
    call roe_averages ( nsc, dl, dr, dltot, drtot, ul, vl, ur, vr, pl, pr, &
                        h0l, h0r, rho_roe, u_roe, v_roe, a_roe, h0_roe )
    un_roe = u_roe*nx + v_roe*ny

    ! Wave-speed estimates (Einfeldt limiting)
    S1 = min(0.d0, (ul*nx + vl*ny) - al, un_ROE - a_ROE)
    S4 = max(0.d0, (ur*nx + vr*ny) + ar, un_ROE + a_ROE)

    ! Left/right inviscid fluxes
    call fluxes_helper ( pl, dl, ul, vl, nx, ny, Frl, Ful, Fvl, Fel )
    call fluxes_helper ( pr, dr, ur, vr, nx, ny, Frr, Fur, Fvr, Fer )

    ! HLL flux (energy via the conservative variable e0 = h0 - p/rho).
    e0l_cons = h0l - pl/dltot
    e0r_cons = h0r - pr/drtot
    inv_denom = 1.d0 / (S4 - S1)

    F_r = ( S4*Frl - S1*Frr + S1*S4*(drtot - dltot)         ) * inv_denom
    F_u = ( S4*Ful - S1*Fur + S1*S4*(drtot*ur - dltot*ul)   ) * inv_denom
    F_v = ( S4*Fvl - S1*Fvr + S1*S4*(drtot*vr - dltot*vl)   ) * inv_denom
    F_E = ( S4*Fel - S1*Fer + S1*S4*(drtot*e0r_cons - dltot*e0l_cons) ) * inv_denom

    ! HLLEM antidiffusion (only in the subsonic range S1 < 0 < S4)
    if (S1 < 0.d0 .and. S4 > 0.d0) then
      drho = drtot - dltot
      dp   = pr - pl
      du   = ur - ul
      dv   = vr - vl

      ! 2D tangent (unique): perpendicular to (nx,ny), unit norm.
      tx = -ny
      ty =  nx

      ! Wave strengths: contact + one shear wave (in 2D the shear is unique)
      alpha_c = drho - dp/(a_ROE*a_ROE)
      alpha_s = rho_ROE*(du*tx + dv*ty)

      ! Antidiffusion coefficient
      delta = a_ROE / (a_ROE + abs(0.5d0*(S1+S4)))
      C     = (S1*S4) * inv_denom

      ! Mass
      F_r = F_r - C*delta*alpha_c
      ! Momentum: contact and shear contributions
      F_u = F_u - C*delta*( alpha_c*u_ROE + alpha_s*tx )
      F_v = F_v - C*delta*( alpha_c*v_ROE + alpha_s*ty )
      ! Energy: contact and shear contributions
      F_E = F_E - C*delta*( alpha_c * 0.5d0 * (u_ROE*u_ROE + v_ROE*v_ROE) &
                          + alpha_s * (u_ROE*tx + v_ROE*ty) )
    end if

    ! Avoid the unused-argument warning on beta (ignored by HLLEM by construction).
    associate (beta_unused => beta); end associate

  end subroutine riemann_HLLEM


  subroutine riemann_HLLEpp ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                              beta, nx, ny, F_r, F_u, F_v, F_E )
    ! HLLE++ in the Tramel/OVERFLOW style (Tramel et al. 2009): HLLE+Roe blend driven
    ! by the shock-detector beta.
    ! 2D formulation: unique tangent, a single shear wave.
    ! Note: `lambda_1 = un_roe` (not max(un_roe,a_roe)) to preserve the shear layers.
    use Q2D_Global_m, only: nsc
    implicit none
    real(R8), intent(in)  :: dl(nsc), ul, vl, pl, al
    real(R8), intent(in)  :: dr(nsc), ur, vr, pr, ar
    real(R8), intent(in)  :: dltot, drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    ! Local variables
    real(R8) :: Sm, Sp, unl, unr
    real(R8) :: Frl, Ful, Fvl, Fel
    real(R8) :: Frr, Fur, Fvr, Fer
    real(R8) :: h0r, h0l
    real(R8) :: inv_denom
    real(R8) :: rho_ROE, u_ROE, v_ROE, a_ROE, un_ROE, h0_roe
    real(R8) :: drho, dp, du, dv, dun
    real(R8) :: tx, ty
    real(R8) :: alpha_c, alpha_s, alpha_m, alpha_p
    real(R8) :: lam1_roe, lam2_roe, lam3_roe, eps
    real(R8) :: lam1_hlle, lam2_hlle, lam3_hlle
    real(R8) :: lam1_pp, lam2_pp, lam3_pp

    unl = ul*nx + vl*ny
    unr = ur*nx + vr*ny

    ! Roe averages
    call roe_averages ( nsc, dl, dr, dltot, drtot, ul, vl, ur, vr, pl, pr, &
                        h0l, h0r, rho_roe, u_roe, v_roe, a_roe, h0_roe )
    un_roe = u_roe*nx + v_roe*ny

    ! HLLE eigenvalues
    Sm = min(0.d0, unl - al, un_ROE - a_ROE)
    Sp = max(0.d0, unr + ar, un_ROE + a_ROE)
    inv_denom = 1.d0 / (Sp - Sm)
    lam1_hlle = abs( (Sp+Sm)*un_roe        *inv_denom - 2d0*Sp*Sm*inv_denom )
    lam2_hlle = abs( (Sp+Sm)*(un_roe-a_roe)*inv_denom - 2d0*Sp*Sm*inv_denom )
    lam3_hlle = abs( (Sp+Sm)*(un_roe+a_roe)*inv_denom - 2d0*Sp*Sm*inv_denom )

    ! Roe eigenvalues (with Harten-Hyman entropy fix on the acoustic waves)
    eps = 0.3d0*(abs(un_ROE) + a_roe)
    lam1_roe = abs(un_ROE)                            ! preserves the shear (no max with a_ROE)
    lam2_roe = hartenHyman(un_roe - a_roe, eps)
    lam3_roe = hartenHyman(un_roe + a_roe, eps)

    ! HLLE++ blend via beta (weight given by the shock detector)
    lam1_pp = beta*lam1_roe + (1.d0 - beta)*lam1_hlle
    lam2_pp = beta*lam2_roe + (1.d0 - beta)*lam2_hlle
    lam3_pp = beta*lam3_roe + (1.d0 - beta)*lam3_hlle

    ! Left/right inviscid fluxes
    call fluxes_helper ( pl, dl, ul, vl, nx, ny, Frl, Ful, Fvl, Fel )
    call fluxes_helper ( pr, dr, ur, vr, nx, ny, Frr, Fur, Fvr, Fer )

    ! Centered average
    F_r = 0.5d0*(Frl + Frr)
    F_u = 0.5d0*(Ful + Fur)
    F_v = 0.5d0*(Fvl + Fvr)
    F_E = 0.5d0*(Fel + Fer)

    ! Diffusive correction (jumps + waves)
    drho = drtot - dltot
    dp   = pr - pl
    du   = ur - ul
    dv   = vr - vl
    dun  = unr - unl

    ! 2D tangent
    tx = -ny
    ty =  nx

    ! Wave strengths: 1 contact + 1 shear + 2 acoustic (in 2D the shear is unique)
    alpha_m = (dp - rho_ROE*a_ROE*dun) / (2d0*a_ROE*a_ROE)
    alpha_c = drho - dp/(a_ROE*a_ROE)
    alpha_s = rho_ROE*(du*tx + dv*ty)
    alpha_p = (dp + rho_ROE*a_ROE*dun) / (2d0*a_ROE*a_ROE)

    ! Mass
    F_r = F_r - 0.5d0*( alpha_c*lam1_pp + alpha_m*lam2_pp + alpha_p*lam3_pp )

    ! Momentum: contact + shear + acoustic
    F_u = F_u - 0.5d0*( alpha_c*lam1_pp*u_ROE + alpha_s*lam1_pp*tx                       &
                      + alpha_m*lam2_pp*(u_roe - a_ROE*nx) + alpha_p*lam3_pp*(u_roe + a_ROE*nx) )
    F_v = F_v - 0.5d0*( alpha_c*lam1_pp*v_ROE + alpha_s*lam1_pp*ty                       &
                      + alpha_m*lam2_pp*(v_roe - a_ROE*ny) + alpha_p*lam3_pp*(v_roe + a_ROE*ny) )

    ! Energy: contact + shear + acoustic
    F_E = F_E - 0.5d0*( alpha_c*lam1_pp*0.5d0*(u_ROE*u_ROE + v_ROE*v_ROE)                &
                      + alpha_s*lam1_pp*(u_ROE*tx + v_ROE*ty)                            &
                      + alpha_m*lam2_pp*(h0_roe - a_ROE*un_ROE)                          &
                      + alpha_p*lam3_pp*(h0_roe + a_ROE*un_ROE) )

  end subroutine riemann_HLLEpp


  subroutine riemann_HLLEMSD ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                               beta, nx, ny, F_r, F_u, F_v, F_E )
    ! Hybrid HLLEM + HLLE with blending via the shock-detector beta.
    use Q2D_Global_m, only: nsc
    implicit none
    real(R8), intent(in)  :: dl(nsc), ul, vl, pl, al
    real(R8), intent(in)  :: dr(nsc), ur, vr, pr, ar
    real(R8), intent(in)  :: dltot, drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: FrE, FuE, FvE, FeE
    real(R8) :: FrC, FuC, FvC, FeC
    real(R8) :: beta_

    beta_ = beta

    if (beta_ == 0d0) then
      call riemann_HLLE  ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                           0d0, nx, ny, F_r, F_u, F_v, F_E )
    else if (beta_ >= 1d0) then
      call riemann_HLLEM ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                           0d0, nx, ny, F_r, F_u, F_v, F_E )
    else
      call riemann_HLLE  ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                           0d0, nx, ny, FrE, FuE, FvE, FeE )
      call riemann_HLLEM ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                           0d0, nx, ny, FrC, FuC, FvC, FeC )
      F_r = beta_*FrC + (1d0-beta_)*FrE
      F_u = beta_*FuC + (1d0-beta_)*FuE
      F_v = beta_*FvC + (1d0-beta_)*FvE
      F_E = beta_*FeC + (1d0-beta_)*FeE
    end if

  end subroutine riemann_HLLEMSD


  subroutine riemann_HLLCSD ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                              beta, nx, ny, F_r, F_u, F_v, F_E )
    ! Hybrid HLLC + HLLE with blending via the shock-detector beta (Tramel HLLC+).
    use Q2D_Global_m, only: nsc
    implicit none
    real(R8), intent(in)  :: dl(nsc), ul, vl, pl, al
    real(R8), intent(in)  :: dr(nsc), ur, vr, pr, ar
    real(R8), intent(in)  :: dltot, drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: FrE, FuE, FvE, FeE
    real(R8) :: FrC, FuC, FvC, FeC
    real(R8) :: beta_

    beta_ = beta

    if (beta_ == 0d0) then
      call riemann_HLLE ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                          0d0, nx, ny, F_r, F_u, F_v, F_E )
    else if (beta_ >= 1d0) then
      call riemann_HLLC ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                          0d0, nx, ny, F_r, F_u, F_v, F_E )
    else
      call riemann_HLLE ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                          0d0, nx, ny, FrE, FuE, FvE, FeE )
      call riemann_HLLC ( dl, ul, vl, pl, al, dltot, dr, ur, vr, pr, ar, drtot, &
                          0d0, nx, ny, FrC, FuC, FvC, FeC )
      F_r = beta_*FrC + (1d0-beta_)*FrE
      F_u = beta_*FuC + (1d0-beta_)*FuE
      F_v = beta_*FvC + (1d0-beta_)*FvE
      F_E = beta_*FeC + (1d0-beta_)*FeE
    end if

  end subroutine riemann_HLLCSD

end module Q2D_Lib_Riemann_HLL
