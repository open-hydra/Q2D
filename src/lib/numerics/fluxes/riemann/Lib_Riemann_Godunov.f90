module Q2D_Lib_Riemann_Godunov
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_exact

contains

  ! Exact Godunov solver (generalized Riemann invariants, Glaister 1988):
  ! the flux energy uses H0 with the full velocity magnitude
  ! sqrt(uf**2+vf**2) (Glaister 1988 §3) through the helper routine `fluxes`,
  ! not with the normal component alone.
  subroutine riemann_exact(dl,ul,vl,p1,a1,dltot,dr,ur,vr,p4,a4,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
    use Q2D_Global_m, only: nsc
    use FLINT_Lib_Thermodynamic, only: f_Rtot, f_gamma
    implicit none
    real(R8), intent(in)  :: dl(nsc),ul,vl,p1,a1
    real(R8), intent(in)  :: dr(nsc),ur,vr,p4,a4
    real(R8), intent(in)  :: dltot,drtot
    real(R8), intent(in)  :: beta, nx, ny
    real(R8), intent(out) :: F_r, F_u, F_v, F_E
    real(R8) :: Rgasl, Rgasr
    integer  :: iter
    real(R8) :: gam1, gam4, gm1_1, gm1_4, delta1, delta4, eta1, eta4, Rgas
    real(R8) :: un1, un4, ut1, ut4, vt1, vt4
    real(R8) :: R2I, R1IV, z, p2, p3, u0, vel, dpr, dpl, a2, a3
    real(R8) :: vw1, vw2, vw3, vw4, dum, sw_l, sw_r
    real(R8) :: af, pf, df(nsc), rhof, uf, vf
    real(R8) :: dummy_beta

    dummy_beta = beta

    Rgasl = f_Rtot(dl)
    Rgasr = f_Rtot(dr)
    gam1 = f_gamma(dl,p1,dltot,Rgasl)
    gam4 = f_gamma(dr,p4,drtot,Rgasr)

    gm1_1 = gam1 - 1.d0
    gm1_4 = gam4 - 1.d0

    delta1 = 0.5d0*gm1_1
    delta4 = 0.5d0*gm1_4

    eta1 = gam1/delta1
    eta4 = gam4/delta4

    un1 = ul*nx + vl*ny
    ut1 = ul - un1*nx; vt1 = vl - un1*ny
    un4 = ur*nx + vr*ny
    ut4 = ur - un4*nx; vt4 = vr - un4*ny

    R2I  = un1 + (2.d0/gm1_1)*a1
    R1IV = un4 - (2.d0/gm1_4)*a4

    if (p1 < p4) then
      dum = gm1_4/(2.d0*gam4)
    else
      dum = gm1_1/(2.d0*gam1)
    end if

    z = (gm1_1/gm1_4)*(a4/a1)*((p1/p4)**dum)
    u0 = (z*R2I + R1IV)/(1.d0 + z)
    vel = u0

    p3 = 1d0
    p2 = 2d0*p3
    iter = 0
    do while ( dabs(1 - (p2/p3)) >= 1d-5 .and. iter <= 1000 )
      if (vel < un1) then
        call shockleft(un1,a1,p1,vel,a2,p2,dpr,sw_l,gam1)
      else
        call rareleft(p1,a1,un1,vel,p2,a2,dpr,vw1,vw2,gam1)
      endif
      if (vel > un4) then
        call shockright(un4,a4,p4,vel,a3,p3,dpl,sw_r,gam4)
      else
        call rareright(p4,a4,un4,vel,p3,a3,dpl,vw3,vw4,gam4)
      endif
      vel = vel - ((p2 - p3)/(dpr - dpl))
      iter = iter + 1
    end do
    if (iter > 1000) write(*,*) "Warning: max NR iter reached (Riemann solver)"

    if (vel > 0.d0) then

      af = a2
      pf = p2
      Rgas = Rgasl

      if (vel <= un1) then
        if (sw_l > 0.d0) then
          af = a1
          vel = un1
          pf = p1
        end if
      else
        if (vw1 > 0.d0) then
          vel = un1
          af = a1
          pf = p1
        else
          if (vw2 > 0.d0) then
            af = (un1 + a1/delta1)/(1.d0 + 1.d0/delta1)
            vel = af
            pf = p1*(af/a1)**eta1
          end if
        end if
      end if

      rhof = gam1*pf/(af**2.d0)
      df = dl/dltot*rhof
      uf = vel*nx + ut1
      vf = vel*ny + vt1

    else

      af = a3
      pf = p3
      Rgas = Rgasr

      if (vel >= un4) then
        if (sw_r < 0.d0) then
          af = a4
          vel = un4
          pf = p4
        end if
      else
        if (vw4 < 0.d0) then
          af = a4
          vel = un4
          pf = p4
        else
          if (vw3 < 0.d0) then
            af = -(un4 - a4/delta4)/(1.d0 + 1.d0/delta4)
            vel = -af
            pf = p4*(af/a4)**eta4
          end if
        end if
      end if

      rhof = gam4*pf/(af**2.d0)
      df = dr/drtot*rhof
      uf = vel*nx + ut4
      vf = vel*ny + vt4

    end if

    ! H0 with full velocity magnitude: H0(...,sqrt(uf**2+vf**2)) (Glaister 1988)
    call fluxes(pf,df,uf,vf,nx,ny,F_r,F_u,F_v,F_E)

  end subroutine riemann_exact


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
    real(R8) :: v_mag

    v_mag = sqrt(u*u + v*v)
    F_r = sum(r)*(u*nx + v*ny)
    F_u = F_r*u + p*nx
    F_v = F_r*v + p*ny
    F_E = F_r*H0(p,r,v_mag)

  end subroutine fluxes


  subroutine rareleft(ps,as,us,ud,pd,ad,dpr,v1,v2,gamma)
    implicit none
    real(R8) :: ps, as, us, ud, pd, ad
    real(R8) :: dum1, dpr
    real(R8) :: gamma, delta
    real(R8) :: v1, v2

    delta = 0.5d0*(gamma - 1)
    dum1 = 2.d0*gamma/(gamma - 1.d0)

    ad = as - delta*(ud - us)
    pd = ps*((ad/as)**dum1)
    dpr = -gamma*pd/ad
    v1 = us - as
    v2 = ud - ad

  end subroutine rareleft


  subroutine rareright(pd,ad,ud,us,ps,as,dpl,v3,v4,gamma)
    implicit none
    real(R8) :: ps, as, us, ud, pd, ad
    real(R8) :: dum1, dpl
    real(R8) :: gamma, delta
    real(R8) :: v3, v4

    delta = 0.5d0*(gamma - 1)
    dum1 = 2.d0*gamma/(gamma - 1.d0)

    as = ad + delta*(us - ud)
    ps = pd*((as/ad)**dum1)
    dpl = gamma*ps/as
    v3 = us + as
    v4 = ud + ad

  end subroutine rareright


  subroutine shockleft(us,as,ps,ud,ad,pd,dpr,ws,gamma)
    implicit none
    real(R8) :: us, as, ps, ud, ad, pd, ws
    real(R8) :: x1, Mlrel, C1, rpd, dpr
    real(R8) :: gamma

    x1 = (gamma + 1.d0)*(ud - us)/(4.d0*as)
    Mlrel = x1 - dsqrt(1.d0 + x1**2.d0)
    C1 = gamma*ps/as
    pd = ps + C1*(ud - us)*Mlrel
    dpr = 2.d0*C1*Mlrel**3.d0/(1.d0 + Mlrel**2.d0)
    rpd = pd/ps
    ws = us + as*Mlrel
    ad = as*dsqrt((gamma + 1.d0 + (gamma - 1.d0)*rpd)/(gamma + 1.d0 + (gamma - 1.d0)/rpd))

  end subroutine shockleft


  subroutine shockright(ud,ad,pd,us,as,ps,dpl,ws,gamma)
    implicit none
    real(R8) :: ud, ad, pd, us, as, ps, ws
    real(R8) :: x4, Mrrel, C4, rps, dpl
    real(R8) :: gamma

    x4 = (gamma + 1.d0)*(us - ud)/(4.d0*ad)
    Mrrel = x4 + dsqrt(1.d0 + x4**2.d0)
    C4 = gamma*pd/ad
    ps = pd + C4*(us - ud)*Mrrel
    dpl = 2.d0*C4*Mrrel**3.d0/(1.d0 + Mrrel**2.d0)
    rps = ps/pd
    ws = ud + ad*Mrrel
    as = ad*dsqrt((gamma + 1.d0 + (gamma - 1.d0)*rps)/(gamma + 1.d0 + (gamma - 1.d0)/rps))

  end subroutine shockright

end module Q2D_Lib_Riemann_Godunov
