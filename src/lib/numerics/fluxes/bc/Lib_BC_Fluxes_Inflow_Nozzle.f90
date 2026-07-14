module Q2D_Lib_BC_Fluxes_Inflow_Nozzle
  ! Router + 4 inner BC
  ! (Supersonic, Subsonic, Subsonic_QS, Borda_Flow). Borda sudden-expansion
  ! injection (Fievisohn 2016 Sec 3.3 / Fievisohn & Yu 2017 Sec II.B Eqs 9-14)
  ! is a Q2D-only feature — sentinel `BC_Mach == -20` with `BC_pAmb`
  ! reinterpreted as A1/A3.
  !
  ! BC_Nozzle_Subsonic keeps the generalised Riemann invariant Glaister 1988
  ! with the fused quadrature Trapezoidal_Riem_PU (a real-gas closure; the MOSE
  ! peer uses the classic closure `u = (Riem-a)/(BC_Sign·(γ-1)·0.5)` valid only
  ! for constant gamma).
  !
  ! Dispatcher router:
  !   Inflow_Outflow >= 0 OR p_b >= p0          -> error 1
  !   BC_Mach == -20                            -> Borda_Flow (g = A1/A3)
  !                                              error==2 (regime c3) -> Supersonic
  !   p_b in [pSub, p0)                         -> Subsonic(BC_Mach=0)
  !   p_b in [pSup, pSub)                       -> Subsonic(BC_Mach=-10)
  !   p_b < pSup AND BC_Mach == 0               -> Supersonic, mdot/A
  !   p_b < pSup AND BC_Mach != 0               -> Supersonic, M_face = BC_Mach
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_Limiters,           only: rlimiter
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index, Trapezoidal, Trapezoidal_Riem_PU
  use Q2D_Mod_Riemann,            only: Riemann
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, f_ss, f_gamma, &
                                        f_tabT, h_tab, cp_tab, dcpi_tab, Ri_tab

  implicit none
  private

  public :: BC_Inflow_Nozzle
  public :: BC_Nozzle_Supersonic
  public :: BC_Nozzle_Subsonic
  public :: BC_Nozzle_Subsonic_QS
  public :: BC_Nozzle_Borda_Flow

contains


  subroutine BC_Inflow_Nozzle ( Im, Jm, Fm, Blk, BC_Mach, BC_T0, BC_p0, BC_rel_fac, &
                                BC_alpha, BC_beta, BC_pAmb, BC_ci, BC_RANS, error )
    implicit none
    integer,             intent(in)    :: Im, Jm, Fm
    real(R8),            intent(in)    :: BC_Mach, BC_T0, BC_p0, BC_rel_fac, BC_alpha, BC_beta
    real(R8),            intent(in)    :: BC_pAmb, BC_ci(nsc), BC_RANS(1:)
    type(Q2D_block_type),intent(inout) :: Blk
    integer,             intent(inout) :: error
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Dir, Face_i, Face_j
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim), Normal(2), int_l, t_Vec(2), t_Mod, BC_Sign, Un
    real(R8) :: Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, Inflow_Outflow, Flux(nprim)
    real(R8) :: pSub, pSup, g, alpha_app, beta_app, BC_supersonic, Bord_mdot


    error         = 0
    BC_supersonic = 0.d0

    pSub      = BC_alpha * 1.d5
    pSup      = BC_beta  * 1.d5
    g         = BC_pAmb
    alpha_app = 0.d0
    beta_app  = 0.d0

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)

    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    Bound_Prim = Blk % P(:,Im,Jm)
    Int_Prim   = Blk % P(:,Int_i,Int_j)

    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l

    select case ( Fm )
      case ( 1 : 2 )
        t_Vec = Blk % node(Face_i,Face_j) % c - Blk % node(Face_i,Face_j-1) % c
      case ( 3 : 4 )
        t_Vec = Blk % node(Face_i,Face_j) % c - Blk % node(Face_i-1,Face_j) % c
    end select

    t_Mod = norm2 ( t_Vec )
    t_Vec = t_Vec / t_Mod

    BC_Sign = real ( modfm2, R8 )

    Un = dot_product ( Bound_Prim(nu:nv), Normal )
    call co_rotot_Rtot ( Bound_Prim(1:nsc), Bound_rho, Bound_Rgas )
    Bound_Sound = f_ss   ( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas )
    Bound_Gamma = f_gamma( Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas )
    Riem        = Bound_Sound + BC_Sign * 0.5d0 * ( Bound_Gamma - 1d0 ) * Un

    Inflow_Outflow = BC_Sign * Un

    if ( Inflow_Outflow >= 0d0 .or. Bound_Prim(np) >= BC_p0 ) then
      error = 1
      return
    elseif ( BC_Mach == -20d0 ) then
      ! Borda sudden-expansion (Fievisohn 2016 Sec 3.3); g = A1/A3.
      call BC_Nozzle_Borda_Flow ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, alpha_app, &
                                  BC_ci, BC_RANS, BC_T0, BC_p0, g, BC_rel_fac,                &
                                  BC_Sign, t_Vec, Normal, int_l, Flux, error, Bord_mdot )
      if ( error == 1 ) return
      if ( error == 2 ) then
        BC_supersonic = 1.d0
        call BC_Nozzle_Supersonic ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2,        &
                                    Bord_mdot, BC_p0, BC_T0, BC_ci, BC_RANS,             &
                                    alpha_app, Normal, int_l, Flux, BC_supersonic )
        error = 0
      endif
    elseif (( Bound_Prim(np) < BC_p0 ) .and. ( Bound_Prim(np) >= pSub )) then
      call BC_Nozzle_Subsonic ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, alpha_app, &
                                0d0, BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign,    &
                                t_Vec, Normal, int_l, Flux )
    elseif (( Bound_Prim(np) < pSub ) .and. ( Bound_Prim(np) >= pSup )) then
      call BC_Nozzle_Subsonic ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, alpha_app, &
                                -10.d0, BC_ci, BC_RANS, BC_T0, g, BC_rel_fac, BC_Sign,     &
                                t_Vec, Normal, int_l, Flux )
    elseif ( Bound_Prim(np) < pSup ) then
      if (BC_Mach == 0.d0) then
        BC_supersonic = 1.d0
        call BC_Nozzle_Supersonic ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, g, BC_p0, &
                                    BC_T0, BC_ci, BC_RANS, alpha_app, Normal, int_l, Flux, BC_supersonic )
      else
        call BC_Nozzle_Supersonic ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p0, &
                                    BC_T0, BC_ci, BC_RANS, alpha_app, Normal, int_l, Flux, BC_supersonic )
      endif
    endif

    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + Modfm2 * Flux

  end subroutine BC_Inflow_Nozzle


  subroutine BC_Nozzle_Supersonic ( Prim2, Int_Prim, modfm, modfm1, modfm2, g, BC_p0, BC_T0, BC_ci, &
                                    BC_RANS, alpha_app, Normal, int_l, Flux, BC_supersonic )
    implicit none
    integer,  intent(in)  :: modfm, modfm1, modfm2
    real(R8), intent(in)  :: Prim2(nprim), Int_Prim(nprim), g, BC_p0, BC_T0, BC_supersonic
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(1:), alpha_app, Normal(2), int_l
    real(R8), intent(out) :: Flux(nprim)
    ! Local
    integer  :: s
    real(R8) :: Sup_Prim(nprim), Sup_T, Sup_Rgas, Sup_rho, Sup_Sound, Prim1(nprim), Prim3(nprim)
    real(R8), dimension(nprim) :: Diff32, Diff21, Slope, Face_Prim, Prim_L, Prim_R
    real(R8) :: rho_L, rho_R, Rgas_L, Rgas_R, Sound_L, Sound_R, F_r, F_u, F_v, F_E
    real(R8) :: su, Sel_L, Sel_R, Fmass
    real(R8) :: h0_, cpSup, dcpSup, hSup, uSup, Sup_h0_, duSup, Sup_p, dpSup, dh0_Sup
    real(R8) :: Sup_gamma, gamSup, dgamSup, daSup
    real(R8) :: fun_T, diff_h0_, alpha


    if ( alpha_app == 0d0 ) then
      alpha = atan ( Normal(2) / ( Normal(1) + 1d-20 ) )
    else
      alpha = alpha_app
    endif

    h0_      = 0d0
    Sup_Rgas = 0d0
    do s = 1, nsc
      h0_      = h0_      + BC_ci(s) * f_tabT( BC_T0, s, h_tab )
      Sup_Rgas = Sup_Rgas + BC_ci(s) * Ri_tab(s)
    enddo
    Sup_T    = 0.5d0 * BC_T0
    diff_h0_ = 1d0

    do while (diff_h0_ >= 1.d-9)
      hSup   = 0d0
      cpSup  = 0d0
      dcpSup = 0d0
      do s = 1, nsc
        hSup   = hSup   + BC_ci(s) * f_tabT( Sup_T, s, h_tab )
        cpSup  = cpSup  + BC_ci(s) * f_tabT( Sup_T, s, cp_tab )
        dcpSup = dcpSup + BC_ci(s) * f_tabT( Sup_T, s, dcpi_tab )
      enddo
      if (BC_supersonic == 1.d0) then
        call Trapezoidal ( sup_T, BC_T0, BC_ci, Sup_Rgas, fun_T )
        Sup_p   = BC_p0 / exp ( fun_T )
        Sup_rho = Sup_p / ( Sup_Rgas * Sup_T )
        uSup    = g / Sup_rho
        Sup_h0_ = hSup + 0.5d0 * uSup**2
        diff_h0_ = abs ( Sup_h0_ - h0_ ) / abs ( h0_ )
        dpSup   = BC_p0 / exp(fun_T) * cpSup / ( Sup_Rgas * Sup_T )
        duSup   = g * Sup_Rgas * ( Sup_p - dpSup * Sup_T ) / Sup_p**2
        dh0_Sup = cpSup + uSup * duSup
      else
        gamSup    = cpSup / ( cpSup - Sup_Rgas )
        Sup_Sound = sqrt( gamSup * Sup_Rgas * Sup_T )
        uSup      = Sup_Sound * g
        Sup_h0_   = hSup + 0.5d0 * uSup**2
        diff_h0_  = abs ( Sup_h0_ - h0_ ) / abs ( h0_ )
        dgamSup   = - Sup_Rgas * dcpSup / ( cpSup - Sup_Rgas )**2
        daSup     = sqrt(Sup_Rgas) * 0.5d0 / sqrt( gamSup * Sup_T ) * ( gamSup + Sup_T * dgamSup )
        duSup     = g * daSup
        dh0_Sup   = cpSup + uSup * duSup
      endif
      Sup_T = Sup_T - ( Sup_h0_ - h0_ ) / dh0_Sup
    enddo

    call Trapezoidal ( sup_T, BC_T0, BC_ci, Sup_Rgas, fun_T )

    Sup_Prim(np) = BC_p0 / exp ( fun_T )
    Sup_rho      = Sup_Prim(np) / ( Sup_Rgas * Sup_T )
    Sup_Prim(1:nsc) = Sup_rho * BC_ci
    Sup_Sound    = f_ss ( Sup_Prim(1:nsc), Sup_Prim(np), Sup_rho, Sup_Rgas )

    if (BC_supersonic == 1.d0) then
      Sup_Prim(nu:nv) = g / Sup_rho * normal
    else
      Sup_Prim(nu:nv) = g * Sup_Sound * normal
    endif
    if (model==2) Sup_Prim(nt:nprim) = Sup_rho * BC_RANS

    Prim1 = modfm * Sup_Prim + modfm1 * Int_Prim
    Prim3 = modfm * Int_Prim + modfm1 * Sup_Prim

    Diff32 = Prim3 - Prim2
    Diff21 = Prim2 - Prim1
    do s = 1, nprim
      Slope(s) = rlimiter ( Diff32(s), Diff21(s) )
    end do
    Face_Prim = Prim2 + 0.5d0 * Slope * modfm2

    Prim_L = modfm * Sup_Prim + modfm1 * Face_Prim
    call co_rotot_Rtot ( Prim_L(1:nsc), rho_L, Rgas_L )
    Sound_L = f_ss ( Prim_L(1:nsc), Prim_L(np), rho_L, Rgas_L )

    Prim_R = modfm * Face_Prim + modfm1 * Sup_Prim
    call co_rotot_Rtot ( Prim_R(1:nsc), rho_R, Rgas_R )
    Sound_R = f_ss ( Prim_R(1:nsc), Prim_R(np), rho_R, Rgas_R )

    ! beta = 1d0 (SD-coupled blending, full HLLC weight). Q2D legacy 0d0 was the
    ! pre-SD default (Borda preceded the SD-coupled solver).
    call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(np), Sound_L, rho_L, &
                   Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(np), Sound_R, rho_R, &
                   1d0, Normal(1), Normal(2), F_r, F_u, F_v, F_E )

    su = sign ( 1d0, F_r )
    Sel_L = ( 1d0 + su ) / 2d0
    Sel_R = ( su  - 1d0 ) / 2d0

    Flux = 0d0   ! mirror MOSE pre-zero: covers the soot/passive band np+1:nt-1
    Fmass = F_r * int_l
    Flux(1:nsc) = Fmass * ( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
    Flux(nu) = F_u * int_l
    Flux(nv) = F_v * int_l
    Flux(np) = F_E * int_l

    if (model==2) Flux(nt:nprim) = Fmass * ( Sel_L*Prim_L(nt:nprim)/rho_L - Sel_R*Prim_R(nt:nprim)/rho_R )

  end subroutine BC_Nozzle_Supersonic


  subroutine BC_Nozzle_Subsonic ( Prim, Un, Sound, Gamma, Riem, alpha_app, BC_Mach, &
                                  BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, &
                                  t_Vec, Normal, int_l, Flux )
    ! Generalised real-gas Riemann invariant Glaister 1988. Fused quadrature
    ! Trapezoidal_Riem_PU.
    implicit none
    real(R8), intent(in)  :: Prim(nprim), Un, Sound, Gamma, Riem, alpha_app, BC_Mach
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(1:), BC_T0, BC_p0, BC_rel_fac, BC_Sign
    real(R8), intent(in)  :: t_Vec(2), Normal(2), int_l
    real(R8), intent(out) :: Flux(nprim)
    ! Local
    integer  :: s
    real(R8) :: b_Vec(3), b_Mod, alpha, XA, XB, XC, XD, XE, XF, check, E1, E2
    real(R8) :: Face_Rgas, Un3, p3, Ut3, Ub3, mdot, Rgas3, cp3, T3, rho3, h3, h0_3
    real(R8) :: dcp_acu, h0_3_iter, p03, h3_iter, u3
    real(R8) :: Face_T, Face_Prim(nprim), Face_rho, Face_Enthalpy, Fmass, Face_Un
    real(R8) :: diff_h0_, du, dT, dp, dh0_, dp0
    real(R8) :: Rgas_acu, T_acu, h_acu, p_acu, u_acu, cp_acu, gamma_acu, a_acu
    real(R8) :: T_bound, Bound_ci(nsc), fun_T, diff_T, riem_K
    real(R8) :: rho_Bound


    b_Vec(1) = 0.d0
    b_Vec(2) = 0.d0
    b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
    b_Mod    = norm2 ( b_Vec )
    b_Vec    = b_Vec / b_Mod

    if ( alpha_app == 0d0 ) then
      alpha = atan ( Normal(2) / ( Normal(1) + 1d-20 ) )
    else
      alpha = alpha_app
    endif

    XA = Normal(2) - Normal(1) * tan(alpha)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha)
    XC =  b_Vec(2) -  b_Vec(1) * tan(alpha)
    XD = 0.d0
    XE = 0.d0
    XF = b_Vec(3)
    check = ( XC * XE - XF * XB )
    E1 = ( XF * XA - XC * XD ) / check
    E2 = ( XD * XB - XE * XA ) / check

    mdot  = BC_p0
    Rgas3 = sum ( BC_ci * Ri_tab )
    call co_rotot_Rtot ( Prim(1:nsc), rho_Bound, Rgas_acu )
    diff_h0_ = 1.d0

    if ( BC_Mach == -10d0 .or. BC_Mach >= 0.d0 ) then
      h0_3 = 0d0
      do s = 1, nsc
        h0_3 = h0_3 + BC_ci(s) * f_tabT( BC_T0, s, h_tab )
      enddo
    endif
    T_Bound  = Prim(np) / rho_Bound / Rgas_acu
    T_acu    = 1.1d0 * T_Bound
    Bound_ci = Prim(1:nsc) / rho_Bound

    do while (diff_h0_ >= 1.d-6)
      h_acu   = 0.d0
      cp_acu  = 0.d0
      dcp_acu = 0.d0
      do s = 1, nsc
        h_acu   = h_acu   + Bound_ci(s) * f_tabT( T_acu, s, h_tab )
        cp_acu  = cp_acu  + Bound_ci(s) * f_tabT( T_acu, s, cp_tab )
        dcp_acu = dcp_acu + Bound_ci(s) * f_tabT( T_acu, s, dcpi_tab )
      enddo
      gamma_acu = cp_acu / ( cp_acu - Rgas_acu )
      a_acu     = sqrt( gamma_acu * Rgas_acu * T_acu )
      ! Glaister 1988 fused quadrature: J = int cp/(R*T) dT (isentropic),
      ! K = int a/((gamma-1)*T) dT (generalised C+/-).
      call Trapezoidal_Riem_PU ( T_Bound, T_acu, Bound_ci, Rgas_acu, fun_T, riem_K )
      p_acu = Prim(np) * exp ( fun_T )
      u_acu = Un - BC_Sign * riem_K
      p3    = p_acu
      u3    = u_acu
      dp    = p_acu * cp_acu / ( Rgas_acu * T_acu )
      du    = - BC_Sign * a_acu / ( ( gamma_acu - 1.d0 ) * T_acu )
      if ( BC_Mach >= 0d0 ) then
        diff_T = 1.d0
        h3     = h0_3 - 0.5d0 * u3**2
        T3     = 0.9 * BC_T0
        do while (diff_T >= 1.d-6)
          cp3     = 0.d0
          h3_iter = 0.d0
          do s = 1, nsc
            h3_iter = h3_iter + BC_ci(s) * f_tabT( T3, s, h_tab )
            cp3     = cp3     + BC_ci(s) * f_tabT( T3, s, cp_tab )
          enddo
          diff_T = abs ( h3 - h3_iter ) / abs ( h3 )
          T3     = T3 - ( h3_iter - h3 ) / cp3
        enddo
        cp3 = 0.d0
        do s = 1, nsc
          cp3 = cp3 + BC_ci(s) * f_tabT( T3, s, cp_tab )
        enddo
        call Trapezoidal ( T3, BC_T0, BC_Ci, Rgas3, fun_T )
        p03      = p3 * exp ( fun_T )
        diff_h0_ = abs ( p03 - BC_p0 ) / BC_p0
        dp0      = dp * exp(fun_T) + p3 * exp(fun_T) / ( Rgas3 * T3 ) * u3 * du
        T_acu    = T_acu - ( p03 - BC_p0 ) / dp0
      else
        dT     = ( 1 / Rgas3 / mdot ) * ( dp * u3 + du * p3 )
        rho3   = mdot / u3
        T3     = p3 / ( Rgas3 * rho3 )
        h3     = 0.d0
        cp3    = 0.d0
        do s = 1, nsc
          h3  = h3  + BC_ci(s) * f_tabT( T3, s, h_tab )
          cp3 = cp3 + BC_ci(s) * f_tabT( T3, s, cp_tab )
        enddo
        if ( BC_Mach == -10d0 ) then
          h0_3_iter = h3 + 0.5d0 * u3**2
          diff_h0_  = abs ( h0_3 - h0_3_iter ) / abs ( h0_3 )
          dh0_      = cp3 * dT + u3 * du
          T_acu     = T_acu - ( h0_3_iter - h0_3 ) / dh0_
        else
          diff_h0_  = abs ( T3 - BC_T0 ) / BC_T0
          T_acu     = T_acu - ( T3 - BC_T0 ) / dT
        endif
      endif
    enddo

    h_acu  = 0.d0
    cp_acu = 0.d0
    do s = 1, nsc
      h_acu  = h_acu  + Bound_ci(s) * f_tabT( T_acu, s, h_tab )
      cp_acu = cp_acu + Bound_ci(s) * f_tabT( T_acu, s, cp_tab )
    enddo
    gamma_acu = cp_acu / ( cp_acu - Rgas_acu )
    a_acu     = sqrt( gamma_acu * Rgas_acu * T_acu )
    call Trapezoidal_Riem_PU ( T_Bound, T_acu, Bound_ci, Rgas_acu, fun_T, riem_K )
    p_acu = Prim(np) * exp ( fun_T )
    u_acu = Un - BC_Sign * riem_K
    p3    = p_acu
    u3    = u_acu

    if ( BC_Mach >= 0.d0 ) then
      h3 = h0_3 - 0.5d0 * u3**2
      T3 = 0.9 * BC_T0
      diff_T = 1.d0
      do while (diff_T >= 1.d-6)
        cp3     = 0.d0
        h3_iter = 0.d0
        do s = 1, nsc
          h3_iter = h3_iter + BC_ci(s) * f_tabT( T3, s, h_tab )
          cp3     = cp3     + BC_ci(s) * f_tabT( T3, s, cp_tab )
        enddo
        diff_T = abs ( h3 - h3_iter ) / abs ( h3 )
        T3     = T3 - ( h3_iter - h3 ) / cp3
      enddo
      rho3 = p3 / ( Rgas3 * T3 )
    else
      rho3 = mdot / u3
    endif

    Un3 = u3
    Ut3 = Un3 * E1
    Ub3 = Un3 * E2
    Face_Rgas = Rgas3
    Face_rho  = rho3
    Face_Prim(1:nsc) = Face_rho * BC_ci
    Face_Prim(nu:nv) = Un3 * Normal + Ut3 * t_Vec
    Face_Prim(np)    = p3
    Face_Un = Un3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    Flux = 0d0   ! mirror MOSE pre-zero: covers the soot/passive band np+1:nt-1
    Face_Enthalpy = 0d0
    Face_T = Face_Prim(np) / ( Face_Rgas * Face_rho )
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT( Face_T, s, h_tab )
      Flux(s) = Face_Prim(s) * Face_Un * int_l
    end do
    Fmass = Face_rho * Face_Un * int_l
    Flux(nu:nv) = Fmass * Face_Prim(nu:nv) + Face_Prim(np) * int_l * Normal
    Flux(np)    = Fmass * ( 0.5d0 * sum( Face_Prim(nu:nv)**2 ) + Face_Enthalpy )
    if (model==2) Flux(nt:nprim) = Face_Prim(nt:nprim) / Face_rho * Fmass

  end subroutine BC_Nozzle_Subsonic


  subroutine BC_Nozzle_Subsonic_QS ( Prim, Un, Sound, Gamma, Riem, alpha_app, BC_Mach, &
                                     BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac, BC_Sign, &
                                     t_Vec, Normal, int_l, Flux )
    ! Quasi-steady variant: no upstream-running characteristic, p3 = Prim(np)
    ! extrapolated. The BC_Mach==-15 branch (HRE tangential injector) was removed
    ! in Q2D legacy; the 3-branch dispatch is preserved.
    implicit none
    real(R8), intent(in)  :: Prim(nprim), Un, Sound, Gamma, Riem, alpha_app, BC_Mach
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(1:), BC_T0, BC_p0, BC_rel_fac, BC_Sign
    real(R8), intent(in)  :: t_Vec(2), Normal(2), int_l
    real(R8), intent(out) :: Flux(nprim)
    ! Local
    integer  :: s
    real(R8) :: alpha
    real(R8) :: Face_Rgas, p3, mdot, Rgas3, cp3, T3, rho3, h3, h0_3, u3
    real(R8) :: Face_T, Face_Prim(nprim), Face_rho, Face_Enthalpy, Fmass, Face_Un
    real(R8) :: diff_h0_, fun_T


    if ( alpha_app == 0d0 ) then
      alpha = atan ( Normal(2) / ( Normal(1) + 1d-20 ) )
    else
      alpha = alpha_app
    endif

    mdot  = BC_p0
    Rgas3 = sum ( BC_ci * Ri_tab )
    diff_h0_ = 1.d0

    if ( BC_Mach == -10d0 .or. BC_Mach >= 0.d0 ) then
      h0_3 = 0d0
      do s = 1, nsc
        h0_3 = h0_3 + BC_ci(s) * f_tabT( BC_T0, s, h_tab )
      enddo
    endif

    p3 = Prim(np)
    Face_Prim(np) = Prim(np)
    T3 = 0.9d0 * BC_T0

    do while (diff_h0_ >= 1.d-6)
      if ( BC_Mach >= 0d0 ) then
        cp3 = 0.d0
        do s = 1, nsc
          cp3 = cp3 + BC_ci(s) * f_tabT( T3, s, cp_tab )
        enddo
        call Trapezoidal ( T3, BC_T0, BC_ci, Rgas3, fun_T )
        diff_h0_ = abs ( fun_T - log(BC_p0/p3) ) / abs ( log(BC_p0/p3) )
        T3       = T3 + ( fun_T - log(BC_p0/p3) ) / ( cp3 / Rgas3 / T3 )
        h3 = 0.d0
        do s = 1, nsc
          h3 = h3 + BC_ci(s) * f_tabT( T3, s, h_tab )
        enddo
        u3 = sqrt( 2 * ( h0_3 - h3 ) )
      elseif ( BC_Mach == -10.d0 ) then
        h3  = 0.d0
        cp3 = 0.d0
        do s = 1, nsc
          h3  = h3  + BC_ci(s) * f_tabT( T3, s, h_tab )
          cp3 = cp3 + BC_ci(s) * f_tabT( T3, s, cp_tab )
        enddo
        u3       = sqrt( 2 * ( h0_3 - h3 ) )
        rho3     = p3 / ( Rgas3 * T3 )
        diff_h0_ = abs ( mdot - rho3 * u3 ) / abs ( mdot )
        T3 = T3 + (rho3 * u3 - mdot) / ( sqrt(2d0) * p3 / Rgas3 *                   &
                                          ( 0.5d0 * cp3 * T3 / sqrt(h0_3 - h3) +    &
                                            sqrt(h0_3 - h3) ) / T3**2 )
        u3 = mdot / ( p3 / Rgas3 / T3 )
      else
        T3       = BC_T0
        rho3     = p3 / ( Rgas3 * T3 )
        u3       = mdot / rho3
        diff_h0_ = 1.d-7
      endif
    enddo

    rho3      = p3 / ( Rgas3 * T3 )
    Face_Rgas = Rgas3
    Face_rho  = rho3
    Face_Prim(1:nsc) = Face_rho * BC_ci
    Face_Prim(nu:nv) = u3 * Normal
    Face_Un          = u3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    Flux = 0d0   ! mirror MOSE pre-zero: covers the soot/passive band np+1:nt-1
    Face_Enthalpy = 0d0
    Face_T = Face_Prim(np) / ( Face_Rgas * Face_rho )
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT( Face_T, s, h_tab )
      Flux(s) = Face_Prim(s) * Face_Un * int_l
    end do
    Fmass       = Face_rho * Face_Un * int_l
    Flux(nu:nv) = Fmass * Face_Prim(nu:nv) + Face_Prim(np) * int_l * Normal
    Flux(np)    = Fmass * ( 0.5d0 * sum( Face_Prim(nu:nv)**2 ) + Face_Enthalpy )
    if (model==2) Flux(nt:nprim) = Face_Prim(nt:nprim) / Face_rho * Fmass

  end subroutine BC_Nozzle_Subsonic_QS


  subroutine BC_Nozzle_Borda_Flow ( Prim, Un, Sound, Gamma, Riem, alpha_app, &
                                    BC_ci, BC_RANS, BC_T0, BC_p0, A1_A3, BC_rel_fac, &
                                    BC_Sign, t_Vec, Normal, int_l, Flux, error, mdot_sup )
    ! Borda sudden-expansion injection (Fievisohn 2016 Sec 3.3; Fievisohn & Yu
    ! 2017 Sec II.B Eqs 9-14). Plenum (BC_T0,BC_p0) -> throat (M1, A1) -> sudden
    ! expansion -> channel face (A3 > A1). r = A1/A3 in (0,1].
    !
    ! Three regimes: (a) choked + sub3 (Borda valid), (b) unchoked (bisection
    ! M1 in (M1_MIN,1) for P1=P2), (c) choked + sup3 (Borda invalid; fallback
    ! isentropic CD-nozzle, sub3 -> Subsonic, sub-sonic post-shock -> Subsonic
    ! mdot_choked, sup3 -> error=2 + caller's Supersonic).
    !
    ! Error codes: 0 = success, 1 = inner Newton not converged, 2 = c3 fallback.
    implicit none
    real(R8), intent(in)  :: Prim(nprim), Un, Sound, Gamma, Riem, alpha_app
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(1:), BC_T0, BC_p0, A1_A3, BC_rel_fac
    real(R8), intent(in)  :: BC_Sign, t_Vec(2), Normal(2), int_l
    real(R8), intent(out) :: Flux(nprim)
    integer,  intent(out) :: error
    real(R8), intent(out) :: mdot_sup

    integer,  parameter :: MAX_NEWTON = 100, MAX_BISECT = 40
    real(R8), parameter :: TOL        = 1.d-6
    real(R8), parameter :: M1_MIN     = 0.01d0
    real(R8), parameter :: P1P2_TOL   = 1.d-4
    real(R8), parameter :: M3_SUP_TH  = 0.99d0

    real(R8) :: Rgas_inj, Rgas_bound, rho_Bound, Bound_ci(nsc), T_Bound, h_plen
    real(R8) :: gamma_plen
    real(R8) :: M1, T1, P1, rho1, v1, a1, cp1, mdot_A3
    real(R8) :: T3, p3, rho3, u3, cp3
    real(R8) :: T_sub_iso, pSub_iso, T_sup_iso, pSup_iso

    real(R8) :: b_Vec(3), b_Mod, alpha, E1, E2
    real(R8) :: Un3, Ut3, Face_Rgas, Face_rho, Face_Prim(nprim), Face_T
    real(R8) :: Face_Un, Face_Enthalpy, Fmass
    integer  :: s
    logical  :: ok, ok_sub, ok_sup


    error    = 0
    mdot_sup = 0.d0

    if ( A1_A3 <= 0.d0 .or. A1_A3 > 1.d0 + 1.d-6 ) then
      write(*,*) 'BC_Nozzle_Borda_Flow: invalid A1_A3 =', A1_A3
      error = 1
      return
    endif

    call set_tangent_coefs ()

    Rgas_inj   = sum ( BC_ci * Ri_tab )
    h_plen     = 0.d0
    gamma_plen = 0.d0
    do s = 1, nsc
      h_plen     = h_plen     + BC_ci(s) * f_tabT ( BC_T0, s, h_tab  )
      gamma_plen = gamma_plen + BC_ci(s) * f_tabT ( BC_T0, s, cp_tab )
    enddo
    gamma_plen = gamma_plen / ( gamma_plen - Rgas_inj )

    call co_rotot_Rtot ( Prim(1:nsc), rho_Bound, Rgas_bound )
    T_Bound  = Prim(np) / ( rho_Bound * Rgas_bound )
    Bound_ci = Prim(1:nsc) / rho_Bound

    if ( A1_A3 >= 1.d0 - 1.d-6 ) then
      call throat_state ( 1.d0, ok )
      if ( .not. ok ) then
        error = 1
        return
      endif

      if ( Prim(np) < P1 ) then
        T3   = T1
        p3   = P1
        rho3 = rho1
        u3   = - BC_Sign * v1
        cp3  = cp1
      else
        call BC_Nozzle_Subsonic ( Prim, Un, Sound, Gamma, Riem, alpha_app, 0.d0, &
                                  BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac,       &
                                  BC_Sign, t_Vec, Normal, int_l, Flux )
        return
      endif
    else
      call solve_at_M1 ( 1.d0, ok )
      if ( .not. ok ) then
        error = 1
        return
      endif

      if ( momentum_P2() > P1 + P1P2_TOL * P1 ) then
        call bisect_M1_unchoked ( ok )
        if ( .not. ok ) then
          error = 1
          return
        endif
      endif

      if ( station3_is_supersonic() ) then
        if ( abs( M1 - 1.d0 ) > 1.d-6 ) then
          call throat_state ( 1.d0, ok )
          if ( .not. ok ) then
            error = 1
            return
          endif
        endif

        call isentropic_face_root ( 0.5d0*(T1 + BC_T0),                                 &
                                    T1    + 1.d-6 *(BC_T0 - T1),                        &
                                    BC_T0 - 1.d-6 *(BC_T0 - T1),                        &
                                    T_sub_iso, pSub_iso, ok_sub )
        call isentropic_face_root ( 0.5d0*T1,                                           &
                                    0.1d0*BC_T0,                                        &
                                    T1    - 1.d-6 * T1,                                 &
                                    T_sup_iso, pSup_iso, ok_sup )

        if ( .not. ok_sub .or. .not. ok_sup ) then
          write(*,*) 'BC_Nozzle_Borda_Flow: isentropic_face_root failed, ', &
                     'falling back to supersonic dispatch'
          mdot_sup = mdot_A3
          error    = 2
          return
        endif

        if ( Prim(np) >= pSub_iso ) then
          call BC_Nozzle_Subsonic ( Prim, Un, Sound, Gamma, Riem, alpha_app, 0.d0, &
                                    BC_ci, BC_RANS, BC_T0, BC_p0, BC_rel_fac,       &
                                    BC_Sign, t_Vec, Normal, int_l, Flux )
          return
        elseif ( Prim(np) >= pSup_iso ) then
          call BC_Nozzle_Subsonic ( Prim, Un, Sound, Gamma, Riem, alpha_app, -10.d0, &
                                    BC_ci, BC_RANS, BC_T0, mdot_A3, BC_rel_fac,        &
                                    BC_Sign, t_Vec, Normal, int_l, Flux )
          return
        else
          mdot_sup = mdot_A3
          error    = 2
          return
        endif
      endif
    endif

    Un3 = u3
    Ut3 = Un3 * E1
    Face_Rgas = Rgas_inj
    Face_rho  = rho3
    Face_Prim(1:nsc) = Face_rho * BC_ci
    Face_Prim(nu:nv) = Un3 * Normal + Ut3 * t_Vec
    Face_Prim(np)    = p3
    Face_Un = Un3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    Flux = 0.d0   ! mirror MOSE pre-zero: covers the soot/passive band np+1:nt-1
    Face_Enthalpy = 0.d0
    Face_T = Face_Prim(np) / ( Face_Rgas * Face_rho )
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT ( Face_T, s, h_tab )
      Flux(s) = Face_Prim(s) * Face_Un * int_l
    enddo
    Fmass       = Face_rho * Face_Un * int_l
    Flux(nu:nv) = Fmass * Face_Prim(nu:nv) + Face_Prim(np) * int_l * Normal
    Flux(np)    = Fmass * ( 0.5d0 * sum( Face_Prim(nu:nv)**2 ) + Face_Enthalpy )
    if (model==2) Flux(nt:nprim) = Face_Prim(nt:nprim) / Face_rho * Fmass

  contains

    subroutine set_tangent_coefs ()
      real(R8) :: XA, XB, XC, XD, XE, XF, check
      b_Vec(1) = 0.d0
      b_Vec(2) = 0.d0
      b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
      b_Mod    = norm2 ( b_Vec )
      b_Vec    = b_Vec / b_Mod
      if ( alpha_app == 0.d0 ) then
        alpha = atan ( Normal(2) / ( Normal(1) + 1.d-20 ) )
      else
        alpha = alpha_app
      endif
      XA = Normal(2) - Normal(1) * tan(alpha)
      XB =  t_Vec(2) -  t_Vec(1) * tan(alpha)
      XC =  b_Vec(2) -  b_Vec(1) * tan(alpha)
      XD = 0.d0
      XE = 0.d0
      XF = b_Vec(3)
      check = ( XC * XE - XF * XB )
      E1 = ( XF * XA - XC * XD ) / check
      E2 = ( XD * XB - XE * XA ) / check
    end subroutine set_tangent_coefs


    subroutine solve_at_M1 ( M1_in, ok_out )
      real(R8), intent(in)  :: M1_in
      logical,  intent(out) :: ok_out
      call throat_state ( M1_in, ok_out )
      if ( .not. ok_out ) return
      call borda_middle_newton ( ok_out )
    end subroutine solve_at_M1


    subroutine throat_state ( M1_in, ok_out )
      real(R8), intent(in)  :: M1_in
      logical,  intent(out) :: ok_out
      real(R8) :: h1, dcp1, gamma1, dgam1, da1, dh_dT1, h_iter, fun_T, diff
      integer  :: s_l, it
      M1 = M1_in
      T1 = BC_T0 / ( 1.d0 + 0.5d0 * ( gamma_plen - 1.d0 ) * M1_in**2 )
      ok_out = .false.
      do it = 1, MAX_NEWTON
        h1 = 0.d0; cp1 = 0.d0; dcp1 = 0.d0
        do s_l = 1, nsc
          h1   = h1   + BC_ci(s_l) * f_tabT ( T1, s_l, h_tab )
          cp1  = cp1  + BC_ci(s_l) * f_tabT ( T1, s_l, cp_tab )
          dcp1 = dcp1 + BC_ci(s_l) * f_tabT ( T1, s_l, dcpi_tab )
        enddo
        gamma1 = cp1 / ( cp1 - Rgas_inj )
        a1 = sqrt ( gamma1 * Rgas_inj * T1 )
        v1 = M1 * a1
        h_iter = h1 + 0.5d0 * v1**2
        diff = abs ( h_iter - h_plen ) / abs ( h_plen )
        if ( diff < 1.d-9 ) then
          ok_out = .true.
          exit
        endif
        dgam1  = - Rgas_inj * dcp1 / ( cp1 - Rgas_inj )**2
        da1    = 0.5d0 * sqrt ( Rgas_inj / ( gamma1 * T1 ) ) * ( gamma1 + T1 * dgam1 )
        dh_dT1 = cp1 + v1 * ( M1 * da1 )
        T1     = T1 - ( h_iter - h_plen ) / dh_dT1
      enddo
      if ( .not. ok_out ) then
        write(*,*) 'BC_Nozzle_Borda_Flow: throat Newton NOT converged, M1=', M1
        return
      endif
      call Trapezoidal ( T1, BC_T0, BC_ci, Rgas_inj, fun_T )
      P1      = BC_p0 / exp ( fun_T )
      rho1    = P1 / ( Rgas_inj * T1 )
      mdot_A3 = rho1 * v1 * A1_A3
    end subroutine throat_state


    subroutine borda_middle_newton ( ok_out )
      logical, intent(out) :: ok_out
      real(R8) :: T_acu, cp_acu, gamma_acu, a_acu
      real(R8) :: dp, du, dT3, drho3, dmdot, res, diff
      integer  :: it
      T_acu = 1.1d0 * T_Bound
      p3 = Prim(np); u3 = 0.d0; T3 = 0.9d0 * BC_T0; rho3 = rho_Bound; cp3 = 0.d0
      ok_out = .false.
      do it = 1, MAX_NEWTON
        call riemann_march ( T_acu, p3, u3, cp_acu, gamma_acu, a_acu, dp, du )
        call invert_h_for_T ( h_plen - 0.5d0 * u3**2, T3, cp3 )
        rho3 = p3 / ( Rgas_inj * T3 )
        res  = rho3 * abs( u3 ) - mdot_A3
        diff = abs( res ) / max( abs( mdot_A3 ), 1.d-30 )
        if ( diff < TOL ) then
          ok_out = .true.
          return
        endif
        dT3   = - u3 * du / cp3
        drho3 = dp / ( Rgas_inj * T3 ) - p3 * dT3 / ( Rgas_inj * T3**2 )
        dmdot = drho3 * abs( u3 ) + rho3 * sign( 1.d0, u3 ) * du
        T_acu = T_acu - res / dmdot
      enddo
      write(*,*) 'BC_Nozzle_Borda_Flow: middle Newton NOT converged'
    end subroutine borda_middle_newton


    subroutine riemann_march ( T_acu, p_out, u_out, cp_out, gam_out, a_out, dp_dT, du_dT )
      real(R8), intent(in)  :: T_acu
      real(R8), intent(out) :: p_out, u_out, cp_out, gam_out, a_out, dp_dT, du_dT
      real(R8) :: fun_T, riem_I
      integer  :: s_l
      cp_out = 0.d0
      do s_l = 1, nsc
        cp_out = cp_out + Bound_ci(s_l) * f_tabT ( T_acu, s_l, cp_tab )
      enddo
      gam_out = cp_out / ( cp_out - Rgas_bound )
      a_out   = sqrt ( gam_out * Rgas_bound * T_acu )
      call Trapezoidal_Riem_PU ( T_Bound, T_acu, Bound_ci, Rgas_bound, fun_T, riem_I )
      p_out = Prim(np) * exp ( fun_T )
      u_out = Un - BC_Sign * riem_I
      dp_dT = p_out * cp_out / ( Rgas_bound * T_acu )
      du_dT = - BC_Sign * a_out / ( ( gam_out - 1.d0 ) * T_acu )
    end subroutine riemann_march


    subroutine invert_h_for_T ( h_target, T_out, cp_out )
      real(R8), intent(in)    :: h_target
      real(R8), intent(inout) :: T_out
      real(R8), intent(out)   :: cp_out
      real(R8) :: h_iter, diff
      integer  :: s_l, it
      do it = 1, MAX_NEWTON
        h_iter = 0.d0; cp_out = 0.d0
        do s_l = 1, nsc
          h_iter = h_iter + BC_ci(s_l) * f_tabT ( T_out, s_l, h_tab )
          cp_out = cp_out + BC_ci(s_l) * f_tabT ( T_out, s_l, cp_tab )
        enddo
        diff = abs ( h_iter - h_target ) / abs ( h_target )
        if ( diff < TOL ) return
        T_out = T_out - ( h_iter - h_target ) / cp_out
      enddo
    end subroutine invert_h_for_T


    subroutine isentropic_face_root ( T_init, T_lo, T_hi, T_root, p_root, ok_out )
      real(R8), intent(in)  :: T_init, T_lo, T_hi
      real(R8), intent(out) :: T_root, p_root
      logical,  intent(out) :: ok_out
      real(R8) :: T, h_T, cp_T, V, p, rho, fun_T, f, drho, dV, df, diff
      real(R8) :: T_step
      integer  :: it, s_l

      T = max( T_lo, min( T_init, T_hi ) )
      ok_out = .false.

      do it = 1, MAX_NEWTON
        h_T  = 0.d0
        cp_T = 0.d0
        do s_l = 1, nsc
          h_T  = h_T  + BC_ci(s_l) * f_tabT ( T, s_l, h_tab  )
          cp_T = cp_T + BC_ci(s_l) * f_tabT ( T, s_l, cp_tab )
        enddo
        V = sqrt ( 2.d0 * ( h_plen - h_T ) )
        call Trapezoidal ( T, BC_T0, BC_ci, Rgas_inj, fun_T )
        p   = BC_p0 / exp ( fun_T )
        rho = p / ( Rgas_inj * T )
        f    = rho * V - mdot_A3
        diff = abs( f ) / max( mdot_A3, 1.d-30 )
        if ( diff < TOL ) then
          ok_out = .true.
          T_root = T
          p_root = p
          return
        endif
        drho = rho * ( cp_T - Rgas_inj ) / ( Rgas_inj * T )
        dV   = - cp_T / V
        df   = drho * V + rho * dV
        T_step = - f / df
        if ( T + T_step <= T_lo .or. T + T_step >= T_hi ) then
          T = 0.5d0 * ( T + 0.5d0 * ( T_lo + T_hi ) )
        else
          T = T + T_step
        endif
      enddo
      write(*,*) 'BC_Nozzle_Borda_Flow: isentropic_face_root NOT converged'
    end subroutine isentropic_face_root


    subroutine bisect_M1_unchoked ( ok_out )
      logical, intent(out) :: ok_out
      real(R8) :: M1_lo, M1_hi, M1_mid, res_lo, res_hi, res_mid
      integer  :: it
      M1_hi  = 1.d0
      res_hi = P1 - momentum_P2()
      call solve_at_M1 ( M1_MIN, ok_out )
      if ( .not. ok_out ) return
      M1_lo  = M1_MIN
      res_lo = P1 - momentum_P2()
      if ( res_lo * res_hi > 0.d0 ) then
        write(*,*) 'BC_Nozzle_Borda_Flow: unchoked but no M1 bracket, keeping choked'
        call solve_at_M1 ( 1.d0, ok_out )
        return
      endif
      do it = 1, MAX_BISECT
        M1_mid = 0.5d0 * ( M1_lo + M1_hi )
        call solve_at_M1 ( M1_mid, ok_out )
        if ( .not. ok_out ) return
        res_mid = P1 - momentum_P2()
        if ( res_mid * res_lo < 0.d0 ) then
          M1_hi = M1_mid;  res_hi = res_mid
        else
          M1_lo = M1_mid;  res_lo = res_mid
        endif
        if ( abs( M1_hi - M1_lo ) < 1.d-5 ) exit
      enddo
      call solve_at_M1 ( 0.5d0 * ( M1_lo + M1_hi ), ok_out )
      if ( ok_out ) then
        if ( abs( P1 - momentum_P2() ) > P1P2_TOL * P1 ) then
          write(*,*) 'BC_Nozzle_Borda_Flow: unchoked bisection residual above tol'
        endif
      endif
    end subroutine bisect_M1_unchoked


    function momentum_P2 () result ( P2 )
      real(R8) :: P2
      P2 = ( mdot_A3 * ( abs( u3 ) - v1 ) - P1 * A1_A3 + p3 ) &
           / ( 1.d0 - A1_A3 )
    end function momentum_P2


    function station3_is_supersonic () result ( supersonic )
      logical :: supersonic
      real(R8) :: gamma3, a3
      gamma3 = cp3 / ( cp3 - Rgas_inj )
      a3 = sqrt ( gamma3 * Rgas_inj * T3 )
      supersonic = ( abs( u3 ) >= M3_SUP_TH * a3 )
    end function station3_is_supersonic

  end subroutine BC_Nozzle_Borda_Flow


end module Q2D_Lib_BC_Fluxes_Inflow_Nozzle
