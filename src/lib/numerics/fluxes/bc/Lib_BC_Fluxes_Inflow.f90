module Q2D_Lib_BC_Fluxes_Inflow
  ! Inflow conditions: one dedicated handler per case code.
  !
  ! 4 handlers (one per BC case):
  !   * BC_Inlet_StagnCond           — BC 401 (T0+p0) + 402 (T0+p0(t)), subsonic stagnation.
  !   * BC_Inlet_MassFlux_T0         — BC 403 (mdot+T0), imposed mass flux + total T, subsonic.
  !   * BC_Inlet_MassFlux_T          — BC 404 (mdot+T_static), imposed mass flux + static T, subsonic.
  !   * BC_Inlet_Supersonic_Static   — BC 405 (M+T_static+p_static), supersonic.
  !
  ! 5 shared private helpers:
  !   * Setup_Inflow_Geometry        — computes modfm/Int_*/Normal/int_l/t_Vec/Bound_Prim/Riem.
  !   * Compute_Inflow_Direction     — interprets the 'normal' sentinel of alpha (huge(R8)).
  !   * Stagn_Inflow                 — iterative discriminant kernel (cases 401/402).
  !   * Supersonic_Inflow            — Riemann kernel with a 3-cell stencil (case 405).
  !   * Compute_Flux_from_Face       — assembles the convective flux from a face state.
  !
  ! Strict 2D reduction (no out-of-plane component w):
  !   * no Km metric, no Int_k, no Face_k.
  !   * Normal(2), t_Vec(2), no out-of-plane binormal vector b_Vec (E2 ≡ 0).
  !   * velocity nu:nv (no nw), no Sup_Prim(nw) / Face_Prim(nw).
  !   * no out-of-plane elevation angle beta: cos(β)=1, sin(β)=0 by construction.
  !   * int_l is the length of the 2D interface (analogue of the 3D face area).
  !
  ! NB: no BC_Mach < 0 sentinel for the imposed mass flux: the dispatch happens
  ! by explicit case code 403/404. The BC_Mach=-20 sentinel (Borda) remains because
  ! it routes to the SEPARATE case 420 (BC_Inflow_Nozzle) and does not interfere
  ! with this module.
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_Limiters,           only: rlimiter
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index
  use Q2D_Lib_RotatingFrame,      only: RF_Convert_BC_Inflow
  use Q2D_Config_Types_m,         only: obj_rot
  use Q2D_Mod_Riemann,            only: Riemann
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, f_ss, f_gamma, f_Rtot, &
                                        f_tabT, h_tab, cp_tab, Ri_tab

  implicit none
  private

  public :: BC_Inlet_StagnCond
  public :: BC_Inlet_MassFlux_T0
  public :: BC_Inlet_MassFlux_T
  public :: BC_Inlet_Supersonic_Static
  ! Public surface — kernels reused by
  ! BC_FarField and possibly by other handlers.
  public :: Stagn_Inflow
  public :: Setup_Inflow_Geometry

contains


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 401 / 402: subsonic inflow — total temperature T0 and total pressure
  ! p0 (stagnation). Flow direction from BC_alpha (or face-normal sentinel).
  !
  ! RotatingFrame conversion APPLIED here: if obj_rot%enabled,
  ! T0_loc/p0_loc/alpha_loc are converted from absolute to relative with
  ! BC_Mach_abs=0d0 fixed (reservoir at rest in the laboratory frame;
  ! the total quantities in the relative frame include the kinetic energy
  ! induced by v_frame).
  subroutine BC_Inlet_StagnCond ( Bm, Im, Jm, Fm, Blk, BC_T0, BC_p0, BC_rel_fac, &
                                   BC_alpha, BC_ci, BC_RANS, BC_mdot, error )
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Fm
    real(R8), intent(in)    :: BC_T0, BC_p0, BC_rel_fac, BC_alpha
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    real(R8), intent(out)   :: BC_mdot
    integer,  intent(inout) :: error
    ! Local variables
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j
    real(R8) :: Normal(2), int_l, t_Vec(2), BC_Sign, alpha
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: T0_loc, p0_loc, alpha_loc, Face_Rgas_loc, Mach_rel_dummy
    real(R8) :: Flux(nprim)

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Normal, int_l, t_Vec, BC_Sign,    &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,&
                                  Bound_Sound, Bound_Gamma, Riem )

    ! Reject if the flow is directed outward (pure-inflow BC).
    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    ! Two-step scheme: `alpha` is the raw output, `alpha_loc` is modifiable
    ! by RF_Convert_BC_Inflow (intent inout).
    call Compute_Inflow_Direction ( Normal, BC_alpha, alpha )

    T0_loc    = BC_T0
    p0_loc    = BC_p0
    alpha_loc = alpha

    ! RotatingFrame conversion: only for the stagnation-condition inflow
    ! (BC_Mach_abs=0d0 fixed: reservoir at rest in the laboratory frame;
    ! the relative total quantities include the kinetic energy of v_frame).
    if (obj_rot%enabled) then
      Face_Rgas_loc = f_Rtot(BC_ci)
      call RF_Convert_BC_Inflow ( Blk, Fm, Im, Jm, 0.0d0, Bound_Gamma, &
                                   Face_Rgas_loc, T0_loc, p0_loc, alpha_loc, &
                                   Mach_rel_dummy )
    end if

    Flux = 0d0
    call Stagn_Inflow ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem, alpha_loc, &
                         BC_ci, BC_RANS, T0_loc, p0_loc, BC_rel_fac, BC_Sign,        &
                         t_Vec, Normal, int_l, Flux, BC_mdot )

    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + modfm2 * Flux


  end subroutine BC_Inlet_StagnCond


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 403: inflow — total temperature T0 + imposed mass flux g [kg/m²/s].
  ! Pressure extrapolated from the interior.
  ! NO RotatingFrame conversion: the relation between imposed mass flux and
  ! absolute/relative T0 is not trivial, and mdot is a frame invariant for
  ! concentric walls.
  subroutine BC_Inlet_MassFlux_T0 ( Bm, Im, Jm, Fm, Blk, BC_T0, BC_g, BC_rel_fac, &
                                     BC_alpha, BC_ci, BC_RANS, error )
    ! No BC_mdot output argument: the output of Compute_Flux_from_Face goes
    ! into the local variable `fmass` (discarded). An intent(out) here would
    ! alias the specification datum `domain%bc(i)%mdot` at the dispatcher:
    ! the output would overwrite the input between iterations, zeroing the
    ! specification down to a division by zero.
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Fm
    real(R8), intent(in)    :: BC_T0, BC_g, BC_rel_fac, BC_alpha
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    integer,  intent(inout) :: error
    ! Local variables
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, s
    real(R8) :: Normal(2), int_l, t_Vec(2), BC_Sign, fmass
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, alpha_loc
    real(R8) :: E1, cp3, Rgas3, T3, rho3, p3, Un3, Ut3
    real(R8) :: Face_Prim(nprim), Face_rho, Face_Rgas, Flux(nprim)
    real(R8) :: XA, XB

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Normal, int_l, t_Vec, BC_Sign,    &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,&
                                  Bound_Sound, Bound_Gamma, Riem )

    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    call Compute_Inflow_Direction ( Normal, BC_alpha, alpha_loc )

    ! Ratio of 2D direction cosines: in 2D the binormal vector is
    ! b_Vec ≡ (0,0,b3), so E2=0 and only E1 = -XA/XB remains.
    XA = Normal(2) - Normal(1) * tan(alpha_loc)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha_loc)
    E1 = -XA / XB

    ! Gas properties from the species composition (T0 used by the solver as total temperature).
    cp3   = 0d0
    Rgas3 = 0d0
    do s = 1, nsc
      cp3   = cp3   + BC_ci(s) * f_tabT(BC_T0, s, cp_tab)
      Rgas3 = Rgas3 + BC_ci(s) * Ri_tab(s)
    end do

    ! T0 is the total temperature: the static T3 is obtained from
    !   T3 = T0 - u²/(2 cp)  with  u = g/ρ = g R T3/p3
    p3 = Bound_Prim(np)   ! pressure extrapolated from the interior
    T3 = (-cp3 + sqrt(cp3**2 + 2d0*cp3*BC_T0*(BC_g*Rgas3/p3)**2)) * (p3 / (BC_g*Rgas3))**2

    rho3 = p3 / (Rgas3 * T3)
    Un3  = BC_g / rho3
    Ut3  = Un3 * E1

    Face_Rgas        = Rgas3
    Face_rho         = rho3
    Face_Prim(1:nsc) = Face_rho * BC_ci
    Face_Prim(nu:nv) = Un3 * Normal + Ut3 * t_Vec
    Face_Prim(np)    = p3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    call Compute_Flux_from_Face ( Face_Prim, Face_rho, Face_Rgas, int_l, Normal, Flux, fmass )
    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + modfm2 * Flux

  end subroutine BC_Inlet_MassFlux_T0


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 404: inflow — static temperature T + imposed mass flux g [kg/m²/s].
  ! Pressure extrapolated from the interior.
  ! NB: here BC_T is the STATIC T (reuse of the bc_type field: bc%T0 holds T_static
  ! for case 404 — contextual interpretation via case code; Mod_BC_Fluxes
  ! passes bc%T0 as BC_T for case 404).
  subroutine BC_Inlet_MassFlux_T ( Bm, Im, Jm, Fm, Blk, BC_T, BC_g, BC_rel_fac, &
                                    BC_alpha, BC_ci, BC_RANS, error, always_in )
    ! Same convention as BC_Inlet_MassFlux_T0: mass flux in the local variable
    ! `fmass`, no BC_mdot output argument (cf. the aliasing note there).
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Fm
    real(R8), intent(in)    :: BC_T, BC_g, BC_rel_fac, BC_alpha
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    integer,  intent(inout) :: error
    logical,  intent(in), optional :: always_in   ! an SRM grain in its inflow mode (502, mode 2): BC_g is the magnitude
                                                  ! of the mass flux into the domain, and the face is not closed when
                                                  ! the gas of the cell is at rest or moves out
    ! Local variables
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, s
    real(R8) :: Normal(2), int_l, t_Vec(2), BC_Sign, fmass
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem, alpha_loc
    real(R8) :: E1, Rgas3, rho3, p3, Un3, Ut3
    real(R8) :: Face_Prim(nprim), Face_rho, Face_Rgas, Flux(nprim)
    real(R8) :: XA, XB

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Normal, int_l, t_Vec, BC_Sign,    &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,&
                                  Bound_Sound, Bound_Gamma, Riem )

    if (BC_Sign * Un >= 0d0) then
      if ( .not. present(always_in) ) then
        error = 1
        return
      else if ( .not. always_in ) then
        error = 1
        return
      end if
    end if

    call Compute_Inflow_Direction ( Normal, BC_alpha, alpha_loc )

    XA = Normal(2) - Normal(1) * tan(alpha_loc)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha_loc)
    E1 = -XA / XB

    Rgas3 = 0d0
    do s = 1, nsc
      Rgas3 = Rgas3 + BC_ci(s) * Ri_tab(s)
    end do

    p3   = Bound_Prim(np)   ! pressure extrapolated from the interior
    rho3 = p3 / (Rgas3 * BC_T)
    Un3  = BC_g / rho3
    if ( present(always_in) ) then
      if ( always_in ) Un3 = - BC_Sign * abs(BC_g) / rho3   ! into the domain on every face
    end if
    Ut3  = Un3 * E1

    Face_Rgas        = Rgas3
    Face_rho         = rho3
    Face_Prim(1:nsc) = Face_rho * BC_ci
    Face_Prim(nu:nv) = Un3 * Normal + Ut3 * t_Vec
    Face_Prim(np)    = p3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    call Compute_Flux_from_Face ( Face_Prim, Face_rho, Face_Rgas, int_l, Normal, Flux, fmass )
    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + modfm2 * Flux


  end subroutine BC_Inlet_MassFlux_T


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 405: supersonic inflow — Mach + static temperature T + static pressure p.
  ! NB: BC_T is T_static, BC_p is p_static (reuse of the bc_type field: bc%T0=T_stat,
  ! bc%pamb=p_stat for case 405 — contextual interpretation via case code).
  subroutine BC_Inlet_Supersonic_Static ( Bm, Im, Jm, Fm, Blk, BC_Mach, BC_T, BC_p, &
                                           BC_rel_fac, BC_alpha, BC_ci, BC_RANS,    &
                                           BC_mdot, error )
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Fm
    real(R8), intent(in)    :: BC_Mach, BC_T, BC_p, BC_rel_fac, BC_alpha
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    real(R8), intent(out)   :: BC_mdot
    integer,  intent(inout) :: error
    ! Local variables
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j
    real(R8) :: Normal(2), int_l, t_Vec(2), BC_Sign
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: Flux(nprim)

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Normal, int_l, t_Vec, BC_Sign,    &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas,&
                                  Bound_Sound, Bound_Gamma, Riem )

    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    ! Supersonic_Inflow treats BC_p/BC_T as static quantities; blending with
    ! the interior + 3-cell stencil + Riemann solver handled in the kernel.
    Flux = 0d0
    call Supersonic_Inflow ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p, &
                              BC_T, BC_ci, BC_RANS, BC_rel_fac, BC_alpha, Normal, int_l,  &
                              Flux, BC_mdot )

    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + modfm2 * Flux


  end subroutine BC_Inlet_Supersonic_Static


  !! Low-level routines called by the high-level BC handlers above.
  !! They operate directly on the primitive arrays; they ignore the block structure.


  !─────────────────────────────────────────────────────────────────────────────
  ! Supersonic_Inflow (low-level kernel) — supersonic inflow for BC 405.
  ! Strict 2D reduction: no nw component, no out-of-plane sin/cos(beta).
  subroutine Supersonic_Inflow ( Prim2, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p0, BC_T0, BC_ci, &
                                  BC_RANS, BC_rel_fac, BC_alpha, Normal, int_l, Flux, Fmass )
    implicit none
    integer,  intent(in)  :: modfm, modfm1, modfm2
    real(R8), intent(in)  :: Prim2(nprim), Int_Prim(nprim), BC_Mach, BC_p0, BC_T0, BC_rel_fac
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(1:), BC_alpha, Normal(2), int_l
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local variables
    integer  :: s
    real(R8) :: Sup_Prim(nprim), Sup_T, Sup_Rgas, Sup_rho, Sup_Sound, Prim1(nprim), Prim3(nprim)
    real(R8), dimension(nprim) :: Diff32, Diff21, Slope, Face_Prim, Prim_L, Prim_R
    real(R8) :: rho_L, rho_R, Rgas_L, Rgas_R, Sound_L, Sound_R, F_r, F_u, F_v, F_E
    real(R8) :: su, Sel_L, Sel_R
    real(R8) :: flow_Rgas, flow_Sound, flow_Mach, BC_Mach_local, BC_p0_local
    real(R8) :: alpha_eff

    ! Flow properties in the interior cell.
    flow_Rgas  = f_Rtot(Prim2(1:nsc))
    flow_Sound = f_ss  (Prim2(1:nsc), Prim2(np), sum(Prim2(1:nsc)), flow_Rgas)
    flow_Mach  = norm2 (Prim2(nu:nv)) / flow_Sound

    ! Mach and static p at the boundary, blended toward the interior via BC_rel_fac.
    BC_Mach_local = BC_rel_fac * BC_Mach + (1 - BC_rel_fac) * flow_Mach
    BC_p0_local   = BC_rel_fac * BC_p0   + (1 - BC_rel_fac) * Prim2(np)

    ! Inflow direction:
    ! 'normal' sentinel huge(R8) → outgoing atan(Normal(2)/(Normal(1)+ε));
    ! otherwise BC_alpha is a literal value (parse_dir_tok in IO_BC produces
    ! huge(R8) for the 'normal' token; bc%alpha=0 from numeric input is a
    ! literal zero angle, NOT a sentinel).
    ! Cf. Compute_Inflow_Direction.
    call Compute_Inflow_Direction(Normal, BC_alpha, alpha_eff)

    ! Supersonic inflow: BC_Mach imposed at the boundary; BC_T0 and BC_p0 are static quantities.
    Sup_Prim(np)     = BC_p0_local
    Sup_T            = BC_T0
    Sup_Rgas         = sum(BC_ci * Ri_Tab)
    Sup_rho          = Sup_Prim(np) / (Sup_Rgas * Sup_T)
    Sup_Prim(1:nsc)  = Sup_rho * BC_ci
    Sup_Sound        = f_ss(Sup_Prim(1:nsc), Sup_Prim(np), Sup_rho, Sup_Rgas)
    Sup_Prim(nu)     = BC_Mach_local * Sup_Sound * cos(alpha_eff)
    Sup_Prim(nv)     = BC_Mach_local * Sup_Sound * sin(alpha_eff)
    if (model==2) Sup_Prim(nt:nprim) = Sup_rho * BC_RANS

    ! Reconstruction with a 3-cell stencil.
    Prim1 = modfm * Sup_Prim + modfm1 * Int_Prim
    Prim3 = modfm * Int_Prim + modfm1 * Sup_Prim

    Diff32 = Prim3 - Prim2
    Diff21 = Prim2 - Prim1
    do s = 1, nprim
      Slope(s) = rlimiter(Diff32(s), Diff21(s))
    end do
    Face_Prim = Prim2 + 0.5d0 * Slope * modfm2

    Prim_L = modfm * Sup_Prim + modfm1 * Face_Prim
    call co_rotot_Rtot(Prim_L(1:nsc), rho_L, Rgas_L)
    Sound_L = f_ss(Prim_L(1:nsc), Prim_L(np), rho_L, Rgas_L)

    Prim_R = modfm * Face_Prim + modfm1 * Sup_Prim
    call co_rotot_Rtot(Prim_R(1:nsc), rho_R, Rgas_R)
    Sound_R = f_ss(Prim_R(1:nsc), Prim_R(np), rho_R, Rgas_R)

    ! beta = 1d0: full HLLC weight in the blending of the SD-aware solvers
    ! (HLLE++/HLLCSD); the classic Riemann solvers ignore the argument.
    call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(np), Sound_L, rho_L, &
                   Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(np), Sound_R, rho_R, &
                   1d0, Normal(1), Normal(2), F_r, F_u, F_v, F_E )

    su    = sign(1d0, F_r)
    Sel_L = (1d0 + su) / 2d0
    Sel_R = (su  - 1d0) / 2d0

    Fmass       = F_r * int_l
    Flux        = 0d0
    Flux(1:nsc) = Fmass * (Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R)
    Flux(nu)    = F_u * int_l
    Flux(nv)    = F_v * int_l
    Flux(np)    = F_E * int_l

    if (model==2) Flux(nt:nprim) = Fmass * (Sel_L*Prim_L(nt:nprim)/rho_L - Sel_R*Prim_R(nt:nprim)/rho_R)

  end subroutine Supersonic_Inflow


  !─────────────────────────────────────────────────────────────────────────────
  ! Stagn_Inflow (low-level kernel) — iterative discriminant solver
  ! for BC 401/402 (subsonic stagnation). Handles only the stagnation-condition
  ! path: the imposed mass flux is handled by the separate handlers
  ! BC_Inlet_MassFlux_T0/T. Strict 2D reduction (no b_Vec/E2/Ub3, no
  ! nw component).
  subroutine Stagn_Inflow ( Prim, Un, Sound, Gamma, Riem, alpha, BC_ci, BC_RANS, &
                             BC_T0, BC_p0, BC_rel_fac, BC_Sign, t_Vec, Normal, int_l, &
                             Flux, Fmass )
    implicit none
    real(R8), intent(in)  :: Prim(nprim), Un, Sound, Gamma, Riem, alpha
    real(R8), intent(in)  :: BC_ci(nsc), BC_RANS(1:), BC_T0, BC_p0, BC_rel_fac, BC_Sign
    real(R8), intent(in)  :: t_Vec(2), Normal(2), int_l
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local variables
    integer  :: s, iter_discr
    real(R8) :: Face_Rgas, Sound0, Coef_D, Coef_C, YA, YB, YC, Discr, Un_Sign, Un3, Sound3
    real(R8) :: p3, Ut3, BC_p0_try
    real(R8) :: Face_T, Face_Sound, Face_Prim(nprim), Face_rho, Face_Enthalpy, Face_Un
    real(R8) :: XA, XB, E1

    ! Ratio of 2D direction cosines: the binormal is out of plane, so E2≡0 and E1 = -XA/XB remains.
    XA = Normal(2) - Normal(1) * tan(alpha)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha)
    E1 = -XA / XB

    ! Iterative discriminant solver (subsonic stagnation): if the discriminant
    ! is negative, p0 is relaxed toward the interior pressure via BC_rel_fac.
    Discr      = -1.d0
    BC_p0_try  = BC_p0
    iter_discr = 0
    do while (discr < 0d0)
      iter_discr = iter_discr + 1
      Face_Rgas = sum(BC_ci * Ri_tab)
      Sound0    = sqrt(Gamma * Face_Rgas * BC_T0)
      Coef_D    = (1d0 + E1**2) * (Gamma - 1d0) / 2d0
      Coef_C    = (BC_p0_try / Prim(np))**((Gamma - 1d0) / Gamma)
      Coef_C    = Coef_C * ((Sound / Sound0)**2)
      YA        = ((Gamma - 1d0) / 2d0)**2 + Coef_C * Coef_D
      YB        = BC_Sign * Riem * ((Gamma - 1d0) / 2d0)
      YC        = Riem**2 - Coef_C * Sound0**2
      Discr     = YB**2 - YA * YC
      if ((discr < 0d0) .and. (abs(discr) < 1d-2)) Discr = 0d0
      if ((discr < 0d0) .and. (abs(discr) > 1d-2)) then
        BC_p0_try = BC_rel_fac * BC_p0_try + (1 - BC_rel_fac) * Prim(np)
      end if
      if (iter_discr > 1000) then
        Discr = 0d0
        exit
      end if
    end do

    Un_Sign         = sign(1d0, Un)
    Un3             = (YB + Un_Sign * sqrt(Discr)) / YA
    Sound3          = sqrt(Sound0**2 - Coef_D * Un3**2)
    p3              = BC_p0_try * (Sound3 / Sound0)**(2d0 * Gamma / (Gamma - 1d0))
    Ut3             = Un3 * E1
    Face_Sound      = Sound3
    Face_Prim(np)   = p3
    Face_rho        = Gamma * Face_Prim(np) / Face_Sound**2
    Face_Prim(1:nsc) = BC_ci * Face_rho
    Face_Prim(nu:nv) = Un3 * Normal + Ut3 * t_Vec
    Face_Un          = Un3
    if (model==2) Face_Prim(nt:nprim) = BC_RANS * Face_rho

    Face_Enthalpy = 0d0
    Face_T = Face_Prim(np) / (Face_Rgas * Face_rho)
    Flux   = 0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT(Face_T, s, h_tab)
      Flux(s)       = Face_Prim(s) * Face_Un * int_l
    end do
    Fmass       = Face_rho * Face_Un * int_l
    Flux(nu:nv) = Fmass * Face_Prim(nu:nv) + Face_Prim(np) * int_l * Normal
    Flux(np)    = Fmass * (0.5d0 * sum(Face_Prim(nu:nv)**2) + Face_Enthalpy)
    if (model==2) Flux(nt:nprim) = Face_Prim(nt:nprim) / Face_rho * Fmass

  end subroutine Stagn_Inflow


  !─────────────────────────────────────────────────────────────────────────────
  ! Setup_Inflow_Geometry — geometric/metric setup common to the inflow handlers.
  ! Strict 2D reduction.
  pure subroutine Setup_Inflow_Geometry ( Blk, Im, Jm, Fm,                  &
                                      modfm, modfm1, modfm2, modfm3,        &
                                      Int_i, Int_j,                         &
                                      Normal, int_l, t_Vec, BC_Sign,        &
                                      Bound_Prim, Int_Prim,                 &
                                      Un, Bound_rho, Bound_Rgas,            &
                                      Bound_Sound, Bound_Gamma, Riem )
    implicit none
    type(Q2D_block_type), intent(in)  :: Blk
    integer,  intent(in)              :: Im, Jm, Fm
    integer,  intent(out)             :: modfm, modfm1, modfm2, modfm3
    integer,  intent(out)             :: Int_i, Int_j
    real(R8), intent(out)             :: Normal(2), int_l, t_Vec(2), BC_Sign
    real(R8), intent(out)             :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8), intent(out)             :: Un, Bound_rho, Bound_Rgas
    real(R8), intent(out)             :: Bound_Sound, Bound_Gamma, Riem
    ! Local variables
    integer  :: Dir, Face_i, Face_j
    real(R8) :: t_Mod

    call Compute_Modfm(Fm, modfm, modfm1, modfm2, modfm3)

    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)

    call Face_Index(Fm, Dir, Im, Jm, Face_i, Face_j)

    Bound_Prim = Blk % P(:,Im,Jm)
    Int_Prim   = Blk % P(:,Int_i,Int_j)

    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % N
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l

    select case (Fm)
      case (1, 2)
        t_Vec = Blk % node(Face_i, Face_j  ) % c - &
                Blk % node(Face_i, Face_j-1) % c
      case (3, 4)
        t_Vec = Blk % node(Face_i,   Face_j) % c - &
                Blk % node(Face_i-1, Face_j) % c
    end select
    t_Mod = norm2(t_Vec)
    t_Vec = t_Vec / t_Mod

    BC_Sign = real(modfm2, R8)   ! -1 for faces 1,3 ; +1 for faces 2,4

    Un          = dot_product(Bound_Prim(nu:nv), Normal)
    call co_rotot_Rtot(Bound_Prim(1:nsc), Bound_rho, Bound_Rgas)
    Bound_Sound = f_ss   (Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas)
    Bound_Gamma = f_gamma(Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas)
    Riem        = Bound_Sound + BC_Sign * 0.5d0 * (Bound_Gamma - 1d0) * Un

  end subroutine Setup_Inflow_Geometry


  !─────────────────────────────────────────────────────────────────────────────
  ! Compute_Inflow_Direction — interprets alpha_in (with the 'normal' sentinel huge(R8))
  ! and returns the resolved angle alpha_out.
  ! 2D reduction: no out-of-plane beta angle.
  pure subroutine Compute_Inflow_Direction ( Normal, alpha_in, alpha_out )
    implicit none
    real(R8), intent(in)  :: Normal(2), alpha_in
    real(R8), intent(out) :: alpha_out

    ! 'normal' sentinel: parse_dir_tok returns huge(R8) to indicate "use the face normal".
    if (alpha_in >= 0.5_R8 * huge(1.0_R8)) then
      alpha_out = atan(Normal(2) / (Normal(1) + 1d-20))
    else
      alpha_out = alpha_in
    end if

  end subroutine Compute_Inflow_Direction


  !─────────────────────────────────────────────────────────────────────────────
  ! Compute_Flux_from_Face — assembles the convective flux vector from a known face state.
  ! Strict 2D reduction (no nw component).
  pure subroutine Compute_Flux_from_Face ( Face_Prim, Face_rho, Face_Rgas, int_l, Normal, Flux, Fmass )
    implicit none
    real(R8), intent(in)  :: Face_Prim(nprim), Face_rho, Face_Rgas, int_l, Normal(2)
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local variables
    integer  :: s
    real(R8) :: Face_Un, Face_T, Face_Enthalpy

    Face_Un       = dot_product(Face_Prim(nu:nv), Normal)
    Face_T        = Face_Prim(np) / (Face_Rgas * Face_rho)
    Fmass         = Face_rho * Face_Un * int_l
    Face_Enthalpy = 0d0
    Flux          = 0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT(Face_T, s, h_tab)
      Flux(s)       = Face_Prim(s) * Face_Un * int_l
    end do
    Flux(nu:nv) = Fmass * Face_Prim(nu:nv) + Face_Prim(np) * int_l * Normal
    Flux(np)    = Fmass * (0.5d0 * sum(Face_Prim(nu:nv)**2) + Face_Enthalpy)
    if (nprim > np) Flux(np+1:nprim) = Face_Prim(np+1:nprim) / Face_rho * Fmass

  end subroutine Compute_Flux_from_Face


end module Q2D_Lib_BC_Fluxes_Inflow
