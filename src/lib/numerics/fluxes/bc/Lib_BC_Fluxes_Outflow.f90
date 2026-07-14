module Q2D_Lib_BC_Fluxes_Outflow
  ! BC_Outflow is a STANDALONE handler (it was previously a non-standalone
  ! helper called by BC_Inflow_Outflow).
  ! Computes Bound_Prim/rho/Sound/Un/Riem in-house from the block; calls the
  ! private Outflow kernel.
  !
  ! Rigid 2D reduction (No-w principle):
  !   * drop Km argument, drop Int_k, drop Face_k.
  !   * Normal(2), no Vel_T 3-component.
  !   * Velocity nu:nv (no nw), no Face_Prim(nw).
  !   * `int_l ↔ Area` rename.
  !
  ! Sentinel BC_pAmb=-3.14d0 (Q2D legacy artifact) replaced by the gate
  ! `BC_pexit > 0`. Functionality preserved: the outflow no-pAmb path → the
  ! Outflow kernel falls into the supersonic-extrapolation else-branch
  ! (extrapolation of all variables).
  !
  ! Q2D-extra preserved features:
  !   * Outflow under-relaxation: Face_Prim(np) = p_exit * BC_rel_fac +
  !     Bound_Prim(np) * (1 - BC_rel_fac) — a Q2D feature for numerical
  !     stability (the MOSE peer code uses p_exit directly without
  !     under-relaxation).
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, f_ss, f_gamma, f_tabT, h_tab

  implicit none
  private

  public :: BC_Outflow
  public :: Outflow   ! re-used by BC_FarField.

contains


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 406: Outflow — prescribed back-pressure BC_pexit (or extrapolation if
  ! BC_pexit ≤ 0 / supersonic).
  subroutine BC_Outflow ( Bm, Im, Jm, Fm, Blk, BC_pexit, BC_rel_fac, BC_mdot, error )
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Fm
    real(R8), intent(in)    :: BC_rel_fac, BC_pexit
    real(R8), intent(out)   :: BC_mdot
    integer,  intent(inout) :: error
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Dir, Face_i, Face_j
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim), Normal(2), int_l, BC_Sign, Un
    real(R8) :: Flux(nprim)

    error = 0

    call Compute_Modfm(Fm, modfm, modfm1, modfm2, modfm3)

    ! (Im,Jm): boundary cell; (Int_i,Int_j): interior cell.
    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)

    call Face_Index(Fm, Dir, Im, Jm, Face_i, Face_j)

    ! Boundary cell primitive variables.
    Bound_Prim = Blk % P(:,Im,Jm)
    Int_Prim   = Blk % P(:,Int_i,Int_j)

    ! Metric stuff.
    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % N
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l

    BC_Sign = real(modfm2, R8)   ! -1 for faces 1,3 ; +1 for faces 2,4
    Un      = dot_product(Bound_Prim(nu:nv), Normal)

    ! Reject if flow is directed inward (pure outflow BC).
    if (BC_Sign * Un < 0d0) then
      error = 1
      return
    end if

    ! BC_rel_fac propagated to the kernel for the Q2D-form under-relaxation.
    call Outflow ( Bound_Prim, BC_Sign, Un, Normal, int_l, BC_pexit, BC_rel_fac, &
                    Flux, BC_mdot )

    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + modfm2 * Flux

  end subroutine BC_Outflow


  !─────────────────────────────────────────────────────────────────────────────
  ! Outflow (private kernel) — BC_Outflow body with rigid 2D reduction
  ! (drop nw, Vel_T(2) instead of Vel_T(3), Area→int_l).
  subroutine Outflow ( Bound_Prim, BC_Sign, Un, Normal, int_l, BC_pexit, BC_rel_fac, &
                        Flux, BC_mdot )
    implicit none
    real(R8), intent(in)  :: Bound_Prim(nprim), BC_Sign, Un, Normal(2), int_l
    real(R8), intent(in)  :: BC_pexit, BC_rel_fac
    real(R8), intent(out) :: Flux(nprim), BC_mdot
    ! Local
    real(R8) :: Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: Vel_N(2), Vel_T(2), Mach_N, Face_Sound, Isentropic
    real(R8) :: Face_Un, Face_Prim(nprim), Face_rho, Face_Rgas, Face_Enthalpy, Face_T
    real(R8) :: Mach_abs, pstar, den, delta, p_exit
    integer :: s

    call co_rotot_Rtot(Bound_Prim(1:nsc), Bound_rho, Bound_Rgas)
    Bound_Sound = f_ss   (Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas)
    Bound_Gamma = f_gamma(Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas)

    ! Riemann invariant `a + BC_Sign·delta·u` directed towards the interface.
    Riem = Bound_Sound + BC_Sign * 0.5d0 * (Bound_Gamma - 1d0) * Un

    ! Interface values computed for the fluxes at the end.
    Vel_N    = Un * Normal
    Vel_T    = Bound_Prim(nu:nv) - Vel_N
    Mach_abs = norm2(Bound_Prim(nu:nv)) / Bound_Sound
    Mach_N   = abs(Un / Bound_Sound)

    if ((Mach_N < 1d0) .and. (BC_pexit > 0d0)) then
      ! Subsonic outflow with prescribed back-pressure.
      p_exit = BC_pexit
      delta  = 0.5d0 * (Bound_Gamma - 1.d0)
      den    = 1.d0 / (1 + delta)**(Bound_Gamma / (Bound_Gamma - 1.d0))
      pstar  = Bound_Prim(np) * den * (1 + delta * Mach_N**2)**(Bound_Gamma / (Bound_Gamma - 1.d0))
      if (pstar >= p_exit) p_exit = pstar
      Isentropic = (Bound_Gamma - 1d0) / (2d0 * Bound_Gamma)

      ! Q2D-extra under-relaxation of the imposed pressure, for transient/
      ! iterative numerical stability (the MOSE peer code uses
      ! Face_Prim(np) = p_exit directly).
      Face_Prim(np) = p_exit * BC_rel_fac + Bound_Prim(np) * (1.d0 - BC_rel_fac)

      ! Isentropic expansion to ambient pressure.
      Face_Sound       = Bound_Sound * ((Face_Prim(np) / Bound_Prim(np))**Isentropic)
      ! Normal velocity at interface using Riemann invariant.
      Face_Un          = (Riem - Face_Sound) / (BC_Sign * 0.5d0 * (Bound_Gamma - 1d0))
      Face_rho         = Bound_Gamma * Face_Prim(np) / Face_Sound**2
      Face_Prim(1:nsc) = Bound_Prim(1:nsc) / Bound_rho * Face_rho   ! mass fractions extrapolated
      Face_Prim(nu:nv) = Vel_T + Face_Un * Normal                   ! tangential extrapolated
      ! Gate `if (nprim > np)`: extrapolate ALL post-pressure variables
      ! (turbulence, soot, passive scalars), mirroring MOSE Outflow (nprim>np).
      ! Q2D legacy gated only the RANS block: `if (RANS) Face_Prim(nt:nprim) = ...`
      ! (Lib_BC.f90 BC_Outflow), missing soot/passive (np+1:nt-1) — a legacy
      ! oversight, not intentional preservation. (Legacy used the `RANS` logical
      ! flag, not `model==2`; `model` did not exist in legacy Lib_BC.)
      if (nprim > np) Face_Prim(np+1:nprim) = Bound_Prim(np+1:nprim) / Bound_rho * Face_rho
      Face_Rgas        = Bound_Rgas

    else
      ! Supersonic exit OR BC_pexit ≤ 0: extrapolation of all variables.
      Face_Prim  = Bound_Prim
      Face_Un    = Un
      Face_Sound = Bound_Sound
      Face_rho   = Bound_rho
      Face_Rgas  = Bound_Rgas
    end if

    Face_Enthalpy = 0d0
    Face_T        = Face_Prim(np) / (Face_Rgas * Face_rho)
    Flux          = 0d0
    do s = 1, nsc
      Face_Enthalpy = Face_Enthalpy + Face_Prim(s) / Face_rho * f_tabT(Face_T, s, h_tab)
      Flux(s)       = Face_Prim(s) * int_l * Face_Un
    end do
    BC_mdot     = Face_rho * int_l * Face_Un
    Flux(nu:nv) = Face_Prim(nu:nv) * BC_mdot + Face_Prim(np) * int_l * Normal
    Flux(np)    = BC_mdot * (0.5d0 * sum(Face_Prim(nu:nv)**2) + Face_Enthalpy)
    ! Gate `if (nprim > np)`: flux-out all post-pressure variables
    ! (turbulence + soot + passive).
    if (nprim > np) Flux(np+1:nprim) = Face_Prim(np+1:nprim) / Face_rho * BC_mdot

    ! NB: Mach_abs computed but not used locally (dead-variable pattern).

  end subroutine Outflow


end module Q2D_Lib_BC_Fluxes_Outflow
