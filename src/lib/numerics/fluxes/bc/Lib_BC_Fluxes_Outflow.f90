module Q2D_Lib_BC_Fluxes_Outflow
  ! BC 406: outflow with prescribed back-pressure.
  ! BC_Outflow is a self-contained handler: it computes Bound_Prim/rho/Sound/Un/Riem
  ! in place from the block and invokes the private kernel Outflow.
  !
  ! 2D formulation:
  !   * 2-component normal, velocity nu:nv (no out-of-plane component);
  !   * 2-component Vel_T;
  !   * int_l: integrated face length (in place of the 3D area).
  !
  ! The `BC_pexit > 0` gate explicitly enables/disables the ambient pressure:
  ! with BC_pexit ≤ 0 the Outflow kernel falls back to the supersonic
  ! extrapolation else branch (extrapolation of all variables).
  !
  ! Under-relaxation of the imposed pressure, for numerical stability:
  !   Face_Prim(np) = p_exit * BC_rel_fac + Bound_Prim(np) * (1 - BC_rel_fac).
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, f_ss, f_gamma, f_tabT, h_tab

  implicit none
  private

  public :: BC_Outflow
  public :: Outflow   ! reused by BC_FarField.

contains


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 406: Outflow — prescribed back-pressure BC_pexit (or extrapolation if
  ! BC_pexit ≤ 0 / supersonic exit).
  subroutine BC_Outflow ( Bm, Im, Jm, Fm, Blk, BC_pexit, BC_rel_fac, BC_mdot, error )
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Fm
    real(R8), intent(in)    :: BC_rel_fac, BC_pexit
    real(R8), intent(out)   :: BC_mdot
    integer,  intent(inout) :: error
    ! Locals
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Dir, Face_i, Face_j
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim), Normal(2), int_l, BC_Sign, Un
    real(R8) :: Flux(nprim)

    error = 0

    call Compute_Modfm(Fm, modfm, modfm1, modfm2, modfm3)

    ! (Im,Jm): boundary cell; (Int_i,Int_j): interior cell.
    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)

    call Face_Index(Fm, Dir, Im, Jm, Face_i, Face_j)

    ! Primitive variables of the boundary cell.
    Bound_Prim = Blk % P(:,Im,Jm)
    Int_Prim   = Blk % P(:,Int_i,Int_j)

    ! Face metrics.
    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % N
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l

    BC_Sign = real(modfm2, R8)   ! -1 for faces 1,3 ; +1 for faces 2,4
    Un      = dot_product(Bound_Prim(nu:nv), Normal)

    ! Reject if the flow is directed inward (pure-outflow BC).
    if (BC_Sign * Un < 0d0) then
      error = 1
      return
    end if

    ! BC_rel_fac propagated to the kernel for under-relaxation of the imposed pressure.
    call Outflow ( Bound_Prim, BC_Sign, Un, Normal, int_l, BC_pexit, BC_rel_fac, &
                    Flux, BC_mdot )

    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + modfm2 * Flux

  end subroutine BC_Outflow


  !─────────────────────────────────────────────────────────────────────────────
  ! Outflow (private kernel) — computational body in 2D formulation
  ! (in-plane velocity nu:nv, tangential Vel_T(2), face length int_l).
  subroutine Outflow ( Bound_Prim, BC_Sign, Un, Normal, int_l, BC_pexit, BC_rel_fac, &
                        Flux, BC_mdot )
    implicit none
    real(R8), intent(in)  :: Bound_Prim(nprim), BC_Sign, Un, Normal(2), int_l
    real(R8), intent(in)  :: BC_pexit, BC_rel_fac
    real(R8), intent(out) :: Flux(nprim), BC_mdot
    ! Locals
    real(R8) :: Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: Vel_N(2), Vel_T(2), Mach_N, Face_Sound, Isentropic
    real(R8) :: Face_Un, Face_Prim(nprim), Face_rho, Face_Rgas, Face_Enthalpy, Face_T
    real(R8) :: Mach_abs, pstar, den, delta, p_exit
    integer :: s

    call co_rotot_Rtot(Bound_Prim(1:nsc), Bound_rho, Bound_Rgas)
    Bound_Sound = f_ss   (Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas)
    Bound_Gamma = f_gamma(Bound_Prim(1:nsc), Bound_Prim(np), Bound_rho, Bound_Rgas)

    ! Riemann invariant `a + BC_Sign·delta·u` directed toward the interface.
    Riem = Bound_Sound + BC_Sign * 0.5d0 * (Bound_Gamma - 1d0) * Un

    ! Interface values used by the fluxes computed at the end.
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

      ! Under-relaxation of the imposed pressure, for numerical stability
      ! in transients and in the iterative process.
      Face_Prim(np) = p_exit * BC_rel_fac + Bound_Prim(np) * (1.d0 - BC_rel_fac)

      ! Espansione isentropica alla pressione ambiente.
      Face_Sound       = Bound_Sound * ((Face_Prim(np) / Bound_Prim(np))**Isentropic)
      ! Velocità normale all'interfaccia dall'invariante di Riemann.
      Face_Un          = (Riem - Face_Sound) / (BC_Sign * 0.5d0 * (Bound_Gamma - 1d0))
      Face_rho         = Bound_Gamma * Face_Prim(np) / Face_Sound**2
      Face_Prim(1:nsc) = Bound_Prim(1:nsc) / Bound_rho * Face_rho   ! frazioni di massa estrapolate
      Face_Prim(nu:nv) = Vel_T + Face_Un * Normal                   ! componente tangenziale estrapolata
      ! Il gate `if (nprim > np)` estrapola TUTTA la banda degli scalari
      ! trasportati oltre la pressione (turbolenza, soot, passivi: np+1:nprim),
      ! non solo quella RANS nt:nprim, per non perdere soot/passivi (np+1:nt-1)
      ! quando il modello RANS è spento.
      if (nprim > np) Face_Prim(np+1:nprim) = Bound_Prim(np+1:nprim) / Bound_rho * Face_rho
      Face_Rgas        = Bound_Rgas

    else
      ! Uscita supersonica O BC_pexit ≤ 0: estrapolazione di tutte le variabili.
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
    ! Il gate `if (nprim > np)` scarica in flusso tutte le variabili oltre la
    ! pressione (turbolenza + soot + passivi).
    if (nprim > np) Flux(np+1:nprim) = Face_Prim(np+1:nprim) / Face_rho * BC_mdot

    ! NB: Mach_abs calcolato ma non usato localmente.

  end subroutine Outflow


end module Q2D_Lib_BC_Fluxes_Outflow
