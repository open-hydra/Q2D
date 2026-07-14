module Q2D_Lib_BC_Fluxes_FarField
  ! BC 407 ambient T0+p0+pback with a
  ! Riemann-based dispatch on the sign of the normal velocity (inflow
  ! Stagn_Inflow vs outflow back-pressure).
  !
  ! Rigid 2D reduction (No-w principle):
  !   * drop Km argument; drop beta usage in the Stagn_Inflow call (Q2D
  !     Stagn_Inflow 17-args vs the 3D peer's 18-args, beta dropped).
  !   * Normal(2), t_Vec(2), Bound_Prim(nprim) with velocity nu:nv only.
  !   * Area→int_l rename.
  !   * rank-3 residual.
  !   * RF_Convert_BC_Inflow Q2D-form 11-args (no Km, no beta).
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm

  implicit none
  private

  public :: BC_FarField

contains


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 407: Ambient BC — T0 + p0 (inflow) or back pressure BC_pback (outflow).
  ! Flow direction is determined from the interior velocity field.
  subroutine BC_FarField ( Bm, Im, Jm, Fm, Blk, BC_T0, BC_p0, BC_pback, &
                            BC_rel_fac, BC_alpha, BC_beta, BC_ci, BC_RANS, BC_mdot )
    use Q2D_Config_Types_m,        only: obj_rot
    use Q2D_Lib_RotatingFrame,     only: RF_Convert_BC_Inflow
    use FLINT_Lib_Thermodynamic,   only: f_Rtot
    use Q2D_Lib_BC_Fluxes_Inflow,  only: Stagn_Inflow, Setup_Inflow_Geometry
    use Q2D_Lib_BC_Fluxes_Outflow, only: Outflow
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Fm
    real(R8), intent(in)    :: BC_T0, BC_p0, BC_pback, BC_rel_fac, BC_alpha, BC_beta
    real(R8), intent(in)    :: BC_ci(nsc), BC_RANS(1:)
    real(R8), intent(out)   :: BC_mdot
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j
    real(R8) :: Normal(2), int_l, t_Vec(2), BC_Sign
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Rgas, Bound_Sound, Bound_Gamma, Riem
    real(R8) :: T0_loc, p0_loc, Mach_loc, alpha_loc, Face_Rgas_loc
    real(R8) :: Flux(nprim)

    ! NB: BC_beta `intent(in)` unused in the body — the No-w principle drops the
    ! out-of-plane elevation angle. The argument is preserved for signature
    ! alignment. ifx remark #7712 (unused intent(in)) is accepted (cf. the
    ! Lib_BC_Fluxes_Rotational Symmetry_Eul stub, same pattern).

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Normal, int_l, t_Vec, BC_Sign,    &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Rgas, &
                                  Bound_Sound, Bound_Gamma, Riem )

    T0_loc    = BC_T0
    p0_loc    = BC_p0
    Mach_loc  = 0.0d0
    alpha_loc = BC_alpha

    Flux = 0d0

    if (BC_Sign * Un >= 0d0) then
      ! Outflow: apply back-pressure outlet treatment (Q2D Outflow kernel takes BC_rel_fac).
      call Outflow ( Bound_Prim, BC_Sign, Un, Normal, int_l, BC_pback, BC_rel_fac, Flux, BC_mdot )
    else
      ! Inflow
      if (obj_rot%enabled) then
        Face_Rgas_loc = f_Rtot(BC_ci)
        call RF_Convert_BC_Inflow ( Blk=Blk, Fm=Fm, Im=Im, Jm=Jm,                &
                                     BC_Mach_abs=0.0d0, gamma=Bound_Gamma,         &
                                     Rgas=Face_Rgas_loc, T0_or_Tstat=T0_loc,       &
                                     p0_or_pstat=p0_loc, alpha=alpha_loc,          &
                                     Mach_rel=Mach_loc )
      end if
      call Stagn_Inflow ( Bound_Prim, Un, Bound_Sound, Bound_Gamma, Riem,    &
                            alpha_loc, BC_ci, BC_RANS, T0_loc, p0_loc,         &
                            BC_rel_fac, BC_Sign, t_Vec, Normal, int_l, Flux, BC_mdot )
    end if

    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + modfm2 * Flux

  end subroutine BC_FarField

end module Q2D_Lib_BC_Fluxes_FarField
