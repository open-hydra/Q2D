!>@brief Module for gas-surface interaction (GSI) subroutines.
!>
!> Provides GSI_Melting (wired in BC 303 Wall_Melting, 1:1 mirror of the MOSE
!> peer) and the GSI_HDPE/PP/HTPB/Carbon pyrolysis routines (wired in BC 304).
!> GSI_Paraffin is a Q2D-extra superset of GSI_Melting (peculiarità #8: hybrid
!> rocket motor paraffin with SEBS-mixed properties + an optional cracking
!> branch); it is exposed here but NOT currently wired into any BC handler, and
!> its SEBS/cracking branches are runtime-dead-by-default (y_SEBS=0,
!> cracking=.false.) so it reduces bit-for-bit to GSI_Melting under every
!> runnable configuration.
!>
!> The material parameters below stay as module-level defaults: there is NO
!> Setup_GSI in Q2D_new (the legacy INI parser of 'paraffin-prop'/'HDPE-prop'/
!> 'PP-prop'/'HTPB-prop' was not re-wired in the refactor). This is a MOSE-shared
!> state — archive-MOSE carries the identical "TODO: Refactor GSI material
!> properties to registry" (Register_Physics.f90) and likewise has no Setup_GSI.
!> Restoring configurability (Setup_GSI + gsi_t + reg%add) is deferred post-v1.0
!> as a coordinated open-hydra decision (TODO B-2/G-7/C-11.e/G-10.q), not a
!> unilateral Q2D wiring (Legge #5 biunivocità).
module Q2D_Lib_GSI
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Global_m,    only: nsc

  implicit none
  private

  public :: GSI_Melting   ! simplified variant (no SEBS, no cracking).
  public :: GSI_Paraffin  ! Q2D-extra — SEBS-mixed + optional cracking branch.
  public :: GSI_Carbon
  public :: GSI_HDPE
  public :: GSI_HTPB
  public :: GSI_PP

  ! Material constants exposed (default-public surface):
  ! consumed by Q2D_Lib_BC_Fluxes_Wall_Pyrolysis to compute the material-specific
  ! q_pyro (HDPE/PP/HTPB) during the Newton-Raphson energy balance.
  public :: cp_HDPE, Dh_HDPE, Ti_HDPE
  public :: cp_PP,   Dh_PP,   Ti_PP
  public :: cp_HTPB, Dh_HTPB, Ti_HTPB

  ! Paraffin / Melting (+ Q2D extras)
  real(R8) :: T_i    = 298.15d0             ! initial temperature of paraffin
  real(R8) :: cp_wax = 1946.0310193884d0    ! NIST: https://webbook.nist.gov/cgi/cbook.cgi?ID=C544854&Units=SI&Mask=6EF
  real(R8) :: Dh_wax = 1.698285789316d+5    ! NIST: same source as above
  real(R8) :: y_SEBS = 0d0                  ! SEBS mass fraction in fuel (default 0 → recovers GSI_Melting)
  integer  :: iFUEL  = 10                   ! standard set for singhWC32 (C32H66)
  integer  :: iSEBS  = 2                    ! C2H4 species index

  ! Q2D-extra (absent from the MOSE peer):
  real(R8) :: f_SEBS  = 1.7d0               ! SEBS cp/Dh = HDPE cp/Dh * f_SEBS
  logical  :: cracking = .false.            ! cracking of paraffin into ethylene+H2 at the surface

  ! Carbon
  real(R8), parameter :: wmc  = 12.0107d0   ! Carbon molecular weight
  real(R8), parameter :: pref = 1.01325d5   ! Atmospheric reference pressure

  ! HDPE
  real(R8) :: cp_HDPE = 1255.2d0
  real(R8) :: Dh_HDPE = 2.72d+6
  real(R8) :: A_HDPE  = 4588.8d3
  real(R8) :: Ea_HDPE = 251039.73d3
  real(R8) :: Ti_HDPE = 298.15d0

  ! PP
  real(R8) :: cp_PP   = 1700d0
  real(R8) :: Dh_PP   = 2.4823d6
  real(R8) :: A_PP    = 2.12d15
  real(R8) :: Ea_PP   = 2.12d8
  real(R8) :: Ti_PP   = 298.15d0

  ! HTPB
  real(R8) :: cp_HTPB  = 1632d0
  real(R8) :: Dh_HTPB  = 1.1d6
  real(R8) :: Ti_HTPB  = 298.15d0
  real(R8) :: A1_HTPB  = 3.965d0
  real(R8) :: Ea1_HTPB = 55.8564d6
  real(R8) :: A2_HTPB  = 11.04d-3
  real(R8) :: Ea2_HTPB = 20.54344d6
  real(R8) :: Ts_HTPB  = 722d0
  real(R8) :: rho_HTPB = 960d0

contains

  !>@brief Energy-balance melting GSI.
  !>
  !> Simplified variant: no SEBS, no cracking.
  !> For the paraffin-superset case (SEBS + cracking) use GSI_Paraffin (Q2D-extra).
  subroutine GSI_Melting ( qw, Tw, omega, mdot )
    implicit none
    real(R8), intent(in)  :: qw, Tw
    real(R8), intent(out) :: omega(nsc), mdot
    real(R8) :: cp_solid, Dh_melting, Dh_conduction, heat_absorbed, mdot_wax

    cp_solid   = cp_wax
    Dh_melting = Dh_wax

    ! energy balance
    Dh_conduction = cp_solid*(Tw - T_i)        ! heat of conduction in solid (J/kg)
    heat_absorbed = Dh_melting + Dh_conduction

    ! mass flux computed as a result of the energy balance at the grain surface
    mdot = qw / heat_absorbed                  ! (kg/(s m^2))

    mdot_wax = (1.0d0 - y_SEBS)*mdot

    !% Species source term
    omega = 0.0d0
    ! Injection of paraffin as C32H66 and ethylene for the SEBS
    if (y_SEBS > 0.0d0) omega(iSEBS) = y_SEBS*mdot   ! C2H4
    omega(iFUEL) = mdot_wax                          ! C32H66

  end subroutine GSI_Melting


  !>@brief Energy-balance melting GSI for paraffin propellant (Q2D-extra).
  !>
  !> Superset of GSI_Melting with two Q2D-only extensions:
  !>  (i)  SEBS-mixed solid properties: cp_solid = (1-y_SEBS)*cp_wax + y_SEBS*cp_SEBS,
  !>       Dh_melting = (1-y_SEBS)*Dh_wax + y_SEBS*Dh_SEBS, where cp_SEBS = f_SEBS*cp_HDPE
  !>       and Dh_SEBS = f_SEBS*Dh_HDPE.
  !>  (ii) optional `cracking` branch: paraffin → ethylene + H2 at the wall
  !>       (C32H66 → 16 C2H4 + H2) for extreme shear/heat at the fuel grain.
  !>
  !> Default y_SEBS=0 + cracking=.false. recovers GSI_Melting exactly.
  !>
  !> WHY the local UB-fix (cp_SEBS_loc, Dh_SEBS_loc): cp_SEBS/Dh_SEBS are declared
  !> but NOT initialised at module level; they are set only by Setup_GSI. Computing
  !> them locally from f_SEBS*cp_HDPE / f_SEBS*Dh_HDPE removes the ordering
  !> dependency and recovers the legacy formulas 1:1.
  subroutine GSI_Paraffin ( qw, Tw, omega, mdot )
    implicit none
    real(R8), intent(in)  :: qw, Tw
    real(R8), intent(out) :: omega(nsc), mdot
    real(R8) :: cp_solid, Dh_melting, Dh_conduction, heat_absorbed, mdot_wax
    real(R8) :: cp_SEBS_loc, Dh_SEBS_loc

    ! UB-fix: cp_SEBS / Dh_SEBS local.
    cp_SEBS_loc = f_SEBS * cp_HDPE
    Dh_SEBS_loc = f_SEBS * Dh_HDPE

    cp_solid   = (1.d0 - y_SEBS)*cp_wax + y_SEBS*cp_SEBS_loc
    Dh_melting = (1.d0 - y_SEBS)*Dh_wax + y_SEBS*Dh_SEBS_loc

    ! energy balance
    Dh_conduction = cp_solid*(Tw - T_i)        ! heat of conduction in solid (J/kg)
    heat_absorbed = Dh_melting + Dh_conduction

    ! mass flux computed as a result of the energy balance at the grain surface
    mdot = qw / heat_absorbed                  ! (kg/(s m^2))

    mdot_wax = (1.d0 - y_SEBS)*mdot

    !% Species source term
    omega = 0.d0
    if (cracking) then
      ! Cracking of paraffin into ethylene and hydrogen at the wall: C32H66 → 16 C2H4 + H2
      omega(6)  = 0.0044d0*mdot_wax                  ! H2
      omega(2)  = 0.9956d0*mdot_wax + y_SEBS*mdot    ! C2H4
    else
      ! Injection of paraffin as C32H66 and ethylene for the SEBS
      omega(2)  = y_SEBS*mdot                        ! C2H4
      omega(10) = mdot_wax                           ! C32H66
    endif

  end subroutine GSI_Paraffin


  subroutine GSI_Carbon ( roi, Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv, Ri_tab, Wm_tab
    implicit none
    real(R8), intent(in)    :: roi(nsc), Tw
    real(R8), intent(inout) :: omegadot(nsc), mdot
    real(R8) :: dotm_H2O, dotm_CO2, dotm_O2, dotm_OH, dotm_O
    real(R8) :: k5, k6, k7, k8, kH2O, kCO2, kOH, kO, Y_term
    real(R8) :: pi(nsc), RuTw, inv_sqrtTw
    integer, parameter :: iH=1, iC2H4=2, iOH=3, iCO=4, iCO2=5
    integer, parameter :: iH2=6, iH2O=7, iO2=8, iO=9

    omegadot = 0d0  ! initialization

    ! Constants
    RuTw = Runiv*Tw
    inv_sqrtTw = 1.0d0/sqrt(Tw)

    ! Partial pressure in atmospheres
    pi = roi*Tw*Ri_tab/pref

    ! Rates
    kH2O = 4.80d5  *exp(-288.00d6/RuTw)  ! H2O
    kCO2 = 9.00d3  *exp(-285.00d6/RuTw)  ! CO2
    kOH  = 3.61d2  *inv_sqrtTw
    kO   = 6.655d2 *inv_sqrtTw

    k5   = 2.40d3  *exp(-125.60d6/RuTw)  ! O2a
    k6   = 2.13d1  *exp(  17.17d6/RuTw)  ! O2b
    k7   = 5.35d-1 *exp( -63.64d6/RuTw)  ! O2c
    k8   = 1.81d7  *exp(-406.10d6/RuTw)  ! O2d

    ! Compute ybig term
    Y_term = 1d0 / ( 1d0 + k8/(k7*pi(iO2)) )

    ! H2O steam-gasification rate (C + H2O -> CO + H2): reads the H2O partial pressure.
    ! Q2D-correct vs legacy/MOSE `pi(iO2)`: legacy & MOSE carry a name/index mismatch
    ! (the slot indexed by their `iO2` parameter actually holds H2O in their species
    ! table); Q2D uses the physically-named `iH2O` consistently. Per user (Sess 36)
    ! this is a Q2D intentional correction, NOT a porting loss -> do not "restore" iO2.
    dotm_H2O = kH2O * sqrt(pi(iH2O))
    dotm_CO2 = kCO2 * sqrt(pi(iCO2))
    dotm_O2  = Y_term * (k5*pi(iO2))/(1.0d0 + k6*pi(iO2)) + k7*pi(iO2)*(1d0-Y_term)
    dotm_OH  = kOH*pi(iOH)
    dotm_O   = kO*pi(iO)
    ! Compute source terms
    omegadot(iO2)  = -(Wm_tab(iO2)/wmc)*(0.5d0*dotm_O2)                                         ! O2
    omegadot(iH2O) = -(Wm_tab(iH2O)/wmc)*(dotm_H2O)                                             ! H2O
    omegadot(iCO)  =  (Wm_tab(iCO)/wmc)*(2.d0*dotm_CO2 + dotm_H2O + dotm_O2 + dotm_OH + dotm_O) ! CO
    omegadot(iCO2) = -(Wm_tab(iCO2)/wmc)*(dotm_CO2)                                             ! CO2
    omegadot(iH2)  =  (Wm_tab(iH2)/wmc)*(dotm_H2O)                                              ! H2
    omegadot(iO)   = -(Wm_tab(iO)/wmc)*(dotm_O)                                                 ! O
    omegadot(iH)   =  (Wm_tab(iH)/wmc)*(dotm_OH)                                                ! H
    omegadot(iOH)  = -(Wm_tab(iOH)/wmc)*(dotm_OH)                                               ! OH

    ! output of dotm
    mdot = sum ( omegadot )

  end subroutine GSI_Carbon


  subroutine GSI_HDPE ( Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv
    implicit none
    real(R8), intent(in)    :: Tw
    real(R8), intent(inout) :: omegadot(nsc), mdot

    omegadot = 0d0
    mdot = A_HDPE * exp(-Ea_HDPE/(2*Runiv*Tw))
    omegadot(2) = mdot ! ethylene, C2H4

  end subroutine GSI_HDPE


  subroutine GSI_HTPB ( Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv
    implicit none
    real(R8), intent(in)    :: Tw
    real(R8), intent(inout) :: omegadot(nsc), mdot

    omegadot = 0d0
    if (Tw <= Ts_HTPB) then
      mdot = rho_HTPB * A1_HTPB * exp(-Ea1_HTPB/(Runiv*Tw))
    else
      mdot = rho_HTPB * A2_HTPB * exp(-Ea2_HTPB/(Runiv*Tw))
    end if
    omegadot(2) = mdot ! C4H6

  end subroutine GSI_HTPB


  subroutine GSI_PP ( Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv
    implicit none
    real(R8), intent(in)    :: Tw
    real(R8), intent(inout) :: omegadot(nsc), mdot
    real(R8) :: Arrh, RR

    omegadot = 0d0
    Arrh = Ea_PP/(Runiv*Tw)
    RR = sqrt ( 4.60517d0 * ( 1d0 - Ti_PP/Tw + Dh_PP/(cp_PP*Tw) ) - Dh_PP/(cp_PP*Tw) )
    mdot = 910d0*sqrt( A_PP * exp(-Arrh) * 1d0/Arrh * 0.2d0/(910d0*cp_PP) * 1/RR )
    omegadot(2) = mdot ! propylene, C3H6

  end subroutine GSI_PP


end module Q2D_Lib_GSI
