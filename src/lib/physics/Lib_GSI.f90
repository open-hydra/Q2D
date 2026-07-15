!>@brief Modulo per le subroutine di interazione gas-superficie (GSI).
!>
!> Fornisce GSI_Melting (collegata alla BC 303 Wall_Melting) e le routine di
!> pirolisi GSI_HDPE/PP/HTPB/Carbon (collegate alla BC 304).
!> GSI_Paraffin è un superset di GSI_Melting (paraffina per motori a razzo
!> ibridi con proprietà miscelate SEBS + ramo opzionale di cracking); è esposta
!> qui ma NON è attualmente collegata ad alcun gestore BC, e i rami
!> SEBS/cracking sono disattivi a runtime per default (y_SEBS=0,
!> cracking=.false.), per cui in ogni configurazione eseguibile si riduce
!> bit per bit a GSI_Melting.
!>
!> I parametri dei materiali sottostanti sono default a livello di modulo:
!> non esiste una routine di setup che li renda configurabili da file INI.
module Q2D_Lib_GSI
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Global_m,    only: nsc

  implicit none
  private

  public :: GSI_Melting   ! variante semplificata (senza SEBS, senza cracking).
  public :: GSI_Paraffin  ! proprietà miscelate SEBS + ramo opzionale di cracking.
  public :: GSI_Carbon
  public :: GSI_HDPE
  public :: GSI_HTPB
  public :: GSI_PP

  ! Costanti dei materiali esposte:
  ! usate da Q2D_Lib_BC_Fluxes_Wall_Pyrolysis per calcolare il q_pyro specifico
  ! del materiale (HDPE/PP/HTPB) nel bilancio energetico Newton-Raphson.
  public :: cp_HDPE, Dh_HDPE, Ti_HDPE
  public :: cp_PP,   Dh_PP,   Ti_PP
  public :: cp_HTPB, Dh_HTPB, Ti_HTPB

  ! Paraffina / Melting
  real(R8) :: T_i    = 298.15d0             ! temperatura iniziale della paraffina
  real(R8) :: cp_wax = 1946.0310193884d0    ! NIST: https://webbook.nist.gov/cgi/cbook.cgi?ID=C544854&Units=SI&Mask=6EF
  real(R8) :: Dh_wax = 1.698285789316d+5    ! NIST: stessa fonte della riga precedente
  real(R8) :: y_SEBS = 0d0                  ! frazione di massa di SEBS nel combustibile (default 0 → recupera GSI_Melting)
  integer  :: iFUEL  = 10                   ! indice del C32H66 nel meccanismo Singh-WC32
  integer  :: iSEBS  = 2                    ! indice della specie C2H4

  ! Parametri delle estensioni SEBS/cracking di GSI_Paraffin:
  real(R8) :: f_SEBS  = 1.7d0               ! cp/Dh del SEBS = cp/Dh dell'HDPE * f_SEBS
  logical  :: cracking = .false.            ! cracking della paraffina in etilene+H2 alla superficie

  ! Carbonio
  real(R8), parameter :: wmc  = 12.0107d0   ! peso molecolare del carbonio
  real(R8), parameter :: pref = 1.01325d5   ! pressione atmosferica di riferimento

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

  !>@brief GSI di fusione a bilancio energetico.
  !>
  !> Variante semplificata: senza SEBS, senza cracking.
  !> Per il caso superset paraffina (SEBS + cracking) usare GSI_Paraffin.
  subroutine GSI_Melting ( qw, Tw, omega, mdot )
    implicit none
    real(R8), intent(in)  :: qw, Tw
    real(R8), intent(out) :: omega(nsc), mdot
    real(R8) :: cp_solid, Dh_melting, Dh_conduction, heat_absorbed, mdot_wax

    cp_solid   = cp_wax
    Dh_melting = Dh_wax

    ! bilancio energetico
    Dh_conduction = cp_solid*(Tw - T_i)        ! calore di conduzione nel solido (J/kg)
    heat_absorbed = Dh_melting + Dh_conduction

    ! flusso di massa risultante dal bilancio energetico alla superficie del grano
    mdot = qw / heat_absorbed                  ! (kg/(s m^2))

    mdot_wax = (1.0d0 - y_SEBS)*mdot

    !% Termine sorgente delle specie
    omega = 0.0d0
    ! iniezione di paraffina come C32H66 e di etilene per il SEBS
    if (y_SEBS > 0.0d0) omega(iSEBS) = y_SEBS*mdot   ! C2H4
    omega(iFUEL) = mdot_wax                          ! C32H66

  end subroutine GSI_Melting


  !>@brief GSI di fusione a bilancio energetico per propellente paraffinico.
  !>
  !> Superset di GSI_Melting con due estensioni:
  !>  (i)  proprietà del solido miscelate con SEBS: cp_solid = (1-y_SEBS)*cp_wax + y_SEBS*cp_SEBS,
  !>       Dh_melting = (1-y_SEBS)*Dh_wax + y_SEBS*Dh_SEBS, con cp_SEBS = f_SEBS*cp_HDPE
  !>       e Dh_SEBS = f_SEBS*Dh_HDPE.
  !>  (ii) ramo opzionale `cracking`: paraffina → etilene + H2 a parete
  !>       (C32H66 → 16 C2H4 + H2) per shear/riscaldamento estremi al grano.
  !>
  !> Il default y_SEBS=0 + cracking=.false. recupera esattamente GSI_Melting.
  !>
  !> Le proprietà del SEBS sono calcolate localmente (cp_SEBS_loc = f_SEBS*cp_HDPE,
  !> Dh_SEBS_loc = f_SEBS*Dh_HDPE): nessuna dipendenza da variabili di modulo
  !> non inizializzate né da un ordine di setup.
  subroutine GSI_Paraffin ( qw, Tw, omega, mdot )
    implicit none
    real(R8), intent(in)  :: qw, Tw
    real(R8), intent(out) :: omega(nsc), mdot
    real(R8) :: cp_solid, Dh_melting, Dh_conduction, heat_absorbed, mdot_wax
    real(R8) :: cp_SEBS_loc, Dh_SEBS_loc

    ! Proprietà SEBS locali, derivate dai parametri HDPE.
    cp_SEBS_loc = f_SEBS * cp_HDPE
    Dh_SEBS_loc = f_SEBS * Dh_HDPE

    cp_solid   = (1.d0 - y_SEBS)*cp_wax + y_SEBS*cp_SEBS_loc
    Dh_melting = (1.d0 - y_SEBS)*Dh_wax + y_SEBS*Dh_SEBS_loc

    ! bilancio energetico
    Dh_conduction = cp_solid*(Tw - T_i)        ! calore di conduzione nel solido (J/kg)
    heat_absorbed = Dh_melting + Dh_conduction

    ! flusso di massa risultante dal bilancio energetico alla superficie del grano
    mdot = qw / heat_absorbed                  ! (kg/(s m^2))

    mdot_wax = (1.d0 - y_SEBS)*mdot

    !% Termine sorgente delle specie
    omega = 0.d0
    if (cracking) then
      ! Cracking della paraffina in etilene e idrogeno a parete: C32H66 → 16 C2H4 + H2
      omega(6)  = 0.0044d0*mdot_wax                  ! H2
      omega(2)  = 0.9956d0*mdot_wax + y_SEBS*mdot    ! C2H4
    else
      ! iniezione di paraffina come C32H66 e di etilene per il SEBS
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

    omegadot = 0d0  ! inizializzazione

    ! Costanti
    RuTw = Runiv*Tw
    inv_sqrtTw = 1.0d0/sqrt(Tw)

    ! Pressioni parziali in atmosfere
    pi = roi*Tw*Ri_tab/pref

    ! Costanti cinetiche
    kH2O = 4.80d5  *exp(-288.00d6/RuTw)  ! H2O
    kCO2 = 9.00d3  *exp(-285.00d6/RuTw)  ! CO2
    kOH  = 3.61d2  *inv_sqrtTw
    kO   = 6.655d2 *inv_sqrtTw

    k5   = 2.40d3  *exp(-125.60d6/RuTw)  ! O2a
    k6   = 2.13d1  *exp(  17.17d6/RuTw)  ! O2b
    k7   = 5.35d-1 *exp( -63.64d6/RuTw)  ! O2c
    k8   = 1.81d7  *exp(-406.10d6/RuTw)  ! O2d

    ! Calcolo del termine Y_term
    Y_term = 1d0 / ( 1d0 + k8/(k7*pi(iO2)) )

    ! Tasso di steam-gasification (C + H2O -> CO + H2): usa deliberatamente la
    ! pressione parziale dell'H2O, indice iH2O, coerente con la tabella locale
    ! delle specie (iH2O=7, iO2=8). Non sostituire l'indice.
    dotm_H2O = kH2O * sqrt(pi(iH2O))
    dotm_CO2 = kCO2 * sqrt(pi(iCO2))
    dotm_O2  = Y_term * (k5*pi(iO2))/(1.0d0 + k6*pi(iO2)) + k7*pi(iO2)*(1d0-Y_term)
    dotm_OH  = kOH*pi(iOH)
    dotm_O   = kO*pi(iO)
    ! Calcolo dei termini sorgente
    omegadot(iO2)  = -(Wm_tab(iO2)/wmc)*(0.5d0*dotm_O2)                                         ! O2
    omegadot(iH2O) = -(Wm_tab(iH2O)/wmc)*(dotm_H2O)                                             ! H2O
    omegadot(iCO)  =  (Wm_tab(iCO)/wmc)*(2.d0*dotm_CO2 + dotm_H2O + dotm_O2 + dotm_OH + dotm_O) ! CO
    omegadot(iCO2) = -(Wm_tab(iCO2)/wmc)*(dotm_CO2)                                             ! CO2
    omegadot(iH2)  =  (Wm_tab(iH2)/wmc)*(dotm_H2O)                                              ! H2
    omegadot(iO)   = -(Wm_tab(iO)/wmc)*(dotm_O)                                                 ! O
    omegadot(iH)   =  (Wm_tab(iH)/wmc)*(dotm_OH)                                                ! H
    omegadot(iOH)  = -(Wm_tab(iOH)/wmc)*(dotm_OH)                                               ! OH

    ! flusso di massa totale in uscita
    mdot = sum ( omegadot )

  end subroutine GSI_Carbon


  subroutine GSI_HDPE ( Tw, omegadot, mdot )
    use FLINT_Lib_Thermodynamic, only: Runiv
    implicit none
    real(R8), intent(in)    :: Tw
    real(R8), intent(inout) :: omegadot(nsc), mdot

    omegadot = 0d0
    mdot = A_HDPE * exp(-Ea_HDPE/(2*Runiv*Tw))
    omegadot(2) = mdot ! etilene, C2H4

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
    omegadot(2) = mdot ! propilene, C3H6

  end subroutine GSI_PP


end module Q2D_Lib_GSI
