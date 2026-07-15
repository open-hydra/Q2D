module Q2D_Mod_BC_Fluxes
  ! Dispatcher tipo BC → subroutine sui codici BC:
  !   101, 103       connessione tra blocchi / accoppiamento multi-solutore
  !   102            chimera
  !   200            periodica rotazionale a settore (non supportata: error stop)
  !   201            periodica traslazionale/rotazionale in piano (via connessione)
  !   300            simmetria
  !   301, 302       parete a flusso termico / a temperatura imposti
  !   303, 304       parete melting / pirolisi
  !   400            estrapolazione
  !   401..405       inflow (T0+p0 / mdot+T0 / mdot+T / supersonico)
  !   406            outflow
  !   407            farfield
  !   410            non supportata: IO_BC.f90 la rifiuta con error stop
  !                  (nessun case nel dispatcher)
  !   420            ugello in choking — Inflow_Nozzle
  !   501            manifold
  !   502            SRM
  use iso_fortran_env,                 only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,            only: Q2D_domain_type
  use Q2D_Global_m,                    only: nsc, nrans, model, nu, nv, np
  use Q2D_Mod_MPI,                     only: is_local_block
  use Q2D_Lib_BC_Fluxes_Connection,    only: BC_Connection_Eul, BC_Connection_Visc
  use Q2D_Lib_BC_Fluxes_Symmetry,      only: BC_Symmetry_Eul, BC_Symmetry_Visc
  use Q2D_Lib_BC_Fluxes_Extrapolation, only: BC_Extrapolation
  use Q2D_Lib_BC_Fluxes_Wall_Heat,     only: BC_Wall_Heat
  use Q2D_Lib_BC_Fluxes_Wall_Temperature, only: BC_Wall_Temperature
  use Q2D_Lib_BC_Fluxes_Wall_Melting,  only: BC_Wall_Melting
  use Q2D_Lib_BC_Fluxes_Wall_Pyrolysis,only: BC_Wall_Pyrolysis
  use Q2D_Lib_BC_Fluxes_Inflow,        only: BC_Inlet_StagnCond, BC_Inlet_MassFlux_T0, &
                                              BC_Inlet_MassFlux_T, BC_Inlet_Supersonic_Static
  use Q2D_Lib_BC_Fluxes_Outflow,       only: BC_Outflow
  use Q2D_Lib_BC_Fluxes_Inflow_Nozzle, only: BC_Inflow_Nozzle
  use Q2D_Lib_BC_Fluxes_FarField,      only: BC_FarField
  use Q2D_Lib_BC_Fluxes_Manifold,      only: BC_Manifold
  use Q2D_Lib_BC_Fluxes_Rotational,    only: BC_Rotational_Periodic_Eul
  use Q2D_Lib_BC_Fluxes_SRM,           only: BC_SRM
  use Q2D_Lib_RotatingFrame,           only: RF_Wall_Velocity, RF_Face_Center, &
                                              RF_Is_Stationary_Face

  implicit none
  private

  public :: BC_Fluxes

contains


  subroutine BC_Fluxes ( domain )
    ! Legge una sola volta i singleton obj_*, dichiara localmente
    ! Sc/Sct/Prt/soot_enabled + SD_limiter/SD_riemann e li propaga a
    ! BC_Connection_Eul/_Visc.
    use Q2D_Config_Types_m, only: obj_rans, obj_soot, obj_space_scheme, obj_riemann, obj_rot
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Locali
    integer  :: f, lower, upper, i, b
    integer  :: error
    integer  :: Bm, Im, Jm, Fm, Bs, Fs
    real(R8) :: T0, g, pstat, BC_mdot
    real(R8) :: Sc, Sct, Prt
    logical  :: SD_limiter, SD_riemann, soot_enabled, RSM
    ! Velocità di parete nel sistema rotante: buffer per iterazione nei case 301/302.
    real(R8) :: r_fc(2), w_wall_rf(2)
    logical  :: rf_enabled

    RSM          = obj_rans%RSM
    SD_limiter   = obj_space_scheme%SD
    SD_riemann   = obj_riemann%SD
    soot_enabled = obj_soot%enabled
    rf_enabled   = obj_rot%enabled
    Sc  = obj_rans%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt

    ! I flussi BC sono calcolati in ordine di blocco e di tipo di faccia.
    upper = 0

    blocks: do b = 1, domain % nb
      faces: do f = 1, 4

        lower = upper + 1                   ! Aggiorna il limite inferiore
        upper = upper + domain % n_bf(b,f)  ! Limite superiore: somma il numero di celle sulla faccia f del blocco b

        !$omp do schedule(dynamic) private(i, Bm, Im, Jm, Fm, error, Bs, Fs, T0, g, pstat, BC_mdot, r_fc, w_wall_rf)
        do i = lower, upper
          Bm = domain % bc(i) % b
          if (.not. is_local_block(Bm)) cycle  ! MPI: salta i blocchi non locali
          Im = domain % bc(i) % i
          Jm = domain % bc(i) % j
          Fm = domain % bc(i) % f

          select case ( domain % bc(i) % type )

            case (103) ! accoppiamento multi-solutore
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              domain % blk(Bm) % R(:,Im,Jm) = domain % blk(Bm) % R(:,Im,Jm) + domain % bc(i) % ext_flux

            case (101, 102, 201) ! connessione (101), chimera (102), periodica traslazionale/rotazionale (201)
              ! L'handler Connection_Eul/Visc copre la connessione tra blocchi
              ! (101), la chimera (102) e la periodicità traslazionale/
              ! rotazionale in piano (201). La BC 201 è smistata via
              ! Fill_BC_Ghost_Connection + BC_Connection_Eul: meccanismo
              ! distinto dalla periodicità rotazionale a settore del case 200
              ! (che richiederebbe facce fuori piano ed è intrinsecamente 3D,
              ! quindi non supportata). La BC 410 (mapped-state) non è
              ! supportata: IO_BC.f90 la rifiuta con error stop e non esiste
              ! alcun case nel dispatcher.
              call BC_Connection_Eul ( Im, Jm, Fm, domain % blk(Bm), SD_limiter, SD_riemann )
              if (model > 0) &
                call BC_Connection_Visc ( Im, Jm, Fm, domain % blk(Bm),       &
                                          domain % bc(i) % Mg(1), domain % bc(i) % Pg, &
                                          Sc, Sct, Prt, soot_enabled )

            case (200) ! periodica rotazionale a settore (non supportata: error stop)
              call BC_Rotational_Periodic_Eul ( Im, Jm, Fm, domain % blk(Bm) )

            case (300) ! simmetria
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) call BC_Symmetry_Visc ( Im, Jm, Fm, domain % blk(Bm), RSM )

            case (301) ! parete: flusso termico imposto
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) then
                ! Patch stazionarie (casing/hub): fisse nel laboratorio →
                ! velocità apparente relativa w_wall = -(Ω × r) ≠ 0, passata alla BC.
                ! Patch non stazionarie (pale): ruotano col riferimento →
                ! w_wall = 0, chiamata senza argomento.
                if ( rf_enabled .and. RF_Is_Stationary_Face(Bm, Fm) ) then
                  call RF_Face_Center( domain % blk(Bm) % node, Im, Jm, Fm, r_fc )
                  w_wall_rf = RF_Wall_Velocity( r_fc )
                  call BC_Wall_Heat ( Im, Jm, Fm, domain % blk(Bm), domain % bc(i) % qw, &
                                      w_wall = w_wall_rf )
                else
                  call BC_Wall_Heat ( Im, Jm, Fm, domain % blk(Bm), domain % bc(i) % qw )
                end if
              end if

            case (302) ! parete: temperatura imposta
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) then
                ! Stesso schema del case 301.
                if ( rf_enabled .and. RF_Is_Stationary_Face(Bm, Fm) ) then
                  call RF_Face_Center( domain % blk(Bm) % node, Im, Jm, Fm, r_fc )
                  w_wall_rf = RF_Wall_Velocity( r_fc )
                  call BC_Wall_Temperature ( Im, Jm, Fm, domain % blk(Bm), domain % bc(i) % Tw, &
                                             w_wall = w_wall_rf )
                else
                  call BC_Wall_Temperature ( Im, Jm, Fm, domain % blk(Bm), domain % bc(i) % Tw )
                end if
              end if

            case (303) ! parete: T_wall + flusso radiativo (melting)
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) &
                call BC_Wall_Melting ( Im, Jm, Fm, domain % blk(Bm), &
                                       domain % bc(i) % Tw, domain % bc(i) % qrad )

            case (304) ! parete: flusso radiativo (pirolisi/ablazione)
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) &
                call BC_Wall_Pyrolysis ( Im, Jm, Fm, domain % blk(Bm), &
                                         domain % bc(i) % qrad, domain % bc(i) % type )

            case (400) ! estrapolazione
              call BC_Extrapolation ( Im, Jm, Fm, domain % blk(Bm) )

            case (401) ! ingresso: condizioni di ristagno T0 + p0
              call BC_Inlet_StagnCond ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                         domain % bc(i) % T0, domain % bc(i) % p0, &
                                         domain % bc(i) % rel_fac, &
                                         domain % bc(i) % alpha, &
                                         domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                         domain % bc(i) % ci(1:nrans), &
                                         domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (402) ! ingresso: T0 + p0 variabile nel tempo (ristagno subsonico)
              if (domain % bc(i) % p0time % exists) &
                domain % bc(i) % p0 = domain % bc(i) % p0time % update(domain % time)
              call BC_Inlet_StagnCond ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                         domain % bc(i) % T0, domain % bc(i) % p0, &
                                         domain % bc(i) % rel_fac, &
                                         domain % bc(i) % alpha, &
                                         domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                         domain % bc(i) % ci(1:nrans), &
                                         domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (403) ! ingresso: temperatura totale T0 + flusso di massa imposto
              ! `domain%bc(i)%mdot` è la specifica di portata in ingresso e non va
              ! passata anche come argomento di output: l'aliasing la
              ! sovrascriverebbe, degradandola fino a divisioni per zero.
              call BC_Inlet_MassFlux_T0 ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                           domain % bc(i) % T0, domain % bc(i) % mdot, &
                                           domain % bc(i) % rel_fac, &
                                           domain % bc(i) % alpha, &
                                           domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                           domain % bc(i) % ci(1:nrans), error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (404) ! ingresso: temperatura statica T + flusso di massa imposto
              ! Come per il case 403: `domain%bc(i)%mdot` è solo input, mai output.
              call BC_Inlet_MassFlux_T ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                          domain % bc(i) % T0, domain % bc(i) % mdot, &
                                          domain % bc(i) % rel_fac, &
                                          domain % bc(i) % alpha, &
                                          domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                          domain % bc(i) % ci(1:nrans), error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (405) ! ingresso supersonico: Mach + temperatura statica + pressione statica
              call BC_Inlet_Supersonic_Static ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                                 domain % bc(i) % mach, &
                                                 domain % bc(i) % T0, domain % bc(i) % pamb, &
                                                 domain % bc(i) % rel_fac, &
                                                 domain % bc(i) % alpha, &
                                                 domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                                 domain % bc(i) % ci(1:nrans), &
                                                 domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (406) ! outflow: back-pressure imposta (o estrapolazione se pAmb ≤ 0)
              call BC_Outflow ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                 domain % bc(i) % pAmb, domain % bc(i) % rel_fac, &
                                 domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (407) ! ambiente: T0 + p0 + back-pressure
              call BC_FarField ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                 domain % bc(i) % T0, domain % bc(i) % p0, &
                                 domain % bc(i) % pamb, &
                                 domain % bc(i) % rel_fac, &
                                 domain % bc(i) % alpha, domain % bc(i) % beta, &
                                 domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                 domain % bc(i) % ci(1:nrans), &
                                 BC_mdot )

            case (420) ! ugello in choking (BC_Mach==-20 → Borda)
              call BC_Inflow_Nozzle ( Im, Jm, Fm, domain % blk(Bm), &
                                      domain % bc(i) % mach, &
                                      domain % bc(i) % T0, domain % bc(i) % p0, &
                                      domain % bc(i) % rel_fac, &
                                      domain % bc(i) % alpha, domain % bc(i) % beta, &
                                      domain % bc(i) % pAmb, &
                                      domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                      domain % bc(i) % ci(1:nrans), &
                                      error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (501) ! connessione manifold a flusso medio
              Bs = domain % bc(i) % bs
              Fs = domain % bc(i) % fs
              call BC_Manifold ( domain % blk(Bm), Fm, domain % blk(Bs), Fs, T0, g, pstat )

            case (502) ! combustione del grain SRM (bc_type esteso con
                        ! aCoeff/n/pRef/Taf/haf).
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              call BC_SRM ( Im, Jm, Fm, domain % blk(Bm),                       &
                            domain % bc(i) % aCoeff, domain % bc(i) % n,         &
                            domain % bc(i) % pRef,   domain % bc(i) % Taf,       &
                            domain % bc(i) % haf,                                &
                            domain % bc(i) % ci(1:nsc) )

          end select
        enddo

      enddo faces
    enddo blocks

  end subroutine BC_Fluxes


end module Q2D_Mod_BC_Fluxes
