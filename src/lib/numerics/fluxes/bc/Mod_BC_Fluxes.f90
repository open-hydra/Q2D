module Q2D_Mod_BC_Fluxes
  ! BC-type → subroutine dispatcher on the BC codes:
  !   101, 103       block-to-block connection / multi-solver coupling
  !   102            chimera
  !   200            sector rotational periodic (not supported: error stop)
  !   201            in-plane translational/rotational periodic (via connection)
  !   300            symmetry
  !   301, 302       wall at imposed heat flux / imposed temperature
  !   303, 304       wall melting / pyrolysis
  !   400            extrapolation
  !   401..405       inflow (T0+p0 / mdot+T0 / mdot+T / supersonic)
  !   406            outflow
  !   407            farfield
  !   410            not supported: IO_BC.f90 rejects it with error stop
  !                  (no case in the dispatcher)
  !   420            choking nozzle — Inflow_Nozzle
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
    ! Reads the obj_* singletons only once, declares locally
    ! Sc/Sct/Prt/soot_enabled + SD_limiter/SD_riemann and propagates them to
    ! BC_Connection_Eul/_Visc.
    use Q2D_Config_Types_m, only: obj_rans, obj_soot, obj_space_scheme, obj_riemann, obj_rot
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Locals
    integer  :: f, lower, upper, i, b
    integer  :: error
    integer  :: Bm, Im, Jm, Fm, Bs, Fs
    real(R8) :: T0, g, pstat, BC_mdot
    real(R8) :: Sc, Sct, Prt
    logical  :: SD_limiter, SD_riemann, soot_enabled, RSM
    ! Wall velocity in the rotating frame: buffer for the iteration in cases 301/302.
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

    ! The BC fluxes are computed in order of block and face type.
    upper = 0

    blocks: do b = 1, domain % nb
      faces: do f = 1, 4

        lower = upper + 1                   ! Update the lower bound
        upper = upper + domain % n_bf(b,f)  ! Upper bound: add the number of cells on face f of block b

        !$omp do schedule(dynamic) private(i, Bm, Im, Jm, Fm, error, Bs, Fs, T0, g, pstat, BC_mdot, r_fc, w_wall_rf)
        do i = lower, upper
          Bm = domain % bc(i) % b
          if (.not. is_local_block(Bm)) cycle  ! MPI: skip non-local blocks
          Im = domain % bc(i) % i
          Jm = domain % bc(i) % j
          Fm = domain % bc(i) % f

          select case ( domain % bc(i) % type )

            case (103) ! multi-solver coupling
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              domain % blk(Bm) % R(:,Im,Jm) = domain % blk(Bm) % R(:,Im,Jm) + domain % bc(i) % ext_flux

            case (101, 102, 201) ! connection (101), chimera (102), translational/rotational periodic (201)
              ! The Connection_Eul/Visc handler covers block-to-block
              ! connection (101), chimera (102) and in-plane translational/
              ! rotational periodicity (201). BC 201 is routed via
              ! Fill_BC_Ghost_Connection + BC_Connection_Eul: a mechanism
              ! distinct from the sector rotational periodicity of case 200
              ! (which would require out-of-plane faces and is intrinsically 3D,
              ! hence not supported). BC 410 (mapped-state) is not
              ! supported: IO_BC.f90 rejects it with error stop and no
              ! case exists in the dispatcher.
              call BC_Connection_Eul ( Im, Jm, Fm, domain % blk(Bm), SD_limiter, SD_riemann )
              if (model > 0) &
                call BC_Connection_Visc ( Im, Jm, Fm, domain % blk(Bm),       &
                                          domain % bc(i) % Mg(1), domain % bc(i) % Pg, &
                                          Sc, Sct, Prt, soot_enabled )

            case (200) ! sector rotational periodic (not supported: error stop)
              call BC_Rotational_Periodic_Eul ( Im, Jm, Fm, domain % blk(Bm) )

            case (300) ! symmetry
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) call BC_Symmetry_Visc ( Im, Jm, Fm, domain % blk(Bm), RSM )

            case (301) ! wall: imposed heat flux
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) then
                ! Stationary patches (casing/hub): fixed in the laboratory →
                ! relative apparent velocity w_wall = -(Ω × r) ≠ 0, passed to the BC.
                ! Non-stationary patches (blades): rotate with the frame →
                ! w_wall = 0, call without the argument.
                if ( rf_enabled .and. RF_Is_Stationary_Face(Bm, Fm) ) then
                  call RF_Face_Center( domain % blk(Bm) % node, Im, Jm, Fm, r_fc )
                  w_wall_rf = RF_Wall_Velocity( r_fc )
                  call BC_Wall_Heat ( Im, Jm, Fm, domain % blk(Bm), domain % bc(i) % qw, &
                                      w_wall = w_wall_rf )
                else
                  call BC_Wall_Heat ( Im, Jm, Fm, domain % blk(Bm), domain % bc(i) % qw )
                end if
              end if

            case (302) ! wall: imposed temperature
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) then
                ! Same scheme as case 301.
                if ( rf_enabled .and. RF_Is_Stationary_Face(Bm, Fm) ) then
                  call RF_Face_Center( domain % blk(Bm) % node, Im, Jm, Fm, r_fc )
                  w_wall_rf = RF_Wall_Velocity( r_fc )
                  call BC_Wall_Temperature ( Im, Jm, Fm, domain % blk(Bm), domain % bc(i) % Tw, &
                                             w_wall = w_wall_rf )
                else
                  call BC_Wall_Temperature ( Im, Jm, Fm, domain % blk(Bm), domain % bc(i) % Tw )
                end if
              end if

            case (303) ! wall: T_wall + radiative flux (melting)
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) &
                call BC_Wall_Melting ( Im, Jm, Fm, domain % blk(Bm), &
                                       domain % bc(i) % Tw, domain % bc(i) % qrad )

            case (304) ! wall: radiative flux (pyrolysis/ablation)
              call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )
              if (model > 0) &
                call BC_Wall_Pyrolysis ( Im, Jm, Fm, domain % blk(Bm), &
                                         domain % bc(i) % qrad, domain % bc(i) % type )

            case (400) ! extrapolation
              call BC_Extrapolation ( Im, Jm, Fm, domain % blk(Bm) )

            case (401) ! inflow: stagnation conditions T0 + p0
              call BC_Inlet_StagnCond ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                         domain % bc(i) % T0, domain % bc(i) % p0, &
                                         domain % bc(i) % rel_fac, &
                                         domain % bc(i) % alpha, &
                                         domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                         domain % bc(i) % ci(1:nrans), &
                                         domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (402) ! inflow: T0 + time-varying p0 (subsonic stagnation)
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

            case (403) ! inflow: total temperature T0 + imposed mass flux
              ! `domain%bc(i)%mdot` is the inflow mass-flow specification and must not
              ! be passed also as an output argument: aliasing would
              ! overwrite it, degrading it down to divisions by zero.
              call BC_Inlet_MassFlux_T0 ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                           domain % bc(i) % T0, domain % bc(i) % mdot, &
                                           domain % bc(i) % rel_fac, &
                                           domain % bc(i) % alpha, &
                                           domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                           domain % bc(i) % ci(1:nrans), error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (404) ! inflow: static temperature T + imposed mass flux
              ! As for case 403: `domain%bc(i)%mdot` is input only, never output.
              call BC_Inlet_MassFlux_T ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                          domain % bc(i) % T0, domain % bc(i) % mdot, &
                                          domain % bc(i) % rel_fac, &
                                          domain % bc(i) % alpha, &
                                          domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                          domain % bc(i) % ci(1:nrans), error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (405) ! supersonic inflow: Mach + static temperature + static pressure
              call BC_Inlet_Supersonic_Static ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                                 domain % bc(i) % mach, &
                                                 domain % bc(i) % T0, domain % bc(i) % pamb, &
                                                 domain % bc(i) % rel_fac, &
                                                 domain % bc(i) % alpha, &
                                                 domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                                 domain % bc(i) % ci(1:nrans), &
                                                 domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (406) ! outflow: imposed back-pressure (or extrapolation if pAmb ≤ 0)
              call BC_Outflow ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                 domain % bc(i) % pAmb, domain % bc(i) % rel_fac, &
                                 domain % bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Bm, Im, Jm, Fm, domain % blk(Bm) )

            case (407) ! ambient: T0 + p0 + back-pressure
              call BC_FarField ( Bm, Im, Jm, Fm, domain % blk(Bm), &
                                 domain % bc(i) % T0, domain % bc(i) % p0, &
                                 domain % bc(i) % pamb, &
                                 domain % bc(i) % rel_fac, &
                                 domain % bc(i) % alpha, domain % bc(i) % beta, &
                                 domain % bc(i) % ci(nrans+1:nrans+nsc), &
                                 domain % bc(i) % ci(1:nrans), &
                                 BC_mdot )

            case (420) ! choking nozzle (BC_Mach==-20 → Borda)
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

            case (501) ! mean-flow manifold connection
              Bs = domain % bc(i) % bs
              Fs = domain % bc(i) % fs
              call BC_Manifold ( domain % blk(Bm), Fm, domain % blk(Bs), Fs, T0, g, pstat )

            case (502) ! SRM grain combustion (bc_type extended with
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
