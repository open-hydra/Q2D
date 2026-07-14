module Q2D_Lib_BC_Fluxes_Wall_Pyrolysis
  ! BC 304 wall_pyrolysis with a
  ! 20-iter Newton-Raphson solve of T_wall (energy balance q_conv = q_pyro +
  ! q_rad) + material-specific GSI dispatch on the `fuel` code.
  !
  ! Fuel codes:
  !   * 9, 304: HDPE / generic ablation (304 = ATLAS qrad-wall type)
  !   * 10:     PP
  !   * 12:     Carbon
  !   * 13:     HTPB
  !
  ! Rigid 2D reduction (No-w principle):
  !   * drop Km argument, Face_k, Ig/Jg/Kg → Ig/Jg.
  !   * Normal(2), M(2,2), Gradient(nprim,2), Stress(2), Blowing(2).
  !   * Area→int_l rename.
  !   * rank-3 residual.
  !
  ! V_dot_G: Stress_Vector called with the 7-arg form (V_dot_G optional,
  ! default 0); consistent with a boundary cell without local area-variation.
  !
  ! GSI_Paraffin: NOT wired in Wall_Pyrolysis (physically incompatible pattern —
  ! Wall_Pyrolysis iterates T_wall while GSI_Paraffin/_Melting compute mdot from
  ! qw at fixed T_wall). GSI_Paraffin (Q2D-extra peculiarità #8) stays exposed in
  ! Q2D_Lib_GSI but is currently wired nowhere; BC 303 Wall_Melting calls
  ! GSI_Melting (MOSE-mirror). See Lib_BC_Fluxes_Wall_Melting.f90 header.
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index, Compute_Wall_Properties
  use Q2D_Lib_Fluid,              only: Stress_Vector
  use Q2D_Lib_GSI,                only: GSI_HDPE, GSI_PP, GSI_Carbon, GSI_HTPB,         &
                                        cp_HDPE, Dh_HDPE, Ti_HDPE,                       &
                                        cp_PP,   Dh_PP,   Ti_PP,                         &
                                        cp_HTPB, Dh_HTPB, Ti_HTPB
  use Q2D_Lib_RANS,               only: RANS_Set_Blowing_Wall, RANS_Diffusive_Flux, RANS_Extrapolate_Wall
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, co_k_mi_lam_Wilke, f_tabT, h_tab

  implicit none
  private

  public :: BC_Wall_Pyrolysis

contains


  subroutine BC_Wall_Pyrolysis ( Im, Jm, Fm, Blk, qrad, fuel, Ovar )
    implicit none
    integer,  intent(in) :: Im, Jm, Fm, fuel
    real(R8), intent(in) :: qrad
    type(Q2D_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(7), intent(inout) :: Ovar
    ! Local
    integer  :: s, modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Ig, Jg, iter
    real(R8) :: Normal(2), int_l, Dist, dl, M(2,2), Prim(nprim), rho, Rgas, T, Gradient(nprim,2)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, omega(nsc), mdot
    real(R8) :: Stress(2), Blowing(2), Sum_Omegai_Hwi, Flux(nprim)
    real(R8) :: T_wall, T_wall_old, balance, balanceold, q_conv, q_pyro, q_rad
    real(R8) :: r_ave


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l
    M      = Blk % M (Im,Jm) % c

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! r_ave for the Stress_Vector wiring (Q2D 2D-axi/area-var aware)
    select case(Fm)
      case(1)
        r_ave = 0.5d0 * ( blk % node(im-1,jm)%c(2) + blk % node(im-1,jm-1)%c(2) )
      case(2)
        r_ave = 0.5d0 * ( blk % node(im,jm)%c(2)   + blk % node(im,jm-1)%c(2)   )
      case(3)
        r_ave = 0.5d0 * ( blk % node(im,jm-1)%c(2) + blk % node(im-1,jm-1)%c(2) )
      case(4)
        r_ave = 0.5d0 * ( blk % node(im,jm)%c(2)   + blk % node(im-1,jm)%c(2)   )
    end select

    ! Initialization
    Gradient = 0d0

    ! Newton-Raphson loop initialization
    T_wall = T                                                   ! First guess: T_wall = T_cell
    rho_wall = Prim(np) / ( Rgas * T_wall )                      ! BL: p(wall)=p(cell), R(w)=R(c), ci(w)=ci(c)
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho              ! Partial densities at wall

    select case ( fuel )
      case ( 9, 304 ) ! HDPE / generic ablation (304 is the ATLAS qrad-wall type)
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 ) ! PP
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 ) ! Carbon
        call GSI_Carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 ) ! HTPB
        call GSI_HTPB(T_wall, omega, mdot)
    end select
    mdot  = mdot  * (-modfm2) ! adjust sign
    omega = omega * (-modfm2) ! adjust sign

    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl ) ! transport

    q_conv = 0d0 ! convective heat flux (initial guess)
    select case ( fuel ) ! pyrolysis heat flux
    case ( 9, 304 )
      q_pyro = mdot * ( Dh_HDPE + cp_HDPE * (T_wall - Ti_HDPE) )
    case ( 10 )
      q_pyro = mdot * ( Dh_PP + cp_PP * (T_wall - Ti_PP) )
    case ( 12 )
      q_pyro = 0d0
      do s = 1, nsc
        q_pyro = q_pyro + omega(s) * f_tabT( T_wall, s, h_tab )
      end do
    case ( 13 )
      q_pyro = mdot * ( Dh_HTPB + cp_HTPB * (T_wall - Ti_HTPB) )
    end select
    q_rad = qrad * modfm2 ! adjust sign for radiation
    balance = q_conv - q_pyro - q_rad

    do iter = 1, 20 ! Newton-Raphson loop
      ! first part of the iteration
      balanceold = balance
      T_wall_old = T_wall
      T_wall = T_wall * (1d0 + 1d-4) ! increment
      Gradient = 0d0
      Gradient(np,Dir) = ( T - T_wall ) * modfm3
      Gradient(np,:)   = matmul ( Gradient(np,:), M )
      rho_wall         = Prim(np) / ( Rgas * T_wall )
      Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho
      select case ( fuel )
      case ( 9, 304 )
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 )
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 )
        call GSI_Carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 )
        call GSI_HTPB(T_wall, omega, mdot)
      end select
      mdot  = mdot  * (-modfm2)
      omega = omega * (-modfm2)
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )
      q_conv = kl * dot_product( Gradient(np,:), Normal )
      select case ( fuel )
      case ( 9, 304 )
        q_pyro = mdot * ( Dh_HDPE + cp_HDPE * (T_wall - Ti_HDPE) )
      case ( 10 )
        q_pyro = mdot * ( Dh_PP + cp_PP * (T_wall - Ti_PP) )
      case ( 12 )
        q_pyro = 0d0
        do s = 1, nsc
          q_pyro = q_pyro + omega(s) * f_tabT( T_wall, s, h_tab )
        end do
      case ( 13 )
        q_pyro = mdot * ( Dh_HTPB + cp_HTPB * (T_wall - Ti_HTPB) )
      end select
      balance = q_conv - q_pyro - q_rad

      ! new wall temperature (Newton step)
      T_wall = T_wall_old - balanceold * (1d-4 * T_wall_old) / (balance - balanceold)

      ! second part of the iteration
      Gradient = 0d0
      Gradient(np,Dir) = ( T - T_wall ) * modfm3
      Gradient(np,:)   = matmul ( Gradient(np,:), M )
      rho_wall         = Prim(np) / ( Rgas * T_wall )
      Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho
      select case ( fuel )
      case ( 9, 304 )
        call GSI_HDPE(T_wall, omega, mdot)
      case ( 10 )
        call GSI_PP(T_wall, omega, mdot)
      case ( 12 )
        call GSI_Carbon(Prim_wall(1:nsc), T_wall, omega, mdot)
      case ( 13 )
        call GSI_HTPB(T_wall, omega, mdot)
      end select
      mdot  = mdot  * (-modfm2)
      omega = omega * (-modfm2)
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )
      q_conv = kl * dot_product( Gradient(np,:), Normal )
      select case ( fuel )
      case ( 9, 304 )
        q_pyro = mdot * ( Dh_HDPE + cp_HDPE * (T_wall - Ti_HDPE) )
      case ( 10 )
        q_pyro = mdot * ( Dh_PP + cp_PP * (T_wall - Ti_PP) )
      case ( 12 )
        q_pyro = 0d0
        do s = 1, nsc
          q_pyro = q_pyro + omega(s) * f_tabT( T_wall, s, h_tab )
        end do
      case ( 13 )
        q_pyro = mdot * ( Dh_HTPB + cp_HTPB * (T_wall - Ti_HTPB) )
      end select
      balance = q_conv - q_pyro - q_rad
    end do

    ! From here on T_wall is known
    Gradient = 0d0
    Gradient(np,Dir) = ( T - T_wall ) * modfm3
    Gradient(np,:)   = matmul ( Gradient(np,:), M )
    q_conv           = kl * dot_product( Gradient(np,:), Normal )
    rho_wall         = Prim(np) / ( Rgas * T_wall )
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho

    Blowing             = mdot / rho_wall * Normal               ! Blowing velocity
    Gradient(nu:nv,Dir) = ( Prim(nu:nv) - Blowing ) * modfm3     ! Velocity gradient due to blowing
    Gradient(nu:nv,:)   = matmul ( Gradient(nu:nv,:), M )        ! Velocity gradient transformation
    Stress = Stress_Vector ( Gradient(nu:nv,:), Normal, mil, 0d0, Prim(nt:), r_ave, 0.d0 )

    Sum_Omegai_Hwi = 0d0
    do s = 1, nsc
      Sum_Omegai_Hwi = Sum_Omegai_Hwi + omega(s) * f_tabT( T_wall, s, h_tab )
    end do

    ! Fluxes
    Flux = 0.0d0
    Flux(1:nsc) = - omega * int_l
    Flux(nu:nv) = ( Stress - mdot * Blowing ) * int_l
    Flux(np)    = int_l * ( q_conv - Sum_Omegai_Hwi - mdot * 0.5d0 * sum( Blowing**2 ) + &
                             dot_product( Stress, Blowing ) )

    if (model==2) then
      dl = blk % yn(im,jm)
      call RANS_Set_Blowing_Wall ( rho=rho_wall, mil=mil,                 &
                                    tau=Stress, rans_variables=Prim_Wall(nt:nprim), &
                                    mdot=mdot, dist=dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
      Gradient(nt:nprim,:)   = matmul ( Gradient(nt:nprim,:), M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim),                  &
                                 rans_variables=Prim_Wall(nt:nprim),   &
                                 vel_gradient=Gradient(nu:nv,:),       &
                                 rans_gradient=Gradient(nt:nprim,:),   &
                                 mul=mil, rho=rho_wall,                &
                                 int_l=int_l, normal=Normal,           &
                                 vr=0.d0, r_ave=r_ave, dist=1d-20 )
      ! Ghost-cell extrapolation of RANS variables
      Ig = Im - guide(Fm,1)
      Jg = Jm - guide(Fm,2)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm) = Blk % r (:,Im,Jm) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties ( stress=Stress, pw=prim(np), qw=q_conv,  &
                                                       mdot=mdot, y=blk%dl(im,jm)%c(dir)*0.5d0, &
                                                       Tw=T_wall, rhow=rho_wall, mu=mil,        &
                                                       exit_array=Ovar )

  end subroutine BC_Wall_Pyrolysis


end module Q2D_Lib_BC_Fluxes_Wall_Pyrolysis
