module Q2D_Lib_BC_Fluxes_Wall_Melting
  ! BC 303 wall_melting: imposed T_wall + radiative flux qrad; GSI for the
  ! ablation/melting mass flux.
  !
  ! 2D formulation:
  !   * ghost-cell indices Ig/Jg; Normal(2), M(2,2), Gradient(nprim,2),
  !     Stress(2), Blowing(2);
  !   * int_l: integrated face length;
  !   * rank-3 residual Blk%r(:,Im,Jm).
  !
  ! V_dot_G: Stress_Vector is invoked without the optional argument V_dot_G
  ! (default 0); for the wall-with-blowing pattern V_dot_G = 0 is consistent
  ! with the absence of local area variation at the wall (boundary cell).
  !
  ! GSI: Wall_Melting invokes GSI_Melting. The extended variant GSI_Paraffin
  ! (mixture with SEBS + optional cracking branch) remains exposed in
  ! Q2D_Lib_GSI but is NOT linked here: its distinctive branches are inactive
  ! with the defaults (y_SEBS=0, cracking=.false.) and cannot be activated in
  ! the absence of an INI parsing of the GSI parameters, not yet implemented.
  ! With the defaults GSI_Paraffin coincides bit-for-bit with GSI_Melting, so
  ! linking the extended variant would give no exercisable benefit.
  !
  ! Compute_Wall_Properties: exit_array dimension(7).
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index, Compute_Wall_Properties
  use Q2D_Lib_Fluid,              only: Stress_Vector
  use Q2D_Lib_GSI,                only: GSI_Melting
  use Q2D_Lib_RANS,               only: RANS_Set_Wall_Values, RANS_Diffusive_Flux, RANS_Extrapolate_Wall
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, co_k_mi_lam_Wilke, f_tabT, h_tab

  implicit none
  private

  public :: BC_Wall_Melting

contains


  subroutine BC_Wall_Melting ( Im, Jm, Fm, Blk, T_wall, q_rad, Ovar )
    implicit none
    integer,             intent(in)    :: Im, Jm, Fm
    real(R8),            intent(in)    :: T_wall, q_rad
    type(Q2D_block_type),intent(inout) :: Blk
    real(R8), optional, dimension(7), intent(inout) :: Ovar
    ! Locals
    integer  :: s, modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Ig, Jg
    real(R8) :: Normal(2), int_l, Dist, dl, M(2,2), Prim(nprim), rho, Rgas, T, Gradient(nprim,2)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, q_conv, omega(nsc), mdot
    real(R8) :: Stress(2), Blowing(2), Sum_Omegai_Hwi, Flux(nprim)
    real(R8) :: r_ave


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    ! Face metrics
    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l
    Dist   = 1d-20
    M      = Blk % M (Im,Jm) % c

    ! Boundary-cell variables
    Prim = Blk % P(:,Im,Jm)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    ! Initialization
    Gradient = 0d0

    Gradient(np,Dir) = ( T - T_wall ) * modfm3                  ! Temperature
    Gradient(np,:)   = matmul ( Gradient(np,:), M )             ! Temperature gradient
    rho_wall         = Prim(np) / ( Rgas * T_wall )             ! Boundary layer: p(wall)=p(cell), R(w)=R(c), ci(w)=ci(c)
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho             ! Partial densities at the wall

    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl ) ! transport properties

    if (model==2) then
      dl = blk % yn(im,jm)
      call RANS_Set_Wall_Values( mil, Prim_Wall(nt:nprim), dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
    end if

    q_conv = kl * dot_product( Gradient(np,:), Normal )         ! Convective heat flux

    call GSI_Melting ( q_conv + q_rad, T_wall, omega, mdot )    ! Gas-surface interaction: paraffin mass flux

    Blowing = mdot / rho_wall * Normal                          ! Blowing velocity
    Gradient(nu:nv,Dir) = ( Prim(nu:nv) - Blowing ) * modfm3    ! Velocity gradient due to blowing

    Gradient(nu:nv,:) = matmul ( Gradient(nu:nv,:), M )         ! Velocity-gradient transformation

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
    Stress = Stress_Vector ( Gradient(nu:nv,:), Normal, mil, 0d0, Prim(nt:), r_ave, 0.d0 )

    Sum_Omegai_Hwi = 0.d0
    do s = 1, nsc
      Sum_Omegai_Hwi = Sum_Omegai_Hwi + omega(s) * f_tabT( T_wall, s, h_tab )
    end do

    ! Fluxes
    Flux = 0.0d0
    Flux(1:nsc) = - omega * int_l
    Flux(nu:nv) = ( Stress - mdot * Blowing ) * int_l
    Flux(np)    = int_l * ( q_conv + q_rad - Sum_Omegai_Hwi - mdot * 0.5d0 * sum( Blowing**2 ) + &
                             dot_product( Stress, Blowing ) )

    if (model==2) then
      Gradient(nt:nprim,:) = matmul ( Gradient(nt:nprim,:), M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim),                  &
                                 rans_variables=Prim_Wall(nt:nprim),   &
                                 vel_gradient=Gradient(nu:nv,:),       &
                                 rans_gradient=Gradient(nt:nprim,:),   &
                                 mul=mil, rho=rho_wall,                &
                                 int_l=int_l, normal=Normal,           &
                                 vr=0.d0, r_ave=r_ave, dist=dist )
      ! Ghost-cell extrapolation of the RANS variables
      Ig = Im - guide(Fm,1)
      Jg = Jm - guide(Fm,2)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm) = Blk % r (:,Im,Jm) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties ( stress=Stress, pw=prim(np), qw=q_conv+q_rad,    &
                                                       mdot=mdot, y=blk%dl(im,jm)%c(dir)*0.5d0,        &
                                                       Tw=T_wall, rhow=rho_wall, mu=mil,               &
                                                       exit_array=Ovar )

  end subroutine BC_Wall_Melting


end module Q2D_Lib_BC_Fluxes_Wall_Melting
