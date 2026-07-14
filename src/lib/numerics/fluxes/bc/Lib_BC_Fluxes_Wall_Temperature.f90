module Q2D_Lib_BC_Fluxes_Wall_Temperature
  ! Twall imposed, mil(wall) computed at Twall.
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index, Compute_Wall_Properties
  use Q2D_Lib_Fluid,              only: Stress_Vector
  use Q2D_Lib_RANS,               only: RANS_Set_Wall_Values, RANS_Diffusive_Flux, RANS_Extrapolate_Wall
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, co_k_mi_lam_Wilke

  implicit none
  private

  public :: BC_Wall_Temperature

contains


  subroutine BC_Wall_Temperature ( Im, Jm, Fm, Blk, T_wall, Ovar, w_wall )
    ! Signature:
    ! w_wall optional dimension(2) — RotatingFrame wall velocity.
    ! Body: present -> Gradient(nu:nv,Dir) = (Prim(nu:nv) - w_wall) * modfm3
    ! (relative to rotating wall); absent -> standard no-slip default.
    ! w_wall(2) vs the 3D peer's w_wall(3).
    implicit none
    integer,             intent(in)    :: Im, Jm, Fm
    real(R8),            intent(in)    :: T_wall
    type(Q2D_block_type),intent(inout) :: Blk
    real(R8), optional, dimension(7), intent(out) :: Ovar
    real(R8), optional, dimension(2), intent(in)  :: w_wall
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, ig, jg
    real(R8) :: Normal(2), int_l, Dist, M(2,2), Prim(nprim), rho, Rgas, T, rho_wall, dl
    real(R8) :: Prim_wall(nprim), mil, kl, Gradient(nprim,2), Stress(2), Flux(nprim)
    real(R8) :: r_ave


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l
    Dist   = 1d-20
    M      = Blk % M (Im,Jm) % c

    Prim = Blk % P(:,Im,Jm)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    Gradient = 0d0

    rho_wall = Prim(np) / ( Rgas * T_wall )
    Prim_wall(1:nsc) = Prim(1:nsc) * rho_wall / rho

    call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, T_wall, mil, kl )

    if ( present(w_wall) ) then
      Gradient(nu:nv,Dir) = ( Prim(nu:nv) - w_wall ) * modfm3
    else
      Gradient(nu:nv,Dir) = Prim(nu:nv) * modfm3
    end if
    Gradient(np,Dir)    = ( T - T_wall ) * modfm3

    if (model==2)  then
      dl = blk % yn(im,jm)
      call RANS_Set_Wall_Values ( mil, Prim_Wall(nt:nprim), dl )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3
    end if

    Gradient = matmul ( Gradient, M )

    select case(Fm)
      case(1)
        r_ave = 0.5d0 * ( blk % node(im-1,jm)%c(2) + blk % node(im-1,jm-1)%c(2) )
      case(2)
        r_ave = 0.5d0 * ( blk % node(im,jm)%c(2) + blk % node(im,jm-1)%c(2) )
      case(3)
        r_ave = 0.5d0 * ( blk % node(im,jm-1)%c(2) + blk % node(im-1,jm-1)%c(2) )
      case(4)
        r_ave = 0.5d0 * ( blk % node(im,jm)%c(2) + blk % node(im-1,jm)%c(2) )
    end select
    Stress = Stress_Vector ( Gradient(nu:nv,:), Normal, mil, 0d0, Prim(nt:), r_ave, 0.d0 )

    Flux = 0d0   ! mirror MOSE full pre-zero (covers soot/passive band np+1:nt-1, not just 1:nsc)
    Flux(nu:nv) = Stress * int_l
    Flux(np)    = int_l * kl * dot_product ( Gradient(np,:), Normal )

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), rans_variables=Prim_Wall(nt:nprim),    &
                                 vel_gradient=Gradient(nu:nv,:), rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, int_l=int_l, normal=normal,          &
                                 vr=0.d0, r_ave=r_ave, dist=dist )
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), rho, rho_wall, Blk % P(nt:nprim,Ig,Jg) )
    endif

    Blk % r (:,Im,Jm) = Blk % r (:,Im,Jm) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties ( stress=Stress, pw=prim(np), Tw=T_wall,        &
                                                      rhow=rho_wall, mu=mil,                        &
                                                      y=blk%dl(im,jm)%c(dir)*0.5d0,                 &
                                                      qw=kl * dot_product ( Gradient(np,:), Normal ), &
                                                      exit_array=Ovar )

  end subroutine BC_Wall_Temperature


end module Q2D_Lib_BC_Fluxes_Wall_Temperature
