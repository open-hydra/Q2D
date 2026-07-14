module Q2D_Lib_BC_Fluxes_Connection
  ! 4-cell stencil across the block boundary (2D reduction: 4 faces, no Km).
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Base_Types_m,           only: Q2D_tensor_2D_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model, gc
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_Reconstruction,     only: Reconstruction
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index
  use Q2D_Mod_Riemann,            only: Riemann
  use Q2D_Lib_Diffusive,          only: Compute_Diffusive_Flux
  use Q2D_Lib_Fluid,              only: V_dot_G_face_avg
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, f_ss, check_gas_state

  implicit none
  private

  public :: BC_Connection_Eul, BC_Connection_Visc

contains


  subroutine BC_Connection_Eul ( Im, Jm, Fm, Blk, SD_limiter, SD_riemann )
    ! 2D-reduced (no Km, no nw). SD_limiter gates beta_ for
    ! Q2D_Lib_Reconstruction; SD_riemann gates beta_ for Riemann (variable
    ! reuse). HLLC+/HLLE++ Tramel SD-aware solvers receive beta_ via
    ! SD_riemann; classic Riemann solvers ignore the arg.
    implicit none
    integer, intent(in) :: Im, Jm, Fm
    logical, intent(in) :: SD_limiter, SD_riemann
    type(Q2D_block_type), intent(inout) :: Blk
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dirm, Face_i, Face_j
    integer :: I1, J1, I2, J2, I3, J3, I4, J4
    integer :: error1, error4
    real(R8), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4
    real(R8), dimension(nprim) :: Prim_L, Prim_R, Flux
    real(R8) :: Normal(2), int_l, Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34
    real(R8) :: rho_L, rho_R, Rtot_L, Rtot_R, a_L, a_R, F_r, F_u, F_v, F_E, Fmass
    real(R8) :: su, sel_L, sel_R, beta_


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dirm, Im, Jm, Face_i, Face_j )

    ! Stencil construction: 4 points for 2nd order reconstruction from both sides at boundary face
    I1 = modfm  * ( Im - guide(Fm,1) * gc ) + modfm1 * ( Im + guide(Fm,1) )
    J1 = modFm  * ( Jm - guide(Fm,2) * gc ) + modFm1 * ( Jm + guide(Fm,2) )
    I2 = modFm1 * Im + modFm  * ( Im - guide(Fm,1) )
    J2 = modFm1 * Jm + modFm  * ( Jm - guide(Fm,2) )
    I3 = modFm1 * ( Im - guide(Fm,1) ) + modFm * Im
    J3 = modFm1 * ( Jm - guide(Fm,2) ) + modFm * Jm
    I4 = modFm1 * ( Im - guide(Fm,1) * gc ) + modFm * ( Im + guide(Fm,1) )
    J4 = modFm1 * ( Jm - guide(Fm,2) * gc ) + modFm * ( Jm + guide(Fm,2) )

    ! Primitive variables for the stencil
    Prim1 = blk % P(:,I1,J1)
    Prim2 = blk % P(:,I2,J2)
    Prim3 = blk % P(:,I3,J3)
    Prim4 = blk % P(:,I4,J4)

    ! SD_limiter gates beta_.
    if (SD_limiter) then
      beta_ = blk % beta(Im,Jm)
    else
      beta_ = 1d0
    end if

    Normal = blk % dir(Dirm) % f(Face_i,Face_j) % n
    int_l  = blk % dir(Dirm) % f(Face_i,Face_j) % int_l
    Dl1 = blk % Dl(I1,J1) % c(Dirm) * 0.5d0
    Dl2 = blk % Dl(I2,J2) % c(Dirm) * 0.5d0
    Dl3 = blk % Dl(I3,J3) % c(Dirm) * 0.5d0
    Dl4 = blk % Dl(I4,J4) % c(Dirm) * 0.5d0

    Dist23 = Dl2 + Dl3
    Dist12 = Dl1 + Dl2
    Dist34 = Dl3 + Dl4

    call Reconstruction ( Prim1, Prim2, Prim3, Prim4, Dist12, Dist23, Dist34, Dl2, Dl3, beta_, Prim_L, Prim_R )

    ! The connection seam is an interior face split across blocks: it must guard the
    ! reconstructed interface state exactly like Q2D_Lib_Convective::Convective_Flux
    ! (floor species + abort on p<0 / T-out-of-table), which the bare seam previously
    ! skipped. Mirror of Q2D_Lib_BC_Fluxes_Connection (coordinated back-port).
    call check_gas_state ( Prim_L(1:nsc), Prim_L(np), error1 )
    call check_gas_state ( Prim_R(1:nsc), Prim_R(np), error4 )
    if (error1 /= 0 .or. error4 /= 0) then
      write(*,*) '[ERROR] Non-physical state detected at the interface.'
      error stop
    end if

    ! Auxiliary variables for state (L)
    call co_rotot_Rtot ( Prim_L(1:nsc), rho_L, Rtot_L )
    a_L = f_ss ( Prim_L(1:nsc), Prim_L(np), rho_L, Rtot_L )

    ! Auxiliary variables for state (R)
    call co_rotot_Rtot ( Prim_R(1:nsc), rho_R, Rtot_R )
    a_R = f_ss ( Prim_R(1:nsc), Prim_R(np), rho_R, Rtot_R )

    ! SD_riemann gates beta_.
    ! SD-aware Riemann (HLLC+/HLLE++) uses beta_ for blending; classic ones ignore it.
    if (SD_riemann) then
      beta_ = blk % beta(Im,Jm)
    else
      beta_ = 1d0
    end if

    call Riemann ( Prim_L(1:nsc), Prim_L(nu), Prim_L(nv), Prim_L(np), a_L, rho_L, &
                   Prim_R(1:nsc), Prim_R(nu), Prim_R(nv), Prim_R(np), a_R, rho_R, &
                   beta_, normal(1), normal(2), F_r, F_u, F_v, F_E )

    su = sign ( 0.5d0, F_r )
    Sel_L = 0.5d0 + su
    Sel_R = su - 0.5d0

    Fmass = F_r * int_l
    Flux(1:nsc) = Fmass * ( Sel_L*Prim_L(1:nsc)/rho_L - Sel_R*Prim_R(1:nsc)/rho_R )
    Flux(nu) = F_u * int_l
    Flux(nv) = F_v * int_l
    Flux(np) = F_E * int_l

    ! Gate `nprim>np` extrapolates ALL post-pressure variables (turbulence +
    ! soot + passive) across the connection. The Q2D legacy gate
    ! `model==2 / nt:nprim` lost soot/passive in non-RANS mode (Q2D legacy
    ! oversight, not intentional preservation).
    if (nprim > np) Flux(np+1:nprim) = Fmass * ( sel_L*Prim_L(np+1:nprim)/rho_L - sel_R*Prim_R(np+1:nprim)/rho_R )

    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + Flux * modfm2

  end subroutine BC_Connection_Eul


  subroutine BC_Connection_Visc ( Im, Jm, Fm, Blk, Mg, Pg, Sc, Sct, Prt, soot_enabled )
    ! Signature: Sc/Sct/Prt/soot_enabled passed as parameters (vs the previous
    ! Q2D version that read obj_rans/obj_soot from a singleton). No Km arg
    ! (2D reduction).
    !
    ! NOTE: parameter Mg is currently unused in the body. Kept for interface
    ! uniformity with BC_Connection_Eul.
    implicit none
    integer, intent(in) :: Im, Jm, Fm
    type(Q2D_tensor_2D_type), intent(in) :: Mg
    real(R8),                 intent(in) :: Pg(nprim,4)
    type(Q2D_block_type),  intent(inout) :: Blk
    real(R8), intent(in) :: Sc, Sct, Prt
    logical,  intent(in) :: soot_enabled
    ! Local
    integer :: dir, modfm2, Face_i, Face_j
    integer :: Ig, Jg
    real(R8) :: Normal(2), int_l, Waldis, M(2,2)
    real(R8), dimension(nprim) :: Prim_loc, Prim_ghost, Visc_loc, Visc_ghost
    real(R8), dimension(nprim) :: Visc_ip, Visc_im, Visc_jp, Visc_jm
    real(R8), dimension(nprim) :: Visc_ghost3, Visc_ghost4
    real(R8) :: Prim(nprim), Gradient_loc(nprim,2), Gradient_ghost(nprim,2)
    real(R8) :: Flux_loc(nprim), Flux_ghost(nprim)
    real(R8) :: r_ave, VG_face, G_loc(2)


    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )
    modfm2 = 1 - 2 * mod (Fm,2)

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Normal = blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l  = blk % dir(Dir) % f(Face_i,Face_j) % int_l
    M = 0.5d0 * ( Blk % M(Im,Jm) % c + Blk % M(Ig,Jg) % c )
    Waldis = 0.5d0 * ( Blk % yn(Im,Jm) + Blk % yn(Ig,Jg) )

    Prim_loc   = Blk % P(:,Im,Jm)
    Prim_ghost = Pg(:,1)

    call Visc_Variables ( Prim_loc,             Visc_loc    )
    call Visc_Variables ( Blk % P(:,Im+1,Jm),   Visc_ip     )
    call Visc_Variables ( Blk % P(:,Im-1,Jm),   Visc_im     )
    call Visc_Variables ( Blk % P(:,Im,Jm+1),   Visc_jp     )
    call Visc_Variables ( Blk % P(:,Im,Jm-1),   Visc_jm     )
    call Visc_Variables ( Prim_ghost,           Visc_ghost  )
    call Visc_Variables ( Pg(:,3),              Visc_ghost3 )
    call Visc_Variables ( Pg(:,4),              Visc_ghost4 )

    select case ( Fm )
      case ( 1 : 2 )
        Gradient_loc   (:,1) = ( Visc_ghost - Visc_loc ) * modfm2
        Gradient_loc   (:,2) = ( Visc_jp - Visc_jm ) * 0.5d0
        Gradient_ghost (:,1) = ( Visc_loc - Visc_ghost ) * (-modfm2)
        Gradient_ghost (:,2) = ( Visc_ghost4 - Visc_ghost3 ) * 0.5d0
      case ( 3 : 4 )
        Gradient_loc   (:,1) = ( Visc_ip - Visc_im ) * 0.5d0
        Gradient_loc   (:,2) = ( Visc_ghost - Visc_loc ) * modfm2
        Gradient_ghost (:,1) = ( Visc_ghost4 - Visc_ghost3 ) * 0.5d0
        Gradient_ghost (:,2) = ( Visc_loc - Visc_ghost ) * (-modfm2)
    end select

    Gradient_loc   = matmul ( Gradient_loc,   M )
    Gradient_ghost = matmul ( Gradient_ghost, M )

    Prim = ( Prim_loc + Prim_ghost ) * 0.5d0

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

    ! Q2D Def A: V·G at the connection face (canonical convention ⟨V⟩·G_loc for
    ! continuity of h(x,y) across coincident connections — G_neigh=G_loc fallback).
    G_loc(1) = Blk%area_geom_term_x(Im,Jm)
    G_loc(2) = Blk%area_geom_term_y(Im,Jm)
    VG_face = V_dot_G_face_avg ( V_loc=Prim_loc(nu:nv), V_neigh=Prim_ghost(nu:nv), &
                                 G_loc=G_loc,           G_neigh=G_loc )

    call Compute_Diffusive_Flux ( Prim, Gradient_loc,   int_l, Normal, Waldis, r_ave, Flux_loc,   &
                                  Sc=Sc, Sct=Sct, Prt=Prt,                                        &
                                  soot_enabled=soot_enabled, V_dot_G=VG_face )
    call Compute_Diffusive_Flux ( Prim, Gradient_ghost, int_l, Normal, Waldis, r_ave, Flux_ghost, &
                                  Sc=Sc, Sct=Sct, Prt=Prt,                                        &
                                  soot_enabled=soot_enabled, V_dot_G=VG_face )

    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) - 0.5d0 * modfm2 * ( Flux_loc + Flux_ghost )

  contains

    subroutine Visc_Variables ( PrimIn, ViscOut )
      implicit none
      real(R8), intent(in)  :: PrimIn (nprim)
      real(R8), intent(out) :: ViscOut(nprim)
      ! Local
      real(R8) :: rho_l, Rgas_l

      call co_rotot_Rtot ( PrimIn(1:nsc), rho_l, Rgas_l )
      ViscOut(1:nsc) = PrimIn(1:nsc) / rho_l
      ViscOut(nu:nv) = PrimIn(nu:nv)
      ViscOut(np)    = PrimIn(np) / ( rho_l * Rgas_l )
      ! Gate `nprim>np` extrapolates ALL post-pressure variables (turbulence +
      ! soot + passive). The Q2D legacy gate `model==2/nt:nprim` lost soot/passive
      ! in non-RANS mode (Q2D legacy oversight, same pattern fixed in
      ! BC_Connection_Eul top-level; the Visc_Variables nested subroutine had not
      ! been updated).
      if (nprim > np) ViscOut(np+1:nprim) = PrimIn(np+1:nprim) / rho_l

    end subroutine Visc_Variables

  end subroutine BC_Connection_Visc


end module Q2D_Lib_BC_Fluxes_Connection
