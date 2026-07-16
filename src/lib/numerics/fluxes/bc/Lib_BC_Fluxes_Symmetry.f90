module Q2D_Lib_BC_Fluxes_Symmetry
  ! Symmetry BC (2D formulation: 4 faces, Normal(2)).
  !
  ! 2D-axi axis r→0: the reflection `Vel2 = Vel1 - 2(V·n)n` leaves ⟨V⟩·n = 0
  ! exact, so the canonical convention for face V·G uses ⟨V⟩·G_loc
  ! (no spurious normal term even if G_n_face != 0).
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, nRANS, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_Limiters,           only: rlimiter
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index
  use Q2D_Mod_Riemann,            only: Riemann
  use Q2D_Lib_Fluid,              only: Stress_Vector, V_dot_G_face_avg
  use Q2D_Lib_RANS,               only: Eddy_Viscosity, RSM_Symmetry, RANS_Diffusive_Flux
  use Q2D_Lib_RotatingFrame,      only: RF_Is_Stationary_Face, RF_Face_Center, RF_Wall_Velocity
  use Q2D_Config_Types_m,         only: obj_rot
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, f_ss, f_laminarViscosity

  implicit none
  private

  public :: BC_Symmetry_Eul, BC_Symmetry_Visc

contains


  subroutine BC_Symmetry_Eul ( Bm, Im, Jm, Fm, Blk )
    ! For a stationary wall in a rotating reference frame
    ! (non-cylindrical casing/shroud): computes
    ! `Un_wall = (Ω × r_face) · Normal`, nonzero when the wall is not
    ! perfectly cylindrical about the rotation axis → enters the Ghost_Prim
    ! mirror and the L/R state of the Riemann problem. For a cylindrical wall
    ! (default, including pure symmetry): Un_wall = 0 → standard mirror.
    implicit none
    integer, intent(in) :: Bm, Im, Jm, Fm
    type(Q2D_block_type), intent(inout) :: Blk
    ! Locals
    integer :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, dir, Face_i, Face_j, s
    real(R8) :: Normal(2), int_l, Bound_Dl, Int_Dl, Un, Un_wall, Velt(2), Dist12, Dist23, Dist
    real(R8), dimension(nprim) :: Bound_Prim, Int_Prim, Ghost_Prim, Prim1, Prim3
    real(R8), dimension(nprim) :: Slope32, Slope21, Slope, Face_Prim
    real(R8) :: rotot, Rtot, Sound, Un1, F_r, F_u, F_v, F_E, Flux(nprim), U1(2), U4(2)
    real(R8) :: r_fc(2)

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)

    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    normal = Blk % dir(dir) % f(Face_i,Face_j) % n
    int_l = Blk % dir(dir) % f(Face_i,Face_j) % int_l
    Bound_Dl = 0.5d0 * Blk % dl(Im,Jm) % c(dir)
    Int_Dl = 0.5d0 * Blk % dl(Int_i,Int_j) % c(dir)

    Bound_Prim = Blk % P(:,Im,Jm)
    Int_Prim = Blk % P(:,Int_i,Int_j)

    ! Normal velocity of the wall face in the rotating reference frame.
    ! Hub rotating with the frame OR pure-symmetry plane: Un_wall = 0.
    ! Non-cylindrical stationary casing/shroud: Un_wall = -(Ω × r_face) · n ≠ 0.
    Un_wall = 0.0_R8
    if ( obj_rot%enabled .and. RF_Is_Stationary_Face(Bm, Fm) ) then
      call RF_Face_Center ( Blk%node, Im, Jm, Fm, r_fc )
      Un_wall = dot_product ( RF_Wall_Velocity(r_fc), Normal )
    end if

    Un = dot_product ( Bound_Prim (nu:nv), Normal )
    Velt = Bound_Prim (nu:nv) - Un * Normal

    ! Mirror velocity of the ghost cell.
    Ghost_Prim = Bound_Prim
    if ( Un_wall /= 0.0_R8 ) then
      Ghost_Prim (nu:nv) = Ghost_Prim (nu:nv) - 2d0 * (Un - Un_wall) * Normal
    else
      Ghost_Prim (nu:nv) = Ghost_Prim (nu:nv) - 2d0 * Un * Normal
    end if

    Prim1 = modfm * Ghost_Prim + modfm1 * Int_Prim
    Dist12 = 2.d0 * Bound_Dl * modfm + ( Bound_Dl + Int_Dl ) * modfm1
    Prim3 = modfm * Int_Prim + modfm1 * Ghost_Prim
    Dist23 = ( Bound_Dl + Int_Dl ) * modFm + 2d0 * Bound_Dl * modfm1

    Slope32 = ( Prim3 - Bound_Prim ) / Dist23
    Slope21 = ( Bound_Prim - Prim1 ) / Dist12

    Slope = 0d0
    do s = 1, np
      Slope(s) = rlimiter ( Slope32(s), Slope21(s) )
    end do

    Dist = Bound_Dl * modfm2
    Face_Prim = Bound_Prim + Slope * Dist

    call co_rotot_Rtot ( Face_Prim(1:nsc), rotot, Rtot )
    Sound = f_ss ( Face_Prim(1:nsc), Face_Prim(np), rotot, Rtot )
    Un  = dot_product ( Face_Prim(nu:nv), Normal )
    Velt = ( Face_Prim(nu:nv) - Un*Normal )
    ! Construction of the Riemann L/R state with normal-velocity shift
    ! for a stationary wall: Un_wall ≠ 0 → shift of L and R by Un_wall*Normal;
    ! standard symmetric mirror otherwise.
    if ( Un_wall /= 0.0_R8 ) then
      Un1 = (Un - Un_wall) * modFm2
      U1  =  Un_wall * Normal + Un1 * Normal + Velt
      U4  =  Un_wall * Normal - Un1 * Normal + Velt
    else
      Un1 = Un * modFm2
      U1  =  Un1 * Normal + Velt
      U4  = -Un1 * Normal + Velt
    end if
    ! beta = 1d0: blending coupled to the SD scheme, full HLLC weight.
    call Riemann ( Face_Prim(1:nsc), U1(1), U1(2), Face_Prim(np), Sound, rotot, &
                   Face_Prim(1:nsc), U4(1), U4(2), Face_Prim(np), Sound, rotot, &
                   1d0, Normal(1), Normal(2), F_r, F_u, F_v, F_e )

    Flux = 0d0
    Flux (nu) = F_u * int_l
    Flux (nv) = F_v * int_l

    Blk % r(nu:nv,Im,Jm) = Blk % r(nu:nv,Im,Jm) + Modfm2 * Flux(nu:nv)

  end subroutine BC_Symmetry_Eul


  subroutine BC_Symmetry_Visc ( Im, Jm, Fm, Blk, RSM )
    ! RSM (gating of the SSGLRR Reynolds-stress model, passed from
    ! Mod_BC_Fluxes): for the RSM the R12 Reynolds stress is ODD
    ! across the symmetry line, hence it carries a nonzero GGDH diffusive
    ! flux that must enter the residual. The SA/SST scalars are even ->
    ! zero flux -> block skipped.
    implicit none
    integer, intent(in) :: Im, Jm, Fm
    logical, intent(in) :: RSM
    type(Q2D_block_type), intent(inout) :: Blk
    ! Locals
    integer :: modfm, modfm1, modfm2, modfm3, dir, Face_i, Face_j, a, b, ig, jg
    real(R8) :: Normal(2), int_l, M(2,2), Prim(nprim), rho, Rgas, mil, mie
    real(R8) :: Prim_Wall(nprim), dist, Gradrij(nRANS,2)
    real(R8) :: Gradient(2,2), Stress(2), Flux(nprim)
    real(R8), dimension(2) :: Vel1, Vel2, Vel3, Vel4, Vel5, Vel6, Vel
    real(R8) :: r_ave, G_loc(2), VG_face

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l = Blk % dir(Dir) % f(Face_i,Face_j) % int_l
    dist = Blk % yn(Im,Jm)
    M = Blk % M (Im,Jm) % c

    Prim = Blk % P(:,Im,Jm)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    mil = f_laminarViscosity ( Prim(1:nsc), Prim(np), rho, Rgas )

    Gradient = 0d0

    select case (Fm)
      case(1:2)
        Vel1  = Blk % P(nu:nv,im,jm)
        Vel3  = Blk % P(nu:nv,im,jm-1)
        Vel4  = Blk % P(nu:nv,im,jm+1)
        a = 1
        b = 2
      case(3:4)
        Vel1  = Blk % P(nu:nv,im,jm)
        Vel3  = Blk % P(nu:nv,im-1,jm)
        Vel4  = Blk % P(nu:nv,im+1,jm)
        a = 2
        b = 1
    end select

    Vel2 = Vel1 - 2d0 * dot_product ( Vel1, Normal ) * Normal
    Vel5 = Vel3 - 2d0 * dot_product ( Vel3, Normal ) * Normal
    Vel6 = Vel4 - 2d0 * dot_product ( Vel4, Normal ) * Normal

    Gradient (:,a) = ( Vel2 - Vel1 ) * modfm2
    Gradient (:,b) = ( Vel4 - Vel3 + Vel6 - Vel5 ) * 0.25d0

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

    Vel = 0.5d0 * ( Vel1 + Vel2 )

    ! V·G at the symmetry face (r→0; canonical convention ⟨V⟩·G_loc —
    ! V_n_face=0 exact for the reflection).
    G_loc(1) = Blk%area_geom_term_x(Im,Jm)
    G_loc(2) = Blk%area_geom_term_y(Im,Jm)
    VG_face = V_dot_G_face_avg ( V_loc=Vel1, V_neigh=Vel2, G_loc=G_loc, G_neigh=G_loc )

    mie = 0d0
    if (model==2) then
      call Eddy_Viscosity ( mut=mie, rans_variables=Prim(nt:nprim),                    &
                            mul=mil, rho=rho, vel_gradient=Gradient,                   &
                            vr=0.5d0*(Prim(nv)+Vel2(2)), r_ave=r_ave, walldist=dist,   &
                            V_dot_G=VG_face                                            )
    end if

    Stress = Stress_Vector ( Gradient, Normal, mil, mie, Prim(nt:), r_ave, &
                             0.5d0*(Prim(nv)+Vel2(2)), V_dot_G=VG_face )

    Flux = 0d0
    Flux(nu:nv) = int_l * Stress
    Flux(np) = int_l * dot_product ( Vel, Stress )

    ! Reynolds-stress diffusive flux at the symmetry plane (RSM only):
    ! RSM_Symmetry imposes the face value (R12=0) and the odd ghost
    ! (-R12); the wall-normal gradient (Prim-Prim_Wall) feeds the
    ! GGDH flux RANS_Diffusive_Flux, integrated into the RANS residual.
    if ( RSM ) then
      Gradrij = 0d0
      ig = Im - guide(Fm,1)
      jg = Jm - guide(Fm,2)
      call RSM_Symmetry ( Prim(nt:nprim), Prim_Wall(nt:nprim), Blk%P(nt:nprim,Ig,Jg), Normal )
      Gradrij(:,dir) = ( Prim(nt:nprim) - Prim_Wall(nt:nprim) ) / rho * modfm3
      Gradrij = matmul ( Gradrij, M )
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), rans_variables=Prim_Wall(nt:nprim),  &
                                 vel_gradient=Gradient, rans_gradient=Gradrij,              &
                                 mul=mil, rho=rho, int_l=int_l, normal=Normal,              &
                                 vr=0.5d0*(Prim(nv)+Vel2(2)), r_ave=r_ave, dist=dist )
      Blk % r (nt:nprim,Im,Jm) = Blk % r (nt:nprim,Im,Jm) - Modfm2 * Flux(nt:nprim)
    endif

    Blk % r (nu:np,Im,Jm) = Blk % r (nu:np,Im,Jm) - Modfm2 * Flux(nu:np)

  end subroutine BC_Symmetry_Visc


end module Q2D_Lib_BC_Fluxes_Symmetry
