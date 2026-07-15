module Q2D_Lib_BC_Fluxes_Symmetry
  ! BC di simmetria (formulazione 2D: 4 facce, Normal(2)).
  !
  ! Asse 2D-axi r→0: la riflessione `Vel2 = Vel1 - 2(V·n)n` lascia ⟨V⟩·n = 0
  ! esatto, quindi la convenzione canonica per V·G di faccia usa ⟨V⟩·G_loc
  ! (nessun termine normale spurio anche se G_n_face != 0).
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
    ! Per una parete stazionaria in un sistema di riferimento rotante
    ! (casing/shroud non cilindrico): calcola
    ! `Un_wall = (Ω × r_face) · Normal`, non nullo quando la parete non è
    ! perfettamente cilindrica attorno all'asse di rotazione → entra nel mirror
    ! di Ghost_Prim e nello stato L/R del problema di Riemann. Per parete
    ! cilindrica (default, inclusa la simmetria pura): Un_wall = 0 → mirror
    ! standard.
    implicit none
    integer, intent(in) :: Bm, Im, Jm, Fm
    type(Q2D_block_type), intent(inout) :: Blk
    ! Locali
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

    ! Velocità normale della faccia di parete nel sistema di riferimento rotante.
    ! Hub rotante col riferimento O piano di pura simmetria: Un_wall = 0.
    ! Casing/shroud stazionario non cilindrico: Un_wall = -(Ω × r_face) · n ≠ 0.
    Un_wall = 0.0_R8
    if ( obj_rot%enabled .and. RF_Is_Stationary_Face(Bm, Fm) ) then
      call RF_Face_Center ( Blk%node, Im, Jm, Fm, r_fc )
      Un_wall = dot_product ( RF_Wall_Velocity(r_fc), Normal )
    end if

    Un = dot_product ( Bound_Prim (nu:nv), Normal )
    Velt = Bound_Prim (nu:nv) - Un * Normal

    ! Velocità mirror della cella ghost.
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
    ! Costruzione dello stato L/R di Riemann con shift della velocità normale
    ! per parete stazionaria: Un_wall ≠ 0 → shift di L e R di Un_wall*Normal;
    ! mirror simmetrico standard altrimenti.
    if ( Un_wall /= 0.0_R8 ) then
      Un1 = (Un - Un_wall) * modFm2
      U1  =  Un_wall * Normal + Un1 * Normal + Velt
      U4  =  Un_wall * Normal - Un1 * Normal + Velt
    else
      Un1 = Un * modFm2
      U1  =  Un1 * Normal + Velt
      U4  = -Un1 * Normal + Velt
    end if
    ! beta = 1d0: blending accoppiato allo schema SD, peso HLLC pieno.
    call Riemann ( Face_Prim(1:nsc), U1(1), U1(2), Face_Prim(np), Sound, rotot, &
                   Face_Prim(1:nsc), U4(1), U4(2), Face_Prim(np), Sound, rotot, &
                   1d0, Normal(1), Normal(2), F_r, F_u, F_v, F_e )

    Flux = 0d0
    Flux (nu) = F_u * int_l
    Flux (nv) = F_v * int_l

    Blk % r(nu:nv,Im,Jm) = Blk % r(nu:nv,Im,Jm) + Modfm2 * Flux(nu:nv)

  end subroutine BC_Symmetry_Eul


  subroutine BC_Symmetry_Visc ( Im, Jm, Fm, Blk, RSM )
    ! RSM (gating del modello agli sforzi di Reynolds SSGLRR, passato da
    ! Mod_BC_Fluxes): per l'RSM lo sforzo di Reynolds nel piano R12 è DISPARI
    ! attraverso la linea di simmetria, quindi porta un flusso diffusivo GGDH
    ! non nullo che deve entrare nel residuo. Gli scalari SA/SST sono pari ->
    ! flusso nullo -> blocco saltato.
    implicit none
    integer, intent(in) :: Im, Jm, Fm
    logical, intent(in) :: RSM
    type(Q2D_block_type), intent(inout) :: Blk
    ! Locali
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

    ! V·G alla faccia di simmetria (r→0; convenzione canonica ⟨V⟩·G_loc —
    ! V_n_face=0 esatto per la riflessione).
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

    ! Flusso diffusivo degli sforzi di Reynolds al piano di simmetria (solo
    ! RSM): RSM_Symmetry impone il valore di faccia (R12=0) e il ghost dispari
    ! (-R12); il gradiente normale alla parete (Prim-Prim_Wall) alimenta il
    ! flusso GGDH RANS_Diffusive_Flux, integrato nel residuo RANS.
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
