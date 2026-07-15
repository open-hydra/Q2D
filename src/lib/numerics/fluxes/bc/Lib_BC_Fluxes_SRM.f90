module Q2D_Lib_BC_Fluxes_SRM
  ! BC 502: combustione del grain SRM. Legge di combustione
  ! `mdot = aCoeff * (P/pRef)^n` (tipo Vieille) + proprietà di fiamma
  ! adiabatica del propellente (Taf, haf, BC_ci).
  !
  ! Formulazione 2D:
  !   * indici di cella ghost Ig/Jg; Normal(2), M(2,2), Gradient(nprim,2),
  !     Stress(2), Blowing(2);
  !   * int_l: lunghezza integrata della faccia;
  !   * residuo di rango 3.
  !
  ! V_dot_G: Stress_Vector è invocata senza l'argomento opzionale V_dot_G
  ! (default 0), coerente con l'assenza di variazione d'area locale a parete.
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Parameters_m,           only: guide
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index, Compute_Wall_Properties
  use Q2D_Lib_Fluid,              only: Stress_Vector
  use Q2D_Lib_RANS,               only: RANS_Set_Blowing_Wall, RANS_Diffusive_Flux, RANS_Extrapolate_Wall
  use FLINT_Lib_Thermodynamic,    only: co_rotot_Rtot, co_k_mi_lam_Wilke

  implicit none
  private

  public :: BC_SRM

contains


  subroutine BC_SRM ( Im, Jm, Fm, Blk, aCoeff, n, pRef, Taf, haf, BC_ci, Ovar )
    implicit none
    integer,  intent(in) :: Im, Jm, Fm
    real(R8), intent(in) :: aCoeff, n, pRef, Taf, haf, BC_ci(nsc)
    type(Q2D_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(7), intent(inout) :: Ovar
    ! Locali
    integer  :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Ig, Jg
    real(R8) :: Normal(2), int_l, Dist, dl, M(2,2), Prim(nprim), rho, Rgas, Gradient(nprim,2)
    real(R8) :: rho_wall, Prim_Wall(nprim), mil, kl, mdot, T
    real(R8) :: Stress(2), Blowing(2), Flux(nprim)
    real(R8) :: q_conv
    real(R8) :: r_ave

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    ! Metriche di faccia
    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l
    Dist   = Blk % yn(Im,Jm)
    M      = Blk % M (Im,Jm) % c

    ! Variabili della cella di bordo
    Prim = Blk % P(:,Im,Jm)
    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    rho_wall = Prim(np) / (Rgas*Taf)

    mdot    = aCoeff * (Prim(np)/pRef)**n * (-modfm2)
    Blowing = (mdot/rho_wall) * Normal ! Velocità di blowing

    ! r_ave per Stress_Vector (termini 2D-assialsimmetrici / variazione d'area)
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

    q_conv   = 0.d0
    Stress   = 0.d0
    Gradient = 0d0
    if (model>0) then
      Prim_wall(1:nsc) = BC_ci * rho_wall ! Densità parziali a parete
      T = Prim(np) / ( rho * Rgas )
      call co_k_mi_lam_Wilke ( Prim_wall(1:nsc), rho_wall, Taf, mil, kl )

      Gradient(np,Dir) = ( T - Taf ) * modfm3
      Gradient(np,:)   = matmul ( Gradient(np,:), M )
      q_conv           = kl * dot_product( Gradient(np,:), Normal )

      Gradient(nu:nv,Dir) = ( Prim(nu:nv) - Blowing ) * modfm3    ! Gradiente di velocità dovuto al blowing
      Gradient(nu:nv,:)   = matmul ( Gradient(nu:nv,:), M )       ! Trasformazione del gradiente di velocità
      Stress = Stress_Vector ( Gradient(nu:nv,:), Normal, mil, 0d0, Prim(nt:), r_ave, 0.d0 )
    endif

    ! Flussi
    Flux(1:nsc) = - mdot * int_l * BC_ci
    Flux(nu:nv) = ( Stress - mdot * Blowing ) * int_l
    Flux(np)    = int_l * ( q_conv - mdot * ( haf + 0.5d0 * sum( Blowing**2 ) ) &
                             + dot_product( Stress, Blowing ) )
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
      ! Estrapolazione in cella ghost delle variabili RANS
      Ig = Im - guide(Fm,1)
      Jg = Jm - guide(Fm,2)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg) )
    endif

    ! Aggiornamento del residuo
    Blk % r (:,Im,Jm) = Blk % r (:,Im,Jm) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties ( stress=Stress, pw=prim(np), qw=q_conv,  &
                                                       mdot=mdot, y=blk%dl(im,jm)%c(dir)*0.5d0, &
                                                       Tw=Taf, rhow=rho_wall, mu=mil,          &
                                                       exit_array=Ovar )

  end subroutine BC_SRM


end module Q2D_Lib_BC_Fluxes_SRM
