module Q2D_Lib_Diffusive
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Diffusive_Flux, Compute_Diffusive_Flux

contains

  subroutine Diffusive_Flux ( normal, int_l, waldis1, waldis2, Prim1, Prim2, Prim3, Prim4, &
                              Prim5, Prim6, M1, M2, Res1, Res2, r2, r1, a, b, &
                              Sc, Sct, Prt, soot_enabled, Gx1, Gy1, Gx2, Gy2 )
    ! Sc, Sct, Prt, soot_enabled are passed as explicit parameters.
    ! `r1, r2`: nodal radial coordinates for the axisymmetric average r_ave;
    ! `Gx1/Gy1/Gx2/Gy2` (optional): geometric terms for the face average V·G.
    ! `int_l` is the interface length (2D analogue of the 3D face area).
    use Q2D_Global_m,            only: nu, nv, np, nt, nsc, nprim, model
    use Q2D_Lib_Fluid,           only: V_dot_G_face_avg
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot
    implicit none
    integer, intent(in)  :: a, b
    logical, intent(in)  :: soot_enabled
    real(R8), intent(in) :: Sc, Sct, Prt
    real(R8), intent(in) :: normal(2), int_l, waldis1, waldis2
    real(R8), intent(in), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4, Prim5, Prim6
    real(R8), intent(in), dimension(2,2) :: M1, M2
    real(R8), intent(in) :: r2, r1
    real(R8), intent(inout), dimension(nprim) :: Res1, Res2
    real(R8), intent(in), optional :: Gx1, Gy1, Gx2, Gy2
    real(R8) :: rho1, Rgas, T1, rho2, T2, Gradient(nprim,2), Prim(nprim), M(2,2), waldis, Flux(nprim), r_ave
    real(R8) :: VG_face

    r_ave = 0.5d0 * ( r1 + r2 )
    call co_rotot_Rtot ( Prim1(1:nsc), rho1, Rgas )
    T1 = Prim1(np) / ( rho1 * Rgas )
    call co_rotot_Rtot ( Prim2(1:nsc), rho2, Rgas )
    T2 = Prim2(np) / ( rho2 * Rgas )

    Gradient ( 1:nsc, a ) = Prim2 ( 1:nsc ) / rho2 - Prim1 ( 1:nsc ) / rho1
    Gradient ( np, a )    = T2 - T1
    Gradient ( nu:nv, a ) = Prim2 ( nu:nv ) - Prim1 ( nu:nv )
    ! The `if (nprim>np)` gate fills the whole band of transported scalars beyond
    ! the pressure (soot+passives+RANS, np+1:nprim), not just the RANS one (nt:nprim).
    ! `waldis` is assigned inside this gate: it is read ONLY by the
    ! Eddy_Viscosity/RANS_Diffusive_Flux calls (active with model==2) in
    ! Compute_Diffusive_Flux, so in laminar cases it may remain uninitialized
    ! without consequences.
    if (nprim>np) then
      Gradient ( np+1:nprim, a ) = Prim2 ( np+1:nprim ) / rho2 - Prim1 ( np+1:nprim ) / rho1
      waldis = 0.5d0 * ( waldis1 + waldis2 )
    end if

    call Tangential_Gradient ( Prim3, Prim4, Prim5, Prim6, Gradient(:,b) )

    M = 0.5d0 * ( M1 + M2 )
    Gradient = matmul ( Gradient, M )

    Prim = 0.5d0 * ( Prim1 + Prim2 )

    ! Convention for the face geometric term: V·G_face = ⟨V⟩_face · ⟨G⟩_face
    ! (V_dot_G_face_avg in Lib_Fluid). Consistent with BC_Symmetry_Visc/BC_Connection_Visc,
    ! which impose G_ghost = G_loc at the boundary. With area_geom_term identically
    ! zero, VG_face = 0 and the term is disabled.
    VG_face = 0.d0
    if (present(Gx1) .and. present(Gy1) .and. present(Gx2) .and. present(Gy2)) then
      VG_face = V_dot_G_face_avg ( V_loc=Prim1(nu:nv), V_neigh=Prim2(nu:nv), &
                                   G_loc=[Gx1,Gy1],   G_neigh=[Gx2,Gy2] )
    end if

    call Compute_Diffusive_Flux ( Prim, Gradient, int_l, normal, waldis, r_ave, Flux, &
                                  Sc=Sc, Sct=Sct, Prt=Prt,                              &
                                  soot_enabled=soot_enabled, V_dot_G=VG_face )

    Res1 = Res1 - Flux
    Res2 = Res2 + Flux

  end subroutine Diffusive_Flux


  subroutine Tangential_Gradient ( Prim1, Prim2, Prim3, Prim4, Gradient )
    use Q2D_Global_m, only: nu, nv, np, nt, nsc, nprim, model
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot
    implicit none
    real(R8), intent(in),  dimension(nprim) :: Prim1, Prim2, Prim3, Prim4
    real(R8), intent(out), dimension(nprim) :: Gradient
    real(R8) :: rho1, rho2, rho3, rho4, T1, T2, T3, T4, Rgas

    call co_rotot_Rtot ( Prim1(1:nsc), rho1, Rgas )
    T1 = Prim1(np) / ( rho1 * Rgas )
    call co_rotot_Rtot ( Prim2(1:nsc), rho2, Rgas )
    T2 = Prim2(np) / ( rho2 * Rgas )
    call co_rotot_Rtot ( Prim3(1:nsc), rho3, Rgas )
    T3 = Prim3(np) / ( rho3 * Rgas )
    call co_rotot_Rtot ( Prim4(1:nsc), rho4, Rgas )
    T4 = Prim4(np) / ( rho4 * Rgas )

    Gradient ( 1:nsc ) = ( Prim2 ( 1:nsc ) / rho2 - Prim1 ( 1:nsc ) / rho1 + &
                           Prim4 ( 1:nsc ) / rho4 - Prim3 ( 1:nsc ) / rho3 ) * 0.25d0

    Gradient ( np ) = ( T2 - T1 + T4 - T3 ) * 0.25d0

    Gradient ( nu:nv ) = ( Prim2 ( nu:nv ) - Prim1 ( nu:nv ) + Prim4 ( nu:nv ) - Prim3 ( nu:nv ) ) * 0.25d0

    ! Riempie l'intera banda degli scalari oltre la pressione (soot+passivi+RANS).
    if (nprim>np) then
      Gradient ( np+1:nprim ) = ( Prim2 ( np+1:nprim ) / rho2 - Prim1 ( np+1:nprim ) / rho1 + &
                                  Prim4 ( np+1:nprim ) / rho4 - Prim3 ( np+1:nprim ) / rho3 ) * 0.25d0
    end if

  end subroutine Tangential_Gradient


  subroutine Compute_Diffusive_Flux ( Prim, Gradient, int_l, normal, waldis, r_ave, Flux, &
                                      Sc, Sct, Prt, soot_enabled, V_dot_G )
    ! Sc, Sct, Prt, soot_enabled sono passati come parametri espliciti.
    ! `r_ave`: raggio medio di faccia (assialsimmetria); `V_dot_G` (opzionale):
    ! termine geometrico mediato sulla faccia.
    ! `int_l` è la lunghezza dell'interfaccia (analogo 2D dell'area di faccia 3D).
    use Q2D_Global_m,            only: nu, nv, np, nt, nsc, nprim, model
    use Q2D_Lib_Fluid,           only: Stress_Vector
    use Q2D_Lib_RANS,            only: Eddy_Viscosity, RANS_Diffusive_Flux
    use Q2D_Mod_Soot,            only: Soot_Diffusive_Flux
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, f_cp_expr, &
                                        co_k_mi_lam_Wilke_expr, f_tabT_expr, h_tab
    implicit none
    real(R8), intent(in)  :: Prim(nprim), Gradient(nprim,2), int_l, normal(2), waldis, r_ave
    real(R8), intent(in)  :: Sc, Sct, Prt
    logical,  intent(in)  :: soot_enabled
    real(R8), intent(out) :: Flux(nprim)
    real(R8), intent(in), optional :: V_dot_G
    integer :: s, T_i, Tint(2)
    real(R8) :: rho, Rgas, T, Tdiff, cp, mil, kl, mie, kappa, h_alpha, h_mix, V_c, VG
    real(R8) :: Dm(nsc), stress(2), DiffHFlux, DmGradYi(nsc)

    VG = 0d0
    if (present(V_dot_G)) VG = V_dot_G

    call co_rotot_Rtot ( Prim(1:nsc), rho, Rgas )
    T = Prim(np) / ( rho * Rgas )

    T_i     = int(T)
    Tdiff   = T - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    cp = f_cp_expr ( Prim(1:nsc), Tint, Tdiff, rho )
    call co_k_mi_lam_Wilke_expr ( Prim(1:nsc), rho, Tint, Tdiff, mil, kl )

    mie = 0d0
    if (model==2) then
      call Eddy_Viscosity ( mut=mie, rans_variables=Prim(nt:nprim),                 &
                            mul=mil, rho=rho, vel_gradient=Gradient(nu:nv,:),       &
                            vr=Prim(nv), r_ave=r_ave, walldist=waldis, V_dot_G=VG   )
    end if

    Dm (1:nsc) = ( mil/Sc + mie/Sct ) / rho

    kappa = kl + mie*cp/Prt

    Stress = Stress_Vector ( Gradient(nu:nv,:), normal, mil, mie, Prim(nt:), r_ave, Prim(nv), V_dot_G=VG )

    DiffHFlux = 0.d0
    h_mix     = 0.d0
    Flux      = 0.d0   ! pre-azzeramento completo: la banda soot/passivi np+1:nt-1
                       ! (e quella RANS quando model/=2) contribuisce 0 a Res,
                       ! evitando valori indefiniti.

    do s = 1, nsc
      DmGradYi(s) = Dm(s) * dot_product ( Gradient(s,:), normal )
      Flux(s)     = rho * int_l * DmGradYi(s)
      h_alpha     = f_tabT_expr ( s, h_tab, Tint, Tdiff )
      DiffHFlux   = DiffHFlux + DmGradYi(s) * h_alpha
      h_mix       = h_mix + ( Prim(s) / rho ) * h_alpha
    enddo

    V_c = sum(DmGradYi(1:nsc))

    Flux(1:nsc) = Flux(1:nsc) - sum(Flux(1:nsc))*Prim(1:nsc)/rho

    Flux(nu:nv) = int_l * Stress

    Flux(np) = int_l * ( dot_product ( Stress, Prim(nu:nv) ) + &
               kappa * dot_product ( Gradient(np,:), normal ) + (DiffHFlux - h_mix*V_c) * rho )

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), rans_variables=Prim(nt:nprim),  &
                                 vel_gradient=Gradient(nu:nv,:),                       &
                                 rans_gradient=Gradient(nt:nprim,:),                   &
                                 mul=mil, rho=rho, int_l=int_l, normal=normal,         &
                                 vr=Prim(nv), r_ave=r_ave, dist=waldis, mut=mie )
    end if

    ! Flussi diffusivi del soot (termoforesi) via Q2D_Mod_Soot::Soot_Diffusive_Flux;
    ! attivazione governata dal parametro `soot_enabled`.
    if (soot_enabled) then
      call Soot_Diffusive_Flux ( flux=Flux, prim=Prim, T_gradient=Gradient(np,:), &
                                 mu=mil+mie, rho=rho, T=T, int_l=int_l, normal=normal )
    end if

  end subroutine Compute_Diffusive_Flux

end module Q2D_Lib_Diffusive
