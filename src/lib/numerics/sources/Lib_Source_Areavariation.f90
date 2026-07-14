module Q2D_Lib_Source_Areavariation
  ! Q2D-extra: the MOSE peer has no numerics/sources/.
  ! (Renamed from Mod_Areavariation to the Lib_* naming convention.)
  !
  ! q2D area-variation source terms:
  !   Eulerian:  Source ~ V·G with G = -∇h/h (geometric gradient from Cross_Section)
  !   Viscous (Q2D Definition A): canonical formulation v7 §N16, Anderson-consistent
  !                               signs — wired via Q2D_Lib_Fluid::Stress_Tensor +
  !                               Q2D_Lib_RANS::{Eddy_Viscosity, RANS_Diffusive_SourceQ2D}.
  !
  ! Mutual exclusion `axisymmetric ⊕ area_variation` enforced in Register_Sources
  ! + double-check in Mod_Sources::Initialize_Sources.
  !
  ! Fallback eps=1.d-12 for Acell_center → 0.
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m, only: Q2D_block_type, Q2D_domain_type
  implicit none
  private

  public :: Precompute_Area_Geometry
  public :: Read_Block_Q2D_Area
  public :: Source_Area_Eul_Blk
  public :: Source_Area_Visc_Blk

contains

  subroutine Precompute_Area_Geometry ( blk )
    ! Precompute G_x = -dA/dx / A_center, G_y = -dA/dy / A_center cell-center.
    ! Uses blk%M(i,j)%c → correct order: called AFTER Compute_Metric_Tensor in
    ! Setup_Metrics. The area_geom_term_x/y fields are allocated by
    ! Mod_Allocate_Data → no self-allocate here.
    implicit none
    type(Q2D_block_type), intent(inout) :: blk
    integer  :: i, j
    real(R8) :: dA_ij(2), dA(2), Acell_center, eps

    eps = 1.d-12
    do j = 1, blk%dim(2)
    do i = 1, blk%dim(1)
      dA_ij(1) = 0.5d0 * ((blk%Cross_Section(i  ,j  ) - blk%Cross_Section(i-1,j  )) + &
                          (blk%Cross_Section(i  ,j-1) - blk%Cross_Section(i-1,j-1)))
      dA_ij(2) = 0.5d0 * ((blk%Cross_Section(i  ,j  ) - blk%Cross_Section(i  ,j-1)) + &
                          (blk%Cross_Section(i-1,j  ) - blk%Cross_Section(i-1,j-1)))
      dA = Matmul ( dA_ij, blk%M(i,j)%c )
      Acell_center = 0.25d0 * ( blk%Cross_Section(i  ,j  ) + blk%Cross_Section(i-1,j  ) + &
                                blk%Cross_Section(i  ,j-1) + blk%Cross_Section(i-1,j-1) )
      blk%area_geom_term_x(i,j) = -dA(1) / max(Acell_center, eps)
      blk%area_geom_term_y(i,j) = -dA(2) / max(Acell_center, eps)
    end do
    end do
  end subroutine Precompute_Area_Geometry


  subroutine Read_Block_Q2D_Area ( domain )
    ! Reads file INPUT/block<N>_area.dat for each block; populates blk%Cross_Section.
    ! Sets blk%area_variation_active=.true. only if the file exists and parsing is OK
    ! (graceful degradation for blocks without a file: warning + cycle).
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer  :: i, j, b, ios, unitfile, ni, nj
    character(len=100)    :: filename
    real(R8), allocatable :: A(:,:)
    logical               :: lexist

    do b = 1, domain%nb
      ni = domain%blk(b)%dim(1)
      nj = domain%blk(b)%dim(2)

      write(filename, '(A,I0,A)') 'INPUT/block', b, '_area.dat'

      inquire(file=filename, exist=lexist)
      if (.not. lexist) then
        write(*,*) ' [Area] Warning: file not found: ', trim(filename), &
                   ' -> disabling area variation on the block'
        cycle
      end if

      open(newunit=unitfile, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
        write(*,*) ' [Area] Warning: cannot open ', trim(filename), &
                   ' (ios=', ios, ') -> disabling area variation on the block'
        cycle
      end if

      domain%blk(b)%area_variation_active = .true.
      allocate(A(0:ni,0:nj))
      read (unitfile, *, iostat=ios) ((A(i,j), i=0, ni), j=0, nj)
      close (unitfile)

      if (ios /= 0) then
        write(*,*) ' [Area] Warning: read error ', trim(filename), &
                   ' (ios=', ios, ') -> disabling area variation on the block'
        domain%blk(b)%area_variation_active = .false.
        deallocate(A)
        cycle
      end if

      domain%blk(b)%Cross_Section = A
      deallocate(A)
    end do
  end subroutine Read_Block_Q2D_Area


  subroutine Source_Area_Eul_Blk ( Prim, Res, sur, n, geometric_term_x, geometric_term_y )
    ! Area-variation Eulerian source terms: V·G accumulated over all equations.
    use Q2D_Global_m,            only: nprim, nsc, nu, nv, np, nt, gc
    use FLINT_Lib_Thermodynamic, only: Co_rotot_Rtot, E0
    implicit none
    integer, intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc),        intent(in)    :: sur
    real(R8), dimension(1:n(1),1:n(2)),                     intent(in)    :: geometric_term_x
    real(R8), dimension(1:n(1),1:n(2)),                     intent(in)    :: geometric_term_y
    ! Local
    integer  :: i, j
    real(R8) :: rho, Rgas, mod_v, Gx, Gy, V_dot_G
    real(R8) :: Source(nprim)

    !$omp do collapse(2) private(rho, Rgas, mod_v, Source, i, j, Gx, Gy, V_dot_G)
    do j = 1, n(2)
    do i = 1, n(1)

      call Co_rotot_Rtot ( Prim(1:nsc,i,j), rho, Rgas )
      mod_v = norm2(Prim(nu:nv,i,j))

      Gx      = geometric_term_x(i,j)
      Gy      = geometric_term_y(i,j)
      V_dot_G = Prim(nu,i,j)*Gx + Prim(nv,i,j)*Gy

      Source(1:nsc) = Prim(1:nsc,i,j) * V_dot_G
      Source(nu)    = rho * Prim(nu,i,j) * V_dot_G
      Source(nv)    = rho * Prim(nv,i,j) * V_dot_G
      Source(np)    = (rho * E0(Prim(np,i,j),Prim(1:nsc,i,j),mod_v) + Prim(np,i,j)) * V_dot_G

      ! Area-variation residual +rho*phi*V_dot_G for ALL post-pressure indices
      ! (soot Nvol, soot rho*Y_s nsoot=2, passive scalars, RANS).
      ! Verified by the term_by_term_sympy audit (LEMMA 5).
      if (nprim > np) then
        Source(np+1:nprim) = Prim(np+1:nprim,i,j) * V_dot_G
      end if

      Res(:,i,j) = Res(:,i,j) - Source * sur(i,j)
    end do; end do
  end subroutine Source_Area_Eul_Blk


  subroutine Source_Area_Visc_Blk ( Prim, Res, M, yn, sur, n, geometric_term_x, geometric_term_y, node )
    ! Area-variation viscous source terms: canonical Q2D Definition A formulation
    ! (S_33 = -V·G), wired via Q2D_Lib_Fluid::Stress_Tensor +
    ! Q2D_Lib_RANS::{Eddy_Viscosity, RANS_Diffusive_SourceQ2D}.
    !
    ! Theory references (consolidated canonical formulation = latest chronological):
    !   - teoria_finale_Q2D.md v6.1 §9.1, §13.7
    !   - teoria_Q2D_integrazione_v7.md §N16 "Final consolidated Q2D URANS equations"
    !
    ! Equations v7 §N16 (Anderson-consistent):
    !   x-mom:    +G_x(τ_33 − τ_11) − G_y·τ_12
    !   y-mom:    −G_x·τ_12 + G_y(τ_33 − τ_22)
    !   energy:   −V·(τ·G) + G·q + Σ_α h_α(G·j_α^*)
    !   species:  +G·j_α^*    [Σ j_α^* = 0 exact via Hirschfelder closure]
    !   RANS:     −D_φ(G·∇φ)
    !
    ! Code convention (Anderson-consistent):
    !   - DmGradYi(α) = D_α·∂Y_α/∂x_i  (positive of the gradient)
    !   - j_α = −ρ·DmGradYi (Fick, sign applied at accumulation time)
    !   - V_c = Σ_α DmGradYi: if Dm uniform over α → V_c = 0 → fix dormant.
    use Q2D_Base_Types_m,        only: Q2D_tensor_2D_type, Q2D_vector_2D_type
    use Q2D_Global_m,            only: nprim, nsc, nu, nv, np, nt, nc, nsoot, gc, model, nRANS
    use Q2D_Config_Types_m,      only: obj_rans, obj_soot
    use Q2D_Lib_Fluid,           only: Stress_Tensor
    use Q2D_Lib_RANS,            only: Eddy_Viscosity, RANS_Diffusive_SourceQ2D
    use Q2D_Mod_Soot,            only: i_CS
    use FLINT_Lib_Thermodynamic, only: Co_rotot_Rtot, f_cp_expr, co_k_mi_lam_Wilke_expr, &
                                        f_tabT_expr, h_tab
    implicit none
    integer, intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc),        intent(in)    :: sur
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc),        intent(in)    :: yn
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: M
    real(R8), dimension(1:n(1),1:n(2)), intent(in) :: geometric_term_x
    real(R8), dimension(1:n(1),1:n(2)), intent(in) :: geometric_term_y
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)), intent(in) :: node
    integer  :: i, j, s, T_i, Tint(2)
    real(R8) :: rho, rho1, rho2, rho3, rho4, Rgas, Rgas1, Rgas2, Rgas3, Rgas4
    real(R8) :: Source(nprim)
    real(R8) :: Stress_t(3,3), Gradient(1:nprim,2)
    real(R8) :: T, T1, T2, T3, T4, Tdiff, cp, mil, kl, mie, Dm(nsc), kappa, h_alpha, h_mix
    real(R8) :: DiffHFlux_x, DiffHFlux_y, DmGradYi_x(nsc), DmGradYi_y(nsc), V_c_x, V_c_y, qx, qy
    real(R8) :: T11, T12, T22, T33, Gx, Gy, V_dot_G, r_ave
    real(R8) :: rans_source(nRANS,2)

    !$omp do collapse (2) private ( rho, rho1, rho2, rho3, rho4, Rgas, Rgas1, Rgas2, Rgas3, Rgas4, Source, i, j, s, Gx, Gy, V_dot_G, r_ave ), &
    !$omp private ( T, T1, T2, T3, T4, T_i, Tdiff, Tint, cp, mil, kl, Gradient, mie, Dm, kappa, Stress_t, h_alpha, h_mix ), &
    !$omp private ( DiffHFlux_x, DiffHFlux_y, DmGradYi_x, DmGradYi_y, V_c_x, V_c_y, qx, qy, T11, T12, T22, T33, rans_source )
    do j = 1, n(2)
    do i = 1, n(1)

      ! Pre-zero the full post-pressure band: nu/nv/np/1:nsc are always assigned, but
      ! the passive band nps:nt-1 and the soot nc band (nsoot==2) are written only under
      ! guards. Mirror the Eulerian sibling's `if(nprim>np) Source(np+1:nprim)=...` coverage.
      Source = 0.d0

      Gx = geometric_term_x(i,j)
      Gy = geometric_term_y(i,j)
      ! Cell-center radius (axi-coherence; in pure area-variation r_ave is unused
      ! by Stress_Tensor_Std/Strain_Tensor_Std but propagated for consistency).
      r_ave = 0.25d0 * ( node(i  ,j  )%c(2) + node(i-1,j  )%c(2) &
                       + node(i  ,j-1)%c(2) + node(i-1,j-1)%c(2) )

      call Co_rotot_Rtot ( Prim(1:nsc,i,j), rho, Rgas )

      call Co_rotot_Rtot ( Prim(1:nsc,i+1,j), rho1, Rgas1 )
      call Co_rotot_Rtot ( Prim(1:nsc,i-1,j), rho2, Rgas2 )
      call Co_rotot_Rtot ( Prim(1:nsc,i,j+1), rho3, Rgas3 )
      call Co_rotot_Rtot ( Prim(1:nsc,i,j-1), rho4, Rgas4 )

      V_dot_G = Prim(nu,i,j)*Gx + Prim(nv,i,j)*Gy

      T = Prim(np,i,j) / ( rho * Rgas )

      T1 = Prim(np,i+1,j) / ( rho1 * Rgas1 )
      T2 = Prim(np,i-1,j) / ( rho2 * Rgas2 )
      T3 = Prim(np,i,j+1) / ( rho3 * Rgas3 )
      T4 = Prim(np,i,j-1) / ( rho4 * Rgas4 )

      T_i     = int(T)
      Tdiff   = T - T_i
      Tint(1) = T_i
      Tint(2) = T_i + 1
      cp = f_cp_expr ( Prim(1:nsc,i,j), Tint, Tdiff, rho )
      call co_k_mi_lam_Wilke_expr ( Prim(1:nsc,i,j), rho, Tint, Tdiff, mil, kl )

      Gradient(1:nsc,1) = ( Prim(1:nsc,i+1,j)/rho1 - Prim(1:nsc,i-1,j)/rho2 ) * 0.5d0
      Gradient(1:nsc,2) = ( Prim(1:nsc,i,j+1)/rho3 - Prim(1:nsc,i,j-1)/rho4 ) * 0.5d0

      Gradient(nu:nv,1) = ( Prim(nu:nv,i+1,j) - Prim(nu:nv,i-1,j) ) * 0.5d0
      Gradient(nu:nv,2) = ( Prim(nu:nv,i,j+1) - Prim(nu:nv,i,j-1) ) * 0.5d0

      Gradient(np,1) = (T1 - T2) * 0.5d0
      Gradient(np,2) = (T3 - T4) * 0.5d0

      if (model==2) then
        Gradient(nt:nprim,1) = ( Prim(nt:nprim,i+1,j)/rho1 - Prim(nt:nprim,i-1,j)/rho2 ) * 0.5d0
        Gradient(nt:nprim,2) = ( Prim(nt:nprim,i,j+1)/rho3 - Prim(nt:nprim,i,j-1)/rho4 ) * 0.5d0
      end if

      Gradient = Matmul ( Gradient, M(i,j) % c )

      ! Eddy_Viscosity (Q2D Def A): V_dot_G + axi-coherent vr/r_ave propagated for
      ! consistency with Stress_Tensor below. For pure area_var (axi=false) vr/r_ave
      ! are not read by Strain_Tensor_Std; passed to avoid fragile coupling.
      mie = 0d0
      if (model==2) then
        call Eddy_Viscosity ( mut=mie, rans_variables=Prim(nt:nprim,i,j), &
                              mul=mil, rho=rho, vel_gradient=Gradient(nu:nv,:), &
                              vr=Prim(nv,i,j), r_ave=r_ave, walldist=yn(i,j), V_dot_G=V_dot_G )
      end if

      ! Sc/Sct scalar (uniform Dm) -> V_c = D*grad(1) = 0 identically (fix dormant).
      Dm (1:nsc) = ( mil/obj_rans%Sc + mie/obj_rans%Sct ) / rho
      kappa = kl + mie*cp/obj_rans%Prt

      ! Stress tensor cell-center with Def A V_dot_G.
      Stress_t = Stress_Tensor ( Gradient(nu:nv,:), mil, mie, Prim(nt:,i,j), r_ave, Prim(nv,i,j), V_dot_G=V_dot_G )
      T11 = Stress_t(1,1); T22 = Stress_t(2,2); T12 = Stress_t(1,2); T33 = Stress_t(3,3)

      qx = -kappa * Gradient(np,1)
      qy = -kappa * Gradient(np,2)

      DiffHFlux_x = 0.d0
      DiffHFlux_y = 0.d0
      h_mix       = 0.d0
      do s = 1, nsc
        DmGradYi_x(s) = Dm(s) * Gradient(s,1)
        DmGradYi_y(s) = Dm(s) * Gradient(s,2)
        h_alpha       = f_tabT_expr ( s, h_tab, Tint, Tdiff )
        DiffHFlux_x   = DiffHFlux_x + DmGradYi_x(s) * h_alpha
        DiffHFlux_y   = DiffHFlux_y + DmGradYi_y(s) * h_alpha
        h_mix         = h_mix + ( Prim(s,i,j) / rho ) * h_alpha
      enddo
      V_c_x = sum(DmGradYi_x(1:nsc))
      V_c_y = sum(DmGradYi_y(1:nsc))

      ! Geometric viscous source terms (v7 §N16 Anderson-consistent).
      Source(nu) = Gx*(T33 - T11) - Gy*T12
      Source(nv) = -Gx*T12 + Gy*(T33 - T22)
      Source(np) = -Prim(nu,i,j)*(T11*Gx + T12*Gy) - Prim(nv,i,j)*(T12*Gx + T22*Gy)   &
                 +  Gx*qx + Gy*qy                                                      &
                 - (DiffHFlux_x - h_mix*V_c_x) * rho * Gx                              &
                 - (DiffHFlux_y - h_mix*V_c_y) * rho * Gy
      Source(1:nsc) = -( rho*DmGradYi_x(1:nsc) - V_c_x*Prim(1:nsc,i,j) ) * Gx        &
                     -( rho*DmGradYi_y(1:nsc) - V_c_y*Prim(1:nsc,i,j) ) * Gy

      ! Soot diffusive flux area-variation G-projection (retroactive).
      ! Verified by the term_by_term_sympy audit (LEMMA 3/6, 7).
      !
      ! Mod_Soot propagates the complete flux to Nvol via assignment:
      !   flux(nc) = flux(i_CS) / (rho * Y_s / N) = flux(i_CS) * Prim(nc)/Prim(i_CS)
      ! where flux(i_CS) post-Hirschfelder + Thermo = (Fick + Hirschfelder + Thermo).
      !
      ! Q2D sign convention area-var: Source(α) chemistry =
      !   -(rho·DmGradYi_x - V_c_x·Prim)·Gx - (rho·DmGradYi_y - V_c_y·Prim)·Gy   (Fick+Hirsch)
      !   -0.55·(mil+mie)·Prim/(rho·T)·(Gx·∂_x T + Gy·∂_y T)                      (Thermo)
      !
      ! For Nvol propagated, scaling Prim(nc)/Prim(i_CS):
      ! - Fick: -rho·(DmGradYi_x·Gx + DmGradYi_y·Gy) · Prim(nc)/Prim(i_CS)        (Y_s>0 guard)
      ! - Hirschfelder: +V_c_x·Prim(nc)·Gx + V_c_y·Prim(nc)·Gy                    (always safe)
      ! - Thermo: -0.55·(mil+mie)·Prim(nc)/(rho·T)·(Gx·∂_x T + Gy·∂_y T)         (always safe)
      !
      ! Skip nsoot=2 (i_CS=-1): soot thermophoretic diffusion is defined only for
      ! nsoot=1 (it is keyed on the soot species index i_CS and propagated to Nvol).
      ! The nsoot=2 heavy-gas model has no soot diffusive face-flux in ANY regime, so
      ! there is no area-projection of it to add here. (The inviscid area-variation
      ! residual for the nsoot=2 soot eqs nc/nc+1 is applied via Source(np+1:nprim)
      ! in Source_Area_Eul_Blk.)
      if (obj_soot%enabled .and. i_CS > 0) then
        ! Y_s thermo + N thermo propagated
        Source(i_CS) = Source(i_CS) - 0.55_R8*(mil + mie)*Prim(i_CS,i,j)/(rho*T) * &
                       (Gx*Gradient(np,1) + Gy*Gradient(np,2))
        Source(nc)   = Source(nc)   - 0.55_R8*(mil + mie)*Prim(nc,  i,j)/(rho*T) * &
                       (Gx*Gradient(np,1) + Gy*Gradient(np,2))
        ! N Fick + Hirschfelder propagated (extension of the area-var residual)
        Source(nc) = Source(nc) + V_c_x * Prim(nc,i,j) * Gx + V_c_y * Prim(nc,i,j) * Gy   ! Hirschfelder
        if (Prim(i_CS,i,j) > 1.0e-10_R8 * rho) then
          Source(nc) = Source(nc) - rho * Prim(nc,i,j) / Prim(i_CS,i,j) * &
                                    (DmGradYi_x(i_CS) * Gx + DmGradYi_y(i_CS) * Gy)
        end if
      end if

      if (model==2) then
        call RANS_Diffusive_SourceQ2D ( source=rans_source,                      &
                                        rans_variables=Prim(nt:nprim,i,j),       &
                                        mul=mil, mut=mie,                        &
                                        vel_gradient=Gradient(nu:nv,:),          &
                                        rans_gradient=Gradient(nt:nprim,:),      &
                                        rho=rho, dist=yn(i,j) )
        Source(nt:nprim) = -(rans_source(:,1)*Gx + rans_source(:,2)*Gy)
      end if

      Res(:,i,j) = Res(:,i,j) - Source * sur(i,j)

    enddo ; enddo

  end subroutine Source_Area_Visc_Blk

end module Q2D_Lib_Source_Areavariation
