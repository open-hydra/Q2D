module Q2D_Lib_Source_Axisymmetry
  ! Q2D axisymmetric source terms:
  !   Eulerian: S(nv) = -ρ·v_r²/r,  S(np) = -(ρ·E0+p)·v_r/r, etc.
  !   viscous: (τ_rr - τ_θθ)/r — realized via Q2D_Lib_Fluid::Stress_Tensor +
  !            Q2D_Lib_RANS::{Eddy_Viscosity, RANS_Diffusive_SourceQ2D}.
  !
  ! Symmetry axis (r → 0): de l'Hôpital fallback with eps=1.d-12.
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m, only: Q2D_block_type
  implicit none
  private

  public :: Precompute_Axisymmetric_Geometry
  public :: Source_Axi_Eul_Blk
  public :: Source_Axi_Visc_Blk

contains

  subroutine Precompute_Axisymmetric_Geometry ( blk )
    ! Precomputes r_ave (cell-center radius) and r_inv = 1/max(r_ave, eps).
    ! The axis_r_ave/axis_inv_r fields are allocated by Mod_Allocate_Data
    ! → no self-allocation here.
    implicit none
    type(Q2D_block_type), intent(inout) :: blk
    integer  :: i, j
    real(R8) :: r_ave, eps

    eps = 1.d-12
    do j = 1, blk%dim(2)
    do i = 1, blk%dim(1)
      r_ave = 0.25d0 * ( blk%node(i  ,j  )%c(2) + blk%node(i-1,j  )%c(2) + &
                         blk%node(i  ,j-1)%c(2) + blk%node(i-1,j-1)%c(2) )
      blk%axis_r_ave(i,j) = r_ave
      blk%axis_inv_r(i,j) = 1.d0 / max(r_ave, eps)
    end do
    end do
  end subroutine Precompute_Axisymmetric_Geometry


  subroutine Source_Axi_Eul_Blk ( Prim, Res, r_ave, r_inv, sur, n )
    ! Eulerian axisymmetric source terms: geometric terms in 1/r.
    use Q2D_Global_m,            only: nprim, nsc, nu, nv, np, nt, gc
    use FLINT_Lib_Thermodynamic, only: Co_rotot_Rtot, E0
    implicit none
    integer, intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc),        intent(in)    :: sur
    real(R8), dimension(1:n(1),1:n(2)),                     intent(in)    :: r_ave
    real(R8), dimension(1:n(1),1:n(2)),                     intent(in)    :: r_inv
    ! Local variables
    integer  :: i, j
    real(R8) :: rho, Rgas, mod_v, invr
    real(R8) :: Source(nprim)

    !$omp do collapse(2) private(rho, Rgas, mod_v, invr, Source, i, j)
    do j = 1, n(2)
    do i = 1, n(1)

      call Co_rotot_Rtot ( Prim(1:nsc,i,j), rho, Rgas )
      mod_v = norm2(Prim(nu:nv,i,j))
      invr  = r_inv(i,j)

      ! Source multiplied by rho: the SA equation is integrated in conservative form
      Source(1:nsc) = -Prim(1:nsc,i,j) * Prim(nv,i,j) * invr
      Source(nu)    = -rho * Prim(nu,i,j) * Prim(nv,i,j) * invr
      Source(nv)    = -rho * Prim(nv,i,j) * Prim(nv,i,j) * invr
      Source(np)    = -(rho * E0(Prim(np,i,j),Prim(1:nsc,i,j),mod_v) + Prim(np,i,j)) * Prim(nv,i,j) * invr

      ! Cylindrical residual -rho*phi*V_r/r for ALL indices beyond the pressure
      ! (soot Nvol, soot rho*Y_s for nsoot=2, passive scalars, RANS).
      if (nprim > np) then
        Source(np+1:nprim) = -Prim(nv,i,j) * Prim(np+1:nprim,i,j) * invr
      end if

      Res(:,i,j) = Res(:,i,j) - Source * sur(i,j)
    end do; end do
  end subroutine Source_Axi_Eul_Blk


  subroutine Source_Axi_Visc_Blk ( Prim, Res, M, yn, r_ave_arr, r_inv, sur, n )
    ! Viscous axisymmetric source terms: geometric diffusive terms
    ! (cylindrical in 1/r).
    !
    ! Structure of the viscous geometric source terms (cylindrical in 1/r):
    !   - Source(nu) = Stress_t(2,1)/r
    !   - Source(nv) = (Stress_t(2,2)-Stress_t(3,3))/r   ! tau_rr - tau_theta_theta
    !   - Source(np) = (-qr + V·tau + (DiffHFlux - h_mix*V_c)*rho)/r
    !   - symmetric Hirschfelder closure for the species (V_c dormant in the current
    !     uniform-Dm regime; it activates with species-dependent Dm via
    !     FLINT co_DS).
    !   - RANS_Diffusive_SourceQ2D for the turbulent variables.
    use Q2D_Base_Types_m,        only: Q2D_tensor_2D_type
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
    real(R8), dimension(1:n(1),1:n(2)), intent(in) :: r_ave_arr
    real(R8), dimension(1:n(1),1:n(2)), intent(in) :: r_inv
    integer  :: i, j, s, Tint(2)
    real(R8) :: rho, rho1, rho2, rho3, rho4, Rgas, Rgas1, Rgas2, Rgas3, Rgas4, r_ave, mod_v, invr
    real(R8) :: Source(nprim)
    real(R8) :: Stress_t(3,3), Gradient(1:nprim,2)
    real(R8) :: T, T1, T2, T3, T4, T_i, Tdiff, cp, mil, kl, mie, Dm(nsc), kappa, DiffHFlux, DmGradYi(nsc), qr
    real(R8) :: h_alpha, h_mix, V_c
    real(R8) :: rans_source(nRANS,2)

    !$omp do collapse (2) private ( rho, rho1, rho2, rho3, rho4, Rgas, Rgas1, Rgas2, Rgas3, Rgas4, r_ave, invr, mod_v, Source, i, j, s ), &
    !$omp private ( T, T1, T2, T3, T4, T_i, Tdiff, Tint, cp, mil, kl, Gradient, mie, Dm, kappa, Stress_t, DiffHFlux, DmGradYi, qr, rans_source, h_alpha, h_mix, V_c )
    do j = 1, n(2)
    do i = 1, n(1)

      ! Pre-zeroing of the whole band beyond the pressure: the passive band nps:nt-1
      ! and the soot band nc (nsoot==2) are written only under a guard; without this
      ! they would reach the read Res(:,i,j) = Res - Source*sur uninitialized.
      ! Coverage mirrors that of Source_Axi_Eul_Blk.
      Source = 0.d0

      call Co_rotot_Rtot ( Prim(1:nsc,i,j), rho, Rgas )

      call Co_rotot_Rtot ( Prim(1:nsc,i+1,j), rho1, Rgas1 )
      call Co_rotot_Rtot ( Prim(1:nsc,i-1,j), rho2, Rgas2 )
      call Co_rotot_Rtot ( Prim(1:nsc,i,j+1), rho3, Rgas3 )
      call Co_rotot_Rtot ( Prim(1:nsc,i,j-1), rho4, Rgas4 )

      mod_v = norm2(Prim(nu:nv,i,j))
      r_ave = r_ave_arr(i,j)
      invr = r_inv(i,j)
      T = Prim(np,i,j) / ( rho * Rgas )

      T1 = Prim(np,i+1,j) / ( rho1 * Rgas1 )
      T2 = Prim(np,i-1,j) / ( rho2 * Rgas2 )
      T3 = Prim(np,i,j+1) / ( rho3 * Rgas3 )
      T4 = Prim(np,i,j-1) / ( rho4 * Rgas4 )

      T_i     = int(T)
      Tdiff   = T - T_i
      Tint(1) = int(T_i)
      Tint(2) = Tint(1) + 1
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

      mie = 0d0
      if (model==2) then
        call Eddy_Viscosity ( mut=mie, rans_variables=Prim(nt:nprim,i,j), &
                              mul=mil, rho=rho, vel_gradient=Gradient(nu:nv,:), &
                              vr=Prim(nv,i,j), r_ave=r_ave, walldist=yn(i,j) )
      end if

      ! Scalar Sc/Sct (uniform Dm) -> V_c = sum Dm*grad(Y) = D*grad(1) = 0 identically.
      ! Underlying Hirschfelder correction dormant (it activates with species-dependent
      ! Dm).
      Dm (1:nsc) = ( mil/obj_rans%Sc + mie/obj_rans%Sct ) / rho
      kappa = kl + mie*cp/obj_rans%Prt
      Stress_t = Stress_Tensor ( Gradient(nu:nv,:), mil, mie, Prim(nt:,i,j), r_ave, Prim(nv,i,j) )

      DiffHFlux = 0.d0
      h_mix     = 0.d0
      do s = 1, nsc
        DmGradYi(s) = Dm(s) * Gradient(s,2)
        h_alpha     = f_tabT_expr ( s, h_tab, Tint, Tdiff )
        DiffHFlux   = DiffHFlux + DmGradYi(s) * h_alpha
        h_mix       = h_mix + ( Prim(s,i,j) / rho ) * h_alpha
      enddo
      V_c = sum(DmGradYi(1:nsc))

      qr = -kappa * Gradient(np,2)

      Source(1:nsc) = (rho * DmGradYi(1:nsc) - V_c * Prim(1:nsc,i,j)) * invr
      Source(nu) = Stress_t(2,1) * invr
      Source(nv) = (Stress_t(2,2)-Stress_t(3,3)) * invr  ! (tau_rr - tau_theta_theta)/r
      Source(np) = (-qr + Prim(nu,i,j) * Stress_t(1,2) + Prim(nv,i,j) * Stress_t(2,2) &
                  + (DiffHFlux - h_mix*V_c) * rho) * invr

      ! Cylindrical residual of the soot diffusive flux, axisymmetric mode.
      !
      ! Mod_Soot propagates the full flux to Nvol via the assignment:
      !   flux(nc) = flux(i_CS) / (rho * Y_s / N) = flux(i_CS) * Prim(nc)/Prim(i_CS)
      ! where flux(i_CS) after Hirschfelder + thermophoresis = (Fick + Hirschfelder + Thermo).
      !
      ! Cell-level cylindrical residual for Nvol =
      !   -j_total,Y_s,r/r · Prim(nc)/Prim(i_CS)
      !
      ! Decomposition for clarity:
      ! - Thermo (i_CS / propagated to N): explicit formula with (mil+mie)·Y_s/T·∂_r T·invr
      ! - Propagated Fick (N): rho·Dm(i_CS)·∂_r Y_s · N/Y_s · invr (requires guard Y_s>0)
      ! - Propagated Hirschfelder (N): -V_c · N · invr (Prim(i_CS) cancels out,
      !   always safe)
      !
      ! Skip nsoot=2 (i_CS=-1): soot thermophoretic diffusion is defined only
      ! for nsoot=1 (it is tied to the soot species index i_CS and propagated to Nvol).
      ! The heavy-gas model nsoot=2 has no face diffusive flux for the soot in
      ! ANY regime, hence there is no curvature correction of it to add
      ! here. (The inviscid cylindrical residual for the soot equations nc/nc+1 of nsoot=2
      ! is applied via Source(np+1:nprim) in Source_Axi_Eul_Blk.)
      if (obj_soot%enabled .and. i_CS > 0) then
        ! Thermophoresis on Y_s + thermophoresis propagated to N
        Source(i_CS) = Source(i_CS) + 0.55_R8*(mil + mie)*Prim(i_CS,i,j)/(rho*T) * Gradient(np,2) * invr
        Source(nc)   = Source(nc)   + 0.55_R8*(mil + mie)*Prim(nc,  i,j)/(rho*T) * Gradient(np,2) * invr
        ! Fick + Hirschfelder propagated to N (extension of the cylindrical residual)
        Source(nc)   = Source(nc)   - V_c * Prim(nc,i,j) * invr   ! propagated Hirschfelder (always safe)
        if (Prim(i_CS,i,j) > 1.0e-10_R8 * rho) then
          Source(nc) = Source(nc)   + rho * DmGradYi(i_CS) * Prim(nc,i,j) / Prim(i_CS,i,j) * invr
        end if
      end if

      if (model==2) then
        call RANS_Diffusive_SourceQ2D ( source=rans_source,                      &
                                        rans_variables=Prim(nt:nprim,i,j),       &
                                        mul=mil, mut=mie,                        &
                                        vel_gradient=Gradient(nu:nv,:),          &
                                        rans_gradient=Gradient(nt:nprim,:),      &
                                        rho=rho, dist=yn(i,j) )
        Source(nt:nprim) = rans_source(:,2) * invr
      end if

      Res(:,i,j) = Res(:,i,j) - Source * sur(i,j)

    enddo ; enddo

  end subroutine Source_Axi_Visc_Blk

end module Q2D_Lib_Source_Axisymmetry
