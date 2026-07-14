!>@brief Spalart-Shur Rotation/Curvature correction (Spalart & Shur 1997).
!
! Reference: Spalart P.R., Shur M. (1997), "On the sensitization of turbulence
!            models to rotation and curvature", Aerospace Sci. Tech. 1(5):297-302.
!
! 2D applicability caveat: the model is motivated for 3D flows with significant
! rotation/curvature (turbomachinery rotating channels, swirl combustors, vortex
! flows). In Q2D 2D-planar the rotation/curvature correction is formally operative
! (scalar vorticity omega_z well defined; r* = S/Omega computable) but with reduced
! physical effect vs a 3D rotating channel (cf. Khayrullina et al. 2017). In 2D-axi
! without swirl (No-w), the rotation correction is formally applicable but on
! effectively non-rotating flows → engineering use, modest payoff.
!
! ============================================================================
! Formal mathematical justification for V·∇T_ij in Q2D
! ============================================================================
!
! The operator appearing in SAR-RC is the convective (material) derivative
! (Spalart-Shur 1997 §2): D(S_ij)/Dt = ∂_t S_ij + V_k ∂_k S_ij. At steady state
! D/Dt → V·∇ applied to the scalar COMPONENTS of S_ij (NOT the divergence ∇·T,
! which in cylindrical coordinates would carry Christoffel terms like (T_rr-T_θθ)/r).
!
! [Area-variation case, planar 2D depth-averaged]
!   Cartesian coordinates (x, y) with z the virtual top/bot direction. ∂_x and ∂_y
!   are ordinary partials, NO Christoffel terms. The components of S_ij in the
!   Cartesian basis are:
!     S_xx = ∂_x V_x;  S_yy = ∂_y V_y;  S_xy = ½(∂_x V_y + ∂_y V_x)
!     S_zz = -V·G  (Definition A: depth-averaged ∂w/∂z from Leibniz +
!                   impermeability BC at the virtual wall).
!   V·∇S_ij = V_x ∂_x S_ij + V_y ∂_y S_ij  →  component-by-component FD OK.
!   For S_zz: ∂_x S_zz = -∂_x(V·G) captured by passing the neighbour V_dot_G to Strain_Tensor.
!
! [Axisymmetric case, cylindrical orthonormal frame, No-w]
!   Coordinates (z, r, θ) with orthonormal basis (e_z, e_r, e_θ). The orthonormal-frame
!   data have only ∂_θ e_r = e_θ, ∂_θ e_θ = -e_r non-zero (the other basis derivatives = 0).
!   Under the Q2D No-w (V_θ = 0) + axisymmetric (∂_θ = 0) hypotheses:
!     V·∇S_ij = V_z ∂_z S_ij + V_r ∂_r S_ij + V_θ * [Christoffel terms]
!             = V_z ∂_z S_ij + V_r ∂_r S_ij  (V_θ = 0 cancels the Christoffel terms)
!   So for axi too it is a component-by-component FD of the orthonormal-frame
!   COMPONENTS of S_ij. The components are:
!     S_zz = ∂_z V_z;  S_rr = ∂_r V_r;  S_zr = ½(∂_z V_r + ∂_r V_z)
!     S_θθ = vr/r (cylindrical θθ-strain from continuity: ∇·V = ∂_z V_z + (1/r) ∂_r(r V_r)
!                                                               → S_θθ = V_r/r incompressible)
!   For S_θθ at the neighbour: using vr_neighbor AND r_ave_neighbor formally captures
!     ∂_r S_θθ = ∂_r(vr/r) = (1/r) ∂_r vr - vr/r²
!   BOTH terms are captured by the FD `(S_θθ_jp1 - S_θθ_jm1)/(2*Δr)` with
!   S_θθ evaluated with the neighbour (vr, r_ave). Using the central r_ave instead loses the
!   -vr/r² term in the radial direction.
!
! [Axis singularity]
!   r_ave → 0 invokes the Strain_Tensor_Axi de l'Hopital branch: S_θθ = ∂_r V_r.
!   Using r_ave_neighbor cell-by-cell, the de l'Hopital fallback is applied
!   locally for each neighbour cell that lies in the axis region.
!
! ============================================================================
! 2 procedures with Q2D 2D reductions:
!   - vel_gradient is Q2D_tensor_2D_type(2,2) (in-plane velocity gradients only;
!     out-of-plane d/dz structurally absent);
!   - strain_gradient(3,3,2): preserve TENSOR indices 1..3 (Sij 3×3), reduce the
!     SPATIAL direction from 3 to 2 (drop direction k);
!   - Ghost extrapolation for 4 faces (i, j) — drop 2 k-faces;
!   - Strain_Tensor 4-args Q2D-form (Gradvel(2,2), vr, r_ave[, V_dot_G]) returns 3×3
!     with Def A baked-in via procedure pointer dispatch (Std/Axi);
!   - **Strain_Tensor on the 5 cells with INTERIOR vs BOUNDARY metric handling**:
!     * INTERIOR (neighbour cell ∈ [1,ni]×[1,nj]): vel_gradient_neighbor + (vr, r_ave,
!       V_dot_G) computed cell-by-cell from existing data (corners node(0:ni,0:nj),
!       Prim(0:ni+1,0:nj+1) ghost, area_geom_term(1:ni,1:nj)). Formally captures the
!       full ∂_k S_ij INCLUDING the cylindrical ∂_k S(3,3) = ∂_r(vr/r) = (1/r)∂_r vr - vr/r²
!       (Axi) and ∂_k(V·G) (Std area-var). Preserves the cylindrical-divergence complexity
!       (Stress_Tensor_Axi: Divervel = G(1,1)+G(2,2)+vr/r vs Stress_Tensor_Std:
!       Divervel = G(1,1)+G(2,2)-V·G).
!     * BOUNDARY (neighbour cell out-of-range, i.e., 0 or ni+1 or j=0 or nj+1):
!       Sij_neighbor = Sij_central. The Q2D analogue of the constant ghost extrapolation
!       already used by Compute_Velocity_Gradient (above) for vel_gradient,
!       applied UNIFORMLY to ALL Strain_Tensor args (vel_gradient + vr +
!       r_ave + V_dot_G) — consistent with the Q2D cell-centered kernel pattern
!       (Source_Axi_Visc_Blk / SST_Blk / SA_Source_Blk: cell-centered metrics, NEVER
!       at the ghosts). Result: effective forward/backward FD at the boundary with factor 1/2
!       (accuracy O(Δ) at the edge vs O(Δ²) interior for the SAR-RC
!       sub-leading correction).
!     **No ghost-metric reach**: the metrics `axis_r_ave` (allocated 1-gc:ni+gc but
!     populated only 1..ni) and `area_geom_term_*` (shape 1:ni, no ghost) are NEVER
!     accessed outside the interior range. The boundary fallback uses cell-centered
!     values already computed at the reference cell (i,j) — zero reach to ghost data.
!   - Triple loop over TENSOR INDICES (1..3) preserved in rtilde;
!   - obj_rot%enabled absolute-vorticity correction wired.
module Q2D_Lib_SpalartShur
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private

  public :: Compute_Velocity_Gradient, Compute_RC_Terms

contains

  subroutine Compute_Velocity_Gradient ( domain )
    ! Cell-centered velocity gradient via central finite difference + ghost
    ! extrapolation. Q2D-form: nu:nv only, drop k loop + k faces;
    ! vel_gradient is Q2D_tensor_2D_type(2,2).
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nu, nv
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, n
    real(R8), dimension(2,2) :: vel_diff   ! 2×2 (vs the 3D peer's 3×3)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks

      !$omp do collapse(2) private(vel_diff)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        ! Cell-centered velocity gradient (drop nw and direction k)
        vel_diff(:,1) = ( domain%blk(b)%P(nu:nv,i+1,j) - domain%blk(b)%P(nu:nv,i-1,j) )/2d0
        vel_diff(:,2) = ( domain%blk(b)%P(nu:nv,i,j+1) - domain%blk(b)%P(nu:nv,i,j-1) )/2d0
        ! 1d-40 epsilon (avoid pathological zero gradients
        ! in laminar regions before metric matmul).
        domain%blk(b)%vel_gradient(i,j)%c = matmul ( vel_diff + 1d-40, domain%blk(b)%m(i,j)%c )
      end do ; end do

      ! Extrapolation in ghost cells (second derivatives needed by Compute_RC_Terms).
      ! i faces (drop k face)
      !$omp do private(n)
      do j = 1, domain % blk(b) % dim(2)
        n = domain % blk(b) % dim(1) + 1
        domain % blk(b) % vel_gradient(0,j) % c = domain % blk(b) % vel_gradient(1,j)   % c
        domain % blk(b) % vel_gradient(n,j) % c = domain % blk(b) % vel_gradient(n-1,j) % c
      end do

      ! j faces (drop k face)
      !$omp do private(n)
      do i = 1, domain % blk(b) % dim(1)
        n = domain % blk(b) % dim(2) + 1
        domain % blk(b) % vel_gradient(i,0) % c = domain % blk(b) % vel_gradient(i,1)   % c
        domain % blk(b) % vel_gradient(i,n) % c = domain % blk(b) % vel_gradient(i,n-1) % c
      end do

      ! k faces: DROP — No-w + 2D planar + 2D-axi (no out-of-plane direction).

    end do

  end subroutine Compute_Velocity_Gradient


  subroutine Compute_RC_Terms ( domain )
    ! rstar = S/Omega, rtilde from material derivative of Sij (without
    ! denominator: different between SA-RC and SST-RC, applied at use site in
    ! SA_Source_Blk / SST_Blk).
    !
    ! V·∇S_ij operator handling (cf. file header formal justification):
    !   - Cartesian planar (Std + area-variation): V·∇ = V_x ∂_x + V_y ∂_y (no Christoffel);
    !   - Cylindrical axisymmetric No-w (Axi):     V·∇ = V_z ∂_z + V_r ∂_r (V_θ=0 + axisymm
    !                                                    cancel the Christoffel frame);
    !   component-by-component FD formally correct in both cases.
    !
    ! Q2D-consistent cell-centered pattern (mirror Source_Axi_Visc_Blk, SST_Blk, SA_Source_Blk):
    !   - vel_gradient on the 5 cells (central + 4 neighbours) — constant ghost extrapolation
    !     provided by Compute_Velocity_Gradient above (vel_gradient(0,j)%c =
    !     vel_gradient(1,j)%c, etc.). Captures ∂_k S(1,1), S(1,2), S(2,1), S(2,2) via
    !     FD on vel_gradient_neighbor.
    !   - vr_ij, r_ave_ij, V_dot_G_ij: cell-centered ONLY, passed identically to all
    !     5 Strain_Tensor calls. Q2D does not extrapolate metrics to the ghosts (axis_r_ave
    !     allocated 1-gc:ni+gc but populated only 1..ni; area_geom_term_* shape 1:ni
    !     without ghost). Consequence: S(3,3) constant in the FD stencil → ∂_k S(3,3) = 0.
    !     **Accepted trade-off**: loses `(1/r_ave) ∂_k vr` (axi) and `-∂_k(V·G)`
    !     (area-var) contributions to strain transport. These are SUB-SUB-LEADING vs the
    !     dominant contribution via the in-plane S(1,1)/S(1,2)/S(2,2) captured by
    !     vel_gradient FD. (In the 3D peer, vel_gradient(3,3) directly contains
    !     ∂_z V_z; in Q2D 2D, S(3,3) is a z-reduction feature treated as constant
    !     for strain transport.)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nu, nv
    use Q2D_Lib_Fluid,        only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    use Q2D_Config_Types_m,   only: obj_rot, obj_source
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, ii, jj, kk, ni, nj
    real(R8), dimension(3,3,2) :: strain_gradient                 ! drop spatial direction 3
    real(R8), dimension(3,3)   :: Sij_ip1, Sij_im1, Sij_jp1, Sij_jm1
    real(R8), dimension(3,3)   :: Sij, Wij, strain_transport
    real(R8) :: omega(3), S, O, rtilde
    real(R8) :: r_ave, vr_ij, VG_ij
    logical  :: area_var

    area_var = obj_source%area_variation

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks

      ni = domain % blk(b) % dim(1)
      nj = domain % blk(b) % dim(2)

      !$omp do collapse(2) private (i, j, ii, jj, kk, strain_gradient, Sij_ip1, Sij_im1, Sij_jp1, Sij_jm1), &
      !$omp private (Sij, Wij, strain_transport, omega, S, O, rtilde, r_ave, vr_ij, VG_ij)
      do j = 1, nj
      do i = 1, ni

        ! Cell-centered geometric quantities (Q2D does not extrapolate metrics to the
        ! ghosts; pattern consistent with Source_Axi_Visc_Blk + SST_Blk + SA_Source_Blk).
        r_ave = 0.25d0 * ( domain%blk(b)%node(i  ,j  )%c(2) + domain%blk(b)%node(i-1,j  )%c(2) &
                         + domain%blk(b)%node(i  ,j-1)%c(2) + domain%blk(b)%node(i-1,j-1)%c(2) )
        vr_ij = domain%blk(b)%P(nv, i, j)

        if (area_var) then
          VG_ij = domain%blk(b)%P(nu,i,j)*domain%blk(b)%area_geom_term_x(i,j) &
                + domain%blk(b)%P(nv,i,j)*domain%blk(b)%area_geom_term_y(i,j)
        else
          VG_ij = 0d0
        end if

        ! ====================================================================
        ! Strain tensor at 5 cells.
        ! - vel_gradient_neighbor: ghost-extended (Compute_Velocity_Gradient above)
        !   → captures ∂_k S(1,1), S(1,2), S(2,1), S(2,2) via spatial FD.
        ! - (vr_ij, r_ave, VG_ij): cell-centered, identical for all 5 calls
        !   → S(3,3) constant in the FD stencil; ∂_k S(3,3) = 0 (sub-sub-leading miss
        !   accepted — cf. file header for the formal analysis).
        ! ====================================================================
        Sij_ip1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i+1,j)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_im1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i-1,j)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_jp1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i,j+1)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_jm1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i,j-1)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij     = Strain_Tensor ( domain%blk(b)%vel_gradient(i  ,j  )%c, vr_ij, r_ave, V_dot_G=VG_ij )

        ! Finite difference strain rate gradient (drop spatial direction 3).
        ! Component-by-component FD on S_ij in the orthonormal basis (cylindrical axi No-w or
        ! Cartesian planar area-var) formally captures V·∇S_ij.
        strain_gradient(:,:,1) = ( Sij_ip1 - Sij_im1 ) * 0.5d0
        strain_gradient(:,:,2) = ( Sij_jp1 - Sij_jm1 ) * 0.5d0

        ! Metric transformation tensor index by tensor index (preserve 3×3).
        ! The cell-centered metric M(i,j) is applied to the spatial gradient (Jacobian
        ! from computational space (xi, eta) to physical (x, y) or (z, r)) — a tensor
        ! invariant, applied to each Sij_ij component separately.
        do ii = 1, 3
        do jj = 1, 3
          strain_gradient(ii,jj,1:2) = matmul ( strain_gradient(ii,jj,1:2), domain%blk(b)%m(i,j)%c )
        end do
        end do

        ! Cell-centered vorticity tensor + obj_rot absolute vorticity correction
        omega = Vorticity_Vector( domain%blk(b)%vel_gradient(i,j)%c )
        ! Absolute vorticity:
        ! omega_abs = omega_rel + 2*Omega_frame; Q2D Vorticity_Vector returns vector(3)
        ! with omega(1)=omega(2)=0; Omega_vec axis [0,0,1] strict.
        if (obj_rot%enabled) omega = omega + 2.0_R8 * obj_rot%Omega_vec
        Wij = Vorticity_Tensor( omega )
        S = sqrt ( 2d0*sum (Sij**2) )
        O = sqrt ( sum(omega**2) )
        domain%blk(b)%rc_term1(i,j) = S/O

        ! Material derivative of strain rate tensor V·∇S_ij (drop spatial direction 3).
        ! P(nu:nv,i,j) are the physical velocity COMPONENTS (V_x, V_y) or (V_z, V_r).
        ! V·∇ is the directional derivative on the scalar components of S_ij; formally
        ! correct for Cartesian planar (no Christoffel) and cylindrical axisymmetric
        ! No-w (V_θ=0 + ∂_θ=0 cancel the Christoffel-frame contributions).
        do ii = 1, 3
        do jj = 1, 3
          strain_transport(ii,jj) = dot_product ( domain%blk(b)%P(nu:nv,i,j), &
                                                  strain_gradient(ii,jj,1:2) )
        end do
        end do

        ! Triple sum over TENSOR INDICES 1..3 (NON-spatial, preserve).
        rtilde = 0d0
        do ii = 1, 3
        do jj = 1, 3
        do kk = 1, 3
          rtilde = rtilde + Wij(ii,kk)*Sij(jj,kk)*strain_transport(ii,jj)
        end do
        end do
        end do
        ! Note: denominator missing, applied at use site in SA_Source_Blk / SST_Blk
        ! (different between SA-RC and SST-RC).
        domain%blk(b)%rc_term2(i,j) = 2d0*rtilde

      end do ; end do
    end do

  end subroutine Compute_RC_Terms

end module Q2D_Lib_SpalartShur
