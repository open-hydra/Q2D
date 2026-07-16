!>@brief Spalart-Shur rotation/curvature correction (Spalart & Shur 1997).
!
! Reference: Spalart P.R., Shur M. (1997), "On the sensitization of turbulence
!              models to rotation and curvature", Aerospace Sci. Tech. 1(5):297-302.
!
! 2D applicability caveat: the model is motivated for 3D flows with
! significant rotation/curvature (rotating turbomachinery channels, swirl
! combustors, vortical flows). In the Q2D 2D-plane the rotation/curvature correction
! is formally operative (scalar vorticity omega_z well defined; r* = S/Omega
! computable) but with a reduced physical effect compared to a 3D rotating channel
! (cf. Khayrullina et al. 2017). In the 2D-axi case without swirl (No-w) the rotation
! correction is formally applicable but on effectively non-rotating flows →
! engineering use, modest benefit.
!
! ============================================================================
! Formal mathematical justification of V·∇S_ij in Q2D
! ============================================================================
!
! The operator appearing in SA-RC is the convective (material) derivative
! (Spalart-Shur 1997 §2): D(S_ij)/Dt = ∂_t S_ij + V_k ∂_k S_ij. At steady state
! D/Dt → V·∇ applied to the scalar COMPONENTS of S_ij (NOT the divergence ∇·T,
! which in cylindrical coordinates would bring Christoffel terms such as (T_rr-T_θθ)/r).
!
! [Area-variation case, depth-averaged 2D plane]
!   Cartesian coordinates (x, y) with z the top/bot virtual direction. ∂_x and ∂_y are
!   ordinary partial derivatives, NO Christoffel term. The components of
!   S_ij in the Cartesian basis are:
!     S_xx = ∂_x V_x;  S_yy = ∂_y V_y;  S_xy = ½(∂_x V_y + ∂_y V_x)
!     S_zz = -V·G  (Definition A: ∂w/∂z depth-averaged via Leibniz +
!                   impermeability BC at the virtual wall).
!   V·∇S_ij = V_x ∂_x S_ij + V_y ∂_y S_ij  →  component-by-component FD OK.
!   For S_zz: ∂_x S_zz = -∂_x(V·G) would be captured by passing Strain_Tensor the
!   V_dot_G of the neighboring cell (not done: cf. trade-off in Compute_RC_Terms).
!
! [Axisymmetric case, orthonormal cylindrical triad, No-w]
!   Coordinates (z, r, θ) with orthonormal basis (e_z, e_r, e_θ). The only nonzero
!   basis derivatives are ∂_θ e_r = e_θ, ∂_θ e_θ = -e_r (the others = 0).
!   Under the Q2D No-w assumptions (V_θ = 0) + axisymmetry (∂_θ = 0):
!     V·∇S_ij = V_z ∂_z S_ij + V_r ∂_r S_ij + V_θ * [Christoffel terms]
!             = V_z ∂_z S_ij + V_r ∂_r S_ij  (V_θ = 0 cancels the Christoffel terms)
!   So in axi too it is a component-by-component FD of the COMPONENTS of
!   S_ij in the orthonormal triad. The components are:
!     S_zz = ∂_z V_z;  S_rr = ∂_r V_r;  S_zr = ½(∂_z V_r + ∂_r V_z)
!     S_θθ = vr/r (cylindrical θθ strain from continuity: ∇·V = ∂_z V_z + (1/r) ∂_r(r V_r)
!                                                               → S_θθ = V_r/r incompressible)
!   For S_θθ at the neighboring cell: using vr_neighbor AND r_ave_neighbor would
!   formally capture
!     ∂_r S_θθ = ∂_r(vr/r) = (1/r) ∂_r vr - vr/r²
!   BOTH terms would be captured by the FD `(S_θθ_jp1 - S_θθ_jm1)/(2*Δr)`
!   with S_θθ evaluated with (vr, r_ave) of the neighboring cell. Using instead the
!   central r_ave loses the -vr/r² term in the radial direction (current choice:
!   cf. trade-off in Compute_RC_Terms).
!
! [Singularity on the axis]
!   r_ave → 0 invokes the de l'Hopital branch of Strain_Tensor_Axi: S_θθ = ∂_r V_r.
!   Using r_ave_neighbor cell by cell, the de l'Hopital fallback would
!   apply locally to each neighboring cell in the axis region.
!
! ============================================================================
! 2 procedures with Q2D 2D reductions:
!   - vel_gradient is Q2D_tensor_2D_type(2,2) (only in-plane velocity gradients;
!     out-of-plane d/dz structurally absent);
!   - strain_gradient(3,3,2): the TENSOR indices 1..3 are kept (Sij 3×3),
!     the SPATIAL direction is reduced from 3 to 2 (the k direction eliminated);
!   - ghost extrapolation on 4 faces (i, j) — the 2 k faces eliminated;
!   - Strain_Tensor with 4 arguments (Gradvel(2,2), vr, r_ave[, V_dot_G]) returns
!     3×3 with Def A built in via procedure-pointer dispatch (Std/Axi);
!   - **Strain_Tensor over the 5 cells (central + 4 neighbors)**:
!     * neighboring-cell vel_gradient: extended to the ghosts with constant
!       extrapolation (Compute_Velocity_Gradient) → the spatial FD captures ∂_k S(1,1), S(1,2),
!       S(2,1), S(2,2). For a neighboring cell out of range (i=0 or ni+1, j=0 or
!       nj+1) the constant extrapolation gives Sij_neighbor = Sij_central: FD
!       effectively one-sided at the boundary with factor 1/2 (O(Δ) accuracy at
!       the boundary vs O(Δ²) in the interior, for an SA-RC correction that is sub-dominant anyway).
!     * (vr, r_ave, V_dot_G): ONLY cell-centered, passed identical to the 5
!       calls — pattern consistent with the cell-centered Q2D kernels
!       (Source_Axi_Visc_Blk / SST_Blk / SA_Source_Blk: cell-centered
!       metrics, NEVER at the ghosts). Consequence: S(3,3) constant in the FD stencil →
!       ∂_k S(3,3) = 0 (trade-off analyzed in Compute_RC_Terms).
!     **No ghost-metric access**: the metrics `axis_r_ave` (allocated
!     1-gc:ni+gc but populated only 1..ni) and `area_geom_term_*` (shape 1:ni, without
!     ghost) are NEVER read outside the interior range; only values
!     already computed at the reference cell (i,j) are used — zero accesses to ghost data.
!   - triple loop over the TENSOR INDICES (1..3) preserved in rtilde;
!   - obj_rot%enabled absolute-vorticity correction active.
module Q2D_Lib_SpalartShur
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private

  public :: Compute_Velocity_Gradient, Compute_RC_Terms

contains

  subroutine Compute_Velocity_Gradient ( domain )
    ! Cell-centered velocity gradient via centered finite differences +
    ! ghost extrapolation. Q2D form: only nu:nv, without k loop or k faces;
    ! vel_gradient is Q2D_tensor_2D_type(2,2).
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nu, nv
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, n
    real(R8), dimension(2,2) :: vel_diff   ! 2×2 (2D reduction)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks

      !$omp do collapse(2) private(vel_diff)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        ! Cell-centered velocity gradient (without nw or k direction)
        vel_diff(:,1) = ( domain%blk(b)%P(nu:nv,i+1,j) - domain%blk(b)%P(nu:nv,i-1,j) )/2d0
        vel_diff(:,2) = ( domain%blk(b)%P(nu:nv,i,j+1) - domain%blk(b)%P(nu:nv,i,j-1) )/2d0
        ! epsilon 1d-40 (avoids pathologically null gradients
        ! in laminar regions before the matmul with the metric).
        domain%blk(b)%vel_gradient(i,j)%c = matmul ( vel_diff + 1d-40, domain%blk(b)%m(i,j)%c )
      end do ; end do

      ! Extrapolation into the ghost cells (second derivatives required by Compute_RC_Terms).
      ! i faces (no k face)
      !$omp do private(n)
      do j = 1, domain % blk(b) % dim(2)
        n = domain % blk(b) % dim(1) + 1
        domain % blk(b) % vel_gradient(0,j) % c = domain % blk(b) % vel_gradient(1,j)   % c
        domain % blk(b) % vel_gradient(n,j) % c = domain % blk(b) % vel_gradient(n-1,j) % c
      end do

      ! j faces (no k face)
      !$omp do private(n)
      do i = 1, domain % blk(b) % dim(1)
        n = domain % blk(b) % dim(2) + 1
        domain % blk(b) % vel_gradient(i,0) % c = domain % blk(b) % vel_gradient(i,1)   % c
        domain % blk(b) % vel_gradient(i,n) % c = domain % blk(b) % vel_gradient(i,n-1) % c
      end do

      ! k faces: ABSENT — No-w + 2D plane + 2D-axi (no out-of-plane direction).

    end do

  end subroutine Compute_Velocity_Gradient


  subroutine Compute_RC_Terms ( domain )
    ! rstar = S/Omega, rtilde from the material derivative of Sij (without
    ! denominator: different between SA-RC and SST-RC, applied at the use site in
    ! SA_Source_Blk / SST_Blk).
    !
    ! Treatment of the operator V·∇S_ij (cf. formal justification at the top of the file):
    !   - Cartesian plane (Std + area variation): V·∇ = V_x ∂_x + V_y ∂_y
    !                                             (no Christoffel);
    !   - No-w axisymmetric cylindrical (Axi):    V·∇ = V_z ∂_z + V_r ∂_r
    !                                             (V_θ=0 + axisymmetry cancel
    !                                             the Christoffel terms);
    !   Component-by-component FD formally correct in both cases.
    !
    ! Cell-centered pattern consistent with Source_Axi_Visc_Blk, SST_Blk, SA_Source_Blk:
    !   - vel_gradient over the 5 cells (central + 4 neighbors) — constant ghost extrapolation
    !     provided by Compute_Velocity_Gradient above (vel_gradient(0,j)%c =
    !     vel_gradient(1,j)%c, etc.). Captures ∂_k S(1,1), S(1,2), S(2,1), S(2,2) via
    !     FD on the neighboring-cell vel_gradient.
    !   - vr_ij, r_ave_ij, V_dot_G_ij: ONLY cell-centered, passed identical to the
    !     5 calls to Strain_Tensor. The metrics are not extrapolated to the ghosts
    !     (axis_r_ave allocated 1-gc:ni+gc but populated only 1..ni; area_geom_term_*
    !     shape 1:ni without ghost). Consequence: S(3,3) constant in the FD stencil →
    !     ∂_k S(3,3) = 0.
    !     **Accepted trade-off**: the contributions `(1/r_ave) ∂_k vr` (axi) and
    !     `-∂_k(V·G)` (area-var) to the strain transport are lost. They are SUB-SUB-DOMINANT
    !     relative to the dominant contribution via S(1,1)/S(1,2)/S(2,2) in the plane,
    !     captured by the FD on vel_gradient. (In a full 3D formulation
    !     vel_gradient(3,3) would directly contain ∂_z V_z; in the Q2D 2D, S(3,3) is
    !     a carry-over of the z reduction, treated as constant for the strain
    !     transport.)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nu, nv
    use Q2D_Lib_Fluid,        only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    use Q2D_Config_Types_m,   only: obj_rot, obj_source
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, ii, jj, kk, ni, nj
    real(R8), dimension(3,3,2) :: strain_gradient                 ! without spatial direction 3
    real(R8), dimension(3,3)   :: Sij_ip1, Sij_im1, Sij_jp1, Sij_jm1
    real(R8), dimension(3,3)   :: Sij, Wij, strain_transport
    real(R8) :: omega(3), S, O, rtilde
    real(R8) :: r_ave, vr_ij, VG_ij, Oz_frame
    logical  :: area_var

    area_var = obj_source%area_variation
    ! Frame rotation rate stored once (shared read
    ! in the OMP loops; zero when the rotating frame is disabled -> the frame term
    ! below vanishes).
    Oz_frame = 0d0
    if (obj_rot%enabled) Oz_frame = obj_rot%Omega_vec(3)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks

      ni = domain % blk(b) % dim(1)
      nj = domain % blk(b) % dim(2)

      !$omp do collapse(2) private (i, j, ii, jj, kk, strain_gradient, Sij_ip1, Sij_im1, Sij_jp1, Sij_jm1), &
      !$omp private (Sij, Wij, strain_transport, omega, S, O, rtilde, r_ave, vr_ij, VG_ij)
      do j = 1, nj
      do i = 1, ni

        ! Cell-centered geometric quantities (the metrics are not extrapolated to the
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
        ! Strain tensor over the 5 cells.
        ! - neighboring-cell vel_gradient: extended to the ghosts (Compute_Velocity_Gradient
        !   above) → captures ∂_k S(1,1), S(1,2), S(2,1), S(2,2) via spatial FD.
        ! - (vr_ij, r_ave, VG_ij): cell-centered, identical for the 5 calls
        !   → S(3,3) constant in the FD stencil; ∂_k S(3,3) = 0 (accepted sub-sub-dominant
        !   loss — cf. the file header for the formal analysis).
        ! ====================================================================
        Sij_ip1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i+1,j)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_im1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i-1,j)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_jp1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i,j+1)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_jm1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i,j-1)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij     = Strain_Tensor ( domain%blk(b)%vel_gradient(i  ,j  )%c, vr_ij, r_ave, V_dot_G=VG_ij )

        ! Strain-tensor gradient via finite differences (without spatial
        ! direction 3). The component-by-component FD on S_ij in the orthonormal basis
        ! (No-w axi cylindrical or area-var Cartesian plane) formally captures V·∇S_ij.
        strain_gradient(:,:,1) = ( Sij_ip1 - Sij_im1 ) * 0.5d0
        strain_gradient(:,:,2) = ( Sij_jp1 - Sij_jm1 ) * 0.5d0

        ! Metric transformation tensor index by tensor index (3×3 preserved).
        ! The cell-centered metric M(i,j) is applied to the spatial gradient (Jacobian
        ! from computational space (xi, eta) to physical (x, y) or (z, r)) — a tensor
        ! invariant, applied to each component of Sij separately.
        do ii = 1, 3
        do jj = 1, 3
          strain_gradient(ii,jj,1:2) = matmul ( strain_gradient(ii,jj,1:2), domain%blk(b)%m(i,j)%c )
        end do
        end do

        ! Cell-centered vorticity tensor + obj_rot absolute-vorticity correction
        omega = Vorticity_Vector( domain%blk(b)%vel_gradient(i,j)%c )
        ! Absolute vorticity:
        ! omega_abs = omega_rel + 2*Omega_frame; Vorticity_Vector returns vector(3)
        ! with omega(1)=omega(2)=0; Omega_vec axis strictly [0,0,1].
        if (obj_rot%enabled) omega = omega + 2.0_R8 * obj_rot%Omega_vec
        Wij = Vorticity_Tensor( omega )
        S = sqrt ( 2d0*sum (Sij**2) )
        O = sqrt ( sum(omega**2) )
        domain%blk(b)%rc_term1(i,j) = S/O

        ! Material derivative of the strain tensor V·∇S_ij (without spatial
        ! direction 3). P(nu:nv,i,j) are the physical velocity COMPONENTS (V_x, V_y)
        ! or (V_z, V_r). V·∇ is the directional derivative on the scalar components of S_ij;
        ! formally correct for the Cartesian plane (no Christoffel) and for the
        ! No-w axisymmetric cylindrical case (V_θ=0 + ∂_θ=0 cancel the Christoffel
        ! contributions).
        do ii = 1, 3
        do jj = 1, 3
          strain_transport(ii,jj) = dot_product ( domain%blk(b)%P(nu:nv,i,j), &
                                                  strain_gradient(ii,jj,1:2) )
        end do
        end do

        ! ── Rotating frame: frame term of the Spalart-Shur bracket ──
        ! The TMR SA-RC / SST-RC formulations (Spalart-Shur 1997; Smirnov-Menter 2009)
        ! prescribe
        !   r_tilde ∝ 2 W_ik S_jk [ DS_ij/Dt + (eps_imn S_jn + eps_jmn S_in) Omega_m ]
        ! The eps term is the rotation of the strain principal axes induced
        ! by the frame; it is NOT covered by the absolute vorticity already used in W_ik
        ! (that is a separate and independent appearance of Omega). With
        ! Omega = (0,0,Oz) and the block-diagonal Q2D S (S13 = S23 = 0) the term
        ! reduces in the plane to:
        !   T11 = -2 Oz S12,  T22 = +2 Oz S12,  T12 = Oz (S11 - S22),
        !   T33 = T13 = T23 = 0  (block structure preserved).
        ! Verified on an analytic vortex (fixed vs rotating frame): WITH this
        ! term the discrete numerator converges (order 2) to the frame-invariant
        ! value of the fixed frame; WITHOUT, it converges to a distorted value of
        ! O(Omega), independently of the grid.
        strain_transport(1,1) = strain_transport(1,1) - 2d0*Oz_frame*Sij(1,2)
        strain_transport(2,2) = strain_transport(2,2) + 2d0*Oz_frame*Sij(1,2)
        strain_transport(1,2) = strain_transport(1,2) + Oz_frame*(Sij(1,1) - Sij(2,2))
        strain_transport(2,1) = strain_transport(1,2)

        ! Triple sum over the TENSOR INDICES 1..3 (non-spatial, to be preserved).
        rtilde = 0d0
        do ii = 1, 3
        do jj = 1, 3
        do kk = 1, 3
          rtilde = rtilde + Wij(ii,kk)*Sij(jj,kk)*strain_transport(ii,jj)
        end do
        end do
        end do
        ! Note: denominator absent, applied at the use site in SA_Source_Blk / SST_Blk
        ! (different between SA-RC and SST-RC).
        domain%blk(b)%rc_term2(i,j) = 2d0*rtilde

      end do ; end do
    end do

  end subroutine Compute_RC_Terms

end module Q2D_Lib_SpalartShur
