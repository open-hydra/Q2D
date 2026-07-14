!> @brief Single Rotating Frame (SRF) source terms — Q2D 2D-planar + q2D area-variation.
!>
!> Rigid 2D reduction (No-w principle:
!> no nw, no Res(nw), no swirl). Scope:
!>
!>   * 2D planar (`obj_source%axisymmetric = .false.`, axis = [0,0,1]) — VALID.
!>   * q2D area-variation (`obj_source%area_variation = .true.`, axis = [0,0,1]) — VALID.
!>   * 2D-axi (`obj_source%axisymmetric = .true.`) — OUT-OF-SCOPE (the only physically
!>     meaningful axis is [1,0,0] and it would require swirl `u_phi = Prim(nw)`, out of
!>     scope under No-w).
!>
!> Setup_RotatingFrame applies a runtime dual-check:
!>   Check #1 → error stop if `obj_source%axisymmetric = .true.` (SKIP-RUNTIME 2D-axi).
!>   Check #2 → error stop if `axis ≠ [0, 0, 1]` (tolerance 1e-12).
!>
!> Governing equations added to the residual:
!>
!>   Momentum (per unit surface, sur is the 2D analogue of the 3D vol):
!>     S_mom = -2 rho (Omega x w)         [Coriolis]
!>           + rho omega^2 r_perp          [Centrifugal]
!>     where `r_perp = r_vec - (r_vec · axis) axis`.
!>     With axis = [0,0,1] and w(3) = 0:
!>       Cor(1) = +2 rho omega w(2)
!>       Cor(2) = -2 rho omega w(1)
!>       Cor(3) = 0  (identically; not representable, dropped from the Res(nw) line)
!>       Cen(3) = 0  (r_vec(3) = 0 in the plane)
!>
!>   Energy (per unit surface):
!>     S_E = w · (rho omega^2 r_perp) = w(1)*Cen(1) + w(2)*Cen(2)
!>     [Coriolis does no work — `w · (Omega x w) = 0` identically.]
!>
!> Rotating-wall BC: `RF_Wall_Velocity(r_face)` produces `w_wall(2)` per face,
!> passed as an optional argument to `BC_Wall_Heat / BC_Wall_Temperature`. Walls
!> listed in `obj_rot%stationary_face` (casing/hub) keep `w_wall = 0`. The Q2D-form
!> symmetry-Eul does not receive `w_wall`: specular reflection is invariant under
!> a rigid-body frame translation with zero normal velocity.
!>
!> Rotating inflow BC: `RF_Convert_BC_Inflow` converts `(T0, p0, alpha)` from the
!> absolute frame to the relative frame BEFORE dispatch to `BC_Inflow_Outflow`.
!> `beta` is enforced to 0 (No-w: no out-of-plane velocity component).
module Q2D_Lib_RotatingFrame
  use iso_fortran_env,      only: I4 => int32, R8 => real64
  use Q2D_Config_Types_m,   only: obj_rot, obj_source

  implicit none
  private

  public :: Setup_RotatingFrame
  public :: RotatingFrame_Source_Terms
  public :: RF_Wall_Velocity
  public :: RF_Is_Stationary_Face
  public :: RF_Convert_BC_Inflow
  public :: RF_Face_Center
  public :: RF_Torque_Walls

contains


  ! ---------------------------------------------------------------------------
  !> @brief Finalise the rotating-frame configuration.
  !>
  !> Includes an additional dual-check enforcing the
  !> Q2D 2D-planar / q2D area-variation scope.
  subroutine Setup_RotatingFrame()
    implicit none
    integer  :: m
    real(R8) :: norm

    if ( obj_rot%model == "" .or. trim(obj_rot%model) == "none" ) then
      obj_rot%enabled = .false.
      return
    elseif ( obj_rot%model == "rigid-body" ) then
      if ( obj_rot%omega < 0.0_R8 ) then
        write(*,'(A)') '[ERROR] Rotational frame model "rigid-body" requires omega >= 0'
        stop
      end if
      obj_rot%enabled = .true.
    end if

    ! Parse axis and origin from registry strings (scalar str).
    read(obj_rot%axis_str, *)   obj_rot%axis
    read(obj_rot%origin_str, *) obj_rot%origin

    ! Normalise axis to unit vector.
    norm = sqrt( sum(obj_rot%axis**2) )
    if ( norm < 1.0e-14_R8 ) then
      write(*,'(A)') '[ERROR] Rotation axis has zero length'
      stop
    end if
    obj_rot%axis = obj_rot%axis / norm

    ! ─── Q2D dual-check (No-w + scope) ───────────────────────────────────────
    !
    ! Check #1: 2D-axi out-of-scope. In axisymmetric geometry the only physically
    ! meaningful rotation axis is [1,0,0] (axis of symmetry) and it requires
    ! swirl `u_phi = Prim(nw)`, out of scope under No-w.
    if ( obj_source%axisymmetric ) then
      write(*,'(A)') '[ERROR] RotatingFrame OUT-OF-SCOPE in 2D-axi (axisymmetric) configuration.'
      write(*,'(A)') '        2D-axi rotation around axis-of-symmetry [1,0,0] requires swirl velocity'
      write(*,'(A)') '        component u_phi (= Prim(nw)), which is OUT-OF-SCOPE under No-w.'
      write(*,'(A)') '        2D-axi rotating frames are reserved for a follow-up activity.'
      stop
    end if
    !
    ! Check #2: 2D planar + q2D area-variation require axis = [0,0,1] (axis
    ! perpendicular to the plane). Coriolis/Centrifugal stay in-plane with no
    ! out-of-plane components (Cor(3) ≡ 0, Cen(3) = 0).
    if ( any( abs(obj_rot%axis - [0d0, 0d0, 1d0]) > 1.0e-12_R8 ) ) then
      write(*,'(A)')           '[ERROR] In Q2D 2D-planar / q2D area-variation, RotatingFrame axis must be [0, 0, 1].'
      write(*,'(A)')           '        Rationale: an axis perpendicular to the (x,y) plane is the only choice that'
      write(*,'(A)')           '        keeps Coriolis/Centrifugal in-plane without requiring out-of-plane swirl.'
      write(*,'(A)')           '        Q2D solves only in the (x,y) plane.'
      write(*,'(A,3ES12.4)')   '        Got axis = ', obj_rot%axis
      stop
    end if
    ! ─────────────────────────────────────────────────────────────────────────

    obj_rot%Omega_vec = obj_rot%omega * obj_rot%axis

    ! Compatibility q2D area-variation + RF: both operate in-plane (the variable
    ! section h(x,y) does not interact with the rotating-frame kinematics).
    ! Informative warning, not an error stop.
    if ( obj_source%area_variation ) then
      obj_rot%warning_message = 'RotatingFrame combined with q2D area-variation: both in-plane (compatible).'
    end if

    ! Parse stationary wall patches from registry strings.
    ! Counter `m` (not `k`): the No-w grep audit `do k =` would false-positive here;
    ! m is a logical counter over RF patches, not a spatial Z index.
    if ( obj_rot%n_stationary > 0 ) then
      allocate( obj_rot%stationary_face(2, obj_rot%n_stationary) )
      do m = 1, obj_rot%n_stationary
        read(obj_rot%stationary_face_str(m), *) &
          obj_rot%stationary_face(1,m), obj_rot%stationary_face(2,m)
      end do
    end if

  end subroutine Setup_RotatingFrame


  ! ---------------------------------------------------------------------------
  !> @brief Add Coriolis and centrifugal source terms to the residual array of
  !>        every block in `domain`.
  subroutine RotatingFrame_Source_Terms ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: b

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      call RF_Source_Blk( domain%blk(b)%P,    &
                          domain%blk(b)%R,    &
                          domain%blk(b)%node, &
                          domain%blk(b)%sur,  &
                          domain%blk(b)%dim   )
    end do

  end subroutine RotatingFrame_Source_Terms


  ! ---------------------------------------------------------------------------
  !> @brief Block-level kernel: accumulate rotating-frame source terms into Res.
  !>
  !> Rigid 2D reduction of the 3D block kernel:
  !>   * drop the k loop (No-w);
  !>   * cell-centre 4-corner average (vs the MOSE peer's 8-corner);
  !>   * Res(nw) line removed (Cor(3) ≡ 0, Cen(3) = 0 with axis [0,0,1] strict);
  !>   * `vol → sur` (the 2D surface metric replaces the 3D volume metric).
  !>
  !> Sign convention: V dU/dt = -R(U), R = CONV - DIFF - SOURCE  ⇒  R -= S * sur.
  subroutine RF_Source_Blk ( Prim, Res, node, sur, n )
    use Q2D_Base_Types_m, only: Q2D_vector_2D_type
    use Q2D_Global_m,     only: nprim, nsc, nu, nv, np, gc

    implicit none
    integer,  intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)),  intent(in) :: node
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc),      intent(in) :: sur
    ! Local
    integer  :: i, j
    real(R8) :: rho, w(2)
    real(R8) :: xc, yc, r_vec(2)
    real(R8) :: Cor(2), Cen(2)
    real(R8) :: omega_z, ox(2), omega2

    ! Cache module-level data in local scalars to avoid repeated indirect lookups.
    ! With axis = [0,0,1] strict (Setup_RotatingFrame Check #2): Omega_vec = (0,0,omega_z),
    ! origin projected in-plane (z component discarded: r_vec is 2D).
    omega_z = obj_rot%Omega_vec(3)
    ox      = obj_rot%origin(1:2)
    omega2  = obj_rot%omega**2

    !$omp do collapse(2) &
    !$omp private( i, j, rho, w, xc, yc, r_vec, Cor, Cen )
    do j = 1, n(2)
    do i = 1, n(1)

      rho = sum( Prim(1:nsc, i, j) )
      w   = Prim(nu:nv, i, j)            ! Relative velocity in the rotating frame (2 components)

      ! ---- Cell centre (4-corner average; vs the MOSE peer's 8-corner for the 2D reduction) ----
      ! Cell (i,j) occupies the quadrilateral i-1..i, j-1..j in rank-2 node coordinates.
      xc = 0.25_R8 * ( node(i-1,j-1)%c(1) + node(i,j-1)%c(1) &
                     + node(i-1,j  )%c(1) + node(i,j  )%c(1) )
      yc = 0.25_R8 * ( node(i-1,j-1)%c(2) + node(i,j-1)%c(2) &
                     + node(i-1,j  )%c(2) + node(i,j  )%c(2) )

      ! ---- Position vector from axis origin (in-plane) ----
      r_vec = [ xc - ox(1), yc - ox(2) ]
      ! r_vec(3) = 0 implicito; r_proj = dot(r_vec, [0,0,1]) = 0; r_perp = r_vec.

      ! ---- Coriolis force per unit surface: -2 rho (Omega x w) -----------
      ! With Omega = (0,0,omega_z), w(3) = 0:
      !   Cor(1) = -2 rho ( Omega(2) w(3) - Omega(3) w(2) ) = +2 rho omega_z w(2)
      !   Cor(2) = -2 rho ( Omega(3) w(1) - Omega(1) w(3) ) = -2 rho omega_z w(1)
      !   Cor(3) ≡ 0 (not representable, no Res(nw)).
      Cor(1) =  2.0_R8 * rho * omega_z * w(2)
      Cor(2) = -2.0_R8 * rho * omega_z * w(1)

      ! ---- Centrifugal force per unit surface: rho omega^2 r_perp -------
      !   Cen(3) = 0 with r_vec(3) = 0.
      Cen(1) = rho * omega2 * r_vec(1)
      Cen(2) = rho * omega2 * r_vec(2)

      ! ---- Accumulate into residuals: R -= source * sur -
      Res(nu, i, j) = Res(nu, i, j) - ( Cor(1) + Cen(1) ) * sur(i,j)
      Res(nv, i, j) = Res(nv, i, j) - ( Cor(2) + Cen(2) ) * sur(i,j)
      ! No Res(nw) line: Cor(3) + Cen(3) ≡ 0 identically (No-w).

      ! Energy: only Cen does work on w; w · Coriolis ≡ 0 identically.
      Res(np, i, j) = Res(np, i, j) - dot_product( w, Cen ) * sur(i,j)

    end do ; end do
    !$omp end do

  end subroutine RF_Source_Blk


  ! ---------------------------------------------------------------------------
  !> @brief Velocity of a lab-stationary wall as seen from the rotating frame.
  !>
  !> 2D reduction. With axis = [0,0,1] strict:
  !>   w_wall(1) = +omega_z * (r_face(2) - origin(2))
  !>   w_wall(2) = -omega_z * (r_face(1) - origin(1))
  !>   w_wall(3) ≡ 0 (No-w: dropped).
  !>
  !> @param[in]  r_face  Cartesian 2D face-centre position [m].
  !> @return     w_wall  Relative velocity of the lab-stationary wall [m/s].
  function RF_Wall_Velocity(r_face) result(w_wall)
    implicit none
    real(R8), intent(in) :: r_face(2)
    real(R8)             :: w_wall(2)
    real(R8) :: r(2), omega_z

    r       = r_face - obj_rot%origin(1:2)
    omega_z = obj_rot%Omega_vec(3)
    w_wall(1) =  omega_z * r(2)
    w_wall(2) = -omega_z * r(1)
  end function RF_Wall_Velocity


  ! ---------------------------------------------------------------------------
  !> @brief 2D face centre as the average of 2 corner nodes (face = edge in 2D).
  !>
  !> 2D reduction. Q2D has 4 faces per cell (no k-faces
  !> 5/6); Fm ∈ {1,2,3,4} only. Fm = 5 or 6 produces an explicit error stop (no
  !> silent fallback).
  subroutine RF_Face_Center(node, Im, Jm, Fm, r_face)
    use Q2D_Base_Types_m, only: Q2D_vector_2D_type
    implicit none
    type(Q2D_vector_2D_type), intent(in)  :: node(0:,0:)
    integer,                  intent(in)  :: Im, Jm, Fm
    real(R8),                 intent(out) :: r_face(2)
    integer :: Fi, Fj, d

    select case (Fm)
    case(1,2)
      Fi = Im - mod(Fm,2) ; Fj = Jm
    case(3,4)
      Fi = Im             ; Fj = Jm - mod(Fm,2)
    case(5,6)
      write(*,'(A,I0,A)') '[ERROR] RF_Face_Center: Fm=', Fm, ' (k-face) invalid in Q2D 2D-planar.'
      stop
    case default
      write(*,'(A,I0,A)') '[ERROR] RF_Face_Center: Fm=', Fm, ' out of range [1..4].'
      stop
    end select

    do d = 1, 2
      select case (Fm)
      case(1,2)   ! i-face: edge between (Fi, Fj-1) and (Fi, Fj)
        r_face(d) = 0.5_R8 * ( node(Fi,Fj-1)%c(d) + node(Fi,Fj)%c(d) )
      case(3,4)   ! j-face: edge between (Fi-1, Fj) and (Fi, Fj)
        r_face(d) = 0.5_R8 * ( node(Fi-1,Fj)%c(d) + node(Fi,Fj)%c(d) )
      end select
    end do

  end subroutine RF_Face_Center


  ! ---------------------------------------------------------------------------
  !> @brief Convert inflow BC `(T0, p0, alpha, Mach_abs) → (T0_rel, p0_rel, alpha_rel, Mach_rel)`.
  !>
  !> Rigid 2D reduction:
  !>   * drop the `Km` argument (no k-direction).
  !>   * drop the `beta` argument (out-of-plane elevation angle); enforced to 0.
  !>   * v_frame(3) ≡ 0 with Omega = (0,0,omega_z), r_face(3) = 0.
  !>
  !> Static conditions are frame-invariant; total conditions and angle must be
  !> converted. Outflow static pressure (pAmb) is frame-independent — the caller
  !> dispatcher does not call this routine for pure outflow.
  subroutine RF_Convert_BC_Inflow ( Blk, Fm, Im, Jm, BC_Mach_abs, gamma, Rgas, &
                                    T0_or_Tstat, p0_or_pstat, alpha, Mach_rel )
    use Q2D_Advanced_Types_m, only: Q2D_block_type
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Fm, Im, Jm
    real(R8), intent(in)    :: BC_Mach_abs, gamma, Rgas
    real(R8), intent(inout) :: T0_or_Tstat, p0_or_pstat, alpha
    real(R8), intent(out)   :: Mach_rel
    ! Local
    integer  :: Face_i, Face_j
    real(R8) :: v_frame(2), r(2), r_fc(2)
    real(R8) :: T_stat, p_stat, a_stat, u_abs(2), w_rel(2), w_norm, cp
    real(R8) :: omega_z

    ! Face index (mirror Q2D Face_Index convention: Fm=1,2 → i-face; Fm=3,4 → j-face).
    select case (Fm)
    case(1,2) ; Face_i = Im - mod(Fm,2) ; Face_j = Jm
    case(3,4) ; Face_i = Im             ; Face_j = Jm - mod(Fm,2)
    case default
      write(*,'(A,I0,A)') '[ERROR] RF_Convert_BC_Inflow: Fm=', Fm, ' invalid in Q2D 2D-planar (expected 1..4).'
      stop
    end select

    ! Face centre (2-node average, 2D).
    select case (Fm)
    case(1,2)   ! i-face: edge between (Face_i, Face_j-1) and (Face_i, Face_j)
      r_fc(1) = 0.5_R8 * ( Blk%node(Face_i,Face_j-1)%c(1) + Blk%node(Face_i,Face_j)%c(1) )
      r_fc(2) = 0.5_R8 * ( Blk%node(Face_i,Face_j-1)%c(2) + Blk%node(Face_i,Face_j)%c(2) )
    case(3,4)   ! j-face: edge between (Face_i-1, Face_j) and (Face_i, Face_j)
      r_fc(1) = 0.5_R8 * ( Blk%node(Face_i-1,Face_j)%c(1) + Blk%node(Face_i,Face_j)%c(1) )
      r_fc(2) = 0.5_R8 * ( Blk%node(Face_i-1,Face_j)%c(2) + Blk%node(Face_i,Face_j)%c(2) )
    end select

    ! Frame velocity at face centre: v_frame = Omega x (r_fc - origin) with Omega = (0,0,omega_z).
    !   v_frame(1) = -omega_z * r(2)
    !   v_frame(2) =  omega_z * r(1)
    !   v_frame(3) = 0  (dropped).
    r          = r_fc - obj_rot%origin(1:2)
    omega_z    = obj_rot%Omega_vec(3)
    v_frame(1) = -omega_z * r(2)
    v_frame(2) =  omega_z * r(1)

    ! Static conditions (frame-invariant).
    if ( BC_Mach_abs > 1.0_R8 ) then
      ! Supersonic: input T0_or_Tstat = T_static, p0_or_pstat = p_static.
      T_stat = T0_or_Tstat
      p_stat = p0_or_pstat
    else
      ! Subsonic: input T0_or_Tstat = T0_abs, p0_or_pstat = p0_abs.
      T_stat = T0_or_Tstat / ( 1.0_R8 + 0.5_R8*(gamma-1.0_R8)*BC_Mach_abs**2 )
      p_stat = p0_or_pstat * ( T_stat / T0_or_Tstat )**( gamma/(gamma-1.0_R8) )
    end if

    ! Speed of sound (frame-invariant).
    a_stat = sqrt( gamma * Rgas * T_stat )

    ! Absolute velocity (no z component, beta enforced 0).
    !   u = M*a * [cos(alpha), sin(alpha)]  with beta = 0 ⇒ cos(beta) = 1, sin(beta) = 0.
    u_abs(1) = BC_Mach_abs * a_stat * cos(alpha)
    u_abs(2) = BC_Mach_abs * a_stat * sin(alpha)

    ! Relative velocity: w = u_abs - v_frame.
    w_rel  = u_abs - v_frame
    w_norm = norm2( w_rel )

    ! Relative Mach number and new flow angle (beta stays 0 by construction in 2D).
    Mach_rel = w_norm / a_stat
    alpha    = atan2( w_rel(2), w_rel(1) )

    if ( BC_Mach_abs > 1.0_R8 ) then
      ! T_static and p_static unchanged (frame-invariant); Mach_rel and alpha updated above.
    else
      ! Relative total conditions from relative velocity.
      cp          = gamma * Rgas / ( gamma - 1.0_R8 )
      T0_or_Tstat = T_stat + 0.5_R8 * w_norm**2 / cp
      p0_or_pstat = p_stat * ( T0_or_Tstat / T_stat )**( gamma/(gamma-1.0_R8) )
    end if

  end subroutine RF_Convert_BC_Inflow


  ! ---------------------------------------------------------------------------
  !> @brief `.true.` if wall patch (block b, face direction f) is lab-stationary.
  !>
  !> Pure logical, regime-independent.
  logical function RF_Is_Stationary_Face(b, f)
    implicit none
    integer, intent(in) :: b, f
    integer :: m

    RF_Is_Stationary_Face = .false.
    if ( .not. allocated(obj_rot%stationary_face) ) return
    do m = 1, size(obj_rot%stationary_face, 2)
      if ( obj_rot%stationary_face(1,m) == b .and. &
           obj_rot%stationary_face(2,m) == f ) then
        RF_Is_Stationary_Face = .true.
        return
      end if
    end do

  end function RF_Is_Stationary_Face


  ! ---------------------------------------------------------------------------
  !> @brief Pressure torque on rotating wall surfaces; writes OUTPUT/torque.dat.
  !>
  !> 2D reduction. In 2D-planar with axis=[0,0,1]:
  !>   torque_vec(1) = r_rel(2)*F_wall(3) - r_rel(3)*F_wall(2) = 0  (F_wall(3)=r_rel(3)=0)
  !>   torque_vec(2) = r_rel(3)*F_wall(1) - r_rel(1)*F_wall(3) = 0
  !>   torque_vec(3) = r_rel(1)*F_wall(2) - r_rel(2)*F_wall(1)  = T_axial in-plane
  !>   T_axial = axis · torque_vec = torque_vec(3)
  !>
  !> Q2D uses `int_l` (interface length) as the 2D analogue of the 3D `A_face`
  !> (cf. `Q2D_f_metrics_type` in `Base_Types_m.f90`). MPI: the partial torque
  !> scalar is reduced to root via `mpi_reduce_sum_r8` (single-rank no-op),
  !> mirror MOSE Lib_RotatingFrame.f90.
  subroutine RF_Torque_Walls(domain, iter, time)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: np
    use Q2D_Lib_Ghost,        only: Is_Wall
    use Q2D_Mod_MPI,          only: is_local_block, mpi_is_root, mpi_reduce_sum_r8
    implicit none
    type(Q2D_domain_type), intent(in) :: domain
    integer,  intent(in) :: iter
    real(R8), intent(in) :: time
    ! Local
    integer,  save :: unit_torque = -1
    logical,  save :: first_call  = .true.
    integer  :: i, Bm, Im, Jm, Fm, Dir, Fi, Fj
    real(R8) :: r_face(2), r_rel(2), F_wall(2)
    real(R8) :: T_axial, p_face, sign_fm
    real(R8) :: N_face(2), int_l_face

    T_axial = 0.0_R8

    do i = 1, domain%nbound

      if ( .not. Is_Wall(domain%bc(i)%type) ) cycle

      Bm = domain%bc(i)%b
      if ( .not. is_local_block(Bm) ) cycle

      Fm = domain%bc(i)%f

      ! Skip lab-stationary faces (casing / hub).
      if ( RF_Is_Stationary_Face(Bm, Fm) ) cycle

      Im = domain%bc(i)%i
      Jm = domain%bc(i)%j

      ! Direction index: 1 for i-faces (Fm=1,2), 2 for j-faces (3,4). No k-faces.
      Dir = (Fm + 1) / 2

      select case (Fm)
      case (1,2) ; Fi = Im - mod(Fm,2) ; Fj = Jm
      case (3,4) ; Fi = Im             ; Fj = Jm - mod(Fm,2)
      end select

      ! Q2D: dir(d)%f%N (unit normal) + dir(d)%f%int_l (face length, 2D analogue of 3D A_face).
      N_face     = domain%blk(Bm)%dir(Dir)%f(Fi,Fj)%N
      int_l_face = domain%blk(Bm)%dir(Dir)%f(Fi,Fj)%int_l

      ! Outward-from-fluid sign: -1 for min faces (Fm odd), +1 for max faces (Fm even).
      sign_fm = real(1 - 2 * mod(Fm, 2), R8)

      ! Pressure at interior cell adjacent to wall (zero-order extrapolation).
      p_face = domain%blk(Bm)%P(np, Im, Jm)

      call RF_Face_Center(domain%blk(Bm)%node, Im, Jm, Fm, r_face)

      ! Pressure force on wall from fluid: F_wall = p * (outward normal) * int_l.
      F_wall = p_face * sign_fm * N_face * int_l_face

      ! In-plane torque about rotation axis origin (scalar, axis=[0,0,1]):
      r_rel   = r_face - obj_rot%origin(1:2)
      T_axial = T_axial + ( r_rel(1)*F_wall(2) - r_rel(2)*F_wall(1) )

    end do

    ! Reduce partial torque sums from all ranks to root (blocking collective;
    ! single-rank no-op). Mirror MOSE Lib_RotatingFrame.f90:572.
    call mpi_reduce_sum_r8(T_axial)

    if ( mpi_is_root ) then
      if ( first_call ) then
        open(newunit=unit_torque, file='OUTPUT/torque.dat', &
             status='replace', form='formatted')
        write(unit_torque, '(A)') &
          '#        iter              time   torque_axial (pressure, non-dim)' // &
          '         power (T*omega, non-dim)'
        first_call = .false.
      end if

      write(unit_torque, '(I9,3E20.10)') iter, time, T_axial, T_axial * obj_rot%omega
      flush(unit_torque)
    end if

  end subroutine RF_Torque_Walls


end module Q2D_Lib_RotatingFrame
