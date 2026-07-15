module Q2D_Lib_BC_Fluxes
  ! Common utilities for the BCs (face indexing, trapezoidal quadratures at
  ! face centers for real-gas nozzle/inflow, packing of wall properties).
  ! 2D reduction: 4 faces, no Km metric.
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Global_m,    only: nsc

  implicit none
  private

  public :: Compute_Modfm
  public :: Face_Index
  public :: Compute_Wall_Properties
  public :: Trapezoidal
  public :: Trapezoidal_RiemU
  public :: Trapezoidal_Riem_PU

contains


  pure subroutine Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    implicit none
    integer, intent(in)  :: fm
    integer, intent(out) :: modfm, modfm1, modfm2, modfm3

    modFm  = mod(Fm,2)
    modfm1 = 1 - modFm
    modfm2 = 1 - 2*modFm
    modfm3 = 2*( 2*modFm - 1 )

  end subroutine Compute_Modfm


  pure subroutine Face_Index ( face, dir, cell_i, cell_j, face_i, face_j )
    implicit none
    integer, intent(in)  :: face, cell_i, cell_j
    integer, intent(out) :: dir, face_i, face_j

    select case ( face )
      case(1:2)
        dir = 1
        face_i = cell_i - mod ( face, 2 )
        face_j = cell_j
      case(3:4)
        dir = 2
        face_i = cell_i
        face_j = cell_j - mod ( face, 2 )
    end select

  end subroutine Face_Index


  subroutine Compute_Wall_Properties ( stress, pw, Tw, rhow, mu, qw, mdot, y, exit_array )
    implicit none
    real(R8), dimension(:), intent(in)    :: stress
    real(R8), intent(in)                  :: pw
    real(R8), dimension(7), intent(inout) :: exit_array
    real(R8), optional :: Tw, qw, rhow, mu, y, mdot
    real(R8) :: tauw, ut

    tauw = sqrt( sum ( Stress**2 ) )
    ut = sqrt ( tauw / rhow )
    exit_array(1) = rhow * y * ut / mu
    exit_array(2:3) = stress
    exit_array(4) = pw
    if (present(Tw))   exit_array(5) = Tw
    if (present(qw))   exit_array(6) = qw
    if (present(mdot)) exit_array(7) = mdot

  end subroutine Compute_Wall_Properties


  subroutine Trapezoidal ( sup_T, BC_T0, BC_ci, Rgas, integral )
    ! Composite trapezoidal quadrature of f(T) = cp(T)/(R*T) on the integer-step
    ! temperature grid (K) between min(sup_T, BC_T0) and max. The result is the
    ! oriented integral
    !   integral = int_{BC_T0}^{sup_T} cp/(R*T) dT
    ! (sign flipped at the end if sup_T > BC_T0, to respect the orientation).
    use FLINT_Lib_Thermodynamic, only: cp_tab, f_tabT
    implicit none
    real(R8), intent(in)  :: sup_T, BC_T0, Rgas, BC_ci(nsc)
    real(R8), intent(out) :: integral
    real(R8), allocatable :: T_vec(:)
    integer :: int_sx, int_dx, n_points, i, s
    real(R8) :: cp_node, f_prev, f_next

    if (sup_T > BC_T0) then
      int_sx = ceiling(BC_T0)
      int_dx = floor(Sup_T)
    else
      int_sx = ceiling(Sup_T)
      int_dx = floor(BC_T0)
    endif

    n_points = int_dx - int_sx + 3
    integral = 0d0
    allocate(T_vec(n_points))

    if (sup_T > BC_T0) then
      T_vec(1) = BC_T0
      T_vec(n_points) = Sup_T
    else
      T_vec(1) = Sup_T
      T_vec(n_points) = BC_T0
    endif
    do i = 2, n_points - 1
      T_vec(i) = int_sx + (i - 2)
    enddo

    cp_node = 0d0
    do s = 1, nsc
      cp_node = cp_node + BC_ci(s) * f_tabT(T_vec(1), s, cp_tab)
    enddo
    f_prev = cp_node / Rgas / T_vec(1)

    do i = 1, n_points - 1
      cp_node = 0d0
      do s = 1, nsc
        cp_node = cp_node + BC_ci(s) * f_tabT(T_vec(i+1), s, cp_tab)
      enddo
      f_next = cp_node / Rgas / T_vec(i+1)
      integral = integral + 0.5d0 * (T_vec(i+1) - T_vec(i)) * (f_prev + f_next)
      f_prev = f_next
    enddo

    if (sup_T > BC_T0) integral = - integral

  end subroutine Trapezoidal


  subroutine Trapezoidal_RiemU ( T_from, T_to, BC_ci, Rgas, integral )
    ! Integral of the generalized real-gas Riemann invariant (Glaister 1988):
    !   integral = int_{T_from}^{T_to} a(T) / ((gamma(T)-1)*T) dT
    ! with a(T) = sqrt(gamma(T)*Rgas*T), gamma(T) = cp(T)/(cp(T)-Rgas).
    ! For constant gamma it reduces exactly to 2(a_to - a_from)/(gamma - 1).
    use FLINT_Lib_Thermodynamic, only: cp_tab, f_tabT
    implicit none
    real(R8), intent(in)  :: T_from, T_to, Rgas, BC_ci(nsc)
    real(R8), intent(out) :: integral
    real(R8), allocatable :: T_vec(:)
    real(R8) :: f_prev, f_next, cp_node, gam_node
    integer  :: i, s, int_sx, int_dx, n_points

    if ( T_to > T_from ) then
      int_sx = ceiling( T_from )
      int_dx = floor  ( T_to   )
    else
      int_sx = ceiling( T_to   )
      int_dx = floor  ( T_from )
    endif
    n_points = int_dx - int_sx + 3
    allocate( T_vec(n_points) )

    if ( T_to > T_from ) then
      T_vec(1)        = T_from
      T_vec(n_points) = T_to
    else
      T_vec(1)        = T_to
      T_vec(n_points) = T_from
    endif
    do i = 2, n_points - 1
      T_vec(i) = int_sx + ( i - 2 )
    enddo

    cp_node = 0.d0
    do s = 1, nsc
      cp_node = cp_node + BC_ci(s) * f_tabT( T_vec(1), s, cp_tab )
    enddo
    gam_node = cp_node / ( cp_node - Rgas )
    f_prev   = sqrt( gam_node * Rgas * T_vec(1) ) / ( ( gam_node - 1.d0 ) * T_vec(1) )

    integral = 0.d0
    do i = 1, n_points - 1
      cp_node = 0.d0
      do s = 1, nsc
        cp_node = cp_node + BC_ci(s) * f_tabT( T_vec(i+1), s, cp_tab )
      enddo
      gam_node = cp_node / ( cp_node - Rgas )
      f_next   = sqrt( gam_node * Rgas * T_vec(i+1) ) / ( ( gam_node - 1.d0 ) * T_vec(i+1) )
      integral = integral + 0.5d0 * ( T_vec(i+1) - T_vec(i) ) * ( f_prev + f_next )
      f_prev   = f_next
    enddo

    if ( T_from > T_to ) integral = - integral

  end subroutine Trapezoidal_RiemU


  subroutine Trapezoidal_Riem_PU ( T_from, T_to, BC_ci, Rgas, J_int, K_int )
    ! Fused trapezoidal quadrature returning BOTH the isentropic pressure
    ! integral J AND the generalized Riemann-invariant integral K in a single
    ! traversal of the grid (one species loop per node: half the cost of
    ! Trapezoidal + Trapezoidal_RiemU). Results are bit-exact with respect to
    ! those two routines.
    use FLINT_Lib_Thermodynamic, only: cp_tab, f_tabT
    implicit none
    real(R8), intent(in)  :: T_from, T_to, Rgas, BC_ci(nsc)
    real(R8), intent(out) :: J_int, K_int
    real(R8), allocatable :: T_vec(:)
    real(R8) :: cp_node, gam_node, fJ_prev, fJ_next, fK_prev, fK_next, dT_seg
    integer  :: i, s, int_sx, int_dx, n_points

    if ( T_to > T_from ) then
      int_sx = ceiling( T_from )
      int_dx = floor  ( T_to   )
    else
      int_sx = ceiling( T_to   )
      int_dx = floor  ( T_from )
    endif
    n_points = int_dx - int_sx + 3
    allocate( T_vec(n_points) )

    if ( T_to > T_from ) then
      T_vec(1)        = T_from
      T_vec(n_points) = T_to
    else
      T_vec(1)        = T_to
      T_vec(n_points) = T_from
    endif
    do i = 2, n_points - 1
      T_vec(i) = int_sx + ( i - 2 )
    enddo

    cp_node = 0.d0
    do s = 1, nsc
      cp_node = cp_node + BC_ci(s) * f_tabT( T_vec(1), s, cp_tab )
    enddo
    gam_node = cp_node / ( cp_node - Rgas )
    fJ_prev  = cp_node / Rgas / T_vec(1)
    fK_prev  = sqrt( gam_node * Rgas * T_vec(1) ) / ( ( gam_node - 1.d0 ) * T_vec(1) )

    J_int = 0.d0
    K_int = 0.d0
    do i = 1, n_points - 1
      cp_node = 0.d0
      do s = 1, nsc
        cp_node = cp_node + BC_ci(s) * f_tabT( T_vec(i+1), s, cp_tab )
      enddo
      gam_node = cp_node / ( cp_node - Rgas )
      fJ_next  = cp_node / Rgas / T_vec(i+1)
      fK_next  = sqrt( gam_node * Rgas * T_vec(i+1) ) / ( ( gam_node - 1.d0 ) * T_vec(i+1) )
      dT_seg   = T_vec(i+1) - T_vec(i)
      J_int    = J_int + 0.5d0 * dT_seg * ( fJ_prev + fJ_next )
      K_int    = K_int + 0.5d0 * dT_seg * ( fK_prev + fK_next )
      fJ_prev  = fJ_next
      fK_prev  = fK_next
    enddo

    if ( T_from > T_to ) then
      J_int = - J_int
      K_int = - K_int
    endif

  end subroutine Trapezoidal_Riem_PU


end module Q2D_Lib_BC_Fluxes
