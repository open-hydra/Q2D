!>@brief Implicit Residual Smoothing with structural 3D->2D reduction:
!>       loop only over directions i,j (ndir=2); rank-3 fields Res(nprim, i, j)
!>       instead of rank-4 Res(nprim, i, j, k); 2-direction Jacobi blending
!>       (njacobi=2).
module Q2D_Lib_IRS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Residual_Smoothing

  integer, parameter, private :: njacobi = 2

contains

  ! Compute the Implicit Residual Smoothing operator: R* = R + beta*LAPLACE(R*)
  ! with Jacobi iterations. 2D form (ndir=2): direction loop su i,j only.
  subroutine Residual_Smoothing(domain, beta)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: ndir
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    real(R8),              intent(in)    :: beta
    ! Local
    integer :: b, d

    do d = 1, ndir

      ! Residual extrapolation in the ghost cells (BC + MPI exchange).
      call Ghost_Residual_Extrapolation(domain)

      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle

        call Residual_Smoothing_Blk(domain%blk(b)%R,   &
                                    domain%blk(b)%RS1, &
                                    domain%blk(b)%RS2, &
                                    domain%blk(b)%dim, d, beta)

      end do ! blocks

    end do ! direction

  end subroutine Residual_Smoothing


  subroutine Residual_Smoothing_Blk(Res, Res_Star, Res_Star_New, n, Dir, beta)
    use Q2D_Global_m, only: nprim, np, gc
    implicit none
    integer,  intent(in) :: n(2), Dir
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res, Res_Star, Res_Star_New
    real(R8), intent(in) :: beta
    ! Local
    integer  :: i, j, jacobi
    integer  :: i1, j1, i2, j2
    real(R8), dimension(nprim) :: r, rm, rp

    ! Res* initialized as Res in the block, including ghost cells.
    !$omp do collapse(2)
    do j = 0, n(2) + 1
    do i = 0, n(1) + 1

      Res_Star(:,i,j) = Res(:,i,j)

    end do; end do

    ! Smoothing with jacobi iterations.
    do jacobi = 1, njacobi

      !$omp do collapse(2) private(r, rm, rp, i1, i2, j1, j2)
      do j = 1, n(2)
      do i = 1, n(1)

        call Jacobi_Stencil(i1, j1, i2, j2, i, j, Dir)
        r (1:np) = Res     (1:np, i,  j )
        rm(1:np) = Res_Star(1:np, i1, j1)
        rp(1:np) = Res_Star(1:np, i2, j2)

        ! Jacobi: r*_New(i) = [ r(i) + e*( r*(i-1) + r*(i+1) ) ] / [ 1 + 2e ]
        ! Only smooth flow equations (1:np); turbulence residuals are left untouched.
        Res_Star_New(1:np, i, j) = (r(1:np) + beta*(rm(1:np) + rp(1:np))) / (1d0 + 2d0*beta)

      end do; end do

      !$omp do collapse(2)
      do j = 1, n(2)
      do i = 1, n(1)

        Res_Star(1:np, i, j) = Res_Star_New(1:np, i, j)

      end do; end do

    end do ! jacobi

    ! R* after directional smoothing becomes the new R.
    !$omp do collapse(2)
    do j = 1, n(2)
    do i = 1, n(1)

      Res(1:np, i, j) = Res_Star(1:np, i, j)

    end do; end do

  contains

    subroutine Jacobi_Stencil(i1, j1, i2, j2, i, j, d)
      implicit none
      integer, intent(in)  :: i, j, d
      integer, intent(out) :: i1, j1, i2, j2

      ! 2D form: the k-direction case (3) is dropped.
      select case (d)
        case (1)
          i1 = i - 1
          j1 = j
          i2 = i + 1
          j2 = j
        case (2)
          i1 = i
          j1 = j - 1
          i2 = i
          j2 = j + 1
      end select

    end subroutine Jacobi_Stencil

  end subroutine Residual_Smoothing_Blk


  subroutine Ghost_Residual_Extrapolation(domain)
    use Q2D_Advanced_Types_m,  only: Q2D_domain_type
    use Q2D_Parameters_m,      only: guide
    use Q2D_Mod_MPI,           only: is_local_block
    use Q2D_Mod_GhostExchange, only: exchange_ghost_R
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, ig, jg
    integer :: bm, im, jm, fm, bs, is, js

    ! MPI: exchange residual data for inter-rank type-1 connections.
    !$omp single
    call exchange_ghost_R(domain)
    !$omp end single

    ! BC array processing — iterate only over local BC entries.
    !$omp do private(ii, i, bm, im, jm, fm, bs, is, js, ig, jg)
    do ii = 1, domain%n_local_bc
      i = domain%local_bc_idx(ii)

      ! Preliminary assignments (k-direction dropped; guide(:,3) preserved from Q2D_Parameters_m).
      bm = domain%bc(i)%b
      im = domain%bc(i)%i
      jm = domain%bc(i)%j
      fm = domain%bc(i)%f

      ! Ghost cell coordinates (2D form: indici i,j only).
      ig = im - guide(fm, 1)
      jg = jm - guide(fm, 2)

      select case (domain%bc(i)%type)

        case (101) ! block connection
          bs = domain%bc(i)%bs
          if (.not. is_local_block(bs)) cycle  ! inter-rank handled by exchange_ghost_R
          is = domain%bc(i)%is
          js = domain%bc(i)%js

          domain%blk(bm)%R(:, ig, jg) = domain%blk(bs)%R(:, is, js)

        case default ! wall, inlet, outlet — zero-gradient ghost
          domain%blk(bm)%R(:, ig, jg) = domain%blk(bm)%R(:, im, jm)

      end select
    end do

  end subroutine Ghost_Residual_Extrapolation

end module Q2D_Lib_IRS
