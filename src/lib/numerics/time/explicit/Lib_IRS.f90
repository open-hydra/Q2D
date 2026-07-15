!>@brief Implicit Residual Smoothing in 2D form: loop over the i,j directions
!>       only (ndir=2); rank-3 fields Res(nprim, i, j); Jacobi blending
!>       over 2 directions (njacobi=2).
module Q2D_Lib_IRS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Residual_Smoothing

  integer, parameter, private :: njacobi = 2

contains

  ! Compute the Implicit Residual Smoothing operator: R* = R + beta*LAPLACE(R*)
  ! with Jacobi iterations. 2D form (ndir=2): direction loop over i,j only.
  subroutine Residual_Smoothing(domain, beta)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: ndir
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    real(R8),              intent(in)    :: beta
    ! Local variables
    integer :: b, d

    do d = 1, ndir

      ! Extrapolate the residual into the ghost cells (BCs + MPI exchange).
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
    ! Local variables
    integer  :: i, j, jacobi
    integer  :: i1, j1, i2, j2
    real(R8), dimension(nprim) :: r, rm, rp

    ! Res* initialized equal to Res in the block, ghost cells included.
    !$omp do collapse(2)
    do j = 0, n(2) + 1
    do i = 0, n(1) + 1

      Res_Star(:,i,j) = Res(:,i,j)

    end do; end do

    ! Smoothing with Jacobi iterations.
    do jacobi = 1, njacobi

      !$omp do collapse(2) private(r, rm, rp, i1, i2, j1, j2)
      do j = 1, n(2)
      do i = 1, n(1)

        call Jacobi_Stencil(i1, j1, i2, j2, i, j, Dir)
        r (1:np) = Res     (1:np, i,  j )
        rm(1:np) = Res_Star(1:np, i1, j1)
        rp(1:np) = Res_Star(1:np, i2, j2)

        ! Jacobi: r*_New(i) = [ r(i) + e*( r*(i-1) + r*(i+1) ) ] / [ 1 + 2e ]
        ! Smoothing acts only on the flow equations (1:np); the turbulence residuals are left unchanged.
        Res_Star_New(1:np, i, j) = (r(1:np) + beta*(rm(1:np) + rp(1:np))) / (1d0 + 2d0*beta)

      end do; end do

      !$omp do collapse(2)
      do j = 1, n(2)
      do i = 1, n(1)

        Res_Star(1:np, i, j) = Res_Star_New(1:np, i, j)

      end do; end do

    end do ! jacobi

    ! R* after the directional smoothing becomes the new R.
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

      ! 2D form: stencil along i (d=1) or along j (d=2).
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
    ! Local variables
    integer :: ii, i, ig, jg
    integer :: bm, im, jm, fm, bs, is, js

    ! MPI: exchange of the residuals for type-1 connections between different ranks.
    !$omp single
    call exchange_ghost_R(domain)
    !$omp end single

    ! Processing of the BC array — only the local BC entries are iterated over.
    !$omp do private(ii, i, bm, im, jm, fm, bs, is, js, ig, jg)
    do ii = 1, domain%n_local_bc
      i = domain%local_bc_idx(ii)

      ! Preliminary assignments (2D form: only columns 1:2 of guide from Q2D_Parameters_m are used).
      bm = domain%bc(i)%b
      im = domain%bc(i)%i
      jm = domain%bc(i)%j
      fm = domain%bc(i)%f

      ! Ghost-cell coordinates (2D form: i,j indices only).
      ig = im - guide(fm, 1)
      jg = jm - guide(fm, 2)

      select case (domain%bc(i)%type)

        case (101) ! block-to-block connection
          bs = domain%bc(i)%bs
          if (.not. is_local_block(bs)) cycle  ! the inter-rank case is handled by exchange_ghost_R
          is = domain%bc(i)%is
          js = domain%bc(i)%js

          domain%blk(bm)%R(:, ig, jg) = domain%blk(bs)%R(:, is, js)

        case default ! wall, inlet, outlet — zero-gradient ghost cell
          domain%blk(bm)%R(:, ig, jg) = domain%blk(bm)%R(:, im, jm)

      end select
    end do

  end subroutine Ghost_Residual_Extrapolation

end module Q2D_Lib_IRS
