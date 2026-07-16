!>@brief 2D multigrid kernel:
!>       4->1 restriction with surface-weighted cell average;
!>       4-weight bilinear prolongation;
!>       rank-3 fields (nprim, i, j);
!>       Coarse_Grid_Blk on rank-2 Q2D_vector_2D_type;
!>       Coarse_IOfield keeps the rank-4 ORION mesh with placeholder Nk=0
!>       (third coordinate mesh(3,...) = 0).
module Q2D_Lib_Multigrid
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Check_Multigrid, Coarse_Grid, Coarse_IOfield
  public :: fine2coarse_prim, coarse2fine_prim
  public :: Coarse_Cross_Section

contains

  subroutine Check_Multigrid(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_multigrid
    use Q2D_Global_m,         only: ndir
    use ir_precision,         only: str
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: b, d, rap, check

    rap = 2**(obj_multigrid%MGL - 1)
    do b = 1, domain%nb
      check = 0
      do d = 1, ndir
        if (Mod(domain%blk(b)%dim(d), rap) == 0) check = check + 1
      end do
      if (check < ndir) then
        write(*,'(A90)') ' [ERROR] in Check_Multigrid, block: '//trim(str(.true., b))
        stop
      end if
    end do

  end subroutine Check_Multigrid


  subroutine fine2coarse_prim(fPrim, cPrim, fSur, cSur, fDim, cDim)
    use Q2D_Global_m,            only: nprim, np, gc
    use FLINT_Lib_Thermodynamic, only: prim2cons, cons2prim, co_rotot_Rtot
    implicit none
    integer,  intent(in),  dimension(2) :: fDim, cDim
    real(R8), intent(in),  dimension(nprim, 1-gc:fDim(1)+gc, 1-gc:fDim(2)+gc) :: fPrim
    real(R8), intent(out), dimension(nprim, 1-gc:cDim(1)+gc, 1-gc:cDim(2)+gc) :: cPrim
    real(R8), intent(in),  dimension(1-gc:fDim(1)+gc, 1-gc:fDim(2)+gc)        :: fSur
    real(R8), intent(in),  dimension(1-gc:cDim(1)+gc, 1-gc:cDim(2)+gc)        :: cSur
    ! Local variables
    integer  :: i, j, i2, j2, i2d, j2d
    real(R8) :: guess(nprim), rho, Rgas, T0
    real(R8), dimension(nprim, fDim(1), fDim(2)) :: fCons
    real(R8), dimension(nprim, cDim(1), cDim(2)) :: cCons

    ! prim -> cons on the fine grid.
    !$omp do collapse(2)
    do j = 1, fDim(2)
    do i = 1, fDim(1)
      fCons(np+1:nprim, i, j) = fPrim(np+1:nprim, i, j)
      fCons(1:np, i, j)       = prim2cons(fPrim(1:np, i, j))
    end do; end do

    ! 4->1 restriction: 4 fine cells averaged into 1 coarse cell with surface
    ! weights (in 2D the cell surface 'sur' plays the role of the volume).
    !$omp do collapse(2)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2  = 2*i
      j2  = 2*j
      i2d = i2 - 1
      j2d = j2 - 1

      cCons(:, i, j) = fCons(:, i2d, j2d) * fSur(i2d, j2d)  &
                     + fCons(:, i2,  j2d) * fSur(i2,  j2d)  &
                     + fCons(:, i2d, j2 ) * fSur(i2d, j2 )  &
                     + fCons(:, i2,  j2 ) * fSur(i2,  j2 )

    end do; end do

    ! cons -> prim on the coarse grid (Newton-Raphson guess from a fine cell).
    !$omp do collapse(2)
    do j = 1, cDim(2)
    do i = 1, cDim(1)
      cCons(:, i, j) = cCons(:, i, j) / cSur(i, j)
      i2 = 2*i
      j2 = 2*j
      guess = fPrim(:, i2, j2)
      call co_rotot_Rtot(guess(:), rho, Rgas)
      T0 = guess(np) / (rho * Rgas)
      cPrim(np+1:nprim, i, j) = cCons(np+1:nprim, i, j) ! turbulence variables
      cPrim(1:np, i, j)       = cons2prim(cCons(1:np, i, j), T0) ! ( rho U p )
    end do; end do

  end subroutine fine2coarse_prim


  subroutine coarse2fine_prim(fPrim, cPrim, fDim, cDim)
    use Q2D_Global_m, only: nprim, gc
    implicit none
    integer,  intent(in),  dimension(2) :: fDim, cDim
    real(R8), intent(out), dimension(nprim, 1-gc:fDim(1)+gc, 1-gc:fDim(2)+gc) :: fPrim
    real(R8), intent(in),  dimension(nprim, 1-gc:cDim(1)+gc, 1-gc:cDim(2)+gc) :: cPrim
    ! Local variables
    integer  :: i, j, rap, i2, j2, i2d, j2d, im, jm, ip, jp
    integer  :: ii, jj, counter, mask(2), id(4)
    real(R8) :: a1, a2, a3, coeffs(4), interp(nprim)

    ! 2D bilinear weights: 9/16, 3/16, 3/16, 1/16.
    a1 = 9.d0/16.d0
    a2 = 3.d0/16.d0
    a3 = 1.d0/16.d0
    coeffs(1:4) = [a1, a2, a2, a3]
    rap = 2

    !$omp do collapse(2)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2  = rap*i
      j2  = rap*j
      i2d = i2 - (rap - 1)
      j2d = j2 - (rap - 1)

      im = Max(1, i - 1)
      jm = Max(1, j - 1)

      ip = Min(cDim(1), i + 1)
      jp = Min(cDim(2), j + 1)

      id(1:4) = [im, jm, ip, jp]

      counter = 1

      do jj = j2d, j2
      do ii = i2d, i2

        ! 2D mask: 4 fine sub-cells per coarse cell.
        if (counter == 1) mask(1:2) = [1, 2]
        if (counter == 2) mask(1:2) = [3, 2]
        if (counter == 3) mask(1:2) = [1, 4]
        if (counter == 4) mask(1:2) = [3, 4]

        interp = coeffs(1) * cPrim(:, i,           j          ) &
               + coeffs(2) * cPrim(:, id(mask(1)), j          ) &
               + coeffs(3) * cPrim(:, i,           id(mask(2))) &
               + coeffs(4) * cPrim(:, id(mask(1)), id(mask(2)))

        fPrim(:, ii, jj) = interp

        counter = counter + 1

      end do; end do

    end do; end do

  end subroutine coarse2fine_prim


  subroutine Coarse_Grid(Fine, Coarse)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    implicit none
    type(Q2D_domain_type), intent(inout) :: Fine, Coarse
    ! Local variables
    integer :: b

    do b = 1, Coarse%nb ! Loop over the blocks
      call Coarse_Grid_Blk(Fine%blk(b)%node, Coarse%blk(b)%node, &
                           Fine%blk(b)%dim,  Coarse%blk(b)%dim)
    end do

  contains

    subroutine Coarse_Grid_Blk(fNode, cNode, fDim, cDim)
      use Q2D_Base_Types_m, only: Q2D_vector_2D_type
      implicit none
      integer,                  intent(in)                                     :: fDim(2), cDim(2)
      type(Q2D_vector_2D_type), intent(in),  dimension(0:fDim(1), 0:fDim(2))   :: fNode
      type(Q2D_vector_2D_type), intent(out), dimension(0:cDim(1), 0:cDim(2))   :: cNode
      ! Local variables
      integer :: i, j, i2, j2

      !$omp do collapse(2)
      do j = 0, fDim(2), 2
      do i = 0, fDim(1), 2

        i2 = i / 2
        j2 = j / 2

        cNode(i2, j2)%c = fNode(i, j)%c

      end do; end do

    end subroutine Coarse_Grid_Blk

  end subroutine Coarse_Grid


  subroutine Coarse_Cross_Section(Fine, Coarse)
    ! Restriction by nodal sub-sampling of blk%Cross_Section:
    ! A_c(I,J) = A_f(2I, 2J) — mirrors Coarse_Grid_Blk for node%c.
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    implicit none
    type(Q2D_domain_type), intent(in)    :: Fine
    type(Q2D_domain_type), intent(inout) :: Coarse
    ! Local variables
    integer :: b, i, j

    do b = 1, Coarse%nb
      if (.not. Fine%blk(b)%area_variation_active) cycle
      Coarse%blk(b)%area_variation_active = .true.
      !$omp do collapse(2)
      do j = 0, Coarse%blk(b)%dim(2)
      do i = 0, Coarse%blk(b)%dim(1)
        Coarse%blk(b)%Cross_Section(i, j) = Fine%blk(b)%Cross_Section(2*i, 2*j)
      end do; end do
    end do

  end subroutine Coarse_Cross_Section


  subroutine Coarse_IOfield(fIOfield, cIOfield)
    use Lib_ORION_data, only: ORION_data
    use Q2D_Global_m,   only: nprim
    implicit none
    type(ORION_data), intent(in)    :: fIOfield
    type(ORION_data), intent(inout) :: cIOfield
    ! Local variables
    integer :: b, nb, Ni, Nj, Onvar

    Onvar = nprim + 3 ! T, gamma, R

    nb = size(fIOfield%block)
    allocate(cIOfield%block(nb))

    do b = 1, nb
      Ni = fIOfield%block(b)%Ni / 2
      Nj = fIOfield%block(b)%Nj / 2
      cIOfield%block(b)%Ni = Ni
      cIOfield%block(b)%Nj = Nj
      cIOfield%block(b)%Nk = 0  ! rank-4 ORION mesh with placeholder Nk=0.
      ! The ORION mesh mesh(coord, i, j, k) is structurally rank-4; in Q2D the third
      ! coordinate mesh(3, ...) stays 0 and the k dimension is degenerate to 0:0.
      allocate(cIOfield%block(b)%mesh(1:3, 0:Ni, 0:Nj, 0:0))
      allocate(cIOfield%block(b)%vars(1:Onvar, 1:Ni, 1:Nj, 1:1))
    end do

  end subroutine Coarse_IOfield

end module Q2D_Lib_Multigrid
