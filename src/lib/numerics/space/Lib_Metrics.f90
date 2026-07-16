module Q2D_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Base_Types_m,     only: Q2D_tensor_2D_type, Q2D_vector_2D_type
  use Q2D_Advanced_Types_m, only: Q2D_block_type, Q2D_bc_type, Q2D_domain_type
  use Q2D_Global_m,         only: ndir, gc
  use Q2D_Parameters_m,     only: guide
  use Q2D_Lib_Ghost,        only: Is_Wall   ! Is_Wall is defined in Lib_Ghost.

  implicit none
  private

  public :: Check_Mesh_Type
  public :: Compute_Norm_Area
  public :: Compute_Metric_Tensor
  public :: Compute_Yn
  public :: Yn_Connection
  public :: BC_Extrapolate_Metrics
  public :: BC_Symmetry_Metrics
  public :: BC_Connect_Metrics
  public :: BC_Connect_Metrics_Periodic

contains

  subroutine Check_Mesh_Type ( domain )
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain

    ! Q2D supports only 1D and 2D computational domains (no out-of-plane
    ! direction). 2D axisymmetry and Q2D area variation are handled by
    ! Lib_Source_* via Cross_Section / axis_r_ave / area_geom_term_*, not
    ! through an out-of-plane angle.
    if (domain%blk(1)%dim(2) == 1) then
      ndir = 1
    else
      ndir = 2
    end if

  end subroutine Check_Mesh_Type


  subroutine Compute_Norm_Area ( b )
    implicit none
    class(Q2D_block_type), intent(inout) :: b
    ! Locals
    real(R8) :: Ai, snix, sniy, Aj, snjx, snjy
    real(R8) :: d1(2), d2(2), d3(2)
    real(R8) :: snixx, sniyy, snjxx, snjyy
    real(R8) :: scal, signi, signj
    integer  :: i, j, im, jm

    im = b%dim(1)
    jm = b%dim(2)

    ! ---- sign of the unit normals to the i faces -----------------------------
    d1 = b%node(1,1)%c - b%node(1,0)%c
    d3(1) =  d1(2)
    d3(2) = -d1(1)
    snixx = d3(1)
    sniyy = d3(2)

    d1 = 0.5d0*( b%node(1,1)%c + b%node(1,0)%c )
    d2 = 0.5d0*( b%node(0,1)%c + b%node(0,0)%c )
    d3 = d1 - d2
    scal  = d3(1)*snixx + d3(2)*sniyy
    signi = sign(1.d0, scal)

    ! ---- sign of the unit normals to the j faces -----------------------------
    d1 = b%node(0,1)%c - b%node(1,1)%c
    d3(1) =  d1(2)
    d3(2) = -d1(1)
    snjxx = d3(1)
    snjyy = d3(2)

    d1 = 0.5d0*( b%node(1,1)%c + b%node(0,1)%c )
    d2 = 0.5d0*( b%node(1,0)%c + b%node(0,0)%c )
    d3 = d1 - d2
    scal  = d3(1)*snjxx + d3(2)*snjyy
    signj = sign(1.d0, scal)

    !$omp parallel default(none) shared(b, im, jm, signi, signj) &
    !$omp          private(i, j, d1, d3, snix, sniy, Ai, snjx, snjy, Aj)

    !$omp do collapse(2)
    do j = 1, jm
      do i = 0, im
        d1 = b%node(i,j)%c - b%node(i,j-1)%c
        d3(1) =  d1(2)
        d3(2) = -d1(1)
        Ai = sqrt( d3(1)**2 + d3(2)**2 )
        if (Ai == 0d0) then
          snix = 0d0
          sniy = 0d0
        else
          snix = d3(1)/Ai*signi
          sniy = d3(2)/Ai*signi
        end if
        b%dir(1)%f(i,j)%int_l = Ai
        b%dir(1)%f(i,j)%N     = [ snix, sniy ]
      end do
    end do

    !$omp do collapse(2)
    do j = 0, jm
      do i = 1, im
        d1 = b%node(i-1,j)%c - b%node(i,j)%c
        d3(1) =  d1(2)
        d3(2) = -d1(1)
        Aj = sqrt( d3(1)**2 + d3(2)**2 )
        if (Aj == 0d0) then
          snjx = 0d0
          snjy = 0d0
        else
          snjx = d3(1)/Aj*signj
          snjy = d3(2)/Aj*signj
        end if
        b%dir(2)%f(i,j)%int_l = Aj
        b%dir(2)%f(i,j)%N     = [ snjx, snjy ]
      end do
    end do
    !$omp end parallel

  end subroutine Compute_Norm_Area


  subroutine Compute_Metric_Tensor ( N1, N2, N3, N4, M, dl, sur )
    implicit none
    type(Q2D_vector_2D_type), intent(in)  :: N1, N2, N3, N4
    type(Q2D_tensor_2D_type), intent(out) :: M
    type(Q2D_vector_2D_type), intent(out) :: dl
    real(R8),                 intent(out) :: sur
    ! Locals
    integer  :: h
    real(R8) :: det, A(2,2), cofactor(2,2)
    real(R8) :: vx(4), vy(4)

    ! Discrete Jacobian A_ij = average of the differences between opposite sides.
    ! Metric tensor M obtained from cofactors/det: the formula produces A^{-T},
    ! which coincides with A^{-1} only for symmetric A (solver convention).
    A(1,:) = N4%c - N2%c + N3%c - N1%c
    A(2,:) = N4%c - N3%c + N2%c - N1%c
    A = 0.5d0 * A

    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    if ( abs(det) == 0d0 ) then
      write(*,'(A)') ' [ERROR] Q2D_Lib_Metrics::Compute_Metric_Tensor: jacobian det=0.'
      error stop
    end if

    cofactor(1,1) =  A(2,2)
    cofactor(1,2) = -A(2,1)
    cofactor(2,1) = -A(1,2)
    cofactor(2,2) =  A(1,1)
    M%c = cofactor / det

    do h = 1, 2
      dl%c(h) = sqrt( A(h,1)**2 + A(h,2)**2 )
    end do

    vx(1) = N1%c(1) ; vy(1) = N1%c(2)
    vx(2) = N2%c(1) ; vy(2) = N2%c(2)
    vx(3) = N3%c(1) ; vy(3) = N3%c(2)
    vx(4) = N4%c(1) ; vy(4) = N4%c(2)

    sur = tsurface(vx, vy, 1, 2, 3, 4)

    if ( sur <= 0d0 ) then
      write(*,'(A)') ' [ERROR] Q2D_Lib_Metrics::Compute_Metric_Tensor: negative cell surface.'
      error stop
    end if

  contains

    pure function tsurface(vx, vy, i1, i2, i3, i4) result(surface)
      implicit none
      real(R8), intent(in) :: vx(4), vy(4)
      integer,  intent(in) :: i1, i2, i3, i4
      real(R8)             :: surface
      real(R8)             :: d1(2), d2(2), cross

      d1(1) = vx(i4) - vx(i1)
      d1(2) = vy(i4) - vy(i1)
      d2(1) = vx(i3) - vx(i2)
      d2(2) = vy(i3) - vy(i2)
      cross = d1(1)*d2(2) - d1(2)*d2(1)
      surface = 0.5d0 * abs(cross)
    end function tsurface

  end subroutine Compute_Metric_Tensor


  subroutine Compute_Yn ( blk, bc, nb, nbound )
    implicit none
    integer,                                intent(in)    :: nb, nbound
    type(Q2D_block_type), dimension(:),     intent(inout) :: blk
    type(Q2D_bc_type),    dimension(nbound), intent(in)   :: bc
    ! Locals
    real(R8), allocatable :: facex(:), facey(:)
    real(R8) :: x1, x2, y1, y2
    integer  :: Fm, Bm, Im, Jm, l, b, nwall
    integer  :: i1, j1, i2, j2

    allocate(facex(nbound), facey(nbound))
    nwall = 0

    do l = 1, nbound
      if ( Is_Wall(bc(l)%type) ) then
        Bm = bc(l)%b
        Im = bc(l)%i
        Jm = bc(l)%j
        Fm = bc(l)%f
        nwall = nwall + 1

        select case (Fm)
          case(1:2)
            i1 = Im - mod(Fm,2) ; j1 = Jm - 1
            i2 = Im - mod(Fm,2) ; j2 = Jm
          case(3:4)
            i1 = Im - 1 ; j1 = Jm - mod(Fm,2)
            i2 = Im     ; j2 = Jm - mod(Fm,2)
        end select

        x1 = blk(Bm)%node(i1,j1)%c(1)
        x2 = blk(Bm)%node(i2,j2)%c(1)
        y1 = blk(Bm)%node(i1,j1)%c(2)
        y2 = blk(Bm)%node(i2,j2)%c(2)

        facex(nwall) = 0.5d0*(x1 + x2)
        facey(nwall) = 0.5d0*(y1 + y2)
      end if
    end do

    do b = 1, nb
      call Wall_Distance_Blk ( blk(b), facex, facey, nwall )
    end do

    deallocate(facex, facey)

  contains

    subroutine Wall_Distance_Blk ( blk, fx, fy, nwall )
      implicit none
      type(Q2D_block_type),       intent(inout) :: blk
      integer,                    intent(in)    :: nwall
      real(R8), dimension(nwall), intent(in)    :: fx, fy
      ! Locals
      integer  :: l, i, j
      real(R8) :: center(2), dummy

      !$omp parallel default(none) shared(blk, fx, fy, nwall) &
      !$omp          private(i, j, l, center, dummy)
      !$omp do collapse(2)
      do j = 1, blk%dim(2)
        do i = 1, blk%dim(1)
          center = 0.25d0*( blk%node(i  ,j  )%c + blk%node(i  ,j-1)%c + &
                            blk%node(i-1,j-1)%c + blk%node(i-1,j  )%c )
          blk%yn(i,j) = 1d8
          do l = 1, nwall
            dummy = sqrt( (center(1) - fx(l))**2 + (center(2) - fy(l))**2 )
            blk%yn(i,j) = min(dummy, blk%yn(i,j))
          end do
        end do
      end do
      !$omp end parallel
    end subroutine Wall_Distance_Blk

  end subroutine Compute_Yn


  subroutine Yn_Connection ( Im, Jm, Fm, blkm, Is, Js, blks )
    implicit none
    integer,              intent(in)    :: Im, Jm, Fm, Is, Js
    type(Q2D_block_type), intent(inout) :: blkm
    type(Q2D_block_type), intent(inout) :: blks
    ! Locals
    integer :: Ig, Jg

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    blkm%yn(Ig,Jg) = blks%yn(Is,Js)

  end subroutine Yn_Connection


  subroutine BC_Extrapolate_Metrics ( Im, Jm, Fm, blk, Mg, dlg, surg )
    implicit none
    integer,                  intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type),     intent(inout) :: blk
    type(Q2D_tensor_2D_type), intent(out)   :: Mg(2)
    type(Q2D_vector_2D_type), intent(out)   :: dlg(2)
    real(R8),                 intent(out)   :: surg(2)
    ! Locals
    integer :: g
    type(Q2D_vector_2D_type), dimension(-gc:gc) :: N1, N2, N3, N4

    select case (Fm)
      case(1)
        do g = 0, gc
          N1(g)%c = blk%node(Im+g-1, Jm-1)%c
          N2(g)%c = blk%node(Im+g-1, Jm  )%c
        end do
        do g = -1, -gc, -1
          N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
          N2(g)%c = 3d0*N2(g+1)%c - 3d0*N2(g+2)%c + N2(g+3)%c
          call Compute_Metric_Tensor( N1(g), N2(g), N1(g+1), N2(g+1), &
                                      blk%M(Im+g, Jm), blk%dl(Im+g, Jm), blk%sur(Im+g, Jm) )
          Mg(-g)   = blk%M  (Im+g, Jm)
          dlg(-g)  = blk%dl (Im+g, Jm)
          surg(-g) = blk%sur(Im+g, Jm)
        end do

      case(2)
        do g = 0, -gc, -1
          N3(g)%c = blk%node(Im+g, Jm-1)%c
          N4(g)%c = blk%node(Im+g, Jm  )%c
        end do
        do g = 1, gc
          N3(g)%c = 3d0*N3(g-1)%c - 3d0*N3(g-2)%c + N3(g-3)%c
          N4(g)%c = 3d0*N4(g-1)%c - 3d0*N4(g-2)%c + N4(g-3)%c
          call Compute_Metric_Tensor( N3(g-1), N4(g-1), N3(g), N4(g), &
                                      blk%M(Im+g, Jm), blk%dl(Im+g, Jm), blk%sur(Im+g, Jm) )
          Mg(g)   = blk%M  (Im+g, Jm)
          dlg(g)  = blk%dl (Im+g, Jm)
          surg(g) = blk%sur(Im+g, Jm)
        end do

      case(3)
        ! dim(2)>2: 3-node quadratic extrapolation; thin blocks in J
        ! (dim(2)<=2, e.g. 1D mesh with J=2 nodes) use the 2-node linear
        ! fallback. The quadratic branch reads node(Im-1,Jm+1), out of the bound
        ! node(0:nj) when nj=1: the guard is therefore necessary.
        if ( blk%dim(2) > 2 ) then
          do g = 0, gc
            N1(g)%c = blk%node(Im-1, Jm+g-1)%c
            N3(g)%c = blk%node(Im  , Jm+g-1)%c
          end do
          do g = -1, -gc, -1
            N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
            N3(g)%c = 3d0*N3(g+1)%c - 3d0*N3(g+2)%c + N3(g+3)%c
            call Compute_Metric_Tensor( N1(g), N1(g+1), N3(g), N3(g+1), &
                                        blk%M(Im, Jm+g), blk%dl(Im, Jm+g), blk%sur(Im, Jm+g) )
            Mg(-g)   = blk%M  (Im, Jm+g)
            dlg(-g)  = blk%dl (Im, Jm+g)
            surg(-g) = blk%sur(Im, Jm+g)
          end do
        else
          do g = 0, 1
            N1(g)%c = blk%node(Im-1, Jm+g-1)%c
            N3(g)%c = blk%node(Im  , Jm+g-1)%c
          end do
          do g = -1, -gc, -1
            N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
            N3(g)%c = 2d0*N3(g+1)%c - N3(g+2)%c
            call Compute_Metric_Tensor( N1(g), N1(g+1), N3(g), N3(g+1), &
                                        blk%M(Im, Jm+g), blk%dl(Im, Jm+g), blk%sur(Im, Jm+g) )
            Mg(-g)   = blk%M  (Im, Jm+g)
            dlg(-g)  = blk%dl (Im, Jm+g)
            surg(-g) = blk%sur(Im, Jm+g)
          end do
        end if

      case(4)
        ! See case(3): same guard with 2-node linear fallback for thin blocks in J.
        if ( blk%dim(2) > 2 ) then
          do g = 0, -gc, -1
            N2(g)%c = blk%node(Im-1, Jm+g)%c
            N4(g)%c = blk%node(Im  , Jm+g)%c
          end do
          do g = 1, gc
            N2(g)%c = 3d0*N2(g-1)%c - 3d0*N2(g-2)%c + N2(g-3)%c
            N4(g)%c = 3d0*N4(g-1)%c - 3d0*N4(g-2)%c + N4(g-3)%c
            call Compute_Metric_Tensor( N2(g-1), N2(g), N4(g-1), N4(g), &
                                        blk%M(Im, Jm+g), blk%dl(Im, Jm+g), blk%sur(Im, Jm+g) )
            Mg(g)   = blk%M  (Im, Jm+g)
            dlg(g)  = blk%dl (Im, Jm+g)
            surg(g) = blk%sur(Im, Jm+g)
          end do
        else
          do g = 0, -1, -1
            N2(g)%c = blk%node(Im-1, Jm+g)%c
            N4(g)%c = blk%node(Im  , Jm+g)%c
          end do
          do g = 1, gc
            N2(g)%c = 2d0*N2(g-1)%c - N2(g-2)%c
            N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
            call Compute_Metric_Tensor( N2(g-1), N2(g), N4(g-1), N4(g), &
                                        blk%M(Im, Jm+g), blk%dl(Im, Jm+g), blk%sur(Im, Jm+g) )
            Mg(g)   = blk%M  (Im, Jm+g)
            dlg(g)  = blk%dl (Im, Jm+g)
            surg(g) = blk%sur(Im, Jm+g)
          end do
        end if
    end select

  end subroutine BC_Extrapolate_Metrics


  subroutine BC_Symmetry_Metrics ( Im, Jm, Fm, blk, Mg, dlg, surg )
    implicit none
    integer,                  intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type),     intent(inout) :: blk
    type(Q2D_tensor_2D_type), intent(out)   :: Mg(2)
    type(Q2D_vector_2D_type), intent(out)   :: dlg(2)
    real(R8),                 intent(out)   :: surg(2)
    ! Locals
    integer :: g

    do g = 1, gc
      select case (Fm)
        case(1)
          blk%M  (Im-g, Jm) = blk%M  (Im+g-1, Jm)
          blk%dl (Im-g, Jm) = blk%dl (Im+g-1, Jm)
          blk%sur(Im-g, Jm) = blk%sur(Im+g-1, Jm)
          Mg(g)   = blk%M  (Im-g, Jm)
          dlg(g)  = blk%dl (Im-g, Jm)
          surg(g) = blk%sur(Im-g, Jm)
        case(2)
          blk%M  (Im+g, Jm) = blk%M  (Im-g+1, Jm)
          blk%dl (Im+g, Jm) = blk%dl (Im-g+1, Jm)
          blk%sur(Im+g, Jm) = blk%sur(Im-g+1, Jm)
          Mg(g)   = blk%M  (Im+g, Jm)
          dlg(g)  = blk%dl (Im+g, Jm)
          surg(g) = blk%sur(Im+g, Jm)
        case(3)
          blk%M  (Im, Jm-g) = blk%M  (Im, Jm+g-1)
          blk%dl (Im, Jm-g) = blk%dl (Im, Jm+g-1)
          blk%sur(Im, Jm-g) = blk%sur(Im, Jm+g-1)
          Mg(g)   = blk%M  (Im, Jm-g)
          dlg(g)  = blk%dl (Im, Jm-g)
          surg(g) = blk%sur(Im, Jm-g)
        case(4)
          blk%M  (Im, Jm+g) = blk%M  (Im, Jm-g+1)
          blk%dl (Im, Jm+g) = blk%dl (Im, Jm-g+1)
          blk%sur(Im, Jm+g) = blk%sur(Im, Jm-g+1)
          Mg(g)   = blk%M  (Im, Jm+g)
          dlg(g)  = blk%dl (Im, Jm+g)
          surg(g) = blk%sur(Im, Jm+g)
      end select
    end do

  end subroutine BC_Symmetry_Metrics


  subroutine BC_Connect_Metrics ( Im, Jm, Fm, blkm, Is, Js, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, surg )
    implicit none
    integer,                  intent(in)    :: Im, Jm, Fm, Is, Js, Fs, d11s, d12s, d21s, d22s
    type(Q2D_block_type),     intent(in)    :: blks
    type(Q2D_block_type),     intent(inout) :: blkm
    type(Q2D_tensor_2D_type), intent(out)   :: Mg(2)
    type(Q2D_vector_2D_type), intent(out)   :: dlg(2)
    real(R8),                 intent(out)   :: surg(2)
    ! Locals
    integer :: g, Is1, Js1, II, JJ, III, JJJ
    integer :: guides(6,3), guides2(6,3)
    type(Q2D_vector_2D_type), dimension(0:gc) :: N1, N2, N3, N4

    ! Stencil offsets for the source block (donor side): 6-face, 3-column tables;
    ! in 2D only columns 1 and 2 are used (column 3 is unused).
    guides(1,1)= 0  ; guides(1,2)= 0  ; guides(1,3)= 0
    guides(2,1)=-1  ; guides(2,2)= 0  ; guides(2,3)= 0
    guides(3,1)= 0  ; guides(3,2)= 0  ; guides(3,3)= 0
    guides(4,1)= 0  ; guides(4,2)=-1  ; guides(4,3)= 0
    guides(5,1)= 0  ; guides(5,2)= 0  ; guides(5,3)= 0
    guides(6,1)= 0  ; guides(6,2)= 0  ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0
    guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0
    guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1
    guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2

    select case (Fm)
      case(1)
        N1(0)%c = blkm%node(Im-1, Jm-1)%c
        N2(0)%c = blkm%node(Im-1, Jm  )%c
      case(2)
        N3(0)%c = blkm%node(Im  , Jm-1)%c
        N4(0)%c = blkm%node(Im  , Jm  )%c
      case(3)
        N1(0)%c = blkm%node(Im-1, Jm-1)%c
        N3(0)%c = blkm%node(Im  , Jm-1)%c
      case(4)
        N2(0)%c = blkm%node(Im-1, Jm  )%c
        N4(0)%c = blkm%node(Im  , Jm  )%c
    end select

    do g = 1, gc
      Is1 = Is + (gc-g)*guides(Fs,1) + (g-1)*guides2(Fs,1)
      Js1 = Js + (gc-g)*guides(Fs,2) + (g-1)*guides2(Fs,2)

      select case (Fs)
        case(1:2)
          II = Is1
          JJ = (2*Js1 + d11s + d21s)/2
        case(3:4)
          II = (2*Is1 + d11s + d21s)/2
          JJ = Js1
      end select

      select case (Fm)
        case(1) ; N2(g)%c = blks%node(II,JJ)%c
        case(2) ; N4(g)%c = blks%node(II,JJ)%c
        case(3) ; N3(g)%c = blks%node(II,JJ)%c
        case(4) ; N4(g)%c = blks%node(II,JJ)%c
      end select

      ! Update of the diagonal corner node (i1-1).
      select case (Fs)
        case(1:2)
          III = II
          JJJ = JJ - d11s
        case(3:4)
          III = II - d11s
          JJJ = JJ
      end select

      select case (Fm)
        case(1) ; N1(g)%c = blks%node(III,JJJ)%c
        case(2) ; N3(g)%c = blks%node(III,JJJ)%c
        case(3) ; N1(g)%c = blks%node(III,JJJ)%c
        case(4) ; N2(g)%c = blks%node(III,JJJ)%c
      end select

      select case (Fm)
        case(1)
          call Compute_Metric_Tensor( N1(g), N2(g), N1(g-1), N2(g-1), &
                                      blkm%M(Im-g, Jm), blkm%dl(Im-g, Jm), blkm%sur(Im-g, Jm) )
          Mg(g)   = blkm%M  (Im-g, Jm)
          dlg(g)  = blkm%dl (Im-g, Jm)
          surg(g) = blkm%sur(Im-g, Jm)
        case(2)
          call Compute_Metric_Tensor( N3(g-1), N4(g-1), N3(g), N4(g), &
                                      blkm%M(Im+g, Jm), blkm%dl(Im+g, Jm), blkm%sur(Im+g, Jm) )
          Mg(g)   = blkm%M  (Im+g, Jm)
          dlg(g)  = blkm%dl (Im+g, Jm)
          surg(g) = blkm%sur(Im+g, Jm)
        case(3)
          call Compute_Metric_Tensor( N1(g), N1(g-1), N3(g), N3(g-1), &
                                      blkm%M(Im, Jm-g), blkm%dl(Im, Jm-g), blkm%sur(Im, Jm-g) )
          Mg(g)   = blkm%M  (Im, Jm-g)
          dlg(g)  = blkm%dl (Im, Jm-g)
          surg(g) = blkm%sur(Im, Jm-g)
        case(4)
          call Compute_Metric_Tensor( N2(g-1), N2(g), N4(g-1), N4(g), &
                                      blkm%M(Im, Jm+g), blkm%dl(Im, Jm+g), blkm%sur(Im, Jm+g) )
          Mg(g)   = blkm%M  (Im, Jm+g)
          dlg(g)  = blkm%dl (Im, Jm+g)
          surg(g) = blkm%sur(Im, Jm+g)
      end select
    end do

  end subroutine BC_Connect_Metrics


  subroutine BC_Connect_Metrics_Periodic ( Im, Jm, Fm, blkm, Is, Js, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, surg )
    ! Ghost metric for a translational periodic connection (BC 201).
    !
    ! WHY a twin of BC_Connect_Metrics rather than its reuse:
    ! BC_Connect_Metrics reads the donor nodes blks%node as ABSOLUTE coordinates,
    ! correct only when the donor is the physically adjacent block (101). For a
    ! periodic 201 seam the donor is the same edge repeated at a constant
    ! Cartesian OFFSET, so its absolute coordinates are translated; mixing them
    ! with the local corner nodes would build a geometrically wrong (folded)
    ! ghost cell. Here the donor nodes are brought back into the LOCAL frame with
    ! a rigid translation T before forming the quadrilateral.
    !
    ! T = (local seam node) - (corresponding donor node), where the corresponding
    ! donor node uses the SAME recurrence Is1/Js1 + (II,JJ) remapping of the g>=1
    ! reads, evaluated at the seam (step g=0).
    ! Compute_Metric_Tensor depends only on the DIFFERENCES between nodes, so
    ! adding the same T to all donor nodes yields the EXACT metric of the periodic
    ! neighbor (not the cubic-extrapolation approximation of
    ! BC_Extrapolate_Metrics). For an adjacent seam (101) T==0 and the routine
    ! reduces to BC_Connect_Metrics.
    !
    ! The body is a copy of BC_Connect_Metrics with: (1) T computed before the
    ! g loop from the pair of nodes corresponding to g=0; (2) +T%c added to
    ! every read of blks%node.
    implicit none
    integer,                  intent(in)    :: Im, Jm, Fm, Is, Js, Fs, d11s, d12s, d21s, d22s
    type(Q2D_block_type),     intent(in)    :: blks
    type(Q2D_block_type),     intent(inout) :: blkm
    type(Q2D_tensor_2D_type), intent(out)   :: Mg(2)
    type(Q2D_vector_2D_type), intent(out)   :: dlg(2)
    real(R8),                 intent(out)   :: surg(2)
    ! Locals
    integer :: g, Is1, Js1, II, JJ, III, JJJ
    integer :: II0, JJ0, III0, JJJ0
    integer :: guides(6,3), guides2(6,3)
    type(Q2D_vector_2D_type), dimension(0:gc) :: N1, N2, N3, N4
    type(Q2D_vector_2D_type) :: T, T2

    guides(1,1)= 0  ; guides(1,2)= 0  ; guides(1,3)= 0
    guides(2,1)=-1  ; guides(2,2)= 0  ; guides(2,3)= 0
    guides(3,1)= 0  ; guides(3,2)= 0  ; guides(3,3)= 0
    guides(4,1)= 0  ; guides(4,2)=-1  ; guides(4,3)= 0
    guides(5,1)= 0  ; guides(5,2)= 0  ; guides(5,3)= 0
    guides(6,1)= 0  ; guides(6,2)= 0  ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0
    guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0
    guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1
    guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2

    select case (Fm)
      case(1)
        N1(0)%c = blkm%node(Im-1, Jm-1)%c
        N2(0)%c = blkm%node(Im-1, Jm  )%c
      case(2)
        N3(0)%c = blkm%node(Im  , Jm-1)%c
        N4(0)%c = blkm%node(Im  , Jm  )%c
      case(3)
        N1(0)%c = blkm%node(Im-1, Jm-1)%c
        N3(0)%c = blkm%node(Im  , Jm-1)%c
      case(4)
        N2(0)%c = blkm%node(Im-1, Jm  )%c
        N4(0)%c = blkm%node(Im  , Jm  )%c
    end select

    ! Rigid periodic offset T from the pair of nodes corresponding to g=0 (seam).
    ! The g=0 step of the same recurrence locates the donor node coincident with the
    ! local seam node; the diagonal step locates the second edge of the seam. T must
    ! come out identical from both (rigid translation) — checked below to intercept a
    ! non-conforming / rotational seam.
    Is1 = Is + gc*guides(Fs,1) - guides2(Fs,1)
    Js1 = Js + gc*guides(Fs,2) - guides2(Fs,2)
    select case (Fs)
      case(1:2)
        II0 = Is1
        JJ0 = (2*Js1 + d11s + d21s)/2
      case(3:4)
        II0 = (2*Is1 + d11s + d21s)/2
        JJ0 = Js1
    end select
    select case (Fs)
      case(1:2)
        III0 = II0
        JJJ0 = JJ0 - d11s
      case(3:4)
        III0 = II0 - d11s
        JJJ0 = JJ0
    end select

    select case (Fm)
      case(1) ; T%c = N2(0)%c - blks%node(II0,JJ0)%c ; T2%c = N1(0)%c - blks%node(III0,JJJ0)%c
      case(2) ; T%c = N4(0)%c - blks%node(II0,JJ0)%c ; T2%c = N3(0)%c - blks%node(III0,JJJ0)%c
      case(3) ; T%c = N3(0)%c - blks%node(II0,JJ0)%c ; T2%c = N1(0)%c - blks%node(III0,JJJ0)%c
      case(4) ; T%c = N4(0)%c - blks%node(II0,JJ0)%c ; T2%c = N2(0)%c - blks%node(III0,JJJ0)%c
    end select

    if ( maxval(abs(T%c - T2%c)) > 1d-10 * (1d0 + maxval(abs(T%c))) ) then
      write(*,'(A)') ' [ERROR] Q2D_Lib_Metrics::BC_Connect_Metrics_Periodic: non-rigid periodic seam'
      write(*,'(A)') '         (BC 201 must be a translational periodic connection; offset is not a constant rigid translation).'
      error stop
    end if

    do g = 1, gc
      Is1 = Is + (gc-g)*guides(Fs,1) + (g-1)*guides2(Fs,1)
      Js1 = Js + (gc-g)*guides(Fs,2) + (g-1)*guides2(Fs,2)

      select case (Fs)
        case(1:2)
          II = Is1
          JJ = (2*Js1 + d11s + d21s)/2
        case(3:4)
          II = (2*Is1 + d11s + d21s)/2
          JJ = Js1
      end select

      select case (Fm)
        case(1) ; N2(g)%c = blks%node(II,JJ)%c + T%c
        case(2) ; N4(g)%c = blks%node(II,JJ)%c + T%c
        case(3) ; N3(g)%c = blks%node(II,JJ)%c + T%c
        case(4) ; N4(g)%c = blks%node(II,JJ)%c + T%c
      end select

      ! Update of the diagonal corner node (i1-1).
      select case (Fs)
        case(1:2)
          III = II
          JJJ = JJ - d11s
        case(3:4)
          III = II - d11s
          JJJ = JJ
      end select

      select case (Fm)
        case(1) ; N1(g)%c = blks%node(III,JJJ)%c + T%c
        case(2) ; N3(g)%c = blks%node(III,JJJ)%c + T%c
        case(3) ; N1(g)%c = blks%node(III,JJJ)%c + T%c
        case(4) ; N2(g)%c = blks%node(III,JJJ)%c + T%c
      end select

      select case (Fm)
        case(1)
          call Compute_Metric_Tensor( N1(g), N2(g), N1(g-1), N2(g-1), &
                                      blkm%M(Im-g, Jm), blkm%dl(Im-g, Jm), blkm%sur(Im-g, Jm) )
          Mg(g)   = blkm%M  (Im-g, Jm)
          dlg(g)  = blkm%dl (Im-g, Jm)
          surg(g) = blkm%sur(Im-g, Jm)
        case(2)
          call Compute_Metric_Tensor( N3(g-1), N4(g-1), N3(g), N4(g), &
                                      blkm%M(Im+g, Jm), blkm%dl(Im+g, Jm), blkm%sur(Im+g, Jm) )
          Mg(g)   = blkm%M  (Im+g, Jm)
          dlg(g)  = blkm%dl (Im+g, Jm)
          surg(g) = blkm%sur(Im+g, Jm)
        case(3)
          call Compute_Metric_Tensor( N1(g), N1(g-1), N3(g), N3(g-1), &
                                      blkm%M(Im, Jm-g), blkm%dl(Im, Jm-g), blkm%sur(Im, Jm-g) )
          Mg(g)   = blkm%M  (Im, Jm-g)
          dlg(g)  = blkm%dl (Im, Jm-g)
          surg(g) = blkm%sur(Im, Jm-g)
        case(4)
          call Compute_Metric_Tensor( N2(g-1), N2(g), N4(g-1), N4(g), &
                                      blkm%M(Im, Jm+g), blkm%dl(Im, Jm+g), blkm%sur(Im, Jm+g) )
          Mg(g)   = blkm%M  (Im, Jm+g)
          dlg(g)  = blkm%dl (Im, Jm+g)
          surg(g) = blkm%sur(Im, Jm+g)
      end select
    end do

  end subroutine BC_Connect_Metrics_Periodic


end module Q2D_Lib_Metrics
