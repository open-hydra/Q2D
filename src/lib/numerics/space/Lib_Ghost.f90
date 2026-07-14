module Q2D_Lib_Ghost
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m, only: Q2D_block_type, Q2D_bc_type, Q2D_domain_type
  use Q2D_Global_m,         only: nprim, nu, nv, np, nsc, gc
  use Q2D_Parameters_m,     only: guide

  implicit none
  private

  public :: Fill_Ghost_Cell
  public :: Ghost_Wall_Extrapolation
  public :: Is_Wall   ! relocated from Q2D_Lib_Metrics.

contains

  subroutine Fill_Ghost_Cell ( domain )
    use Q2D_Mod_MPI,           only: is_local_block, mpi_size_
    use Q2D_Mod_GhostExchange, only: exchange_ghost_P_post_recv, exchange_ghost_P_pack, &
                                     exchange_ghost_P_post_send, exchange_ghost_P_wait_unpack, &
                                     exchange_ghost_P_wait_send, &
                                     Ghost_Interrank, exchange_ghost_Pg, ghost_sched
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, Bm, Im, Jm, Fm, Bs, Is, Js, Fs, d11s, d12s, d21s, d22s
    integer :: fg

    ! 2D-reduced: 4 faces, no Km/Ks.
    ! MPI exchange interleaved with local computation: pack/post send during
    ! local processing, then wait/unpack inter-rank entries.
    !
    ! NB: this routine does NOT open `!$omp parallel`. The parallel region is
    ! opened by the callers (Mod_Explicit in the main loop, Wrap_Setup
    ! single-threaded setup context); the inner `!$omp single`/`!$omp do` act as
    ! sequential no-ops when called outside an `!$omp parallel`.

    ! MPI: post persistent receives, pack buffer in parallel, then post sends
    !$omp single
    call exchange_ghost_P_post_recv(domain)
    !$omp end single

    ! Pack send buffer in parallel over face groups
    !$omp do schedule(static) private(fg)
    do fg = 1, ghost_sched%n_send_faces
      call exchange_ghost_P_pack(domain, fg, fg)
    end do

    ! Post sends (must wait for all packing to complete — implicit barrier from !$omp do)
    !$omp single
    call exchange_ghost_P_post_send(domain)
    !$omp end single nowait

    ! Process LOCAL BC entries while MPI communication is in flight.
    ! Uses pre-filtered local_bc_idx to avoid scanning all nbound entries.
    !$omp do schedule(dynamic) private(ii, i, Bm, Im, Jm, Fm, Bs, Is, Js, Fs)
    do ii = 1, domain%n_local_bc
      i  = domain%local_bc_idx(ii)
      Bm = domain%bc(i)%b
      Im = domain%bc(i)%i
      Jm = domain%bc(i)%j
      Fm = domain%bc(i)%f
      select case (domain%bc(i)%type)
        case(101, 201) ! block connection (101) + translational/rotational periodic (201)
          Bs = domain%bc(i)%bs
          if (.not. is_local_block(Bs)) cycle  ! inter-rank handled after MPI completes
          Is = domain%bc(i)%is
          Js = domain%bc(i)%js
          Fs = domain%bc(i)%fs
          call Ghost_Connection ( Im, Jm, Fm, domain%blk(Bm), &
                                  Is, Js, Fs, domain%blk(Bs), &
                                  domain%bc(i)%Pg )
        case(300) ! symmetry
          call Ghost_Symmetry ( Im, Jm, Fm, domain%blk(Bm) )
        case(0, 401:407, 420) ! inflow / outflow / extrapolation: zero gradient
          call Ghost_ZG_Extrapolate ( Im, Jm, Fm, domain%blk(Bm) )
        case(102) ! chimera
          call Ghost_Chimera ( domain%nb, domain%blk, domain%bc(i) )
        case default ! 2nd order extrapolation
          call Ghost_Extrapolate ( Im, Jm, Fm, domain%blk(Bm) )
      end select
    end do

    ! Compute Pg(:,3:4) for type-1 connections where source block is local.
    ! This only reads local blk%P data, so it can run before the MPI exchange.
    !$omp do schedule(dynamic) private(ii, i, Bs, Is, Js, Fs, d11s, d12s, d21s, d22s)
    do ii = 1, domain%n_local_bs
      i    = domain%local_bs_idx(ii)
      Bs   = domain%bc(i)%bs
      Is   = domain%bc(i)%is
      Js   = domain%bc(i)%js
      Fs   = domain%bc(i)%fs
      d11s = domain%bc(i)%d11
      d12s = domain%bc(i)%d12
      d21s = domain%bc(i)%d21
      d22s = domain%bc(i)%d22
      call Fill_BC_Ghost_Connection ( Is, Js, Fs, d11s, d12s, d21s, d22s, &
                                      domain%blk(Bs), domain%bc(i)%Pg )
    end do

    ! Compute Pg for chimera (102) connections where Bm is local.
    !$omp do schedule(dynamic) private(ii, i, Bm, Im, Jm, Fm)
    do ii = 1, domain%n_local_bc
      i  = domain%local_bc_idx(ii)
      Bm = domain%bc(i)%b
      Im = domain%bc(i)%i
      Jm = domain%bc(i)%j
      Fm = domain%bc(i)%f
      select case (domain%bc(i)%type)
        case(102) ! chimera
          call Fill_BC_Ghost_Chimera ( Im, Jm, Fm, domain%blk(Bm), domain%bc(i)%Pg )
      end select
    end do

    ! MPI: wait for P receives to complete
    !$omp single
    call exchange_ghost_P_wait_unpack(domain)
    !$omp end single

    ! Process INTER-RANK type-1 entries (Bm local, Bs remote).
    ! Ghost_Interrank fills Pg(:,1:gc) ghost layer using P data received via MPI.
    !$omp do schedule(dynamic) private(ii, i, Bm, Im, Jm, Fm, Bs)
    do ii = 1, domain%n_local_bc
      i  = domain%local_bc_idx(ii)
      if (domain%bc(i)%type /= 101 .and. domain%bc(i)%type /= 201) cycle  ! 201 cross-rank halo = pure copy like 101 (§1 div.)
      Bs = domain%bc(i)%bs
      if (is_local_block(Bs)) cycle  ! already processed above
      Bm = domain%bc(i)%b
      Im = domain%bc(i)%i
      Jm = domain%bc(i)%j
      Fm = domain%bc(i)%f
      call Ghost_Interrank(Im, Jm, Fm, domain%blk(Bm), domain%bc(i)%Pg)
    end do

    ! Wait for P sends to complete before reusing buffers; exchange Pg(:,3:4) stencil.
    !$omp single
    call exchange_ghost_P_wait_send(domain)
    call exchange_ghost_Pg(domain)
    !$omp end single

  end subroutine Fill_Ghost_Cell


  subroutine Ghost_Connection ( Im, Jm, Fm, blkm, Is, Js, Fs, blks, Pg )
    implicit none
    integer,              intent(in)    :: Im, Jm, Fm, Is, Js, Fs
    type(Q2D_block_type), intent(in)    :: blks
    type(Q2D_block_type), intent(inout) :: blkm
    real(R8),             intent(inout) :: Pg(:,:)
    ! Local
    integer :: g, Ig, Jg

    do g = 1, gc
      Ig = Im - guide(Fm,1)*g
      Jg = Jm - guide(Fm,2)*g
      blkm%P(:,Ig,Jg) = blks%P(:, Is + guide(Fs,1)*(g-1), &
                                  Js + guide(Fs,2)*(g-1))
      Pg(:,g) = blkm%P(:,Ig,Jg)
    end do

  end subroutine Ghost_Connection


  subroutine Ghost_Symmetry ( Im, Jm, Fm, blk )
    implicit none
    integer,              intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type), intent(inout) :: blk
    ! Local
    integer  :: Ig, Jg, Ic, Jc
    real(R8) :: normal(2), vel(2), udotn

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)

    if (Fm <= 2) then
      Ic = Im - mod(Fm,2)
      normal = blk%dir(1)%f(Ic, Jm)%N
    else
      Jc = Jm - mod(Fm,2)
      normal = blk%dir(2)%f(Im, Jc)%N
    end if

    vel = blk%P(nu:nv, Im, Jm)
    udotn = dot_product(vel, normal)
    vel = vel - 2d0*udotn*normal

    blk%P(nu:nv,    Ig, Jg) = vel
    blk%P(1:nsc,    Ig, Jg) = blk%P(1:nsc,    Im, Jm)
    blk%P(np:nprim, Ig, Jg) = blk%P(np:nprim, Im, Jm)

    ! RSM (SSGLRR): np:nprim stays zero-gradient here, mirroring MOSE
    ! Ghost_Symmetry. The odd in-plane shear R12 ghost (-R12) is set in the
    ! viscous BC flux Lib_BC_Fluxes_Symmetry::BC_Symmetry_Visc (RSM_Symmetry),
    ! where the symmetry Reynolds-stress diffusive flux is also assembled —
    ! canonical MOSE placement (setting it here too would diverge from MOSE).

  end subroutine Ghost_Symmetry


  subroutine Ghost_ZG_Extrapolate ( Im, Jm, Fm, blk )
    implicit none
    integer,              intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type), intent(inout) :: blk
    ! Local
    integer :: g

    do g = 1, gc
      select case (Fm)
        case(1)
          blk%P(:, Im-g, Jm) = blk%P(:, Im, Jm)
        case(2)
          blk%P(:, Im+g, Jm) = blk%P(:, Im, Jm)
        case(3)
          blk%P(:, Im, Jm-g) = blk%P(:, Im, Jm)
        case(4)
          blk%P(:, Im, Jm+g) = blk%P(:, Im, Jm)
      end select
    end do

  end subroutine Ghost_ZG_Extrapolate


  subroutine Ghost_Chimera ( nb, blk, bc )
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, EOS, prim2cons, cons2prim
    implicit none
    integer,              intent(in)    :: nb
    type(Q2D_block_type), intent(inout) :: blk(nb)
    type(Q2D_bc_type),    intent(inout) :: bc
    ! Local
    integer  :: Bm, Im, Jm, Fm, Ig, Jg, Ig2, Jg2, Bs, Is, Js, c
    real(R8) :: rho, Rgas, temperature, consi(nprim), consg(nprim)

    Bm = bc%b
    Im = bc%i
    Jm = bc%j
    Fm = bc%f

    Ig  = Im - guide(Fm,1)
    Jg  = Jm - guide(Fm,2)
    Ig2 = Im - guide(Fm,1)*2
    Jg2 = Jm - guide(Fm,2)*2

    call co_rotot_Rtot ( blk(Bm)%P(1:nsc, Im, Jm), rho, Rgas )
    temperature = EOS( blk(Bm)%P(np, Im, Jm), rho, Rgas )

    consg = 0.d0
    do c = 1, bc%ni(1)
      Bs = bc%donorID(c, 1)
      Is = bc%donorID(c, 2)
      Js = bc%donorID(c, 3)
      consi(1:np) = prim2cons ( blk(Bs)%P(1:np, Is, Js) )
      consg = consg + consi*bc%volume_fraction(c)
    end do
    blk(Bm)%P(1:np, Ig, Jg) = cons2prim ( consg(1:np), temperature )
    bc%Pg(:,1) = blk(Bm)%P(:, Ig, Jg)

    consg = 0.d0
    do c = bc%ni(1)+1, sum(bc%ni)
      Bs = bc%donorID(c, 1)
      Is = bc%donorID(c, 2)
      Js = bc%donorID(c, 3)
      consi(1:np) = prim2cons ( blk(Bs)%P(1:np, Is, Js) )
      consg = consg + consi*bc%volume_fraction(c)
    end do
    blk(Bm)%P(1:np, Ig2, Jg2) = cons2prim ( consg(1:np), temperature )
    bc%Pg(:,2) = blk(Bm)%P(:, Ig2, Jg2)

  end subroutine Ghost_Chimera


  subroutine Ghost_Extrapolate ( Im, Jm, Fm, blk )
    implicit none
    integer,              intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type), intent(inout) :: blk
    ! Local
    integer :: g

    ! 2nd-order Taylor extrapolation: x4 = 3*x3 - 3*x2 + x1.
    do g = 1, gc
      select case (Fm)
        case(1)
          blk%P(:, Im-g, Jm) = 3d0*blk%P(:, Im-g+1, Jm) &
                              - 3d0*blk%P(:, Im-g+2, Jm) &
                              +     blk%P(:, Im-g+3, Jm)
        case(2)
          blk%P(:, Im+g, Jm) = 3d0*blk%P(:, Im+g-1, Jm) &
                              - 3d0*blk%P(:, Im+g-2, Jm) &
                              +     blk%P(:, Im+g-3, Jm)
        case(3)
          blk%P(:, Im, Jm-g) = 3d0*blk%P(:, Im, Jm-g+1) &
                              - 3d0*blk%P(:, Im, Jm-g+2) &
                              +     blk%P(:, Im, Jm-g+3)
        case(4)
          blk%P(:, Im, Jm+g) = 3d0*blk%P(:, Im, Jm+g-1) &
                              - 3d0*blk%P(:, Im, Jm+g-2) &
                              +     blk%P(:, Im, Jm+g-3)
      end select
    end do

  end subroutine Ghost_Extrapolate


  subroutine Fill_BC_Ghost_Connection ( Is, Js, Fs, d11s, d12s, d21s, d22s, blks, Pg )
    implicit none
    integer,              intent(in)    :: Is, Js, Fs, d11s, d12s, d21s, d22s
    type(Q2D_block_type), intent(in)    :: blks
    real(R8),             intent(inout) :: Pg(:,:)
    ! Local
    integer :: i1, j1, i2, j2

    select case (Fs)
      case(1:2)
        i1 = Is
        j1 = Js - d11s
        i2 = Is
        j2 = Js + d11s
      case(3:4)
        i1 = Is - d11s
        j1 = Js
        i2 = Is + d11s
        j2 = Js
    end select

    Pg(:,3) = blks%P(:, i1, j1)
    Pg(:,4) = blks%P(:, i2, j2)

  end subroutine Fill_BC_Ghost_Connection


  subroutine Fill_BC_Ghost_Chimera ( Im, Jm, Fm, blk, Pg )
    implicit none
    integer,              intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type), intent(inout) :: blk
    real(R8),             intent(inout) :: Pg(:,:)
    ! Local
    integer :: i1, j1, i2, j2
    integer :: dim2(2), Ig, Jg

    dim2 = blk%dim
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)

    select case (Fm)
      case(1:2)
        if (Jm == 1)        blk%P(:, Ig, Jm-1) = blk%P(:, Ig, Jm)
        if (Jm == dim2(2))  blk%P(:, Ig, Jm+1) = blk%P(:, Ig, Jm)
      case(3:4)
        if (Im == 1)        blk%P(:, Im-1, Jg) = blk%P(:, Im, Jg)
        if (Im == dim2(1))  blk%P(:, Im+1, Jg) = blk%P(:, Im, Jg)
    end select

    select case (Fm)
      case(1:2)
        i1 = Im ; j1 = Jm - 1
        i2 = Im ; j2 = Jm + 1
      case(3:4)
        i1 = Im - 1 ; j1 = Jm
        i2 = Im + 1 ; j2 = Jm
    end select

    Pg(:,3) = blk%P(:, i1, j1)
    Pg(:,4) = blk%P(:, i2, j2)

  end subroutine Fill_BC_Ghost_Chimera


  subroutine Ghost_Wall_Extrapolation ( domain )
    ! 2D-reduced (drop km/kg/nw). MPI: iterates only over local BC via
    ! domain%n_local_bc / domain%local_bc_idx (populated by build_local_bc_index).
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, ig, jg
    integer :: bm, im, jm, fm

    !$omp do private(ii, i, bm, im, jm, fm, ig, jg)
    ! BC array processing — iterate only over local BC entries
    do ii = 1, domain%n_local_bc
      i = domain%local_bc_idx(ii)
      bm = domain%bc(i)%b
      im = domain%bc(i)%i
      jm = domain%bc(i)%j
      fm = domain%bc(i)%f
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)

      if ( Is_Wall(domain%bc(i)%type) ) then
        ! Invert velocity; copy species + pressure + scalar tail (soot/passive/RANS) by
        ! zero-gradient, mirroring MOSE P(np:nprim). The RANS sub-band nt:nprim is then
        ! overwritten by RANS_Extrapolate_Wall in the wall BC routines (as in MOSE); the
        ! soot/passive sub-band np+1:nt-1 is set ONLY here — copying only P(np) (the
        ! pre-Sess-35 form) left it stale when nsoot/npass>0 (Q2D-only porting loss).
        domain%blk(bm)%P(nu:nv,    ig, jg) = -domain%blk(bm)%P(nu:nv,    im, jm)
        domain%blk(bm)%P(1:nsc,    ig, jg) =  domain%blk(bm)%P(1:nsc,    im, jm)
        domain%blk(bm)%P(np:nprim, ig, jg) =  domain%blk(bm)%P(np:nprim, im, jm)
      end if
    end do

  end subroutine Ghost_Wall_Extrapolation


  ! Wall BC type test (single definition; Q2D_Lib_Metrics imports from here).
  pure function Is_Wall(type) result(ans)
    implicit none
    integer, intent(in) :: type
    logical             :: ans

    select case (type)
      case (301, 302, 303, 304)
        ans = .true.
      case default
        ans = .false.
    end select
  end function Is_Wall

end module Q2D_Lib_Ghost
