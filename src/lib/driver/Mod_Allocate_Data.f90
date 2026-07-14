module Q2D_Mod_Allocate_Data
  use, intrinsic :: iso_fortran_env, only: iostat_end

  implicit none
  private

  public :: Setup_Data_Structure
  public :: Allocate_Block
  public :: deallocate_remote_computation_data

contains

  subroutine Setup_Data_Structure(domain, IOfield)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_io, obj_chemistry
    use Q2D_Global_m,         only: nu, nv, np, nt, nc, nps, &
                                    nsc, nsoot, npass, nrans, &
                                    nprim, nres
    use Lib_ORION_data,       only: ORION_data
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    type(ORION_data),      intent(inout) :: IOfield
    ! Local
    integer :: b, d, nblocks

    ! Primitive variable indices (no nw — V=[u,v] on 2 components).
    nu = nsc + 1
    nv = nu + 1
    np = nv + 1
    nprim = np + nsoot + npass + nrans

    if (nsoot > 0) nc  = np + 1
    if (npass > 0) nps = np + nsoot + 1
    if (nrans > 0) then
      nt = np + nsoot + npass + 1
    else
      nt = nprim
    end if

    nres = nprim - nsc + 1

    ! Domain is an alias of IOfield.
    nblocks = size(IOfield%block)
    allocate(domain%blk(1:nblocks))
    domain%nb = nblocks

    ! i,j dimensions of each block (2 dimensions; no Nk).
    do b = 1, nblocks
      domain%blk(b)%dim(1) = IOfield%block(b)%Ni
      domain%blk(b)%dim(2) = IOfield%block(b)%Nj
    end do

    if (size(IOfield%block(1)%vars, 1) < nprim) then
      write(*,'(A)')         '[ERROR] Number of variables in IOfield does not match Q2D expectation.'
      write(*,'(A,I0,A,I0)') '        Expected: ', nprim, ', Found: ', size(IOfield%block(1)%vars, 1)
      stop
    end if

    do b = 1, nblocks
      call Allocate_Block(domain%blk(b), domain%blk(b)%dim)

      ! Import mesh nodes (d=1..2, no z component).
      do d = 1, 2
        domain%blk(b)%node(0:IOfield%block(b)%Ni, 0:IOfield%block(b)%Nj)%c(d) &
          = IOfield%block(b)%mesh(d, 0:, 0:, 0)
      end do

      ! Import primitives (k=1 fixed from ORION).
      do d = 1, nprim
        domain%blk(b)%P(d, 1:IOfield%block(b)%Ni, 1:IOfield%block(b)%Nj) &
          = IOfield%block(b)%vars(d, 1:, 1:, 1)
      end do
      domain%time = obj_io%IOtime

      if (allocated(obj_chemistry%no_chem_list)) then
        if (any(obj_chemistry%no_chem_list == b)) then
          domain%blk(b)%no_chem = .true.
        end if
      end if
    end do

  end subroutine Setup_Data_Structure


  subroutine Allocate_Block(blk, nij)
    use iso_fortran_env,      only: R8 => real64
    use Q2D_Advanced_Types_m, only: Q2D_block_type
    use Q2D_Config_Types_m,   only: obj_irs, obj_rans, obj_source, obj_multigrid, obj_shock_detector
    use Q2D_Global_m,         only: gc, nprim
    implicit none
    integer,              intent(in)    :: nij(2)
    type(Q2D_block_type), intent(inout) :: blk
    ! Local
    integer :: ni, nj

    ni = nij(1)
    nj = nij(2)

    ! Metrics (2D shape without k; sur replaces the 3D vol).
    allocate(blk%node(0:ni, 0:nj))
    allocate(blk%M   (1-gc:ni+gc, 1-gc:nj+gc))
    allocate(blk%dl  (1-gc:ni+gc, 1-gc:nj+gc))
    allocate(blk%sur (1-gc:ni+gc, 1-gc:nj+gc))
    allocate(blk%dir(1)%f(0:ni, nj))
    allocate(blk%dir(2)%f(ni, 0:nj))
    ! yn always allocated:
    ! the legacy Q2D `if (VISC) allocate(yn)` was defensive — in laminar yn was
    ! valid but unused (downstream gated `if (model==2)`). Allocate unconditionally
    ! + run Compute_Yn only for `model>1` (yn=garbage in laminar but discarded).
    allocate(blk%yn(1-gc:ni+gc, 1-gc:nj+gc))

    ! Cross_Section for area-variation. source=0 is defense-in-depth: with the
    ! per-block gate (Mod_Metrics) a no-file block never reaches Precompute_Area_
    ! Geometry, but zero-init guarantees no uninitialised read even if a future
    ! consumer drops the guard.
    if (obj_source%area_variation) then
      allocate(blk%Cross_Section(0:ni, 0:nj), source=0.0_R8)
    end if

    ! Primitive and residuals (rank-3 with first dim = nprim).
    allocate(blk%P(nprim, 1-gc:ni+gc, 1-gc:nj+gc))
    allocate(blk%PO, blk%R, mold = blk%P)

    ! Local time-step (cell-center, no ghost).
    allocate(blk%dtlocal(1:ni, 1:nj))

    ! Velocity gradient + Spalart-Shur corrections (conditional).
    if (obj_rans%SpalartShur) then
      allocate(blk%vel_gradient(1-gc:ni+gc, 1-gc:nj+gc))
      allocate(blk%rc_term1(ni, nj))
      allocate(blk%rc_term2(ni, nj))
    end if

    ! Implicit residual smoothing storage.
    if (obj_irs%enabled) then
      allocate(blk%RS1, blk%RS2, mold = blk%R)
    end if

    ! Shock detector flag (rank-3 reduced to rank-2 for 2D).
    ! beta is allocated unconditionally (no SD gate). Write_Diagnostic and
    ! Mod_GhostExchange::gather_diagnostic_to_root read `blk%beta(i,j)`
    ! regardless of SD: the `if (SD)` gate was a divergence that broke the 2D
    ! steady-state regression tests (chang-nozzle + prandtl-meyer).
    allocate(blk%beta(1:ni, 1:nj))

    ! Q2D-extra precomputed axisymmetric/area_variation fields.
    ! axis_r_ave/axis_inv_r allocated `(1:ni, 1:nj)` no-ghost mirror the Q2D legacy
    ! Precompute_Axisymmetric_Geometry alloc — the dummy in Source_Axi_*_Blk is
    ! `r_ave(1:n(1), 1:n(2))` no-ghost. When allocated with ghost bounds the
    ! assumed-shape slice passed from the lower-bound (1-gc) made dummy r_ave(1,1)
    ! read from the uninitialized ghost cell axis_r_ave(1-gc, 1-gc) → garbage →
    ! axisym source diverges.
    if (obj_source%axisymmetric) then
      allocate(blk%axis_r_ave(1:ni, 1:nj))
      allocate(blk%axis_inv_r(1:ni, 1:nj))
    end if
    ! area_geom_term_x/y allocated UNCONDITIONALLY with zero-fill (Definition B
    ! default; with V_dot_G=0 it reduces to the standard form). When conditional, an unallocated
    ! allocatable was passed as an explicit-shape dummy in
    ! Mod_Fluxes::Fluxes_blk_omp_do for non-area_variation NS tests → SIGSEGV
    ! (latent bug masked by Euler-only ctest). Definition A (V_dot_G != 0) is
    ! activated by Precompute_Area_Geometry only when obj_source%area_variation=
    ! .true.. Design intent documented in Lib_Diffusive + Lib_SST + Lib_Wilcox2006
    ! + Lib_Source_Areavariation ("universal zero-fill" pattern). Memory cost
    ! ~64 KB/block negligible: universal zero-fill = 5 LOC vs conditional dispatch
    ! refactor = 50 LOC over 7 files, same observable behavior (IEEE 0*x=0 exact).
    allocate(blk%area_geom_term_x(1:ni, 1:nj), source=0.0_R8)
    allocate(blk%area_geom_term_y(1:ni, 1:nj), source=0.0_R8)
    ! area_variation_active stays driven SOLELY by Read_Block_Q2D_Area (.true. iff
    ! the block's area file exists+parses) — no force-override here. A block without
    ! its file keeps the flag .false. so Mod_Metrics skips Precompute_Area_Geometry
    ! and the universal zero-fill (above) is preserved. Mirrors the legacy per-block
    ! graceful degradation; MG propagates the per-block flag (Coarse_Cross_Section).

  end subroutine Allocate_Block


  subroutine deallocate_remote_computation_data(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block, mpi_is_root
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, d, i, c
    logical, allocatable :: needs_remote_P(:)

    ! In a serial build is_local_block(b)=.true. and mpi_is_root=.true., so this
    ! routine is essentially a no-op.
    if (.not. allocated(domain%bc)) return  ! BC not yet populated.

    allocate(needs_remote_P(domain%nb))
    needs_remote_P = .false.

    do i = 1, domain%nbound
      select case (domain%bc(i)%type)
        case (102) ! chimera
          if (allocated(domain%bc(i)%donorID)) then
            do c = 1, size(domain%bc(i)%donorID, 1)
              b = domain%bc(i)%donorID(c, 1)
              if (.not. is_local_block(b)) needs_remote_P(b) = .true.
            end do
          end if
        case (501) ! manifold
          b = domain%bc(i)%bs
          if (b > 0 .and. .not. is_local_block(b)) needs_remote_P(b) = .true.
      end select
    end do

    do b = 1, domain%nb
      if (is_local_block(b)) cycle

      if (allocated(domain%blk(b)%PO))               deallocate(domain%blk(b)%PO)
      if (allocated(domain%blk(b)%R))                deallocate(domain%blk(b)%R)
      if (allocated(domain%blk(b)%RS1))              deallocate(domain%blk(b)%RS1)
      if (allocated(domain%blk(b)%RS2))              deallocate(domain%blk(b)%RS2)
      if (allocated(domain%blk(b)%dtlocal))          deallocate(domain%blk(b)%dtlocal)
      if (allocated(domain%blk(b)%beta))             deallocate(domain%blk(b)%beta)
      if (allocated(domain%blk(b)%vel_gradient))     deallocate(domain%blk(b)%vel_gradient)
      if (allocated(domain%blk(b)%rc_term1))         deallocate(domain%blk(b)%rc_term1)
      if (allocated(domain%blk(b)%rc_term2))         deallocate(domain%blk(b)%rc_term2)
      if (allocated(domain%blk(b)%axis_r_ave))       deallocate(domain%blk(b)%axis_r_ave)
      if (allocated(domain%blk(b)%axis_inv_r))       deallocate(domain%blk(b)%axis_inv_r)
      if (allocated(domain%blk(b)%area_geom_term_x)) deallocate(domain%blk(b)%area_geom_term_x)
      if (allocated(domain%blk(b)%area_geom_term_y)) deallocate(domain%blk(b)%area_geom_term_y)

      if (.not. needs_remote_P(b)) then
        if (allocated(domain%blk(b)%P)) deallocate(domain%blk(b)%P)
      end if

      if (.not. mpi_is_root) then
        if (allocated(domain%blk(b)%node))          deallocate(domain%blk(b)%node)
        if (allocated(domain%blk(b)%M))             deallocate(domain%blk(b)%M)
        if (allocated(domain%blk(b)%dl))            deallocate(domain%blk(b)%dl)
        if (allocated(domain%blk(b)%sur))           deallocate(domain%blk(b)%sur)
        if (allocated(domain%blk(b)%yn))            deallocate(domain%blk(b)%yn)
        if (allocated(domain%blk(b)%Cross_Section)) deallocate(domain%blk(b)%Cross_Section)
        if (.not. needs_remote_P(b)) then
          do d = 1, 2
            if (allocated(domain%blk(b)%dir(d)%f)) deallocate(domain%blk(b)%dir(d)%f)
          end do
        end if
      end if
    end do

    deallocate(needs_remote_P)

  end subroutine deallocate_remote_computation_data

end module Q2D_Mod_Allocate_Data
