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
    ! Locals
    integer :: b, d, nblocks

    ! Primitive-variable indices: the velocity has only two components
    ! V=[u,v] (indices nu, nv), consistently with the 2D formulation.
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

    ! The domain mirrors the block structure of IOfield.
    nblocks = size(IOfield%block)
    allocate(domain%blk(1:nblocks))
    domain%nb = nblocks

    ! i,j dimensions of each block (2D grid).
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

      ! Import of the mesh nodes (d=1..2: x,y coordinates).
      do d = 1, 2
        domain%blk(b)%node(0:IOfield%block(b)%Ni, 0:IOfield%block(b)%Nj)%c(d) &
          = IOfield%block(b)%mesh(d, 0:, 0:, 0)
      end do

      ! Import of the primitives (k index of the ORION structure fixed to 1).
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
    ! Locals
    integer :: ni, nj

    ni = nij(1)
    nj = nij(2)

    ! Metrics on the 2D grid: sur is the cell area and plays the role of the volume.
    allocate(blk%node(0:ni, 0:nj))
    allocate(blk%M   (1-gc:ni+gc, 1-gc:nj+gc))
    allocate(blk%dl  (1-gc:ni+gc, 1-gc:nj+gc))
    allocate(blk%sur (1-gc:ni+gc, 1-gc:nj+gc))
    allocate(blk%dir(1)%f(0:ni, nj))
    allocate(blk%dir(2)%f(ni, 0:nj))
    ! yn (wall distance) allocated unconditionally; Compute_Yn is executed
    ! only for `model>1`, so in laminar cases yn remains uninitialized but is
    ! never read (downstream consumers are guarded by `if (model==2)`).
    allocate(blk%yn(1-gc:ni+gc, 1-gc:nj+gc))

    ! Cross_Section for the area-variation source. source=0 is defense in
    ! depth: with the per-block gate (Mod_Metrics) a block without an area
    ! file never reaches Precompute_Area_Geometry, but the zero-initialization
    ! guarantees the absence of uninitialized reads even if a future consumer
    ! omitted the gate.
    if (obj_source%area_variation) then
      allocate(blk%Cross_Section(0:ni, 0:nj), source=0.0_R8)
    end if

    ! Primitives and residuals (rank-3 arrays with first dimension nprim).
    allocate(blk%P(nprim, 1-gc:ni+gc, 1-gc:nj+gc))
    allocate(blk%PO, blk%R, mold = blk%P)

    ! Local time step (cell centers, no ghosts).
    allocate(blk%dtlocal(1:ni, 1:nj))

    ! Velocity gradient and Spalart-Shur correction terms
    ! (conditional allocation).
    if (obj_rans%SpalartShur) then
      allocate(blk%vel_gradient(1-gc:ni+gc, 1-gc:nj+gc))
      allocate(blk%rc_term1(ni, nj))
      allocate(blk%rc_term2(ni, nj))
    end if

    ! Storage for the implicit residual smoothing.
    if (obj_irs%enabled) then
      allocate(blk%RS1, blk%RS2, mold = blk%R)
    end if

    ! Shock-sensor flag (rank-2 field on the 2D grid).
    ! beta is allocated unconditionally, with no gate on the sensor: Write_Diagnostic
    ! and Mod_GhostExchange::gather_diagnostic_to_root read `blk%beta(i,j)` even
    ! with the sensor disabled, so a conditional allocation would produce reads
    ! of an unallocated array (2D steady-state regressions: chang-nozzle and
    ! prandtl-meyer).
    allocate(blk%beta(1:ni, 1:nj))

    ! Precomputed geometric fields for the axisymmetric source.
    ! axis_r_ave/axis_inv_r are allocated `(1:ni, 1:nj)` WITHOUT ghosts: the dummy in
    ! Source_Axi_*_Blk is the explicit-shape `r_ave(1:n(1), 1:n(2))`, so with
    ! ghost bounds the actual-argument association would start from the lower bound
    ! (1-gc) and r_ave(1,1) would read the uninitialized ghost cell
    ! axis_r_ave(1-gc, 1-gc) → spurious values → diverging axisymmetric source.
    if (obj_source%axisymmetric) then
      allocate(blk%axis_r_ave(1:ni, 1:nj))
      allocate(blk%axis_inv_r(1:ni, 1:nj))
    end if
    ! area_geom_term_x/y allocated UNCONDITIONALLY and zero-initialized (default
    ! Definition B: with V_dot_G=0 the formulation reduces to the standard form).
    ! With a conditional allocation an unallocated allocatable would be passed
    ! as an explicit-shape dummy to Mod_Fluxes::Fluxes_blk_omp_do in NS cases
    ! without area_variation → SIGSEGV. Definition A (V_dot_G != 0) is activated
    ! by Precompute_Area_Geometry only with obj_source%area_variation=.true..
    ! Same universal zero-initialization pattern in Lib_Diffusive, Lib_SST,
    ! Lib_Wilcox2006 and Lib_Source_Areavariation. Memory cost ~64 KB/block,
    ! negligible: the universal zero-initialization avoids a conditional dispatch
    ! spread over 7 files with identical observable behavior
    ! (in IEEE arithmetic 0*x=0 is exact).
    allocate(blk%area_geom_term_x(1:ni, 1:nj), source=0.0_R8)
    allocate(blk%area_geom_term_y(1:ni, 1:nj), source=0.0_R8)
    ! area_variation_active is governed EXCLUSIVELY by Read_Block_Q2D_Area
    ! (.true. if and only if the block's area file exists and is readable):
    ! no override here. A block without a file keeps the flag .false., so
    ! Mod_Metrics skips Precompute_Area_Geometry and the universal
    ! zero-initialization (above) remains valid — controlled per-block
    ! degradation. The multigrid propagates the per-block flag
    ! (Coarse_Cross_Section).

  end subroutine Allocate_Block


  subroutine deallocate_remote_computation_data(domain)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block, mpi_is_root
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Locals
    integer :: b, d, i, c
    logical, allocatable :: needs_remote_P(:)

    ! In a serial build is_local_block(b)=.true. and mpi_is_root=.true., so
    ! this routine is effectively a no-op.
    if (.not. allocated(domain%bc)) return  ! BCs not populated yet.

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
