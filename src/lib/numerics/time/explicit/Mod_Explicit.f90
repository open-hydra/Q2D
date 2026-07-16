module Q2D_Mod_Explicit
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Explicit_Step

contains

  subroutine Explicit_Step ( domain, External_Function )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_sim_param, obj_time_scheme, obj_multigrid, &
                                    obj_chemistry, obj_io, obj_soot, obj_rot
    use Q2D_Global_m,         only: nres, model, SOURCE_TERM
    use Q2D_Mod_dt,           only: Set_Global_dt, Compute_dt
    use Q2D_Lib_Ghost,        only: Fill_Ghost_Cell, Ghost_Wall_Extrapolation
    use Q2D_Mod_Fluxes,       only: Fluxes
    use Q2D_Mod_BC_Fluxes,    only: BC_Fluxes
    use Q2D_Mod_Newstate,     only: RK_Newstate
    use Q2D_Mod_Sources,      only: Compute_Sources
    use Q2D_Mod_Diagnostic,   only: Compute_Residual
    use Q2D_Lib_RANS,         only: RANS_Source_Terms
    use Q2D_Lib_Chemistry,    only: Chemistry_Newstate, Chemistry_Equilibrium
    use Q2D_Lib_RotatingFrame,only: RotatingFrame_Source_Terms
    use Q2D_Mod_Soot,         only: Soot_Source_Terms
    use Q2D_Mod_MPI,          only: is_local_block, mpi_reduce_sum_r8_array, &
                                    mpi_is_root, mpi_bcast_logical, mpi_bcast_integer
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    external :: External_Function
    ! Local variables
    logical  :: endsim, iosim, endmg
    integer  :: i_rk, i_strang, b, level
    real(R8) :: average(nres)

    level = obj_multigrid%MG_level
    obj_multigrid%change_MG = .false.
    domain(level) % iter = domain(level) % iter + 1
    obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
    obj_sim_param%iter_general   = obj_sim_param%iter_general + 1

    ! Time-step computation. The HYDRA_* flags are COUPLING flags: in standalone
    ! Q2D builds they stay at the default `.false.`; the high-level ATLAS
    ! orchestrator can set them to .true. for the HYDRA coupling.
    if (obj_sim_param%HYDRA_time_accurate) then
      call Set_Global_dt ( domain(level) )
      domain(level) % time = domain(level) % time + domain(level) % dtglobal
    else
      domain(level) % dtglobal = 1d5
      call Compute_dt ( domain(level), obj_time_scheme%cfl, obj_time_scheme%vnn, obj_time_scheme%rampa_cfl_iter )
      if ( obj_time_scheme%time_accurate ) then
        call Set_Global_dt ( domain(level) )
        domain(level) % time = domain(level) % time + domain(level) % dtglobal
      endif
    endif

    !$omp parallel
    strang: do i_strang = 1, obj_chemistry%N_strang ! Strang splitting: L = L_rk(dt/2) * L_chem(dt) * L_rk(dt/2)

      if (obj_chemistry%imodel == 2) then
        do b = 1, domain(level) % nb
          if (.not. is_local_block(b)) cycle
          call Chemistry_Equilibrium ( domain(level)%blk(b)%P, domain(level)%blk(b)%dim )
        enddo
      end if

      call Copy_State ( domain(level) )

      rk: do i_rk = 1, obj_time_scheme%n_RK

        call Fill_Ghost_Cell ( domain(level) )               ! ghost-cell filling
        call Fluxes ( domain(level) )                        ! convective and diffusive fluxes
        call BC_Fluxes ( domain(level) )                     ! boundary fluxes

        call External_Function ( domain(level) )             ! external function (e.g. source terms)

        ! The wall extrapolation into the ghost cells MUST precede every source
        ! term: the geometric source terms (viscous, from area variation) and the
        ! RANS ones read the wall ghost cells through the gradient stencil
        ! Prim(i±1,j±1). The gate includes SOURCE_TERM in addition to model==2
        ! (RANS ⇔ model==2) so as to also cover the laminar case with area-
        ! variation sources.
        if (model==2 .or. SOURCE_TERM) call Ghost_Wall_Extrapolation ( domain(level) )

        ! Dispatch of the geometric source terms (axisymmetric /
        ! area_variation). SOURCE_TERM gate: initialized by Initialize_Sources,
        ! invoked in Read_Inifile / Wrap_Setup; default .false. → no
        ! operation.
        if (SOURCE_TERM) call Compute_Sources ( domain(level) )

        ! Coriolis and centrifugal source terms in the rotating frame of
        ! reference (gate obj_rot%enabled — Setup_RotatingFrame applies a double
        ! check: axis strictly [0,0,1] and axisymmetric=.false.).
        if (obj_rot%enabled) call RotatingFrame_Source_Terms ( domain(level) )

        if (model==2) then
          call RANS_Source_Terms ( domain(level) )           ! RHS of the turbulence equations
        endif

        if (obj_soot%enabled) call Soot_Source_Terms ( domain(level) )

        call RK_Newstate ( domain(level), i_rk )             ! state update

      enddo rk

      ! Chemistry, only if i_strang=1 and N_strang=2
      if ( i_strang + obj_chemistry%N_strang == 3 ) then
        do b = 1, domain(level) % nb
          if (.not. is_local_block(b)) cycle
          if (domain(level)%blk(b)%no_chem) cycle
          call Chemistry_Newstate ( domain(level)%blk(b)%P, domain(level)%blk(b)%dtlocal, &
                                    domain(level)%blk(b)%dim )
        enddo
      endif

    enddo strang
    !$omp end parallel

    ! Saving of the residuals
    if ( (obj_sim_param%iter_from_call == 1) .or. (mod (domain(level) % iter, obj_io%res_diter) == 0) ) then
      obj_sim_param%residuotot = 0.d0
      do b = 1, domain(level) % nb
        if (.not. is_local_block(b)) cycle
        call Compute_Residual ( new=domain(level)%blk(b)%P, &
                                old=domain(level)%blk(b)%PO, &
                                dt=domain(level)%blk(b)%dtlocal, &
                                n=domain(level)%blk(b)%dim, &
                                average=average, &
                                total=obj_sim_param%residuotot )
      enddo
      call mpi_reduce_sum_r8_array(obj_sim_param%residuotot, nres)
      if (mpi_is_root) obj_sim_param%residuotot = sqrt ( obj_sim_param%residuotot ) ! L2 norm of the time derivative
    endif

    ! Determine the simulation-control flags on the root rank, then broadcast
    ! them to all ranks. endsim/endmg are local variables: at each call only one
    ! of them is assigned (depending on level), but both are read further below —
    ! they must be zeroed here, otherwise the unassigned one holds undefined
    ! stack values (without zeroing, the coarse multigrid levels can terminate
    ! the whole simulation after 1 iteration, as observed with TESTING-type
    ! gfortran builds).
    endsim = .false.
    endmg  = .false.
    if (mpi_is_root) then
      if (level == 1) then
        endsim = ( obj_sim_param%iter_from_call >= domain(1) % itermax ) &
            .or. ( obj_sim_param%residuotot(1) <= obj_sim_param%res_threshold ) &
            .or. ( domain(1) % time >= obj_sim_param%time_threshold )
      else
        endmg = ( obj_sim_param%iter_from_call >= domain(level) % itermax ) &
           .or. ( obj_sim_param%residuotot(1) <= obj_sim_param%res_threshold ) &
           .or. ( domain(level) % time >= obj_sim_param%time_threshold )
      endif

      iosim  = ( mod (domain(level) % iter, obj_io%sol_diter) == 0) &
          .or. ( domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime )

      if ( endsim ) then
        obj_sim_param%TODO = 3
      elseif ( iosim .or. endmg ) then
        obj_sim_param%TODO = 2
        if (endmg .and. obj_multigrid%MGL > 1) then
          obj_multigrid%change_MG = .true.
        endif
      else
        obj_sim_param%TODO = 1
      endif
    end if

    ! Broadcast of the simulation control from the root to all ranks
    call mpi_bcast_integer(obj_sim_param%TODO)
    call mpi_bcast_logical(obj_multigrid%change_MG)

  end subroutine Explicit_Step


  subroutine Copy_State ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: i, j, b

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(2)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)

        domain % blk(b) % PO (:,i,j) = domain % blk(b) % P (:,i,j)

      enddo; enddo
      !$omp end do
    enddo

  end subroutine Copy_State

end module Q2D_Mod_Explicit
