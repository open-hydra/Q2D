module Q2D_Lib_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Newstate_Conservative, Newstate_Primitive
  ! Check_And_Fix_State is private.

contains

  ! ===========================================================================
  !  CONSERVATIVE NEWSTATE
  ! ===========================================================================
  subroutine Newstate_Conservative ( domain, irk )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_time_scheme, obj_irs, obj_chemistry
    use Q2D_Global_m,         only: nprim, np
    use Q2D_Lib_IRS,          only: Residual_Smoothing
    use Q2D_Lib_RK,           only: RK_stage
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: irk
    ! Local
    integer  :: i, j, b
    integer  :: n_rk
    logical  :: irs_enabled
    real(R8) :: irs_beta, strangcoeff

    n_rk        = obj_time_scheme%n_RK
    irs_enabled = obj_irs%enabled
    irs_beta    = obj_irs%beta
    strangcoeff = obj_chemistry%strangcoeff

    ! Branch on IRS once outside all loops.
    if ( irs_enabled ) then

      ! PHASE 1: scale residuals only (no state update yet)
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(2) schedule(static) private(i,j)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Scale_Residual_Cons( domain%blk(b)%R(:,i,j),     &
                                    domain%blk(b)%dtlocal(i,j), &
                                    domain%blk(b)%sur(i,j), strangcoeff )
        enddo; enddo
        !$omp end do
      enddo

      ! PHASE 2: smooth residuals (serial across blocks by design)
      call Residual_Smoothing( domain, irs_beta )

      ! PHASE 3: update state after smoothing
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(2) schedule(static) private(i,j)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Update_State_Cons_IRS( domain%blk(b)%P(:,i,j),  &
                                      domain%blk(b)%PO(:,i,j), &
                                      domain%blk(b)%R(:,i,j),  &
                                      irk, n_rk, b, i, j )
        enddo; enddo
        !$omp end do
      enddo

    else

      ! No IRS: one pass — compute and update in same loop
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(2) schedule(static) private(i,j)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_Cons_NoIRS( domain%blk(b)%P(:,i,j),    &
                                   domain%blk(b)%PO(:,i,j),   &
                                   domain%blk(b)%R(:,i,j),    &
                                   domain%blk(b)%dtlocal(i,j),&
                                   domain%blk(b)%sur(i,j),    &
                                   irk, n_rk, b, i, j, strangcoeff )
        enddo; enddo
        !$omp end do
      enddo

    endif

  contains

    ! Scale residual in-place for conservative integration (IRS path).
    ! Q2D divides by the cell surface (the 2D analogue of the 3D volume).
    subroutine Scale_Residual_Cons( residual, dt, surface, strangcoeff )
      use Q2D_Global_m, only: nprim
      implicit none
      real(R8), intent(inout) :: residual(nprim)
      real(R8), intent(in)    :: dt, surface, strangcoeff
      residual = -residual / surface * dt * strangcoeff
    end subroutine Scale_Residual_Cons

    ! Apply RK stage + cons2prim after IRS smoothing (residual already scaled).
    subroutine Update_State_Cons_IRS( prim, primold, residual, irk, n_rk, b, i, j )
      use Q2D_Global_m,            only: nsc, np, nprim
      use FLINT_Lib_Thermodynamic, only: prim2cons, cons2prim, co_rotot_Rtot
      implicit none
      real(R8), intent(inout) :: prim(nprim), residual(nprim)
      real(R8), intent(in)    :: primold(nprim)
      integer,  intent(in)    :: irk, n_rk, b, i, j
      real(R8) :: rho, Rgas, temperature, consold(nprim), cons(nprim)

      consold(1:np) = prim2cons( primold(1:np) )
      cons   (1:np) = prim2cons( prim   (1:np) )
      if (nprim > np) then
        consold(np+1:nprim) = primold(np+1:nprim)
        cons   (np+1:nprim) = prim   (np+1:nprim)
      endif

      cons = RK_stage( irk, n_rk, cons, consold, residual )

      ! Reuse current prim pressure/density for temperature guess (cheap)
      call co_rotot_Rtot( prim(1:nsc), rho, Rgas )
      temperature = prim(np) / ( rho * Rgas )

      prim(1:np) = cons2prim( cons(1:np), temperature )
      if (nprim > np) prim(np+1:nprim) = cons(np+1:nprim)

      call Check_And_Fix_State( prim, b, i, j )

    end subroutine Update_State_Cons_IRS

    ! Full compute+update without IRS: prim => cons => RK stage => prim => check.
    subroutine Compute_Cons_NoIRS( prim, primold, residual, dt, surface, irk, n_rk, b, i, j, strangcoeff )
      use Q2D_Global_m,            only: nsc, np, nprim
      use FLINT_Lib_Thermodynamic, only: prim2cons, cons2prim, co_rotot_Rtot
      implicit none
      real(R8), intent(inout) :: prim(nprim), residual(nprim)
      real(R8), intent(in)    :: primold(nprim), dt, surface, strangcoeff
      integer,  intent(in)    :: irk, n_rk, b, i, j
      real(R8) :: rho, Rgas, temperature, consold(nprim), cons(nprim)

      residual = -residual / surface * dt * strangcoeff

      consold(1:np) = prim2cons( primold(1:np) )
      cons   (1:np) = prim2cons( prim   (1:np) )
      if (nprim > np) then
        consold(np+1:nprim) = primold(np+1:nprim)
        cons   (np+1:nprim) = prim   (np+1:nprim)
      endif

      cons = RK_stage( irk, n_rk, cons, consold, residual )

      call co_rotot_Rtot( prim(1:nsc), rho, Rgas )
      temperature = prim(np) / ( rho * Rgas )

      prim(1:np) = cons2prim( cons(1:np), temperature )
      if (nprim > np) prim(np+1:nprim) = cons(np+1:nprim)

      call Check_And_Fix_State( prim, b, i, j )

    end subroutine Compute_Cons_NoIRS

  end subroutine Newstate_Conservative


  ! ===========================================================================
  !  PRIMITIVE NEWSTATE
  ! ===========================================================================
  subroutine Newstate_Primitive ( domain, irk )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_time_scheme, obj_irs, obj_chemistry
    use Q2D_Global_m,         only: nprim
    use Q2D_Lib_IRS,          only: Residual_Smoothing
    use Q2D_Lib_RK,           only: RK_stage
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: irk
    ! Local
    integer  :: i, j, b
    integer  :: n_rk
    logical  :: irs_enabled
    real(R8) :: irs_beta, strangcoeff

    n_rk        = obj_time_scheme%n_RK
    irs_enabled = obj_irs%enabled
    irs_beta    = obj_irs%beta
    strangcoeff = obj_chemistry%strangcoeff

    ! Branch on IRS once outside all loops.
    if ( irs_enabled ) then

      ! PHASE 1: convert residual to primitive space + scale (no state update yet)
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(2) schedule(static) private(i,j)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_Prim_Residual( domain%blk(b)%P(:,i,j),    &
                                      domain%blk(b)%R(:,i,j),    &
                                      domain%blk(b)%dtlocal(i,j),&
                                      domain%blk(b)%sur(i,j), strangcoeff )
        enddo; enddo
        !$omp end do
      enddo

      ! PHASE 2: smooth primitive residuals
      call Residual_Smoothing( domain, irs_beta )

      ! PHASE 3: apply RK stage after smoothing
      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(2) schedule(static) private(i,j)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          domain%blk(b)%P(:,i,j) = RK_stage( irk, n_rk,                  &
                                             domain%blk(b)%P(:,i,j),     &
                                             domain%blk(b)%PO(:,i,j),    &
                                             domain%blk(b)%R(:,i,j) )
          call Check_And_Fix_State( domain%blk(b)%P(:,i,j), b, i, j )
        enddo; enddo
        !$omp end do
      enddo

    else

      do b = 1, domain%nb
        if (.not. is_local_block(b)) cycle
        !$omp do collapse(2) schedule(static) private(i,j)
        do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          call Compute_And_Update_Prim( domain%blk(b)%P(:,i,j),    &
                                        domain%blk(b)%PO(:,i,j),   &
                                        domain%blk(b)%R(:,i,j),    &
                                        domain%blk(b)%dtlocal(i,j),&
                                        domain%blk(b)%sur(i,j),    &
                                        irk, n_rk, b, i, j, strangcoeff )
        enddo; enddo
        !$omp end do
      enddo

    endif

  contains

    ! Converts conservative residual to primitive residual in-place
    subroutine Compute_Prim_Residual( prim, residual, dt, surface, strangcoeff )
      use Q2D_Global_m,            only: nsc, nu, nv, np, nprim
      use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, f_tabT_expr, h_tab, cp_tab, Ri_tab
      implicit none
      real(R8), intent(in)    :: prim(nprim), dt, surface, strangcoeff
      real(R8), intent(inout) :: residual(nprim)
      integer  :: T_i, Tint(2), s
      real(R8) :: Tdiff, inv_surface, inv_rhosurface, rho, Rgas, temperature
      real(R8) :: dwrotot, dwE, dwroe, dwTemp, cv, eiroi
      real(R8) :: velocity(2)

      call co_rotot_Rtot( prim(1:nsc), rho, Rgas )
      temperature = prim(np) / ( rho * Rgas )

      T_i     = int( temperature )
      Tdiff   = temperature - T_i
      Tint(1) = T_i
      Tint(2) = T_i + 1

      inv_surface    = 1d0 / surface
      inv_rhosurface = inv_surface / rho
      velocity       = prim(nu:nv)

      ! Species densities
      residual(1:nsc) = -residual(1:nsc) * inv_surface
      dwrotot = sum( residual(1:nsc) )

      ! Velocity (2 components, no nw)
      residual(nu:nv) = -residual(nu:nv) * inv_rhosurface - velocity / rho * dwrotot

      ! Pressure via energy
      dwE   = -residual(np) * inv_surface
      dwroe = dwE - rho * dot_product( velocity, residual(nu:nv) ) &
                  - 0.5d0 * dwrotot * sum( velocity**2 )
      eiroi = 0d0; cv = 0d0
      do s = 1, nsc
        eiroi = eiroi + ( f_tabT_expr(s,h_tab,Tint,Tdiff) - Ri_tab(s)*temperature ) * residual(s)
        cv    = cv    + prim(s) * ( f_tabT_expr(s,cp_tab,Tint,Tdiff) - Ri_tab(s) )
      enddo
      dwTemp = ( dwroe - eiroi ) / cv

      residual(np) = sum( Ri_tab(1:nsc) * ( temperature*residual(1:nsc) + prim(1:nsc)*dwTemp ) )

      ! Turbulence + soot + passive scalars (everything beyond np)
      if (nprim > np) residual(np+1:nprim) = -residual(np+1:nprim) * inv_surface

      ! Time scaling (Strang operator splitting)
      residual = residual * dt * strangcoeff

    end subroutine Compute_Prim_Residual

    ! Full primitive residual + RK update in one shot (no IRS)
    subroutine Compute_And_Update_Prim( prim, primold, residual, dt, surface, irk, n_rk, b, i, j, strangcoeff )
      use Q2D_Global_m, only: nprim
      implicit none
      real(R8), intent(inout) :: prim(nprim), residual(nprim)
      real(R8), intent(in)    :: primold(nprim), dt, surface, strangcoeff
      integer,  intent(in)    :: irk, n_rk, b, i, j

      call Compute_Prim_Residual( prim, residual, dt, surface, strangcoeff )
      prim = RK_Stage( irk, n_rk, prim, primold, residual )
      call Check_And_Fix_State( prim, b, i, j )

    end subroutine Compute_And_Update_Prim

  end subroutine Newstate_Primitive


  ! ===========================================================================
  !  STATE SANITY CHECK AND FLOOR ENFORCEMENT
  ! ===========================================================================
  subroutine Check_And_Fix_State ( prim, b, i, j )
    use Q2D_Global_m, only: nsc, nt, np, nc, nsoot, nprim, model
    use Q2D_Lib_RANS, only: RANS_Enforce_Realizability

    implicit none
    real(R8), dimension(:), intent(inout) :: prim
    integer, intent(in) :: b, i, j
    integer :: s

    ! Bail out on NaN or negative pressure (any(isnan) + a separate pressure
    ! check is cheaper than the product trick).
    if ( any( isnan(prim(1:nprim)) ) .or. prim(np) < 0d0 ) then
      write(*,'(A,3I4)') "Integration failed at b, i, j:", b, i, j
      write(*,*) prim(1:nprim)
      stop "NaN or p<0 detected"
    endif

    ! Floor species densities
    do s = 1, nsc
      if ( prim(s) < 1d-20 ) prim(s) = 1d-20
    enddo

    ! Realizability for turbulence (gate `model==2`; formally
    ! equivalent to the Q2D legacy `if (RANS)` — RANS .iff. model==2).
    if (model==2) call RANS_Enforce_Realizability ( prim(nt:nprim) )

    ! Floor soot variables.
    if (nsoot == 2) then
      if ( prim(nc)   < 1d-20 ) prim(nc)   = 1d-20
      if ( prim(nc+1) < 1d+10 ) prim(nc+1) = 1d+10
    elseif (nsoot == 1) then
      if ( prim(nc)   < 1d+10 ) prim(nc)   = 1d+10
    endif

  end subroutine Check_And_Fix_State

end module Q2D_Lib_Newstate
