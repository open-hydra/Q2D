module Q2D_Lib_Chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Chemistry
  public :: Chemistry_Newstate
  public :: Chemistry_Equilibrium

contains

  subroutine Setup_Chemistry()
    use Q2D_Config_Types_m,       only: obj_chemistry
    use Q2D_Global_m,             only: nsc
    use FLINT_Lib_Chemistry_wdot, only: Assign_Mechanism
    use FLINT_CEA_setup,          only: CEA_initialize_global
    use oslo,                     only: Setup_ODEsolver
    implicit none

    ! Initialize chemistry model
    if (obj_chemistry%model=='equilibrium') then

      ! Master thread initializes; other threads will lazy-init
      ! inside CEA_solve on first call (threadprivate variables).
      call CEA_initialize_global()
      obj_chemistry%use_strang = .false.
      obj_chemistry%imodel = 2

    elseif (obj_chemistry%model=='finite-rate') then

      call Assign_Mechanism(obj_chemistry%mechanism_name)
      obj_chemistry%iopt = 0
      obj_chemistry%iopt(1) = obj_chemistry%max_ode_steps
      call Setup_ODEsolver(N=(nsc+1), solver=obj_chemistry%ode_name, &
                           RT=obj_chemistry%RT, AT=obj_chemistry%AT, &
                           iopt=obj_chemistry%iopt)
      obj_chemistry%use_strang = .true.
      obj_chemistry%imodel = 1

    else

      obj_chemistry%use_strang = .false.
      obj_chemistry%imodel = 0
      return

    endif

    ! exclude-blocks (no-chemistry block list) is read upstream with no length cap
    ! by Scan_NoChem (Q2D_Backend_INI), invoked from both INI read paths before
    ! Setup_Chemistry runs. Guarantee an allocated (possibly empty) list here so the
    ! consumers (Mod_Allocate_Data, Mod_Diagnostic) never see an unallocated array.
    if (.not. allocated(obj_chemistry%no_chem_list)) allocate(obj_chemistry%no_chem_list(0))

  end subroutine Setup_Chemistry


  subroutine Chemistry_Newstate ( P, dt, n )
    use Q2D_Global_m,            only: nprim, nsc, np, gc
    use FLINT_Lib_Chemistry_rhs, only: rhs_native
    use FLINT_Lib_Thermodynamic, only: f_Rtot
    use oslo,                    only: Run_ODESolver

    implicit none
    integer,  intent(in) :: n(2)
    real(R8), intent(in) :: dt(n(1),n(2))
    real(R8) :: P(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc)
    ! Local
    integer  :: lz
    real(R8) :: Z(nsc+1)
    real(R8) :: dtchem, t0, rho
    integer  :: i, j, s, err

    lz = nsc + 1 ! Z vector length: number of species + temperature.

    !$omp do collapse (2) schedule (dynamic) private(i,j,s,Z,dtchem,t0,err,rho)
    do j = 1, n(2); do i = 1, n(1)

      rho = sum( P(1:nsc,i,j) )

      Z(1:nsc) = P(1:nsc,i,j)
      Z(lz)    = P(np,i,j)/(rho*f_Rtot(Z(1:nsc)))
      dtchem   = dt(i,j)
      t0 = 0d0

      ! ODE integration according to selected method
      call Run_ODESolver ( lz, t0, dtchem, Z, rhs_native, err )

      ! New state from ODE integration in dtchem
      do s = 1, nsc
        if (Z(s)<=1d-20) Z(s) = 1d-20
        ! mass fractions
        P(s,i,j) = Z(s)
      end do
      ! pressure
      P(np,i,j) = rho*f_Rtot(Z(1:nsc))*Z(lz)

    enddo; enddo

  end subroutine Chemistry_Newstate


  subroutine Chemistry_Equilibrium ( P, n )
    use Q2D_Global_m,            only: nprim, nsc, np, gc
    use FLINT_CEA_solver,        only: CEA_solve
    use FLINT_Lib_Thermodynamic, only: f_Rtot
    implicit none
    integer,  intent(in) :: n(2)
    real(R8) :: P(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc)
    ! Local
    real(R8) :: rho, temp, teq, yeq(nsc), departure
    integer  :: i, j, s

    !$omp do collapse (2) schedule (dynamic) &
    !$omp private(i,j,s,rho,temp,teq,yeq,departure)
    do j = 1, n(2); do i = 1, n(1)

      rho  = sum( P(1:nsc,i,j) )
      temp = P(np,i,j)/(rho*f_Rtot(P(1:nsc,i,j)))

      call CEA_solve(temp, P(1:nsc,i,j), teq, yeq)

      ! Measure departure from equilibrium (max relative change in Y)
      departure = 0.0_R8
      do s = 1, nsc
        if (yeq(s) > 1d-20) then
          departure = max(departure, abs(yeq(s) - P(s,i,j)/rho) / yeq(s))
        end if
      end do
      departure = max(departure, abs(teq - temp)/teq)

      ! Skip correction if state is already at equilibrium to avoid noise accumulation
      if (departure < 1d-2) cycle

      ! Replace composition with equilibrium
      do s = 1, nsc
        if (yeq(s) <= 1d-20) yeq(s) = 1d-20
        P(s,i,j) = yeq(s)*rho
      end do

      ! Pressure from flow-solver EOS
      P(np,i,j) = rho*f_Rtot(P(1:nsc,i,j))*teq

    enddo; enddo

  end subroutine Chemistry_Equilibrium


end module Q2D_Lib_Chemistry
