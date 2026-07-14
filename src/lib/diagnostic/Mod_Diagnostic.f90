module Q2D_Mod_Diagnostic
  use iso_fortran_env,      only: R8 => real64
  use Q2D_Advanced_Types_m, only: Q2D_simulation_type, Q2D_domain_type
  use Q2D_Parameters_m,     only: llen, clen
  use Lib_ORION_data,       only: ORION_data

  implicit none
  private

  public :: Print_Header
  public :: Check_Tables
  public :: Check_Input
  public :: Check_IC
  public :: Check_BC
  public :: Print_Warnings
  public :: Stop_Simulation
  public :: Print_Shell_Info
  public :: Write_Diagnostic
  public :: Compute_Residual

  character(len=llen), private :: Dvarnames = '"rho" "rhou" "rhov" "rhoe" "dt" "beta"'

contains

  subroutine Print_Header
    ! Q2D ASCII art header.
    ! Exposed in Mod_Diagnostic for the driver/diagnostic separation.
    write(*,*)
    write(*,'(A44)') "    ,----..          ,----,      ,---,     "
    write(*,'(A44)') "   /   /   \       .'   .' \   .'  .' `\   "
    write(*,'(A44)') "  /   .     :    ,----,'    |,---.'     \  "
    write(*,'(A44)') " .   /   ;.  \   |    :  .  ;|   |  .`\  | "
    write(*,'(A44)') ".   ;   /  ` ;   ;    |.'  / :   : |  '  | "
    write(*,'(A44)') ";   |  ; \ ; |   `----'/  ;  |   ' '  ;  : "
    write(*,'(A44)') "|   :  | ; | '     /  ;  /   '   | ;  .  | "
    write(*,'(A44)') ".   |  ' ' ' :    ;  /  /-,  |   | :  |  ' "
    write(*,'(A44)') "'   ;  \; /  |   /  /  /.`|  '   : | /  ;  "
    write(*,'(A44)') " \   \  ',  . \./__;      :  |   | '` ,/   "
    write(*,'(A44)') "  ;   :      ; |   :    .'   ;   :  .'     "
    write(*,'(A44)') "   \   \ .'`--';   | .'      |   ,.'       "
    write(*,'(A44)') "    `---`      `---'         '---'         "
    write(*,*)
  end subroutine Print_Header


  subroutine Check_Tables
    ! Verifica obj_thermo/transport/chemistry%error_message after Load_*; stop on error.
    use Q2D_Config_Types_m, only: obj_thermo, obj_transport, obj_chemistry
    implicit none
    logical :: has_error

    has_error = .false.

    write(*,'(A)') ' ========================================================================================='
    write(*,'(A)') ' Loading'
    write(*,'(A)') ' ========================================================================================='

    ! Thermodynamics
    if (index(obj_thermo%error_message,'ERROR') > 0) then
      write(*,'(A,T35,A)') '   Thermodynamics', 'FAIL'
      write(*,'(4X,A)') trim(obj_thermo%error_message)
      has_error = .true.
    else
      write(*,'(A,T35,A)') '   Thermodynamics', 'OK'
    endif

    ! Transport
    if (index(obj_transport%error_message,'ERROR') > 0) then
      write(*,'(A,T35,A)') '   Transport', 'FAIL'
      write(*,'(4X,A)') trim(obj_transport%error_message)
      has_error = .true.
    else
      write(*,'(A,T35,A)') '   Transport', 'OK'
    endif

    ! Chemistry
    if (index(obj_chemistry%error_message,'ERROR') > 0) then
      write(*,'(A,T35,A)') '   Chemistry', 'FAIL'
      write(*,'(4X,A)') trim(obj_chemistry%error_message)
      has_error = .true.
    else
      write(*,'(A,T35,A)') '   Chemistry', 'OK'
    endif

    if (has_error) stop
  end subroutine Check_Tables


  subroutine Check_Input
    ! Validate_Registry is Q2D_Registry::Validate_Registry.
    use Q2D_Registry_m,   only: Validate_Registry
    use Q2D_Parameters_m, only: hlen
    implicit none
    character(len=hlen) :: out

    out = Validate_Registry()
    if (index(out,'ERROR') > 0) then
      write(*,'(A,T35,A)') '   Input file', 'FAIL'
      write(*,'(4X,A)') trim(out)
      stop
    else
      write(*,'(A,T35,A)') '   Input file', 'OK'
    endif
  end subroutine Check_Input


  subroutine Check_IC
    use Q2D_Config_Types_m, only: obj_io
    implicit none

    if (index(obj_io%error_message,'ERROR') > 0) then
      write(*,'(A,T35,A)') '   Initial conditions', 'FAIL'
      write(*,'(4X,A)') trim(obj_io%error_message)
      stop
    else
      write(*,'(A,T35,A)') '   Initial conditions', 'OK'
    endif
  end subroutine Check_IC


  subroutine Check_BC
    use Q2D_Config_Types_m, only: obj_io_bc
    implicit none

    if (index(obj_io_bc%error_message,'ERROR') > 0) then
      write(*,'(A,T35,A)') '   Boundary conditions', 'FAIL'
      write(*,'(4X,A)') trim(obj_io_bc%error_message)
      stop
    else
      write(*,'(A,T35,A)') '   Boundary conditions', 'OK'
    endif
    write(*,'(A)') ' ========================================================================================='
  end subroutine Check_BC


  subroutine Print_Warnings
    ! Stampa obj_*%warning_message of all singletons.
    use Q2D_Config_Types_m, only: obj_io, obj_thermo, obj_transport, obj_rans, obj_rot, &
                                  obj_chemistry, obj_soot, obj_time_scheme, obj_riemann, &
                                  obj_space_scheme, obj_irs
    implicit none

    ! IO warnings
    if (index(obj_io%warning_message,'WARNING') > 0) write(*,'(A)') obj_io%warning_message
    ! Physics warnings
    if (index(obj_thermo%warning_message,'WARNING') > 0)    write(*,'(A)') obj_thermo%warning_message
    if (index(obj_transport%warning_message,'WARNING') > 0) write(*,'(A)') obj_transport%warning_message
    if (index(obj_rans%warning_message,'WARNING') > 0)      write(*,'(A)') obj_rans%warning_message
    if (index(obj_rot%warning_message,'WARNING') > 0)       write(*,'(A)') obj_rot%warning_message
    if (index(obj_chemistry%warning_message,'WARNING') > 0) write(*,'(A)') obj_chemistry%warning_message
    if (index(obj_soot%warning_message,'WARNING') > 0)      write(*,'(A)') obj_soot%warning_message
    ! Numerical scheme warnings
    if (index(obj_time_scheme%warning_message,'WARNING') > 0)  write(*,'(A)') obj_time_scheme%warning_message
    if (index(obj_riemann%warning_message,'WARNING') > 0)      write(*,'(A)') obj_riemann%warning_message
    if (index(obj_space_scheme%warning_message,'WARNING') > 0) write(*,'(A)') obj_space_scheme%warning_message
    if (index(obj_irs%warning_message,'WARNING') > 0)          write(*,'(A)') obj_irs%warning_message
  end subroutine Print_Warnings


  subroutine Stop_Simulation
    ! Aggrega obj_*%error_message di rans/rot/soot/time_scheme/riemann/space_scheme/irs
    ! e stop su error.
    use Q2D_Config_Types_m, only: obj_rans, obj_rot, obj_soot, obj_time_scheme, &
                                  obj_riemann, obj_space_scheme, obj_irs
    implicit none
    logical :: has_error

    has_error = .false.

    ! Physics errors
    if (index(obj_rans%error_message,'ERROR') > 0) then
      write(*,'(A)') obj_rans%error_message;    has_error = .true.
    endif
    if (index(obj_rot%error_message,'ERROR') > 0) then
      write(*,'(A)') obj_rot%error_message;     has_error = .true.
    endif
    if (index(obj_soot%error_message,'ERROR') > 0) then
      write(*,'(A)') obj_soot%error_message;    has_error = .true.
    endif
    ! Numerical scheme errors
    if (index(obj_time_scheme%error_message,'ERROR') > 0) then
      write(*,'(A)') obj_time_scheme%error_message;  has_error = .true.
    endif
    if (index(obj_riemann%error_message,'ERROR') > 0) then
      write(*,'(A)') obj_riemann%error_message;      has_error = .true.
    endif
    if (index(obj_space_scheme%error_message,'ERROR') > 0) then
      write(*,'(A)') obj_space_scheme%error_message; has_error = .true.
    endif
    if (index(obj_irs%error_message,'ERROR') > 0) then
      write(*,'(A)') obj_irs%error_message;          has_error = .true.
    endif

    if (has_error) stop
  end subroutine Stop_Simulation


  subroutine Print_Shell_Info(simulation)
    ! 2D-reduced: `dim(1:3)` -> `dim(1:2)`. Domain topology, BC summary,
    ! IO, gas model, chemistry, soot, rotating frame, numerical scheme.
    use Q2D_Global_m,       only: model, nsc, npass, nsoot
    use Q2D_Config_Types_m, only: obj_io, obj_thermo, obj_chemistry, obj_soot, obj_rans, &
                                  obj_rot, obj_sim_param, obj_time_scheme, obj_riemann, &
                                  obj_space_scheme, obj_irs, obj_multigrid, obj_source
    use Q2D_IO_BC,          only: Print_BC_Summary
    use Q2D_IO_Probes,      only: nprobes
    use IR_Precision,       only: str
    implicit none
    type(Q2D_simulation_type), intent(in) :: simulation
    ! Local
    character(llen) :: eosword
    integer :: b, total_cells, m

    ! ----- Domain topology -----
    total_cells = 0
    do b = 1, simulation%domain(1)%nb
      total_cells = total_cells + product(simulation%domain(1)%blk(b)%dim(1:2))
    end do

    write(*,*)
    write(*,'(A)') ' ========================================================================================='
    write(*,'(A)') ' Set-up'
    write(*,'(A)') ' ========================================================================================='
    write(*,'(A)') ' Domain'
    write(*,'(A,T35,I0)') '   Blocks',         simulation%domain(1)%nb
    write(*,'(A,T35,I0)') '   Cells',          total_cells
    write(*,'(A,T35,I0)') '   Boundary faces', simulation%domain(1)%nbound
    if ( obj_multigrid%MGL > 1 ) &
      write(*,'(A,T35,I0)') '   Multigrid levels', obj_multigrid%MGL

    ! ----- Boundary conditions -----
    call Print_BC_Summary()

    ! ----- IO -----
    write(*,*)
    write(*,'(A)') ' Input/Output'
    write(*,'(A,T35,A)') '   Initial conditions file', trim(obj_io%nameinit)
    write(*,'(A,T35,A)') '   Solution format',         trim(obj_io%sol_format)
    write(*,'(A,T35,A)') '   Probes number',           str(.true., nprobes)
    write(*,'(A)') ' ========================================================================================='

    ! ----- Physical model -----
    eosword = 'Ideal'

    write(*,*)
    write(*,'(A)') ' ========================================================================================='
    write(*,'(A)') ' Physical model'
    write(*,'(A)') ' ========================================================================================='
    write(*,'(A)') ' Gas model'
    write(*,'(A,T35,A)') '   Equations', trim(obj_sim_param%description)
    if (model==1) then
      write(*,'(A,T35,A)') '   Viscous model', 'Laminar'
    elseif (model==2) then
      write(*,'(A,T35,A)') '   Viscous model', 'Turbulent, '//trim(obj_rans%model)
    end if
    write(*,'(A,T35,A)')  '   Equation of state', trim(eosword)
    write(*,'(A,T35,A)')  '   Thermodynamics',    trim(obj_thermo%description)
    write(*,'(A,T35,I0)') '   Species',           nsc

    if (obj_chemistry%model /= 'frozen') then
      write(*,*)
      write(*,'(A)') ' Chemistry'
      if (obj_chemistry%model == 'finite-rate') then
        write(*,'(A)')       '   Finite-rate '
        write(*,'(A,T35,A)') '   - Mechanism', trim(obj_chemistry%mechanism_name)
        write(*,'(A,T35,A)') '   - Splitting', 'Strang operator'
        write(*,'(A,T35,A)') '   - ODE solver', trim(obj_chemistry%ode_name)
      elseif (obj_chemistry%model == 'equilibrium') then
        write(*,'(A)')       '   CEA equilibrium '
      endif
      if ( size(obj_chemistry%no_chem_list) == 0 ) then
        write(*,'(A,T35,A)') '   Active blocks', 'all'
      else
        write(*,'(A,T35,*(I0,:," "))') '   Inactive blocks', obj_chemistry%no_chem_list
      endif
    endif

    if (obj_soot%enabled) then
      write(*,*)
      write(*,'(A)') ' Soot'
      if (trim(obj_soot%model) == 'LL91') then
        write(*,'(A,T35,A)') '   Mechanism', 'Leung (1991)'
      elseif (trim(obj_soot%model) == 'LIN') then
        write(*,'(A,T35,A)') '   Mechanism', 'Lindsted (1994)'
      end if
      if (nsoot == 1) write(*,'(A,T35,A)') '   Model', 'C(gr) => Heavy Gas'
    end if

    if (npass > 0) write(*,'(A,T35,I0)') '   Passive scalars', npass

    if (obj_rot%enabled) then
      write(*,*)
      write(*,'(A)')         ' Rotating Frame '
      write(*,'(A,F14.4,A)') '   omega  = ', obj_rot%omega,  ' rad/s'
      write(*,'(A,3F11.5)')  '   axis   = ', obj_rot%axis
      write(*,'(A,3F11.5)')  '   origin = ', obj_rot%origin
      if ( allocated(obj_rot%stationary_face) ) then
        write(*,'(A,I0,A)') '   stationary patches = ', size(obj_rot%stationary_face,2), ':'
        ! Loop counter `m` (not `k`): this is a logical counter (stationary RF
        ! patches), not a spatial z-direction index.
        do m = 1, size(obj_rot%stationary_face, 2)
          write(*,'(A,I0,A,I0)') '     block ', obj_rot%stationary_face(1,m), &
                                  '  face-dir ', obj_rot%stationary_face(2,m)
        end do
      else
        write(*,'(A)') '   all type-5/6 walls rotate with the frame'
      end if
    end if

    ! Geometric source terms (Q2D-extra, peculiarità #2/#3: area-variation V·G /
    ! axisymmetric 1/r — the MOSE peer has no numerics/sources/). Reported so the
    ! active Q2D-specific regime is visible in the run header.
    if (obj_source%axisymmetric .or. obj_source%area_variation) then
      write(*,*)
      write(*,'(A)') ' Geometric sources'
      if (obj_source%axisymmetric) &
        write(*,'(A,T35,A)') '   Model', 'Axisymmetric (1/r)'
      if (obj_source%area_variation) &
        write(*,'(A,T35,A)') '   Model', 'Area-variation (V.G, Definition A)'
    end if

    write(*,'(A)') ' ========================================================================================='

    ! ----- Numerical scheme -----
    write(*,*)
    write(*,'(A)') ' ========================================================================================='
    write(*,'(A)') ' Numerical scheme'
    write(*,'(A)') ' ========================================================================================='
    write(*,'(A,T35,A)') ' Space'
    write(*,'(A,T35,A)') '   Reconstruction', trim(obj_space_scheme%description)
    if (len_trim(obj_space_scheme%flux_limiter) > 0) &
      write(*,'(A,T35,A)') '   Flux limiter', trim(obj_space_scheme%flux_limiter)
    write(*,*)
    write(*,'(A,T35,A)') ' Time'
    write(*,'(A,T35,A)') '   Scheme', trim(obj_time_scheme%description)
    write(*,'(A,T35,A)') '   Integration variables', trim(obj_time_scheme%integration_variables)
    if (len_trim(obj_irs%description) > 0) &
      write(*,'(A,T35,A)') '   Implicit residual smoothing', trim(obj_irs%description)
    write(*,*)
    write(*,'(A,T35,A)') ' Fluxes'
    write(*,'(A,T35,A)') '   Riemann solver', trim(obj_riemann%description)
    ! OMP dispatcher (Q2D-extra, peculiarità #6: Assign_Fluxes('do'|'task')
    ! runtime-selectable — the MOSE peer has only !$omp do). Reported so the
    ! active Q2D-specific parallel strategy is visible in the run header.
    write(*,'(A,T35,A)') '   OMP parallelism', trim(obj_time_scheme%omp_strategy)
    write(*,'(A)') ' ========================================================================================='
    write(*,*)
  end subroutine Print_Shell_Info


  subroutine Write_Diagnostic(domain, IOfield, file)
    ! 2D-reduced (no k-loop / Nk; IOfield k-index = 1 fixed).
    use Q2D_Advanced_Types_m,  only: Q2D_domain_type
    use Q2D_Config_Types_m,    only: obj_io
    use Q2D_Global_m,          only: nsc, nu, np, nres
    use Q2D_Mod_MPI,           only: mpi_is_root
    use Q2D_Mod_GhostExchange, only: gather_diagnostic_to_root, mpi_io_barrier
    use Lib_VTK,               only: vtk_write_structured_multiblock
    use Lib_Tecplot,           only: tec_write_structured_multiblock
    use strings,               only: parse
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    type(ORION_data),      intent(inout) :: IOfield
    character(llen),       intent(in)    :: file
    ! Local
    character(len=llen) :: path, localpath_vtk
    integer  :: E_IO, b, i, j
    character(len=clen) :: format(2)
    character(len=6)    :: extension

    ! Gather R, dtlocal, beta from all ranks to root (collective).
    call gather_diagnostic_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'

      format = ''
      call parse(obj_io%sol_format, ' ', format)

      ! Diagnostic residual field — EXACT mirror of the MOSE peer
      ! (archive-MOSE/src/lib/diagnostic/Mod_Diagnostic.f90:84-96), 2D-reduced:
      ! MOSE writes vars(2:5)=r(nu:np) (4-wide, 3D mean-flow incl. rhow); Q2D drops
      ! nw -> vars(2:4)=R(nu:np) (3-wide). dt at slot 5 (MOSE 6), beta at slot 6
      ! (MOSE 7). The RANS turbulence block is commented out exactly as in MOSE.
      ! Nvars = nres+2 mirrors MOSE :113. The SOTA-coherent form (write the full
      ! np+1:nprim tail + dynamic Dvarnames/Nvars) is a back-port to do in lockstep
      ! on BOTH repos -> TODO (MOSE READ-ONLY, no unilateral Q2D divergence).
      do b = 1, size(IOfield%block)
        IOfield%block(b)%vars(1, :, :, 1) = sum( &
          domain%blk(b)%R(1:nsc, 1:IOfield%block(b)%Ni, 1:IOfield%block(b)%Nj), dim=1 )
        IOfield%block(b)%vars(2:4, :, :, 1) = &
          domain%blk(b)%R(nu:np, 1:IOfield%block(b)%Ni, 1:IOfield%block(b)%Nj)
        !if ( nrans > 0 ) then ! turbulence variables
        !  IOfield%block(b)%vars(5:5+nrans-1, :, :, 1) = &
        !    domain%blk(b)%R(nt:nprim, 1:IOfield%block(b)%Ni, 1:IOfield%block(b)%Nj)
        !end if
        ! Auxiliary variables: dt
        do j = 1, IOfield%block(b)%Nj
        do i = 1, IOfield%block(b)%Ni
          IOfield%block(b)%vars(5, i, j, 1) = domain%blk(b)%dtlocal(i, j)
        enddo; enddo
        ! Auxiliary variables: beta (shock detector)
        do j = 1, IOfield%block(b)%Nj
        do i = 1, IOfield%block(b)%Ni
          IOfield%block(b)%vars(6, i, j, 1) = domain%blk(b)%beta(i, j)
        enddo; enddo
      enddo

      ! Write the IOfield in the chosen solution format
      select case(trim(format(1)))
      case('vtk')
        IOfield%vtk%format = trim(format(2))
        localpath_vtk = trim(path)//'/vtk'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield, &
                                                vtspath=trim(localpath_vtk)//trim(file), &
                                                vtmpath=trim(path)//trim(file), &
                                                varnames=Dvarnames, time=domain%time)
      case('tecplot')
        IOfield%tec%format = trim(format(2))
        ! Extension conditional: `.szplt` for binary tecplot SZPLT format,
        ! `.tec` for ASCII text.
        if (trim(format(2)) == 'binary') then
          extension = '.szplt'
        else
          extension = '.tec'
        end if
        E_IO = tec_write_structured_multiblock(Nvars=nres+2, orion=IOfield, &
                                                 varnames=Dvarnames, &
                                                 filename=trim(path)//trim(file)//trim(extension))
      end select
    end if

    ! Synchronize all ranks after I/O.
    call mpi_io_barrier()
  end subroutine Write_Diagnostic


  subroutine Compute_Residual(new, old, dt, n, average, total)
    ! 2D-reduced (no k-index/Nk; collapse(2) instead of collapse(3)).
    use Q2D_Global_m, only: nsc, nu, nprim, nres, gc
    implicit none
    integer,  intent(in)    :: n(2)
    real(R8), intent(in)    :: new(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc)
    real(R8), intent(in)    :: old(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc)
    real(R8), intent(in)    :: dt(n(1), n(2))
    real(R8), intent(out)   :: average(nres)
    real(R8), intent(inout) :: total(nres)
    ! Local
    integer  :: i, j
    real(R8) :: resn(nres), residuo(nres), reslocal(nres, n(1), n(2))

    residuo = 0d0

    !$omp parallel
    !$omp do private(resn) reduction(+:residuo) collapse(2)
    do j = 1, n(2)
    do i = 1, n(1)
      resn(1)      = abs( sum(new(1:nsc, i, j) - old(1:nsc, i, j)) )    ! delta(rho)
      resn(2:nres) = abs( new(nu:nprim, i, j) - old(nu:nprim, i, j) )    ! delta(...)
      reslocal(:, i, j) = resn
      residuo = residuo + resn*resn
    enddo; enddo
    !$omp end parallel

    ! Average residual of this block
    average = sqrt( residuo / real(n(1)*n(2), R8) )
    ! Overall residual
    total = total + residuo
  end subroutine Compute_Residual

end module Q2D_Mod_Diagnostic
