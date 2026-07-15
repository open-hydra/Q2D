module Q2D_IO_Solution
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Input_Solution
  public :: Read_IC
  public :: Setup_Output_Solution
  public :: Write_solution

  !> Puntatore a procedura verso il writer di soluzione selezionato
  procedure(w_solution_if), pointer :: Write_solution => null()

  abstract interface
    subroutine w_solution_if ( domain, IOfield, file )
      use Q2D_Advanced_Types_m, only: Q2D_domain_type
      use Q2D_Parameters_m,     only: llen
      use Lib_ORION_data,       only: ORION_data
      implicit none
      type(Q2D_domain_type), intent(inout) :: domain
      type(ORION_data),      intent(inout) :: IOfield
      character(llen),       intent(in)    :: file
    end subroutine w_solution_if
  end interface

contains

  ! Wrapper vuoto (no-op), mantenuto per completezza dell'interfaccia pubblica:
  ! il chiamante lo invoca prima di Read_IC, ma tutta la logica di ricerca del
  ! file e di selezione del formato del reader (ini_format su nuova run,
  ! sol_format su restart) è già incorporata nella subroutine concreta Read_IC
  ! (vedi sotto).
  subroutine Setup_Input_Solution()
    implicit none
  end subroutine Setup_Input_Solution


  subroutine Read_IC ( IOfield )
    use Q2D_Config_Types_m, only: obj_sim_param, obj_io
    use Q2D_Global_m,       only: Q2D_phase_prefix
    use Q2D_Parameters_m,   only: llen, clen
    use Lib_ORION_data,     only: ORION_data
    use Lib_Tecplot,        only: tec_read_structured_multiblock
    use Lib_VTK,            only: vtk_read_structured_multiblock
    use IR_Precision,       only: str
    use strings,            only: parse
    implicit none
    type(ORION_data), intent(inout) :: IOfield
    ! Variabili locali
    logical             :: present
    integer             :: i, error
    character(len=clen) :: format(2)
    character(len=6)    :: extension
    character(len=llen) :: try

    ! Selezione del formato del reader: ini_format su nuova run, sol_format su
    ! restart. Il ramo restart più sotto deduce già l'estensione del file da
    ! sol_format; ricavare il dispatch del reader dalla stessa stringa mantiene
    ! coerenti estensione e reader quando ini_format /= sol_format.
    format = ''
    if (obj_sim_param%newrun) then
      call parse(obj_io%ini_format, ' ', format)
    else
      call parse(obj_io%sol_format, ' ', format)
    end if

    if (obj_sim_param%newrun) then
      if (index(obj_io%ini_format, 'vtk') > 0) then
        extension = '.vtm'
      else
        if (index(obj_io%ini_format, 'ascii') > 0) then
          extension = '.tec'
        else
          extension = '.szplt'
        end if
      end if
      obj_io%nameinit = 'INPUT/'//trim(Q2D_phase_prefix)//'ic'//trim(extension)
      inquire(file=obj_io%nameinit, exist=present)
    else
      if (index(obj_io%sol_format, 'vtk') > 0) then
        extension = '.vtm'
      else
        if (index(obj_io%sol_format, 'ascii') > 0) then
          extension = '.tec'
        else
          extension = '.szplt'
        end if
      end if
      obj_io%nameinit = 'OUTPUT/'//trim(Q2D_phase_prefix)//'field'//trim(extension)
      inquire(file=obj_io%nameinit, exist=present)
      if (.not. present) then
        i = 0
        do
          i = i + 1
          try = 'OUTPUT/'//trim(Q2D_phase_prefix)//'field'//trim(str(.true.,i))//trim(extension)
          inquire(file=try, exist=present)
          if (present) then
            obj_io%nameinit = try
          else
            exit
          end if
        end do
      end if
    end if

    obj_io%IOtime = 0.0_R8
    error = 0

    if (.not. present) then
      obj_io%error_message = '[ERROR] reading input file '//trim(obj_io%nameinit)
      return
    end if

    select case (trim(format(1)))
      case ('tecplot')
        IOfield%tec%format = trim(format(2))
        error = tec_read_structured_multiblock(orion=IOfield, filename=trim(obj_io%nameinit))
      case ('vtk')
        IOfield%vtk%format = trim(format(2))
        error = vtk_read_structured_multiblock( &
                  orion   = IOfield, &
                  vtmpath = obj_io%nameinit(1:len(trim(obj_io%nameinit))-4), &
                  vtspath = 'INPUT/vtk/field', &
                  time    = obj_io%IOtime )
    end select

    if (error /= 0) obj_io%error_message = '[ERROR] reading input file '//trim(obj_io%nameinit)

  end subroutine Read_IC


  subroutine Setup_Output_Solution ( IOfield )
    use Q2D_Config_Types_m, only: obj_multigrid, obj_io
    use Q2D_Global_m,       only: nsc, nrans, nprim, npass, nsoot
    use Q2D_Parameters_m,   only: llen
    use Lib_ORION_data,     only: ORION_data
    use IR_Precision,       only: str
    implicit none
    type(ORION_data), intent(inout) :: IOfield(:)
    ! Variabili locali
    integer :: b, i, m

    ! Specifica delle variabili di IO — nessuna "w" (nessuna velocità fuori dal piano).
    obj_io%Ovarnames = ' '
    obj_io%ORANSname = ' '
    if (nrans == 1) then
      obj_io%ORANSname = ' "mi_tilde"'
    else if (nrans == 2) then
      obj_io%ORANSname = ' "kappa" "omega"'
    else if (nrans == 5) then
      ! RSM SSG/LRR-omega: 4 tensioni di Reynolds mantenute + omega (R13 e R23
      ! eliminate: da 7 a 5 variabili). Nomi senza apostrofi (ruu = r u'u', ecc.)
      ! per compatibilità con i parser Tecplot/ORION.
      obj_io%ORANSname = ' "ruu" "rvv" "rww" "ruv" "omega"'
    end if
    do i = 1, nsc
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"rho('//trim(str(.true.,i))//')"'
    end do
    obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "u" "v" "p"'
    if (nsoot == 1) then
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "rho*Np"'
    end if
    if (nsoot == 2) then
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "rho*Yp" "rho*Np"'
    end if
    if (npass > 0) then
      do i = 1, npass
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//'"Pass('//trim(str(.true.,i))//')"'
      end do
    end if
    obj_io%Ovarnames = trim(obj_io%Ovarnames)//trim(obj_io%ORANSname)

    ! Variabili ausiliarie: T è scritta di default.
    obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "T"'
    obj_io%Onvar = nprim + 1
    obj_io%write_thermo = .false.
    if (index(obj_io%sol_variables, 'thermo') > 0) then
      obj_io%write_thermo = .true.
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "g" "R"'
      obj_io%Onvar = obj_io%Onvar + 2
    end if
    obj_io%write_transport = .false.
    if (index(obj_io%sol_variables, 'transport') > 0) then
      obj_io%write_transport = .true.
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "mil" "kl"'
      obj_io%Onvar = obj_io%Onvar + 2
    end if
    if (obj_io%write_transport .and. (nrans > 0)) then
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "mit"'
      obj_io%Onvar = obj_io%Onvar + 1
    end if

    do m = 1, obj_multigrid%MGL
      IOfield(m)%vtk%node = .false.
      IOfield(m)%tec%node = .false.
      IOfield(m)%tec%bc   = .false.
      if (index(obj_io%sol_format, 'ascii') > 0) then
        IOfield(m)%tec%extension = '.tec'
        IOfield(m)%tec%format    = 'ascii'
        IOfield(m)%vtk%format    = 'ascii'
      else
        IOfield(m)%tec%extension = '.szplt'
        IOfield(m)%tec%format    = 'binary'
        IOfield(m)%vtk%format    = 'binary'
      end if
    end do

    ! Concretizzazione del writer
    Write_solution => Write_vtk_tec

    do m = 1, obj_multigrid%MGL
      do b = 1, size(IOfield(m)%block)
        IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
      end do
    end do

    ! Rialloca le vars di IOfield se Onvar è cambiato (es. backup scritto con un set di variabili diverso).
    do m = 1, obj_multigrid%MGL
      if (obj_io%Onvar /= size(IOfield(m)%block(1)%vars, 1)) then
        do b = 1, size(IOfield(m)%block)
          IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
          deallocate(IOfield(m)%block(b)%vars)
          allocate(IOfield(m)%block(b)%vars(1:obj_io%Onvar, &
                                            1:IOfield(m)%block(b)%Ni, &
                                            1:IOfield(m)%block(b)%Nj, &
                                            1))   ! Nk = 1
        end do
      end if
    end do

  end subroutine Setup_Output_Solution


  subroutine Write_vtk_tec ( domain, IOfield, file )
    use Q2D_Advanced_Types_m,    only: Q2D_domain_type
    use Q2D_Config_Types_m,      only: obj_io
    use Q2D_Global_m,            only: nsc, nu, nv, np, nt, nprim, nrans
    use Q2D_Parameters_m,        only: llen
    use Q2D_Lib_RANS,            only: Eddy_Viscosity
    use Q2D_Mod_MPI,             only: mpi_is_root
    use Q2D_Mod_GhostExchange,   only: gather_P_to_root, mpi_io_barrier
    use Lib_ORION_data,          only: ORION_data
    use Lib_VTK,                 only: vtk_write_structured_multiblock
    use Lib_Tecplot,             only: tec_write_structured_multiblock
    use FLINT_Lib_Thermodynamic, only: f_Rtot, f_gamma, eos, co_k_mi_lam_Wilke
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    type(ORION_data),      intent(inout) :: IOfield
    character(llen),       intent(in)    :: file
    ! Variabili locali
    character(len=llen) :: path, localpath_vtk
    integer  :: E_IO, b, i, j, pv
    real(R8) :: p, rho, Rtot, rhoi(nsc), T, mu, kl, mut, r_ave
    real(R8) :: vel_gradient(2,2)

    ! Raccolta di P da tutti i rank sul root (chiamata collettiva — tutti i rank devono eseguirla).
    call gather_P_to_root(domain)

    if (mpi_is_root) then
    path = 'OUTPUT/'

    ! Aggiorna le variabili di IOfield con le primitive del dominio + le variabili ausiliarie.
    do b = 1, size(IOfield%block)
      IOfield%block(b)%vars(1:nprim_local(),:,:,1) = &
        domain%blk(b)%P(:, 1:IOfield%block(b)%Ni, 1:IOfield%block(b)%Nj)
      do j = 1, IOfield%block(b)%Nj
        do i = 1, IOfield%block(b)%Ni
          rhoi = domain%blk(b)%P(1:nsc, i, j)
          p    = domain%blk(b)%P(np,    i, j)
          rho  = sum(rhoi)
          Rtot = f_Rtot(rhoi)
          T    = eos(p=p, rho=rho, R=Rtot)
          pv = 1
          IOfield%block(b)%vars(nprim_local()+1, i, j, 1) = T
          if (obj_io%write_thermo) then
            pv = pv + 1
            IOfield%block(b)%vars(nprim_local()+pv, i, j, 1) = f_gamma(rhoi, p, rho, Rtot)
            pv = pv + 1
            IOfield%block(b)%vars(nprim_local()+pv, i, j, 1) = Rtot
          end if
          if (obj_io%write_transport) then
            call co_k_mi_lam_Wilke(rhoi, rho, T, mu, kl)
            pv = pv + 1
            IOfield%block(b)%vars(nprim_local()+pv, i, j, 1) = mu
            pv = pv + 1
            IOfield%block(b)%vars(nprim_local()+pv, i, j, 1) = kl
            if (nrans > 0) then
              vel_gradient(:,1) = ( domain%blk(b)%P(nu:nv,i+1,j) - domain%blk(b)%P(nu:nv,i-1,j) ) / 2d0
              vel_gradient(:,2) = ( domain%blk(b)%P(nu:nv,i,j+1) - domain%blk(b)%P(nu:nv,i,j-1) ) / 2d0
              vel_gradient = matmul ( vel_gradient + 1d-40, domain%blk(b)%M(i,j)%c )
              r_ave = 0.25d0 * ( domain%blk(b)%node(i,  j  )%c(2) + domain%blk(b)%node(i-1,j  )%c(2) + &
                                 domain%blk(b)%node(i,  j-1)%c(2) + domain%blk(b)%node(i-1,j-1)%c(2) )
              call Eddy_Viscosity ( mut=mut, rans_variables=domain%blk(b)%P(nt:nprim,i,j),     &
                                    mul=mu, rho=rho, vel_gradient=vel_gradient,                  &
                                    vr=domain%blk(b)%P(nv,i,j), r_ave=r_ave,                     &
                                    walldist=domain%blk(b)%yn(i,j) )
              pv = pv + 1
              IOfield%block(b)%vars(nprim_local()+pv, i, j, 1) = mut
            end if
          end if
        end do
      end do
    end do

    ! Scrittura di IOfield nel formato di output selezionato.
    if (index(obj_io%sol_format, 'vtk') > 0) then
      localpath_vtk = trim(path)//'vtk/'
      call execute_command_line('mkdir -p '//trim(localpath_vtk))
      E_IO = vtk_write_structured_multiblock( &
               orion    = IOfield, &
               vtspath  = trim(localpath_vtk)//trim(file), &
               vtmpath  = trim(path)//trim(file), &
               varnames = obj_io%Ovarnames, &
               time     = domain%time )
    else
      E_IO = tec_write_structured_multiblock( &
               orion    = IOfield, &
               varnames = obj_io%Ovarnames, &
               filename = trim(path)//trim(file)//trim(IOfield%tec%extension) )
    end if
    end if  ! mpi_is_root

    ! Sincronizzazione di tutti i rank dopo l'I/O.
    call mpi_io_barrier()

  end subroutine Write_vtk_tec


  ! Helper locale per evitare lo shadowing di `use Q2D_Global_m, only: nprim`
  ! dentro Write_vtk_tec (che già usa nu/nv/np). nprim è letto qui direttamente
  ! dal modulo globale.
  pure function nprim_local() result(n)
    use Q2D_Global_m, only: nprim
    integer :: n
    n = nprim
  end function nprim_local

end module Q2D_IO_Solution
