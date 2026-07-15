module Q2D_IO_Wall
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Lib_ORION_data,  only: ORION_data

  implicit none
  private
  public :: IOwall
  public :: Initialize_Wall_File
  public :: Write_Wall_Solution
  public :: Extract_Properties_From_BC

  type(ORION_data), save :: IOwall

contains

  subroutine Initialize_Wall_File ( domain, extension )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_io, obj_io_bc
    use IR_Precision,         only: str
    implicit none
    type(Q2D_domain_type), intent(in) :: domain
    character(len=*),      intent(in) :: extension
    ! Variabili locali
    integer :: b, f, bb, i, j

    IOwall%tec%extension = extension

    obj_io%write_wall_mechanical = .false.
    obj_io%write_wall_thermal    = .false.
    obj_io%write_wall_mass       = .false.

    ! Intestazione: nessuna "tauZ" (nessuna componente di sforzo fuori dal piano).
    obj_io%Onwall = 1
    obj_io%Owallnames = '"y+"'
    if (index(obj_io%wall_variables, 'mech') > 0) then
      obj_io%Onwall = obj_io%Onwall + 3
      obj_io%Owallnames = trim(obj_io%Owallnames)//' "tauX" "tauY" "pw"'
      obj_io%write_wall_mechanical = .true.
    end if
    if (index(obj_io%wall_variables, 'thermal') > 0) then
      obj_io%Onwall = obj_io%Onwall + 2
      obj_io%Owallnames = trim(obj_io%Owallnames)//' "Tw" "qw"'
      obj_io%write_wall_thermal = .true.
    end if
    if (index(obj_io%wall_variables, 'mass') > 0) then
      obj_io%Onwall = obj_io%Onwall + 1
      obj_io%Owallnames = trim(obj_io%Owallnames)//' "mdot"'
      obj_io%write_wall_mass = .true.
    end if

    ! Alloca un blocco IOwall per ogni faccia viscosa (4 facce per blocco nel 2D).
    allocate( IOwall%block(1:count(obj_io_bc%viscous_flag)) )

    bb = 0
    do b = 1, domain%nb
      do f = 1, 4
        if (obj_io_bc%viscous_flag(b, f)) then
          bb = bb + 1
          IOwall%block(bb)%name = 'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
          select case (f)
            case (1)
              allocate( IOwall%block(bb)%mesh(1:2, 0:0, 0:domain%blk(b)%dim(2), 0:0) )
              allocate( IOwall%block(bb)%vars(1:obj_io%Onwall, 1, 1:domain%blk(b)%dim(2), 1) )
              do j = 0, domain%blk(b)%dim(2)
                IOwall%block(bb)%mesh(1, 0, j, 0) = domain%blk(b)%node(0, j)%c(1)
                IOwall%block(bb)%mesh(2, 0, j, 0) = domain%blk(b)%node(0, j)%c(2)
              end do
            case (2)
              allocate( IOwall%block(bb)%mesh(1:2, 0:0, 0:domain%blk(b)%dim(2), 0:0) )
              allocate( IOwall%block(bb)%vars(1:obj_io%Onwall, 1, 1:domain%blk(b)%dim(2), 1) )
              do j = 0, domain%blk(b)%dim(2)
                IOwall%block(bb)%mesh(:, 0, j, 0) = domain%blk(b)%node(domain%blk(b)%dim(1), j)%c
              end do
            case (3)
              allocate( IOwall%block(bb)%mesh(1:2, 0:domain%blk(b)%dim(1), 0:0, 0:0) )
              allocate( IOwall%block(bb)%vars(1:obj_io%Onwall, 1:domain%blk(b)%dim(1), 1, 1) )
              do i = 0, domain%blk(b)%dim(1)
                IOwall%block(bb)%mesh(:, i, 0, 0) = domain%blk(b)%node(i, 0)%c
              end do
            case (4)
              allocate( IOwall%block(bb)%mesh(1:2, 0:domain%blk(b)%dim(1), 0:0, 0:0) )
              allocate( IOwall%block(bb)%vars(1:obj_io%Onwall, 1:domain%blk(b)%dim(1), 1, 1) )
              do i = 0, domain%blk(b)%dim(1)
                IOwall%block(bb)%mesh(:, i, 0, 0) = domain%blk(b)%node(i, domain%blk(b)%dim(2))%c
              end do
          end select
        end if
      end do
    end do

  end subroutine Initialize_Wall_File


  subroutine Write_Wall_Solution ( domain, file )
    use Q2D_Advanced_Types_m,  only: Q2D_domain_type
    use Q2D_Config_Types_m,    only: obj_io
    use Q2D_Parameters_m,      only: llen, clen
    use Q2D_Mod_MPI,           only: mpi_is_root
    use Q2D_Mod_GhostExchange, only: mpi_io_barrier
    use Lib_VTK,               only: vtk_write_structured_multiblock
    use Lib_Tecplot,           only: tec_write_structured_multiblock
    use strings,               only: parse
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    character(len=llen),   intent(in)    :: file
    ! Variabili locali
    character(len=llen) :: path, localpath_vtk
    integer             :: E_IO
    character(len=clen) :: format(2)

    ! Ogni rank calcola le proprietà a parete per i propri blocchi locali e le
    ! riduce sul root. Chiamata collettiva — tutti i rank la eseguono.
    call Extract_Properties_From_BC(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'
      format = ''
      call parse(obj_io%sol_format, ' ', format)

      select case (trim(format(1)))
        case ('vtk')
          IOwall%vtk%format = trim(format(2))
          localpath_vtk = trim(path)//'vtk/'
          call execute_command_line('mkdir -p '//trim(localpath_vtk))
          E_IO = vtk_write_structured_multiblock( &
                   orion    = IOwall, &
                   vtspath  = trim(localpath_vtk)//trim(file), &
                   vtmpath  = trim(path)//trim(file), &
                   varnames = obj_io%Owallnames, &
                   time     = domain%time )
        case ('tecplot')
          IOwall%tec%format = trim(format(2))
          E_IO = tec_write_structured_multiblock( &
                   orion    = IOwall, &
                   varnames = obj_io%Owallnames, &
                   filename = trim(path)//trim(file)//trim(IOwall%tec%extension) )
      end select
    end if

    ! Sincronizzazione di tutti i rank dopo l'I/O.
    call mpi_io_barrier()

  end subroutine Write_Wall_Solution


  subroutine Extract_Properties_From_BC ( domain )
    !> Calcola le proprietà delle BC di parete su tutti i rank (solo blocchi
    !> locali), riduce sul root e riempie IOwall sul root. Chiamata collettiva —
    !> tutti i rank devono eseguirla. Schema in 5 fasi (7 variabili per cella,
    !> 4 facce per blocco nel 2D).
    use Q2D_Advanced_Types_m,               only: Q2D_domain_type
    use Q2D_Config_Types_m,                 only: obj_io, obj_io_bc
    use Q2D_Mod_MPI,                        only: mpi_is_root, is_local_block, mpi_reduce_sum_r8_array
    use Q2D_Lib_BC_Fluxes_Wall_Heat,        only: BC_Wall_Heat
    use Q2D_Lib_BC_Fluxes_Wall_Temperature, only: BC_Wall_Temperature
    use Q2D_Lib_BC_Fluxes_Wall_Melting,     only: BC_Wall_Melting
    use Q2D_Lib_BC_Fluxes_Wall_Pyrolysis,   only: BC_Wall_Pyrolysis
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Variabili locali
    real(R8), dimension(7)  :: Ovar
    integer :: f, lower, upper, i, b, bb, pv, Bm, Im, Jm, Fm, Is, Js
    integer :: n_wall_total, wall_i
    real(R8), allocatable :: wall_buf(:)

    ! FASE 1 — Conteggio delle celle BC totali su tutte le facce viscose
    n_wall_total = 0
    do b = 1, domain%nb
      do f = 1, 4
        if (obj_io_bc%viscous_flag(b,f)) n_wall_total = n_wall_total + domain%n_bf(b,f)
      end do
    end do

    ! FASE 2 — Allocazione del buffer piatto 7*n_wall_total (7 variabili per cella)
    allocate(wall_buf(7 * n_wall_total))
    wall_buf = 0d0

    ! FASE 3 — Ogni rank calcola Ovar per i propri blocchi locali; il salto dei
    ! blocchi remoti lascia zero in wall_buf, così la riduzione SUM dà valori
    ! corretti (ogni cella appartiene esattamente a un rank).
    upper = 0; wall_i = 0
    do b = 1, domain%nb
      do f = 1, 4
        lower = upper + 1
        upper = upper + domain%n_bf(b, f)
        if (.not. obj_io_bc%viscous_flag(b, f)) cycle

        do i = lower, upper
          wall_i = wall_i + 1
          Bm = domain%bc(i)%b
          if (.not. is_local_block(Bm)) cycle  ! blocco remoto: lascia zero nel buffer

          Im = domain%bc(i)%i
          Jm = domain%bc(i)%j
          Fm = domain%bc(i)%f
          Ovar = 0.0_R8
          select case (domain%bc(i)%type)
            case (301) ! parete, qw assegnato
              call BC_Wall_Heat ( Im, Jm, Fm, domain%blk(Bm), domain%bc(i)%qw, Ovar )
            case (302) ! parete, T assegnata
              call BC_Wall_Temperature ( Im, Jm, Fm, domain%blk(Bm), domain%bc(i)%Tw, Ovar )
            case (303) ! parete, T_wall + qrad (fusione)
              call BC_Wall_Melting ( Im, Jm, Fm, domain%blk(Bm),         &
                                      domain%bc(i)%Tw, domain%bc(i)%qrad, &
                                      Ovar )
            case (304) ! parete, qrad (pirolisi/ablazione); fuel = bc%type → ramo HDPE
              call BC_Wall_Pyrolysis ( Im, Jm, Fm, domain%blk(Bm), &
                                        domain%bc(i)%qrad,         &
                                        domain%bc(i)%type, Ovar )
          end select
          wall_buf(7*(wall_i-1)+1 : 7*wall_i) = Ovar
        end do
      end do
    end do

    ! FASE 4 — Riduzione di wall_buf sul root. Ogni cella appartiene esattamente
    ! a un rank, quindi SUM dà valori corretti. Sul root, wall_buf è aggiornato
    ! in place da mpi_reduce_sum_r8_array.
    call mpi_reduce_sum_r8_array(wall_buf, 7 * n_wall_total)

    ! FASE 5 — Il root spacchetta wall_buf in IOwall. L'ordine dei cicli DEVE
    ! coincidere col riempimento della Fase 3 (stessa sequenza b/f/i) per una
    ! mappatura corretta.
    if (mpi_is_root) then
      upper = 0; wall_i = 0; bb = 0
      do b = 1, domain%nb
        do f = 1, 4
          lower = upper + 1
          upper = upper + domain%n_bf(b, f)
          if (.not. obj_io_bc%viscous_flag(b, f)) cycle
          bb = bb + 1

          do i = lower, upper
            wall_i = wall_i + 1
            Ovar = wall_buf(7*(wall_i-1)+1 : 7*wall_i)

            select case (f)
              case (1, 2); Is = 1               ; Js = domain%bc(i)%j
              case (3, 4); Is = domain%bc(i)%i  ; Js = 1
            end select
            IOwall%block(bb)%vars(1, Is, Js, 1) = Ovar(1)
            pv = 1
            if (obj_io%write_wall_mechanical) then
              IOwall%block(bb)%vars(pv+1:pv+3, Is, Js, 1) = Ovar(2:4)
              pv = pv + 3
            end if
            if (obj_io%write_wall_thermal) then
              IOwall%block(bb)%vars(pv+1:pv+2, Is, Js, 1) = Ovar(5:6)
              pv = pv + 2
            end if
            if (obj_io%write_wall_mass) then
              IOwall%block(bb)%vars(pv+1, Is, Js, 1) = Ovar(7)
              pv = pv + 1  ! incremento di 1: una sola variabile assegnata
              ! (pv non è più usato dopo questo ramo).
            end if
          end do
        end do
      end do
    end if

    deallocate(wall_buf)

  end subroutine Extract_Properties_From_BC

end module Q2D_IO_Wall
