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
    ! Local variables
    integer :: b, f, bb, i, j

    IOwall%tec%extension = extension

    obj_io%write_wall_mechanical = .false.
    obj_io%write_wall_thermal    = .false.
    obj_io%write_wall_mass       = .false.

    ! Header: no "tauZ" (no out-of-plane stress component).
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

    ! Allocate one IOwall block for each viscous face (4 faces per block in 2D).
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
    ! Local variables
    character(len=llen) :: path, localpath_vtk
    integer             :: E_IO
    character(len=clen) :: format(2)

    ! Each rank computes the wall properties for its own local blocks and
    ! reduces them onto the root. Collective call — all ranks execute it.
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

    ! Synchronize all ranks after the I/O.
    call mpi_io_barrier()

  end subroutine Write_Wall_Solution


  subroutine Extract_Properties_From_BC ( domain )
    !> Computes the wall BC properties on all ranks (local blocks only),
    !> reduces onto the root and fills IOwall on the root. Collective call —
    !> all ranks must execute it. 5-phase scheme (7 variables per cell,
    !> 4 faces per block in 2D).
    use Q2D_Advanced_Types_m,               only: Q2D_domain_type
    use Q2D_Config_Types_m,                 only: obj_io, obj_io_bc
    use Q2D_Mod_MPI,                        only: mpi_is_root, is_local_block, mpi_reduce_sum_r8_array
    use Q2D_Lib_BC_Fluxes_Wall_Heat,        only: BC_Wall_Heat
    use Q2D_Lib_BC_Fluxes_Wall_Temperature, only: BC_Wall_Temperature
    use Q2D_Lib_BC_Fluxes_Wall_Melting,     only: BC_Wall_Melting
    use Q2D_Lib_BC_Fluxes_Wall_Pyrolysis,   only: BC_Wall_Pyrolysis
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    real(R8), dimension(7)  :: Ovar
    integer :: f, lower, upper, i, b, bb, pv, Bm, Im, Jm, Fm, Is, Js
    integer :: n_wall_total, wall_i
    real(R8), allocatable :: wall_buf(:)

    ! PHASE 1 — Count the total BC cells over all viscous faces
    n_wall_total = 0
    do b = 1, domain%nb
      do f = 1, 4
        if (obj_io_bc%viscous_flag(b,f)) n_wall_total = n_wall_total + domain%n_bf(b,f)
      end do
    end do

    ! PHASE 2 — Allocate the flat buffer 7*n_wall_total (7 variables per cell)
    allocate(wall_buf(7 * n_wall_total))
    wall_buf = 0d0

    ! PHASE 3 — Each rank computes Ovar for its own local blocks; skipping the
    ! remote blocks leaves zero in wall_buf, so the SUM reduction gives correct
    ! values (each cell belongs to exactly one rank).
    upper = 0; wall_i = 0
    do b = 1, domain%nb
      do f = 1, 4
        lower = upper + 1
        upper = upper + domain%n_bf(b, f)
        if (.not. obj_io_bc%viscous_flag(b, f)) cycle

        do i = lower, upper
          wall_i = wall_i + 1
          Bm = domain%bc(i)%b
          if (.not. is_local_block(Bm)) cycle  ! remote block: leave zero in the buffer

          Im = domain%bc(i)%i
          Jm = domain%bc(i)%j
          Fm = domain%bc(i)%f
          Ovar = 0.0_R8
          select case (domain%bc(i)%type)
            case (301) ! wall, qw prescribed
              call BC_Wall_Heat ( Im, Jm, Fm, domain%blk(Bm), domain%bc(i)%qw, Ovar )
            case (302) ! wall, T prescribed
              call BC_Wall_Temperature ( Im, Jm, Fm, domain%blk(Bm), domain%bc(i)%Tw, Ovar )
            case (303) ! wall, T_wall + qrad (melting)
              call BC_Wall_Melting ( Im, Jm, Fm, domain%blk(Bm),         &
                                      domain%bc(i)%Tw, domain%bc(i)%qrad, &
                                      Ovar )
            case (304) ! wall, qrad (pyrolysis/ablation); fuel = bc%type → HDPE branch
              call BC_Wall_Pyrolysis ( Im, Jm, Fm, domain%blk(Bm), &
                                        domain%bc(i)%qrad,         &
                                        domain%bc(i)%type, Ovar )
          end select
          wall_buf(7*(wall_i-1)+1 : 7*wall_i) = Ovar
        end do
      end do
    end do

    ! PHASE 4 — Reduce wall_buf onto the root. Each cell belongs to exactly one
    ! rank, so SUM gives correct values. On the root, wall_buf is updated in
    ! place by mpi_reduce_sum_r8_array.
    call mpi_reduce_sum_r8_array(wall_buf, 7 * n_wall_total)

    ! PHASE 5 — The root unpacks wall_buf into IOwall. The loop order MUST
    ! match the Phase 3 filling (same b/f/i sequence) for a correct mapping.
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
              pv = pv + 1  ! increment by 1: a single variable assigned
              ! (pv is no longer used after this branch).
            end if
          end do
        end do
      end do
    end if

    deallocate(wall_buf)

  end subroutine Extract_Properties_From_BC

end module Q2D_IO_Wall
