module Q2D_Read_IO
  use iso_fortran_env,    only: I4 => int32, R8 => real64
  use Q2D_Registry_m,     only: reg
  use Q2D_Config_Types_m, only: obj_io, obj_io_probes, obj_io_lineprobes

  implicit none
  private
  public :: Register_IO_Fields, Register_Probes, Register_LineProbes

contains

  subroutine Register_IO_Fields()
    use Q2D_Parameters_m, only: codename
    implicit none
    character(len=:), allocatable :: section

    section = trim(codename)//'-IO'

    obj_io%warning_message = 'none'
    obj_io%error_message   = 'none'
    obj_io%description     = 'none'

    call reg%add( section, 'ini-format', obj_io%ini_format, 'tecplot ascii', 'Initial condition (INPUT/ic.*) format', 'tecplot ascii, tecplot binary, vtk ascii, vtk raw', .false. )
    call reg%add( section, 'sol-format', obj_io%sol_format, 'tecplot ascii', 'Solution (OUTPUT/field.*) format', 'tecplot ascii, tecplot binary, vtk ascii, vtk raw', .false. )

    call reg%add( section, 'sol-diter', obj_io%sol_diter, '1000000000', 'Solution output iter frequency', '> 0', .false. )
    call reg%add( section, 'sol-dtime', obj_io%sol_dtime, '1e30', 'Solution output time frequency', '> 0', .false. )
    call reg%add( section, 'sol-overwrite', obj_io%sol_overwrite, 'true', 'Overwrite solution files', 'true, false', .false. )

    ! Definizione delle variabili da scrivere
    call reg%add( section, 'sol-variables', obj_io%sol_variables, 'thermo', 'Solution variables to write', '', .false. )
    call reg%add( section, 'wall-variables', obj_io%wall_variables, 'mech thermal', 'Wall variables to write', '', .false. )

    ! File di storia dei residui
    call reg%add( section, 'res-diter', obj_io%res_diter, '1', 'Residual history iter frequency', '> 0', .false. )

    ! Opzioni della shell
    call reg%add( section, 'shell-diter', obj_io%shell_diter, '1', 'Shell update iter frequency', '> 0', .false. )
    call reg%add( section, 'ini-diter', obj_io%ini_diter, '10000', 'input.ini update iter frequency', '> 0', .false. )

  end subroutine Register_IO_Fields


  subroutine Register_Probes(n, probes_name)
    use Q2D_Parameters_m, only: codename
    use IR_Precision,     only: str
    implicit none
    integer, intent(in) :: n
    character(*), intent(in) :: probes_name(n)
    integer :: p

    allocate(obj_io_probes(n))
    do p = 1, n
      call reg%add( trim(codename)//'-Probes', 'probe'//trim(str(.true.,p)), obj_io_probes(p)%file, 'probe-placeholder', 'Probe file name', '', .false. )
      call Register_One_Probe(p, probes_name(p))
    end do

  end subroutine Register_Probes


  subroutine Register_One_Probe(p, probe_name)
    use IR_Precision, only: str
    implicit none
    integer, intent(in) :: p
    character(*), intent(in) :: probe_name

    obj_io_probes(p)%warning_message = 'none'
    obj_io_probes(p)%error_message   = 'none'
    obj_io_probes(p)%description     = 'none'

    ! Nomi delle variabili
    call reg%add( trim(probe_name), 'variables', obj_io_probes(p)%varnames, 'none', 'Probe variables to write', '', .false. )

    ! Frequenza
    call reg%add( trim(probe_name), 'dtime', obj_io_probes(p)%dtime, '1e30', 'Probe output time frequency', '> 0', .false. )
    call reg%add( trim(probe_name), 'diter', obj_io_probes(p)%diter, '1000000000', 'Probe output iter frequency', '> 0', .false. )

    ! Posizione — indici 2D (b,i,j) e coordinate (x,y); nessun indice k né coordinata z.
    call reg%add( trim(probe_name), 'index-position', obj_io_probes(p)%iloc, '0 0 0', 'Probe location by index (b i j)', '>= 0', .false. )
    call reg%add( trim(probe_name), 'position', obj_io_probes(p)%loc, '0.0 0.0', 'Probe location by 2D coordinates', '', .false. )

  end subroutine Register_One_Probe


  subroutine Register_LineProbes(n, lineprobes_names)
    ! Line-probe: campionamento lungo una linea a indice i o j fissato.
    ! Gestione via registry (stessa struttura di Register_Probes): la
    ! configurazione di ogni line-probe è legata a obj_io_lineprobes(p) e
    ! Setup_LineProbes usa i valori caricati — nessun accesso diretto fini%get.
    use Q2D_Parameters_m, only: codename
    use IR_Precision,     only: str
    implicit none
    integer, intent(in) :: n
    character(*), intent(in) :: lineprobes_names(n)
    integer :: p

    if (allocated(obj_io_lineprobes)) deallocate(obj_io_lineprobes)
    allocate(obj_io_lineprobes(n))
    do p = 1, n
      call reg%add( trim(codename)//'-LineProbes', 'lineprobe'//trim(str(.true.,p)), &
                    obj_io_lineprobes(p)%name, 'lineprobe-placeholder', &
                    'Line-probe section name', '', .false. )
      call Register_One_LineProbe(p, lineprobes_names(p))
    end do

  end subroutine Register_LineProbes


  subroutine Register_One_LineProbe(p, lineprobe_name)
    implicit none
    integer, intent(in) :: p
    character(*), intent(in) :: lineprobe_name

    obj_io_lineprobes(p)%warning_message = 'none'
    obj_io_lineprobes(p)%error_message   = 'none'
    obj_io_lineprobes(p)%description     = 'none'

    call reg%add( trim(lineprobe_name), 'variables', obj_io_lineprobes(p)%varnames, 'none', &
                  'Variables to record along the line (whitespace-separated)', '', .false. )

    ! Frequenza
    call reg%add( trim(lineprobe_name), 'dtime', obj_io_lineprobes(p)%dtime, '1e30', &
                  'Line-probe time frequency', '> 0', .false. )
    call reg%add( trim(lineprobe_name), 'diter', obj_io_lineprobes(p)%diter, '1000000000', &
                  'Line-probe iter frequency', '> 0', .false. )

    ! Posizione — `[block, direction]` a 2 componenti (direction=1 linea a i
    ! fissato, direction=2 a j fissato) + `fixed-index` scalare. Entrambi
    ! OBBLIGATORI al momento del solve: i default sentinella `0` permettono a
    ! Setup_LineProbes di emettere `error stop` se l'utente li omette.
    call reg%add( trim(lineprobe_name), 'index-position', obj_io_lineprobes(p)%iloc, '0 0', &
                  'Line-probe location [block, direction]: direction=1 (i-fixed) | 2 (j-fixed)', &
                  '>= 0', .false. )
    call reg%add( trim(lineprobe_name), 'fixed-index', obj_io_lineprobes(p)%fixed_index, '0', &
                  'Fixed-axis sample index along the line', '>= 0', .false. )

  end subroutine Register_One_LineProbe

end module Q2D_Read_IO
