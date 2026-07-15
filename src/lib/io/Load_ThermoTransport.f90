! Wrapper Q2D per il caricamento delle tabelle termodinamiche e di trasporto via FLINT.
! Strato sottile attorno a
! FLINT_Load_ThermoTransport::{read_idealgas_thermo, read_idealgas_transport}
! che riporta gli esiti nei descrittori obj_thermo / obj_transport e imposta
! Q2D_Global_m::nsc dal numero di specie di FLINT.
module Q2D_Load_ThermoTransport
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_ThermoTransport

contains

  subroutine Load_ThermoTransport()
    use Q2D_Config_Types_m,         only: obj_thermo, obj_transport
    use Q2D_Global_m,               only: nsc, Q2D_phase_prefix
    use FLINT_Load_ThermoTransport, only: read_idealgas_thermo, read_idealgas_transport
    use FLINT_Lib_Thermodynamic,    only: FLINT_phase_prefix, ns
    implicit none
    integer :: error

    !! Termodinamica
    obj_thermo%warning_message = 'none'
    obj_thermo%error_message   = 'none'
    obj_thermo%description     = 'none'

    FLINT_phase_prefix = Q2D_phase_prefix
    error = read_idealgas_thermo('INPUT')

    select case (error)
      case (0)
        obj_thermo%description = 'Thermally perfect gas'
      case (5)
        obj_thermo%description = 'Thermally perfect gas'
      case (1)
        obj_thermo%error_message = '[ERROR] Phase file (phase.txt) not found'
      case (2)
        obj_thermo%error_message = '[ERROR] Phase file (phase.txt) found but could not be read'
      case (3)
        obj_thermo%error_message = '[ERROR] Thermo table file not found'
      case (4)
        obj_thermo%error_message = '[ERROR] Thermo table file found but could not be read'
      case (6)
        obj_thermo%error_message = '[ERROR] Composition file found but could not be read'
      case default
        obj_thermo%error_message = '[ERROR] Unknown error loading thermodynamic data'
    end select

    !! Trasporto
    obj_transport%warning_message = 'none'
    obj_transport%error_message   = 'none'
    obj_transport%description     = 'none'

    error = read_idealgas_transport('INPUT')
    select case (error)
      case (0)
        obj_transport%description = 'Mixture-averaged transport properties'
      case (1)
        obj_transport%description = 'Unavailable'
      case (2)
        obj_transport%error_message = '[ERROR] Transport table file found but could not be read'
      case (3)
        obj_transport%error_message = '[ERROR] Transport data temperature range does not match thermo data range'
      case default
        obj_transport%error_message = '[ERROR] Unknown error loading transport data'
    end select

    nsc = ns   ! numero di specie di Q2D_Global_m impostato dal conteggio di FLINT

  end subroutine Load_ThermoTransport

end module Q2D_Load_ThermoTransport
