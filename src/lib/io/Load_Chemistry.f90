! Wrapper Q2D per il caricamento dei dati di chimica via FLINT.
! Strato sottile attorno a FLINT_Load_Chemistry::read_chemistry che riporta
! l'esito tramite obj_chemistry%description / %error_message.
! La logica del solutore ODE / splitting di Strang NON è qui: risiede in
! Q2D_Assign_Setup (caricamento dati e setup ode/strang sono separati).
module Q2D_Load_Chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_Chemistry

contains

  subroutine Load_Chemistry()
    use Q2D_Config_Types_m,  only: obj_chemistry
    use FLINT_Load_Chemistry, only: read_chemistry
    implicit none
    integer :: error

    error = read_chemistry(mech_name = obj_chemistry%mechanism_name)
    select case (error)
      case (0)
        obj_chemistry%description = 'See next section for details'
      case (1)
        obj_chemistry%description = 'Frozen mixture'
      case (2)
        obj_chemistry%error_message = '[ERROR] Chemistry info file found but could not be read'
      case (3)
        obj_chemistry%error_message = '[ERROR] Chemistry rate table file not found'
      case (4)
        obj_chemistry%error_message = '[ERROR] Chemistry rate table file found but could not be read'
      case default
        obj_chemistry%error_message = '[ERROR] Unknown error loading chemistry'
    end select

  end subroutine Load_Chemistry

end module Q2D_Load_Chemistry
