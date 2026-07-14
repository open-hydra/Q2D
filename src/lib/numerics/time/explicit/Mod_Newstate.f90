module Q2D_Mod_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Assign_Integration_Variables

  !> Concrete procedure pointing to one of the function realizations
  procedure(newstate_if), pointer, public :: RK_Newstate

  !> Abstract interface relative to the "Newstate" procedure
  abstract interface
    subroutine Newstate_if ( domain, irk )
      use Q2D_Advanced_Types_m, only: Q2D_domain_type
      implicit none
      type(Q2D_domain_type), intent(inout) :: domain
      integer, intent(in)                   :: irk
    end subroutine Newstate_if
  end interface

contains

  subroutine Assign_Integration_Variables ()
    use Q2D_Config_Types_m, only: obj_time_scheme
    use Q2D_Lib_Newstate,   only: Newstate_Primitive, Newstate_Conservative
    implicit none

    select case (obj_time_scheme%integration_variables)
    case ('prim')
      RK_Newstate => Newstate_Primitive
      obj_time_scheme%integration_variables = 'Primitive'
    case ('cons')
      RK_Newstate => Newstate_Conservative
      obj_time_scheme%integration_variables = 'Conservative'
    end select

  end subroutine Assign_Integration_Variables

end module Q2D_Mod_Newstate
