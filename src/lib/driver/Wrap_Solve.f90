module Q2D_Wrap_Solve
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Q2D_solve

contains

  subroutine Q2D_solve ( simulation, External_Function )
    use Q2D_Advanced_Types_m, only: Q2D_simulation_type
    use Q2D_Config_Types_m,   only: obj_time_scheme, obj_multigrid
    use Q2D_Mod_Multigrid,    only: Prolongation
    use Q2D_Mod_Explicit,     only: Explicit_Step
    implicit none
    type(Q2D_simulation_type), intent(inout) :: simulation
    external :: External_Function

    select case (trim(obj_time_scheme%solver_type))

      case ('euler', 'RK2', 'RK3')

        ! 1. Calculate solution on local grid
        call Explicit_Step( simulation%domain, External_Function )

        ! 2. Interpolate solution on finer grid (multigrid V-cycle)
        if (obj_multigrid%change_MG) then
          call Prolongation ( Fine=simulation%domain(obj_multigrid%MG_level-1), &
                              Coarse=simulation%domain(obj_multigrid%MG_level) )
        endif

    end select

  end subroutine Q2D_solve

end module Q2D_Wrap_Solve
