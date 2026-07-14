module Q2D_Read_Sources
  ! Q2D-specific module (no analogue in the MOSE peer) — entry point of the config sources system.
  ! Absorbs the Q2D legacy IO_Sources content.

  use iso_fortran_env,    only: I4 => int32, R8 => real64
  use Q2D_Parameters_m,   only: llen, codename
  use Q2D_Config_Types_m, only: obj_source
  use Q2D_Registry_m,     only: reg

  implicit none
  private
  public :: Register_Sources

contains

  subroutine Register_Sources()
    implicit none
    character(len=:), allocatable :: section

    section = trim(codename)//'-Sources'

    obj_source%warning_message = 'none'
    obj_source%error_message   = 'none'
    obj_source%description     = 'none'

    ! Active geometric source (mutually exclusive — Mod_Sources and Mod_RANS
    ! enforce the exclusion explicitly with `error stop`):
    !   axisymmetric    → Lib_Source_Axisymmetry  (1/r terms)
    !   area_variation  → Lib_Source_Areavariation (G = ∇h, V·G)
    !
    ! A whitespace-separated list is parsed into `models(:)` for future
    ! extensibility; the mutual exclusion is verified at Setup time. The registry
    ! allowed-values list is scalar.
    call reg%add( section, 'active', obj_source%active_str, 'none', &
                  'Active geometric source model (mutually exclusive)', &
                  'none, axisymmetric, area_variation', .false. )

  end subroutine Register_Sources

end module Q2D_Read_Sources
