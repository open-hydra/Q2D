module Q2D_Read_Sources
  ! Punto di ingresso della configurazione dei termini sorgente geometrici.

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

    ! Sorgente geometrica attiva (mutuamente esclusive — Mod_Sources e Mod_RANS
    ! impongono l'esclusione in modo esplicito con `error stop`):
    !   axisymmetric    → Lib_Source_Axisymmetry  (termini 1/r)
    !   area_variation  → Lib_Source_Areavariation (G = ∇h, V·G)
    !
    ! Una lista separata da spazi viene analizzata in `models(:)` per estensioni
    ! future; la mutua esclusione è verificata in fase di Setup. La lista dei
    ! valori ammessi nel registry è scalare.
    call reg%add( section, 'active', obj_source%active_str, 'none', &
                  'Active geometric source model (mutually exclusive)', &
                  'none, axisymmetric, area_variation', .false. )

  end subroutine Register_Sources

end module Q2D_Read_Sources
