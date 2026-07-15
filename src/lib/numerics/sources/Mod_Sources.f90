module Q2D_Mod_Sources
  ! Punto di ingresso del sistema dei termini sorgente: dispatcher polimorfico su
  ! container_sources(:) % source_realized.
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m, only: Q2D_domain_type
  use Q2D_Lib_Sources,      only: SourceArray, AxisymmetricSource, AreaVariationSource

  implicit none
  private

  type(SourceArray), allocatable :: container_sources(:)
  public :: Initialize_Sources, Compute_Sources

contains

  subroutine Initialize_Sources()
    ! Popola container_sources(:) leggendo obj_source%models(:).
    ! Mutua esclusione `axisymmetric ⊕ area_variation` imposta con error stop
    ! (difesa in profondità — verificata anche in Register_Sources).
    !
    ! Impostare SOURCE_TERM=.true. abilita la chiamata a Compute_Sources nel ciclo
    ! principale (gate `if (SOURCE_TERM)` in Mod_Explicit::Explicit_Step).
    !
    ! Parser della stringa attiva: interpreta `obj_source%active_str` come lista
    ! separata da spazi → `obj_source%models(:)` e imposta i flag booleani derivati
    ! `obj_source%axisymmetric` / `area_variation`.
    use Q2D_Global_m,       only: SOURCE_TERM
    use Q2D_Config_Types_m, only: obj_source
    use Q2D_Parameters_m,   only: llen
    use strings,            only: parse
    implicit none
    integer         :: i
    character(llen) :: args(10)
    character(len=llen) :: clean_str

    ! Interpreta active_str in models(:) (token separati da spazi).
    args = ''
    clean_str = adjustl(obj_source%active_str)
    if (len_trim(clean_str) == 0 .or. trim(clean_str) == 'none') then
      obj_source%n_models = 0
      return
    end if
    call parse(clean_str, ' ', args)
    obj_source%n_models = count(args /= '')
    if (allocated(obj_source%models)) deallocate(obj_source%models)
    allocate(obj_source%models(obj_source%n_models))
    do i = 1, obj_source%n_models
      obj_source%models(i) = trim(adjustl(args(i)))
    end do

    ! Imposta i flag booleani derivati consumati a valle (Mod_Allocate_Data,
    ! Mod_Metrics, Wrap_Setup): ciascun flag è vero se il token corrispondente
    ! compare in models(:).
    obj_source%axisymmetric   = any(index(adjustl(obj_source%models), 'axisymmetric')   > 0)
    obj_source%area_variation = any(index(adjustl(obj_source%models), 'area_variation') > 0)

    if (obj_source%axisymmetric .and. obj_source%area_variation) then
      error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
    end if

    if (allocated(container_sources)) deallocate(container_sources)
    allocate(container_sources(obj_source%n_models))

    do i = 1, obj_source%n_models
      select case (trim(obj_source%models(i)))
        case ('axisymmetric')
          allocate(AxisymmetricSource :: container_sources(i)%source_realized)
        case ('area_variation')
          allocate(AreaVariationSource :: container_sources(i)%source_realized)
        case default
          write(*,*) 'Warning: source type "', trim(obj_source%models(i)), '" not recognized.'
          stop
      end select
    end do
    SOURCE_TERM = .true.
  end subroutine Initialize_Sources


  subroutine Compute_Sources(domain)
    ! Dispatcher polimorfico: invoca compute() su ogni container_sources(i).
    ! Chiamato da Q2D_Mod_Explicit::Explicit_Step nel ciclo RK principale
    ! (gate `if (SOURCE_TERM)`).
    use Q2D_Config_Types_m, only: obj_source
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: i

    do i = 1, obj_source%n_models
      call container_sources(i)%source_realized%compute(domain)
    end do
  end subroutine Compute_Sources

end module Q2D_Mod_Sources
