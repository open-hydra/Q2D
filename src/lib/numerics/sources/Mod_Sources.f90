module Q2D_Mod_Sources
  ! Q2D-extra: MOSE has no numerics/sources/.
  ! Q2D-only entry point of the sources system: polymorphic dispatcher over
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
    ! Populates container_sources(:) reading obj_source%models(:).
    ! Mutual exclusion `axisymmetric ⊕ area_variation` enforced with error stop
    ! (defense in depth — also verified in Register_Sources).
    !
    ! Setting SOURCE_TERM=.true. enables the Compute_Sources call in the main loop
    ! (Mod_Explicit::Explicit_Step gate `if (SOURCE_TERM)`).
    !
    ! Active-string parser: parses `obj_source%active_str` whitespace-separated
    ! list → `obj_source%models(:)` + sets derived boolean flags
    ! `obj_source%axisymmetric` / `area_variation`.
    use Q2D_Global_m,       only: SOURCE_TERM
    use Q2D_Config_Types_m, only: obj_source
    use Q2D_Parameters_m,   only: llen
    use strings,            only: parse
    implicit none
    integer         :: i
    character(llen) :: args(10)
    character(len=llen) :: clean_str

    ! Parse active_str into models(:) (whitespace-separated tokens).
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

    ! Set derived boolean flags consumed downstream (Mod_Allocate_Data,
    ! Mod_Metrics, Wrap_Setup):
    ! `axis_source_active = any(index(source_names, 'axisymmetric') > 0)`.
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
    ! Polymorphic dispatcher: invokes compute() on each container_sources(i).
    ! Called by Q2D_Mod_Explicit::Explicit_Step in the main RK loop
    ! (gate `if (SOURCE_TERM)`). Q2D-extra — MOSE has no such call.
    use Q2D_Config_Types_m, only: obj_source
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: i

    do i = 1, obj_source%n_models
      call container_sources(i)%source_realized%compute(domain)
    end do
  end subroutine Compute_Sources

end module Q2D_Mod_Sources
