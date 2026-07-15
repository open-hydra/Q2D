module Q2D_Lib_Sources
  ! Pattern a oggetti per i termini sorgente:
  !   type, abstract :: SourceTerm + class(SourceTerm), allocatable :: source_realized
  !   2 estensioni di tipo: AxisymmetricSource, AreaVariationSource
  !   container SourceArray per lo storage polimorfico in Q2D_Mod_Sources.
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m, only: Q2D_domain_type
  implicit none
  private

  public :: SourceTerm
  public :: AxisymmetricSource
  public :: AreaVariationSource
  public :: SourceArray

  ! Classe base astratta
  type, abstract :: SourceTerm
  contains
    procedure(source_if), deferred, pass :: compute
  end type SourceTerm

  abstract interface
    subroutine source_if(this, domain)
      import :: Q2D_domain_type, SourceTerm
      class(SourceTerm),     intent(inout) :: this
      type(Q2D_domain_type), intent(inout) :: domain
    end subroutine source_if
  end interface

  ! Estensioni di tipo
  type, extends(SourceTerm) :: AxisymmetricSource
  contains
    procedure, pass(this) :: compute => source_axisymmetric
  end type AxisymmetricSource

  type, extends(SourceTerm) :: AreaVariationSource
  contains
    procedure, pass(this) :: compute => source_area_variation
  end type AreaVariationSource

  ! Container polimorfico (wrapper per un array di class allocatable in Mod_Sources)
  type :: SourceArray
    class(SourceTerm), allocatable :: source_realized
  end type SourceArray

contains

  subroutine source_axisymmetric ( this, domain )
    use Q2D_Global_m,                only: model
    use Q2D_Lib_Source_Axisymmetry,  only: Source_Axi_Eul_Blk, Source_Axi_Visc_Blk
    use Q2D_Mod_MPI,                 only: is_local_block
    implicit none
    class(AxisymmetricSource), intent(inout) :: this
    type(Q2D_domain_type),     intent(inout) :: domain
    ! Variabili locali
    integer :: b

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle

      call Source_Axi_Eul_Blk ( domain%blk(b)%P,            &
                                domain%blk(b)%R,            &
                                domain%blk(b)%axis_r_ave,   &
                                domain%blk(b)%axis_inv_r,   &
                                domain%blk(b)%sur,          &
                                domain%blk(b)%dim           )

      if (model>0) then
        call Source_Axi_Visc_Blk ( domain%blk(b)%P,            &
                                   domain%blk(b)%R,            &
                                   domain%blk(b)%M,            &
                                   domain%blk(b)%yn,           &
                                   domain%blk(b)%axis_r_ave,   &
                                   domain%blk(b)%axis_inv_r,   &
                                   domain%blk(b)%sur,          &
                                   domain%blk(b)%dim           )
      end if

    end do
    associate(dummy => this); end associate
  end subroutine source_axisymmetric


  subroutine source_area_variation ( this, domain )
    use Q2D_Global_m,                 only: model
    use Q2D_Lib_Source_Areavariation, only: Source_Area_Eul_Blk, Source_Area_Visc_Blk
    use Q2D_Mod_MPI,                  only: is_local_block
    implicit none
    class(AreaVariationSource), intent(inout) :: this
    type(Q2D_domain_type),      intent(inout) :: domain
    ! Variabili locali
    integer :: b

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      if (.not. domain%blk(b)%area_variation_active) cycle

      call Source_Area_Eul_Blk ( domain%blk(b)%P,                &
                                 domain%blk(b)%R,                &
                                 domain%blk(b)%sur,              &
                                 domain%blk(b)%dim,              &
                                 domain%blk(b)%area_geom_term_x, &
                                 domain%blk(b)%area_geom_term_y  )

      if (model>0) then
        call Source_Area_Visc_Blk ( domain%blk(b)%P,                &
                                    domain%blk(b)%R,                &
                                    domain%blk(b)%M,                &
                                    domain%blk(b)%yn,               &
                                    domain%blk(b)%sur,              &
                                    domain%blk(b)%dim,              &
                                    domain%blk(b)%area_geom_term_x, &
                                    domain%blk(b)%area_geom_term_y, &
                                    domain%blk(b)%node              )
      end if

    end do
    associate(dummy => this); end associate
  end subroutine source_area_variation

end module Q2D_Lib_Sources
