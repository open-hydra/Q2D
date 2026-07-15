module Q2D_Base_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: Q2D_tensor_2D_type
  public :: Q2D_tensor_3D_type
  public :: Q2D_vector_2D_type
  public :: Q2D_tensor_2D_R2_type
  public :: Q2D_f_metrics_type
  public :: Q2D_d_metrics_type

  type :: Q2D_tensor_2D_type
    real(R8) :: c(2,2)                        !> Metric tensor components.
  end type Q2D_tensor_2D_type

  type :: Q2D_tensor_3D_type
    real(R8) :: c(3,3)                        !> 3x3 stress/strain tensor (S_33 = -V·G from Definition A).
  end type Q2D_tensor_3D_type

  type :: Q2D_vector_2D_type
    real(R8) :: c(2)                          !> Mean cell-length components.
  end type Q2D_vector_2D_type

  type :: Q2D_tensor_2D_R2_type
    real(R8) :: c(2,2,2)                      !> Metric tensor components.
  end type Q2D_tensor_2D_R2_type

  type :: Q2D_f_metrics_type
    real(R8) :: N(2)                          !> Unit normal vector.
    real(R8) :: int_l                         !> Interface length.
  end type Q2D_f_metrics_type

  type :: Q2D_d_metrics_type
    type(Q2D_f_metrics_type), allocatable :: f(:,:)
  end type Q2D_d_metrics_type

end module Q2D_Base_Types_m
