module Q2D_Procedures_m
  use Q2D_Wrap_Setup,       only: Q2D_setup
  use Q2D_Wrap_Solve,       only: Q2D_solve
  use Q2D_Wrap_Postprocess, only: Q2D_postprocess

  implicit none
  private
  public :: Q2D_type

  type :: Q2D_type
  contains
    procedure, nopass :: setup       => Q2D_setup
    procedure, nopass :: solve       => Q2D_solve
    procedure, nopass :: postprocess => Q2D_postprocess
  end type Q2D_type

end module Q2D_Procedures_m
