module Q2D_Lib_Shock_Detector
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: SD_Tramel

  real(R8), parameter, private :: DELTA = 20.d0

contains

  ! Jameson-type pressure-based shock detector (Tramel):
  ! 2D formulation with stencil pmatrix(0:2,0:2) along the two grid directions.
  pure function SD_Tramel(pmatrix) result(beta)
    implicit none
    real(R8), intent(in) :: pmatrix(0:2,0:2)
    real(R8)             :: beta
    real(R8) :: kp1, kp2, sw, phi
    real(R8) :: p1, p2, p3, p4, p

    p  = pmatrix(1,1)
    p1 = pmatrix(0,1)
    p2 = pmatrix(2,1)
    p3 = pmatrix(1,0)
    p4 = pmatrix(1,2)

    kp1 = abs( (p2 - 2d0*p + p1) / (p2 + 2d0*p + p1) )
    kp2 = abs( (p4 - 2d0*p + p3) / (p4 + 2d0*p + p3) )
    sw = max(kp1, kp2)

    if (sw < 1d0/DELTA) then
      phi = max(sw/DELTA, 0d0)
      beta = 1d0 - tanh(10d0*phi*phi*phi)
    else
      beta = 0d0
    end if

  end function SD_Tramel

end module Q2D_Lib_Shock_Detector
