module Q2D_Parameters_m

  implicit none
  private

  public :: hlen, llen, clen, codename, guide

  integer, parameter             :: hlen = 1024
  integer, parameter             :: llen = 256
  integer, parameter             :: clen = 16
  character(len=clen), parameter :: codename = 'Q2D'          ! `parameter` constant (immutable): the codename does not change at runtime.

  integer, dimension(6,3), parameter ::  guide  = reshape( [ 1, 0, 0, &
                                                            -1, 0, 0, &
                                                             0, 1, 0, &
                                                             0,-1, 0, &
                                                             0, 0, 1, &
                                                             0, 0,-1  ] , shape(guide), order=[2,1] )

end module Q2D_Parameters_m
