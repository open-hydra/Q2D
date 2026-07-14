module Q2D_Global_m
  use Q2D_Parameters_m, only: clen

  implicit none
  private

  public :: Q2D_phase_prefix
  public :: nsc, nu, nv, np, nt, nrans
  public :: nc, nsoot, npass, nps
  public :: nprim, nres, ndir, gc, model
  public :: SOURCE_TERM   ! Q2D-extra — gate Compute_Sources in Mod_Explicit

  character(len=clen) :: Q2D_phase_prefix = ''

  integer :: nsc           ! Number of species
  integer :: nu, nv        ! Velocity/momentum components indexes in prim/residuals vector (no nw)
  integer :: np            ! Pressure/energy index in prim/residuals vector
  integer :: nt, nrans     ! First RANS variable index in prim/residuals vector and number of RANS variables
  integer :: nc, nsoot     ! Index of first soot variable, number of soot variables
  integer :: npass = 0     ! Number of passive scalar variables
  integer :: nps           ! First passive scalar index in prim/residuals vector
  integer :: nprim         ! Number of primitive variables
  integer :: nres          ! Number of variables whose residuals are printed
  integer :: ndir = 2      ! Number of dimensions of computational frame (2D reduction)
  integer :: gc = 2        ! Ghost cells
  integer :: model         ! Model index (e.g., 0: Euler, 1: laminar NS, 2: RANS)

  logical :: SOURCE_TERM = .false.   ! Set true in Initialize_Sources

end module Q2D_Global_m
