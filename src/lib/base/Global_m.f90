module Q2D_Global_m
  use Q2D_Parameters_m, only: clen

  implicit none
  private

  public :: Q2D_phase_prefix
  public :: nsc, nu, nv, np, nt, nrans
  public :: nc, nsoot, npass, nps
  public :: nprim, nres, ndir, gc, model
  public :: SOURCE_TERM   ! Gate enabling Compute_Sources in Mod_Explicit

  character(len=clen) :: Q2D_phase_prefix = ''

  integer :: nsc           ! Number of species
  integer :: nu, nv        ! Indices of the velocity/momentum components in the primitive/residual vector (no nw in 2D)
  integer :: np            ! Index of pressure/energy in the primitive/residual vector
  integer :: nt, nrans     ! Index of the first RANS variable in the primitive/residual vector and number of RANS variables
  integer :: nc, nsoot     ! Index of the first soot variable, number of soot variables
  integer :: npass = 0     ! Number of passive scalars
  integer :: nps           ! Index of the first passive scalar in the primitive/residual vector
  integer :: nprim         ! Number of primitive variables
  integer :: nres          ! Number of variables whose residuals are printed
  integer :: ndir = 2      ! Number of dimensions of the computational domain (2D formulation)
  integer :: gc = 2        ! Number of ghost cells
  integer :: model         ! Model index (e.g. 0: Euler, 1: laminar NS, 2: RANS)

  logical :: SOURCE_TERM = .false.   ! Set to .true. in Initialize_Sources

end module Q2D_Global_m
