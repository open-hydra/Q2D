!>@brief: Soot reaction-rate kernels (LL91, LIN). 2D-reduced.
!
! 2 public routines: LL91_soot_rates, LIN_soot_rates.
! Scalar kernels (no 2D/3D dependence) working on
! (rho, T, C_C2H2, C_O2, Y_soot, Nvol).
!
! A HACA_soot_rates dispatcher (commented out in the MOSE peer) is NOT included —
! out of scope (no innovation).
module Q2D_Lib_Soot
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Parameters_m

  implicit none

  real(R8), parameter :: Mc   = 12.011_R8     ! kg/kmol — soot carbon molar mass
  real(R8), parameter :: rhos = 1800.0_R8     ! kg/m^3  — soot density

contains

  subroutine LL91_soot_rates ( rho, T, C_C2H2, C_O2, Y_soot, Nvol, R1, R2, R3, R4 )
    ! Leung-Lindstedt 1991 soot formation/oxidation kernel
    ! (nucleation + growth + oxidation + coagulation).
    implicit none
    real(R8), intent(in)  :: rho, T, C_C2H2, C_O2, Y_soot, Nvol
    real(R8), intent(out) :: R1, R2, R3, R4

    ! Constants
    real(R8), parameter :: pi   = 3.141592653589793_R8
    real(R8), parameter :: NA   = 6.02214076e26_R8
    real(R8), parameter :: kB   = 1.380649e-23_R8
    integer,  parameter :: Cmin = 100
    real(R8), parameter :: Ca   = 9.0_R8

    ! Locals
    real(R8) :: C_mol, d_p, S_area, f_Area
    real(R8) :: k1, k2, k3, k4

    ! Molar concentration of soot carbon
    C_mol = rho*max(Y_soot,0.0_R8)/Mc  ! kmol/m^3

    ! Mean particle diameter from volume-based number density
    d_p = (6.0_R8*Mc*C_mol / (pi*rhos*Nvol))**(1.0_R8/3.0_R8)

    ! Surface area per unit volume
    S_area = pi * d_p**2 * Nvol
    f_Area = sqrt(pi*((6.0_R8*Mc)/(pi*rhos))**(2.0_R8/3.0_R8)) &
            * C_mol**(1.0_R8/3.0_R8) * Nvol**(1.0_R8/6.0_R8)

    ! --- Rate constants (LL91) ---
    k1 = 0.1e5_R8 * exp(-21100.0_R8/T)
    k2 = 0.6e4_R8 * exp(-12100.0_R8/T)
    k3 = 0.1e5_R8 * sqrt(T) * exp(-19680.0_R8/T)
    k4 = 2.0_R8 * Ca * ((6.0_R8*Mc)/(pi*rhos))**(1.0_R8/6.0_R8) &
                     * sqrt((6.0_R8*kB*T)/rhos)

    ! --- Reaction rates ---
    R1 = k1 * C_C2H2                   ! nucleation [kmol/m^3/s]
    R2 = k2 * C_C2H2 * f_Area          ! surface growth [kmol/m^3/s]
    R3 = k3 * S_area * C_O2            ! oxidation [kmol/m^3/s]
    R4 = (2.0_R8 * NA / Cmin) * R1 &   ! Nucleation-Coagulation [#/m^3/s]
         - k4 * C_mol**(1.0_R8/6.0_R8) * Nvol**(11.0_R8/6.0_R8)
  end subroutine LL91_soot_rates


  subroutine LIN_soot_rates ( rho, T, C_C2H2, C_O2, Y_soot, Nvol, R1, R2, R3, R4 )
    ! Lindstedt 2005 nucleation + Lindstedt 1994 surface growth + LL91
    ! oxidation/coagulation.
    implicit none
    real(R8), intent(in)  :: rho, T, C_C2H2, C_O2, Y_soot, Nvol
    real(R8), intent(out) :: R1, R2, R3, R4

    ! Constants
    real(R8), parameter :: pi   = 3.141592653589793_R8
    real(R8), parameter :: NA   = 6.02214076e26_R8
    real(R8), parameter :: kB   = 1.380649e-23_R8
    integer,  parameter :: Cmin = 100
    real(R8), parameter :: Ca   = 9.0_R8

    ! Locals
    real(R8) :: C_mol, d_p, S_area
    real(R8) :: k1, k2, k3, k4

    ! Molar concentration of soot carbon
    C_mol = rho * max(Y_soot, 0.0_R8) / Mc   ! kmol/m^3

    ! Mean particle diameter from number density
    d_p = (6.0_R8 * Mc * C_mol / (pi * rhos * Nvol))**(1.0_R8/3.0_R8)

    ! Surface area per unit volume
    S_area = pi * d_p**2 * Nvol

    ! --- Rate constants ---
    ! Lindstedt (2005) nucleation
    k1 = 0.63e4_R8 * exp(-21100.0_R8 / T)
    ! Lindstedt (1994) surface growth (Bockhorn eq. 27.36)
    k2 = 0.1e-11_R8 * exp(-12100.0_R8 / T)
    ! LL91 oxidation
    k3 = 0.1e5_R8 * sqrt(T) * exp(-19680.0_R8 / T)
    ! LL91 coagulation
    k4 = 2.0_R8 * Ca * ((6.0_R8*Mc)/(pi*rhos))**(1.0_R8/6.0_R8) &
                     * sqrt((6.0_R8*kB*T)/rhos)

    ! --- Reaction rates ---
    R1 = k1 * C_C2H2                   ! nucleation [kmol/m^3/s]
    R2 = k2 * C_C2H2 * Nvol            ! surface growth [kmol/m^3/s]
    R3 = k3 * S_area * C_O2            ! oxidation
    R4 = (2.0_R8 * NA / Cmin) * R1 &   ! Nucleation-Coagulation [#/m^3/s]
         - k4 * C_mol**(1.0_R8/6.0_R8) * Nvol**(11.0_R8/6.0_R8)
  end subroutine LIN_soot_rates

end module Q2D_Lib_Soot
