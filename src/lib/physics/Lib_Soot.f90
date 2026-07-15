!>@brief: Kernel dei tassi di reazione del soot (LL91, LIN).
!
! 2 routine pubbliche: LL91_soot_rates, LIN_soot_rates.
! Kernel scalari puntuali (indipendenti dalla dimensionalità) che operano su
! (rho, T, C_C2H2, C_O2, Y_soot, Nvol).
module Q2D_Lib_Soot
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Q2D_Parameters_m

  implicit none

  real(R8), parameter :: Mc   = 12.011_R8     ! kg/kmol — massa molare del carbonio del soot
  real(R8), parameter :: rhos = 1800.0_R8     ! kg/m^3  — densità del soot

contains

  subroutine LL91_soot_rates ( rho, T, C_C2H2, C_O2, Y_soot, Nvol, R1, R2, R3, R4 )
    ! Kernel di formazione/ossidazione del soot di Leung-Lindstedt 1991
    ! (nucleazione + crescita superficiale + ossidazione + coagulazione).
    implicit none
    real(R8), intent(in)  :: rho, T, C_C2H2, C_O2, Y_soot, Nvol
    real(R8), intent(out) :: R1, R2, R3, R4

    ! Costanti
    real(R8), parameter :: pi   = 3.141592653589793_R8
    real(R8), parameter :: NA   = 6.02214076e26_R8
    real(R8), parameter :: kB   = 1.380649e-23_R8
    integer,  parameter :: Cmin = 100
    real(R8), parameter :: Ca   = 9.0_R8

    ! Variabili locali
    real(R8) :: C_mol, d_p, S_area, f_Area
    real(R8) :: k1, k2, k3, k4

    ! Concentrazione molare del carbonio del soot
    C_mol = rho*max(Y_soot,0.0_R8)/Mc  ! kmol/m^3

    ! Diametro medio delle particelle dalla densità numerica per volume
    d_p = (6.0_R8*Mc*C_mol / (pi*rhos*Nvol))**(1.0_R8/3.0_R8)

    ! Area superficiale per unità di volume
    S_area = pi * d_p**2 * Nvol
    f_Area = sqrt(pi*((6.0_R8*Mc)/(pi*rhos))**(2.0_R8/3.0_R8)) &
            * C_mol**(1.0_R8/3.0_R8) * Nvol**(1.0_R8/6.0_R8)

    ! --- Costanti cinetiche (LL91) ---
    k1 = 0.1e5_R8 * exp(-21100.0_R8/T)
    k2 = 0.6e4_R8 * exp(-12100.0_R8/T)
    k3 = 0.1e5_R8 * sqrt(T) * exp(-19680.0_R8/T)
    k4 = 2.0_R8 * Ca * ((6.0_R8*Mc)/(pi*rhos))**(1.0_R8/6.0_R8) &
                     * sqrt((6.0_R8*kB*T)/rhos)

    ! --- Tassi di reazione ---
    R1 = k1 * C_C2H2                   ! nucleazione [kmol/m^3/s]
    R2 = k2 * C_C2H2 * f_Area          ! crescita superficiale [kmol/m^3/s]
    R3 = k3 * S_area * C_O2            ! ossidazione [kmol/m^3/s]
    R4 = (2.0_R8 * NA / Cmin) * R1 &   ! nucleazione-coagulazione [#/m^3/s]
         - k4 * C_mol**(1.0_R8/6.0_R8) * Nvol**(11.0_R8/6.0_R8)
  end subroutine LL91_soot_rates


  subroutine LIN_soot_rates ( rho, T, C_C2H2, C_O2, Y_soot, Nvol, R1, R2, R3, R4 )
    ! Nucleazione Lindstedt 2005 + crescita superficiale Lindstedt 1994 +
    ! ossidazione/coagulazione LL91.
    implicit none
    real(R8), intent(in)  :: rho, T, C_C2H2, C_O2, Y_soot, Nvol
    real(R8), intent(out) :: R1, R2, R3, R4

    ! Costanti
    real(R8), parameter :: pi   = 3.141592653589793_R8
    real(R8), parameter :: NA   = 6.02214076e26_R8
    real(R8), parameter :: kB   = 1.380649e-23_R8
    integer,  parameter :: Cmin = 100
    real(R8), parameter :: Ca   = 9.0_R8

    ! Variabili locali
    real(R8) :: C_mol, d_p, S_area
    real(R8) :: k1, k2, k3, k4

    ! Concentrazione molare del carbonio del soot
    C_mol = rho * max(Y_soot, 0.0_R8) / Mc   ! kmol/m^3

    ! Diametro medio delle particelle dalla densità numerica
    d_p = (6.0_R8 * Mc * C_mol / (pi * rhos * Nvol))**(1.0_R8/3.0_R8)

    ! Area superficiale per unità di volume
    S_area = pi * d_p**2 * Nvol

    ! --- Costanti cinetiche ---
    ! Nucleazione Lindstedt (2005)
    k1 = 0.63e4_R8 * exp(-21100.0_R8 / T)
    ! Crescita superficiale Lindstedt (1994) (Bockhorn eq. 27.36)
    k2 = 0.1e-11_R8 * exp(-12100.0_R8 / T)
    ! Ossidazione LL91
    k3 = 0.1e5_R8 * sqrt(T) * exp(-19680.0_R8 / T)
    ! Coagulazione LL91
    k4 = 2.0_R8 * Ca * ((6.0_R8*Mc)/(pi*rhos))**(1.0_R8/6.0_R8) &
                     * sqrt((6.0_R8*kB*T)/rhos)

    ! --- Tassi di reazione ---
    R1 = k1 * C_C2H2                   ! nucleazione [kmol/m^3/s]
    R2 = k2 * C_C2H2 * Nvol            ! crescita superficiale [kmol/m^3/s]
    R3 = k3 * S_area * C_O2            ! ossidazione
    R4 = (2.0_R8 * NA / Cmin) * R1 &   ! nucleazione-coagulazione [#/m^3/s]
         - k4 * C_mol**(1.0_R8/6.0_R8) * Nvol**(11.0_R8/6.0_R8)
  end subroutine LIN_soot_rates

end module Q2D_Lib_Soot
