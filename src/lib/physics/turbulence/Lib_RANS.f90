!>@brief: Interfacce astratte e puntatori a procedura RANS (riduzione 2D Q2D).
!
! 8 puntatori a procedura + 8 interfacce astratte, tra cui:
!   RANS_Diffusive_SourceQ2D / RANS_Diffusive_SourceQ2D_if: sorgente diffusiva
!     per la variazione d'area;
!   Eddy_Viscosity_if a 9 argomenti: include vr/r_ave (assialsimmetrico) e
!     V_dot_G (opzionale, Definizione A);
!   RANS_Diffusive_Flux_if a 12 argomenti: include int_l (analogo 2D dell'area
!     di faccia) + vr/r_ave + mut opzionale.
!
! `RSM_Symmetry` è consumata dal blocco `if(RSM)` di
! Lib_BC_Fluxes_Symmetry::BC_Symmetry_Visc (docs/theory/q2d-rans-closures.md §6.7).
module Q2D_Lib_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: Eddy_Viscosity, RANS_Diffusive_Flux, RANS_Diffusive_SourceQ2D
  public :: RANS_Source_Terms, RANS_Set_Wall_Values, RANS_Set_Blowing_Wall
  public :: RANS_Extrapolate_Wall, RANS_Enforce_Realizability
  public :: RSM_Symmetry

  ! Puntatori a procedura (8)
  procedure(Eddy_Viscosity_if),             pointer :: Eddy_Viscosity
  procedure(RANS_Diffusive_Flux_if),        pointer :: RANS_Diffusive_Flux
  procedure(RANS_Diffusive_SourceQ2D_if),   pointer :: RANS_Diffusive_SourceQ2D => null()
  procedure(RANS_Source_Terms_if),          pointer :: RANS_Source_Terms
  procedure(RANS_Set_Wall_Values_if),       pointer :: RANS_Set_Wall_Values
  procedure(RANS_Set_Blowing_Wall_if),      pointer :: RANS_Set_Blowing_Wall
  procedure(RANS_Extrapolate_Wall_if),      pointer :: RANS_Extrapolate_Wall
  procedure(RANS_Enforce_Realizability_if), pointer :: RANS_Enforce_Realizability

  !> Eddy_Viscosity_if a 9 argomenti.
  !> V_dot_G opzionale: presente -> propagato a Strain_Tensor per la correzione S_33 = -V·G
  !> (modifica |S|^2 = |S|^2_2D + 2(V·G)^2 nel Boussinesq); assente -> Definizione B.
  abstract interface
    subroutine Eddy_Viscosity_if ( mut, rans_variables, mul, rho, vel_gradient, vr, r_ave, walldist, V_dot_G )
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(in), dimension(nRANS) :: rans_variables
      real(R8), intent(in)                   :: mul
      real(R8), intent(in)                   :: rho
      real(R8), intent(in), dimension(2,2)   :: vel_gradient
      real(R8), intent(in)                   :: walldist
      real(R8), intent(in)                   :: vr, r_ave
      real(R8), intent(out)                  :: mut
      real(R8), intent(in), optional         :: V_dot_G
    end subroutine Eddy_Viscosity_if
  end interface

  !> RANS_Diffusive_Flux_if a 12 argomenti.
  !> mut opzionale: se presente, la viscosità turbolenta del chiamante è riusata (evita
  !> il ricalcolo in SST e garantisce consistenza con la mu_t dello sforzo di quantità
  !> di moto in presenza di variazione d'area).
  !> int_l è l'analogo 2D dell'area di faccia 3D.
  abstract interface
    subroutine RANS_Diffusive_Flux_if(flux, rans_variables, mul, vel_gradient, rans_gradient, rho, &
                                      int_l, normal, vr, r_ave, dist, mut)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(in), dimension(nRANS)   :: rans_variables
      real(R8), intent(in)                     :: mul
      real(R8), intent(in), dimension(2,2)     :: vel_gradient
      real(R8), intent(in), dimension(nRANS,2) :: rans_gradient
      real(R8), intent(in)                     :: rho
      real(R8), intent(in)                     :: int_l, normal(2), dist
      real(R8), intent(in)                     :: vr, r_ave
      real(R8), intent(out), dimension(nRANS)  :: flux
      real(R8), intent(in), optional           :: mut
    end subroutine RANS_Diffusive_Flux_if
  end interface

  !> RANS_Diffusive_SourceQ2D_if (variazione d'area).
  !> Restituisce source(nrans,2); il chiamante proietta sulla geometria (G o 1/r).
  abstract interface
    subroutine RANS_Diffusive_SourceQ2D_if(source, rans_variables, mul, mut, vel_gradient, rans_gradient, rho, dist)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(in), dimension(nRANS)    :: rans_variables
      real(R8), intent(in)                      :: mul
      real(R8), intent(in)                      :: mut
      real(R8), intent(in), dimension(2,2)      :: vel_gradient
      real(R8), intent(in), dimension(nRANS,2)  :: rans_gradient
      real(R8), intent(in)                      :: rho
      real(R8), intent(in)                      :: dist
      real(R8), intent(out), dimension(nRANS,2) :: source
    end subroutine RANS_Diffusive_SourceQ2D_if
  end interface

  abstract interface
    subroutine RANS_Source_Terms_if(domain)
      use Q2D_Advanced_Types_m, only: Q2D_domain_type
      implicit none
      type(Q2D_domain_type), intent(inout) :: domain
    end subroutine RANS_Source_Terms_if
  end interface

  abstract interface
    subroutine RANS_Set_Wall_Values_if(mil, rans_variables, dist)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(out), dimension(nRANS)  :: rans_variables
      real(R8), intent(in)                     :: mil
      real(R8), intent(in)                     :: dist
    end subroutine RANS_Set_Wall_Values_if
  end interface

  !> tau(2): vettore sforzo a parete ridotto 2D.
  abstract interface
    subroutine RANS_Set_Blowing_Wall_if(rho, mil, tau, rans_variables, mdot, dist)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(out), dimension(nRANS)  :: rans_variables
      real(R8), intent(in)                     :: rho, mil, mdot, dist
      real(R8), intent(in), dimension(2)       :: tau
    end subroutine RANS_Set_Blowing_Wall_if
  end interface

  abstract interface
    subroutine RANS_Extrapolate_Wall_if ( prim, wall_rans, rho, wall_rho, ghost_rans )
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nRANS
      implicit none
      real(R8), intent(in),  dimension(nRANS) :: prim, wall_rans
      real(R8), intent(in) :: rho, wall_rho
      real(R8), intent(out), dimension(nRANS) :: ghost_rans
    end subroutine RANS_Extrapolate_Wall_if
  end interface

  abstract interface
    subroutine RANS_Enforce_Realizability_if ( var )
      use iso_fortran_env, only: R8 => real64
      implicit none
      real(R8), intent(inout),  dimension(:) :: var
    end subroutine RANS_Enforce_Realizability_if
  end interface

contains

  ! Valori degli sforzi di Reynolds all'interfaccia di simmetria e nella cella ghost
  ! (docs/theory/q2d-rans-closures.md §6.7). Riduzione 2D Q2D (4 sforzi + omega;
  ! R13=R23 eliminati, 7->5): il piano di simmetria è un piano x o y (nessun piano z
  ! in 2D). Gli sforzi normali R11,R22,R33 e lo scalare omega sono PARI (ghost a
  ! gradiente nullo); lo shear nel piano R12 è DISPARI su qualunque piano di
  ! simmetria nel piano (avendo un indice lungo la normale al piano), quindi il suo
  ! valore di faccia è nullo e il ghost cambia segno. È la restrizione 2D della
  ! regola di parità 3D (piano x: R12,R13 dispari; piano y: R12,R23 dispari) con gli
  ! sforzi incrociati fuori piano assenti.
  ! var = (rho*R11, rho*R22, rho*R33, rho*R12, rho*omega); n è la normale di faccia.
  subroutine RSM_Symmetry ( var, var_w, var_g, n )
    use Q2D_Global_m, only: nRANS
    implicit none
    real(R8), intent(in),  dimension(nRANS) :: var
    real(R8), intent(out), dimension(nRANS) :: var_w, var_g
    real(R8), intent(in),  dimension(2)     :: n

    ! Sforzi normali + omega: gradiente nullo (pari)
    var_w(1:3) = var(1:3)
    var_w(5)   = var(5)
    var_g(1:3) = var(1:3)
    var_g(5)   = var(5)

    ! Shear nel piano R12: dispari su ogni piano di simmetria nel piano
    var_w(4) = 0d0
    var_g(4) = - var(4)

  end subroutine RSM_Symmetry

end module Q2D_Lib_RANS
