module Q2D_Mod_Riemann
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Assign_Riemann_Solver

  procedure(riemann_if), pointer, public :: Riemann

  abstract interface
    subroutine riemann_if(dl,ul,vl,pl,al,dltot,dr,ur,vr,pr,ar,drtot,beta,nx,ny,F_r,F_u,F_v,F_E)
      use iso_fortran_env, only: R8 => real64
      use Q2D_Global_m, only: nsc
      implicit none
      real(R8), intent(in)  :: dl(nsc),ul,vl,pl,al
      real(R8), intent(in)  :: dr(nsc),ur,vr,pr,ar
      real(R8), intent(in)  :: dltot,drtot
      real(R8), intent(in)  :: beta, nx, ny
      real(R8), intent(out) :: F_r, F_u, F_v, F_E
    end subroutine riemann_if
  end interface

contains

  subroutine Assign_Riemann_Solver()
    use Q2D_Config_Types_m,    only: obj_riemann
    use Q2D_Lib_Riemann_AUSM,    only: riemann_Hanel, riemann_AUSMp, &
                                       riemann_AUSMp_up, riemann_AUSMp_up2
    use Q2D_Lib_Riemann_Godunov, only: riemann_exact
    use Q2D_Lib_Riemann_HLL,     only: riemann_HLLE, riemann_HLLEM, riemann_HLLC, &
                                       riemann_HLLEpp, riemann_HLLEMSD, riemann_HLLCSD, &
                                       riemann_HLLCHLLE
    use Q2D_Lib_Riemann_LF,      only: riemann_LLF, Riemann_PLLF
    use Q2D_Lib_Riemann_SLAU,    only: riemann_SLAU, riemann_SLAU2
    implicit none

    nullify(Riemann)

    select case (obj_riemann%description)

      ! Solutori di tipo AUSM
      case ('Hanel')
        Riemann => riemann_Hanel
        obj_riemann%description = 'Hanel'
      case ('AUSM+')
        Riemann => riemann_AUSMp
        obj_riemann%description = 'AUSM+'
      case ('AUSM+-up')
        Riemann => riemann_AUSMp_up
        obj_riemann%description = 'AUSM+-up'
      case ('AUSM+-up2')
        Riemann => riemann_AUSMp_up2
        obj_riemann%description = 'AUSM+-up2'

      ! Solutori di Godunov
      case ('exact','Exact')
        Riemann => riemann_exact
        obj_riemann%description = 'Exact Solver'

      ! Solutori di tipo HLL
      case ('HLLE')
        Riemann => riemann_HLLE
        obj_riemann%description = 'HLLE'
      case ('HLLEM')
        Riemann => riemann_HLLEM
        obj_riemann%description = 'HLLEM'
      case ('HLLC')
        ! Unica parola chiave di input accettata: 'HLLC' (la descrizione riportata
        ! è 'HLLC Batten'); la stringa 'HLLC Batten' non è una parola chiave valida.
        Riemann => riemann_HLLC
        obj_riemann%description = 'HLLC Batten'
      case ('HLLC+')
        ! HLLC+ di Tramel — usa il rilevatore d'urto: blend HLLC+HLLE tramite beta.
        obj_riemann%SD = .true.
        Riemann => riemann_HLLCSD
        obj_riemann%description = 'Tramel HLLC+'
      case ('HLLE++')
        ! HLLE++ di Tramel — usa il rilevatore d'urto: blend HLLE+Roe tramite beta.
        obj_riemann%SD = .true.
        Riemann => riemann_HLLEpp
        obj_riemann%description = 'Tramel HLLE++'
      case ('HLLC Rotated')
        Riemann => riemann_HLLCHLLE
        obj_riemann%description = 'Rotated HLLC Batten / HLLE'

      ! Solutori di tipo Lax-Friedrichs
      case ('LLF','Rusanov')
        Riemann => riemann_LLF
        obj_riemann%description = 'Local Lax-Friedrichs (Rusanov)'
      case ('PLLF')
        Riemann => Riemann_PLLF
        obj_riemann%description = 'Preconditioned Local Lax-Friedrichs'

      ! Solutori di tipo SLAU
      case ('SLAU')
        Riemann => riemann_SLAU
        obj_riemann%description = 'SLAU'
      case ('SLAU2')
        Riemann => riemann_SLAU2
        obj_riemann%description = 'SLAU2'

      ! Caso di default
      case default
        Riemann => riemann_HLLC
        obj_riemann%description = 'HLLC Batten'

    end select

  end subroutine Assign_Riemann_Solver

end module Q2D_Mod_Riemann
