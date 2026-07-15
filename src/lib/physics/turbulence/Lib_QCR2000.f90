!>@brief Quadratic Constitutive Relations (Spalart 2000) — riduzione 2D per Q2D.
!
! Riferimento: Spalart P.R. (2000), "Strategies for turbulence modelling and simulations",
!              Int. J. Heat Fluid Flow 21(3):252-263.
!
! Limite di applicabilità in 2D: la QCR nasce per le separazioni d'angolo 3D
! (giunzioni ala-fusoliera, interazioni al mozzo delle palettature). In Q2D
! (Tau_xz=Tau_yz=0 a priori perché i gradienti di velocità fuori piano sono
! assenti) la QCR agisce solo sull'anisotropia nel piano (Tau_xx, Tau_yy,
! Tau_xy): il beneficio fisico è quindi ridotto rispetto al caso d'uso 3D
! d'origine. Per la RANS 2D ingegneristica i modelli primari restano SST/SA;
! la QCR è raccomandata come correzione opzionale per casi con anisotropia
! significativa.
!
! Q2D fornisce 2 realizzazioni (Stress_Vector_QCR2000_Std/Axi) secondo lo
! schema di commutazione geometrica Std/Axi di Stress_Vector_Std/Axi in
! Lib_Fluid. L'interfaccia astratta Stress_Vector_if di Q2D richiede 8
! argomenti (Gradvel(2,2), Stress(2)), quindi l'eventuale forma unica a 5
! argomenti (Gradvel(3,3), Normal(3), Stress(3)) non è assegnabile al
! procedure pointer — da qui la separazione Std/Axi.
!
! Dispatch: QCR2000 si aggancia solo nel ramo SA; per SST/Wilcox, quando
! k_energy_coupling è attivo, si usa il percorso Stress_Vector_2eq, quindi la
! QCR è mutuamente esclusiva con esso.
module Q2D_Lib_QCR2000
  use iso_fortran_env,    only: I4 => int32, R8 => real64
  use Q2D_Config_Types_m, only: obj_rot   ! acceduto in sola lettura nelle funzioni pure sottostanti

  implicit none
  private

  real(R8), parameter, private :: C_cr1             = 0.3d0
  real(R8), parameter, private :: negativetwothird  = -2.d0/3.d0

  public :: Stress_Vector_QCR2000_Std, Stress_Vector_QCR2000_Axi
  public :: Stress_Tensor_QCR2000_Std, Stress_Tensor_QCR2000_Axi

contains

  pure function Stress_Vector_QCR2000_Std ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! Firma a 8 argomenti conforme a Stress_Vector_if (Lib_Fluid).
    ! Definizione A di Q2D: Divervel = G(1,1) + G(2,2) - VG (incorporata via VG).
    ! 2D piano: Tau_xz_turb = Tau_yz_turb = 0 a priori
    ! (nessun gradiente di velocità fuori piano) -> correzione QCR solo su
    ! Tau_xx, Tau_yy, Tau_xy, e Tau_zz_QCR = Tau_zz (nessun mescolamento fuori piano).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t                ! turbolento (deviatorico di Boussinesq)
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR          ! turbolento corretto QCR
    real(R8) :: Tau_xx, Tau_yy, Tau_xy                      ! totale laminare + turbolento
    real(R8) :: lam_l, Divervel, lamDivervel, VG
    real(R8) :: omega_z, Den, O_off

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    ! Divergenza secondo la Definizione A di Q2D
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG

    ! Sforzo turbolento (deviatorico di Boussinesq, con mi_t). 2D piano.
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! Componente z della vorticità (omega(1)=omega(2)=0)
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den = sqrt(sum(Gradvel_full**2)) include la componente (3,3)=-VG.
    ! Evita la divisione per zero nelle regioni laminari.
    Den = sqrt ( sum ( Gradvel**2 ) + VG**2 ) + 1d-40
    ! ── Oggettività in sistema rotante (verificata su un vortice analitico,
    ! assi fissi vs assi rotanti): O_ik = 2 W_ik/||grad u|| va costruito con
    ! grandezze in assi ASSOLUTI — omega_z traslata di 2*Oz (come ogni altro
    ! consumatore di vorticità) E Den dal gradiente assoluto = relativo +
    ! gradiente di trascinamento d(Omega x r)_a/dx_b = [[0,-Oz],[Oz,0]] nel
    ! piano (la componente (3,3) è invariante di sistema). La sola traslazione
    ! di omega NON basta: un Den dal gradiente relativo lascia tau_QCR
    ! dipendente dal sistema di riferimento e sovracorregge; con entrambe le
    ! correzioni tau_QCR è esattamente oggettivo (zero macchina).
    if (obj_rot%enabled) then
      omega_z = omega_z + 2.0_R8*obj_rot%Omega_vec(3)
      Den = sqrt ( Gradvel(1,1)**2 + (Gradvel(1,2)-obj_rot%Omega_vec(3))**2 &
                 + (Gradvel(2,1)+obj_rot%Omega_vec(3))**2 + Gradvel(2,2)**2 &
                 + VG**2 ) + 1d-40
    end if

    ! O_ij(1,2) = 2*W_ij(1,2)/Den = -omega_z/Den; O_ij(2,1) = +omega_z/Den.
    ! Modulo fuori diagonale in senso di componente z (uniche entrate non nulle in 2D).
    O_off = omega_z / Den   ! = -O_ij(1,2) = +O_ij(2,1)

    ! Correzione QCR (le 6 componenti si riducono alle 3 nel piano):
    ! T_ij_QCR = T_ij - C_cr1*[O_ik*T_jk + O_jk*T_ik], somme su k = {1,2,3}.
    ! In 2D piano: O_ik = 0 tranne (1,2),(2,1); T_ik = 0 per k=3 nello sforzo turbolento.
    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1   ! = Tau_xx_t + 2*O_off*Tau_xy_t*C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1
                                                                ! = Tau_xy_t + O_off*(Tau_yy_t - Tau_xx_t)*C_cr1

    ! Laminare (newtoniano lineare, senza correzione QCR). Speculare a Stress_Vector_Std.
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_QCR2000_Std


  pure function Stress_Vector_QCR2000_Axi ( Gradvel, Normal, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress )
    ! Firma a 8 argomenti conforme a Stress_Vector_if. Divergenza cilindrica.
    ! Argomento opzionale V_dot_G per uniformità di interfaccia (NON usato in
    ! Axi: i flussi assialsimmetrici non hanno parete virtuale top/bot nel
    ! senso della Definizione A).
    ! 2D-axi No-w (senza swirl): Tau_xz_turb = Tau_yz_turb = 0, e
    ! Tau_zz_turb (theta-theta) NON entra nello Stress(1:2) nel piano.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), Normal(2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress(2)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR
    real(R8) :: Tau_xx, Tau_yy, Tau_xy
    real(R8) :: lam_l, Divervel, lamDivervel, G_zz
    real(R8) :: omega_z, Den, O_off

    ! Divergenza in coordinate cilindriche (fallback di de l'Hopital all'asse r_ave -> 0)
    if (r_ave .le. 1.d-6) then
      Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(2,2)  ! de l'Hopital vr/r --> dv/dr
      G_zz     = Gradvel(2,2)
    else
      Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave
      G_zz     = vr/r_ave
    end if

    ! Sforzo turbolento
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! Componente z della vorticità
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den include la componente (3,3) = vr/r_ave (axi)
    Den = sqrt ( sum ( Gradvel**2 ) + G_zz**2 ) + 1d-40
    ! ── Oggettività in sistema rotante — vedi il commento della variante Std
    ! (oggettività verificata su un vortice analitico): omega assoluta E
    ! gradiente assoluto in Den; la componente cilindrica (3,3) = vr/r è
    ! invariante di sistema. (La combinazione 2D-axi + sistema rotante è oggi
    ! bloccata dal controllo di compatibilità in Setup_RotatingFrame; questo
    ! ramo mantiene la realizzazione Axi coerente col sistema di riferimento
    ! in vista dell'estensione 2D-axi con swirl.)
    if (obj_rot%enabled) then
      omega_z = omega_z + 2.0_R8*obj_rot%Omega_vec(3)
      Den = sqrt ( Gradvel(1,1)**2 + (Gradvel(1,2)-obj_rot%Omega_vec(3))**2 &
                 + (Gradvel(2,1)+obj_rot%Omega_vec(3))**2 + Gradvel(2,2)**2 &
                 + G_zz**2 ) + 1d-40
    end if
    O_off = omega_z / Den

    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1

    ! Laminare
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    Stress(1) = dot_product ( [Tau_xx, Tau_xy], Normal )
    Stress(2) = dot_product ( [Tau_xy, Tau_yy], Normal )
  end function Stress_Vector_QCR2000_Axi


  pure function Stress_Tensor_QCR2000_Std ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! Tensore degli sforzi 3x3 con correzione QCR2000 (Spalart 2000), variante Std.
    ! Speculare allo schema di Stress_Vector_QCR2000_Std con:
    !   - firma a 7 argomenti conforme a Stress_Tensor_if (niente Normal, uscita 3x3
    !     invece di vector(2)·n).
    !   - Definizione A di Q2D: Divervel = G(1,1) + G(2,2) - VG (incorporata via VG).
    !   - correzione QCR sulle componenti nel piano (xx, yy, xy); Tau_QCR_zz = Tau_zz
    !     INVARIATO (No-w → O_θ·=0 → nessuna QCR sulla componente (3,3)).
    !   - Tau_xz, Tau_yz = 0 (No-w).
    ! Verifica di coerenza: O_off=0 → Tau_QCR_xx/yy/xy = Tau_xx/yy/xy non-QCR ✓.
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t          ! turbolento (deviatorico di Boussinesq, mi_t)
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR    ! turbolento corretto QCR
    real(R8) :: Tau_xx, Tau_yy, Tau_xy, Tau_zz        ! completo (laminare + turbolento_QCR)
    real(R8) :: lam_l, Divervel, lamDivervel, VG
    real(R8) :: omega_z, Den, O_off
    real(R8) :: mi_combined

    VG = 0.d0
    if (present(V_dot_G)) VG = V_dot_G

    ! Divergenza secondo la Definizione A di Q2D (incorporata via VG)
    Divervel = Gradvel(1,1) + Gradvel(2,2) - VG

    ! Sforzo turbolento (deviatorico di Boussinesq, con mi_t). 2D piano.
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! Componente z della vorticità (omega(1)=omega(2)=0)
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den include la componente (3,3) = -VG (Definizione A con variazione d'area).
    Den = sqrt ( sum ( Gradvel**2 ) + VG**2 ) + 1d-40
    ! ── Oggettività in sistema rotante (verificata su un vortice analitico,
    ! assi fissi vs assi rotanti): O_ik = 2 W_ik/||grad u|| va costruito con
    ! grandezze in assi ASSOLUTI — omega_z traslata di 2*Oz (come ogni altro
    ! consumatore di vorticità) E Den dal gradiente assoluto = relativo +
    ! gradiente di trascinamento d(Omega x r)_a/dx_b = [[0,-Oz],[Oz,0]] nel
    ! piano (la componente (3,3) è invariante di sistema). La sola traslazione
    ! di omega NON basta: un Den dal gradiente relativo lascia tau_QCR
    ! dipendente dal sistema di riferimento e sovracorregge; con entrambe le
    ! correzioni tau_QCR è esattamente oggettivo (zero macchina).
    if (obj_rot%enabled) then
      omega_z = omega_z + 2.0_R8*obj_rot%Omega_vec(3)
      Den = sqrt ( Gradvel(1,1)**2 + (Gradvel(1,2)-obj_rot%Omega_vec(3))**2 &
                 + (Gradvel(2,1)+obj_rot%Omega_vec(3))**2 + Gradvel(2,2)**2 &
                 + VG**2 ) + 1d-40
    end if
    O_off = omega_z / Den   ! = -O_ij(1,2) = +O_ij(2,1)

    ! Correzione QCR (le 6 componenti si riducono alle 3 nel piano):
    ! T_ij_QCR = T_ij - C_cr1*[O_ik*T_jk + O_jk*T_ik], somme su k.
    ! O_θ·=0 (No-w) → nessuna correzione QCR sulla componente (3,3).
    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1

    ! Laminare (newtoniano lineare, senza correzione QCR).
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    ! Nel piano: laminare + turbolento corretto QCR
    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    ! Tau_zz: V·G incorporato + NESSUNA correzione QCR (No-w → O_θ·=0).
    ! Speculare a Stress_Tensor_Std: Tau_zz = lamDivervel_combined + 2*mi*(-VG).
    mi_combined = mi_l + mi_t
    Tau_zz = (negativetwothird * mi_combined) * Divervel + 2d0 * mi_combined * (-VG)

    ! Assemblaggio del tensore 3x3 (τ_xz = τ_yz = 0)
    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz

    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy

    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_QCR2000_Std


  pure function Stress_Tensor_QCR2000_Axi ( Gradvel, mi_l, mi_t, var, r_ave, vr, V_dot_G ) result ( Stress_t )
    ! Tensore degli sforzi 3x3 con correzione QCR2000, variante Axi.
    ! Speculare a Stress_Tensor_Axi (Lib_Fluid) + QCR sulle componenti nel piano
    ! (xx, yy, xy); Tau_QCR_θθ = Tau_θθ INVARIATO (No-w → O_θ·=0).
    ! Argomento opzionale V_dot_G per uniformità di interfaccia (NON usato in Axi).
    ! Tensore valutato al centro cella (r_ave > 0 per costruzione della mesh) →
    ! nessun fallback di de l'Hopital (stessa convenzione di Stress_Tensor_Axi).
    implicit none
    real(R8), intent(in) :: Gradvel(2,2), mi_t, mi_l, var(:), r_ave, vr
    real(R8), intent(in), optional :: V_dot_G
    real(R8) :: Stress_t(3,3)
    ! Local
    real(R8) :: Tau_xx_t, Tau_xy_t, Tau_yy_t
    real(R8) :: Tau_xx_QCR, Tau_yy_QCR, Tau_xy_QCR
    real(R8) :: Tau_xx, Tau_yy, Tau_xy, Tau_zz
    real(R8) :: lam_l, Divervel, lamDivervel, G_zz
    real(R8) :: omega_z, Den, O_off
    real(R8) :: mi_combined

    ! Divergenza in coordinate cilindriche (V_r/r nella componente (3,3)).
    Divervel = Gradvel(1,1) + Gradvel(2,2) + vr/r_ave
    G_zz     = vr/r_ave

    ! Sforzo turbolento (deviatorico di Boussinesq, con mi_t).
    Tau_xx_t = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy_t = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_xy_t =       mi_t * ( Gradvel(1,2) + Gradvel(2,1)        )

    ! Componente z della vorticità
    omega_z = Gradvel(2,1) - Gradvel(1,2)

    ! Den include la componente (3,3) = vr/r_ave (termine di Christoffel cilindrico).
    Den = sqrt ( sum ( Gradvel**2 ) + G_zz**2 ) + 1d-40
    ! ── Oggettività in sistema rotante — vedi il commento della variante Std
    ! (oggettività verificata su un vortice analitico): omega assoluta E
    ! gradiente assoluto in Den; la componente cilindrica (3,3) = vr/r è
    ! invariante di sistema. (La combinazione 2D-axi + sistema rotante è oggi
    ! bloccata dal controllo di compatibilità in Setup_RotatingFrame; questo
    ! ramo mantiene la realizzazione Axi coerente col sistema di riferimento
    ! in vista dell'estensione 2D-axi con swirl.)
    if (obj_rot%enabled) then
      omega_z = omega_z + 2.0_R8*obj_rot%Omega_vec(3)
      Den = sqrt ( Gradvel(1,1)**2 + (Gradvel(1,2)-obj_rot%Omega_vec(3))**2 &
                 + (Gradvel(2,1)+obj_rot%Omega_vec(3))**2 + Gradvel(2,2)**2 &
                 + G_zz**2 ) + 1d-40
    end if
    O_off = omega_z / Den

    ! Correzione QCR (solo nel piano, O_θ·=0)
    Tau_xx_QCR = Tau_xx_t - 2d0 * (-O_off) * Tau_xy_t * C_cr1
    Tau_yy_QCR = Tau_yy_t - 2d0 * ( O_off) * Tau_xy_t * C_cr1
    Tau_xy_QCR = Tau_xy_t - ( (-O_off) * Tau_yy_t + ( O_off) * Tau_xx_t ) * C_cr1

    ! Laminare
    lam_l = negativetwothird * mi_l
    lamDivervel = lam_l * Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_xy =               mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR

    ! Tau_zz (= τ_θθ in lettura cilindrica): vr/r incorporato + NESSUNA QCR (No-w).
    ! Speculare a Stress_Tensor_Axi: Tau_zz = 2*mi*vr/r_ave + lamDivervel_combined.
    mi_combined = mi_l + mi_t
    Tau_zz = (negativetwothird * mi_combined) * Divervel + 2d0 * mi_combined * G_zz

    ! Assemblaggio del tensore 3x3 (fuori diagonale τ_zr/rz nulli)
    Stress_t(1,1) = Tau_xx
    Stress_t(2,2) = Tau_yy
    Stress_t(3,3) = Tau_zz

    Stress_t(1,2) = Tau_xy
    Stress_t(2,1) = Tau_xy

    Stress_t(1,3) = 0.d0
    Stress_t(3,1) = 0.d0
    Stress_t(2,3) = 0.d0
    Stress_t(3,2) = 0.d0
  end function Stress_Tensor_QCR2000_Axi


end module Q2D_Lib_QCR2000
