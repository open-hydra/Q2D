!>@brief: Dispatcher RANS (Setup_RANS_Model).
!
! Caratteristiche:
!   - match parziale `index(model, 'SA') > 0` per gestire le varianti del modello
!     (suffissi -R/-RC/-QCR2000/-blowcorr/-SD);
!   - binding dei puntatori Stress_Vector/Stress_Tensor/Strain_Tensor per ogni
!     regime (Std/Axi; laminare, 2eq, RSM, QCR2000);
!   - binding di RANS_Diffusive_SourceQ2D (sorgente diffusiva di variazione d'area);
!   - guardie `error stop` di difesa in profondità contro axisymmetric+area_variation
!     in ogni ramo viscoso;
!   - keyword `'none'` per il caso laminare, coerente con Assign_Setup e con il
!     default di Config_Types `obj_rans%model`.
!
! Inizializzazione dei flag pre-dispatch:
!   `obj_rans%SpalartShur/SAR/SAcomp/QCR2000/blowing_corr` ricavati dalle sottostringhe
!   `'-R'`/`'-RC'`/`'-QCR2000'`/`'-blowcorr'` di `obj_rans%model`.
!   I flag NON sono registrati in Register_Physics; la forma canonica è il parsing
!   (override rigido: tutti azzerati a false, poi ricavati dalla stringa).
module Q2D_Mod_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_RANS_Model

contains

  subroutine Setup_RANS_Model ()
    use Q2D_Config_Types_m, only: obj_rans, obj_source
    use Q2D_Global_m,       only: nRANS
    use Q2D_Lib_RANS,       only: Eddy_Viscosity, RANS_Diffusive_Flux, RANS_Diffusive_SourceQ2D, &
                                  RANS_Source_Terms, RANS_Set_Wall_Values, RANS_Set_Blowing_Wall, &
                                  RANS_Extrapolate_Wall, RANS_Enforce_Realizability
    use Q2D_Lib_Fluid,      only: Stress_Vector, Stress_Tensor, Strain_Tensor, &
                                  Stress_Vector_Std, Stress_Vector_Axi, &
                                  Stress_Vector_2eq_Std, Stress_Vector_2eq_Axi, &
                                  Stress_Vector_RSM_Std, Stress_Vector_RSM_Axi, &
                                  Stress_Tensor_Std, Stress_Tensor_Axi, &
                                  Stress_Tensor_2eq_Std, Stress_Tensor_2eq_Axi, &
                                  Stress_Tensor_RSM_Std, Stress_Tensor_RSM_Axi, &
                                  Strain_Tensor_Std, Strain_Tensor_Axi
    use Q2D_Lib_SSGLRR,     only: SSGLRR_Eddy_Viscosity, SSGLRR_RANS_Diffusive_Flux, &
                                  SSGLRR_SD_RANS_Diffusive_Flux, SSGLRR_Source_Diffusive_Q2D, &
                                  SSGLRR_Source_Terms, SSGLRR_Set_Wall_Values, &
                                  SSGLRR_Extrapolate_Wall, SSGLRR_Enforce_Realizability
    use Q2D_Lib_QCR2000,    only: Stress_Vector_QCR2000_Std, Stress_Vector_QCR2000_Axi, &
                                  Stress_Tensor_QCR2000_Std, Stress_Tensor_QCR2000_Axi
    use Q2D_Lib_Spalart,    only: Spalart_Eddy_Viscosity, Spalart_RANS_Diffusive_Flux, &
                                  Spalart_Source_Diffusive_Q2D, Spalart_Source_Terms, &
                                  Spalart_Set_Wall_Values, Spalart_Set_Blowing_Wall, &
                                  Spalart_Extrapolate_Wall, Spalart_Enforce_Realizability
    use Q2D_Lib_SST,        only: SST_Eddy_Viscosity, SST_RANS_Diffusive_Flux, &
                                  SST_Source_Diffusive_Q2D, SST_Source_Terms, &
                                  SST_Set_Wall_Values, SST_Blowing_Correction, SST_Blowing_noCorrection, &
                                  SST_Extrapolate_Wall, SST_Enforce_Realizability
    use Q2D_Lib_Wilcox2006, only: Wilcox2006_Eddy_Viscosity, Wilcox2006_RANS_Diffusive_Flux, &
                                  Wilcox2006_Source_Diffusive_Q2D, Wilcox2006_Source_Terms, &
                                  Wilcox2006_Set_Wall_Values, &
                                  Wilcox2006_Extrapolate_Wall, Wilcox2006_Enforce_Realizability

    implicit none

    ! ====================================================================
    ! Inizializzazione dei flag pre-dispatch.
    ! Override rigido: tutti azzerati a false, poi parsing della stringa modello.
    ! ====================================================================
    obj_rans%SpalartShur  = .false.
    obj_rans%SAR          = .false.
    obj_rans%SAcomp       = .false.
    obj_rans%QCR2000      = .false.
    obj_rans%blowing_corr = .false.

    if (trim(obj_rans%model) == 'SAcomp') obj_rans%SAcomp = .true.

    if (index(trim(obj_rans%model), '-R') > 0) then
      if (index(trim(obj_rans%model), '-RC') > 0) then
        obj_rans%SpalartShur = .true.
      else
        obj_rans%SAR = .true.
      end if
    end if

    if (index(trim(obj_rans%model), '-QCR2000')  > 0) obj_rans%QCR2000      = .true.
    if (index(trim(obj_rans%model), '-blowcorr') > 0) obj_rans%blowing_corr = .true.

    obj_rans%description = 'RANS model: '//trim(obj_rans%model)

    ! ====================================================================
    ! Dispatcher del modello selezionato
    ! ====================================================================
    if ( index ( trim(obj_rans%model), 'none' ) > 0 ) then
      ! Navier-Stokes laminare (nessun modello RANS). Keyword 'none', coerente con
      ! Assign_Setup e con il default %model di Config_Types.
      nRANS = 0
      Eddy_Viscosity              => null()
      RANS_Diffusive_Flux         => null()
      RANS_Diffusive_SourceQ2D    => null()
      if ( obj_source%axisymmetric ) then
        Stress_Vector => Stress_Vector_Axi
        Stress_Tensor => Stress_Tensor_Axi
        Strain_Tensor => Strain_Tensor_Axi
      else
        Stress_Vector => Stress_Vector_Std
        Stress_Tensor => Stress_Tensor_Std
        Strain_Tensor => Strain_Tensor_Std
      end if
      RANS_Enforce_Realizability => null()

    elseif ( index ( trim(obj_rans%model), 'SSGLRR' ) > 0 ) then
      ! Trasporto degli sforzi di Reynolds SSG/LRR-omega (4 sforzi + omega, nRANS=5).
      ! Derivazione formale per regime geometrico (2D-piano/2D-axi/variazione d'area):
      ! docs/theory/q2d-rans-closures.md §6 + Appendice A.
      ! Niente Boussinesq: accoppiamento con la quantità di moto via Stress_Vector_RSM
      ! (= mu_l laminare meno rho*R_ij, ibid. §6.8). Le sorgenti geometriche convettive
      ! sono aggiunte genericamente da Source_Axi/Area_Eul_Blk (indici post-pressione);
      ! la sorgente geometrica diffusiva da RANS_Diffusive_SourceQ2D (ibid. §6.1).
      nRANS = 5
      obj_rans%RSM = .true.
      obj_rans%SD  = index ( trim(obj_rans%model), '-SD' ) > 0
      Eddy_Viscosity => SSGLRR_Eddy_Viscosity
      if ( obj_rans%SD ) then
        RANS_Diffusive_Flux => SSGLRR_SD_RANS_Diffusive_Flux
      else
        RANS_Diffusive_Flux => SSGLRR_RANS_Diffusive_Flux
      end if
      ! Mutua esclusione axi+area_var (difesa in profondità).
      if ( obj_source%axisymmetric .and. obj_source%area_variation ) then
        error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
      end if
      RANS_Diffusive_SourceQ2D => SSGLRR_Source_Diffusive_Q2D ! sorgente diffusiva area-var
      ! Niente Boussinesq: sforzo totale = laminare(mu_l) - rho*R_ij (Stress_*_RSM).
      ! Strain_Tensor resta _Std/_Axi: fornisce la S_33 diagonale a blocchi
      ! (0 / vr/r / -V·G) consumata da SSGLRR_Blk (docs/theory/q2d-rans-closures.md
      ! §6.1, Appendice A.1).
      if ( obj_source%axisymmetric ) then
        Stress_Vector => Stress_Vector_RSM_Axi
        Stress_Tensor => Stress_Tensor_RSM_Axi
        Strain_Tensor => Strain_Tensor_Axi
      else
        Stress_Vector => Stress_Vector_RSM_Std
        Stress_Tensor => Stress_Tensor_RSM_Std
        Strain_Tensor => Strain_Tensor_Std
      end if
      RANS_Source_Terms          => SSGLRR_Source_Terms
      RANS_Set_Wall_Values       => SSGLRR_Set_Wall_Values
      RANS_Extrapolate_Wall      => SSGLRR_Extrapolate_Wall
      RANS_Enforce_Realizability => SSGLRR_Enforce_Realizability

    elseif ( index ( trim(obj_rans%model), 'SA' ) > 0 ) then
      ! Spalart-Allmaras (1 equazione). Varianti: SA, SA-R, SA-RC, SAcomp, SA-QCR2000.
      nRANS = 1
      Eddy_Viscosity      => Spalart_Eddy_Viscosity
      RANS_Diffusive_Flux => Spalart_RANS_Diffusive_Flux
      ! Mutua esclusione axi+area_var (difesa in profondità).
      if ( obj_source%axisymmetric .and. obj_source%area_variation ) then
        error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
      end if
      RANS_Diffusive_SourceQ2D => Spalart_Source_Diffusive_Q2D ! sorgente diffusiva area-var
      ! Dispatch degli sforzi: QCR2000 attivo solo nel ramo SA.
      ! Anche Stress_Tensor va legato a una variante QCR: le sorgenti geometriche
      ! cilindriche/di variazione d'area consumano il tensore 3x3 completo.
      ! Strain_Tensor resta Std/Axi (invariante cinematico: QCR è una correzione
      ! costitutiva, non modifica lo strain).
      if ( obj_rans%QCR2000 ) then
        if ( obj_source%axisymmetric ) then
          Stress_Vector => Stress_Vector_QCR2000_Axi
          Stress_Tensor => Stress_Tensor_QCR2000_Axi   ! cilindrico
          Strain_Tensor => Strain_Tensor_Axi           ! invariante cinematico (nessuna variante QCR)
        else
          Stress_Vector => Stress_Vector_QCR2000_Std
          Stress_Tensor => Stress_Tensor_QCR2000_Std   ! area-var
          Strain_Tensor => Strain_Tensor_Std           ! invariante cinematico (nessuna variante QCR)
        end if
      else
        if ( obj_source%axisymmetric ) then
          Stress_Vector => Stress_Vector_Axi
          Stress_Tensor => Stress_Tensor_Axi
          Strain_Tensor => Strain_Tensor_Axi
        else
          Stress_Vector => Stress_Vector_Std
          Stress_Tensor => Stress_Tensor_Std
          Strain_Tensor => Strain_Tensor_Std
        end if
      end if
      RANS_Source_Terms          => Spalart_Source_Terms
      RANS_Set_Wall_Values       => Spalart_Set_Wall_Values
      RANS_Set_Blowing_Wall      => Spalart_Set_Blowing_Wall
      RANS_Extrapolate_Wall      => Spalart_Extrapolate_Wall
      RANS_Enforce_Realizability => Spalart_Enforce_Realizability

    elseif ( index ( trim(obj_rans%model), 'Wilcox2006' ) > 0 ) then
      ! Wilcox 2006 k-omega (2 equazioni).
      nRANS = 2
      Eddy_Viscosity      => Wilcox2006_Eddy_Viscosity
      RANS_Diffusive_Flux => Wilcox2006_RANS_Diffusive_Flux
      if ( obj_source%axisymmetric .and. obj_source%area_variation ) then
        error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
      end if
      RANS_Diffusive_SourceQ2D => Wilcox2006_Source_Diffusive_Q2D ! sorgente diffusiva area-var
      if (obj_rans%k_energy_coupling) then
        if ( obj_source%axisymmetric ) then
          Stress_Vector => Stress_Vector_2eq_Axi
          Stress_Tensor => Stress_Tensor_2eq_Axi
          Strain_Tensor => Strain_Tensor_Axi
        else
          Stress_Vector => Stress_Vector_2eq_Std
          Stress_Tensor => Stress_Tensor_2eq_Std
          Strain_Tensor => Strain_Tensor_Std
        end if
      else
        if ( obj_source%axisymmetric ) then
          Stress_Vector => Stress_Vector_Axi
          Stress_Tensor => Stress_Tensor_Axi
          Strain_Tensor => Strain_Tensor_Axi
        else
          Stress_Vector => Stress_Vector_Std
          Stress_Tensor => Stress_Tensor_Std
          Strain_Tensor => Strain_Tensor_Std
        end if
      end if
      RANS_Source_Terms    => Wilcox2006_Source_Terms
      RANS_Set_Wall_Values => Wilcox2006_Set_Wall_Values
      ! Wilcox riusa SST_Blowing_Correction/noCorrection.
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall      => Wilcox2006_Extrapolate_Wall
      RANS_Enforce_Realizability => Wilcox2006_Enforce_Realizability

    elseif ( index ( trim(obj_rans%model), 'SST' ) > 0 ) then
      ! Menter k-omega SST (2 equazioni). Varianti: SST, SST-RC, SST-blowcorr.
      nRANS = 2
      Eddy_Viscosity      => SST_Eddy_Viscosity
      RANS_Diffusive_Flux => SST_RANS_Diffusive_Flux
      if ( obj_source%axisymmetric .and. obj_source%area_variation ) then
        error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
      end if
      RANS_Diffusive_SourceQ2D => SST_Source_Diffusive_Q2D ! sorgente diffusiva area-var
      if (obj_rans%k_energy_coupling) then
        if ( obj_source%axisymmetric ) then
          Stress_Vector => Stress_Vector_2eq_Axi
          Stress_Tensor => Stress_Tensor_2eq_Axi
          Strain_Tensor => Strain_Tensor_Axi
        else
          Stress_Vector => Stress_Vector_2eq_Std
          Stress_Tensor => Stress_Tensor_2eq_Std
          Strain_Tensor => Strain_Tensor_Std
        end if
      else
        if ( obj_source%axisymmetric ) then
          Stress_Vector => Stress_Vector_Axi
          Stress_Tensor => Stress_Tensor_Axi
          Strain_Tensor => Strain_Tensor_Axi
        else
          Stress_Vector => Stress_Vector_Std
          Stress_Tensor => Stress_Tensor_Std
          Strain_Tensor => Strain_Tensor_Std
        end if
      end if
      RANS_Source_Terms    => SST_Source_Terms
      RANS_Set_Wall_Values => SST_Set_Wall_Values
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall      => SST_Extrapolate_Wall
      RANS_Enforce_Realizability => SST_Enforce_Realizability

    else
      ! Flusso non viscoso (Eulero).
      nRANS = 0
      Eddy_Viscosity           => null()
      RANS_Diffusive_Flux      => null()
      RANS_Diffusive_SourceQ2D => null()
      Stress_Vector            => null()
      obj_rans%model           = '  Inviscid flow'

    endif

  end subroutine Setup_RANS_Model

end module Q2D_Mod_RANS
