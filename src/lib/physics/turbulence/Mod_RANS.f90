!>@brief: RANS dispatcher (Setup_RANS_Model).
!
! Features:
!   - partial match `index(model, 'SA') > 0` to handle the model variants
!     (suffixes -R/-RC/-QCR2000/-blowcorr/-SD);
!   - binding of the Stress_Vector/Stress_Tensor/Strain_Tensor pointers for each
!     regime (Std/Axi; laminar, 2eq, RSM, QCR2000);
!   - binding of RANS_Diffusive_SourceQ2D (area-variation diffusive source);
!   - defense-in-depth `error stop` guards against axisymmetric+area_variation
!     in every viscous branch;
!   - keyword `'none'` for the laminar case, consistent with Assign_Setup and with
!     the Config_Types default `obj_rans%model`.
!
! Pre-dispatch flag initialization:
!   `obj_rans%SpalartShur/SAR/SAcomp/QCR2000/blowing_corr` derived from the substrings
!   `'-R'`/`'-RC'`/`'-QCR2000'`/`'-blowcorr'` of `obj_rans%model`.
!   The flags are NOT registered in Register_Physics; the canonical form is the parsing
!   (hard override: all reset to false, then derived from the string).
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
    ! Pre-dispatch flag initialization.
    ! Hard override: all reset to false, then parse the model string.
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
    ! Dispatcher for the selected model
    ! ====================================================================
    if ( index ( trim(obj_rans%model), 'none' ) > 0 ) then
      ! Laminar Navier-Stokes (no RANS model). Keyword 'none', consistent with
      ! Assign_Setup and with the Config_Types %model default.
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
      ! SSG/LRR-omega Reynolds-stress transport (4 stresses + omega, nRANS=5).
      ! Formal derivation per geometric regime (2D-plane/2D-axi/area variation):
      ! docs/theory/q2d-rans-closures.md §6 + Appendix A.
      ! No Boussinesq: coupling with momentum via Stress_Vector_RSM
      ! (= laminar mu_l minus rho*R_ij, ibid. §6.8). The convective geometric sources
      ! are added generically by Source_Axi/Area_Eul_Blk (post-pressure indices);
      ! the diffusive geometric source by RANS_Diffusive_SourceQ2D (ibid. §6.1).
      nRANS = 5
      obj_rans%RSM = .true.
      obj_rans%SD  = index ( trim(obj_rans%model), '-SD' ) > 0
      Eddy_Viscosity => SSGLRR_Eddy_Viscosity
      if ( obj_rans%SD ) then
        RANS_Diffusive_Flux => SSGLRR_SD_RANS_Diffusive_Flux
      else
        RANS_Diffusive_Flux => SSGLRR_RANS_Diffusive_Flux
      end if
      ! Mutual exclusion axi+area_var (defense in depth).
      if ( obj_source%axisymmetric .and. obj_source%area_variation ) then
        error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
      end if
      RANS_Diffusive_SourceQ2D => SSGLRR_Source_Diffusive_Q2D ! area-var diffusive source
      ! No Boussinesq: total stress = laminar(mu_l) - rho*R_ij (Stress_*_RSM).
      ! Strain_Tensor stays _Std/_Axi: it provides the block-diagonal S_33
      ! (0 / vr/r / -V·G) consumed by SSGLRR_Blk (docs/theory/q2d-rans-closures.md
      ! §6.1, Appendix A.1).
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
      ! Spalart-Allmaras (1 equation). Variants: SA, SA-R, SA-RC, SAcomp, SA-QCR2000.
      nRANS = 1
      Eddy_Viscosity      => Spalart_Eddy_Viscosity
      RANS_Diffusive_Flux => Spalart_RANS_Diffusive_Flux
      ! Mutual exclusion axi+area_var (defense in depth).
      if ( obj_source%axisymmetric .and. obj_source%area_variation ) then
        error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
      end if
      RANS_Diffusive_SourceQ2D => Spalart_Source_Diffusive_Q2D ! area-var diffusive source
      ! Stress dispatch: QCR2000 active only in the SA branch.
      ! Stress_Tensor must also be bound to a QCR variant: the cylindrical/area-
      ! variation geometric sources consume the full 3x3 tensor.
      ! Strain_Tensor stays Std/Axi (kinematic invariant: QCR is a constitutive
      ! correction, it does not modify the strain).
      if ( obj_rans%QCR2000 ) then
        if ( obj_source%axisymmetric ) then
          Stress_Vector => Stress_Vector_QCR2000_Axi
          Stress_Tensor => Stress_Tensor_QCR2000_Axi   ! cylindrical
          Strain_Tensor => Strain_Tensor_Axi           ! kinematic invariant (no QCR variant)
        else
          Stress_Vector => Stress_Vector_QCR2000_Std
          Stress_Tensor => Stress_Tensor_QCR2000_Std   ! area-var
          Strain_Tensor => Strain_Tensor_Std           ! kinematic invariant (no QCR variant)
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
      ! Wilcox 2006 k-omega (2 equations).
      nRANS = 2
      Eddy_Viscosity      => Wilcox2006_Eddy_Viscosity
      RANS_Diffusive_Flux => Wilcox2006_RANS_Diffusive_Flux
      if ( obj_source%axisymmetric .and. obj_source%area_variation ) then
        error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
      end if
      RANS_Diffusive_SourceQ2D => Wilcox2006_Source_Diffusive_Q2D ! area-var diffusive source
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
      ! Wilcox reuses SST_Blowing_Correction/noCorrection.
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall      => Wilcox2006_Extrapolate_Wall
      RANS_Enforce_Realizability => Wilcox2006_Enforce_Realizability

    elseif ( index ( trim(obj_rans%model), 'SST' ) > 0 ) then
      ! Menter k-omega SST (2 equations). Variants: SST, SST-RC, SST-blowcorr.
      nRANS = 2
      Eddy_Viscosity      => SST_Eddy_Viscosity
      RANS_Diffusive_Flux => SST_RANS_Diffusive_Flux
      if ( obj_source%axisymmetric .and. obj_source%area_variation ) then
        error stop 'Q2D: "area_variation" and "axisymmetric" cannot be present together.'
      end if
      RANS_Diffusive_SourceQ2D => SST_Source_Diffusive_Q2D ! area-var diffusive source
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
      ! Inviscid flow (Euler).
      nRANS = 0
      Eddy_Viscosity           => null()
      RANS_Diffusive_Flux      => null()
      RANS_Diffusive_SourceQ2D => null()
      Stress_Vector            => null()
      obj_rans%model           = '  Inviscid flow'

    endif

  end subroutine Setup_RANS_Model

end module Q2D_Mod_RANS
