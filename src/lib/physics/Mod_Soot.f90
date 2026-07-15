!>@brief: Dispatcher del soot — Setup_Soot, Soot_Source_Terms, Soot_Diffusive_Flux.
!
! 3 routine pubbliche. Convenzioni 2D:
!   - loop di blocco sui soli indici (i, j), dimensioni n(1), n(2);
!   - `blk%sur(i,j)` è l'elemento d'area 2D che pesa i termini sorgente;
!   - Soot_Diffusive_Flux: T_gradient(2) e normal(2) al centro faccia; `int_l`
!     è la lunghezza dell'interfaccia 2D.
!
! Setup_Soot è il dispatcher logico (LL91/LIN/none) più il conteggio nsoot.
module Q2D_Mod_Soot
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: Setup_Soot, Soot_Source_Terms, Soot_Diffusive_Flux
  public :: i_CS  ! esposto per i residui soot dei termini sorgente assialsimmetrici/area-variation

  integer, private :: i_C2H2, i_CO2, i_CO, i_H2, i_O2, i_C6H6
  integer          :: i_CS  ! indice della massa di soot in 1:nsc (nsoot=1) oppure -1 (nsoot=2)

contains

  subroutine Setup_Soot()
    ! Dispatcher logico + ricerca degli indici in species_names da FLINT_Lib_Thermodynamic.
    ! nsoot impostato a 1 se C(gr) è presente nel set di specie (soot come "gas
    ! pesante"), 2 altrimenti (scalari soot separati).
    use Q2D_Config_Types_m,      only: obj_soot
    use Q2D_Global_m,            only: nsc, nsoot
    use FLINT_Lib_Thermodynamic, only: species_names
    implicit none
    integer :: i

    if (trim(obj_soot%model) == 'none') then
      obj_soot%enabled     = .false.
      obj_soot%description = 'Soot model disabled'
    else if (index(obj_soot%model, 'LL91') > 0) then
      obj_soot%enabled     = .true.
      obj_soot%LL91        = .true.
      obj_soot%LIN         = .false.
      obj_soot%description = 'Leung and Lindstedt (1991) model for soot formation and oxidation'
    else if (index(obj_soot%model, 'LIN') > 0) then
      obj_soot%enabled     = .true.
      obj_soot%LL91        = .false.
      obj_soot%LIN         = .true.
      obj_soot%description = 'Linear model for soot formation and oxidation'
    else
      obj_soot%error_message = 'Unknown soot model: '//trim(obj_soot%model)
    end if

    nsoot = 0
    if (.not. obj_soot%enabled) return

    ! Inizializza gli indici di specie
    i_O2   = -1
    i_C2H2 = -1
    i_CO2  = -1
    i_CO   = -1
    i_H2   = -1
    i_CS   = -1
    i_C6H6 = -1

    do i = 1, nsc
      select case (trim(adjustl(species_names(i))))
      case ('O2')
        i_O2 = i
      case ('C2H2')
        i_C2H2 = i
      case ('CO2')
        i_CO2 = i
      case ('CO')
        i_CO = i
      case ('H2')
        i_H2 = i
      case ('C(gr)')
        i_CS = i
      case ('C6H6')
        i_C6H6 = i
      end select
    end do

    if (i_CS == -1) then
      nsoot = 2
    else
      nsoot = 1
    end if
  end subroutine Setup_Soot


  subroutine Soot_Source_Terms ( domain )
    ! Driver con loop sui blocchi e gating MPI (elabora solo i blocchi locali).
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_soot
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Variabili locali
    integer :: b
    logical :: LL91, LIN

    LL91 = obj_soot%LL91
    LIN  = obj_soot%LIN

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle

      call Soot_Source_blk ( domain%blk(b)%P,    &
                             domain%blk(b)%R,    &
                             domain%blk(b)%sur,  &
                             domain%blk(b)%dim,  &
                             LL91, LIN )
    enddo

  contains

    subroutine Soot_Source_blk ( prim, res, sur, n, LL91, LIN )
      ! Kernel di blocco: prim/res di forma (1:nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc),
      ! sur di forma (1-gc:n(1)+gc, 1-gc:n(2)+gc); sur è l'elemento d'area 2D
      ! che pesa i termini sorgente.
      use Q2D_Global_m,            only: nsc, np, nc, nsoot, nprim, gc
      use Q2D_Lib_Soot,            only: Mc, LL91_soot_rates, LIN_soot_rates
      use FLINT_Lib_Thermodynamic, only: Co_rotot_Rtot, Wm_tab, f_tabT, h_tab
      implicit none
      integer,  intent(in)    :: n(2)
      logical,  intent(in)    :: LL91, LIN
      real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: prim
      real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: res
      real(R8), dimension(       1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: sur
      ! Variabili locali
      integer  :: i, j, s
      real(R8) :: rho, Rgas, T
      real(R8) :: Yc2h2, Yo2, C_c2h2, C_o2, Y_soot, Nvol
      real(R8) :: R1, R2, R3, R4, S_rhoYs, S_rhoN, omega(nsc), energy

      !$omp do collapse(2) private(rho, Rgas, T, Yc2h2, Yo2, Y_soot, Nvol, C_c2h2, C_o2) &
      !$omp private(R1, R2, R3, R4, S_rhoYs, S_rhoN, omega, energy, i, j, s)
      do j = 1, n(2)
      do i = 1, n(1)

        call Co_rotot_Rtot ( Prim(1:nsc, i, j), rho, Rgas )
        T = prim(np, i, j) / (rho * Rgas)

        ! --- dati specie / soot ---
        Yc2h2 = prim(i_C2H2, i, j) / rho
        Yo2   = prim(i_O2,   i, j) / rho
        if (nsoot == 1) then
          Y_soot = prim(i_CS, i, j) / rho      ! frazione di massa del soot
          Nvol   = prim(nc,   i, j)            ! densità numerica [#/m^3]
        else
          Y_soot = prim(nc,   i, j) / rho      ! frazione di massa del soot
          Nvol   = prim(nc+1, i, j)            ! densità numerica [#/m^3]
        end if

        ! concentrazioni molari (kmol/m^3)
        C_c2h2 = rho * Yc2h2 / Wm_tab(i_C2H2)
        C_o2   = rho * Yo2   / Wm_tab(i_O2)

        ! --- tassi di reazione (LL91 / LIN) ---
        if (LL91) then
          call LL91_soot_rates ( rho, T, C_c2h2, C_o2, Y_soot, Nvol, R1, R2, R3, R4 )
        elseif (LIN) then
          call LIN_soot_rates  ( rho, T, C_c2h2, C_o2, Y_soot, Nvol, R1, R2, R3, R4 )
        endif

        ! --- sorgenti per i residui ---
        ! massa soot: d(rho*Ys)/dt = Mc*(2R1 + R2 - R3)  [kg/m^3/s]
        S_rhoYs = Mc * ( 2.0d0*R1 + R2 - R3 )
        ! densità numerica di particelle: d(rhoN)/dt = R4  [#/m^3/s]
        S_rhoN  = R4

        ! Accumulo nei residui (sur = elemento d'area 2D)
        if (nsoot == 2) then
          res(nc,   i, j) = res(nc,   i, j) - S_rhoYs * sur(i, j)
          res(nc+1, i, j) = res(nc+1, i, j) - S_rhoN  * sur(i, j)
        else
          omega = 0.0_R8
          ! sorgenti di massa, scritte come produzione (negative per consumo)
          omega(i_C2H2) = -(R1+R2) * Wm_tab(i_C2H2)
          omega(i_O2)   = -0.5d0*R3 * Wm_tab(i_O2)
          omega(i_H2)   = +(R1+R2) * Wm_tab(i_H2)
          omega(i_CO)   = +R3      * Wm_tab(i_CO)
          omega(i_CS)   = -sum(omega)  ! impone la conservazione della massa
          ! aggiornamento dei residui delle specie
          res(1:nsc, i, j) = res(1:nsc, i, j) - omega(1:nsc) * sur(i, j)
          ! aggiornamento del residuo della densità numerica
          res(nc,    i, j) = res(nc,    i, j) - S_rhoN * sur(i, j)
          ! Termine sorgente di energia
          energy = 0.0_R8
          do s = 1, nsc
            energy = energy + omega(s) * f_tabT(T, s, h_tab)
          enddo
          res(np, i, j) = res(np, i, j) - energy * sur(i, j)
        endif

      enddo; enddo
    end subroutine Soot_Source_blk

  end subroutine Soot_Source_Terms


  subroutine Soot_Diffusive_Flux ( flux, prim, T_gradient, mu, rho, T, int_l, normal )
    ! Flusso termoforetico del soot attraverso una faccia: T_gradient(2) e
    ! normal(2) al centro faccia, int_l lunghezza dell'interfaccia 2D.
    use Q2D_Global_m, only: nc, nprim
    implicit none
    ! Ingressi
    real(R8), intent(in)    :: prim(nprim)
    real(R8), intent(in)    :: T_gradient(2)
    real(R8), intent(in)    :: mu, rho, T, int_l
    real(R8), intent(in)    :: normal(2)
    ! Uscite
    real(R8), intent(inout) :: flux(nprim)
    ! Costanti / parametri del modello
    real(R8), parameter :: CT = 0.55_R8        ! coefficiente termoforetico
    ! Variabili locali
    real(R8) :: Ys, Nvol, Thermo

    ! Estrazione degli scalari del soot
    Ys   = prim(i_CS) / rho     ! adimensionale
    Nvol = prim(nc)             ! [#/m^3]

    ! Flusso di massa termoforetico (la diffusione molecolare è già calcolata a monte)
    Thermo = -CT * (mu * Ys) * dot_product(T_gradient, normal) / T * int_l   ! [kg/s]
    flux(i_CS) = flux(i_CS) + Thermo

    ! Flusso equivalente su Nvol
    if (Ys > 1.0e-10_R8) flux(nc) = flux(i_CS) / (rho * Ys / Nvol)
  end subroutine Soot_Diffusive_Flux

end module Q2D_Mod_Soot
