!> @brief Termini sorgente Single Rotating Frame (SRF) — 2D-planare e area-variation.
!>
!> Formulazione rigidamente 2D (principio No-w: nessuna componente di velocità
!> fuori piano, nessun residuo associato, nessuno swirl). Ambito:
!>
!>   * 2D planare (`obj_source%axisymmetric = .false.`, axis = [0,0,1]) — VALIDO.
!>   * Q2D area-variation (`obj_source%area_variation = .true.`, axis = [0,0,1]) — VALIDO.
!>   * 2D-axi (`obj_source%axisymmetric = .true.`) — FUORI AMBITO (l'unico asse
!>     fisicamente sensato è [1,0,0] e richiederebbe una componente di swirl u_phi
!>     assente dal vettore di stato, fuori ambito sotto No-w).
!>
!> Setup_RotatingFrame applica un doppio controllo a runtime:
!>   Controllo #1 → error stop se `obj_source%axisymmetric = .true.` (2D-axi fuori ambito).
!>   Controllo #2 → error stop se `axis ≠ [0, 0, 1]` (tolleranza 1e-12).
!>
!> Equazioni aggiunte al residuo:
!>
!>   Quantità di moto (per unità di superficie; sur è l'analogo 2D del volume di cella):
!>     S_mom = -2 rho (Omega x w)         [Coriolis]
!>           + rho omega^2 r_perp          [centrifuga]
!>     con `r_perp = r_vec - (r_vec · axis) axis`.
!>     Con axis = [0,0,1] e w(3) = 0:
!>       Cor(1) = +2 rho omega w(2)
!>       Cor(2) = -2 rho omega w(1)
!>       Cor(3) = 0  (identicamente; non rappresentabile, nessuna riga di residuo associata)
!>       Cen(3) = 0  (r_vec(3) = 0 nel piano)
!>
!>   Energia (per unità di superficie):
!>     S_E = w · (rho omega^2 r_perp) = w(1)*Cen(1) + w(2)*Cen(2)
!>     [Coriolis non compie lavoro — `w · (Omega x w) = 0` identicamente.]
!>
!> BC di parete rotante: `RF_Wall_Velocity(r_face)` produce `w_wall(2)` per faccia,
!> passato come argomento opzionale a `BC_Wall_Heat / BC_Wall_Temperature`. Le pareti
!> elencate in `obj_rot%stationary_face` (cassa/mozzo) mantengono `w_wall = 0`. La BC
!> di simmetria (Eulero) non riceve `w_wall`: la riflessione speculare è invariante
!> rispetto a una traslazione rigida del sistema con velocità normale nulla.
!>
!> BC di inflow rotante: `RF_Convert_BC_Inflow` converte `(T0, p0, alpha)` dal
!> sistema assoluto a quello relativo PRIMA del dispatch a `BC_Inflow_Outflow`.
!> `beta` è imposto a 0 (No-w: nessuna componente di velocità fuori piano).
module Q2D_Lib_RotatingFrame
  use iso_fortran_env,      only: I4 => int32, R8 => real64
  use Q2D_Config_Types_m,   only: obj_rot, obj_source

  implicit none
  private

  public :: Setup_RotatingFrame
  public :: RotatingFrame_Source_Terms
  public :: RF_Wall_Velocity
  public :: RF_Is_Stationary_Face
  public :: RF_Convert_BC_Inflow
  public :: RF_Face_Center
  public :: RF_Torque_Walls

contains


  ! ---------------------------------------------------------------------------
  !> @brief Finalizza la configurazione del sistema rotante.
  !>
  !> Include il doppio controllo che impone l'ambito
  !> 2D-planare / Q2D area-variation.
  subroutine Setup_RotatingFrame()
    implicit none
    integer  :: m
    real(R8) :: norm

    if ( obj_rot%model == "" .or. trim(obj_rot%model) == "none" ) then
      obj_rot%enabled = .false.
      return
    elseif ( obj_rot%model == "rigid-body" ) then
      if ( obj_rot%omega < 0.0_R8 ) then
        write(*,'(A)') '[ERROR] Rotational frame model "rigid-body" requires omega >= 0'
        stop
      end if
      obj_rot%enabled = .true.
    end if

    ! Legge asse e origine dalle stringhe del registry (str scalari).
    read(obj_rot%axis_str, *)   obj_rot%axis
    read(obj_rot%origin_str, *) obj_rot%origin

    ! Normalizza l'asse a vettore unitario.
    norm = sqrt( sum(obj_rot%axis**2) )
    if ( norm < 1.0e-14_R8 ) then
      write(*,'(A)') '[ERROR] Rotation axis has zero length'
      stop
    end if
    obj_rot%axis = obj_rot%axis / norm

    ! ─── Doppio controllo Q2D (No-w + ambito) ────────────────────────────────
    !
    ! Controllo #1: 2D-axi fuori ambito. In geometria assialsimmetrica l'unico
    ! asse di rotazione fisicamente sensato è [1,0,0] (asse di simmetria) e
    ! richiederebbe una componente di swirl u_phi, fuori ambito sotto No-w.
    if ( obj_source%axisymmetric ) then
      write(*,'(A)') '[ERROR] RotatingFrame OUT-OF-SCOPE in 2D-axi (axisymmetric) configuration.'
      write(*,'(A)') '        2D-axi rotation around axis-of-symmetry [1,0,0] requires swirl velocity'
      write(*,'(A)') '        component u_phi (= Prim(nw)), which is OUT-OF-SCOPE under No-w.'
      write(*,'(A)') '        2D-axi rotating frames are reserved for a follow-up activity.'
      stop
    end if
    !
    ! Controllo #2: 2D planare + Q2D area-variation richiedono axis = [0,0,1]
    ! (asse perpendicolare al piano). Coriolis/centrifuga restano nel piano,
    ! senza componenti fuori piano (Cor(3) ≡ 0, Cen(3) = 0).
    if ( any( abs(obj_rot%axis - [0d0, 0d0, 1d0]) > 1.0e-12_R8 ) ) then
      write(*,'(A)')           '[ERROR] In Q2D 2D-planar / Q2D area-variation, RotatingFrame axis must be [0, 0, 1].'
      write(*,'(A)')           '        Rationale: an axis perpendicular to the (x,y) plane is the only choice that'
      write(*,'(A)')           '        keeps Coriolis/Centrifugal in-plane without requiring out-of-plane swirl.'
      write(*,'(A)')           '        Q2D solves only in the (x,y) plane.'
      write(*,'(A,3ES12.4)')   '        Got axis = ', obj_rot%axis
      stop
    end if
    ! ─────────────────────────────────────────────────────────────────────────

    obj_rot%Omega_vec = obj_rot%omega * obj_rot%axis

    ! Compatibilità Q2D area-variation + RF: entrambi operano nel piano (la
    ! sezione variabile h(x,y) non interagisce con la cinematica del sistema
    ! rotante). Avviso informativo, non un error stop.
    if ( obj_source%area_variation ) then
      obj_rot%warning_message = 'RotatingFrame combined with Q2D area-variation: both in-plane (compatible).'
    end if

    ! Legge i patch di parete fissa dalle stringhe del registry.
    ! Il contatore `m` è un indice logico sui patch RF, non un indice spaziale.
    if ( obj_rot%n_stationary > 0 ) then
      allocate( obj_rot%stationary_face(2, obj_rot%n_stationary) )
      do m = 1, obj_rot%n_stationary
        read(obj_rot%stationary_face_str(m), *) &
          obj_rot%stationary_face(1,m), obj_rot%stationary_face(2,m)
      end do
    end if

  end subroutine Setup_RotatingFrame


  ! ---------------------------------------------------------------------------
  !> @brief Aggiunge i termini sorgente di Coriolis e centrifugo all'array dei
  !>        residui di ogni blocco in `domain`.
  subroutine RotatingFrame_Source_Terms ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer :: b

    do b = 1, domain%nb
      if ( .not. is_local_block(b) ) cycle
      call RF_Source_Blk( domain%blk(b)%P,    &
                          domain%blk(b)%R,    &
                          domain%blk(b)%node, &
                          domain%blk(b)%sur,  &
                          domain%blk(b)%dim   )
    end do

  end subroutine RotatingFrame_Source_Terms


  ! ---------------------------------------------------------------------------
  !> @brief Kernel a livello di blocco: accumula in Res i termini sorgente del sistema rotante.
  !>
  !> Formulazione 2D:
  !>   * centro cella come media dei 4 nodi d'angolo;
  !>   * nessuna riga di residuo fuori piano (Cor(3) ≡ 0, Cen(3) = 0 con axis = [0,0,1] stretto);
  !>   * `sur` (metrica di superficie 2D, analogo del volume di cella) pesa i termini sorgente.
  !>
  !> Convenzione di segno: V dU/dt = -R(U), R = CONV - DIFF - SOURCE  ⇒  R -= S * sur.
  subroutine RF_Source_Blk ( Prim, Res, node, sur, n )
    use Q2D_Base_Types_m, only: Q2D_vector_2D_type
    use Q2D_Global_m,     only: nprim, nsc, nu, nv, np, gc

    implicit none
    integer,  intent(in) :: n(2)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)    :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(inout) :: Res
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)),  intent(in) :: node
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc),      intent(in) :: sur
    ! Variabili locali
    integer  :: i, j
    real(R8) :: rho, w(2)
    real(R8) :: xc, yc, r_vec(2)
    real(R8) :: Cor(2), Cen(2)
    real(R8) :: omega_z, ox(2), omega2

    ! Copia i dati di modulo in scalari locali per evitare accessi indiretti ripetuti.
    ! Con axis = [0,0,1] stretto (Controllo #2 di Setup_RotatingFrame): Omega_vec = (0,0,omega_z),
    ! origine proiettata nel piano (componente z scartata: r_vec è 2D).
    omega_z = obj_rot%Omega_vec(3)
    ox      = obj_rot%origin(1:2)
    omega2  = obj_rot%omega**2

    !$omp do collapse(2) &
    !$omp private( i, j, rho, w, xc, yc, r_vec, Cor, Cen )
    do j = 1, n(2)
    do i = 1, n(1)

      rho = sum( Prim(1:nsc, i, j) )
      w   = Prim(nu:nv, i, j)            ! velocità relativa nel sistema rotante (2 componenti)

      ! ---- Centro cella (media dei 4 nodi d'angolo) ----
      ! La cella (i,j) occupa il quadrilatero i-1..i, j-1..j nelle coordinate dei nodi.
      xc = 0.25_R8 * ( node(i-1,j-1)%c(1) + node(i,j-1)%c(1) &
                     + node(i-1,j  )%c(1) + node(i,j  )%c(1) )
      yc = 0.25_R8 * ( node(i-1,j-1)%c(2) + node(i,j-1)%c(2) &
                     + node(i-1,j  )%c(2) + node(i,j  )%c(2) )

      ! ---- Vettore posizione dall'origine dell'asse (nel piano) ----
      r_vec = [ xc - ox(1), yc - ox(2) ]
      ! r_vec(3) = 0 implicito; r_proj = dot(r_vec, [0,0,1]) = 0; r_perp = r_vec.

      ! ---- Forza di Coriolis per unità di superficie: -2 rho (Omega x w) --
      ! Con Omega = (0,0,omega_z), w(3) = 0:
      !   Cor(1) = -2 rho ( Omega(2) w(3) - Omega(3) w(2) ) = +2 rho omega_z w(2)
      !   Cor(2) = -2 rho ( Omega(3) w(1) - Omega(1) w(3) ) = -2 rho omega_z w(1)
      !   Cor(3) ≡ 0 (non rappresentabile nel piano; nessuna riga di residuo).
      Cor(1) =  2.0_R8 * rho * omega_z * w(2)
      Cor(2) = -2.0_R8 * rho * omega_z * w(1)

      ! ---- Forza centrifuga per unità di superficie: rho omega^2 r_perp --
      !   Cen(3) = 0 con r_vec(3) = 0.
      Cen(1) = rho * omega2 * r_vec(1)
      Cen(2) = rho * omega2 * r_vec(2)

      ! ---- Accumulo nei residui: R -= sorgente * sur ----
      Res(nu, i, j) = Res(nu, i, j) - ( Cor(1) + Cen(1) ) * sur(i,j)
      Res(nv, i, j) = Res(nv, i, j) - ( Cor(2) + Cen(2) ) * sur(i,j)
      ! Nessun residuo fuori piano: Cor(3) + Cen(3) ≡ 0 identicamente (No-w).

      ! Energia: solo la centrifuga compie lavoro su w; w · Coriolis ≡ 0 identicamente.
      Res(np, i, j) = Res(np, i, j) - dot_product( w, Cen ) * sur(i,j)

    end do ; end do
    !$omp end do

  end subroutine RF_Source_Blk


  ! ---------------------------------------------------------------------------
  !> @brief Velocità di una parete fissa nel laboratorio vista dal sistema rotante.
  !>
  !> Con axis = [0,0,1] stretto:
  !>   w_wall(1) = +omega_z * (r_face(2) - origin(2))
  !>   w_wall(2) = -omega_z * (r_face(1) - origin(1))
  !>   w_wall(3) ≡ 0 (No-w: non rappresentata).
  !>
  !> @param[in]  r_face  Posizione cartesiana 2D del centro faccia [m].
  !> @return     w_wall  Velocità relativa della parete fissa nel laboratorio [m/s].
  function RF_Wall_Velocity(r_face) result(w_wall)
    implicit none
    real(R8), intent(in) :: r_face(2)
    real(R8)             :: w_wall(2)
    real(R8) :: r(2), omega_z

    r       = r_face - obj_rot%origin(1:2)
    omega_z = obj_rot%Omega_vec(3)
    w_wall(1) =  omega_z * r(2)
    w_wall(2) = -omega_z * r(1)
  end function RF_Wall_Velocity


  ! ---------------------------------------------------------------------------
  !> @brief Centro faccia 2D come media dei 2 nodi estremi (in 2D la faccia è uno spigolo).
  !>
  !> Q2D ha 4 facce per cella; Fm ∈ {1,2,3,4}. Fm = 5 o 6 (facce k) produce un
  !> error stop esplicito (nessun fallback silenzioso).
  subroutine RF_Face_Center(node, Im, Jm, Fm, r_face)
    use Q2D_Base_Types_m, only: Q2D_vector_2D_type
    implicit none
    type(Q2D_vector_2D_type), intent(in)  :: node(0:,0:)
    integer,                  intent(in)  :: Im, Jm, Fm
    real(R8),                 intent(out) :: r_face(2)
    integer :: Fi, Fj, d

    select case (Fm)
    case(1,2)
      Fi = Im - mod(Fm,2) ; Fj = Jm
    case(3,4)
      Fi = Im             ; Fj = Jm - mod(Fm,2)
    case(5,6)
      write(*,'(A,I0,A)') '[ERROR] RF_Face_Center: Fm=', Fm, ' (k-face) invalid in Q2D 2D-planar.'
      stop
    case default
      write(*,'(A,I0,A)') '[ERROR] RF_Face_Center: Fm=', Fm, ' out of range [1..4].'
      stop
    end select

    do d = 1, 2
      select case (Fm)
      case(1,2)   ! faccia i: spigolo tra (Fi, Fj-1) e (Fi, Fj)
        r_face(d) = 0.5_R8 * ( node(Fi,Fj-1)%c(d) + node(Fi,Fj)%c(d) )
      case(3,4)   ! faccia j: spigolo tra (Fi-1, Fj) e (Fi, Fj)
        r_face(d) = 0.5_R8 * ( node(Fi-1,Fj)%c(d) + node(Fi,Fj)%c(d) )
      end select
    end do

  end subroutine RF_Face_Center


  ! ---------------------------------------------------------------------------
  !> @brief Converte la BC di inflow `(T0, p0, alpha, Mach_abs) → (T0_rel, p0_rel, alpha_rel, Mach_rel)`.
  !>
  !> Formulazione 2D:
  !>   * nessun angolo di elevazione fuori piano: `beta` imposto a 0;
  !>   * v_frame(3) ≡ 0 con Omega = (0,0,omega_z), r_face(3) = 0.
  !>
  !> Le condizioni statiche sono invarianti per cambio di sistema; le condizioni
  !> totali e l'angolo vanno convertiti. La pressione statica di outflow (pAmb) è
  !> indipendente dal sistema di riferimento — il dispatcher chiamante non invoca
  !> questa routine per outflow puro.
  subroutine RF_Convert_BC_Inflow ( Blk, Fm, Im, Jm, BC_Mach_abs, gamma, Rgas, &
                                    T0_or_Tstat, p0_or_pstat, alpha, Mach_rel )
    use Q2D_Advanced_Types_m, only: Q2D_block_type
    implicit none
    type(Q2D_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Fm, Im, Jm
    real(R8), intent(in)    :: BC_Mach_abs, gamma, Rgas
    real(R8), intent(inout) :: T0_or_Tstat, p0_or_pstat, alpha
    real(R8), intent(out)   :: Mach_rel
    ! Variabili locali
    integer  :: Face_i, Face_j
    real(R8) :: v_frame(2), r(2), r_fc(2)
    real(R8) :: T_stat, p_stat, a_stat, u_abs(2), w_rel(2), w_norm, cp
    real(R8) :: omega_z

    ! Indice di faccia (convenzione Face_Index: Fm=1,2 → faccia i; Fm=3,4 → faccia j).
    select case (Fm)
    case(1,2) ; Face_i = Im - mod(Fm,2) ; Face_j = Jm
    case(3,4) ; Face_i = Im             ; Face_j = Jm - mod(Fm,2)
    case default
      write(*,'(A,I0,A)') '[ERROR] RF_Convert_BC_Inflow: Fm=', Fm, ' invalid in Q2D 2D-planar (expected 1..4).'
      stop
    end select

    ! Centro faccia (media di 2 nodi, 2D).
    select case (Fm)
    case(1,2)   ! faccia i: spigolo tra (Face_i, Face_j-1) e (Face_i, Face_j)
      r_fc(1) = 0.5_R8 * ( Blk%node(Face_i,Face_j-1)%c(1) + Blk%node(Face_i,Face_j)%c(1) )
      r_fc(2) = 0.5_R8 * ( Blk%node(Face_i,Face_j-1)%c(2) + Blk%node(Face_i,Face_j)%c(2) )
    case(3,4)   ! faccia j: spigolo tra (Face_i-1, Face_j) e (Face_i, Face_j)
      r_fc(1) = 0.5_R8 * ( Blk%node(Face_i-1,Face_j)%c(1) + Blk%node(Face_i,Face_j)%c(1) )
      r_fc(2) = 0.5_R8 * ( Blk%node(Face_i-1,Face_j)%c(2) + Blk%node(Face_i,Face_j)%c(2) )
    end select

    ! Velocità di trascinamento al centro faccia: v_frame = Omega x (r_fc - origin), Omega = (0,0,omega_z).
    !   v_frame(1) = -omega_z * r(2)
    !   v_frame(2) =  omega_z * r(1)
    !   v_frame(3) = 0  (non rappresentata).
    r          = r_fc - obj_rot%origin(1:2)
    omega_z    = obj_rot%Omega_vec(3)
    v_frame(1) = -omega_z * r(2)
    v_frame(2) =  omega_z * r(1)

    ! Condizioni statiche (invarianti per cambio di sistema).
    if ( BC_Mach_abs > 1.0_R8 ) then
      ! Supersonico: in ingresso T0_or_Tstat = T statica, p0_or_pstat = p statica.
      T_stat = T0_or_Tstat
      p_stat = p0_or_pstat
    else
      ! Subsonico: in ingresso T0_or_Tstat = T0 assoluta, p0_or_pstat = p0 assoluta.
      T_stat = T0_or_Tstat / ( 1.0_R8 + 0.5_R8*(gamma-1.0_R8)*BC_Mach_abs**2 )
      p_stat = p0_or_pstat * ( T_stat / T0_or_Tstat )**( gamma/(gamma-1.0_R8) )
    end if

    ! Velocità del suono (invariante per cambio di sistema).
    a_stat = sqrt( gamma * Rgas * T_stat )

    ! Velocità assoluta (nessuna componente z, beta imposto a 0).
    !   u = M*a * [cos(alpha), sin(alpha)]  con beta = 0 ⇒ cos(beta) = 1, sin(beta) = 0.
    u_abs(1) = BC_Mach_abs * a_stat * cos(alpha)
    u_abs(2) = BC_Mach_abs * a_stat * sin(alpha)

    ! Velocità relativa: w = u_abs - v_frame.
    w_rel  = u_abs - v_frame
    w_norm = norm2( w_rel )

    ! Mach relativo e nuovo angolo di flusso (beta resta 0 per costruzione in 2D).
    Mach_rel = w_norm / a_stat
    alpha    = atan2( w_rel(2), w_rel(1) )

    if ( BC_Mach_abs > 1.0_R8 ) then
      ! T statica e p statica invariate (invarianti); Mach_rel e alpha già aggiornati sopra.
    else
      ! Condizioni totali relative dalla velocità relativa.
      cp          = gamma * Rgas / ( gamma - 1.0_R8 )
      T0_or_Tstat = T_stat + 0.5_R8 * w_norm**2 / cp
      p0_or_pstat = p_stat * ( T0_or_Tstat / T_stat )**( gamma/(gamma-1.0_R8) )
    end if

  end subroutine RF_Convert_BC_Inflow


  ! ---------------------------------------------------------------------------
  !> @brief `.true.` se il patch di parete (blocco b, direzione di faccia f) è fisso nel laboratorio.
  !>
  !> Funzione puramente logica, indipendente dal regime.
  logical function RF_Is_Stationary_Face(b, f)
    implicit none
    integer, intent(in) :: b, f
    integer :: m

    RF_Is_Stationary_Face = .false.
    if ( .not. allocated(obj_rot%stationary_face) ) return
    do m = 1, size(obj_rot%stationary_face, 2)
      if ( obj_rot%stationary_face(1,m) == b .and. &
           obj_rot%stationary_face(2,m) == f ) then
        RF_Is_Stationary_Face = .true.
        return
      end if
    end do

  end function RF_Is_Stationary_Face


  ! ---------------------------------------------------------------------------
  !> @brief Coppia di pressione sulle pareti rotanti; scrive OUTPUT/torque.dat.
  !>
  !> In 2D-planare con axis=[0,0,1]:
  !>   torque_vec(1) = r_rel(2)*F_wall(3) - r_rel(3)*F_wall(2) = 0  (F_wall(3)=r_rel(3)=0)
  !>   torque_vec(2) = r_rel(3)*F_wall(1) - r_rel(1)*F_wall(3) = 0
  !>   torque_vec(3) = r_rel(1)*F_wall(2) - r_rel(2)*F_wall(1)  = T_axial nel piano
  !>   T_axial = axis · torque_vec = torque_vec(3)
  !>
  !> Q2D usa `int_l` (lunghezza dell'interfaccia) come analogo 2D dell'area di
  !> faccia (cf. `Q2D_f_metrics_type` in `Base_Types_m.f90`). MPI: la coppia
  !> parziale scalare è ridotta su root via `mpi_reduce_sum_r8` (no-op a rank
  !> singolo).
  subroutine RF_Torque_Walls(domain, iter, time)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: np
    use Q2D_Lib_Ghost,        only: Is_Wall
    use Q2D_Mod_MPI,          only: is_local_block, mpi_is_root, mpi_reduce_sum_r8
    implicit none
    type(Q2D_domain_type), intent(in) :: domain
    integer,  intent(in) :: iter
    real(R8), intent(in) :: time
    ! Variabili locali
    integer,  save :: unit_torque = -1
    logical,  save :: first_call  = .true.
    integer  :: i, Bm, Im, Jm, Fm, Dir, Fi, Fj
    real(R8) :: r_face(2), r_rel(2), F_wall(2)
    real(R8) :: T_axial, p_face, sign_fm
    real(R8) :: N_face(2), int_l_face

    T_axial = 0.0_R8

    do i = 1, domain%nbound

      if ( .not. Is_Wall(domain%bc(i)%type) ) cycle

      Bm = domain%bc(i)%b
      if ( .not. is_local_block(Bm) ) cycle

      Fm = domain%bc(i)%f

      ! Salta le facce fisse nel laboratorio (cassa / mozzo).
      if ( RF_Is_Stationary_Face(Bm, Fm) ) cycle

      Im = domain%bc(i)%i
      Jm = domain%bc(i)%j

      ! Indice di direzione: 1 per le facce i (Fm=1,2), 2 per le facce j (3,4).
      Dir = (Fm + 1) / 2

      select case (Fm)
      case (1,2) ; Fi = Im - mod(Fm,2) ; Fj = Jm
      case (3,4) ; Fi = Im             ; Fj = Jm - mod(Fm,2)
      end select

      ! dir(Dir)%f%N (normale unitaria) + dir(Dir)%f%int_l (lunghezza della faccia 2D).
      N_face     = domain%blk(Bm)%dir(Dir)%f(Fi,Fj)%N
      int_l_face = domain%blk(Bm)%dir(Dir)%f(Fi,Fj)%int_l

      ! Segno uscente dal fluido: -1 per le facce min (Fm dispari), +1 per le facce max (Fm pari).
      sign_fm = real(1 - 2 * mod(Fm, 2), R8)

      ! Pressione nella cella interna adiacente alla parete (estrapolazione di ordine zero).
      p_face = domain%blk(Bm)%P(np, Im, Jm)

      call RF_Face_Center(domain%blk(Bm)%node, Im, Jm, Fm, r_face)

      ! Forza di pressione del fluido sulla parete: F_wall = p * (normale uscente) * int_l.
      F_wall = p_face * sign_fm * N_face * int_l_face

      ! Coppia nel piano rispetto all'origine dell'asse di rotazione (scalare, axis=[0,0,1]):
      r_rel   = r_face - obj_rot%origin(1:2)
      T_axial = T_axial + ( r_rel(1)*F_wall(2) - r_rel(2)*F_wall(1) )

    end do

    ! Riduce su root le coppie parziali di tutti i rank (collettiva bloccante;
    ! no-op a rank singolo).
    call mpi_reduce_sum_r8(T_axial)

    if ( mpi_is_root ) then
      if ( first_call ) then
        open(newunit=unit_torque, file='OUTPUT/torque.dat', &
             status='replace', form='formatted')
        write(unit_torque, '(A)') &
          '#        iter              time   torque_axial (pressure, non-dim)' // &
          '         power (T*omega, non-dim)'
        first_call = .false.
      end if

      write(unit_torque, '(I9,3E20.10)') iter, time, T_axial, T_axial * obj_rot%omega
      flush(unit_torque)
    end if

  end subroutine RF_Torque_Walls


end module Q2D_Lib_RotatingFrame
