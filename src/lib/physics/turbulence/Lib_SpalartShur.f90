!>@brief Correzione di rotazione/curvatura Spalart-Shur (Spalart & Shur 1997).
!
! Riferimento: Spalart P.R., Shur M. (1997), "On the sensitization of turbulence
!              models to rotation and curvature", Aerospace Sci. Tech. 1(5):297-302.
!
! Avvertenza di applicabilità 2D: il modello è motivato per flussi 3D con
! rotazione/curvatura significative (canali rotanti di turbomacchine, combustori
! con swirl, flussi vorticosi). Nel 2D-piano Q2D la correzione di rotazione/curvatura
! è formalmente operativa (vorticità scalare omega_z ben definita; r* = S/Omega
! calcolabile) ma con effetto fisico ridotto rispetto a un canale rotante 3D
! (cfr. Khayrullina et al. 2017). Nel 2D-axi senza swirl (No-w) la correzione di
! rotazione è formalmente applicabile ma su flussi effettivamente non rotanti →
! uso ingegneristico, beneficio modesto.
!
! ============================================================================
! Giustificazione matematica formale di V·∇S_ij in Q2D
! ============================================================================
!
! L'operatore che compare in SA-RC è la derivata convettiva (materiale)
! (Spalart-Shur 1997 §2): D(S_ij)/Dt = ∂_t S_ij + V_k ∂_k S_ij. A regime
! D/Dt → V·∇ applicato alle COMPONENTI scalari di S_ij (NON la divergenza ∇·T,
! che in coordinate cilindriche porterebbe termini di Christoffel come (T_rr-T_θθ)/r).
!
! [Caso variazione d'area, 2D piano mediato in profondità]
!   Coordinate cartesiane (x, y) con z direzione virtuale top/bot. ∂_x e ∂_y sono
!   derivate parziali ordinarie, NESSUN termine di Christoffel. Le componenti di
!   S_ij nella base cartesiana sono:
!     S_xx = ∂_x V_x;  S_yy = ∂_y V_y;  S_xy = ½(∂_x V_y + ∂_y V_x)
!     S_zz = -V·G  (Definizione A: ∂w/∂z mediato in profondità da Leibniz +
!                   BC di impermeabilità alla parete virtuale).
!   V·∇S_ij = V_x ∂_x S_ij + V_y ∂_y S_ij  →  FD componente per componente OK.
!   Per S_zz: ∂_x S_zz = -∂_x(V·G) si catturerebbe passando a Strain_Tensor il
!   V_dot_G della cella vicina (non fatto: cfr. trade-off in Compute_RC_Terms).
!
! [Caso assialsimmetrico, terna ortonormale cilindrica, No-w]
!   Coordinate (z, r, θ) con base ortonormale (e_z, e_r, e_θ). Le uniche derivate
!   di base non nulle sono ∂_θ e_r = e_θ, ∂_θ e_θ = -e_r (le altre = 0).
!   Sotto le ipotesi Q2D No-w (V_θ = 0) + assialsimmetria (∂_θ = 0):
!     V·∇S_ij = V_z ∂_z S_ij + V_r ∂_r S_ij + V_θ * [termini di Christoffel]
!             = V_z ∂_z S_ij + V_r ∂_r S_ij  (V_θ = 0 annulla i termini di Christoffel)
!   Quindi anche in axi è una FD componente per componente delle COMPONENTI di
!   S_ij nella terna ortonormale. Le componenti sono:
!     S_zz = ∂_z V_z;  S_rr = ∂_r V_r;  S_zr = ½(∂_z V_r + ∂_r V_z)
!     S_θθ = vr/r (strain θθ cilindrico dalla continuità: ∇·V = ∂_z V_z + (1/r) ∂_r(r V_r)
!                                                               → S_θθ = V_r/r incomprimibile)
!   Per S_θθ alla cella vicina: usare vr_neighbor E r_ave_neighbor catturerebbe
!   formalmente
!     ∂_r S_θθ = ∂_r(vr/r) = (1/r) ∂_r vr - vr/r²
!   ENTRAMBI i termini sarebbero catturati dalla FD `(S_θθ_jp1 - S_θθ_jm1)/(2*Δr)`
!   con S_θθ valutato con (vr, r_ave) della cella vicina. Usare invece l'r_ave
!   centrale perde il termine -vr/r² nella direzione radiale (scelta attuale:
!   cfr. trade-off in Compute_RC_Terms).
!
! [Singolarità sull'asse]
!   r_ave → 0 invoca il ramo de l'Hopital di Strain_Tensor_Axi: S_θθ = ∂_r V_r.
!   Usando r_ave_neighbor cella per cella, il fallback de l'Hopital si
!   applicherebbe localmente a ciascuna cella vicina nella regione dell'asse.
!
! ============================================================================
! 2 procedure con riduzioni 2D Q2D:
!   - vel_gradient è Q2D_tensor_2D_type(2,2) (soli gradienti di velocità nel piano;
!     d/dz fuori piano strutturalmente assente);
!   - strain_gradient(3,3,2): si conservano gli indici TENSORIALI 1..3 (Sij 3×3),
!     si riduce la direzione SPAZIALE da 3 a 2 (eliminata la direzione k);
!   - estrapolazione ghost su 4 facce (i, j) — eliminate le 2 facce k;
!   - Strain_Tensor a 4 argomenti (Gradvel(2,2), vr, r_ave[, V_dot_G]) restituisce
!     3×3 con la Def A incorporata via dispatch del puntatore a procedura (Std/Axi);
!   - **Strain_Tensor sulle 5 celle (centrale + 4 vicine)**:
!     * vel_gradient della cella vicina: esteso ai ghost con estrapolazione costante
!       (Compute_Velocity_Gradient) → la FD spaziale cattura ∂_k S(1,1), S(1,2),
!       S(2,1), S(2,2). Per una cella vicina fuori intervallo (i=0 o ni+1, j=0 o
!       nj+1) l'estrapolazione costante dà Sij_neighbor = Sij_central: FD
!       effettivamente decentrata al bordo con fattore 1/2 (accuratezza O(Δ) al
!       bordo vs O(Δ²) all'interno, per una correzione SA-RC comunque sub-dominante).
!     * (vr, r_ave, V_dot_G): SOLO centrati in cella, passati identici alle 5
!       chiamate — pattern coerente con i kernel Q2D centrati in cella
!       (Source_Axi_Visc_Blk / SST_Blk / SA_Source_Blk: metriche centrate in
!       cella, MAI ai ghost). Conseguenza: S(3,3) costante nello stencil FD →
!       ∂_k S(3,3) = 0 (trade-off analizzato in Compute_RC_Terms).
!     **Nessun accesso a metriche ghost**: le metriche `axis_r_ave` (allocata
!     1-gc:ni+gc ma popolata solo 1..ni) e `area_geom_term_*` (shape 1:ni, senza
!     ghost) non sono MAI lette fuori dall'intervallo interno; si usano solo valori
!     già calcolati alla cella di riferimento (i,j) — zero accessi a dati ghost.
!   - triplo loop sugli INDICI TENSORIALI (1..3) conservato in rtilde;
!   - correzione di vorticità assoluta obj_rot%enabled attiva.
module Q2D_Lib_SpalartShur
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private

  public :: Compute_Velocity_Gradient, Compute_RC_Terms

contains

  subroutine Compute_Velocity_Gradient ( domain )
    ! Gradiente di velocità centrato in cella via differenze finite centrate +
    ! estrapolazione ai ghost. Forma Q2D: solo nu:nv, senza loop k né facce k;
    ! vel_gradient è Q2D_tensor_2D_type(2,2).
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nu, nv
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, n
    real(R8), dimension(2,2) :: vel_diff   ! 2×2 (riduzione 2D)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: salta i blocchi non locali

      !$omp do collapse(2) private(vel_diff)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        ! Gradiente di velocità centrato in cella (senza nw né direzione k)
        vel_diff(:,1) = ( domain%blk(b)%P(nu:nv,i+1,j) - domain%blk(b)%P(nu:nv,i-1,j) )/2d0
        vel_diff(:,2) = ( domain%blk(b)%P(nu:nv,i,j+1) - domain%blk(b)%P(nu:nv,i,j-1) )/2d0
        ! epsilon 1d-40 (evita gradienti patologicamente nulli
        ! nelle regioni laminari prima del matmul con la metrica).
        domain%blk(b)%vel_gradient(i,j)%c = matmul ( vel_diff + 1d-40, domain%blk(b)%m(i,j)%c )
      end do ; end do

      ! Estrapolazione nelle celle ghost (derivate seconde richieste da Compute_RC_Terms).
      ! Facce i (nessuna faccia k)
      !$omp do private(n)
      do j = 1, domain % blk(b) % dim(2)
        n = domain % blk(b) % dim(1) + 1
        domain % blk(b) % vel_gradient(0,j) % c = domain % blk(b) % vel_gradient(1,j)   % c
        domain % blk(b) % vel_gradient(n,j) % c = domain % blk(b) % vel_gradient(n-1,j) % c
      end do

      ! Facce j (nessuna faccia k)
      !$omp do private(n)
      do i = 1, domain % blk(b) % dim(1)
        n = domain % blk(b) % dim(2) + 1
        domain % blk(b) % vel_gradient(i,0) % c = domain % blk(b) % vel_gradient(i,1)   % c
        domain % blk(b) % vel_gradient(i,n) % c = domain % blk(b) % vel_gradient(i,n-1) % c
      end do

      ! Facce k: ASSENTI — No-w + 2D piano + 2D-axi (nessuna direzione fuori piano).

    end do

  end subroutine Compute_Velocity_Gradient


  subroutine Compute_RC_Terms ( domain )
    ! rstar = S/Omega, rtilde dalla derivata materiale di Sij (senza
    ! denominatore: diverso tra SA-RC e SST-RC, applicato al sito d'uso in
    ! SA_Source_Blk / SST_Blk).
    !
    ! Trattamento dell'operatore V·∇S_ij (cfr. giustificazione formale in testa al file):
    !   - Cartesiano piano (Std + variazione d'area): V·∇ = V_x ∂_x + V_y ∂_y
    !                                                 (niente Christoffel);
    !   - Cilindrico assialsimmetrico No-w (Axi):     V·∇ = V_z ∂_z + V_r ∂_r
    !                                                 (V_θ=0 + assialsimmetria annullano
    !                                                 i termini di Christoffel);
    !   FD componente per componente formalmente corretta in entrambi i casi.
    !
    ! Pattern centrato in cella coerente con Source_Axi_Visc_Blk, SST_Blk, SA_Source_Blk:
    !   - vel_gradient sulle 5 celle (centrale + 4 vicine) — estrapolazione ghost costante
    !     fornita da Compute_Velocity_Gradient sopra (vel_gradient(0,j)%c =
    !     vel_gradient(1,j)%c, ecc.). Cattura ∂_k S(1,1), S(1,2), S(2,1), S(2,2) via
    !     FD sul vel_gradient della cella vicina.
    !   - vr_ij, r_ave_ij, V_dot_G_ij: SOLO centrati in cella, passati identici alle
    !     5 chiamate a Strain_Tensor. Le metriche non sono estrapolate ai ghost
    !     (axis_r_ave allocata 1-gc:ni+gc ma popolata solo 1..ni; area_geom_term_*
    !     shape 1:ni senza ghost). Conseguenza: S(3,3) costante nello stencil FD →
    !     ∂_k S(3,3) = 0.
    !     **Trade-off accettato**: si perdono i contributi `(1/r_ave) ∂_k vr` (axi) e
    !     `-∂_k(V·G)` (area-var) al trasporto dello strain. Sono SUB-SUB-DOMINANTI
    !     rispetto al contributo dominante via S(1,1)/S(1,2)/S(2,2) nel piano,
    !     catturato dalla FD su vel_gradient. (In una formulazione 3D completa
    !     vel_gradient(3,3) conterrebbe direttamente ∂_z V_z; nel 2D Q2D, S(3,3) è
    !     un portato della riduzione in z, trattato come costante per il trasporto
    !     dello strain.)
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nu, nv
    use Q2D_Lib_Fluid,        only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    use Q2D_Config_Types_m,   only: obj_rot, obj_source
    use Q2D_Mod_MPI,          only: is_local_block

    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, ii, jj, kk, ni, nj
    real(R8), dimension(3,3,2) :: strain_gradient                 ! senza direzione spaziale 3
    real(R8), dimension(3,3)   :: Sij_ip1, Sij_im1, Sij_jp1, Sij_jm1
    real(R8), dimension(3,3)   :: Sij, Wij, strain_transport
    real(R8) :: omega(3), S, O, rtilde
    real(R8) :: r_ave, vr_ij, VG_ij, Oz_frame
    logical  :: area_var

    area_var = obj_source%area_variation
    ! Velocità di rotazione del riferimento memorizzata una volta (lettura condivisa
    ! nei loop OMP; nulla a riferimento rotante disabilitato -> il termine di frame
    ! più sotto si annulla).
    Oz_frame = 0d0
    if (obj_rot%enabled) Oz_frame = obj_rot%Omega_vec(3)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle  ! MPI: salta i blocchi non locali

      ni = domain % blk(b) % dim(1)
      nj = domain % blk(b) % dim(2)

      !$omp do collapse(2) private (i, j, ii, jj, kk, strain_gradient, Sij_ip1, Sij_im1, Sij_jp1, Sij_jm1), &
      !$omp private (Sij, Wij, strain_transport, omega, S, O, rtilde, r_ave, vr_ij, VG_ij)
      do j = 1, nj
      do i = 1, ni

        ! Quantità geometriche centrate in cella (le metriche non sono estrapolate ai
        ! ghost; pattern coerente con Source_Axi_Visc_Blk + SST_Blk + SA_Source_Blk).
        r_ave = 0.25d0 * ( domain%blk(b)%node(i  ,j  )%c(2) + domain%blk(b)%node(i-1,j  )%c(2) &
                         + domain%blk(b)%node(i  ,j-1)%c(2) + domain%blk(b)%node(i-1,j-1)%c(2) )
        vr_ij = domain%blk(b)%P(nv, i, j)

        if (area_var) then
          VG_ij = domain%blk(b)%P(nu,i,j)*domain%blk(b)%area_geom_term_x(i,j) &
                + domain%blk(b)%P(nv,i,j)*domain%blk(b)%area_geom_term_y(i,j)
        else
          VG_ij = 0d0
        end if

        ! ====================================================================
        ! Tensore di deformazione sulle 5 celle.
        ! - vel_gradient della cella vicina: esteso ai ghost (Compute_Velocity_Gradient
        !   sopra) → cattura ∂_k S(1,1), S(1,2), S(2,1), S(2,2) via FD spaziale.
        ! - (vr_ij, r_ave, VG_ij): centrati in cella, identici per le 5 chiamate
        !   → S(3,3) costante nello stencil FD; ∂_k S(3,3) = 0 (perdita sub-sub-dominante
        !   accettata — cfr. intestazione del file per l'analisi formale).
        ! ====================================================================
        Sij_ip1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i+1,j)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_im1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i-1,j)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_jp1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i,j+1)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij_jm1 = Strain_Tensor ( domain%blk(b)%vel_gradient(i,j-1)%c, vr_ij, r_ave, V_dot_G=VG_ij )
        Sij     = Strain_Tensor ( domain%blk(b)%vel_gradient(i  ,j  )%c, vr_ij, r_ave, V_dot_G=VG_ij )

        ! Gradiente del tensore di deformazione via differenze finite (senza direzione
        ! spaziale 3). La FD componente per componente su S_ij nella base ortonormale
        ! (cilindrica axi No-w o cartesiana piana area-var) cattura formalmente V·∇S_ij.
        strain_gradient(:,:,1) = ( Sij_ip1 - Sij_im1 ) * 0.5d0
        strain_gradient(:,:,2) = ( Sij_jp1 - Sij_jm1 ) * 0.5d0

        ! Trasformazione metrica indice tensoriale per indice tensoriale (3×3 conservato).
        ! La metrica centrata in cella M(i,j) è applicata al gradiente spaziale (jacobiano
        ! dallo spazio computazionale (xi, eta) al fisico (x, y) o (z, r)) — invariante
        ! tensoriale, applicata a ciascuna componente di Sij separatamente.
        do ii = 1, 3
        do jj = 1, 3
          strain_gradient(ii,jj,1:2) = matmul ( strain_gradient(ii,jj,1:2), domain%blk(b)%m(i,j)%c )
        end do
        end do

        ! Tensore di vorticità centrato in cella + correzione di vorticità assoluta obj_rot
        omega = Vorticity_Vector( domain%blk(b)%vel_gradient(i,j)%c )
        ! Vorticità assoluta:
        ! omega_abs = omega_rel + 2*Omega_frame; Vorticity_Vector restituisce vector(3)
        ! con omega(1)=omega(2)=0; asse di Omega_vec [0,0,1] stretto.
        if (obj_rot%enabled) omega = omega + 2.0_R8 * obj_rot%Omega_vec
        Wij = Vorticity_Tensor( omega )
        S = sqrt ( 2d0*sum (Sij**2) )
        O = sqrt ( sum(omega**2) )
        domain%blk(b)%rc_term1(i,j) = S/O

        ! Derivata materiale del tensore di deformazione V·∇S_ij (senza direzione
        ! spaziale 3). P(nu:nv,i,j) sono le COMPONENTI fisiche di velocità (V_x, V_y)
        ! o (V_z, V_r). V·∇ è la derivata direzionale sulle componenti scalari di S_ij;
        ! formalmente corretta per il cartesiano piano (niente Christoffel) e per il
        ! cilindrico assialsimmetrico No-w (V_θ=0 + ∂_θ=0 annullano i contributi di
        ! Christoffel).
        do ii = 1, 3
        do jj = 1, 3
          strain_transport(ii,jj) = dot_product ( domain%blk(b)%P(nu:nv,i,j), &
                                                  strain_gradient(ii,jj,1:2) )
        end do
        end do

        ! ── Riferimento rotante: termine di frame della parentesi di Spalart-Shur ──
        ! Le formulazioni TMR SA-RC / SST-RC (Spalart-Shur 1997; Smirnov-Menter 2009)
        ! prescrivono
        !   r_tilde ∝ 2 W_ik S_jk [ DS_ij/Dt + (eps_imn S_jn + eps_jmn S_in) Omega_m ]
        ! Il termine in eps è la rotazione degli assi principali dello strain indotta
        ! dal riferimento; NON è coperto dalla vorticità assoluta già usata in W_ik
        ! (quella è una comparsa separata e indipendente di Omega). Con
        ! Omega = (0,0,Oz) e la S Q2D diagonale a blocchi (S13 = S23 = 0) il termine
        ! si riduce nel piano a:
        !   T11 = -2 Oz S12,  T22 = +2 Oz S12,  T12 = Oz (S11 - S22),
        !   T33 = T13 = T23 = 0  (struttura a blocchi conservata).
        ! Verificato su un vortice analitico (riferimento fisso vs rotante): CON questo
        ! termine il numeratore discreto converge (ordine 2) al valore invariante di
        ! frame del riferimento fisso; SENZA, converge a un valore distorto di
        ! O(Omega), indipendentemente dalla griglia.
        strain_transport(1,1) = strain_transport(1,1) - 2d0*Oz_frame*Sij(1,2)
        strain_transport(2,2) = strain_transport(2,2) + 2d0*Oz_frame*Sij(1,2)
        strain_transport(1,2) = strain_transport(1,2) + Oz_frame*(Sij(1,1) - Sij(2,2))
        strain_transport(2,1) = strain_transport(1,2)

        ! Tripla somma sugli INDICI TENSORIALI 1..3 (non spaziali, da conservare).
        rtilde = 0d0
        do ii = 1, 3
        do jj = 1, 3
        do kk = 1, 3
          rtilde = rtilde + Wij(ii,kk)*Sij(jj,kk)*strain_transport(ii,jj)
        end do
        end do
        end do
        ! Nota: denominatore assente, applicato al sito d'uso in SA_Source_Blk / SST_Blk
        ! (diverso tra SA-RC e SST-RC).
        domain%blk(b)%rc_term2(i,j) = 2d0*rtilde

      end do ; end do
    end do

  end subroutine Compute_RC_Terms

end module Q2D_Lib_SpalartShur
