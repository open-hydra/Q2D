module Q2D_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: Setup_Metrics

contains

  subroutine Setup_Metrics ( domain )
    use Q2D_Base_Types_m,             only: Q2D_vector_2D_type
    use Q2D_Advanced_Types_m,         only: Q2D_domain_type
    use Q2D_Global_m,                 only: model
    use Q2D_Config_Types_m,           only: obj_source
    use Q2D_Mod_MPI,                  only: is_local_block
    use Q2D_Lib_Metrics,              only: Check_Mesh_Type, Compute_Norm_Area, Compute_Metric_Tensor, &
                                             Compute_Yn, Yn_Connection, &
                                             BC_Extrapolate_Metrics, BC_Symmetry_Metrics, BC_Connect_Metrics, &
                                             BC_Connect_Metrics_Periodic
    use Q2D_Lib_Source_Axisymmetry,   only: Precompute_Axisymmetric_Geometry
    use Q2D_Lib_Source_Areavariation, only: Precompute_Area_Geometry
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Locali
    integer :: b, i, j
    integer :: Bm, Im, Jm, Fm, Bs, Is, Js, Fs, d11s, d12s, d21s, d22s
    type(Q2D_vector_2D_type) :: N1, N2, N3, N4

    call Check_Mesh_Type ( domain )

    ! area_variation: Cross_Section è un input geometrico nodale, letto una sola
    ! volta in Wrap_Setup su domain(1) e ristretto fine→coarse tramite
    ! Coarse_Cross_Section in Setup_Multigrid. Nessuna rilettura per livello qui:
    ! produrrebbe una re-indicizzazione errata per MGL>1.

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle  ! MPI: salta i blocchi non locali
      !$omp parallel default(none) shared(domain, b) private(i, j, N1, N2, N3, N4)
      !$omp do collapse(2)
      do j = 1, domain%blk(b)%dim(2)
        do i = 1, domain%blk(b)%dim(1)
          N1%c = domain%blk(b)%node(i-1, j-1)%c
          N2%c = domain%blk(b)%node(i-1, j  )%c
          N3%c = domain%blk(b)%node(i  , j-1)%c
          N4%c = domain%blk(b)%node(i  , j  )%c
          call Compute_Metric_Tensor( N1, N2, N3, N4, &
                                      domain%blk(b)%M(i,j), domain%blk(b)%dl(i,j), domain%blk(b)%sur(i,j) )
        end do
      end do
      !$omp end parallel

      call Compute_Norm_Area ( domain%blk(b) )
    end do

    ! Precalcolo della geometria dei termini sorgente dopo il popolamento di M.
    ! axisymmetric e area_variation sono mutuamente esclusivi (imposto in
    ! Register_Sources + Mod_Sources::Initialize_Sources).
    do b = 1, domain%nb
      if (obj_source%axisymmetric)   call Precompute_Axisymmetric_Geometry(domain%blk(b))
      ! Gate per blocco: un blocco privo del file INPUT/block<b>_area.dat mantiene
      ! area_variation_active=.false. (Read_Block_Q2D_Area lo salta) e NON deve avere
      ! area_geom_term calcolato da una Cross_Section non inizializzata; conserva
      ! invece il riempimento a zero (V_dot_G=0). Gestione robusta del caso
      ! multiblocco parziale.
      if (obj_source%area_variation .and. domain%blk(b)%area_variation_active) &
        call Precompute_Area_Geometry(domain%blk(b))
    end do

    ! Calcolo di yn condizionato a `model>1`: gli usi a valle di yn sono
    ! condizionati da `if (model==2)` in Lib_Diffusive/Lib_Source_*Visc, quindi
    ! nel caso laminare yn non serve. yn è comunque sempre allocato; eseguire
    ! Compute_Yn solo per `model>1` evita il relativo ciclo nel caso laminare.
    if (model>1) then
      call Compute_Yn ( domain%blk(:), domain%bc, domain%nb, domain%nbound )

      !$omp parallel default(none) shared(domain) private(i, Bm, Im, Jm, Fm, Bs, Is, Js)
      !$omp do schedule(dynamic)
      do i = 1, domain%nbound
        select case (domain%bc(i)%type)
          case(101, 201)
            ! Una 201 (periodica traslazionale) è una giunzione di connessione: il suo
            ! yn ghost è la distanza da parete del donatore (uno scalare, invariante per
            ! traslazione/rotazione). Nessuna guardia is_local_block(Bs): Compute_Yn qui
            ! sopra è globale (tutti i blocchi, senza guardia) e i nodi del donatore sono
            ! residenti per tutta la durata di Setup_Metrics, quindi blks%yn è valido
            ! anche per un donatore su altro rank.
            Bm = domain%bc(i)%b
            Im = domain%bc(i)%i
            Jm = domain%bc(i)%j
            Fm = domain%bc(i)%f
            Bs = domain%bc(i)%bs
            Is = domain%bc(i)%is
            Js = domain%bc(i)%js
            call Yn_Connection ( Im, Jm, Fm, domain%blk(Bm), Is, Js, domain%blk(Bs) )
        end select
      end do
      !$omp end parallel
    end if

    !$omp parallel default(none) shared(domain) &
    !$omp          private(i, Bm, Im, Jm, Fm, Bs, Is, Js, Fs, d11s, d12s, d21s, d22s)
    !$omp do schedule(dynamic)
    do i = 1, domain%nbound
      Bm = domain%bc(i)%b
      Im = domain%bc(i)%i
      Jm = domain%bc(i)%j
      Fm = domain%bc(i)%f
      select case (domain%bc(i)%type)
        case(101)
          Bs   = domain%bc(i)%bs
          Is   = domain%bc(i)%is
          Js   = domain%bc(i)%js
          Fs   = domain%bc(i)%fs
          d11s = domain%bc(i)%d11
          d12s = domain%bc(i)%d12
          d21s = domain%bc(i)%d21
          d22s = domain%bc(i)%d22
          call BC_Connect_Metrics ( Im, Jm, Fm, domain%blk(Bm), &
                                    Is, Js, Fs, domain%blk(Bs), d11s, d12s, d21s, d22s, &
                                    domain%bc(i)%Mg, domain%bc(i)%dlg, domain%bc(i)%surg )
        case(201)
          ! Una 201 (periodica traslazionale) è una connessione il cui donatore giace a
          ! un offset cartesiano costante; BC_Connect_Metrics_Periodic ricostruisce la
          ! metrica ghost dai nodi del donatore traslati nel riferimento locale (ESATTA,
          ! contro l'estrapolazione cubica che una 201 riceverebbe altrimenti dal case
          ! default). Le 101 restano su BC_Connect_Metrics.
          Bs   = domain%bc(i)%bs
          Is   = domain%bc(i)%is
          Js   = domain%bc(i)%js
          Fs   = domain%bc(i)%fs
          d11s = domain%bc(i)%d11
          d12s = domain%bc(i)%d12
          d21s = domain%bc(i)%d21
          d22s = domain%bc(i)%d22
          call BC_Connect_Metrics_Periodic ( Im, Jm, Fm, domain%blk(Bm), &
                                    Is, Js, Fs, domain%blk(Bs), d11s, d12s, d21s, d22s, &
                                    domain%bc(i)%Mg, domain%bc(i)%dlg, domain%bc(i)%surg )
        case(300)
          call BC_Symmetry_Metrics ( Im, Jm, Fm, domain%blk(Bm), &
                                     domain%bc(i)%Mg, domain%bc(i)%dlg, domain%bc(i)%surg )
        case default
          call BC_Extrapolate_Metrics ( Im, Jm, Fm, domain%blk(Bm), &
                                        domain%bc(i)%Mg, domain%bc(i)%dlg, domain%bc(i)%surg )
      end select
    end do
    !$omp end parallel

  end subroutine Setup_Metrics

end module Q2D_Mod_Metrics
