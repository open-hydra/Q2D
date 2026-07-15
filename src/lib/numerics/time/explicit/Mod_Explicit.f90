module Q2D_Mod_Explicit
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Explicit_Step

contains

  subroutine Explicit_Step ( domain, External_Function )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_sim_param, obj_time_scheme, obj_multigrid, &
                                    obj_chemistry, obj_io, obj_soot, obj_rot
    use Q2D_Global_m,         only: nres, model, SOURCE_TERM
    use Q2D_Mod_dt,           only: Set_Global_dt, Compute_dt
    use Q2D_Lib_Ghost,        only: Fill_Ghost_Cell, Ghost_Wall_Extrapolation
    use Q2D_Mod_Fluxes,       only: Fluxes
    use Q2D_Mod_BC_Fluxes,    only: BC_Fluxes
    use Q2D_Mod_Newstate,     only: RK_Newstate
    use Q2D_Mod_Sources,      only: Compute_Sources
    use Q2D_Mod_Diagnostic,   only: Compute_Residual
    use Q2D_Lib_RANS,         only: RANS_Source_Terms
    use Q2D_Lib_Chemistry,    only: Chemistry_Newstate, Chemistry_Equilibrium
    use Q2D_Lib_RotatingFrame,only: RotatingFrame_Source_Terms
    use Q2D_Mod_Soot,         only: Soot_Source_Terms
    use Q2D_Mod_MPI,          only: is_local_block, mpi_reduce_sum_r8_array, &
                                    mpi_is_root, mpi_bcast_logical, mpi_bcast_integer
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    external :: External_Function
    ! Variabili locali
    logical  :: endsim, iosim, endmg
    integer  :: i_rk, i_strang, b, level
    real(R8) :: average(nres)

    level = obj_multigrid%MG_level
    obj_multigrid%change_MG = .false.
    domain(level) % iter = domain(level) % iter + 1
    obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
    obj_sim_param%iter_general   = obj_sim_param%iter_general + 1

    ! Calcolo del passo temporale. I flag HYDRA_* sono flag di ACCOPPIAMENTO:
    ! nelle build standalone di Q2D restano al default `.false.`; l'orchestratore
    ! di alto livello ATLAS può impostarli a .true. per l'accoppiamento HYDRA.
    if (obj_sim_param%HYDRA_time_accurate) then
      call Set_Global_dt ( domain(level) )
      domain(level) % time = domain(level) % time + domain(level) % dtglobal
    else
      domain(level) % dtglobal = 1d5
      call Compute_dt ( domain(level), obj_time_scheme%cfl, obj_time_scheme%vnn, obj_time_scheme%rampa_cfl_iter )
      if ( obj_time_scheme%time_accurate ) then
        call Set_Global_dt ( domain(level) )
        domain(level) % time = domain(level) % time + domain(level) % dtglobal
      endif
    endif

    !$omp parallel
    strang: do i_strang = 1, obj_chemistry%N_strang ! splitting di Strang: L = L_rk(dt/2) * L_chem(dt) * L_rk(dt/2)

      if (obj_chemistry%imodel == 2) then
        do b = 1, domain(level) % nb
          if (.not. is_local_block(b)) cycle
          call Chemistry_Equilibrium ( domain(level)%blk(b)%P, domain(level)%blk(b)%dim )
        enddo
      end if

      call Copy_State ( domain(level) )

      rk: do i_rk = 1, obj_time_scheme%n_RK

        call Fill_Ghost_Cell ( domain(level) )               ! riempimento delle celle ghost
        call Fluxes ( domain(level) )                        ! flussi convettivi e diffusivi
        call BC_Fluxes ( domain(level) )                     ! flussi al contorno

        call External_Function ( domain(level) )             ! funzione esterna (es. termini sorgente)

        ! L'estrapolazione di parete nelle celle ghost DEVE precedere ogni termine
        ! sorgente: i termini sorgente geometrici (viscosi da variazione d'area) e
        ! quelli RANS leggono le celle ghost di parete tramite lo stencil dei
        ! gradienti Prim(i±1,j±1). Il gate include SOURCE_TERM oltre a model==2
        ! (RANS ⇔ model==2) per coprire anche il caso laminare con sorgenti da
        ! variazione d'area.
        if (model==2 .or. SOURCE_TERM) call Ghost_Wall_Extrapolation ( domain(level) )

        ! Dispatch dei termini sorgente geometrici (assialsimmetrico /
        ! area_variation). Gate SOURCE_TERM: inizializzato da Initialize_Sources,
        ! invocata in Read_Inifile / Wrap_Setup; default .false. → nessuna
        ! operazione.
        if (SOURCE_TERM) call Compute_Sources ( domain(level) )

        ! Termini sorgente di Coriolis e centrifughi in sistema di riferimento
        ! rotante (gate obj_rot%enabled — Setup_RotatingFrame applica un doppio
        ! controllo: asse strettamente [0,0,1] e axisymmetric=.false.).
        if (obj_rot%enabled) call RotatingFrame_Source_Terms ( domain(level) )

        if (model==2) then
          call RANS_Source_Terms ( domain(level) )           ! RHS delle equazioni di turbolenza
        endif

        if (obj_soot%enabled) call Soot_Source_Terms ( domain(level) )

        call RK_Newstate ( domain(level), i_rk )             ! aggiornamento dello stato

      enddo rk

      ! Chimica, solo se i_strang=1 e N_strang=2
      if ( i_strang + obj_chemistry%N_strang == 3 ) then
        do b = 1, domain(level) % nb
          if (.not. is_local_block(b)) cycle
          if (domain(level)%blk(b)%no_chem) cycle
          call Chemistry_Newstate ( domain(level)%blk(b)%P, domain(level)%blk(b)%dtlocal, &
                                    domain(level)%blk(b)%dim )
        enddo
      endif

    enddo strang
    !$omp end parallel

    ! Salvataggio dei residui
    if ( (obj_sim_param%iter_from_call == 1) .or. (mod (domain(level) % iter, obj_io%res_diter) == 0) ) then
      obj_sim_param%residuotot = 0.d0
      do b = 1, domain(level) % nb
        if (.not. is_local_block(b)) cycle
        call Compute_Residual ( new=domain(level)%blk(b)%P, &
                                old=domain(level)%blk(b)%PO, &
                                dt=domain(level)%blk(b)%dtlocal, &
                                n=domain(level)%blk(b)%dim, &
                                average=average, &
                                total=obj_sim_param%residuotot )
      enddo
      call mpi_reduce_sum_r8_array(obj_sim_param%residuotot, nres)
      if (mpi_is_root) obj_sim_param%residuotot = sqrt ( obj_sim_param%residuotot ) ! norma L2 della derivata temporale
    endif

    ! Determina i flag di controllo della simulazione sul rank root, poi li
    ! propaga in broadcast a tutti i rank. endsim/endmg sono variabili locali:
    ! a ogni chiamata ne viene assegnata una sola (in base a level), ma entrambe
    ! sono lette più sotto — vanno azzerate qui, altrimenti quella non assegnata
    ! contiene valori indefiniti dello stack (senza azzeramento, i livelli
    ! multigrid coarse possono terminare l'intera simulazione dopo 1 iterazione,
    ! come osservato con build gfortran di tipo TESTING).
    endsim = .false.
    endmg  = .false.
    if (mpi_is_root) then
      if (level == 1) then
        endsim = ( obj_sim_param%iter_from_call >= domain(1) % itermax ) &
            .or. ( obj_sim_param%residuotot(1) <= obj_sim_param%res_threshold ) &
            .or. ( domain(1) % time >= obj_sim_param%time_threshold )
      else
        endmg = ( obj_sim_param%iter_from_call >= domain(level) % itermax ) &
           .or. ( obj_sim_param%residuotot(1) <= obj_sim_param%res_threshold ) &
           .or. ( domain(level) % time >= obj_sim_param%time_threshold )
      endif

      iosim  = ( mod (domain(level) % iter, obj_io%sol_diter) == 0) &
          .or. ( domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime )

      if ( endsim ) then
        obj_sim_param%TODO = 3
      elseif ( iosim .or. endmg ) then
        obj_sim_param%TODO = 2
        if (endmg .and. obj_multigrid%MGL > 1) then
          obj_multigrid%change_MG = .true.
        endif
      else
        obj_sim_param%TODO = 1
      endif
    end if

    ! Broadcast del controllo della simulazione dal root a tutti i rank
    call mpi_bcast_integer(obj_sim_param%TODO)
    call mpi_bcast_logical(obj_multigrid%change_MG)

  end subroutine Explicit_Step


  subroutine Copy_State ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Variabili locali
    integer :: i, j, b

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(2)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)

        domain % blk(b) % PO (:,i,j) = domain % blk(b) % P (:,i,j)

      enddo; enddo
      !$omp end do
    enddo

  end subroutine Copy_State

end module Q2D_Mod_Explicit
