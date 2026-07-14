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
    ! Local
    integer :: b, i, j
    integer :: Bm, Im, Jm, Fm, Bs, Is, Js, Fs, d11s, d12s, d21s, d22s
    type(Q2D_vector_2D_type) :: N1, N2, N3, N4

    call Check_Mesh_Type ( domain )

    ! area_variation: Cross_Section is a nodal geometric input, read once in
    ! Wrap_Setup on domain(1) and restricted fine→coarse via Coarse_Cross_Section
    ! in Setup_Multigrid. No per-level re-read here (it caused scrambled
    ! re-indexing for MGL>1).

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks
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

    ! Precompute the source geometry after M is populated. axisymmetric and
    ! area_variation are mutually exclusive (enforced in Register_Sources +
    ! Mod_Sources::Initialize_Sources).
    do b = 1, domain%nb
      if (obj_source%axisymmetric)   call Precompute_Axisymmetric_Geometry(domain%blk(b))
      ! Per-block gate (mirror legacy Wrap_Setup): a block whose INPUT/block<b>_area.dat
      ! is absent keeps area_variation_active=.false. (Read_Block_Q2D_Area cycles) and
      ! must NOT have its area_geom_term computed from an uninitialised Cross_Section;
      ! it keeps the universal zero-fill (V_dot_G=0) instead. Graceful partial-multiblock.
      if (obj_source%area_variation .and. domain%blk(b)%area_variation_active) &
        call Precompute_Area_Geometry(domain%blk(b))
    end do

    ! yn computation gated `model>1`:
    ! the legacy Q2D `if (VISC) Compute_Yn` was defensive (yn computed even in
    ! laminar but unused — downstream gated `if (model==2)` in
    ! Lib_Diffusive/Lib_Source_*Visc). yn is always allocated; Compute_Yn runs
    ! only for `model>1`, saving the Compute_Yn loop in laminar.
    if (model>1) then
      call Compute_Yn ( domain%blk(:), domain%bc, domain%nb, domain%nbound )

      !$omp parallel default(none) shared(domain) private(i, Bm, Im, Jm, Fm, Bs, Is, Js)
      !$omp do schedule(dynamic)
      do i = 1, domain%nbound
        select case (domain%bc(i)%type)
          case(101, 201)
            ! 201 (translational periodic) is a connection seam: its ghost yn is the
            ! donor's wall distance (a scalar, invariant under translation/rotation).
            ! No is_local_block(Bs) guard: Compute_Yn above is global (all blocks, no
            ! guard) and donor nodes are resident through Setup_Metrics, so blks%yn is
            ! valid for a cross-rank donor too. Cf. piano §1 divergenza #21.
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
          ! 201 (translational periodic) is a connection whose donor lives at a
          ! constant Cartesian offset; BC_Connect_Metrics_Periodic rebuilds the ghost
          ! metric from the donor nodes shifted into the local frame (EXACT, vs the
          ! cubic extrapolation a 201 would otherwise get via case default). 101 stays
          ! on BC_Connect_Metrics (MOSE-mirror). Cf. piano §1 divergenza #21.
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
