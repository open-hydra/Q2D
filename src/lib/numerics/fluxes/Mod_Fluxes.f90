module Q2D_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Assign_Fluxes

  procedure(Fluxes_if), pointer, public :: Fluxes

  abstract interface
    subroutine Fluxes_if(domain)
      use Q2D_Advanced_Types_m, only: Q2D_domain_type
      implicit none
      type(Q2D_domain_type), intent(inout) :: domain
    end subroutine Fluxes_if
  end interface

contains

  ! Q2D-extra: dispatcher for the OMP strategy
  ! ('do' = classic !$omp do, 'task' = !$omp taskloop+taskgroup intra-block).
  ! Intel ifx workaround (explicit taskwait) preserved in Fluxes_blk_omp_task.
  subroutine Assign_Fluxes(Fluxes_word)
    implicit none
    character(len=*), intent(inout) :: Fluxes_word

    nullify(Fluxes)

    select case (Fluxes_word)
    case ('do')
      Fluxes => Fluxes_omp_do
      Fluxes_word = 'OMP parallelism, omp do'
    case ('task')
      Fluxes => Fluxes_omp_task
      Fluxes_word = 'OMP parallelism, omp task'
    case default
      Fluxes => Fluxes_omp_do
      Fluxes_word = 'OMP parallelism, omp do'
    end select

  end subroutine Assign_Fluxes


  subroutine Fluxes_omp_do ( domain )
    ! Reads obj_rans/obj_soot/obj_space_scheme from singletons, declares local
    ! Sc/Sct/Prt/soot_enabled/SD_limiter, propagates them down to block-level.
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_shock_detector, obj_space_scheme, obj_rans, obj_soot
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer  :: b
    logical  :: soot_enabled, SD_limiter
    real(R8) :: Sc, Sct, Prt

    soot_enabled = obj_soot%enabled
    Sc  = obj_rans%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    SD_limiter = obj_space_scheme%SD

    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks
      if (obj_shock_detector%SD) then
        call Fluxes_blk_omp_do ( domain % blk(b) % P,                  &
                                 domain % blk(b) % r,                  &
                                 domain % blk(b) % dir,                &
                                 domain % blk(b) % dl,                 &
                                 domain % blk(b) % yn,                 &
                                 domain % blk(b) % m,                  &
                                 domain % blk(b) % node,               &
                                 domain % blk(b) % dim,                &
                                 domain % blk(b) % area_geom_term_x,   &
                                 domain % blk(b) % area_geom_term_y,   &
                                 Sc, Sct, Prt, soot_enabled,           &
                                 SD_limiter,                           &
                                 beta = domain % blk(b) % beta         )
      else
        call Fluxes_blk_omp_do ( domain % blk(b) % P,                  &
                                 domain % blk(b) % r,                  &
                                 domain % blk(b) % dir,                &
                                 domain % blk(b) % dl,                 &
                                 domain % blk(b) % yn,                 &
                                 domain % blk(b) % m,                  &
                                 domain % blk(b) % node,               &
                                 domain % blk(b) % dim,                &
                                 domain % blk(b) % area_geom_term_x,   &
                                 domain % blk(b) % area_geom_term_y,   &
                                 Sc, Sct, Prt, soot_enabled,           &
                                 SD_limiter                            )
      end if
    end do

  end subroutine Fluxes_omp_do


  subroutine Fluxes_blk_omp_do ( Prim, Res, Dir, dl, yn, M, node, n, geom_x, geom_y, &
                                 Sc, Sct, Prt, soot_enabled, SD_limiter, beta )
    ! Signature: Sc/Sct/Prt/soot_enabled/SD_limiter passed as parameters.
    ! geom_x/geom_y + optional beta are Q2D-extra (V_dot_G + SD).
    use Q2D_Global_m,            only: nprim, np, gc, model
    use Q2D_Base_Types_m,        only: Q2D_d_metrics_type, Q2D_vector_2D_type, Q2D_tensor_2D_type
    use Q2D_Lib_Convective,      only: Convective_Flux
    use Q2D_Lib_Diffusive,       only: Diffusive_Flux
    use Q2D_Lib_Shock_Detector,  only: SD_Tramel
    implicit none
    integer, dimension(2), intent(in) :: n
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)  :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(out) :: Res
    type(Q2D_d_metrics_type), dimension(2), intent(in) :: Dir
    type(Q2D_vector_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: dl
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: yn
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: M
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)), intent(in) :: node
    real(R8), dimension(1:n(1), 1:n(2)), intent(in) :: geom_x, geom_y
    real(R8), intent(in) :: Sc, Sct, Prt
    logical,  intent(in) :: soot_enabled, SD_limiter
    real(R8), dimension(1:n(1), 1:n(2)), intent(inout), optional :: beta
    integer :: i, j

    !$omp do collapse(2)
    do j = 1, n(2)
    do i = 1, n(1)
      Res(:,i,j) = 0d0
    end do; end do

    ! Shock detector (2D with a 3x3 stencil on
    ! Prim(np, i-1:i+1, j-1:j+1)). beta(i,j) is active only if
    ! obj_shock_detector%SD (dispatcher gate Fluxes_omp_do); passed as a scalar
    ! to Convective_Flux for Roe-Rusanov blending.
    if (present(beta)) then
      !$omp do collapse(2)
      do j = 1, n(2)
      do i = 1, n(1)
        beta(i,j) = SD_Tramel( Prim(np, i-1:i+1, j-1:j+1) )
      end do; end do
    end if

    !$omp do
    do j = 1, n(2)
    do i = 1, n(1) - 1

      if (present(beta)) then
        call Convective_Flux ( dl(i-1:i+2,j) % c(1), &
                               Dir(1) % f(i,j) % N,  &
                               Dir(1) % f(i,j) % int_l, &
                               Prim(:,i-1:i+2,j),    &
                               Res (:,i:i+1,j),      &
                               beta=beta(i,j), SD_limiter=SD_limiter )
      else
        call Convective_Flux ( dl(i-1:i+2,j) % c(1), &
                               Dir(1) % f(i,j) % N,  &
                               Dir(1) % f(i,j) % int_l, &
                               Prim(:,i-1:i+2,j),    &
                               Res (:,i:i+1,j),      &
                               beta=1d0, SD_limiter=SD_limiter )
      end if

      if (model>0) then
        call Diffusive_Flux ( Dir(1) % f(i,j) % n,    &
                              Dir(1) % f(i,j) % int_l,&
                              yn(i  ,j),              &
                              yn(i+1,j),              &
                              Prim(:,i  ,j),          &
                              Prim(:,i+1,j),          &
                              Prim(:,i  ,j-1),        &
                              Prim(:,i  ,j+1),        &
                              Prim(:,i+1,j-1),        &
                              Prim(:,i+1,j+1),        &
                              M(i  ,j) % c,           &
                              M(i+1,j) % c,           &
                              Res (:,i  ,j),          &
                              Res (:,i+1,j),          &
                              node (i,j  ) % c(2),    &
                              node (i,j-1) % c(2),    &
                              1, 2,                   &
                              Sc, Sct, Prt, soot_enabled, &
                              Gx1=geom_x(i  ,j),      &
                              Gy1=geom_y(i  ,j),      &
                              Gx2=geom_x(i+1,j),      &
                              Gy2=geom_y(i+1,j)       )
      end if

    end do; end do

    !$omp do
    do i = 1, n(1)
    do j = 1, n(2) - 1

      if (present(beta)) then
        call Convective_Flux ( dl(i,j-1:j+2) % c(2), &
                               Dir(2) % f(i,j) % n,  &
                               Dir(2) % f(i,j) % int_l, &
                               Prim(:,i,j-1:j+2),    &
                               Res (:,i,j:j+1),      &
                               beta=beta(i,j), SD_limiter=SD_limiter )
      else
        call Convective_Flux ( dl(i,j-1:j+2) % c(2), &
                               Dir(2) % f(i,j) % n,  &
                               Dir(2) % f(i,j) % int_l, &
                               Prim(:,i,j-1:j+2),    &
                               Res (:,i,j:j+1),      &
                               beta=1d0, SD_limiter=SD_limiter )
      end if

      if (model>0) then
        call Diffusive_Flux ( Dir(2) % f(i,j) % n,    &
                              Dir(2) % f(i,j) % int_l,&
                              yn(i,j  ),              &
                              yn(i,j+1),              &
                              Prim(:,i,j  ),          &
                              Prim(:,i,j+1),          &
                              Prim(:,i-1,j  ),        &
                              Prim(:,i+1,j  ),        &
                              Prim(:,i-1,j+1),        &
                              Prim(:,i+1,j+1),        &
                              M(i,j  ) % c,           &
                              M(i,j+1) % c,           &
                              Res (:,i,j  ),          &
                              Res (:,i,j+1),          &
                              node (i  ,j) % c(2),    &
                              node (i-1,j) % c(2),    &
                              2, 1,                   &
                              Sc, Sct, Prt, soot_enabled, &
                              Gx1=geom_x(i,j  ),      &
                              Gy1=geom_y(i,j  ),      &
                              Gx2=geom_x(i,j+1),      &
                              Gy2=geom_y(i,j+1)       )
      end if

    end do; end do

  end subroutine Fluxes_blk_omp_do


  subroutine Fluxes_omp_task ( domain )
    ! Analogous to Fluxes_omp_do: reads obj_rans/obj_soot/obj_space_scheme,
    ! propagates Sc/Sct/Prt/soot_enabled/SD_limiter to block-level. Q2D-extra HPC
    ! variant (OMP-task; the MOSE peer has no such variant — Q2D-only).
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_shock_detector, obj_space_scheme, obj_rans, obj_soot
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    integer  :: b
    logical  :: soot_enabled, SD_limiter
    real(R8) :: Sc, Sct, Prt

    soot_enabled = obj_soot%enabled
    Sc  = obj_rans%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    SD_limiter = obj_space_scheme%SD

    !$omp single
    !$omp taskgroup
    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle  ! MPI: skip non-local blocks
      if (obj_shock_detector%SD) then
        !$omp task firstprivate(b) default(shared)
        call Fluxes_blk_omp_task ( domain % blk(b) % P,                  &
                                   domain % blk(b) % r,                  &
                                   domain % blk(b) % dir,                &
                                   domain % blk(b) % dl,                 &
                                   domain % blk(b) % yn,                 &
                                   domain % blk(b) % m,                  &
                                   domain % blk(b) % node,               &
                                   domain % blk(b) % dim,                &
                                   domain % blk(b) % area_geom_term_x,   &
                                   domain % blk(b) % area_geom_term_y,   &
                                   Sc, Sct, Prt, soot_enabled,           &
                                   SD_limiter,                           &
                                   beta = domain % blk(b) % beta         )
        !$omp end task
      else
        !$omp task firstprivate(b) default(shared)
        call Fluxes_blk_omp_task ( domain % blk(b) % P,                  &
                                   domain % blk(b) % r,                  &
                                   domain % blk(b) % dir,                &
                                   domain % blk(b) % dl,                 &
                                   domain % blk(b) % yn,                 &
                                   domain % blk(b) % m,                  &
                                   domain % blk(b) % node,               &
                                   domain % blk(b) % dim,                &
                                   domain % blk(b) % area_geom_term_x,   &
                                   domain % blk(b) % area_geom_term_y,   &
                                   Sc, Sct, Prt, soot_enabled,           &
                                   SD_limiter                            )
        !$omp end task
      end if
    end do
    !$omp end taskgroup
    !$omp end single

  end subroutine Fluxes_omp_task


  subroutine Fluxes_blk_omp_task ( Prim, Res, Dir, dl, yn, M, node, n, geom_x, geom_y, &
                                   Sc, Sct, Prt, soot_enabled, SD_limiter, beta )
    ! Signature: Sc/Sct/Prt/soot_enabled/SD_limiter passed as parameters
    ! (analogous to Fluxes_blk_omp_do). HPC OMP-task variant Q2D-extra.
    use Q2D_Global_m,            only: nprim, np, gc, model
    use Q2D_Base_Types_m,        only: Q2D_d_metrics_type, Q2D_vector_2D_type, Q2D_tensor_2D_type
    use Q2D_Lib_Convective,      only: Convective_Flux
    use Q2D_Lib_Diffusive,       only: Diffusive_Flux
    use Q2D_Lib_Shock_Detector,  only: SD_Tramel
    implicit none
    integer, dimension(2), intent(in) :: n
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in)  :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(out) :: Res
    type(Q2D_d_metrics_type), dimension(2), intent(in) :: Dir
    type(Q2D_vector_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: dl
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: yn
    type(Q2D_tensor_2D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc), intent(in) :: M
    type(Q2D_vector_2D_type), dimension(0:n(1), 0:n(2)), intent(in) :: node
    real(R8), dimension(1:n(1), 1:n(2)), intent(in) :: geom_x, geom_y
    real(R8), intent(in) :: Sc, Sct, Prt
    logical,  intent(in) :: soot_enabled, SD_limiter
    real(R8), dimension(1:n(1), 1:n(2)), intent(inout), optional :: beta
    integer :: i, j

    ! NOTE: taskwait directives below are theoretically redundant (OMP 5.0
    ! Sec. 2.10.2: taskloop without nogroup includes an implicit taskgroup).
    ! However, Intel Fortran's runtime does not correctly synchronize the
    ! implicit taskgroup of taskloop when nested inside an explicit task
    ! region. Without these taskwait, SIGSEGV occurs. Tested with ifx/ifort.
    ! Keep these until the compiler bug is fixed.
    !
    ! NOTE: explicit `private(i,j)` clauses added per
    ! OMP 5.0 §2.19.1 rule 2 ("A loop iteration variable for a sequential
    ! loop in a parallel or task generating construct is private in the
    ! innermost such construct that encloses the loop"). Empirical: gfortran
    ! 7.5.0 `!$omp taskloop default(shared)` does NOT auto-privatize inner
    ! sequential loop iteration variables (test_taskloop_private.f90: 7750
    ! race errors / 50000 elements without private(i); 0 errors with).
    ! `!$omp do` auto-privatizes correctly. ifx/ifort honor rule 2 for both.
    ! Explicit private(i,j) makes the spec-prescribed semantics enforceable
    ! across compilers.
    !$omp taskgroup
    !$omp taskloop collapse(2) default(shared) private(i,j)
    do j = 1, n(2)
    do i = 1, n(1)
      Res(:,i,j) = 0d0
    end do; end do
    !$omp end taskloop
    !$omp taskwait

    ! Shock detector (task strategy; mirror Fluxes_blk_omp_do).
    if (present(beta)) then
      !$omp taskloop collapse(2) default(shared) private(i,j)
      do j = 1, n(2)
      do i = 1, n(1)
        beta(i,j) = SD_Tramel( Prim(np, i-1:i+1, j-1:j+1) )
      end do; end do
      !$omp end taskloop
      !$omp taskwait
    end if

    !$omp taskloop default(shared) private(i)
    do j = 1, n(2)
      do i = 1, n(1) - 1
        if (present(beta)) then
          call Convective_Flux ( dl(i-1:i+2,j) % c(1), &
                                 Dir(1) % f(i,j) % n,  &
                                 Dir(1) % f(i,j) % int_l, &
                                 Prim(:,i-1:i+2,j),    &
                                 Res (:,i:i+1,j),      &
                                 beta=beta(i,j), SD_limiter=SD_limiter )
        else
          call Convective_Flux ( dl(i-1:i+2,j) % c(1), &
                                 Dir(1) % f(i,j) % n,  &
                                 Dir(1) % f(i,j) % int_l, &
                                 Prim(:,i-1:i+2,j),    &
                                 Res (:,i:i+1,j),      &
                                 beta=1d0, SD_limiter=SD_limiter )
        end if

        if (model>0) then
          call Diffusive_Flux ( Dir(1) % f(i,j) % n,    &
                                Dir(1) % f(i,j) % int_l,&
                                yn(i  ,j),              &
                                yn(i+1,j),              &
                                Prim(:,i  ,j),          &
                                Prim(:,i+1,j),          &
                                Prim(:,i  ,j-1),        &
                                Prim(:,i  ,j+1),        &
                                Prim(:,i+1,j-1),        &
                                Prim(:,i+1,j+1),        &
                                M(i  ,j) % c,           &
                                M(i+1,j) % c,           &
                                Res (:,i  ,j),          &
                                Res (:,i+1,j),          &
                                node (i,j  ) % c(2),    &
                                node (i,j-1) % c(2),    &
                                1, 2,                   &
                                Sc, Sct, Prt, soot_enabled, &
                                Gx1=geom_x(i  ,j),      &
                                Gy1=geom_y(i  ,j),      &
                                Gx2=geom_x(i+1,j),      &
                                Gy2=geom_y(i+1,j)       )
        end if
      end do
    end do
    !$omp end taskloop
    !$omp taskwait

    !$omp taskloop default(shared) private(j)
    do i = 1, n(1)
      do j = 1, n(2) - 1
        if (present(beta)) then
          call Convective_Flux ( dl(i,j-1:j+2) % c(2), &
                                 Dir(2) % f(i,j) % n,  &
                                 Dir(2) % f(i,j) % int_l, &
                                 Prim(:,i,j-1:j+2),    &
                                 Res (:,i,j:j+1),      &
                                 beta=beta(i,j), SD_limiter=SD_limiter )
        else
          call Convective_Flux ( dl(i,j-1:j+2) % c(2), &
                                 Dir(2) % f(i,j) % n,  &
                                 Dir(2) % f(i,j) % int_l, &
                                 Prim(:,i,j-1:j+2),    &
                                 Res (:,i,j:j+1),      &
                                 beta=1d0, SD_limiter=SD_limiter )
        end if

        if (model>0) then
          call Diffusive_Flux ( Dir(2) % f(i,j) % n,    &
                                Dir(2) % f(i,j) % int_l,&
                                yn(i,j  ),              &
                                yn(i,j+1),              &
                                Prim(:,i,j  ),          &
                                Prim(:,i,j+1),          &
                                Prim(:,i-1,j  ),        &
                                Prim(:,i+1,j  ),        &
                                Prim(:,i-1,j+1),        &
                                Prim(:,i+1,j+1),        &
                                M(i,j  ) % c,           &
                                M(i,j+1) % c,           &
                                Res (:,i,j  ),          &
                                Res (:,i,j+1),          &
                                node (i  ,j) % c(2),    &
                                node (i-1,j) % c(2),    &
                                2, 1,                   &
                                Sc, Sct, Prt, soot_enabled, &
                                Gx1=geom_x(i,j  ),      &
                                Gy1=geom_y(i,j  ),      &
                                Gx2=geom_x(i,j+1),      &
                                Gy2=geom_y(i,j+1)       )
        end if
      end do
    end do
    !$omp end taskloop
    !$omp taskwait
    !$omp end taskgroup

  end subroutine Fluxes_blk_omp_task

end module Q2D_Mod_Fluxes
