! WHY: SD_limiter scope. SD_limiter gates ONLY the MUSCL
!   Reconstruction (beta_ = beta if SD_limiter, else 1d0). The Riemann solver
!   always receives raw `beta`, because HLLE++/HLLEMSD/HLLCSD are SD-aware
!   *by design* (cf. Lib_Riemann_HLL header "Modified HLLE with SD/blending");
!   disabling it collapses the solver to plain HLL and removes that behaviour.
module Q2D_Lib_Convective
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Convective_Flux

contains

  subroutine Convective_Flux ( dl, normal, int_l, Prim, Res, beta, SD_limiter )
    use Q2D_Global_m, only: nsc, nu, nv, np, nprim
    use Q2D_Lib_Reconstruction, only: Reconstruction
    use Q2D_Mod_Riemann, only: Riemann
    use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, f_ss, check_gas_state
    implicit none
    real(R8), intent(in), dimension(-1:2) :: dl
    real(R8), intent(in) :: normal(2), int_l
    real(R8), intent(in), dimension(nprim,-1:2) :: Prim
    real(R8), intent(inout), dimension(nprim,0:1) :: Res
    real(R8), intent(in) :: beta
    logical,  intent(in) :: SD_limiter
    real(R8) :: l0, l1, l2, lm, lp
    real(R8), dimension(nprim) :: Prim1, Prim4
    real(R8) :: rotot1, Rtot1, rotot4, Rtot4, a1, a4
    real(R8) :: beta_
    real(R8) :: F_r, F_u, F_v, F_E, su, sp, sm, Flux(nprim)
    integer  :: error1, error4

    l0 = 0.5d0 * ( dl(-1) + dl(0) )
    l1 = 0.5d0 * ( dl(0)  + dl(1) )
    l2 = 0.5d0 * ( dl(1)  + dl(2) )
    lm = 0.5d0 * dl(0)
    lp = 0.5d0 * dl(1)

    if (SD_limiter) then
      beta_ = beta
    else
      beta_ = 1d0
    endif

    call Reconstruction ( Prim(:,-1), Prim(:,0), Prim(:,1), Prim(:,2), &
                          l0, l1, l2, lm, lp, beta_, Prim1, Prim4 )

    call check_gas_state ( Prim1(1:nsc), Prim1(np), error1 )
    call check_gas_state ( Prim4(1:nsc), Prim4(np), error4 )
    if (error1 /= 0 .or. error4 /= 0) then
      write(*,*) '[ERROR] Non-physical state detected at the interface.'
      error stop
    end if

    call co_rotot_Rtot ( Prim1(1:nsc), rotot1, Rtot1 )
    call co_rotot_Rtot ( Prim4(1:nsc), rotot4, Rtot4 )

    a1 = f_ss ( Prim1(1:nsc), Prim1(np), rotot1, Rtot1 )
    a4 = f_ss ( Prim4(1:nsc), Prim4(np), rotot4, Rtot4 )

    call Riemann ( Prim1(1:nsc), Prim1(nu), Prim1(nv), Prim1(np), a1, rotot1, &
                   Prim4(1:nsc), Prim4(nu), Prim4(nv), Prim4(np), a4, rotot4, &
                   beta, normal(1), normal(2), F_r, F_u, F_v, F_E )

    su = sign ( 0.5d0, F_r )
    sp = 0.5d0 + su
    sm = su - 0.5d0

    Flux(1:nsc) = F_r * int_l * ( sp*(Prim1(1:nsc)/rotot1) - sm*(Prim4(1:nsc)/rotot4) )
    Flux(nu)    = F_u * int_l
    Flux(nv)    = F_v * int_l
    Flux(np)    = F_E * int_l

    if (nprim > np) then
      Flux(np+1:nprim) = F_r * int_l * ( sp * Prim1(np+1:nprim)/rotot1 - sm * Prim4(np+1:nprim)/rotot4 )
    end if

    Res(:,0) = Res(:,0) + Flux
    Res(:,1) = Res(:,1) - Flux

  end subroutine Convective_Flux

end module Q2D_Lib_Convective
