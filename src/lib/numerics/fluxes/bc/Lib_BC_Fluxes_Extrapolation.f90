module Q2D_Lib_BC_Fluxes_Extrapolation
  ! Zero-gradient extrapolation with
  ! the flux built from the boundary state (no Riemann, no reconstruction).
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim, model
  use Q2D_Lib_BC_Fluxes,          only: Face_Index
  use FLINT_Lib_Thermodynamic,    only: f_Rtot, H0

  implicit none
  private

  public :: BC_Extrapolation

contains


  subroutine BC_Extrapolation ( Im, Jm, Fm, Blk )
    implicit none
    integer, intent(in) :: Im, Jm, Fm
    type(Q2D_block_type), intent(inout) :: Blk
    ! Local
    integer :: Dir, Face_i, Face_j, FluxSign
    real(R8) :: Normal(2), int_l, Prim(nprim), Un, Rgas, Fmass, Flux(nprim)


    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    Normal = Blk % dir(Dir) % f(Face_i,Face_j) % n
    int_l  = Blk % dir(Dir) % f(Face_i,Face_j) % int_l

    Prim = Blk % P(:,Im,Jm)
    Un   = dot_product ( Prim(nu:nv), Normal )

    Rgas        = f_Rtot ( Prim(1:nsc) )
    Flux(1:nsc) = Prim(1:nsc) * Un * int_l
    Fmass       = sum ( Flux(1:nsc) )
    Flux(nu:nv) = Fmass * Prim(nu:nv) + Prim(np) * int_l * Normal
    Flux(np)    = Fmass * H0 ( Prim(np), Prim(1:nsc), norm2 ( Prim(nu:nv) ) )
    ! Gate `nprim>np` extrapolates ALL post-pressure variables (turbulence + soot
    ! + passive). The Q2D legacy gate `model==2 / nt:nprim` lost soot/passive in
    ! non-RANS mode (Q2D legacy oversight, not intentional preservation).
    if (nprim > np) Flux(np+1:nprim) = Prim(np+1:nprim) * Un * int_l

    FluxSign = 1 - 2 * mod(Fm,2)
    Blk % r(:,Im,Jm) = Blk % r(:,Im,Jm) + FluxSign * Flux

  end subroutine BC_Extrapolation


end module Q2D_Lib_BC_Fluxes_Extrapolation
