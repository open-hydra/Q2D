module Q2D_Lib_BC_Fluxes_Rotational
  ! Rotational symmetry/periodicity BC.
  !
  !   * BC_Rotational_Symmetry_Eul: inviscid symmetric mirror for the y=0 plane
  !     in 2D axisymmetric: `Flux(nv) = Prim(np) * int_l * Normal(2)`.
  !     Compatible with 2D planar (axis perpendicular to the plane). Currently
  !     without a call-site in the dispatcher: symmetry (BC 300) uses
  !     BC_Symmetry_Eul from Lib_BC_Fluxes_Symmetry.
  !
  !   * BC_Rotational_Periodic_Eul (BC 200): NOT supported. Rotational
  !     periodicity between 2 faces of a 3D sector mesh requires the
  !     delta-theta angle and the swirl velocity component (3 components),
  !     outside the scope of a 2D solver: the routine produces an explicit
  !     `error stop` (no silent fallback). In 2D planar the periodicity
  !     collapses into the translational one, handled by BC 101 (Connection)
  !     with a Cartesian offset.
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nv, np, nprim
  use Q2D_Lib_BC_Fluxes,          only: Compute_Modfm, Face_Index

  implicit none
  private

  public :: BC_Rotational_Symmetry_Eul
  public :: BC_Rotational_Periodic_Eul

contains


  subroutine BC_Rotational_Symmetry_Eul ( Im, Jm, Fm, Blk )
    ! Inviscid mirror: pressure-only contribution on the nv component.
    implicit none
    integer,             intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type),intent(inout) :: Blk
    ! Locals
    integer  :: modfm, modfm1, modfm2, modfm3, dir, Face_i, Face_j
    real(R8) :: Normal(2), int_l
    real(R8) :: Prim(nprim), Flux(nprim)

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Face_i, Face_j )

    Normal = Blk % dir(dir) % f(Face_i,Face_j) % n
    int_l  = Blk % dir(dir) % f(Face_i,Face_j) % int_l
    Prim   = Blk % P(:,Im,Jm)

    Flux = 0d0
    Flux(nv) = Prim(np) * int_l * Normal(2)
    ! Residual update (momentum component nv only)
    Blk % r(nv,Im,Jm) = Blk % r(nv,Im,Jm) + modfm2 * Flux(nv)

  end subroutine BC_Rotational_Symmetry_Eul


  subroutine BC_Rotational_Periodic_Eul ( Im, Jm, Fm, Blk )
    ! Not supported in this 2D solver: sector rotational periodicity
    ! (delta-theta + swirl) requires a 3D mesh and the out-of-plane velocity
    ! component. For 2D planar periodicity use BC 101 (Connection) with a
    ! translational offset.
    implicit none
    integer,             intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type),intent(inout) :: Blk

    write(*,'(A)') '[ERROR] BC_Rotational_Periodic_Eul OUT-OF-SCOPE in 2D planar Q2D.'
    write(*,'(A)') '        Sector-mesh rotational periodicity requires delta-theta angle + swirl component.'
    write(*,'(A)') '        RotatingFrame 2D-axi periodicity is out of scope for this solver.'
    write(*,'(A)') '        Use BC 101 (Connection) with translational offset for 2D planar periodicity.'
    write(*,'(A,3(I0,1X))') '        Caller context (Im, Jm, Fm) = ', Im, Jm, Fm
    if (allocated(Blk % P)) write(*,'(A)') '        Block reference is allocated — runtime trip, not init-time.'
    error stop 'BC_Rotational_Periodic_Eul: out-of-scope (use BC 101 connection)'

  end subroutine BC_Rotational_Periodic_Eul


end module Q2D_Lib_BC_Fluxes_Rotational
