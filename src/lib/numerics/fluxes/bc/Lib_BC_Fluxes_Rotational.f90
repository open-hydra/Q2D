module Q2D_Lib_BC_Fluxes_Rotational
  ! BC di simmetria/periodicità rotazionale.
  !
  !   * BC_Rotational_Symmetry_Eul: mirror simmetrico inviscido per il piano
  !     y=0 in 2D assialsimmetrico: `Flux(nv) = Prim(np) * int_l * Normal(2)`.
  !     Compatibile con il 2D piano (asse perpendicolare al piano). Attualmente
  !     senza call-site nel dispatcher: la simmetria (BC 300) usa
  !     BC_Symmetry_Eul di Lib_BC_Fluxes_Symmetry.
  !
  !   * BC_Rotational_Periodic_Eul (BC 200): NON supportata. La periodicità
  !     rotazionale tra 2 facce di una mesh 3D a settore richiede l'angolo
  !     delta-theta e la componente di swirl della velocità (3 componenti),
  !     fuori dallo scope di un solutore 2D: la routine produce un `error stop`
  !     esplicito (nessun fallback silenzioso). In 2D piano la periodicità
  !     collassa in quella traslazionale, gestita dalla BC 101 (Connection)
  !     con offset cartesiano.
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
    ! Mirror inviscido: contributo di sola pressione sulla componente nv.
    implicit none
    integer,             intent(in)    :: Im, Jm, Fm
    type(Q2D_block_type),intent(inout) :: Blk
    ! Locali
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
    ! Aggiornamento del residuo (sola componente di quantità di moto nv)
    Blk % r(nv,Im,Jm) = Blk % r(nv,Im,Jm) + modfm2 * Flux(nv)

  end subroutine BC_Rotational_Symmetry_Eul


  subroutine BC_Rotational_Periodic_Eul ( Im, Jm, Fm, Blk )
    ! Non supportata in questo solutore 2D: la periodicità rotazionale a
    ! settore (delta-theta + swirl) richiede una mesh 3D e la componente di
    ! velocità fuori piano. Per la periodicità 2D piana usare la BC 101
    ! (Connection) con offset traslazionale.
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
