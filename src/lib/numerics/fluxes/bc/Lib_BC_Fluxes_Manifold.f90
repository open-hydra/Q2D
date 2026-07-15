module Q2D_Lib_BC_Fluxes_Manifold
  ! BC 501, connessione manifold a flusso medio:
  ! integra le proprietà medie (T0, flusso di massa, p_static) di un blocco
  ! "Slave" connesso a una faccia del blocco "Master" attraverso una superficie
  ! potenzialmente non conforme (area_M ≠ area_S).
  !
  ! Riduzione 2D rigorosa (nessuna componente fuori piano w):
  !   * il dispatch sulle facce copre le sole 4 facce 2D (nessuna faccia k).
  !   * nessun loop sul terzo indice (`dim(3)` non esiste).
  !   * int_l è la lunghezza della faccia 2D (analogo dell'area 3D).
  !   * velocità nel solo piano: norm2(meancons_S(nu:nv)).
  use iso_fortran_env,            only: I4 => int32, R8 => real64
  use Q2D_Advanced_Types_m,       only: Q2D_block_type
  use Q2D_Global_m,               only: nsc, nu, nv, np, nt, nprim
  use Q2D_Lib_BC_Fluxes,          only: Face_Index

  implicit none
  private

  public :: BC_Manifold

contains


  subroutine BC_Manifold ( Blkm, Fm, Blks, Fs, T0, g, pstat )
    use FLINT_Lib_Thermodynamic, only: prim2cons, cons2prim, co_rotot_Rtot, f_tabT_expr, cp_tab

    implicit none
    type(Q2D_block_type), intent(in) :: Blkm, Blks
    integer,  intent(in)  :: Fm, Fs
    real(R8), intent(out) :: T0, g, pstat
    ! Variabili locali
    integer  :: im_min, jm_min, im_max, jm_max, dirm
    integer  :: is_min, js_min, is_max, js_max, dirs
    integer  :: im, jm, is, js, Face_i, Face_j
    real(R8) :: int_l, area_M, area_S, Trozza, meanTrozza, meanTrozza_S
    real(R8), dimension(1:np) :: cons, meancons, meancons_S, meanprim_S
    real(R8) :: rho, Rgas, v, T, Tdiff, cp
    integer  :: T_i, Tint(2)

    ! Calcolo della lunghezza della faccia del blocco locale (Master)
    im_min = 1
    jm_min = 1
    im_max = Blkm % dim(1)
    jm_max = Blkm % dim(2)
    select case(Fm)
      case(1)
        im_max = 1
        dirm = 1
      case(2)
        im_min = Blkm % dim(1)
        dirm = 1
      case(3)
        jm_max = 1
        dirm = 2
      case(4)
        jm_min = Blkm % dim(2)
        dirm = 2
    end select

    area_M = 0.0d0

    do jm = jm_min, jm_max; do im = im_min, im_max

      call Face_Index ( Fm, dirm, im, jm, Face_i, Face_j )

      int_l = Blkm % dir(dirm) % f(Face_i,Face_j) % int_l
      area_M = area_M + int_l

    enddo; enddo

    ! Calcolo delle proprietà medie del blocco connesso (Slave)
    is_min = 1
    js_min = 1
    is_max = Blks % dim(1)
    js_max = Blks % dim(2)
    select case(Fs)
      case(1)
        is_max = 1
        dirs = 1
      case(2)
        is_min = Blks % dim(1)
        dirs = 1
      case(3)
        js_max = 1
        dirs = 2
      case(4)
        js_min = Blks % dim(2)
        dirs = 2
    end select

    area_S     = 0.0d0
    meancons   = 0.0d0
    meanTrozza = 0.0d0
    do js = js_min, js_max; do is = is_min, is_max

      call Face_Index ( Fs, dirs, is, js, Face_i, Face_j )

      int_l  = Blks % dir(dirs) % f(Face_i,Face_j) % int_l
      area_S = area_S + int_l

      cons = prim2cons( Blks % P(1:np,is,js) ) * int_l

      meancons = meancons + cons

      ! Stima iniziale di temperatura per cons2prim
      call co_rotot_Rtot ( Blks % P(1,is,js), rho, Rgas )
      Trozza = Blks % P(np,is,js) / ( rho * Rgas ) * int_l
      meanTrozza = meanTrozza + Trozza

    enddo; enddo

    meanTrozza_S = meanTrozza / area_S
    meancons_S   = meancons   / area_S
    meanprim_S   = cons2prim ( meancons_S, meanTrozza_S )

    ! Flusso di massa
    g = norm2(meancons_S(nu:nv)) * area_S / area_M

    ! Temperatura totale
    v = norm2(meanprim_S (nu:nv))
    call co_rotot_Rtot ( meanprim_S(1), rho, Rgas )
    T = meanprim_S(np) / (rho * Rgas)

    T_i     = int ( T )
    Tdiff   = T - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    cp = f_tabT_expr(1, cp_tab, Tint, Tdiff)
    T0 = T + 0.5d0 * v**2 / cp

    ! Pressione statica
    pstat = meanprim_S (np)

  end subroutine BC_Manifold


end module Q2D_Lib_BC_Fluxes_Manifold
