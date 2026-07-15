!>@brief Driver multigrid 2D:
!>       Setup_Multigrid alloca i blocchi coarse simulation%domain(2..MGL),
!>       chiama Coarse_Grid/Coarse_IOfield e popola la mesh IOfield con un
!>       doppio ciclo (i, j);
!>       Restriction/Prolongation iterano sui blocchi con gating MPI
!>       (is_local_block).
module Q2D_Mod_Multigrid
  use Q2D_Advanced_Types_m, only: Q2D_domain_type, Q2D_simulation_type

  implicit none
  private

  public :: Setup_Multigrid
  public :: Restriction
  public :: Prolongation

contains

  subroutine Setup_Multigrid(simulation)
    use Q2D_Config_Types_m,    only: obj_multigrid, obj_source
    use Q2D_Lib_Multigrid,     only: Check_Multigrid, Coarse_Grid, Coarse_IOfield, &
                                     Coarse_Cross_Section
    use Q2D_Mod_Allocate_Data, only: Allocate_Block
    implicit none
    type(Q2D_simulation_type), intent(inout) :: simulation
    ! Variabili locali
    integer :: b, m, i, j, nb, rap

    call Check_Multigrid(simulation%domain(1))

    nb = simulation%domain(1)%nb

    do m = 2, obj_multigrid%MGL

      allocate(simulation%domain(m)%blk(nb))
      simulation%domain(m)%nb   = nb
      simulation%domain(m)%time = simulation%domain(1)%time

      rap = 2
      do b = 1, nb
        simulation%domain(m)%blk(b)%dim = simulation%domain(m-1)%blk(b)%dim / rap
        call Allocate_Block(simulation%domain(m)%blk(b), simulation%domain(m)%blk(b)%dim)
      end do

      call Coarse_Grid(simulation%domain(m-1), simulation%domain(m))

      ! Restrizione di Cross_Section fine→coarse (sotto-campionamento nodale,
      ! speculare a Coarse_Grid_Blk per node%c). Necessaria per MGL>1: senza,
      ! la Cross_Section dei livelli coarse resterebbe non inizializzata.
      if (obj_source%area_variation) &
        call Coarse_Cross_Section(simulation%domain(m-1), simulation%domain(m))

      call Coarse_IOfield(simulation%IOfield(m-1), simulation%IOfield(m))

      ! Popola la mesh IOfield dai Node%c coarse (ciclo 2D, nessun indice k).
      ! Mesh ORION mesh(coord, i, j, k) rank-4 con placeholder k=0; Q2D_vector_2D_type%c
      ! ha 2 componenti (1:2), la terza coordinata mesh(3,...) resta 0.
      do b = 1, nb
        do j = 0, simulation%domain(m)%blk(b)%dim(2)
        do i = 0, simulation%domain(m)%blk(b)%dim(1)
          simulation%IOfield(m)%block(b)%mesh(1:2, i, j, 0) = simulation%domain(m)%blk(b)%node(i, j)%c(1:2)
          simulation%IOfield(m)%block(b)%mesh(3,   i, j, 0) = 0d0
        end do; end do
      end do

    end do

    ! Legge il numero di iterazioni per ciascun dominio.
    do m = 1, obj_multigrid%MGL
      simulation%domain(m)%itermax = obj_multigrid%iter_threshold(m)
    end do

  end subroutine Setup_Multigrid


  subroutine Restriction(Fine, Coarse)
    use Q2D_Lib_Multigrid, only: fine2coarse_prim
    use Q2D_Mod_MPI,       only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: Fine, Coarse
    ! Variabili locali
    integer :: b

    do b = 1, Coarse%nb ! Ciclo sui blocchi
      if (.not. is_local_block(b)) cycle
      ! Restrizione pesata sulla superficie di cella 'sur' (in 2D svolge il ruolo
      ! del volume).
      call fine2coarse_prim(Fine%blk(b)%P,   Coarse%blk(b)%P,   &
                            Fine%blk(b)%sur, Coarse%blk(b)%sur, &
                            Fine%blk(b)%dim, Coarse%blk(b)%dim)
    end do

  end subroutine Restriction


  subroutine Prolongation(Fine, Coarse)
    use Q2D_Lib_Multigrid, only: coarse2fine_prim
    use Q2D_Mod_MPI,       only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: Fine, Coarse
    ! Variabili locali
    integer :: b

    do b = 1, Coarse%nb ! Ciclo sui blocchi
      if (.not. is_local_block(b)) cycle
      call coarse2fine_prim(Fine%blk(b)%P, Coarse%blk(b)%P, &
                            Fine%blk(b)%dim, Coarse%blk(b)%dim)
    end do

  end subroutine Prolongation

end module Q2D_Mod_Multigrid
