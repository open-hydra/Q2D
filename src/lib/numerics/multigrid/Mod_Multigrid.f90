!>@brief Multigrid driver — port of MOSE Mod_Multigrid with structural 3D->2D
!>       reduction:
!>       Setup_Multigrid allocates simulation%domain(2..MGL) coarse blocks +
!>       calls Coarse_Grid/Coarse_IOfield + populates the IOfield mesh with a
!>       double (i, j) loop only (k dropped);
!>       Restriction/Prolongation iterate over the blocks with MPI gating
!>       (is_local_block).
!>       The 3D branch is dropped (Dim(3)=Max(1,Dim(3)) in MOSE Mod_Multigrid
!>       is stripped).
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
    ! Local
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

      ! R1 restriction of Cross_Section fine→coarse (nodal sub-sampling, mirror
      ! Coarse_Grid_Blk for node%c). Fixes Cross_Section corruption for MGL>1.
      if (obj_source%area_variation) &
        call Coarse_Cross_Section(simulation%domain(m-1), simulation%domain(m))

      call Coarse_IOfield(simulation%IOfield(m-1), simulation%IOfield(m))

      ! Populate IOfield mesh from coarse Node%c (2D loop, no k).
      ! ORION mesh(coord, i, j, k) rank-4 with k=0 placeholder; Q2D_vector_2D_type%c
      ! has 2 components (1:2), the 3rd coordinate mesh(3,...) stays 0.
      do b = 1, nb
        do j = 0, simulation%domain(m)%blk(b)%dim(2)
        do i = 0, simulation%domain(m)%blk(b)%dim(1)
          simulation%IOfield(m)%block(b)%mesh(1:2, i, j, 0) = simulation%domain(m)%blk(b)%node(i, j)%c(1:2)
          simulation%IOfield(m)%block(b)%mesh(3,   i, j, 0) = 0d0
        end do; end do
      end do

    end do

    ! Read number of iterations for each domain.
    do m = 1, obj_multigrid%MGL
      simulation%domain(m)%itermax = obj_multigrid%iter_threshold(m)
    end do

  end subroutine Setup_Multigrid


  subroutine Restriction(Fine, Coarse)
    use Q2D_Lib_Multigrid, only: fine2coarse_prim
    use Q2D_Mod_MPI,       only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b

    do b = 1, Coarse%nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle
      ! 'sur' Q2D-form replaces the MOSE 3D 'vol'.
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
    ! Local
    integer :: b

    do b = 1, Coarse%nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle
      call coarse2fine_prim(Fine%blk(b)%P, Coarse%blk(b)%P, &
                            Fine%blk(b)%dim, Coarse%blk(b)%dim)
    end do

  end subroutine Prolongation

end module Q2D_Mod_Multigrid
