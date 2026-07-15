module Q2D_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_dt, Set_Global_dt

contains

  subroutine Compute_dt ( domain, cfl, vnn, rampa_iter )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Global_m,         only: nsc, nu, nv, np, nt, nprim
    use Q2D_Mod_MPI,          only: is_local_block, mpi_allreduce_min_r8
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: cfl, vnn
    integer,  intent(in) :: rampa_iter
    ! Local variables
    integer  :: i, j, b
    real(R8) :: dtcell, dtglobal, dtglobal_mpi

    dtglobal = domain % dtglobal

    do b = 1, domain % nb ! loop over blocks
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(2) private ( dtcell ), reduction ( min : dtglobal )
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)

        ! Compute the local cell dt according to the CFL and VNN numbers
        call compute ( rhoi = domain%blk(b)%P(1:nsc,i,j), &
                       vel  = domain%blk(b)%P(nu:nv,i,j), &
                       p    = domain%blk(b)%P(np,i,j), &
                       rans_= domain%blk(b)%P(nt:nprim,i,j), &
                       met  = domain%blk(b)%M(i,j)%c, &
                       dl   = domain%blk(b)%dl(i,j)%c, &
                       dtmin = dtcell, &
                       cfl  = cfl, &
                       vnn  = vnn )

        ! Apply the CFL ramp if requested
        if ( domain%iter < rampa_iter ) dtcell = dtcell * domain%iter / rampa_iter

        ! Update the local cell dt and the global minimum
        domain % blk(b) % dtlocal(i,j) = dtcell
        dtglobal = min ( dtcell, dtglobal )

      enddo; enddo
      !$omp end parallel
    enddo ! end of the loop over blocks

    ! MPI: global minimum over all ranks
    call mpi_allreduce_min_r8(dtglobal, dtglobal_mpi)
    domain % dtglobal = dtglobal_mpi

    contains

      subroutine compute ( rhoi, vel, p, rans_, met, dl, dtmin, cfl, vnn )
        use Q2D_Global_m,            only: nsc, ndir, model
        use Q2D_Lib_RANS,            only: Eddy_Viscosity
        use FLINT_Lib_Thermodynamic, only: co_rotot_Rtot, f_ss, f_laminarViscosity
        implicit none
        real(R8), intent(in)  :: rhoi(nsc), vel(2), p, rans_(:)
        real(R8), intent(in)  :: met(2,2), dl(2), cfl, vnn
        real(R8), intent(out) :: dtmin
        ! Local variables
        integer  :: d
        real(R8) :: rho, Rgas, Sound, dt, versor(2), lambda, mie, mil, mi
        real(R8) :: vel_, dummy(2,2)=0d0

        call co_rotot_Rtot ( rhoi, rho, Rgas )
        Sound = f_ss ( rhoi, p, rho, Rgas )

        dtmin = 1d8
        do d = 1, ndir
          ! CFL condition along direction d
          versor = met(d,:) / norm2 ( met(d,:) )
          vel_ = abs ( dot_product ( vel, versor ) )

          lambda  = vel_ + Sound
          dt = dl(d) / lambda * cfl
          dtmin = min ( dt, dtmin )
        enddo

        if ( model==0 ) return

        mie = 0.d0
        mil = f_laminarViscosity ( rhoi, p, rho, Rgas )
        if (model==2) then
          ! Approximation: vel_gradient=0 and small walldist.
          ! vr=0 zeroes the strain also in the axisymmetric case → dummy.
          call Eddy_Viscosity ( mut=mie, rans_variables=rans_, mul=mil, rho=rho, &
                                vel_gradient=dummy, vr=0.d0, r_ave=1.d0, walldist=1d-6 )
        endif
        mi = mie + mil

        do d = 1, ndir
          ! VNN condition along direction d
          dt = ( rho * dl(d)**2 * vnn ) / mi
          dtmin = min ( dt, dtmin )
        enddo

      end subroutine compute

  end subroutine Compute_dt


  subroutine Set_Global_dt ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Mod_MPI,          only: is_local_block
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain
    ! Local variables
    integer :: b, i, j

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(2)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        domain % blk(b) % dtlocal(i,j) = domain % dtglobal
      enddo; enddo
      !$omp end parallel
    enddo

  end subroutine Set_Global_dt

end module Q2D_Mod_dt
