module Q2D_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_BC
  public :: Print_BC_Summary

  ! Counters per BC type — populated by Read_BCfile, consumed by Print_BC_Summary
  integer :: nconnect, nwall, ngsi, nSRM, nio, nsym, nper, nchimera, ncoupled, next, nmanifold, nstate
  logical :: has_tdep_bc

contains

  subroutine Setup_BC ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    use Q2D_Config_Types_m,   only: obj_multigrid, obj_io_bc
    use Q2D_Global_m,         only: model
    implicit none
    type(Q2D_domain_type), intent(inout) :: domain(:)
    ! Local variables
    integer :: m, error

    do m = 1, obj_multigrid%MGL

      !! Phase 1: allocation and validation
      call Allocate_BC ( domain(m) )

      ! 2D solver: 4 faces per block. The flags are allocated as (nb*6, 4):
      ! the factor 6 in the first dimension is kept for layout symmetry with
      ! the 3D six-face case, but only 4 faces are actually used.
      if (allocated(obj_io_bc%viscous_flag))  deallocate(obj_io_bc%viscous_flag)
      allocate ( obj_io_bc%viscous_flag( 6 * domain(m) % nb, 4 ) )
      obj_io_bc%viscous_flag = .false.

      if (allocated(obj_io_bc%coupling_flag)) deallocate(obj_io_bc%coupling_flag)
      allocate ( obj_io_bc%coupling_flag( 6 * domain(m) % nb, 4 ) )
      obj_io_bc%coupling_flag = .false.

      error = Check_BC ( domain(m) % nbound, m )
      if (error /= 0) cycle

      !! Phase 2: reading
      call Read_BCfile ( domain(m) % bc, domain(m) % n_bf, m )

      if ((any(obj_io_bc%viscous_flag)) .and. (model==0)) then
        obj_io_bc%error_message = '[ERROR] Viscous wall BCs assigned with Eulerian model'
      end if

    end do

  end subroutine Setup_BC


  subroutine Allocate_BC ( domain )
    use Q2D_Advanced_Types_m, only: Q2D_domain_type
    implicit none
    class(Q2D_domain_type) :: domain
    ! Local variables
    integer :: b, n, ni, nj

    ! 4 2D faces per block — n = sum_b (2*nj + 2*ni).
    n = 0
    do b = 1, domain % nb
      ni = domain % blk(b) % dim(1)
      nj = domain % blk(b) % dim(2)
      n = n + 2*nj + 2*ni
    end do

    domain % nbound = n
    allocate ( domain % bc ( domain % nbound ) )
    allocate ( domain % n_bf ( domain % nb, 4 ) )

  end subroutine Allocate_BC


  function Check_BC ( n, level ) result(ios)
    use Q2D_Config_Types_m, only: obj_io_bc
    use Q2D_Global_m,       only: Q2D_phase_prefix
    use IR_Precision,       only: str
    implicit none
    integer, intent(in) :: n, level
    integer             :: ios
    integer             :: di(5), ti, unitfile, n_proof
    integer             :: c, ci, cii

    ios = 0

    ! Open the file
    if (level == 1) then
      open(newunit=unitfile, file='INPUT/'//trim(Q2D_phase_prefix)//'bc.txt', &
           status='old', iostat=ios, action='read')
    else
      open(newunit=unitfile, file='INPUT/'//trim(Q2D_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt', &
           status='old', iostat=ios, action='read')
    end if
    if (ios /= 0) then
      obj_io_bc%error_message = '[ERROR] Boundary condition file not found for grid '//trim(str(.true.,level))
      return
    end if

    ! Consistency check of the BC file.
    ! The header (ATLAS format, meshType==-2) is ALWAYS 5 integers
    ! (b, i, j, f, type); no k field, no autodetection of alternative
    ! formats. bc%type is read DIRECTLY: the codes are exactly those
    ! emitted by ATLAS, with no translation.
    ! Initialization `n_proof = -1` with no exit guard: the EOF iteration
    ! performs an empty increment that compensates the initial -1, so the
    ! final value of `n_proof` equals the number of headers read successfully.
    ios = 0; n_proof = -1
    do while (ios == 0)
      read(unitfile, *, iostat=ios) di(1), di(2), di(3), di(4), ti
      ! BC type 410 (mapped-state) is not supported in Q2D: it is rejected
      ! already at validation with an explicit message, before Read_BCfile
      ! terminates with error stop.
      if (ti == 410) then
        obj_io_bc%error_message = '[ERROR] BC type 410 (mapped-state) is not supported. ' // &
                                  'Use BC 401-407 (inflow) or BC 400 (extrapolation) instead.'
        ! The nonzero return value makes Setup_BC skip Read_BCfile
        ! (cycle): the rejection surfaces through Mod_Diagnostic::Check_BC as
        ! a clean '[ERROR] ...' diagnostic, instead of the raw error stop in
        ! Read_BCfile (kept as a defensive safeguard).
        ios = -1
        close(unitfile)
        return
      end if
      ! Dispatch to skip the correct number of additional lines.
      ! IDs with ONE property line:
      !   101/103 = connection/coupled | 201 = periodic | 301-309 = wall |
      !   401-407,420 = inflow/outflow | 501 = manifold | 502 = srm
      ! IDs WITHOUT a property line:
      !   200 = axisymmetry | 300 = symmetry | 400 = extrapolation
      ! IDs with VARIABLE-length property lines:
      !   102 = chimera
      select case (ti)
      case (101, 103, 201, 301:309, 401:407, 420, 501:502)
        read(unitfile, *, iostat=ios)
      case (102) ! chimera: variable-length property block
        read(unitfile, *, iostat=ios) ci, cii
        do c = 1, ci+cii
          read(unitfile, *, iostat=ios)
        end do
      case (200, 300, 400)
        ! no property line
      end select
      n_proof = n_proof + 1
    end do

    if (n_proof /= n) then
      obj_io_bc%error_message = '[ERROR] Boundary conditions number ('//str(.true.,n_proof) &
                              //') is different than the one of the initial conditions ('//str(.true.,n)//')'
      close(unitfile)
      return
    end if

    ! Validation passed: reset ios (nonzero because of the EOF) to signal success.
    ios = 0
    close(unitfile)

  end function Check_BC


  subroutine Read_BCfile ( bc, n_bf, level )
    use Q2D_Advanced_Types_m, only: Q2D_bc_type
    use Q2D_Config_Types_m,   only: obj_io, obj_io_bc
    use Q2D_Global_m,         only: nsc, nrans, nprim, Q2D_phase_prefix
    use IR_Precision,         only: str
    implicit none
    type(Q2D_bc_type), dimension(:), intent(inout) :: bc
    integer,                          intent(in)    :: level
    integer, dimension(1:,1:),        intent(inout) :: n_bf
    ! Local variables
    integer :: cc, i, s
    integer :: unitfile, ios
    character(len=256) :: p0file
    character(len=32)  :: alpha_tok, beta_tok

    ! Open the file
    if (level == 1) then
      open(newunit=unitfile, file='INPUT/'//trim(Q2D_phase_prefix)//'bc.txt', &
           status='old', iostat=ios, action='read')
    else
      open(newunit=unitfile, file='INPUT/'//trim(Q2D_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt', &
           status='old', iostat=ios, action='read')
    end if
    if (ios /= 0) then
      obj_io_bc%error_message = '[ERROR] Boundary condition file not found for grid '//trim(str(.true.,level))
      return
    end if

    ! Counters (reset only at grid level 1)
    if (level == 1) then
      nconnect = 0; nwall = 0; ngsi = 0; nSRM = 0
      nio = 0; nsym = 0; nper = 0; nchimera = 0
      ncoupled = 0; nmanifold = 0; next = 0; nstate = 0
      has_tdep_bc = .false.
    end if

    ! Counter of the number of cells per face in each block
    n_bf = 0

    ! Read the file
    do i = 1, size(bc)

      ! ── First line: block, ij, face, ATLAS BC ID ─────────────────────────────
      ! Header (ATLAS format, meshType==-2) of 5 integers (b, i, j, f, type).
      ! No k field: ks=1 is imposed. The codes are exactly those emitted by
      ! ATLAS — bc%type is read DIRECTLY.
      read(unitfile, *, iostat=ios) bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%f, bc(i)%type
      bc(i)%ks = 1
      if (ios /= 0) write(*,'(A,I0)') '  Error in BC file at record ', i

      ! Update n_bf
      n_bf( bc(i)%b, bc(i)%f ) = n_bf( bc(i)%b, bc(i)%f ) + 1

      ! ── Second line: property data, depending on the BC type ─────────────────
      ! 2D form: the connection/coupled property line is 8 integers
      ! (bs,is,js,fs + d11,d12,d21,d22); no ks. The 2nd dimension of Pg is 4
      ! (4 faces in 2D).
      select case ( bc(i)%type )

        ! ─────────────────────────────────────────────────────────────────────
        ! Connection and periodic BCs
        case (101, 201)
          if (level == 1) nconnect = nconnect + 1
          read(unitfile, *, iostat=ios) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          allocate( bc(i)%Pg(nprim, 4) )

        ! ─────────────────────────────────────────────────────────────────────
        ! Chimera overset BC
        case (102)
          if (level == 1) nchimera = nchimera + 1
          read(unitfile, *, iostat=ios) (bc(i)%ni(cc), cc=1,2)
          allocate( bc(i)%donorID(1:sum(bc(i)%ni), 1:4) )
          allocate( bc(i)%volume_fraction(1:sum(bc(i)%ni)) )
          do s = 1, bc(i)%ni(1)
            read(unitfile, *, iostat=ios) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
          end do
          do s = bc(i)%ni(1)+1, bc(i)%ni(1)+bc(i)%ni(2)
            read(unitfile, *, iostat=ios) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
          end do
          allocate( bc(i)%Pg(nprim, 4) )

        ! ─────────────────────────────────────────────────────────────────────
        ! Axisymmetry
        case (200)
          if (level == 1) nper = nper + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Eulerian symmetry
        case (300)
          if (level == 1) nsym = nsym + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Extrapolation
        case (400)
          if (level == 1) next = next + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, prescribed heat flux
        ! Second line: q, roughness_ks, emissivity_eps
        case (301)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b, bc(i)%f ) = .true.
          read(unitfile, *, iostat=ios) bc(i)%qw, bc(i)%k_rough, bc(i)%eps_wall

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, prescribed temperature
        ! Second line: T, roughness_ks, emissivity_eps
        case (302)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b, bc(i)%f ) = .true.
          read(unitfile, *, iostat=ios) bc(i)%Tw, bc(i)%k_rough, bc(i)%eps_wall

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, temperature + radiative flux (melting)
        ! Second line: T, qrad, roughness_ks
        case (303)
          if (level == 1) nwall = nwall + 1
          if (level == 1) ngsi  = ngsi  + 1
          obj_io_bc%viscous_flag( bc(i)%b, bc(i)%f ) = .true.
          read(unitfile, *, iostat=ios) bc(i)%Tw, bc(i)%qrad, bc(i)%k_rough

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, radiative flux (pyrolysis/ablation)
        ! Second line: qrad, roughness_ks
        case (304)
          if (level == 1) nwall = nwall + 1
          if (level == 1) ngsi  = ngsi  + 1
          obj_io_bc%viscous_flag( bc(i)%b, bc(i)%f ) = .true.
          read(unitfile, *, iostat=ios) bc(i)%qrad, bc(i)%k_rough

        ! ─────────────────────────────────────────────────────────────────────
        ! Inflow, stagnation conditions (T0, p0)
        ! Second line: T0, p0, alpha, beta, rel_fac, massf(1:nsc), turb(1:nrans)
        case (401)
          if (level == 1) nio = nio + 1
          allocate( bc(i)%ci(1 : nsc+nrans) )
          read(unitfile, *, iostat=ios) &
            bc(i)%T0, bc(i)%p0, alpha_tok, beta_tok, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Inflow, time-varying p0
        ! Second line: T0, p0_time_file, alpha, beta, rel_fac, massf, turb
        case (402)
          if (level == 1) nio = nio + 1
          has_tdep_bc = .true.
          allocate( bc(i)%ci(1 : nsc+nrans) )
          read(unitfile, *, iostat=ios) &
            bc(i)%T0, p0file, alpha_tok, beta_tok, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          call bc(i)%p0time%initialize(file=p0file, bar=.true.)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Inflow, mass flux g + T0
        ! Second line: T0, g, alpha, beta, rel_fac, massf, turb
        case (403)
          if (level == 1) nio = nio + 1
          allocate( bc(i)%ci(1 : nsc+nrans) )
          read(unitfile, *, iostat=ios) &
            bc(i)%T0, bc(i)%mdot, alpha_tok, beta_tok, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Inflow, mass flux g + T (static)
        ! Second line: T, g, alpha, beta, rel_fac, massf, turb
        case (404)
          if (level == 1) nio = nio + 1
          allocate( bc(i)%ci(1 : nsc+nrans) )
          read(unitfile, *, iostat=ios) &
            bc(i)%T0, bc(i)%mdot, alpha_tok, beta_tok, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Supersonic inflow, M + T (static) + p (static)
        ! Second line: mach, T, p, alpha, beta, rel_fac, massf, turb
        case (405)
          if (level == 1) nio = nio + 1
          allocate( bc(i)%ci(1 : nsc+nrans) )
          read(unitfile, *, iostat=ios) &
            bc(i)%mach, bc(i)%T0, bc(i)%pamb, alpha_tok, beta_tok, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Outflow, prescribed back pressure (extrapolation if zero)
        ! Second line: p_back, rel_fac
        case (406)
          if (level == 1) nio = nio + 1
          read(unitfile, *, iostat=ios) bc(i)%pamb, bc(i)%rel_fac

        ! ─────────────────────────────────────────────────────────────────────
        ! Ambient, T0 + p0 + back pressure
        ! Second line: T0, p0, p_back, alpha, beta, rel_fac, massf, turb
        case (407)
          if (level == 1) nio = nio + 1
          allocate( bc(i)%ci(1 : nsc+nrans) )
          read(unitfile, *, iostat=ios) &
            bc(i)%T0, bc(i)%p0, bc(i)%pamb, alpha_tok, beta_tok, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Mapped state (BC 410) — not supported in Q2D: rejected with an
        ! explicit error stop (no silent fallback, no partial implementation).
        ! The mirrored, earlier rejection is in Check_BC.
        case (410)
          obj_io_bc%error_message = '[ERROR] BC type 410 (mapped-state) is not supported. '           // &
                                    'Use BC 401-407 (inflow) or BC 400 (extrapolation) instead.'
          error stop '[FATAL] BC type 410 is not supported (see error_message)'

        ! ─────────────────────────────────────────────────────────────────────
        ! Choked nozzle (Borda flow encoded with BC_Mach==-20 in the mach
        ! field). pamb carries the Borda area ratio A1/A3, consumed as A1_A3
        ! by BC_Nozzle_Borda_Flow.
        ! Second line: mach, T0, p0, alpha, beta, pamb, rel_fac, massf, turb
        case (420)
          if (level == 1) nio = nio + 1
          allocate( bc(i)%ci(1 : nsc+nrans) )
          read(unitfile, *, iostat=ios) &
            bc(i)%mach, bc(i)%T0, bc(i)%p0, alpha_tok, beta_tok, bc(i)%pamb, bc(i)%rel_fac, &
            (bc(i)%ci(s), s = nrans+1, nrans+nsc), (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Mean-flow manifold connection
        case (501)
          if (level == 1) nmanifold = nmanifold + 1
          read(unitfile, *, iostat=ios) bc(i)%bs, bc(i)%fs

        ! ─────────────────────────────────────────────────────────────────────
        ! SRM grain combustion (BC 502).
        ! Reads Taf, aCoeff, n, pRef, rhoGrain, SF_geo, massf(1:nsc);
        ! haf forced to 0 — not present in the ATLAS output format.
        case (502)
          if (level == 1) nSRM = nSRM + 1
          allocate( bc(i)%ci(1 : nsc) )
          bc(i)%haf = 0.0_R8   ! haf absent in ATLAS output; set to zero
          read(unitfile, *, iostat=ios) bc(i)%Taf, bc(i)%aCoeff, bc(i)%n, bc(i)%pRef, &
            bc(i)%rhoGrain, bc(i)%SF_geo, (bc(i)%ci(s), s = 1, nsc)

        ! ─────────────────────────────────────────────────────────────────────
        ! Multi-solver coupled wall
        case (103)
          if (level == 1) ncoupled = ncoupled + 1
          obj_io_bc%coupling_flag( bc(i)%b, bc(i)%f ) = .true.
          read(unitfile, *, iostat=ios) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          allocate( bc(i)%ext_flux(nprim) )
          bc(i)%ext_flux = 0.0_R8
          allocate( bc(i)%Pg(1, 4) )

        ! Unhandled codes (e.g. 0 = zero-gradient extrapolation along the
        ! collapsed direction, used for the trivial j direction in the 1D
        ! shock-tube tests, and the range 305:309, which has no body) fall
        ! through with no effect — no property line, no allocation. The
        ! absence of a `case default` is intentional.

      end select

    end do

    if (nwall > 0) obj_io%write_wall = .true.

    close(unitfile)

  end subroutine Read_BCfile


  subroutine Print_BC_Summary ()
    implicit none

    write(*,*)
    write(*,'(A)') ' Boundary conditions'
    if (nconnect  > 0) write(*,'(A,T35,I0)') '   Connection',         nconnect
    if (nwall     > 0) write(*,'(A,T35,I0)') '   Viscous wall',       nwall
    if (ngsi      > 0) write(*,'(A,T35,I0)') '   GSI (inc. wall)',    ngsi
    if (nio       > 0) write(*,'(A,T35,I0)') '   Inflow/outflow',     nio
    if (nsym      > 0) write(*,'(A,T35,I0)') '   Symmetry',           nsym
    if (nper      > 0) write(*,'(A,T35,I0)') '   Periodicity',        nper
    if (next      > 0) write(*,'(A,T35,I0)') '   Extrapolation',      next
    if (nchimera  > 0) write(*,'(A,T35,I0)') '   Chimera',            nchimera
    if (ncoupled  > 0) write(*,'(A,T35,I0)') '   Coupled wall',       ncoupled
    if (nmanifold > 0) write(*,'(A,T35,I0)') '   Manifold',           nmanifold
    if (nSRM      > 0) write(*,'(A,T35,I0)') '   SRM grain',          nSRM
    if (nstate    > 0) write(*,'(A,T35,I0)') '   Q2D mapped (410)',   nstate
    if (has_tdep_bc)   write(*,'(A)')        '   Time-dependent BC detected'

  end subroutine Print_BC_Summary


  ! ---------------------------------------------------------------------------
  ! parse_dir_tok: converts a direction token (read as character) into a real.
  !   If the token contains 'normal', returns huge(R8) (flag for a direction
  !   normal to the face). Otherwise it is interpreted as a real number.
  ! ATLAS writes the literal token 'normal' for alpha/beta in types
  ! 401/402/403/404/405/407/420 when the inflow direction is normal to the
  ! face.
  ! ---------------------------------------------------------------------------
  pure function parse_dir_tok(tok) result(val)
    implicit none
    character(len=*), intent(in) :: tok
    real(R8)                     :: val
    integer                      :: ios_loc
    character(len=len(tok))      :: tok_

    tok_ = adjustl(tok)
    if (index(trim(tok_), 'normal') > 0) then
      val = huge(1.0_R8)
    else
      read(tok_, *, iostat=ios_loc) val
      if (ios_loc /= 0) val = 0.0_R8
    end if
  end function parse_dir_tok

end module Q2D_IO_BC
