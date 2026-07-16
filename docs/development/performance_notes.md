# Q2D Performance Optimisation Notes

> Target file: `src/lib/numerics/fluxes/Mod_Fluxes.f90` and related
> flux routines.

This page collects performance-tuning notes for the Q2D solver.  The
priorities mirror the MOSE list, adapted for the 2D Q2D loops and the
OMP `do`/`task` runtime dispatcher.

!!! warning "Some notes below were written early and not
    verified against the current source; corrected / flagged where found wrong:"

    - **Priority 6** (AoS→SoA): the example `dl(d,i,j)` was WRONG for Fortran column-major — the
      component index must be LAST (`dl(i,j,d)`); and the layout is byte-shared with MOSE, so SoA is
      a parity-coupled joint change, not Q2D-unilateral. (Both fixed in §Priority 6.)
    - **Priority 3** (precompute): correct as written (`co_rotot_Rtot` per-cell precompute). NOTE the
      thermo property `f_cp` is **already a 1-K lookup table** (`FLINT_Lib_Thermodynamic`), not a
      polynomial — so this note was rewritten to this per-cell framing, NOT a new LUT.
    - **Priorities 1/5** (rlimiter/Riemann static dispatch): the dispatch id is a **runtime** config
      value → a module `integer limiter_id`/`riemann_id` (saved variable), **NOT** `integer,
      parameter` (impossible) and **not** constant-folded; the win is inlining/SIMD.
    - **Priority 4** (fuse loops): the code already runs under ONE `!$omp parallel` region with
      orphaned `!$omp do` — there are no separate `!$omp parallel do` regions to fuse; the only
      lever is `nowait`, unsafe for the `Res`/`beta`-accumulating loops.
    - **Priority 8** (compiler flags): `-march=native` is listed as "currently used" but is NOT in
      `cmake/SetFortranFlags.cmake` (stale); `-flto` is deliberately excluded for MPI +
      mixed-language static-lib fragility, and adding it is a parity-coupled cmake change.

---

## Priority 1 — Replace `rlimiter` procedure pointer with static dispatch

`rlimiter` (in `Q2D_Lib_Limiters`) is a procedure pointer called
inside the `elemental` subroutine `Reconstruction`, which is called
**once per face per RK stage**.  Since the limiter choice does not
change during a run, the pointer prevents:

- **Inlining** — the limiter is a one-liner (e.g., Van Leer: 3
  arithmetic ops), but the indirect call overhead dominates.
- **SIMD vectorisation** — the compiler cannot vectorise through an
  indirect call.

### Suggested fix

Option A — Use a module-level integer `limiter_id` and a `select case`
inside `Reconstruction`:

```fortran
select case (limiter_id)
  case (LIM_MINMOD)  ; slopel = rlimiter_MINMOD(slope1, slope0)
  case (LIM_VANLEER) ; slopel = rlimiter_VANLEER(slope1, slope0)
  case (LIM_LIMO3)   ; slopel = rlimiter_LIMO3(slope1, slope0)   ! LIMO3
  ...
end select
```

Option B — Duplicate `Reconstruction` for each limiter using
preprocessor macros or separate routines, and dispatch at the loop
level in `Fluxes_blk`.

**Expected gain:** 10–20 % on the full solver (limiter is in the
innermost loop).

---

## Priority 2 — Eliminate array temporaries in `Convective_Flux` calls

Each call to `Convective_Flux` typically passes:

| Argument | Example | Issue |
|---|---|---|
| `Prim(:,i-1:i+2,j)` | 4-column slice, dim-2 stride | Non-contiguous → runtime copy |
| `Res(:,i:i+1,j)` | 2-column slice, dim-2 stride | Non-contiguous → runtime copy back |
| `Prim(:,i,j-1:j+2)` | Dir-2: stride across dim-3 | Worse: stride across two dims |

### Suggested fix

Refactor `Convective_Flux` to receive the full `Prim` array plus
index coordinates, or use local stencil variables:

```fortran
real(R8) :: P_stencil(nprim, -1:2), R_face(nprim, 0:1)
! Manually copy the stencil before the call
P_stencil(:,-1) = Prim(:,i-1,j)
P_stencil(:, 0) = Prim(:,i  ,j)
P_stencil(:, 1) = Prim(:,i+1,j)
P_stencil(:, 2) = Prim(:,i+2,j)
```

**Expected gain:** 10–15 %.

---

## Priority 3 — Precompute thermodynamic quantities per cell

`co_rotot_Rtot` (density + gas constant from species mass fractions)
is called redundantly:

- 2× per face in `Convective_Flux` (left + right reconstructed
  states).
- 4× per face in `Diffusive_Flux` (2 in the main body + 2 in
  `Tangential_Gradient`; only one tangential direction in 2D).
- 1× more in `Compute_Diffusive_Flux`.

For a 2D cell with 4 neighbours, that's up to ~12 redundant
evaluations.

### Suggested fix

Add per-cell precomputation before the face loops:

```fortran
real(R8) :: rho_cell(1-gc:n(1)+gc, 1-gc:n(2)+gc)
real(R8) :: T_cell  (1-gc:n(1)+gc, 1-gc:n(2)+gc)
real(R8) :: Rgas_cell(1-gc:n(1)+gc, 1-gc:n(2)+gc)

!$omp do collapse(2)
do j = 1-gc, n(2)+gc
  do i = 1-gc, n(1)+gc
    call co_rotot_Rtot(Prim(1:nsc,i,j), rho_cell(i,j), Rgas_cell(i,j))
    T_cell(i,j) = Prim(np,i,j) / (rho_cell(i,j) * Rgas_cell(i,j))
  end do
end do
```

**Expected gain:** 5–15 %.

---

## Priority 4 — Fuse convective + diffusive loops

Currently, convective and diffusive fluxes for each direction are
computed in separate double loops.  Each loop pass streams through
the entire `Prim` and `Res` arrays, causing redundant cache misses.

### Suggested fix

Merge into a single loop per direction:

```fortran
!$omp do collapse(2)
do j = 1, n(2)
  do i = 1, n(1) - 1
    call Convective_Flux(...)   ! data already in L1/L2 cache
    if (model > 0) call Diffusive_Flux(...)  ! reuse the same data
  end do
end do
```

**Expected gain:** 5–10 % (cache locality).

---

## Priority 5 — Static dispatch for `Riemann`

Same issue as `rlimiter`.  The `Riemann` solver is a procedure
pointer called per face.  While each call is heavier than `rlimiter`
(so the relative overhead is smaller), the indirect call still
blocks inlining and compiler optimisation.

### Suggested fix

Dispatch at the loop level:

```fortran
select case (riemann_id)
  case (RIEMANN_HLLC)
    do j = 1, n(2); do i = 1, n(1)-1
      call riemann_HLLC(...)
    enddo; enddo
  case (RIEMANN_AUSM)
    ...
end select
```

Or use a generic `Riemann` interface with a discriminating integer
argument.

**Expected gain:** 5–10 %.

---

## Priority 6 — AoS → SoA for metric types

The derived types used for mesh metrics store small arrays inside
structs:

```fortran
type :: Q2D_vector_2D_type
  real(R8) :: c(2)         ! dl(i,j)%c(d) — Array of Structures
end type
```

This layout causes stride-2 access when sweeping over cells for a
single component.  The optimal layout for SIMD is Structure of Arrays.

!!! warning "Fortran is column-major — put the component index LAST"

    Fortran stores arrays column-major (first index varies fastest), so
    the SoA component index must be the **trailing** dimension to make a
    fixed-component cell sweep contiguous (stride-1):

    ```fortran
    real(R8) :: dl(1-gc:n1+gc, 1-gc:n2+gc, 2)  ! dl(i,j,d) — d LAST
    ```

    Putting `d` first (`dl(d,i,j)`) reproduces the AoS byte-pattern in
    column-major memory and yields **no** locality benefit — a loop over
    cells with `d` fixed would still stride by 2.

Same applies to `Q2D_tensor_2D_type` and `Q2D_tensor_3D_type` (the
latter is preserved for the Definition A stress storage,
`Stress_t(3,3)`) and `Q2D_f_metrics_type` (`N(2)`, `A`).

**Expected gain:** 5–10 % (better vectorisation).

!!! danger "Parity-coupled — NOT a Q2D-unilateral change"

    The metric derived types (`Q2D_vector_2D_type %c(2)`,
    `Q2D_tensor_2D_type %c(2,2)`, `Q2D_f_metrics_type %N(2)+int_l`,
    `Q2D_d_metrics_type`) are **byte-identical to MOSE** (`MOSE_*_type`)
    modulo the sanctioned 2D↔3D dimensional reduction.
    The AoS layout is a **shared open-hydra architectural choice**, not a
    Q2D inefficiency. An SoA migration would create a Q2D-only divergence
    from byte-shared MOSE → it is **parity-coupled**, applicable only
    as a joint MOSE+Q2D change (upstream back-port), never
    Q2D-unilaterally.

---

## Priority 7 — OMP `do` vs `task` strategy

Q2D ships **two** OpenMP loop strategies in `Mod_Fluxes`, selected
at run time via `[Q2D-Numerics] OMP-parallel = do | task`:

| Strategy | Pragma pattern | When preferred |
|---|---|---|
| `Fluxes_omp_do` (default) | `!$omp parallel + !$omp do collapse(2)` | Diagnostic-friendly |
| `Fluxes_omp_task` | `!$omp single + !$omp taskloop collapse(2) + taskwait` | Imbalanced block decomposition (HPC tuning) |

The `task` variant is a Q2D-only addition.  For balanced block sizes
the `do` variant is usually faster; the `task` variant becomes
beneficial when block sizes vary by more than ≈ 30 % across the
decomposition.

!!! warning "`taskloop` index data-sharing"

    The `taskloop` body declares its loop indices explicitly as
    `private(i, j)`.  Under the OpenMP task data-sharing rules a
    variable that is `private` in the enclosing `parallel` region is
    **firstprivate** (shared by reference) inside a task unless
    re-declared, so omitting the explicit `private(i, j)` lets
    concurrent tasks race on the loop indices.  Keep the explicit
    `private` clause when editing the `task` variant.

Profile on representative cases before changing the default.

---

## Priority 8 — Compiler flags

### Currently used (Release)

`-O3`, `-funroll-loops`, `-ffast-math`, `-finline-functions`,
`-fomit-frame-pointer` (per `cmake/SetFortranFlags.cmake` RELEASE block).
(Note: `-march=native` is NOT actually set — it was listed here in error.)

### Additions to consider

| Flag (gfortran) | Effect |
|---|---|
| `-flto` | Link-time optimisation — enables cross-module inlining. |
| `-fopt-info-vec-missed` | Diagnostic: shows what the compiler fails to vectorise. |
| `-fopenmp-simd` | Enables `!$omp simd` directives for explicit vectorisation hints. |

**Expected gain:** 5–10 % from `-flto` alone (cross-module inlining
of small routines like limiters).

---

## Summary

| # | Action | Effort | Expected gain |
|---|---|---|---|
| 1 | Static dispatch for `rlimiter` | small | 10–20 % |
| 2 | Eliminate `Prim`/`Res` array temporaries | medium | 10–15 % |
| 3 | Precompute `rho`, `T`, `a` per cell | medium | 5–15 % |
| 4 | Fuse convective + diffusive loops | small | 5–10 % |
| 5 | Static dispatch for `Riemann` | small | 5–10 % |
| 6 | AoS → SoA for metric types | large | 5–10 % |
| 7 | OMP `do`/`task` strategy tuning | trivial (run-time flag) | case-dependent |
| 8 | Enable `-flto` | trivial | 5–10 % |
