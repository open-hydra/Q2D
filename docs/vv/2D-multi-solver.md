# Multi-Solver Coupling BC 103 — Dispatcher Smoke Test

A quiescent 2D box whose faces use **BC 103**, the multi-solver
coupling wall (Q2D legacy code 999).  BC 103 applies the symmetry/wall
condition and then adds an externally supplied coupling flux:

```
BC_Symmetry_Eul + R += ext_flux      (`case (103)` in `src/lib/numerics/fluxes/bc/Mod_BC_Fluxes.f90`)
```

With the default `ext_flux = 0` initialization (BC reader defaults, `src/lib/io/IO_BC.f90`),
BC 103 behaves bit-identically to BC 300 (symmetry/wall).

**Feature exercised:** **BC `103`** multi-solver coupling dispatcher
(Q2D legacy code 999) on all faces, 2D-plane regime, quiescent-IC
equilibrium; MUSCL/vanleer + SLAU.  Dispatcher smoke test only — see
the warning below.

!!! warning "This is a smoke test, not a validation case"

    The ctest `multi-solver-103` verifies **only** that the BC 103
    dispatcher path executes without crashing and preserves the
    quiescent-IC equilibrium with zero external flux.  It does **not**
    verify a coupled two-solver exchange: no partner solver supplies a
    non-zero `ext_flux`, so the coupling data path itself (filling
    `ext_flux` from another code and recovering the correct interface
    flux) is *not* exercised.  What it protects against is a
    dispatcher regression — a broken case branch, index error, or
    uninitialized-flux NaN in the BC 103 code path.

## Problem setup

| Parameter | Value |
|---|---|
| Geometry | single block, 20-cell faces, all faces BC `103` |
| Initial condition | quiescent air, $p_0 = 1.0\times10^{5}$ Pa |
| Iterations | 500 |
| Time scheme | RK3, CFL 0.5 |
| Space reconstruction / limiter | MUSCL / vanleer |
| Riemann solver | SLAU |

## Acceptance gates (`check.py`, ctest `multi-solver-103`)

All hard, NaN-safe (a NaN anywhere is a dispatcher failure):

- the simulation completes 500 iterations with no NaN/Inf anywhere in
  the field;
- min($\rho$), min($p$), min($T$) > 0;
- velocity remains bounded: $\max|u|$, $\max|v|$ < **1 m/s**;
- pressure stays near $p_0$: relative deviation < **1 %**.

## How to run

```bash
cd test/2D/euler/multi-solver-103
./Q2D.sh solve
python3 -B check.py
```
