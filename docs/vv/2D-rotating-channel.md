# Rotating Channel — Coriolis Geostrophic Balance

A straight 2D plane channel with uniform axial throughflow in a frame
rotating about the out-of-plane axis.  The case verifies the Coriolis
source-term coupling of [the rotating-frame theory
page](../theory/rotating-frame.md) in isolation: it is the MOSE 3D
`coriolis-channel` mirror reduced to the 2D plane
(axis $= [0,0,1]$, strict No-$w$).

**Feature exercised:** rigid-body **rotating frame**
(`rotational-frame = rigid-body`, `[Q2D-Rotating-Frame]`
$\omega = 10$ rad/s, axis $[0,0,1]$) — Coriolis source term
(`Lib_RotatingFrame.f90`) in the 2D-plane regime; **BC `201`**
periodic block connection (`Lib_BC_Fluxes_Connection`) in $x$ +
BC `300` walls in $y$; MUSCL/vanleer + SLAU.

## Problem setup

**Geometry and grid.** Channel $x \in [-L_x/2, L_x/2]$,
$y \in [-L_y/2, L_y/2]$ with $L_x = 0.4$ m, $L_y = 0.2$ m, centered
at the origin.

**Rotating frame** (from `input.ini` `[Q2D-Rotating-Frame]`):
$\omega = 10$ rad/s, axis $[0, 0, 1]$, origin $(0,0,0)$.

**Initial/inflow state.** Uniform axial flow $u = U_\infty = 100$ m/s
(Mach 0.288) in quiescent air at $T_0 = 300$ K,
$p_0 = 1.0\times10^{5}$ Pa ($\rho_0 = 1.16144$ kg/m³, $R = 287$
J/kg/K, $\gamma = 1.4$).

**Boundary conditions.** Periodic in $x$ (BC `201`, the
block-connection path `Lib_BC_Fluxes_Connection`) and inviscid walls
in $y$ (BC `300`) — the throughflow recirculates indefinitely while
the walls support the transverse pressure gradient.

**Expected behaviour.** With $\boldsymbol{\Omega} = \omega\hat{z}$
and $\boldsymbol{u} = U_\infty\hat{x}$ the Coriolis acceleration is

$$\mathrm{Cor}_x = +2\omega v = 0, \qquad
  \mathrm{Cor}_y = -2\omega u = -2\omega U_\infty,$$

so the steady isothermal equilibrium is the geostrophic balance

$$\frac{\partial p}{\partial y} = -2\,\omega\,\rho\,U_\infty,
  \qquad p(y) = p_0\, e^{-2\omega U_\infty y /(R T_0)} .$$

The analytic slope for this configuration is
$\mathrm{d}p/\mathrm{d}y = -2\omega\rho_0 U_\infty = -2322.88$ Pa/m
(committed in `reference/dpdy_analytical.dat`).  The centrifugal
contribution is negligible by construction: Cen/Cor = 1.118 % at the
domain corner, i.e. a Coriolis-dominant regime.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3, local time stepping (2000 iterations) |
| CFL | 0.1 (rise threshold 100) |
| IRS β | 0.0 |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | SLAU |

## Results

<figure>
  {% include "vv/images/RC-dpdy.svg" %}
</figure>

## Acceptance gates (`check.py`, ctest `rotating-channel`)

All hard, NaN-safe:

- sanity: min($\rho$), min($p$), min($T$) > 0;
- **pressure-gradient balance**: the fitted CFD slope
  $\mathrm{d}p/\mathrm{d}y$ (linear fit of the $x$-averaged pressure
  vs $y$) vs the analytic $-2\omega\rho_0 U_\infty$:
  $\varepsilon_\text{rel}$ < **5e-2**;
- **drift gates**: the IC ($u = U_\infty$, $v = 0$) is an exact fixed
  point of the equilibrium, so the solved field must not drift from
  it:
    - $\max|u - U_\infty|$ < **1.5e-2 m/s** (measured 6.1e-3 on the
      committed 2000-iteration run),
    - $\max|v|$ < **8.0e-2 m/s** (measured 3.7e-2),

  both set at ~2× the measured discretization-level settling of the
  committed run.

## How to run

```bash
cd test/2D/euler/rotating-channel
./Q2D.sh solve
python3 -B check.py
```
