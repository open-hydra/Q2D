# SRM Rocket Motor — Vieille Burn-Rate BC

A 2D axisymmetric solid rocket motor whose closed-end grain face
injects combustion products through the **BC 502 SRM boundary
condition** with the Vieille (Saint-Robert) burn-rate law
$r_b = a\,(P/P_\text{ref})^n$.  The converged chamber pressure is
validated against the analytic Vieille steady-state balance with a
choked throat, and the products against their total temperature.  Q2D
applies the grain as MOSE does (MOSE: `AblatingWall`, `BCInput`); the
ctest `SRM-rocket-motor-inflow` runs the same chamber with the grain in
mode 2 (an inflow), and `SRM-grain-state` checks the state of the gas a
grain lets in over one step, in both modes, on a low face (face 1) and on
a high face (face 4, along the contour of the nozzle), and the pressure on
the grain face, that of the cell, with the gas moving away from it.

**Reference:** analytic Vieille (Saint-Robert) steady-state
chamber-pressure balance — the full derivation and every input
constant are committed in the header of
`reference/P_eq_vieille.dat`.

**Feature exercised:** **BC `502`** SRM grain-face injection with the
Vieille burn-rate law, the products at the total temperature $T_{af}$; **2D axisymmetric** geometric source
(`[Q2D-Sources] active = axisymmetric`); MUSCL/vanleer + SLAU;
chemistry-labelled ctest (`2D;chemistry;validation;fast`).

## Problem setup

The geometry reuses the chang-nozzle 2D-axi mesh (190 × 40 cells,
$R_\text{chamber} = 0.0635$ m, $R_\text{throat} = 0.0203$ m,
contraction ratio $K_n = 9.75$).  BC 502 parameters (from
`INPUT/bc.txt`):

| Parameter | Value |
|---|---|
| Mass flux of the propellant gas $\dot m_\text{ref}$ (`mdot_ref`) | 50 kg/(m²·s) at $P_\text{ref}$ |
| Burn-rate exponent $n$ | 0.3 |
| Reference pressure $P_\text{ref}$ | $1.0\times10^{6}$ Pa |
| Flame temperature $T_{af}$, the total temperature of the products | 3000 K |
| Mode | 1, a wall that blows the products in (2 in `SRM-rocket-motor-inflow`) |
| Gas | air-like, $R = 287.05$ J/kg/K, $\gamma = 1.4$ |

`mdot_ref` is the mass flux of the propellant gas at $P_\text{ref}$, the
linear burn rate times the grain density (and its scale factor), as the
boundary-condition builder of ATLAS writes it.  The contraction ratio is
small for a real SRM, so the mass flux is boosted (50 kg/m²/s vs ~9 for
real AP/HTPB at $P_\text{ref} = 1$ MPa) to bring the choked-throat
Vieille equilibrium within reach of a ctest budget (8000 iterations at
CFL 0.5, the steady state of 40000 at CFL 0.1).

**Analytic reference.** The Vieille steady-state balance with a
choked throat,

$$a\,(P/P_\text{ref})^n\,A_\text{chamber} = \frac{P\,A_\text{throat}}{c^*}
  \;\Rightarrow\;
  P_\text{eq} = \left[\frac{a\,c^*\,(A_\text{chamber}/A_\text{throat})}{P_\text{ref}^{\,n}}\right]^{1/(1-n)},$$

evaluates to $P_\text{eq} = 5.532\times10^{5}$ Pa
($c^* = 1355.25$ m/s), committed in `reference/P_eq_vieille.dat`
together with the full input list.

## Numerical setup

| Parameter | Value |
|---|---|
| Regime | 2D axisymmetric (`[Q2D-Sources] active = axisymmetric`) |
| Time scheme | RK3, local time stepping (8000 iterations) |
| CFL | 0.5 (rise threshold 100) |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | SLAU |

## Results

<figure>
  {% include "vv/images/SRM-pressure.svg" %}
</figure>

## Acceptance gates (`check.py`, ctest `SRM-rocket-motor`)

All hard:

1. **Vieille chamber pressure**: the chamber-averaged pressure
   (cells with $x < x_\text{throat} - 0.005$ m) vs the analytic
   $P_\text{eq}$: $|P_\text{chamber} - P_\text{eq}|/P_\text{eq}$ <
   **0.15**.  The 15 % tolerance (per the script's own rationale) is
   the chang-nozzle-pattern ±5 % avg / ±15 % max acceptance relaxed
   for (a) the spatial-vs-lumped discretization residual of comparing
   a 2D field against a 0D balance: the mean over the chamber cells
   takes the lower pressure of the converging part — measured
   $\varepsilon_\text{rel}$ = **3.9e-2** at the steady state (gfortran
   RELEASE build).
2. **Convergence sentinel**: final density residual < **0.5**
   (NaN-safe): measured 3.7e-5 at the steady state; the threshold
   catches divergence / NaN blow-up.
3. **NaN-safe positivity**: $\rho$, $p$, $T$, $g$, $R$ > 0; $u$, $v$
   finite.
4. **Total temperature of the products**: the grain lets the gas in
   with the total enthalpy $h(T_{af})$ of the tables and the steady
   Euler flow keeps it, so that $T + |\mathbf{u}|^2/(2c_p) = T_{af}$ in
   every cell: largest deviation < **5e-3** of $T_{af}$ (measured
   1.8e-3).  The former BC 502 injected the products with $h = 0$, the
   gas next to the grain at 15 K: the chamber showed a quasi-steady
   pressure at 5000 iterations only because the cold, dense gas barely
   moved.

The same gates hold for the grain in mode 2 (`SRM-rocket-motor-inflow`,
`run-inflow.sh`).

## How to run

```bash
cd test/2D/euler/SRM-rocket-motor
./Q2D.sh solve
python3 -B check.py
```
