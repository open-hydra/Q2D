# SRM Rocket Motor — Vieille Burn-Rate BC

A 2D axisymmetric solid rocket motor whose closed-end grain face
injects combustion products through the **BC 502 SRM boundary
condition** with the Vieille (Saint-Robert) burn-rate law
$r_b = a\,(P/P_\text{ref})^n$.  The converged chamber pressure is
validated against the analytic Vieille steady-state balance with a
choked throat.  MOSE has `Lib_BC_Fluxes_SRM` but no
SRM test-case mirror.

**Reference:** analytic Vieille (Saint-Robert) steady-state
chamber-pressure balance — the full derivation and every input
constant are committed in the header of
`reference/P_eq_vieille.dat`.

**Feature exercised:** **BC `502`** SRM grain-face injection with the
Vieille burn-rate law; **2D axisymmetric** geometric source
(`[Q2D-Sources] active = axisymmetric`); MUSCL/vanleer + SLAU;
chemistry-labelled ctest (`2D;chemistry;validation;fast`).

## Problem setup

The geometry reuses the chang-nozzle 2D-axi mesh (190 × 40 cells,
$R_\text{chamber} = 0.0635$ m, $R_\text{throat} = 0.0203$ m,
contraction ratio $K_n = 9.75$).  BC 502 parameters (from
`INPUT/bc.txt`):

| Parameter | Value |
|---|---|
| Flame temperature $T_{af}$ | 3000 K |
| Burn-rate coefficient $a$ | 50 kg/(m²·s) at $P_\text{ref}$ |
| Burn-rate exponent $n$ | 0.3 |
| Reference pressure $P_\text{ref}$ | $1.0\times10^{6}$ Pa |
| Grain density $\rho_G$ | 1800 kg/m³ |
| Gas | air-like, $R = 287.05$ J/kg/K, $\gamma = 1.4$ |

The contraction ratio is small for a real SRM, so the burn-rate
coefficient is boosted (50 kg/m²/s vs ~9 for real AP/HTPB at
$P_\text{ref} = 1$ MPa) to bring the choked-throat Vieille
equilibrium within reach of a ctest budget (5000 iterations).

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
| Time scheme | RK3, local time stepping (5000 iterations) |
| CFL | 0.1 (rise threshold 100) |
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
   a 2D field against a 0D balance, and (b) the residual
   pressurization transient: at the 5000-iteration ctest budget the
   chamber is still slowly filling toward the asymptote — measured
   $\varepsilon_\text{rel}$ = **9.4e-2** (2026-07-15, gfortran
   TESTING build).
2. **Convergence sentinel**: final density residual < **0.5**
   (NaN-safe).  The run is quasi-steady, not deep-converged: the
   measured final value is 2.6e-1 (slow monotonic chamber
   pressurization); the threshold is ~2× the measured final and
   catches divergence / NaN blow-up.
3. **NaN-safe positivity**: $\rho$, $p$, $T$, $g$, $R$ > 0; $u$, $v$
   finite.

## How to run

```bash
cd test/2D/euler/SRM-rocket-motor
./Q2D.sh solve
python3 -B check.py
```
