# Noh Cylindrical Implosion — Axisymmetric Source Verification

The cylindrical Noh (1987) implosion in the 2D axisymmetric $(x, r)$
regime: a uniform cold gas moves radially toward the symmetry axis at
constant speed and an infinite-strength shock reflects from the axis.
This is the only axisymmetric case in the suite with an **analytic
exact solution** — it verifies the axisymmetric geometric source terms
of [the axisymmetric-derivation theory
page](../theory/axisymmetric-derivation.md) in isolation (the other
axi cases validate against digitized published data or committed
reference fields).

**Feature exercised:** **axisymmetric geometric source terms**
(`[Q2D-Sources] active = axisymmetric`,
`Lib_Source_Axisymmetry.f90` — every Euler $1/r$ term, including the
$r \to 0$ axis treatment via $1/\max(r, \varepsilon)$) with an
analytic-exact reference; BC `300` on the axis row (the chang-nozzle
convention) and as $x$-symmetry; BC `400` outer extrapolation;
first-order + HLLE near vacuum ($p_0 = 10^{-9}$).

## Problem setup

**Geometry and grid.** Thin axisymmetric strip: $x \in [0, 0.01]$ m
(4 cells, $\Delta x = \Delta r$), $r \in [0, 2.5]$ m (1000 cells,
$\Delta r = 2.5\times10^{-3}$).  The exact solution is $x$-invariant,
so the strip doubles as a negative control for spurious $x$-fluxes.
Everything is regenerated deterministically by `gen_case.py`
(`INPUT/ic.tec`, `INPUT/bc.txt`, `MESH/mesh.tec`).

**Initial state** (mirror of the 1D [Noh87](1D.md#noh-implosion-problem)
choices; same thermally-perfect air table, constant $C_p$ →
$\gamma = 1.4$): $\rho_0 = 1$ kg/m³, $u_r = -U_0 = -1$ m/s everywhere,
$u_x = 0$, $p_0 = 10^{-9}$ Pa.  The 1D case uses $p_0 = 10^{-15}$,
exact there because the planar pre-shock state is uniform; in
cylindrical geometry the pre-shock compression is non-uniform and the
truncation/roundoff noise in $p$ sits at the $10^{-14}$ scale, so the
floor is raised to $10^{-9}$ (still dynamically negligible:
$a_0 \approx 4\times10^{-5} \ll U_0$, $p_0/p_\text{post} \sim 10^{-10}$).

**Boundary conditions.** Axis row $r = 0$: BC `300` (symmetry — how
chang-nozzle marks its axis; BC `200` is the rotational-periodic code,
not an axis BC).  $x$-faces: BC `300` (mirror ghosts are exact for an
$x$-invariant solution).  Outer $r = 2.5$: BC `400` extrapolation —
the outer flow is a pressureless supersonic *inflow*, so the
zero-gradient ghost mis-states the incoming mass flux of the
converging profile at an $O(1)$ rate; no time-dependent
supersonic-inflow BC exists in the dispatcher, so the boundary error
is *outrun* instead: it advects inward at $|u_r| + a \approx U_0$,
reaching only $r \approx 1.5$ by $t_\text{final}$, and all gates are
restricted to the validated region $r \le 1$ (measured density error
$1.8\times10^{-4}$ at $r = 1.4$ vs $1.8\times10^{-2}$ at $r = 2.0$).

**Reference:** W. F. Noh, _J. Comput. Phys._ **72** (1987) 78–120.

**Exact solution** (Noh 1987, $m = 1$;
`test/common/pycheck/q2d_exact1d.py::noh_exact`, cross-checked against
the classic $\gamma = 5/3$ values $D = 1/3$, $\rho_2 = 16$,
$p_2 = 16/3$).  The reflected shock moves outward at

$$D = \frac{\gamma - 1}{2}\, U_0 = 0.2 \text{ m/s}, \qquad
  r_s(t) = D\,t ,$$

with the uniform stagnation state behind it and the geometric
pre-shock compression ahead of it:

$$r < r_s:\quad
  \rho = \rho_0 \left(\frac{\gamma+1}{\gamma-1}\right)^{2} = 36, \quad
  u_r = 0, \quad
  p = \frac{\gamma+1}{2}\,\rho_0 U_0^2\,\frac{\gamma+1}{\gamma-1} = 7.2 ,$$

$$r > r_s:\quad
  \rho = \rho_0 \left(1 + \frac{U_0\, t}{r}\right), \quad
  u_r = -U_0, \quad
  p \simeq p_0 .$$

The run is time-accurate to $t_\text{final} = 1$ s
(`time-threshold = 1.0`), so $r_s = 0.2$ m; the check evaluates the
exact solution at the actual output `SOLUTIONTIME` (< 0.1 %
overshoot).

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3, time-accurate to $t = 1.0$ s |
| CFL | 0.50 |
| Space reconstruction | first-order |
| Riemann solver | HLLE |
| Sources | `active = axisymmetric` |

(Scheme settings mirror the 1D Noh87 case.)

## Results

<figure>
  {% include "vv/images/NohCyl.svg" %}
</figure>

The density dip toward the axis is the classic Noh *wall-heating*
artifact (Noh 1987, §5), common to shock-capturing schemes; the gates
below account for it.

## Acceptance gates (`check.py`, ctest `noh-cylindrical`)

All hard, NaN-safe, evaluated on the $x$-averaged radial profile for
$r \le 1$ (calibrated on the 1000-cell run, 2026-07-15):

- sanity: min($\rho$), min($p$), min($T$) > 0; $u$, $v$ finite;
- **whole-region rel L1($\rho$)** vs exact: < **8e-2** (measured
  4.77e-2, dominated by first-order shock smearing + the axis
  wall-heating dip; ~1.7× headroom);
- **post-shock plateau density**, mean over $r \in (0.4, 0.9)\,r_s$
  (window excludes the wall-heating dip and the smeared shock), vs
  the analytic $\rho_0((\gamma+1)/(\gamma-1))^2$:
  $\varepsilon_\text{rel}$ < **5e-2** (measured 2.53e-2, ~2×);
- **shock position**: steepest density-gradient face vs $r_s = D\,t$,
  within **4 cells** (measured offset +1.93 cells);
- **$x$-uniformity** (negative control for spurious $x$-fluxes through
  the axisymmetric source): max relative $x$-span of $\rho$/$v$/$p$
  < **1e-8** (measured 3.2e-12 / 4.8e-12 / 1.9e-12 — roundoff).

## How to run

```bash
cd test/2D/euler/noh-cylindrical
./Q2D.sh solve
python3 -B check.py
```
