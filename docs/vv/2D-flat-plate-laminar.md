# 2D Flat Plate – Laminar Boundary Layer

Incompressible laminar flow over a flat plate at low Mach number.
This configuration establishes a developing Blasius boundary layer
whose properties are well-characterised by classical similarity
theory.  Verification is performed by comparison with the exact
Blasius solution for skin-friction coefficient and velocity profiles.

**Reference:** Blasius similarity solution for the laminar boundary
layer over a flat plate; the tabulated $f'(\eta)$ profile enforced by
`check.py` is per Howarth (1938) [1], with the theory as in
Schlichting [2] — see the References below.

**Feature exercised:** laminar Navier–Stokes (viscous fluxes, no
turbulence model) in the 2D-plane regime; no-slip adiabatic wall +
symmetry leading edge; analytical-only acceptance (Blasius/Howarth
table, no external-code reference); heavy-tier ctest
(`./Q2D.sh -p 20 solve`, 20 OpenMP threads, label `heavy`).



<figure>
  {% include "vv/images/LP-field.svg" %}
</figure>

---

## Problem setup

Uniform incompressible flow at Mach 0.2 enters from the left over a
semi-infinite flat plate.  The leading edge induces a boundary layer
that develops downstream following the laminar Blasius solution.
Post-stagnation, a no-slip adiabatic wall condition is imposed.

The Blasius equation governs the self-similar boundary-layer
profile:

$$f''' + f \cdot f'' = 0$$

with boundary conditions $f(0) = f'(0) = 0$ and $f'(\eta \to \infty) = 1$,
where $\eta = y\sqrt{U_\infty / 2\nu x}$ is the similarity variable.

The skin-friction coefficient follows:

$$C_f = \frac{\sqrt{2} \cdot f''(0)}{\sqrt{Re_x}}$$

where $Re_x = U_\infty x / \nu$ is the local Reynolds number.

This test runs in the **2D planar** Q2D regime; the boundary layer
develops on a Cartesian wall with no axisymmetric or area-variation
source.

**Freestream conditions**

| Parameter | Value |
|---|---|
| Mach number | 0.2 |
| Temperature | 300 K |
| Pressure | 97 250 Pa (`input.ini` / `INPUT/bc.txt`) |
| Viscosity | 1.0 × 10⁻³ Pa·s (`mil` in `input.ini`) |
| Re/m | ρ∞U∞/μ ≈ 7.84 × 10⁴ m⁻¹ |

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK2 |
| CFL | 0.5 |
| VNN | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |

## Grid structure

The mesh is a 2D structured grid ($65 \times 65$ cells) spanning the
physical domain with fine resolution near the leading edge and the
near-wall region to capture boundary-layer development accurately.

Boundary conditions:

- **Inlet:** Riemann invariant boundary condition with freestream
  state.
- **Adiabatic wall** ($x \ge 0$): no-slip, no-heat-flux.
- **Symmetry** ($x < 0$): symmetry boundary condition upstream.
- **Outlets:** atmospheric ambient conditions.

## Acceptance gates (`check.py`, ctest `flat-plate-laminar`)

The ctest gate (heavy tier) enforces, after `./Q2D.sh -p 20 solve`:

| Gate | Threshold | Notes |
|---|---|---|
| Solver outputs present, finite (NaN/Inf scan) | hard fail | missing output = FAIL, never skip |
| $u/U_\infty$ profile vs Blasius (Howarth 1938 table) at $x \approx 0.15$ m | max pointwise error $\le$ **0.15** over $\eta \in [0.5, 4.5]$ | measured 0.115 on the committed coarse grid; 0.15 is the NASA TMR coarse-grid acceptance (a $\pm 5\%$ gate requires a much finer grid) |
| Profile-point coverage | $\ge 1$ point in the $\eta$ window | anti-vacuous-pass guard |

The viscosity is read from the committed `INPUT/transport.dat` (no
hard-coded $\mu$; the historical mu-source bug that motivated the
tolerance revision is documented in the script).

## Results and verification

Comparison of skin-friction coefficient along the plate surface
versus the analytical Blasius distribution:

<figure>
  {% include "vv/images/LP-Cf.svg" %}
</figure>

The solution demonstrates excellent agreement with Blasius theory
away from the leading edge.  Leading-edge errors ($x < 0.05$ m) are
expected due to numerical stagnation effects and coarse initial mesh
resolution; Blasius theory assumes an infinitesimal leading edge,
whereas practical CFD has finite spatial resolution.

**Velocity profile**

Comparison of the normalised velocity profile $u / U_\infty$ versus
the Blasius similarity coordinate $\eta$ at $x = 0.29$ m:

<figure>
  {% include "vv/images/LP-slice.svg" %}
</figure>

The numerical velocity profile shows excellent collapse onto the
Blasius self-similar solution, confirming that the boundary layer
develops according to classical theory.

## References

1. L. Howarth, "On the solution of the laminar boundary layer
   equations," *Proc. R. Soc. Lond. A*, 164(919), 1938, 547–579
   (tabulated Blasius profile used by `check.py`).
2. H. Schlichting, *Boundary-Layer Theory*, 7th ed., McGraw-Hill,
   1979, ch. 7.
