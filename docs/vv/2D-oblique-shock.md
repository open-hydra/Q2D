# Oblique Shock

A Mach 4.0 supersonic flow encountering a wedge ramp at 30°
deflection angle.  This configuration generates a single oblique
shock whose properties are governed by the oblique shock relations
from inviscid shock theory.  The test verifies that numerical
solutions across different Riemann solvers accurately satisfy these
theoretical shock jump conditions.

**Reference:** standard oblique shock relations (_e.g._ J. D.
Anderson, _Modern Compressible Flow with Historical Perspective_,
3rd ed., McGraw-Hill, 2003).

**Feature exercised:** 2D-plane Euler on a 2-block grid (multi-block
BC wiring); attached oblique shock vs θ–β–M theory; Riemann-solver
robustness sweep (HLLC, HLLC+, HLLE, SLAU) via `./Q2D.sh test`.

!!! success "V&V status"

    Wired into CTest as `oblique-shock` and passes on the default
    mesh.  `check.py` gates (all hard, NaN-safe):

    - sanity: min(ρ), min(p), min(T) > 0 in every block;
    - convergence sentinel: final density residual < **1e-3**
      (measured 1.5e-5 on the committed 3000-iteration run);
    - probe-window guard: the hard-coded probe index windows must lie
      inside the block cell dimensions (a re-meshed case fails
      loudly, not silently);
    - **θ–β–M theory**: p₂/p₁, T₂/T₁, ρ₂/ρ₁ and M₂ within **1.5 %**
      of the oblique-shock relations, with the deflection δ measured
      from the Block-2 wall nodes (HLLE-only single-Riemann ctest
      run, per the MOSE `verify.py` probe layout).

    Quantitative shock-jump tables across HLLC / HLLC+ / HLLE / SLAU
    are produced via `verify.py` matplotlib from Q2D's
    `OUTPUT/field-<scheme>.tec`.

---

## Problem setup

Uniform Mach 4 flow enters from the left at unit pressure and
temperature.  The flow encounters a sharp ramp deflection angle
$\delta$.  The oblique shock is uniquely determined by the shock
angle $\beta$ through the theta-beta-Mach relation:

$$\tan(\delta) = 2\cot(\beta) \frac{\sin^2(\beta) - M_1^{-2} \sin^2(\delta)}{\sin(2\beta) + M_1^2 \sin^2(\delta)}$$

Post-shock properties (pressure ratio, temperature ratio, density
ratio, exit Mach) are uniform.

This test runs in the **2D planar** Q2D regime.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.8 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| **Riemann solvers tested** | **HLLC, HLLC+, HLLE, SLAU** |

## Grid structure

The mesh is composed of two structured blocks:

- **Block 1** (pre-shock): $48 \times 16$ cells over
  $[0, 0.6] \times [0, 0.2]$.
- **Block 2** (post-shock): $240 \times 64$ cells over
  $[0, 3.0] \times [0.2, 1.0]$.

!!! note
    The grid is built to trigger the unstable behaviour of some
    Riemann solvers.

Boundary conditions:

- Inlet: Riemann invariant extrapolation with Mach 4 state.
- Wall: slip wall (symmetry).
- Upper boundary: zero-gradient extrapolation.
- Exit: zero-gradient extrapolation.

## Results and verification

Comparison of numerical shock jump conditions across all Riemann
solvers:

<figure>
  {% include "vv/images/OS-fields.svg" %}
</figure>

Theoretically, a uniform flow should exist behind the oblique shock
wave.  However, most Riemann solvers except HLLE exhibit
instabilities that produce an undesirable zebra-like oscillatory
pattern in the post-shock region.  The HLLC+ solver, which is
explicitly designed to suppress such numerical artefacts, delivers
the cleanest solution.

The transversal density profile normalised with the analytical
solution is plotted along a vertical line in the post-shock region:

<figure>
  {% include "vv/images/OS-slice.svg" %}
</figure>

The HLLE solver produces a quite uniform flow with density matching
the analytical solution very closely.  In contrast, HLLC exhibits
errors up to 10 % compared to the theoretical value, while HLLC+
remains slightly unstable with oscillations of approximately 2 %
around the analytical density.

!!! note "V&V acceptance"

    The ctest acceptance is the analytical θ–β–M gate set above
    (1.5 % on the shock-jump ratios and post-shock Mach). Field data
    from all four solvers (HLLE, HLLC, HLLC+, SLAU) demonstrate the
    shock-jump conditions and post-shock density profiles.
