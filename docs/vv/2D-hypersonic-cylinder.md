# Hypersonic Cylinder

A Mach 8.1 supersonic flow around a cylindrical body, generating a
detached bow shock upstream of the cylinder stagnation point.  This
test case verifies the accuracy of numerical solutions against
theoretical normal shock relations at the stagnation line, where the
bow shock becomes approximately normal.  The test is particularly
valuable for validating Riemann-solver behaviour in high-speed flows
known to be prone to the generation of numerical instabilities such
as carbuncle phenomena.

**Reference:** normal shock relations (Anderson [1]) and the
Billig (1967) bow-shock standoff correlation [2] — see the
References below.

**Feature exercised:** 2D-plane Euler on an O-grid (wall at j=0,
farfield at j=jmax); carbuncle-prone M=8.1 bow shock captured with
the carbuncle-stable HLLE solver; no geometric source
(`[Q2D-Sources]` omitted).

!!! success "V&V status"

    Wired into CTest as `hypersonic-cylinder` and passes on the
    default mesh.  `check.py` gates (all hard, NaN-safe;
    min(ρ), min(p), min(T) > 0):

    1. **Normal-shock + isentropic-stagnation theory** at the
       stagnation-line probe (M₁ = 8.1, γ = 1.4): total pressure p₀
       within **5 %** and total temperature T₀ within **1 %** of the
       analytical values.
    2. **Symmetry**: the case is mirror-symmetric about the
       stagnation line; the ρ asymmetry between mirrored O-grid cells
       (scaled by the local max) is gated max < **4e-1**, mean <
       **1e-2** (measured 1.9e-1 max — a single shock-band cell — and
       1.8e-3 mean on the committed run).  A carbuncle-type breakdown
       produces O(1) asymmetry over wide regions and blows the mean
       gate immediately.
    3. **Bow-shock standoff sanity**: the standoff distance of the
       density jump (ρ crossing 2ρ∞) along the stagnation line must
       fall within **±40 %** of the Billig (1967) cylinder
       correlation Δ/R = 0.386·exp(4.67/M²) = 0.414 at M = 8.1
       (measured 0.423 on the committed run).

---

## Problem setup

Uniform Mach 8.1 flow impinges on a circular cylinder.  At the
stagnation point on the cylinder surface, a detached bow shock forms
approximately normal to the free-stream direction.  The post-shock
flow decelerates isentropically to the stagnation conditions of the
cylinder.

The analytical solution at the stagnation line is governed by normal
shock relations:

$$M_2 = \sqrt{\frac{M_1^2 + 2/(\gamma-1)}{2\gamma M_1^2/(\gamma-1) - 1}}$$

$$\frac{p_2}{p_1} = 1 + \frac{2\gamma}{\gamma+1}(M_1^2 - 1)$$

$$\frac{T_2}{T_1} = \frac{p_2/p_1}{\rho_2/\rho_1}$$

The test verifies that numerical predictions of total pressure $p_0$
and total temperature $T_0$ behind the shock match theoretical
values calculated from the post-shock Mach number.

This test runs in the **2D planar** Q2D regime (the cylinder
cross-section is 2D Cartesian, no axisymmetric source needed).

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Albada |
| **Riemann solvers tested** | **HLLC, HLLC+, HLLE, SLAU** |

## Grid structure

The mesh is a structured grid around the cylindrical body with 160
cells in the circumferential direction and 80 cells in the radial
direction, extending from the cylinder surface to a far-field
boundary located at 1.5 cylinder radii.

!!! note
    The grid is built to trigger the unstable behaviour of some
    Riemann solvers.  In particular, there is no refinement in the
    shock layer or use of cells with a high aspect ratio to avoid the
    carbuncle phenomenon.

Boundary conditions:

- Inlet: Mach 8.1 free-stream state.
- Cylinder surface: slip wall (symmetry).
- Far-field boundaries: zero-gradient extrapolation.

The test probes are positioned at the stagnation streamline to
measure post-shock total quantities.

## Results and verification

Q2D solutions compared across four Riemann solvers in terms of
non-dimensional density contours flow field and wall trend:

<figure>
  {% include "vv/images/HC-fields.svg" %}
</figure>

<figure>
  {% include "vv/images/HC-slice.svg" %}
</figure>

!!! note "Figure generation"

    HC-fields is a 1×4 multipanel showing density contours from the four
    Riemann solvers (HLLE, HLLC, HLLC+, SLAU). HC-slice overlays the
    stagnation-line density profiles for the four solvers on a single 1D
    plot, demonstrating the carbuncle behaviour discussed below. The
    underlying field data was cross-checked against the MOSE twin
    solver during development (agreement to 7e-13 relative L∞ error for
    each of the four solvers) — a one-off code-to-code study requiring
    the MOSE repository, not gated by this repo's ctest (the runnable
    gates are those listed above).

All Riemann solvers capture the fundamental shock structure, but
exhibit varying levels of accuracy in preserving a smooth post-shock
region.  In particular, the HLLC solver shows significant
oscillations in the shock layer, while HLLC+ and HLLE provide
smoother solutions.  The SLAU solver demonstrates an intermediate
behaviour with moderate oscillations.

The density profile along the cylinder wall confirms the overall
behaviour previously observed.

Quantitative error metrics at the stagnation point compared to the
analytical normal-shock solution:

| Solver | $p_0$ [MPa] | Error $p_0$ [%] | $T_0$ [K] | Error $T_0$ [%] |
|---|---|---|---|---|
| HLLC | 8.53825 | -0.79 | 4081.36 | +0.30 |
| HLLC+ | 8.67777 | +0.83 | 4079.80 | +0.26 |
| HLLE | 8.74553 | +1.62 | 4060.20 | -0.22 |
| SLAU | 8.34013 | -3.09 | 4079.70 | +0.26 |
| **Analytical** | **8.60638** | **—** | **4069.11** | **—** |

Notwithstanding only HLLE and HLLC+ provide smooth post-shock
solutions, all solvers exhibit errors well within engineering
tolerances.

## References

1. J. D. Anderson, *Modern Compressible Flow with Historical
   Perspective*, 3rd ed., McGraw-Hill, 2003 (normal shock relations).
2. F. S. Billig, "Shock-wave shapes around spherical- and
   cylindrical-nosed bodies," *J. Spacecraft Rockets*, 4(6), 1967,
   822–823 (bow-shock standoff correlation).
