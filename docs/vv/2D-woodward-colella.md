# Woodward-Colella Supersonic Step

A Mach 3.0 supersonic flow over a forward step geometry.  This is a
classic test case for compressible flow solvers that exercises the
scheme's ability to capture oblique shocks, expansion fans, and their
interactions in a multi-dimensional setting.  The solution features
a reflected shock, expansion fan, and contact surface, providing a
stringent validation of shock-capturing capabilities.

**Reference:** P. Woodward & P. Colella, _J. Comput. Phys._ **54**
(1984) 115–173.

**Feature exercised:** 2D-plane Euler, **time-accurate** RK3 run to a
deterministic endpoint ($t = 4.0$); 2-block grid (multi-block BC
wiring); HLLE shock capture with shock/expansion interactions.

!!! success "V&V status"

    Wired into CTest as `supersonic-forward-step` and passes on the
    default mesh.  `check.py` gates (all hard):

    1. **Density trace vs committed baseline**: ρ along block-2 cell
       row j=23 (nearest y=0.5, through the bow-shock/Mach-stem
       region) vs `reference/trace_baseline.dat` (provenance in
       `reference/README.md`: gfortran TESTING build, main @ 8c60e91,
       2026-07-15).  Max relative difference per point < **5 %** —
       single-thread reruns of the same binary are bit-identical
       (measured 0.0), so 5 % is pure cross-compiler headroom.  The
       trace's x cell centers must match the baseline's to 1e-9
       (mesh guard).
    2. **Bow-shock position**: the shock x-position along the same
       row, located by the steepest density gradient, within
       **3 cells** (|dx| ≤ 0.0375) of the baseline's (measured
       x_shock = 0.425 at t = 4.0).
    3. **NaN-safe positivity** on every block: ρ, p, T, g, R > 0;
       u, v finite.

    The thermo-data setup is via the FLINT `phase.txt`+`thermo.dat`
    symlinks to `test/common/nondim-Air/`.

---

## Problem setup

A uniform Mach 3 flow enters from the left boundary and faces a
forward step obstacle at the lower boundary.  This produces a
classical shock-deflection pattern:

- Primary oblique shock deflecting around the step corner.
- Second shock reflecting from the upper boundary.
- Expansion fan originating from the step corner.
- Contact surface separating regions of different entropy.

This test runs in the **2D planar** Q2D regime (no axisymmetric or
area-variation source); it provides a direct comparison with MOSE
(2D plane of a 3D simulation) for verification of the inviscid Euler
core.

## Numerical setup

| Parameter | Value |
|---|---|
| Domain | $[0, 3] \times [0, 1]$ |
| Final time | $t = 4.0$ |
| Mach number | 3.0 |
| Time scheme | RK3 |
| CFL | 0.8 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLE |

## Grid structure

The mesh is composed of two structured blocks:

- **Block 1** (step region): $48 \times 16$ cells over
  $[0, 0.6] \times [0, 0.2]$.
- **Block 2** (downstream): $240 \times 64$ cells over
  $[0, 3.0] \times [0.2, 1.0]$.

Boundary conditions:

- Inlet: Riemann invariant extrapolation with Mach 3 state.
- Step wall: slip wall (symmetry).
- Upper boundary: symmetry.
- Exit: zero-gradient extrapolation.

## Results and verification

Q2D solutions compared against an OpenFOAM (`rhoCentralFoam`)
reference on the baseline mesh (the OpenFOAM data lives in
`test/2D/euler/supersonic-forward-step/reference/OpenFOAM/`; this
cross-code comparison is **manual-only** — the ctest gate is the
committed density-trace baseline described above):

<figure>
  {% include "vv/images/WC-fields.svg" %}
</figure>

Density profile cut at $y = 0.5$ demonstrating sharp resolution of
shock structures:

<figure>
  {% include "vv/images/WC-slice.svg" %}
</figure>

The solution accurately captures all major flow features including
the oblique shock angle, the shock interaction with the upper
boundary, and the expansion fan structure.  Pressure and density
jumps are well-resolved with minimal oscillation.

!!! note "Figure provenance"

    The WC-fields (3-panel p, T, u contours) and WC-slice density-line
    SVG figures are Q2D-side renderings. The ctest acceptance record
    is the committed density-trace baseline + bow-shock-position gate
    described in the V&V status box above.
