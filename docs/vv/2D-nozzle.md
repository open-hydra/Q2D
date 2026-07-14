# Rocket Nozzle – Multi-Chemistry

A 2D axisymmetric Euler simulation of a rocket nozzle operating on
RP-1 / LOX propellants mocking the Falcon 9 Merlin1D.  Three
configurations are tested: **frozen composition** (fixed species
frozen at chamber state), **equilibrium** (species adjust
instantaneously to maintain thermochemical equilibrium), and
**finite-rate** (species composition evolves according to reaction
kinetics).  Area-weighted quasi-1D averages extracted from 2D CFD
solutions are compared against 1D CEA (Chemical Equilibrium with
Applications) references at discrete area-ratio stations.  The tests
verify multi-species solver accuracy across different thermochemical
assumptions.

**Reference:** CEA (Chemical Equilibrium with Applications), NASA
Glenn Research Center.

!!! success "V&V status"

    All three sub-cases (frozen, equilibrium, finite-rate) are wired
    into CTest as `rocket-nozzle-{frozen,equilibrium,finite-rate}` and
    pass on the default mesh.  Per-case `check.py` validates the 2D-axi
    solution against MOSE (3D-wedge K=0 layer, multi-zone cross-check
    comparison via `tools/diff_mose_q2d_field.py`) and against the CEA
    1D equilibrium reference for the equilibrium and finite-rate
    configurations.

    Residuals at termination:

    - frozen: 3185 iter, 54 s, residual 9.99e-11 (full convergence)
    - equilibrium: 10000 iter, 5.7 min, residual ~6e-2 (chemistry-source
      noise floor; same behaviour observed in MOSE)
    - finite-rate: 2000 iter, 15.7 min, residual ~8.3e-5

    `Linf_rel` vs MOSE (K=0 layer) ≤ 2 % over the nozzle interior; ~4-6 %
    at the axis row, where the 1/r source term dominates.

!!! note "Q2D regime"

    This test runs in the **2D axisymmetric** Q2D regime, exercising
    the cylindric-residual source `Lib_Source_Axisymmetry`
    (peculiarity #2) and the `Stress_Vector_Axi` viscous closure
    (peculiarity #3).  The axis of symmetry coincides with the lower
    block boundary (face 3) where the BC `200` (axisymmetry) is
    applied.  The Borda-flow choked-nozzle BC (peculiarity #4,
    sentinel `BC_Mach == -20`) is also available for inflow boundary
    treatment in this geometry.

<figure>
  {% include "vv/images/Merlin-field.svg" %}
  Contours of Mach number from the 2D frozen-flow simulation of the Merlin nozzle.
</figure>

---

## Problem setup

Combustion products of RP-1 burned with liquid oxygen at an
oxidiser-to-fuel ratio of 2.35 are expanded from a stagnation
pressure of 108 bar and a stagnation temperature of 3644 K through
the nozzle contour.

### Grid structure

A single structured block with $120 \times 40$ cells (axial × radial)
covers the full nozzle from the converging inlet to the nozzle exit.
The lower boundary is the axis of symmetry and the upper boundary
follows the Merlin ideal contour.

Boundary conditions (same for all configurations):

- Inlet: CEA chamber equilibrium state (subsonic inflow) with 13
  species from RP-1 / LOX combustion products.
- Exit: zero-gradient supersonic extrapolation.
- Axis: axisymmetry (BC type `200`).
- Wall: symmetry (inviscid).

## Test configurations

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.9 |
| IRS β | 0.3 |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | SLAU |

- **Frozen composition** — species mass fractions are fixed at the
  chamber-equilibrium state and do not react during the expansion.
- **Equilibrium chemistry** — species composition adjusts
  instantaneously at each point to maintain thermochemical
  equilibrium.
- **Finite-rate chemistry** — species composition evolves according
  to reaction kinetics using the CORIA reaction mechanism with the
  NASA9 thermodynamic database.

## Results and verification

The 2D averages from each configuration are compared against 1D CEA
references.  The frozen-composition solver reproduces isentropic
expansion closely, validating the multi-species Euler solver.  The
equilibrium configuration verifies that thermochemistry updates are
correctly integrated with the flow solver across the nozzle
expansion.  Finite-rate chemistry results confirm that the solver
can handle non-equilibrium effects, although the high-temperature,
high-pressure conditions in the nozzle lead to near-equilibrium
behaviour.  Overall, the tests demonstrate accurate multi-species
flow predictions across different thermochemical assumptions in a
rocket nozzle geometry.

<figure>
  {% include "vv/images/Merlin-XY.svg" %}
  Area-weighted 1D averages of Mach number, pressure ratio, and temperature ratio compared against CEA reference stations.
</figure>

!!! note "Figure provenance"

    The Merlin-field and Merlin-XY SVG figures above are Q2D-side
    renderings produced by `tools/gen_vv_plots.py --all-2d`:

    - **Merlin-field**: density contour from the **frozen** sub-case
      (gnuplot pm3d).
    - **Merlin-XY**: 1×3 multiplot of area-weighted axial averages of
      Mach, p/p_in, T/T_in for all three sub-cases (frozen,
      equilibrium, finite-rate). Area weighting is the axi-cylindric
      `r dr dx` integral over the radial column at each axial station,
      enforced by the 2D-axisymmetric geometry of the test case
      (peculiarity #2). Throat Mach numbers extracted:
      frozen ~3.70, equilibrium ~3.49, finite-rate ~3.51.

    The rocket-nozzle 2D-axi solution agrees with MOSE
    (3D-wedge K=0 layer) at `Linf_rel ≤ 2 %` over the interior; ~4-6 %
    at the axis row.
    These 2D-axi figures provide the strict V&V record;
    visual comparison with MOSE figures may show minor differences due to solver specifics.
