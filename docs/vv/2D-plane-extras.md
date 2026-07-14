# 2D Plane Euler — Q2D-extra cases

This page documents three 2D plane Euler validation cases that exist
in the Q2D V&V suite but have no direct counterpart in MOSE
`docs/vv/`. They are wired into CTest at the default mesh and were
validated as part of the v1.2.0 validation suite.

| Test                 | CTest target            | Riemann | Status |
|---|---|---|---|
| Chang Nozzle         | `chang-nozzle`          | SLAU    | PASS, Linf_rel ≤ 7e-13 vs MOSE (K=0 layer; axis row 4-6 %) |
| Prandtl–Meyer        | `prandtl-meyer`         | SLAU    | PASS, Linf_rel ≤ 7e-13 vs MOSE (K=0 layer) |
| Ramp Channel         | `ramp-channel`          | SLAU    | PASS, Linf_rel ≤ 7e-13 vs MOSE (K=0 layer) |

---

## Chang Nozzle

A 2D axisymmetric Euler simulation of the Chang (1996) supersonic
nozzle with single-species Air. The test exercises the Q2D
peculiarity #2 (`Lib_Source_Axisymmetry`) on a pure-air gas
without chemistry coupling, providing an isolated regression for the
axisymmetric source-term wiring before the full multi-species
`rocket-nozzle` cases.

**Reference**: T. K. Chang, _J. Spacecraft Rockets_ **33** (1996)
856–865.

### Setup

- Mesh: extracted from MOSE 3D-wedge K=0 layer
  (`tools/extract_mose_K0.py`).
- Domain: nozzle contour from converging inlet through diverging
  exit.
- BCs: axisymmetry on the axis (face 3, BC `200`); inviscid wall
  on the nozzle contour (face 4, BC `301`); supersonic
  inflow/extrapolation at the inlet and exit.

### Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.9 |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | SLAU |

### Validation

The 2D-axi solution agrees with MOSE (3D-wedge K=0 layer)
at relative tolerance (`Linf_rel ≤ 7e-13`) over the nozzle interior. At the
axis row (j=0) the agreement is relaxed to ~4-6 % `eps_rel`: this is
the **T-1.2D.axis discrepancy** — 
Q2D applies the 1/r axisymmetric source with an `eps = 1e-12` floor,
while MOSE-3D-wedge implicitly takes `dz → 0` at `r → 0`. The two
asymptotic regimes converge under mesh refinement.

---

## Prandtl–Meyer Expansion

A 2D plane Euler simulation of a supersonic flow encountering a
convex corner, producing a centred expansion fan whose properties
are fully characterised by the analytical Prandtl–Meyer function
$\nu(M)$. The case verifies the solver's ability to capture smooth
supersonic expansions away from any shock structure.

**Reference**: L. Prandtl & T. Meyer (1908), classical
expansion-fan theory; see also Anderson, _Modern Compressible
Flow_, McGraw-Hill.

### Setup

- Inlet Mach: 2.0, $p_\infty = 1$, $T_\infty = 1$ (non-dimensional Air).
- Corner deflection: 15° expansion.
- Block: structured single-block mesh covering the corner region
  on each side of the expansion fan.

### Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.9 |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | SLAU |

### Validation

`ctest -R prandtl-meyer` runs the case and `check.py` cross-checks the
2D plane solution against MOSE (3D-wedge K=0 layer),
`Linf_rel ≤ 7e-13` on the full domain. The smooth
expansion-fan structure provides a sanity check that the high-order
MUSCL reconstruction does not introduce spurious oscillations in the
absence of shocks (a complementary regression to the shock-dominated
oblique-shock case).

---

## Ramp Channel

A 2D plane Mach 3 supersonic channel flow over an internal ramp.
The configuration generates an oblique shock that reflects off the
upper wall, producing a regular reflection pattern. The test
verifies the consistency of the oblique-shock relations under
multiple shock encounters within a confined geometry.

**Reference**: classical 2D oblique shock theory and channel-flow
references (Anderson, _Modern Compressible Flow_, McGraw-Hill).

### Setup

- Inlet Mach: 3.0.
- Ramp deflection angle: 10° at the lower wall.
- Channel: rectangular bounded above and below.

### Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.9 |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | SLAU |

### Validation

`ctest -R ramp-channel` runs the case (~3 min wall-time on default
mesh) and `check.py` cross-checks the 2D plane solution against MOSE
(3D-wedge K=0 layer), `Linf_rel ≤ 7e-13` on the full
domain. Q2D and MOSE agree on the reflected shock
structure.

`./Q2D.sh test` runs the multi-Riemann-solver sweep across HLLE, HLLC,
HLLC+, SLAU producing `OUTPUT/field-<scheme>.tec` per solver. `verify.py`
then overlays the four-solver density fields interactively (`--plot`); the
ramp-channel case is verified by interactive inspection rather than a
shipped figure. Field visualization is produced offline with `verify.py`
when matplotlib is available.

---

## Visual figures

The 2D contour figures (Mach contours, pressure contours, density
slices along reference lines) for these three cases are produced
offline with matplotlib + ORION `read_TEC` and are inserted by
the release-engineer environment when matplotlib is available. The
underlying field data is archived per-mesh in
`test/regression_baseline/<case>/`; the V&V
acceptance is enforced by the `check.py` regression
(analytical reference + MOSE cross-check), independent
of matplotlib availability.
