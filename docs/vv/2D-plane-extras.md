# 2D Euler — Q2D-extra cases

This page documents five 2D Euler validation cases that exist in the
Q2D V&V suite but have no direct counterpart in MOSE `docs/vv/`.
Despite the historical page name, they are not all 2D-plane: the
Chang nozzle pair runs in the **2D axisymmetric** regime; the
Prandtl–Meyer and ramp-channel cases are 2D plane.  All are wired
into CTest at the default mesh.

| Test | CTest target | Regime | Primary gate |
|---|---|---|---|
| Chang Nozzle | `chang-nozzle` | 2D axi | wall/axis Mach vs Chang 1996 digitized data (avg < 5e-2, max < 1e-1) |
| Chang Nozzle + MG2 | `chang-nozzle-mg2` | 2D axi | same Chang gates + coarse-restriction consistency |
| Prandtl–Meyer | `prandtl-meyer` | 2D plane | exit-probe Mach vs analytic PM function (< 0.5 %) |
| Prandtl–Meyer line probe | `prandtl-meyer-lineprobe` | 2D plane | line-probe output reproduces the field (atol 1e-6 + rtol 1e-3) |
| Ramp Channel | `ramp-channel` | 2D plane | θ–β–M shock-jump ratios (1–2 %) + shock-detector audit |

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

**Feature exercised:** **2D axisymmetric** geometric source
(`Lib_Source_Axisymmetry`, peculiarity #2) on single-species air;
BC `200` axisymmetry axis + BC `301` inviscid wall; MUSCL/vanleer +
SLAU.

### Setup

- Mesh: 190 × 40 structured block ($j$ is the radial direction in the
  axisymmetric regime; K=1 collapsed).
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

`check.py` gates (all hard, NaN-safe):

- min($\rho$), min($p$), min($T$) > 0;
- convergence sentinel: final density residual < **1e-6** (measured
  1.0e-10 on the committed run — the solver's own stopping criterion
  drives it to ~1e-10);
- Mach distribution along the **wall** vs the Chang 1996 digitized
  reference `reference/ref2.dat` and along the **axis** vs
  `reference/ref1.dat`: avg $\varepsilon_\text{rel}$ < **5e-2**,
  max $\varepsilon_\text{rel}$ < **1e-1**.

The relaxed max tolerance covers the axis-row cells: Q2D applies the
$1/r$ axisymmetric source with an $\varepsilon = 10^{-12}$ floor,
while the MOSE 3D-wedge equivalent implicitly takes $dz \to 0$ at
$r \to 0$ — two distinct approximations of axisymmetric physics that
differ ~4 % at the axis row (the **T-1.2D.axis discrepancy**) and
converge under mesh refinement.

### Multigrid variant (`chang-nozzle-mg2`)

Same Chang nozzle run with 2-level multigrid
(`[Q2D-Multigrid] levels = 2`); the coarse-grid `bc2.txt` is
generated via `tools/gen_coarse_bc.py` (factor-2 restriction).

**Feature exercised:** **2-level multigrid** ×  the 2D axisymmetric
source — the coarse level must carry a consistent restriction of the
axisymmetric geometry.

Mirror of the planar-nozzle-mg2 three-layer acceptance
(see [Planar Nozzle](2D-planar-nozzle.md)):

1. **Convergence sentinel** (hard): final density residual <
   **1e-6** (measured 1.0e-10 on the committed run).  A broken
   coarse-level operator stalls/diverges instead of converging.
2. **Physics**: the converged fine field matches the *same*
   Chang 1996 digitized references within the *same* tolerances as
   `chang-nozzle` — multigrid must not corrupt the physics.
3. **Multigrid consistency**: `OUTPUT/field-level2.tec` must exist
   (hard: the MGL=2 path did not run otherwise), coarse NI/NJ must be
   exactly fine NI/2, NJ/2, and the coarse-level Mach must reproduce
   the 2×2-pairwise-restricted fine Mach: max rel diff < **2e-1**,
   mean < **5e-2** (measured 8.1e-2 max — at the axis/throat
   singularity row — and 5.2e-3 mean on the committed run; a
   scrambled coarse restriction is off by O(1) everywhere).

<figure>
  {% include "vv/images/CN-mg2.svg" %}
</figure>

---

## Prandtl–Meyer Expansion

A 2D plane Euler simulation of a supersonic flow encountering a
convex corner, producing a centred expansion fan whose properties
are fully characterised by the analytical Prandtl–Meyer function
$\nu(M)$. The case verifies the solver's ability to capture smooth
supersonic expansions away from any shock structure.

**Reference**: the classical expansion-fan theory of Prandtl and
Meyer (1908); modern treatment in J. D. Anderson, _Modern
Compressible Flow with Historical Perspective_, 3rd ed.,
McGraw-Hill, 2003, ch. 4.

**Feature exercised:** 2D-plane smooth supersonic expansion;
point-probe output (`[Q2D-Probes]` exit probe, `OUTPUT/exit.txt`);
MUSCL/vanleer + SLAU.  The `-lineprobe` variant additionally
exercises the Q2D-only **line-probe** registry path
(`Q2D_IO_LineProbes`).

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

`ctest -R prandtl-meyer` runs the case; `check.py` reads the
converged Mach at the exit-probe location (`OUTPUT/exit.txt`,
`[Q2D-Probes]` in `input.ini`) and gates it against the analytical
isentropic expansion result:
$|M_\text{CFD} - M_\text{analytical}| / M_\text{analytical}$ <
**0.5 %**, plus the field sanity min($\rho$, $p$, $T$) > 0
everywhere (hard, NaN-safe).  The smooth expansion-fan structure
provides a sanity check that the high-order MUSCL reconstruction does
not introduce spurious oscillations in the absence of shocks (a
complementary regression to the shock-dominated oblique-shock case).

### Line-probe variant (`prandtl-meyer-lineprobe`)

The same expansion fan with a single $j$-fixed line probe
(`[midrow]`, `fixed-index = 25`) sampling cell row $j=25$ along $i$.
It exercises the Q2D-only line-probe feature (`Q2D_IO_LineProbes`;
MOSE has no line probes) end-to-end through the registry-driven path:
`Scan_Ini` count → `Register_LineProbes` bind → `Setup_LineProbes`
read → `Write_LineProbe_Tecplot` output.

The reference is **self-consistency** (no invented thresholds): a
correct line probe must reproduce, at the sampled grid line, the same
primitive field the main writer puts in `OUTPUT/field.tec` — both
write the identical `P(:, i, j)` slice.  `check.py` asserts:

- **(struct)** `OUTPUT/midrow.tec` exists with the expected header
  (`rho(1), u, v, p, T, M`) and degenerate-strip zone geometry (J=2);
- **(geometry)** the written strip nodes are bit-exactly the mesh
  node rows $j = 24/25$ read from the case's own `INPUT/ic.tec`
  (tolerance 1e-12) — the (x, y) really trace that grid line;
- **(physical)** along the sampled line the expansion fan
  accelerates the flow: Mach anchored at the inflow $M = 2$ (±0.1) at
  the line start, ending more than 1.0 above it; $p$, $\rho$, $T$
  all dropping;
- **(self-consistency)** line-probe values vs field cells at row
  $j=25$: $|lp - \text{field}| < 10^{-6} + 10^{-3}\,|\text{field}|$
  per variable ($\rho, u, v, p, T$).  The residual offset between the
  probe's last write and the field write is bounded by the converged
  density residual (< 1e-8 here); a real extraction bug — wrong row,
  transposed index, wrong variable map — is off by O(1)–O(100 %).

<figure>
  {% include "vv/images/PM-lineprobe.svg" %}
</figure>

---

## Ramp Channel

A 2D plane Mach 3 supersonic channel flow over an internal ramp.
The oblique shock from the ramp tip (at $x = 0.1$) reflects off the
upper boundary; the ramp shoulder ($x \approx 0.187$) sheds an
expansion fan. The test verifies the oblique-shock relations inside a
confined geometry and audits the shock-detector field.

**Reference**: classical 2D oblique shock theory (J. D. Anderson,
_Modern Compressible Flow with Historical Perspective_, 3rd ed.,
McGraw-Hill, 2003).

**Feature exercised:** 2D-plane oblique-shock capture with the
**shock detector** (`SD_Tramel`) coupled to IRS
(`irs-beta = 0.5`); MUSCL/vanleer + SLAU; multi-Riemann-solver sweep
(HLLE, HLLC, HLLC+, SLAU) via `./Q2D.sh test`.

### Setup

- Inlet Mach: 3.0.
- Ramp deflection angle: 30° ramp geometry at the lower wall
  (the deflection δ used by the θ–β–M gate is measured from the
  wall nodes themselves).
- Channel: rectangular bounded above and below.

### Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.9 |
| IRS β | 0.5 (with shock detector) |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | SLAU |

### Validation

`ctest -R ramp-channel` runs the case (~3 min wall-time on the
default mesh); `check.py` gates (all hard):

1. **Oblique-shock θ–β–M theory** ($M_1 = 3$, $\gamma = 1.4$, ramp δ
   measured from the wall nodes): post/pre-shock ratios averaged over
   a pre-shock probe and a post-shock plateau probe (inside the
   uniform wedge region, downstream of the shock foot and upstream of
   the shoulder expansion).  Tolerances calibrated at ~2× the
   measured errors (2026-07-15, gfortran TESTING build: $p_2/p_1$
   −0.002 %, $T_2/T_1$ +0.95 %, $\rho_2/\rho_1$ −0.94 %, $M_2$
   −0.94 %): $|p_2/p_1|$ err < **1 %**; $|T_2/T_1|$, $|\rho_2/\rho_1|$,
   $|M_2|$ err < **2 %**.
2. **Shock-detector audit** (`SD_Tramel`, β=0 in shock cells, β=1 in
   smooth flow): the β field differentiates; ≥ 1 % of cells detected;
   detected cells localized downstream of the ramp tip (≥ 90 %;
   observed 99.3 %); smooth pre-shock block β mean > 0.99; and a
   negative control — an "expected shock zone" placed upstream of the
   tip must be ~empty (≤ 2 %).
3. **Convergence + IRS stability** (`irs-beta = 0.5` active with the
   shock detector): residual history finite and final density
   residual < **1e-3** (measured final 2.49e-4 after 10000
   iterations, peak 7.98e-3 at shock formation ~iteration 99; the
   threshold is ~4× the measured final — a destabilized IRS diverges
   instead).
4. **NaN-safe sanity**: $\rho$, $p$, $T$, $g$, $R$ > 0; $u$, $v$
   finite.

`./Q2D.sh test` additionally runs the multi-Riemann-solver sweep
across HLLE, HLLC, HLLC+, SLAU producing `OUTPUT/field-<scheme>.tec`
per solver; `verify.py` overlays the four-solver density fields
interactively (`--plot`).  The sweep is inspected interactively
rather than gated.

---

## Visual figures

The 2D contour figures (Mach contours, pressure contours, density
slices along reference lines) for these cases are produced offline
with matplotlib + ORION `read_TEC` via each case's `verify.py` when
matplotlib is available.  The V&V acceptance is enforced by the
`check.py` gates described above (analytical / digitized-reference
comparisons), independent of matplotlib availability.
