# Planar Nozzle — Q2D Area Variation

A converging–diverging nozzle run in the **Q2D area-variation**
regime: the out-of-plane thickness $h(x)$ (`Cross_Section`,
`[Q2D-Sources] active = area_variation`) plays the role of the
quasi-1D area $A(x)$.  For $y$-uniform flow the $h$-weighted
continuity reduces to $\mathrm{d}(h\rho u)/\mathrm{d}x = 0$, so the
converged inviscid solution must match the analytic isentropic
area law (Shapiro).  This is the flagship validation of the
area-variation source term (see
[Q2D Area-Variation Derivation](../theory/area-variation-derivation.md)).

Two ctests share the configuration: `planar-nozzle` (single grid) and
`planar-nozzle-mg2` (identical case plus 2-level multigrid).

**Feature exercised:** Q2D area-variation geometric source
(`[Q2D-Sources] active = area_variation`, Definition A
$\boldsymbol{V}\cdot\boldsymbol{G}$ coupling with
$\boldsymbol{G} = -\nabla A / A$); BC `401` stagnation inflow,
BC `400` extrapolation exit, BC `300` symmetry walls; MUSCL/vanleer +
HLLC with IRS.  The `-mg2` variant additionally exercises
`[Q2D-Multigrid] levels = 2` with the `Cross_Section` fine→coarse R1
restriction (`Coarse_Cross_Section`).

**Reference:** A. H. Shapiro, _The Dynamics and Thermodynamics of
Compressible Fluid Flow_, Ronald Press, 1953 (quasi-1D isentropic
area law).

## Problem setup

Constant-$\gamma$ air ($\gamma = 1.4$).  The area law (identical in
`gen_case.py` and `check.py`) is a smooth cosine blend with
$\mathrm{d}A/\mathrm{d}x = 0$ at the throat:

| Parameter | Value |
|---|---|
| Length $L$ | 1.0 |
| Throat position $x_t$ | 0.40 |
| $A_t$ / $A_\text{in}$ / $A_\text{ex}$ | 1.0 / 2.0 / 2.5 |
| Grid | 200 × 4 cells (the $y$-direction is passive) |
| BCs | face 1 inflow `401` ($T_0$, $p_0$); face 2 extrapolation `400`; faces 3–4 symmetry `300` |

The subsonic inflow chokes at the throat and expands supersonically
to the exit: the analytic Mach distribution follows
$A(x)/A^*$ through the isentropic area–Mach relation (subsonic branch
upstream of $x_t$, supersonic downstream).

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3, local time stepping |
| CFL | 0.5 (rise threshold 50) |
| IRS β | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | HLLC |

## Results

<figure>
  {% include "vv/images/PN-Mach.svg" %}
</figure>

## Acceptance gates — `planar-nozzle`

`check.py` compares the converged Mach distribution
($y$-column means of `OUTPUT/field.tec`) against the analytic
area law, using $\gamma$ and $R$ read from the field itself
(isolates discretisation error):

- positivity: min($\rho$), min($p$), min($T$) > 0 (hard, NaN-safe);
- exit Mach vs analytic: rel err < **2.0e-3**;
- smooth region ($|x-x_t| > 0.05$): max rel err < **5.0e-3**;
- whole domain (incl. throat corner): max rel err < **1.5e-2**.

The full grid-convergence study (Richardson/GCI) is in the case's
`verify.py`.

## Acceptance gates — `planar-nozzle-mg2` (multigrid)

Same nozzle with `[Q2D-Multigrid] levels = 2`.  The coarse level
receives `Cross_Section` restricted fine→coarse via R1 nodal
sub-sampling (`Coarse_Cross_Section`,
$A_c(I,J) = A_f(2I,2J)$).  `check.py` is a three-layer acceptance
that verifies the *interpretation* of the feature (integration test
for the B-10 coarse-restriction fix):

1. **Convergence sentinel** — final density residual < **1.0e-6**
   (measured 9.9e-11 on the committed run, the solver's own 1e-10
   stopping floor).  With the original B-10 bug (`Cross_Section`
   re-read with coarse dims → scrambled re-indexing) the coarse level
   carries a meaningless geometric vector $G = -\nabla A / A$ and the
   solver stalls or diverges.
2. **Physics** — the converged fine field matches the *same* analytic
   area law within the *same* tolerances as the MGL=1 baseline above
   (exit 2.0e-3 / smooth 5.0e-3 / whole 1.5e-2): multigrid must not
   corrupt the area-variation physics.
3. **Multigrid consistency** — `OUTPUT/field-level2.tec` must exist,
   coarse NI must be exactly fine NI/2, and the coarse Mach must
   reproduce the pairwise-restricted fine Mach to
   max rel diff < **5.0e-3** (observed ~4.5e-4 on the committed run;
   one decade of headroom that a scrambled coarse `Cross_Section`
   blows through immediately).

## How to run

```bash
cd test/2D/euler/planar-nozzle        # or planar-nozzle-mg2
./Q2D.sh solve
python3 -B check.py

# via ctest
ctest --test-dir build -R 'planar-nozzle'
```
