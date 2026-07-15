# Verification & Validation

This section documents the Verification & Validation (V&V) test
suite for Q2D.  The purpose is to demonstrate that the solver
produces correct results across a range of physical regimes, by
comparison against analytical solutions, published/experimental
references, or committed deterministic baselines, and to protect
Q2D-specific features (geometric-source regimes, special boundary
conditions, multigrid, probes) with hard regression gates.

## Test suite

**34 tests** are wired into CTest (`test/CMakeLists.txt`), split into
a **fast** tier (run on every push/PR: GitHub CI runs
`ctest -LE heavy`) and a **heavy** tier (RANS/viscous cases needing
many threads and/or hours of wall clock; excluded from hosted CI and
run locally: `ctest -L heavy`).

### 1D cases (fast)

| Test | Physics | Verification | Feature exercised | Ref |
|---|---|---|---|---|
| [Sod79](1D.md#sod-shock-tube) | Shock, contact, rarefaction | Analytical (exact Riemann) | 1D strip layout; MUSCL/MC; SLAU | Sod (1978) |
| [Ein91](1D.md#einfeldt-double-rarefaction) | Near-vacuum double rarefaction | Analytical (exact Riemann) | positivity near vacuum; vanleer; HLLE | Einfeldt et al. (1991) |
| [Tor99](1D.md#toro-test-case-3) | Transonic rarefaction + shock | Analytical (exact Riemann) | entropy fix; vanleer; HLLE++ | Toro (1999) |
| [Noh87](1D.md#noh-implosion-problem) | Infinite-strength shocks | Analytical (Noh) | first-order; HLLE; $p \sim 10^{-15}$ IC | Noh (1987) |
| [ShuOsher-LIMO3](1D-shuosher.md) | Shock/entropy-wave interaction | Committed limiter baselines (TV gates + negative control) | **LIMO3 third-order limiter** (Q2D-only); HLLC | Shu & Osher (1989) |
| [Sod79-convergence](convergence-table-1D.md#convergence-gates-ctest) | ×1–×8 refinement ladder | Analytic exact solution, observed L1 order ≥ 0.70 | `verify.py --gate` convergence gating | Toro (2009) |
| [Ein91-convergence](convergence-table-1D.md#convergence-gates-ctest) | ×1–×8 refinement ladder | observed L1 order ≥ 0.78 | `verify.py --gate` convergence gating | Toro (2009) |
| [Tor99-convergence](convergence-table-1D.md#convergence-gates-ctest) | ×1–×8 refinement ladder | observed L1 order ≥ 0.75 | `verify.py --gate` convergence gating | Toro (2009) |
| [Noh87-convergence](convergence-table-1D.md#convergence-gates-ctest) | ×1–×8 refinement ladder | observed L1 order ≥ 0.85 | `verify.py --gate` convergence gating | Noh (1987) |

### 2D Euler cases (fast)

| Test | Regime | Physics | Verification | Feature exercised | Ref |
|---|---|---|---|---|---|
| [chang-nozzle](2D-plane-extras.md#chang-nozzle) | 2D axi | transonic nozzle expansion | Chang 1996 digitized wall/axis Mach | axisymmetric source (`Lib_Source_Axisymmetry`); BC 200/301 | Chang (1996) |
| [chang-nozzle-mg2](2D-plane-extras.md#multigrid-variant-chang-nozzle-mg2) | 2D axi | same + 2-level multigrid | same Chang gates + coarse-restriction consistency | **multigrid** × axisymmetric source; `gen_coarse_bc.py` | Chang (1996) |
| [prandtl-meyer](2D-plane-extras.md#prandtlmeyer-expansion) | 2D plane | M=2 expansion fan | Analytical PM function (< 0.5 %) | smooth-expansion MUSCL; point probes | Prandtl–Meyer (1908) |
| [prandtl-meyer-lineprobe](2D-plane-extras.md#line-probe-variant-prandtl-meyer-lineprobe) | 2D plane | same + line probe | field self-consistency (atol 1e-6 + rtol 1e-3) | **line probes** (`Q2D_IO_LineProbes`, registry path) | Prandtl–Meyer (1908) |
| [planar-nozzle](2D-planar-nozzle.md) | Q2D area-var | choked C-D nozzle | quasi-1D isentropic area law | **area-variation source** (Definition A V·G); BC 401 | Shapiro (1953) |
| [planar-nozzle-mg2](2D-planar-nozzle.md#acceptance-gates-planar-nozzle-mg2-multigrid) | Q2D area-var | same + 2-level multigrid | area law + coarse-restriction consistency | **multigrid** × area variation (`Coarse_Cross_Section` R1) | Shapiro (1953) |
| [borda-channel](2D-borda-channel.md) | 2D plane | choked Borda injection | Fievisohn 2016 eqs. 3.26–3.31 | **BC 420 Borda injector** (`BC_Mach == -20`) | Fievisohn (2016) |
| [borda-channel-negctl](2D-borda-channel.md#negative-control-ctest-borda-channel-negctl) | 2D plane | corrupted-expectation control | exit 0 iff the corrupted gate FAILS | **negative control** (false-green tripwire) | Fievisohn (2016) |
| [hypersonic-cylinder](2D-hypersonic-cylinder.md) | 2D plane | M=8.1 bow shock, carbuncle | normal-shock/stagnation theory + symmetry + Billig standoff | carbuncle-stable HLLE; O-grid | Billig (1967) |
| [noh-cylindrical](2D-noh-cylindrical.md) | 2D axi | cylindrical implosion, infinite-strength shock | Analytical (Noh m=1: D=(γ−1)U₀/2, ρ₂=((γ+1)/(γ−1))², 8 %) + x-uniformity control | **axisymmetric source** (`Lib_Source_Axisymmetry`) vs analytic exact — the only axi-exact case | Noh (1987) |
| [ramp-channel](2D-plane-extras.md#ramp-channel) | 2D plane | M=3 oblique shock reflection | θ–β–M theory (1–2 %) + shock-detector audit | **shock detector** (`SD_Tramel`) + IRS | Anderson |
| [oblique-shock](2D-oblique-shock.md) | 2D plane | M=4 attached oblique shock | θ–β–M theory (1.5 %) | 2-block grid; Riemann-solver sweep | Anderson |
| [supersonic-forward-step](2D-woodward-colella.md) | 2D plane | M=3 forward step, time-accurate | committed density-trace baseline (5 %) + bow-shock position | time-accurate RK3; 2-block | Woodward & Colella (1984) |
| [rocket-nozzle-frozen](2D-nozzle.md) | 2D axi | 13-species frozen chemistry | committed reference field (Linf_rel < 5e-2, L2_rel < 1e-2) + gas-law consistency | multi-species frozen mode | CEA / NASA |
| [rocket-nozzle-equilibrium](2D-nozzle.md) | 2D axi | 13-species CEA equilibrium | same committed-reference gates | per-cell **CEA equilibrium** | CEA / NASA |
| [rocket-nozzle-finite-rate](2D-nozzle.md) | 2D axi | 19-species finite-rate (CORIA) | same committed-reference gates | **finite-rate kinetics** (sdirk4b) | CEA / CORIA |
| [SRM-rocket-motor](2D-srm-rocket-motor.md) | 2D axi | solid-motor chamber pressure | analytic Vieille balance (15 %) | **BC 502 SRM** Vieille burn rate | Vieille law |
| [rotating-channel](2D-rotating-channel.md) | 2D plane | Coriolis geostrophic balance | analytic dp/dy = −2ωρU∞ (5 %) + drift gates | **rotating frame** (Coriolis); **BC 201 periodic** | MOSE coriolis-channel mirror |
| [rotating-channel-rsm](2D-rotating-channel-rsm.md) | 2D plane | frame-invariance of the SSGLRR rotating-frame closures | homogeneous ODE (RK4 reference, per-cell) 2 % + no-G_ij discriminance ×5 | **rotating frame + SSGLRR** (exact Coriolis production G_ij, W_abs C5); time-accurate RK3 | Q2D-original (frame-invariance ODE) |
| [multi-solver-103](2D-multi-solver.md) | 2D plane | quiescent equilibrium (smoke) | bounded drift from IC | **BC 103** coupling dispatcher (smoke test) | — |
| [ssglrr-planar-nozzle-turbulent](2D-ssglrr.md) | Q2D area-var | RSM vs SST pair | code-verification gates G1–G5 | **SSG/LRR-ω RSM** × area variation | — (code verification) |

### Viscous / RANS cases (heavy tier)

| Test | Regime | Physics | Verification | Feature exercised | Ref |
|---|---|---|---|---|---|
| [flat-plate-laminar](2D-flat-plate-laminar.md) | 2D plane | laminar boundary layer | Blasius similarity (Howarth table) | laminar NS; adiabatic wall | Blasius (1908) |
| [flat-plate-turbulent-coarse-kw](2D-flat-plate-turbulent.md) | 2D plane | turbulent boundary layer | NASA TMR SST Cf (`cf-sst.dat`, FUN3D zone), physical units (15 %; measured 4.7 %) | **Menter SST**; RANS 12-arg dispatch | NASA TMR |
| [swbli](2D-swbli.md) | 2D plane | M=5 shock/boundary-layer interaction | Schulein/SU2/Wind-US separation envelope + MOSE-SA Cf | **Spalart–Allmaras**; skip semantics (exit 77) | Schulein (2002) |
| [ssglrr-chang-nozzle-turbulent](2D-ssglrr.md) | 2D axi | RSM vs SST pair | code-verification gates G1–G5 | **SSG/LRR-ω RSM** × axisymmetric source | — (code verification) |

## Coverage matrix

The complete feature x test x reference x docs x figure map — including
the reference-type hierarchy and the honest list of known gaps
(untested Riemann solvers, turbulence variants, BC codes, ...) — is
maintained in the [V&V Coverage Matrix](coverage-matrix.md).

## Running the tests

Every case directory under `test/` contains a `check.py` — the
CTest-driven hard gate (exit 0 pass / exit 1 fail / exit 77 skip
where a legitimate skip condition exists).  In addition, the four 1D
Riemann/Noh cases and the two ssglrr family roots carry a
`verify.py` (standalone study / family gate); for the 1D cases
`verify.py --gate` is itself a ctest (the convergence ladder).

```bash
# single case, by hand
cd test/2D/euler/planar-nozzle
./Q2D.sh solve
python3 -B check.py

# whole fast tier (what CI runs)
ctest --test-dir build -LE heavy --output-on-failure

# heavy tier (local, many-core node)
ctest --test-dir build -L heavy --output-on-failure
```

---

## Execution status

The suite is tiered (see
[Testing](../development/testing.md)):

- **Fast tier** (30 tests, labels without `heavy`): runs on every
  push/PR in GitHub CI via `ctest -LE heavy`.  About **45 min**
  single-threaded serially, ~25 min with `ctest -j`; the longest
  fast test is `rocket-nozzle-finite-rate` at ~22 min single-thread.
- **Heavy tier** (4 tests, label `heavy`): RANS/viscous cases run
  with 20 OpenMP threads and wall-clock budgets up to hours
  (`swbli` has an 18000 s timeout); excluded from hosted CI, run
  locally on many-core nodes with `ctest -L heavy`.

Every test uses `SKIP_RETURN_CODE 77`: a check script that cannot run
its quantitative gate exits 77 and shows up as **Skipped** — never as
a silent "Passed" (false green).  Tests sharing a case directory
(e.g. `borda-channel` + its negative control) are serialized via a
CTest `RESOURCE_LOCK`.

For each case, `check.py` validates the Q2D solution against an
analytical reference (1D Riemann/Noh, area law, θ–β–M,
Coriolis/Vieille balances), published or experimental data
(Chang 1996, Blasius, NASA TMR, Schulein/SU2/Wind-US), or a committed
deterministic baseline (rocket-nozzle reference fields,
forward-step density trace, ShuOsher limiter baselines).  The gate
values are quoted verbatim on each case page and in the check-script
docstrings — they are the single source of truth.

1D mesh-refinement convergence (×1/×2/×4/×8) is gated in CTest and
summarised in [the 1D convergence page](convergence-table-1D.md).

---

## Verification Harness

The V&V suite is driven by CTest (`test/CMakeLists.txt`).  Each test
runs the solver in its own case directory and then executes the
case's `check.py` (or `verify.py --gate` / family `verify.py`), which
computes the errors against the reference and returns a hard exit
code.  Plotting (`verify.py` without `--gate`) is optional and never
part of the acceptance.
