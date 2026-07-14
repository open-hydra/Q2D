# Verification & Validation

This section documents the Verification & Validation (V&V) test
suite for Q2D.  The purpose is to demonstrate that the solver
produces correct results across a range of physical regimes, by
comparison against analytical or well-established reference
solutions, and to specify the regression tolerance budget for
bit-per-bit reproducibility against Q2D legacy outputs.

## Test suite

27 cases are wired into CTest and pass on the default mesh. The suite
includes inviscid Riemann/canonical 2D (13 cases) and viscous regimes
(flat-plate-laminar, flat-plate-turbulent, shock-wave-boundary-layer-interaction).

| Test | Dim | Mach | Physics | Verification | Solver(s) | Ref |
|---|---|---|---|---|---|---|
| [Sod Shock Tube](1D.md#sod-shock-tube) | 1D | subsonic–supersonic | Shock, contact, rarefaction | Analytical (Toro Riemann) | SLAU | Sod (1978) |
| [Einfeldt Double Rarefaction](1D.md#einfeldt-double-rarefaction) | 1D | subsonic | Near-vacuum, two rarefactions | Analytical (Toro Riemann) | HLLE | Einfeldt et al. (1991) |
| [Toro Test 3](1D.md#toro-test-case-3) | 1D | subsonic–supersonic | Compound wave (shock + contact + shock) | Analytical (Toro Riemann) | HLLE++ | Toro (1999) |
| [Noh Implosion](1D.md#noh-implosion-problem) | 1D | supersonic | Strong shock, density jump | Analytical (Noh) | HLLE | Noh (1987) |
| [Chang Nozzle](2D-plane-extras.md#chang-nozzle) | 2D axi | transonic | Single-species axisymmetric expansion (peculiarity #2) | MOSE (3D-wedge K=0) | SLAU | Chang (1996) |
| [Prandtl–Meyer](2D-plane-extras.md#prandtlmeyer-expansion) | 2D planar | 2.0 | Smooth supersonic expansion fan | Analytical (PM theory) + MOSE | SLAU | Prandtl–Meyer (1908) |
| [Hypersonic Cylinder](2D-hypersonic-cylinder.md) | 2D planar | 8.1 | Bow shock, carbuncle | Analytical normal shock + MOSE | HLLE (carbuncle-stable) | Anderson |
| [Ramp Channel](2D-plane-extras.md#ramp-channel) | 2D planar | 3.0 | Oblique shock reflection | Analytical θ–β–M + MOSE | SLAU | classic |
| [Oblique Shock](2D-oblique-shock.md) | 2D planar | 4.0 | Oblique shock, post-shock instabilities | Analytical θ–β–M + MOSE | HLLC, HLLC+, HLLE, SLAU | oblique shock theory |
| [Supersonic Forward Step (WC)](2D-woodward-colella.md) | 2D planar | 3.0 | Shock + expansion + interactions | Reference + MOSE | HLLE | Woodward & Colella (1984) |
| [Rocket Nozzle — frozen](2D-nozzle.md) | 2D axi | transonic–supersonic | 13-species frozen multi-fluid (peculiarity #2 axisymmetric) | CEA 1D + MOSE | SLAU | CEA / NASA |
| [Rocket Nozzle — equilibrium](2D-nozzle.md) | 2D axi | transonic–supersonic | 13-species instantaneous equilibrium | CEA 1D + MOSE | SLAU | CEA / NASA |
| [Rocket Nozzle — finite-rate](2D-nozzle.md) | 2D axi | transonic–supersonic | 19-species CORIA reaction mechanism | CEA 1D + MOSE | SLAU | CEA / NASA + CORIA |

**Viscous and RANS regimes** (wired and passing):

| Test | Dim | Mach | Physics | Verification | Solver(s) | Ref |
|---|---|---|---|---|---|---|
| [Laminar Flat Plate](2D-flat-plate-laminar.md) | 2D planar | 0.2 | Laminar boundary layer, viscous effects | Blasius similarity | HLLC | Blasius (1908) |
| [Turbulent Flat Plate](2D-flat-plate-turbulent.md) | 2D planar | 0.2 | Turbulent boundary layer, RANS closures | NASA TMR (CFL3D, FUN3D) | HLLC | NASA |
| [Shock Wave-Boundary Layer Interaction](2D-swbli.md) | 2D | 5.0 | Shock-boundary-layer interaction, separation/reattachment | Schulein/SU2/Wind-US comparison | HLLC (SA) | Schulein + code-to-code |

| Document | Scope |
|---|---|
| [Regression Tolerance](regression_tolerance.md) | Numerical tolerance for bit-per-bit / floating-point regression of Q2D solver outputs |

## Running the tests

Each test case lives under `test/` and contains a `verify.py`
script.  The script reads the solver output, compares it against the
exact (or reference) solution, and writes output figures.

```bash
cd test/1D/<TestName>
./Q2D.sh solve
python verify.py          # check errors and export figure
python verify.py --plot   # as above, also display figure interactively
```

---

## Regression Strategy

Q2D aims at **bit-per-bit reproducibility** of Q2D legacy outputs
in `RELEASE -O0` (relative error $< 10^{-12}$) and
**tolerance-bounded reproducibility** in `RELEASE -O3` (relative
error $< 10^{-9}$, accommodating floating-point reordering).  The
tolerance is documented per-test in
[`regression_tolerance.md`](regression_tolerance.md).

The benchmark cases above share the test grid layout with MOSE where
applicable, with adjustments for Q2D's three geometric regimes (2D
planar, 2D axisymmetric, q2D area-variation).

!!! success "V&V execution status (v1.2.0)"

    The full V&V suite is **closed and passing**. All 27 test cases (4×1D Riemann,
    6×2D planar, 3×2D axisymmetric chemistry, 3×viscous/RANS, 11 additional inviscid)
    pass CTest on the default mesh: `ctest --test-dir build --output-on-failure`
    produces 27/27 PASS in ~18 min (rocket-nozzle-finite-rate is the longest at 957 s).

    For each case, the per-test `verify.py` script validates the Q2D solution
    against either an analytical reference (1D Riemann + Noh), experimental data
    (Blasius, NASA TMR, Schulein), or MOSE comparison (3D-wedge K=0 layer).
    Q2D vs MOSE agrees at `Linf_rel ≤ 1e-12` for all 2D planar cases; the 2D-axi cases
    (chang-nozzle, rocket-nozzle ×3) exhibit a 4–6 % `eps_rel` discrepancy along the
    axis row (documented in closure archives).

    1D mesh-refinement convergence (×1/×2/×4/×8, 16 runs) is summarised
    in [convergence-table-1D.md](convergence-table-1D.md): all four Riemann + Noh
    cases meet the 1–1.5 % pointwise relative error acceptance criterion on probe
    points at ×8 refinement.

---

## Verification Harness

The V&V suite is driven by CTest. Each test case includes a `verify.py` script that
computes point-wise errors, compares against reference solutions, and generates
output plots. Regression tolerance parameters are documented in
[`regression_tolerance.md`](regression_tolerance.md).
