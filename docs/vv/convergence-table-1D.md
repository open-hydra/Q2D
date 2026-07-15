# 1D mesh-refinement convergence

Each 1D Riemann/Noh case is run on a ×1/×2/×4/×8 refinement ladder
(`./Q2D.sh test` produces `OUTPUT/field_x{1,2,4,8}.tec`) and verified
against the exact solution. Two complementary metrics:

- **L1 relative error** (integral over the whole domain): convergence-friendly metric for shock-capturing schemes. Expected `p_obs` ≈ 1.0 for first-order shocks, ≈ 1.5-2.0 in smooth regions. Insensitive to single-cell shock-position mismatch between numerical and exact.
- **Linf relative error on probe points** (cells at fixed fractions 0.20/0.50/0.60/0.80/0.90 of the domain): the 1% acceptance threshold evaluated away from discontinuities.

**Feature exercised:** observed-order-of-convergence verification of
the 1D Euler core on the ×1/×2/×4/×8 refinement ladder
(`./Q2D.sh test`): MUSCL reconstruction with the per-case limiter and
Riemann solver (Sod79: MC+SLAU; Ein91: vanleer+HLLE; Tor99:
vanleer+HLLE++; Noh87: first-order+HLLE), gated against the analytic
exact solutions in `test/common/pycheck/q2d_exact1d.py` via
`verify.py --gate`.

<figure>
  {% include "vv/images/convergence-1D.svg" %}
</figure>

## Convergence gates (ctest)

The ladder is wired into CTest as four dedicated tests —
`Sod79-convergence`, `Ein91-convergence`, `Tor99-convergence`,
`Noh87-convergence` (labels `1D;convergence;fast`, 900 s timeout).
Each runs

```bash
./Q2D.sh test && python3 -B verify.py --gate
```

In `--gate` mode `verify.py` is a numpy-only hard gate (exit 0/1, no
matplotlib/ExactPack): it evaluates the exact solution **analytically
at every cell center** of every ladder level via
`test/common/pycheck/q2d_exact1d.py` — an exact ideal-gas Riemann
solver (E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid
Dynamics*, 3rd ed., Springer, 2009, ch. 4: Newton iteration for $p^*$,
positivity-guarded
for the Ein91 near-vacuum star state $p^* \sim 1.9\times10^{-3}$, then
vectorized similarity sampling) and the planar Noh (1987) closed form
(shocks at $x_0 \pm Dt$, $D = (\gamma-1)u_0/2$, stagnation
$\rho = \rho_0(\gamma+1)/(\gamma-1)$).  The committed 100-point
`reference/exact.txt` tables are *not* used here: linear interpolation
of a table smears every discontinuity over one table spacing — an
O(1e-3) L1 error floor that would dominate the ×8 level and destroy
the observed order.

Per-case asserts (all hard):

- NaN-safe positivity of $\rho$/$p$ on every level;
- whole-domain relative $L_1(\rho)$ errors monotonically
  non-increasing across the ladder;
- observed order $p_\text{obs}$ (least-squares slope of
  $\log L_1(\rho)$ vs $\log h$) above a per-case threshold set at
  ~0.85× the measured value (calibration 2026-07):

| ctest | Gate $p_\text{obs} \geq$ | Measured $p_\text{obs}$ | Measured L1(ρ) ×1/×2/×4/×8 |
|---|---:|---:|---|
| `Sod79-convergence` | 0.70 | 0.85 | 8.97e-3 / 4.72e-3 / 2.78e-3 / 1.50e-3 |
| `Ein91-convergence` | 0.78 | 0.92 | 2.20e-2 / 1.16e-2 / 6.94e-3 / 3.10e-3 |
| `Tor99-convergence` | 0.75 | 0.89 | 1.37e-2 / 7.37e-3 / 3.86e-3 / 2.18e-3 |
| `Noh87-convergence` | 0.85 | 1.01 | 1.63e-2 / 7.98e-3 / 3.81e-3 / 2.03e-3 |

These are shock-dominated problems: the formal order collapses to ~1
at best (for Noh87, with infinite-strength shocks, L1 order ~1 is the
theoretical ceiling), which is what the gates encode.

Without `--gate`, `verify.py` is the original offline study —
ExactPack exact solution, 2×2 matplotlib SVG figure, ×8 probe-point
gates, L2/Linf norm tables. The tables below are the measured
refinement-ladder record; the L1(ρ) columns are the calibration
values quoted in each case's `verify.py` docstring.

## Sod79  (riemann)

### L1 relative error

| Level | Cells | L1_rel rho | L1_rel u | L1_rel p | p_obs rho | p_obs u | p_obs p |
|------:|------:|-----------:|---------:|---------:|----------:|--------:|--------:|
| x1 | 100 | 8.965e-03 | 2.202e-02 | 7.042e-03 | -- | -- | -- |
| x2 | 200 | 4.719e-03 | 1.162e-02 | 3.674e-03 | 0.93 | 0.92 | 0.94 |
| x4 | 400 | 2.784e-03 | 6.942e-03 | 2.123e-03 | 0.76 | 0.74 | 0.79 |
| x8 | 800 | 1.497e-03 | 3.373e-03 | 9.873e-04 | 0.90 | 1.04 | 1.10 |

### Linf probe-point relative error (1% acceptance)

| Level | Cells | Linf rho [%] | Linf u [%] | Linf p [%] |
|------:|------:|-------------:|-----------:|-----------:|
| x1 | 100 | 0.179 | 0.474 | 0.015 |
| x2 | 200 | 0.040 | 0.125 | 0.012 |
| x4 | 400 | 0.009 | 0.028 | 0.004 |
| x8 | 800 | 0.018 | 0.046 | 0.018 |

## Ein91  (riemann)

### L1 relative error

| Level | Cells | L1_rel rho | L1_rel u | L1_rel p | p_obs rho | p_obs u | p_obs p |
|------:|------:|-----------:|---------:|---------:|----------:|--------:|--------:|
| x1 | 100 | 2.203e-02 | 2.035e-02 | 3.265e-02 | -- | -- | -- |
| x2 | 200 | 1.159e-02 | 1.213e-02 | 1.592e-02 | 0.93 | 0.75 | 1.04 |
| x4 | 400 | 6.936e-03 | 1.120e-02 | 1.043e-02 | 0.74 | 0.12 | 0.61 |
| x8 | 800 | 3.097e-03 | 8.622e-03 | 4.089e-03 | 1.16 | 0.38 | 1.35 |

### Linf probe-point relative error (1% acceptance)

| Level | Cells | Linf rho [%] | Linf u [%] | Linf p [%] |
|------:|------:|-------------:|-----------:|-----------:|
| x1 | 100 | 5.043 | 1.999 | 6.846 |
| x2 | 200 | 3.046 | 1.219 | 4.118 |
| x4 | 400 | 2.100 | 0.846 | 2.833 |
| x8 | 800 | 0.930 | 0.374 | 1.256 |

## Tor99  (riemann)

### L1 relative error

| Level | Cells | L1_rel rho | L1_rel u | L1_rel p | p_obs rho | p_obs u | p_obs p |
|------:|------:|-----------:|---------:|---------:|----------:|--------:|--------:|
| x1 | 100 | 1.372e-02 | 2.148e-02 | 1.024e-02 | -- | -- | -- |
| x2 | 200 | 7.372e-03 | 9.703e-03 | 4.713e-03 | 0.90 | 1.15 | 1.12 |
| x4 | 400 | 3.864e-03 | 4.309e-03 | 2.146e-03 | 0.93 | 1.17 | 1.14 |
| x8 | 800 | 2.175e-03 | 1.722e-03 | 1.015e-03 | 0.83 | 1.32 | 1.08 |

### Linf probe-point relative error (1% acceptance)

| Level | Cells | Linf rho [%] | Linf u [%] | Linf p [%] |
|------:|------:|-------------:|-----------:|-----------:|
| x1 | 100 | 1.097 | 0.958 | 1.532 |
| x2 | 200 | 0.469 | 0.088 | 0.141 |
| x4 | 400 | 0.057 | 0.003 | 0.002 |
| x8 | 800 | 0.002 | 0.002 | 0.001 |

## Noh87  (noh)

### L1 relative error

| Level | Cells | L1_rel rho | L1_rel u | L1_rel p | p_obs rho | p_obs u | p_obs p |
|------:|------:|-----------:|---------:|---------:|----------:|--------:|--------:|
| x1 | 100 | 1.633e-02 | 2.061e-02 | 2.048e-02 | -- | -- | -- |
| x2 | 200 | 7.977e-03 | 1.012e-02 | 1.001e-02 | 1.03 | 1.03 | 1.03 |
| x4 | 400 | 3.815e-03 | 4.872e-03 | 4.809e-03 | 1.06 | 1.05 | 1.06 |
| x8 | 800 | 2.029e-03 | 2.561e-03 | 2.542e-03 | 0.91 | 0.93 | 0.92 |

### Linf probe-point relative error (1% acceptance)

| Level | Cells | Linf rho [%] | Linf u [%] | Linf p [%] |
|------:|------:|-------------:|-----------:|-----------:|
| x1 | 100 | 0.986 | 0.001 | 0.002 |
| x2 | 200 | 0.666 | 0.000 | 0.001 |
| x4 | 400 | 0.459 | 0.000 | 0.000 |
| x8 | 800 | 0.321 | 0.000 | 0.000 |

