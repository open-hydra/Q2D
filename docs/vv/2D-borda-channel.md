# Borda Channel — Choked Sudden-Expansion Injector BC

A straight plane channel fed by the Q2D-only **Borda
sudden-expansion injector** inflow boundary condition
(`BC_Nozzle_Borda_Flow`; Fievisohn 2016 thesis, Sec. 3.3,
eqs. 3.26–3.31).  The BC models a Borda (re-entrant) injector: a
choked throat (station 1) followed by a sudden expansion into the
channel (station 3), with the associated Borda–Carnot
stagnation-pressure loss.  MOSE has no nozzle/Borda BC — this is a
Q2D-specific feature regression.

The BC is selected with code **420** and the sentinel
`BC_Mach == -20`; the area ratio $A_1/A_3$ is carried in the `pamb`
field of `INPUT/bc.txt`.

**Reference:** R. T. Fievisohn, *Development and Application of
Theoretical Models for Rotating Detonation Engine Flowfields*, PhD
dissertation, University of Maryland, College Park, 2016 — Sec. 3.3,
eqs. 3.26–3.31.

**Feature exercised:** Borda injector inflow **BC `420`**
(`BC_Nozzle_Borda_Flow`, peculiarity #4, sentinel `BC_Mach == -20`,
$A_1/A_3$ in the `pamb` field of `bc.txt`) in the 2D-plane regime;
BC `400` extrapolation exit; MUSCL/vanleer + HLLC with IRS; the
`-negctl` sibling exercises the **negative-control** path
(`BORDA_NEG_CONTROL=1`) and the CTest `RESOURCE_LOCK` serialization
of tests sharing one case directory.

## Problem setup

| Parameter | Value |
|---|---|
| Geometry | straight plane channel, 80 × 8 cells |
| Gas | constant-$\gamma$ air (`common/Air`, $\gamma = 1.4$, $R = 287$) |
| Plenum | $P_0 = 1.0\times10^{6}$ Pa, $T_0 = 300$ K |
| Area ratio $A_1/A_3$ | 0.2 |
| Interior back pressure | $P_3 = 0.70\,P_0$ |
| BCs | face 1 Borda inflow `420`; face 2 extrapolation `400`; faces 3–4 symmetry |

With $P_3 = 0.70\,P_0$ the case sits in **regime (A)**: choked throat
($M_1 = 1$), subsonic station 3 — the regime where the Fievisohn jump
algebra is valid.  The injected face state is checked against the
paper's own equations evaluated at this fully-specified gas (the
paper tabulates no numeric case; no thresholds are invented).

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3, local time stepping (5000 iterations) |
| CFL | 0.5 (rise threshold 50) |
| IRS β | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | HLLC |

## Results

<figure>
  {% include "vv/images/Borda-state.svg" %}
</figure>

## Acceptance gates (`check.py`, ctest `borda-channel`)

- positivity: min($\rho$, $p$, $T$) > 0 (hard, NaN-safe);
- convergence sentinel (hard): final density residual < **1e-6**
  (measured 1.0e-10 on the committed run);
- **guard — Borda branch actually entered**: the injected mass flux
  must exceed 10 % of the analytic choked value.  A silent wall
  fallback would give $u \approx 0$ at the inlet and pass a naive
  check trivially; this guard is the false-green prevention;
- **(A) choked mass flux**: inlet $\rho u$ vs the analytic
  $\dot m_{A_3} = \rho_1 v_1 A_1/A_3$ (eq. 3.29 at $M_1 = 1$),
  rel err < **3e-2**;
- **(B) isentropic choke state**: the stagnation temperature
  recovered from the injected $(\rho, u, p, T)$ via
  isentropic + energy relations must match the plenum $T_0$,
  rel err < **3e-2**;
- **(C) stagnation-pressure-loss interpretation** (the Borda–Carnot
  hallmark): the injected stagnation pressure $p_{0,\text{inj}}$ must
  be strictly **below** the plenum $P_0$ (a real loss, entropy
  non-decreasing).

## Negative control (ctest `borda-channel-negctl`)

A second ctest re-runs `check.py` with `BORDA_NEG_CONTROL=1`: the
check corrupts its own expectation (the choked $\dot m_{A_3}$ is
multiplied by 1.5) and must **detect the mismatch** — it exits 0 iff
the corrupted gate (A) FAILS.  This institutionalizes the false-green
tripwire in CI.  The two tests share the case directory (and thus
`OUTPUT/`), so CTest serializes them via a `RESOURCE_LOCK`.

## How to run

```bash
cd test/2D/euler/borda-channel
./Q2D.sh solve
python3 -B check.py                     # validation
BORDA_NEG_CONTROL=1 python3 -B check.py # negative control
```
