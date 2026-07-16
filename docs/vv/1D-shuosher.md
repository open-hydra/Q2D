# Shu–Osher Problem — LIMO3 Limiter

The Shu–Osher shock/entropy-wave interaction (Shu & Osher, _J.
Comput. Phys._ **83** (1989) 32, Example 8): a Mach-3 shock runs into
a sinusoidal density field.  Downstream of the shock the entropy
waves steepen into high-frequency **smooth** density oscillations —
the textbook "shock + regular waves" configuration and the canonical
demonstrator of a high-order limiter.

The case exists in the suite to regression-protect the Q2D-only
third-order limiter **LIMO3** (Čada & Torrilhon 2009;
`rlimiter_LIMO3` in `src/lib/numerics/space/Lib_Limiters.f90`,
dispatched by the uppercase case label `LIMO3` in
`Mod_Space::Assign_Limiter`).  MOSE has no LIMO3.  A second-order TVD
limiter (MC / van Leer / minmod) clips the smooth extrema back to
first order ($\phi(r<0)=0$); LIMO3 does not ($\phi(r<0)>0$, third
order at smooth extrema), so it preserves a larger post-shock
oscillation amplitude while keeping the main shock monotone.  That
no-clip-at-smooth-extrema property — not a bare value match — is what
`check.py` verifies.

**Reference:** C.-W. Shu & S. Osher, _J. Comput. Phys._ **83** (1989)
32–78; M. Čada & M. Torrilhon, _J. Comput. Phys._ **228** (2009)
4118–4145.

**Feature exercised:** 1D problem on the 2D-plane strip layout;
**LIMO3 third-order limiter** (Q2D-only) with MUSCL reconstruction; HLLC Riemann solver;
time-accurate RK3; BC `400` (extrapolation) on the streamwise ends,
BC `300` (symmetry) spanwise; committed-baseline TV gates including a
built-in first-order negative control.

## Initial conditions

Domain $x \in [-5, 5]$, $\gamma = 1.4$, final time $t = 1.8$:

| Region | $\rho$ | $u$ | $p$ |
|---|---|---|---|
| $x < -4$ | 3.857143 | 2.629369 | 10.33333 |
| $x \geq -4$ | $1 + 0.2\sin(5x)$ | 0 | 1 |

At $t = 1.8$ the main shock sits near $x \approx 2.4$ and the
boundary cells stay in constant regions, so plain extrapolation
(BC `400`) on the streamwise ends is exact; the spanwise faces are
symmetry (BC `300`), mirroring the Ein91/Tor99 1D strip layout.

## Numerical setup

| Parameter | Value |
|---|---|
| Cells | 400 (ctest mesh) |
| Final time | $t = 1.8$ (time-accurate) |
| Time scheme | RK3 |
| CFL | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | **LIMO3** |
| Riemann solver | HLLC |

## Results

<figure>
  {% include "vv/images/ShuOsher.svg" %}
</figure>

## Acceptance gates (`check.py`, ctest `ShuOsher-LIMO3`)

All gates operate on the total variation (TV) of density in the
post-shock smooth-wave window $x \in [-3.0, 2.3]$, compared against
committed baselines (`reference/baseline_{iord,minmod,vanleer,mc}.dat`:
$x$, $\rho$ at $N=400$, $t=1.8$, generated once with the same binary
and the respective limiter — deterministic committed data, not re-run
per ctest):

1. **Negative control** (the regression that matters).  The
   documented failure mode is a *silent* fallback to first-order
   (the historical lowercase-`limo3` case-label bug,
   `Mod_Space.f90:69-72`).  The check asserts LIMO3's post-shock
   TV ≥ **1.5×** the committed first-order/IORD baseline; the
   measured margin is **2.12×**, so the 1.5× gate fails loudly if
   LIMO3 is silently degraded.
2. **Third-order ≥ lower-order TVD ordering**:
   TV(LIMO3) ≥ TV(vanleer) ≥ TV(minmod) — the "less clipping"
   property, which holds with margin.
3. **Comparable-to-MC band**: MC is the least-dissipative
   second-order TVD limiter; on Shu–Osher its raw TV is essentially
   tied with LIMO3.  The check does not claim LIMO3 beats MC — it
   asserts LIMO3 stays within **±10 %** of the MC TV (the
   high-resolution class).
4. **Physics sanity**: NaN-safe positivity of $\rho$, $p$, $T$;
   density within the physical Shu–Osher bounds $[0.7, 5.0]$; and the
   main shock located at the analytic position $x \approx 2.4$
   (travelled from $x=-4$ over $t=1.8$) within a ±0.4 band.

Run it standalone:

```bash
cd test/1D/ShuOsher-LIMO3
./Q2D.sh solve
python3 -B check.py
```
