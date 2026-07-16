"""ShuOsher-LIMO3 regression check — third-order LIMO3 limiter.

Activates the LIMO3 limiter (Cada & Torrilhon 2009; rlimiter_LIMO3 in
numerics/space/Lib_Limiters.f90, dispatched by the uppercase case label 'LIMO3' in
Mod_Space::Assign_Limiter). MOSE has no LIMO3.

The test runs the Shu-Osher problem (Mach-3 shock into a sinusoidal density field;
Shu & Osher, JCP 83 (1989) 32). Downstream of the shock the entropy waves steepen into
high-frequency SMOOTH density oscillations — the textbook "shock + regular waves" setup.

WHAT THIS CHECK VERIFIES (the *interpretation* of LIMO3, not a bare value match):

  1. NEG-CONTROL (the regression that matters). LIMO3 is third-order at smooth extrema
     (φ(r<0)>0, no clipping), so it preserves much more total variation in the smooth-wave
     region than a first-order scheme. The documented failure mode is a
     SILENT fallback to first-order IORD (the historical lowercase-'limo3' case-label bug).
     If that regressed, LIMO3 output would collapse to the IORD
     baseline. We assert LIMO3's post-shock TV >= 1.5x the committed first-order/IORD
     baseline. Measured margin is 2.12x, so 1.5x is a comfortable, cross-platform-stable gate
     that fails loudly if LIMO3 is silently degraded.

  2. THIRD-ORDER >= LOWER-ORDER TVD ordering. LIMO3 resolves the smooth oscillations at
     least as well as the more dissipative TVD limiters: TV(LIMO3) >= TV(vanleer) >= TV(minmod)
     (committed baselines). This is the "less clipping" property, holds with margin.

  3. COMPARABLE-TO-MC band. MC is the least-dissipative 2nd-order TVD limiter; on Shu-Osher
     its raw TV is essentially tied with LIMO3 (the HLLC dissipation dominates the limiter
     choice for the bulk amplitude). We do NOT claim LIMO3 beats MC; we assert LIMO3 stays in
     MC's high-resolution class (within +/-10% of MC TV), proving it is not degraded.

  4. PHYSICS sanity: NaN-safe positivity (rho,p,T>0 via q2d_check_lib.assert_positive — a plain
     ``min() <= 0`` is False for an all-NaN field), physical density bounds, and the main shock
     sits at the analytic location x≈2.4 (it travels right from x=-4 over t=1.8) within a
     tolerance. Every output/baseline read goes through q2d_check_lib.require_file.

Reference baselines (reference/baseline_*.dat: x, rho at N=400, t=1.8) were generated once
with the same binary using the respective limiters (mc/vanleer/minmod) and first-order (iord),
mirror the gen_case.py setup. They are committed deterministic data, not re-run per ctest.
"""
import sys
import os
import re
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent

sys.path.insert(0, str(ROOT.parents[1] / "common" / "pycheck"))
from q2d_check_lib import fail, require_file, assert_positive, assert_finite

# Post-shock smooth-wave window (entropy waves are upstream of the shock at x≈2.4).
WIN_LO, WIN_HI = -3.0, 2.3
# Analytic main-shock location at t=1.8 (travels right from x=-4); tolerant band.
SHOCK_X, SHOCK_TOL = 2.4, 0.4


def read_field_rho(path):
    """Minimal Q2D 1D BLOCK-tec reader -> (xc, rho) cell-centered."""
    lines = Path(path).read_text().splitlines()
    m = re.search(r'I=\s*(\d+).*J=\s*(\d+).*K=\s*(\d+)', lines[1])
    ni, nj, nk = int(m.group(1)), int(m.group(2)), int(m.group(3))
    nn = ni * nj * nk
    names = re.findall(r'"([^"]+)"', lines[0])
    data = np.array([float(s) for s in (ln.strip() for ln in lines[2:]) if s])
    off = 2 * nn                       # skip x-nodes and y-nodes
    nc = (ni - 1) * (nj - 1)
    cells = {names[2 + k]: data[off + k * nc: off + (k + 1) * nc]
             for k in range(len(names) - 2)}
    xn = data[0:nn].reshape(nj, ni)[0]
    xc = 0.5 * (xn[:-1] + xn[1:])
    return xc, cells, names


def total_variation(x, r, lo, hi):
    seg = r[(x >= lo) & (x <= hi)]
    return float(np.sum(np.abs(np.diff(seg))))


# --- LIMO3 run output ---
require_file(ROOT / "OUTPUT" / "field.tec", "run './Q2D.sh solve' first")
xc, cells, names = read_field_rho(ROOT / "OUTPUT" / "field.tec")
rho = cells[names[2]]
p   = cells[names[5]]
T   = cells[names[6]]

# 4. Physics sanity (NaN-safe: assert_positive checks finiteness first) ------
assert_positive("rho", rho)
assert_positive("p", p)
assert_positive("T", T)
if not (0.7 <= rho.min() and rho.max() <= 5.0):
    fail(f"Density out of physical Shu-Osher bounds: [{rho.min():.3f}, {rho.max():.3f}]")

ishock = int(np.argmax(np.abs(np.diff(rho))))
xshock = xc[ishock]
if abs(xshock - SHOCK_X) > SHOCK_TOL:
    fail(f"Main shock at x={xshock:.3f}, expected {SHOCK_X}+/-{SHOCK_TOL}")

# Total variation of the LIMO3 solution in the smooth-wave window
tv_limo3 = total_variation(xc, rho, WIN_LO, WIN_HI)

# --- committed baselines ---
def baseline_tv(tag):
    path = ROOT / "reference" / f"baseline_{tag}.dat"
    require_file(path, "committed limiter baseline missing")
    arr = np.loadtxt(path)
    assert_finite(f"baseline_{tag}", arr)
    return total_variation(arr[:, 0], arr[:, 1], WIN_LO, WIN_HI)

tv_iord    = baseline_tv("iord")
tv_minmod  = baseline_tv("minmod")
tv_vanleer = baseline_tv("vanleer")
tv_mc      = baseline_tv("mc")

# 1. NEG-CONTROL: LIMO3 must be far above the first-order floor (not silently degraded).
if tv_limo3 < 1.5 * tv_iord:
    fail(f"NEG-CONTROL FAIL: LIMO3 TV={tv_limo3:.3f} not >= 1.5x first-order TV={tv_iord:.3f} "
         f"(ratio {tv_limo3 / tv_iord:.2f}). LIMO3 may have silently fallen back to IORD "
         f"(case-label regression in Mod_Space::Assign_Limiter).")

# 2. THIRD-ORDER >= LOWER-ORDER TVD ordering (less clipping than the diffusive limiters).
if tv_limo3 < tv_vanleer:
    fail(f"ORDERING FAIL: LIMO3 TV={tv_limo3:.3f} < vanleer TV={tv_vanleer:.3f}")
if tv_vanleer < tv_minmod:
    fail(f"ORDERING FAIL: vanleer TV={tv_vanleer:.3f} < minmod TV={tv_minmod:.3f}")

# 3. COMPARABLE-TO-MC: LIMO3 stays in the high-resolution class (within +/-10% of MC).
if not (0.90 * tv_mc <= tv_limo3 <= 1.10 * tv_mc):
    fail(f"MC-BAND FAIL: LIMO3 TV={tv_limo3:.3f} outside [0.90, 1.10]x MC TV={tv_mc:.3f} "
         f"(ratio {tv_limo3 / tv_mc:.3f})")

print(f"ShuOsher-LIMO3: shock x={xshock:.3f}  rho in [{rho.min():.3f},{rho.max():.3f}]")
print(f"  post-shock TV  LIMO3={tv_limo3:.3f}  iord={tv_iord:.3f} (ratio {tv_limo3/tv_iord:.2f}x)"
      f"  minmod={tv_minmod:.3f}  vanleer={tv_vanleer:.3f}  mc={tv_mc:.3f}")
print("Correct!")

# Explicit success exit: keeps the ctest contract (0 = all gates passed)
# aligned with the other 1D checks rather than relying on fall-through.
sys.exit(0)
