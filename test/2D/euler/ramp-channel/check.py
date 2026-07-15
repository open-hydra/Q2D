"""ramp-channel 2D steady-state verification check.

M=3 plane Euler flow over a 30-deg ramp inside a channel (oblique shock from
the ramp tip at x=0.1; the ramp shoulder at x~0.187 sheds an expansion fan
and the shock reflects off the upper boundary further downstream).

Gates (all hard):
  (1) Oblique-shock theta-beta-Mach theory (M1=3, gamma=1.4, ramp delta
      measured from wall nodes i=20..50 j=0): post/pre shock ratios averaged
      over a pre-shock probe (cells i=5:15, j=70:90) and a post-shock plateau
      probe (cells i=36:48, j=0:12 -- inside the uniform wedge region,
      downstream of the shock foot and upstream of the shoulder expansion).
      Tolerances calibrated at ~2x the measured errors (2026-07-15, gfortran
      TESTING build: p2/p1 -0.002%, T2/T1 +0.95%, rho2/rho1 -0.94%,
      M2 -0.94%): |p2/p1 err| < 1%, |T2/T1|, |rho2/rho1|, |M2 err| < 2%.
  (2) Shock-detector beta field (SD_Tramel, beta=0 in shock cells, beta=1 in
      smooth flow): beta differentiates, >=1% of cells detected, detected
      cells localized downstream of the ramp tip (>=90%; observed 99.3%),
      smooth pre-shock block beta mean > 0.99, and a negative control (an
      "expected shock zone" placed upstream of the tip must be ~empty, <=2%).
  (3) Convergence + IRS stability (irs-beta=0.5 active with SD): residual
      history finite and final density residual < 1e-3 (measured final
      2.49e-4 after 10000 iters, peak 7.98e-3 at shock formation ~iter 99;
      threshold ~4x the measured final -- a destabilized IRS diverges
      instead, cf. the hypersonic-cylinder M=8.1 + IRS experiment).
  (4) NaN-safe sanity: rho, p, T, g, R > 0; u, v finite.
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (assert_converged, assert_finite, assert_positive,
                           fail, gate, read_q2d_field, require_file)

require_file("OUTPUT/field.tec", "solver did not write the field")
require_file("OUTPUT/diagnostic.tec", "solver did not write the SD diagnostic")


# ---------------------------------------------------------------------------
# Oblique-shock theta-beta-Mach theory
# ---------------------------------------------------------------------------
def tbm(M, beta, gamma):
    c1 = M * M * math.sin(beta) ** 2 - 1.0
    c2 = M * M * (gamma + math.cos(2.0 * beta)) + 2.0
    return 2.0 * (1.0 / math.tan(beta)) * (c1 / c2)


def sckobl(mach1, delta, gamma, tol):
    """Weak-solution oblique shock for wedge angle delta: returns
    (M2, beta, p2/p1, T2/T1, rho2/rho1)."""
    theta1 = 1.2 * delta
    delt1 = tbm(mach1, theta1, gamma)
    theta2 = 1.02 * theta1
    while True:
        delt2 = tbm(mach1, theta2, gamma)
        if abs((delt2 - math.tan(delta)) / math.tan(delta)) <= tol:
            break
        theta = theta2 - (math.tan(delta) - delt2) * (theta2 - theta1) / (delt1 - delt2)
        theta1 = theta2
        delt1 = delt2
        theta2 = theta
    theta = theta2
    rmn1 = mach1 * math.sin(theta)
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0
    r2r1 = gp1 * rmn1 ** 2 / (gm1 * rmn1 ** 2 + 2.0)
    p2p1 = 1.0 + 2.0 * gamma * (rmn1 ** 2 - 1.0) / gp1
    rmn2 = math.sqrt((rmn1 ** 2 + 2.0 / gm1) / (2.0 * gamma * rmn1 ** 2 / gm1 - 1.0))
    t2t1 = p2p1 / r2r1
    mach2 = rmn2 / math.sin(theta - delta)
    return mach2, theta, p2p1, t2t1, r2r1


M1 = 3.0
GAMMA = 1.4

x_, y_, var_, names = read_q2d_field("OUTPUT/field.tec")
cc = [n.strip() for n in names[3:]]
missing = [n for n in ("rho(1)", "u", "v", "p", "T", "g", "R") if n not in cc]
if missing:
    fail("missing variables in OUTPUT/field.tec: %s (found %s)" % (missing, cc))


def field(name):
    return np.asarray(var_[0][cc.index(name)], float)[:, :, 0]


xn = x_[0][:, :, 0]
yn = y_[0][:, :, 0]
rho, u, v = field("rho(1)"), field("u"), field("v")
p, T, g, R = field("p"), field("T"), field("g"), field("R")

# (4) NaN-safe sanity
for name, arr in (("rho", rho), ("p", p), ("T", T), ("g", g), ("R", R)):
    assert_positive(name, arr)
for name, arr in (("u", u), ("v", v)):
    assert_finite(name, arr)
print("(4) sanity OK: rho, p, T, g, R > 0; u, v finite")

# (1) theta-beta-Mach hard gates
delta_deg = float(np.degrees(np.arctan2(yn[50, 0] - yn[20, 0], xn[50, 0] - xn[20, 0])))
delta = math.radians(delta_deg)
mach2_ref, beta_th, p2p1_ref, t2t1_ref, r2r1_ref = sckobl(M1, delta, GAMMA, 1e-6)
print("(1) oblique-shock theory (M1=%g gamma=%g):" % (M1, GAMMA))
print("    ramp delta=%.4f deg  beta=%.4f deg" % (delta_deg, math.degrees(beta_th)))
print("    theory: p2/p1=%.6f T2/T1=%.6f r2/r1=%.6f M2=%.6f"
      % (p2p1_ref, t2t1_ref, r2r1_ref, mach2_ref))

# Pre-shock probe: far upstream, away from the wall.
i_pre, j_pre = slice(5, 15), slice(70, 90)
# Post-shock plateau probe: on the ramp, behind the shock foot, upstream of
# the shoulder (i=50) expansion fan.
i_post, j_post = slice(36, 48), slice(0, 12)


def probe(arr, si, sj):
    return float(np.mean(arr[si, sj]))


p1_, rho1, T1_ = probe(p, i_pre, j_pre), probe(rho, i_pre, j_pre), probe(T, i_pre, j_pre)
p2, rho2, T2 = probe(p, i_post, j_post), probe(rho, i_post, j_post), probe(T, i_post, j_post)
u2, v2 = probe(u, i_post, j_post), probe(v, i_post, j_post)
g2, R2 = probe(g, i_post, j_post), probe(R, i_post, j_post)
m2 = math.sqrt(u2 ** 2 + v2 ** 2) / math.sqrt(g2 * R2 * T2)

print("    probe pre  (i=5:15, j=70:90): p1=%.6f rho1=%.6f T1=%.6f" % (p1_, rho1, T1_))
print("    probe post (i=36:48, j=0:12): p2=%.6f rho2=%.6f T2=%.6f M2=%.4f"
      % (p2, rho2, T2, m2))

gate("|p2/p1 - theory| / theory", abs(p2 / p1_ - p2p1_ref) / p2p1_ref, 0.01)
gate("|T2/T1 - theory| / theory", abs(T2 / T1_ - t2t1_ref) / t2t1_ref, 0.02)
gate("|rho2/rho1 - theory| / theory", abs(rho2 / rho1 - r2r1_ref) / r2r1_ref, 0.02)
gate("|M2 - theory| / theory", abs(m2 - mach2_ref) / mach2_ref, 0.02)

# ---------------------------------------------------------------------------
# (2) Shock-detector beta gates -- piano §5 Fase 8a criterion C1.
# Convention: SD_Tramel returns beta=0 in shock cells, beta=1 in smooth flow
# (per src/lib/numerics/fluxes/Lib_Shock_Detector.f90).
# ---------------------------------------------------------------------------
print("(2) shock-detector beta verification:")
xd, yd, vard, named = read_q2d_field("OUTPUT/diagnostic.tec")
ccd = [n.strip().lower() for n in named[3:]]
if "beta" not in ccd:
    fail("beta column not in OUTPUT/diagnostic.tec VARIABLES (%s)" % ccd)
beta2d = np.asarray(vard[0][ccd.index("beta")], float)[:, :, 0]
assert_finite("beta", beta2d)
nc = beta2d.size
n_shock = int((beta2d < 0.99).sum())
print("    beta_min=%.6f beta_max=%.6f" % (beta2d.min(), beta2d.max()))
print("    shock-detected cells (beta<0.99): %d/%d (%.2f%%)"
      % (n_shock, nc, 100.0 * n_shock / nc))
if beta2d.max() == beta2d.min():
    fail("beta field is uniform (SD not differentiating)")
if n_shock < nc // 100:
    fail("shock cells < 1%% (%.4f%%) -- SD not detecting shocks" % (100.0 * n_shock / nc))
print("    PASS: SD beta differentiates (shock + smooth regions both present)")

# (2b) SD spatial localization + negative control. The sensor must fire WHERE
# the oblique shock is (downstream of the ramp tip) and stay quiet in the
# smooth pre-shock region -- not everywhere, not nowhere. Thresholds read from
# the committed steady run: >=90% of beta<0.5 cells downstream of the tip
# (observed 99.3%); far-upstream smooth block beta mean > 0.99 (observed
# 1.000). Negative control: an "expected shock zone" placed UPSTREAM of the
# tip must contain (almost) no detected cells.
xnd = xd[0][:, :, 0]
xc = 0.25 * (xnd[:-1, :-1] + xnd[1:, :-1] + xnd[:-1, 1:] + xnd[1:, 1:])
x_tip = float(xnd[20, 0])
shock_mask = beta2d < 0.5
n_sh = int(shock_mask.sum())
print("(2b) SD spatial localization + negative control (ramp tip x ~ %.4f):" % x_tip)
if n_sh == 0:
    fail("no beta<0.5 cells to localize")
frac_downstream = float((xc[shock_mask] > x_tip).mean())
print("    beta<0.5 cells downstream of tip: %.1f%% (%d cells, x in [%.3f,%.3f])"
      % (100 * frac_downstream, n_sh, xc[shock_mask].min(), xc[shock_mask].max()))
if frac_downstream < 0.90:
    fail("only %.1f%% of detected cells are downstream of the ramp tip (< 90%%) "
         "-- sensor not localized to the shock" % (100 * frac_downstream))
smooth_mean = float(beta2d[2:12, 40:60].mean())
print("    smooth pre-shock block (i=2:12 j=40:60) beta mean = %.4f" % smooth_mean)
if smooth_mean < 0.99:
    fail("pre-shock smooth region beta mean %.4f < 0.99 -- sensor firing in "
         "smooth flow (false positives)" % smooth_mean)
upstream_cells = xc < x_tip
n_sh_upstream = int((shock_mask & upstream_cells).sum())
n_up = int(upstream_cells.sum())
frac_up = n_sh_upstream / max(n_up, 1)
print("    [neg-control] beta<0.5 cells UPSTREAM of tip: %d/%d (%.2f%%) -- must be ~0"
      % (n_sh_upstream, n_up, 100 * frac_up))
if frac_up > 0.02:
    fail("negative control -- %.2f%% of upstream cells flagged as shock (> 2%%); "
         "sensor mislocalized" % (100 * frac_up))
print("    PASS: SD localized to the shock (downstream of tip), quiet upstream")

# ---------------------------------------------------------------------------
# (3) Convergence + IRS stability (irs-beta=0.5 active with SD): the density
# residual (column 2) must be finite throughout and settle below 1e-3
# (measured final 2.49e-4, shock-formation peak 7.98e-3 at ~iter 99).
# ---------------------------------------------------------------------------
print("(3) convergence + IRS stability (irs-beta=0.5 active with SD):")
assert_converged(1e-3, column=2)

print("Correct!")
