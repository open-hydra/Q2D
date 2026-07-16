"""Noh87 1D strong-shock regression check (Noh 1987 planar Cartesian).

Compares ``OUTPUT/field.tec`` (Q2D tecplot ASCII BLOCK packing on a 2D mesh --
x,y nodal + rho/u/v/p/T/g/R cell-centered) against the analytical
strong-shock solution in ``reference/exact.txt`` (5 cols: x, rho, u, p, sie),
generated offline via ``tools/gen_noh_exact.py`` (Rankine-Hugoniot post-shock
stagnation: rho=6, u=0, p=1.2 for gamma=1.4, u0=1).

Noh87 setup (mirror MOSE):
  rho0=1.0, u0=+/-1.0 (symmetric collapse), p0=1e-15
  gamma=1.4, x0=0.5, t_final=1.0, domain=[0,1]
Wave structure at t=1.0:
  pre-shock left  (x < 0.3): rho=1, u=+1, p=0
  stagnation     (0.3<x<0.7): rho=6, u=0,  p=1.2
  pre-shock right (x > 0.7): rho=1, u=-1, p=0

Tec reading (ORION + z=0 patch) and all guards go through the shared helper
library ``test/common/pycheck/q2d_check_lib.py``. Gates:

  1. NaN-safe finiteness + positivity of rho/p/T (a plain ``min() <= 0`` is
     False for an all-NaN field, so positivity MUST use ``assert_positive``);
  2. probe-point rho/u/p/sie errors (1% rho/u/sie, 1.5% p) away from the two
     shock fronts (historical MOSE gate, unchanged);
  3. whole-domain relative L1(rho) and L1(p) vs exact.txt — measured on the
     x8 solve mesh (Noh87 Q2D.sh solve links ic_x8; calibration 2026-07):
     L1(rho)=7.15e-3, L1(p)=1.00e-2 (dominated by the fixed interpolation
     smearing of the 100-point table at the two rho=1->6 jumps, not by the
     scheme); tolerances set ~1.75x the measured values;
  4. shock positions: steepest density-gradient face on each half-domain vs
     the analytic shocks x = 0.5 -/+ D*t = 0.30 / 0.70 (D = (g-1)u0/2 = 0.2),
     each within 3 cells.
"""
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_finite, assert_positive,
                           rel_l1, gate, read_q2d_field)

X_SHOCK_L, X_SHOCK_R = 0.30, 0.70  # x0 -/+ D*t, D = (gamma-1)*u0/2 = 0.2

require_file('reference/exact.txt')
ref = np.loadtxt('reference/exact.txt')
xe, rhoe, ue, pe, siee = ref[:, 0], ref[:, 1], ref[:, 2], ref[:, 3], ref[:, 4]

require_file('OUTPUT/field.tec', "run './Q2D.sh solve' first")
[x_, y_, var_, vnames] = read_q2d_field('OUTPUT/field.tec')

xc  = 0.5 * (x_[0][:-1] + x_[0][1:])[:, 0, 0]
# Q2D cell-centered order: rho(1), u, v, p, T, g, R
rho = var_[0][0][:, 0, 0]
u   = var_[0][1][:, 0, 0]
p   = var_[0][3][:, 0, 0]
T   = var_[0][4][:, 0, 0]
gam = var_[0][5][:, 0, 0]
R   = var_[0][6][:, 0, 0]
cv  = R / (gam - 1)
sie = cv * T

rhoe_i = np.interp(xc, xe, rhoe)
ue_i   = np.interp(xc, xe, ue)
pe_i   = np.interp(xc, xe, pe)
siee_i = np.interp(xc, xe, siee)

# Probes positioned away from the two shock fronts (x=0.3 and x=0.7). Sample
# pre-shock left, stagnation core, pre-shock right (avoid shock smearing).
ncells = len(rho)
iprobe = [int(ncells * f) for f in [0.10, 0.50, 0.90]]
threshold_rho  = 1.0   # percent
threshold_u    = 1.0
threshold_p    = 1.5   # mirror MOSE Noh87 verify.py (p harder near vacuum)
threshold_sie  = 1.0

# --- 1. NaN-safe sanity (finiteness is checked inside assert_positive) ------
assert_positive('rho', rho)
assert_positive('p', p)
assert_positive('T', T)
assert_finite('u', u)

# --- 2. probe-point gates (historical MOSE regression gate) -----------------
for i in iprobe:
    if rhoe_i[i] == 0.0:
        continue
    perr = abs((rho[i] - rhoe_i[i]) / rhoe_i[i]) * 100
    if perr > threshold_rho:
        fail(f'Density at cell {i}:  CFD={rho[i]:.6f}  Exact={rhoe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold_rho}%')

for i in iprobe:
    if abs(ue_i[i]) < 1e-10:
        # Skip stagnation region (u=0 makes relative error undefined; sanity
        # check above already validates physical positivity).
        continue
    perr = abs((u[i] - ue_i[i]) / ue_i[i]) * 100
    if perr > threshold_u:
        fail(f'Velocity at cell {i}:  CFD={u[i]:.6f}  Exact={ue_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold_u}%')

for i in iprobe:
    if pe_i[i] < 1e-12:
        # Pre-shock vacuum: absolute floor instead of relative.
        if p[i] > 1.0e-3:
            fail(f'Pressure at cell {i}:  CFD={p[i]:.6e} >> vacuum reference {pe_i[i]:.3e}')
        continue
    perr = abs((p[i] - pe_i[i]) / pe_i[i]) * 100
    if perr > threshold_p:
        fail(f'Pressure at cell {i}:  CFD={p[i]:.6f}  Exact={pe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold_p}%')

for i in iprobe:
    if abs(siee_i[i]) < 1e-12:
        continue
    perr = abs((sie[i] - siee_i[i]) / siee_i[i]) * 100
    if perr > threshold_sie:
        fail(f'SIE at cell {i}:  CFD={sie[i]:.6f}  Exact={siee_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold_sie}%')

# --- 3. whole-domain relative L1 gates vs exact.txt --------------------------
# Measured on the x8 (800-cell) solve mesh: L1(rho)=7.15e-03, L1(p)=1.00e-02.
gate('whole-domain rel L1(rho) vs exact', rel_l1(rho, rhoe_i), 1.25e-2)
gate('whole-domain rel L1(p) vs exact',   rel_l1(p, pe_i),     1.75e-2)

# --- 4. shock-position gates (both fronts) -----------------------------------
# Steepest density gradient per half-domain (rho jumps 1 -> 6 at each shock;
# the solution is symmetric about x0=0.5). Measured offset on x8: 0.0 cells.
h = (xc[-1] - xc[0]) / (ncells - 1)
mid = ncells // 2
for tag, sl, x_exact in (('left', slice(0, mid), X_SHOCK_L),
                         ('right', slice(mid, ncells), X_SHOCK_R)):
    seg_x, seg_r = xc[sl], rho[sl]
    isk = int(np.argmax(np.abs(np.diff(seg_r))))
    x_shock = 0.5 * (seg_x[isk] + seg_x[isk + 1])
    if abs(x_shock - x_exact) > 3.0 * h:
        fail(f'{tag} shock at x={x_shock:.4f}, expected {x_exact}+/-{3.0*h:.4f} (3 cells)')
    print(f'  [PASS] {tag} shock position: x={x_shock:.4f} (exact {x_exact}, tol 3h={3.0*h:.4f})')

print('Correct!')
