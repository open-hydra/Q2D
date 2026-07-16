"""Sod79 1D regression check.

Compares ``OUTPUT/field.tec`` (Q2D tecplot ASCII BLOCK packing, 2D no-w layout:
x,y nodal + rho/u/v/p/T/g/R cell-centered) against the analytical
Riemann solution in ``reference/exact.txt`` (5 cols: x, rho, u, p, sie).

Tec reading (ORION + z=0 patch) and all guards go through the shared helper
library ``test/common/pycheck/q2d_check_lib.py``. Gates:

  1. NaN-safe finiteness + positivity of rho/p/T (a plain ``min() <= 0`` is
     False for an all-NaN field, so positivity MUST use ``assert_positive``);
  2. probe-point rho/u/p/sie errors < 1% at smooth-region cells (historical
     MOSE/test/1D/Sod79 gate, unchanged);
  3. whole-domain relative L1(rho) and L1(p) vs exact.txt — measured on the
     x4 solve mesh (calibration 2026-07): L1(rho)=2.94e-3, L1(p)=2.08e-3;
     tolerances set ~1.7x the measured values;
  4. shock position: steepest pressure-gradient face vs the analytic shock
     location x = x0 + S*t = 0.85043 (S = 1.75216), within 3 cells.
"""
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_finite, assert_positive,
                           rel_l1, gate, read_q2d_field)

X_SHOCK_EXACT = 0.85043  # x0 + S*t = 0.5 + 1.75216*0.20 (exact Riemann solver)

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

ncells = len(rho)
iprobe = [int(ncells * f) for f in [0.20, 0.50, 0.60, 0.80, 0.90]]
threshold = 1.0  # percent

# --- 1. NaN-safe sanity (finiteness is checked inside assert_positive) ------
assert_positive('rho', rho)
assert_positive('p', p)
assert_positive('T', T)
assert_finite('u', u)

# --- 2. probe-point gates (historical 1% regression gate) -------------------
for i in iprobe:
    perr = abs((rho[i] - rhoe_i[i]) / rhoe_i[i]) * 100
    if perr > threshold:
        fail(f'Density at cell {i}:  CFD={rho[i]:.6f}  Exact={rhoe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')

for i in iprobe:
    perr = abs((u[i] - ue_i[i]) / (ue_i[i] + 1e-10)) * 100
    if perr > threshold:
        fail(f'Velocity at cell {i}:  CFD={u[i]:.6f}  Exact={ue_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')

for i in iprobe:
    perr = abs((p[i] - pe_i[i]) / pe_i[i]) * 100
    if perr > threshold:
        fail(f'Pressure at cell {i}:  CFD={p[i]:.6f}  Exact={pe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')

for i in iprobe:
    perr = abs((sie[i] - siee_i[i]) / siee_i[i]) * 100
    if perr > threshold:
        fail(f'SIE at cell {i}:  CFD={sie[i]:.6f}  Exact={siee_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')

# --- 3. whole-domain relative L1 gates vs exact.txt --------------------------
# Measured on the x4 (400-cell) solve mesh: L1(rho)=2.94e-03, L1(p)=2.08e-03.
gate('whole-domain rel L1(rho) vs exact', rel_l1(rho, rhoe_i), 5.0e-3)
gate('whole-domain rel L1(p) vs exact',   rel_l1(p, pe_i),     3.5e-3)

# --- 4. shock-position gate ---------------------------------------------------
# Steepest pressure gradient isolates the shock (the contact carries no p jump
# and the rarefaction is smooth). Measured offset on x4: 0.0021 (~0.8 cells).
h = (xc[-1] - xc[0]) / (ncells - 1)
isk = int(np.argmax(np.abs(np.diff(p))))
x_shock = 0.5 * (xc[isk] + xc[isk + 1])
if abs(x_shock - X_SHOCK_EXACT) > 3.0 * h:
    fail(f'shock at x={x_shock:.4f}, expected {X_SHOCK_EXACT}+/-{3.0*h:.4f} (3 cells)')
print(f'  [PASS] shock position: x={x_shock:.4f} (exact {X_SHOCK_EXACT}, tol 3h={3.0*h:.4f})')

print('Correct!')
