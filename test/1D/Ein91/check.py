"""Einfeldt 1991 (Ein91) strong-rarefaction 1D regression check.

Compares ``OUTPUT/field.tec`` (Q2D_new tecplot ASCII BLOCK packing, peculiarita'
#1: 2D mesh) against the analytical Riemann solution in ``reference/exact.txt``
(5 cols: x, rho, u, p, sie). Mirror of Sod79/check.py with tolerance widened (5%)
to account for the harder Ein91 problem (low-density rarefaction near contact).

Tec reading (ORION + z=0 patch) and all guards go through the shared helper
library ``test/common/pycheck/q2d_check_lib.py``. Gates:

  1. NaN-safe finiteness + positivity of rho/p/T (a plain ``min() <= 0`` is
     False for an all-NaN field, so positivity MUST use ``assert_positive``);
  2. probe-point rho/u/p errors < 5% on the smooth rarefaction tails
     (historical gate, unchanged);
  3. whole-domain relative L1(rho) and L1(p) vs exact.txt — measured on the
     x4 solve mesh (calibration 2026-07): L1(rho)=6.85e-3, L1(p)=1.03e-2;
     tolerances set ~1.75x the measured values.

No shock-position gate: the 1-2-0-3 solution is two rarefactions with no
discontinuity, so a steepest-gradient locator is not meaningful here.
"""
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_finite, assert_positive,
                           rel_l1, gate, read_q2d_field)

require_file('reference/exact.txt')
ref = np.loadtxt('reference/exact.txt')
xe, rhoe, ue, pe, siee = ref[:, 0], ref[:, 1], ref[:, 2], ref[:, 3], ref[:, 4]

require_file('OUTPUT/field.tec', "run './Q2D.sh solve' first")
[x_, y_, var_, vnames] = read_q2d_field('OUTPUT/field.tec')

xc  = 0.5 * (x_[0][:-1] + x_[0][1:])[:, 0, 0]
# Q2D_new cell-centered order: rho(1), u, v, p, T, g, R
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

ncells = len(rho)
# Probe on the smooth rarefaction tails away from the central trough.
iprobe = [int(ncells * f) for f in [0.05, 0.15, 0.85, 0.95]]
threshold = 5.0  # percent — looser than Sod79 due to strong rarefaction

# --- 1. NaN-safe sanity (finiteness is checked inside assert_positive) ------
assert_positive('rho', rho)
assert_positive('p', p)
assert_positive('T', T)
assert_finite('u', u)

# --- 2. probe-point gates (historical 5% regression gate) -------------------
for i in iprobe:
    perr = abs((rho[i] - rhoe_i[i]) / rhoe_i[i]) * 100
    if perr > threshold:
        fail(f'Density at cell {i}:  CFD={rho[i]:.6f}  Exact={rhoe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')

for i in iprobe:
    # Ein91 has |u| ≈ 2 across the field — pure relative check is safe.
    perr = abs((u[i] - ue_i[i]) / (abs(ue_i[i]) + 1e-10)) * 100
    if perr > threshold:
        fail(f'Velocity at cell {i}:  CFD={u[i]:.6f}  Exact={ue_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')

for i in iprobe:
    perr = abs((p[i] - pe_i[i]) / pe_i[i]) * 100
    if perr > threshold:
        fail(f'Pressure at cell {i}:  CFD={p[i]:.6f}  Exact={pe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')

# --- 3. whole-domain relative L1 gates vs exact.txt --------------------------
# Measured on the x4 (400-cell) solve mesh: L1(rho)=6.85e-03, L1(p)=1.03e-02.
gate('whole-domain rel L1(rho) vs exact', rel_l1(rho, rhoe_i), 1.2e-2)
gate('whole-domain rel L1(p) vs exact',   rel_l1(p, pe_i),     1.8e-2)

print('Correct!')
