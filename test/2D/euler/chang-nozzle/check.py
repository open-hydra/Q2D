"""chang-nozzle 2D steady-state regression check.

Validates the converged solution of the Chang 1996 axisymmetric nozzle (cf.
Chang 1996 J. Spacecraft Rockets). Mirror MOSE/test/2D/euler/chang-nozzle/
verify.py adapted to the Q2D 2D no-w layout (j is the radial direction in
axisymmetric mode; single collapsed K strip) with the axisymmetric source
active in input.ini [Q2D-Sources].

Tolerances per ``docs/vv/regression_tolerance.md`` Section 3.2:
- ``min(rho), min(p), min(T) > 0`` (hard, NaN-safe via q2d_check_lib).
- convergence sentinel (hard): final density residual < 1e-6
  (measured 1.0e-10 on the committed run; the solver's own stopping
  criterion drives it to ~1e-10).
- Mach distribution on the wall vs Chang 1996 reference ref2.dat and on
  the axis vs ref1.dat: avg eps_rel < 5e-2, max eps_rel < 1e-1.
"""
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_positive, assert_finite,
                           assert_converged, read_q2d_field)

require_file('OUTPUT/field.tec', "did ./Q2D.sh solve run?")
x_, y_, var_, vnames = read_q2d_field('OUTPUT/field.tec')

xn  = x_[0][:, :, 0]
yn  = y_[0][:, :, 0]
rho = var_[0][0][:, :, 0]
u   = var_[0][1][:, :, 0]
v   = var_[0][2][:, :, 0]
p   = var_[0][3][:, :, 0]
T   = var_[0][4][:, :, 0]
gam = var_[0][5][:, :, 0]
R   = var_[0][6][:, :, 0]

# Sanity (hard, NaN-safe): positivity of rho, p, T; finiteness of u, v.
assert_positive("rho", rho)
assert_positive("p", p)
assert_positive("T", T)
assert_finite("u", u)
assert_finite("v", v)

a = np.sqrt(gam * p / rho)
M = np.sqrt(u**2 + v**2) / a
assert_finite("Mach", M)

# Convergence sentinel (hard): final density residual (column 2 of
# OUTPUT/residual-history.dat) below 1e-6. Measured 1.0e-10 on the
# committed run -- 1e-6 is a regression band a stalled run blows through.
assert_converged(1.0e-6)

# Mach along wall (j=last index, far from axis) vs Chang 1996 reference ref2.dat
x_wall = xn[: M.shape[0], -1]
M_wall = M[:, -1]
ref2 = np.loadtxt(require_file('reference/ref2.dat'))
M_ref_wall = np.interp(x_wall, ref2[:, 0], ref2[:, 1])

# Mach along axis (j=0) vs Chang 1996 reference ref1.dat
x_axis = xn[: M.shape[0], 0]
M_axis = M[:, 0]
ref1 = np.loadtxt(require_file('reference/ref1.dat'))
M_ref_axis = np.interp(x_axis, ref1[:, 0], ref1[:, 1])

# Tolerance: avg eps_rel < 5e-2 (regression_tolerance.md Section 3.2 wall avg).
# Max eps_rel checked vs separate threshold 1e-1, because Q2D 2D axisym source
# (1/r with eps=1e-12 floor) and MOSE 3D wedge slab (dz→0 at axis) are two
# distinct approximations of axisymmetric physics that differ ~4% at axis-row
# cells. Inherent to the Q2D 2D axisymmetric vs MOSE 3D wedge
# discretization.
def report_eps(label, sim, ref):
    denom = np.maximum(np.abs(ref), 1e-3)
    eps_rel = np.abs(sim - ref) / denom
    eps_max = eps_rel.max()
    eps_avg = eps_rel.mean()
    print(f'Mach {label}: eps_rel max={eps_max:.3e} avg={eps_avg:.3e}')
    return eps_max, eps_avg

eps_wall_max, eps_wall_avg = report_eps('wall', M_wall, M_ref_wall)
eps_axis_max, eps_axis_avg = report_eps('axis', M_axis, M_ref_axis)

thr_avg = 5.0e-2   # avg-cell tolerance
thr_max = 1.0e-1   # max-cell tolerance (relaxed for axis-singularity row)
if eps_wall_avg > thr_avg:
    fail(f'Mach wall avg mismatch: {eps_wall_avg:.3e} > {thr_avg}')
if eps_axis_avg > thr_avg:
    fail(f'Mach axis avg mismatch: {eps_axis_avg:.3e} > {thr_avg}')
if eps_wall_max > thr_max:
    fail(f'Mach wall max mismatch: {eps_wall_max:.3e} > {thr_max}')
if eps_axis_max > thr_max:
    fail(f'Mach axis max mismatch: {eps_axis_max:.3e} > {thr_max}')

print('Correct!')
