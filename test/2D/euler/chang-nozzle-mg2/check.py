"""chang-nozzle-mg2 2D steady-state regression check (Multigrid, 2 levels).

Same Chang 1996 axisymmetric nozzle as ``test/2D/euler/chang-nozzle`` (cf.
Chang 1996 J. Spacecraft Rockets; Q2D peculiarita' #1 j-radial + #2
axisymmetric source), but run with 2-level Multigrid
(``[Q2D-Multigrid] levels=2``). Mirror of the planar-nozzle-mg2 3-layer
acceptance: the physics gates are the SAME as the MGL=1 baseline, plus a
multigrid-consistency layer.

  1. (convergence)  hard sentinel: final density residual < 1e-6
     (measured 1.0e-10 on the committed run). A broken coarse-level
     operator stalls/diverges instead of converging.
  2. (physical)     the converged fine field matches the SAME Chang 1996
     digitized references (ref1.dat axis, ref2.dat wall) within the SAME
     tolerances as chang-nozzle -> multigrid does not corrupt the physics.
  3. (multigrid)    OUTPUT/field-level2.tec exists (hard), coarse NI/NJ are
     exactly fine NI/2, NJ/2, and the coarse-level Mach reproduces the
     2x2-pairwise-restricted fine Mach: max rel diff < 2e-1, mean < 5e-2
     (measured 8.1e-2 max -- at the axis/throat singularity row -- and
     5.2e-3 mean on the committed run; a scrambled coarse restriction is
     off by O(1) everywhere).

Field sanity per ``docs/vv/regression_tolerance.md`` Section 3.2:
``min(rho), min(p), min(T) > 0`` (hard, NaN-safe via q2d_check_lib).
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

# ---------- layer 1: convergence sentinel (hard) ----------
# Final density residual (column 2 of OUTPUT/residual-history.dat) < 1e-6.
# Measured 1.0e-10 on the committed run; a broken coarse level (e.g. a
# scrambled coarse geometric source) stalls or diverges instead.
assert_converged(1.0e-6)

# ---------- layer 2: physical -- same Chang 1996 gates as chang-nozzle ----------
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
# cells (cf. chang-nozzle check.py, note T-1.2D.axis-discrepancy).
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

# ---------- layer 3: multigrid -- coarse level reproduces the fine restriction ----------
# field-level2.tec must exist (hard: MGL=2 path did not run otherwise) and be
# an exact 2:1 coarsening of the fine grid in both directions.
require_file('OUTPUT/field-level2.tec',
             "the coarse level was not written; MGL=2 path did not run")
xc_, yc_, varc_, vnamesc = read_q2d_field('OUTPUT/field-level2.tec')
rho_c = varc_[0][0][:, :, 0]
u_c   = varc_[0][1][:, :, 0]
v_c   = varc_[0][2][:, :, 0]
p_c   = varc_[0][3][:, :, 0]
T_c   = varc_[0][4][:, :, 0]
gam_c = varc_[0][5][:, :, 0]

NI, NJ = M.shape
NIc, NJc = rho_c.shape[:2]
if NIc != NI // 2 or NJc != NJ // 2:
    fail(f"coarse cells ({NIc},{NJc}) are not fine cells/2 ({NI//2},{NJ//2}): "
         "2:1 coarsening broken")

assert_positive("coarse rho", rho_c)
assert_positive("coarse p", p_c)
assert_positive("coarse T", T_c)
M_c = np.sqrt(u_c**2 + v_c**2) / np.sqrt(gam_c * p_c / rho_c)
assert_finite("coarse Mach", M_c)

# 2x2 pairwise restriction of the fine Mach onto the coarse cells. A correct
# coarse level tracks this closely; a scrambled restriction (B-10-type bug)
# is off by O(1). Measured on the committed run: max 8.1e-2 (axis/throat
# row, the axisymmetric-singularity cells), mean 5.2e-3. Tolerances:
# max 2.0e-1 (~2.5x measured -- the axis row caps how tight this can be),
# mean 5.0e-2 (~10x measured, same headroom policy as planar-nozzle-mg2).
M_restr = M.reshape(NIc, 2, NJc, 2).mean(axis=(1, 3))
cf_rel = np.abs(M_c - M_restr) / np.maximum(np.abs(M_restr), 1e-30)
cf_max = float(cf_rel.max())
cf_mean = float(cf_rel.mean())
print(f'Multigrid coarse/fine Mach: rel diff max={cf_max:.3e} mean={cf_mean:.3e} '
      f'(coarse {NIc}x{NJc} vs fine {NI}x{NJ})')
if cf_max > 2.0e-1:
    fail(f'coarse/fine Mach max rel {cf_max:.3e} > 2.0e-1 '
         '(coarse level does not reproduce the fine restriction)')
if cf_mean > 5.0e-2:
    fail(f'coarse/fine Mach mean rel {cf_mean:.3e} > 5.0e-2 '
         '(coarse level does not reproduce the fine restriction)')

print('Correct!')
