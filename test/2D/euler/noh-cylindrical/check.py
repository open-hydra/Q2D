"""noh-cylindrical verification check (Noh 1987 cylindrical implosion, m=1).

The ONLY axisymmetric case in the suite gated against an ANALYTIC exact
solution: it verifies the axisymmetric geometric source terms
(Lib_Source_Axisymmetry.f90, [Q2D-Sources] active = axisymmetric) in
isolation, on a Cartesian (x, r) strip where every 1/r term is exercised
and the answer is known in closed form.

Setup (gen_case.py; mirror of the 1D Noh87 choices): thermally-perfect
air (test/common/Air, constant Cp -> gamma = 1.40000000229 as written by
the solver itself in the 'g' output variable -- read back here so exact
and CFD use the SAME gamma), rho0 = 1, u_r = -U0 = -1, p0 = 1e-9,
r in [0, 2.5] resolved by 1000 cells (dr = 2.5e-3), 4-cell x-strip,
t_final = 1.0 (time-threshold; the exact solution is evaluated at the
actual output SOLUTIONTIME, which overshoots by < 0.1 %).

Exact solution (test/common/pycheck/q2d_exact1d.py::noh_exact, m = 1;
cross-checked there against gamma = 5/3: D = 1/3, rho2 = 16, p2 = 16/3):

  shock radius  r_s = D t,   D = (gamma-1)/2 U0 = 0.2
  r < r_s:  rho = rho0 ((g+1)/(g-1))^2 = 36,  u_r = 0,
            p = rho0 U0^2 (g+1)/2 (g+1)/(g-1) = 7.2
  r > r_s:  rho = rho0 (1 + U0 t / r),  u_r = -U0,  p ~ p0

Gates are evaluated on the x-averaged radial profile, restricted to the
validated region r <= R_VALID = 1.0: the zero-gradient outer boundary at
r = 2.5 mis-feeds the supersonic inflow (see gen_case.py) and its error
band only reaches r ~ 2.5 - U0 t = 1.5 by t_final (measured rel err
1.8e-4 at r = 1.4 vs 1.8e-2 at r = 2.0).

Calibration (NJ = 1000 run, 2026-07-15; first-order HLLE):
  1. NaN-safe positivity of rho/p/T, finiteness of u/v;
  2. whole-region rel L1(rho) vs exact on r <= 1: measured 4.77e-2
     (dominated by first-order shock smearing + the classic Noh axis
     wall-heating dip) -> tol 8e-2 (~1.7x headroom);
  3. post-shock plateau density, mean over r in (0.4, 0.9) r_s (window
     excludes both the axis wall-heating dip and the smeared shock):
     measured 35.09 vs 36, rel err 2.53e-2 -> tol 5e-2 (~2x);
  4. shock position (steepest density-gradient face) vs r_s = D t:
     measured offset +1.93 cells -> tol 4 cells;
  5. x-uniformity (negative control for spurious x-fluxes through the
     axisymmetric source): max over r-rows of the x-span of rho/v/p,
     relative to each field's max-abs: measured 3.2e-12 / 4.8e-12 /
     1.9e-12 -> tol 1e-8 (catches any real x-flux while allowing
     roundoff-level noise).
"""
import re
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_finite, assert_positive,
                           rel_l1, gate, read_q2d_field)
from q2d_exact1d import noh_exact

# --- Noh m=1 setup: keep in sync with gen_case.py ---------------------------
RHO0, U0, P0 = 1.0, 1.0, 1.0e-9
T_FINAL = 1.0          # input.ini time-threshold (SOLUTIONTIME read below)
R_VALID = 1.0          # gated region (outer-BC error band stays at r > 1.4)

TOL_L1_RHO   = 8.0e-2  # measured 4.77e-2 (2026-07-15) -> ~1.7x headroom
TOL_PLATEAU  = 5.0e-2  # measured 2.53e-2 -> ~2x headroom
TOL_SHOCK    = 4.0     # cells; measured offset +1.93 cells
TOL_XUNIF    = 1.0e-8  # measured 3.2e-12 (roundoff); real x-fluxes are >> this
PLATEAU_WIN  = (0.4, 0.9)  # plateau window in units of r_s (see docstring)

path = require_file('OUTPUT/field.tec', "run './Q2D.sh solve' first")

# Actual output time: the time-accurate loop stops at the first step with
# t >= time-threshold, so SOLUTIONTIME overshoots T_FINAL by up to one dt
# (< 0.1 %); evaluating the exact solution there removes that bias.
zone = open(str(path)).readlines()[1]
m = re.search(r'SOLUTIONTIME=\s*([0-9.Ee+-]+)', zone)
if m is None:
    fail('OUTPUT/field.tec zone record has no SOLUTIONTIME')
t_out = float(m.group(1))
if not (T_FINAL <= t_out < 1.01 * T_FINAL):
    fail('output time %.6f not in [%g, %g) -- run did not reach t_final?'
         % (t_out, T_FINAL, 1.01 * T_FINAL))

[x_, y_, var_, vnames] = read_q2d_field('OUTPUT/field.tec')
yn = y_[0][:, :, 0]                     # nodes, [i, j]
# Q2D_new cell-centered order: rho(1), u, v, p, T, g, R  (v = radial velocity)
rho = var_[0][0][:, :, 0]
u   = var_[0][1][:, :, 0]
v   = var_[0][2][:, :, 0]
p   = var_[0][3][:, :, 0]
T   = var_[0][4][:, :, 0]
gam = var_[0][5][:, :, 0]

# --- 1. NaN-safe sanity ------------------------------------------------------
assert_positive('rho', rho)
assert_positive('p', p)
assert_positive('T', T)
assert_finite('u', u)
assert_finite('v', v)

# gamma as the solver wrote it (constant-Cp air => uniform field)
assert_finite('g', gam)
if gam.max() - gam.min() > 1e-9:
    fail('gamma field not uniform (spread %.3e): constant-Cp air expected'
         % (gam.max() - gam.min()))
g = float(gam.mean())

# --- exact solution on the x-averaged radial profile -------------------------
rc = 0.5 * (yn[0, :-1] + yn[0, 1:])     # cell-center radii (uniform in i)
h  = rc[1] - rc[0]
rho_prof = rho.mean(axis=0)             # x-average (axis 0 = i)
ex = noh_exact(rc, t_out, x0=0.0, rho0=RHO0, u0=U0, p0=P0, gamma=g, m=1)
r_s = ex['shock_right']                 # = D * t_out, D = (g-1)/2 * U0
rho_post = float(ex['rho'][0])          # = RHO0 ((g+1)/(g-1))^2
print('exact (gamma=%.6f, t=%.6f): D=%.6f, r_s=%.6f, rho_post=%.4f, p_post=%.4f'
      % (g, t_out, ex['shock_speed'], r_s, rho_post, float(ex['p'][0])))

valid = rc <= R_VALID

# --- 2. whole-region rel L1(rho) vs exact (r <= R_VALID) ---------------------
gate('rel L1(rho) vs exact, r <= %g' % R_VALID,
     rel_l1(rho_prof[valid], ex['rho'][valid]), TOL_L1_RHO)

# --- 3. post-shock plateau density vs analytic value -------------------------
win = (rc > PLATEAU_WIN[0] * r_s) & (rc < PLATEAU_WIN[1] * r_s)
if not win.any():
    fail('empty plateau window (%g, %g) r_s -- shock not where expected?'
         % PLATEAU_WIN)
rho_plateau = float(rho_prof[win].mean())
gate('post-shock plateau rho (window (%.1f, %.1f) r_s) vs %.4f' %
     (PLATEAU_WIN[0], PLATEAU_WIN[1], rho_post),
     abs(rho_plateau - rho_post) / rho_post, TOL_PLATEAU)

# --- 4. shock position vs D t -------------------------------------------------
isk = int(np.argmax(np.abs(np.diff(rho_prof))))
r_shock = 0.5 * (rc[isk] + rc[isk + 1])
if abs(r_shock - r_s) > TOL_SHOCK * h:
    fail('shock at r=%.5f, expected %.5f +/- %g cells (h=%.4g)'
         % (r_shock, r_s, TOL_SHOCK, h))
print('  [PASS] shock position: r=%.5f (exact %.5f, offset %.2f cells, tol %g)'
      % (r_shock, r_s, (r_shock - r_s) / h, TOL_SHOCK))

# --- 5. x-uniformity (spurious x-flux negative control) ----------------------
for name, f in (('rho', rho), ('v', v), ('p', p)):
    span = float((f.max(axis=0) - f.min(axis=0)).max() / np.abs(f).max())
    gate('x-uniformity rel span(%s)' % name, span, TOL_XUNIF)

print('Correct!')
