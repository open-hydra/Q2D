"""prandtl-meyer 2D steady-state regression check.

Validates the converged solution of the Prandtl-Meyer expansion fan against
the analytical isentropic expansion result (cf. Prandtl-Meyer 1908). Mirror
MOSE/test/2D/euler/prandtl-meyer/check.py.

Reads the converged Mach at the exit-probe location from ``OUTPUT/exit.txt``
(Q2D-Probes section in input.ini). Tolerance: |M_CFD - M_analytical| / M_analytical
< 0.5%.

Also enforces sanity check on the field (regression_tolerance.md Section 3.2):
``min(rho), min(p), min(T) > 0`` everywhere (hard, NaN-safe via q2d_check_lib).
"""
import math
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_positive, assert_finite,
                           read_q2d_field)

gamma = 1.4


def nufunc(M, gamma):
    return math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(
        math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) - math.atan(math.sqrt(M ** 2 - 1))


def mach_from_nu(nu, gamma, a=1.01, b=30.0, tol=1e-6, max_iter=1000000):
    for _ in range(max_iter):
        M_mid = (a + b) / 2
        nu_err = nufunc(M_mid, gamma) - nu
        if abs(nu_err) < tol:
            return M_mid
        if np.sign(nu_err) < 0:
            a = M_mid
        else:
            b = M_mid
    raise RuntimeError("Bisection did not converge.")


def prandtl_meyer_expansion(theta, M):
    nu1 = nufunc(M, gamma)
    nu2 = nu1 - theta * np.pi / 180 if theta < 0 else nu1 + theta * np.pi / 180
    return mach_from_nu(nu2, gamma)


# Field-level sanity check (hard, NaN-safe)
require_file('OUTPUT/field.tec', "did ./Q2D.sh solve run?")
x_, y_, var_, vnames = read_q2d_field('OUTPUT/field.tec')

rho = var_[0][0]
u   = var_[0][1]
v   = var_[0][2]
p   = var_[0][3]
T   = var_[0][4]
assert_positive("rho", rho)
assert_positive("p", p)
assert_positive("T", T)
assert_finite("u", u)
assert_finite("v", v)

# Analytical Prandtl-Meyer Mach at the exit
M1 = 2.0
theta = 34.971
M2_analytical = prandtl_meyer_expansion(theta, M1)

# CFD Mach at exit probe (last line of probe output). Guarded: the probe
# file must exist, be non-empty and its last line must parse to
# ``iter  Mach`` -- an empty/truncated probe file is a hard failure, not
# an IndexError.
require_file('OUTPUT/exit.txt', "exit probe was not written (Q2D-Probes)")
with open('OUTPUT/exit.txt', 'r') as file:
    lines = [ln for ln in file.readlines() if ln.strip()]
if not lines:
    fail("OUTPUT/exit.txt contains no data lines (probe never sampled)")
columns = lines[-1].strip().split()
if len(columns) < 2:
    fail(f"OUTPUT/exit.txt last line has {len(columns)} column(s), expected >= 2: "
         f"'{lines[-1].strip()}'")
try:
    CFDval = float(columns[1])
except ValueError:
    fail(f"OUTPUT/exit.txt last line is not numeric: '{lines[-1].strip()}'")
if not math.isfinite(CFDval):
    fail(f"exit-probe Mach is non-finite: {CFDval}")

perr = abs(CFDval - M2_analytical) / M2_analytical * 100
if perr > 0.5:
    print(f'CFD={CFDval}  Analytical={M2_analytical}')
    print(f'Error = {perr}%, threshold = 0.5%')
    fail('Mach number mismatch (CFD vs analytical)')
print(f'Natural reference: M_CFD={CFDval}, M_analytical={M2_analytical:.6f}, '
      f'eps_rel={perr:.4f}%')

print('Correct!')
