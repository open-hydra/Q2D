"""Sod79 x1/x2/x4/x8 refinement-ladder verification.

Two modes:

* ``verify.py --gate`` (ctest, after ``./Q2D.sh test``): numpy-only hard gate.
  Reads OUTPUT/field_x{1,2,4,8}.tec, evaluates the exact Riemann solution
  analytically at every cell center (test/common/pycheck/q2d_exact1d.py — the
  100-point reference table would smear the discontinuities by one table
  spacing and stall the x8 error), and asserts:
    - NaN-safe positivity of rho/p on every level;
    - whole-domain relative L1(rho) errors monotonically non-increasing;
    - observed order p_obs (least-squares slope of log L1(rho) vs log h)
      >= 0.70. Measured p_obs = 0.85 (calibration 2026-07, errors
      8.97e-3 / 4.72e-3 / 2.78e-3 / 1.50e-3); threshold ~0.85x measured —
      shock-dominated problem, formal order collapses to ~1 at best.
  No plotting and no matplotlib/ExactPack import in this mode.

* without ``--gate``: original offline study — ExactPack exact solution,
  2x2 matplotlib SVG figure, x8 probe-point gates, L2/Linf norm tables.
"""
import argparse
import sys
import os
import warnings
import numpy as np
from numpy import linspace
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_positive, rel_l1,
                           read_q2d_field)
from q2d_exact1d import riemann_exact, observed_order

# Define command-line arguments
parser = argparse.ArgumentParser(description='Plot CFD and Exact data')
parser.add_argument('--plot', action='store_true', help='Specify whether to plot the data or not')
parser.add_argument('--gate', action='store_true',
                    help='ctest mode: hard PASS/FAIL on the observed L1(rho) convergence order (no plotting)')

# Parse the command-line arguments
args = parser.parse_args()
output_dir = "OUTPUT"
folder_name = os.path.basename(os.getcwd())

# Sod (1978) shock-tube setup (matches input.ini / reference/exact.txt)
XMIN, XD0, XMAX, T_FINAL = 0.0, 0.5, 1.0, 0.20
RHOL, UL, PL = 1.0, 0.0, 1.0
RHOR, UR, PR = 0.125, 0.0, 0.1
LEVELS = (1, 2, 4, 8)
P_OBS_MIN = 0.70  # measured p_obs = 0.85 (2026-07); gate at ~0.85x measured


def _load_level(level):
    """Cell centers + named cell-centered fields of one ladder level."""
    path = 'OUTPUT/field_x%d.tec' % level
    require_file(path, "run './Q2D.sh test' first")
    x_, y_, var_, vnames = read_q2d_field(path)
    xc = 0.5 * (x_[0][:-1] + x_[0][1:])[:, 0, 0]
    cc = [n.strip() for n in vnames][3:]  # after nodal x, y, z
    fields = {name: np.asarray(var_[0][k], float)[:, 0, 0]
              for k, name in enumerate(cc)}
    return xc, fields


def run_gate():
    print('Sod79 convergence gate: L1(rho) vs analytic Riemann solution')
    errs, hs = [], []
    for level in LEVELS:
        xc, f = _load_level(level)
        rho, p = f['rho(1)'], f['p']
        assert_positive('rho (x%d)' % level, rho)
        assert_positive('p (x%d)' % level, p)
        ex = riemann_exact(xc, T_FINAL, XD0, RHOL, UL, PL, RHOR, UR, PR)
        err = rel_l1(rho, ex['rho'])
        h = (xc[-1] - xc[0]) / (len(xc) - 1)
        errs.append(err)
        hs.append(h)
        print('  x%-2d  n=%4d  h=%.6g  rel L1(rho) = %.5e' % (level, len(xc), h, err))
    for i in range(len(errs) - 1):
        if errs[i + 1] > errs[i]:
            print('  [FAIL] L1(rho) not monotonically non-increasing: '
                  'x%d -> x%d : %.5e -> %.5e' % (LEVELS[i], LEVELS[i + 1], errs[i], errs[i + 1]))
            fail('refinement did not reduce the L1(rho) error')
    print('  [PASS] L1(rho) monotonically non-increasing across the ladder')
    p_obs = observed_order(hs, errs)
    status = 'PASS' if p_obs >= P_OBS_MIN else 'FAIL'
    print('  [%s] observed order p_obs = %.3f (threshold %.2f, measured 0.85 at calibration)'
          % (status, p_obs, P_OBS_MIN))
    if p_obs < P_OBS_MIN:
        fail('observed L1(rho) order %.3f below threshold %.2f' % (p_obs, P_OBS_MIN))
    print('PASS: Sod79 convergence gate')
    return 0


if args.gate:
    sys.exit(run_gate())

# ---------------------------------------------------------------------------
# Offline study (original behavior): ExactPack + matplotlib, x8 probe gates.
# ---------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
from exactpack.solvers.riemann.ep_riemann import IGEOS_Solver

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

# Try to make the script work even when conda is not activated
root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break

try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "ORION package not found. Run 'conda activate base' or install ORION in the current Python environment."
    )

# Evaluate exact solution via ExactPack
xmin, xd0, xmax, t = 0.0, 0.5, 1.0, 0.20
rl, ul, pl, gl = 1.0, 0.0, 1.0, 1.4
rr, ur, pr, gr = 0.125, 0.0, 0.1, 1.4
solver_ig = IGEOS_Solver(xmin=xmin, xd0=xd0, xmax=xmax, t=t,
                            rl=rl, ul=ul, pl=pl, gl=gl,
                            rr=rr, ur=ur, pr=pr, gr=gr)
t = 0.20
nx = 100
xe = linspace(0.5*1.0/nx, 1.0-0.5*1.0/nx, nx)
soln1 = solver_ig(xe, t)
nx = 200
xe = linspace(0.5*1.0/nx, 1.0-0.5*1.0/nx, nx)
soln2 = solver_ig(xe, t)
nx = 400
xe = linspace(0.5*1.0/nx, 1.0-0.5*1.0/nx, nx)
soln4 = solver_ig(xe, t)
nx = 800
xe = linspace(0.5*1.0/nx, 1.0-0.5*1.0/nx, nx)
soln8 = solver_ig(xe, t)

# Read the CFD solution
# Name-based CC index dispatch (robust Q2D 2D ndir=2 no-w vs the 3D peer ndir=3).
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field_x1.tec')
_ndir = 2 if x_[0].shape[2] == 1 else 3
_cc = vnames[_ndir:]
def _i(name):
    return _cc.index(name if name in _cc else ('rho(1)' if name == 'rho' else name))
xc = 0.5 * (x_[0][:-1] + x_[0][1:])
x1 = xc; var1 = var_[0]
cv1 = var1[_i('R')]/(var1[_i('g')]-1)
T1 = var1[_i('T')]
ei1 = cv1*var1[_i('T')]
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field_x2.tec')
xc = 0.5 * (x_[0][:-1] + x_[0][1:])
x2 = xc; var2 = var_[0]
cv2 = var2[_i('R')]/(var2[_i('g')]-1)
T2 = var2[_i('T')]
ei2 = cv2*var2[_i('T')]
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field_x4.tec')
xc = 0.5 * (x_[0][:-1] + x_[0][1:])
x4 = xc; var4 = var_[0]
cv4 = var4[_i('R')]/(var4[_i('g')]-1)
T4 = var4[_i('T')]
ei4 = cv4*var4[_i('T')]
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field_x8.tec')
xc = 0.5 * (x_[0][:-1] + x_[0][1:])
x8 = xc; var8 = var_[0]
cv8 = var8[_i('R')]/(var8[_i('g')]-1)
T8 = var8[_i('T')]
ei8 = cv8*var8[_i('T')]

# Create a figure with 2x2 subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

# Density
axs[0, 0].plot(xe, soln8['density'], label='Exact', color='black')
axs[0, 0].plot(x1[:, 0, 0], var1[_i('rho(1)')][:, 0, 0], label='CFD x1', color='red')
axs[0, 0].plot(x8[:, 0, 0], var8[_i('rho(1)')][:, 0, 0], label='CFD x8', color='green')
axs[0, 0].set_xlabel('x')
axs[0, 0].set_ylabel('Density')
axs[0, 0].grid(True, alpha=0.3)
axs[0, 0].legend()

# Velocity
axs[0, 1].plot(xe, soln8['velocity'], label='Exact', color='black')
axs[0, 1].plot(x1[:, 0, 0], var1[_i('u')][:, 0, 0], label='CFD x1', color='red')
axs[0, 1].plot(x8[:, 0, 0], var8[_i('u')][:, 0, 0], label='CFD x8', color='green')
axs[0, 1].set_xlabel('x')
axs[0, 1].set_ylabel('Velocity')
axs[0, 1].grid(True, alpha=0.3)
axs[0, 1].legend()

# Pressure
axs[1, 0].plot(xe, soln8['pressure'], label='Exact', color='black')
axs[1, 0].plot(x1[:, 0, 0], var1[_i('p')][:, 0, 0], label='CFD x1', color='red')
axs[1, 0].plot(x8[:, 0, 0], var8[_i('p')][:, 0, 0], label='CFD x8', color='green')
axs[1, 0].set_xlabel('x')
axs[1, 0].set_ylabel('Pressure')
axs[1, 0].grid(True, alpha=0.3)
axs[1, 0].legend()

# SIE
axs[1, 1].plot(xe, soln8['specific_internal_energy'], label='Exact', color='black')
axs[1, 1].plot(x1[:, 0, 0], cv1[:,0,0]*T1[:, 0, 0], label='CFD x1', color='red')
axs[1, 1].plot(x8[:, 0, 0], cv8[:,0,0]*T8[:, 0, 0], label='CFD x8', color='green')
axs[1, 1].set_xlabel('x')
axs[1, 1].set_ylabel('SIE')
axs[1, 1].grid(True, alpha=0.3)
axs[1, 1].legend()

# Adjust layout to prevent overlap of subplots
plt.tight_layout()

out_file = os.path.join(output_dir, f"{folder_name}.svg")
plt.savefig(out_file, bbox_inches="tight", transparent=True)

if args.plot:
    # Display the figure with four distinct plots
    plt.show()

iprobe = [8*20,8*50,8*60,8*80,8*90]
l2_norm = np.zeros((4, 4))
linf_norm = np.zeros((4, 4))

# Density
l2_norm[0,0] = np.linalg.norm(var1[_i('rho(1)')][:,0,0] - soln1['density'])
l2_norm[0,1] = np.linalg.norm(var2[_i('rho(1)')][:,0,0] - soln2['density'])
l2_norm[0,2] = np.linalg.norm(var4[_i('rho(1)')][:,0,0] - soln4['density'])
l2_norm[0,3] = np.linalg.norm(var8[_i('rho(1)')][:,0,0] - soln8['density'])
linf_norm[0,0] = np.max(np.abs(var1[_i('rho(1)')][:,0,0] - soln1['density']))
linf_norm[0,1] = np.max(np.abs(var2[_i('rho(1)')][:,0,0] - soln2['density']))
linf_norm[0,2] = np.max(np.abs(var4[_i('rho(1)')][:,0,0] - soln4['density']))
linf_norm[0,3] = np.max(np.abs(var8[_i('rho(1)')][:,0,0] - soln8['density']))
for i in range(len(iprobe)):
    perr = np.abs((var8[_i('rho(1)')][iprobe[i],0,0] - soln8['density'][iprobe[i]]) / soln8['density'][iprobe[i]]) * 100
    if (perr > 1.0):
        print('Density')
        print(iprobe[i])
        print(var8[_i('rho(1)')][iprobe[i],0,0], soln8['density'][iprobe[i]])
        print('Error = ',perr, 'threshold = ',1.0)
        sys.exit(1)

# Velocity
l2_norm[1,0] = np.linalg.norm(var1[_i('u')][:,0,0] - soln1['velocity'])
l2_norm[1,1] = np.linalg.norm(var2[_i('u')][:,0,0] - soln2['velocity'])
l2_norm[1,2] = np.linalg.norm(var4[_i('u')][:,0,0] - soln4['velocity'])
l2_norm[1,3] = np.linalg.norm(var8[_i('u')][:,0,0] - soln8['velocity'])
linf_norm[1,0] = np.max(np.abs(var1[_i('u')][:,0,0] - soln1['velocity']))
linf_norm[1,1] = np.max(np.abs(var2[_i('u')][:,0,0] - soln2['velocity']))
linf_norm[1,2] = np.max(np.abs(var4[_i('u')][:,0,0] - soln4['velocity']))
linf_norm[1,3] = np.max(np.abs(var8[_i('u')][:,0,0] - soln8['velocity']))
for i in range(len(iprobe)):
    perr = np.abs((var8[_i('u')][iprobe[i],0,0] - soln8['velocity'][iprobe[i]]) / (soln8['velocity'][iprobe[i]]+1e-10)) * 100
    if (perr > 1.0):
        print('Velocity')
        print(iprobe[i])
        print(var8[_i('u')][iprobe[i],0,0], soln8['velocity'][iprobe[i]])
        print('Error = ',perr, 'threshold = ',1.0)
        sys.exit(1)

# Pressure
l2_norm[2,0] = np.linalg.norm(var1[_i('p')][:,0,0] - soln1['pressure'])
l2_norm[2,1] = np.linalg.norm(var2[_i('p')][:,0,0] - soln2['pressure'])
l2_norm[2,2] = np.linalg.norm(var4[_i('p')][:,0,0] - soln4['pressure'])
l2_norm[2,3] = np.linalg.norm(var8[_i('p')][:,0,0] - soln8['pressure'])
linf_norm[2,0] = np.max(np.abs(var1[_i('p')][:,0,0] - soln1['pressure']))
linf_norm[2,1] = np.max(np.abs(var2[_i('p')][:,0,0] - soln2['pressure']))
linf_norm[2,2] = np.max(np.abs(var4[_i('p')][:,0,0] - soln4['pressure']))
linf_norm[2,3] = np.max(np.abs(var8[_i('p')][:,0,0] - soln8['pressure']))
for i in range(len(iprobe)):
    perr = np.abs((var8[_i('p')][iprobe[i],0,0] - soln8['pressure'][iprobe[i]]) / soln8['pressure'][iprobe[i]]) * 100
    if (perr > 1.0):
        print('Pressure')
        print(iprobe[i])
        print(var8[_i('p')][iprobe[i],0,0], soln8['pressure'][iprobe[i]])
        print('Error = ',perr, 'threshold = ',1.0)
        sys.exit(1)

# SIE
l2_norm[3,0] = np.linalg.norm(ei1[:,0,0] - soln1['specific_internal_energy'])
l2_norm[3,1] = np.linalg.norm(ei2[:,0,0] - soln2['specific_internal_energy'])
l2_norm[3,2] = np.linalg.norm(ei4[:,0,0] - soln4['specific_internal_energy'])
l2_norm[3,3] = np.linalg.norm(ei8[:,0,0] - soln8['specific_internal_energy'])
linf_norm[3,0] = np.max(np.abs(ei1[:,0,0] - soln1['specific_internal_energy']))
linf_norm[3,1] = np.max(np.abs(ei2[:,0,0] - soln2['specific_internal_energy']))
linf_norm[3,2] = np.max(np.abs(ei4[:,0,0] - soln4['specific_internal_energy']))
linf_norm[3,3] = np.max(np.abs(ei8[:,0,0] - soln8['specific_internal_energy']))
for i in range(len(iprobe)):
    perr = np.abs((ei8[iprobe[i],0,0] - soln8['specific_internal_energy'][iprobe[i]]) / soln8['specific_internal_energy'][iprobe[i]]) * 100
    if (perr > 1.0):
        print('SIE')
        print(iprobe[i])
        print(ei8[iprobe[i],0,0], soln8['specific_internal_energy'][iprobe[i]])
        print('Error = ',perr, 'threshold = ',1.0)
        sys.exit(1)
