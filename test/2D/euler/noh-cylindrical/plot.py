"""noh-cylindrical V&V figure -> docs/vv/images/NohCyl.svg.

x-averaged radial profiles of rho, u_r and p from the converged
noh-cylindrical run (OUTPUT/field.tec) against the analytic cylindrical
Noh (1987) m=1 exact solution (test/common/pycheck/q2d_exact1d.py::
noh_exact) at the actual output time -- the same data pipeline as
check.py.  Profiles are shown on the validated region r <= 1 (the
extrapolation outer-boundary band lives at r > 1.4, see gen_case.py).

Run './Q2D.sh solve' first, then 'python3 plot.py'.
"""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

# House style: transparent, theme-aware SVG (same rcParams as the repo's
# committed plot/verify scripts); fixed hashsalt for reproducible SVG ids.
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
    "svg.hashsalt": "Q2D-vv",
})

HERE = Path(__file__).resolve().parent           # <repo>/test/2D/euler/noh-cylindrical
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "NohCyl.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file, read_q2d_field
from q2d_exact1d import noh_exact

# Setup constants: keep in sync with gen_case.py / check.py.
RHO0, U0, P0 = 1.0, 1.0, 1.0e-9
R_VALID = 1.0

path = require_file(HERE / "OUTPUT" / "field.tec", "run './Q2D.sh solve' first")
zone = open(str(path)).readlines()[1]
t_out = float(re.search(r'SOLUTIONTIME=\s*([0-9.Ee+-]+)', zone).group(1))

x_, y_, var_, vnames = read_q2d_field(path)
yn = y_[0][:, :, 0]
rho = var_[0][0][:, :, 0]
v = var_[0][2][:, :, 0]                          # v = radial velocity
p = var_[0][3][:, :, 0]
g = float(var_[0][5][:, :, 0].mean())            # uniform constant-Cp gamma

# same reduction as check.py: x-average (axis 0 = i), cell-center radii
rc = 0.5 * (yn[0, :-1] + yn[0, 1:])
sel = rc <= R_VALID
rc = rc[sel]
rho_prof = rho.mean(axis=0)[sel]
v_prof = v.mean(axis=0)[sel]
p_prof = p.mean(axis=0)[sel]

r_fine = np.linspace(rc[0], R_VALID, 2000)
ex = noh_exact(r_fine, t_out, x0=0.0, rho0=RHO0, u0=U0, p0=P0, gamma=g, m=1)
r_s = ex["shock_right"]

fig, axs = plt.subplots(1, 3, figsize=(12, 4.0))
panels = (
    (rho_prof, ex["rho"], r"$\rho$ [kg/m$^3$]", "density"),
    (v_prof, ex["u"], r"$u_r$ [m/s]", "radial velocity"),
    (p_prof, ex["p"], r"$p$ [Pa]", "pressure"),
)
for ax, (q_cfd, q_ex, ylab, title) in zip(axs, panels):
    ax.plot(r_fine, q_ex, "k--", lw=1.5, label="exact (Noh 1987, m=1)")
    ax.plot(rc, q_cfd, "o", color="tab:red", ms=2.5, mfc="none",
            label="Q2D, x-averaged")
    ax.axvline(r_s, color="tab:blue", ls=":", lw=1,
               label=r"shock $r_s = D\,t$ = %.3f" % r_s)
    ax.set_xlabel("r [m]")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")

fig.suptitle("Cylindrical Noh implosion (axisymmetric source verification): "
             r"Q2D vs exact at $t$ = %.3f s ($\gamma$ = %.4f)" % (t_out, g))
plt.tight_layout(rect=(0, 0, 1, 0.93))
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("shock at r=%.4f, post-shock rho=%.4f (exact %.4f)"
      % (r_s, rho_prof[rc < 0.5 * r_s].mean(), ex["rho"][0]))
print("Figure saved to %s" % OUT_SVG)
