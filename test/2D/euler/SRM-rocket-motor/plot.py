"""SRM-rocket-motor V&V figure -> docs/vv/images/SRM-pressure.svg.

Left: axial pressure profile (j-averaged) of the steady SRM field
(OUTPUT/field.tec), the chamber-average pressure gated by check.py, and the
analytic Vieille equilibrium level from reference/P_eq_vieille.dat with the
15% acceptance band shaded. Throat located from MESH/profile.dat (as in
check.py). The solver writes no pressure-vs-iteration history, so the
time-side evidence (right) is the density-residual history: the decay of
the chamber-pressurisation transient toward the Vieille
equilibrium, with the check.py convergence sentinel.

Run './Q2D.sh solve' first, then 'python3 plot.py'.
"""
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

HERE = Path(__file__).resolve().parent           # <repo>/test/2D/euler/SRM-rocket-motor
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "SRM-pressure.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file, read_q2d_field, read_residual_history

# check.py gates: Vieille acceptance band and convergence sentinel.
VIEILLE_TOL = 0.15
RESID_THR = 0.5
CHAMBER_MARGIN = 0.005                           # chamber = x < x_throat - 0.005 m

x_, y_, var_, vnames = read_q2d_field(require_file(HERE / "OUTPUT" / "field.tec",
                                                   "run './Q2D.sh solve' first"))
xn = x_[0][:, :, 0]
cc = [n.strip() for n in vnames][3:]             # after nodal x, y, z
p = np.asarray(var_[0][cc.index("p")], float)[:, :, 0]

xc = 0.5 * (xn[:-1, :-1] + xn[1:, 1:])[:, 0]

profile = np.loadtxt(require_file(HERE / "MESH" / "profile.dat"))
x_throat = float(profile[np.argmin(profile[:, 1]), 0])

mask_chamber = xc < (x_throat - CHAMBER_MARGIN)
p_chamber = float(p[mask_chamber, :].mean())

ref = np.loadtxt(require_file(HERE / "reference" / "P_eq_vieille.dat"),
                 comments="#")
P_eq = float(ref if ref.ndim == 0 else ref[0])
eps = abs(p_chamber - P_eq) / P_eq

res = read_residual_history(str(HERE / "OUTPUT" / "residual-history.dat"))

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))

ax = axs[0]
ax.axhspan(P_eq * (1 - VIEILLE_TOL) / 1e6, P_eq * (1 + VIEILLE_TOL) / 1e6,
           color="tab:green", alpha=0.15, label="acceptance band (15%)")
ax.axhline(P_eq / 1e6, color="k", ls="--", lw=1.2,
           label="Vieille equilibrium P_eq")
ax.plot(xc, p.mean(axis=1) / 1e6, "-", color="tab:red", lw=1.5,
        label="Q2D p(x), radially averaged")
ax.axhline(p_chamber / 1e6, color="tab:blue", ls="-.", lw=1.2,
           label="chamber average (rel err %.1e)" % eps)
ax.axvline(x_throat, color="0.5", ls=":", lw=1)
ax.text(x_throat, 0.03, " throat", fontsize=8, color="0.4",
        transform=ax.get_xaxis_transform())
ax.set_xlabel("x [m]")
ax.set_ylabel("p [MPa]")
ax.set_title("Chamber pressure vs Vieille equilibrium")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="center left")

ax = axs[1]
ax.semilogy(res[:, 0], res[:, 2], "-", color="tab:red", lw=1.2,
            label="density residual")
ax.axhline(RESID_THR, color="k", ls="--", lw=1,
           label="convergence sentinel (check.py)")
ax.set_xlabel("iteration")
ax.set_ylabel("density residual")
ax.set_title("Pressurisation transient (to the steady state)")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=8, loc="lower right")

fig.suptitle("SRM rocket motor (BC 502, Vieille burn rate): Q2D vs analytic "
             "equilibrium")
plt.tight_layout(rect=(0, 0, 1, 0.94))
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("P_chamber %.4e vs P_eq %.4e (rel err %.3f)" % (p_chamber, P_eq, eps))
print("Figure saved to %s" % OUT_SVG)
