"""chang-nozzle-mg2 V&V figure -> docs/vv/images/CN-mg2.svg.

Left: converged fine-grid Mach along the nozzle axis (j = 0) and wall
(j = last) against the Chang 1996 digitized reference points
(reference/ref1.dat axis, ref2.dat wall) — the check.py layer-2 physics gate.
Right: the check.py layer-3 multigrid gate — coarse-level Mach
(OUTPUT/field-level2.tec) vs the 2x2 pairwise restriction of the fine Mach,
compared per j-averaged axial profile and quantified by the cellwise
relative-difference statistics against the acceptance thresholds.

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

HERE = Path(__file__).resolve().parent           # <repo>/test/2D/euler/chang-nozzle-mg2
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "CN-mg2.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file, read_q2d_field

# check.py layer-3 acceptance on the coarse/fine Mach restriction.
CF_TOL_MAX, CF_TOL_MEAN = 2.0e-1, 5.0e-2


def load(path):
    """(nodal x [i,j], cell-center x [i,j], Mach [i,j]) of a Q2D field file."""
    x_, y_, var_, vnames = read_q2d_field(require_file(path))
    xn = x_[0][:, :, 0]
    rho = var_[0][0][:, :, 0]
    u = var_[0][1][:, :, 0]
    v = var_[0][2][:, :, 0]
    p = var_[0][3][:, :, 0]
    gam = var_[0][5][:, :, 0]
    M = np.sqrt(u**2 + v**2) / np.sqrt(gam * p / rho)
    xc = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
    return xn, xc, M


xn, xc, M = load(HERE / "OUTPUT" / "field.tec")
_, xc_c, M_c = load(HERE / "OUTPUT" / "field-level2.tec")

# 2x2 pairwise restriction of the fine Mach (as in check.py layer 3)
NIc, NJc = M_c.shape
M_restr = M.reshape(NIc, 2, NJc, 2).mean(axis=(1, 3))
cf_rel = np.abs(M_c - M_restr) / np.maximum(np.abs(M_restr), 1e-30)

ref1 = np.loadtxt(require_file(HERE / "reference" / "ref1.dat"))  # axis
ref2 = np.loadtxt(require_file(HERE / "reference" / "ref2.dat"))  # wall

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))

ax = axs[0]
ax.plot(xc[:, 0], M[:, 0], "-", color="tab:red", lw=1.5, label="Q2D axis")
ax.plot(ref1[::4, 0], ref1[::4, 1], "o", color="k", ms=4, mfc="none",
        label="Chang 1996 axis (digitized)")
ax.plot(xc[:, -1], M[:, -1], "-", color="tab:blue", lw=1.5, label="Q2D wall")
ax.plot(ref2[::4, 0], ref2[::4, 1], "s", color="k", ms=4, mfc="none",
        label="Chang 1996 wall (digitized)")
ax.set_xlabel("x [m]")
ax.set_ylabel("Mach")
ax.set_title("Fine grid (MG 2-level) vs Chang 1996")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper left")

ax = axs[1]
x_prof = xc_c.mean(axis=1)
ax.plot(x_prof, M_restr.mean(axis=1), "-", color="tab:red", lw=1.5,
        label="fine Mach, 2x2 restricted (j-avg)")
ax.plot(x_prof[::2], M_c.mean(axis=1)[::2], "o", color="tab:blue", ms=4,
        mfc="none", label="coarse level (j-avg)")
ax.set_xlabel("x [m]")
ax.set_ylabel("Mach")
ax.set_title("MG consistency: coarse level vs restricted fine")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper left")
stats = ("cellwise rel diff:\n"
         "max  %.2e (tol %.0e)\n"
         "mean %.2e (tol %.0e)"
         % (cf_rel.max(), CF_TOL_MAX, cf_rel.mean(), CF_TOL_MEAN))
ax.text(0.97, 0.05, stats, transform=ax.transAxes, fontsize=8,
        family="monospace", va="bottom", ha="right",
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8))

fig.suptitle("Chang nozzle (2D-axi), 2-level multigrid: Q2D vs reference")
plt.tight_layout(rect=(0, 0, 1, 0.94))
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("coarse/fine rel diff: max %.3e mean %.3e" % (cf_rel.max(), cf_rel.mean()))
print("Figure saved to %s" % OUT_SVG)
