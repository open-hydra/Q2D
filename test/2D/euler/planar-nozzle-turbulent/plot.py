"""planar-nozzle-turbulent V&V figure -> docs/vv/images/SSGLRR-axis.svg.

RSM (SSG/LRR-omega) vs SST centreline-pressure overlay (the verify.py G3
mean-flow-consistency content) plus the Reynolds-stress anisotropy that
develops along the centreline from the exactly isotropic IC (the G5
discriminator vs the Boussinesq sibling). Data extraction is reused from
verify.py (read_field) — no duplicated parsing or constants.

Run './Q2D.sh solve' in kw/ and rsm/ first, then 'python3 plot.py'.
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

HERE = Path(__file__).resolve().parent  # <repo>/test/2D/euler/planar-nozzle-turbulent
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "SSGLRR-axis.svg"

# Reuse the family verify.py parsing (its module body only defines helpers;
# main() is __main__-guarded). It also sets up the pycheck import path.
sys.path.insert(0, str(HERE))
from verify import read_field

(kw, _, _), (rsm, ni, nj) = read_field("kw"), read_field("rsm")

jc = kw["p"].shape[0] // 2                       # centreline row (as in G3)
x_ax = rsm["xc"][jc, :]
p_kw, p_rsm = kw["p"][jc, :], rsm["p"][jc, :]
rel_linf = float(np.max(np.abs(p_rsm - p_kw) / np.abs(p_kw)))

# anisotropy tensor a_ij = R_ij/k - 2/3 delta_ij (as in G5)
R11, R22, R33, R12 = rsm["ruu"], rsm["rvv"], rsm["rww"], rsm["ruv"]
k = 0.5 * (R11 + R22 + R33)
a11 = R11 / k - 2.0 / 3.0
a22 = R22 / k - 2.0 / 3.0
a33 = R33 / k - 2.0 / 3.0
a12 = R12 / k

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))

ax = axs[0]
ax.plot(x_ax, p_kw / 1e5, "-", color="tab:blue", lw=1.5, label="SST (k-omega)")
ax.plot(x_ax, p_rsm / 1e5, "--", color="tab:red", lw=1.5,
        label="SSG/LRR-omega (RSM)")
ax.set_xlabel("x [m]")
ax.set_ylabel("p [bar] (centreline)")
ax.set_title("Mean flow: centreline pressure "
             "(max rel diff %.1e)" % rel_linf)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper right")

ax = axs[1]
for arr, lbl, color in ((a11, "a11", "tab:blue"),
                        (a22, "a22", "tab:orange"),
                        (a33, "a33 (depth)", "tab:green"),
                        (a12, "a12 (shear)", "tab:red")):
    ax.plot(x_ax, arr[jc, :], "-", color=color, lw=1.2, label=lbl)
ax.axhline(0.0, color="0.5", ls=":", lw=1)
ax.set_xlabel("x [m]")
ax.set_ylabel("a_ij = R_ij/k - 2/3 d_ij (centreline)")
ax.set_title("RSM anisotropy (IC was exactly isotropic)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper left")

fig.suptitle("Q2D area-variation turbulent nozzle: SSG/LRR-omega vs SST")
plt.tight_layout(rect=(0, 0, 1, 0.94))
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("centreline p max rel diff (rsm vs kw) = %.3e" % rel_linf)
print("Figure saved to %s" % OUT_SVG)
