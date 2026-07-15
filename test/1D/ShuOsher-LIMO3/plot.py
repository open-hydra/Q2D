"""ShuOsher-LIMO3 V&V figure -> docs/vv/images/ShuOsher.svg.

Post-shock density profile of the fresh LIMO3 run (OUTPUT/field.tec) against
the four committed limiter baselines (reference/baseline_{mc,vanleer,minmod,
iord}.dat), with a zoom on the resolved entropy-wave (sine) packet and the
post-shock total-variation (TV) values computed exactly as in check.py.

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

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]                      # <repo>/test/1D/ShuOsher-LIMO3
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "ShuOsher.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file, read_q2d_field

# Post-shock smooth-wave TV window and analytic shock location:
# keep in sync with check.py (WIN_LO/WIN_HI, SHOCK_X).
WIN_LO, WIN_HI = -3.0, 2.3
SHOCK_X = 2.4


def total_variation(x, r, lo, hi):
    """Same TV metric as check.py: sum |diff(rho)| inside [lo, hi]."""
    seg = r[(x >= lo) & (x <= hi)]
    return float(np.sum(np.abs(np.diff(seg))))


# --- fresh LIMO3 run ---------------------------------------------------------
field = HERE / "OUTPUT" / "field.tec"
require_file(field, "run './Q2D.sh solve' first")
x_, y_, var_, vnames = read_q2d_field(field)
xc = 0.5 * (x_[0][:-1] + x_[0][1:])[:, 0, 0]
cc = [n.strip() for n in vnames][3:]             # after nodal x, y, z
rho = np.asarray(var_[0][cc.index("rho(1)")], float)[:, 0, 0]

# --- committed limiter baselines (x, rho at N=400, t=1.8) --------------------
baselines = {}
for tag in ("mc", "vanleer", "minmod", "iord"):
    arr = np.loadtxt(require_file(HERE / "reference" / ("baseline_%s.dat" % tag)))
    baselines[tag] = (arr[:, 0], arr[:, 1])

tv = {"LIMO3": total_variation(xc, rho, WIN_LO, WIN_HI)}
for tag, (xb, rb) in baselines.items():
    tv[tag] = total_variation(xb, rb, WIN_LO, WIN_HI)

# --- figure ------------------------------------------------------------------
# Fixed series order/colors: LIMO3 (the run) in red; baselines by decreasing
# resolution class: MC blue, van Leer green, minmod orange, 1st-order gray.
SERIES = [
    ("mc",      "MC (baseline)",        "tab:blue",   "-",  1.0),
    ("vanleer", "van Leer (baseline)",  "tab:green",  "-",  1.0),
    ("minmod",  "minmod (baseline)",    "tab:orange", "-",  1.0),
    ("iord",    "1st order (baseline)", "0.5",        "-",  1.0),
]

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))

for ax in axs:
    for tag, label, color, ls, lw in SERIES:
        xb, rb = baselines[tag]
        ax.plot(xb, rb, ls, color=color, lw=lw, label=label)
    ax.plot(xc, rho, "-", color="tab:red", lw=1.8, label="Q2D LIMO3 (this run)")
    ax.set_xlabel("x")
    ax.set_ylabel("rho")
    ax.grid(True, alpha=0.3)

# left: full domain, TV window shaded
axs[0].axvspan(WIN_LO, WIN_HI, color="0.5", alpha=0.12,
               label="post-shock TV window")
axs[0].set_xlim(-5.0, 5.0)
axs[0].set_title("Shu-Osher, N=400, t=1.8")
axs[0].legend(fontsize=8, loc="lower left")

# right: zoom on the resolved entropy-wave (sine) packet behind the shock
axs[1].set_xlim(0.3, SHOCK_X + 0.15)
axs[1].set_ylim(2.6, 4.8)
axs[1].set_title("Resolved sine packet (zoom)")
tv_text = "post-shock TV (x in [%g, %g]):\n" % (WIN_LO, WIN_HI) + "\n".join(
    "%-9s %.3f" % (k, tv[k])
    for k in ("LIMO3", "mc", "vanleer", "minmod", "iord"))
axs[1].text(0.03, 0.03, tv_text, transform=axs[1].transAxes,
            fontsize=8, family="monospace", va="bottom", ha="left",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8))

plt.tight_layout()
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("Figure saved to %s" % OUT_SVG)
