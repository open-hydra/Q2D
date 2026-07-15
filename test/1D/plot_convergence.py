"""1D refinement-ladder convergence V&V figure -> docs/vv/images/convergence-1D.svg.

Log-log L1(rho) error vs cell size h for the four 1D Euler ladders
(Sod79, Ein91, Tor99, Noh87), each with its least-squares observed order.
Exactly the same data pipeline as the per-case ``verify.py --gate`` scripts:
OUTPUT/field_x{1,2,4,8}.tec cell centers vs the analytic exact solutions in
test/common/pycheck/q2d_exact1d.py, relative L1 via q2d_check_lib.rel_l1.

Run './Q2D.sh test' in each of the four case directories first, then
'python3 plot_convergence.py'.
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

HERE = Path(__file__).resolve().parent           # <repo>/test/1D
REPO_ROOT = HERE.parents[1]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "convergence-1D.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file, read_q2d_field, rel_l1
from q2d_exact1d import riemann_exact, noh_exact_planar, observed_order

LEVELS = (1, 2, 4, 8)

# Case setups: keep in sync with each <case>/verify.py (which matches the
# case's input.ini / reference/exact.txt).
RIEMANN_CASES = {
    #        x0 , t_final, rho_l, u_l , p_l, rho_r, u_r, p_r
    "Sod79": (0.5, 0.20, 1.0, 0.00, 1.0, 0.125, 0.0, 0.1),
    "Ein91": (0.5, 0.15, 1.0, -2.00, 0.4, 1.000, 2.0, 0.4),
    "Tor99": (0.3, 0.20, 1.0, 0.75, 1.0, 0.125, 0.0, 0.1),
}
NOH_CASE = ("Noh87", 0.5, 1.0, 1.0, 1.0, 1e-15, 1.4)  # x0, t, rho0, u0, p0, g

# Fixed series order/colors/markers (one per case).
STYLE = {
    "Sod79": ("tab:blue",   "o"),
    "Ein91": ("tab:orange", "s"),
    "Tor99": ("tab:green",  "^"),
    "Noh87": ("tab:red",    "D"),
}


def load_level(case, level):
    """Cell centers + rho of one ladder level (same as verify.py _load_level)."""
    path = HERE / case / "OUTPUT" / ("field_x%d.tec" % level)
    require_file(path, "run './Q2D.sh test' in test/1D/%s first" % case)
    x_, y_, var_, vnames = read_q2d_field(path)
    xc = 0.5 * (x_[0][:-1] + x_[0][1:])[:, 0, 0]
    cc = [n.strip() for n in vnames][3:]         # after nodal x, y, z
    rho = np.asarray(var_[0][cc.index("rho(1)")], float)[:, 0, 0]
    return xc, rho


def ladder(case):
    """(h, rel L1(rho) error) arrays over the x1..x8 ladder of one case."""
    hs, errs = [], []
    for level in LEVELS:
        xc, rho = load_level(case, level)
        if case == "Noh87":
            _, x0, t, rho0, u0, p0, g = NOH_CASE
            ex = noh_exact_planar(xc, t, x0=x0, rho0=rho0, u0=u0, p0=p0, gamma=g)
        else:
            x0, t, rl, ul, pl, rr, ur, pr = RIEMANN_CASES[case]
            ex = riemann_exact(xc, t, x0, rl, ul, pl, rr, ur, pr)
        errs.append(rel_l1(rho, ex["rho"]))
        hs.append((xc[-1] - xc[0]) / (len(xc) - 1))
    return np.array(hs), np.array(errs)


fig, ax = plt.subplots(1, 1, figsize=(6.4, 5.0))

h_all = []
for case in ("Sod79", "Ein91", "Tor99", "Noh87"):
    hs, errs = ladder(case)
    p_obs = observed_order(hs, errs)
    color, marker = STYLE[case]
    ax.loglog(hs, errs, marker=marker, ls="-", color=color, ms=6,
              label="%s  (observed order %.2f)" % (case, p_obs))
    h_all.append(hs)
    print("%s: p_obs = %.3f, L1(rho) = %s" % (case, p_obs, errs))

# order-1 reference slope anchored at the coarsest Sod79 point
h_ref = np.array([h_all[0].max(), h_all[0].min()])
e_anchor = 2.0e-2
ax.loglog(h_ref, e_anchor * h_ref / h_ref[0], "k--", lw=1,
          label="order 1 (reference slope)")

ax.set_xlabel("h (cell size)")
ax.set_ylabel("relative L1(rho) error")
ax.set_title("Q2D 1D Euler ladders: grid convergence (x1 - x8)")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=9, loc="lower right")

plt.tight_layout()
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("Figure saved to %s" % OUT_SVG)
