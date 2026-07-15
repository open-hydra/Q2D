"""planar-nozzle V&V figure -> docs/vv/images/PN-Mach.svg.

Left: y-averaged Mach(x) of the converged Q2D area-variation nozzle
(OUTPUT/field.tec) against the quasi-1D isentropic area-law analytic curve
(Shapiro), with the multigrid sibling's coarse level overlaid.
Right: planar-nozzle-mg2 multigrid consistency — relative difference between
the coarse-level Mach (field-level2.tec) and the 2:1 restriction of the fine
Mach, against the check.py acceptance threshold.

Run './Q2D.sh solve' here and in ../planar-nozzle-mg2 first, then
'python3 plot.py'.
"""
import math
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

HERE = Path(__file__).resolve().parent           # <repo>/test/2D/euler/planar-nozzle
MG2 = HERE.parent / "planar-nozzle-mg2"
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "PN-Mach.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file, read_q2d_field

# --- nozzle area law (identical to gen_case.py / check.py) -------------------
L, XT, AT, AIN, AEX = 1.0, 0.40, 1.0, 2.0, 2.5
def area(x):
    return (AT + (AIN-AT)*0.5*(1+math.cos(math.pi*x/XT)) if x <= XT
            else AT + (AEX-AT)*0.5*(1-math.cos(math.pi*(x-XT)/(L-XT))))
def AR_of_M(M, g): return (1.0/M)*((2.0/(g+1))*(1+0.5*(g-1)*M*M))**((g+1)/(2*(g-1)))
def M_of_AR(ar, g, sup):
    if abs(ar-1) < 1e-13: return 1.0
    lo, hi = (1+1e-12, 80.0) if sup else (1e-9, 1-1e-12)
    for _ in range(300):
        m = 0.5*(lo+hi); f = AR_of_M(m, g)-ar
        if sup: hi, lo = (m, lo) if f > 0 else (hi, m)
        else:   lo, hi = (m, hi) if f > 0 else (lo, m)
    return 0.5*(lo+hi)

# check.py layer-3 acceptance on the coarse/fine Mach restriction (mg2)
CF_TOL = 5.0e-3


def load_mach(path):
    """(x cell centers along i, y-averaged Mach(x)) from a Q2D field file."""
    x_, y_, var_, vnames = read_q2d_field(require_file(path))
    xn = x_[0][:, :, 0]                          # nodal x, axes = (i, j)
    cc = [n.strip() for n in vnames][3:]         # after nodal x, y, z
    def f(name):
        return np.asarray(var_[0][cc.index(name)], float)[:, :, 0]
    u, v = f("u"), f("v")
    T, g, R = f("T"), f("g"), f("R")
    M = np.sqrt(u**2 + v**2) / np.sqrt(g * R * T)
    xc = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
    return xc[:, 0], M.mean(axis=1), float(np.median(g))


x_f, M_f, g = load_mach(HERE / "OUTPUT" / "field.tec")
x_mg, M_mg, _ = load_mach(MG2 / "OUTPUT" / "field.tec")
x_c, M_c, _ = load_mach(MG2 / "OUTPUT" / "field-level2.tec")

M_a = np.array([M_of_AR(area(x) / AT, g, x > XT) for x in x_f])

# 2:1 restriction of the mg2 fine Mach onto the coarse cells (as in
# planar-nozzle-mg2 check.py layer 3)
M_restr = M_mg.reshape(len(x_c), 2).mean(axis=1)
cf_rel = np.abs(M_c - M_restr) / np.maximum(np.abs(M_restr), 1e-30)

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))

ax = axs[0]
ax.plot(x_f, M_a, "k-", lw=1.5, label="quasi-1D area law (analytic)")
ax.plot(x_f, M_f, "-", color="tab:red", lw=1.2, label="Q2D (single grid)")
ax.plot(x_c[::4], M_c[::4], "o", color="tab:blue", ms=4, mfc="none",
        label="Q2D MG 2-level, coarse level")
ax.axvline(XT, color="0.5", ls=":", lw=1)
ax.text(XT, ax.get_ylim()[0] + 0.05, " throat", fontsize=8, color="0.4")
ax.set_xlabel("x [m]")
ax.set_ylabel("Mach (y-averaged)")
ax.set_title("C-D nozzle: Mach vs isentropic area law")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper left")

ax = axs[1]
ax.semilogy(x_c, cf_rel, "o-", color="tab:blue", ms=4, lw=1,
            label="|M_coarse / M_restricted - 1|")
ax.axhline(CF_TOL, color="k", ls="--", lw=1,
           label="acceptance (check.py, 5e-3)")
ax.set_xlabel("x [m]")
ax.set_ylabel("coarse/fine Mach relative difference")
ax.set_title("MG consistency: coarse level vs 2:1 restricted fine")
ax.set_ylim(top=3e-2)                            # headroom above the tolerance
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=8, loc="upper right")

plt.tight_layout()
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("max coarse/fine rel diff = %.3e (tol %.0e)" % (cf_rel.max(), CF_TOL))
print("Figure saved to %s" % OUT_SVG)
