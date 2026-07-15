"""borda-channel V&V figure -> docs/vv/images/Borda-state.svg.

Borda sudden-expansion injected state (first interior cell column adjacent to
the inflow face, profiles across the channel height y) against the Fievisohn
2016 analytic expectations (thesis Sec 3.3, eqs 3.26-3.31) evaluated at the
case's plenum/gas — the same anchors and tolerances as check.py:

  (A) injected mass flux rho*u  vs  the choked value mdot_A3 = rho1*a1*A1/A3
      (rel tol 3e-2, shaded);
  (B) recovered stagnation temperature T0 (isentropic + energy) vs the plenum
      T0 (rel tol 3e-2, shaded);
  (C) recovered stagnation pressure p0 vs the plenum P0: a real Borda-Carnot
      stagnation-pressure LOSS (p0 < P0).

Run './Q2D.sh solve' first, then 'python3 plot.py'.
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

HERE = Path(__file__).resolve().parent           # <repo>/test/2D/euler/borda-channel
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "Borda-state.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file, read_q2d_field

# Borda plenum + area ratio (identical to gen_case.py / check.py) and the
# check.py acceptance tolerance on anchors (A) and (B).
P0, T0, R_AREA = 1.0e6, 300.0, 0.2
TOL = 3.0e-2

x_, y_, var_, vnames = read_q2d_field(require_file(HERE / "OUTPUT" / "field.tec",
                                                   "run './Q2D.sh solve' first"))
yn = y_[0][:, :, 0]                              # nodal y, axes = (i, j)
cc = [n.strip() for n in vnames][3:]             # after nodal x, y, z
def f(name):
    return np.asarray(var_[0][cc.index(name)], float)[:, :, 0]

rho, u, p, T = f("rho(1)"), f("u"), f("p"), f("T")
g = float(np.median(f("g")))
R = float(np.median(f("R")))

# inlet column = first interior cells adjacent to the Borda face (i = 0)
yc = 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:])[0, :]
rho_in, u_in, p_in, T_in = rho[0, :], u[0, :], p[0, :], T[0, :]

mflux = rho_in * u_in                            # (A) conserved jump quantity
M_in = u_in / np.sqrt(g * R * T_in)
tf = 1.0 + 0.5 * (g - 1.0) * M_in**2
T0_inj = T_in * tf                               # (B) recovered T0
p0_inj = p_in * tf ** (g / (g - 1.0))            # (C) recovered p0

# analytic choked anchors at this gas (Fievisohn eq 3.29 + isentropic choke)
T1 = T0 * 2.0 / (g + 1.0)
P1 = P0 * (2.0 / (g + 1.0)) ** (g / (g - 1.0))
rho1 = P1 / (R * T1)
v1 = math.sqrt(g * R * T1)                       # a1 (M1 = 1)
mdot_A3 = rho1 * v1 * R_AREA                     # per unit A3

fig, axs = plt.subplots(1, 3, figsize=(12, 4.0), sharey=True)

def band(ax, level, tol, unit, label):
    ax.axvspan(level * (1 - tol), level * (1 + tol), color="tab:green",
               alpha=0.15, label="tolerance band (3%)")
    ax.axvline(level, color="k", ls="--", lw=1.2, label=label)

ax = axs[0]
band(ax, mdot_A3, TOL, "kg/m2/s", "choked mdot/A3 (Fievisohn eq 3.29)")
ax.plot(mflux, yc, "-", color="tab:red", lw=1.5, label="Q2D injected rho*u")
ax.set_xlabel("mass flux rho*u [kg/(m^2 s)]")
ax.set_ylabel("y [m]")
ax.set_title("(A) choked mass flux")
ax.set_xlim(mdot_A3 * (1 - 3 * TOL), mdot_A3 * (1 + 3 * TOL))
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="lower left")

ax = axs[1]
band(ax, T0, TOL, "K", "plenum T0 = %.0f K" % T0)
ax.plot(T0_inj, yc, "-", color="tab:red", lw=1.5, label="Q2D recovered T0")
ax.set_xlabel("recovered T0 [K]")
ax.set_title("(B) stagnation temperature")
ax.set_xlim(T0 * (1 - 3 * TOL), T0 * (1 + 3 * TOL))
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="lower left")

ax = axs[2]
ax.axvline(P0 / 1e6, color="k", ls="--", lw=1.2, label="plenum P0")
ax.plot(p0_inj / 1e6, yc, "-", color="tab:red", lw=1.5,
        label="Q2D recovered p0")
ax.set_xlabel("recovered p0 [MPa]")
ax.set_title("(C) Borda-Carnot p0 loss (p0 < P0)")
ax.set_xlim(0.6, 1.1)
ax.annotate("", xy=(np.mean(p0_inj) / 1e6, yc[len(yc) // 2]),
            xytext=(P0 / 1e6, yc[len(yc) // 2]),
            arrowprops=dict(arrowstyle="->", color="0.4"))
ax.text(0.5 * (np.mean(p0_inj) + P0) / 1e6, yc[len(yc) // 2], " loss",
        fontsize=8, color="0.4", va="bottom", ha="center")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="lower left")

fig.suptitle("Borda injector: Q2D injected inlet state vs Fievisohn 2016 "
             "(A1/A3 = %.1f, gamma = %.1f)" % (R_AREA, g))
plt.tight_layout(rect=(0, 0, 1, 0.94))
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("mflux %.2f vs %.2f | T0 %.2f vs %.0f | p0 %.0f < %.0f"
      % (mflux.mean(), mdot_A3, T0_inj.mean(), T0, p0_inj.mean(), P0))
print("Figure saved to %s" % OUT_SVG)
