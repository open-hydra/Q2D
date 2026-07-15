"""prandtl-meyer-lineprobe V&V figure -> docs/vv/images/PM-lineprobe.svg.

The Q2D-only line-probe (OUTPUT/midrow.tec, last zone, cell row j = 25)
overlaid on the same row extracted from the main field writer
(OUTPUT/field.tec) — the check.py self-consistency gate — plus the analytic
post-expansion Prandtl-Meyer Mach level (M1 = 2, theta = 34.971 deg, same
functions as ../prandtl-meyer/check.py). Left: Mach; right: pressure.

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

HERE = Path(__file__).resolve().parent  # <repo>/test/2D/euler/prandtl-meyer-lineprobe
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "PM-lineprobe.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import fail, require_file

# Line-probe layout constants: keep in sync with check.py (input.ini [midrow]).
FIXED_J = 25
LP_VARS = ["rho", "u", "v", "p", "T", "M"]
FIELD_VARS = ["rho", "u", "v", "p", "T", "g", "R"]

# Analytic Prandtl-Meyer expansion (same setup and functions as the sibling
# ../prandtl-meyer/check.py: inflow M1 = 2 from bc.txt BC 405, wall deflection
# theta = 34.971 deg; INPUT/ic.tec and bc.txt are identical in both cases).
GAMMA = 1.4
M1, THETA_DEG = 2.0, 34.971


def nufunc(M, gamma):
    a = np.sqrt((gamma + 1.0) / (gamma - 1.0))
    return (a * np.arctan(np.sqrt(M**2 - 1.0) / a) -
            np.arctan(np.sqrt(M**2 - 1.0)))


def mach_from_nu(nu, gamma, a=1.01, b=30.0, tol=1e-6, max_iter=1000000):
    for _ in range(max_iter):
        M_mid = 0.5 * (a + b)
        nu_err = nufunc(M_mid, gamma) - nu
        if abs(nu_err) < tol:
            return M_mid
        if np.sign(nu_err) < 0:
            a = M_mid
        else:
            b = M_mid
    return 0.5 * (a + b)


M2_analytic = mach_from_nu(nufunc(M1, GAMMA) + THETA_DEG * np.pi / 180.0, GAMMA)


def read_last_zone(path):
    """(ni, nj, nk, flat token list) of the LAST Tecplot zone (as check.py)."""
    lines = Path(path).read_text().splitlines()
    zone_idx = [i for i, l in enumerate(lines) if "ZONE" in l]
    if not zone_idx:
        fail("no ZONE record in %s" % path)
    z = zone_idx[-1]
    m = re.search(r"I=\s*(\d+).*J=\s*(\d+).*K=\s*(\d+)", lines[z])
    ni, nj, nk = int(m.group(1)), int(m.group(2)), int(m.group(3))
    toks = [float(t) for ln in lines[z + 1:] for t in ln.split()]
    return ni, nj, nk, toks


# --- line probe (last zone): strip nodes + cell values along the line --------
lp_ni, lp_nj, _, lp_toks = read_last_zone(
    require_file(HERE / "OUTPUT" / "midrow.tec", "run './Q2D.sh solve' first"))
Ni = lp_ni - 1
nxy = 2 * (Ni + 1)
lp_x, lp_y = lp_toks[:nxy], lp_toks[nxy:2 * nxy]
lp = {}
o = 2 * nxy
for name in LP_VARS:
    lp[name] = np.array(lp_toks[o:o + Ni])
    o += Ni

# cell-center x along the sampled line (mean of the 4 strip nodes per cell)
xlo, xup = np.array(lp_x[:Ni + 1]), np.array(lp_x[Ni + 1:])
x_line = 0.25 * (xlo[:-1] + xlo[1:] + xup[:-1] + xup[1:])

# --- same row from the main field writer (as check.py self-consistency) ------
f_ni, f_nj, _, f_toks = read_last_zone(
    require_file(HERE / "OUTPUT" / "field.tec", "run './Q2D.sh solve' first"))
nnode, ncell = f_ni * f_nj, (f_ni - 1) * (f_nj - 1)
o = 2 * nnode
fvars = {}
for name in FIELD_VARS:
    fvars[name] = np.array(f_toks[o:o + ncell]).reshape(f_nj - 1, f_ni - 1)
    o += ncell
row = {name: fvars[name][FIXED_J - 1, :] for name in FIELD_VARS}
g_row, R_row = row["g"], row["R"]
M_field = (np.sqrt(row["u"]**2 + row["v"]**2) /
           np.sqrt(g_row * R_row * row["T"]))

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))

ax = axs[0]
ax.axhline(M2_analytic, color="k", ls="--", lw=1.2,
           label="analytic post-expansion M = %.3f" % M2_analytic)
ax.plot(x_line, M_field, "-", color="tab:red", lw=1.5,
        label="Q2D field row j = %d" % FIXED_J)
ax.plot(x_line[::3], lp["M"][::3], "o", color="tab:blue", ms=4, mfc="none",
        label="Q2D line probe (midrow.tec)")
ax.set_xlabel("x")
ax.set_ylabel("Mach")
ax.set_title("Mach along the sampled row (max |lp-field| = %.1e)"
             % float(np.abs(lp["M"] - M_field).max()))
ax.set_ylim(top=M2_analytic + 0.2)               # headroom above the analytic level
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="lower right")

ax = axs[1]
ax.plot(x_line, row["p"], "-", color="tab:red", lw=1.5,
        label="Q2D field row j = %d" % FIXED_J)
ax.plot(x_line[::3], lp["p"][::3], "o", color="tab:blue", ms=4, mfc="none",
        label="Q2D line probe (midrow.tec)")
ax.set_xlabel("x")
ax.set_ylabel("p")
ax.set_title("Pressure along the sampled row (max |lp-field| = %.1e)"
             % float(np.abs(lp["p"] - row["p"]).max()))
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper right")

fig.suptitle("Prandtl-Meyer expansion: line-probe self-consistency "
             "(M1 = %.0f, theta = %.3f deg)" % (M1, THETA_DEG))
plt.tight_layout(rect=(0, 0, 1, 0.94))
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("M: line %.4f -> %.4f (analytic post-expansion %.4f)"
      % (lp["M"][0], lp["M"][-1], M2_analytic))
print("Figure saved to %s" % OUT_SVG)
