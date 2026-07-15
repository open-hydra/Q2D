"""flat-plate-turbulent/coarse-kw V&V figure -> docs/vv/images/FPT-Cf.svg.

Physical skin-friction coefficient Cf(x) = tauX / (0.5 rho_inf U_inf^2) of
the converged Q2D SST run (OUTPUT/wall.tec) against the committed NASA
Turbulence Modeling Resource SST reference curve (reference/cf-sst.dat,
Tecplot zone T="FUN3D, triangles"), with the check.py primary-gate +/-15%
band over the gated plate region x in [0.05, 1.0] m. The Wilcox2006 curve
(reference/cf-wilcox.dat, two zones: CFL3D and FUN3D) is overlaid as the
informative cross-check only (different model from the SST run — see
check.py); its zones are split where the reference x resets so the zone
boundary is not drawn as a spurious segment.

Everything numerical is imported from check.py (wall.tec parser, reference
reader, freestream state rho_inf/U_inf/q_inf, gate function and tolerance):
this script plots exactly what the ctest gates, it does not re-derive it.

Run './Q2D.sh -p 20 solve' here first (heavy tier), then 'python3 plot.py'.
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

HERE = Path(__file__).resolve().parent  # <repo>/test/2D/viscous/flat-plate-turbulent/coarse-kw
REPO_ROOT = HERE.parents[4]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "FPT-Cf.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file

# Reuse the ctest's own parsing + gate machinery (no duplicated numbers).
sys.path.insert(0, str(HERE))
import check  # noqa: E402  (freestream state + parsers + gates)

# check.py primary physical-Cf gate: tol=0.15 over x in [0.05, 1.0]
GATE_TOL, GATE_XLO, GATE_XHI = 0.15, 0.05, 1.0

wall_path = HERE / "OUTPUT" / "wall.tec"
ref_sst = HERE / "reference" / "cf-sst.dat"
ref_wilcox = HERE / "reference" / "cf-wilcox.dat"

_, x_nodal, tauX_cc = check._parse_tecplot_wall(str(require_file(str(wall_path))))
sst_x, sst_cf = check._read_cf_reference(str(require_file(str(ref_sst))))
wx_x, wx_cf = check._read_cf_reference(str(require_file(str(ref_wilcox))))

# Q2D physical Cf at wall-cell centers (identical to check.cf_physical_compare)
x_nodal = np.asarray(x_nodal, float)
xc = 0.5 * (x_nodal[:-1] + x_nodal[1:])
cf_q2d = np.asarray(tauX_cc, float) / check.Q_INF

# Reproduce the PRIMARY gate number for the figure annotation.
pass_phys, status_phys = check.cf_physical_compare(
    list(x_nodal), list(tauX_cc), sst_x, sst_cf, tol=GATE_TOL,
    label="physical Cf vs cf-sst.dat", x_lo=GATE_XLO, x_hi=GATE_XHI)
print(status_phys)
if not pass_phys:
    raise SystemExit("gate FAILED on committed OUTPUT — refusing to plot a "
                     "passing-looking figure: " + status_phys)

sst_x = np.asarray(sst_x, float)
sst_cf = np.asarray(sst_cf, float)
band = (sst_x >= GATE_XLO) & (sst_x <= GATE_XHI)

fig, ax = plt.subplots(figsize=(6.4, 4.2))

ax.fill_between(sst_x[band], (1 - GATE_TOL) * sst_cf[band],
                (1 + GATE_TOL) * sst_cf[band], color="0.6", alpha=0.25, lw=0,
                label=u"±15% gate band (check.py, x ∈ [0.05, 1.0])")
ax.plot(sst_x, sst_cf, "k-", lw=1.5,
        label="NASA TMR FUN3D SST reference (cf-sst.dat)")
# cf-wilcox.dat concatenates a CFL3D and a FUN3D zone: split where x resets
# so the zone boundary is not drawn as a spurious segment.
wx_x = np.asarray(wx_x, float)
wx_cf = np.asarray(wx_cf, float)
zone_starts = [0] + [i + 1 for i in np.nonzero(np.diff(wx_x) < 0)[0]] + [len(wx_x)]
for k in range(len(zone_starts) - 1):
    s = slice(zone_starts[k], zone_starts[k + 1])
    ax.plot(wx_x[s], wx_cf[s], ":", color="tab:blue", lw=1.2,
            label="NASA TMR Wilcox2006, CFL3D+FUN3D zones (informative)"
            if k == 0 else "_nolegend_")
on_plate = xc >= 0.0
ax.plot(xc[on_plate], cf_q2d[on_plate], "-", color="tab:red", lw=1.2,
        label="Q2D SST (coarse-kw, wall.tec)")
for xg in (GATE_XLO, GATE_XHI):
    ax.axvline(xg, color="0.5", ls=":", lw=1)

ax.set_xlim(0.0, 2.0)
ax.set_ylim(0.0, 8.0e-3)
ax.set_xlabel("x [m]")
ax.set_ylabel(r"$C_f = \tau_x / (\frac{1}{2}\rho_\infty U_\infty^2)$")
ax.set_title("Turbulent flat plate: physical $C_f(x)$ vs NASA TMR (SST)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper right")
ax.annotate(status_phys.split(":", 1)[1].split("cells_compared")[0].strip()
            .replace("max_rel_diff=", "max rel diff (gated region) = "),
            xy=(0.03, 0.05), xycoords="axes fraction", fontsize=8, color="0.35")

plt.tight_layout()
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("Figure saved to %s" % OUT_SVG)
