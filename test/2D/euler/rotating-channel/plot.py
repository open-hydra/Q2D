"""rotating-channel V&V figure -> docs/vv/images/RC-dpdy.svg.

Left: x-averaged pressure profile p(y) of the converged rotating-channel run
(OUTPUT/field.tec), its linear fit, and the analytic geostrophic (Coriolis)
slope dp/dy = -2*omega*rho_0*U_inf from reference/dpdy_analytical.dat —
the same data pipeline as check.py.
Right: residual drift of the velocity field from the equilibrium IC
(u = U_inf, v = 0), per-y-column maxima, against the check.py drift gates.

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

HERE = Path(__file__).resolve().parent           # <repo>/test/2D/euler/rotating-channel
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "RC-dpdy.svg"

sys.path.insert(0, str(REPO_ROOT / "test" / "common" / "pycheck"))
from q2d_check_lib import require_file, read_q2d_field

# IC axial flow and drift gates: keep in sync with check.py.
U_INF = 100.0     # m/s (reference/dpdy_analytical.dat header)
THR_U = 1.5e-2    # check.py gate on max|u - U_inf| [m/s]
THR_V = 8.0e-2    # check.py gate on max|v| [m/s]

x_, y_, var_, vnames = read_q2d_field(require_file(HERE / "OUTPUT" / "field.tec",
                                                   "run './Q2D.sh solve' first"))
yn = y_[0][:, :, 0]
rho = var_[0][0][:, :, 0]
u = var_[0][1][:, :, 0]
v = var_[0][2][:, :, 0]
p = var_[0][3][:, :, 0]

# same reduction as check.py: mean across x (axis 0 = i), y at cell centers
yc = 0.5 * (yn[:-1, :-1] + yn[1:, 1:])
y_avg = yc[0, :]
p_avg = p.mean(axis=0)

slope, intercept = np.polyfit(y_avg, p_avg, 1)
ref = np.loadtxt(require_file(HERE / "reference" / "dpdy_analytical.dat"),
                 comments="#")
dpdy_ref = float(ref if ref.ndim == 0 else ref[0])
eps_rel = abs(slope - dpdy_ref) / abs(dpdy_ref)

# analytic profile with the geostrophic slope, anchored at the fitted mean
p_fit = slope * y_avg + intercept
p_ana = dpdy_ref * (y_avg - y_avg.mean()) + p_avg.mean()

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))

ax = axs[0]
ax.plot(y_avg, p_ana / 1e3, "k--", lw=1.5,
        label="analytic slope -2*omega*rho_0*U_inf\n(%.1f Pa/m)" % dpdy_ref)
ax.plot(y_avg, p_avg / 1e3, "o", color="tab:red", ms=4, mfc="none",
        label="Q2D p(y), x-averaged")
ax.plot(y_avg, p_fit / 1e3, "-", color="tab:blue", lw=1.2,
        label="linear fit (%.1f Pa/m, rel err %.1e)" % (slope, eps_rel))
ax.set_xlabel("y [m]")
ax.set_ylabel("p [kPa]")
ax.set_title("Coriolis equilibrium: pressure gradient")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, loc="upper right")

ax = axs[1]
u_drift = np.abs(u - U_INF).max(axis=0)          # per-y-column max over x
v_drift = np.abs(v).max(axis=0)
ax.semilogy(y_avg, u_drift, "-", color="tab:blue", lw=1.2,
            label="max_x |u - U_inf|")
ax.axhline(THR_U, color="tab:blue", ls="--", lw=1,
           label="gate on |u - U_inf| (check.py)")
ax.semilogy(y_avg, v_drift, "-", color="tab:orange", lw=1.2,
            label="max_x |v|")
ax.axhline(THR_V, color="tab:orange", ls="--", lw=1,
           label="gate on |v| (check.py)")
ax.set_xlabel("y [m]")
ax.set_ylabel("velocity drift [m/s]")
ax.set_title("Residual drift from the equilibrium IC")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=8, loc="lower right")

fig.suptitle("Rotating channel (omega = 10 rad/s, U_inf = %.0f m/s): "
             "Q2D vs geostrophic balance" % U_INF)
plt.tight_layout(rect=(0, 0, 1, 0.94))
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(str(OUT_SVG), bbox_inches="tight", transparent=True)
print("slope %.4e vs analytic %.4e (rel err %.2e)" % (slope, dpdy_ref, eps_rel))
print("Figure saved to %s" % OUT_SVG)
