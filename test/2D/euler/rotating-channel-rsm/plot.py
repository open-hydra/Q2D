"""rotating-channel-rsm V&V figure -> docs/vv/images/RC-RSM-frame.svg.

Left: homogeneous-ODE trajectories of the Reynolds stresses R11, R22, R33, R12
(RK4, the same reference integrated by verify.py) with the converged Q2D
mid-channel values at t = T overlaid as markers; the no-G_ij trajectory of R12
is dashed — the visual form of the verify.py discriminance gate G3 (a binary
without the exact Coriolis production lands on the dashed curve, three orders
of magnitude away from the solid one).
Right: y-profile of R12(T) across the channel, Q2D vs the per-cell ODE
prediction (each cell has its own IC because rho(y) follows the Coriolis
geostrophic equilibrium — verify.py gate G1); the symmetry walls pin R12 = 0.

Run './Q2D.sh solve' first, then 'python3 plot.py'.
Constants and readers: keep in sync with verify.py.
"""
import re
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

HERE = Path(__file__).resolve().parent           # <repo>/test/2D/euler/rotating-channel-rsm
REPO_ROOT = HERE.parents[3]
OUT_SVG = REPO_ROOT / "docs" / "vv" / "images" / "RC-RSM-frame.svg"

# ---- model constants and ODE: keep in sync with verify.py -------------------
OMEGA, T = 10.0, 0.05
Cmu = 0.09
C1e, C2e, C5e, be = 1.7, 1.05, 0.2, 0.0828

def rhs(y, withG):
    R = np.array([[y[0], y[3], 0], [y[3], y[1], 0], [0, 0, y[2]]]); om = y[4]
    k = 0.5*np.trace(R); eps = Cmu*k*om
    a = R/k - 2.0/3.0*np.eye(3)
    aa = np.sum(a*a)
    W = np.array([[0, -OMEGA, 0], [OMEGA, 0, 0], [0, 0, 0.0]])
    Pi = (-(C1e*eps)*a + C2e*eps*(a@a - aa/3*np.eye(3))
          + C5e*k*(a@W.T + (a@W.T).T))
    dR = Pi - 2.0/3.0*eps*np.eye(3)
    if withG:
        G = np.zeros((3, 3))
        G[0, 0] = 4*OMEGA*R[0, 1]; G[1, 1] = -4*OMEGA*R[0, 1]
        G[0, 1] = G[1, 0] = 2*OMEGA*(R[1, 1]-R[0, 0])
        dR = dR + G
    return np.array([dR[0, 0], dR[1, 1], dR[2, 2], dR[0, 1], -be*om*om])

def rk4_traj(y0, withG, n=2000):
    y, dt = np.array(y0, float), T/n
    ts, ys = [0.0], [y.copy()]
    for _ in range(n):
        k1 = rhs(y, withG); k2 = rhs(y+dt/2*k1, withG)
        k3 = rhs(y+dt/2*k2, withG); k4 = rhs(y+dt*k3, withG)
        y = y + dt/6*(k1+2*k2+2*k3+k4)
        ts.append(ts[-1]+dt); ys.append(y.copy())
    return np.array(ts), np.array(ys)

def read_field(path):
    txt = open(path).read().split('\n')
    m = re.search(r'I=\s*(\d+), J=\s*(\d+)', txt[1]); I, J = int(m.group(1)), int(m.group(2))
    nn, nc = I*J, (I-1)*(J-1)
    vals = [float(t) for ln in txt[2:] for t in ln.split()]
    nvar = (len(vals)-2*nn)//nc
    nodes = np.array(vals[:2*nn]).reshape(2, I, J, order='C')
    cells = np.array(vals[2*nn:2*nn+nvar*nc]).reshape(nvar, I-1, J-1, order='C')
    return nodes, cells

# ---- data --------------------------------------------------------------------
_, c = read_field(HERE / "OUTPUT" / "field.tec")
nodes, ic = read_field(HERE / "INPUT" / "ic.tec")
rho, rho0 = c[0], ic[0]
Rcode = {k: c[i]/rho for k, i in (('11', 4), ('22', 5), ('33', 6), ('12', 7))}
ni, nj = rho.shape
mid = (ni//2, nj//2)
yc = 0.25*(nodes[1, :-1, :-1] + nodes[1, 1:, :-1] + nodes[1, :-1, 1:] + nodes[1, 1:, 1:])

y0 = [ic[v][mid]/rho0[mid] for v in (4, 5, 6, 7, 8)]
ts, yF = rk4_traj(y0, True)
_, yN = rk4_traj(y0, False)

# ---- figure ------------------------------------------------------------------
fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.5, 4.2))

labels = [r"$R_{11}$", r"$R_{22}$", r"$R_{33}$", r"$R_{12}$"]
colors = ["#1f77b4", "#2ca02c", "#9467bd", "#d62728"]
for q, (lab, col) in enumerate(zip(labels, colors)):
    axL.plot(ts, yF[:, q], color=col, lw=1.6, label=lab + " ODE (with $G_{ij}$)")
    axL.plot(T, Rcode[('11', '22', '33', '12')[q]][mid], "o", color=col, ms=7,
             mfc="none", mew=1.8, label=lab + " Q2D" if q == 3 else None)
axL.plot(ts, yN[:, 3], "--", color="#d62728", lw=1.4,
         label=r"$R_{12}$ ODE without $G_{ij}$")
axL.set_xlabel(r"$t$ [s]")
axL.set_ylabel(r"$R_{ij}$ [m$^2$/s$^2$]")
axL.set_title("Homogeneous evolution, mid-channel cell")
axL.legend(fontsize=7.5, ncol=2, frameon=False)
axL.grid(alpha=0.3)

# right: R12(T) profile, per-cell ODE (rho(y) geostrophic) vs Q2D
i0 = mid[0]
r12_ode = []
for jj in range(nj):
    r0 = rho0[i0, jj]
    y0c = [ic[v][i0, jj]/r0 for v in (4, 5, 6, 7, 8)]
    r12_ode.append(rk4_traj(y0c, True, n=800)[1][-1, 3])
axR.plot(np.array(r12_ode), yc[i0, :], "-", color="#1f77b4", lw=1.6,
         label="per-cell ODE (with $G_{ij}$)")
axR.plot(Rcode['12'][i0, :], yc[i0, :], "o", ms=4, mfc="none", mew=1.4,
         color="#d62728", label="Q2D")
axR.axvline(yN[-1, 3], ls="--", color="#7f7f7f", lw=1.2,
            label="no-$G_{ij}$ prediction")
axR.set_xlabel(r"$R_{12}(T)$ [m$^2$/s$^2$]")
axR.set_ylabel(r"$y$ [m]")
axR.set_title(r"$R_{12}$ profile at $t = T$ (walls pin $R_{12}=0$)")
axR.legend(fontsize=8, frameon=False)
axR.grid(alpha=0.3)

fig.tight_layout()
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_SVG)
print(f"figure: {OUT_SVG}")
