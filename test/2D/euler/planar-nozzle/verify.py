# -*- coding: utf-8 -*-
#
# verify.py — planar-nozzle (q2D area-variation): publication-grade grid-convergence
# verification against the analytic quasi-1D area-law (Shapiro) solution, constant-gamma
# air. Produces two Journal-of-Propulsion-and-Power-style figures:
#
#   nozzle_solution.{png,svg}    : nozzle geometry, Mach & p/p0 vs x/L (CFD vs area law),
#                                  and the local relative error (the throat spike).
#   nozzle_convergence.{png,svg} : L1/L2/Linf Mach-error vs dx (log-log) with observed
#                                  order, slope-1/2 reference triangles, smooth-region
#                                  order, and Richardson-extrapolation + GCI.
#
# The throat is a weak singularity: A/A*-1 ~ (x-x_t)^2 and M-1 ~ sqrt(A/A*-1) => M-1 ~
# |x-x_t| (a corner). The kink signature of a 2nd-order scheme is Linf~O(dx), L2~O(dx^1.5),
# L1~O(dx^2); the smooth region and Richardson recover full 2nd order. gamma and R are
# read FROM the solver field so the comparison isolates discretisation error.
#
# Usage:  python3 verify.py            (runs the convergence sweep in ./.conv and plots)

import os, sys, math, shutil, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

HERE   = os.path.dirname(os.path.abspath(__file__))
GEN    = os.path.join(HERE, "gen_case.py")
ROOT   = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
SOLVER = os.path.join(ROOT, "bin", "Q2D")
SCRATCH = os.path.join(HERE, ".conv")
GRIDS  = [50, 100, 200, 400, 800]

# ---- nozzle area law (identical to gen_case.py) ----
L, XT, AT, AIN, AEX = 1.0, 0.40, 1.0, 2.0, 2.5
def area(x):
    x = np.asarray(x, float); A = np.empty_like(x); c = x <= XT
    A[c] = AT + (AIN-AT)*0.5*(1+np.cos(np.pi*x[c]/XT))
    dz = (x[~c]-XT)/(L-XT); A[~c] = AT + (AEX-AT)*0.5*(1-np.cos(np.pi*dz))
    return A
def AR_of_M(M, g): return (1.0/M)*((2.0/(g+1))*(1+0.5*(g-1)*M*M))**((g+1)/(2*(g-1)))
def M_of_AR(ar, g, sup):
    if abs(ar-1) < 1e-13: return 1.0
    lo, hi = (1+1e-12, 80.0) if sup else (1e-9, 1-1e-12)
    for _ in range(300):
        m = 0.5*(lo+hi); f = AR_of_M(m, g)-ar
        if sup: hi, lo = (m, lo) if f > 0 else (hi, m)
        else:   lo, hi = (m, hi) if f > 0 else (lo, m)
    return 0.5*(lo+hi)
def mach_curve(xq, g):
    ar = area(np.atleast_1d(xq)*L)/AT
    return np.array([M_of_AR(a, g, x > XT/L) for a, x in zip(ar, np.atleast_1d(xq))])

def parse_field(path):
    Ls = open(path).readlines()
    zi = next(i for i, l in enumerate(Ls) if l.strip().startswith("ZONE"))
    z = Ls[zi].replace(",", " ").split()
    NIp = int(next(t for t in z if t.startswith("I="))[2:])
    NJp = int(next(t for t in z if t.startswith("J="))[2:])
    NI, NJ = NIp-1, NJp-1
    v = []
    for l in Ls[zi+1:]:
        for t in l.replace(",", " ").split():
            try: v.append(float(t))
            except ValueError: pass
    v = np.array(v); off = [0]
    def take(n): a = v[off[0]:off[0]+n]; off[0] += n; return a
    nn, nc = NIp*NJp, NI*NJ
    xn = take(nn).reshape(NJp, NIp); take(nn)
    rho, u = take(nc).reshape(NJ, NI), take(nc).reshape(NJ, NI)
    take(nc); p = take(nc).reshape(NJ, NI)
    T, g, R = (take(nc).reshape(NJ, NI) for _ in range(3))
    m = lambda a: a.mean(axis=0)
    return dict(NI=NI, xc=0.5*(xn[0, :-1]+xn[0, 1:]),
                rho=m(rho), u=m(u), p=m(p), T=m(T), g=float(np.median(g)), R=float(np.median(R)))

def run_grid(NI):
    subprocess.run([sys.executable, GEN, str(NI), "4", SCRATCH], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    shutil.copy(SOLVER, os.path.join(SCRATCH, "bin", "Q2D"))
    for f in os.listdir(os.path.join(SCRATCH, "OUTPUT")):
        os.remove(os.path.join(SCRATCH, "OUTPUT", f))
    subprocess.run(["./bin/Q2D"], cwd=SCRATCH, env=dict(os.environ, OMP_NUM_THREADS="1"),
                   stdout=open(os.path.join(SCRATCH, "run.log"), "w"), stderr=subprocess.STDOUT)
    return parse_field(os.path.join(SCRATCH, "OUTPUT", "field.tec"))

# ======================= run sweep =======================
os.makedirs(SCRATCH, exist_ok=True)
print("# planar-nozzle grid-convergence verification (q2D area-variation vs quasi-1D area law)")
data = []
for NI in GRIDS:
    d = run_grid(NI); g0 = d["g"]
    xc = d["xc"]; M_q = d["u"]/np.sqrt(d["g"]*d["R"]*d["T"]); M_a = mach_curve(xc/L, g0)
    e = np.abs(M_q-M_a)
    interior = (xc > 2*L/NI) & (xc < L-2*L/NI)
    smooth = interior & (np.abs(xc-XT) > 0.05)
    rec = dict(NI=NI, dx=L/NI, g=g0, R=d["R"], xc=xc, M_q=M_q, M_a=M_a, p=d["p"], ar=area(xc)/AT,
               relerr=e/M_a,
               L1=float(np.mean(e[interior])), L2=float(np.sqrt(np.mean(e[interior]**2))),
               Linf=float(np.max(e[interior])), xLinf=float(xc[interior][np.argmax(e[interior])]),
               Linf_s=float(np.max(e[smooth])), L2_s=float(np.sqrt(np.mean(e[smooth]**2))),
               M_exit=float(M_q[-1]))
    data.append(rec)
    print(f"  N={NI:4d} dx={rec['dx']:.5f} g={g0:.5f} R={rec['R']:.3f}  "
          f"L1={rec['L1']:.2e} L2={rec['L2']:.2e} Linf={rec['Linf']:.2e}@x={rec['xLinf']:.3f}  "
          f"Linf_s={rec['Linf_s']:.2e} M_exit={rec['M_exit']:.5f}")

dxs = np.array([r["dx"] for r in data])
fitp = lambda key: float(np.polyfit(np.log(dxs), np.log([r[key] for r in data]), 1)[0])
orders = {k: fitp(k) for k in ("L1", "L2", "Linf", "Linf_s", "L2_s")}
print("# observed order (log-log fit): " + "  ".join(f"{k}={v:.2f}" for k, v in orders.items()))

# Richardson + GCI on exit Mach (3 finest, r=2)
f3, f2, f1 = data[-3]["M_exit"], data[-2]["M_exit"], data[-1]["M_exit"]
r = 2.0
p_RE = math.log(abs((f3-f2)/(f2-f1)))/math.log(r)
f_ext = f1 + (f1-f2)/(r**p_RE-1); Fs = 1.25
GCI_f = Fs*abs((f2-f1)/f1)/(r**p_RE-1); GCI_c = Fs*abs((f3-f2)/f2)/(r**p_RE-1)
asym = GCI_c/(r**p_RE*GCI_f); M_true = M_of_AR(AEX/AT, data[-1]["g"], True)
print(f"# Richardson(exit Mach): p={p_RE:.3f}  M_extrap={f_ext:.6f}  M_true={M_true:.6f}  "
      f"|d|={abs(f_ext-M_true):.1e}  GCI_fine={100*GCI_f:.4f}%  asympt={asym:.3f}")

# ======================= publication style (mathtext, not full LaTeX) =======================
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm",
    "font.size": 11, "axes.labelsize": 12.5, "axes.titlesize": 12,
    "legend.fontsize": 9.5, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "axes.linewidth": 0.9, "lines.linewidth": 1.9, "lines.markersize": 5.5,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "figure.dpi": 140, "savefig.dpi": 300, "savefig.bbox": "tight",
    "legend.framealpha": 0.96, "legend.edgecolor": "0.6", "legend.fancybox": False,
})
NAVY, CRIM, FOREST, ORANGE, PURPLE = "#1f3b73", "#b2182b", "#1b7837", "#d8731a", "#762a83"

# ---------- FIGURE 1 : nozzle solution ----------
sol = next(r for r in data if r["NI"] == 200)
xa = np.linspace(0, 1, 500); ga = sol["g"]
Ma = mach_curve(xa, ga); tf = 1+0.5*(ga-1)*Ma**2; pa = tf**(-ga/(ga-1)); ARa = area(xa*L)/AT
ev = max(1, sol["NI"]//26); xf = sol["xc"]/L

fig = plt.figure(figsize=(6.7, 7.3))
gs = fig.add_gridspec(3, 1, height_ratios=[0.85, 2.2, 1.25], hspace=0.14)
axg = fig.add_subplot(gs[0])
axg.fill_between(xa, 0.5*ARa, -0.5*ARa, color="#d9e1ee", lw=0)
axg.plot(xa, 0.5*ARa, color=NAVY, lw=1.4); axg.plot(xa, -0.5*ARa, color=NAVY, lw=1.4)
axg.axvline(XT/L, color="0.45", ls=(0, (4, 3)), lw=1.0)
axg.text(XT/L+0.015, 1.0, "throat ($M=1$)", fontsize=9, color="0.3", va="top")
axg.set_ylabel(r"$\pm\,A/A^{*}$"); axg.set_xlim(0, 1); axg.set_yticks([])
axg.set_title(r"Planar q2D area-variation nozzle: CFD vs quasi-1D area law  ($\gamma=%.3f$)" % ga)
axg.tick_params(labelbottom=False)

axm = fig.add_subplot(gs[1], sharex=axg)
axm.plot(xa, Ma, "-", color=NAVY, label=r"Mach — area law")
axm.plot(xf[::ev], sol["M_q"][::ev], "o", color=NAVY, mfc="white", mew=1.2,
         label=r"Mach — Q2D ($N=200$)")
axm.axvline(XT/L, color="0.45", ls=(0, (4, 3)), lw=1.0)
axm.set_ylabel("Mach number"); axm.tick_params(labelbottom=False)
axp = axm.twinx()
axp.plot(xa, pa, "--", color=CRIM, label=r"$p/p_0$ — area law")
axp.plot(xf[::ev], sol["p"][::ev]/1.0e6, "s", color=CRIM, mfc="white", mew=1.2,
         label=r"$p/p_0$ — Q2D")
axp.set_ylabel(r"$p/p_0$", color=CRIM); axp.tick_params(axis="y", colors=CRIM, direction="in")
axp.spines["right"].set_color(CRIM); axp.minorticks_on()
h1, l1 = axm.get_legend_handles_labels(); h2, l2 = axp.get_legend_handles_labels()
axm.legend(h1+h2, l1+l2, loc="center right", bbox_to_anchor=(0.99, 0.46))

axe = fig.add_subplot(gs[2], sharex=axg)
for rr, c, lbl in ((data[2], NAVY, r"$N=200$"), (data[4], ORANGE, r"$N=800$")):
    axe.semilogy(rr["xc"]/L, 100*rr["relerr"], "-", color=c, lw=1.6, label=lbl)
axe.axvline(XT/L, color="0.45", ls=(0, (4, 3)), lw=1.0)
axe.set_xlabel(r"$x/L$"); axe.set_ylabel(r"$|M_\mathrm{CFD}-M_\mathrm{exact}|\,/\,M$   [%]")
axe.set_xlim(0, 1); axe.legend(loc="upper right", ncol=2)
axe.annotate(r"throat: $1^\mathrm{st}$ order", xy=(XT/L, 100*data[2]["relerr"].max()),
             xytext=(XT/L+0.07, 100*data[2]["relerr"].max()*0.45),
             fontsize=9, color="0.3",
             arrowprops=dict(arrowstyle="->", color="0.55", lw=0.9))
fig.savefig(os.path.join(HERE, "nozzle_solution.png"))
fig.savefig(os.path.join(HERE, "nozzle_solution.svg"))

# ---------- FIGURE 2 : grid convergence ----------
fig2, ax = plt.subplots(figsize=(6.7, 5.6))
series = [("Linf",   r"$L_\infty$ (incl. throat)", CRIM,   "o"),
          ("Linf_s", r"$L_\infty$ smooth region",  FOREST, "^"),
          ("L2",     r"$L_2$",                     NAVY,   "s"),
          ("L1",     r"$L_1$",                     PURPLE, "D")]
for key, lbl, c, mk in series:
    e = np.array([r[key] for r in data])
    ax.loglog(dxs, e, mk+"-", color=c, mfc="white", mew=1.3, label=f"{lbl},  $p={orders[key]:.2f}$")

def slope_tri(x0, y0, slope, c, lab):
    x1 = x0*2.0; y1 = y0*2.0**slope
    ax.plot([x0, x1, x1, x0], [y0, y0, y1, y0], color=c, lw=1.0)
    ax.text(x1*1.05, math.sqrt(y0*y1), lab, color=c, fontsize=10, va="center")
LinfA = np.array([r["Linf"] for r in data]); L1A = np.array([r["L1"] for r in data])
slope_tri(dxs[2], LinfA[2]*1.45, 1.0, "0.35", "1")
slope_tri(dxs[1], L1A[1]*0.32,  2.0, "0.35", "2")

ax.set_xlabel(r"grid spacing  $\Delta x = L/N$")
ax.set_ylabel(r"Mach error norm vs area law")
ax.set_title(r"Grid convergence — MUSCL $2^\mathrm{nd}$ order;  throat corner $\Rightarrow L_\infty\to 1$")
ax.grid(True, which="major", ls="-", lw=0.4, alpha=0.35)
ax.grid(True, which="minor", ls=":", lw=0.3, alpha=0.25)
ax.legend(loc="upper left", title="Mach-error convergence")
box = ("Richardson / GCI  (exit Mach)\n"
       f"observed order  $p = {p_RE:.2f}$\n"
       f"$M_\\mathrm{{extrap}} = {f_ext:.5f}$  vs exact ${M_true:.5f}$\n"
       f"$\\mathrm{{GCI}}_\\mathrm{{fine}} = {100*GCI_f:.4f}$ %,  asympt. ${asym:.2f}$\n"
       r"kink: $L_\infty\sim\Delta x,\ L_2\sim\Delta x^{3/2},\ L_1\sim\Delta x^{2}$")
ax.text(0.975, 0.04, box, transform=ax.transAxes, fontsize=8.4, va="bottom", ha="right",
        bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="0.6", lw=0.8))
fig2.savefig(os.path.join(HERE, "nozzle_convergence.png"))
fig2.savefig(os.path.join(HERE, "nozzle_convergence.svg"))
print(f"# figures: {HERE}/nozzle_solution.png  {HERE}/nozzle_convergence.png  (+ .svg)")
shutil.rmtree(SCRATCH, ignore_errors=True)
