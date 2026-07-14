"""planar-nozzle 2D steady-state acceptance check (ctest).

q2D area-variation converging-diverging nozzle, constant-gamma air (oblique-shock
thermo table, gamma=1.4). The out-of-plane thickness h(x)=Cross_Section plays the role
of the quasi-1D area A(x): h-weighted continuity for y-uniform flow is d/dx(h*rho*u)=0,
so the converged inviscid solution must match the analytic isentropic area-law (Shapiro).

Validates the converged OUTPUT/field.tec against the analytic quasi-1D area law,
using gamma and R read from the field itself (isolates discretisation error):
  - positivity: min(rho), min(p), min(T) > 0           (hard)
  - exit Mach vs analytic                  rel err < 2.0e-3
  - smooth-region (|x-x_t|>0.05) Mach       max rel err < 5.0e-3
  - whole domain (incl. throat corner)      max rel err < 1.5e-2

The full grid-convergence study (Richardson/GCI, JPP plots) is in verify.py.
"""
import os
import sys
import math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIELD = os.path.join(HERE, "OUTPUT", "field.tec")

# --- nozzle area law (identical to gen_case.py) ---
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
    return (NI, 0.5*(xn[0, :-1]+xn[0, 1:]),
            m(rho), m(u), m(p), m(T), float(np.median(g)), float(np.median(R)))

def fail(msg):
    print(f"[planar-nozzle] FAIL: {msg}")
    sys.exit(1)

if not os.path.isfile(FIELD):
    fail(f"missing {FIELD} (did ./Q2D.sh solve run?)")

NI, xc, rho, u, p, T, g, R = parse_field(FIELD)

# positivity (hard)
for nm, a in (("rho", rho), ("p", p), ("T", T)):
    if np.min(a) <= 0.0:
        fail(f"non-physical min({nm}) = {np.min(a):.3e}")

M_q = u/np.sqrt(g*R*T)
M_a = np.array([M_of_AR(area(x)/AT, g, x > XT) for x in xc])
rel = np.abs(M_q-M_a)/M_a
interior = (xc > 2*L/NI) & (xc < L-2*L/NI)
smooth = interior & (np.abs(xc-XT) > 0.05)

exit_err = abs(M_q[-1]-M_a[-1])/M_a[-1]
smooth_max = float(np.max(rel[smooth]))
whole_max = float(np.max(rel[interior]))

print(f"[planar-nozzle] N={NI} gamma={g:.5f} R={R:.3f} | "
      f"M_exit={M_q[-1]:.5f} (analytic {M_a[-1]:.5f}) | "
      f"exit_relerr={exit_err:.2e} smooth_max={smooth_max:.2e} whole_max={whole_max:.2e}")

ok = True
if exit_err   > 2.0e-3: ok = False; print(f"  exit Mach rel err {exit_err:.2e} > 2.0e-3")
if smooth_max > 5.0e-3: ok = False; print(f"  smooth-region max rel err {smooth_max:.2e} > 5.0e-3")
if whole_max  > 1.5e-2: ok = False; print(f"  whole-domain max rel err {whole_max:.2e} > 1.5e-2")
if not ok:
    fail("Mach distribution does not match the quasi-1D area law within tolerance")

print("[planar-nozzle] PASS (q2D area-variation matches quasi-1D area law within tolerance)")
sys.exit(0)
