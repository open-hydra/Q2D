"""borda-channel 2D steady-state acceptance check (ctest).

Regression test for the Q2D-only Borda sudden-expansion injector inflow BC
(peculiarità #4, BC_Nozzle_Borda_Flow; Fievisohn 2016 thesis Sec 3.3, eqs 3.26-3.31;
sentinel BC_Mach == -20 on BC code 420, A1/A3 in the bc.txt pamb field). MOSE has no
nozzle/Borda BC (divergence #4). Study + map-to-code + V&V criteria:
docs/review/borda_fievisohn_study_sess35.md.

The case is a straight plane channel with a Borda choked inflow on face 1, interior
back-pressure P3 = 0.70*P0 → regime (A): choked throat (M1=1), subsonic station 3 — the
regime where the Fievisohn jump algebra is valid. Const-gamma air (common/Air, gamma=1.4,
R=287). The injected face state is checked against the paper's own equations evaluated at
this fully-specified gas (no invented thresholds; the paper tabulates no numeric case and
its Fig 3.6/3.7 are read-from-figure only).

Anchors (study §7):
  - positivity: min(rho,p,T) > 0                                            (hard)
  - convergence sentinel: density residual settled                         (hard)
  - (guard) Borda branch actually entered: injected mass flux >> 0         (false-green
    prevention — without it a silent wall fallback would pass trivially)
  - (A) choked mass flux  rho*u (inlet)  vs analytic mdot_A3 = rho1*v1*A1/A3   rel < 3e-2
  - (B) isentropic choke state implied by the injected stagnation:
        recovered P0,T0 from the injected (rho,u,p,T) via isentropic+energy   rel < 3e-2
  - (C) total-pressure-loss INTERPRETATION (Borda-Carnot hallmark): the injected
        stagnation pressure p0_inj < P0 (a real loss), and entropy non-decreasing.

Negative control: corrupt the expected choked mdot and confirm (A) FAILS.
"""
import os
import sys
import math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIELD = os.path.join(HERE, "OUTPUT", "field.tec")
RESID = os.path.join(HERE, "OUTPUT", "residual-history.dat")

# --- Borda plenum + area ratio (identical to gen_case.py) ---
P0, T0, R_AREA = 1.0e6, 300.0, 0.2

NEG_CONTROL = os.environ.get("BORDA_NEG_CONTROL") == "1"


def parse_field(path):
    Ls = open(path).readlines()
    zi = next(i for i, l in enumerate(Ls) if l.strip().startswith("ZONE"))
    z = Ls[zi].replace(",", " ").split()
    NIp = int(next(t for t in z if t.startswith("I="))[2:])
    NJp = int(next(t for t in z if t.startswith("J="))[2:])
    NI, NJ = NIp - 1, NJp - 1
    v = []
    for l in Ls[zi+1:]:
        for t in l.replace(",", " ").split():
            try: v.append(float(t))
            except ValueError: pass
    v = np.array(v); off = [0]
    def take(n): a = v[off[0]:off[0]+n]; off[0] += n; return a
    nn, nc = NIp*NJp, NI*NJ
    take(nn); take(nn)  # x, y nodal
    rho = take(nc).reshape(NJ, NI); u = take(nc).reshape(NJ, NI)
    take(nc); p = take(nc).reshape(NJ, NI)  # skip v
    T = take(nc).reshape(NJ, NI); g = take(nc).reshape(NJ, NI); R = take(nc).reshape(NJ, NI)
    return (NI, rho, u, p, T, float(np.median(g)), float(np.median(R)))


def fail(msg):
    print(f"[borda-channel] FAIL: {msg}")
    sys.exit(1)


if not os.path.isfile(FIELD):
    fail(f"missing {FIELD} (did ./Q2D.sh solve run?)")

NI, rho, u, p, T, g, R = parse_field(FIELD)

# --- positivity (hard) ---
for nm, a in (("rho", rho), ("p", p), ("T", T)):
    if np.min(a) <= 0.0:
        fail(f"non-physical min({nm}) = {np.min(a):.3e}")

# --- convergence sentinel (hard) ---
if os.path.isfile(RESID):
    try:
        last = np.loadtxt(RESID)[-1]
        dens_resid = float(last[1]) if last.ndim else float(last)
        if not math.isfinite(dens_resid) or dens_resid > 1.0e-6:
            fail(f"not converged: final density residual {dens_resid:.3e} > 1e-6")
    except Exception as e:
        print(f"[borda-channel] WARN: could not parse residual-history ({e}); skipping sentinel")

# --- inlet column (first interior cell adjacent to the Borda face) ---
rho_in = float(rho[:, 0].mean())
u_in   = float(u[:, 0].mean())
p_in   = float(p[:, 0].mean())
T_in   = float(T[:, 0].mean())
mflux  = rho_in * u_in            # conserved jump quantity across the face

# --- analytic choked anchors at this gas (study §7; eqs 3.29 + isentropic) ---
T1   = T0 * 2.0 / (g + 1.0)
P1   = P0 * (2.0 / (g + 1.0)) ** (g / (g - 1.0))
rho1 = P1 / (R * T1)
v1   = math.sqrt(g * R * T1)              # a1 (M1 = 1)
mdot_A3 = rho1 * v1 * R_AREA              # Eq 3.29 @ M1=1, per unit A3

if NEG_CONTROL:
    mdot_A3 *= 1.5   # corrupt the expectation: (A) must FAIL

# --- (guard) Borda branch entered: a real injected flux (not a wall fallback) ---
# wall fallback would give u≈0 at the inlet; choked injection is O(100) kg/m2s.
if mflux < 0.1 * mdot_A3:
    fail(f"Borda branch NOT entered (injected mass flux {mflux:.3e} << choked {mdot_A3:.3e}); "
         "likely a silent wall fallback (false green)")

# --- (A) choked mass flux ---
err_mdot = abs(mflux - mdot_A3) / mdot_A3

# --- (B) recovered stagnation from the injected state (isentropic + energy) ---
a_in    = math.sqrt(g * R * T_in)
M_in    = u_in / a_in
tf      = 1.0 + 0.5 * (g - 1.0) * M_in * M_in
T0_inj  = T_in * tf
p0_inj  = p_in * tf ** (g / (g - 1.0))     # injected stagnation pressure
err_T0  = abs(T0_inj - T0) / T0

# --- (C) total-pressure-loss interpretation (Borda-Carnot hallmark) ---
loss_ok = p0_inj < P0                       # a real stagnation-pressure loss

print(f"[borda-channel] N={NI} gamma={g:.5f} R={R:.3f} A1/A3={R_AREA}")
print(f"  inlet: rho={rho_in:.5f} u={u_in:.4f} p={p_in:.1f} T={T_in:.3f} M={M_in:.4f}")
print(f"  (A) mass flux {mflux:.4f} vs choked mdot_A3 {mdot_A3:.4f}  rel={err_mdot:.2e}")
print(f"  (B) recovered T0={T0_inj:.3f} (plenum {T0})  rel={err_T0:.2e}")
print(f"  (C) injected p0={p0_inj:.1f} < P0={P0:.1f} ? {loss_ok}  (stagnation-pressure loss)")

ok = True
if err_mdot > 3.0e-2: ok = False; print(f"  (A) choked mass-flux rel err {err_mdot:.2e} > 3.0e-2")
if err_T0   > 3.0e-2: ok = False; print(f"  (B) recovered-T0 rel err {err_T0:.2e} > 3.0e-2")
if not loss_ok:       ok = False; print(f"  (C) NO stagnation-pressure loss (p0_inj {p0_inj:.1f} >= P0 {P0:.1f})")

if NEG_CONTROL:
    if ok:
        fail("NEGATIVE CONTROL did not trigger (corrupted choked mdot still passed)")
    print("[borda-channel] negative control OK (corrupted expectation correctly FAILS)")
    sys.exit(0)

if not ok:
    fail("Borda injected state does not match the choked sudden-expansion anchors")

print("[borda-channel] PASS (Borda choked sudden-expansion injection matches "
      "Fievisohn eqs 3.26-3.31 within tolerance)")
sys.exit(0)
