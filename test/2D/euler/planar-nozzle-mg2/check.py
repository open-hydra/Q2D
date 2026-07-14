"""planar-nozzle-mg2 acceptance check (ctest) — peculiarita' #3: area-variation x Multigrid.

Same q2D area-variation C-D nozzle as ``test/2D/euler/planar-nozzle`` (the out-of-plane
thickness h(x)=Cross_Section plays the role of the quasi-1D area A(x); the converged
inviscid solution must match the analytic isentropic area-law, Shapiro), but run with
2-level Multigrid (``[Q2D-Multigrid] levels=2``).

This is the integration test deferred as T-10 from the B-10 fix
(``docs/review/multigrid_area_variation_restriction.md`` Section 4.2). It verifies the
*interpretation* of the feature, not just a value-match:

  1. (convergence)  the MGL=2 area-variation run reaches steady state. With the original
     B-10 bug (Cross_Section re-read with coarse dims -> scrambled re-indexing, Section 2.2)
     the coarse level carries a meaningless geometric vector G = -grad(A)/A and the solver
     stalls / diverges; convergence to the residual floor is the sentinel that the R1
     restriction (Coarse_Cross_Section) is correct.
  2. (physical)    the converged fine field matches the SAME analytic quasi-1D area law as
     the MGL=1 baseline, within the SAME tolerances -> multigrid does not corrupt the
     area-variation physics.
  3. (multigrid)   the coarse level (field-level2.tec) is a consistent 2:1 restriction of
     the fine field: the coarse Mach distribution reproduces the fine Mach averaged over
     each coarse cell pair to O(h_c^2). This is the direct evidence that Cross_Section was
     sub-sampled (A_c(I,J)=A_f(2I,2J), Section 1.2), not scrambled -- a wrong restriction
     would break the coarse residual operator and the coarse Mach would not track the fine.

Thresholds are the baseline planar-nozzle area-law tolerances (physical reference, Shapiro)
plus an O(h_c^2)-scale coarse/fine restriction tolerance anchored to the committed run; none
are invented for this case.
"""
import os
import sys
import math
import numpy as np

HERE  = os.path.dirname(os.path.abspath(__file__))
FIELD = os.path.join(HERE, "OUTPUT", "field.tec")
COARSE = os.path.join(HERE, "OUTPUT", "field-level2.tec")
RESID = os.path.join(HERE, "OUTPUT", "residual-history.dat")

# --- nozzle area law (identical to planar-nozzle gen_case.py) ---
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

def parse_field(path, nnodal):
    """Parse a Q2D field.tec. nnodal = number of NODAL leading vars (fine=2 [x,y];
    coarse field-levelN=3 [x,y,z]). Returns NI, cell-x-centres, and per-i column means
    of rho,u,v,p,T plus median(g),median(R)."""
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
    xn = take(nn).reshape(NJp, NIp)
    for _ in range(nnodal-1): take(nn)            # skip y (and z when present)
    rho, u = take(nc).reshape(NJ, NI), take(nc).reshape(NJ, NI)
    v_ = take(nc).reshape(NJ, NI); p = take(nc).reshape(NJ, NI)
    T = take(nc).reshape(NJ, NI)
    g = take(nc).reshape(NJ, NI); R = take(nc).reshape(NJ, NI)
    m = lambda a: a.mean(axis=0)
    M = np.sqrt(u*u + v_*v_)/np.sqrt(g*R*T)
    return (NI, 0.5*(xn[0, :-1]+xn[0, 1:]),
            m(rho), m(u), m(p), m(T), m(M),
            float(np.median(g)), float(np.median(R)))

def fail(msg):
    print(f"[planar-nozzle-mg2] FAIL: {msg}")
    sys.exit(1)

if not os.path.isfile(FIELD):
    fail(f"missing {FIELD} (did ./Q2D.sh solve run?)")
if not os.path.isfile(COARSE):
    fail(f"missing {COARSE} -- the coarse level was not written; MGL=2 path did not run")

# ---------- layer 1: convergence (B-10 sentinel) ----------
# The density residual must drop to steady state. With the B-10 bug the coarse-level
# source is meaningless and the run stalls / diverges instead of converging.
if not os.path.isfile(RESID):
    fail(f"missing {RESID}")
res = np.loadtxt(RESID)
# col 0 = iter index; remaining columns are residuals (cf. chang-nozzle-mg2 check.py).
res_norm = np.linalg.norm(res[:, 1:], axis=1)
first = res_norm[0] if res_norm[0] > 0 else 1.0
drop = res_norm[-1] / first
dens_resid = float(np.min(np.abs(res[-1, 1:])))   # density residual at last iter

# ---------- layer 2: physical -- fine field matches the analytic area law ----------
NI, xc, rho, u, p, T, M_q, g, R = parse_field(FIELD, 2)

for nm, a in (("rho", rho), ("p", p), ("T", T)):
    if np.min(a) <= 0.0:
        fail(f"non-physical min({nm}) = {np.min(a):.3e}")

M_a = np.array([M_of_AR(area(x)/AT, g, x > XT) for x in xc])
rel = np.abs(M_q-M_a)/M_a
interior = (xc > 2*L/NI) & (xc < L-2*L/NI)
smooth = interior & (np.abs(xc-XT) > 0.05)
exit_err = abs(M_q[-1]-M_a[-1])/M_a[-1]
smooth_max = float(np.max(rel[smooth]))
whole_max = float(np.max(rel[interior]))

# ---------- layer 3: multigrid -- coarse level reproduces the fine restriction ----------
# field-level2.tec carries a 3rd nodal column (z); coarse NI = fine NI / 2.
NIc, xcc, _, _, _, _, M_c, _, _ = parse_field(COARSE, 3)
if NIc != NI // 2:
    fail(f"coarse NI={NIc} is not fine NI/2={NI//2}: 2:1 coarsening broken")
# R1 restriction of the fine Mach column means onto the coarse cells: average each
# fine cell pair (2I-1,2I). A correct Cross_Section restriction makes the coarse Mach
# track this to O(h_c^2); a scrambled one (the B-10 bug) does not.
M_q_restr = M_q.reshape(NIc, 2).mean(axis=1)
cf_rel = np.abs(M_c - M_q_restr)/np.maximum(np.abs(M_q_restr), 1e-30)
cf_max = float(np.max(cf_rel))

print(f"[planar-nozzle-mg2] N={NI} (coarse {NIc}) gamma={g:.5f} R={R:.3f} | "
      f"resid_drop={drop:.2e} dens_resid={dens_resid:.2e} | "
      f"M_exit={M_q[-1]:.5f} (analytic {M_a[-1]:.5f}) | "
      f"exit_relerr={exit_err:.2e} smooth_max={smooth_max:.2e} whole_max={whole_max:.2e} | "
      f"coarse/fine Mach max_rel={cf_max:.2e}")

ok = True
# layer 1 -- convergence sentinel (B-10: must converge, not stall/diverge).
if not np.isfinite(dens_resid) or dens_resid > 1.0e-6:
    ok = False; print(f"  MGL=2 not converged: density residual {dens_resid:.2e} > 1.0e-6 "
                       f"(B-10 regression? coarse Cross_Section may be scrambled)")
# layer 2 -- same area-law tolerances as the MGL=1 baseline (planar-nozzle check.py).
if exit_err   > 2.0e-3: ok = False; print(f"  exit Mach rel err {exit_err:.2e} > 2.0e-3")
if smooth_max > 5.0e-3: ok = False; print(f"  smooth-region max rel err {smooth_max:.2e} > 5.0e-3")
if whole_max  > 1.5e-2: ok = False; print(f"  whole-domain max rel err {whole_max:.2e} > 1.5e-2")
# layer 3 -- coarse level reproduces fine within O(h_c^2). Observed ~4.5e-4 on the
# committed run; 5e-3 is a regression band (one decade of headroom) that a scrambled
# coarse Cross_Section (the B-10 failure) blows through immediately.
if not np.isfinite(cf_max) or cf_max > 5.0e-3:
    ok = False; print(f"  coarse/fine Mach max rel {cf_max:.2e} > 5.0e-3 "
                      f"(coarse level does not reproduce the fine restriction)")

if not ok:
    fail("MGL=2 area-variation run does not satisfy convergence / area-law / "
         "coarse-restriction acceptance")

print("[planar-nozzle-mg2] PASS (area-variation x multigrid: converged, matches the "
      "quasi-1D area law, coarse level reproduces the fine restriction)")
sys.exit(0)
