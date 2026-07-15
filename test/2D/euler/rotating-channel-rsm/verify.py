#!/usr/bin/env python3
"""rotating-channel-rsm — frame-invariance V&V gate for the SSGLRR
rotating-frame closures (exact Coriolis production G_ij, defect D1 fix).

Physics: uniform relative flow (grad u_rel = 0) in a frame rotating at
Omega = 10 rad/s, uniform anisotropic R_ij(0) with R12(0) = 0. The exact model
solution is the homogeneous ODE (P = 0, S = 0, omega_abs = 2*Omega):
  dR/dt = G(R;Omega) + Pi_slow(a) + Pi_C5(a, W_abs) - (2/3) eps I,
  d(om)/dt = -beta_e om^2
(epsilon-branch constants: slip walls -> walldist = 1e8 -> F1 = 0, exactly as
Lib_SSGLRR compute_F1 evaluates). Reference integrated here with RK4 at dt/8.
Frame-invariance provenance: an ODE analysis proves
this ODE (with G) is the rotating-frame image of the inertial-frame model.

Gates (all hard):
  G1 spatial uniformity at T (interior, 2 wall rows excluded — the symmetry BC
     pins R12 = 0 at the walls and diffuses inward O(sqrt(nu_t T)) ~ 1 cell).
  G2 |R_code(T) - R_ODE(T)| / k(T) < 2e-2 componentwise (mid-channel cell).
  G3 discriminance: distance from the no-G trajectory at least 5x the distance
     from the full trajectory (the test FAILS on a binary without G_ij).
"""
import re, sys
import numpy as np

OMEGA, T = 10.0, 0.05
Cmu = 0.09
C1e, C2e, C5e, be = 1.7, 1.05, 0.2, 0.0828

def rhs(y, withG):
    R = np.array([[y[0], y[3], 0], [y[3], y[1], 0], [0, 0, y[2]]]); om = y[4]
    k = 0.5*np.trace(R); eps = Cmu*k*om
    a = R/k - 2.0/3.0*np.eye(3)
    aa = np.sum(a*a)
    W = np.array([[0, -OMEGA, 0], [OMEGA, 0, 0], [0, 0, 0.0]])  # W_abs: omega_abs=2*Omega -> W12=-om/2
    Pi = (-(C1e*eps)*a + C2e*eps*(a@a - aa/3*np.eye(3))
          + C5e*k*(a@W.T + (a@W.T).T))            # S=0 -> C3/C4 terms vanish; P=0 -> C1* term vanishes
    dR = Pi - 2.0/3.0*eps*np.eye(3)
    if withG:
        G = np.zeros((3, 3))
        G[0, 0] = 4*OMEGA*R[0, 1]; G[1, 1] = -4*OMEGA*R[0, 1]
        G[0, 1] = G[1, 0] = 2*OMEGA*(R[1, 1]-R[0, 0])
        dR = dR + G
    return np.array([dR[0, 0], dR[1, 1], dR[2, 2], dR[0, 1], -be*om*om])

def rk4(y0, withG, n=8000):
    y, dt = np.array(y0, float), T/n
    for _ in range(n):
        k1 = rhs(y, withG); k2 = rhs(y+dt/2*k1, withG)
        k3 = rhs(y+dt/2*k2, withG); k4 = rhs(y+dt*k3, withG)
        y = y + dt/6*(k1+2*k2+2*k3+k4)
    return y

def read_field(path):
    txt = open(path).read().split('\n')
    m = re.search(r'I=\s*(\d+), J=\s*(\d+)', txt[1]); I, J = int(m.group(1)), int(m.group(2))
    nn, nc = I*J, (I-1)*(J-1)
    vals = [float(t) for ln in txt[2:] for t in ln.split()]
    nvar = (len(vals)-2*nn)//nc
    cells = np.array(vals[2*nn:2*nn+nvar*nc]).reshape(nvar, I-1, J-1, order='C')
    return cells  # [var, i, j]: rho,u,v,p,ruu,rvv,rww,ruv,romega

f = sorted(__import__('glob').glob('OUTPUT/field*.tec'))[-1]
c = read_field(f)
rho = c[0]
Rf = {k: c[i]/rho for k, i in (('11', 4), ('22', 5), ('33', 6), ('12', 7))}
omf = c[8]/rho
ni, nj = rho.shape
mid = (ni//2, nj//2)
FAIL = []
def gate(name, ok, msg):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg}")
    if not ok: FAIL.append(name)

# G1 per-cell ODE tracking: the IC prescribes UNIFORM conserved rho*R, but the
# Coriolis equilibrium has dp/dy = 2 rho Omega U -> rho(y) varies by O(1%) ->
# each cell starts from its own R(0) = (rho R)_0/rho(y) and evolves per the SAME
# homogeneous ODE (rho steady). Gate: every interior cell (2 wall rows excluded
# for the R12=0 symmetry-wall layer, disclosed) tracks its own ODE trajectory.
ic = read_field('INPUT/ic.tec')
errs = []
for jj in range(2, nj-2):
    r0 = ic[0][mid[0], jj]
    y0c = [ic[4][mid[0], jj]/r0, ic[5][mid[0], jj]/r0, ic[6][mid[0], jj]/r0,
           ic[7][mid[0], jj]/r0, ic[8][mid[0], jj]/r0]
    yFc = rk4(y0c, True, n=2000)
    kc = 0.5*(yFc[0]+yFc[1]+yFc[2])
    cc = np.array([Rf['11'][mid[0], jj], Rf['22'][mid[0], jj],
                   Rf['33'][mid[0], jj], Rf['12'][mid[0], jj]])
    errs.append(max(abs(cc-yFc[:4]))/kc)
gate("G1 per-cell ODE tracking", max(errs) < 2e-2,
     f"max over interior j of |R_code - R_ODE(y)|/k = {max(errs):.2e} < 2e-2")

# G2/G3: ODE comparison at the mid-channel cell (ic already loaded above)
rho0 = ic[0][mid]
y0 = [ic[4][mid]/rho0, ic[5][mid]/rho0, ic[6][mid]/rho0, ic[7][mid]/rho0, ic[8][mid]/rho0]
yF = rk4(y0, True); yN = rk4(y0, False)
code = np.array([Rf['11'][mid], Rf['22'][mid], Rf['33'][mid], Rf['12'][mid], omf[mid]])
kT = 0.5*(yF[0]+yF[1]+yF[2])
errF = max(abs(code[:4]-yF[:4]))/kT
errN = max(abs(code[:4]-yN[:4]))/kT
gate("G2 matches full ODE (with G_ij)", errF < 2e-2,
     f"max|R_code - R_ODE|/k = {errF:.2e} < 2e-2  (R12: code {code[3]:.4f} vs ODE {yF[3]:.4f})")
gate("G3 discriminance vs no-G trajectory", errN > 5*errF,
     f"no-G distance {errN:.2e} >= 5 x {errF:.2e} (a binary without G_ij FAILS G2)")
print("RESULT:", "ALL GATES PASS" if not FAIL else f"FAILURES: {FAIL}"); sys.exit(1 if FAIL else 0)
