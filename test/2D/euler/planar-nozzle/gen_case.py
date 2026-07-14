# -*- coding: utf-8 -*-
#
# planar-nozzle — generate a q2D area-variation converging-diverging nozzle case
# (constant-gamma air). The q2D out-of-plane thickness h(x)=Cross_Section plays the
# role of the quasi-1D area A(x): h-weighted continuity for y-uniform flow gives
# h*rho*u = mdot, so the converged inviscid solution must match the quasi-1D isentropic
# area-law (Shapiro). Used by check.py (acceptance) and verify.py (grid convergence).
#
# Usage:  python3 gen_case.py [NI] [NJ] [target_dir]
#         defaults: NI=200, NJ=4, target_dir = this script's directory
#
# Gas: the VERIFIED constant-Cp air of test/2D/euler/oblique-shock (Cp=1004.5 -> gamma=1.4).

import os, sys, math, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
NI   = int(sys.argv[1]) if len(sys.argv) > 1 else 200
NJ   = int(sys.argv[2]) if len(sys.argv) > 2 else 4
TGT  = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else HERE
OBLIQUE = os.path.join(HERE, "..", "oblique-shock", "INPUT")

# --- gas (match the table exactly) ---
RU, WMOL, CP = 8314.462618, 28.970418, 1004.5
RGAS  = RU / WMOL
GAMMA = CP / (CP - RGAS)
P0, T0 = 1.0e6, 300.0
RHO0 = P0 / (RGAS * T0)

# --- C-D nozzle area law A(x)=h(x): smooth, dA/dx=0 at the throat (sonic point) ---
L, XT, AT, AIN, AEX, W = 1.0, 0.40, 1.0, 2.0, 2.5, 0.10
def area(x):
    if x <= XT:
        return AT + (AIN-AT)*0.5*(1.0 + math.cos(math.pi*x/XT))
    return AT + (AEX-AT)*0.5*(1.0 - math.cos(math.pi*(x-XT)/(L-XT)))

def AR_of_M(M):
    g = GAMMA
    return (1.0/M)*((2.0/(g+1))*(1.0+0.5*(g-1)*M*M))**((g+1)/(2*(g-1)))
def M_of_AR(ar, sup):
    if abs(ar-1.0) < 1e-12: return 1.0
    lo, hi = (1.0+1e-9, 60.0) if sup else (1e-7, 1.0-1e-9)
    for _ in range(200):
        mid = 0.5*(lo+hi); f = AR_of_M(mid)-ar
        if sup: hi, lo = (mid, lo) if f > 0 else (hi, mid)
        else:   lo, hi = (mid, hi) if f > 0 else (lo, mid)
    return 0.5*(lo+hi)
def prim(x):
    ar = area(x)/AT; M = M_of_AR(ar, x > XT)
    tf = 1.0 + 0.5*(GAMMA-1)*M*M
    p = P0*tf**(-GAMMA/(GAMMA-1)); T = T0/tf
    return p/(RGAS*T), M*math.sqrt(GAMMA*RGAS*T), 0.0, p     # rho,u,v,p

os.makedirs(os.path.join(TGT, "INPUT"), exist_ok=True)
os.makedirs(os.path.join(TGT, "OUTPUT"), exist_ok=True)
os.makedirs(os.path.join(TGT, "bin"), exist_ok=True)

xn = [i*L/NI for i in range(NI+1)]
yn = [j*W/NJ for j in range(NJ+1)]
xc = [0.5*(xn[i]+xn[i+1]) for i in range(NI)]
fm = lambda v: f"{v: .15E}"

# ic.tec  (BLOCK; x,y NODAL ; rho1,u,v,p CELLCENTERED ; i fastest)
with open(os.path.join(TGT, "INPUT", "ic.tec"), "w") as f:
    f.write(' VARIABLES ="x"  "y"  "rho1"  "u"  "v"  "p"\n')
    f.write(f' ZONE  T = B1-IG, I={NI+1}, J={NJ+1}, K=1, DATAPACKING=BLOCK, '
            'VARLOCATION=([1-2]=NODAL,[3-6]=CELLCENTERED), SOLUTIONTIME=0.000000000000000E+000\n')
    for j in range(NJ+1):
        for i in range(NI+1): f.write(fm(xn[i]) + "\n")
    for j in range(NJ+1):
        for i in range(NI+1): f.write(fm(yn[j]) + "\n")
    for vidx in range(4):
        for j in range(NJ):
            for i in range(NI): f.write(fm(prim(xc[i])[vidx]) + "\n")

# block1_area.dat  (nodal h, ((A(i,j),i=0,ni),j=0,nj), i fastest)
with open(os.path.join(TGT, "INPUT", "block1_area.dat"), "w") as f:
    for j in range(NJ+1):
        for i in range(NI+1): f.write(f"  {area(xn[i]):.15E}\n")

# bc.txt  (b i j k f code) : face1 inflow(401, T0/p0) ; face2 extrap(400) ; face3,4 sym(300)
with open(os.path.join(TGT, "INPUT", "bc.txt"), "w") as f:
    def rec(b, i, j, k, fc, code, prop=None):
        f.write(f"{b:8d}{i:8d}{j:8d}{k:8d}{fc:8d}{code:8d}\n")
        if prop is not None: f.write(prop + "\n")
    inflow = (f"  {T0:.6E},  {P0:.6E},         normal,         normal,"
              "    0.100000E+01,    0.100000E+01,")
    for j in range(1, NJ+1): rec(1, 1,  j, 1, 1, 401, inflow)
    for j in range(1, NJ+1): rec(1, NI, j, 1, 2, 400)
    for i in range(1, NI+1): rec(1, i,  1, 1, 3, 300)
    for i in range(1, NI+1): rec(1, i, NJ, 1, 4, 300)

# thermo / phase : the VERIFIED constant-gamma air table
shutil.copy(os.path.join(OBLIQUE, "thermo.dat"), os.path.join(TGT, "INPUT", "thermo.dat"))
shutil.copy(os.path.join(OBLIQUE, "phase.txt"),  os.path.join(TGT, "INPUT", "phase.txt"))

# input.ini  (euler + area_variation ; numerics mirrored from working euler cases)
with open(os.path.join(TGT, "input.ini"), "w") as f:
    f.write("[Q2D-Parameters]\niter-threshold = 40000\n\n"
            "[Q2D-Physics]\nequations = euler\n\n"
            "[Q2D-Numerics]\ncfl = 0.5\ncfl-rise-threshold = 50\ntime-scheme = RK3\n"
            "time-accurate = false\nirs-beta = 0.5\nintegration-variables = cons\n"
            "space-reconstruction = MUSCL\nflux-limiter = vanleer\nriemann-solver = HLLC\n\n"
            "[Q2D-IO]\nshell-diter = 500\n\n"
            "[Q2D-Sources]\nactive = area_variation\n")

if __name__ == "__main__":
    print(f"planar-nozzle: NI={NI} NJ={NJ} -> {TGT}")
    print(f"  gamma={GAMMA:.6f} R={RGAS:.4f} P0={P0:.3e} T0={T0}")
    print(f"  throat x={XT} A*={AT} ; inlet AR={AIN} (M={M_of_AR(AIN,False):.4f}) ; exit AR={AEX} (M={M_of_AR(AEX,True):.4f})")
