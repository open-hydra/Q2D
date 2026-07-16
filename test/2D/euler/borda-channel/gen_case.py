# -*- coding: utf-8 -*-
#
# borda-channel — regression case for the Borda sudden-expansion injector
# inflow BC (BC_Nozzle_Borda_Flow; Fievisohn 2016 thesis Sec 3.3,
# eqs 3.26-3.31). MOSE has no nozzle/Borda BC. The BC is a quasi-1D
# analytic jump: plenum (P0,T0) -> isentropic nozzle to throat (station 1, area A1,
# choked M1=1) -> sudden expansion (A2=A3-A1) -> channel face (station 3, area A3),
# with r = A1/A3 the area ratio carried by the bc.txt `pamb` field. Activated by the
# sentinel BC_Mach == -20 on BC code 420.
#
# Geometry: a straight 2D plane channel (NO area-variation source). Face 1 = Borda
# choked inflow; face 2 = supersonic extrapolation outflow; faces 3,4 = symmetry.
# Interior is initialised at a subsonic back-pressure P3 = 0.70*P0 so the inlet runs
# in regime (A): CHOKED throat (M1=1) with a SUBSONIC station 3 — the regime in which
# the Fievisohn jump algebra is valid and the analytic anchors below are exact.
#
# Gas: the common constant-Cp air table (test/common/Air, Cp=1004.5 -> gamma=1.4,
# R=287.0), symlinked into INPUT by Q2D.sh (mirror ramp-channel). The Fievisohn test
# case (Fig 3.6, P0=10 atm, T0=300 K, A1/A3=0.2) uses H2-air whose gamma,R are never
# numerically given in the thesis; eqs 3.26-3.31 are gamma-parametric, so evaluating
# them at a fully-specified gas gives EXACT analytic references while reproducing the
# paper's regime structure (A1/A3=0.2 kept). See docs/vv/2D-borda-channel.md.
#
# Analytic choked anchors (M1=1, used by check.py — derived, not invented):
#   T1/T0 = 2/(g+1);  P1/P0 = (2/(g+1))^(g/(g-1));  v1 = a1 = sqrt(g*R*T1)
#   choked mass flux (the conserved jump quantity, Eq 3.29 @ M1=1, per unit A3):
#     mdot_A3 = rho1 * v1 * (A1/A3)
#
# Usage:  python3 gen_case.py [NI] [NJ] [target_dir]
#         defaults: NI=80, NJ=8, target_dir = this script's directory

import os, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
NI   = int(sys.argv[1]) if len(sys.argv) > 1 else 80
NJ   = int(sys.argv[2]) if len(sys.argv) > 2 else 8
TGT  = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else HERE

# --- gas (must match test/common/Air: Cp=1004.5, air 28.970418) ---
RU, WMOL, CP = 8314.462618, 28.970418, 1004.5
RGAS  = RU / WMOL
GAMMA = CP / (CP - RGAS)

# --- Borda plenum + area ratio (paper r = A1/A3 = 0.2) ---
P0, T0, R_AREA = 1.0e6, 300.0, 0.2

# --- straight channel geometry ---
L, W = 1.0, 0.10

# --- interior IC: subsonic back-pressure P3 = 0.70*P0 -> regime (A) choked+subsonic-3 ---
PINT, TINT, UINT = 0.70e6, 285.0, 50.0
RHO_INT = PINT / (RGAS * TINT)

os.makedirs(os.path.join(TGT, "INPUT"), exist_ok=True)
os.makedirs(os.path.join(TGT, "OUTPUT"), exist_ok=True)
os.makedirs(os.path.join(TGT, "bin"), exist_ok=True)

xn = [i*L/NI for i in range(NI+1)]
yn = [j*W/NJ for j in range(NJ+1)]
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
    for val in (RHO_INT, UINT, 0.0, PINT):
        for j in range(NJ):
            for i in range(NI): f.write(fm(val) + "\n")

# bc.txt  (b i j k f code) : face1 Borda(420) ; face2 extrap(400) ; face3,4 sym(300)
# BC 420 property line (reader order): mach, T0, p0, alpha, beta, pamb(=A1/A3), rel_fac, ci
with open(os.path.join(TGT, "INPUT", "bc.txt"), "w") as f:
    def rec(b, i, j, k, fc, code, prop=None):
        f.write(f"{b:8d}{i:8d}{j:8d}{k:8d}{fc:8d}{code:8d}\n")
        if prop is not None: f.write(prop + "\n")
    borda = (f"   -0.200000E+02,  {T0:.6E},  {P0:.6E},         normal,         normal,"
             f"  {R_AREA:.6E},    0.100000E+01,    0.100000E+01,")
    for j in range(1, NJ+1): rec(1, 1,  j, 1, 1, 420, borda)
    for j in range(1, NJ+1): rec(1, NI, j, 1, 2, 400)
    for i in range(1, NI+1): rec(1, i,  1, 1, 3, 300)
    for i in range(1, NI+1): rec(1, i, NJ, 1, 4, 300)

# input.ini  (plane euler ; NO area_variation ; numerics mirrored from working euler cases)
with open(os.path.join(TGT, "input.ini"), "w") as f:
    f.write("[Q2D-Parameters]\niter-threshold = 5000\n\n"
            "[Q2D-Physics]\nequations = euler\n\n"
            "[Q2D-Numerics]\ncfl = 0.5\ncfl-rise-threshold = 50\ntime-scheme = RK3\n"
            "time-accurate = false\nirs-beta = 0.5\nintegration-variables = cons\n"
            "space-reconstruction = MUSCL\nflux-limiter = vanleer\nriemann-solver = HLLC\n\n"
            "[Q2D-IO]\nshell-diter = 500\n")

if __name__ == "__main__":
    T1 = T0*2/(GAMMA+1); P1 = P0*(2/(GAMMA+1))**(GAMMA/(GAMMA-1))
    rho1 = P1/(RGAS*T1); a1 = math.sqrt(GAMMA*RGAS*T1)
    print(f"borda-channel: NI={NI} NJ={NJ} -> {TGT}")
    print(f"  gamma={GAMMA:.6f} R={RGAS:.4f} P0={P0:.3e} T0={T0} A1/A3={R_AREA}")
    print(f"  choke anchors: T1={T1:.3f} P1={P1:.1f} v1=a1={a1:.4f} mdot_A3={rho1*a1*R_AREA:.4f}")
