# -*- coding: utf-8 -*-
#
# ShuOsher-LIMO3 — regression case for the third-order LIMO3 limiter
# (rlimiter_LIMO3, Cada & Torrilhon 2009). MOSE has no LIMO3
# (the limiter set there stops at the 6 TVD limiters); LIMO3 lives in
# numerics/space/Lib_Limiters.f90 and is dispatched by the uppercase case label
# 'LIMO3' in Mod_Space::Assign_Limiter.
#
# WHY Shu-Osher: it is THE textbook "shock + regular waves" problem and the
# canonical demonstrator of a high-order limiter. A Mach-3 shock runs into a
# sinusoidal density field; downstream of the shock the entropy waves steepen
# into high-frequency SMOOTH density oscillations. A 2nd-order TVD limiter
# (MC/van Leer/minmod) CLIPS those smooth extrema back to first order
# (φ(r<0)=0); LIMO3 does NOT clip them (φ(r<0)>0, third order at extrema), so
# it preserves a larger oscillation amplitude while keeping the main shock
# monotone. THIS no-clip-at-smooth-extrema property IS the behaviour that
# check.py verifies — not a bare value match.
#
# Shu & Osher, J. Comput. Phys. 83 (1989) 32, problem of §6 (Example 8):
#   domain x in [-5, 5], gamma = 1.4
#   x <  -4 : (rho,u,p) = (3.857143, 2.629369, 10.33333)   [post-Mach-3 shock]
#   x >= -4 : (rho,u,p) = (1 + 0.2*sin(5x), 0, 1)
#   final time t = 1.8
# At t=1.8 the shock sits near x≈2.4 and the boundary cells stay in constant
# regions, so plain extrapolation (BC 400) on the streamwise ends is exact; the
# spanwise faces are symmetry (BC 300), mirror Ein91/Tor99 1D layout.
#
# Gas: the common constant-Cp air table (test/common/Air, Cp=1004.5 -> gamma=1.4,
# R=287.06), symlinked into INPUT by Q2D.sh (mirror Sod79/Ein91). The Shu-Osher
# IC states are the standard non-dimensional values evaluated at gamma=1.4, which
# this table reproduces (gamma=Cp/(Cp-R)=1.4006, within rounding of the canonical
# 1.4 — the limiter discriminator is insensitive to that 4th-digit difference).
#
# Usage:  python3 gen_case.py [NI] [target_dir]
#         defaults: NI=400 (the ctest mesh), target_dir = INPUT here.
#   The committed limiter baselines (reference/baseline_*.dat) were generated
#   once with the same binary at NI=400 using the mc/vanleer/minmod/iord limiters.

import os, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
NI   = int(sys.argv[1]) if len(sys.argv) > 1 else 400
TGT  = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(HERE, "INPUT")
NJ   = 1                      # pseudo-1D: a single interior cell row (J=2 nodes)

# --- domain (Shu-Osher) ---
XL, XR = -5.0, 5.0
XSHOCK = -4.0

def ic_state(xc):
    """Cell-centered Shu-Osher initial (rho,u,v,p) at cell-center xc."""
    if xc < XSHOCK:
        return 3.857143, 2.629369, 0.0, 10.33333
    else:
        return 1.0 + 0.2*math.sin(5.0*xc), 0.0, 0.0, 1.0

def write_case(target, ni):
    os.makedirs(target, exist_ok=True)
    dx = (XR - XL) / ni
    xn = [XL + i*dx for i in range(ni+1)]          # node coords (i fastest)
    yn = [0.0, 1.0]                                 # J=2 nodes, unit span
    xc = [0.5*(xn[i] + xn[i+1]) for i in range(ni)] # cell centers
    fm = lambda v: f"{v: .15E}"

    # ic.tec  (BLOCK; x,y NODAL ; rho1,u,v,p CELLCENTERED ; i fastest)
    with open(os.path.join(target, "ic.tec"), "w") as f:
        f.write(' VARIABLES ="x" "y"  "rho1" "u" "v" "p"\n')
        f.write(f' ZONE  T = B1-IG, I={ni+1}, J={NJ+1}, K=1, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-2]=NODAL,[3-6]=CELLCENTERED), '
                'SOLUTIONTIME=0.000000000000000E+000\n')
        for j in range(NJ+1):
            for i in range(ni+1): f.write(fm(xn[i]) + "\n")
        for j in range(NJ+1):
            for i in range(ni+1): f.write(fm(yn[j]) + "\n")
        # cell data blocks: rho, u, v, p (each ni cells, j-row single)
        states = [ic_state(x) for x in xc]
        for comp in range(4):
            for s in states: f.write(fm(s[comp]) + "\n")

    # bc.txt : face1,2 extrapolation(400) streamwise ends ; face3,4 symmetry(300)
    #          mirror Ein91/Tor99 1D layout (1D = symmetry top/bottom).
    # Q2D bc.txt head is ALWAYS 5 ints (b, i, j, f, type) -- no K column
    # (see IO_BC.f90 Check_BC); the 6-column MOSE 3D layout crashes the
    # metrics setup with misparsed face codes.
    with open(os.path.join(target, "bc.txt"), "w") as f:
        def rec(b, i, j, fc, code):
            f.write(f"{b:8d}{i:8d}{j:8d}{fc:8d}{code:8d}\n")
        rec(1, 1,  1, 1, 400)
        rec(1, ni, 1, 2, 400)
        for i in range(1, ni+1): rec(1, i, 1, 3, 300)
        for i in range(1, ni+1): rec(1, i, 1, 4, 300)

if __name__ == "__main__":
    write_case(TGT, NI)
    g = 1.4
    print(f"ShuOsher-LIMO3: NI={NI} NJ={NJ} domain [{XL},{XR}] -> {TGT}")
    print(f"  IC: x<{XSHOCK}: (3.857143, 2.629369, 10.33333) ; "
          f"x>= : (1+0.2 sin5x, 0, 1)  t_final=1.8")
