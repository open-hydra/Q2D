# -*- coding: utf-8 -*-
#
# noh-cylindrical — cylindrical Noh (1987) implosion, the axisymmetric-source
# verification case with an ANALYTIC exact solution (the only axi case in the
# suite not gated against digitized published data or a committed reference
# field).  A uniform cold gas (rho0 = 1, near-vacuum p0 — see the P0 note
# below) moves radially toward the symmetry axis at u_r = -U0; an
# infinite-strength shock reflects from the axis and travels outward at
# D = (gamma-1)/2 * U0 (Noh, JCP 72 (1987) 78, m=1 relations — see
# test/common/pycheck/q2d_exact1d.py::noh_exact, cross-checked there against
# the classic gamma=5/3 values D = 1/3, rho_post = 16, p_post = 16/3).
#
# Q2D axisymmetric convention (docs/theory/axisymmetric-derivation.md +
# Lib_Source_Axisymmetry.f90 Precompute_Axisymmetric_Geometry): the SECOND
# grid coordinate y is the radius r (cell radius r_ave from node%c(2)) and
# v = Prim(nv) is the radial velocity.  So the domain is thin in x (i) and
# resolves r in [0, R_MAX] along j; the imploding IC is v = -U0, u = 0.
#
# Geometry / resolution:
#   dr = 2.5e-3 (NJ = 1000 cells on r in [0, R_MAX = 2.5]): the t_final shock
#        radius D*t = 0.2 (air gamma ~ 1.4) spans 80 cells — a clean
#        post-shock plateau — while the single-thread bounds-checked solve
#        stays in the ~10 s range;
#   NI = 4 cells in x with dx = dr (LX = NI*dr): the exact solution is
#        x-invariant, so a few cells suffice and they double as the negative
#        control for spurious x-fluxes (check.py x-uniformity gate).
#
# Boundary conditions (bc.txt, 5-int records b,i,j,f,type — IO_BC.f90):
#   faces 1/2 (x ends) : 300 symmetry — the exact solution has u = 0 and
#        d/dx = 0, so mirror ghosts reproduce the interior state exactly and
#        the x-fluxes cancel identically (preserves x-uniformity; the 1D
#        cases use the same trick on their transverse faces);
#   face 3 (j=1, r=0)  : 300 symmetry — the axis; VERIFIED against
#        chang-nozzle (the reference axi case), which marks its axis row
#        with 300 (BC 200 is the ROTATIONAL-PERIODIC code, not an axis BC:
#        Mod_BC_Fluxes.f90 'case (200) ! periodic').  The 1/r geometric
#        source handles the axis via r_inv = 1/max(r_ave, eps);
#   face 4 (j=NJ, r=R_MAX) : 400 extrapolation — mirror of the 1D Noh87
#        outer boundary, PLUS a buffer band so its error cannot reach the
#        validated region.  The outer flow is a pressureless supersonic
#        INFLOW (u_r = -U0, a ~ 0): all characteristics enter the domain,
#        so the zero-gradient ghost mis-states the incoming mass flux of
#        the converging profile rho = rho0 (1 + U0 t/r) at an O(1) RATE
#        (ghost error (dr/2)|drho/dr| divided by the cell size dr gives
#        d(rho_err)/dt ~ U0 |drho/dr|/2, independent of dr — measured
#        rho = 2.72 vs exact 2.00 at the boundary cell of an R_MAX = 1 trial).
#        No time-dependent supersonic-inflow BC exists in the dispatcher
#        (405 imposes a CONSTANT Mach/T/p — wrong by the same order), so the
#        boundary error is instead OUTRUN: it advects inward at the pre-shock
#        signal speed |u_r| + a ~ U0, reaching r = R_MAX - U0 t_final = 1.5
#        at t_final, and check.py gates only r <= R_VALID = 1.0 (0.5 margin).
#
# Usage:  python3 gen_case.py [NJ] [target_dir]
#         defaults: NJ = 1000 (the ctest mesh), target_dir = this case dir.

import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
NJ   = int(sys.argv[1]) if len(sys.argv) > 1 else 1000  # cells along r
TGT  = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else HERE
NI   = 4                        # thin x-strip (x-uniformity negative control)

# --- Noh (1987) cylindrical implosion IC (mirror the 1D Noh87 choices) ---
RHO0 = 1.0                      # uniform initial density
U0   = 1.0                      # implosion speed: u_r = -U0 everywhere
# Near-zero pressure floor.  The 1D Noh87 case uses 1e-15, which works there
# because the planar pre-shock state is UNIFORM (zero truncation error);
# in cylindrical geometry the pre-shock compression rho = rho0 (1 + U0 t/r)
# is non-uniform and the residual truncation/roundoff noise in p is at the
# ~1e-14 scale (measured: with P0 = 1e-15 the run aborts at t ~ 0.29 with
# p = -3e-27 near the outer boundary).  1e-9 keeps ~1e5 headroom over that
# noise while remaining dynamically negligible: pre-shock sound speed
# a0 = sqrt(g P0/RHO0) ~ 4e-5 << U0, and P0/p_post ~ 1.4e-10.
P0   = 1.0e-9
R_MAX = 2.5                     # radial extent (validated region: r <= 1.0)

DR = R_MAX / NJ                 # radial spacing
LX = NI * DR                    # x extent chosen so dx = dr (isotropic cells)


def write_case(target, nj):
    inp = os.path.join(target, "INPUT")
    msh = os.path.join(target, "MESH")
    os.makedirs(inp, exist_ok=True)
    os.makedirs(msh, exist_ok=True)
    dr = R_MAX / nj
    dx = LX / NI
    xn = [i * dx for i in range(NI + 1)]      # x nodes
    rn = [j * dr for j in range(nj + 1)]      # r (= y) nodes, axis at r = 0
    fm = lambda v: "%20.15E" % v

    # MESH/mesh.tec (BLOCK; x,y NODAL; i fastest — chang-nozzle header format)
    with open(os.path.join(msh, "mesh.tec"), "w") as f:
        f.write(' VARIABLES ="x" "y"\n')
        f.write(' ZONE  T = Block1, I=%d, J=%d, K=1, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-2]=NODAL)\n' % (NI + 1, nj + 1))
        for j in range(nj + 1):
            for i in range(NI + 1):
                f.write(fm(xn[i]) + "\n")
        for j in range(nj + 1):
            for i in range(NI + 1):
                f.write(fm(rn[j]) + "\n")

    # INPUT/ic.tec (BLOCK; x,y NODAL; rho1,u,v,p CELLCENTERED; i fastest)
    with open(os.path.join(inp, "ic.tec"), "w") as f:
        f.write(' VARIABLES ="x" "y"  "rho1" "u" "v" "p"\n')
        f.write(' ZONE  T = B1-IG, I=%d, J=%d, K=1, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-2]=NODAL,[3-6]=CELLCENTERED), '
                'SOLUTIONTIME=0.000000000000000E+000\n' % (NI + 1, nj + 1))
        for j in range(nj + 1):
            for i in range(NI + 1):
                f.write(fm(xn[i]) + "\n")
        for j in range(nj + 1):
            for i in range(NI + 1):
                f.write(fm(rn[j]) + "\n")
        # cell blocks (NI*nj cells each): rho, u (axial), v (radial), p
        ncell = NI * nj
        for value in (RHO0, 0.0, -U0, P0):
            for _ in range(ncell):
                f.write(fm(value) + "\n")

    # INPUT/bc.txt : 5-int records (b, i, j, f, type) — see header comment.
    with open(os.path.join(inp, "bc.txt"), "w") as f:
        def rec(b, i, j, fc, code):
            f.write("%8d%8d%8d%8d%8d\n" % (b, i, j, fc, code))
        for j in range(1, nj + 1):
            rec(1, 1, j, 1, 300)          # x-min: symmetry (x-invariance)
        for j in range(1, nj + 1):
            rec(1, NI, j, 2, 300)         # x-max: symmetry (x-invariance)
        for i in range(1, NI + 1):
            rec(1, i, 1, 3, 300)          # r = 0: axis (chang-nozzle mirror)
        for i in range(1, NI + 1):
            rec(1, i, nj, 4, 400)         # r = R_MAX: extrapolation (Noh87 mirror)


if __name__ == "__main__":
    write_case(TGT, NJ)
    print("noh-cylindrical: NI=%d NJ=%d, x in [0,%g], r in [0,%g] -> %s"
          % (NI, NJ, LX, R_MAX, TGT))
    print("  IC: rho=%g, u_r=%g, p=%g (Noh m=1); shock D*t = (g-1)/2*t at t_final"
          % (RHO0, -U0, P0))
