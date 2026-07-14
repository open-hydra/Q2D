#!/usr/bin/env python3
"""Generate the axisymmetric turbulent chang-nozzle inputs (kw = Menter SST and
rsm = SSG/LRR-omega) from the committed Euler chang-nozzle case.

PURPOSE (SSGLRR axi V&V, Phase C): a *running* 2D-axisymmetric RSM case that
exercises every axi-specific RSM path (theory q2d-rans-closures.md §6):
  - P_theta_theta = -2 rho R33 u_r/r  (azimuthal-stretching production, §6.2)
  - convective  -rho u_r R_ij / r     (Source_Axi_Eul_Blk, §6.1)
  - diffusive   +F^(ij)_r / r         (SSGLRR_Source_Diffusive_Q2D, §6.4)
  - axis r->0 (face 3, y=1e-8): RSM symmetry BC + de l'Hopital branch (§6.7)
  - momentum coupling tau_lam(mu_l) - rho R_ij (Stress_*_RSM_Axi, §6.8)

DESIGN (no invented numbers — every value derives from committed sources):
  - geometry/mesh/IC/numerics: VERBATIM the committed `test/2D/euler/chang-nozzle`
    (NASA-Chang converging-diverging nozzle, T0=555.56 K, p0=1.0342e6 Pa,
    BC: face1=401 inflow, face2=400 extrapolation, face3=300 axis,
    face4=300 slip contour).
  - The contour stays SLIP (300): the committed Euler mesh has no wall
    clustering (first cell-centre y+ would be O(10^3), the boundary layer is
    not resolvable), so the case is a turbulent CORE flow. With no viscous
    wall the wall distance is yn=1e8 -> F1=0 -> the model runs on its
    free-shear (SSG / k-eps-branch) coefficient set, which is the correct
    model behaviour for a wall-free flow. Both BC-300 faces (axis + contour)
    exercise the RSM symmetry diffusive-flux block of BC_Symmetry_Visc.
  - turbulence inflow/IC: the committed flat-plate `coarse-kw` freestream
    levels (k=2.25e-3 m^2/s^2, omega=6250 1/s at T=555 K, M=0.2, p~1 atm,
    input.ini ICB-Block1) expressed NON-DIMENSIONALLY and rescaled to the
    chang inflow state (computed from the committed ic.tec inflow column):
        k_in     = (k/U^2)_kw * U_in^2          [same turbulence intensity]
        omega_in = rho_in k_in / ((mut/mu)_kw * mu)   [same mut/mu ratio]
    with mu = 1.1858685985e-5 (common/Air-viscous transport.dat, constant).
  - RSM initialization: ISOTROPIC from k (standard RSM init, identical rule
    to flat-plate gen_rsm_inputs.py): R11=R22=R33=(2/3)k, R12=0, omega same.

Column conventions read in-session from src/lib/io/IO_BC.f90 (BC 401 property
line: T0, p0, alpha, beta, rel_fac, ci(species), ci(1:nrans)) and
src/lib/io/IO_Solution.f90 (ic.tec VARIABLES naming: kappa/omega for nRANS=2,
ruu/rvv/rww/ruv/omega for nRANS=5).
"""
import re
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent.parent / "euler" / "chang-nozzle"

MU_LAM = 1.1858685985e-5          # common/Air-viscous transport.dat (constant)
GAMMA, RGAS = 1.4, 287.0          # air (phase/thermo: cp=1004.5, M=28.97)

# committed flat-plate coarse-kw freestream (input.ini ICB-Block1 + bc.txt 403)
KW_T, KW_M, KW_P0 = 555.0, 0.2, 102828.0
KW_K, KW_OMEGA = 2.25e-3, 6250.0

TWO_THIRDS = 2.0 / 3.0


def kw_anchors():
    """Non-dimensional turbulence anchors of the committed coarse-kw case."""
    a = math.sqrt(GAMMA * RGAS * KW_T)
    U = KW_M * a
    p = KW_P0 / (1.0 + 0.2 * KW_M**2) ** 3.5
    rho = p / (RGAS * KW_T)
    k_over_U2 = KW_K / U**2
    mut_over_mu = rho * KW_K / (KW_OMEGA * MU_LAM)
    return k_over_U2, mut_over_mu


def read_ic():
    src = (SRC / "INPUT" / "ic.tec").read_text().splitlines()
    hdr_vars, hdr_zone = src[0], src[1]
    m = re.search(r"I\s*=\s*(\d+).*?J\s*=\s*(\d+)", hdr_zone)
    ni, nj = int(m.group(1)), int(m.group(2))
    n_nodal, n_cc = ni * nj, (ni - 1) * (nj - 1)
    data = src[2:]
    assert len(data) == 2 * n_nodal + 4 * n_cc, \
        f"ic.tec data lines {len(data)} != {2*n_nodal + 4*n_cc}"
    o = 0
    def take(n):
        nonlocal o
        chunk = data[o:o + n]; o += n; return chunk
    x = take(n_nodal); y = take(n_nodal)
    rho = take(n_cc); u = take(n_cc); v = take(n_cc); p = take(n_cc)
    return hdr_zone, ni, nj, x, y, rho, u, v, p


def inflow_state(ni, nj, rho, u, v):
    """Mean inflow-column (i=1 cells) state from the committed IC."""
    nci = ni - 1
    rho_in, U_in, n = 0.0, 0.0, 0
    for j in range(nj - 1):
        idx = j * nci          # i=1 cell of row j (BLOCK data, i fastest)
        r = float(rho[idx]); uu = float(u[idx]); vv = float(v[idx])
        rho_in += r; U_in += math.hypot(uu, vv); n += 1
    return rho_in / n, U_in / n


def gen_ic(variant, k_in, omega_in):
    hdr_zone, ni, nj, x, y, rho, u, v, p = read_ic()
    n_cc = (ni - 1) * (nj - 1)
    kap = [f" {k_in:.15E}"] * n_cc
    ome = [f" {omega_in:.15E}"] * n_cc
    if variant == "kw":
        new_vars = ' VARIABLES ="x" "y"  "rho1" "u" "v" "p"  "kappa"  "omega"'
        new_zone = hdr_zone.replace("[3-6]=CELLCENTERED", "[3-8]=CELLCENTERED")
        blocks = (x, y, rho, u, v, p, kap, ome)
    else:
        rss = [f" {TWO_THIRDS * k_in:.15E}"] * n_cc
        zero = [f" {0.0:.15E}"] * n_cc
        new_vars = (' VARIABLES ="x" "y"  "rho1" "u" "v" "p"'
                    '  "ruu"  "rvv"  "rww"  "ruv"  "omega"')
        new_zone = hdr_zone.replace("[3-6]=CELLCENTERED", "[3-11]=CELLCENTERED")
        blocks = (x, y, rho, u, v, p, rss, rss, rss, zero, ome)
    out = [new_vars, new_zone]
    for b in blocks:
        out.extend(b)
    dst = HERE / variant / "INPUT" / "ic.tec"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(out) + "\n")
    print(f"{variant}/INPUT/ic.tec: {ni}x{nj}, {len(blocks)} vars")


def gen_bc(variant, k_in, omega_in):
    """Append the nRANS trailing tokens to every BC-401 property line; all
    other records (300/400, no property line) pass through verbatim."""
    lines = (SRC / "INPUT" / "bc.txt").read_text().splitlines()
    if variant == "kw":
        rans = [f"{k_in:.6E}", f"{omega_in:.6E}"]
    else:
        rss = TWO_THIRDS * k_in
        rans = [f"{rss:.6E}"] * 3 + [f"{0.0:.6E}", f"{omega_in:.6E}"]
    out, transform_next = [], False
    idx_re = re.compile(r"^\s*\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s*$")
    for ln in lines:
        m = idx_re.match(ln)
        if m:
            out.append(ln)
            transform_next = int(m.group(1)) == 401
            continue
        if transform_next:
            toks = [t for t in re.split(r"[,\s]+", ln.strip()) if t]
            out.append("    " + ",    ".join(toks + rans) + ",")
            transform_next = False
        else:
            out.append(ln)
    dst = HERE / variant / "INPUT" / "bc.txt"
    dst.write_text("\n".join(out) + "\n")
    print(f"{variant}/INPUT/bc.txt: {len(out)} lines")


def copy_mesh(variant):
    dst = HERE / variant / "MESH"
    dst.mkdir(parents=True, exist_ok=True)
    for f in ("mesh.tec", "profile.dat"):
        (dst / f).write_text((SRC / "MESH" / f).read_text())
    print(f"{variant}/MESH: mesh.tec + profile.dat copied")


if __name__ == "__main__":
    k_over_U2, mut_over_mu = kw_anchors()
    _, ni, nj, x, y, rho, u, v, p = read_ic()
    rho_in, U_in = inflow_state(ni, nj, rho, u, v)
    k_in = k_over_U2 * U_in**2
    omega_in = rho_in * k_in / (mut_over_mu * MU_LAM)
    print(f"anchors (coarse-kw): k/U^2={k_over_U2:.6E}  mut/mu={mut_over_mu:.6E}")
    print(f"chang inflow (from committed ic.tec): rho={rho_in:.6f} kg/m^3, "
          f"|V|={U_in:.4f} m/s")
    print(f"derived: k_in={k_in:.6E} m^2/s^2  omega_in={omega_in:.6E} 1/s  "
          f"(R_iso={TWO_THIRDS*k_in:.6E})")
    for variant in ("kw", "rsm"):
        gen_ic(variant, k_in, omega_in)
        gen_bc(variant, k_in, omega_in)
        copy_mesh(variant)
