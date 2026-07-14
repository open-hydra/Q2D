#!/usr/bin/env python3
"""Generate the q2D AREA-VARIATION turbulent planar-nozzle inputs (kw = Menter SST
and rsm = SSG/LRR-omega) from the committed-gas planar-nozzle C-D area-variation
case.

PURPOSE (SSGLRR area-variation V&V, Phase C — the *running* sibling of the 2D-axi
chang-nozzle-turbulent case): a running q2D area-variation RSM case that exercises
every area-variation-specific RSM path (theory q2d-rans-closures.md §6, Appendix
A.5):
  - P_33 = +2 rho R33 (V.G)          depth-compression production (§6.2 / A.5),
                                     GRADIENT-FREE -> the running-case evidence metric
  - convective  +rho R_ij (V.G)      Source_Area_Eul_Blk generic post-pressure line
  - diffusive   -F^(ij).G            SSGLRR_Source_Diffusive_Q2D, Source_Area_Visc_Blk
  - momentum coupling tau_lam(mu_l) - rho R_ij, Definition A tau_33 = lamDiv
                                     - 2 mu_l (V.G) - rho R33 (Stress_Tensor_RSM_Std, §6.8)
  - slip contour (BC 300): both top/bottom faces exercise the RSM symmetry
    diffusive-flux block of BC_Symmetry_Visc with the area-var V.G face average.

DESIGN (no invented numbers — every value derives from the planar-nozzle base
and the committed coarse-kw anchors, identical rule to the axi gen_inputs.py):
  - geometry/mesh/area-law/IC/numerics: VERBATIM the planar-nozzle base case
    (C-D nozzle, gamma=1.4 air, P0=1e6 Pa, T0=300 K, throat XT=0.40, AR_in=2.0,
    AR_ex=2.5; BC face1=401 inflow, face2=400 extrapolation, face3/4=300 slip;
    block1_area.dat = h(x) area law) + equations=navier-stokes, vnn=1.0,
    turbulence=SST|SSGLRR (the ONLY diff between the two inputs is the
    turbulence keyword).
  - The contour stays SLIP (300): the inviscid mesh has no wall clustering, so
    the case is a turbulent CORE flow (no boundary layer). With no viscous wall
    the wall distance is yn=1e8 -> F1=0 -> free-shear (SSG/eps-branch)
    coefficients, the correct model regime for a wall-free flow.
  - turbulence inflow/IC: the committed flat-plate `coarse-kw` freestream
    levels (k=2.25e-3 m^2/s^2, omega=6250 1/s at T=555 K, M=0.2) expressed
    NON-DIMENSIONALLY (k/U^2, mut/mu) and rescaled to the planar-nozzle inflow
    state (x=0, AR=2.0, computed from the analytic isentropic area-law):
        k_in     = (k/U^2)_kw * U_in^2          [same turbulence intensity]
        omega_in = rho_in k_in / ((mut/mu)_kw * mu)   [same mut/mu ratio]
    with mu = 1.1858685985e-5 (the constant-viscosity air transport.dat).
  - RSM initialization: ISOTROPIC from k (standard RSM init, identical rule to
    the axi/flat-plate gen_rsm_inputs.py): R11=R22=R33=(2/3)k, R12=0, omega same.

Column conventions (read in-session, identical to the axi case):
  ic.tec VARIABLES kappa/omega for nRANS=2, ruu/rvv/rww/ruv/omega for nRANS=5;
  BC 401 property line trailing tokens ci(1:nrans).
"""
import re
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "planar-nozzle"          # the area-variation base case (user drift)

MU_LAM = 1.1858685985e-5          # constant-viscosity air transport.dat
GAMMA, RGAS = 1.4, 287.0          # planar-nozzle gas (cp=1004.5, M=28.970418)
P0, T0 = 1.0e6, 300.0             # planar-nozzle stagnation state (gen_case.py)
AT, AIN = 1.0, 2.0                # throat / inlet area ratio (gen_case.py)

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
    return KW_K / U**2, rho * KW_K / (KW_OMEGA * MU_LAM)


def AR_of_M(M):
    g = GAMMA
    return (1.0 / M) * ((2.0 / (g + 1)) * (1.0 + 0.5 * (g - 1) * M * M)) ** ((g + 1) / (2 * (g - 1)))


def M_of_AR_sub(ar):
    lo, hi = 1e-7, 1.0 - 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if AR_of_M(mid) - ar > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def inflow_state():
    """Analytic isentropic inflow (x=0, AR=AIN/AT) — the IC the base ic.tec carries."""
    M = M_of_AR_sub(AIN / AT)
    tf = 1.0 + 0.5 * (GAMMA - 1) * M * M
    p = P0 * tf ** (-GAMMA / (GAMMA - 1))
    T = T0 / tf
    rho = p / (RGAS * T)
    U = M * math.sqrt(GAMMA * RGAS * T)
    return rho, U


def read_ic():
    src = (BASE / "INPUT" / "ic.tec").read_text().splitlines()
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
        chunk = data[o:o + n]
        o += n
        return chunk

    x = take(n_nodal); y = take(n_nodal)
    rho = take(n_cc); u = take(n_cc); v = take(n_cc); p = take(n_cc)
    return hdr_zone, ni, nj, x, y, rho, u, v, p


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
    """Append the nRANS trailing tokens to every BC-401 property line; all other
    records (300/400, no property line) pass through verbatim."""
    lines = (BASE / "INPUT" / "bc.txt").read_text().splitlines()
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


# constant-viscosity air transport (mirror the axi turbulent case values); the
# table T-range MUST match the planar-nozzle thermo.dat range (1..NT_THERMO) or
# the solver rejects it ("Transport data temperature range does not match thermo
# data range"). The base thermo.dat has NT=10000, so the transport table is
# generated at NT=10000 too (the axi case happened to share NT=5000 thermo).
MU_CONST, K_CONST = 1.1858686e-5, 0.016544513989


def gen_transport(dst, n_thermo):
    rows = ["TITLE = \"Transport Properties\"",
            "VARIABLES = \"Temperature\", \"Viscosity\", \"Conductivity\"",
            "ZONE T=\"Air\"", f"I={n_thermo}, F=POINT"]
    for i in range(1, n_thermo + 1):
        rows.append(f"{float(i):.1f} {MU_CONST:.12f} {K_CONST:.12f} ")
    (dst / "transport.dat").write_text("\n".join(rows) + "\n")


def thermo_n(text):
    import re as _re
    m = _re.search(r"I=\s*(\d+)", text)
    return int(m.group(1))


def copy_aux(variant):
    """area law + thermo/phase + matching transport (self-contained: the new case
    does not reference the drift base at solve time)."""
    dst = HERE / variant / "INPUT"
    dst.mkdir(parents=True, exist_ok=True)
    thermo_text = (BASE / "INPUT" / "thermo.dat").read_text()
    for f in ("block1_area.dat", "thermo.dat", "phase.txt"):
        (dst / f).write_text((BASE / "INPUT" / f).read_text())
    gen_transport(dst, thermo_n(thermo_text))   # range-matched to thermo.dat
    (HERE / variant / "OUTPUT").mkdir(parents=True, exist_ok=True)
    (HERE / variant / "bin").mkdir(parents=True, exist_ok=True)


def gen_ini(variant):
    turb = "SST" if variant == "kw" else "SSGLRR"
    ini = (
        "[Q2D-Parameters]\niter-threshold = 50000\n\n"
        "[Q2D-Physics]\nequations  = navier-stokes\n"
        f"turbulence = {turb}\n\n"
        "[Q2D-Numerics]\ncfl                   = 0.5\ncfl-rise-threshold    = 50\n"
        "vnn                   = 1.0\ntime-scheme           = RK3\ntime-accurate         = false\n"
        "irs-beta              = 0.5\nintegration-variables = cons\n"
        "space-reconstruction  = MUSCL\nflux-limiter          = vanleer\n"
        "riemann-solver        = HLLC\n\n"
        "[Q2D-IO]\nsol-diter   = 5000\nshell-diter = 500\n\n"
        "[Q2D-Sources]\nactive = area_variation\n"
    )
    (HERE / variant / "input.ini").write_text(ini)
    print(f"{variant}/input.ini: navier-stokes + {turb} + area_variation")


if __name__ == "__main__":
    k_over_U2, mut_over_mu = kw_anchors()
    rho_in, U_in = inflow_state()
    k_in = k_over_U2 * U_in**2
    omega_in = rho_in * k_in / (mut_over_mu * MU_LAM)
    print(f"anchors (coarse-kw): k/U^2={k_over_U2:.6E}  mut/mu={mut_over_mu:.6E}")
    print(f"planar-nozzle inflow (x=0, AR={AIN/AT:.1f}): rho={rho_in:.6f} kg/m^3, "
          f"|V|={U_in:.4f} m/s")
    print(f"derived: k_in={k_in:.6E} m^2/s^2  omega_in={omega_in:.6E} 1/s  "
          f"(R_iso={TWO_THIRDS*k_in:.6E})")
    for variant in ("kw", "rsm"):
        copy_aux(variant)
        gen_ic(variant, k_in, omega_in)
        gen_bc(variant, k_in, omega_in)
        gen_ini(variant)
