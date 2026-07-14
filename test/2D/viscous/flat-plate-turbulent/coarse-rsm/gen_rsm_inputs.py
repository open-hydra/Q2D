#!/usr/bin/env python3
"""Generate the SSG/LRR-omega RSM inputs (ic.tec, bc.txt) for the NASA TMR
flat-plate from the existing coarse-kw (Menter SST) inputs.

The RSM transports 5 variables at index nt (theory q2d-rans-closures.md §6.1):
    rho*R11, rho*R22, rho*R33, rho*R12, rho*omega
(the two cross-stresses R13,R23 vanish by the regime mirror symmetry, 7->5).

INITIALIZATION (user-decided, non-arbitrary): the Reynolds stress is initialized
ISOTROPIC from the existing kw turbulent kinetic energy k:
    R11 = R22 = R33 = (2/3) k ,   R12 = 0 ,   omega = omega_kw
This is the standard RSM initialization (isotropic stress from the freestream k)
and derives 1:1 from the committed coarse-kw IC/BC — no invented numbers
(feedback_no_random_params_check_mose). The kw freestream is k=2.25e-3, omega=6250
(input.ini ICB-Block1), so R11=R22=R33=1.5e-3, R12=0.

It transforms:
  - INPUT/ic.tec : Tecplot BLOCK; the kappa block -> ruu,rvv,rww (=2/3 kappa) +
                   ruv (=0) blocks; omega block unchanged; header VARIABLES +
                   VARLOCATION updated (8 -> 11 variables).
  - INPUT/bc.txt : every BC-403/407 data line carries the nRANS trailing values
                   (kappa, omega for SST). They are replaced by the 5 RSM values
                   (2/3 kappa, 2/3 kappa, 2/3 kappa, 0, omega). All other lines
                   (index lines, 300/301 wall lines) pass through unchanged.

Source of truth for the column conventions (read in-session):
  - src/lib/io/IO_BC.f90 (BC 403/407: trailing ci(1:nrans) = RANS, mass-specific)
  - src/lib/io/IO_Solution.f90 (ic.tec VARIABLES/VARLOCATION block layout)
"""
import re
from pathlib import Path

KW  = Path(__file__).resolve().parent.parent / "coarse-kw"
OUT = Path(__file__).resolve().parent

TWO_THIRDS = 2.0 / 3.0


def gen_ic():
    # Preserve the EXACT original line formatting (leading-space ' %.15E') of the
    # kept blocks — the ORION Tecplot reader is whitespace/format sensitive, so we
    # slice the original VERBATIM lines by block rather than re-stringifying them.
    src = (KW / "INPUT" / "ic.tec").read_text().splitlines()
    hdr_vars, hdr_zone = src[0], src[1]
    data = src[2:]

    m = re.search(r"I\s*=\s*(\d+).*?J\s*=\s*(\d+)", hdr_zone)
    ni, nj = int(m.group(1)), int(m.group(2))
    n_nodal = ni * nj
    n_cc = (ni - 1) * (nj - 1)
    assert len(data) == 2 * n_nodal + 6 * n_cc, \
        f"ic.tec data lines {len(data)} != {2*n_nodal + 6*n_cc}"

    o = 0
    def take(n):
        nonlocal o
        chunk = data[o:o + n]; o += n; return chunk
    x  = take(n_nodal); y = take(n_nodal)         # verbatim original lines
    rho = take(n_cc); u = take(n_cc); v = take(n_cc); p = take(n_cc)
    kap = take(n_cc); ome = take(n_cc)

    # Isotropic stress from k: R11=R22=R33=2/3 k, R12=0. Match the kept-block
    # format exactly: a literal leading space + %.15E (as in the kw ic.tec).
    ruu = [f" {TWO_THIRDS * float(s):.15E}" for s in kap]
    ruv = [f" {0.0:.15E}" for _ in kap]

    new_vars = ' VARIABLES ="x"  "y"  "rho1"  "u"  "v"  "p"  "ruu"  "rvv"  "rww"  "ruv"  "omega"'
    new_zone = hdr_zone.replace("[3-8]=CELLCENTERED", "[3-11]=CELLCENTERED")

    out = [new_vars, new_zone]
    for block in (x, y, rho, u, v, p, ruu, ruu, ruu, ruv, ome):
        out.extend(block)
    (OUT / "INPUT" / "ic.tec").write_text("\n".join(out) + "\n")
    print(f"ic.tec: {ni}x{nj}  nodal={n_nodal} cc={n_cc}  11 vars (format-preserved)")


def gen_bc():
    lines = (KW / "INPUT" / "bc.txt").read_text().splitlines()
    out, transform_next = [], False
    idx_re = re.compile(r"^\s*\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s*$")
    for ln in lines:
        m = idx_re.match(ln)
        if m:
            code = int(m.group(1))
            out.append(ln)
            transform_next = code in (401, 402, 403, 404, 405, 407, 420)
            continue
        if transform_next:
            # data line: trailing nRANS=2 tokens (kappa, omega) -> 5 RSM tokens
            toks = [t for t in re.split(r"[,\s]+", ln.strip()) if t]
            head, kap, ome = toks[:-2], float(toks[-2]), float(toks[-1])
            rss = TWO_THIRDS * kap
            new = head + [f"{rss:.6E}", f"{rss:.6E}", f"{rss:.6E}", f"{0.0:.6E}", f"{ome:.6E}"]
            out.append("    " + ",    ".join(new) + ",")
            transform_next = False
        else:
            out.append(ln)
    (OUT / "INPUT" / "bc.txt").write_text("\n".join(out) + "\n")
    print(f"bc.txt: {len(out)} lines written")


if __name__ == "__main__":
    gen_ic()
    gen_bc()
