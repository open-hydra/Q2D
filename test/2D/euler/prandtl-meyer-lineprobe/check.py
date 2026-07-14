"""prandtl-meyer-lineprobe regression check (Q2D-extra peculiarità #5).

Exercises the Q2D-only line-probe feature (``Q2D_IO_LineProbes``; MOSE has no
line-probes) end-to-end through the registry-driven path migrated in commit
7d47234 (Scan_Ini count -> Register_LineProbes bind -> Setup_LineProbes read ->
Write_LineProbe_Tecplot output). The case is the prandtl-meyer planar expansion
fan with a single j-fixed line-probe ``[midrow]`` (``index-position = 1 2``,
``fixed-index = 25``) sampling the cell row j=25 along i.

Reference (no invented thresholds, self-consistency): a correct line-probe MUST
reproduce, at the sampled grid line, the same primitive field the main solution
writer (``OUTPUT/field.tec``) writes there. Both write the identical
``domain%blk(b)%P(:, i, j)`` slice (IO_Solution.f90:261 vs
IO_LineProbes.f90:127), so the line-probe values at row j=25 equal the field
cell values at (i, 25). The line-probe writes per ``diter`` (a converged-time
series); we compare its LAST zone against the field. The residual offset between
the line-probe's last write and the field write is bounded by the converged
density residual (< 1e-8 here, the field being steady between the two iters), so
the acceptance is |lp - field| < atol + rtol*|field| with atol=1e-6, rtol=1e-3
(a real extraction bug -- wrong row / transposed index / wrong variable map --
would be off by O(1)..O(100%), far outside this band).

Three further checks confirm the output makes sense *as a line-probe* (not merely
that values match the field):
  (struct)   the file/header/geometry are written (registry binding + writer ran);
  (geometry) the written degenerate strip nodes ARE the mesh node rows the
             line-probe interpretation prescribes -- lower row == mesh j=FIXED_J-1,
             upper row == mesh j=FIXED_J, bit-exact, read from the case's own
             ic.tec mesh (Setup_Nodes dir=2: node(ip,idx-1)/node(ip,idx)); i.e. the
             (x,y) really trace that geometric line in the domain;
  (physical) along the sampled line the Prandtl-Meyer expansion fan accelerates
             the flow -- Mach up (anchored to inflow M=2 at the line start), p/rho/T
             down -- the qualitative invariant of what the line-probe samples.
"""
import os
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

LP_FILE = "OUTPUT/midrow.tec"
FIELD_FILE = "OUTPUT/field.tec"
FIXED_J = 25                       # the line-probe fixed-index (cell row)
LP_VARS = ["rho", "u", "v", "p", "T", "M"]
FIELD_VARS = ["rho", "u", "v", "p", "T", "g", "R"]   # field.tec column order
ATOL = 1.0e-6
RTOL = 1.0e-3
NEARZERO = 1.0e-6                   # below this magnitude, compare on atol only


def read_last_zone(path):
    """Return (ni, nj, nk, flat_token_list) for the LAST Tecplot zone."""
    with open(path) as f:
        lines = f.read().splitlines()
    zone_idx = [i for i, l in enumerate(lines) if "ZONE" in l]
    if not zone_idx:
        raise ValueError(f"no ZONE record in {path}")
    z = zone_idx[-1]
    zone = lines[z]
    m = re.search(r"I=\s*(\d+).*J=\s*(\d+).*K=\s*(\d+)", zone)
    ni, nj, nk = int(m.group(1)), int(m.group(2)), int(m.group(3))
    toks = [float(t) for ln in lines[z + 1:] for t in ln.split()]
    return ni, nj, nk, toks


def fail(msg):
    print(f"FAIL: {msg}")
    sys.exit(1)


# --- structural check: the line-probe file exists with the expected header ----
if not os.path.exists(LP_FILE):
    fail(f"{LP_FILE} not written -- line-probe registry binding/setup did not run")

with open(LP_FILE) as f:
    header = f.readline()
for v in ['"rho(1)"', '"u"', '"v"', '"p"', '"T"', '"M"']:
    if v not in header:
        fail(f"line-probe header missing variable {v}: {header.strip()}")

lp_ni, lp_nj, lp_nk, lp_toks = read_last_zone(LP_FILE)
if lp_nj != 2 or lp_nk != 1:
    fail(f"line-probe zone geometry unexpected: I={lp_ni} J={lp_nj} K={lp_nk} (expected J=2 K=1)")
Ni = lp_ni - 1                      # npoints (sampled cells along i)

# line-probe BLOCK layout: x then y, each = lower-row(0..Ni) ++ upper-row(0..Ni)
# = 2*(Ni+1) values (the degenerate j=0/j=1 strip bracketing the sampled cells);
# then each var = Ni cell values along the line (IO_LineProbes.f90:301-305).
nxy = 2 * (Ni + 1)
idx = 0
lp_x = lp_toks[idx:idx + nxy]; idx += nxy
lp_y = lp_toks[idx:idx + nxy]; idx += nxy
lp_xlo, lp_xup = lp_x[:Ni + 1], lp_x[Ni + 1:]
lp_ylo, lp_yup = lp_y[:Ni + 1], lp_y[Ni + 1:]
lp = {}
for name in LP_VARS:
    lp[name] = lp_toks[idx:idx + Ni]
    idx += Ni
if idx != len(lp_toks):
    fail(f"line-probe token count {len(lp_toks)} != expected {idx} (Ni={Ni}, {len(LP_VARS)} vars)")

print(f"(struct) line-probe OUTPUT/midrow.tec OK: {len(LP_VARS)} vars, Ni={Ni} sampled cells, J=2")

# --- geometry check: the strip nodes ARE the mesh node rows the line-probe means
# A j-fixed line-probe at fixed-index=FIXED_J writes a degenerate strip whose lower
# row = mesh node row j=FIXED_J-1 and upper row = mesh node row j=FIXED_J
# (Setup_Nodes dir=2: node(ip, idx-1) / node(ip, idx)). Confirm the written (x,y)
# really is that geometric line in the domain, read from the case's own ic.tec mesh.
with open("INPUT/ic.tec") as f:
    ic_lines = f.read().splitlines()
m = re.search(r"I=\s*(\d+).*J=\s*(\d+)", ic_lines[1])
ic_ni, ic_nj = int(m.group(1)), int(m.group(2))
ic_toks = [float(t) for ln in ic_lines[2:] for t in ln.split()]
ic_nnode = ic_ni * ic_nj
ic_X = ic_toks[:ic_nnode]
ic_Y = ic_toks[ic_nnode:2 * ic_nnode]


def node(arr, i, j):                # 0-based node indices; BLOCK = i fastest within j
    return arr[j * ic_ni + i]


err_lo = max(abs(lp_xlo[i] - node(ic_X, i, FIXED_J - 1)) +
             abs(lp_ylo[i] - node(ic_Y, i, FIXED_J - 1)) for i in range(Ni + 1))
err_up = max(abs(lp_xup[i] - node(ic_X, i, FIXED_J)) +
             abs(lp_yup[i] - node(ic_Y, i, FIXED_J)) for i in range(Ni + 1))
if err_lo > 1e-12 or err_up > 1e-12:
    fail(f"line-probe strip geometry mismatch vs mesh: "
         f"|lower-row - mesh j={FIXED_J-1}|={err_lo:.3e}, |upper-row - mesh j={FIXED_J}|={err_up:.3e}")
print(f"(geometry) strip nodes == mesh rows j={FIXED_J-1}/{FIXED_J} exactly "
      f"(|d_lower|={err_lo:.1e}, |d_upper|={err_up:.1e}); "
      f"line spans x=[{lp_xlo[0]:.3f},{lp_xlo[-1]:.3f}]")

# --- physical-sense check: along a Prandtl-Meyer expansion fan the flow accelerates
# (Mach up; p, rho, T down). Anchor: inflow Mach = 2 (bc.txt BC 405 M=2) at the line
# start. These are the qualitative invariants of the line-probe's *interpretation*
# (it samples the real flow along the line), not invented numeric thresholds.
M = lp["M"]
if not (1.9 < M[0] < 2.1):
    fail(f"line start Mach {M[0]:.3f} not at inflow M=2 (line-probe not anchored to the inflow edge)")
if not (M[-1] > M[0] + 1.0):
    fail(f"line end Mach {M[-1]:.3f} did not accelerate vs start {M[0]:.3f} (no expansion sampled)")
for name, arr in (("p", lp["p"]), ("rho", lp["rho"]), ("T", lp["T"])):
    if not (arr[-1] < arr[0]):
        fail(f"{name} did not drop along the expansion: start {arr[0]:.4e} end {arr[-1]:.4e}")
inc = sum(1 for i in range(1, Ni) if M[i] >= M[i - 1] - 1e-9)
print(f"(physical) Prandtl-Meyer expansion along line OK: M {M[0]:.3f}->{M[-1]:.3f} "
      f"(increasing {inc}/{Ni-1}), p {lp['p'][0]:.3f}->{lp['p'][-1]:.3f}, "
      f"T {lp['T'][0]:.3f}->{lp['T'][-1]:.3f} (all dropping)")

# --- self-consistency check: line-probe vs field at row j=25 -------------------
if not os.path.exists(FIELD_FILE):
    fail(f"{FIELD_FILE} not written")

f_ni, f_nj, f_nk, f_toks = read_last_zone(FIELD_FILE)
nnode = f_ni * f_nj
ncell = (f_ni - 1) * (f_nj - 1)
if (f_ni - 1) != Ni:
    fail(f"field i-cells {(f_ni-1)} != line-probe npoints {Ni}")

idx = 2 * nnode
fvars = {}
for name in FIELD_VARS:
    fvars[name] = f_toks[idx:idx + ncell]
    idx += ncell


def fcell(name, i, j):             # 1-based cell indices; BLOCK = i fastest within j
    return fvars[name][(j - 1) * (f_ni - 1) + (i - 1)]


worst = {}
for name in ["rho", "u", "v", "p", "T"]:
    me = 0.0
    ai = 0
    for i in range(1, Ni + 1):
        a = lp[name][i - 1]
        b = fcell(name, i, FIXED_J)
        diff = abs(a - b)
        tol = ATOL + RTOL * abs(b)
        if diff > tol:
            metric = diff - tol
            if metric > me:
                me = metric
                ai = i
    worst[name] = (me, ai)
    if me > 0.0:
        i = ai
        fail(f"{name}: line-probe {lp[name][i-1]:.8e} vs field {fcell(name,i,FIXED_J):.8e} "
             f"at i={i} exceeds atol+rtol*|field| (|diff|={abs(lp[name][i-1]-fcell(name,i,FIXED_J)):.3e})")

for name in ["rho", "u", "v", "p", "T"]:
    mx = max(abs(lp[name][i] - fcell(name, i + 1, FIXED_J)) for i in range(Ni))
    print(f"(self-consistency) {name}: max|lp-field@j={FIXED_J}| = {mx:.3e}  (atol {ATOL:.0e} + rtol {RTOL:.0e})")

print("Correct!")
