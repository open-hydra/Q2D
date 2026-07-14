"""hypersonic-cylinder 2D steady-state regression check (DOUBLE reference).

Validates the converged 2D plane Euler solution for hypersonic flow M=8.1
over a cylinder cross-section (bow shock at stagnation line). The check
applies TWO independent references (CLAUDE.md Legge #5 biunivocita' Q2D≡MOSE
+ user direction 2026-05-15 "doppia ref"):

  (1) Natural physical reference: normal-shock + isentropic-stagnation theory
      for inflow M1=8.1, gamma=1.4 (mirror MOSE verify.py probe at I=80 J=0).

  (2) MOSE bit-exact reference: Q2D field K=0 vs MOSE field K=0 cell-by-cell
      (Linf_rel < 1e-9, L2_rel < 1e-12 — RELEASE -O3 numerical-precision class).
      Reference cache at test/regression_baseline/MOSE_K0/hypersonic-cylinder/.

Geometric classification (rigorous, z-coord analysis of MOSE wedge mesh):
  2D PLANE — z=+-5e-4 const per K-layer, X/Y identical between K-layers,
             half-angle wedge tilt = 0 deg => no axisymmetric source.
  Q2D config: [Q2D-Sources] active=none (omitted from input.ini).

Tolerances per ``docs/vv/regression_tolerance.md`` Section 3.2:
  Sanity: min(rho), min(p), min(T) > 0 (hard, any violation = FAIL).
  Natural: total p0 within 5% of analytical; total T0 within 1% (conservative).
  Bit-exact MOSE: Linf_rel < 1e-9 (RELEASE-O3), L2_rel < 1e-12.
"""
import math
import os
import re
import subprocess
import sys
import tempfile
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break

try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "ORION package not found. Run 'conda activate base' or install ORION in the current Python environment."
    )


def patch_q2d_tec_to_3d(in_path):
    """Patch Q2D 2D-strip field.tec (K=1) to ORION-compatible 3D-wedge layout
    (K=2 with fake z=0). Mirror chang-nozzle/prandtl-meyer helper."""
    with open(in_path) as f:
        lines = f.readlines()
    header = lines[0]
    zone   = lines[1]
    m = re.search(r'I=\s*(\d+).*J=\s*(\d+).*K=\s*(\d+)', zone)
    ni, nj, nk = int(m.group(1)), int(m.group(2)), int(m.group(3))
    nnodes = ni * nj * nk
    new_header = header.replace('"y"', '"y" "z"', 1)
    m2 = re.search(r'VARLOCATION=\(\[(\d+)-(\d+)\]=NODAL,\[(\d+)-(\d+)\]=CELLCENTERED', zone)
    if m2:
        na, nb_, ca, cb = (int(m2.group(i)) for i in (1, 2, 3, 4))
        new_zone = zone.replace(
            f'[{na}-{nb_}]=NODAL,[{ca}-{cb}]=CELLCENTERED',
            f'[{na}-{nb_+1}]=NODAL,[{ca+1}-{cb+1}]=CELLCENTERED',
        )
    else:
        new_zone = zone
    z_block = ' 0.000000000000000E+000\n' * nnodes
    out = [new_header, new_zone]
    out.extend(lines[2 : 2 + nnodes])
    out.extend(lines[2 + nnodes : 2 + 2 * nnodes])
    out.append(z_block)
    out.extend(lines[2 + 2 * nnodes :])
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.tec', delete=False)
    tmp.writelines(out)
    tmp.close()
    return tmp.name


# ---------------------------------------------------------------------------
# (1) NATURAL REFERENCE: normal-shock + isentropic stagnation
# ---------------------------------------------------------------------------
def normal_shock(M1, gamma):
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0
    M1sq = M1 * M1
    M2sq = (M1sq + 2.0/gm1) / (2.0*gamma/gm1 * M1sq - 1.0)
    M2   = math.sqrt(M2sq)
    p2p1 = 1.0 + 2.0*gamma/gp1 * (M1sq - 1.0)
    r2r1 = gp1*M1sq / (gm1*M1sq + 2.0)
    t2t1 = p2p1 / r2r1
    return M2, p2p1, t2t1, r2r1


def total_pressure(p, M, gamma):
    return p * (1.0 + 0.5*(gamma - 1.0)*M*M) ** (gamma / (gamma - 1.0))


def total_temperature(T, M, gamma):
    return T * (1.0 + 0.5*(gamma - 1.0)*M*M)


M1    = 8.1
gamma = 1.4
p1    = 101325.0
T1    = 288.14

M2_ref, p2p1_ref, t2t1_ref, r2r1_ref = normal_shock(M1, gamma)
p2_ref = p2p1_ref * p1
T2_ref = t2t1_ref * T1
p0_ref = total_pressure(p2_ref, M2_ref, gamma)
T0_ref = total_temperature(T2_ref, M2_ref, gamma)

print(f"(1) Natural reference (normal shock M1={M1} gamma={gamma}):")
print(f"    M2={M2_ref:.6f}  p0_post={p0_ref:.4f} Pa  T0_post={T0_ref:.4f} K")

patched = patch_q2d_tec_to_3d('OUTPUT/field.tec')
try:
    [x_, y_, z_, var_, vnames] = read_TEC(patched)
finally:
    os.unlink(patched)

rho = var_[0][0]; u = var_[0][1]; v = var_[0][2]
p   = var_[0][3]; T = var_[0][4]; g = var_[0][5]; R = var_[0][6]

if rho.min() <= 0.0:
    print(f'Sanity FAIL: min(rho)={rho.min()} <= 0'); sys.exit(1)
if p.min() <= 0.0:
    print(f'Sanity FAIL: min(p)={p.min()} <= 0'); sys.exit(1)
if T.min() <= 0.0:
    print(f'Sanity FAIL: min(T)={T.min()} <= 0'); sys.exit(1)

I_PROBE = 80
J_PROBE = 0
p_pr = float(p[I_PROBE, J_PROBE, 0]); T_pr = float(T[I_PROBE, J_PROBE, 0])
u_pr = float(u[I_PROBE, J_PROBE, 0]); v_pr = float(v[I_PROBE, J_PROBE, 0])
g_pr = float(g[I_PROBE, J_PROBE, 0]); R_pr = float(R[I_PROBE, J_PROBE, 0])
speed = math.sqrt(u_pr**2 + v_pr**2)
a_pr  = math.sqrt(g_pr * R_pr * T_pr)
M_pr  = speed / a_pr
p0_pr = total_pressure(p_pr, M_pr, g_pr)
T0_pr = total_temperature(T_pr, M_pr, g_pr)
err_p0 = 100.0 * (p0_pr - p0_ref) / p0_ref
err_T0 = 100.0 * (T0_pr - T0_ref) / T0_ref
print(f"    Probe (I={I_PROBE} J={J_PROBE}): p={p_pr:.4e}  T={T_pr:.4f}  M={M_pr:.6f}")
print(f"    Q2D: p0={p0_pr:.4f} Pa err={err_p0:+.4f}%  T0={T0_pr:.4f} K err={err_T0:+.4f}%")

thr_p0_natural = 5.0
thr_T0_natural = 1.0
if abs(err_p0) > thr_p0_natural:
    print(f'(1) FAIL: p0 natural mismatch {err_p0:+.4f}% > {thr_p0_natural}%'); sys.exit(1)
if abs(err_T0) > thr_T0_natural:
    print(f'(1) FAIL: T0 natural mismatch {err_T0:+.4f}% > {thr_T0_natural}%'); sys.exit(1)
print('(1) Natural reference: OK')

# ---------------------------------------------------------------------------
# (2) MOSE BIT-EXACT REFERENCE: Q2D field K=0 cell-by-cell vs MOSE field K=0
# ---------------------------------------------------------------------------
mose_ref = root.parents[3] / 'test' / 'regression_baseline' / 'MOSE_K0' / 'hypersonic-cylinder' / 'field.tec'
diff_tool = root.parents[3] / 'tools' / 'diff_mose_q2d_field.py'
if not mose_ref.exists():
    print(f'(2) WARN: MOSE bit-exact reference not found at {mose_ref}; skipping bit-exact check')
elif not diff_tool.exists():
    print(f'(2) WARN: diff utility not found at {diff_tool}; skipping bit-exact check')
else:
    print(f"(2) MOSE bit-exact reference vs Q2D field (Linf_rel < 1e-9, L2_rel < 1e-12):")
    res = subprocess.run([sys.executable, str(diff_tool), str(mose_ref), 'OUTPUT/field.tec',
                          '--tol-linf', '1e-9', '--tol-l2', '1e-12'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    print(res.stdout)
    if res.returncode != 0:
        print(res.stderr)
        print('(2) FAIL: MOSE bit-exact mismatch'); sys.exit(1)
    print('(2) MOSE bit-exact reference: OK')

print('Correct!')
