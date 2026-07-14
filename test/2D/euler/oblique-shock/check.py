"""oblique-shock 2D steady-state regression check (DOUBLE reference).

Validates the converged 2D plane Euler solution for M=4 flow over a ramp
with attached oblique shock (2-block grid: Block1 pre-ramp flat 21x219,
Block2 post-ramp 39x219). Mirror MOSE verify.py probe layout (i_pre0/j_pre0
in Block1, i_post1/j_post1 in Block2).

References:
  (1) Natural: oblique-shock theta-beta-Mach theory with delta measured from
      Block2 wall nodes (i=0..end, j=0). Tolerance 1.5% per MOSE verify.py.
  (2) MOSE bit-exact: Q2D 2-zone field K=0 vs MOSE 2-zone field K=0.
      Reference cache test/regression_baseline/MOSE_K0/oblique-shock/field.tec.

Geometric classification (z-coord analysis): 2D PLANE — z=+-5e-4 const per
K-layer per zone, X/Y identical between K-layers => no axisymmetric source.

Tolerances:
  Sanity: min(rho), min(p), min(T) > 0 (hard).
  Natural: 1.5% per MOSE verify.py (HLLE-only single-Riemann run).
  Bit-exact MOSE: Linf_rel < 1e-9, L2_rel < 1e-12.
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
    raise ModuleNotFoundError("ORION package not found. Run 'conda activate base'.")


def patch_q2d_multizone_tec_to_3d(in_path):
    """Patch Q2D 2D-strip multi-zone field.tec (K=1) to ORION-readable 3D-wedge
    (K=2 with fake z=0). Multi-zone variant (mirror chang-nozzle helper)."""
    with open(in_path) as f:
        lines = f.readlines()
    # Header
    out = []
    header = lines[0]
    new_header = header.replace('"y"', '"y" "z"', 1)
    out.append(new_header)
    i = 1
    # walk zones
    while i < len(lines):
        l = lines[i]
        if 'ZONE' in l.upper() and 'I=' in l.upper():
            zone = l
            m = re.search(r'I=\s*(\d+).*J=\s*(\d+).*K=\s*(\d+)', zone)
            if m is None:
                out.append(l); i += 1; continue
            ni, nj, nk = int(m.group(1)), int(m.group(2)), int(m.group(3))
            nnodes = ni * nj * nk
            m2 = re.search(r'VARLOCATION=\(\[(\d+)-(\d+)\]=NODAL,\[(\d+)-(\d+)\]=CELLCENTERED', zone)
            if m2:
                na, nb_, ca, cb = (int(m2.group(j)) for j in (1, 2, 3, 4))
                new_zone = zone.replace(
                    f'[{na}-{nb_}]=NODAL,[{ca}-{cb}]=CELLCENTERED',
                    f'[{na}-{nb_+1}]=NODAL,[{ca+1}-{cb+1}]=CELLCENTERED',
                )
            else:
                new_zone = zone
            out.append(new_zone)
            i += 1
            # Read nodal vars: na .. nb_ — original n_nodal = nb_ - na + 1
            n_nodal_orig = nb_ - na + 1 if m2 else 2
            n_cc_orig    = cb - ca + 1 if m2 else 0
            # Copy nodal arrays (each nnodes lines)
            for vi in range(n_nodal_orig):
                # one var = nnodes lines
                out.extend(lines[i : i + nnodes])
                i += nnodes
            # Inject fake z column: nnodes lines of '0.0'
            out.append(' 0.000000000000000E+000\n' * nnodes)
            # Copy cellcentered (each ncells lines)
            ncells = (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1)
            for vi in range(n_cc_orig):
                out.extend(lines[i : i + ncells])
                i += ncells
        else:
            out.append(l)
            i += 1
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.tec', delete=False)
    tmp.writelines(out); tmp.close()
    return tmp.name


# (1) Natural ref: oblique-shock theta-beta-Mach (mirror MOSE verify.py sckobl)
def tbm(M, beta, gamma):
    c1 = M*M * math.sin(beta)**2 - 1.0
    c2 = M*M * (gamma + math.cos(2.0*beta)) + 2.0
    return 2.0 * (1.0 / math.tan(beta)) * (c1 / c2)


def sckobl(mach1, delta, gamma, tol):
    theta1 = 1.2 * delta
    delt1 = tbm(mach1, theta1, gamma)
    theta2 = 1.02 * theta1
    while True:
        delt2 = tbm(mach1, theta2, gamma)
        if abs((delt2 - math.tan(delta)) / math.tan(delta)) <= tol:
            break
        theta = theta2 - (math.tan(delta) - delt2) * (theta2 - theta1) / (delt1 - delt2)
        theta1 = theta2; delt1 = delt2; theta2 = theta
    theta = theta2
    rmn1 = mach1 * math.sin(theta)
    gp1 = gamma + 1.0; gm1 = gamma - 1.0
    r2r1 = gp1 * rmn1**2 / (gm1 * rmn1**2 + 2.0)
    p2p1 = 1.0 + 2.0 * gamma * (rmn1**2 - 1.0) / gp1
    t1 = (rmn1**2 + 2.0 / gm1)
    t2 = (2.0 * gamma * rmn1**2 / gm1 - 1.0)
    rmn2 = math.sqrt(t1 / t2)
    t2t1 = p2p1 / r2r1
    mach2 = rmn2 / math.sin(theta - delta)
    return mach2, theta, p2p1, t2t1, r2r1


M1 = 4.0; gamma = 1.4; tol = 1e-6

patched = patch_q2d_multizone_tec_to_3d('OUTPUT/field.tec')
try:
    [x_, y_, z_, var_, vnames] = read_TEC(patched)
finally:
    os.unlink(patched)

# Sanity: all blocks
for b in range(len(var_)):
    if var_[b][0].min() <= 0: print(f'Sanity FAIL b{b} min(rho)={var_[b][0].min()}'); sys.exit(1)
    if var_[b][3].min() <= 0: print(f'Sanity FAIL b{b} min(p)={var_[b][3].min()}'); sys.exit(1)
    if var_[b][4].min() <= 0: print(f'Sanity FAIL b{b} min(T)={var_[b][4].min()}'); sys.exit(1)

# Measure ramp angle from Block2 wall nodes (j=0): MOSE verify.py L120
xw = x_[1][:, 0, 0]
yw = y_[1][:, 0, 0]
delta_deg = float(np.degrees(np.arctan2(yw[-1] - yw[0], xw[-1] - xw[0])))
delta = math.radians(delta_deg)
mach2_ref, beta, p2p1_ref, t2t1_ref, r2r1_ref = sckobl(M1, delta, gamma, tol)
print(f"(1) Natural reference (oblique shock M1={M1} gamma={gamma}):")
print(f"    delta={delta_deg:.4f} deg beta={math.degrees(beta):.4f} deg")
print(f"    p2/p1={p2p1_ref:.6f} T2/T1={t2t1_ref:.6f} r2/r1={r2r1_ref:.6f} M2={mach2_ref:.6f}")

# Pre-shock probe Block1 (mirror MOSE verify.py: i_pre0=5:15, j_pre0=5:15)
f0 = var_[0]
i_pre = slice(5, 15); j_pre = slice(5, 15)
p1_  = float(np.mean(f0[3][i_pre, j_pre, 0]))
rho1 = float(np.mean(f0[0][i_pre, j_pre, 0]))
T1_  = float(np.mean(f0[4][i_pre, j_pre, 0]))

# Post-shock probe Block2 (mirror MOSE verify.py: i_post1=14:35, j_post1=5:50)
f1 = var_[1]
i_post = slice(14, 35); j_post = slice(5, 50)
p2  = float(np.mean(f1[3][i_post, j_post, 0]))
rho2 = float(np.mean(f1[0][i_post, j_post, 0]))
T2  = float(np.mean(f1[4][i_post, j_post, 0]))
u2  = float(np.mean(f1[1][i_post, j_post, 0]))
v2  = float(np.mean(f1[2][i_post, j_post, 0]))
g2  = float(np.mean(f1[5][i_post, j_post, 0]))
R2  = float(np.mean(f1[6][i_post, j_post, 0]))

p2p1 = p2 / p1_
r2r1 = rho2 / rho1
t2t1 = T2 / T1_
a2 = math.sqrt(g2 * R2 * T2)
m2 = math.sqrt(u2**2 + v2**2) / a2

e_p = (p2p1 - p2p1_ref) / p2p1_ref * 100.0
e_T = (t2t1 - t2t1_ref) / t2t1_ref * 100.0
e_r = (r2r1 - r2r1_ref) / r2r1_ref * 100.0
e_m = (m2   - mach2_ref) / mach2_ref * 100.0

print(f"    Probe pre Block1:  p1={p1_:.6f} rho1={rho1:.6f} T1={T1_:.6f}")
print(f"    Probe post Block2: p2={p2:.6f} rho2={rho2:.6f} T2={T2:.6f} M2={m2:.4f}")
print(f"    Errors: p2/p1={e_p:+.2f}% T2/T1={e_T:+.2f}% r2/r1={e_r:+.2f}% M2={e_m:+.2f}%")

tol_pct = 1.5
fail = False
for label, err in [('p2/p1', e_p), ('T2/T1', e_T), ('rho2/rho1', e_r), ('M2', e_m)]:
    if abs(err) > tol_pct:
        print(f'(1) FAIL: {label} err {err:+.4f}% > {tol_pct}%'); fail = True
if fail: sys.exit(1)
print('(1) Natural reference: OK')

# (2) MOSE bit-exact reference
mose_ref = root.parents[3] / 'test' / 'regression_baseline' / 'MOSE_K0' / 'oblique-shock' / 'field.tec'
diff_tool = root.parents[3] / 'tools' / 'diff_mose_q2d_field.py'
if not mose_ref.exists():
    print(f'(2) WARN: MOSE bit-exact reference not found at {mose_ref}; skipping')
elif not diff_tool.exists():
    print(f'(2) WARN: diff utility not found at {diff_tool}; skipping')
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
