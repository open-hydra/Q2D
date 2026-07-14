"""prandtl-meyer 2D steady-state regression check.

Validates the converged solution of the Prandtl-Meyer expansion fan against
the analytical isentropic expansion result (cf. Prandtl-Meyer 1908). Mirror
MOSE/test/2D/euler/prandtl-meyer/check.py.

Reads the converged Mach at the exit-probe location from ``OUTPUT/exit.txt``
(Q2D-Probes section in input.ini). Tolerance: |M_CFD - M_analytical| / M_analytical
< 0.5%.

Also enforces sanity check on the field (regression_tolerance.md Section 3.2):
``min(rho), min(p), min(T) > 0`` everywhere.
"""
import math
import os
import re
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


gamma = 1.4


def nufunc(M, gamma):
    return math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(
        math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) - math.atan(math.sqrt(M ** 2 - 1))


def mach_from_nu(nu, gamma, a=1.01, b=30.0, tol=1e-6, max_iter=1000000):
    for _ in range(max_iter):
        M_mid = (a + b) / 2
        nu_err = nufunc(M_mid, gamma) - nu
        if abs(nu_err) < tol:
            return M_mid
        if np.sign(nu_err) < 0:
            a = M_mid
        else:
            b = M_mid
    raise RuntimeError("Bisection did not converge.")


def prandtl_meyer_expansion(theta, M):
    nu1 = nufunc(M, gamma)
    nu2 = nu1 - theta * np.pi / 180 if theta < 0 else nu1 + theta * np.pi / 180
    return mach_from_nu(nu2, gamma)


# Field-level sanity check
patched = patch_q2d_tec_to_3d('OUTPUT/field.tec')
try:
    [x_, y_, z_, var_, vnames] = read_TEC(patched)
finally:
    os.unlink(patched)

rho = var_[0][0]
p   = var_[0][3]
T   = var_[0][4]
if rho.min() <= 0.0:
    print(f'Sanity: min(rho)={rho.min()} <= 0'); sys.exit(1)
if p.min() <= 0.0:
    print(f'Sanity: min(p)={p.min()} <= 0'); sys.exit(1)
if T.min() <= 0.0:
    print(f'Sanity: min(T)={T.min()} <= 0'); sys.exit(1)

# Analytical Prandtl-Meyer Mach at the exit
M1 = 2.0
theta = 34.971
M2_analytical = prandtl_meyer_expansion(theta, M1)

# CFD Mach at exit probe (last line of probe output)
with open('OUTPUT/exit.txt', 'r') as file:
    lines = file.readlines()
    last_line = lines[-1]
    columns = last_line.strip().split()
    CFDval = columns[1]

perr = abs(float(CFDval) - M2_analytical) / M2_analytical * 100
if perr > 0.5:
    print('Mach number mismatch (CFD vs analytical)')
    print(f'CFD={CFDval}  Analytical={M2_analytical}')
    print(f'Error = {perr}%, threshold = 0.5%')
    sys.exit(1)
else:
    print(f'(1) Natural reference: M_CFD={CFDval}, M_analytical={M2_analytical:.6f}, eps_rel={perr:.4f}%')

# (2) MOSE bit-exact reference: Q2D field K=0 vs cached MOSE field (2D plane case;
# bit-exact equivalence verified iter-by-iter from Sess 5b-bis 2026-05-15).
import subprocess
mose_ref = Path(__file__).resolve().parents[4] / 'test' / 'regression_baseline' / 'MOSE_K0' / 'prandtl-meyer' / 'field.tec'
diff_tool = Path(__file__).resolve().parents[4] / 'tools' / 'diff_mose_q2d_field.py'
if not mose_ref.exists():
    print(f'(2) WARN: MOSE bit-exact reference not found at {mose_ref}; skipping')
elif not diff_tool.exists():
    print(f'(2) WARN: diff utility not found at {diff_tool}; skipping')
else:
    print(f"(2) MOSE bit-exact reference vs Q2D field (Linf_rel < 1e-9, L2_rel < 1e-12, atol 1e-10):")
    res = subprocess.run([sys.executable, str(diff_tool), str(mose_ref), 'OUTPUT/field.tec',
                          '--tol-linf', '1e-9', '--tol-l2', '1e-12'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    print(res.stdout)
    if res.returncode != 0:
        print(res.stderr)
        print('(2) FAIL: MOSE bit-exact mismatch'); sys.exit(1)
    print('(2) MOSE bit-exact reference: OK')

print('Correct!')
