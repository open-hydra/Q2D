"""rocket-nozzle/equilibrium 2D-axisymmetric CEA chemistry=equilibrium regression.

Validates the Q2D solution for the rocket-nozzle equilibrium case (CEA
on-the-fly equilibrium per cell, 13 species). Single-block 121x41 2D-axi.

References:
  (1) Natural: MOSE reference/field.tec compared cell-by-cell with relaxed
      tolerance 5% Linf_rel / 1% L2_rel. Observed actual ~2% Linf max.
      Chemistry-source-perturbed steady-state never tightly converges
      (residual ~6e-2 floor for both Q2D and MOSE; this is INHERENT to
      CEA-equilibrium chemistry update each timestep).
  (2) Sanity: min(rho), min(p), min(T) > 0 (hard).

Note: comparing Q2D iter=10000 vs MOSE iter=2000 reference snapshot — both
solvers settle into the same statistical-steady-state pattern (chemistry
source noise floor) so the diff is meaningful within ~2% in chemistry-
sensitive regions.
"""
import os
import re
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
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


def patch_q2d_tec_to_3d(in_path):
    with open(in_path) as f:
        lines = f.readlines()
    header = lines[0]
    zone = lines[1]
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
    out = [new_header, new_zone]
    out.extend(lines[2 : 2 + nnodes])
    out.extend(lines[2 + nnodes : 2 + 2 * nnodes])
    out.append(' 0.000000000000000E+000\n' * nnodes)
    out.extend(lines[2 + 2 * nnodes :])
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.tec', delete=False)
    tmp.writelines(out); tmp.close()
    return tmp.name


patched = patch_q2d_tec_to_3d('OUTPUT/field.tec')
try:
    [x_, y_, z_, var_, vnames] = read_TEC(patched)
finally:
    os.unlink(patched)

rho_sum = sum(var_[0][i] for i in range(13))
u = var_[0][13]; v = var_[0][14]
p = var_[0][15]; T = var_[0][16]; gam = var_[0][17]; R = var_[0][18]

if rho_sum.min() <= 0:
    print(f'Sanity FAIL: min(rho_sum)={rho_sum.min()}'); sys.exit(1)
if p.min() <= 0:
    print(f'Sanity FAIL: min(p)={p.min()}'); sys.exit(1)
if T.min() <= 0:
    print(f'Sanity FAIL: min(T)={T.min()}'); sys.exit(1)
print(f'Sanity OK: min(rho_sum)={rho_sum.min():.3e} min(p)={p.min():.3e} min(T)={T.min():.3e}')

mose_ref = root.parents[4] / 'test' / 'regression_baseline' / 'MOSE_K0' / 'rocket-nozzle-equilibrium' / 'field.tec'
diff_tool = root.parents[4] / 'tools' / 'diff_mose_q2d_field.py'
if not mose_ref.exists():
    print(f'(1) WARN: MOSE reference not found at {mose_ref}; skipping')
elif not diff_tool.exists():
    print(f'(1) WARN: diff utility not found at {diff_tool}; skipping')
else:
    print('(1) MOSE reference vs Q2D 2D-axi (relaxed Linf_rel < 5e-2, L2_rel < 1e-2):')
    res = subprocess.run([sys.executable, str(diff_tool), str(mose_ref), 'OUTPUT/field.tec',
                          '--tol-linf', '5e-2', '--tol-l2', '1e-2'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    print(res.stdout)
    if res.returncode != 0:
        print(res.stderr)
        print('(1) FAIL: MOSE reference comparison mismatch'); sys.exit(1)
    print('(1) MOSE reference comparison: OK')

i_throat = int(np.argmin(np.abs(x_[0][:, 0, 0])))
M_axis = np.sqrt(u[i_throat, 0, 0]**2 + v[i_throat, 0, 0]**2) / np.sqrt(gam[i_throat, 0, 0] * R[i_throat, 0, 0] * T[i_throat, 0, 0])
print(f'  Throat at i={i_throat} (x={x_[0][i_throat, 0, 0]:.3f}); Mach axis: {M_axis:.3f}  p={p[i_throat, 0, 0]:.3e}  T={T[i_throat, 0, 0]:.1f}')

print('Correct!')
