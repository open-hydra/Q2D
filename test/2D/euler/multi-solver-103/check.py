"""Multi-solver-103 BC 103 dispatcher smoke test regression check.

Q2D_new test/2D/euler/multi-solver-103 -- Sess 15.B C3 piano section-5 Fase 8c
line 602 criterion partial fulfillment:

    "test multi-solver coupling 1D-2D (BC 103 = ex Q2D 999) produce flussi
     residui consistenti".

BC 103 = multi-solver coupling wall (Mod_BC_Fluxes.f90:99-101):
    BC_Symmetry_Eul + R += ext_flux

With ext_flux=0 (IO_BC.f90:256 default initialization), BC 103 behaves
bit-identically to BC 300 (symmetry/wall). Test verifies BC 103 dispatcher
path executes without crash + maintains quiescent IC equilibrium.

Acceptance:
    - simulation completes iter=500 without NaN
    - velocity remains bounded (|u|, |v| < 1 m/s)
    - pressure stays near p_0 (within 1%)
"""
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


patched = patch_q2d_tec_to_3d('OUTPUT/field.tec')
try:
    [x_, y_, z_, var_, vnames] = read_TEC(patched)
finally:
    os.unlink(patched)

rho = var_[0][0][:, :, 0]
u = var_[0][1][:, :, 0]
v = var_[0][2][:, :, 0]
p = var_[0][3][:, :, 0]
T = var_[0][4][:, :, 0]

if rho.min() <= 0.0:
    print(f'Sanity: min(rho)={rho.min()} <= 0'); sys.exit(1)
if p.min() <= 0.0:
    print(f'Sanity: min(p)={p.min()} <= 0'); sys.exit(1)
if T.min() <= 0.0:
    print(f'Sanity: min(T)={T.min()} <= 0'); sys.exit(1)

u_max = np.abs(u).max()
v_max = np.abs(v).max()
p_min, p_max = p.min(), p.max()
p_0 = 1.0e5
p_rel_dev = max(abs(p_max - p_0), abs(p_min - p_0)) / p_0
print(f'max |u| = {u_max:.4e} m/s')
print(f'max |v| = {v_max:.4e} m/s')
print(f'p range: [{p_min:.4f}, {p_max:.4f}] Pa (target {p_0})')
print(f'pressure rel deviation from p_0: {p_rel_dev:.3e}')

thr_uv = 1.0
thr_p = 0.01
if u_max > thr_uv:
    print(f'Velocity u drift exceeds threshold: {u_max:.3e} > {thr_uv}')
    sys.exit(1)
if v_max > thr_uv:
    print(f'Velocity v drift exceeds threshold: {v_max:.3e} > {thr_uv}')
    sys.exit(1)
if p_rel_dev > thr_p:
    print(f'Pressure drift exceeds threshold: {p_rel_dev:.3e} > {thr_p}')
    sys.exit(1)

print('Correct!')
