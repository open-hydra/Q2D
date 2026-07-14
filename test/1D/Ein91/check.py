"""Einfeldt 1991 (Ein91) strong-rarefaction 1D regression check.

Compares ``OUTPUT/field.tec`` (Q2D_new tecplot ASCII BLOCK packing, peculiarita'
#1: 2D mesh) against the analytical Riemann solution in ``reference/exact.txt``
(5 cols: x, rho, u, p, sie). Mirror of Sod79/check.py with tolerance widened (5%)
to account for the harder Ein91 problem (low-density rarefaction near contact).

ORION's read_TEC assumes 3D file (x,y,z nodal); for the Q2D 2D output we insert
a z=0 nodal block on the fly + patch VARIABLES/VARLOCATION header so ORION reads
cleanly.
"""
import sys
import os
import re
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
    """Insert z=0 nodal block + update VARIABLES/VARLOCATION header so ORION
    (which assumes 3 nodal vars x,y,z) can read a Q2D 2D file. Returns the
    /tmp patched-file path."""
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


ref = np.loadtxt('reference/exact.txt')
xe, rhoe, ue, pe, siee = ref[:, 0], ref[:, 1], ref[:, 2], ref[:, 3], ref[:, 4]

patched = patch_q2d_tec_to_3d('OUTPUT/field.tec')
try:
    [x_, y_, z_, var_, vnames] = read_TEC(patched)
finally:
    os.unlink(patched)

xc  = 0.5 * (x_[0][:-1] + x_[0][1:])
rho = var_[0][0][:, 0, 0]
u   = var_[0][1][:, 0, 0]
p   = var_[0][3][:, 0, 0]
T   = var_[0][4][:, 0, 0]
gam = var_[0][5][:, 0, 0]
R   = var_[0][6][:, 0, 0]
cv  = R / (gam - 1)
sie = cv * T

rhoe_i = np.interp(xc[:, 0, 0], xe, rhoe)
ue_i   = np.interp(xc[:, 0, 0], xe, ue)
pe_i   = np.interp(xc[:, 0, 0], xe, pe)

ncells = len(rho)
# Probe on the smooth rarefaction tails away from the central trough.
iprobe = [int(ncells * f) for f in [0.05, 0.15, 0.85, 0.95]]
threshold = 5.0  # percent — looser than Sod79 due to strong rarefaction

if rho.min() <= 0.0:
    print(f'Sanity: min(rho)={rho.min()} <= 0'); sys.exit(1)
if p.min() <= 0.0:
    print(f'Sanity: min(p)={p.min()} <= 0'); sys.exit(1)
if T.min() <= 0.0:
    print(f'Sanity: min(T)={T.min()} <= 0'); sys.exit(1)

for i in iprobe:
    perr = abs((rho[i] - rhoe_i[i]) / rhoe_i[i]) * 100
    if perr > threshold:
        print(f'Density at cell {i}:  CFD={rho[i]:.6f}  Exact={rhoe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')
        sys.exit(1)

for i in iprobe:
    # Ein91 has |u| ≈ 2 across the field — pure relative check is safe.
    perr = abs((u[i] - ue_i[i]) / (abs(ue_i[i]) + 1e-10)) * 100
    if perr > threshold:
        print(f'Velocity at cell {i}:  CFD={u[i]:.6f}  Exact={ue_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')
        sys.exit(1)

for i in iprobe:
    perr = abs((p[i] - pe_i[i]) / pe_i[i]) * 100
    if perr > threshold:
        print(f'Pressure at cell {i}:  CFD={p[i]:.6f}  Exact={pe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')
        sys.exit(1)

print('Correct!')
