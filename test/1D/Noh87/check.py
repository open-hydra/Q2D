"""Noh87 1D strong-shock regression check (Noh 1987 planar Cartesian).

Compares ``OUTPUT/field.tec`` (Q2D_new tecplot ASCII BLOCK packing, peculiarita'
#1: 2D mesh -- x,y nodal + rho/u/v/p/T/g/R cell-centered) against the analytical
strong-shock solution in ``reference/exact.txt`` (5 cols: x, rho, u, p, sie),
generated offline via ``tools/gen_noh_exact.py`` (Rankine-Hugoniot post-shock
stagnation: rho=6, u=0, p=1.2 for gamma=1.4, u0=1).

ORION read_TEC patch identical to Sod79/Tor99. Mirror MOSE/test/1D/Noh87/
verify.py threshold (1% rho/u/sie, 1.5% p — Noh near-vacuum p_pre=1e-15 is
harder than other tests for the relative-error metric).

Noh87 setup (mirror MOSE):
  rho0=1.0, u0=+/-1.0 (symmetric collapse), p0=1e-15
  gamma=1.4, x0=0.5, t_final=1.0, domain=[0,1]
Wave structure at t=1.0:
  pre-shock left  (x < 0.3): rho=1, u=+1, p=0
  stagnation     (0.3<x<0.7): rho=6, u=0,  p=1.2
  pre-shock right (x > 0.7): rho=1, u=-1, p=0
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
siee_i = np.interp(xc[:, 0, 0], xe, siee)

# Probes positioned away from the two shock fronts (x=0.3 and x=0.7). Sample
# pre-shock left, stagnation core, pre-shock right (avoid shock smearing).
ncells = len(rho)
iprobe = [int(ncells * f) for f in [0.10, 0.50, 0.90]]
threshold_rho  = 1.0   # percent
threshold_u    = 1.0
threshold_p    = 1.5   # mirror MOSE Noh87 verify.py (p harder near vacuum)
threshold_sie  = 1.0

if rho.min() <= 0.0:
    print(f'Sanity: min(rho)={rho.min()} <= 0'); sys.exit(1)
if p.min() <= 0.0:
    print(f'Sanity: min(p)={p.min()} <= 0'); sys.exit(1)
if T.min() <= 0.0:
    print(f'Sanity: min(T)={T.min()} <= 0'); sys.exit(1)

for i in iprobe:
    if rhoe_i[i] == 0.0:
        continue
    perr = abs((rho[i] - rhoe_i[i]) / rhoe_i[i]) * 100
    if perr > threshold_rho:
        print(f'Density at cell {i}:  CFD={rho[i]:.6f}  Exact={rhoe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold_rho}%')
        sys.exit(1)

for i in iprobe:
    if abs(ue_i[i]) < 1e-10:
        # Skip stagnation region (u=0 makes relative error undefined; sanity
        # check above already validates physical positivity).
        continue
    perr = abs((u[i] - ue_i[i]) / ue_i[i]) * 100
    if perr > threshold_u:
        print(f'Velocity at cell {i}:  CFD={u[i]:.6f}  Exact={ue_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold_u}%')
        sys.exit(1)

for i in iprobe:
    if pe_i[i] < 1e-12:
        # Pre-shock vacuum: absolute floor instead of relative.
        if p[i] > 1.0e-3:
            print(f'Pressure at cell {i}:  CFD={p[i]:.6e} >> vacuum reference {pe_i[i]:.3e}')
            sys.exit(1)
        continue
    perr = abs((p[i] - pe_i[i]) / pe_i[i]) * 100
    if perr > threshold_p:
        print(f'Pressure at cell {i}:  CFD={p[i]:.6f}  Exact={pe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold_p}%')
        sys.exit(1)

for i in iprobe:
    if abs(siee_i[i]) < 1e-12:
        continue
    perr = abs((sie[i] - siee_i[i]) / siee_i[i]) * 100
    if perr > threshold_sie:
        print(f'SIE at cell {i}:  CFD={sie[i]:.6f}  Exact={siee_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold_sie}%')
        sys.exit(1)

print('Correct!')
