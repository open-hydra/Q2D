"""Rotating-channel 2D plane Euler Coriolis equilibrium regression check.

Mirror MOSE topology pattern: periodic in x (BC 201) + walls in y (BC 300).
Sess 15.B 2026-05-25 enabled via BC 201 port to Q2D Lib_Ghost.f90:63 +
Mod_BC_Fluxes.f90:103 (mirror MOSE Lib_Ghost.f90:51 + Mod_BC_Fluxes.f90:72).

Piano section-5 Fase 8e line 619 criterion C5 fulfillment:

    "test 2D piano canale rotante (Coriolis significativo) confronta con
     soluzione analitica".

Geometry: 2D channel x in [-L_x/2, L_x/2], y in [-L_y/2, L_y/2]; periodic-x via
BC 201 (Q2D-Lib_BC_Fluxes_Connection mirror MOSE) + walls-y via BC 300.
Rotation z-axis omega=10 rad/s; uniform axial flow u=U_inf=100 m/s.

Lib_RotatingFrame.f90:175-242 Coriolis balance in inviscid Euler (peculiarita 1
No-w, axis=[0,0,1] strict):

    Cor_x = +2*omega*v       (=0 for v=0 IC)
    Cor_y = -2*omega*u       (=-2*omega*U_inf, force in -y direction)

Steady-state isothermal equilibrium:
    dp/dy = -2*omega*rho*u
    p(y) = p_0 * exp(-2*omega*U_inf*y/(R*T_0))     [isothermal: T_0 const]
    rho(y) = p(y)/(R*T_0)

Tolerance: gradient relative error < 5% (Coriolis-dominant regime, Cen/Cor < 5%
at corner). Drift test: IC is exact fixed-point; check CFD does NOT drift.
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

xn = x_[0][:, :, 0]
yn = y_[0][:, :, 0]
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

yc = 0.5 * (yn[:-1, :-1] + yn[1:, 1:])
p_avg_y = p.mean(axis=0)
y_avg = yc[0, :]

slope, intercept = np.polyfit(y_avg, p_avg_y, 1)
print(f'CFD slope dp/dy (mean across x): {slope:.4e} Pa/m')

ref = np.loadtxt('reference/dpdy_analytical.dat', comments='#')
dpdy_ref = float(ref if ref.ndim == 0 else ref[0])
print(f'Analytical dp/dy = -2*omega*rho_0*U_inf: {dpdy_ref:.4e} Pa/m')

eps_rel = abs(slope - dpdy_ref) / abs(dpdy_ref)
print(f'Coriolis gradient rel error: eps_rel = {eps_rel:.3e}')

u_mean = u.mean()
u_std = u.std()
v_mean = v.mean()
v_max = abs(v).max()
print(f'u mean = {u_mean:.4f} m/s (target 100.0), std = {u_std:.4e}')
print(f'v mean = {v_mean:.4e} m/s (target 0), max |v| = {v_max:.4e}')

thr = 5.0e-2  # piano section-5 tolerance: 5%
if eps_rel > thr:
    print(f'dp/dy mismatch vs Coriolis analytical: {eps_rel:.3e} > {thr}')
    sys.exit(1)

print('Correct!')
