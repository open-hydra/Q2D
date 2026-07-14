"""SRM-rocket-motor 2D-axi steady-state regression check.

Validates the converged chamber pressure of a solid rocket motor with Vieille
burn-rate law at the closed-end grain face (BC 502 SRM). Q2D-extra test
(MOSE has Lib_BC_Fluxes_SRM but no SRM test case mirror); piano §5 Fase 8c
line 602 acceptance criterion C2:

    "test SRM 2D-axi con BC 502 produce evolution di pressione coerente con
     regressione Vieille".

Geometry: chang-nozzle 2D-axi mesh (190x40, R_chamber=0.0635 m, R_throat=0.0203 m,
Kn=9.75). The contraction ratio is small for a real SRM, so the burn-rate
coefficient aCoeff is boosted (50 kg/m^2/s vs ~9 for real AP/HTPB at pRef=1 MPa)
to guarantee a choked-throat steady-state Vieille equilibrium within reach of a
ctest budget.

Vieille steady-state balance with choked throat:

    aCoeff * (P/pRef)^n * A_chamber = P * A_throat / c_star

    => P_eq = [aCoeff * c_star * (A_chamber/A_throat) / pRef^n]^(1/(1-n))

Acceptance tolerance per piano §5 Fase 8c (analogous chang-nozzle pattern ±5%
avg, ±15% max relaxed for spatial-vs-lumped discretization residual):

    |P_chamber_avg - P_eq_analytical| / P_eq_analytical < 0.15

Reads reference/P_eq_vieille.dat for the analytical reference value.
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
u   = var_[0][1][:, :, 0]
v   = var_[0][2][:, :, 0]
p   = var_[0][3][:, :, 0]
T   = var_[0][4][:, :, 0]

if rho.min() <= 0.0:
    print(f'Sanity: min(rho)={rho.min()} <= 0'); sys.exit(1)
if p.min() <= 0.0:
    print(f'Sanity: min(p)={p.min()} <= 0'); sys.exit(1)
if T.min() <= 0.0:
    print(f'Sanity: min(T)={T.min()} <= 0'); sys.exit(1)

res_path = Path('OUTPUT/residual-history.dat')
if res_path.exists():
    res = np.loadtxt(res_path, comments=('#',))
    if res.ndim == 2 and res.shape[0] > 10:
        res_norm = np.linalg.norm(res[:, 1:], axis=1)
        first_res = res_norm[0] if res_norm[0] > 0 else 1.0
        last_res = res_norm[-1]
        rel_drop = last_res / first_res
        print(f'Residual drop (L2 norm last/first): {rel_drop:.3e}')

profile = np.loadtxt('MESH/profile.dat')
x_throat = profile[np.argmin(profile[:, 1]), 0]
print(f'Throat location: x={x_throat:.4f} m')

xc = 0.5 * (xn[:-1, :-1] + xn[1:, 1:])
mask_chamber = xc[:, 0] < (x_throat - 0.005)
n_cells_chamber = int(mask_chamber.sum())
print(f'Chamber cells (x < x_throat - 0.005): {n_cells_chamber} / {xc.shape[0]}')
if n_cells_chamber < 5:
    print('ERROR: fewer than 5 chamber cells; check geometry'); sys.exit(1)

p_chamber = p[mask_chamber, :].mean()
p_chamber_std = p[mask_chamber, :].std()
print(f'P_chamber (mean over chamber cells): {p_chamber:.4e} Pa = {p_chamber/1e6:.4f} MPa')
print(f'P_chamber std-dev: {p_chamber_std:.4e} Pa ({100*p_chamber_std/p_chamber:.2f} %)')

ref_path = Path('reference/P_eq_vieille.dat')
ref_data = np.loadtxt(ref_path, comments='#')
P_eq_ref = float(ref_data if ref_data.ndim == 0 else ref_data[0])
print(f'P_eq Vieille analytical: {P_eq_ref:.4e} Pa = {P_eq_ref/1e6:.4f} MPa')

eps_rel = abs(p_chamber - P_eq_ref) / P_eq_ref
print(f'Vieille rel error: eps_rel = {eps_rel:.3e}')

thr = 0.15
if eps_rel > thr:
    print(f'P_chamber mismatch vs Vieille analytical: {eps_rel:.3e} > {thr}')
    sys.exit(1)

print('Correct!')
