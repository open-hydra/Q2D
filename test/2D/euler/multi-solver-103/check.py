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

Acceptance (hard, NaN-safe via q2d_check_lib):
    - simulation completes iter=500 without NaN/Inf anywhere in the field
    - min(rho), min(p), min(T) > 0
    - velocity remains bounded (|u|, |v| < 1 m/s)
    - pressure stays near p_0 (within 1%)
"""
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_positive, assert_finite,
                           read_q2d_field)

require_file('OUTPUT/field.tec', "did ./Q2D.sh solve run?")
x_, y_, var_, vnames = read_q2d_field('OUTPUT/field.tec')

rho = var_[0][0][:, :, 0]
u = var_[0][1][:, :, 0]
v = var_[0][2][:, :, 0]
p = var_[0][3][:, :, 0]
T = var_[0][4][:, :, 0]

# Sanity (hard, NaN-safe): a NaN anywhere is a dispatcher failure.
assert_positive("rho", rho)
assert_positive("p", p)
assert_positive("T", T)
assert_finite("u", u)
assert_finite("v", v)

u_max = float(np.abs(u).max())
v_max = float(np.abs(v).max())
p_min, p_max = float(p.min()), float(p.max())
p_0 = 1.0e5
p_rel_dev = max(abs(p_max - p_0), abs(p_min - p_0)) / p_0
print(f'max |u| = {u_max:.4e} m/s')
print(f'max |v| = {v_max:.4e} m/s')
print(f'p range: [{p_min:.4f}, {p_max:.4f}] Pa (target {p_0})')
print(f'pressure rel deviation from p_0: {p_rel_dev:.3e}')

thr_uv = 1.0
thr_p = 0.01
if u_max > thr_uv:
    fail(f'Velocity u drift exceeds threshold: {u_max:.3e} > {thr_uv}')
if v_max > thr_uv:
    fail(f'Velocity v drift exceeds threshold: {v_max:.3e} > {thr_uv}')
if p_rel_dev > thr_p:
    fail(f'Pressure drift exceeds threshold: {p_rel_dev:.3e} > {thr_p}')

print('Correct!')
