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

Gates (all hard):
  - sanity: min(rho), min(p), min(T) > 0 (NaN-safe via q2d_check_lib);
  - dp/dy vs Coriolis analytic: eps_rel < 5e-2 (piano section-5 tolerance;
    Coriolis-dominant regime, Cen/Cor < 5% at corner);
  - drift: the IC (u=U_inf, v=0) is an exact fixed point of the equilibrium;
    the solved field must not drift from it:
        max|u - U_inf| < 1.5e-2 m/s   (measured 6.1e-3 on the committed run)
        max|v|         < 8.0e-2 m/s   (measured 3.7e-2 on the committed run)
    both ~2x the measured drift of the committed 2000-iteration run.
"""
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_positive, assert_finite,
                           read_q2d_field)

U_INF = 100.0    # m/s, IC axial flow (reference/dpdy_analytical.dat header)

require_file('OUTPUT/field.tec', "did ./Q2D.sh solve run?")
x_, y_, var_, vnames = read_q2d_field('OUTPUT/field.tec')

xn = x_[0][:, :, 0]
yn = y_[0][:, :, 0]
rho = var_[0][0][:, :, 0]
u = var_[0][1][:, :, 0]
v = var_[0][2][:, :, 0]
p = var_[0][3][:, :, 0]
T = var_[0][4][:, :, 0]

# Sanity (hard, NaN-safe)
assert_positive("rho", rho)
assert_positive("p", p)
assert_positive("T", T)
assert_finite("u", u)
assert_finite("v", v)

yc = 0.5 * (yn[:-1, :-1] + yn[1:, 1:])
p_avg_y = p.mean(axis=0)
y_avg = yc[0, :]

slope, intercept = np.polyfit(y_avg, p_avg_y, 1)
print(f'CFD slope dp/dy (mean across x): {slope:.4e} Pa/m')

ref = np.loadtxt(require_file('reference/dpdy_analytical.dat'), comments='#')
dpdy_ref = float(ref if ref.ndim == 0 else ref[0])
print(f'Analytical dp/dy = -2*omega*rho_0*U_inf: {dpdy_ref:.4e} Pa/m')

eps_rel = abs(slope - dpdy_ref) / abs(dpdy_ref)
print(f'Coriolis gradient rel error: eps_rel = {eps_rel:.3e}')

u_drift = float(np.abs(u - U_INF).max())
v_max = float(np.abs(v).max())
print(f'u mean = {u.mean():.4f} m/s (target {U_INF}), max |u-U_inf| = {u_drift:.4e}')
print(f'v mean = {v.mean():.4e} m/s (target 0), max |v| = {v_max:.4e}')

thr = 5.0e-2  # piano section-5 tolerance: 5%
if eps_rel > thr:
    fail(f'dp/dy mismatch vs Coriolis analytical: {eps_rel:.3e} > {thr}')

# Drift gates (hard): IC is an exact fixed point of the Coriolis
# equilibrium; the committed 2000-iteration run measures
# max|u-U_inf| = 6.1e-3 m/s and max|v| = 3.7e-2 m/s (discretization-level
# settling, not growth). Tolerances ~2x measured: 1.5e-2 and 8.0e-2 m/s.
thr_u = 1.5e-2
thr_v = 8.0e-2
if u_drift > thr_u:
    fail(f'u drift from equilibrium: max|u-U_inf| = {u_drift:.3e} m/s > {thr_u}')
if v_max > thr_v:
    fail(f'v drift from equilibrium: max|v| = {v_max:.3e} m/s > {thr_v}')

print('Correct!')
