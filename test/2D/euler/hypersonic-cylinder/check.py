"""hypersonic-cylinder 2D steady-state regression check.

Validates the converged 2D plane Euler solution for hypersonic flow M=8.1
over a cylinder cross-section (bow shock at stagnation line).

  (1) Natural physical reference: normal-shock + isentropic-stagnation theory
      for inflow M1=8.1, gamma=1.4 (mirror MOSE verify.py probe at I=80 J=0):
      total p0 within 5% of analytical; total T0 within 1%.

  (2) Symmetry: the case is mirror-symmetric about the stagnation line
      (i <-> ni-1-i on the O-grid). Asymmetry of rho between mirrored cells
      (scaled by the local max) is gated: max < 4e-1, mean < 1e-2 (measured
      1.9e-1 max -- a single shock-band cell -- and 1.8e-3 mean on the
      committed run). A carbuncle-type breakdown produces O(1) asymmetry
      over wide regions and blows the mean gate immediately.

  (3) Bow-shock standoff sanity: the standoff distance of the density jump
      (rho crossing 2*rho_inf) along the stagnation line must fall within
      +-40% of the Billig (1967) cylinder correlation
      Delta/R = 0.386*exp(4.67/M^2) = 0.414 at M=8.1
      (measured 0.423 on the committed run).

Geometric classification (rigorous, z-coord analysis of MOSE wedge mesh):
  2D PLANE — z=+-5e-4 const per K-layer, X/Y identical between K-layers,
             half-angle wedge tilt = 0 deg => no axisymmetric source.
  Q2D config: [Q2D-Sources] active=none (omitted from input.ini).
  O-grid: wall at j=0 (cylinder R=0.5), farfield at j=jmax, stagnation
  point at the i-midline (flow from -x).

Sanity per ``docs/vv/regression_tolerance.md`` Section 3.2:
  min(rho), min(p), min(T) > 0 (hard, NaN-safe via q2d_check_lib).
"""
import math
import sys
import warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (fail, require_file, assert_positive, assert_finite,
                           read_q2d_field)


# ---------------------------------------------------------------------------
# (1) NATURAL REFERENCE: normal-shock + isentropic stagnation
# ---------------------------------------------------------------------------
def normal_shock(M1, gamma):
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0
    M1sq = M1 * M1
    M2sq = (M1sq + 2.0/gm1) / (2.0*gamma/gm1 * M1sq - 1.0)
    M2   = math.sqrt(M2sq)
    p2p1 = 1.0 + 2.0*gamma/gp1 * (M1sq - 1.0)
    r2r1 = gp1*M1sq / (gm1*M1sq + 2.0)
    t2t1 = p2p1 / r2r1
    return M2, p2p1, t2t1, r2r1


def total_pressure(p, M, gamma):
    return p * (1.0 + 0.5*(gamma - 1.0)*M*M) ** (gamma / (gamma - 1.0))


def total_temperature(T, M, gamma):
    return T * (1.0 + 0.5*(gamma - 1.0)*M*M)


M1    = 8.1
gamma = 1.4
p1    = 101325.0
T1    = 288.14

M2_ref, p2p1_ref, t2t1_ref, r2r1_ref = normal_shock(M1, gamma)
p2_ref = p2p1_ref * p1
T2_ref = t2t1_ref * T1
p0_ref = total_pressure(p2_ref, M2_ref, gamma)
T0_ref = total_temperature(T2_ref, M2_ref, gamma)

print(f"(1) Natural reference (normal shock M1={M1} gamma={gamma}):")
print(f"    M2={M2_ref:.6f}  p0_post={p0_ref:.4f} Pa  T0_post={T0_ref:.4f} K")

require_file('OUTPUT/field.tec', "did ./Q2D.sh solve run?")
x_, y_, var_, vnames = read_q2d_field('OUTPUT/field.tec')

xn = x_[0][:, :, 0]
yn = y_[0][:, :, 0]
rho = var_[0][0]; u = var_[0][1]; v = var_[0][2]
p   = var_[0][3]; T = var_[0][4]; g = var_[0][5]; R = var_[0][6]

# Sanity (hard, NaN-safe)
assert_positive("rho", rho)
assert_positive("p", p)
assert_positive("T", T)
assert_finite("u", u)
assert_finite("v", v)

ni, nj = rho.shape[0], rho.shape[1]

I_PROBE = 80
J_PROBE = 0
if I_PROBE >= ni or J_PROBE >= nj:
    fail(f"probe (I={I_PROBE},J={J_PROBE}) outside cell dims ({ni},{nj}) -- "
         "mesh changed? re-calibrate the probe")
p_pr = float(p[I_PROBE, J_PROBE, 0]); T_pr = float(T[I_PROBE, J_PROBE, 0])
u_pr = float(u[I_PROBE, J_PROBE, 0]); v_pr = float(v[I_PROBE, J_PROBE, 0])
g_pr = float(g[I_PROBE, J_PROBE, 0]); R_pr = float(R[I_PROBE, J_PROBE, 0])
speed = math.sqrt(u_pr**2 + v_pr**2)
a_pr  = math.sqrt(g_pr * R_pr * T_pr)
M_pr  = speed / a_pr
p0_pr = total_pressure(p_pr, M_pr, g_pr)
T0_pr = total_temperature(T_pr, M_pr, g_pr)
err_p0 = 100.0 * (p0_pr - p0_ref) / p0_ref
err_T0 = 100.0 * (T0_pr - T0_ref) / T0_ref
print(f"    Probe (I={I_PROBE} J={J_PROBE}): p={p_pr:.4e}  T={T_pr:.4f}  M={M_pr:.6f}")
print(f"    Q2D: p0={p0_pr:.4f} Pa err={err_p0:+.4f}%  T0={T0_pr:.4f} K err={err_T0:+.4f}%")

thr_p0_natural = 5.0
thr_T0_natural = 1.0
if not abs(err_p0) <= thr_p0_natural:
    fail(f'(1) p0 natural mismatch {err_p0:+.4f}% > {thr_p0_natural}%')
if not abs(err_T0) <= thr_T0_natural:
    fail(f'(1) T0 natural mismatch {err_T0:+.4f}% > {thr_T0_natural}%')
print('(1) Natural reference: OK')

# ---------------------------------------------------------------------------
# (2) SYMMETRY about the stagnation line (carbuncle sentinel)
# ---------------------------------------------------------------------------
# The O-grid runs i=0..ni-1 around the half-cylinder with the stagnation
# line at the i-midline; the exact solution satisfies rho(i,j)=rho(ni-1-i,j).
rho2d = rho[:, :, 0]
mir = rho2d[::-1, :]
asym = np.abs(rho2d - mir) / np.maximum(np.abs(rho2d), np.abs(mir))
asym_max = float(asym.max())
asym_mean = float(asym.mean())
print(f"(2) Symmetry: rho mirror asymmetry max={asym_max:.3e} mean={asym_mean:.3e}")
# Measured on the committed run: max 1.9e-1 (one shock-band cell: a half-cell
# shock offset between mirrored halves), mean 1.8e-3. Tolerances: max 4.0e-1
# (~2x measured; still < the O(1) of a broken shock), mean 1.0e-2 (~5x
# measured). A carbuncle breaks the mean gate by orders of magnitude.
thr_asym_max, thr_asym_mean = 4.0e-1, 1.0e-2
if not asym_max <= thr_asym_max:
    fail(f'(2) rho mirror asymmetry max {asym_max:.3e} > {thr_asym_max} '
         '(carbuncle-type shock asymmetry?)')
if not asym_mean <= thr_asym_mean:
    fail(f'(2) rho mirror asymmetry mean {asym_mean:.3e} > {thr_asym_mean} '
         '(carbuncle-type shock asymmetry?)')
print('(2) Symmetry: OK')

# ---------------------------------------------------------------------------
# (3) BOW-SHOCK STANDOFF sanity (Billig 1967 cylinder correlation +-40%)
# ---------------------------------------------------------------------------
# Stagnation line = the two cell columns straddling the i-midline. All
# geometry is derived from the mesh itself: wall at j=0; the wall's first
# and last nodes are diametrically opposite => R_cyl = half their distance.
R_cyl = 0.5 * math.hypot(xn[-1, 0] - xn[0, 0], yn[-1, 0] - yn[0, 0])
i_lo, i_hi = ni // 2 - 1, ni // 2
rho_st = 0.5 * (rho2d[i_lo, :] + rho2d[i_hi, :])
# cell centres along the stagnation line, distance from the wall point
xc = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
yc = 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:])
xs = 0.5 * (xc[i_lo, :] + xc[i_hi, :])
ys = 0.5 * (yc[i_lo, :] + yc[i_hi, :])
xw, yw = xn[ni // 2, 0], yn[ni // 2, 0]     # wall node on the stagnation line
d = np.sqrt((xs - xw) ** 2 + (ys - yw) ** 2)

R_gas = float(np.median(R))
rho_inf = p1 / (R_gas * T1)
# Density-jump detector level: 2*rho_inf sits between the freestream (1x)
# and the normal-shock post-state (r2/r1 = 5.5 at M=8.1), so the crossing
# unambiguously locates the captured shock regardless of its smearing.
thr_rho = 2.0 * rho_inf
j_cross = None
for j in range(nj - 1, 0, -1):               # walk inward from the farfield
    if rho_st[j] < thr_rho and rho_st[j - 1] >= thr_rho:
        j_cross = j
        break
if j_cross is None:
    fail('(3) no bow-shock density jump found along the stagnation line '
         f'(rho never crosses 2*rho_inf={thr_rho:.4f}): shock missing or detached '
         'outside the domain')
frac = (thr_rho - rho_st[j_cross]) / (rho_st[j_cross - 1] - rho_st[j_cross])
d_shock = float(d[j_cross] + frac * (d[j_cross - 1] - d[j_cross]))
delta_R = d_shock / R_cyl
billig = 0.386 * math.exp(4.67 / (M1 * M1))  # Billig 1967, cylinder
lo, hi = 0.6 * billig, 1.4 * billig
print(f"(3) Bow-shock standoff: Delta={d_shock:.4f} (R_cyl={R_cyl:.4f}) "
      f"Delta/R={delta_R:.4f}  Billig={billig:.4f}  window=[{lo:.4f},{hi:.4f}]")
if not (lo <= delta_R <= hi):
    fail(f'(3) standoff Delta/R={delta_R:.4f} outside Billig +-40% window '
         f'[{lo:.4f},{hi:.4f}]')
print('(3) Bow-shock standoff: OK')

print('Correct!')
