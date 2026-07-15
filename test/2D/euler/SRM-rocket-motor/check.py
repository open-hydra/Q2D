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

Gates (all hard):
  (1) Vieille chamber pressure: |P_chamber_avg - P_eq_analytical| /
      P_eq_analytical < 0.15. Acceptance tolerance per piano §5 Fase 8c
      (analogous chang-nozzle pattern ±5% avg, ±15% max relaxed for
      spatial-vs-lumped discretization residual AND for the residual
      pressurization transient: at the 5000-iteration ctest budget the
      chamber is still slowly filling toward the asymptote -- measured
      eps_rel = 9.4e-2 on 2026-07-15, gfortran TESTING build).
      Reads reference/P_eq_vieille.dat for the analytical reference value.
  (2) Convergence sentinel: final density residual < 0.5 (L2 norm of the
      per-iteration density change, residual-history.dat column 2; NaN-safe).
      The run is quasi-steady, not deep-converged: the measured final value
      is 2.6e-1 (slow monotonic chamber pressurization). The threshold is
      ~2x the measured final and catches divergence/NaN blow-up.
  (3) NaN-safe positivity: rho, p, T, g, R > 0; u, v finite.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (assert_converged, assert_finite, assert_positive,
                           fail, gate, read_q2d_field, require_file)

require_file("OUTPUT/field.tec", "solver did not write the field")
require_file("MESH/profile.dat", "mesh profile missing")
require_file("reference/P_eq_vieille.dat", "analytical Vieille reference missing")

x, y, var, names = read_q2d_field("OUTPUT/field.tec")
cc = [n.strip() for n in names[3:]]  # cell-centered variables (after x, y, z)
missing = [n for n in ("rho(1)", "u", "v", "p", "T", "g", "R") if n not in cc]
if missing:
    fail("missing variables in OUTPUT/field.tec: %s (found %s)" % (missing, cc))


def field(name):
    return np.asarray(var[0][cc.index(name)], float)[:, :, 0]


xn = x[0][:, :, 0]
p = field("p")

# (3) NaN-safe positivity / finiteness
for name in ("rho(1)", "p", "T", "g", "R"):
    assert_positive(name, field(name))
for name in ("u", "v"):
    assert_finite(name, field(name))
print("(3) positivity/finiteness OK: rho, p, T, g, R > 0; u, v finite")

# (2) convergence sentinel (density residual, column 2; finite + below threshold)
print("(2) convergence sentinel:")
assert_converged(0.5, column=2)

# (1) Vieille chamber-pressure balance
print("(1) Vieille chamber pressure:")
profile = np.loadtxt("MESH/profile.dat")
x_throat = profile[np.argmin(profile[:, 1]), 0]
print("  throat location: x=%.4f m" % x_throat)

xc = 0.5 * (xn[:-1, :-1] + xn[1:, 1:])
mask_chamber = xc[:, 0] < (x_throat - 0.005)
n_cells_chamber = int(mask_chamber.sum())
print("  chamber cells (x < x_throat - 0.005): %d / %d" % (n_cells_chamber, xc.shape[0]))
if n_cells_chamber < 5:
    fail("fewer than 5 chamber cells; check geometry")

p_chamber = float(p[mask_chamber, :].mean())
p_chamber_std = float(p[mask_chamber, :].std())
print("  P_chamber (mean over chamber cells): %.4e Pa = %.4f MPa"
      % (p_chamber, p_chamber / 1e6))
print("  P_chamber std-dev: %.4e Pa (%.2f %%)"
      % (p_chamber_std, 100 * p_chamber_std / p_chamber))

ref_data = np.loadtxt("reference/P_eq_vieille.dat", comments="#")
P_eq_ref = float(ref_data if ref_data.ndim == 0 else ref_data[0])
print("  P_eq Vieille analytical: %.4e Pa = %.4f MPa" % (P_eq_ref, P_eq_ref / 1e6))

gate("Vieille rel error |P_chamber - P_eq| / P_eq",
     abs(p_chamber - P_eq_ref) / P_eq_ref, 0.15)

print("Correct!")
