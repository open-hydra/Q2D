"""rocket-nozzle/frozen 2D-axisymmetric multi-species (13) Euler steady-state regression.

Validates the Q2D converged solution for the rocket-nozzle frozen-chemistry
case (CEA-equilibrium upstream stagnation, frozen evolution through nozzle).
Single-block 121x41 2D-axi (axisymmetric source active). 13 species from
chamber composition.

Gates (all hard):
  (1) Cell-by-cell field comparison vs the committed reference/field.tec
      (3D-wedge K=2 snapshot of the same physics; the Q2D 2D-axi K=1 strip
      is compared on the single common cell layer, 120x40). Per-variable
      tolerances Linf_rel < 5e-2, L2_rel < 1e-2. Measured 2026-07-15
      (gfortran TESTING build): worst Linf_rel = 1.36e-2 (v near axis),
      worst L2_rel = 1.63e-3; all species ~1.0e-3 Linf_rel.
  (2) Species/ideal-gas consistency: max_cells |p - sum_s rho_s*R*T| / p
      < 1e-8 (no total-rho variable is written, so the mixture gas law is
      the available closure; measured 9.3e-15 -- tec ASCII roundoff only).
  (3) NaN-safe positivity: every species density, p, T, g, R > 0 strictly;
      u, v finite (v is legitimately negative toward the axis).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "common" / "pycheck"))
from q2d_check_lib import (assert_finite, assert_positive, compare_fields,
                           fail, gate, read_q2d_field, require_file)

N_SPECIES = 13
OUT = "OUTPUT/field.tec"
REF = "reference/field.tec"
require_file(OUT, "solver did not write the field")
require_file(REF, "committed reference field missing")

x, y, var, names = read_q2d_field(OUT)
cc = [n.strip() for n in names[3:]]  # cell-centered variables (after x, y, z)
species = ["rho(%d)" % (s + 1) for s in range(N_SPECIES)]
missing = [n for n in species + ["u", "v", "p", "T", "g", "R"] if n not in cc]
if missing:
    fail("missing variables in %s: %s (found %s)" % (OUT, missing, cc))


def field(name):
    return np.asarray(var[0][cc.index(name)], float)


# (3) NaN-safe positivity / finiteness for ALL variables
for name in species + ["p", "T", "g", "R"]:
    assert_positive(name, field(name))
for name in ("u", "v"):
    assert_finite(name, field(name))
print("(3) positivity/finiteness OK: %d species densities, p, T, g, R > 0; u, v finite"
      % N_SPECIES)

# (2) species/ideal-gas consistency: p == (sum_s rho_s) * R_mix * T
rho_sum = sum(field(s) for s in species)
p, T, R = field("p"), field("T"), field("R")
cons = float(np.max(np.abs(p - rho_sum * R * T) / p))
print("(2) species consistency:")
gate("max |p - sum(rho_s)*R*T| / p", cons, 1e-8)

# (1) committed reference comparison (hard per-variable gates)
print("(1) reference comparison:")
compare_fields(OUT, REF, tol_linf=5e-2, tol_l2=1e-2)

# Diagnostic report (not a gate): Mach on the axis at the throat
u, v, g = field("u"), field("v"), field("g")
i_throat = int(np.argmin(np.abs(x[0][:, 0, 0])))
M_axis = (np.sqrt(u[i_throat, 0, 0] ** 2 + v[i_throat, 0, 0] ** 2)
          / np.sqrt(g[i_throat, 0, 0] * R[i_throat, 0, 0] * T[i_throat, 0, 0]))
print("  throat at i=%d (x=%.3f); Mach axis: %.3f  p=%.3e  T=%.1f"
      % (i_throat, x[0][i_throat, 0, 0], M_axis, p[i_throat, 0, 0], T[i_throat, 0, 0]))

print("Correct!")
