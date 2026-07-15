"""supersonic-forward-step (Woodward-Colella 1984) 2D time-accurate regression.

Validates the Q2D solution at the deterministic endpoint t=4.0 for the
Woodward-Colella forward-step benchmark (M=3 inflow into a channel with a
0.2-high step at x=0.6; 2-block grid: block 1 pre-step 48x16 cells, block 2
above-step 240x64 cells).

Gates (all hard):
  (1) Density trace vs the committed baseline reference/trace_baseline.dat
      (rho along block-2 cell row j=23, nearest y=0.5, through the
      bow-shock/Mach-stem region; provenance in reference/README.md:
      gfortran TESTING build, main @ 8c60e91, 2026-07-15). Max relative
      difference per point < 5% -- single-thread reruns of the same binary
      are bit-identical (measured 0.0), so 5% is pure cross-compiler
      headroom. The trace's x cell centers must match the baseline's to
      1e-9 (mesh guard).
  (2) Bow-shock x-position along the same row, located by the steepest
      density gradient: within 3 cells (|dx| <= 0.0375) of the baseline's
      (measured x_shock = 0.425 at t=4.0).
  (3) NaN-safe positivity on every block: rho, p, T, g, R > 0; u, v finite.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import (assert_finite, assert_positive, fail, gate,
                           read_q2d_field, require_file)

OUT = "OUTPUT/field.tec"
BASELINE = "reference/trace_baseline.dat"
TRACE_ZONE = 1        # block 2: above-step, 240x64 cells
TRACE_Y = 0.5         # extract the cell row nearest this height
TOL_TRACE = 0.05      # max relative difference per point
TOL_SHOCK_X = 0.0375  # 3 cells of dx=0.0125

require_file(OUT, "solver did not write the field")
require_file(BASELINE, "committed trace baseline missing")

x_, y_, var_, names = read_q2d_field(OUT)
cc = [n.strip() for n in names[3:]]
missing = [n for n in ("rho(1)", "u", "v", "p", "T", "g", "R") if n not in cc]
if missing:
    fail("missing variables in %s: %s (found %s)" % (OUT, missing, cc))
if len(var_) < 2:
    fail("expected 2 zones in %s, found %d" % (OUT, len(var_)))

# (3) NaN-safe positivity / finiteness on every block
for b in range(len(var_)):
    for name in ("rho(1)", "p", "T", "g", "R"):
        assert_positive("block%d %s" % (b + 1, name),
                        np.asarray(var_[b][cc.index(name)], float))
    for name in ("u", "v"):
        assert_finite("block%d %s" % (b + 1, name),
                      np.asarray(var_[b][cc.index(name)], float))
print("(3) positivity/finiteness OK on all %d blocks: rho, p, T, g, R > 0; u, v finite"
      % len(var_))

# ---------------------------------------------------------------------------
# Extract the density trace: block-2 cell row nearest y=TRACE_Y
# ---------------------------------------------------------------------------
xn = x_[TRACE_ZONE][:, :, 0]
yn = y_[TRACE_ZONE][:, :, 0]
rho = np.asarray(var_[TRACE_ZONE][cc.index("rho(1)")], float)[:, :, 0]
xc = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
yc = 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:])
j = int(np.argmin(np.abs(yc[0, :] - TRACE_Y)))
trace_x = xc[:, j]
trace_rho = rho[:, j]

base = np.loadtxt(BASELINE)
if base.ndim != 2 or base.shape[1] != 2:
    fail("baseline '%s' is not a 2-column (x, rho) table" % BASELINE)
if base.shape[0] != trace_x.size:
    fail("trace length %d != baseline length %d (mesh or extraction row changed)"
         % (trace_x.size, base.shape[0]))
dx_mesh = float(np.abs(np.diff(base[:, 0])).mean())
if float(np.abs(trace_x - base[:, 0]).max()) > 1e-9:
    fail("trace x cell centers differ from baseline (mesh changed?); "
         "max |dx| = %.3e" % float(np.abs(trace_x - base[:, 0]).max()))

# (1) trace comparison
print("(1) density trace vs committed baseline (zone 2, row j=%d, y=%.4f):"
      % (j, yc[0, j]))
rel = np.abs(trace_rho - base[:, 1]) / np.abs(base[:, 1])
print("    trace points: %d  mean rel diff: %.3e" % (rel.size, float(rel.mean())))
gate("max rel density diff along trace", float(rel.max()), TOL_TRACE)

# (2) bow-shock x-position by steepest density gradient
i_run = int(np.argmax(np.abs(np.diff(trace_rho))))
i_ref = int(np.argmax(np.abs(np.diff(base[:, 1]))))
x_run = 0.5 * (trace_x[i_run] + trace_x[i_run + 1])
x_ref = 0.5 * (base[i_ref, 0] + base[i_ref + 1, 0])
print("(2) bow-shock x-position (steepest density gradient along the row):")
print("    run x_shock=%.4f  baseline x_shock=%.4f  (dx_cell=%.4f)"
      % (x_run, x_ref, dx_mesh))
gate("|x_shock - baseline|", abs(x_run - x_ref), TOL_SHOCK_X)

print("Correct!")
