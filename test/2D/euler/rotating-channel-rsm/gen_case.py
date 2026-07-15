#!/usr/bin/env python3
"""Generate rotating-channel-rsm inputs from the committed Euler parent.

Derivation of every value (no invented numbers):
  - mesh / bc.txt / phase.txt / thermo.dat: VERBATIM the committed
    `test/2D/euler/rotating-channel` (periodic 201 in x, symmetry 300 in y,
    omega = 10 rad/s, U_inf = 100 m/s);
  - transport.dat: the committed `planar-nozzle-turbulent/rsm` file (same
    single-species air phase) — SSGLRR activates navier-stokes and the
    transport reader;
  - turbulence IC (uniform in the CONSERVED variables, appended as 5 extra
    cell-centred columns to the parent ic.tec):
      rho*R11 = 2.0, rho*R22 = 1.0, rho*R33 = 1.5  [k/U_inf^2 ~ 2e-4:
        one-way coupling to the mean flow],
      rho*R12 = 0                                   [compatible with the
        symmetry walls, which pin R12 = 0; the discriminating signal R12(t)
        is then generated purely by G12 = 2*Omega*(R22 - R11)],
      rho*omega_t = 400                             [Cmu*omega_t*C1*T ~ 3:
        the slow pressure-strain acts visibly on the T = 0.05 s window].

Gates and reference ODE: verify.py (frame-invariance ODE proof). Figure: plot.py.
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "rotating-channel"

R11, R22, R33, R12, OMEGA_T = 2.0, 1.0, 1.5, 0.0, 400.0

lines = (SRC / "INPUT" / "ic.tec").read_text().split("\n")
m = re.search(r"I=(\d+), J=(\d+)", lines[1])
I, J = int(m.group(1)), int(m.group(2))
nnod, ncell = I*J, (I-1)*(J-1)
vals = [float(t) for ln in lines[2:] for t in ln.split()]
assert len(vals) == 2*nnod + 4*ncell, (len(vals), 2*nnod + 4*ncell)

with open(HERE / "INPUT" / "ic.tec", "w") as f:
    f.write(' VARIABLES ="x" "y"  "rho1" "u" "v" "p"  "ruu"  "rvv"  "rww"  "ruv"  "omega"\n')
    f.write(' ZONE  T = B1-IG, I=%d, J=%d, K=1, DATAPACKING=BLOCK, '
            'VARLOCATION=([1-2]=NODAL,[3-11]=CELLCENTERED), '
            'SOLUTIONTIME=0.000000000000000E+000\n' % (I, J))
    for v in vals:
        f.write(' %.15E\n' % v)
    for c in (R11, R22, R33, R12, OMEGA_T):
        for _ in range(ncell):
            f.write(' %.15E\n' % c)

(HERE / "INPUT" / "bc.txt").write_text((SRC / "INPUT" / "bc.txt").read_text())
for name in ("phase.txt", "thermo.dat"):
    (HERE / "INPUT" / name).write_text((SRC / "INPUT" / name).read_text())
(HERE / "INPUT" / "transport.dat").write_text(
    (HERE.parent / "planar-nozzle-turbulent" / "rsm" / "INPUT" / "transport.dat").read_text())
(HERE / "MESH").mkdir(exist_ok=True)
(HERE / "MESH" / "mesh.tec").write_text((SRC / "MESH" / "mesh.tec").read_text())
print("rotating-channel-rsm inputs generated: ic.tec %dx%d (9 cell vars), "
      "bc/phase/thermo/transport/mesh copied" % (I, J))
