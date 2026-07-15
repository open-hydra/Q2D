#!/usr/bin/env python3
"""Standalone verification of the Q2D AREA-VARIATION turbulent planar-nozzle pair
(kw = Menter SST, rsm = SSG/LRR-omega) — SSGLRR area-variation running V&V
(Phase C). The area-variation sibling of the 2D-axi chang-nozzle-turbulent
verify.py.

Mirror of the test-family standalone verify.py pattern (matplotlib, console
gates, PNG figure). Code-verification evidence, NOT validation (no independent
area-var RSM reference exists; model/term/discretisation correctness is
established by {ssglrr_regime_reduction_sympy, ssglrr_term_level_
check, mms_ssglrr_areavar} — this script verifies the RUNNING case exercises
the area-var-RSM machinery and behaves physically):

  G1  both runs converged (density residual < 1e-8)
  G2  realizability everywhere on the converged RSM field
      (R_aa >= 0 ; |R12| <= sqrt(R11 R22))
  G3  mean-flow consistency RSM vs SST: with the derived inflow turbulence
      (mut/mu ~ 2e-2) the mean flow is turbulence-insensitive, so the centreline
      pressure profiles must agree closely (rel. L_inf < 1%)
  G4  the depth-compression production P_33 = +2 rho R33 (V.G) (the area-var
      specific RSM term, theory §6.2 / A.5, GRADIENT-FREE) is dynamically
      significant: max |P_33| / (rho eps) = O(1) somewhere in the field
      (rho eps = rho C_mu k omega, the local dissipation scale; both
      gradient-free, computed exactly from cell values). G = -grad(h)/h from the
      committed area law block1_area.dat (the same h(x) the solver uses).
  G5  Reynolds-stress anisotropy develops from the exactly isotropic IC
      (a_ij = R_ij/k - 2/3 delta_ij != 0 where the nozzle strains the
      turbulence) — the RSM discriminator vs the Boussinesq sibling.

Run:  python3 verify.py   (after both kw/ and rsm/ have been solved)

Exit semantics: every gate G1-G5 is HARD — any violation -> exit 1; exit 0
only when all five gates ran and passed. Missing/empty required inputs
({kw,rsm}/OUTPUT/residual-history.dat, {kw,rsm}/OUTPUT/field.tec,
rsm/INPUT/block1_area.dat) -> exit 1 with an explicit FAIL message (never a
silent pass). All parsed fields are NaN/Inf-guarded via the shared
q2d_check_lib.assert_finite before the gates run. The PNG figure is the only
optional artifact: it is best-effort and never affects the exit code.
"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import fail, require_file, assert_finite

HERE = Path(__file__).resolve().parent
C_MU = 0.09


def read_field(case):
    p = HERE / case / "OUTPUT" / "field.tec"
    require_file(p, "solver did not write the {0} field output".format(case))
    lines = p.read_text().splitlines()
    names = re.findall(r'"([^"]+)"', lines[0])
    if len(lines) < 2:
        fail("{0}: truncated (no ZONE record)".format(p))
    m = re.search(r"I\s*=\s*(\d+).*?J\s*=\s*(\d+)", lines[1])
    if m is None:
        fail("{0}: no I=/J= zone record (truncated output?)".format(p))
    ni, nj = int(m.group(1)), int(m.group(2))
    n_nodal, n_cc = ni * nj, (ni - 1) * (nj - 1)
    vals = []
    for ln in lines[2:]:
        vals.extend(float(t) for t in ln.split())
    out, o = {}, 0
    if len(vals) < 2 * n_nodal + (len(names) - 2) * n_cc:
        fail("{0}: {1} values parsed, expected {2} (truncated output?)".format(
            p, len(vals), 2 * n_nodal + (len(names) - 2) * n_cc))
    for k, name in enumerate(names):
        n = n_nodal if k < 2 else n_cc
        a = np.array(vals[o:o + n]); o += n
        assert_finite("{0} '{1}'".format(p, name), a)
        if k < 2:
            out[name] = a.reshape(nj, ni)
        else:
            out[name] = a.reshape(nj - 1, ni - 1)
    xc = 0.25 * (out["x"][:-1, :-1] + out["x"][1:, :-1] +
                 out["x"][:-1, 1:] + out["x"][1:, 1:])
    yc = 0.25 * (out["y"][:-1, :-1] + out["y"][1:, :-1] +
                 out["y"][:-1, 1:] + out["y"][1:, 1:])
    out["xc"], out["yc"] = xc, yc
    return out, ni, nj


def read_area(ni, nj):
    """Nodal area law h(x,y) from block1_area.dat ((A(i,j),i=0,ni),j=0,nj)."""
    p = HERE / "rsm" / "INPUT" / "block1_area.dat"
    require_file(p, "committed area-law input missing (required for G4)")
    raw = p.read_text().split()
    if len(raw) != ni * nj:
        fail("{0}: {1} values, expected ni*nj = {2}".format(p, len(raw), ni * nj))
    A = np.array([float(t) for t in raw]).reshape(nj, ni)   # j outer, i inner
    assert_finite(str(p), A)
    return A


def area_G(field, A):
    """G = -grad(h)/h at cell centres, central difference on the nodal area law
    projected to cell centres (h_cc = 4-node average), matching the code's
    Precompute_Area_Geometry sign convention G_x = -dA/dx / A_center."""
    # cell-centre h and its central-difference gradient on the cell grid
    h_cc = 0.25 * (A[:-1, :-1] + A[1:, :-1] + A[:-1, 1:] + A[1:, 1:])
    xc, yc = field["xc"], field["yc"]
    Gx = np.zeros_like(h_cc); Gy = np.zeros_like(h_cc)
    Gx[:, 1:-1] = -(h_cc[:, 2:] - h_cc[:, :-2]) / (xc[:, 2:] - xc[:, :-2]) / h_cc[:, 1:-1]
    Gy[1:-1, :] = -(h_cc[2:, :] - h_cc[:-2, :]) / (yc[2:, :] - yc[:-2, :]) / h_cc[1:-1, :]
    return Gx, Gy


def main():
    ok_all = True

    def gate(name, ok, detail):
        nonlocal ok_all
        ok_all &= ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    print("=" * 74)
    print("planar-nozzle-turbulent — SSGLRR area-variation running V&V (kw vs rsm)")
    print("=" * 74)

    res = {}
    for case in ("kw", "rsm"):
        hp = HERE / case / "OUTPUT" / "residual-history.dat"
        require_file(hp, "solver did not write the {0} residual history".format(case))
        h = np.loadtxt(hp)
        if h.ndim == 1:
            h = h.reshape(1, -1)
        assert_finite("{0} residual history".format(case), h)
        res[case] = h
        gate(f"G1 {case} converged", h[-1, 2] < 1e-8,
             f"iter {int(h[-1,0])}, density residual {h[-1,2]:.2e}")

    (kw, _, _), (rsm, ni, nj) = read_field("kw"), read_field("rsm")

    rho = rsm["rho(1)"]
    R11, R22, R33, R12 = rsm["ruu"], rsm["rvv"], rsm["rww"], rsm["ruv"]
    omg = rsm["omega"]
    k = 0.5 * (R11 + R22 + R33)

    min_norm = min(R11.min(), R22.min(), R33.min())
    schwarz = np.max(np.abs(R12) - np.sqrt(np.maximum(R11 * R22, 0)))
    gate("G2 realizability", (min_norm >= 0) and (schwarz <= 1e-12),
         f"min R_aa = {min_norm:.3e} >= 0 ; max(|R12|-sqrt(R11 R22)) = {schwarz:.3e} <= 0")

    # G3 mean-flow consistency along the centreline row (mid-j)
    jc = kw["p"].shape[0] // 2
    p_kw, p_rsm = kw["p"][jc, :], rsm["p"][jc, :]
    rel = np.max(np.abs(p_rsm - p_kw) / np.abs(p_kw))
    gate("G3 mean flow rsm == kw (centreline p)", rel < 0.01,
         f"max rel diff = {rel:.3e} < 1e-2")

    sig = k > 1e-2 * k.max()

    # G4 depth-compression production P_33 = +2 rho R33 (V.G), V.G = u Gx + v Gy
    A = read_area(ni, nj)
    Gx, Gy = area_G(rsm, A)
    u, v = rsm["u"], rsm["v"]
    VG = u * Gx + v * Gy
    P33 = 2.0 * rho * R33 * VG
    rhoeps = rho * C_MU * k * omg
    ratio = np.abs(P33) / np.maximum(rhoeps, 1e-300)
    rs = ratio[sig]
    gate("G4 P_33 exercised", rs.max() > 0.5,
         f"max |P_33|/(rho eps) = {rs.max():.2f} on the significant-k region "
         f"(unmasked max {ratio.max():.0f} in floored-k cells, disclosed)")

    a11 = R11 / k - 2.0 / 3.0
    a22 = R22 / k - 2.0 / 3.0
    a33 = R33 / k - 2.0 / 3.0
    a12 = R12 / k
    amax = max(np.abs(a11[sig]).max(), np.abs(a22[sig]).max(),
               np.abs(a33[sig]).max())
    amax_full = max(np.abs(a11).max(), np.abs(a22).max(), np.abs(a33).max())
    gate("G5 anisotropy develops", amax > 0.05,
         f"max |a_ii| = {amax:.3f} (significant-k region; unmasked "
         f"{amax_full:.3f}), max |a12| = {np.abs(a12[sig]).max():.3f}; "
         "IC was exactly isotropic")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2, 2, figsize=(12, 8))
        x_ax = rsm["xc"][jc, :]

        ax = axs[0, 0]
        ax.plot(x_ax, p_kw / 1e5, "b-", label="SST (kw)")
        ax.plot(x_ax, p_rsm / 1e5, "r--", label="SSG/LRR-omega (rsm)")
        ax.set_xlabel("x [m]"); ax.set_ylabel("p [bar] (centreline)")
        ax.set_title("G3 mean flow: centreline pressure (overlay)")
        ax.legend(); ax.grid(alpha=0.3)

        ax = axs[0, 1]
        for arr, lbl in ((a11, "a11"), (a22, "a22"),
                         (a33, "a33 (depth)"), (a12, "a12 (shear)")):
            ax.plot(rsm["xc"][jc, :], arr[jc, :], label=lbl)
        ax.set_xlabel("x [m]"); ax.set_ylabel("a_ij = R_ij/k - 2/3 d_ij")
        ax.set_title("G5 anisotropy at centreline (IC was isotropic)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axs[1, 0]
        pc = ax.pcolormesh(rsm["xc"], rsm["yc"], np.log10(np.maximum(ratio, 1e-12)),
                           shading="auto", cmap="magma")
        fig.colorbar(pc, ax=ax, label="log10 |P_33| / (rho eps)")
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        ax.set_title("G4 depth-compression production +2 rho R33 (V.G)")

        ax = axs[1, 1]
        ax.semilogy(res["kw"][:, 0], res["kw"][:, 2], "b-", label="SST")
        ax.semilogy(res["rsm"][:, 0], res["rsm"][:, 2], "r--", label="SSGLRR")
        ax.set_xlabel("iteration"); ax.set_ylabel("density residual")
        ax.set_title("G1 convergence (identical numerics)")
        ax.legend(); ax.grid(alpha=0.3)

        fig.suptitle("planar-nozzle-turbulent: SSG/LRR-omega area-variation running V&V")
        fig.tight_layout()
        out = HERE / "rsm" / "OUTPUT" / "verify_RSM_areavar.png"
        fig.savefig(out, dpi=130)
        print(f"  figure: {out}")
    except Exception as e:
        print(f"  (figure skipped: {e})")

    print("=" * 74)
    print("RESULT:", "ALL GATES PASS" if ok_all else "FAILURES PRESENT")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
