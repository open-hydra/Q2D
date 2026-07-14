#!/usr/bin/env python3
"""Standalone verification of the axisymmetric turbulent chang-nozzle pair
(kw = Menter SST, rsm = SSG/LRR-omega) — SSGLRR 2D-axi running V&V (Phase C).

Mirror of the test-family standalone verify.py pattern (matplotlib, console
gates, PNG figure). Code-verification evidence, NOT validation (no independent
axi RSM reference exists; the model/term/discretisation correctness is
established by tools/audit/{ssglrr_regime_reduction_sympy, ssglrr_term_level_
check, mms_ssglrr_axi} — this script verifies the RUNNING case exercises the
axi-RSM machinery and behaves physically):

  G1  both runs converged (density residual < 1e-8)
  G2  realizability everywhere on the converged RSM field
      (R_aa >= 0 ; |R12| <= sqrt(R11 R22))
  G3  mean-flow consistency RSM vs SST: with the derived inflow turbulence
      (mut/mu ~ 2e-2) the mean flow is turbulence-insensitive, so the axis
      pressure profiles must agree closely (rel. L_inf < 1%)
  G4  the azimuthal-stretching production P_tt = -2 rho R33 u_r/r (the
      axi-specific RSM term, theory §6.2) is dynamically significant:
      max |P_tt| / (rho eps) = O(1) somewhere in the field
      (rho eps = rho C_mu k omega, the local dissipation scale; both
      gradient-free, computed exactly from cell values)
  G5  Reynolds-stress anisotropy develops from the exactly isotropic IC
      (a_ij = R_ij/k - 2/3 delta_ij != 0 where the nozzle strains the
      turbulence) — the RSM discriminator vs the Boussinesq sibling.

Run:  python3 verify.py   (after both kw/ and rsm/ have been solved)
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
C_MU = 0.09


def read_field(case):
    p = HERE / case / "OUTPUT" / "field.tec"
    lines = p.read_text().splitlines()
    names = re.findall(r'"([^"]+)"', lines[0])
    m = re.search(r"I\s*=\s*(\d+).*?J\s*=\s*(\d+)", lines[1])
    ni, nj = int(m.group(1)), int(m.group(2))
    n_nodal, n_cc = ni * nj, (ni - 1) * (nj - 1)
    vals = []
    for ln in lines[2:]:
        vals.extend(float(t) for t in ln.split())
    out, o = {}, 0
    for k, name in enumerate(names):
        n = n_nodal if k < 2 else n_cc
        a = np.array(vals[o:o + n]); o += n
        if k < 2:
            out[name] = a.reshape(nj, ni)
        else:
            out[name] = a.reshape(nj - 1, ni - 1)
    # cell-centre coordinates from the nodal x,y
    xc = 0.25 * (out["x"][:-1, :-1] + out["x"][1:, :-1] +
                 out["x"][:-1, 1:] + out["x"][1:, 1:])
    rc = 0.25 * (out["y"][:-1, :-1] + out["y"][1:, :-1] +
                 out["y"][:-1, 1:] + out["y"][1:, 1:])
    out["xc"], out["rc"] = xc, rc
    return out


def main():
    ok_all = True

    def gate(name, ok, detail):
        nonlocal ok_all
        ok_all &= ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    print("=" * 74)
    print("chang-nozzle-turbulent — SSGLRR 2D-axi running V&V (kw vs rsm)")
    print("=" * 74)

    # G1 convergence
    res = {}
    for case in ("kw", "rsm"):
        h = np.loadtxt(HERE / case / "OUTPUT" / "residual-history.dat")
        res[case] = h
        gate(f"G1 {case} converged", h[-1, 2] < 1e-8,
             f"iter {int(h[-1,0])}, density residual {h[-1,2]:.2e}")

    kw, rsm = read_field("kw"), read_field("rsm")

    # RSM state (mass-specific R_ij and omega; field.tec writes them /rho? ->
    # the solver outputs primitive RANS vars rho*Rij / rho = mass-specific:
    # names ruu rvv rww ruv omega; density is rho(1))
    rho = rsm["rho(1)"]
    R11, R22, R33, R12 = rsm["ruu"], rsm["rvv"], rsm["rww"], rsm["ruv"]
    omg = rsm["omega"]
    k = 0.5 * (R11 + R22 + R33)

    # G2 realizability
    min_norm = min(R11.min(), R22.min(), R33.min())
    schwarz = np.max(np.abs(R12) - np.sqrt(np.maximum(R11 * R22, 0)))
    gate("G2 realizability", (min_norm >= 0) and (schwarz <= 1e-12),
         f"min R_aa = {min_norm:.3e} >= 0 ; max(|R12|-sqrt(R11 R22)) = {schwarz:.3e} <= 0")

    # G3 mean-flow consistency along the axis row (j=0)
    p_kw, p_rsm = kw["p"][0, :], rsm["p"][0, :]
    rel = np.max(np.abs(p_rsm - p_kw) / np.abs(p_kw))
    gate("G3 mean flow rsm == kw (axis p)", rel < 0.01,
         f"max rel diff = {rel:.3e} < 1e-2")

    # Significance mask: cells whose k is dynamically meaningful. Where the
    # supersonic expansion has decayed the stresses toward the realizability
    # floor (1e-16), a_ij = R_ij/k - 2/3 is a 0/0 of floored values and the
    # production/dissipation ratio is noise — measured on the significant
    # region AND reported unmasked, not hidden.
    sig = k > 1e-2 * k.max()

    # G4 azimuthal-stretching production P_tt = -2 rho R33 v/r vs rho*eps
    v = rsm["v"]
    rc = rsm["rc"]
    P33 = -2.0 * rho * R33 * v / rc
    rhoeps = rho * C_MU * k * omg
    ratio = np.abs(P33) / np.maximum(rhoeps, 1e-300)
    rs = ratio[sig]
    gate("G4 P_tt exercised", rs.max() > 0.5,
         f"max |P_tt|/(rho eps) = {rs.max():.2f} on the significant-k region "
         f"(unmasked max {ratio.max():.0f} in floored-k cells, disclosed)")

    # G5 anisotropy development from the exactly isotropic IC
    a11 = R11 / k - 2.0 / 3.0
    a22 = R22 / k - 2.0 / 3.0
    a33 = R33 / k - 2.0 / 3.0
    a12 = R12 / k
    amax = max(np.abs(a11[sig]).max(), np.abs(a22[sig]).max(),
               np.abs(a33[sig]).max())
    amax_full = max(np.abs(a11).max(), np.abs(a22).max(), np.abs(a33).max())
    gate("G5 anisotropy develops", amax > 0.05,
         f"max |a_ii| = {amax:.3f} (significant-k region; unmasked "
         f"{amax_full:.3f} -> one-component limit reached only in floored-k "
         f"cells), max |a12| = {np.abs(a12[sig]).max():.3f}; IC was exactly isotropic")

    # ------------------------------------------------------------------
    # figure
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2, 2, figsize=(12, 8))
        x_ax = rsm["xc"][0, :]

        ax = axs[0, 0]
        ax.plot(x_ax, p_kw / 1e5, "b-", label="SST (kw)")
        ax.plot(x_ax, p_rsm / 1e5, "r--", label="SSG/LRR-omega (rsm)")
        ax.set_xlabel("x [m]"); ax.set_ylabel("p [bar] (axis)")
        ax.set_title("G3 mean flow: axis pressure (overlay)")
        ax.legend(); ax.grid(alpha=0.3)

        ax = axs[0, 1]
        jmid = a11.shape[0] // 2
        for arr, lbl in ((a11, "a11 (axial)"), (a22, "a22 (radial)"),
                         (a33, "a33 (azimuthal)"), (a12, "a12 (shear)")):
            ax.plot(rsm["xc"][jmid, :], arr[jmid, :], label=lbl)
        ax.set_xlabel("x [m]"); ax.set_ylabel("a_ij = R_ij/k - 2/3 d_ij")
        ax.set_title("G5 anisotropy at mid-radius (IC was isotropic)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axs[1, 0]
        pc = ax.pcolormesh(rsm["xc"], rsm["rc"], np.log10(np.maximum(ratio, 1e-12)),
                           shading="auto", cmap="magma")
        fig.colorbar(pc, ax=ax, label="log10 |P_tt| / (rho eps)")
        ax.set_xlabel("x [m]"); ax.set_ylabel("r [m]")
        ax.set_title("G4 azimuthal-stretching production -2 rho R33 v/r")

        ax = axs[1, 1]
        ax.semilogy(res["kw"][:, 0], res["kw"][:, 2], "b-", label="SST")
        ax.semilogy(res["rsm"][:, 0], res["rsm"][:, 2], "r--", label="SSGLRR")
        ax.set_xlabel("iteration"); ax.set_ylabel("density residual")
        ax.set_title("G1 convergence (identical numerics, CFL 0.8)")
        ax.legend(); ax.grid(alpha=0.3)

        fig.suptitle("chang-nozzle-turbulent: SSG/LRR-omega 2D-axi running V&V")
        fig.tight_layout()
        out = HERE / "rsm" / "OUTPUT" / "verify_RSM_axi.png"
        fig.savefig(out, dpi=130)
        print(f"  figure: {out}")
    except Exception as e:  # figure is best-effort; gates are the record
        print(f"  (figure skipped: {e})")

    print("=" * 74)
    print("RESULT:", "ALL GATES PASS" if ok_all else "FAILURES PRESENT")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
