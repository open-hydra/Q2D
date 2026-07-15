"""swbli 2D viscous RANS-SA acceptance regression check.

History:
  Sess 8 (2026-05-11): smoke infrastructure (positivity + NaN check).
  Sess 12.A (2026-05-22): Cf wall extraction + comparison vs MOSE-SA.tec
    reference (±10% engineering V&V tolerance).
  Sess 39 (2026-06-15, F13): independent-reference V&V on the separation
    bubble — the shock/boundary-layer interaction's defining physics. The
    separation and reattachment locations (where the wall shear tauX changes
    sign) are a nondimensional, scale-invariant V&V metric (independent of the
    tauX->Cf normalization). They are bracketed against the independent
    references already in reference/ (Schulein-2002 experiment + SU2-SA CFD +
    Wind-US tunnel), turning the case from consistency-only (vs MOSE) into an
    independent-reference gate (closes F13).

Acceptance criteria:
  - Smoke: OUTPUT/field*.tec exists, no NaN/Inf, n_tokens > 0.
  - Cf consistency (vs MOSE-SA.tec, Sess 12.A): max rel diff of tauX <= 0.1%
    over the developed-plate range (MOSE-shared SA -> essentially bit-tight).
  - Separation-bubble independent-reference V&V (F13, Sess 39): the solver's
    separation x_sep and reattachment x_reatt (tauX sign changes in the
    interaction zone) must fall inside the envelope of the independent
    references within a 1.5% engineering band. The envelope (separation
    [0.3337, 0.3375], reattachment [0.3430, 0.3458]) is derived from the
    zero-Cf crossings of schulein.dat / SU2-SA.dat / wind-SA.dat (no invented
    numbers). The 1.5% band is set by the inter-reference spread itself (the
    references disagree by ~1.1% on x_sep — the irreducible physical
    uncertainty), NOT tuned to pass; SA is known to slightly delay separation
    and Q2D lands 0.44% from the nearest CFD-SA peer (SU2). A run that fails to
    separate, or separates outside the band, FAILS (built-in negative control).
    Committed converged run (SOLUTIONTIME=1e5): x_sep=0.3390, x_reatt=0.3435
    (== MOSE bit-coincident).

Exit semantics (2026-07 hardening; replaces the old "SA-INFORMATIVE mode"
that silently PASSED smoke-only when the quantitative gates could not run):
  - exit 0  : smoke + BOTH quantitative gates ran and passed.
  - exit 1  : any failure — missing OUTPUT/field*.tec, missing OUTPUT/wall.tec,
    missing committed reference files, NaN/Inf in the field, parse errors, or
    a quantitative-gate violation. A missing solver output is always a
    failure (Q2D.sh wipes OUTPUT/ before each solve).
  - exit 77 : SKIP (CTest SKIP_RETURN_CODE — reported "Skipped", never
    "Passed") ONLY when wall.tec exists but SOLUTIONTIME < 95000, i.e. the
    run is genuinely under-converged for the engineering V&V gates.

References:
  - MOSE-SA.tec: extracted from MOSE upstream commit 6d5a3a7 Sess 11.A
    (consistency reference, same SA model).
  - schulein.dat: Schulein 2002 experimental skin-friction (independent).
  - SU2-SA.dat: SU2 Spalart-Allmaras CFD (independent code).
  - wind-SA.dat: Wind-US tunnel/CFD reference (independent).
"""
import os
import sys
import glob
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common" / "pycheck"))
from q2d_check_lib import fail, skip, require_file, assert_finite


def smoke_check():
    """Field-output smoke gate: missing output or NaN/Inf -> fail() (exit 1)."""
    out_dir = "OUTPUT"
    if not os.path.isdir(out_dir):
        fail("OUTPUT/ missing (solver did not run)")

    candidates = glob.glob(os.path.join(out_dir, "field*.tec"))
    if not candidates:
        fail("no field*.tec in OUTPUT/ (solver produced no field output)")

    field_file = max(candidates, key=os.path.getmtime)
    require_file(field_file)
    print("[swbli check] inspecting {0}".format(field_file))

    values = []
    with open(field_file, "r") as fh:
        for line in fh:
            tok = line.strip()
            if not tok or tok.startswith("VARIABLES") or tok.startswith("ZONE"):
                continue
            for piece in tok.split():
                try:
                    values.append(float(piece))
                except ValueError:
                    continue

    print("[swbli check] tokens parsed: {0}".format(len(values)))
    if not values:
        fail("no numeric data parsed from {0} (truncated output?)".format(field_file))
    assert_finite("{0} values".format(field_file), values)


def _parse_tecplot_block(path, nodal_var_count, cc_var_count):
    """Parse Tecplot ASCII DATAPACKING=BLOCK file.

    Returns (I, J, K, nodal_arrays, cc_arrays) where:
      nodal_arrays[v_idx] = list of n_nodes floats per nodal var.
      cc_arrays[v_idx] = list of n_cells floats per cc var.
    """
    import re

    with open(path) as fh:
        lines = fh.read().split('\n')

    # find ZONE line + parse I, J, K
    i = 0
    while i < len(lines) and 'ZONE' not in lines[i].upper():
        i += 1
    if i >= len(lines):
        raise RuntimeError("no ZONE line in {0}".format(path))
    zone = lines[i]
    m = re.search(r'I=\s*(\d+).*?J=\s*(\d+).*?K=\s*(\d+)', zone)
    if not m:
        raise RuntimeError("bad ZONE line: {0}".format(zone))
    I, J, K = int(m.group(1)), int(m.group(2)), int(m.group(3))
    i += 1  # skip past ZONE line

    n_nodes = I * J * K
    n_cells = max(I - 1, 1) * max(J - 1, 1) * max(K - 1, 1)

    nodal_arrays = []
    cc_arrays = []

    # read floats sequentially: nodal_var_count blocks of n_nodes, then cc_var_count blocks of n_cells
    floats_iter = _float_stream(lines, i)
    for v in range(nodal_var_count):
        arr = []
        for _ in range(n_nodes):
            arr.append(next(floats_iter))
        nodal_arrays.append(arr)
    for v in range(cc_var_count):
        arr = []
        for _ in range(n_cells):
            arr.append(next(floats_iter))
        cc_arrays.append(arr)
    return I, J, K, nodal_arrays, cc_arrays


def _float_stream(lines, start_idx):
    for ln in lines[start_idx:]:
        s = ln.replace(',', ' ').strip()
        if not s:
            continue
        for tok in s.split():
            try:
                yield float(tok.replace('D', 'E').replace('d', 'e'))
            except ValueError:
                continue


def _cell_centers_x(x_nodal, I, J, K):
    """Compute cell-center x from nodal x (BLOCK k-fastest ordering).

    For 2D with K=1 (Q2D) or K=2 (MOSE K=1-slice): cells along i are between
    consecutive nodal i; for J=1 cell layer, simply average i and i+1.
    """
    # for Q2D K=1 J=1: nodal x[i, j=0, k=0] at index i + 0*I + 0*I*J = i (where i in [0, I-1])
    # We need wall row cells. Wall is at j=0 plate row (I-1 cells).
    # For simplicity: take pairs (x[i], x[i+1])/2 along i for the wall.
    cc = []
    for i in range(I - 1):
        cc.append(0.5 * (x_nodal[i] + x_nodal[i + 1]))
    return cc


def _wall_tec_solutiontime(path):
    """Extract SOLUTIONTIME from Tecplot ZONE line; returns float or None."""
    import re
    try:
        with open(path) as fh:
            for line in fh:
                if 'ZONE' in line.upper():
                    m = re.search(r'SOLUTIONTIME\s*=\s*([\dE\.\+\-]+)', line, re.IGNORECASE)
                    if m:
                        return float(m.group(1).replace('D', 'E').replace('d', 'e'))
                    return None
    except Exception:
        return None
    return None


def cf_check():
    """Compare Q2D wall.tec tauX vs MOSE-SA.tec wall tauX (K=1 slice).

    Returns: (pass_bool, status_msg). The caller has already verified that
    wall.tec + the reference exist and the run is converged; any problem here
    is a hard failure.
    """
    q2d_wall = "OUTPUT/wall.tec"
    mose_ref = "reference/MOSE-SA.tec"

    # Q2D wall.tec: NODAL (x, y) + CC (y+, tauX, tauY, pw, Tw, qw)
    try:
        Iq, Jq, Kq, nq, cq = _parse_tecplot_block(q2d_wall,
                                                  nodal_var_count=2,
                                                  cc_var_count=6)
    except Exception as e:
        return False, "Q2D wall.tec parse error: {0}".format(e)
    # MOSE wall.tec: NODAL (x, y, z) + CC (y+, tauX, tauY, tauZ, pw, Tw, qw)
    try:
        Im, Jm, Km, nm, cm = _parse_tecplot_block(mose_ref,
                                                  nodal_var_count=3,
                                                  cc_var_count=7)
    except Exception as e:
        return False, "MOSE-SA.tec parse error: {0}".format(e)

    if Iq != Im or Jq != Jm:
        return False, "I/J mismatch Q2D ({0},{1}) vs MOSE ({2},{3})".format(
            Iq, Jq, Im, Jm)

    # Q2D nodal x = nq[0], MOSE nodal x = nm[0] first I*J=850 values (K=1 slice)
    q2d_x_nodal = nq[0]
    mose_x_slice = nm[0][:Iq * Jq]
    q2d_tauX = cq[1]  # cc var index 1 = tauX (after y+)
    mose_tauX = cm[1]  # cc var index 1 = tauX (after y+, before tauY, tauZ)

    if len(q2d_tauX) != len(mose_tauX):
        return False, "tauX length mismatch Q2D {0} vs MOSE {1}".format(
            len(q2d_tauX), len(mose_tauX))

    # Skip leading-edge + trailing-edge (boundary effects)
    n_skip = 5  # boundary cells excluded (renamed: do not shadow lib skip())
    n_cells = len(q2d_tauX)
    cell_x = _cell_centers_x(q2d_x_nodal, Iq, Jq, 1)

    max_rel_diff = 0.0
    max_abs_diff = 0.0
    n_compared = 0
    n_diverge = 0
    tol = 0.001  # 0.1% tight tolerance — Q2D + MOSE both SA model + same 2D plane case
                  # should produce essentially bit-tight match per user-direction Sess 13
                  # "stesso modello → identici, molto più basse". Empirically Q2D achieves
                  # max_rel_diff ≈ 9e-04 (0.09%) vs MOSE-SA.tec reference; ±0.1% is
                  # comfortable headroom while signalling Q2D ≠ MOSE divergence if it
                  # ever arises. Reference: docs/vv/regression_tolerance.md §2.3.
    for k in range(n_skip, n_cells - n_skip):
        tq = q2d_tauX[k]
        tm = mose_tauX[k]
        if abs(tm) < 1e-30:
            continue
        rel = abs(tq - tm) / abs(tm)
        max_rel_diff = max(max_rel_diff, rel)
        max_abs_diff = max(max_abs_diff, abs(tq - tm))
        n_compared += 1
        if rel > tol:
            n_diverge += 1

    status = "max_rel_diff={0:.3e} max_abs_diff={1:.3e} cells_compared={2} cells_above_tol={3}".format(
        max_rel_diff, max_abs_diff, n_compared, n_diverge)
    return (max_rel_diff < tol), status


def _sign_crossings(cx, tau, lo, hi):
    """Interpolated x where tau changes sign within (lo, hi). Returns sorted list."""
    zc = []
    prev = None
    px = None
    for k in range(len(tau)):
        if lo < cx[k] < hi:
            if prev is not None and (prev < 0.0) != (tau[k] < 0.0):
                # linear interpolation to the zero
                xz = px - prev * (cx[k] - px) / (tau[k] - prev)
                zc.append(xz)
            prev = tau[k]
            px = cx[k]
    return zc


def _ref_zero_crossings(path):
    """Zero-Cf crossings of an (x, Cf) reference dataset."""
    pts = []
    try:
        with open(path) as fh:
            for ln in fh:
                s = ln.split()
                if len(s) >= 2:
                    try:
                        pts.append((float(s[0]), float(s[1])))
                    except ValueError:
                        continue
    except Exception:
        return []
    zc = []
    for k in range(1, len(pts)):
        if (pts[k - 1][1] < 0.0) != (pts[k][1] < 0.0):
            x0, y0 = pts[k - 1]
            x1, y1 = pts[k]
            zc.append(x0 - y0 * (x1 - x0) / (y1 - y0))
    return zc


def separation_check():
    """Independent-reference V&V on the separation bubble (F13, Sess 39).

    The shock/boundary-layer interaction separates the wall boundary layer; the
    separation x_sep and reattachment x_reatt (tauX sign changes) are a
    scale-invariant physical V&V metric. Bracket the solver's bubble against the
    envelope of the independent references (Schulein exp + SU2 + Wind).

    Returns (pass_bool, status_msg). The caller has already verified that
    wall.tec exists and the run is converged; any problem here (including
    unparseable committed references) is a hard failure.
    """
    q2d_wall = "OUTPUT/wall.tec"

    try:
        Iq, Jq, Kq, nq, cq = _parse_tecplot_block(q2d_wall,
                                                  nodal_var_count=2,
                                                  cc_var_count=6)
    except Exception as e:
        return False, "Q2D wall.tec parse error: {0}".format(e)

    x_nodal = nq[0]
    tauX = cq[1]
    cell_x = _cell_centers_x(x_nodal, Iq, Jq, 1)

    # Interaction zone (well inside the plate; brackets all references' crossings).
    lo, hi = 0.25, 0.45
    q_cross = _sign_crossings(cell_x, tauX, lo, hi)
    if len(q_cross) < 2:
        return False, ("solver shows no separation bubble in x[{0},{1}] "
                       "(found {2} sign change(s); SA must separate here)").format(
            lo, hi, len(q_cross))
    q_sep, q_reatt = q_cross[0], q_cross[-1]

    # Build the independent-reference envelope (no invented numbers).
    seps, reatts = [], []
    for ref in ("schulein.dat", "SU2-SA.dat", "wind-SA.dat"):
        zc = _ref_zero_crossings(os.path.join("reference", ref))
        if len(zc) >= 2:
            seps.append(zc[0])
            reatts.append(zc[-1])
    if len(seps) < 2:
        # The references are committed files: failing to parse >= 2 of them is
        # a repository/parse regression, not an optional artifact -> FAIL.
        return False, ("insufficient independent references parsed "
                       "({0}/3 usable; committed reference/*.dat corrupt or "
                       "missing)").format(len(seps))

    # Engineering V&V band around the independent envelope. The references
    # themselves disagree by ~1.1% on x_sep (Schulein 0.3337 vs SU2 0.3375) —
    # that inter-reference spread is the irreducible physical uncertainty of
    # "where separation is", so a band tighter than it would be unphysical. SA
    # is additionally known to slightly delay/shorten SWBLI separation, and Q2D
    # (== MOSE, same SA) sits 0.44% from the nearest CFD-SA peer (SU2) and 1.59%
    # from the experiment. A 1.5% relative band (covers the inter-reference
    # spread + margin; standard for RANS separation-location V&V) is the
    # defensible tolerance — NOT tuned to pass. (At 1 cell ~0.3% the band would
    # be tighter than the references disagree with each other.)
    TOL = 0.015
    sep_lo, sep_hi = min(seps) * (1.0 - TOL), max(seps) * (1.0 + TOL)
    reatt_lo, reatt_hi = min(reatts) * (1.0 - TOL), max(reatts) * (1.0 + TOL)

    ok_sep = sep_lo <= q_sep <= sep_hi
    ok_reatt = reatt_lo <= q_reatt <= reatt_hi

    status = ("x_sep={0:.4f} in [{1:.4f},{2:.4f}]={3} ; "
              "x_reatt={4:.4f} in [{5:.4f},{6:.4f}]={7} ; "
              "refs sep{8} reatt{9}").format(
        q_sep, sep_lo, sep_hi, ok_sep,
        q_reatt, reatt_lo, reatt_hi, ok_reatt,
        ["{0:.4f}".format(s) for s in seps],
        ["{0:.4f}".format(r) for r in reatts])
    return (ok_sep and ok_reatt), status


def main():
    smoke_check()  # fail()s (exit 1) on any smoke violation

    # The quantitative gates need the wall output and the committed MOSE
    # reference. Q2D.sh wipes OUTPUT/ before each solve, so a missing wall.tec
    # means the fresh run did not produce it -> hard FAIL (never smoke-pass).
    q2d_wall = "OUTPUT/wall.tec"
    require_file(q2d_wall, "solver did not write the wall output "
                 "(required for the Cf + separation V&V gates)")
    require_file("reference/MOSE-SA.tec",
                 "committed MOSE consistency reference missing from the repo")

    # Under-converged run: the engineering V&V gates are meaningless below
    # SOLUTIONTIME 95k (target 100k per input.ini) -> SKIP (exit 77, reported
    # "Skipped" by CTest, never "Passed").
    q2d_st = _wall_tec_solutiontime(q2d_wall)
    if q2d_st is not None and q2d_st < 95000:
        skip("Q2D SOLUTIONTIME={0:.0f} < 95000: insufficient convergence for "
             "the Cf/separation engineering V&V gates".format(q2d_st))

    failed = False

    # (1) Cf consistency vs MOSE-SA.tec (Sess 12.A) — MOSE-shared SA, bit-tight.
    cf_pass, cf_status = cf_check()
    print("[swbli check] Cf vs MOSE-SA.tec: {0}".format(cf_status))
    if cf_pass:
        print("[swbli check] Cf consistency PASS (tolerance 0.1%)")
    else:
        print("[swbli check] Cf consistency FAIL (tolerance 0.1% exceeded)")
        failed = True

    # (2) Separation-bubble independent-reference V&V (F13, Sess 39) —
    #     Schulein exp + SU2 + Wind envelope.
    sep_pass, sep_status = separation_check()
    print("[swbli check] separation vs independent refs: {0}".format(sep_status))
    if sep_pass:
        print("[swbli check] separation-bubble V&V PASS (within exp/CFD envelope)")
    else:
        print("[swbli check] separation-bubble V&V FAIL (outside exp/CFD envelope)")
        failed = True

    if failed:
        return 1
    print("[swbli check] ALL GATES PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
