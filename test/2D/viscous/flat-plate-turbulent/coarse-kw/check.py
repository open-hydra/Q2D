"""flat-plate-turbulent/coarse-kw 2D viscous regression check (NASA TMR reference enforcement).

History:
  - smoke infrastructure (positivity + NaN check).
  - cf-sst.dat NASA TMR reference enforcement (matching turbulence = SST
    model per input.ini); cf-wilcox.dat preserved as informative
    cross-check (Wilcox2006 model — different from SST run).
  - hardening: removed the two silent "SMOKE PASS only" pass-through
    paths (missing wall.tec / missing cf-sst.dat now FAIL); added a PRIMARY
    de-normalized physical-Cf gate (the old shape gate normalizes both curves
    by their own maxima, which cancels any constant multiplicative error in
    tauX — e.g. a wrong mu or a factor-2 in the wall-shear extraction would
    still "PASS" the shape gate).

Acceptance criteria:
  - Smoke: OUTPUT/field*.tec exists, no NaN/Inf, n_tokens > 0. Missing ->
    FAIL (exit 1) — Q2D.sh wipes OUTPUT/ before each solve, so a missing
    output is always a fresh-run defect, never a stale leftover.
  - PRIMARY physical-Cf gate: Cf(x) = tauX / (0.5 rho_inf U_inf^2) vs the
    NASA TMR SST curve (cf-sst.dat, Tecplot zone "FUN3D, triangles") in
    PHYSICAL units; max rel diff
    <= 15% over x in [0.05, 1.0]. The freestream state is derived
    unambiguously from the committed inputs (see FREESTREAM DERIVATION in
    the code): T_inf = 555 K, p_inf = 1e5 Pa, M_inf = 0.2, gamma = 1.4,
    R = 287 -> rho_inf = 0.62781 kg/m^3, U_inf = 94.446 m/s. Cross-checks:
    rho_inf*U_inf = 59.2937 kg/(m^2 s) == the BC-403 mass flux 59.2934 in
    INPUT/bc.txt; Re/m = rho_inf U_inf / mil = 5.00e6 = the NASA TMR
    flat-plate unit Reynolds number. Measured on the committed converged run
    (SOLUTIONTIME=1e5): max rel diff 4.7% — the 15% tolerance is coarse-grid
    engineering headroom, not tuning.
  - SECONDARY Cf shape gate vs cf-sst.dat: both curves normalized by their
    plate-region maxima, max rel shape diff <= 10% over x in [0.05, 1.0]
    (scale-free curve-shape tripwire).
  - cf-wilcox.dat cross-check INFORMATIVE only (different model — Wilcox2006
    vs SST).

Exit semantics: exit 0 only when smoke + BOTH Cf gates pass; exit 1 on any
missing output, missing committed reference, NaN/Inf, parse error, or gate
violation. This script never exits 77 (no legitimate skip condition).

References:
  - cf-sst.dat: copied from MOSE/test/2D/viscous/flat-plate-turbulent/reference/
    (NASA TMR for the SST model; the file's own Tecplot zone header reads
    T="FUN3D, triangles" — FUN3D, not CFL3D); matches Q2D turbulence = SST
    per input.ini.
  - cf-wilcox.dat: NASA TMR for the Wilcox2006 model (two Tecplot zones:
    T="CFL3D" and T="FUN3D"); informative cross-check only — Q2D ran SST,
    a different model, so a larger relative diff than cf-sst.dat is
    expected.
"""
import os
import sys
import glob
import math
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "common" / "pycheck"))
from q2d_check_lib import fail, require_file, assert_finite

# ---------------------------------------------------------------------------
# FREESTREAM DERIVATION (from committed inputs only; no invented numbers):
#   input.ini [GPB-Phase1]  : gamma = 1.4, R = 287.0, mil = 1.1858685985e-5
#   input.ini [ICB-Block1]  : T = 555 K, mach = 0.2
#   INPUT/bc.txt BC 407 exit: p_back = 1.0e5 Pa (subsonic plate -> p_inf)
#   INPUT/bc.txt BC 403 inlet: T0 = 555 K, mass flux g = 59.2934 kg/(m^2 s)
# Consistency: rho_inf = p_inf/(R T_inf) = 0.627806 kg/m^3;
#   U_inf = M sqrt(gamma R T_inf) = 94.4456 m/s;
#   rho_inf*U_inf = 59.2937 == bc.txt g = 59.2934;
#   Re/m = rho_inf U_inf / mil = 5.00e6 == NASA TMR unit Reynolds number.
# ---------------------------------------------------------------------------
GAMMA = 1.4
R_GAS = 287.0
T_INF = 555.0
P_INF = 1.0e5
MACH_INF = 0.2
RHO_INF = P_INF / (R_GAS * T_INF)
U_INF = MACH_INF * math.sqrt(GAMMA * R_GAS * T_INF)
Q_INF = 0.5 * RHO_INF * U_INF ** 2  # dynamic pressure ~ 2799.7 Pa


def smoke_check():
    """Field-output smoke gate: missing output or NaN/Inf -> fail() (exit 1)."""
    out_dir = "OUTPUT"
    if not os.path.isdir(out_dir):
        fail("OUTPUT/ missing (solver did not run)")

    candidates = sorted(glob.glob(os.path.join(out_dir, "field*.tec")),
                        key=os.path.getmtime)
    if not candidates:
        fail("no field*.tec in OUTPUT/ (solver produced no field output)")

    field_file = candidates[-1]
    require_file(field_file)
    print("[FP-turb-kw check] inspecting {0}".format(field_file))

    values = []
    with open(field_file, "r") as fh:
        for line in fh:
            tok = line.strip()
            if not tok or tok.startswith("VARIABLES") or tok.startswith("ZONE"):
                continue
            for x in tok.split():
                try:
                    values.append(float(x))
                except ValueError:
                    continue

    print("[FP-turb-kw check] smoke: n_tot={0}".format(len(values)))
    if not values:
        fail("no numeric data parsed from {0} (truncated output?)".format(field_file))
    assert_finite("{0} values".format(field_file), values)


def _parse_tecplot_wall(path):
    """Parse Q2D wall.tec NODAL(x,y) + CC(y+,tauX,tauY,pw,Tw,qw) BLOCK format.

    Returns: (I, x_nodal_list, tauX_cc_list).
    Q2D wall.tec is 2D with K=1, J=1 (single wall row).
    """
    with open(path) as fh:
        lines = fh.readlines()
    I = J = K = None
    nvar = nodal_var_count = cc_var_count = 0
    data_start = 0
    for i, ln in enumerate(lines):
        u = ln.upper()
        if 'VARIABLES' in u:
            nvar = len(re.findall(r'"[^"]*"', ln))
        if 'ZONE' in u:
            m = re.search(r'I\s*=\s*(\d+)', ln, re.IGNORECASE)
            if m:
                I = int(m.group(1))
            m = re.search(r'J\s*=\s*(\d+)', ln, re.IGNORECASE)
            if m:
                J = int(m.group(1))
            m = re.search(r'K\s*=\s*(\d+)', ln, re.IGNORECASE)
            if m:
                K = int(m.group(1))
            m = re.search(r'VARLOCATION\s*=\s*\(\s*\[\s*1\s*-\s*(\d+)\s*\]\s*=\s*NODAL',
                          ln, re.IGNORECASE)
            if m:
                nodal_var_count = int(m.group(1))
                cc_var_count = nvar - nodal_var_count
            data_start = i + 1
            break
    if I is None or J is None or K is None or nodal_var_count == 0:
        raise RuntimeError("wall.tec ZONE parse failed: {0}".format(path))

    n_nodes = I * J * K
    if K == 1 and J == 1:
        n_cells = max(1, I - 1)
    else:
        n_cells = max(1, I - 1) * max(1, J - 1) * max(1, K - 1)

    floats = []
    for ln in lines[data_start:]:
        s = ln.replace(',', ' ').strip()
        if not s:
            continue
        for t in s.split():
            try:
                floats.append(float(t.replace('D', 'E').replace('d', 'e')))
            except ValueError:
                continue

    idx = 0
    nodal_arrays = []
    for _ in range(nodal_var_count):
        nodal_arrays.append(floats[idx:idx + n_nodes])
        idx += n_nodes
    cc_arrays = []
    for _ in range(cc_var_count):
        cc_arrays.append(floats[idx:idx + n_cells])
        idx += n_cells

    x_nodal = nodal_arrays[0]
    tauX_cc = cc_arrays[1] if len(cc_arrays) > 1 else []
    return I, x_nodal, tauX_cc


def _read_cf_reference(path):
    """Read cf-*.dat reference (VARIABLES = x, cf; ZONE T=...). Returns
    (x_list, cf_list) skipping header lines."""
    x_list = []
    cf_list = []
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.upper().startswith('VARIABLES') or s.upper().startswith('ZONE'):
                continue
            tok = s.split()
            if len(tok) < 2:
                continue
            try:
                xv = float(tok[0])
                cv = float(tok[1])
            except ValueError:
                continue
            x_list.append(xv)
            cf_list.append(cv)
    return x_list, cf_list


def _interp_linear(x_grid, y_grid, x_query):
    if not x_grid:
        return None
    if x_query <= x_grid[0]:
        return y_grid[0]
    if x_query >= x_grid[-1]:
        return y_grid[-1]
    lo, hi = 0, len(x_grid) - 1
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if x_grid[mid] <= x_query:
            lo = mid
        else:
            hi = mid
    dx = x_grid[hi] - x_grid[lo]
    if dx == 0:
        return y_grid[lo]
    return y_grid[lo] + (y_grid[hi] - y_grid[lo]) * (x_query - x_grid[lo]) / dx


def cf_shape_compare(q2d_x_cells, q2d_tauX, ref_x, ref_cf, tol, label,
                      x_lo=0.05, x_hi=1.0):
    """Compare normalized cf-like shape: tauX/max(plate) vs cf/max(plate).

    Plate region [x_lo, x_hi] (skip leading-edge BC effects + trailing-edge).
    Returns (pass_bool, status_msg).
    """
    q2d_xc = []
    for i in range(len(q2d_x_cells) - 1):
        q2d_xc.append(0.5 * (q2d_x_cells[i] + q2d_x_cells[i + 1]))
    if len(q2d_xc) != len(q2d_tauX):
        return False, "{0}: x_cc len {1} != tauX len {2}".format(
            label, len(q2d_xc), len(q2d_tauX))

    q2d_plate = [t for x, t in zip(q2d_xc, q2d_tauX)
                  if x_lo <= x <= x_hi and t > 0]
    if not q2d_plate:
        return False, "{0}: no positive Q2D tauX in plate region [{1},{2}]".format(
            label, x_lo, x_hi)
    q2d_max = max(q2d_plate)

    ref_plate = [(x, c) for x, c in zip(ref_x, ref_cf)
                  if x_lo <= x <= x_hi and c > 0]
    if not ref_plate:
        return False, "{0}: no positive ref cf in plate region [{1},{2}]".format(
            label, x_lo, x_hi)
    ref_max = max(c for _, c in ref_plate)
    ref_x_plate = [x for x, _ in ref_plate]
    ref_cf_plate = [c for _, c in ref_plate]

    max_rel_diff = 0.0
    n_compared = 0
    n_diverge = 0
    for x, t in zip(q2d_xc, q2d_tauX):
        if not (x_lo <= x <= x_hi) or t <= 0:
            continue
        q2d_norm = t / q2d_max
        ref_at_x = _interp_linear(ref_x_plate, ref_cf_plate, x)
        if ref_at_x is None or ref_at_x <= 0:
            continue
        ref_norm = ref_at_x / ref_max
        rel = abs(q2d_norm - ref_norm) / max(abs(ref_norm), 1e-30)
        max_rel_diff = max(max_rel_diff, rel)
        n_compared += 1
        if rel > tol:
            n_diverge += 1

    if n_compared == 0:
        return False, "{0}: no overlap Q2D x_cc vs ref x in plate region".format(label)
    status = "{0}: max_rel_diff={1:.3e} cells_compared={2} cells_above_tol={3}".format(
        label, max_rel_diff, n_compared, n_diverge)
    return (max_rel_diff < tol), status


def cf_physical_compare(q2d_x_nodal, q2d_tauX, ref_x, ref_cf, tol, label,
                        x_lo=0.05, x_hi=1.0):
    """PRIMARY de-normalized gate: physical Cf = tauX / (0.5 rho_inf U_inf^2)
    vs the NASA TMR reference cf(x) in physical units (no self-normalization,
    so a constant multiplicative error in tauX — wrong mu, factor-2 in the
    wall-shear extraction, wrong reference state — trips this gate even though
    it cancels in the shape gate). Max rel diff <= tol over [x_lo, x_hi].

    Returns (pass_bool, status_msg).
    """
    q2d_xc = []
    for i in range(len(q2d_x_nodal) - 1):
        q2d_xc.append(0.5 * (q2d_x_nodal[i] + q2d_x_nodal[i + 1]))
    if len(q2d_xc) != len(q2d_tauX):
        return False, "{0}: x_cc len {1} != tauX len {2}".format(
            label, len(q2d_xc), len(q2d_tauX))

    ref_plate = [(x, c) for x, c in zip(ref_x, ref_cf)
                 if x_lo <= x <= x_hi and c > 0]
    if not ref_plate:
        return False, "{0}: no positive ref cf in plate region [{1},{2}]".format(
            label, x_lo, x_hi)
    ref_x_plate = [x for x, _ in ref_plate]
    ref_cf_plate = [c for _, c in ref_plate]

    max_rel_diff = 0.0
    n_compared = 0
    n_diverge = 0
    for x, t in zip(q2d_xc, q2d_tauX):
        if not (x_lo <= x <= x_hi):
            continue
        cf_q2d = t / Q_INF
        ref_at_x = _interp_linear(ref_x_plate, ref_cf_plate, x)
        if ref_at_x is None or ref_at_x <= 0:
            continue
        rel = abs(cf_q2d - ref_at_x) / abs(ref_at_x)
        max_rel_diff = max(max_rel_diff, rel)
        n_compared += 1
        if rel > tol:
            n_diverge += 1

    if n_compared == 0:
        return False, "{0}: no overlap Q2D x_cc vs ref x in plate region".format(label)
    status = ("{0}: max_rel_diff={1:.3e} cells_compared={2} cells_above_tol={3} "
              "(q_inf={4:.4g} Pa, rho_inf={5:.5g}, U_inf={6:.5g})").format(
        label, max_rel_diff, n_compared, n_diverge, Q_INF, RHO_INF, U_INF)
    return (max_rel_diff < tol), status


def main():
    smoke_check()  # fail()s (exit 1) on any smoke violation

    wall_path = "OUTPUT/wall.tec"
    ref_sst_path = "reference/cf-sst.dat"
    ref_wilcox_path = "reference/cf-wilcox.dat"

    # Q2D.sh wipes OUTPUT/ before each solve: missing wall.tec means the fresh
    # run did not write it -> hard FAIL (was: silent "SMOKE PASS only").
    require_file(wall_path, "solver did not write the wall output "
                 "(required for the Cf gates)")
    # cf-sst.dat is a committed reference: missing means repo damage -> FAIL
    # (was: silent "SMOKE PASS only").
    require_file(ref_sst_path,
                 "committed NASA TMR SST reference missing from the repo")

    try:
        I, x_nodal, tauX_cc = _parse_tecplot_wall(wall_path)
    except Exception as e:
        fail("wall.tec parse error: {0}".format(e))
    print("[FP-turb-kw check] wall.tec I={0} cells={1}".format(I, len(tauX_cc)))
    assert_finite("wall.tec tauX", tauX_cc)
    assert_finite("wall.tec x", x_nodal)

    sst_x, sst_cf = _read_cf_reference(ref_sst_path)

    # PRIMARY: de-normalized physical Cf vs NASA TMR SST (15% tolerance:
    # coarse 68x48 grid engineering V&V; catches multiplicative errors).
    pass_phys, status_phys = cf_physical_compare(
        x_nodal, tauX_cc, sst_x, sst_cf, tol=0.15,
        label="physical Cf vs cf-sst.dat (NASA TMR SST)")
    print("[FP-turb-kw check] {0}".format(status_phys))

    # SECONDARY: scale-free shape gate, kept as a curve-shape tripwire.
    pass_sst, status_sst = cf_shape_compare(
        x_nodal, tauX_cc, sst_x, sst_cf, tol=0.10,
        label="Cf shape vs cf-sst.dat (NASA TMR SST)")
    print("[FP-turb-kw check] {0}".format(status_sst))

    # Wilcox2006 informative cross-check.
    if os.path.isfile(ref_wilcox_path):
        wx_x, wx_cf = _read_cf_reference(ref_wilcox_path)
        _, status_wx = cf_shape_compare(
            x_nodal, tauX_cc, wx_x, wx_cf, tol=0.10,
            label="Cf shape vs cf-wilcox.dat (NASA TMR Wilcox2006)")
        print("[FP-turb-kw check] {0} (INFORMATIVE cross-check; Q2D ran SST)".format(status_wx))

    failed = False
    if not pass_phys:
        print("[FP-turb-kw check] FAIL: physical Cf vs cf-sst.dat +/-15% acceptance NOT met")
        failed = True
    if not pass_sst:
        print("[FP-turb-kw check] FAIL: Cf shape vs cf-sst.dat +/-10% acceptance NOT met")
        failed = True
    if failed:
        return 1
    print("[FP-turb-kw check] PASS physical Cf +/-15% + shape +/-10% + smoke acceptance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
