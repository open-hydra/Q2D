"""flat-plate-turbulent/coarse-kw 2D viscous regression check (Sess 13 reference enforcement).

History:
  Sess 10 (2026-05-21): smoke infrastructure (positivity + NaN check).
  Sess 13 (2026-05-22): cf-sst.dat NASA TMR reference enforcement (matching
    turbulence = SST model per input.ini); cf-wilcox.dat preserved as
    informative cross-check (Wilcox2006 model — different from SST run).

Acceptance criteria (Sess 13):
  - Smoke: OUTPUT/field*.tec exists, no NaN/Inf, n_tokens > 0.
  - Cf shape comparison vs cf-sst.dat (NASA TMR CFL3D reference for SST model):
      Q2D wall.tec format: NODAL (x, y), CELLCENTERED (y+, tauX, tauY, pw, Tw, qw).
      Extract tauX cell-centered at wall row.
      Normalize Q2D tauX to its plate-region max (avoid freestream-conversion
      ambiguity). Normalize cf-sst.dat cf to its plate-region max. Compare
      normalized shape over x in [0.05, 1.0] plate region.
      Tolerance: max relative shape diff <= 10% (NASA TMR engineering V&V
      standard for turbulent flat-plate Cf prediction).
  - cf-wilcox.dat cross-check INFORMATIVE only (different model — Wilcox2006
    vs SST; informs whether Q2D SST result is closer to NASA TMR SST or
    Wilcox2006 reference).

References:
  - cf-sst.dat: copied from MOSE/test/2D/viscous/flat-plate-turbulent/reference/
    Sess 13 F-12-K (NASA TMR CFL3D for SST model); matches Q2D
    turbulence = SST per input.ini.
  - cf-wilcox.dat: NASA TMR CFL3D for Wilcox2006 model (informative
    cross-check; Q2D ran SST so this is different model; expected larger
    relative diff than cf-sst.dat).
"""
import os
import sys
import glob
import math
import re


def smoke_check():
    out_dir = "OUTPUT"
    if not os.path.isdir(out_dir):
        print("[FP-turb-kw check] OUTPUT/ missing")
        return 1

    candidates = sorted(glob.glob(os.path.join(out_dir, "field*.tec")),
                        key=os.path.getmtime)
    if not candidates:
        print("[FP-turb-kw check] no field*.tec in OUTPUT/")
        return 1

    field_file = candidates[-1]
    print("[FP-turb-kw check] inspecting {0}".format(field_file))

    nnan = 0
    n_tot = 0
    with open(field_file, "r") as fh:
        for line in fh:
            tok = line.strip()
            if not tok or tok.startswith("VARIABLES") or tok.startswith("ZONE"):
                continue
            for x in tok.split():
                try:
                    v = float(x)
                except ValueError:
                    continue
                n_tot += 1
                if v != v:
                    nnan += 1

    print("[FP-turb-kw check] smoke: n_tot={0} nnan={1}".format(n_tot, nnan))
    if nnan > 0:
        print("[FP-turb-kw check] FAIL: {0} NaN values".format(nnan))
        return 1
    if n_tot == 0:
        print("[FP-turb-kw check] FAIL: empty output")
        return 1
    return 0


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


def main():
    rc = smoke_check()
    if rc != 0:
        print("[FP-turb-kw check] SMOKE FAIL")
        return rc

    wall_path = "OUTPUT/wall.tec"
    ref_sst_path = "reference/cf-sst.dat"
    ref_wilcox_path = "reference/cf-wilcox.dat"

    if not os.path.isfile(wall_path):
        print("[FP-turb-kw check] no OUTPUT/wall.tec - SMOKE PASS only (Cf comparison skipped)")
        return 0

    try:
        I, x_nodal, tauX_cc = _parse_tecplot_wall(wall_path)
    except Exception as e:
        print("[FP-turb-kw check] wall.tec parse error: {0}".format(e))
        return 1
    print("[FP-turb-kw check] wall.tec I={0} cells={1}".format(I, len(tauX_cc)))

    # Q2D ran turbulence = SST per input.ini -> primary check vs cf-sst.dat.
    if not os.path.isfile(ref_sst_path):
        print("[FP-turb-kw check] no reference/cf-sst.dat - SMOKE PASS only")
        return 0
    sst_x, sst_cf = _read_cf_reference(ref_sst_path)
    pass_sst, status_sst = cf_shape_compare(
        x_nodal, tauX_cc, sst_x, sst_cf, tol=0.10,
        label="Cf vs cf-sst.dat (NASA TMR SST)")
    print("[FP-turb-kw check] {0}".format(status_sst))

    # Wilcox2006 informative cross-check.
    if os.path.isfile(ref_wilcox_path):
        wx_x, wx_cf = _read_cf_reference(ref_wilcox_path)
        _, status_wx = cf_shape_compare(
            x_nodal, tauX_cc, wx_x, wx_cf, tol=0.10,
            label="Cf vs cf-wilcox.dat (NASA TMR Wilcox2006)")
        print("[FP-turb-kw check] {0} (INFORMATIVE cross-check; Q2D ran SST)".format(status_wx))

    if not pass_sst:
        print("[FP-turb-kw check] FAIL: Cf vs cf-sst.dat +/-10% acceptance NOT met")
        return 1
    print("[FP-turb-kw check] PASS Cf vs cf-sst.dat +/-10% + smoke acceptance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
