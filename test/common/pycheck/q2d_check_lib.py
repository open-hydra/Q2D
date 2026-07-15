"""Shared verification helpers for Q2D ctest check scripts.

Centralizes the utilities that were previously copy-pasted across the
per-test ``check.py`` scripts:

- ORION Tecplot reader path setup + the 2D->3D tec patch helper;
- NaN-safe positivity / finiteness guards (``x.min() <= 0`` is False for
  NaN fields: every positivity gate MUST go through ``assert_positive``);
- relative L1/L2/Linf error norms;
- cell-by-cell field comparison against a committed reference field;
- uniform exit conventions: ``fail()`` -> exit 1, ``skip()`` -> exit 77
  (CTest SKIP_RETURN_CODE), so a skipped gate is visibly "Skipped",
  never silently "Passed".

Import boilerplate for a test living at ``test/<family>/<case>``::

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[<n>] / "common" / "pycheck"))
    from q2d_check_lib import *
"""
import math
import os
import re
import sys
import tempfile
from pathlib import Path

import numpy as np

SKIP_RC = 77  # keep in sync with SKIP_RETURN_CODE in test/CMakeLists.txt

_here = Path(__file__).resolve()
REPO_ROOT = _here.parents[3]  # <repo>/test/common/pycheck/q2d_check_lib.py


def _setup_orion_path():
    candidate = REPO_ROOT / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
    else:
        fail(f"ORION python package not found at {candidate}")


def fail(msg):
    """Print a FAIL message and exit 1 (test failure)."""
    print(f"FAIL: {msg}")
    sys.exit(1)


def skip(msg):
    """Print a SKIP message and exit 77 (CTest 'Skipped', never 'Passed')."""
    print(f"SKIP: {msg}")
    sys.exit(SKIP_RC)


def require_file(path, hint=""):
    """Hard-fail if an expected solver output is missing.

    Q2D.sh wipes OUTPUT/ before each solve, so a missing file always means
    the fresh run did not produce it (crash, writer regression, or
    under-convergence) -- never a stale leftover.
    """
    p = Path(path)
    if not p.is_file():
        fail(f"expected output '{path}' not found{': ' + hint if hint else ''}")
    if p.stat().st_size == 0:
        fail(f"expected output '{path}' is empty{': ' + hint if hint else ''}")
    return p


def assert_finite(name, arr):
    arr = np.asarray(arr, dtype=float)
    if not np.isfinite(arr).all():
        n_bad = int((~np.isfinite(arr)).sum())
        fail(f"{name} contains {n_bad} non-finite values (NaN/Inf)")


def assert_positive(name, arr, strict=True):
    """NaN-safe positivity gate (finiteness checked first)."""
    assert_finite(name, arr)
    arr = np.asarray(arr, dtype=float)
    lo = arr.min()
    if strict and lo <= 0.0:
        fail(f"min({name}) = {lo:g} <= 0")
    if not strict and lo < 0.0:
        fail(f"min({name}) = {lo:g} < 0")


def rel_l1(computed, exact, floor=1e-300):
    c, e = np.asarray(computed, float), np.asarray(exact, float)
    return float(np.abs(c - e).sum() / max(np.abs(e).sum(), floor))


def rel_l2(computed, exact, floor=1e-300):
    c, e = np.asarray(computed, float), np.asarray(exact, float)
    return float(np.sqrt(((c - e) ** 2).sum()) / max(np.sqrt((e ** 2).sum()), floor))


def rel_linf(computed, exact, floor=1e-300):
    c, e = np.asarray(computed, float), np.asarray(exact, float)
    return float(np.abs(c - e).max() / max(np.abs(e).max(), floor))


def gate(label, value, tol):
    """Print a PASS/FAIL line for a scalar gate; exit 1 on violation."""
    status = "PASS" if value <= tol else "FAIL"
    print(f"  [{status}] {label}: {value:.6g} (tol {tol:g})")
    if value > tol:
        fail(f"{label} = {value:.6g} exceeds tolerance {tol:g}")
    return value


def patch_q2d_tec_to_3d(in_path):
    """Rewrite a Q2D 2D tecplot ASCII file as 3D (adds a zero 'z' block per
    zone) so the ORION reader can parse it. Handles single- and multi-zone
    files. Returns the path of a temporary patched copy."""
    require_file(in_path)
    with open(in_path) as f:
        lines = f.readlines()
    out = [lines[0].replace('"y"', '"y" "z"', 1)]
    i = 1
    n_zones = 0
    while i < len(lines):
        line = lines[i]
        if 'ZONE' in line.upper() and 'I=' in line.upper():
            m = re.search(r'I=\s*(\d+).*J=\s*(\d+).*K=\s*(\d+)', line)
            if m is None:
                fail(f"'{in_path}' zone record has no I=/J=/K= (truncated output?)")
            ni, nj, nk = int(m.group(1)), int(m.group(2)), int(m.group(3))
            nnodes = ni * nj * nk
            ncells = max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1)
            m2 = re.search(r'VARLOCATION=\(\[(\d+)-(\d+)\]=NODAL,\[(\d+)-(\d+)\]=CELLCENTERED', line)
            if m2 is None:
                fail(f"'{in_path}' zone record has no VARLOCATION (unsupported layout)")
            na, nb_, ca, cb = (int(m2.group(g)) for g in (1, 2, 3, 4))
            out.append(line.replace(
                f'[{na}-{nb_}]=NODAL,[{ca}-{cb}]=CELLCENTERED',
                f'[{na}-{nb_ + 1}]=NODAL,[{ca + 1}-{cb + 1}]=CELLCENTERED',
            ))
            i += 1
            n_nodal = nb_ - na + 1
            n_cc = cb - ca + 1
            out.extend(lines[i: i + n_nodal * nnodes])   # nodal vars (x, y)
            i += n_nodal * nnodes
            out.append(' 0.000000000000000E+000\n' * nnodes)  # fake z
            out.extend(lines[i: i + n_cc * ncells])      # cell-centered vars
            i += n_cc * ncells
            n_zones += 1
        else:
            out.append(line)
            i += 1
    if n_zones == 0:
        fail(f"'{in_path}' has no ZONE record (truncated output?)")
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.tec', delete=False)
    tmp.writelines(out)
    tmp.close()
    return tmp.name


def read_q2d_field(path):
    """Read a Q2D tecplot field file. Returns (x, y, var, varnames) where
    ``var`` follows the ORION read_TEC convention (var[block][ivar] arrays,
    cell-centered variables only) and ``varnames`` lists ALL variables as
    stored (nodal coordinates first: varnames[3:] map onto var[block][:]).

    2D Q2D outputs (no "z" variable) are patched to 3D on the fly; files
    that already carry a "z" coordinate (e.g. committed 3D-wedge reference
    fields) are read as-is."""
    _setup_orion_path()
    from ORION import read_TEC
    require_file(path)
    if 'z' in [n.strip().lower() for n in field_variable_names(path)[:3]]:
        x, y, _z, var, vnames = read_TEC(str(path))
        return x, y, var, vnames
    patched = patch_q2d_tec_to_3d(path)
    try:
        x, y, _z, var, vnames = read_TEC(patched)
    finally:
        os.unlink(patched)
    return x, y, var, vnames


def field_variable_names(path):
    """Variable names from a Q2D tec header (order as stored)."""
    with open(path) as f:
        header = f.readline()
    return re.findall(r'"([^"]+)"', header)


def compare_fields(out_path, ref_path, tol_linf, tol_l2, skip_vars=("x", "y", "z"),
                   floor=1e-30):
    """Cell-by-cell comparison of a solver field against a committed reference.

    For every common cell-centered variable, computes relative Linf and L2
    (normalized by the reference's max-abs, floored) and gates both.
    """
    xo, yo, vo, no = read_q2d_field(out_path)
    xr, yr, vr, nr = read_q2d_field(ref_path)
    # ORION returns ALL names but only the cell-centered arrays: the first
    # 3 names are the nodal coordinates, so var[block][i] <-> names[3 + i].
    cc_o = [n.strip() for n in no[3:]]
    cc_r = [n.strip() for n in nr[3:]]
    common = [n for n in cc_r if n in cc_o and n not in skip_vars]
    if not common:
        fail(f"no common variables between '{out_path}' and '{ref_path}'")
    if len(vo) != len(vr):
        fail(f"block count mismatch: {len(vo)} in '{out_path}' vs {len(vr)} in '{ref_path}'")
    worst = (None, 0.0)
    print(f"  field comparison {out_path} vs {ref_path} "
          f"(tol Linf {tol_linf:g}, L2 {tol_l2:g}):")
    for name in common:
        a = np.concatenate([np.asarray(vb[cc_o.index(name)], float).ravel() for vb in vo])
        b = np.concatenate([np.asarray(vb[cc_r.index(name)], float).ravel() for vb in vr])
        if a.shape != b.shape:
            fail(f"variable '{name}' shape mismatch {a.shape} vs {b.shape}")
        assert_finite(f"output '{name}'", a)
        scale = max(np.abs(b).max(), floor)
        linf = float(np.abs(a - b).max() / scale)
        l2 = float(np.sqrt(((a - b) ** 2).mean()) / scale)
        ok = linf <= tol_linf and l2 <= tol_l2
        print(f"    [{'PASS' if ok else 'FAIL'}] {name:16s} Linf_rel={linf:.3e}  L2_rel={l2:.3e}")
        if linf > worst[1]:
            worst = (name, linf)
        if not ok:
            fail(f"variable '{name}' exceeds tolerance (Linf_rel={linf:.3e}, L2_rel={l2:.3e})")
    print(f"  worst variable: {worst[0]} (Linf_rel={worst[1]:.3e})")


def read_residual_history(path="OUTPUT/residual-history.dat"):
    """Return the residual history as a 2D float array (iters x columns).

    Q2D layout: column 0 = iteration, column 1 = time (constant -1.0 for
    steady runs), column 2 = density residual, columns 3.. = the remaining
    equation residuals.
    """
    require_file(path, "solver did not write the residual history")
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def assert_converged(threshold, column=2, path="OUTPUT/residual-history.dat"):
    """Hard convergence sentinel: final residual (given column) < threshold.

    ``column`` defaults to 2 = the density residual (column 1 is the time
    stamp, constant -1.0 for steady runs -- gating on it would always pass).
    """
    data = read_residual_history(path)
    final = float(data[-1, column])
    assert_finite("residual history", data)
    status = "PASS" if final < threshold else "FAIL"
    print(f"  [{status}] convergence: final residual {final:.3e} (threshold {threshold:g})")
    if not final < threshold:
        fail(f"final residual {final:.3e} did not reach threshold {threshold:g}")
    return final
