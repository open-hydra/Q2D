"""supersonic-forward-step (Woodward-Colella 1984) 2D time-accurate regression.

Validates the Q2D solution at t=4.0 for the Woodward-Colella forward-step
benchmark (M=3 inflow into a channel with a 0.2-high step at x=0.6,
2-block grid: Block1 pre-step 49x17, Block2 above-step 241x65). Mirror
MOSE verify.py probe layout adapted for Q2D 2D-plane K=0 strip.

References:
  (1) Natural: OpenFOAM EnSight binary reference data shipped under
      reference/OpenFOAM/ (forwardStep.0000.{mesh,p,T,U}). REPORTED-only
      (mirror MOSE verify.py which produces a plot without hard threshold).
  (2) MOSE bit-exact: Q2D 2-zone field K=0 at t=4 vs MOSE 2-zone field K=0
      at t=4. Reference cache regression_baseline/MOSE_K0/
      supersonic-forward-step/field.tec.

Geometric classification (z-coord analysis): 2D PLANE (MOSE z=+-5e-4 const
per K-layer per zone, X/Y identical between K-layers => no axisymmetric).

Tolerances:
  Sanity:        min(rho), min(p), min(T) > 0 (hard).
  Natural OF:    REPORTED (no hard threshold; OpenFOAM is a different solver
                 with different numerics, mirror MOSE pattern).
  Bit-exact MOSE: Linf_rel < 1e-9, L2_rel < 1e-12 (atol 1e-10).
"""
import os
import re
import struct
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
warnings.filterwarnings("ignore")

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break

try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError("ORION package not found. Run 'conda activate base'.")


def patch_q2d_multizone_tec_to_3d(in_path):
    """Patch Q2D 2D-strip multi-zone field.tec (K=1) to ORION-readable 3D-wedge
    (K=2 with fake z=0). Multi-zone variant (mirror oblique-shock helper)."""
    with open(in_path) as f:
        lines = f.readlines()
    out = []
    header = lines[0]
    new_header = header.replace('"y"', '"y" "z"', 1)
    out.append(new_header)
    i = 1
    while i < len(lines):
        l = lines[i]
        if 'ZONE' in l.upper() and 'I=' in l.upper():
            zone = l
            m = re.search(r'I=\s*(\d+).*J=\s*(\d+).*K=\s*(\d+)', zone)
            if m is None:
                out.append(l); i += 1; continue
            ni, nj, nk = int(m.group(1)), int(m.group(2)), int(m.group(3))
            nnodes = ni * nj * nk
            m2 = re.search(r'VARLOCATION=\(\[(\d+)-(\d+)\]=NODAL,\[(\d+)-(\d+)\]=CELLCENTERED', zone)
            if m2:
                na, nb_, ca, cb = (int(m2.group(j)) for j in (1, 2, 3, 4))
                new_zone = zone.replace(
                    f'[{na}-{nb_}]=NODAL,[{ca}-{cb}]=CELLCENTERED',
                    f'[{na}-{nb_+1}]=NODAL,[{ca+1}-{cb+1}]=CELLCENTERED',
                )
            else:
                new_zone = zone
            out.append(new_zone)
            i += 1
            n_nodal_orig = nb_ - na + 1 if m2 else 2
            n_cc_orig    = cb - ca + 1 if m2 else 0
            for vi in range(n_nodal_orig):
                out.extend(lines[i : i + nnodes])
                i += nnodes
            out.append(' 0.000000000000000E+000\n' * nnodes)
            ncells = (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1)
            for vi in range(n_cc_orig):
                out.extend(lines[i : i + ncells])
                i += ncells
        else:
            out.append(l)
            i += 1
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.tec', delete=False)
    tmp.writelines(out); tmp.close()
    return tmp.name


# =============================================================
# Load Q2D solution
# =============================================================
patched = patch_q2d_multizone_tec_to_3d('OUTPUT/field.tec')
try:
    [x_, y_, z_, var_, vnames] = read_TEC(patched)
finally:
    os.unlink(patched)

# Sanity: all blocks; vars order rho1, u, v, p, T, g, R (CC indices 0..6)
for b in range(len(var_)):
    if var_[b][0].min() <= 0:
        print(f'Sanity FAIL b{b} min(rho)={var_[b][0].min()}'); sys.exit(1)
    if var_[b][3].min() <= 0:
        print(f'Sanity FAIL b{b} min(p)={var_[b][3].min()}'); sys.exit(1)
    if var_[b][4].min() <= 0:
        print(f'Sanity FAIL b{b} min(T)={var_[b][4].min()}'); sys.exit(1)
print('Sanity: all blocks min(rho)/min(p)/min(T) > 0 OK')

# =============================================================
# (1) Natural reference: OpenFOAM EnSight binary (REPORTED only)
# =============================================================
REF_DIR = root / 'reference' / 'OpenFOAM'


def _es_str(f):
    return f.read(80).decode('ascii').strip('\x00').strip()


def read_ensight_mesh(path):
    with open(path, 'rb') as f:
        for _ in range(5):
            _es_str(f)
        _es_str(f)
        struct.unpack('i', f.read(4))[0]
        _es_str(f)
        _es_str(f)
        nn = struct.unpack('i', f.read(4))[0]
        xn = np.frombuffer(f.read(4 * nn), dtype=np.float32)
        yn = np.frombuffer(f.read(4 * nn), dtype=np.float32)
        f.read(4 * nn)
        _es_str(f)
        ne = struct.unpack('i', f.read(4))[0]
        conn = np.frombuffer(f.read(8 * 4 * ne), dtype=np.int32).reshape(ne, 8)
    return xn, yn, conn, ne


def read_ensight_scalar(path, ne):
    with open(path, 'rb') as f:
        _es_str(f); _es_str(f)
        struct.unpack('i', f.read(4))[0]
        _es_str(f)
        return np.frombuffer(f.read(4 * ne), dtype=np.float32).astype(float)


def read_ensight_vector(path, ne):
    with open(path, 'rb') as f:
        _es_str(f); _es_str(f)
        struct.unpack('i', f.read(4))[0]
        _es_str(f)
        return np.frombuffer(f.read(3 * 4 * ne), dtype=np.float32).reshape(3, ne).astype(float)


if REF_DIR.exists() and (REF_DIR / 'forwardStep.0000.mesh').exists():
    xn_of, yn_of, conn_of, ne_of = read_ensight_mesh(REF_DIR / 'forwardStep.0000.mesh')
    of_p = read_ensight_scalar(REF_DIR / 'forwardStep.0000.p', ne_of)
    of_T = read_ensight_scalar(REF_DIR / 'forwardStep.0000.T', ne_of)
    of_U = read_ensight_vector(REF_DIR / 'forwardStep.0000.U', ne_of)
    of_xc = xn_of[conn_of - 1].mean(axis=1)
    of_yc = yn_of[conn_of - 1].mean(axis=1)
    of_umag = np.sqrt(of_U[0]**2 + of_U[1]**2 + of_U[2]**2)

    # Q2D 2D cell centres from nodal corners (z fake=0 layer; index 0)
    q2d_xc_parts, q2d_yc_parts = [], []
    q2d_p_parts, q2d_T_parts, q2d_umag_parts = [], [], []
    for zone in range(len(x_)):
        xn = np.asarray(x_[zone])
        yn = np.asarray(y_[zone])
        vn = np.asarray(var_[zone])
        xc = 0.25 * (xn[:-1, :-1, 0] + xn[1:, :-1, 0] +
                     xn[:-1, 1:, 0]  + xn[1:, 1:, 0])
        yc = 0.25 * (yn[:-1, :-1, 0] + yn[1:, :-1, 0] +
                     yn[:-1, 1:, 0]  + yn[1:, 1:, 0])
        q2d_xc_parts.append(xc.ravel())
        q2d_yc_parts.append(yc.ravel())
        # cellcentered indices: 0=rho1, 1=u, 2=v, 3=p, 4=T (Q2D 2D plane: no w)
        q2d_p_parts.append(vn[3, :, :, 0].ravel())
        q2d_T_parts.append(vn[4, :, :, 0].ravel())
        q2d_umag_parts.append(np.sqrt(vn[1, :, :, 0]**2 + vn[2, :, :, 0]**2).ravel())

    q2d_xc = np.concatenate(q2d_xc_parts)
    q2d_yc = np.concatenate(q2d_yc_parts)
    q2d_p = np.concatenate(q2d_p_parts)
    q2d_T = np.concatenate(q2d_T_parts)
    q2d_umag = np.concatenate(q2d_umag_parts)

    try:
        from scipy.interpolate import NearestNDInterpolator
        of_pts = np.column_stack([of_xc, of_yc])
        q2d_pts = np.column_stack([q2d_xc, q2d_yc])
        of_p_at = NearestNDInterpolator(of_pts, of_p)(q2d_pts)
        of_T_at = NearestNDInterpolator(of_pts, of_T)(q2d_pts)
        of_umag_at = NearestNDInterpolator(of_pts, of_umag)(q2d_pts)

        def l2_rel(a, b):
            return np.sqrt(np.mean((a - b)**2)) / np.sqrt(np.mean(b**2))

        def linf_rel(a, b):
            return np.max(np.abs(a - b)) / np.sqrt(np.mean(b**2))

        print('(1) Natural reference (OpenFOAM EnSight, REPORTED-only):')
        print(f'    Q2D cells: {len(q2d_xc)}  OF cells: {ne_of}')
        for name, qv, ov in [('p', q2d_p, of_p_at),
                             ('T', q2d_T, of_T_at),
                             ('|U|', q2d_umag, of_umag_at)]:
            print(f'    {name:<5} L2_rel={100*l2_rel(qv, ov):.2f}%  Linf_rel={100*linf_rel(qv, ov):.2f}%')
    except ImportError:
        print('(1) WARN: scipy not available; skipping OpenFOAM comparison')
else:
    print('(1) WARN: reference/OpenFOAM/ not found; skipping')

# =============================================================
# (2) MOSE bit-exact reference
# =============================================================
mose_ref = root.parents[3] / 'test' / 'regression_baseline' / 'MOSE_K0' / 'supersonic-forward-step' / 'field.tec'
diff_tool = root.parents[3] / 'tools' / 'diff_mose_q2d_field.py'
if not mose_ref.exists():
    print(f'(2) WARN: MOSE bit-exact reference not found at {mose_ref}; skipping')
elif not diff_tool.exists():
    print(f'(2) WARN: diff utility not found at {diff_tool}; skipping')
else:
    print('(2) MOSE bit-exact reference vs Q2D field (Linf_rel < 1e-9, L2_rel < 1e-12):')
    res = subprocess.run([sys.executable, str(diff_tool), str(mose_ref), 'OUTPUT/field.tec',
                          '--tol-linf', '1e-9', '--tol-l2', '1e-12'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    print(res.stdout)
    if res.returncode != 0:
        print(res.stderr)
        print('(2) FAIL: MOSE bit-exact mismatch'); sys.exit(1)
    print('(2) MOSE bit-exact reference: OK')

print('Correct!')
