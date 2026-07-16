"""
Woodward-Colella forward step — Q2D vs OpenFOAM comparison

Loads:
  - Q2D output : OUTPUT/field.tec  (Tecplot BLOCK, 2 structured zones)
  - OpenFOAM ref: reference/OpenFOAM/forwardStep.0000.{mesh,p,T,U}
                  (EnSight Gold binary, single-layer hexa8 mesh)

Interpolates the OpenFOAM fields onto the Q2D cell centres with a
nearest-neighbour scheme, then computes relative L2 and L-inf errors
for pressure, temperature, and velocity magnitude.

Usage:
    python verify.py           # prints error table, saves figure to OUTPUT/verify.png
    python verify.py --plot    # shows interactive figure instead
"""

import argparse
import struct
import sys
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import NearestNDInterpolator
import warnings
warnings.filterwarnings("ignore")

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

output_dir = "OUTPUT"

# -------------------------------------------------------------
# Locate ORION Python package
# -------------------------------------------------------------
root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break

try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "ORION package not found. Run 'conda activate base' or "
        "install ORION in the current Python environment."
    )

# -------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Woodward-Colella forward step — Q2D vs OpenFOAM"
)
parser.add_argument("--plot", action="store_true", help="Show interactive plots")
args = parser.parse_args()

# =============================================================
# 1.  Load Q2D solution (OUTPUT/field.tec)
# =============================================================
# Variable layout: name-based via vnames (robust to Q2D 2D 'rho,u,v,p,T,g,R'
# vs the 3D peer's 'rho,u,v,w,p,T,g,R'). ORION read_TEC + VARLOCATION → ndir.
x_, y_, z_, var_, vnames = read_TEC("OUTPUT/field.tec")
_ndir = 2 if x_[0].shape[2] == 1 else 3
_cc_names = vnames[_ndir:]
_iu = _cc_names.index('u')
_iv = _cc_names.index('v')
_iw = _cc_names.index('w') if 'w' in _cc_names else None   # Q2D 2D: no w (2D no-w layout)
_ip = _cc_names.index('p')
_iT = _cc_names.index('T')

q2d_xc_parts, q2d_yc_parts = [], []
q2d_p_parts, q2d_T_parts, q2d_umag_parts = [], [], []

for zone in range(len(x_)):
    xn = np.asarray(x_[zone])      # (I, J, K) – nodal coordinates
    yn = np.asarray(y_[zone])
    vn = np.asarray(var_[zone])    # (Nvar, I-1, J-1, max(K-1,1)) – cell-centred

    # 2-D cell centres (face-averaged from k=0 nodal layer)
    xc = 0.25 * (xn[:-1, :-1, 0] + xn[1:, :-1, 0] +
                 xn[:-1,  1:, 0] + xn[1:,  1:, 0])   # (I-1, J-1)
    yc = 0.25 * (yn[:-1, :-1, 0] + yn[1:, :-1, 0] +
                 yn[:-1,  1:, 0] + yn[1:,  1:, 0])

    q2d_xc_parts.append(xc.ravel())
    q2d_yc_parts.append(yc.ravel())
    q2d_p_parts.append(vn[_ip, :, :, 0].ravel())
    q2d_T_parts.append(vn[_iT, :, :, 0].ravel())
    _u2 = vn[_iu, :, :, 0]**2 + vn[_iv, :, :, 0]**2
    if _iw is not None:
        _u2 = _u2 + vn[_iw, :, :, 0]**2
    q2d_umag_parts.append(np.sqrt(_u2).ravel())

q2d_xc   = np.concatenate(q2d_xc_parts)
q2d_yc   = np.concatenate(q2d_yc_parts)
q2d_p    = np.concatenate(q2d_p_parts)
q2d_T    = np.concatenate(q2d_T_parts)
q2d_umag = np.concatenate(q2d_umag_parts)

# =============================================================
# 2.  Load OpenFOAM EnSight Gold reference
# =============================================================
REF_DIR = Path("reference/OpenFOAM")


def _es_str(f):
    """Read one 80-byte EnSight Gold string record."""
    return f.read(80).decode("ascii").strip("\x00").strip()


def read_ensight_mesh(path):
    """Parse EnSight Gold binary geometry -> node coords + hexa8 connectivity."""
    with open(path, "rb") as f:
        for _ in range(5):
            _es_str(f)                          # C binary / desc / author / node_id / elem_id
        _es_str(f)                              # "part"
        struct.unpack("i", f.read(4))[0]        # part number
        _es_str(f)                              # part description
        _es_str(f)                              # "coordinates"
        nn = struct.unpack("i", f.read(4))[0]
        xn   = np.frombuffer(f.read(4 * nn), dtype=np.float32)
        yn   = np.frombuffer(f.read(4 * nn), dtype=np.float32)
        f.read(4 * nn)                          # z coords — ignored (2-D problem)
        _es_str(f)                              # "hexa8"
        ne   = struct.unpack("i", f.read(4))[0]
        conn = np.frombuffer(f.read(8 * 4 * ne), dtype=np.int32).reshape(ne, 8)
    return xn, yn, conn, ne


def read_ensight_scalar(path, ne):
    """Parse EnSight Gold per-element scalar field."""
    with open(path, "rb") as f:
        _es_str(f)                              # field description
        _es_str(f)                              # "part"
        struct.unpack("i", f.read(4))[0]        # part number
        _es_str(f)                              # element type
        data = np.frombuffer(f.read(4 * ne), dtype=np.float32)
    return data.astype(float)


def read_ensight_vector(path, ne):
    """Parse EnSight Gold per-element vector field -> (3, ne)."""
    with open(path, "rb") as f:
        _es_str(f)
        _es_str(f)
        struct.unpack("i", f.read(4))[0]
        _es_str(f)
        data = np.frombuffer(f.read(3 * 4 * ne), dtype=np.float32).reshape(3, ne)
    return data.astype(float)


xn_of, yn_of, conn_of, ne_of = read_ensight_mesh(REF_DIR / "forwardStep.0000.mesh")

of_p    = read_ensight_scalar(REF_DIR / "forwardStep.0000.p", ne_of)
of_T    = read_ensight_scalar(REF_DIR / "forwardStep.0000.T", ne_of)
of_U    = read_ensight_vector(REF_DIR / "forwardStep.0000.U", ne_of)

# Element centroids (connectivity is 1-indexed)
of_xc   = xn_of[conn_of - 1].mean(axis=1)
of_yc   = yn_of[conn_of - 1].mean(axis=1)
of_umag = np.sqrt(of_U[0]**2 + of_U[1]**2 + of_U[2]**2)

# =============================================================
# 3.  Interpolate OpenFOAM fields onto Q2D cell centres
# =============================================================
of_pts = np.column_stack([of_xc, of_yc])

interp_p    = NearestNDInterpolator(of_pts, of_p)
interp_T    = NearestNDInterpolator(of_pts, of_T)
interp_umag = NearestNDInterpolator(of_pts, of_umag)

q2d_pts = np.column_stack([q2d_xc, q2d_yc])
of_p_at_q2d    = interp_p(q2d_pts)
of_T_at_q2d    = interp_T(q2d_pts)
of_umag_at_q2d = interp_umag(q2d_pts)

# =============================================================
# 4.  Error metrics
# =============================================================
def l2_rel(a, b):
    """Relative L2 error: ||a-b||_2 / ||b||_2."""
    return np.sqrt(np.mean((a - b) ** 2)) / np.sqrt(np.mean(b ** 2))


def linf_rel(a, b):
    """Relative L-inf error: max|a-b| / RMS(b)."""
    return np.max(np.abs(a - b)) / np.sqrt(np.mean(b ** 2))


print("Woodward-Colella forward step — Q2D vs OpenFOAM")
print("=" * 52)
print(f"  Q2D cells    : {len(q2d_xc):>6d}")
print(f"  OpenFOAM cells: {ne_of:>5d}")
print()
print(f"  {'Variable':<10}  {'L2 rel [%]':>12}  {'Linf rel [%]':>14}")
print("  " + "-" * 40)
for vname, mv, ov in [
    ("p",   q2d_p,    of_p_at_q2d),
    ("T",   q2d_T,    of_T_at_q2d),
    ("|U|", q2d_umag, of_umag_at_q2d),
]:
    print(f"  {vname:<10}  {100*l2_rel(mv,ov):>12.2f}  {100*linf_rel(mv,ov):>14.2f}")
print()

# =============================================================
# 5.  Plots
# =============================================================
fig, axes = plt.subplots(3, 2, figsize=(14, 9))

# Helper: rebuild structured arrays for one Q2D zone
def _zone_arrays(zone_idx, var_idx):
    xn = np.asarray(x_[zone_idx])
    yn = np.asarray(y_[zone_idx])
    vn = np.asarray(var_[zone_idx])
    return xn[:, :, 0], yn[:, :, 0], vn[var_idx, :, :, 0]


STEP_X = [0.0, 0.6, 0.6]
STEP_Y = [0.2, 0.2, 0.0]

for row, (vname, var_idx, of_vals) in enumerate([
    ("p",   _ip, of_p),
    ("T",   _iT, of_T),
    ("|U|", None, of_umag),
]):
    # Collect Q2D cell centres and values for scatter plot
    q2d_xc_all, q2d_yc_all, q2d_vals_all = [], [], []
    for zone in range(len(x_)):
        xn = np.asarray(x_[zone])
        yn = np.asarray(y_[zone])
        vn = np.asarray(var_[zone])
        xc = 0.25 * (xn[:-1, :-1, 0] + xn[1:, :-1, 0] +
                     xn[:-1,  1:, 0] + xn[1:,  1:, 0])
        yc = 0.25 * (yn[:-1, :-1, 0] + yn[1:, :-1, 0] +
                     yn[:-1,  1:, 0] + yn[1:,  1:, 0])
        if var_idx is not None:
            vals = vn[var_idx, :, :, 0]
        else:
            _u2 = vn[_iu,:,:,0]**2 + vn[_iv,:,:,0]**2
            if _iw is not None:
                _u2 = _u2 + vn[_iw,:,:,0]**2
            vals = np.sqrt(_u2)
        q2d_xc_all.append(xc.ravel())
        q2d_yc_all.append(yc.ravel())
        q2d_vals_all.append(vals.ravel())
    
    q2d_xc_plot = np.concatenate(q2d_xc_all)
    q2d_yc_plot = np.concatenate(q2d_yc_all)
    q2d_vals_plot = np.concatenate(q2d_vals_all)
    
    # shared colour scale
    vmin = min(q2d_vals_plot.min(), of_vals.min())
    vmax = max(q2d_vals_plot.max(), of_vals.max())

    # ---- Q2D scatter plot ----
    ax = axes[row, 0]
    sc = ax.scatter(q2d_xc_plot, q2d_yc_plot, c=q2d_vals_plot, s=4, cmap="turbo", vmin=vmin, vmax=vmax,
                    linewidths=0, rasterized=True)
    fig.colorbar(sc, ax=ax, label=vname)
    ax.axvline(0.6, color="gray", lw=0.8, ls=":", alpha=0.7)
    ax.set_title(f"Q2D — {vname}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1)

    # ---- OpenFOAM scatter plot ----
    ax = axes[row, 1]
    sc = ax.scatter(of_xc, of_yc, c=of_vals, s=4, cmap="turbo", vmin=vmin, vmax=vmax,
                    linewidths=0,rasterized=True)
    fig.colorbar(sc, ax=ax, label=vname)
    ax.axvline(0.6, color="gray", lw=0.8, ls=":", alpha=0.7)
    ax.set_title(f"OpenFOAM — {vname}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1)

plt.tight_layout()

if args.plot:
    plt.show()
else:
    out_path = os.path.join(output_dir, f"WC-fields.svg")
    plt.savefig(out_path, bbox_inches="tight", transparent=True)
    print(f"  Figure saved to {out_path}")

# =============================================================
# 6.  Density profile along y = 0.5 horizontal line
# =============================================================
Y_SLICE = 0.5

# --- Q2D: collect cells from both zones nearest to y = Y_SLICE ---
q2d_x_slice, q2d_rho_slice = [], []
for zone in range(len(x_)):
    xn = np.asarray(x_[zone])
    yn = np.asarray(y_[zone])
    vn = np.asarray(var_[zone])
    xc = 0.25 * (xn[:-1, :-1, 0] + xn[1:, :-1, 0] +
                 xn[:-1,  1:, 0] + xn[1:,  1:, 0])   # (I-1, J-1)
    yc = 0.25 * (yn[:-1, :-1, 0] + yn[1:, :-1, 0] +
                 yn[:-1,  1:, 0] + yn[1:,  1:, 0])
    # For each column i, pick the row j closest to Y_SLICE
    j_best = np.argmin(np.abs(yc[0, :] - Y_SLICE))
    # Only include if y is within half a cell height
    dy = np.abs(yc[0, j_best] - Y_SLICE)
    dy_cell = np.abs(yc[0, min(j_best + 1, yc.shape[1] - 1)] - yc[0, j_best])
    if dy <= dy_cell:
        q2d_x_slice.append(xc[:, j_best])
        _irho = _cc_names.index('rho(1)') if 'rho(1)' in _cc_names else _cc_names.index('rho')
        q2d_rho_slice.append(vn[_irho, :, j_best, 0])

q2d_x_slice   = np.concatenate(q2d_x_slice)
q2d_rho_slice = np.concatenate(q2d_rho_slice)
sort_idx = np.argsort(q2d_x_slice)
q2d_x_slice   = q2d_x_slice[sort_idx]
q2d_rho_slice = q2d_rho_slice[sort_idx]

# --- OpenFOAM: cells in a band around Y_SLICE, then sort by x ---
R_gas = 1.0 / 1.4          # same non-dimensional R used by both solvers
of_rho = of_p / (R_gas * of_T)

# use a band of half the typical cell height (~1/128 of domain height)
band = 1.0 / 128.0
mask_of = np.abs(of_yc - Y_SLICE) <= band
# average duplicate x positions (cells at same x, different z layer)
of_x_raw   = of_xc[mask_of]
of_rho_raw = of_rho[mask_of]
sort_of = np.argsort(of_x_raw)
of_x_slice   = of_x_raw[sort_of]
of_rho_slice = of_rho_raw[sort_of]

fig2, ax2 = plt.subplots(figsize=(10, 4))
ax2.plot(q2d_x_slice, q2d_rho_slice, label="Q2D",     lw=2.5)
ax2.plot(of_x_slice,   of_rho_slice,   label="OpenFOAM", lw=2.5, ls="--")
ax2.set_xlabel("x")
ax2.set_ylabel(r"$\rho$")
ax2.set_xlim(0.0, 3.0)
ax2.legend()
ax2.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()

if args.plot:
    plt.show()
else:
    out_path2 = os.path.join(output_dir, f"WC-slice.svg")
    fig2.savefig(out_path2, bbox_inches="tight", transparent=True)
    print(f"  Figure saved to {out_path2}")