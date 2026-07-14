"""Flat-plate turbulent boundary layer — SSG/LRR-omega RSM validation.

MODEL-CONSISTENT reference: this compares the Q2D SSG/LRR-omega result against the
NASA Turbulence Modeling Resource (TMR) CFL3D solution for the SAME model
(SSG/LRR-RSM-w2012 baseline), file reference/turb_contours_cfl3d_ssglrrrsm.dat
(downloaded from the NASA TMR; variables x, z, R11, R22, R33, R13, omega).

  Model identity verified (Sess 26) against the authoritative Eisfeld & Rumsey,
  "Length-Scale Correction for Reynolds-Stress Modeling", AIAA J. 2020, Appendix A
  Eqs 17-37 + Table 2: the Q2D implementation is the SSG/LRR-omega *baseline*
  (w2012); the w2019 length-scale correction is inactive on attached flows, so the
  w2012 CFL3D field is the matching reference. NOT compared against SA/SST data
  (different models).

NASA uses the (x,z) plane with z=wall-normal; Q2D uses (x,y) with y=wall-normal.
Index permutation:  NASA(R11, R22, R33, R13)  <->  Q2D(R11, R33, R22, R12).

The discriminating RSM quantity is the near-wall Reynolds-stress ANISOTROPY
R_ij/k (non-dimensional, grid/normalization independent) — a Boussinesq 2-eq model
gives isotropic normal stresses; the RSM reproduces the measured anisotropy.

Also checks the RSM physical invariants: realizability and the 7->5 reduction.

Usage:  python verify.py [--plot]
"""
import argparse
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    c = parent / "lib" / "ORION" / "src" / "python"
    if c.exists():
        sys.path.insert(0, str(c)); break
from ORION import read_TEC

ap = argparse.ArgumentParser()
ap.add_argument('--plot', action='store_true')
args = ap.parse_args()

XSTATION = 1.50   # developed station (both Q2D and NASA fully developed BL)

# ---------------- NASA CFL3D SSG/LRR reference (index-mapped to Q2D) ----------
ref = np.load(root / "reference" / "nasa_ssglrr_station.npz")
# NASA z=wall-normal; map to Q2D: R22q<-R33z(wall-normal), R33q<-R22y(spanwise), R12q<-R13
nasa = dict(z=ref['z'], R11=ref['R11'], R22=ref['R33z'], R33=ref['R22y'], R12=ref['R13'])
nasa_k = 0.5 * (nasa['R11'] + nasa['R22'] + nasa['R33'])
jp = int(np.argmax(nasa['R11']))
nasa_aniso = {c: float(nasa[c][jp] / nasa_k[jp]) for c in ('R11', 'R22', 'R33', 'R12')}
print("NASA CFL3D SSG/LRR-RSM near-wall anisotropy R_ij/k (Q2D convention):")
print(f"  R11/k={nasa_aniso['R11']:.3f}  R22/k={nasa_aniso['R22']:.3f}  "
      f"R33/k={nasa_aniso['R33']:.3f}  R12/k={nasa_aniso['R12']:.3f}")

# ---------------- Q2D field ---------------------------------------------------
x_, y_, z_, var_, vnames = read_TEC("OUTPUT/field.tec")
nd = 2 if x_[0].shape[2] == 1 else 3
cc = vnames[nd:]
def C(n): return var_[0][cc.index(n)]

for n in ('ruu', 'rvv', 'rww', 'ruv', 'omega'):
    assert n in cc, f"missing RSM var {n} in field.tec ({cc})"
assert 'ruw' not in cc and 'rvw' not in cc, "cross-stresses R13/R23 present — 7->5 violated"
print("\n7->5 reduction: only {ruu,rvv,rww,ruv,omega} transported  [PASS]")

rho = C('rho(1)') if 'rho(1)' in cc else C('rho1')
Ruu, Rvv, Rww, Ruv = C('ruu')/rho, C('rvv')/rho, C('rww')/rho, C('ruv')/rho

# realizability
neg = float(min(Ruu.min(), Rvv.min(), Rww.min()))
schwarz = float(np.max(np.abs(Ruv) - np.sqrt(np.maximum(Ruu*Rvv, 0.0))))
realizable = (neg >= -1e-12) and (schwarz <= 1e-8)
print(f"realizability: min(R_aa)={neg:.2e} (>=0), max(|R12|-sqrt(R11 R22))={schwarz:.2e} (<=0)  "
      f"[{'PASS' if realizable else 'FAIL'}]")

# Q2D anisotropy at the matched station, near-wall R11 peak
xc = 0.5 * (x_[0][:-1, 0, 0] + x_[0][1:, 0, 0]) if x_[0].shape[0] > Ruu.shape[0] else x_[0][:, 0, 0]
icol = int(np.argmin(np.abs(xc - XSTATION)))
# wall-normal profiles at the matched x-station (cell-centered, k=0 plane)
ruu_p, rvv_p, rww_p, ruv_p = (Ruu[icol, :, 0], Rvv[icol, :, 0], Rww[icol, :, 0], Ruv[icol, :, 0])
kq = 0.5 * (ruu_p + rvv_p + rww_p)
jpk = int(np.argmax(ruu_p))
q2d_aniso = {'R11': float(ruu_p[jpk]/kq[jpk]), 'R22': float(rvv_p[jpk]/kq[jpk]),
             'R33': float(rww_p[jpk]/kq[jpk]), 'R12': float(ruv_p[jpk]/kq[jpk])}
print(f"\nQ2D SSG/LRR near-wall anisotropy R_ij/k at x={xc[icol]:.3f}:")
print(f"  R11/k={q2d_aniso['R11']:.3f}  R22/k={q2d_aniso['R22']:.3f}  "
      f"R33/k={q2d_aniso['R33']:.3f}  R12/k={q2d_aniso['R12']:.3f}")

err = {c: abs(q2d_aniso[c] - nasa_aniso[c]) for c in q2d_aniso}
maxerr = max(err.values())
print(f"\nAnisotropy abs error vs NASA CFL3D: "
      + "  ".join(f"{c}={err[c]:.3f}" for c in err))
status = "PASS" if maxerr < 0.15 else "FAIL"
print(f"  max abs anisotropy error = {maxerr:.3f}   [{status}] (tol 0.15)")

# ---------------- SOTA V&V figure: Q2D overlaid on NASA CFL3D ------------------
# Wall-normal coordinate of the Q2D cell-centre profile at the station (physical y).
yq_n = y_[0][icol, :, 0]
if yq_n.shape[0] > kq.shape[0]:
    yq = 0.5 * (yq_n[:-1] + yq_n[1:])          # node -> cell centre
else:
    yq = yq_n
nz = nasa['z']
comps = [('R11', 'C0', 'o'), ('R22', 'C1', 's'), ('R33', 'C2', '^'), ('R12', 'C3', 'D')]
labels = {'R11': "$R_{11}/k$ (streamwise)", 'R22': "$R_{22}/k$ (wall-normal)",
          'R33': "$R_{33}/k$ (spanwise)", 'R12': "$R_{12}/k$ (shear)"}
q2d_prof = {'R11': ruu_p, 'R22': rvv_p, 'R33': rww_p, 'R12': ruv_p}

fig, ax = plt.subplots(1, 2, figsize=(14, 6))
# Panel 1: Reynolds-stress anisotropy profiles, Q2D (markers) over NASA (lines)
for c, col, mk in comps:
    ax[0].plot(nasa[c] / nasa_k, nz, col + '-', lw=2.4, label=f'NASA CFL3D {labels[c]}')
    ax[0].plot(q2d_prof[c] / kq, yq, col + mk, ms=5, mfc='none', mew=1.4,
               label=f'Q2D SSG/LRR-$\\omega$  {c}/k')
ax[0].set_ylim(0, 0.04)
ax[0].set_xlabel('$R_{ij}/k$  (Reynolds-stress anisotropy)')
ax[0].set_ylabel('wall-normal distance  $y$  [m]')
ax[0].legend(fontsize=7, ncol=2, loc='upper right')
ax[0].grid(alpha=0.3)
ax[0].set_title(f'SSG/LRR-$\\omega$ RSM flat plate, $x={xc[icol]:.2f}$ m — Q2D vs NASA TMR CFL3D')

# Panel 2: near-wall-peak anisotropy bar comparison (the RSM signature)
xb = np.arange(4); w = 0.36
nasa_b = [nasa_aniso[c] for c in ('R11', 'R22', 'R33', 'R12')]
q2d_b = [q2d_aniso[c] for c in ('R11', 'R22', 'R33', 'R12')]
ax[1].bar(xb - w/2, nasa_b, w, label='NASA CFL3D', color='0.5')
ax[1].bar(xb + w/2, q2d_b, w, label='Q2D SSG/LRR-$\\omega$', color='C0')
ax[1].set_xticks(xb)
ax[1].set_xticklabels(['$R_{11}/k$', '$R_{22}/k$', '$R_{33}/k$', '$R_{12}/k$'])
ax[1].axhline(0, color='k', lw=0.8)
ax[1].set_ylabel('near-wall anisotropy $R_{ij}/k$')
ax[1].legend(fontsize=9)
ax[1].grid(alpha=0.3, axis='y')
ax[1].set_title(f'Near-wall anisotropy — max |Δ| = {maxerr:.3f} ({status}, tol 0.15)')
for k_, (nb, qb) in enumerate(zip(nasa_b, q2d_b)):
    ax[1].text(k_, max(nb, qb) + 0.03, f'{abs(qb-nb):.3f}', ha='center', fontsize=8)

plt.tight_layout()
plt.savefig("OUTPUT/verify_RSM.png", dpi=130)
print("\nfigure -> OUTPUT/verify_RSM.png")
if args.plot:
    plt.show()

print(f"\nOVERALL: 7->5 PASS  realizability {'PASS' if realizable else 'FAIL'}  "
      f"anisotropy_vs_NASA {status}")
sys.exit(0 if (realizable and status == 'PASS') else 1)
