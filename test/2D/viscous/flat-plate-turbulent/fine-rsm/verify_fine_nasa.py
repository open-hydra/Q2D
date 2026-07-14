"""Fine-grid (exact NASA 545x385) SSG/LRR-omega RSM V&V vs NASA TMR CFL3D.

Observables extracted from the converged Q2D run and plotted against the
respective NASA Turbulence Modeling Resource (tmbwg.github.io/turbmodels)
SSG/LRR-RSM-w2012 references:
  1. Skin friction Cf vs Re_x   -> NASA CFL3D verification point Cf(x=0.97008)=2.6899e-3
                                    (FlatPlate/SSGLRRRSM/cf_convergence_ssglrrrsm.dat, 545x385)
  2. Reynolds-stress anisotropy R_ij/k profiles at x-stations
                                 -> turb_contours_cfl3d_ssglrrrsm.dat (NASA CFL3D field)
  3. u+ vs y+ (law of the wall) -> theory (viscous sublayer u+=y+, log law 1/k ln y+ + B)

Run after the fine-rsm field.tec is converged:  python verify_fine_nasa.py
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    c = parent / "lib" / "ORION" / "src" / "python"
    if c.exists():
        sys.path.insert(0, str(c)); break
from ORION import read_TEC

REF = root / "reference"
MU = 1.1858685985e-5            # molecular viscosity (input.ini / transport.dat freestream)

# ---------- Q2D converged fine field ----------
x_, y_, _, var_, vn = read_TEC(str(root / "OUTPUT" / "field.tec"))
cc = vn[2:]
def C(n): return var_[0][cc.index('rho(1)' if n == 'rho1' else n)][:, :, 0]
xn = x_[0][:, :, 0]; yn = y_[0][:, :, 0]
rho = C('rho1'); u = C('u'); v = C('v'); p = C('p')
Ruu, Rvv, Rww, Ruv = C('ruu')/rho, C('rvv')/rho, C('rww')/rho, C('ruv')/rho
Ni, Nj = u.shape
# cell-centre coords
xc = 0.25*(xn[:-1, :-1]+xn[1:, :-1]+xn[:-1, 1:]+xn[1:, 1:])
yc = 0.25*(yn[:-1, :-1]+yn[1:, :-1]+yn[:-1, 1:]+yn[1:, 1:])
# freestream reference (top boundary, downstream)
i_fs, j_fs = int(Ni*0.7), Nj-2
rho_inf = rho[i_fs, j_fs]; u_inf = u[i_fs, j_fs]
q_inf = 0.5*rho_inf*u_inf**2
print(f"freestream: rho={rho_inf:.4f}  u={u_inf:.3f}  q_inf={q_inf:.2f}  Re/L={rho_inf*u_inf/MU:.3e}")

# ---------- (1) Cf vs Re_x ----------
xw = xc[:, 0]                       # streamwise (first cell row)
duwall = u[:, 0] / yc[:, 0]         # du/dy at wall (u_wall=0)
tau_w = MU * duwall
Cf = tau_w / q_inf
Re_x = rho_inf*u_inf*np.maximum(xw, 1e-9)/MU
plate = xw > 0.0
# NASA verification point
nasa_cf_x097 = 2.6898663500e-3
ix097 = np.argmin(np.abs(xw-0.970084071))
print(f"\n(1) Cf @ x=0.97008:  Q2D={Cf[ix097]:.5e}   NASA CFL3D(545x385)={nasa_cf_x097:.5e}"
      f"   rel.err={abs(Cf[ix097]-nasa_cf_x097)/nasa_cf_x097*100:.2f}%")

# ---------- (2) R_ij/k profiles vs NASA contour ----------
vals = []
for ln in open(REF/"turb_contours_cfl3d_ssglrrrsm.dat"):
    s = ln.strip()
    if not s or s[0] in '#"' or s.upper().startswith(('VARIABLES', 'ZONE', 'I=', 'DATAPACKING', 'DT=')):
        continue
    vals += [float(t) for t in s.split()]
I, J = 544, 384; N = I*J; vv = np.array(vals[:7*N]).reshape(7, N)
xN = vv[0].reshape(J, I); zN = vv[1].reshape(J, I)
R11N, R22N, R33N = vv[2].reshape(J, I), vv[3].reshape(J, I), vv[4].reshape(J, I)
# NASA z=wall-normal: map to Q2D R22q<-R33z(wall-normal), R33q<-R22y(spanwise)
stations = [0.5, 0.97008, 1.5]

# ---------- (3) u+ vs y+ at x=0.97008 ----------
utau = np.sqrt(tau_w[ix097]/rho[ix097, 0])
yp = yc[ix097, :]*utau*rho[ix097, :]/MU
up = u[ix097, :]/utau
yt = np.logspace(-1, 3.3, 200)
ut_law = np.where(yt < 11.0, yt, 1.0/0.41*np.log(yt)+5.0)

# ================= FIGURE =================
fig = plt.figure(figsize=(16, 5))
# Panel 1: Cf vs Re_x
ax1 = fig.add_subplot(1, 3, 1)
ax1.semilogx(Re_x[plate], Cf[plate], 'C0-', lw=2, label='Q2D SSG/LRR-$\\omega$ (545$\\times$385)')
Re097 = rho_inf*u_inf*0.970084071/MU
ax1.semilogx([Re097], [nasa_cf_x097], 'rs', ms=11, mfc='none', mew=2.2,
             label=f'NASA CFL3D x=0.97 ({nasa_cf_x097:.4e})')
ax1.set_xlabel('$Re_x$'); ax1.set_ylabel('$C_f$'); ax1.grid(alpha=0.3, which='both')
ax1.set_xlim(1e4, 6e6); ax1.set_ylim(0, 6e-3); ax1.legend(fontsize=8)
ax1.set_title(f'Skin friction — Q2D vs NASA CFL3D\n@x=0.97: {abs(Cf[ix097]-nasa_cf_x097)/nasa_cf_x097*100:.2f}% err')
# Panel 2: R_ij/k profiles at x=0.97008
ax2 = fig.add_subplot(1, 3, 2)
xs = 0.970084071
iq = np.argmin(np.abs(xc[:, 0]-xs)); ic = np.argmin(np.abs(xN[0, :]-xs))
kq = 0.5*(Ruu[iq, :]+Rvv[iq, :]+Rww[iq, :]); yq = yc[iq, :]
kN = 0.5*(R11N[:, ic]+R22N[:, ic]+R33N[:, ic]); zc_ = zN[:, ic]
for (q, n, col, mk, lab) in [(Ruu[iq, :]/kq, R11N[:, ic]/kN, 'C0', 'o', '$R_{11}/k$'),
                             (Rvv[iq, :]/kq, R33N[:, ic]/kN, 'C1', 's', '$R_{22}/k$ (wall-normal)'),
                             (Rww[iq, :]/kq, R22N[:, ic]/kN, 'C2', '^', '$R_{33}/k$ (spanwise)'),
                             (Ruv[iq, :]/kq, vv[5].reshape(J, I)[:, ic]/kN, 'C3', 'D', '$R_{12}/k$')]:
    ax2.plot(n, zc_, col+'-', lw=2)
    ax2.plot(q, yq, col+mk, ms=4, mfc='none', label='Q2D '+lab)
ax2.set_ylim(0, 0.03); ax2.set_xlabel('$R_{ij}/k$'); ax2.set_ylabel('y [m]')
ax2.grid(alpha=0.3); ax2.legend(fontsize=7, loc='upper right')
ax2.set_title('Reynolds-stress anisotropy x=0.97\nlines=NASA CFL3D, markers=Q2D')
# Panel 3: u+ vs y+
ax3 = fig.add_subplot(1, 3, 3)
ax3.semilogx(yt, ut_law, 'k--', lw=1.5, label='law of the wall ($u^+=y^+$; $\\frac{1}{0.41}\\ln y^+ +5$)')
m = yp > 0.1
ax3.semilogx(yp[m], up[m], 'C0o-', ms=4, mfc='none', label='Q2D SSG/LRR-$\\omega$')
ax3.set_xlabel('$y^+$'); ax3.set_ylabel('$u^+$'); ax3.grid(alpha=0.3, which='both')
ax3.set_xlim(0.1, 2e3); ax3.set_ylim(0, 28); ax3.legend(fontsize=8)
ax3.set_title(f'Law of the wall x=0.97  ($u_\\tau$={utau:.3f})')
plt.tight_layout()
out = root/"OUTPUT"/"verify_fine_nasa.png"
plt.savefig(out, dpi=140)
print(f"\nfigure -> {out}")
