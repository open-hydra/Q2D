# SSG/LRR-ω RSM — Per-Regime Code Verification

Two ctests exercise the SSG/LRR-ω Reynolds-stress model in the two
Q2D geometric-source regimes, each as a **kw/rsm sibling pair** (the
only input difference is the turbulence keyword: `kw` = Menter SST,
`rsm` = SSG/LRR-ω) sharing one family-root `verify.py`:

| ctest | Regime | Case directory | Tier |
|---|---|---|---|
| `ssglrr-chang-nozzle-turbulent` | 2D axisymmetric | `test/2D/viscous/chang-nozzle-turbulent/{kw,rsm}` | heavy (20 threads) |
| `ssglrr-planar-nozzle-turbulent` | Q2D area variation | `test/2D/euler/planar-nozzle-turbulent/{kw,rsm}` | fast (~2 min total) |

!!! warning "Code verification, not validation"

    These are **code-verification** cases, not validation: no
    independent axisymmetric / area-variation RSM reference exists.
    The model/term/discretisation correctness was established
    separately during development with symbolic (SymPy) regime-reduction
    and MMS term-level audits; the running cases verify
    that the regime-specific RSM machinery is actually **exercised**
    by a converged flow and behaves physically.  The validated RSM
    reference case remains the [2D planar turbulent flat
    plate](2D-flat-plate-turbulent.md) vs NASA TMR CFL3D.

**Feature exercised:** **SSG/LRR-ω Reynolds-stress model** (5-equation
RSM transport: $R_{11}, R_{22}, R_{33}, R_{12}, \omega$) vs the
Menter SST sibling in the two Q2D geometric-source regimes — the
axisymmetric pair exercises the azimuthal-stretching RSM production
term on the 2D-axi source, the planar pair the depth-compression term
on the Q2D area-variation source (Definition A
$\boldsymbol{V}\cdot\boldsymbol{G}$); RANS 12-arg diffusive-flux
dispatch; family-root `verify.py` gating (G1–G5).

## Results

<figure>
  {% include "vv/images/SSGLRR-axis.svg" %}
</figure>

## Gates G1–G5 (family-root `verify.py`)

Every gate is **hard** — any violation → exit 1; exit 0 only when all
five ran and passed.  Missing/empty required inputs
(`{kw,rsm}/OUTPUT/residual-history.dat`, `{kw,rsm}/OUTPUT/field.tec`,
and for the area-variation pair `rsm/INPUT/block1_area.dat`) → exit 1
with an explicit FAIL, never a silent pass.  All parsed fields are
NaN/Inf-guarded before the gates run; the PNG figure is best-effort
and never affects the exit code.

- **G1 — convergence**: both runs converged, density residual
  < **1e-8**.
- **G2 — realizability** everywhere on the converged RSM field:
  $R_{\alpha\alpha} \geq 0$ and
  $|R_{12}| \leq \sqrt{R_{11} R_{22}}$ (Schwarz inequality).
- **G3 — mean-flow consistency RSM vs SST**: with the derived inflow
  turbulence ($\mu_t/\mu \sim 2\times10^{-2}$) the mean flow is
  turbulence-insensitive, so the axis/centreline pressure profiles of
  the two siblings must agree: rel. $L_\infty$ < **1 %**.
- **G4 — regime-specific production term dynamically significant**:
  $\max |P| / (\rho\varepsilon) >$ **0.5** somewhere in the
  significant-$k$ region ($\rho\varepsilon = \rho C_\mu k \omega$,
  gradient-free, computed exactly from cell values):
    - axisymmetric pair: the azimuthal-stretching production
      $P_{\theta\theta} = -2\rho R_{33}\, u_r / r$;
    - area-variation pair: the depth-compression production
      $P_{33} = +2\rho R_{33} (\boldsymbol{V}\cdot\boldsymbol{G})$
      with $\boldsymbol{G} = -\nabla h / h$ from the committed area
      law `block1_area.dat` (the same $h(x)$ the solver uses).
- **G5 — anisotropy develops** from the exactly isotropic IC:
  $\max |a_{ii}| = \max |R_{ii}/k - 2/3| >$ **0.05** in the
  significant-$k$ region — the RSM discriminator vs the Boussinesq
  sibling, whose normal stresses stay isotropic.

The significance mask ($k > 10^{-2}\,k_\text{max}$) exists because
cells where the supersonic expansion has decayed the stresses toward
the realizability floor (1e-16) make $a_{ij}$ a 0/0 of floored
values; the unmasked extrema are still printed (disclosed, not
hidden).

## How to run

```bash
# area-variation pair (fast)
cd test/2D/euler/planar-nozzle-turbulent
(cd kw && ./Q2D.sh solve) && (cd rsm && ./Q2D.sh solve)
python3 -B verify.py

# axisymmetric pair (heavy: 20 threads)
cd test/2D/viscous/chang-nozzle-turbulent
(cd kw && ./Q2D.sh -p 20 solve) && (cd rsm && ./Q2D.sh -p 20 solve)
python3 -B verify.py
```

`verify.py` also writes a best-effort 2×2 diagnostic figure
(`rsm/OUTPUT/verify_RSM_*.png`: G3 pressure overlay, G5 anisotropy
profiles, G4 production-ratio map, G1 residual histories).
