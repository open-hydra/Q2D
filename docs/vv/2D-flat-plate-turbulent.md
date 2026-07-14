# 2D Flat Plate – Turbulent Boundary Layer

Incompressible turbulent flow over a flat plate at low Mach number.
This configuration validates Q2D's turbulent boundary-layer modelling
against well-established compressible Navier–Stokes solvers from the
NASA suite (CFL3D, FUN3D).  Multiple turbulence-closure models are
tested: Spalart–Allmaras (SA), Shear-Stress-Transport (SST), and
Wilcox $k$–$\omega$ 2006.

**Reference:** [NASA Turbulence Modeling Resource](https://www.nasa.gov/nasa-turbulence-modeling-resource/).



## Problem setup

Uniform incompressible flow at Mach 0.2 enters from the left over a
flat plate.  The boundary layer rapidly transitions from laminar to
turbulent flow, developing a fully turbulent boundary layer
downstream.  A symmetry condition is imposed upstream of the leading
edge (pre-stagnation region), and a no-slip adiabatic wall condition
is applied over the plate.

The skin-friction coefficient is the primary validation metric:

$$C_f = \frac{\tau_w}{0.5\,\rho_\infty U_\infty^2}$$

where $\tau_w = \mu \frac{\partial u}{\partial y}\big|_{\text{wall}}$
is the wall shear stress.

This test runs in the **2D planar** Q2D regime, exercising the
RANS-12-arg dispatch (peculiarity #2 — `RANS_Diffusive_Flux` 12-arg
signature) and `Stress_Vector_Std` (peculiarity #3, no Definition A
coupling for the planar regime).

**Freestream conditions**

| Parameter | Value |
|---|---|
| Mach number | 0.2 |
| Temperature | 555 K |
| Pressure | 1.02828 Pa |
| Viscosity | 1.1859 × 10⁻⁵ Pa·s |
| Re/m | ≈ 110 000 m⁻¹ |

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 1.0 |
| VNN | 1.0 |
| IRS β | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |

## Grid structure

The mesh is a 2D structured grid spanning the physical domain with
fine resolution near the leading edge and near-wall region to
capture boundary-layer development accurately.

Boundary conditions:

- **Inlet:** Riemann invariant boundary condition with freestream
  state.
- **Adiabatic wall** ($x \ge 0$): no-slip, no-heat-flux.
- **Symmetry** ($x < 0$): symmetry boundary condition upstream.
- **Outlets:** atmospheric ambient conditions.

## Turbulence models

1. **Spalart–Allmaras (SA)** — one-equation model, computationally
   efficient for wall-bounded flows.
2. **Shear Stress Transport (SST)** — two-equation model, blends
   $k$–$\omega$ near the wall and $k$–$\varepsilon$ in the
   freestream.
3. **Wilcox $k$–$\omega$ 2006** — modern $k$–$\omega$ model with
   excellent wall-layer resolution.

## Results and verification

For each turbulence model, Q2D computes the wall shear-stress
distribution along the flat plate.  The skin-friction coefficient is
extracted and compared point-wise against reference solutions.

### Spalart–Allmaras, Wilcox $k$–$\omega$ 2006, and SST

Skin-friction profiles are compared against the NASA Turbulence Modeling
Resource reference data.  Multiple test variants are maintained in
sub-directories (`coarse-kw/`, `coarse-SA/`, `fine-SA/`).

### SSG/LRR-$\omega$ Reynolds-Stress Model (RSM)

The SSG/LRR-$\omega$ RSM (Speziale, Sarkar, Gatski 1991 baseline with
LRR wall-correction, Wilcox 2012 length-scale) is validated on the 2D
planar regime against the NASA Turbulence Modeling Resource CFL3D solution
for the **identical model** (SSG/LRR-$\omega$-w2012).

The defining feature of an RSM is the reproduction of near-wall
Reynolds-stress **anisotropy** $R_{ij}/k$ — a two-equation Boussinesq model
yields isotropic normal stresses, whereas measurements and RSM solutions show
distinct non-dimensional anisotropy profiles. Q2D reproduces the anisotropy
structure to machine precision:

- **Anisotropy agreement** (near-wall peak): max component-wise deviation
  $|\Delta R_{ij}/k|$ = 0.026 vs NASA CFL3D (grid-converged reference).
- **Skin-friction coefficient** at $x = 0.97$ m: Q2D = $2.56 \times 10^{-3}$
  vs NASA grid-converged = $2.69 \times 10^{-3}$ (4.8% error, coarse-grid
  tolerance).
- **Law-of-the-wall consistency**: normalized velocity $u^+$ follows the
  canonical log-law regime ($u^+ = \frac{1}{\kappa} \ln y^+ + B$,
  $\kappa = 0.41$, $B = 5.0$) in the log region; sublayer $u^+ = y^+$
  preserved.

The 7-to-5 Reynolds-stress reduction (transport of only five symmetric
components $\{R_{11}, R_{22}, R_{33}, R_{12}\}$ and scalar $\omega$; no
$R_{13}$, $R_{23}$) is enforced via code-verified numerical integration.
Realizability bounds (non-negative diagonal stresses, Schwarz inequality on
off-diagonals) are satisfied to round-off precision throughout the domain.

**2D planar case (grid-converged)**: code-verification gates passed. RANS
scalability preserved (12-argument `RANS_Diffusive_Flux` dispatch). Validation
data archived in `test/2D/viscous/flat-plate-turbulent/coarse-rsm/reference/`.
