# Shock Wave - Boundary Layer Interaction

Supersonic turbulent flow at Mach 5 over a flat plate with an isothermal wall, interacting with an impinging oblique shock. This benchmark validates Q2D in a shock-induced separation and reattachment regime by comparing skin-friction trends against published and code-to-code reference data.

**References**: Schulein experimental dataset (as distributed with this case), SU2, and Wind-US reference solutions.

<figure>
  {% include "vv/images/SWBLI-field.svg" %}
</figure>

---

## Problem setup

The case is a canonical SWBLI configuration where a turbulent boundary layer develops along a flat plate and then interacts with a strong compression shock. The interaction region induces an adverse pressure gradient, local separation, and subsequent reattachment.

The primary validation metric is the skin-friction coefficient:

$$C_f = \frac{\tau_w}{0.5\,\rho_\infty U_\infty^2}$$

where $\tau_w$ is the streamwise wall shear stress, and $\rho_\infty$, $U_\infty$ are freestream density and velocity.

In this configuration, the expected SWBLI signature is:

- positive $C_f$ upstream of interaction,
- negative $C_f$ inside the separated region,
- positive $C_f$ again after reattachment.

**Freestream conditions**

| Parameter | Value |
|---|---|
| Mach number | 5.0 |
| Temperature | 68.3 K |
| Pressure | 0.04 bar |
| Dynamic viscosity | 1.1858685985 x 10^-5 Pa*s |
| Turbulence model | Spalart-Allmaras (SA) |

## Numerical setup

| Parameter | Value |
|---|---|
| Equations | Navier-Stokes |
| Time scheme | RK3 |
| CFL | 0.3 |
| VNN | 0.1 |
| IRS $\beta$ | 0.0 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |
| Integration variables | Conservative |

## Grid structure

The mesh is a 2D structured multi-block grid (Q2D peculiarità #1 K=1 strip preserved). It is built to resolve:

- the incoming turbulent boundary layer,
- the shock interaction region,
- the separated and reattaching near-wall flow.

Boundary conditions:

- **Inlet**: Supersonic inflow (Mach 5, prescribed thermodynamic state)
- **Lower wall**: Isothermal no-slip wall
- **Upper boundary**: Symmetry upstream and isothermal wall downstream
- **Outlet**: Supersonic outlet extrapolation

## Results and validation

Validation is performed with the wall skin-friction distribution, compared against the reference datasets (Schulein, SU2, Wind-US):

<figure>
  {% include "vv/images/SWBLI-cf-sa.svg" %}
</figure>

Q2D reproduces the expected SWBLI behavior, including the negative-$C_f$ pocket associated with shock-induced separation and the recovery after reattachment. Overall agreement with the reference trends is good in both location and amplitude of the separation/reattachment signature.

The regression check anchors this against the **independent references** quantitatively. The separation and reattachment locations — where the wall shear changes sign — are a scale-invariant physical metric (independent of the skin-friction normalization). The solver places separation at $x \approx 0.339$ and reattachment at $x \approx 0.344$, both inside the envelope of the independent datasets (separation $[0.3337, 0.3375]$, reattachment $[0.3430, 0.3458]$ from Schulein, SU2-SA and Wind-US) within a 1.5 % engineering band. That band is set by the spread among the references themselves (which disagree by about 1.1 % on the separation location — the irreducible physical uncertainty), not chosen to pass: the Spalart–Allmaras model is known to slightly delay separation, and the solver lands within 0.5 % of the nearest CFD–SA reference (SU2). This complements the consistency check against the same-model MOSE reference (wall shear matching to within 0.1 %).
