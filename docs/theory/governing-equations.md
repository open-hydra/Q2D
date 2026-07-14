# Governing Equations

Q2D solves the compressible Navier–Stokes equations for multi-species
ideal-gas mixtures in conservative form on structured multi-block grids,
specialised for **three geometric regimes**: 2D planar, 2D axisymmetric
(without swirl), and q2D area-variation.  Both the inviscid (Euler) and
the viscous (Navier–Stokes) limits are supported.

The governing system is written in conservative form as

$$
\frac{\partial \mathbf{U}}{\partial t}
+ \nabla\!\cdot\!\mathbf{F}^{c}(\mathbf{U})
= \nabla\!\cdot\!\mathbf{F}^{v}(\mathbf{U})
+ \mathbf{S}
$$

where $\mathbf{U}$ is the vector of conservative variables,
$\mathbf{F}^{c}$ the convective (inviscid) flux,
$\mathbf{F}^{v}$ the viscous (diffusive) flux,
and $\mathbf{S}$ a source-term vector that aggregates chemistry,
turbulence, soot, and — in the axisymmetric and area-variation
regimes — geometric source terms.  The two non-planar regimes
are described in detail in [Geometric Source Terms](geometric-sources.md).

---

## State Vectors

The flow state is described by two complementary representations:
**primitive variables** for physical interpretation and I/O, and
**conservative variables** for the temporal integration and flux
evaluation in the solver.

### Primitive variables

The Q2D primitive state vector is

$$
\mathbf{P} = [\rho_1, \rho_2, \ldots, \rho_{N_s}, u, v, p, f_\text{soot}, f_\text{pass}, \mathbf{q}_\text{rans}]
$$

where the components are:

| Component | Symbol | Description |
|:-----------:|:------:|----------|
| Species densities | $\rho_s$ ($s=1,\ldots,N_s$) | Partial mass density of each species |
| Velocity | $u, v$ | In-plane Cartesian velocity components |
| Pressure | $p$ | Thermodynamic pressure |
| Soot variable | $f_\text{soot}$ | Soot volume fraction or total soot density (optional) |
| Passive scalar | $f_\text{pass}$ | Inert tracer field for mixing studies (optional) |
| RANS variables | $\mathbf{q}_\text{rans}$ | Turbulence model variables (optional) |

!!! note "Q2D is two-dimensional natively"

    Q2D is a two-dimensional solver: the primitive vector carries
    **only** $(u, v)$.  There is no out-of-plane velocity component
    as a primitive variable, no $k$-direction in the mesh, and
    no $\hat n_z$ in face metrics.  All routines are written
    natively in 2D — the formulation is **not** obtained by
    setting $w = 0$ in 3D code.  The dimensionality
    $\text{ndir} = 2$ is fixed at module level.

**Soot Model Variables (if active):**

- $f_\text{soot} = \rho_\text{soot} / \rho$: volume fraction of soot particles
- Additional soot moments or precursor scalars may be tracked

**RANS Model Variables (if active):**

- **Spalart–Allmaras**: $\tilde{\nu}$ (modified eddy viscosity)
- **$k$–$\omega$ (SST, Wilcox 2006)**: $k$ (turbulence kinetic energy),
  $\omega$ (specific dissipation rate)

### Conservative variables

Conservation laws require the flux-divergence form, which motivates the
conservative variable vector:

$$
\mathbf{U} = [\rho_1, \rho_2, \ldots, \rho_{N_s}, \rho\,u, \rho\,v, \rho\,E_0, f_\text{soot}, f_\text{pass}, \mathbf{q}_\text{rans}]
$$

where the first $N_s + 3$ components correspond to mass-momentum-energy
conservation:

$$
\rho = \sum_{s=1}^{N_s} \rho_s \quad \text{(total density)}
$$

$$
E_0 = e + \frac{1}{2}(u^2 + v^2) \quad \text{(total specific energy)}
$$

$$
e = h - \frac{p}{\rho} \quad \text{(internal specific energy)}
$$

The mixture specific enthalpy is computed as a mass-fraction-weighted
sum:

$$
h = \sum_{s=1}^{N_s} Y_s\, h_s(T), \quad Y_s = \frac{\rho_s}{\rho}
$$

where species enthalpies $h_s(T)$ are obtained from NASA polynomial
correlations (see [Thermodynamic and Transport Properties](thermo.md)).

**Passive and Turbulence Variables:**

Soot, passive scalars, and RANS variables do not enter the energy
equation or the equation of state.  They are advected with the bulk
flow and obey scalar transport equations with their own source terms.

### Transformation between conservative and primitive variables

Forward transformation ($\mathbf{P} \to \mathbf{U}$): species densities
remain unchanged, momentum is obtained as $\rho u, \rho v$, and the
conservative energy is computed by combining internal and kinetic
energy.

Inverse transformation ($\mathbf{U} \to \mathbf{P}$): recovering
primitive variables from conservative ones requires solving for
temperature implicitly via Newton–Raphson on the calorically imperfect
mixture EOS:

1. Extract densities and momentum:
   $\rho_s = U_s$ for each species; $\rho = \sum_s \rho_s$
2. Recover velocity:
   $(u, v) = (U_u, U_v)/\rho$
3. Compute specific internal energy:
   $e = U_E/\rho - \tfrac{1}{2}(u^2 + v^2)$
4. Solve for temperature iteratively:
   $T^{(k+1)} = T^{(k)} + (e - e(T^{(k)}))/c_v(T^{(k)})$
   where $e(T)$ and $c_v(T)$ are from NASA polynomial fits.
   Iterate until the relative error drops below $10^{-6}$.
5. Compute pressure from the ideal-gas EOS:
   $p = \rho\, R_\text{mix}\, T$
   where $R_\text{mix} = \sum_s Y_s\, R_s$ is the mass-fraction-weighted
   mixture gas constant.

---

## Convective Fluxes

The convective flux through a face of outward unit normal $\hat{\mathbf{n}}$
in the 2D plane is

$$
\mathbf{F}^{c}\!\cdot\!\hat{\mathbf{n}} = [
  \rho_1\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}),
  \ldots,
  \rho_{N_s}\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}),
  \rho\,u\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}) + p\,\hat{n}_x,
  \rho\,v\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}) + p\,\hat{n}_y,
  \rho\,H_0\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}})
]
$$

where $\mathbf{v} = (u, v)$, $\hat{\mathbf{n}} = (\hat n_x, \hat n_y)$,
and $H_0 = h + \tfrac{1}{2}\lvert\mathbf{v}\rvert^2$ is the total
specific enthalpy.  The normal mass flux is
$\dot{m} = \rho\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}})$.

---

## Viscous Fluxes

### Stress tensor

The viscous stress tensor follows the Newtonian fluid model with the
Stokes hypothesis $\lambda = -\tfrac{2}{3}\mu$:

$$
\tau_{ij} =
\mu\!\left(\frac{\partial v_i}{\partial x_j}
+ \frac{\partial v_j}{\partial x_i}\right)
- \frac{2}{3}\,\mu\,(\nabla\!\cdot\!\mathbf{v})\,\delta_{ij}
$$

where $\mu = \mu_\ell + \mu_t$ is the sum of laminar and eddy
viscosity.

For two-equation RANS models the isotropic part of the Reynolds stress
is included:

$$
\tau_{ij}^{\text{2eq}} = \tau_{ij} - \tfrac{2}{3}\,\rho\,k\,\delta_{ij}
$$

!!! info "Out-of-plane stress in the area-variation regime"

    In the area-variation regime the stress tensor is stored as a
    $3\times 3$ array even though the velocity vector is
    two-dimensional, because depth-averaging produces a non-zero
    out-of-plane component

    $$
    \tau_{zz} = -\tfrac{4}{3}\,\mu\,(\mathbf{V}\!\cdot\!\mathbf{G})
                - \tfrac{2}{3}\,\mu\,(\nabla_\text{2D}\!\cdot\!\mathbf{V})
    $$

    with $\mathbf{G} = -\nabla h/h$ the height-gradient vector.
    The full derivation is in
    [q2D Area-Variation Derivation](area-variation-derivation.md#63-theorem-63-depth-averaged-newton-stress).

### Species diffusion

Each species has a diffusive mass flux

$$
\mathbf{j}_s = -\rho\,D_{m,s}\,\nabla Y_s
$$

where the mixture-averaged diffusion coefficient is

$$
D_{m,s} = \frac{\mu_\ell}{\rho\,\text{Sc}}
         + \frac{\mu_t}{\rho\,\text{Sc}_t}
$$

A Hirschfelder correction velocity
$\mathbf{V}_c = \sum_s D_{m,s}\,\nabla Y_s$ enforces
$\sum_s \mathbf{j}_s^\ast = 0$ exactly; the corrected fluxes
$\mathbf{j}_s^\ast = \mathbf{j}_s + \rho Y_s\,\mathbf{V}_c$ are
used everywhere downstream, including the
$\sum_s h_s(\mathbf{j}_s^\ast)$ enthalpy diffusion term in the
energy equation.

### Heat conduction

The effective thermal conductivity is

$$
\kappa = k_\ell + \frac{\mu_t\,c_p}{\text{Pr}_t}
$$

### Diffusive flux vector

The complete viscous flux projected onto $\hat{\mathbf{n}}$ is

$$
\mathbf{F}^{v}\!\cdot\!\hat{\mathbf{n}} = [
  -\mathbf{j}_1\!\cdot\!\hat{\mathbf{n}},
  \ldots,
  -\mathbf{j}_{N_s}\!\cdot\!\hat{\mathbf{n}},
  \boldsymbol{\tau}\!\cdot\!\hat{\mathbf{n}},
  (\boldsymbol{\tau}\!\cdot\!\mathbf{v})\!\cdot\!\hat{\mathbf{n}}
  + \kappa\,\nabla T\!\cdot\!\hat{\mathbf{n}}
  + \rho\!\sum_s D_{m,s}\,(\nabla Y_s\!\cdot\!\hat{\mathbf{n}})\,h_s
]
$$

The energy flux therefore comprises viscous work, heat conduction, and
species enthalpy diffusion.

---

## Source Terms

Beyond convection and diffusion, the flow may be driven by source
terms aggregated into a single right-hand-side vector $\mathbf{S}$:

| Physical process | Affected variables | Integration strategy |
|---|---|---|
| Chemical kinetics | Species partial densities $\rho_s$ | Strang operator splitting; finite-rate or equilibrium |
| Turbulence closure | RANS variables $k$, $\omega$, $\tilde\nu$ | Explicit production-dissipation in spatial residual |
| Soot transport | Soot volume fraction $f_\text{soot}$ | Explicit source; nucleation, surface reactions, coagulation |
| Axisymmetric metric | All transport equations ($-1/r$ contributions) | See [Geometric Source Terms](geometric-sources.md#2-two-dimensional-axisymmetric-regime) |
| Area-variation metric | All transport equations ($\mathbf{V}\!\cdot\!\mathbf{G}$ contributions) | See [Geometric Source Terms](geometric-sources.md#3-quasi-two-dimensional-area-variation-regime) |

The chemical, turbulence, and soot sources are detailed in
[Reactive Kinetics](kinetics.md), [Turbulence Modelling](turbulence.md),
and the soot modelling section of the same.  The two
geometric-metric sources — present only in the non-planar regimes —
are derived from first principles in
[Geometric Source Terms](geometric-sources.md).
