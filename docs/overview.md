---
title: Overview
---

# Overview

Q2D (Quasi-2D solver for compressible flow — 2D planar, 2D axisymmetric, and Q2D area-variation) is an open-source compressible-flow solver for the Euler and Navier-Stokes equations on multi-block structured grids, written in modern Fortran. It targets a wide range of inviscid and viscous compressible flow problems — from classical shock-tube benchmarks to turbulent reacting flows.

---

## Hydra CFD Suite

Q2D is the **ideal-gas compressible-flow solver** of the **Hydra** CFD ecosystem — an integrated suite of tools for multi-physics simulation of complex systems.

| Component | Role | Status |
|-----------|------|--------|
| [**ATLAS**](https://github.com/open-hydra/ATLAS) | Pre-processor: mesh prep, initial & boundary conditions, thermodynamic/chemical data | Separate package |
| **Q2D** | Solver: compressible Euler/Navier–Stokes with finite-rate chemistry | This package |

!!! info "Using Q2D without ATLAS"
    The input files required by Q2D (initial conditions, boundary conditions, thermodynamic tables, transport and chemistry data) are **typically produced by ATLAS**. If ATLAS is not available, all input files can be prepared manually; see the [User Guide](user/using.md) for the expected formats.

---

## Q2D Capabilities

Q2D solves the compressible Navier–Stokes equations for multi-species ideal-gas mixtures in conservative form:

$$
\frac{\partial \mathbf{U}}{\partial t}
+ \nabla\!\cdot\!\mathbf{F}^{c}(\mathbf{U})
= \nabla\!\cdot\!\mathbf{F}^{v}(\mathbf{U})
+ \mathbf{S}
$$

where $\mathbf{U}$ is the conservative-variable vector, $\mathbf{F}^c$ convective fluxes, $\mathbf{F}^v$ viscous fluxes, and $\mathbf{S}$ source terms (chemistry, turbulence, rotating frame). Both the inviscid limit (Euler equations) and the full viscous form (Navier–Stokes) are supported.

---

## Physical Models

### Thermally perfect gas

Q2D models the gas mixture as a collection of $N_s$ **thermally perfect** species, each obeying the ideal-gas equation of state $p = \rho R_\text{mix} T$. All thermodynamic and transport properties are **temperature-dependent** and evaluated from tabulated data (NASA 7/9-coefficient polynomials) via the [FLINT](https://github.com/MarcoGrossi92/FLINT) library. Mixture quantities are formed by mass-weighted averaging.

| Property | Model |
|----------|-------|
| Equation of state | Thermally perfect gas ($p = \rho R_\text{mix} T$) |
| Specific heat $c_p(T)$, enthalpy $h(T)$, entropy $s(T)$ | NASA polynomial tables (FLINT) |
| Dynamic viscosity $\mu(T)$ | Species-level kinetic-theory fits; Wilke mixture rule |
| Thermal conductivity $\lambda(T)$ | Species-level kinetic-theory fits; Wilke mixture rule |
| Species diffusion | Fick's law with constant Schmidt number |
| Back-end (optional) | [Cantera](https://cantera.org) integration for extended evaluation |

### Chemistry

Q2D supports three chemistry modes, selected via
`[Q2D-Physics] chemistry = frozen | equilibrium | finite-rate`
(default `frozen`):

**Frozen composition** — the species mass fractions are transported
but do not react (the default; multi-species thermodynamics only).

**Finite-rate kinetics** — species mass fractions evolve by elementary, three-body, and pressure-dependent reactions following the modified Arrhenius law $k_f = A\,T^b\,\exp(-E_a/R_u T)$. The kinetics framework is provided by [FLINT](https://github.com/MarcoGrossi92/FLINT). Stiff source terms are integrated with Strang operator splitting via the [OSlo](https://github.com/MarcoGrossi92/OSlo) ODE library.

**Chemical equilibrium** — NASA CEA integration for flows where the chemistry time scale is much shorter than the flow time scale (e.g. rocket nozzle expansions).

**Soot formation** — Leung–Lindstedt (1991) formation/oxidation model (`LL91`) and a linear formation/oxidation model (`LIN`), selected via `[Q2D-Physics] soot-generation = LL91 | LIN | none`.

| Capability | Details |
|------------|---------|
| Reaction types | Elementary, three-body, pressure-dependent (Lindemann/Troe) |
| Rate law | Modified Arrhenius; forward–backward via equilibrium constants |
| Stiff integration | Strang splitting + [OSlo](https://github.com/MarcoGrossi92/OSlo) / SUNDIALS |
| Equilibrium | NASA CEA back-end |
| Soot | `LL91` (Leung–Lindstedt 1991) and `LIN` (linear) models |

### Turbulence

All turbulence models are RANS closures selected at run time. The eddy-viscosity models link the Reynolds stress tensor to the mean strain rate via the Boussinesq hypothesis; the SSG/LRR-$\omega$ Reynolds-stress model transports the stresses directly.

| Model | Type | Key feature |
|-------|------|-------------|
| Spalart–Allmaras (SA) | 1-equation | Variants: `SA-R`, `SA-RC` (rotation/curvature), `SAcomp` |
| Menter SST | 2-equation | $k$–$\omega$ / $k$–$\varepsilon$ blending; robust in adverse pressure gradients |
| Wilcox k-$\omega$ 2006 | 2-equation | Updated cross-diffusion; improved freestream sensitivity |
| SSG/LRR-$\omega$ (`SSGLRR`, `SSGLRR-SD`) | Reynolds-stress | Speziale–Sarkar–Gatski / Launder–Reece–Rodi blend; `-SD` shock-detector diffusive flux |
| QCR2000 | Algebraic correction | Non-linear Boussinesq (Quadratic Constitutive Relation), `-QCR2000` model suffix |

The `-QCR2000` and `-blowcorr` (wall-blowing correction) suffixes are
parsed as substrings of the `[Q2D-Physics] turbulence` model string
(e.g. `turbulence = SA-RC-QCR2000`); the
[registry](user/registry.md#q2d-physics) lists the base model values.

### Gas-Surface Interaction (GSI)

The GSI model enables **reactive wall boundary conditions** for high-enthalpy flows where the gas composition at the wall changes due to heterogeneous reactions (catalysis, ablation). Species mass fluxes at the surface are computed from user-specified catalytic efficiencies or recombination probabilities.

GSI is **BC-code driven, not registry-exposed**: it is activated through the wall BC codes `303`/`304` in `INPUT/bc.txt` (see [Boundary Conditions](user/boundary-conditions.md)), has no `[Q2D-*]` INI keys yet, and is not exercised by the ctest suite (see the [Coverage Matrix](vv/coverage-matrix.md)).

| Capability | Details |
|------------|---------|
| Catalytic walls | Fully catalytic, partially catalytic, or non-catalytic |
| Species fluxes | Computed from surface reaction rates and diffusion |
| Coupling | Implicit coupling with the flow solver at each time step |

### Rotating Frame (SRF)

The **Single Rotating Frame** model re-casts the Navier–Stokes equations in a reference frame rotating rigidly at constant angular velocity $\boldsymbol{\Omega}$. The solver stores and evolves the **relative velocity** $\boldsymbol{w} = \boldsymbol{u}_\text{abs} - \boldsymbol{\Omega}\times\boldsymbol{r}$, so the convective and diffusive flux modules require no modification.

Additional source terms account for the inertial forces:

$$
\boldsymbol{S}_\text{mom} = -2\rho(\boldsymbol{\Omega}\times\boldsymbol{w}) + \rho\,\omega^2\boldsymbol{r}_\perp
$$

| Term | Expression |
|------|-----------|
| Coriolis | $-2\rho(\boldsymbol{\Omega}\times\boldsymbol{w})$ |
| Centrifugal | $+\rho\,\omega^2\boldsymbol{r}_\perp$ |
| Energy source | Modified rothalpy conservation |

Typical applications: turbomachinery (impellers, rotors), swirling jets, and propeller aerodynamics.

---

## Numerical Methods

### Spatial discretisation

| | Details |
|-|---------|
| Framework | Cell-centred finite volume on structured multi-block hexahedral grids |
| Reconstruction | `MUSCL` (second-order), `MUSCL-SD` (MUSCL with shock-detector blending), `first-order` — `space-reconstruction` key |
| Flux limiters | van Leer · van Albada · minmod · MC · superbee · LIMO3 |
| Riemann solvers | HLLC · HLLC+ · HLLE · HLLE++ · HLLEM · SLAU · SLAU2 · AUSM⁺ · AUSM⁺-up · AUSM⁺-up2 · exact |
| Diffusive fluxes | 10-point stencil with face-metric tensor mapping |

### Time integration

| | Details |
|-|---------|
| Explicit | Forward Euler, RK2, RK3 (`time-scheme = euler | RK2 | RK3`) |
| Implicit | Implicit Residual Smoothing (IRS) |
| Chemistry coupling | Strang operator splitting for stiff source terms (via [OSlo](https://github.com/MarcoGrossi92/OSlo)) |

### Parallel computing

| Mode | Details |
|------|---------|
| Shared memory | OpenMP thread-level parallelism within a block |
| Distributed memory | MPI block distribution (compile-time option `USE_MPI`); greedy load balancing by cell count |
| Hybrid | OpenMP + MPI combined runs (`MPI_THREAD_FUNNELED`) |

---

## Code Dependencies

### Required libraries

| Library | Role | Source |
|---------|------|--------|
| [FLINT](https://github.com/MarcoGrossi92/FLINT) | Thermodynamic & chemical kinetics database | Bundled submodule |
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O — Tecplot, VTK, Plot3D | Bundled submodule |
| [OSlo](https://github.com/MarcoGrossi92/OSlo) | Stiff ODE solver (chemistry integration) | Bundled submodule |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser | Bundled submodule |

### Optional libraries

| Library | Role |
|---------|------|
| OpenMP | Shared-memory thread parallelism |
| MPI | Distributed-memory parallelism |
| SUNDIALS | Alternative stiff ODE integrator |
| Cantera | Extended thermochemistry back-end |
| TecIO | Binary Tecplot output |

### Build toolchain

| Tool | Minimum version |
|------|----------------|
| CMake | 3.23 (`cmake_minimum_required` in `CMakeLists.txt`) |
| Fortran compiler | GNU gfortran (CI-tested 12–14) or Intel ifx/ifort |
| C compiler | GCC or ICC (for C bindings) |
| Python | 3.9+ (for validation scripts) |

---

## Documentation Guide

| Section | What you'll find |
|---------|-----------------|
| [**Getting Started**](getting-started/index.md) | Installation, prerequisites, and first run |
| [**User Guide**](user/index.md) | Running simulations, configuring input files, boundary conditions, output |
| [**Input File Reference**](user/input.md) | `input.ini` structure, all sections, auto-generated parameter registry |
| [**Theory Guide**](theory/index.md) | Governing equations, numerical methods, turbulence and chemistry models |
| [**Verification & Validation**](vv/index.md) | Worked examples validated against exact or reference solutions |
| [**Developer Guide**](development/index.md) | Repository architecture, testing framework, contribution guidelines |
| [**About**](about/index.md) | License, acknowledgements, and contributors |

---

## License

Q2D is free and open-source software released under the **[GNU General Public License v3.0](about/license.md)** (GPL-3.0).

| Permission | |
|------------|-|
| :white_check_mark: Use freely | For any purpose, including commercial |
| :white_check_mark: Modify | Change the source code as needed |
| :white_check_mark: Distribute | Share original or modified versions |
| :white_check_mark: Patent grant | Contributors grant patent rights |
| :warning: Share-alike | Derivative works must use GPL-3.0 |
| :warning: Disclose source | Source code must be provided when distributing |

Full license text: [`LICENSE`](about/license.md)
