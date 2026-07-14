# Theoretical Guide

This section provides the theoretical foundations for the physical
models and numerical methods implemented in Q2D.  The material spans
the governing equations of compressible multi-species flow, the
finite-volume spatial discretization, Riemann solvers, time
integration, turbulence closures, the rotating-frame formulation,
the geometric source terms of the axisymmetric and area-variation
regimes, the thermochemistry framework, and the equilibrium solver.

<div class="grid cards" markdown>

-   :material-math-integral:{ .lg .middle } __Governing Equations__

    ---

    Compressible Navier–Stokes for multi-species ideal-gas mixtures: conservative form, state vectors, convective and viscous fluxes, source terms

    [:octicons-arrow-right-24: Governing equations](governing-equations.md)

-   :material-grid:{ .lg .middle } __Spatial Discretization__

    ---

    Cell-centred finite volume method, MUSCL reconstruction, flux limiters, Jameson shock detection

    [:octicons-arrow-right-24: Spatial discretization](numerics.md)

-   :material-waves:{ .lg .middle } __Riemann Solvers__

    ---

    Approximate and exact solvers: AUSM family, HLL family, Godunov, Lax–Friedrichs, SLAU

    [:octicons-arrow-right-24: Riemann solvers](riemann-solvers.md)

-   :material-timer-outline:{ .lg .middle } __Time Integration__

    ---

    SSP Runge–Kutta, CFL/VNN stability, Strang splitting, implicit residual smoothing, geometric multigrid

    [:octicons-arrow-right-24: Time integration](time-integration.md)

-   :material-weather-windy:{ .lg .middle } __Turbulence Modelling__

    ---

    Spalart–Allmaras (+ SA-RC, SAcomp, SAR), Menter SST $k$–$\omega$, Wilcox 2006 $k$–$\omega$

    [:octicons-arrow-right-24: Turbulence models](turbulence.md)

-   :material-rotate-3d:{ .lg .middle } __Rotating Frame__

    ---

    Coriolis and centrifugal forces, relative velocity formulation, turbomachinery applications, Spalart–Shur correction

    [:octicons-arrow-right-24: Rotating frame](rotating-frame.md)

-   :material-axis-arrow:{ .lg .middle } __Geometric Source Terms__

    ---

    Overview hub for the axisymmetric ($1/r$ cylindric residual) and area-variation ($\mathbf{V}\!\cdot\!\mathbf{G}$ depth-averaging) regimes, with full formal derivations in [2D Axisymmetric](axisymmetric-derivation.md), [q2D Area-Variation](area-variation-derivation.md), and [Q2D RANS Closures](q2d-rans-closures.md)

    [:octicons-arrow-right-24: Geometric source terms](geometric-sources.md)

-   :material-thermometer:{ .lg .middle } __Thermodynamic and Transport Properties__

    ---

    Equation of state, mixture rules, property evaluation, and Wilke transport models

    [:octicons-arrow-right-24: Thermodynamic models](thermo.md)

-   :material-chart-bell-curve:{ .lg .middle } __Finite-Rate Kinetics__

    ---

    Arrhenius reactions, three-body collisions, Lindemann falloff, Troe formulation, surface kinetics

    [:octicons-arrow-right-24: Kinetics models](kinetics.md)

-   :material-scale-balance:{ .lg .middle } __Chemical Equilibrium__

    ---

    NASA CEA algorithm for constant-UV equilibrium problems

    [:octicons-arrow-right-24: Equilibrium solver](equilibrium.md)

</div>

---

## Overview

The Q2D solver advances the compressible Navier–Stokes equations in
conservative form on structured multi-block grids in **two
dimensions**, with three supported geometric regimes:

| Regime | Geometry | Geometric source |
|--------|----------|------------------|
| **2D planar** | Standard Cartesian $(x, y)$ | None |
| **2D axisymmetric** | Cylindrical $(x, r)$, $u_\theta \equiv 0$ | $\sim 1/r$ cylindric residual |
| **q2D area-variation** | Depth-averaged channel with height $h(x, y)$ | $\sim \mathbf{V}\!\cdot\!\mathbf{G}$ depth-averaging |

The numerical pipeline can be summarised as follows:

| Stage | Method | Page |
|-------|--------|------|
| **Governing system** | Multi-species Euler / Navier–Stokes (two-dimensional) | [Governing Equations](governing-equations.md) |
| **Spatial discretization** | Cell-centred FVM + MUSCL-limited reconstruction | [Spatial Discretization](numerics.md) |
| **Interface fluxes** | Approximate Riemann solvers (AUSM, HLL, Godunov, LF, SLAU) | [Riemann Solvers](riemann-solvers.md) |
| **Time marching** | 3-stage SSP RK3, IRS, geometric multigrid | [Time Integration](time-integration.md) |
| **Turbulence closure** | SA, SST $k$–$\omega$, Wilcox 2006 $k$–$\omega$ | [Turbulence Modelling](turbulence.md) |
| **Geometric sources** | Cylindrical axisymmetric or depth-averaged area-variation | [Geometric Source Terms](geometric-sources.md) |
| **Thermochemistry** | Ideal-gas mixture, NASA polynomials | [Thermodynamic Properties](thermo.md) |
| **Chemistry** | Finite-rate (Arrhenius/Troe) or equilibrium (CEA) | [Kinetics](kinetics.md), [Equilibrium](equilibrium.md) |

