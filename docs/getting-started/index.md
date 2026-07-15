# Getting Started

Welcome to the Q2D getting-started guide! This section will help you install Q2D and run your first examples.

!!! info "Q2D within Hydra"
    Q2D is the quasi-2D compressible flow solver of **Hydra**. The pre-processor **ATLAS** (distributed separately) is used to generate the input files that Q2D requires: initial conditions, boundary conditions, and thermodynamic/transport/chemical data. See the [overview page](../overview.md#hydra-cfd-suite) for an overview of the ecosystem.

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } __Installation__

    ---

    Build the API and the executable from source.

    [:octicons-arrow-right-24: Install Q2D](installation.md)

-   :material-rocket-launch:{ .lg .middle } __Quick Start__

    ---

    Get up and running in minutes (Prandtl–Meyer expansion).

    [:octicons-arrow-right-24: Quick start tutorial](quick-start.md)

</div>

## Scope of This Section

This getting-started guide covers:

1. **Installation**
   - Step-by-step build instructions.
   - Build options and customisation.
   - Activation of the in-tree git hooks
     (No-$w$ pre-commit + CTest pre-push).

2. **Quick Start**
   - Run your first case (Prandtl–Meyer expansion in 2D planar).
   - Verify the result against the analytic Prandtl–Meyer
     prediction.

## Prerequisites at a Glance

Q2D requires:

- A modern Fortran 2008+ compiler (GNU `gfortran`, CI-tested with
  versions 12–14, or Intel `ifx`).
- CMake ≥ 3.23.
- Companion libraries [FLINT](https://github.com/MarcoGrossi92/FLINT),
  [ORION](https://github.com/MarcoGrossi92/ORION),
  [OSlo](https://github.com/MarcoGrossi92/OSlo),
  [FiNeR](https://github.com/szaghi/FiNeR) (managed as git
  submodules in `lib/`).

Optional: OpenMP runtime, MPI implementation (e.g., Intel MPI),
SUNDIALS, Cantera, TecIO.

For the full installation procedure see
[Installation](installation.md).
