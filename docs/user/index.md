# User Guide

The User Guide describes how to set up, configure, and run simulations with Q2D.

<div class="grid cards" markdown>

-   :material-play-circle-outline:{ .lg .middle } __Running Q2D__

    ---

    Simulation workflow, directory structure, running the solver, output, and restart.

    [:octicons-arrow-right-24: How to run a simulation](using.md)

-   :material-file-document-outline:{ .lg .middle } __Input File__

    ---

    `input.ini` structure, sections overview, and the auto-generated parameter registry.

    [:octicons-arrow-right-24: Input file reference](input.md)

-   :material-database-outline:{ .lg .middle } __Initial Conditions__

    ---

    Supported file format and variables requirements.

    [:octicons-arrow-right-24: Initial conditions reference](initial-conditions.md)

-   :material-arrow-left-right:{ .lg .middle } __Boundary Conditions__

    ---

    File format and complete type reference.

    [:octicons-arrow-right-24: Boundary conditions reference](boundary-conditions.md)


</div>

This section is intended for:

- Users running Q2D simulations for the first time
- Experienced users looking for a quick reference

If you are new to Q2D, start from [Getting Started](../getting-started/index.md).

---

## Q2D-specific Geometric Regimes

Q2D supports **three geometric regimes**, selected at setup time via
the `[Q2D-Sources] active` INI key:

- **2D planar** (`active = none`, default): standard Cartesian flow.
- **2D axisymmetric** (`active = axisymmetric`): rotational symmetry around an axis (without swirl). Q2D adds geometric source terms `Source_Axi` for cylindric Cartesian-FV residuals.
- **Q2D area-variation** (`active = area_variation`): three-dimensional flow with prescribed out-of-plane height $h(x, y)$, projected to a 2D mesh with geometric source terms `Source_Area` (cf. [Q2D Area-Variation Derivation](../theory/area-variation-derivation.md)).

The geometric regime selection drives the dispatch of stress/strain procedure pointers `Stress_Vector_Std/_Axi`, `Stress_Tensor_Std/_Axi`, `Strain_Tensor_Std/_Axi`.
