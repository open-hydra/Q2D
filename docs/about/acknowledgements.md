# Acknowledgements

Q2D is built upon several open-source projects shared with its
companion solver MOSE in the `open-hydra` release.  The Q2D codebase
itself originated as a fork from an early version of MOSE
(pre-modernisation) and evolved independently with focus on the three
Q2D geometric regimes (2D planar, 2D axisymmetric, Q2D
area-variation), Q2D-specific numerical features and applications
to hybrid rocket motors and Borda-flow choked nozzles.

## FiNeR

**Fortran INI ParseR and generator**

- **Repository:** [github.com/szaghi/FiNeR](https://github.com/szaghi/FiNeR)
- **License:** GPL v3.0

FiNeR is a pure Fortran 2003+ OOP library for reading and writing
INI configuration files.

## ORION

**I/O Library for Fortran**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** GPL v3.0

ORION provides built-in functions to read files in different formats.
Q2D uses ORION to import thermodynamic, transport, and chemistry
tabulated data.

## FLINT

**Fortran Library for INtegrated Thermochemistry**

- **Repository:** [github.com/MarcoGrossi92/FLINT](https://github.com/MarcoGrossi92/FLINT)
- **License:** GPL v3.0

FLINT is a lightweight, modular Fortran library for thermochemistry
in reactive flow simulations.

## OSLO

**Fortran ODE Solver Library**

- **Repository:** [github.com/MarcoGrossi92/OSLO](https://github.com/MarcoGrossi92/OSLO)
- **License:** GPL v3.0

OSlo provides stiff ODE solvers (Radau IIA, SDIRK, Rosenbrock) used
by Q2D for the chemistry sub-step under Strang splitting.

## Documentation Tools

### MkDocs

**Static Site Generator**

- **Website:** [mkdocs.org](https://www.mkdocs.org)
- **License:** BSD-2-Clause
- **Contribution:** Documentation site generation.

MkDocs transforms Q2D's documentation into a beautiful, searchable
website.

### Material for MkDocs

**Modern Documentation Theme**

- **Website:** [squidfunk.github.io/mkdocs-material](https://squidfunk.github.io/mkdocs-material/)
- **License:** MIT
- **Contribution:** Professional documentation theme with advanced
  features.

Material for MkDocs provides the sleek, modern interface for Q2D's
documentation.

## Institutional Support

While Q2D is an independent open-source project, it has benefited
from:

- Academic research environments.
- Access to various CFD datasets for validation.

## License Compliance

Q2D respects all licenses of dependencies and foundations:

- **GPL v3.0** — FLINT, ORION, OSLO, FiNeR.

See the [License](license.md) page for Q2D's full license text.
