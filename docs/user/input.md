# Input File

Q2D is configured through a single **INI-format** file called `input.ini`, located in the case root directory.

The file is organized into parameter blocks — each block controls a specific aspect of the simulation (flow model, numerics, output, turbulence, chemistry, etc.). Parameters not specified take their default values.

## File Structure

```ini
[SECTION-NAME]
parameter = value
```

Example:

```ini
[Q2D-Parameters]
res-threshold = 1e-8
iter-threshold = 500000

[Q2D-Physics]
equations = navier-stokes
turbulence = SST

[Q2D-Numerics]
cfl = 0.9
time-scheme = RK3
space-reconstruction = MUSCL
```

!!! warning "Default values"
    If a parameter is not specified, the **default value** is used silently.

!!! warning "Parameter names are case-sensitive"
    Incorrect names or values are ignored and the default is used instead. Check the [Parameter Registry](registry.md) for the exact spelling.

---

## Section namespacing

All Q2D configuration sections share the prefix `[Q2D-...]`.  The
prefix comes from the solver **codename**, a compile-time constant
(`Q2D_Parameters_m`): nothing needs to be set in `input.ini` to
activate the namespace.

```ini
[Q2D-Parameters]
res-threshold = 1.0e-8
```

The same convention applies to optional feature sections such as
line probes (see [Line probes](../development/structure.md#line-probes)).

---

## Sections

| Section | Description | Reference |
|---------|-------------|-----------|
| `[Q2D-Parameters]` | Simulation control: convergence thresholds, restart | [→](registry.md#q2d-parameters) |
| `[Q2D-Physics]` | Equations, turbulence, chemistry, soot, rotating frame | [→](registry.md#q2d-physics) |
| `[Q2D-IO]` | Output formats, frequencies, and variables | [→](registry.md#q2d-io) |
| `[Q2D-Probes]` | Point probes for time-history recording | [→](registry.md#q2d-probes) |
| `[Q2D-LineProbes]` | Line probes: 1D profiles along fixed-index grid lines | [→](registry.md#q2d-lineprobes) |
| `[Q2D-Numerics]` | Spatial/temporal discretization, Riemann solver, CFL | [→](registry.md#q2d-numerics) |
| `[Q2D-Multigrid]` | Multigrid acceleration settings | [→](registry.md#q2d-multigrid) |
| `[Q2D-Chemistry]` | ODE solver and tolerances for finite-rate chemistry | [→](registry.md#q2d-chemistry) |
| `[Q2D-Turbulence]` | Turbulent Prandtl/Schmidt numbers | [→](registry.md#q2d-turbulence) |
| `[Q2D-Rotating-Frame]` | Rotating reference frame parameters | [→](registry.md#q2d-rotating-frame) |
| `[Q2D-Sources]` | Geometric source selection (axisymmetric / area-variation) | [→](registry.md#q2d-sources) |

---

## Parameter Registry

The complete list of all parameters, defaults, and allowed values is in the **[Parameter Registry](registry.md)**. It is generated automatically from the source code at every release — do not edit it manually.
