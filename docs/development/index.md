# Developer Guide

Architecture, testing, and contribution guidelines for Q2D.

## Code Structure

Q2D develops with the following architectural template:

```
src/
├── lib/
│   ├── base/         # Q2D_*_m types, parameters, advanced types, series data
│   ├── config/       # Q2D_Config_Types_m, Registry, IO_Setup
│   ├── driver/       # Wrap_Setup, Wrap_Solve, Wrap_Postprocess, Procedures_m
│   ├── numerics/
│   │   ├── space/    # Limiters, Reconstruction, Metrics, Mod_Space
│   │   ├── fluxes/   # Mod_Fluxes, Lib_Convective, Lib_Diffusive, Lib_Shock_Detector
│   │   │   ├── bc/    # 16 BC handlers
│   │   │   └── riemann/  # 12 Riemann solvers
│   │   ├── time/     # Mod_dt, time/explicit/{Mod_Explicit, Lib_Newstate, Lib_RK, Lib_IRS, Lib_Strang}
│   │   ├── multigrid/# Mod_Multigrid, Lib_Multigrid
│   │   └── sources/  # Q2D-extra: Source_Axisymmetry + Source_Areavariation
│   ├── physics/      # Lib_Chemistry, Lib_GSI, Lib_Soot, Lib_RotatingFrame
│   │   └── turbulence/  # Lib_Fluid, Mod_RANS, Lib_RANS, Lib_Spalart, Lib_SST, Lib_Wilcox2006, Lib_QCR2000, Lib_SpalartShur
│   ├── parallel/     # Mod_MPI, Mod_GhostExchange
│   ├── io/           # IO_Setup, IO_BC, IO_Solution, IO_Wall, IO_Probes, IO_LineProbes (Q2D-extra), Load_Chemistry, Load_ThermoTransport
│   └── diagnostic/   # Mod_Diagnostic
└── app/
    └── main.f90      # Q2D entry-point
```

## Testing

Q2D integrates with CTest for validation:

```bash
cd build
ctest --output-on-failure
```

Test cases:

- 1D: Sod 1979, Einfeldt 1991.
- 2D steady: Prandtl-Meyer expansion, Chang 1996.
- 2D unsteady: shock diffraction, supersonic forward step.

Reference solutions (Q2D legacy `OUTREF/`) are the ground truth for bit-per-bit (`-O0`) and tolerance-bounded (`-O3`) regression. See [`docs/vv/regression_tolerance.md`](../vv/regression_tolerance.md).

## Pre-Push Hooks

Activate the in-tree git hooks:

```bash
git config core.hooksPath .githooks
```

The hooks run:

- **`pre-commit`**: No-w invariants check (`do k =`, `Normal(3)` violations) + orphan TODO check on staged sources.
- **`pre-push`**: CTest validation.

## Contributing

Contributing to Q2D follows the conventions in the [Contributing](contributing.md) guide: the **five immutable laws** (don't lose anything; don't hallucinate; stay on scope; don't deceive; one-to-one structural parity with MOSE) plus the formal review protocol.
