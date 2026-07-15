# supersonic-forward-step reference data

## trace_baseline.dat

Committed regression baseline for `check.py`: density along a fixed cell row
of the above-step block at the deterministic time-accurate endpoint t=4.0 of
the Woodward-Colella (1984) M=3 forward-step benchmark.

- **Quantity**: cell-centered `rho(1)` along zone 2 (above-step block,
  240x64 cells), cell row `j=23` (the row nearest y=0.5, cutting through the
  bow-shock / Mach-stem region). Columns: `x_cell_center  rho`.
- **Generated**: 2026-07-15, from a fresh `./Q2D.sh solve` run with
  `OMP_NUM_THREADS=1` (single-thread, deterministic).
- **Solver build**: `bin/Q2D`, gfortran TESTING build, repository `main`
  branch at commit `8c60e91` ("FiNeR v2.0.6").
- **Case setup**: `input.ini` as committed (time-accurate RK3, CFL 0.8,
  HLLE + MUSCL/vanleer, `time-threshold = 4.0`); meshes under `MESH/`.

`check.py` gates every future run against this trace at max relative
difference 5 percent per point (single-thread reruns of the same binary are
bit-identical; the 5 percent headroom absorbs compiler/optimization-level
differences away from a shock crossing) and gates the bow-shock x-position
(steepest density gradient along the row) within 3 cells (0.0375) of the
baseline's.

Regenerate (only after an intentional physics/numerics change, with
justification in the commit message): run `./Q2D.sh solve`, then extract the
same row with the snippet documented in `check.py` and update this README's
build/commit provenance.

## OpenFOAM/

EnSight export of an OpenFOAM rhoCentralFoam run of the same benchmark
(different solver/numerics). Kept for manual cross-solver comparison only;
not used by `check.py`.
