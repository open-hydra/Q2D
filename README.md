<p align="center">
  <h1 align="center">Q2D</h1>
  <p align="center"><b>A Fortran solver for 2D planar, 2D axisymmetric, and quasi-2D area-variation compressible flows</b></p>
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.21394015"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.21394015.svg" alt="DOI"></a>
  <a href="https://github.com/open-hydra/Q2D/actions/workflows/ctest.yml"><img src="https://github.com/open-hydra/Q2D/actions/workflows/ctest.yml/badge.svg" alt="CTest"></a>
  <a href="https://github.com/open-hydra/Q2D/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
  <a href="https://open-hydra.github.io/Q2D/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran-734f96.svg" alt="Language: Fortran">
</p>

---

Q2D is an open-source compressible flow solver for the Euler and Navier–Stokes equations on structured grids, specialised for **three reduced-dimensional geometric regimes**: 2D planar, 2D axisymmetric, and quasi-2D area-variation. It is the companion of [MOSE](https://github.com/open-hydra/MOSE) (the full-3D solver) in the `open-hydra` release, sharing the same Fortran library stack and architectural philosophy.

## Features

- **Compressible Euler & Navier–Stokes** in three geometric regimes (2D planar, 2D axisymmetric, Q2D area-variation), steady-state and time-accurate.
- **High-order numerics** — MUSCL reconstruction with multiple flux limiters and Riemann solvers.
- **Turbulence modeling** — RANS closures: Spalart–Allmaras (with rotation/curvature corrections via Spalart–Shur), Menter SST, Wilcox k-ω 2006, QCR2000.
- **Finite-rate chemistry & soot** — multi-species reacting flows with finite-rate kinetics (FLINT) or chemical equilibrium (NASA CEA). Soot transport included.
- **Borda-flow inflow** — Fievisohn 2016/2017 generalised Riemann invariant for choked sudden-expansion injection (PDE/RDE chamber inflow).
- **Parallel execution** — shared-memory parallelism via OpenMP with a runtime-selectable dispatcher (`!$omp do` vs `!$omp task/taskloop`). MPI support for multi-block distributed-memory runs.
- **Flexible I/O** — Tecplot (ASCII and binary via TecIO), VTK, Plot3D via ORION. Line probes and point probes for time-history recording. Restart capability.

## Q2D 2D modelling features

Q2D is structurally a mirror of MOSE, specialised for reduced-dimensional flows. Highlights:

- **No-w principle** — physical removal of all `w` (out-of-plane velocity) residuals; tensors retain the 3×3 layout.
- **Area-variation and axisymmetric source closure**.
- **Line probes** — sampling along fixed-index i / j lines.
- **OMP do/task dispatcher** — runtime-selectable OpenMP parallelisation strategy.
- **Bootstrap-guard input.ini** — `[Q2D-*]` section prefix matched at runtime via the `codename` mechanism.

## Quick Start

### Prerequisites

| Requirement | Details |
|---|---|
| **CMake** | ≥ 3.23 |
| **Fortran compiler** | GNU (`gfortran`) or Intel/oneAPI (`ifort` / `ifx`) |
| **C/C++ compiler** | Required only for optional components (Cantera, TecIO, SUNDIALS) |
| **Python** | ≥ 3.6 with `numpy` (required by `check.py` regression scripts); `scipy` + `matplotlib` only for the standalone `verify.py' |

### Build

```bash
git clone --recurse-submodules https://github.com/open-hydra/Q2D.git
cd Q2D

# Build with GNU compilers and OpenMP
./install.sh build --compiler=gnu --use-openmp

# — or with Intel compilers and full feature set —
./install.sh build --compiler=intel --use-openmp --use-mpi \
                    --use-tecio --use-sundials --use-cantera
```

The executable is placed in `bin/Q2D`.

See the [Installation Guide](https://open-hydra.github.io/Q2D/getting-started/installation/) for all build options, CMake presets, and troubleshooting.

### Run the Sod Shock Tube

```bash
cd test/1D/Sod79
./Q2D.sh test     # mesh-refinement on x1/x2/x4/x8
python3 -B verify.py   # standalone visualization + ExactPack analytical reference
```

For the CTest-driven regression workflow:

```bash
ctest --test-dir build --output-on-failure -R Sod79   # one case
ctest --test-dir build -LE heavy --output-on-failure  # full fast tier (CI set)
```

See the [Quick Start](https://open-hydra.github.io/Q2D/getting-started/quick-start/) for a full walkthrough.

## Dependencies

Q2D is built on top of the same companion libraries as MOSE, included as Git submodules:

| Library | Role |
|---|---|
| [FLINT](https://github.com/MarcoGrossi92/FLINT) | Thermodynamic and chemical kinetics database |
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O (Tecplot, VTK, Plot3D) |
| [OSLO](https://github.com/MarcoGrossi92/OSLO) | ODE solver library (stiff chemistry integration) |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser |
| [ExactPack](https://github.com/lanl/ExactPack) | LANL analytical reference solutions for V&V |

Optional external libraries: **OpenMP**, **MPI**, **SUNDIALS**, **Cantera**, **TecIO**.

## Project Structure

```
Q2D/
├── src/
│   ├── app/           # Main application
│   └── lib/           # Solver library sources
│       ├── numerics/  # Space, time, fluxes, sources, BCs
│       ├── physics/   # Turbulence, chemistry, soot, GSI
│       ├── io/        # IO_BC, IO_Solution, IO_Probes, IO_LineProbes
│       ├── driver/    # Wrap_Setup, Wrap_Solve, Wrap_Postprocess
│       ├── diagnostic # Residuals, headers, checks
│       └── config/    # Parameters, Registry, Backend_INI
├── lib/               # Git submodule dependencies
├── test/              # Verification & validation cases
│   ├── 1D/            #   1-D shock tubes (Sod79, Ein91, Tor99, Noh87)
│   └── 2D/
│       ├── euler/     #   2-D Euler (chang-nozzle, hypersonic-cylinder,
│       │              #              ramp-channel, oblique-shock, prandtl-meyer,
│       │              #              rocket-nozzle ×3, supersonic-forward-step)
│       └── viscous/   #   2-D viscous (flat-plate-laminar, flat-plate-turbulent ×3, swbli)
├── docs/              # MkDocs documentation source
├── cmake/             # CMake modules
├── install.sh         # Build helper script
└── CMakeLists.txt
```

## Documentation

Full documentation is available at **[open-hydra.github.io/Q2D](https://open-hydra.github.io/Q2D/)**, covering:

- Installation & quick start
- User guide & input file reference
- Verification & validation results (1D mesh-refinement convergence + 2D bit-exact MOSE comparison)
- Theory guide (governing equations, three geometric regimes, area-variation closures)

## Citation

If you use Q2D in academic work, please cite the metadata in [`CITATION.cff`](CITATION.cff).

## License

Q2D is free and open-source software released under the [GNU General Public License v3.0](LICENSE).
