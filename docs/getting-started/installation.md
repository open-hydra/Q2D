# Installation

This document describes how to obtain and build **Q2D**.  The
instructions describe the current `install.sh` script, CMake
configuration and Git submodule layout.

!!! note
    Q2D has a dual nature: it is both a library and an executable.
    The installation process described here produces both the
    static library `libQ2D.a` and the main executable `bin/Q2D`.
    If you are only interested in using Q2D as a library, you can
    link against `libQ2D.a` without caring about the executable.

## Prerequisites

Before attempting to build Q2D make sure your system provides the
following external tools and compilers:

- **CMake** — 3.23 or newer.
- **Fortran compiler** — either the GNU toolchain (`gfortran`) or
  Intel/oneAPI (`ifort`/`ifx`) are supported.
- **C++ compiler** — required only if you wish to enable optional
  components that depend on C++ (Cantera, TecIO).
- **C compiler** — needed when building SUNDIALS (optional).
- **OpenMP** — needed for optional parallelisation support.

### Git submodules

Q2D depends on the same companion library stack as MOSE — included as
Git submodules:

| Path                               | Repository URL                                        | Purpose                                  |
|------------------------------------|-------------------------------------------------------|------------------------------------------|
| `lib/third_party/ExactPack`        | `https://github.com/lanl/ExactPack.git`              | Exact solutions database for reference   |
| `lib/ORION`                        | `https://github.com/MarcoGrossi92/ORION.git`         | I/O routines (TecIO, VTK, Plot3D, ...)   |
| `lib/third_party/FiNeR`            | `https://github.com/szaghi/FiNeR.git`                | INI file parser                           |
| `build/lib/OSlo` (installed by FLINT) | `https://github.com/MarcoGrossi92/OSlo.git`       | ODE solvers library                      |
| `lib/FLINT`                        | `https://github.com/MarcoGrossi92/FLINT.git`         | Thermodynamic database utilities          |

## Build methods

First clone the repository with submodules:

```bash
git clone https://github.com/open-hydra/Q2D.git
cd Q2D
# initialise submodules
git submodule update --init --recursive
```

To fully install Q2D, you may either use the bundled install script
or invoke CMake manually.  The script is convenient and is the
preferred route for most users.

### Build with `install.sh` (recommended)

The script exposes three commands: `build`, `compile`, and `update`.
It also maintains a `CMakePresets.json` file that records the
configuration used for the most recent `build` invocation.

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**Global options**

* `-v`, `--verbose` — enable verbose logging.
* `-h`, `--help` — print the usage summary and exit.

**`build` command**

Performs a clean configure + build cycle.  Example usage:

```bash
# minimal GNU build with OpenMP enabled
./install.sh build --compilers=gnu --use-openmp

# full configuration with Intel compilers and all optional features
./install.sh build --compilers=intel \
                  --use-openmp --use-mpi \
                  --use-tecio --use-sundials --use-cantera
```

Options accepted by `build`:

* `--compilers=<gnu|intel>` — select the compiler family.  If
  omitted, no `FC`/`CC`/`CXX` exports are made and CMake auto-detects
  the system compilers.
* `--use-openmp` — enable OpenMP parallelisation.
* `--use-mpi` — enable MPI parallelisation.
* `--use-tecio` — enable TecIO support (requires a C++ compiler).
* `--use-sundials` — build and link the bundled SUNDIALS library
  (requires a C compiler).
* `--use-cantera` — build and link Cantera (requires a C++ compiler;
  if enabled, SUNDIALS should also be enabled).
* `--include-orion=PATH` — use an external ORION tree instead of the
  submodule.
* `--include-flint=PATH` — same for FLINT.
* `--include-oslo=PATH` — same for OSlo.
* `--include-finer=PATH` — same for FiNeR.

The script sets up environment variables for the chosen compilers
and then invokes CMake.  After a successful build, a
`CMakePresets.json` file is written in the source root so that
subsequent compilations can reuse the configuration.

**`compile` command**

Re-runs CMake using the previously generated preset and rebuilds
the project without clearing the build directory.  Useful during
development when only the source has changed:

```bash
./install.sh compile
```

**`update` command**

Synchronises the Git submodules.  By default it checks out the
commit recorded in `.gitmodules`; passing `--remote` will fetch the
latest commit from each remote branch.

```bash
./install.sh update            # sync to recorded commit
./install.sh update --remote   # update to newest remote commit
```

### Build with CMake

If you prefer fine-grained control, perform the configuration
yourself.  This is essentially what `install.sh` does under the
hood.

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \  # or ifx
    -DUSE_OPENMP=ON \                    # or OFF
    -DUSE_TECIO=ON \                     # optional (needs C++ compiler)
    -DUSE_SUNDIALS=ON \                  # optional (needs C compiler)
    -DUSE_CANTERA=ON \                   # optional (needs C++ compiler)
    -DORION_PATH=/path/to/ORION \        # optional
    -DFLINT_PATH=/path/to/FLINT \        # optional
    -DOSLO_PATH=/path/to/OSlo \          # optional
    -DFINER_PATH=/path/to/FiNeR          # optional
cmake --build . --parallel
```

CMake build artefacts (objects, modules, the static `Q2DL` library)
stay inside `build/`, but the main executable is emitted to the
**repository-root** `bin/Q2D` (`src/app/CMakeLists.txt` sets
`CMAKE_RUNTIME_OUTPUT_DIRECTORY` to `<source-root>/bin`) — this is
the path the test drivers (`Q2D.sh`) expect.

## CMake presets

The file `CMakePresets.json` produced by the install script records
all of the cache variables that were used during configuration.
You can build the project later simply with

```bash
cmake --preset default
cmake --build build
```

or using the `compile` command of the install script as described
above.

## Optional components

### TecIO

Q2D can be built with support for TecIO, a library for writing
Tecplot files.  This is an optional feature shipped with ORION.  If
enabled, Q2D will be able to read and write output in Tecplot binary
format, which can be useful for post-processing and visualisation
of large datasets.  Enabling TecIO requires a working C++ compiler.

### SUNDIALS

SUNDIALS is a suite of solvers for ODEs.  The library is included
as a submodule in OSLo and will be built automatically if you
select the `--use-sundials` option during installation.  Enabling
SUNDIALS requires a working C compiler.

### Cantera

Q2D's core thermodynamic library does **not** depend directly on
Cantera, but only on FLINT.  The same warnings and observations that
apply to MOSE are applicable to Q2D: Cantera should be used for
small projects when advanced chemistry not already included in FLINT
is required, but for large-scale simulations with parallel support
the FLINT Fortran interface should be preferred.  See the
[FLINT documentation](https://marcogrossi92.github.io/FLINT/) for
details.

!!! note
    If Cantera is enabled, a working C++ compiler is mandatory and
    `--use-sundials` should also be selected because Cantera itself
    relies on the SUNDIALS library.

## Library linking (advanced)

To use Q2D from an external Fortran program you can compile as
follows:

```bash
gfortran -I/path/to/Q2D/include \
         -L/path/to/Q2D/lib \
         -lQ2D \
         your_program.f90 -o your_program
```

or, from a CMake project:

```cmake
find_package(Q2D REQUIRED)
add_executable(myapp main.f90)
target_link_libraries(myapp Q2D::Q2D)
```

Installation prefix and other details may be customised via standard
CMake variables such as `CMAKE_INSTALL_PREFIX` and by running
`cmake --install`.

## Activate in-tree git hooks (recommended for contributors)

After cloning, contributors should activate the in-tree hooks that
enforce the No-$w$ principle and run the validation tests before
each push:

```bash
git config core.hooksPath .githooks
```

The pre-commit hook runs grep-based invariant checks (No-$w$
principle, orphan TODO markers) on staged sources; the pre-push hook
runs the CTest validation suite against the `build/` directory.  See
[Testing](../development/testing.md#pre-push-hooks).

## Next steps

* **[Quick Start Tutorial](quick-start.md)** — build and verify the
  installation.
