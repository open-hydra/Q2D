# Running Q2D

This page describes the complete workflow for running a Q2D
simulation: preparing the case directory, defining the needed files,
launching the solver, and inspecting the output.

All configuration is driven by:

- `input.ini` — the main INI file that defines the simulation
  parameters, including the numerical settings and physical models.
- `INPUT/` — a directory containing the supporting data files for the
  simulation, such as thermodynamic properties, initial conditions,
  and boundary conditions.

!!! tip
    All files within `INPUT/` can be built using the Q2D pre-processing
    tools provided by ATLAS (shared with MOSE).

!!! warning
    ATLAS is a powerful tool, but it is not required to run Q2D.  All
    files can be created manually or with custom scripts as long as
    they follow the expected format.

For the detailed format of thermodynamic, initial and boundary
conditions see the dedicated reference pages:

- [Thermodynamic Data (FLINT library)](../theory/thermo.md)
- [Initial Conditions](initial-conditions.md)
- [Boundary Conditions](boundary-conditions.md)

---

## Simulation Workflow

A typical Q2D run follows the following phases (using ATLAS):

1. **Setup** — define the mesh, thermodynamic data, initial and
   boundary conditions.
2. **Preprocess** — run ATLAS pre-processing tools to generate the
   necessary input files in `INPUT/`.
3. **Solve** — advance the solution in time (explicit time
   integration).
4. **Postprocess** — write solution files, wall data, probes,
   line probes, and diagnostics.

All configuration is driven by the INI file (`input.ini`).

---

## Geometric regime selection

Q2D supports three geometric regimes selected at setup time through
the `active` key of the `[Q2D-Sources]` section:

```ini
[Q2D-Sources]
active = none   ; none | axisymmetric | area_variation
```

| `active` | Resulting regime |
|:-:|---|
| `none` | **2D planar** (default) — Cartesian flow |
| `axisymmetric` | **2D axisymmetric** (no swirl) — cylindric residual −1/r·(…) |
| `area_variation` | **Q2D area-variation** — depth-averaged residual ∇h/h·(…) + Definition A V·G; `INPUT/height.tec` required |

The two geometric source models are mutually exclusive; the solver
stops with an error if both are requested.

The selection drives the dispatch of the geometric source procedure
pointers and the `Stress_Vector_*Std` / `_Axi` viscous closures; see
[Theory — Geometric Source Terms](../theory/geometric-sources.md).

## Rotating frame

A rigid-body rotating reference frame is activated through the
`rotational-frame` key of `[Q2D-Physics]` plus the
`[Q2D-Rotating-Frame]` section:

```ini
[Q2D-Physics]
rotational-frame = rigid-body   ; rigid-body | none (default none)

[Q2D-Rotating-Frame]
omega  = 10.0            ; rotation rate [rad/s]
axis   = 0.0 0.0 1.0     ; rotation axis (default z)
origin = 0.0 0.0 0.0     ; a point on the axis [m]
```

See `test/2D/euler/rotating-channel` for a complete working example
(gated against the analytic Coriolis pressure balance — see the
[rotating channel V&V page](../vv/2D-rotating-channel.md)) and
[Theory — Rotating Frame](../theory/rotating-frame.md) for the
source-term formulation.

## Chemistry selection

The chemistry mode is selected through the `chemistry` key of
`[Q2D-Physics]`; the stiff-ODE integrator used by finite-rate
kinetics is configured in `[Q2D-Chemistry]`:

```ini
[Q2D-Physics]
chemistry = finite-rate   ; frozen (default) | equilibrium | finite-rate

[Q2D-Chemistry]
ode-solver    = sdirk4b   ; H-radau5 (default) | sdirk4b | ros4
ode-max-steps = 100000    ; default
; per-equation tolerances (all default 1e-5):
;   ode-relative-tol-species, ode-relative-tol-temperature,
;   ode-absolute-tol-species, ode-absolute-tol-temperature
```

Finite-rate runs additionally need the mechanism files in `INPUT/`
(`chemistry-info.txt`, `chemistry-Arrhenius.dat`,
`chemistry-Troe.dat`) alongside the multi-species `phase.txt` /
`thermo.dat`.  See the three
`test/2D/euler/rocket-nozzle/{frozen,equilibrium,finite-rate}`
working examples (gated in ctest — see the
[rocket-nozzle V&V page](../vv/2D-nozzle.md)) and
[Theory — Reactive Kinetics](../theory/kinetics.md).

---

## Case Directory Structure

Every Q2D case follows a standard layout:

```
my_case/
├── input.ini              ← case configuration
├── Q2D.sh                 ← convenience run script (optional)
├── INPUT/
│   ├── phase.txt          ← species and mixture definition
│   ├── thermo.dat         ← thermodynamic property tables
│   ├── ...                ← other files (chemistry, transport, etc.)
│   ├── ic.tec             ← initial condition grid + flow field
│   ├── height.tec         ← (Q2D area-variation only) prescribed h(x,y)
│   └── bc.txt             ← face-by-face boundary condition data
└── OUTPUT/                ← created at runtime
    ├── field.tec          ← solution fields (ρ, u, v, w=0, p, T, …)
    ├── field_0001.tec     ← solution snapshots (when sol-overwrite = false)
    ├── ...
    ├── wall.tec           ← wall quantities (viscous cases only)
    ├── residual-history.dat  ← residual norm vs. iteration
    ├── <name>.txt         ← time-history at point probe <name>
    └── <name>.tec         ← profile history at line probe <name>
```

!!! warning
    Q2D reads the grid geometry directly from the initial condition
    file (the `x`, `y`, `z` node coordinates).  There is no separate
    mesh file.

---

## Running the Solver

**Using `Q2D.sh`**

Each test case includes a convenience script:

```bash
./Q2D.sh solve                  # single-threaded
./Q2D.sh -p 4 solve             # OpenMP with 4 threads
./Q2D.sh -b solve               # background execution
./Q2D.sh -b -p 8 solve          # background + 8 threads
```

**Directly**

```bash
export OMP_NUM_THREADS=4
ulimit -s unlimited
./bin/Q2D
```

!!! warning
    Q2D may require a large stack size.  Always run
    `ulimit -s unlimited` (or set `export KMP_STACKSIZE=100M` for
    Intel compilers) before launching the solver.

!!! info "OpenMP strategy"

    The OpenMP loop dispatcher in `Mod_Fluxes` can be selected at
    run time via `[Q2D-Numerics] OMP-parallel`:

    - `OMP-parallel = do` (default) — `!$omp do collapse(2)` pattern.
      Diagnostic-friendly.
    - `OMP-parallel = task` — `!$omp taskloop collapse(2)` pattern with
      explicit `taskwait`.  Beneficial for imbalanced block
      decomposition.

    See [Theory — Spatial Discretization](../theory/numerics.md#openmp-dotask-dispatcher).

### Running with MPI

Q2D can also be built with distributed-memory parallelism
(`src/lib/parallel/Mod_MPI.f90`).  MPI support is a compile-time
option: configure with `-DUSE_MPI=ON` (or `./install.sh build
--use-mpi`) and launch through the MPI runtime:

```bash
mpirun -np 4 ./bin/Q2D
```

Parallelisation is *block-distributed*: the grid blocks are
partitioned across the ranks with a greedy load-balancing strategy
(each block is assigned to the least-loaded rank, weighted by its
cell count).  The number of ranks should therefore not exceed the
number of blocks.  When compiled without `USE_MPI`, all MPI wrappers
fall back to serial no-ops.

---

## Multigrid

Geometric multigrid acceleration is activated through the
`[Q2D-Multigrid]` section:

```ini
[Q2D-Multigrid]
levels      = 2       ; number of grid levels (1 = single grid, default)
level1-iter = 40000   ; iterations on the fine level
level2-iter = 200     ; iterations on the coarse level
```

Note that `levels` is read directly by the INI backend
(`Backend_INI.f90`) rather than through the parameter-registry
mechanism, so it does not appear in the auto-generated
[Parameter Registry](registry.md) table for `[Q2D-Multigrid]`.

Each coarse level requires its own boundary-condition file in
`INPUT/`: `bc2.txt` for level 2 (and so on), obtained by restricting
the fine-grid `bc.txt`.  The helper script `tools/gen_coarse_bc.py`
generates these files automatically.  See
`test/2D/euler/planar-nozzle-mg2` for a complete working example.

---

## Logging and Diagnostics

Q2D prints diagnostic information to the console during loading and
iteration:

- loading status of the input files (thermodynamics, transport,
  chemistry, initial and boundary conditions);
- a summary of the case set-up (blocks, cells, boundary faces and
  their types, physical model, numerical scheme choices);
- iteration count and density residual during the run.

!!! tip "Loading checks"

    If any issue occurs during loading (missing files, format errors),
    Q2D prints an `[ERROR] ...` diagnostic and exits — a clean
    `Loading` block with all `OK` entries is a good sign that the case
    is set up correctly.

Console output of the `prandtl-meyer` validation case
(`test/2D/euler/prandtl-meyer`, captured from a real run):

```
 =========================================================================================
 Loading
 =========================================================================================
   Thermodynamics                 OK
   Transport                      OK
   Chemistry                      OK
   Input file                     OK
   Initial conditions             OK
   Boundary conditions            OK
 =========================================================================================

 =========================================================================================
 Set-up
 =========================================================================================
 Domain
   Blocks                         1
   Cells                          2500
   Boundary faces                 200

 Boundary conditions
   Inflow/outflow                 50
   Symmetry                       50
   Extrapolation                  100

 Input/Output
   Initial conditions file        INPUT/ic.tec
   Solution format                tecplot ascii
   Probes number                  1
 =========================================================================================

 =========================================================================================
 Physical model
 =========================================================================================
 Gas model
   Equations                      Euler
   Equation of state              Ideal
   Thermodynamics                 Thermally perfect gas
   Species                        1
 =========================================================================================

 =========================================================================================
 Numerical scheme
 =========================================================================================
 Time
   Scheme                         Third-order Runge-Kutta
   Integration variables          Conservative
   Implicit residual smoothing    Beta set to 0.500000E+00

 Fluxes
   Riemann solver                 HLLE
 =========================================================================================

Q2D | Iter =      100 | Global iter =      100 | Density residual = 0.125027E-03
Q2D | Iter =      200 | Global iter =      200 | Density residual = 0.231340E-06
Q2D | Iter =      300 | Global iter =      300 | Density residual = 0.221740E-08
Q2D | Iter =      384 | Global iter =      384 | Density residual = 0.998459E-10

   Time of operation was   3.1026850000000002E-002 min
```

The `Density residual` is the L2 norm of the density time derivative;
steady runs stop when it reaches `res-threshold`
(`[Q2D-Parameters]`), when `iter-threshold` iterations are exceeded,
or when `time-threshold` is reached. With multigrid the iteration
line carries the grid level (`Q2D Grid Level 2 | Iter = ...`).

---

## Output

After a run completes, the `OUTPUT/` directory contains:

| File | Description |
|------|-------------|
| `field.tec` | Solution fields ($\rho$, $u$, $v$, $w \equiv 0$, $p$, $T$, …) |
| `residual-history.dat` | Iteration-by-iteration residual norms |
| `wall.tec` | Wall quantities (skin friction, heat flux) — viscous cases only |
| `<name>.txt` | Time-history data at point probe `<name>` |
| `<name>.tec` | Profile history at line probe `<name>` |

### Output Formats

Controlled by `sol-format` in `[Q2D-IO]`:

| Value | Format |
|-------|--------|
| `tecplot ascii` | Tecplot ASCII (`.tec`) — default |
| `tecplot binary` | Tecplot binary (`.szplt`) — requires TecIO |
| `vtk ascii` | VTK ASCII (`.vtm`) |
| `vtk raw` | VTK binary (`.vtm`) |

### Output Variables

By default, Q2D writes primitive variables ($\rho$, $u$, $v$, $w$,
$p$, $T$).  The $w$ column is always zero but is preserved in the
output for tooling compatibility.  Additional
variable groups can be enabled in `[Q2D-IO]` via `sol-variables`:

- `thermo` — specific heat ratio $\gamma$, gas constant $R$.
- `transport` — laminar viscosity $\mu_l$, thermal conductivity
  $k_l$, turbulent viscosity $\mu_t$.

### Probes

**Point probes** — point-measurement time histories are configured
with probe sections:

```ini
[Q2D-Probes]
probe1 = exit

[exit]
variables      = p T u
dtime          = 1e-4
index-position = 1 50 10   ; block, i-index, j-index
```

Each `probeN = <name>` entry points to a `[<name>]` section.  The
probe location is given either by index (`index-position = b i j`)
or by physical coordinates (`position = x y`).  Each point probe
writes a text file `OUTPUT/exit.txt` with columns for time and the
requested variables.

**Line probes** — Q2D additionally provides *line-sampling* probes
that record a 1D profile along an $i$- or $j$-fixed grid line at the
requested cadence:

```ini
[Q2D-LineProbes]
lineprobe1 = midrow

[midrow]
variables      = rho(1) u v p T M
diter          = 100    ; cadence in iterations
index-position = 1 2    ; block, direction (1 = i fixed, 2 = j fixed)
fixed-index    = 25     ; the constant i or j index defining the line
```

Each `lineprobeN = <name>` entry points to a `[<name>]` section.
`index-position` takes the block index and the line direction:
`1` samples along $j$ at fixed $i$, `2` samples along $i$ at fixed
$j$; `fixed-index` is the constant index of the sampled line.  Both
keys are mandatory.  Each line probe writes a Tecplot file
`OUTPUT/midrow.tec` (one zone per output instant).  Line probes are
implemented in `io/IO_LineProbes.f90`; see
`test/2D/euler/prandtl-meyer-lineprobe` for a working example.

---

## Restart

To resume from a previous run, set `newrun = false` in
`[Q2D-Parameters]`.  Q2D will read the latest solution file in
`OUTPUT/` and restore the flow field, physical time, and iteration
counter.  Multi-grid restart files follow the same `field_x*.tec`
naming convention as the IC files.

---

## Smoke-test execution without `INPUT/`

Q2D can be launched without an `INPUT/input.ini` file for CI smoke
testing.  When `Read_Ini` finds that `input.ini` is missing, the
solver prints its banner and exits cleanly (bootstrap-guard):

```bash
cd /tmp/q2d-empty-case      # no INPUT/ here
/path/to/bin/Q2D            # exits 0 with banner + no-op message
```

This capability is useful for build-system verification and
continuous integration without setting up a full test case.  See
[Code Structure — Solver Pipeline](../development/structure.md#solver-pipeline)
for implementation details.
