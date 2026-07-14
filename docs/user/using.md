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
the `[Q2D-Source]` section:

```ini
[Q2D-Source]
axisymmetric   = .false.   ; cylindric residual −1/r·(…)
area_variation = .false.   ; q2D depth-averaged residual ∇h/h·(…) + Definition A V·G
```

| `axisymmetric` | `area_variation` | Resulting regime |
|:-:|:-:|---|
| `.false.` | `.false.` | **2D planar** (default) — Cartesian flow |
| `.true.`  | `.false.` | **2D axisymmetric** (no swirl) — cylindric residual |
| `.false.` | `.true.`  | **q2D area-variation** — `INPUT/height.tec` required |

The selection drives the dispatch of the geometric source procedure
pointers and the `Stress_Vector_*Std` / `_Axi` viscous closures; see
[Theory — Geometric Source Terms](../theory/geometric-sources.md).

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
│   ├── height.tec         ← (q2D area-variation only) prescribed h(x,y)
│   └── bc.txt             ← face-by-face boundary condition data
└── OUTPUT/                ← created at runtime
    ├── field.tec          ← solution fields (ρ, u, v, w=0, p, T, …)
    ├── field_0001.tec     ← solution snapshots (when sol-overwrite = false)
    ├── ...
    ├── wall.tec           ← wall quantities (viscous cases only)
    ├── residual-history.dat  ← residual norm vs. iteration
    ├── probe_<name>.txt   ← time-history at point probes
    └── lprobe_<name>.txt  ← time-history at line probes
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
    run time via `[Q2D-Numerics] omp_name`:

    - `omp_name = do` (default) — `!$omp do collapse(2)` pattern.
      Diagnostic-friendly.
    - `omp_name = task` — `!$omp taskloop collapse(2)` pattern with
      explicit `taskwait`.  Beneficial for imbalanced block
      decomposition.

    See [Theory — Spatial Discretization](../theory/numerics.md#openmp-dotask-dispatcher).

---

## Output

After a run completes, the `OUTPUT/` directory contains:

| File | Description |
|------|-------------|
| `field.tec` | Solution fields ($\rho$, $u$, $v$, $w \equiv 0$, $p$, $T$, …) |
| `residual-history.dat` | Iteration-by-iteration residual norms |
| `wall.tec` | Wall quantities (skin friction, heat flux) — viscous cases only |
| `probe_*.txt` | Time-history data at point probes |
| `lprobe_*.txt` | Time-history data at line probes |

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
variables = p T u
dtime = 1e-4
```

Each point probe writes a text file `OUTPUT/probe_exit.txt` with
columns for time and the requested variables.

**Line probes** — Q2D additionally provides *line-sampling* probes
that record a 1D profile along an $i$- or $j$-fixed line at the
requested cadence:

```ini
[Q2D-LineProbes]
lprobe1 = exit_line

[exit_line]
variables = p T u v
direction = j           ; sample along the j-line of fixed i
i         = 50          ; i-index
diter     = 100         ; cadence in iterations
```

Each line probe writes a text file `OUTPUT/lprobe_exit_line.txt`.
Line probes are implemented in `io/IO_LineProbes.f90`.

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
