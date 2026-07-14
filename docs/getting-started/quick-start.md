# Quick Start: Prandtl–Meyer Expansion

This guide walks you through running the **Prandtl–Meyer expansion**
test case — a canonical 2-D supersonic flow verification that
demonstrates corner-based isentropic expansion with analytic
validation.  The case runs in the **2D planar** Q2D regime (no
axisymmetric or area-variation source).

!!! note "Prerequisites"
    - Q2D built and available at `bin/Q2D` (see [Installation](installation.md))
    - Python 3 with NumPy (for post-run verification)

---

## Physics Overview

The Prandtl–Meyer expansion describes the isentropic turning of a
supersonic flow over a smooth convex corner.  As the flow expands,
the Mach number increases and pressure decreases while entropy
remains constant.

**Canonical parameters:**

| Parameter | Value | Notes |
|-----------|-------|-------|
| Upstream Mach | $M_1 = 2.0$ | Supersonic inlet condition |
| Expansion angle | $\theta \approx 35°$ | Corner deflection |
| Specific heat ratio | $\gamma = 1.4$ | Diatomic gas (air) |
| Upstream T, p | $T_1 = 1$ K, $p_1 = 1$ Pa | Normalised reference state |

**Expected result:** the downstream Mach $M_2$ computed by Q2D
matches the analytic Prandtl–Meyer prediction to within **0.5 %**
relative error.

The analytic downstream Mach is found by solving:

$$\nu(M_2) = \nu(M_1) + \theta$$

where the Prandtl–Meyer function is:

$$\nu(M) = \sqrt{\frac{\gamma + 1}{\gamma - 1}} \arctan\sqrt{\frac{\gamma - 1}{\gamma + 1}(M^2 - 1)} - \arctan\sqrt{M^2 - 1}$$

---

## Case Directory Structure

The case is located at `test/2D/euler/prandtl-meyer/`:

```
test/2D/euler/prandtl-meyer/
├── input.ini           # Solver configuration
├── Q2D.sh              # Run script
├── check.py            # Verification script
├── INPUT/
│   ├── ic.tec          # Initial condition (mesh + field data)
│   └── bc.txt          # Boundary condition table
├── MESH/
│   └── mesh.tec        # Grid in Tecplot format
└── OUTPUT/             # (Created at runtime)
```

Shared thermodynamic data (phase table, transport properties) lives
in `test/common/Air/` and is symlinked at runtime.

### Input Files Description

- **`input.ini`** — solver configuration specifying physics (Euler
  equations), numerics (time scheme, Riemann solver,
  reconstruction), and I/O preferences.
- **`ic.tec`** — initial condition file in Tecplot format containing
  both the computational mesh and field data (density, pressure,
  velocity components $(u, v)$ — no $w$ stored).
  The mesh is embedded here, so the separate `MESH/` directory is
  optional.
- **`bc.txt`** — boundary condition table defining edge or face-wise
  boundary treatment (e.g., inlet, exit, slip wall) referenced by
  the solver.  The Q2D `bc.txt` reader transparently accepts both
  legacy 5-field and MOSE-form 6-field record widths via
  `IO_BC::detect_record_width`.
- **Thermodynamic data (via FLINT)** — pre-compiled phase table,
  transport coefficients, and species properties from FLINT
  (located in `test/common/Air/`) provide material properties needed
  by the solver at runtime.

---

## Solver Configuration

The `input.ini` file specifies the physics and numerics:

```ini
[Q2D-Physics]
equations = euler            # Inviscid compressible Euler

[Q2D-Numerics]
cfl = 1.0                    # Courant–Friedrichs–Lewy number
time-scheme = RK3            # 3rd-order Runge–Kutta temporal scheme
time-accurate = false        # Steady-state integration (no physical time)
irs-beta = 0.5               # Implicit residual smoothing parameter
riemann-solver = HLLE        # HLLE numerical flux
space-reconstruction = MUSCL # 2nd-order spatial reconstruction
flux-limiter = vanleer       # van Leer slope limiter
integration-variables = cons # Conserved variables (ρ, ρu, ρv, E)
omp_name = do                # 'do' or 'task' OMP loop dispatcher

[Q2D-IO]
shell-diter = 100            # Diagnostic output every 100 iterations
sol-format = tecplot ascii   # Tecplot ASCII output format

[Q2D-Probes]
probe1 = exit                # Probe at the exit plane

[exit]
index-position = 1 50 33 1   # Probe location (b, i, j, k=1) — k always 1 in Q2D
variables = M                # Mach number output
diter = 100                  # Probe output frequency
```

Key choices:

- **Steady-state integration** (`time-accurate = false`) — physical
  time is not needed; we iterate to steady state.
- **HLLE solver** — robust, entropy-satisfying choice for expansion
  waves.
- **van Leer limiter** — smooth and non-oscillatory.
- **Probe output** — captures Mach number at the exit plane for
  verification.

---

## Running the Simulation

### Basic Execution

Navigate to the case directory and run:

```bash
cd test/2D/euler/prandtl-meyer
./Q2D.sh solve
```

The `Q2D.sh` script:

1. Symlinks shared thermodynamic data (`phase.txt`, `thermo.dat`)
   into `INPUT/`.
2. Creates `OUTPUT/` and `bin/` directories if needed.
3. Launches the solver executable.

**Typical runtime:** 5–30 seconds (depending on mesh size and
hardware).

### Parallel Execution

To run with N OpenMP threads:

```bash
./Q2D.sh solve -p 4    # Run with 4 threads
```

For best performance, use `N = (number of physical cores)` on your
machine.

---

## Post-Run Verification

After the run completes, verify that the CFD result matches the
analytic prediction:

```bash
python3 check.py
```

This script:

1. Reads `OUTPUT/exit.txt` (probe data from the exit plane).
2. Computes the expected Mach using the Prandtl–Meyer function with
   $M_1 = 2.0$ and the expansion angle.
3. Compares the CFD Mach (`M_CFD`) against the analytic Mach
   (`M_analytical`).
4. Reports **success** (exit code 0) if relative error < 0.5 %,
   **failure** (exit code 1) otherwise.

**Example output:**

```
Correct!
```

or on failure:

```
Mach number mismatch (CFD vs analytical)
3.42 3.50
Error = 2.3%, threshold = 0.5%
```

### Inspecting the Solution

The flow field is written to `OUTPUT/field.tec` (Tecplot ASCII
format) at the end of the run.  You can visualise this with
Paraview.

Key variables in the output: density $\rho$, velocity $(u, v)$
(plus `w = 0` placeholder for tooling compatibility), pressure $p$,
temperature $T$.

---

## Understanding the Output

### Diagnostics on screen

During execution, Q2D prints iteration count, residuals, and solver
status.  Example:

```
  Iteration:    100  Time:  0.00000E+00  RMS Res:  1.23E-02
  Iteration:    200  Time:  0.00000E+00  RMS Res:  4.56E-03
  Iteration:    300  Time:  0.00000E+00  RMS Res:  1.89E-03
  ...
  Converged: RMS Res < 1E-6
```

Watch for the RMS residual decreasing smoothly.  If it stalls or
grows, the simulation is not stable.

### Output Files

| File | Content |
|------|---------|
| `OUTPUT/field.tec` | Full 2-D field data (ρ, u, v, p, T, M, ...) at convergence |
| `OUTPUT/exit.txt` | Probe data (Mach vs. probe location) from the exit boundary |
| `OUTPUT/diagnostic.tec` | Convergence history and diagnostic quantities |

---

## Next Steps

Now that you understand the basic workflow:

- **Run other 2-D cases:** see `test/2D/euler/` for additional
  verification problems (oblique shock, hypersonic cylinder, rocket
  nozzle, ...).
- **Try the q2D area-variation regime:** use a case with an
  `INPUT/height.tec` height field to exercise the geometric source
  terms and Definition A V·G coupling.
- **Study the code:** refer to [Code Structure](../development/structure.md)
  and [Contributing](../development/contributing.md) for developer
  documentation.

---

## References

- **Prandtl–Meyer expansion:** [NACA Report 1135](https://ntrs.nasa.gov/citations/19930090976).
- **Riemann solvers:** Toro, E. F. *Riemann Solvers and Numerical
  Methods for Fluid Dynamics* (3rd ed.).
- **MUSCL reconstruction:** van Leer, B. "Towards the Ultimate
  Conservative Difference Scheme. II. Monotonicity and Conservation
  Combined in a Second-Order Scheme" (1974).
