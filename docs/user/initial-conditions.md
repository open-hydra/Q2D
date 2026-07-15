# Initial Conditions

Q2D supports two ways to specify the initial flow field:

1. **Inline** — piecewise-constant state(s) defined directly in
   `input.ini` via `[ICB-Block*]` sections.
2. **File-based** — a Tecplot ASCII (`.tec`) file containing the full
   grid geometry and an initial flow field.

The two methods can coexist: the grid is always read from
`INPUT/ic.tec` (or equivalent), while the flow field may be
overwritten by the inline `[ICB-Block*]` definitions.

---

## Multi-grid Naming Convention

When multigrid is enabled, Q2D expects one IC file per grid level:

| File | Level |
|------|-------|
| `INPUT/ic_x1.tec` | Finest grid |
| `INPUT/ic_x2.tec` | 2× coarser |
| `INPUT/ic_x4.tec` | 4× coarser |
| `INPUT/ic_x8.tec` | 8× coarser |

For single-grid cases the file is simply `INPUT/ic.tec`.

---

## File-Based Initial Conditions

### Tecplot ASCII Format (`.tec`)

The IC file is a standard Tecplot ASCII structured file.  It stores
both the **grid node coordinates** and the **cell-centred flow
variables** for each block.

**Variable list:**

```
VARIABLES = "x"  "y"  "z"  "rho1" ["rho2" ...] "u"  "v"  "w"  "p"
```

- `x`, `y`, `z` — Cartesian node coordinates (NODAL).
- `rho1`, `rho2`, … — partial densities for each species
  $[\text{kg/m}^3]$ (CELLCENTERED).
- `u`, `v`, `w` — velocity components $[\text{m/s}]$
  (CELLCENTERED).
- `p` — static pressure $[\text{Pa}]$ (CELLCENTERED).

For a single-species gas there is only `rho1`.  For a mixture of
$n_s$ species the density variables are `rho1 rho2 ... rho{ns}`,
ordered as in `phase.txt`.

!!! note "Q2D 2D format compatibility"

    Q2D writes and reads `.tec` files with the same MOSE-form variable
    list including `z` and `w` for tooling compatibility.  However:

    - The `z` coordinate is uniform within a block (degenerate
      out-of-plane coordinate) and used purely for the structured-grid
      shape.
    - The `w` velocity column is **always zero** and is ignored by the
      solver (No-$w$ principle).  Tecplot/VTK
      visualisations that read `w` will see zero out-of-plane
      velocity.
    - In Q2D area-variation cases, the local height $h(x, y)$ is
      provided as a separate `INPUT/height.tec` file, **not** as the
      `z` coordinate.

**Zone header:**

Each block is stored as one Tecplot zone with the following
attributes:

```
ZONE T = "B{m}-IG", I={ni+1}, J={nj+1}, K={nk+1},
     DATAPACKING=BLOCK,
     VARLOCATION=([1-3]=NODAL, [4-{4+ns+3}]=CELLCENTERED),
     SOLUTIONTIME={t}
```

| Field | Description |
|-------|-------------|
| `T` | Zone title — must be `B{m}-IG` where `{m}` is the 1-based block index |
| `I`, `J`, `K` | Number of **nodes** in each direction (cells + 1 per dimension); `K = 2` for Q2D 2D blocks |
| `DATAPACKING` | Must be `BLOCK` |
| `VARLOCATION` | Coordinates are NODAL; all flow variables are CELLCENTERED |
| `SOLUTIONTIME` | Physical time of the snapshot (used for restart) |

**Data layout:**

Data are written in block order — all values of variable 1, then all
values of variable 2, etc.:

```
x(1,1,1)  x(2,1,1)  ...  x(ni+1,nj+1,nk+1)    ← node coords, (ni+1)×(nj+1)×(nk+1) values
y(1,1,1)  ...
z(1,1,1)  ...
rho1(1,1,1)  rho1(2,1,1)  ...                   ← cell-centred, ni×nj×nk values
u(1,1,1)  ...
v(1,1,1)  ...
w(1,1,1)  ...                                   ← always zero in Q2D
p(1,1,1)  ...
```

**Minimal 1-D example (100 cells, single block):**

```
 VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"
 ZONE  T = B1-IG, I=101, J=2, K=2, DATAPACKING=BLOCK,
       VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0
 0.000000E+000
 1.000000E-002
 ...
 1.000000E+000
 0.000000E+000    ← y coordinates (all zeros for 1D)
 ...
 0.000000E+000    ← z coordinates
 ...
 1.225000E+000    ← rho1 (100 values)
 ...
 0.000000E+000    ← u
 ...
 0.000000E+000    ← v
 0.000000E+000    ← w  (zero in Q2D)
 1.01325E+005     ← p
 ...
```

### Tecplot Binary (`.szplt`)

Compact binary format produced by TecIO, selected with
`ini-format = tecplot binary` in `[Q2D-IO]`. Recommended for large
multi-block cases where ASCII I/O becomes a bottleneck; requires
ORION built with TecIO support (`USE_TECIO=ON`, the default local
build — the CI builds use `USE_TECIO=OFF` and are ASCII-only). The
variable list and zone structure are identical to the ASCII format.

### VTK Structured (`.vts` / `.vtm`)

Open-standard XML format readable by ParaView and VisIt, selected
with `ini-format = vtk ascii` or `vtk raw` in `[Q2D-IO]`. Multi-block
cases use a `.vtm` container referencing individual `.vts` structured
pieces under `INPUT/vtk/` (Q2D reads them through ORION's
`vtk_read_structured_multiblock`; see `src/lib/io/IO_Solution.f90`).
Variable names and layout follow the same conventions as the Tecplot
formats.

## Required Variables

Every file-based IC must contain node coordinates and the
cell-centred flow state:

| Variable | Description | Units |
|----------|-------------|-------|
| `x`, `y` | Node coordinates | m |
| `rho1` … `rhoN` | Partial density of each species | kg/m³ |
| `u`, `v` | Velocity components | m/s |
| `p` | Pressure | Pa |

If a turbulence model is active, the IC file must also carry its
transported variables, appended after `p` (the reader maps variables
positionally; the names below are the ones Q2D's own solution writer
emits, `src/lib/io/IO_Solution.f90`):

| Model (`turbulence =`) | Additional variables | Example case |
|---|---|---|
| `SA` family (1 eq.) | `mi_tilde` | `test/2D/viscous/swbli` |
| `SST`, `Wilcox2006` (2 eq.) | `kappa`, `omega` | `test/2D/euler/planar-nozzle-turbulent/kw` |
| `SSGLRR(-SD)` (5 eq.) | `ruu`, `rvv`, `rww`, `ruv`, `omega` | `test/2D/euler/planar-nozzle-turbulent/rsm` |

Because the mapping is positional, an IC written for a different
model family (wrong number of transported variables) is not a valid
restart file for the active model.

---

## Inline Initial Conditions (`[ICB-Block*]`)

For simple cases the flow field can be specified directly in
`input.ini` using `[ICB-Block*]` sections.  Q2D still needs the IC
file for the grid geometry, but the cell-centred data are overwritten
by the inline values.

### Homogeneous State

A single uniform state over an entire block:

```ini
[ICB-Block1]
type  = homogeneous
p     = 97250.0      ; static pressure [Pa]
T     = 300.0        ; static temperature [K]
mach  = 0.2          ; Mach number
```

| Key | Description |
|-----|-------------|
| `type` | Must be `homogeneous` (default if omitted) |
| `p` | Static pressure $p$ [Pa] |
| `p0` | Total pressure $p_0$ [Pa] (alternative to `p`) |
| `T` | Static temperature $T$ [K] |
| `T0` | Total temperature $T_0$ [K] (alternative to `T`) |
| `mach` | Mach number (used to derive velocity from $T$ and $\gamma$) |
| `u`, `v` | In-plane velocity components $[\text{m/s}]$ (alternative to `mach`) — `w` is **not** accepted |
| `kappa` | Turbulent kinetic energy $k$ (for RANS, if turbulence model active) |
| `omega` | Specific dissipation $\omega$ (for SST/k-ω, if turbulence model active) |

**Example — viscous subsonic flow:**

```ini
[ICB-Block1]
type  = homogeneous
p     = 97250.0
T     = 300.0
mach  = 0.2
kappa = 2.25e-3
omega = 6250.0
```

### Piecewise (Riemann-problem) State

For shock-tube or Riemann-problem configurations, the domain can be
split into two regions along an in-plane coordinate direction:

```ini
[ICB-Block1]
direction = x
zone1     = left
zone2     = right
range1    = 0.0 0.5
range2    = 0.5 1.0

[left]
p = 1.0e5
T = 348.43

[right]
p = 0.1e5
T = 278.60
```

| Key | Description |
|-----|-------------|
| `direction` | Splitting axis: `x` or `y` (no `z` in Q2D) |
| `zone1`, `zone2` | Names of sections that define each region |
| `range1`, `range2` | Coordinate range `[xmin xmax]` for each region |

Each named region (`[left]`, `[right]`) accepts the same keys as the
homogeneous state above.

### Nozzle State

A special initialisation for nozzle flows that sets a stagnation-based
uniform state:

```ini
[ICB-Block1]
type = nozzle
p0   = 10.342e5
T0   = 555.56
```

### Wildcard Assignment

To apply the same state to all blocks without listing them
individually, use the wildcard `[ICB-Block*]`:

```ini
[ICB-Block*]
p    = 1.0
T    = 1.0
mach = 4.0
```

This sets the same homogeneous state for every block in the domain.

---

## Geometric-regime considerations

The three Q2D geometric regimes (cf. [User Guide](index.md))
introduce specific requirements on the initial condition:

- **2D planar** — no additional requirements; the standard `(u, v)`
  primitives carry the full velocity field.
- **2D axisymmetric** (no swirl) — the IC file must use the
  $(x, r)$ Cartesian-style coordinates, with $r \ge 0$ the radial
  distance from the symmetry axis.  The radial velocity $v$ should
  satisfy $v(r = 0) = 0$ for regularity; near the axis Q2D applies a
  $r_\text{ave} \to 10^{-6}$ de l'Hôpital fallback in the source-term
  evaluation.
- **Q2D area-variation** — additional `INPUT/height.tec` file
  (Tecplot ASCII, same grid as `ic.tec`) prescribing the local
  out-of-plane height $h(x, y)$.  The depth-averaged primitive state
  in `ic.tec` is interpreted as the volume-averaged state over the
  height $h$.

---

## Restart

When resuming from a previous run (`newrun = false`), the initial
condition is read from the last solution file in `OUTPUT/` rather
than from `INPUT/`.  Q2D automatically detects the latest file and
restores the flow field, physical time, and iteration counter.  See
[Running Q2D — Restart](using.md#restart) for details.
