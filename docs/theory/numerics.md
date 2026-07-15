# Spatial Discretization

Q2D uses a cell-centred finite volume method (FVM) on structured
multi-block grids.  This page describes the discretization framework,
the MUSCL reconstruction that provides second-order accuracy, the
slope limiters that ensure monotonicity, the shock-detection sensor
that adaptively reduces the order near discontinuities, and the
Q2D-specific extensions wired in `numerics/sources/` and in
`Mod_Fluxes`.

---

## Finite Volume Framework

The integral form of the conservation law over a 2D cell $\Omega_i$ with
boundary $\partial\Omega_i$ reads

$$
\frac{\mathrm{d}}{\mathrm{d}t}
\int_{\Omega_i}\!\mathbf{U}\,\mathrm{d}V
+ \oint_{\partial\Omega_i}\!\mathbf{F}^{c}\!\cdot\!\hat{\mathbf{n}}\,\mathrm{d}A
= \oint_{\partial\Omega_i}\!\mathbf{F}^{v}\!\cdot\!\hat{\mathbf{n}}\,\mathrm{d}A
+ \int_{\Omega_i}\!\mathbf{S}\,\mathrm{d}V
$$

Approximating the integrals with midpoint quadrature yields the
semi-discrete update

$$
\frac{\mathrm{d}\mathbf{U}_i}{\mathrm{d}t}
= -\frac{1}{V_i}\sum_{f}\bigl(\mathbf{F}^{c}_f
  - \mathbf{F}^{v}_f\bigr)\,A_f
+ \mathbf{S}_i
$$

where $V_i$ is the cell **area** in 2D planar and 2D axisymmetric, or
$V_i = A_i\,h_i$ in Q2D area-variation (with $A_i$ the 2D Jacobian and
$h_i$ the local out-of-plane height).  The face sum runs over the
**four** faces $f$ of each quadrilateral cell — there are no top/bottom
faces in Q2D.

The numerical algorithm consists of the following steps:

- **Convective fluxes**: left and right states at each face are
  reconstructed via MUSCL; an approximate or exact Riemann solver
  returns the numerical convective flux $\mathbf{F}^{c}_f$ from these
  states and the face normal vector (see
  [Riemann Solvers](riemann-solvers.md)).
- **Diffusive fluxes**: velocity, temperature, and species gradients
  are computed at each face with a 6-point stencil (2 normal + 4
  tangential — only one tangential direction in 2D)
  and mapped to Cartesian coordinates using the face metric tensor.
- **Geometric sources** (Q2D): axisymmetric and
  area-variation residuals are added cell-wise after the face-flux
  divergence, dispatched caller-side in `Mod_Explicit` through the
  `obj_source%models(:)` polymorphic registry.

---

## MUSCL Reconstruction

First- and second-order spatial accuracy is available.  Second order
is obtained through piecewise-linear reconstruction (MUSCL — Monotone
Upstream-centred Schemes for Conservation Laws).

At each cell interface the algorithm uses a four-point stencil
$\{i{-}1,\;i,\;i{+}1,\;i{+}2\}$ along the face-normal direction and
proceeds as follows:

1. Compute slopes from consecutive cell values and physical spacings
   $\Delta l_0,\,\Delta l_1,\,\Delta l_2$:

$$
s_0 = \frac{P_i - P_{i-1}}{\Delta l_0},\qquad
s_1 = \frac{P_{i+1} - P_i}{\Delta l_1},\qquad
s_2 = \frac{P_{i+2} - P_{i+1}}{\Delta l_2}
$$

2. Limit the slopes to prevent spurious oscillations near
   discontinuities:

$$
\bar{s}_L = \phi(s_1,\, s_0),\qquad
\bar{s}_R = \phi(s_2,\, s_1)
$$

where $\phi$ is the chosen limiter function (see below).

3. Reconstruct interface values:

$$
P_L = P_i     + \beta\,\bar{s}_L\,\delta l_L,\qquad
P_R = P_{i+1} - \beta\,\bar{s}_R\,\delta l_R
$$

where $\delta l_L,\,\delta l_R$ are the distances from the cell
centres to the interface and $\beta$ is a blending parameter
controlled by the shock detector.

!!! note "Role of β"

    $\beta = 1$ recovers the full second-order MUSCL scheme.
    $\beta = 0$ reduces reconstruction to first order (donor-cell),
    used in the immediate vicinity of shocks.  The shock-detection
    sensor (see below) smoothly modulates $\beta$.

---

## Flux Limiters

Several limiter functions $\phi(a,b)$ are available.  All are
expressed in terms of the slope ratio $r = b\,/\,a$.

**MinMod**

$$
\phi(a,b) =
\begin{cases}
\operatorname{sign}(a)\,\min(|a|,\,|b|) & \text{if } a\,b > 0 \\
0 & \text{otherwise}
\end{cases}
$$

Most diffusive TVD limiter; strictly enforces monotonicity.

**Van Leer**

$$
\phi(a,b) = \frac{r + |r|}{1 + |r|}\;a
$$

Smooth, differentiable limiter.  Good balance of dissipation and
accuracy.

**Van Albada**

$$
\phi(a,b) = \frac{r^2 + r}{1 + r^2}\;a
$$

Smooth, less dissipative than van Leer near local extrema.

**MC (Monotonicity-Centred)**

$$
\phi(a,b) =
\max\!\Bigl(0,\;\min\!\bigl(2r,\;\tfrac{1}{2}(1+r),\;2\bigr)\Bigr)\;a
$$

Symmetric TVD limiter that permits the steepest gradients within the
TVD region.

**Superbee**

$$
\phi(a,b) =
\max\!\Bigl(0,\;\min(2r,\,1),\;\min(r,\,2)\Bigr)\;a
$$

Sharpest TVD limiter; excellent at capturing steep gradients but may
introduce mild staircase artefacts.

!!! info "LIMO3 — high-resolution third-order limiter"

    Q2D ships an **additional** limiter, LIMO3, beyond the five
    entries above.  LIMO3 enables a high-resolution third-order
    monotonic reconstruction suitable for moderately smooth flows
    where standard TVD limiters are excessively diffusive.  Selected
    at run time via the input file (`flux-limiter = LIMO3`).
    Implementation: `Lib_Limiters.f90` (entry 7).

### TVD Diagram

| Limiter | Dissipation | Smoothness |
|---------|:-----------:|:----------:|
| MinMod | High | $C^0$ |
| Van Leer | Medium | $C^\infty$ |
| Van Albada | Medium-low | $C^\infty$ |
| MC | Low | $C^0$ |
| Superbee | Lowest | $C^0$ |
| **LIMO3** (Q2D) | Low (3rd-order) | $C^0$ |

### Selection guidance

- **`vanleer`** — the workhorse default: smooth ($C^\infty$),
  well-balanced dissipation; used by most of the V&V suite
  (24 of the committed test inputs).
- **`mc`** — sharper gradients than van Leer at slightly higher
  oscillation risk; gated on [Sod79](../vv/1D.md) and the
  [Chang nozzle](../vv/2D-plane-extras.md) cases.
- **`minmod`** — the most dissipative TVD choice; use when robustness
  matters more than resolution.
- **`LIMO3`** — third-order accuracy in smooth regions with shock
  capturing; the [Shu–Osher case](../vv/1D-shuosher.md) gates its
  total-variation ordering against minmod/van Leer/MC baselines.
- **`superbee`** — the sharpest TVD limiter, but prone to
  over-steepening of smooth profiles; **not yet exercised by the
  ctest suite** (see the
  [Coverage Matrix known gaps](../vv/coverage-matrix.md#known-gaps-prioritized)).

---

## Shock Detection

Near strong shocks the MUSCL reconstruction can produce oscillations
even with TVD limiters (carbuncle phenomenon).  Q2D uses a
Jameson-type pressure sensor to detect shocks and smoothly degrade to
first order in their vicinity.

### Sensor formulation

A $3\times 3$ stencil of pressures centred on cell $(i,j)$ is used to
compute second-order undivided differences in each in-plane coordinate
direction (two directions only, no $\kappa_3$):

$$
\kappa_d =
\frac{|p_{+} - 2p_0 + p_{-}|}{p_{+} + 2p_0 + p_{-}},
\qquad d = 1, 2
$$

The shock strength is the maximum over the two in-plane directions:

$$
\sigma = \max(\kappa_1,\;\kappa_2)
$$

### Blending function

A smooth blending is applied using a tuneable threshold $\Delta = 20$:

$$
\beta =
\begin{cases}
1 - \tanh\!\bigl(10\,(\sigma\,\Delta)^3\bigr)
  & \text{if } \sigma < 1/\Delta \\[4pt]
0 & \text{otherwise}
\end{cases}
$$

- $\beta \to 1$ in smooth flow — full second-order reconstruction
- $\beta \to 0$ at shocks — first-order (no reconstruction)

This $\beta$ is passed directly to the MUSCL reconstruction step,
where it multiplies the limited slope.

!!! note "Solver-level integration"
    Some Riemann solvers also use the shock-detection flag internally
    to tune their numerical dissipation, providing an additional layer
    of robustness (see [Riemann Solvers](riemann-solvers.md)).

---

## Gradient Computation for Diffusive Fluxes

Velocity, temperature, and species gradients needed for the viscous
flux are computed with a **6-point stencil** at each cell interface
(the 2D reduction of a full 3D 10-point stencil):

- 2 points in the face-normal direction (the two cells sharing the
  face)
- 4 points in the single in-plane tangential direction (two per side)

The gradient in computational space $(\xi, \eta)$ is mapped to
Cartesian coordinates $(x, y)$ using the face metric tensor
$M_{2\times 2}$, which is computed from the grid geometry and stored
at each face.  The same metric is used to compute the physical
spacing $\Delta l$ for the MUSCL reconstruction.

In Q2D area-variation, the metric tensor additionally carries the
height-gradient field $\mathbf{G} = -\nabla h / h$ used for Definition
A coupling (see
[Q2D Area-Variation Derivation](area-variation-derivation.md)).  The
face-averaged contraction $\mathbf{V}\!\cdot\!\mathbf{G}$ is computed
by `V_dot_G_face_avg` and passed as an optional argument to
`Stress_Vector` and `Strain_Tensor` along the Definition A path.

---

## Q2D-specific extensions

### Geometric source dispatch

The axisymmetric and area-variation source terms are evaluated
caller-side in `Mod_Explicit` (a Q2D-only state machine):

```fortran
! Mod_Explicit.f90 — pseudo-code
do source_kind = 1, obj_source%n_models
  call obj_source%models(source_kind)%compute( &
       Pg, R, obj_metrics, obj_chem, ...)
end do
```

The `obj_source%models(:)` registry holds polymorphic
`class(SourceTerm)` extensions
(`AxisymmetricSource`, `AreaVariationSource`); registration is done at
setup time via `Register_Sources.f90`.

### OpenMP do/task dispatcher

The Q2D `Mod_Fluxes` module exposes **two** flux-loop variants that
coexist in the binary and are selected at run time via the INI key
`[Q2D-Numerics] OMP-parallel = do | task`:

| Variant | Pragma pattern | When preferred |
|---|---|---|
| `Fluxes_omp_do` | `!$omp parallel` outer + `!$omp do collapse(2)` intra-block | Default — diagnostic-friendly |
| `Fluxes_omp_task` | `!$omp single + taskgroup` block-level + `!$omp taskloop collapse(2)` intra-block + explicit `taskwait` | Imbalanced block decomposition; HPC tuning |

The dispatcher `Assign_Fluxes('do' | 'task')` selects the variant at
run time.

### Borda-flow inflow nozzle BC

For choked-nozzle inflow conditions following Fievisohn (2016, 2017),
Q2D adds a **Borda-flow** generalised Riemann-invariant treatment
selected by the sentinel `BC_Mach == -20` in the input.  This gates a
dedicated branch in `Lib_BC_Fluxes_Inflow_Nozzle.f90` that applies the
Glaister 1988 generalised Riemann invariant with full velocity
magnitude $H_0$.

---

## OpenMP and loop structure

All Q2D loops are 2D-collapsed.  The canonical pattern is

```fortran
!$omp parallel do collapse(2) private(i, j, ...) &
!$omp&         reduction(+:residuotot) default(none)
do j = 1, Nj
  do i = 1, Ni
    ! body — no `do k`, no (:, :, 1) slice
  end do
end do
```

There is no `Nk` in `Global_m`, no `do k` in the source, and no
`(:, :, 1)` slice — these are structural invariants of the 2D
(No-$w$) reduction.

---

## Known limitations

!!! note "Soot thermophoretic enthalpy in the energy budget"
    The diffusive energy flux (`Compute_Diffusive_Flux` in
    `Lib_Diffusive.f90`) accumulates the species enthalpy contribution
    $\sum_\alpha h_\alpha\, j_\alpha^\text{Fick}$ only over the
    gas-phase species loop $s = 1\ldots n_{sc}$.  The soot
    thermophoretic mass flux
    $j_\text{thermo} = -C_T\,\mu\,Y_s\,\nabla T / T$ is added to the
    soot mass flux in `Mod_Soot::Soot_Diffusive_Flux`, but the
    corresponding enthalpy contribution $h_\text{soot}\,j_\text{thermo}$
    is **not** added to the energy face-flux `Flux(np)`, nor to the
    cylindric residual `Source(np)` (axisymmetric mode), nor to the
    area-variation residual `Source(np)` (Definition A).  The omission
    activates only when soot is enabled **and** $Y_s > 0$ **and**
    $\nabla T \neq 0$; for runs without soot it is inactive.
    Engineering significance: $|h_\text{soot}\,j_\text{thermo}|$ is of
    the same order as the Fourier heat flux $q = -k\,\nabla T$ in
    soot + thermal-gradient regions (flames, boundary layers), so it
    is **not sub-leading** for soot-active configurations.

---

## References

1. B. van Leer, "Towards the ultimate conservative difference scheme.
   V. A second-order sequel to Godunov's method," *J. Comput. Phys.*,
   32(1), 1979.
2. P. K. Sweby, "High resolution schemes using flux limiters for
   hyperbolic conservation laws," *SIAM J. Numer. Anal.*, 21(5), 1984.
3. A. Jameson, W. Schmidt, E. Turkel, "Numerical solution of the Euler
   equations by finite volume methods using Runge–Kutta time-stepping
   schemes," AIAA-81-1259, 1981.
4. R. Tramel, R. Nichols, P. Buning, "Addition of improved
   shock-capturing schemes to OVERFLOW 2.1," 19th AIAA Computational
   Fluid Dynamics, 2009. AIAA Paper 2009-3988.
5. P. K. Glaister, "An approximate linearised Riemann solver for the
   Euler equations for real gases," *J. Comput. Phys.*, 74(2), 1988
   — generalised Riemann invariant used by the Borda-flow nozzle BC.
6. R. T. Fievisohn, K. H. Yu, "Steady-state analysis of rotating
   detonation engine flowfields with the method of characteristics,"
   *J. Propul. Power*, 33(1), 2017, 89–99.
7. R. T. Fievisohn, *Development and Application of Theoretical
   Models for Rotating Detonation Engine Flowfields*, PhD
   dissertation, University of Maryland, College Park, 2016.
