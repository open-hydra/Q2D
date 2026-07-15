# Single Rotating Frame (SRF) Theory

## Overview

The Single Rotating Frame (SRF) model allows simulation of the
compressible Navier–Stokes equations in a reference frame rotating
rigidly at constant angular velocity $\boldsymbol{\Omega}$ [rad/s]
about a fixed axis.  This is particularly useful for turbomachinery
applications where the computational domain rotates with an impeller
or rotor.

!!! warning "Q2D 2D-planar-only scope"

    Q2D supports rotating frames **only in 2D planar mode** with the
    rotation axis fixed to $\boldsymbol{e} = [0, 0, 1]$
    (perpendicular to the simulation plane).  In this configuration
    the Coriolis force collapses to a scalar
    $-2\Omega\,(\partial v_y/\partial x - \partial v_x/\partial y)$
    and the centrifugal acceleration acts in-plane.

    **2D-axi-swirl** (rotating frame about the axisymmetric axis with
    a non-zero $u_\theta$) is **out of scope**; it is reserved for a
    follow-up activity, requires dedicated pre-processing not included
    here, and is structurally incompatible with the 2D No-$w$ reduction.

## Velocity Definition

The velocity field stored in the solution vector is the **relative
(rotating-frame) velocity**:

$$\boldsymbol{w} = \boldsymbol{u}_{\text{abs}} - \boldsymbol{\Omega} \times \boldsymbol{r}$$

where:

- $\boldsymbol{u}_{\text{abs}}$ is the absolute (inertial-frame)
  velocity.
- $\boldsymbol{\Omega}$ is the angular velocity vector
  ($\boldsymbol{\Omega} = \Omega\,\boldsymbol{e}_z$ in Q2D).
- $\boldsymbol{r}$ is the position vector from the axis of rotation.

This formulation ensures that no modifications are required to the
standard convective or diffusive flux modules — they operate directly
on the relative velocity without change.

In Q2D the in-plane velocity stored in the primitives is
$(u, v) = (w_x, w_y)$, i.e. the Cartesian components of $\boldsymbol{w}$.

## Governing Equations in the Rotating Frame

### Continuity and Species Equations

The continuity and species mass-fraction equations are **unchanged**
from the inertial frame:

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{w}) = 0$$

$$\frac{\partial (\rho Y_k)}{\partial t} + \nabla \cdot (\rho Y_k \boldsymbol{w}) = -\nabla \cdot \boldsymbol{J}_k + \dot{m}_k$$

### Momentum Equation

The momentum equation includes additional source terms accounting for
Coriolis and centrifugal accelerations (per unit volume):

$$\frac{\partial (\rho \boldsymbol{w})}{\partial t} + \nabla \cdot (\rho \boldsymbol{w} \boldsymbol{w}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \boldsymbol{S}_{\text{mom}}$$

where the source term is:

$$\boldsymbol{S}_{\text{mom}} = -2\rho (\boldsymbol{\Omega} \times \boldsymbol{w}) + \rho \omega^2 \boldsymbol{r}_{\perp}$$

**Components:**

- **Coriolis force:** $-2\rho (\boldsymbol{\Omega} \times \boldsymbol{w})$ — in the Q2D 2D plane with $\boldsymbol{\Omega} = \Omega\,\boldsymbol{e}_z$ this collapses to $(-2\Omega\,\rho\,(-w_y,\,w_x))$.
- **Centrifugal force:** $\rho \omega^2 \boldsymbol{r}_{\perp}$ — purely in-plane in Q2D.

where $\boldsymbol{r}_{\perp} = \boldsymbol{r} - (\boldsymbol{r} \cdot \boldsymbol{e}) \boldsymbol{e}$ is the position vector perpendicular to the rotation axis (with $\boldsymbol{e}$ the unit vector along the axis), and $\omega = |\boldsymbol{\Omega}|$.

### Energy Equation

The energy equation (with specific total energy $E = c_v T + \frac{1}{2} |\boldsymbol{w}|^2$) receives a source term:

$$\frac{\partial (\rho E)}{\partial t} + \nabla \cdot (\rho E \boldsymbol{w}) = -\nabla \cdot p\boldsymbol{w} + \nabla \cdot (\boldsymbol{\tau} \cdot \boldsymbol{w}) - \nabla \cdot \boldsymbol{q} + \boldsymbol{w} \cdot \boldsymbol{S}_{\text{mom}}$$

The work done by body forces simplifies to:

$$\boldsymbol{S}_{E} = \boldsymbol{w} \cdot (\rho \omega^2 \boldsymbol{r}_{\perp})$$

**Note:** the Coriolis force does no work since
$\boldsymbol{w} \cdot (\boldsymbol{\Omega} \times \boldsymbol{w}) = 0$
identically.

### RANS Turbulence Equations

For **scalar** transported quantities ($\tilde\nu$, $k$, $\omega$) the
transport equations themselves are **unchanged** from the inertial
frame: the exact Coriolis contribution to the second-moment budget is
trace-free, so its contraction onto $k$ (and hence $\omega$, $\tilde\nu$)
vanishes identically (sympy-proved). For the
**tensor-valued** Reynolds-stress transport (SSGLRR) this is *not* true —
see the exact Coriolis production $G_{ij}$ below. In addition, every
closure that consumes the vorticity (rather than the strain rate) must
use the **absolute** vorticity. For constant
$\boldsymbol{\Omega}$ the kinematic identity

$$\boldsymbol{\omega}_{\rm abs} = \boldsymbol{\omega}_{\rm rel} + 2\boldsymbol{\Omega}$$

follows from $\nabla\times(\boldsymbol{\Omega}\times\boldsymbol{r}) = 2\boldsymbol{\Omega}$,
and Q2D applies it wherever vorticity enters a closure:

| Model | Rotating-frame terms |
|---|---|
| Spalart–Allmaras (baseline) | $\tilde S = \|\boldsymbol{\omega}_{\rm abs}\|$ in the production term |
| SA + Spalart–Shur (SA-RC) | $\boldsymbol{\omega}_{\rm abs}$ in $W_{ik}$, $r^*$, $D^2$ **and** the frame term $(\varepsilon_{imn}S_{jn}+\varepsilon_{jmn}S_{in})\Omega_m$ in the $\tilde r$ bracket (TMR SA-RC) |
| Menter SST $k$–$\omega$ (SST-RC) | $\boldsymbol{\omega}_{\rm abs}$ in $O = \sqrt{2 W_{ij}W_{ij}}$ **and** the same $\tilde r$-bracket frame term (Smirnov–Menter 2009) |
| Wilcox 2006 $k$–$\omega$ | $f_\beta$ vortex-stretching limiter $W_{ij}W_{jk}\hat S_{ki}$ with $\boldsymbol{\omega}_{\rm abs}$ |
| SA-QCR2000 | $O_{ik} = 2W_{ik}/\|\nabla\boldsymbol{u}\|$ built from **absolute** vorticity *and* absolute gradient norm (objectivity requires both) |
| SSG/LRR-$\omega$ RSM | $W_{ij}^{\rm abs}$ in the $C_5$ rapid pressure-strain **and** the exact Coriolis production $G_{ij}$ (below) |

In the Q2D 2D plane with $\boldsymbol{\Omega} = \Omega\,\boldsymbol{e}_z$
the vorticity is a scalar $\omega_z = \partial v/\partial x - \partial u/\partial y$
and the shift collapses to $\omega_{z,{\rm abs}} = \omega_{z,{\rm rel}} + 2\Omega$.

**Strain-rate-based production is frame-invariant.** The $S_{ij}$-driven
production used by SST and Wilcox 2006 needs no correction — only the
vorticity-dependent auxiliary terms above receive the $2\boldsymbol{\Omega}$
shift.

**Caveat for baseline SA at high $\Omega$.** Because $|\boldsymbol{\omega}_{\rm abs}| = 2|\boldsymbol{\Omega}|$
in a region of pure rigid-body co-rotation, baseline SA produces a
non-zero $\tilde S$ there and will generate eddy viscosity from the
frame rotation alone. This is the correct inertial-frame answer (the
absolute vorticity really is non-zero in that region), but in flows
where streamline curvature does not enhance turbulence the production
is spurious. SA-RC's $f_{r1}$ is designed to attenuate exactly this
contribution — so the recommended pairing in a rotating frame is
**SA-RC**, not baseline SA.

**SSG/LRR-$\omega$ Reynolds-stress model** (`turbulence = SSGLRR`): two
distinct rotating-frame terms are required, with no double counting
(Speziale–Sarkar–Gatski 1991 noninertial form):

1. the in-plane vorticity entering the $C_5$ rapid pressure-strain is
   shifted to the absolute frame,
   $\boldsymbol{\omega}_{\rm abs} = \boldsymbol{\omega}_{\rm rel} + 2\boldsymbol{\Omega}$
   (a *modelled* term);
2. the **exact Coriolis production** — obtained by taking moments of the
   fluctuating momentum equation, no closure involved —

   $$G_{ij} = -2\,\Omega_k\left(\varepsilon_{ikl}R_{lj} + \varepsilon_{jkl}R_{li}\right)$$

   which with axis $[0,0,1]$ and the Q2D 5-variable set reduces to
   $G_{11} = +4\Omega R_{12}$, $G_{22} = -4\Omega R_{12}$, $G_{33} = 0$,
   $G_{12} = 2\Omega(R_{22}-R_{11})$, with $G_{13}=G_{23}=0$ (the 7→5
   reduction stays closed under rotation). $G_{ij}$ is trace-free ($k$
   and $\omega$ rows unaffected) and **not absorbable** in the production
   built from relative gradients: $P(\nabla u_{\rm abs}) -
   P(\nabla u_{\rm rel}) = G/2$ exactly. It is the leading-order
   mechanism for rotation-induced stabilization/destabilization
   (e.g. the spanwise-rotating channel asymmetry).

Verification: `a sympy audit of the exact rotating-frame second-moment algebra` (exact
derivation + transcription vs the published Table-2 constants) and the
frame-invariance analyses: the rotating-frame model with $G_{ij}$
reproduces the rotated inertial-frame solution to the integrator's own
order (RK4 order 4.0), while dropping $G_{ij}$ leaves an $O(1)$,
$dt$-independent error; analogous gates hold for the $\tilde r$
bracket and QCR objectivity. The solver-level regression is the
[rotating-channel-rsm](../vv/2D-rotating-channel-rsm.md) V&V case.

The frame itself is configured in the `[Q2D-Rotating-Frame]`
section: `omega` (angular speed, rad/s), plus the registry
parameters `axis` (rotation-axis direction, default `0.0 0.0 1.0`)
and `origin` (a point on the axis, default `0.0 0.0 0.0`).

**References.**

- Speziale, C.G. (1989) *"Turbulence modeling in noninertial frames of
  reference,"* Theoretical and Computational Fluid Dynamics **1**,
  3–19. — General statement of why vorticity-dependent closures require
  the absolute vorticity in a rotating frame.
- Spalart, P.R. and Shur, M. (1997) *"On the sensitization of
  turbulence models to rotation and curvature,"* Aerospace Science and
  Technology **1**(5), 297–302. — SA-RC; absolute vorticity in
  $\tilde r$.
- Smirnov, P.E. and Menter, F.R. (2009) *"Sensitization of the SST
  turbulence model to rotation and curvature by applying the
  Spalart–Shur correction term,"* ASME Journal of Turbomachinery
  **131**, 041010. — SST-RC.
- Wilcox, D.C. (2006) *Turbulence Modeling for CFD*, 3rd ed., DCW
  Industries. — Frame indifference and rotating-frame closures, §2.4.

## Viscous Fluxes

The stress tensor in the rotating frame uses the relative velocity:

$$\boldsymbol{\tau} = \mu \left( \nabla\boldsymbol{w} + (\nabla\boldsymbol{w})^T - \frac{2}{3}(\nabla \cdot \boldsymbol{w}) \mathbf{I} \right)$$

The viscous energy flux is computed as
$\boldsymbol{\tau} \cdot \boldsymbol{w}$, and the heat conduction flux
remains $\boldsymbol{q} = -k \nabla T$.

Since all quantities already use the relative velocity stored in the
solution, **no modifications are required** to standard diffusive flux
routines.

## Boundary Conditions

### Standard Boundaries

Inflow/outflow and other standard boundary conditions are applied
directly to the relative velocity field $\boldsymbol{w}$ without
modification. When defining the inflow properties it is necessary to
specify them in the **laboratory frame of reference**. The correction
in the rotating frame is performed automatically by
`RF_Convert_BC_Inflow`.

### Wall Boundaries

On rotating walls, the no-slip condition applies to the **absolute
velocity**.  When the wall rotates with the frame at angular velocity
$\boldsymbol{\Omega}$, the relative velocity of fluid at a wall point
is:

$$\boldsymbol{w}_{\text{wall}} = \boldsymbol{u}_{\text{abs, wall}} - (\boldsymbol{\Omega} \times \boldsymbol{r}_{\text{wall}})$$

For a perfectly rotating solid body
$\boldsymbol{u}_{\text{abs, wall}} = \boldsymbol{\Omega} \times \boldsymbol{r}_{\text{wall}}$,
hence:

$$\boldsymbol{w}_{\text{wall}} = \boldsymbol{0}$$

(i.e., zero relative velocity at a wall moving with the frame).

## Physical Interpretation

The rotating-frame formulation is equivalent to solving the full
inertial-frame problem but with fictitious body forces representing:

1. **Coriolis acceleration** — deflects moving fluid perpendicular to
   both the rotation axis and the velocity vector.
2. **Centrifugal acceleration** — acts radially outward from the
   rotation axis.

For steady-state turbomachinery flows, a steady rotating-frame
solution avoids the need to track time-periodic structures in the
inertial frame, enabling efficient steady-state analysis.

## Implementation Notes

- The computational domain is assumed to rotate as a **rigid body** at
  constant $\boldsymbol{\Omega}$.
- In Q2D the axis is fixed to $\boldsymbol{e} = [0, 0, 1]$
  (perpendicular to the simulation plane); only $\Omega = |\boldsymbol{\Omega}|$
  is configurable.
- The origin of the coordinate system (where $\boldsymbol{r} = \boldsymbol{0}$)
  should coincide with the axis of rotation.
- All flux routines (convective and viscous) operate on the relative
  velocity without modification.
- The rotating-frame source terms are added cell-wise via the
  `obj_source%models(:)` registry mechanism (caller-side dispatch in
  `Mod_Explicit`).
