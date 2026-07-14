# 2D Axisymmetric Derivation

This page derives the geometric source terms of the two-dimensional
axisymmetric regime from the three-dimensional Navier–Stokes equations
written in cylindrical coordinates $(x, r, \theta)$, reduced under
azimuthal independence and zero swirl.  It gives the full strain, stress
and vorticity tensors, the cylindrical-divergence decomposition that
isolates the source, the complete set of axisymmetric conservation
equations with their right-hand-side sources, a proof that the metric
pressure contribution cancels in the radial momentum equation, and the
treatment of the $1/r$ factor on the axis.

---

## 1. Three-dimensional Navier–Stokes in cylindrical coordinates

The compressible multi-species Navier–Stokes system in conservative form
is

$$
\frac{\partial \mathbf{U}}{\partial t} + \nabla\!\cdot\!\mathbf{F} = 0,
\qquad
\mathbf{U} = \bigl(\rho,\; \rho u_x,\; \rho u_r,\; \rho E_0,\; \rho_s,\; \rho_\phi\bigr)^\mathsf{T},
$$

with $\rho_s = \rho Y_s$ the species partial densities, $\rho_\phi$ the
passive-scalar / turbulence partial densities, $E_0$ the total specific
energy, and $\mathbf{F}$ the sum of convective and viscous fluxes.  The
velocity is $\mathbf{V} = (u_x, u_r, u_\theta)$.

The Newtonian strain-rate tensor $S_{ij} = \tfrac12(\partial_i V_j + \partial_j V_i)$
has, in the cylindrical basis,

$$
\begin{aligned}
S_{xx} &= \frac{\partial u_x}{\partial x}, &
S_{rr} &= \frac{\partial u_r}{\partial r}, &
S_{\theta\theta} &= \frac{u_r}{r} + \frac{1}{r}\frac{\partial u_\theta}{\partial\theta}, \\[4pt]
S_{xr} &= \tfrac12\!\left(\frac{\partial u_x}{\partial r} + \frac{\partial u_r}{\partial x}\right), &
S_{x\theta} &= \tfrac12\!\left(\frac{1}{r}\frac{\partial u_x}{\partial\theta} + \frac{\partial u_\theta}{\partial x}\right), &
S_{r\theta} &= \tfrac12\!\left(\frac{1}{r}\frac{\partial u_r}{\partial\theta} + \frac{\partial u_\theta}{\partial r} - \frac{u_\theta}{r}\right).
\end{aligned}
$$

The metric terms $u_r/r$ in $S_{\theta\theta}$ and $-u_\theta/r$ in
$S_{r\theta}$ arise from the curvature of the cylindrical basis vectors,
not from spatial velocity gradients.  The constitutive law with the
Stokes hypothesis $\lambda = -\tfrac23\mu$ is
$\tau_{ij} = 2\mu\,S_{ij} + \lambda\,(\nabla\!\cdot\!\mathbf{V})\,\delta_{ij}$.

---

## 2. Reduction under $\partial/\partial\theta \equiv 0$ and $u_\theta \equiv 0$

The axisymmetric regime models a flow **independent of the azimuthal
angle** and with **no swirl**, $\partial/\partial\theta \equiv 0$,
$u_\theta \equiv 0$.  The strain tensor collapses to a planar block plus
one surviving metric component:

$$
S_{xx} = \frac{\partial u_x}{\partial x}, \quad
S_{rr} = \frac{\partial u_r}{\partial r}, \quad
\boxed{\,S_{\theta\theta} = \frac{u_r}{r}\,}, \quad
S_{xr} = \tfrac12\!\left(\frac{\partial u_x}{\partial r} + \frac{\partial u_r}{\partial x}\right),
$$

with $S_{x\theta} = S_{r\theta} = 0$.  The $\theta\theta$ component is
**non-zero even at uniform radial velocity**: a purely radial motion
stretches every annular fluid ring.  The cylindrical divergence is

$$
\nabla\!\cdot\!\mathbf{V}
= S_{xx} + S_{rr} + S_{\theta\theta}
= \frac{\partial u_x}{\partial x} + \frac{\partial u_r}{\partial r} + \frac{u_r}{r}.
$$

The stress components follow with the cylindrical divergence inside the
dilatational term:

$$
\begin{aligned}
\tau_{xx} &= 2\mu\,\frac{\partial u_x}{\partial x} + \lambda\,(\nabla\!\cdot\!\mathbf{V}), &
\tau_{rr} &= 2\mu\,\frac{\partial u_r}{\partial r} + \lambda\,(\nabla\!\cdot\!\mathbf{V}), \\[4pt]
\tau_{\theta\theta} &= 2\mu\,\frac{u_r}{r} + \lambda\,(\nabla\!\cdot\!\mathbf{V}), &
\tau_{xr} &= \mu\!\left(\frac{\partial u_x}{\partial r} + \frac{\partial u_r}{\partial x}\right),
\end{aligned}
$$

with $\lambda = -\tfrac23\mu$ and $\tau_{x\theta} = \tau_{r\theta} = 0$.
Of the vorticity, only the azimuthal component
$\omega_\theta = \partial_x u_r - \partial_r u_x$ survives (the in-plane
shear); $\omega_x$ and $\omega_r$ vanish identically with $u_\theta$ and
$\partial_\theta$.

---

## 3. The cylindrical-divergence decomposition

The finite-volume solver discretises the conservation laws on the
meridional grid $(x, r)$ as if the divergence were planar.  For a vector
flux $\mathbf{F} = (F_x, F_r)$ with $\partial/\partial\theta \equiv 0$,

$$
\nabla\!\cdot\!\mathbf{F}
= \frac{\partial F_x}{\partial x} + \frac{1}{r}\frac{\partial(r\,F_r)}{\partial r}
= \underbrace{\frac{\partial F_x}{\partial x} + \frac{\partial F_r}{\partial r}}_{\text{planar residual}}
+ \frac{F_r}{r}.
$$

The planar reconstruction over the four cell edges captures the first
two terms; the third, $F_r/r$, has **no counterpart** and becomes a
geometric source on the right-hand side.  For the radial momentum
equation the source carries an additional metric term, derived in §5.2.

---

## 4. The complete set of axisymmetric conservation equations

Writing each equation as (planar flux divergence) $=$ (geometric source),
with the planar divergence $\mathcal{D}_{\text{2D}}(F) := \partial_x F_x + \partial_r F_r$:

$$
\begin{aligned}
\textbf{Continuity:}\quad &
\partial_t\rho + \mathcal{D}_{\text{2D}}(\rho\mathbf{V})
= -\frac{\rho\,u_r}{r}, \\[4pt]
\textbf{Species:}\quad &
\partial_t\rho_s + \mathcal{D}_{\text{2D}}(\rho_s\mathbf{V} + \mathbf{j}_s)
= -\frac{\rho_s\,u_r}{r} + \frac{\rho D_s\,\partial_r Y_s - V_c\,\rho_s}{r}, \\[4pt]
\textbf{$x$-momentum:}\quad &
\partial_t(\rho u_x) + \mathcal{D}_{\text{2D}}(\rho u_x\mathbf{V} + p\,\mathbf{e}_x - \boldsymbol\tau_x)
= -\frac{\rho u_x u_r}{r} + \frac{\tau_{xr}}{r}, \\[4pt]
\textbf{$r$-momentum:}\quad &
\partial_t(\rho u_r) + \mathcal{D}_{\text{2D}}(\rho u_r\mathbf{V} + p\,\mathbf{e}_r - \boldsymbol\tau_r)
= -\frac{\rho u_r^2}{r} + \frac{\tau_{rr} - \tau_{\theta\theta}}{r}, \\[4pt]
\textbf{Energy:}\quad &
\partial_t(\rho E_0) + \mathcal{D}_{\text{2D}}(\mathbf{F}^E)
= -\frac{(\rho E_0 + p)u_r}{r} + \frac{-q_r + u_x\tau_{xr} + u_r\tau_{rr} + (\mathcal{H}_d - h_{\text{mix}}V_c)\rho}{r}, \\[4pt]
\textbf{Scalar/RANS $\phi$:}\quad &
\partial_t\rho_\phi + \mathcal{D}_{\text{2D}}(\rho_\phi\mathbf{V} - D_\phi\nabla\phi)
= -\frac{\rho_\phi\,u_r}{r} + \frac{D_\phi\,\partial_r\phi}{r}.
\end{aligned}
$$

Here $\mathbf{j}_s = -\rho D_s\nabla Y_s + \rho_s\mathbf{V}_c$ is the
Hirschfelder-corrected species flux,
$q_r = -\kappa\,\partial_r T$ the conductive heat flux,
$\mathcal{H}_d = \sum_s D_s\,(\partial_r Y_s)\,h_s$ the diffusive
enthalpy flux, $h_{\text{mix}} = \sum_s Y_s h_s$ the mixture enthalpy,
$V_c = \sum_s D_s\,\partial_r Y_s$ the Hirschfelder correction velocity
(enforcing $\sum_s \mathbf{j}_s = 0$), and $D_\phi$ the effective
diffusion coefficient of each turbulence variable.  The two non-trivial
sources — the $r$-momentum metric term and the energy viscous work — are
derived in §5.

---

## 5. Assembly of the geometric sources

Collecting the right-hand sides above into an Eulerian (convective) part
and a viscous part:

$$
\mathbf{S}_{\text{Eul}}^{\text{axi}} = -\frac{1}{r}
\begin{pmatrix}
\rho_s\,u_r \\[2pt] \rho\,u_x\,u_r \\[2pt] \rho\,u_r^2 \\[2pt]
(\rho E_0 + p)\,u_r \\[2pt] \rho_\phi\,u_r
\end{pmatrix},
\qquad
\mathbf{S}_{\text{visc}}^{\text{axi}} = \frac{1}{r}
\begin{pmatrix}
\rho D_s\,\partial_r Y_s - V_c\,\rho_s \\[2pt]
\tau_{xr} \\[2pt]
\boxed{\,\tau_{rr} - \tau_{\theta\theta}\,} \\[2pt]
-q_r + u_x\tau_{xr} + u_r\tau_{rr} + (\mathcal{H}_d - h_{\text{mix}}V_c)\rho \\[2pt]
D_\phi\,\partial_r\phi
\end{pmatrix}.
$$

### 5.1 The scalar entries

For continuity, species, axial momentum and energy the radial flux is a
plain advective (plus diffusive) term, so the Eulerian source is simply
$-F_r/r$ and the viscous source is the corresponding $+ (\text{diffusive flux})/r$.
The species viscous entry $\rho D_s\,\partial_r Y_s - V_c\,\rho_s$ is the
radial component of $-\mathbf{j}_s$, so the source equals $-j_{s,r}/r$;
the energy viscous entry is the radial component of the viscous work and
heat/diffusion fluxes, $-(q_r - u_x\tau_{xr} - u_r\tau_{rr} + \cdots)/r$.

### 5.2 Radial momentum: the pressure-cancellation proof

The radial momentum balance is governed by the divergence of the
**momentum-flux tensor** $\Pi_{ij} = \rho V_i V_j + p\,\delta_{ij} - \tau_{ij}$,
whose radial cylindrical divergence carries the metric term
$-\Pi_{\theta\theta}/r$:

$$
[\nabla\!\cdot\!\boldsymbol{\Pi}]_r
= \frac{\partial \Pi_{xr}}{\partial x} + \frac{1}{r}\frac{\partial(r\,\Pi_{rr})}{\partial r} - \frac{\Pi_{\theta\theta}}{r}
= \underbrace{\frac{\partial \Pi_{xr}}{\partial x} + \frac{\partial \Pi_{rr}}{\partial r}}_{\text{planar residual}}
+ \frac{\Pi_{rr} - \Pi_{\theta\theta}}{r}.
$$

With $u_\theta \equiv 0$,
$\Pi_{rr} = \rho u_r^2 + p - \tau_{rr}$ and
$\Pi_{\theta\theta} = p - \tau_{\theta\theta}$, so

$$
\Pi_{rr} - \Pi_{\theta\theta} = \rho u_r^2 - \tau_{rr} + \tau_{\theta\theta}.
$$

The thermodynamic pressure $p$ appears in **both** $\Pi_{rr}$ and
$\Pi_{\theta\theta}$ and **cancels identically**: the radial geometric
source carries no isolated pressure term.  Equivalently, the cylindrical
divergence of the isotropic part is
$[\nabla\!\cdot\!(p\mathbf I)]_r = (1/r)\partial_r(r p) - p/r = \partial_r p$,
the bare pressure gradient — the planar $\partial_r p$ already evaluated
by the solver, with no $1/r$ residue.  $\square$

The source therefore splits as $-\rho u_r^2/r$ (Eulerian) and
$+(\tau_{rr} - \tau_{\theta\theta})/r$ (viscous).  The viscous term is
**non-zero even in pure radial expansion**, because $\tau_{\theta\theta}$
carries the irreducible metric strain $2\mu\,u_r/r$.

### 5.3 Turbulence variables

The turbulence transport variables $\phi \in \{k, \omega, \tilde\nu\}$
follow the same passive-scalar pattern.  Their radial flux is
$F_r = \rho_\phi u_r - D_\phi\,\partial_r\phi$, so the geometric source
$-F_r/r$ splits into an Eulerian metric source $-\rho_\phi u_r/r$ and a
viscous geometric source $+D_\phi\,\partial_r\phi/r$ — the **same sign**
as the conductive ($-q_r/r = +\kappa\,\partial_r T/r$) and species
diffusive terms, since all share the $-F_r^{\text{diff}}/r$ structure
with $F_r^{\text{diff}}$ a down-gradient flux.  The effective diffusivity
$D_\phi$ is the model-specific value ($\mu_l + \sigma_k\mu_t$ for $k$,
$\mu_l + \sigma_\omega\mu_t$ for $\omega$,
$(\mu_l + \rho\tilde\nu)/\sigma$ for Spalart–Allmaras).  The vorticity
magnitude that drives the Spalart–Allmaras production reduces to the
azimuthal component $\omega_\theta$ and receives no extra metric
contribution.

---

## 6. The $1/r$ factor on the axis

The geometric sources carry a factor $1/r$.  By the primal/dual mesh
construction the cell centres lie at
$r_{\text{ave}} = \tfrac14\sum_{\text{nodes}} r_{\text{node}} > 0$,
strictly **off** the axis, so the cell-centred evaluation of $u_r/r$ and
$F_r/r$ is regular.  Two safeguards are nonetheless applied, with
distinct thresholds:

- **Geometric source.** The inverse radius is regularised once per grid
  as $r^{-1}_{\text{reg}} = 1/\max(r_{\text{ave}}, \epsilon_0)$ with
  $\epsilon_0 = 10^{-12}$, and used for every $1/r$ source term.
- **Strain $\theta\theta$ component.** The constructor applies an
  explicit de l'Hôpital fallback,

  $$
  S_{\theta\theta} =
  \begin{cases}
  \dfrac{u_r}{r_{\text{ave}}}, & r_{\text{ave}} > 10^{-6}, \\[6pt]
  \dfrac{\partial u_r}{\partial r}, & r_{\text{ave}} \le 10^{-6},
  \end{cases}
  $$

  consistent with $\displaystyle\lim_{r\to0}\frac{u_r}{r} = \frac{\partial u_r}{\partial r}$
  for a smooth flow with $u_r(x, 0) = 0$ (radial velocity odd about the
  axis).  The **stress tensor**, evaluated only at cell centres where
  $r_{\text{ave}} > 0$, uses $u_r/r_{\text{ave}}$ directly with **no
  fallback**.

!!! note "Two thresholds, two roles"

    $\epsilon_0 = 10^{-12}$ merely floors the precomputed inverse radius
    of the geometric source; it never activates on a physical grid.  The
    $10^{-6}$ de l'Hôpital branch is a defensive limit inside the strain
    constructor for the limiting case $r_{\text{ave}} \to 0$, where the
    $\theta\theta$ strain is the analytic derivative rather than the
    quotient.

---

## 7. Implementation reference

The Eulerian and viscous cylindrical sources are assembled per block by
`Source_Axi_Eul_Blk` and `Source_Axi_Visc_Blk`; the regularised inverse
radius is precomputed by `Precompute_Axisymmetric_Geometry`.  The
cylindrical stress tensor — in particular
$\tau_{\theta\theta} = 2\mu\,u_r/r + \lambda\,(\nabla\!\cdot\!\mathbf{V})$
— and the $\theta\theta$ strain with its de l'Hôpital branch are produced
by the same constitutive routines used for the planar face fluxes,
evaluated in their cylindrical branch.

---

## References

1. R. B. Bird, W. E. Stewart, E. N. Lightfoot, *Transport Phenomena*,
   2nd ed., Wiley, 2002 (chapter 3, cylindrical Newton–Stokes stress
   tensor).
2. J. D. Anderson, *Modern Compressible Flow with Historical
   Perspective*, 3rd ed., McGraw-Hill, 2003 (chapter 6).
3. F. M. White, *Viscous Fluid Flow*, 3rd ed., McGraw-Hill, 2006
   (curvilinear forms of the viscous stress tensor).
