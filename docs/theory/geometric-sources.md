# Geometric Source Terms

Q2D supports three geometric regimes — two-dimensional planar,
two-dimensional axisymmetric (cylindrical), and quasi-two-dimensional
area-variation (depth-averaging of a thin three-dimensional layer between
virtual walls).  All three share a common two-dimensional primitive
vector $(u, v)$ and identical convective fluxes in the $(x, y)$ plane.
The two non-planar regimes add **geometric source terms** that capture
the implicit out-of-plane physics — the divergence of the cylindrical
$r$-component and the depth-averaged wall trace, respectively.

This page is the **overview hub**: it states the common framework and the
essentials of each regime, and points to the full formal derivations.

- [2D Axisymmetric Derivation](axisymmetric-derivation.md) — the complete
  cylindrical Navier–Stokes reduction, the radial-momentum
  pressure-cancellation proof, and the axis treatment.
- [Q2D Area-Variation Derivation](area-variation-derivation.md) — the
  full depth-averaging derivation: hypotheses, asymptotic scaling with the
  order of every term, the algebraic properties of the mean operator, the
  transverse-velocity profile, the tensor machinery, and the
  step-by-step $z$-integration of each conservation law.
- [Q2D RANS Closures](q2d-rans-closures.md) — the URANS $k$–$\omega$ and
  Spalart–Allmaras equations in the area-variation regime.

---

## 1. Common framework

The cell-wise residual update in the two-dimensional planar regime reads

$$
V_i\,\frac{\mathrm{d}\mathbf{U}_i}{\mathrm{d}t}
= -\sum_{f=1}^{4}\bigl(\mathbf{F}\!\cdot\!\hat{\mathbf{n}}\bigr)_f\,A_f
$$

with $\mathbf{U} = (\rho, \rho u, \rho v, \rho E_0, \rho Y_s, \rho\phi)^\mathsf{T}$
the conservative state, $\mathbf{F}$ the inviscid + viscous flux, and $f$
running over the four edges of the quadrilateral cell.  In the non-planar
regimes the residual gains an additive source term

$$
V_i\,\frac{\mathrm{d}\mathbf{U}_i}{\mathrm{d}t}
= -\sum_{f=1}^{4}(\mathbf{F}\!\cdot\!\hat{\mathbf{n}})_f\,A_f
  + \mathbf{S}_\text{geom}(\mathbf{U}_i, \ldots)\,V_i
$$

carrying the out-of-plane physics no longer represented by face fluxes
after the dimensional reduction.  A single regime is active,
selected via the `active` parameter in `[Q2D-Sources]`
(`none | axisymmetric | area_variation`).  The two
regimes yield **structurally different** sources, with **opposite signs**
for the same physical flux:

| | 2D axisymmetric (residual $-F_r/r$) | Q2D area-variation (residue $+\mathbf{G}\!\cdot\!\langle\mathbf{F}\rangle$) |
|---|---|---|
| convective | $-\rho_\phi u_r/r$ | $+\rho_\phi(\mathbf{V}\!\cdot\!\mathbf{G})$ |
| diffusive ($\phi$) | $+D_\phi\,\partial_r\phi/r$ | $-D_\phi(\mathbf{G}\!\cdot\!\nabla\phi)$ |

The difference in sign is intrinsic: the axisymmetric source is the
missing cylindrical residual $-F_r/r$, whereas the area-variation source
is the algebraic residue $+\mathbf{G}\!\cdot\!\langle\mathbf{F}\rangle$ of
writing the globally conservative $h$-weighted form per unit area.

---

## 2. Two-dimensional axisymmetric regime

The axisymmetric regime models a fully three-dimensional flow independent
of the azimuthal angle ($\partial/\partial\theta \equiv 0$) and without
swirl ($u_\theta \equiv 0$).  The primitive vector $(u_x, u_r)$ maps onto
the meridional grid $(x, r)$.

**Cylindrical divergence.** For $\mathbf{F} = (F_x, F_r)$,

$$
\nabla\!\cdot\!\mathbf{F}
= \frac{\partial F_x}{\partial x} + \frac{1}{r}\frac{\partial(r F_r)}{\partial r}
= \underbrace{\frac{\partial F_x}{\partial x} + \frac{\partial F_r}{\partial r}}_{\text{planar (face fluxes)}}
+ \frac{F_r}{r}.
$$

The residual $F_r/r$ becomes the geometric source $-F_r/r$.  The radial
momentum equation is the only non-trivial case: the pressure cancels in
the metric residual $(\Pi_{rr} - \Pi_{\theta\theta})/r$, leaving an
Eulerian source $-\rho u_r^2/r$ and a viscous source
$(\tau_{rr} - \tau_{\theta\theta})/r$.

**$\theta\theta$ strain non-zero.** The cylindrical metric makes
$S_{\theta\theta} = u_r/r$ non-zero even at uniform radial velocity, so
$\tau_{\theta\theta} = 2\mu\,u_r/r + \lambda\,(\nabla\!\cdot\!\mathbf{V})$
is non-zero — the origin of the radial viscous source.

**Axis.** The $1/r$ factor is regular at cell centres (which sit at
$r_{\text{ave}} = \Delta r/2 > 0$).  The geometric source uses a floored
inverse radius $1/\max(r_{\text{ave}}, 10^{-12})$; the $\theta\theta$
**strain** constructor carries a de l'Hôpital fallback
$S_{\theta\theta} = \partial_r u_r$ for $r_{\text{ave}} \le 10^{-6}$,
while the cell-centred **stress** uses $u_r/r$ directly.

For the full cylindrical Navier–Stokes reduction, the pressure-cancellation
proof, the complete equation set, and the axis treatment, see the
[2D Axisymmetric Derivation](axisymmetric-derivation.md).

---

## 3. Quasi-two-dimensional area-variation regime

The area-variation regime models a thin three-dimensional layer bounded
by two **virtual walls** $z = \pm h(x, y)/2$ of variable local thickness,
obtained by depth-averaging the 3D conservation equations across $h$.

**Height-gradient coupling.** The fluid domain has $h \in C^2(\Omega_{\text{2D}})$,
$h_{\min} > 0$, and the height-gradient vector and its flux contraction
are

$$
\mathbf{G} := -\frac{\nabla h}{h},
\qquad
\mathbf{V}\!\cdot\!\mathbf{G} := u\,G_x + v\,G_y .
$$

The asymptotic parameter $\varepsilon := \sup|\nabla h|/L_\text{ref} \ll 1$
controls the thin-layer ratio.

**Depth-averaging operator.** For $f(x, y, z, t)$,
$\langle f\rangle_z := \tfrac1h\int_{-h/2}^{+h/2} f\,\mathrm{d}z$, a
linear operator acting on functions.

**Key kinematic identity.** Under wall impermeability the transverse
velocity satisfies, exactly,

$$
\boxed{\;\langle\partial_z w\rangle_z = -\,\mathbf{V}\!\cdot\!\mathbf{G}\;}
$$

(the wall traces of $w$ are fixed by impermeability, and the fundamental
theorem of calculus turns the $z$-integral into the jump).  This is the
**Definition A** out-of-plane strain $S_{33} = -\mathbf{V}\!\cdot\!\mathbf{G}$:
it propagates through the Newtonian stress to give
$T_{33} = -\tfrac43\mu(\mathbf{V}\!\cdot\!\mathbf{G}) - \tfrac23\mu\,\nabla_{\text{2D}}\!\cdot\!\mathbf{V}$,
the structural fingerprint of the depth-averaged tensor.

**Geometric sources.** Each depth-averaged conservation law gains a
convective source $+\rho_\phi(\mathbf{V}\!\cdot\!\mathbf{G})$ and, for the
viscous/diffusive equations, a source built from $\mathbf{G}$ acting on
the depth-averaged fluxes — the momentum sources
$G_x(T_{33} - T_{11}) - G_y T_{12}$ and
$-G_x T_{12} + G_y(T_{33} - T_{22})$, the energy source
$-\mathbf{V}\!\cdot\!(\boldsymbol\tau\!\cdot\!\mathbf{G}) + \mathbf{G}\!\cdot\!\mathbf{q} + \sum_\alpha h_\alpha(\mathbf{G}\!\cdot\!\mathbf{j}^{\text{eff}}_\alpha)$,
and the species source $\mathbf{G}\!\cdot\!\mathbf{j}^{\text{eff}}_\alpha$
with Hirschfelder closure $\sum_\alpha\mathbf{j}^{\text{eff}}_\alpha = 0$.

!!! note "Definition A is mandatory"

    $S_{33} := \langle\partial_z w\rangle_z = -\mathbf{V}\!\cdot\!\mathbf{G}$
    (Definition A) is the only reduction that preserves momentum and
    energy conservation under depth-averaging.  The alternative
    $S_{33} = 0$ introduces an $O(\mathbf{V}\!\cdot\!\mathbf{G}) = O(U/L)$
    error in $T_{33}$ — comparable to the retained terms, **not**
    subleading.

For the hypotheses, the asymptotic scaling with the order of every term,
the algebraic properties of the mean operator, the transverse-velocity
profile, the complete tensor machinery, and the step-by-step
$z$-integration, see the
[Q2D Area-Variation Derivation](area-variation-derivation.md).  For the
turbulence closures, see [Q2D RANS Closures](q2d-rans-closures.md).

---

## 4. Implementation reference

The two geometric regimes are dispatched through the source registry at
residual-assembly time; at most one is active per block, and the planar
regime sets $\mathbf{S}_\text{geom} \equiv 0$.  The depth-averaging
coupling $\mathbf{V}\!\cdot\!\mathbf{G}$ is computed at faces and
propagated to the stress-tensor constructor through an optional argument,
so the same constitutive routine serves both face fluxes and geometric
sources.

---

## References

1. R. B. Bird, W. E. Stewart, E. N. Lightfoot, *Transport Phenomena*,
   2nd ed., Wiley, 2002 (chapters 1–3).
2. J. D. Anderson, *Modern Compressible Flow with Historical
   Perspective*, 3rd ed., McGraw-Hill, 2003.
3. A. H. Shapiro, *The Dynamics and Thermodynamics of Compressible Fluid
   Flow*, Vol. I, Ronald Press, 1953.
4. J.-F. Gerbeau, B. Perthame, "Derivation of viscous Saint-Venant system
   for laminar shallow water; numerical validation", *Discrete and
   Continuous Dynamical Systems — Series B*, 1(1):89–102, 2001.
