# Q2D RANS Closures

This page gives the RANS closures of the Q2D geometric regimes in their
**consolidated, sign-consistent form**.  Sections 1–5 cover the
eddy-viscosity closures of the Q2D area-variation regime — the Boussinesq
two-equation stress tensor with the turbulent kinetic energy in the
diagonal, the Menter SST $k$–$\omega$ (2003m) equations, the
Spalart–Allmaras model, and the energy-coupling term — each inheriting the
depth-averaging geometric sources of the
[area-variation derivation](area-variation-derivation.md): a convective
source $\rho\phi(\mathbf{V}\!\cdot\!\mathbf{G})$ and a diffusive
geometric source $-D_\phi(\mathbf{V}\!\cdot\!\mathbf{G}\text{-projection})$,
with the **Anderson-coherent negative sign** that the $+\mathbf{G}\!\cdot\!\langle\mathbf{F}\rangle$
structure of the area-variation regime dictates, and the production
correction $2(\mathbf{V}\!\cdot\!\mathbf{G})^2$.

Section 6 gives the **full second-moment closure** — the SSG/LRR-$\omega$
Reynolds-stress model (RSM), which transports the Reynolds-stress tensor
$\rho R_{ij}$ itself (no Boussinesq hypothesis) — reduced to **all three**
Q2D regimes (2D-plane, 2D-axi, Q2D area-variation).  Its rigorous
top-down derivation (coordinate-free model → cylindrical projection /
$z$-averaging → the $7\!\to\!5$ reduction) is given in
[Appendix A](#appendix-a-coordinate-free-rsm-derivation),
and is confirmed symbolically term-by-term by a SymPy reduction
audit (internal verification artefact, referenced below as §S$n$).

Throughout, the depth-averaged primitives are written without an
overbar, $\mu_{\text{eff}} = \mu_l + \mu_t$, the planar divergence is
$\nabla_{\text{2D}}\!\cdot\!\mathbf{V} = \partial_x u + \partial_y v$, and
the geometric contraction is $\mathbf{V}\!\cdot\!\mathbf{G} = uG_x + vG_y$.

---

## 1. Boussinesq two-equation stress tensor

The turbulent stress is the Boussinesq closure added to the
depth-averaged Newtonian tensor (Theorem 2 of the area-variation
derivation), with the eddy viscosity $\mu_t$ replacing $\mu$ and the
isotropic turbulent-kinetic-energy term $-\tfrac23\rho k\,\delta_{ij}$
on the diagonal:

$$
\begin{aligned}
\tau_{\text{tot},11} &= 2\mu_{\text{eff}}\,\partial_x u - \tfrac23\mu_{\text{eff}}\bigl(\nabla_{\text{2D}}\!\cdot\!\mathbf{V} - \mathbf{V}\!\cdot\!\mathbf{G}\bigr) - \tfrac23\rho k, \\
\tau_{\text{tot},22} &= 2\mu_{\text{eff}}\,\partial_y v - \tfrac23\mu_{\text{eff}}\bigl(\nabla_{\text{2D}}\!\cdot\!\mathbf{V} - \mathbf{V}\!\cdot\!\mathbf{G}\bigr) - \tfrac23\rho k, \\
\tau_{\text{tot},12} &= \mu_{\text{eff}}\bigl(\partial_y u + \partial_x v\bigr), \\
\tau_{\text{tot},33} &= -\tfrac43\mu_{\text{eff}}\,(\mathbf{V}\!\cdot\!\mathbf{G}) - \tfrac23\mu_{\text{eff}}\,\nabla_{\text{2D}}\!\cdot\!\mathbf{V} - \tfrac23\rho k.
\end{aligned}
$$

The out-of-plane component $\tau_{\text{tot},33}$ retains the
$-\tfrac43\mu_{\text{eff}}(\mathbf{V}\!\cdot\!\mathbf{G})$ depth-averaging
fingerprint and the $-\tfrac23\rho k$ diagonal term.  This diagonal
$\rho k$ is the device that makes the energy equation consistent without
an explicit effective pressure: the combination of the convective source
$(\rho E_{\text{res}} + p)(\mathbf{V}\!\cdot\!\mathbf{G})$ and the viscous
geometric source $-\mathbf{V}\!\cdot\!(\boldsymbol\tau_{\text{tot}}\!\cdot\!\mathbf{G})$
produces, through the diagonal $-\tfrac23\rho k$, exactly
$+\tfrac23\rho k\,(\mathbf{V}\!\cdot\!\mathbf{G})$ — equivalent to an
effective pressure $P^\ast = p + \tfrac23\rho k$ **without** writing
$P^\ast$, and with no double-counting.

!!! note "Why $p$ and not $P^\ast$"

    Writing the thermodynamic pressure $p$ in the convective fluxes and
    carrying $-\tfrac23\rho k$ only in the diagonal of
    $\boldsymbol\tau_{\text{tot}}$ is the convention that avoids counting
    $\tfrac23\rho k\,(\mathbf{V}\!\cdot\!\mathbf{G})$ twice — once in a
    $P^\ast$ convective source and once in the viscous geometric source.
    The effective-pressure contribution emerges once, from the viscous
    geometric term alone.

The corresponding axisymmetric tensor replaces
$-\mathbf{V}\!\cdot\!\mathbf{G}$ in the dilatational term with the
cylindrical metric strain $u_r/r$ and sets
$\tau_{\theta\theta} = 2\mu_{\text{eff}}\,u_r/r + \lambda\,(\nabla\!\cdot\!\mathbf{V}) - \tfrac23\rho k$
(see the [axisymmetric derivation](axisymmetric-derivation.md)).

---

## 2. Menter SST $k$–$\omega$ (2003m)

### 2.1 Turbulent kinetic energy $k$

$$
\boxed{\;
\partial_t(\rho k) + \partial_x(\rho u k - D_k\,\partial_x k) + \partial_y(\rho v k - D_k\,\partial_y k)
= \rho k(\mathbf{V}\!\cdot\!\mathbf{G}) - D_k(\mathbf{G}\!\cdot\!\nabla k) + P_k^{\text{Q2D}} - \beta^\ast\rho\omega k
\;}
$$

with $D_k = \mu_l + \sigma_k\mu_t$.  The two area-variation sources are
the convective $\rho k(\mathbf{V}\!\cdot\!\mathbf{G})$ and the diffusive
geometric $-D_k(\mathbf{G}\!\cdot\!\nabla k)$.  The **negative sign** of
the diffusive geometric source is the leading-order consequence of the
$+\mathbf{G}\!\cdot\!\langle\mathbf{F}^k\rangle$ structure with the
diffusive flux $\mathbf{F}^k_{\text{diff}} = -D_k\nabla k$, so
$\mathbf{G}\!\cdot\!\mathbf{F}^k_{\text{diff}} = -D_k(\mathbf{G}\!\cdot\!\nabla k)$.

### 2.2 Production with the area-variation correction

$$
\boxed{\;
P_k^{\text{Q2D}} = \mu_t\bigl[\,|S|^2_{\text{2D}} + 2(\mathbf{V}\!\cdot\!\mathbf{G})^2\bigr],
\qquad
|S|^2_{\text{2D}} = 2\bigl[(\partial_x u)^2 + (\partial_y v)^2\bigr] + (\partial_y u + \partial_x v)^2.
\;}
$$

The term $2(\mathbf{V}\!\cdot\!\mathbf{G})^2$ is the depth-average of the
out-of-plane strain contribution $2S_{33}^2 = 2(\mathbf{V}\!\cdot\!\mathbf{G})^2$
to $|S|^2$, and is **leading-order** — not a subleading correction:
$\langle|S|^2\rangle_z = |S|^2_{\text{2D}} + 2(\mathbf{V}\!\cdot\!\mathbf{G})^2 + O(\varepsilon^2)$.
The production is limited by $P_k \leftarrow \min(P_k, 10\,\beta^\ast\rho\omega k)$.

### 2.3 Specific dissipation rate $\omega$

$$
\boxed{\;
\partial_t(\rho\omega) + \partial_x(\rho u\omega - D_\omega\partial_x\omega) + \partial_y(\rho v\omega - D_\omega\partial_y\omega)
= \rho\omega(\mathbf{V}\!\cdot\!\mathbf{G}) - D_\omega(\mathbf{G}\!\cdot\!\nabla\omega) + \frac{\gamma\rho}{\mu_t}P_k^{\text{Q2D}} - \beta\rho\omega^2 + \langle\mathrm{CD}_{k\omega}\rangle
\;}
$$

with $D_\omega = \mu_l + \sigma_\omega\mu_t$ and the cross-diffusion term

$$
\langle\mathrm{CD}_{k\omega}\rangle = 2\rho\sigma_{\omega2}(1 - F_1)\,\frac{1}{\omega}\bigl(\partial_x k\,\partial_x\omega + \partial_y k\,\partial_y\omega\bigr).
$$

The cross-diffusion gradient product $\nabla k\!\cdot\!\nabla\omega$ is a
scalar invariant and carries no geometric metric term.  The blending
functions $F_1, F_2$ and the blended constants
($\sigma_k = \sigma_{k1}F_1 + \sigma_{k2}(1-F_1)$, etc.) are scalar
functions of scalar arguments and are unaffected by the geometric
reduction.

### 2.4 Wall boundary condition

At a solid wall the standard Menter condition applies,
$\omega_{\text{wall}} = \dfrac{60\,\nu}{\beta_1\,y^2}$, with $y$ the
distance to the first cell centre; $k_{\text{wall}} = 0$.

---

## 3. Spalart–Allmaras

The single transport equation for the modified eddy viscosity
$\tilde\nu$ is

$$
\boxed{\;
\partial_t(\rho\tilde\nu) + \partial_x\!\Bigl(\rho u\tilde\nu - \tfrac{1}{\sigma}(\mu_l + \rho\tilde\nu)\partial_x\tilde\nu\Bigr) + \partial_y\!\Bigl(\rho v\tilde\nu - \tfrac{1}{\sigma}(\mu_l + \rho\tilde\nu)\partial_y\tilde\nu\Bigr)
= \rho\tilde\nu(\mathbf{V}\!\cdot\!\mathbf{G}) - \tfrac{1}{\sigma}(\mu_l + \rho\tilde\nu)(\mathbf{G}\!\cdot\!\nabla\tilde\nu) + P_{\tilde\nu} + D^{\,\text{cb2}}_{\tilde\nu} - Y_{\tilde\nu}
\;}
$$

with production $P_{\tilde\nu} = c_{b1}\tilde{S}\rho\tilde\nu$, the
non-conservative diffusion $D^{\,\text{cb2}}_{\tilde\nu} = (c_{b2}/\sigma)\rho|\nabla\tilde\nu|^2$,
and destruction $Y_{\tilde\nu} = c_{w1}f_w\rho(\tilde\nu/d)^2$ ($d$ the
wall distance).  The diffusive geometric source
$-\tfrac{1}{\sigma}(\mu_l + \rho\tilde\nu)(\mathbf{G}\!\cdot\!\nabla\tilde\nu)$
carries the same **negative** Anderson-coherent sign as the SST sources.

The vorticity magnitude that drives $\tilde{S}$ reduces to the in-plane
component and receives **no** leading-order area-variation correction,
$\langle|\boldsymbol\omega|^2\rangle_z = \omega_3^2 + O(\varepsilon^2)$;
$|\nabla\tilde\nu|^2$ is the scalar invariant
$(\partial_x\tilde\nu)^2 + (\partial_y\tilde\nu)^2$.  The model variants
are scalar modifications of $\tilde{S}$ that compose with the geometric
sources unchanged:

- **SA-RC** (rotation/curvature, Spalart–Shur): multiplies the production
  by $f_{r1}(r^\ast, \tilde r)$ with $r^\ast = S/\Omega$ a scalar ratio.
- **SAcomp** (compressibility, Paciorri–Sabetta): a scalar correction
  factor.
- **SA-R** (with curvature via the strain magnitude): uses
  $|S|^2 = |S|^2_{\text{2D}} + 2(\mathbf{V}\!\cdot\!\mathbf{G})^2$ and so
  automatically inherits the area-variation production correction.

---

## 4. Energy coupling $\Phi_k$ (Strategy B)

With the total energy defined **without** the turbulent kinetic energy,
$E_{\text{res}} = e + \tfrac12|\mathbf{V}|^2$, the energy equation
(area-variation derivation §8.5) carries the coupling source

$$
\Phi_k = \bigl(-P_k + \beta^\ast\rho\omega k\bigr)\cdot\mathbb{1}_{\text{k-coupling}},
$$

active only when the $k$-energy coupling is enabled.  It returns to the
mean-flow energy the production extracted by, and the dissipation
returned from, the turbulence field, so that the total (mean +
turbulent) energy budget is consistent under the $E_{\text{res}}$
convention.

---

## 5. Boundary conditions

- **Walls:** $k = 0$, $\omega_{\text{wall}} = 60\nu/(\beta_1 y^2)$ (SST);
  $\tilde\nu = 0$ (SA).  Under the free-slip virtual walls of the
  area-variation regime, the turbulence diffusive flux normal to the wall
  vanishes at leading order, which is what makes the wall traction of the
  $k$, $\omega$, $\tilde\nu$ equations zero (area-variation derivation
  §8, master identity).
- **Inflow:** $k$, $\omega$ (or $\tilde\nu$) from a prescribed turbulence
  intensity and length scale.
- **Outflow:** zero-gradient extrapolation.

---

## 6. Reynolds-stress transport (SSG/LRR-$\omega$ RSM)

The SSG/LRR-$\omega$ model (SSG/LRR-RSM-w2019, NASA Turbulence Modeling
Resource) is a **full second-moment closure**: rather than modelling the
Reynolds stress through an eddy viscosity, it transports the six
independent components of the symmetric stress tensor
$\rho R_{ij} = \rho\langle u_i' u_j'\rangle$ together with the specific
dissipation $\rho\omega$.  There is **no Boussinesq hypothesis** — the
transported $\rho R_{ij}$ enter the mean-momentum flux directly (§6.8).
Through the SST blending function $F_1$ the model interpolates a near-wall
LRR-$\omega$ constant set ($\cdot_w$) and a free-stream SSG-$\varepsilon$
set ($\cdot_e$); $C_\mu = 0.09$ and the turbulent kinetic energy is
$k = \tfrac12(R_{11}+R_{22}+R_{33})$.

The results below are stated in their consolidated per-regime form; every
curvature/dilatation and geometric-source term is **derived top-down** in
[Appendix A](#appendix-a-coordinate-free-rsm-derivation)
and confirmed symbolically by the SymPy reduction audit (§S$n$).

### 6.1 Transported variables and the unified regime structure

In 3D the model carries seven variables
$(\rho R_{11},\rho R_{22},\rho R_{33},\rho R_{12},\rho R_{13},\rho R_{23},\rho\omega)$.
Each Q2D regime **drops the two cross-stresses that vanish under a mirror
symmetry of the regime** (Appendix A.6), leaving **five** transported
variables $(\rho R_{11},\rho R_{22},\rho R_{33},\rho R_{12},\rho\omega)$ —
"four-stress $+\,\omega$".  The component map:

| Q2D variable | 2D-plane | 2D-axi (no swirl) | area-variation |
|---|---|---|---|
| $R_{11}$ | $\langle u'u'\rangle$ | $\langle u_x' u_x'\rangle$ (axial) | $\langle u'u'\rangle$ |
| $R_{22}$ | $\langle v'v'\rangle$ | $\langle u_r' u_r'\rangle$ (radial) | $\langle v'v'\rangle$ |
| $R_{33}$ | $\langle w'w'\rangle$ (spanwise) | $\langle u_\theta' u_\theta'\rangle$ (azimuthal) | $\langle w'w'\rangle$ (depth) |
| $R_{12}$ | $\langle u'v'\rangle$ | $\langle u_x' u_r'\rangle$ | $\langle u'v'\rangle$ |
| dropped | $R_{13},R_{23}$ ($\partial_z\!\equiv\!0$) | $R_{x\theta},R_{r\theta}$ ($\theta$-reflection) | $R_{13},R_{23}$ ($z$-parity) |

$R_{33}$ is **retained in every regime**: it is a non-zero normal stress
fed by pressure-strain redistribution even where its direct production
vanishes (hence "four-stress", not three).

All three regimes share one structure (Appendix A; sympy §S2–S6): the
**same** SSG+LRR pressure-strain, dissipation and realizability *3×3*
tensor algebra (§6.3); the **same** exact production (§6.2) fed a
*block-diagonal* mean-velocity gradient whose only out-of-plane entry is
the regime dilatation

$$
(\nabla\mathbf{U})_{33} = S_{33} =
\begin{cases}
0, & \text{2D-plane},\\[2pt]
u_r/r, & \text{2D-axi},\\[2pt]
-\,\mathbf{V}\!\cdot\!\mathbf{G}, & \text{area-variation},
\end{cases}
\qquad (\nabla\mathbf{U})_{3k}=(\nabla\mathbf{U})_{k3}=0\ (k\neq3);
$$

and geometric transport/diffusion sources **identical in form to the
scalar** mean-flow RANS sources of each regime (the
[axisymmetric](axisymmetric-derivation.md) $-F_r/r$ and the
[area-variation](area-variation-derivation.md) $+\mathbf{G}\!\cdot\!\langle\mathbf{F}\rangle$
structures), now applied per stress component $\rho R_{ij}$:

| regime | convective geometric source (per component) | diffusive geometric source |
|---|---|---|
| 2D-plane | $0$ | $0$ |
| 2D-axi | $-\,\rho\,u_r\,R_{ij}/r$ | $+\,\mathsf{F}^{(ij)}_r/r$ |
| area-variation | $+\,\rho\,R_{ij}\,(\mathbf{V}\!\cdot\!\mathbf{G})$ | $-\,\mathsf{F}^{(ij)}\!\cdot\!\mathbf{G}$ |

with $\mathsf{F}^{(ij)}$ the GGDH diffusive flux of $\rho R_{ij}$ (§6.4).
The out-of-plane strain $S_{33}$ is exactly the value already produced by
`Strain_Tensor_Axi` ($u_r/r$, de l'Hôpital at the axis) and
`Strain_Tensor_Std` ($-\mathbf{V}\!\cdot\!\mathbf{G}$, Definition A) — so
no widening of the stored $2\times2$ velocity gradient is required
(Appendix A.1).

### 6.2 Production (exact, with the curvature/dilatation term)

The production requires no modelling:

$$
\boxed{\;P_{ij} = -\rho\bigl(R_{ik}\,\partial_k U_j + R_{jk}\,\partial_k U_i\bigr)\;}
$$

evaluated with the block-diagonal regime gradient.  The out-of-plane
stress therefore acquires the **curvature/dilatation production** (sympy
§S3)

$$
\boxed{\;P_{33} = -2\rho R_{33}\,(\nabla\mathbf{U})_{33} =
\begin{cases}
0, & \text{plane},\\[2pt]
-2\rho R_{33}\,u_r/r, & \text{2D-axi (azimuthal stretching)},\\[2pt]
+2\rho R_{33}\,(\mathbf{V}\!\cdot\!\mathbf{G}), & \text{area-variation (depth compression)}.
\end{cases}\;}
$$

The in-plane productions reduce to the planar $2\times2$ block,
$P_{11} = -2\rho(R_{11}\partial_x U + R_{12}\partial_y U)$,
$P_{22} = -2\rho(R_{22}\partial_y V + R_{12}\partial_x V)$,
$P_{12} = -\rho(R_{11}\partial_x V + R_{22}\partial_y U + R_{12}(\partial_x U+\partial_y V))$
(axial/radial indices in 2D-axi), and the trace
$P_{kk} = P_{11}+P_{22}+P_{33}$ — **including** the out-of-plane term —
drives the $\omega$-equation (§6.5).

### 6.3 Pressure-strain, dissipation, realizability (regime-invariant)

These are pure $3\times3$ tensor algebra in the anisotropy
$a_{ij} = R_{ij}/k - \tfrac23\delta_{ij}$, the strain $S_{ij}$, its
deviator $S^\ast_{ij} = S_{ij} - \tfrac13 S_{kk}\delta_{ij}$ and the
vorticity $W_{ij}$; they are **identical in all three regimes**, provided
the contracted traces $a_{kl}a_{kl}$, $a_{kl}S_{kl}$, $P_{kk}$ run over the
full $3\times3$ (the $33$-component included).  With $\rho\varepsilon = \rho C_\mu k\omega$,

$$
\Pi_{ij} = \underbrace{-\bigl(C_1\rho\varepsilon + \tfrac12 C_1^\ast P_{kk}\bigr)a_{ij}}_{C_1}
+ \underbrace{C_2\rho\varepsilon\bigl(a_{ik}a_{kj} - \tfrac13 a_{kl}a_{kl}\,\delta_{ij}\bigr)}_{C_2}
+ \underbrace{\bigl(C_3 - C_3^\ast\sqrt{a_{kl}a_{kl}}\bigr)\rho k\,S^\ast_{ij}}_{C_3}
$$
$$
+ \underbrace{C_4\rho k\bigl(a_{ik}S_{jk} + a_{jk}S_{ik} - \tfrac23 a_{kl}S_{kl}\,\delta_{ij}\bigr)}_{C_4}
+ \underbrace{C_5\rho k\bigl(a_{ik}W_{jk} + a_{jk}W_{ik}\bigr)}_{C_5},
$$

$$
\rho\varepsilon_{ij} = \tfrac23\rho\varepsilon\,\delta_{ij}
\quad(\text{isotropic dissipation}),
\qquad
\boxed{\;\text{Source}_{ij} = P_{ij} + \Pi_{ij} - \rho\varepsilon_{ij}\;}
$$

The blended constants $C_1\!=\!F_1C_{1w}+F_2C_{1e}$, … (LRR-$\omega$:
$C_{1w}\!=\!1.8$, $C_{1w}^\ast\!=\!0$, $C_{2w}\!=\!0$, $C_{3w}\!=\!0.8$,
$C_{3w}^\ast\!=\!0$, $C_{4w}\!=\!\tfrac12(18C_{2}^{\text{LRR}}+12)/11$,
$C_{5w}\!=\!\tfrac12(-14C_{2}^{\text{LRR}}+20)/11$ with $C_2^{\text{LRR}}\!=\!0.52$;
SSG-$\varepsilon$: $C_{1e}\!=\!1.7$, $C_{1e}^\ast\!=\!0.9$, $C_{2e}\!=\!1.05$,
$C_{3e}\!=\!0.8$, $C_{3e}^\ast\!=\!0.65$, $C_{4e}\!=\!0.625$, $C_{5e}\!=\!0.2$)
are scalar functions of scalar arguments — unaffected by the geometric
reduction.  The $C_5$ rapid term uses the **absolute** vorticity
$W_{ij}^{\text{abs}}$ (the relative vorticity plus $2\Omega_{\text{frame}}$
in a rotating frame); the rotating-frame 2D-axi extension is out of scope
(project decision), so in the base regimes $W_{ij}=W_{ij}^{\text{rel}}$
with the single surviving in-plane component
$\omega_3 = \partial_x V - \partial_y U$ (azimuthal in 2D-axi), exactly the
`Vorticity_Vector`/`Vorticity_Tensor` 2D reduction.

!!! note "Why the pressure-strain is regime-blind"

    $\Pi_{ij}$, $\rho\varepsilon_{ij}$ and the realizability limiter are
    **true tensor equations** — frame-independent algebraic functions of
    $a$, $S$, $W$.  Once those tensors carry the correct $33$-component per
    regime ($S_{33}=(\nabla\mathbf{U})_{33}$, $a_{33}=R_{33}/k-\tfrac23$,
    $W_{3\cdot}=0$), the same code evaluates them in every regime.  The
    sympy audit (§S2, §S6) confirms this by reducing the full 3D operator
    and matching the plane restriction term-by-term.

**Realizability** (applied after each update, regime-invariant): the
normal stresses are floored, $R_{\alpha\alpha} \ge 10^{-16}$ (no sum), and
the shear stresses are Schwarz-limited,
$|R_{\alpha\beta}| \le \sqrt{R_{\alpha\alpha}R_{\beta\beta}}$.

### 6.4 Diffusion: generalized gradient (Daly–Harlow GGDH)

The molecular-plus-turbulent diffusive flux of each stress component is
the anisotropic generalized-gradient (Daly–Harlow) form — the turbulent
part uses the **stress tensor itself** as the diffusivity:

$$
\mathsf{F}^{(ij)}_m = \mu_l\,\partial_m R_{ij}
+ D_\ast\,\rho R_{ml}\,\partial_l R_{ij},
\qquad
D_\ast = \frac{D_w F_1 + D_e F_2}{C_\mu\,\omega},
$$

with $D_w = 0.75\,C_\mu$, $D_e = 0.22$ (the shock-detector `-SD` variant
replaces the GGDH tensor by the isotropic
$(\mu_l + \tfrac{D_\ast}{C_\mu}\mu_t)\partial_m R_{ij}$ with
$\mu_t=\rho k/\omega$).  Its regime divergence produces the geometric
diffusive source of §6.1: $+\mathsf{F}^{(ij)}_r/r$ in 2D-axi (same
$-F_r/r$ structure as the scalar conductive/species fluxes,
[axisymmetric derivation §5.3](axisymmetric-derivation.md)) and
$-\mathsf{F}^{(ij)}\!\cdot\!\mathbf{G}$ in area-variation (the
**negative** Anderson-coherent sign of the
$+\mathbf{G}\!\cdot\!\langle\mathbf{F}\rangle$ structure,
[area-variation derivation §8.5](area-variation-derivation.md)),
precisely as the SST/SA diffusive geometric sources of §2–§3.

### 6.5 The $\omega$-equation

$\omega$ is a scalar specific-dissipation rate and transports exactly like
the SST $\omega$ of §2.3, with the same scalar geometric sources
($-\rho\omega u_r/r$ convective $+\,D_\omega\partial_r\omega/r$ diffusive in
2D-axi; $+\rho\omega(\mathbf{V}\!\cdot\!\mathbf{G}) - D_\omega(\mathbf{G}\!\cdot\!\nabla\omega)$
in area-variation):

$$
\boxed{\;\text{Source}_\omega
= \alpha_\omega\,\frac{\omega}{k}\cdot\tfrac12 P_{kk}
- \beta_\omega\,\rho\omega^2
+ \sigma_d\,\frac{\rho}{\omega}\max\!\bigl(\nabla k\!\cdot\!\nabla\omega,\,0\bigr)\;}
$$

with $D_\omega = \mu_l + \sigma_\omega\rho k/\omega$ and the blended
$\alpha_\omega,\beta_\omega,\sigma_\omega,\sigma_d$ (LRR-$\omega$:
$0.5556,\,0.075,\,0.5,\,0$; SSG-$\varepsilon$:
$0.44,\,0.0828,\,0.856,\,1.712$).  The production half-trace $\tfrac12 P_{kk}$
inherits the regime curvature/dilatation through $P_{33}$ (§6.2); the
cross-diffusion gradient product $\nabla k\!\cdot\!\nabla\omega$ is a scalar
invariant (in-plane $\partial_x k\,\partial_x\omega + \partial_y k\,\partial_y\omega$),
carrying no extra metric term.  The blending $F_1 = \tanh(\zeta^4)$ uses
the SST argument $\zeta = \min\!\bigl(\max(\sqrt k/(C_\mu\omega d),\,500\nu/(\rho\omega d^2)),\,4\sigma_{\omega e}\rho k/(\mathrm{CD}\,d^2)\bigr)$,
$\mathrm{CD} = \sigma_{de}(\rho/\omega)\max(\nabla k\!\cdot\!\nabla\omega,10^{-20})$,
$d$ the wall distance — all scalar, regime-blind.

### 6.6 The $7\to5$ reduction

$$
\boxed{\;
\text{With } R_{13}=R_{23}=0 \text{ and the block-diagonal regime gradient,}\;
\text{Source}_{13}\equiv\text{Source}_{23}\equiv 0.
\;}
$$

Production, pressure-strain and dissipation **never regenerate** the two
dropped cross-stresses, so the five-equation system is closed.  Proved in
Appendix A.6 from the regime mirror symmetry, and confirmed symbolically
(sympy §S2: $\text{Source}_{5}=\text{Source}_{6}=0$ identically).

### 6.7 Wall, axis and symmetry boundary conditions

- **Solid wall** (regime-invariant): $R_{ij}=0$ (no-slip kills all
  fluctuations) and the Menter smooth-wall value
  $\omega_{\text{wall}} = 800\,\nu/d^2$, $d$ the first cell-centre
  distance; the ghost enforces $R_{ij,\text{ghost}} = -R_{ij}$ (so the
  face value is zero) and extrapolates $\omega$.
- **Axis of symmetry** (2D-axi, $r\to0$): the axis is a reflection line,
  $r\to-r$, $u_r\to-u_r$, $u_x\to u_x$.  The normal stresses
  $R_{xx},R_{rr},R_{\theta\theta}$ are **even** (zero-gradient ghost); the
  shear $R_{xr}=\langle u_x'u_r'\rangle$ is **odd**, so $R_{xr}=0$ on the
  axis (sign-flipped ghost).  This is the RSM analogue of the mean-flow
  axis treatment; the $1/r$ source and the $\theta\theta$ strain use the
  same regularisation / de l'Hôpital branch as the mean flow (§6.1,
  [axisymmetric derivation §6](axisymmetric-derivation.md)).
- **Symmetry plane** (Cartesian): the off-diagonal stresses normal to the
  plane sign-flip, the diagonals take a zero-gradient ghost — the standard
  `RSM_Symmetry` rule.
- **Area-variation virtual walls**: the top/bottom $z = \pm h/2$ walls are
  free-slip and impermeable; they impose **no** $R_{ij}$ condition — their
  effect is entirely the geometric source of §6.1 (Appendix A.5), exactly
  as for the scalar turbulence variables (§5).

### 6.8 Mean-momentum coupling (no Boussinesq)

The Reynolds stress enters the mean-momentum and energy fluxes **directly**
— the total stress is the laminar Newtonian tensor (molecular viscosity
$\mu_l$ only) **minus** the transported Reynolds stress:

$$
\boxed{\;\tau^{\text{tot}}_{ij} = \tau^{\text{lam}}_{ij}(\mu_l) - \rho R_{ij}\;}
$$

i.e. $\tau^{\text{tot}}_{11} = \lambda_l(\nabla\!\cdot\!\mathbf{V}) + 2\mu_l\partial_x u - \rho R_{11}$,
$\tau^{\text{tot}}_{12} = \mu_l(\partial_y u + \partial_x v) - \rho R_{12}$,
and the out-of-plane $\tau^{\text{tot}}_{33} = \tau^{\text{lam}}_{33}(\mu_l) - \rho R_{33}$
with $\tau^{\text{lam}}_{33}$ the regime laminar value ($2\mu_l u_r/r + \lambda_l\nabla\!\cdot\!\mathbf{V}$
in 2D-axi; $-\tfrac43\mu_l(\mathbf{V}\!\cdot\!\mathbf{G}) - \tfrac23\mu_l\nabla_{\text{2D}}\!\cdot\!\mathbf{V}$
under Definition A).  There is **no** $\mu_t$, **no** $-\tfrac23\rho k$
Boussinesq diagonal: the actual $\rho R_{ij}$ — carrying the full
anisotropy — replace the modelled isotropic-plus-eddy stress of §1.  The
regime laminar tensor is the existing `Stress_Tensor_Std`/`_Axi` evaluated
with $\mu_t=0$, so the same constitutive path again serves the face fluxes
and the geometric sources (§6.1), now with the $-\rho R_{ij}$ subtraction.

---

## 7. Implementation reference

The Boussinesq two-equation tensor is assembled by
`Stress_Tensor_2eq_Std`/`_Axi` (the $-\tfrac23\rho k$ diagonal entering
through `rhokdelta`); the SST and SA transport terms and their diffusive
geometric sources are produced by `Lib_SST` and `Lib_Spalart` and routed
into the area-variation / axisymmetric source assembly through
`RANS_Diffusive_SourceQ2D`.  The eddy viscosity $\mu_t$ is propagated to
the stress constructor with the same depth-averaging coupling
$\mathbf{V}\!\cdot\!\mathbf{G}$ used for the laminar tensor, so a single
constitutive path serves the face fluxes and the geometric sources.

The SSG/LRR-$\omega$ RSM (§6) maps onto the same regime machinery.  Its
wiring mirrors `Lib_Fluid`'s `_Std`/`_Axi` procedure-pointer pattern:
`Lib_SSGLRR` provides the regime realizations supplying the
block-diagonal gradient (the out-of-plane $S_{33}$ from
`Strain_Tensor_Std`/`_Axi`), the regime-invariant
pressure-strain/dissipation/realizability algebra (§6.3), the GGDH
diffusive flux (§6.4) routed through the same per-regime divergence as the
viscous flux, the $\omega$-equation (§6.5), the wall/axis/symmetry BC
(§6.7), and the `Stress_Vector_RSM` momentum coupling (§6.8) subtracting
the retained $\rho R_{ij}$ from the laminar tensor.  The state layout adds
five variables at `nt` ($\rho R_{11},\rho R_{22},\rho R_{33},\rho R_{12},\rho\omega$,
$\texttt{nrans}=5$); `Mod_RANS` dispatches `SSGLRR` (and the
stress-dependent `SSGLRR-SD` variant) directly.  The model is verified
per regime: 2D-plane against the NASA Turbulence Modeling Resource, and
2D-axisymmetric and Q2D area-variation by per-regime wiring sweeps,
method-of-manufactured-solutions instruments, and running
verification-and-validation cases.

---

## Appendix A — Coordinate-free RSM derivation

This appendix derives the per-regime Reynolds-stress transport of §6
**top-down**, in the manner of the
[axisymmetric](axisymmetric-derivation.md) and
[area-variation](area-variation-derivation.md) derivations: start from the
general model in coordinate-free (basis-free) differential form, then
project into cylindrical coordinates (2D-axi) and depth-average
(area-variation); the 3D-Cartesian MOSE code is only the **cross-check of
the Cartesian/plane limit**, never the starting point.  The reduction is
verified term-by-term by the SymPy reduction audit
(referenced as §S$n$).

!!! note "Derivation panel (multi-perspective)"

    The derivation is cross-checked from four expert standpoints and an
    adversarial reviewer (the review arc is recorded in A.8): **(RSM
    theory)** the pressure-strain/$\omega$ closure and realizability;
    **(differential geometry)** the covariant transport of a rank-2 tensor
    in cylindrical and depth-averaged metrics; **(Q2D peculiarità)** the
    No-$w$ reduction, the geometric sources and the $\mathbf{V}\!\cdot\!\mathbf{G}$
    Definition A; **(numerics/V&V)** stiffness, BC and NASA-TMR validation.
    Every contested term is refuted, then confirmed symbolically.

### A.1 The coordinate-free transport model

The SSG/LRR-$\omega$ model transports the symmetric Reynolds-stress tensor
$\rho\mathbf{R}$ (components $\rho R_{ij}=\rho\langle u_i'u_j'\rangle$) and
the scalar $\rho\omega$.  In intrinsic (basis-free) form,

$$
\partial_t(\rho\mathbf{R}) + \nabla\!\cdot\!(\rho\mathbf{U}\otimes\mathbf{R})
= \mathbf{P} + \boldsymbol\Pi - \rho\boldsymbol\varepsilon + \nabla\!\cdot\!\mathsf{D},
$$

with the **production** $\mathbf{P} = -\rho\bigl(\mathbf{R}\!\cdot\!(\nabla\mathbf{U})^{\!\top} + (\nabla\mathbf{U})\!\cdot\!\mathbf{R}\bigr)$
(components $P_{ij} = -\rho(R_{ik}\partial_k U_j + R_{jk}\partial_k U_i)$),
the **pressure-strain** $\boldsymbol\Pi$ (SSG slow + LRR rapid, §6.3), the
**isotropic dissipation** $\rho\boldsymbol\varepsilon = \tfrac23\rho\varepsilon\,\mathbf{I}$
($\rho\varepsilon = \rho C_\mu k\omega$), and the **GGDH diffusion** flux
$\mathsf{D}$ (a rank-3 object, §6.4).  $\boldsymbol\Pi$, $\boldsymbol\varepsilon$
and the realizability limiter are **algebraic tensor functions** of
$\mathbf{a} = \mathbf{R}/k - \tfrac23\mathbf{I}$, the strain
$\mathbf{S} = \mathrm{sym}\,\nabla\mathbf{U}$ and the vorticity
$\mathbf{W} = \mathrm{skew}\,\nabla\mathbf{U}$ — hence frame-independent and
identical in every regime once $\nabla\mathbf{U}$ is the correct regime
tensor.

**The velocity-gradient embedding.**  In each Q2D regime the mean velocity
has no third (out-of-plane) component varying in-plane, and the mean fields
do not vary out-of-plane; the physical velocity-gradient tensor is
therefore **block-diagonal**,

$$
\nabla\mathbf{U} =
\begin{pmatrix} \partial_1 U_1 & \partial_2 U_1 & 0\\ \partial_1 U_2 & \partial_2 U_2 & 0\\ 0 & 0 & (\nabla\mathbf{U})_{33}\end{pmatrix},
\qquad
(\nabla\mathbf{U})_{33} = \begin{cases} 0 & \text{plane}\\ u_r/r & \text{axi}\\ -\mathbf{V}\!\cdot\!\mathbf{G} & \text{area-var}\end{cases}
$$

The lone out-of-plane entry is a *scalar* — exactly what `Strain_Tensor`'s
$E_{33}$ already carries — so the stored $2\times2$ gradient never needs
widening (Q2D peculiarità #1/#3).  This is the structural pivot: it makes
the *production* and *pressure-strain* regime-correct through one scalar,
and leaves only the *transport* (convection + diffusion) divergence to be
projected per regime (A.3, A.5).

### A.2 The transport divergence is the only genuinely geometric piece

Production, pressure-strain, dissipation and realizability are **pointwise
algebra** (no spatial derivatives of $\mathbf{R}$), so their regime form
follows entirely from the gradient embedding A.1 and the value of $S_{33}$
— nothing else to project.  The genuinely coordinate-dependent work is the
divergence of the two tensor fluxes: the convective $\rho\mathbf{U}\otimes\mathbf{R}$
and the diffusive $\mathsf{D}$.  We project these in A.3 (cylindrical) and
A.5 ($z$-average).

### A.3 2D-axi: cylindrical projection of the transport divergence

Coordinates $(x,r,\theta)$, physical components, scale factors $(1,1,r)$;
axisymmetry $\partial_\theta\equiv0$, no swirl $u_\theta\equiv0$.  The
convective flux is the rank-3 tensor $T_{kij} = \rho U_k R_{ij}$, contracted
on the first (flux) index $k$.  Writing $\mathbf{T} = T_{kij}\,\hat{\mathbf e}_k\hat{\mathbf e}_i\hat{\mathbf e}_j$
and using $\partial_\theta\hat{\mathbf e}_r = \hat{\mathbf e}_\theta$,
$\partial_\theta\hat{\mathbf e}_\theta = -\hat{\mathbf e}_r$,
$\partial_\theta\hat{\mathbf e}_x = 0$ (the only non-zero basis
derivatives),

$$
(\nabla\!\cdot\!\mathbf{T})_{ij}
= \sum_m \frac{1}{h_m}\,\hat{\mathbf e}_m\!\cdot\!\partial_{q_m}\!\bigl(T_{kij}\hat{\mathbf e}_k\hat{\mathbf e}_i\hat{\mathbf e}_j\bigr).
$$

Only the $m=\theta$ term carries basis derivatives (prefactor $1/r$); with
$\partial_\theta T_{kij}=0$ it splits into three pieces:

- **(A)** $\hat{\mathbf e}_\theta\!\cdot\!(\partial_\theta\hat{\mathbf e}_k)$
  picks $k=r$ (since $\hat{\mathbf e}_\theta\!\cdot\!\hat{\mathbf e}_\theta=1$):
  contributes $T_{rij}/r$, completing the annular flux $\tfrac1r\partial_r(rT_{rij})$.
- **(B), (C)** the basis derivatives on the *free* indices $i,j$ require
  $\hat{\mathbf e}_\theta\!\cdot\!\hat{\mathbf e}_k=\delta_{k\theta}$, i.e.
  the flux component $T_{\theta ij} = \rho\,u_\theta\,R_{ij} = 0$.  **They
  vanish identically — there is no swirl to advect the rotating basis.**

Hence

$$
\boxed{\;(\nabla\!\cdot\!\mathbf{T})_{ij}
= \partial_x(\rho u_x R_{ij}) + \frac1r\partial_r(r\,\rho u_r R_{ij})
= \underbrace{\partial_x(\rho u_x R_{ij}) + \partial_r(\rho u_r R_{ij})}_{\text{planar residual (FV)}}
+ \frac{\rho u_r R_{ij}}{r}\;}
$$

so the apparent source is the **pure scalar** $-\rho u_r R_{ij}/r$ per
component — with **no inter-component Christoffel coupling** (sympy §S4).

!!! note "Contrast with the momentum equation (the reviewer's objection)"

    One might expect a $-R_{\theta\theta}/r$-type apparent term in the
    $R_{rr}$ equation, by analogy with the $-\sigma_{\theta\theta}/r$ of
    the [radial momentum balance](axisymmetric-derivation.md#52-radial-momentum-the-pressure-cancellation-proof).
    It does **not** arise.  There, the stress $\boldsymbol\sigma$ is the
    *flux* — a rank-2 tensor with a non-zero $\sigma_{\theta\theta}$ in the
    flux slot — so its divergence keeps $-\sigma_{\theta\theta}/r$.  Here
    $\mathbf{R}$ is the *transported density* and the flux is
    $\rho\mathbf{U}\otimes\mathbf{R}$, whose $\theta$-component vanishes
    with $u_\theta$.  The sympy audit reproduces **both**: the RSM
    no-coupling result *and* the momentum $-\sigma_{\theta\theta}/r$
    contrast (§S4).

The diffusive flux $\mathsf{D}^{(ij)}$ (a vector per component) divergence
is identical scalar-type: $\tfrac1r\partial_r(r\,\mathsf{F}^{(ij)}_r)$ gives
the geometric source $+\mathsf{F}^{(ij)}_r/r$ (the down-gradient
$-F_r/r$ structure of the conductive/species fluxes,
[axisymmetric §5.3](axisymmetric-derivation.md)).

**Production curvature emerges from the embedding.**  With the cylindrical
velocity gradient $(\nabla\mathbf{U})_{\theta\theta}=u_r/r$ (the metric
strain of azimuthal stretching, identical to $S_{\theta\theta}$), the exact
production gives $P_{\theta\theta} = -2\rho R_{\theta\theta}\,u_r/r$ — the
curvature production **emerges** from the projected gradient, not bolted on
(sympy §S3).  At the axis $r\to0$, $u_r/r\to\partial_r u_r$ (de l'Hôpital,
$u_r$ odd), exactly the `Strain_Tensor_Axi` branch.

### A.4 The exact axisymmetric system (per component)

Collecting A.1–A.3, each retained component obeys

$$
\partial_t(\rho R_{ij}) + \mathcal D_{\text{2D}}\!\bigl(\rho\mathbf{U}R_{ij} - \mathsf{D}^{(ij)}\bigr)
= P_{ij} + \Pi_{ij} - \tfrac23\rho\varepsilon\delta_{ij}
\;\underbrace{-\,\frac{\rho u_r R_{ij}}{r} + \frac{\mathsf{F}^{(ij)}_r}{r}}_{\text{geometric sources}},
$$

with $\mathcal D_{\text{2D}}(F)=\partial_x F_x + \partial_r F_r$ the planar
divergence the solver reconstructs, $P_{\theta\theta}=-2\rho R_{\theta\theta}u_r/r$,
and the $33\equiv\theta\theta$ row retained.  The Eulerian metric source
$-\rho u_r R_{ij}/r$ is the literal stress-component analogue of the scalar
$-\rho_\phi u_r/r$ of every passive/RANS scalar
([axisymmetric §4](axisymmetric-derivation.md)).

### A.5 area-variation: $z$-averaging of the fluctuation-tensor transport

The regime depth-averages the 3D equations across a thin layer between two
impermeable free-slip virtual walls $z=\pm h/2$, $\mathbf{G}=-\nabla h/h$.
Apply the $h$-weighted master identity
([area-variation §5.6, §8](area-variation-derivation.md)) to each stress
flux $\mathbf{F}^{(ij)} = \rho\mathbf{U}R_{ij} - \mathsf{D}^{(ij)}$.  The
wall traction vanishes (convection: $\mathbf{U}\!\cdot\!\hat{\mathbf n}=0$
impermeability; diffusion: the turbulence flux normal to the free-slip wall
is zero at leading order, [§5](#5-boundary-conditions)), leaving the
**geometric source** $+\mathbf{G}\!\cdot\!\langle\mathbf{F}^{(ij)}\rangle_z$:

$$
\boxed{\;\partial_t(\rho R_{ij}) + \nabla_{\text{2D}}\!\cdot\!\langle\mathbf{F}^{(ij)}\rangle_z
= P_{ij} + \Pi_{ij} - \tfrac23\rho\varepsilon\delta_{ij}
\;\underbrace{+\,\rho R_{ij}(\mathbf{V}\!\cdot\!\mathbf{G}) - \mathsf{F}^{(ij)}\!\cdot\!\mathbf{G}}_{\text{geometric sources}}\;}
$$

— convective $+\rho R_{ij}(\mathbf{V}\!\cdot\!\mathbf{G})$ (the species/$k$
structure of [§8.2](area-variation-derivation.md)) and diffusive
$-\mathsf{F}^{(ij)}\!\cdot\!\mathbf{G}$ (the negative Anderson-coherent sign
of §2.1).

**The out-of-plane production via Lemma 3.**  The depth-averaged mean
strain has $S_{33} = \langle\partial_z w\rangle_z = -\mathbf{V}\!\cdot\!\mathbf{G}$
([area-variation derivation](area-variation-derivation.md) Lemma 3, §6.4,
*exact*, $w=-(\mathbf{V}\!\cdot\!\mathbf{G})z$), and $\partial_z u,\partial_z v$
average to zero ($z$-odd) — the block-diagonal embedding A.1 with
$(\nabla\mathbf{U})_{33}=-\mathbf{V}\!\cdot\!\mathbf{G}$.  Hence

$$
P_{33} = -2\rho R_{33}(\nabla\mathbf{U})_{33} = +2\rho R_{33}\,(\mathbf{V}\!\cdot\!\mathbf{G}),
$$

the depth-compression production of $\langle w'w'\rangle$ (sympy §S3, §S5).

!!! note "Well-posedness of the depth-averaged $\langle w'w'\rangle$ (honest scope)"

    The transport of $R_{33}=\langle w'w'\rangle$ in the depth-averaged
    regime is closed at the **same leading order** as the mean-flow
    Definition A: the turbulence statistics $R_{ij}$ are treated as
    $z$-even at $O(\varepsilon^2)$ (like the primary fields,
    [area-variation §6.1](area-variation-derivation.md)), so
    $\langle R_{ij}\rangle_z = R_{ij}(x,y)$ and the system closes with the
    production $P_{33}=+2\rho R_{33}(\mathbf{V}\!\cdot\!\mathbf{G})$ and the
    pressure-strain redistribution.  This is a **stated leading-order
    modelling assumption**, not an exact identity: a $z$-dependent
    $R_{33}$ profile would add $O(\varepsilon^2)$ corrections outside the
    closure — exactly the status of every other turbulence quantity in the
    Q2D regime, and consistent with the $O(\varepsilon^2)$ model error of
    the area-variation derivation §1.4.

### A.6 The $7\to5$ reduction — mirror-symmetry proof

In each regime a discrete reflection leaves the mean flow invariant while
flipping the out-of-plane fluctuation:

| regime | reflection | $u_3'$ parity | $\Rightarrow R_{13},R_{23}$ |
|---|---|---|---|
| 2D-plane | spanwise invariance $\partial_z\equiv0$, $w\equiv0$ | — | $\langle u'w'\rangle=\langle v'w'\rangle=0$ |
| 2D-axi | $\theta\to-\theta$ | $u_\theta'\to-u_\theta'$ | $\langle u_x'u_\theta'\rangle=\langle u_r'u_\theta'\rangle=0$ |
| area-var | $z\to-z$ (H9 $z$-even data) | $w'$ $z$-odd | $\langle u'w'\rangle,\langle v'w'\rangle$ $z$-odd $\Rightarrow\langle\cdot\rangle_z=0$ |

The cross-stresses $R_{13},R_{23}$ therefore vanish.  It remains to show
they are **never regenerated**: with $R_{13}=R_{23}=0$ and the
block-diagonal gradient, the source of the $13$ and $23$ equations is
identically zero.  *Proof (term-by-term, sympy §S2).*  Production:
$P_{13} = -\rho(R_{1k}G_{3k}+R_{3k}G_{1k})$; $G_{3k}=0$ except $k=3$ (and
$R_{13}=0$), $R_{3k}=0$ except $k=3$ (and $G_{13}=0$), so $P_{13}=0$.
Pressure-strain: $a_{13}=R_{13}/k=0$ kills $C_1$; $a_{1k}a_{k3}=0$ kills
$C_2$; $S^\ast_{13}=S_{13}=0$ kills $C_3$; $a_{1k}S_{3k}+a_{3k}S_{1k}=0$
kills $C_4$; $W_{3k}=0$ (in-plane vorticity only) kills $C_5$.  Dissipation
$\tfrac23\rho\varepsilon\delta_{13}=0$.  Diffusion of a zero field with a
block-diagonal GGDH tensor stays zero.  Hence $\text{Source}_{13}\equiv\text{Source}_{23}\equiv0$:
the four-stress $+\,\omega$ system is **closed**. $\square$

### A.7 2D-plane: the Cartesian limit (MOSE cross-check)

Setting $(\nabla\mathbf{U})_{33}=0$ and dropping the geometric sources, the
reduced five-equation system is the 3D-Cartesian MOSE model restricted to
$\partial_z(\cdot)=0$, $w=0$, $R_{13}=R_{23}=0$ — verified term-by-term
(sympy §S6: the reduced sources equal MOSE's plane restriction).  $R_{33}$
remains a transported variable: $P_{33}=0$ but the pressure-strain
redistributes energy into it (it is **not** algebraically slaved).  This
closes the loop demanded by the methodology — the projected forms reproduce
the MOSE code in the plane limit, and only there is MOSE used as a check.

### A.8 Review arc and formal assertions

**Adversarial review (both sides).**

1. *Reviewer:* "the $R_{rr}$ equation should carry a $-2R_{\theta\theta}u_r/r$
   curvature source like the literature axisymmetric RSM."  *Resolution:*
   that term appears in formulations **with swirl or in a rotating frame**;
   for the inertial no-swirl reduction the convective transport has no
   inter-component coupling (A.3, sympy §S4).  The curvature enters only
   $P_{\theta\theta}$ (production) and the pressure-strain (through
   $S_{\theta\theta}$, $a_{\theta\theta}$) — confirmed.
2. *Reviewer:* "dropping $R_{13},R_{23}$ may break the pressure-strain
   trace $a_{kl}a_{kl}$."  *Resolution:* the traces still run over the full
   $3\times3$; only the off-diagonal $a_{13},a_{23}$ are zero, which is
   consistent (sympy §S2 keeps the full traces).
3. *Reviewer:* "the depth-averaged $\langle w'w'\rangle$ is not rigorous."
   *Resolution:* conceded as a **leading-order** closure (A.5 note), with
   the same $O(\varepsilon^2)$ status as every Q2D turbulence quantity —
   stated, not hidden (Legge #4).
4. *Reviewer:* "absolute vorticity in $C_5$ needs the frame term."
   *Resolution:* the in-plane vorticity entering $C_5$ is shifted to
   the absolute frame ($+2\boldsymbol\Omega$) whenever the rotating
   frame is enabled (`Lib_SSGLRR`); without frame rotation the base
   regimes use the relative vorticity.

| # | Assertion | Type | Source |
|---|---|---|---|
| A1 | $\boldsymbol\Pi,\boldsymbol\varepsilon$, realizability are frame-invariant $3\times3$ algebra | identity | §6.3, A.1 |
| A2 | $\nabla\mathbf{U}$ block-diagonal; $(\nabla\mathbf{U})_{33}\in\{0,u_r/r,-\mathbf{V}\!\cdot\!\mathbf{G}\}$ | construction | A.1 |
| A3 | axi convective transport: only $-\rho u_r R_{ij}/r$, no Christoffel coupling | **exact** | A.3, §S4 |
| A4 | $P_{33} = -2\rho R_{33}(\nabla\mathbf{U})_{33}$ (curvature/dilatation) | identity | A.3–A.5, §S3 |
| A5 | area-var geometric sources $+\rho R_{ij}(\mathbf{V}\!\cdot\!\mathbf{G})$, $-\mathsf{F}^{(ij)}\!\cdot\!\mathbf{G}$ | asymptotic (LO) | A.5, §S5 |
| A6 | $R_{13}=R_{23}=0$ never regenerated ($\text{Source}_{13}\equiv\text{Source}_{23}\equiv0$) | **exact** | A.6, §S2 |
| A7 | plane limit of the reduced system $=$ MOSE plane restriction | identity | A.7, §S6 |

Assertions A3, A6, A7 are exact identities; A4 follows from the exact
embedding; A5 is leading-order asymptotic (the only modelled step, A.5
note).  The derivation is rigorous at leading order with the controlled
$O(\varepsilon^2)$ status of the host regime.

---

## References

1. F. R. Menter, M. Kuntz, R. Langtry, "Ten years of industrial
   experience with the SST turbulence model", in *Turbulence, Heat and
   Mass Transfer 4*, Begell House, 2003.
2. D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries,
   2006.
3. P. R. Spalart, S. R. Allmaras, "A one-equation turbulence model for
   aerodynamic flows", AIAA-92-0439, 1992.
4. P. R. Spalart, M. Shur, "On the sensitization of turbulence models to
   rotation and curvature", *Aerospace Science and Technology*,
   1(5):297–302, 1997.
5. S. B. Pope, *Turbulent Flows*, Cambridge University Press, 2000.
6. C. G. Speziale, S. Sarkar, T. B. Gatski, "Modelling the pressure–strain
   correlation of turbulence: an invariant dynamical systems approach",
   *Journal of Fluid Mechanics*, 227:245–272, 1991 (SSG).
7. B. E. Launder, G. J. Reece, W. Rodi, "Progress in the development of a
   Reynolds-stress turbulence closure", *Journal of Fluid Mechanics*,
   68(3):537–566, 1975 (LRR).
8. B. J. Daly, F. H. Harlow, "Transport equations in turbulence",
   *Physics of Fluids*, 13(11):2634–2649, 1970 (GGDH diffusion).
9. B. Eisfeld, "Reynolds stress modelling for aeronautical flows", in
   *Progress in Hybrid RANS-LES Modelling*, Springer, 2010 (SSG/LRR-$\omega$
   blending).
10. A.-D. Cécora, R. Radespiel, B. Eisfeld, A. Probst, "Differential
    Reynolds-stress modeling for aeronautics", *AIAA Journal*,
    53(3):739–755, 2015.
11. NASA Turbulence Modeling Resource, "SSG/LRR-RSM-w2019",
    `https://turbmodels.larc.nasa.gov/rsm-ssglrr.html`.
12. R. B. Bird, W. E. Stewart, E. N. Lightfoot, *Transport Phenomena*,
    2nd ed., Wiley, 2002 (cylindrical tensor divergence).
