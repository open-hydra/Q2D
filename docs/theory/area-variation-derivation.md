# Q2D Area-Variation Model — Formal Foundation

This page is the **single authoritative derivation** of the quasi-two-dimensional
area-variation model: hypotheses, exact identities, the complete order analysis,
the field structure in both admissible Reynolds regimes, the depth-averaged
closure, the truncation-admissibility theorem, the itemised loss budget, the
validity envelope, and the registry of legitimate objections.  It supersedes,
as the formulation of record, the historical documents listed in Appendix C
(which remain unmodified as audit trail) and absorbs all their errata.

Every claim in this page is labelled **EXACT** (identity, no asymptotics),
**LEADING** (leading-order along the declared path, controlled remainder),
**POSTULATE** (assumed, with declared role), or **CONJECTURE**.  Every step
cites the hypotheses it uses.  Every algebraic claim has a machine-verifiable
twin in `tools/audit/` (Appendix B).

---

## 1. Declaration of the model

### 1.1 What the model is

Let the *pointwise* flow be governed by the 3D compressible (URANS) Navier–
Stokes equations (P1) in a thin symmetric channel
$\Omega = \{(x,y,z): (x,y)\in\Omega_{2\mathrm D},\ |z| < h(x,y)/2\}$
bounded by two impermeable free-slip walls (P2–P4).  The **Q2D model** is:

> the truncation of the depth-averaged 3D system at **leading order in
> $\varepsilon$ along the distinguished path $\Pi$ of §2.4**, retaining
> *every* term of order $\varepsilon^0$ — whatever its Reynolds/Prandtl/
> Schmidt prefactor — and discarding *only* terms of order
> $O(\varepsilon^2)$.

Two clauses of this declaration are load-bearing and historically left
implicit; they are stated here once and used everywhere.

**(Truncation rule M3) — truncate in $\varepsilon$ only, never in
$\mathrm{Re}$.**  A term of magnitude $\varepsilon^0 \cdot g(\mathrm{Re},
\mathrm{Pr}, \mathrm{Sc}, \mu_t/\mu_l)$ is retained for **any** $g$; a term of
magnitude $\varepsilon^2\cdot g$ is discarded for any $g$.  Rationale: the
effective Reynolds number $\mathrm{Re}_{\rm eff} = \rho U L/(\mu_l+\mu_t)$
varies by 4–6 orders of magnitude across a turbulent reacting domain, so no
fixed ordering between $1/\mathrm{Re}_{\rm eff}$ and $\varepsilon^2$ holds
uniformly in space; truncating in $\mathrm{Re}$ would discard physics that is
leading elsewhere in the same domain.  The model is therefore the
$\varepsilon$-leading-order system **uniformly over the admitted Reynolds
regimes** of (M2).  Consequence (stated for honesty): in the
inviscid-dominant regime the model retains viscous geometric sources that are
*smaller* than the discarded $O(\varepsilon^2)$ terms; they are inert there
and become leading where $\mathrm{Re}_{\rm eff}$ drops locally (§13, Q7).

**(First discarded order is the second).**  Under the symmetry hypotheses
(G1)+(D1) all odd-order coefficients of the $\varepsilon$-expansion vanish
identically (Theorem 5.1).  The gap between what is kept ($\varepsilon^0$) and
what is discarded ($\varepsilon^2$) is therefore **two powers of
$\varepsilon$ by theorem, not by choice**; if the $z$-symmetry is broken
(curved mid-surface, $z$-asymmetric injection) the model degrades to
first-order accuracy (§9.3).

### 1.2 What the model is not

- Not a shallow-water model: no free surface, no gravity driver, walls are
  rigid, $\partial_t h \equiv 0$; the viscous no-slip shallow-water regime
  $\mathrm{Re}\,\varepsilon^2 = \Theta(1)$ is excluded by (M2).
- Not a homogenised rough-channel model: wall corrugation on scales shorter
  than $L$ is excluded by (G1); see the counterexample in §2.6.
- Not a low-Mach or hypersonic asymptotic model: the path $\Pi$ fixes
  $\mathrm{Ma} = \Theta(1)$ (§2.4); other Mach paths require re-scaled
  pressure and are outside the certified envelope (§10.4).
- Not a model of gap-normal wave dynamics: transverse acoustics and
  wall-to-wall shock reflection live at the discarded order (§9.2).

### 1.3 Certified configuration class $\mathcal{C}$ (target of this document)

The validity statements of §8–§10 are asserted for the configuration class
the model is built for:

> **$\mathcal{C}$:** an **unwrapped RDE annulus** (2D plane: $x$ azimuthal-
> periodic, $y$ axial) or, equivalently, a **convergent–divergent nozzle with
> free-slip out-of-plane walls**, in which
>
> 1. the gap $h(x,y)$ varies **gently and smoothly**, like a small C–D
>    nozzle: single-scale, no corrugation, slope $\sigma = O(\varepsilon)$
>    — i.e. (G1) holds with moderate constants $c_0, M$;
> 2. the **initial data carry no out-of-plane gradients**:
>    $\partial_z f|_{t=0} \equiv 0$ for all primary fields and
>    $w|_{t=0} = -(\mathbf V\!\cdot\!\mathbf G)\,z$ (or $0$ where
>    $\mathbf V_0 = 0$) — hypothesis (D2) in its **strong form**, which makes
>    the second-order data contribution vanish identically,
>    $f_2(\cdot,0) \equiv \langle f_2\rangle(\cdot,0)$;
> 3. in-plane boundary data are $z$-invariant (D1);
> 4. the flow is smooth URANS-mean except for captured fronts, whose local
>    violation of (M1) is declared and bounded in §10.5.

Everything in §3–§7 is derived for general data satisfying (D1)–(D2-weak);
class $\mathcal C$ is the case in which every constant of the error budget is
sharpest and the regime-(b) propagation argument (§5.3) is cleanest.

---

## 2. Hypotheses

Conventions: indices $1=x$, $2=y$, $3=z$; $\mathbf V_{3\mathrm D} = (u,v,w)$;
$h_x = \partial_x h$;

$$
\mathbf G := -\frac{\nabla h}{h},\qquad
\mathbf V\!\cdot\!\mathbf G = u G_x + v G_y,\qquad
\langle f\rangle_z := \frac1h\int_{-h/2}^{h/2} f\,\mathrm dz,\qquad
f' := f - \langle f\rangle_z,\qquad [f] := f|_+ - f|_-,
$$

with $f|_\pm = f(x,y,\pm h/2,t)$.  Two dimensionless geometric measures are
used throughout and **must not be conflated**:

$$
\varepsilon := \frac{h_{\max}}{L}\ \ \text{(thinness)},\qquad\qquad
\sigma := \lVert\nabla h\rVert_{\infty}\ \ \text{(wall slope)} .
$$

They are linked by $\sigma = \varepsilon\,(L/L_h)$ where
$L_h := \inf_{\Omega_{2\mathrm D}} h/|\nabla h|$ is the geometry's own
variation length.  Thinness does **not** imply gentle slope (§2.6,
counterexample CE-1); the model needs both.

### 2.1 Geometry

**(G1) Admissible geometry** — *(replaces legacy (H2)+(H4); incorporates the
previously implicit slope hypothesis).*  There exist a length $L>0$ — the
**problem scale**, fixed by the geometry and the boundary/initial data, *not*
defined through the solution — a parameter $\varepsilon\in(0,\varepsilon_0]$,
and constants $c_0\in(0,1]$, $M\ge 1$ independent of $\varepsilon$, such that

$$
h(x,y) \;=\; \varepsilon\,L\,\hat h\!\Big(\frac{x}{L},\frac{y}{L}\Big),
\qquad
\hat h\in C^2(\overline{\Omega^*_{2\mathrm D}}),\qquad
c_0\le \hat h\le 1,\qquad
\lVert\hat h\rVert_{C^2}\le M,
$$

and $\partial_t h\equiv 0$.

**Lemma 2.1 (slope lemma).**  Under (G1):

$$
\frac{h_{\max}}{L}\le\varepsilon,\qquad
\boxed{\;\sigma=\lVert\nabla h\rVert_\infty\le M\varepsilon\;},\qquad
|\mathbf G|\le \frac{M}{c_0 L},\qquad
|\nabla\mathbf G|\le \frac{C(M,c_0)}{L^2},\qquad
L_h\ \ge\ \frac{c_0}{M}\,L .
$$

*Proof.* Chain rule on $h=\varepsilon L\hat h(x/L,y/L)$:
$\nabla h = \varepsilon\nabla_{*}\hat h$, so $\sigma\le M\varepsilon$;
$\mathbf G = -\nabla_*\hat h/(L\hat h)$ gives $|\mathbf G|\le M/(c_0L)$ and,
differentiating once more with $\hat h\in C^2$, $|\nabla\mathbf G| =
O(1/L^2)$. $\square$  **EXACT.**

The statement "$h_x = O(\varepsilon)$", used at every ordering step below, is
from here on a **theorem** (Lemma 2.1), not an unstated assumption.  Its
independent necessity is shown by counterexample CE-1 in §2.6 and by the
kinematic compatibility Lemma 5.4.

### 2.2 Physics

**(P1) Governing system** — *(legacy (H1))*.  3D compressible multi-species
(U)RANS in conservative form with Newton–Stokes stress
$\tau_{ij}=2\mu S_{ij}-\tfrac23\mu(\partial_\ell u_\ell)\delta_{ij}$
(Stokes hypothesis $\lambda=-\tfrac23\mu$; bulk-viscosity caveat in §13 Q10),
Fourier flux $\mathbf q=-\kappa\nabla T$, Fick diffusion with symmetric
Hirschfelder correction $\mathbf j_k=-\rho D_k\nabla Y_k+\rho Y_k\mathbf V_c$,
$\mathbf V_c=\sum_\beta D_\beta\nabla Y_\beta$, total energy
$E = e(T,Y)+\tfrac12|\mathbf V|^2$, chemical source $\omega_k(T,Y,\rho)\in
C^2$ with the Damköhler/Zel'dovich bound of (M2-iv).  For RANS the same
structure holds with $\mu_{\rm eff}=\mu_l+\mu_t$, Boussinesq
$R_{ij}=2\mu_tS_{ij}-\tfrac23\rho k\,\delta_{ij}$, and $k,\omega$ (or
$\tilde\nu$) as primary fields.

**(P2) Impermeable free-slip walls** — *(legacy (H5))*.
$\mathbf V\!\cdot\!\hat{\mathbf n}=0$ and
$(\boldsymbol\tau_{\rm tot}\!\cdot\!\hat{\mathbf n})\!\cdot\!\mathbf t=0$
pointwise on $z=\pm h/2$ ($\boldsymbol\tau_{\rm tot}=\boldsymbol\tau+R$ for
RANS).

**(P3) Adiabatic walls** — *(legacy (H6))*. $\mathbf q_{\rm tot}\!\cdot\!
\hat{\mathbf n}=0$.

**(P4) Non-catalytic walls** — *(legacy (H7))*. $\mathbf j_k\!\cdot\!
\hat{\mathbf n}=0$ for every $k$; for the turbulence variables
$\partial_n k=\partial_n\omega=\partial_n\tilde\nu=0$ (mirror-extension
consistency; legacy §N12).

*Physical scope note (declared, quantified in §10.6):* (P2)–(P4) are an
**idealisation** of the real no-slip, conducting walls of an RDE gap.  The
$O(\varepsilon^2)$ error statements of this page are **relative to the
idealised 3D reference**; the idealisation error itself (boundary-layer
momentum deficit, wall heat flux) is a separate, additive model-form error.

### 2.3 Data

**(D1) $z$-symmetric data and BC** — *(legacy (H9))*.  Initial data and
in-plane boundary data are even in $z$ for
$\rho,u,v,T,p,Y_k,k,\omega,\tilde\nu$ and odd in $z$ for $w$.

**(D2) Prepared data** — *(legacy (G6), strengthened for class
$\mathcal C$).*
*Weak form:* the initial data do not depend on $\varepsilon$ explicitly;
hence by Taylor expansion in $z=\varepsilon L z^*$ their $z^*$-dependence is
$O(\varepsilon^2)$ (Lemma 5.2).
*Strong form (class $\mathcal C$):* $\partial_z f|_{t=0}\equiv 0$ for the
primary fields and $w|_{t=0}=-(\mathbf V\!\cdot\!\mathbf G)z$; then
$f_2(\cdot,0)$ is $z^*$-independent and the entire initial fluctuation
vanishes.

### 2.4 Method: expansion and the distinguished path $\Pi$

**(M1) Regular asymptotic expansion** — *(legacy (H8))*. **POSTULATE.**
Each primary field admits
$f(x,y,z,t;\varepsilon)=\sum_{n\ge0}\varepsilon^n f_n(x,y,z^*,t)$,
$z^*=z/(\varepsilon L)$, truncatable at finite order with remainder
$O(\varepsilon^{N+1})$ uniformly on $\overline\Omega\times[0,T_{\rm fin}]$,
coefficients $f_n$ with bounded $C^2$ norms **in the slow variables
$x^*=x/L$** (this clause carries the statement "the fields vary on the scale
$L$", removing the circular operational definition of $L$ used historically).

**(M2) Distinguished path $\Pi$.**  The model is the leading-order truncation
along the parameter path

$$
\Pi:\quad \varepsilon\to0\ \ \text{with}\ \
\begin{cases}
\text{(i)}& \sigma\le M\varepsilon\quad\text{(automatic from (G1))},\\[2pt]
\text{(ii)}& \mathrm{Re}_{\rm eff}\,\varepsilon^2\notin[\Lambda^{-1},\Lambda]
\ \ \text{pointwise-zonally, i.e. either}\\
&\quad\text{(a) viscous-transverse-dominant: } \mathrm{Re}_{\rm eff}\,
\varepsilon^2\le\Lambda^{-1},\ \ \text{or}\\
&\quad\text{(b) inviscid-dominant: } \mathrm{Re}_{\rm eff}\,\varepsilon^2\ge
\Lambda,\\[2pt]
\text{(iii)}& \mathrm{Ma}=\Theta(1)\quad\text{(pressure scale }\rho_{\rm
ref}U^2\text{)},\\[2pt]
\text{(iv)}& \mathrm{Da}=\Theta(1)\ \text{and}\ \beta^2\varepsilon^4\ll1,
\ \ \beta := E_a/(R\,T_0)\ \text{(Zel'dovich)} .
\end{cases}
$$

Departures from $\Pi$ are **different models**: $\sigma/\varepsilon\to\infty$
$\Rightarrow$ homogenisation; $\mathrm{Re}\,\varepsilon^2=\Theta(1)$
$\Rightarrow$ viscous shallow water (Gerbeau–Perthame); $\mathrm{Ma}\to0$
$\Rightarrow$ low-Mach rescaling; $\mathrm{Da}\to\infty$ $\Rightarrow$ front
inner problems ((M1) violated locally).

**(M3) Truncation rule.**  As declared in §1.1: keep all $\varepsilon^0$
terms uniformly over regimes (a)/(b); discard all $O(\varepsilon^2)$ terms.

### 2.5 Technical regularity

**(T1)** — *(legacy (H11))*. A $C^2$ solution of (P1)–(P4) exists on
$\overline\Omega\times[0,T_{\rm fin}]$.  Local-in-time existence and
uniqueness of classical solutions hold: for the viscous system by
Matsumura–Nishida (1980) / Valli (1982); for the inviscid (regime-(b))
leading system by Kato (1975) / Majda (1984), and in bounded domains with
slip walls by Schochet (1986).  Global regularity is open in general;
$T_{\rm fin}$ is the smooth horizon. **POSTULATE** beyond the local theory.

**(T2)** — *(legacy (H12))*. $\mu,\kappa,D_k,c_p\in C^1$ with bounded
derivatives on the admissible state set; $\mu_{\min}=\inf\mu>0$.

**(T3)** — *(legacy (H13))*. The expansion coefficients $f_n$ are uniformly
bounded in $z^*\in[-\tfrac12,\tfrac12]$ together with their first
$z^*$-derivatives.

**(T4)** — *(legacy (H14))*. The $\varepsilon$-dependence is $C^\infty$ with
uniform Taylor remainders (no $\sqrt\varepsilon$ or
$\varepsilon\log\varepsilon$ layers).

**(T5) Uniform leading bounds** — *(formal definition of the historical,
previously undefined label "(H_b)").* **POSTULATE.**  On $[0,T_{\rm fin}]$
the leading-order fields satisfy
$\lVert(\rho_0,u_0,v_0,T_0,Y_{k,0})\rVert_{C^2(\overline\Omega^*_{2\mathrm
D})}\le K_0$ and $\rho_0\ge\rho_{\min}>0$, $T_0\ge T_{\min}>0$, with
$K_0,\rho_{\min},T_{\min}$ independent of $\varepsilon$ and of
$\mathrm{Re}_{\rm eff}$ within the admitted regimes.  Role: makes every
order-constant of §4–§9 finite and $\varepsilon$-uniform; in regime (b) it
substitutes for the unavailable global inviscid theory.

### 2.6 Necessity table and counterexamples

| Hypothesis | First place used | What fails without it | Counterexample / failure mode |
|---|---|---|---|
| (G1) thinness | §4 (all tables) | no small parameter | — |
| (G1) **slope clause** ($\sigma\le M\varepsilon$) | §4.3 BC ordering; §5.4 | wall terms promote to $O(1)$; $w$-scales incompatible | **CE-1**: $h=\varepsilon L\,(1+\tfrac12\sin\frac{x}{\varepsilon L})$ is thin, $C^2$, $h\ge h_{\min}$, yet $\sigma=O(1)$, $L_h=O(\varepsilon L)$: in-plane geometric derivatives become as large as transverse ones; correct limit is a **homogenised** model, not Q2D |
| (P2) impermeability | §3.4 | $[w]$ undetermined; no $\mathbf V\!\cdot\!\mathbf G$ | porous wall: mass source appears |
| (P2) slip | §5.2; §7.3 | wall shear layers; plug profile lost | no-slip: Gerbeau–Perthame regime, shape factor 6/5, friction |
| (P3)+(P4) | §7.5–7.6 | energy/species wall fluxes survive | conducting wall: Nusselt closure needed |
| (D1) | §5.1 | odd orders survive | asymmetric injection: error degrades to $O(\varepsilon)$ |
| (D2) | §5.3 | $f_2(0)=O(1)$; in regime (b) **never relaxes** | $z$-structured injectors: permanent $O(1)$ $z$-fluctuation transported by Euler |
| (M1) | §5 | no hierarchy | shock interior: gradients on $\delta_s\ll L$ |
| (M2-ii) | §5.2–5.3 | neither closure mechanism valid | $\mathrm{Re}\varepsilon^2=\Theta(1)$: $z$-profile genuinely 3D |
| (M2-iii) | §4.1 | pressure scale wrong | $\mathrm{Ma}\to0$: $p=O(\rho U^2/\mathrm{Ma}^2)$ |
| (M2-iv) | §7.7 | chemistry closure constant blows up | $\beta^2\varepsilon^4=O(1)$: $\langle\omega\rangle\ne\omega(\langle\cdot\rangle)$ at leading order |
| (T1)–(T5) | everywhere | constants not $\varepsilon$-uniform | blow-up, vacuum, $\mu\to0$ |

No hypothesis is implied by the others.  In particular CE-1 proves that the
slope clause of (G1) has **exclusion content** beyond thinness — it is the
formal answer to the question "is gentle slope implicit in thinness?":
*no*; it was the intended but never-formalised hypothesis of the historical
documents (Appendix C, item C-1).

---

## 3. The exact kinematic layer (no asymptotics)

Everything in this section is an **identity**, valid for any $h\in
C^1$, any slope, any Reynolds number.  At this level there is no model yet:
the wall traces are not expressible through the averages alone (the closure
problem; §13 Q1).

**3.1 Linearity and time commutation.** $\langle\cdot\rangle_z$ is linear and,
since $\partial_t h=0$,
$\langle\partial_t f\rangle_z=\partial_t\langle f\rangle_z$. **EXACT.**

**3.2 Lemma (Leibniz commutation).** For $f\in C^1(\overline\Omega)$,

$$
\langle\partial_x f\rangle_z=\partial_x\langle f\rangle_z-G_x\langle
f\rangle_z+\frac{G_x}{2}\big(f|_++f|_-\big),
$$

and likewise for $\partial_y$.  *Proof:* Leibniz rule for the moving limits
$\pm h/2$, divide by $h$, use $h_x/h=-G_x$. $\square$ **EXACT.**

**3.3 Fundamental theorem of calculus.**
$\langle\partial_z f\rangle_z=[f]/h$. **EXACT.**

**3.4 Lemma (wall traces of $w$).**  Impermeability (P2) with
$\hat{\mathbf n}_\pm\,\mathrm dA=(-h_x/2,-h_y/2,\pm1)\,\mathrm dx\,\mathrm dy$
gives, pointwise and exactly,

$$
w|_\pm=\pm\tfrac12\big(h_x\,u|_\pm+h_y\,v|_\pm\big)
=\mp\tfrac{h}{2}\,(\mathbf V\!\cdot\!\mathbf G)\big|_\pm ,
\qquad
[w]=\tfrac{h_x}{2}\big(u|_++u|_-\big)+\tfrac{h_y}{2}\big(v|_++v|_-\big).
$$

**EXACT** (BC identity, independent of any expansion).

**3.5 Lemma (key kinematic identity).** Combining 3.3 and 3.4,

$$
\boxed{\;\langle\partial_z w\rangle_z=\frac{[w]}{h}
=-\Big(G_x\,\bar u_{\rm w}+G_y\,\bar v_{\rm w}\Big)\;}\qquad
\bar u_{\rm w}:=\tfrac12(u|_++u|_-),
$$

**EXACT in the wall-trace averages**; replacing traces by depth averages,
$\langle\partial_z w\rangle_z=-\mathbf V\!\cdot\!\mathbf G+O(\varepsilon^2)
\cdot U/L$ under (D1)+(M1) (**LEADING**; the historical "exact" label refers
to the trace form).  The identity is purely kinematic: a definite integral
equal to a BC-determined jump — no profile information is used.

**3.6 Master identity.** For any flux $\mathbf F$,

$$
\langle\nabla_{3\mathrm D}\!\cdot\!\mathbf F\rangle_z
=\nabla_{2\mathrm D}\!\cdot\!\langle\mathbf F\rangle_z
-\mathbf G\!\cdot\!\langle\mathbf F\rangle_z
+\frac1h\Big(\mathbf F\!\cdot\!\hat{\mathbf n}\big|_+ +\mathbf
F\!\cdot\!\hat{\mathbf n}\big|_-\Big).
$$

*Proof:* sum 3.2 (twice) and 3.3; the trace groups assemble into
$\mathbf F\!\cdot\!\hat{\mathbf n}|_\pm$. $\square$ **EXACT.**  The middle
term is the geometric source; the last term is the physical wall flux,
evaluated under (P2)–(P4) in §7.

---

## 4. General two-parameter order analysis

This section is **general bookkeeping**: every term of every equation
receives its magnitude as a function of the groups
$(\varepsilon,\sigma,\mathrm{Re},\mathrm{Pr},\mathrm{Sc},\mathrm{Ma},
\mathrm{Da})$, *before* any truncation.  The model emerges only in §5–§8 when
the path $\Pi$ is invoked.  This is the formal sense in which "an order
analysis is a general analysis".

**4.1 Scales.** $x,y\sim L$; $z\sim\varepsilon L$ (from (G1)); $u,v\sim U$;
$t\sim L/U$; $\rho\sim\rho_{\rm ref}$; $p\sim\rho_{\rm ref}U^2$ (valid along
$\Pi$-(iii): for $\mathrm{Ma}\to0$ the pressure scale is $\rho U^2/
\mathrm{Ma}^2$ and the tables below must be re-derived);
$\mu\sim\mu_{\rm ref}$; geometry: $|\nabla h|\le\sigma$,
$|\mathbf G|\le M/(c_0L)$ (Lemma 2.1).

**4.2 Scale of $w$ (deduced).**  Non-dimensional continuity reads
$\partial_{x^*}(\rho^*u^*)+\partial_{y^*}(\rho^*v^*)+\varepsilon^{-1}
\partial_{z^*}(\rho^*\hat w)=0$ with $\hat w=w/U$; a non-degenerate balance —
*given that the fields vary on $L$, clause (M1)* — forces
$\hat w=O(\varepsilon)$, i.e. $w=O(\varepsilon U)$.  Independently, Lemma 3.4
gives $|w|_{\rm wall}\le\sigma U$.  The consistency of these two
determinations is **exactly** the slope clause $\sigma=O(\varepsilon)$ —
made precise in Lemma 5.4.

**4.3 $x$-momentum.** Dividing by the convective scale $\rho_{\rm ref}U^2/L$:

| Term | order (general) | along $\Pi$ ($\sigma\le M\varepsilon$) |
|---|---|---|
| $\rho\,\partial_t u,\ \rho u\partial_xu,\ \rho v\partial_yu,\ \rho w\partial_zu$ | $1$ | $1$ |
| $-\partial_x p$ | $1$ | $1$ |
| $\frac1{\mathrm{Re}}\partial^2_{x}u,\ \frac1{\mathrm{Re}}\partial^2_{y}u$, mixed $\lambda$-terms | $1/\mathrm{Re}$ | $1/\mathrm{Re}$ |
| $\boxed{\frac{1}{\mathrm{Re}\,\varepsilon^2}\partial^2_{z^*}u}$ | $1/(\mathrm{Re}\,\varepsilon^2)$ | dominant in regime (a) |
| $\frac1{\mathrm{Re}}\partial_z(\mu\,\partial_x w)$ | $1/\mathrm{Re}$ | $1/\mathrm{Re}$ |

The **only** term carrying $\varepsilon^{-2}$ is the transverse viscous one;
substituting the expansion, $\varepsilon^{-2}$ acts on $f_0$ and
$\varepsilon^{-1}$ on $f_1$ — the structural origin of the hierarchy
(§5.2).  The enumeration of all 14 sub-terms of
$\partial_j\tau_{1j}$ with variable $\mu(T)$, confirming that none other
produces $\varepsilon^{-1}$ or $\varepsilon^{-2}$, is inherited verbatim from
the audited historical analysis (Appendix C, item C-2) and machine-checked
(Appendix B).

**4.4 $z$-momentum.** Dividing by the dominant scale
$\rho_{\rm ref}U^2/(\varepsilon L)$:

| Term | order (general) |
|---|---|
| $\rho\,\mathrm Dw/\mathrm Dt$ | $\varepsilon^2$ |
| $-\partial_z p$ | $1$ (dominant) |
| $\partial_x\tau_{13},\ \partial_y\tau_{23}$ | $\varepsilon^2/\mathrm{Re}$ |
| $\partial_z\tau_{33}$ | $1/\mathrm{Re}$ |

**Caution (corrects the historical argument).**  In regime (a)
($\mathrm{Re}$ fixed) the entry $1/\mathrm{Re}$ is **not** asymptotically
subleading in $\varepsilon$: at order $\varepsilon^0$ the balance is
$-\partial_{z^*}p_0+\mathrm{Re}^{-1}\partial_{z^*}(\tau_{33})_0=0$, and the
conclusion $\partial_{z^*}p_0=0$ requires the additional fact that
$(\tau_{33})_0$ is $z^*$-independent — which is a consequence of the
velocity/temperature hierarchy and is therefore proved **after** it
(Lemma 5.5), not before.  The historical presentation inverted this order
(Appendix C, item C-3).

**4.5 Energy, species, RANS.** Identical structure with transverse-dominant
terms $\frac{1}{\mathrm{Re}\mathrm{Pr}\,\varepsilon^2}\partial^2_{z^*}T$,
$\frac{1}{\mathrm{Re}\mathrm{Sc}\,\varepsilon^2}\partial^2_{z^*}Y_k$,
$\frac{D_\phi}{\rho UL\,\varepsilon^2}\partial^2_{z^*}\phi$; chemical source
$\omega_k$: order $\mathrm{Da}$, with $z^*$-coupling controlled by (M2-iv)
(§7.7).

**4.6 The slip condition, exactly, and its $\sigma$-ordering.**  On
$z=+h/2$, with tangents $\mathbf t_1=(1,0,h_x/2)$, $\mathbf
t_2=(0,1,h_y/2)$ and $\hat{\mathbf n}=(-h_x/2,-h_y/2,1)$, the condition
$(\boldsymbol\tau\!\cdot\!\hat{\mathbf n})\!\cdot\!\mathbf t_1=0$ reads
**exactly**

$$
\tau_{13}\Big(1-\frac{h_x^2}{4}\Big)
=\frac{h_x}{2}\,(\tau_{11}-\tau_{33})+\frac{h_y}{2}\,\tau_{12}
+\frac{h_xh_y}{4}\,\tau_{23},
$$

(and the $x\!\leftrightarrow\!y$ analogue for $\tau_{23}$).  Ordering with
the two parameters: the retained right-hand terms are $O(\sigma)\cdot
\tau_{\rm in\text{-}plane}$; the discarded corrections are **relatively
$O(\sigma^2)$** (the $(1-h_x^2/4)$ factor and the $\tau_{23}$ cross-term).
Along $\Pi$, $\sigma^2\le M^2\varepsilon^2$ and the truncation

$$
\tau_{13}\big|_\pm=\pm\Big[\frac{h_x}{2}\,(\tau_{11}-\tau_{33})
+\frac{h_y}{2}\,\tau_{12}\Big]+O(\sigma^3\,\mu U/L)
$$

is second-order consistent. **LEADING.**  Off the path
($\sigma=\Theta(1)$, e.g. 45° walls) the "discarded" terms are $O(1)$
relative: the truncation is illegitimate — a (L7) failure, visible only in
the two-parameter bookkeeping.

---

## 5. Field structure along $\Pi$

### 5.1 Parity: all odd orders vanish

**Theorem 5.1.** Under (G1)+(P1)–(P4)+(D1)+(D2)+(T1), on the smooth horizon:
the solution is even in $z$ for the primary fields and odd for $w$; moreover
every odd-order coefficient of the $\varepsilon$-expansion vanishes:
$f_n\equiv0$ for odd $n$ (primary fields), $\hat w_n\equiv0$ for even $n$.

*Proof.* **(i) Solution parity.** The reflected field
$(\rho,u,v,-w,T,Y_k,p)(x,y,-z,t)$ solves the same system: each equation of
(P1) is invariant under $(z\to-z,\,w\to-w)$ — verified term by term,
including the transformation of all nine $\tau_{ij}$ and the closures
$\mu(T),\kappa(T),D_k(T)$, which depend only on even scalars.  The geometry
(G1) and the BCs (P2)–(P4) are invariant.  The data (D1) are invariant.
Local-in-time uniqueness of classical solutions — Matsumura–Nishida (1980)/
Valli (1982) for the viscous system; Kato (1975)/Majda (1984)/Schochet
(1986, slip-bounded domain) for the inviscid leading system of regime (b) —
forces the solution to coincide with its reflection on $[0,T_{\rm loc}]$;
iteration extends this over the smooth horizon under (T1). **(ii)
$\varepsilon$-parity of the hierarchy.** In stretched variables
$(x^*,y^*,z^*,t^*)$ with $\hat h$ $\varepsilon$-free (G1), the
non-dimensional system depends on $\varepsilon$ **only through
$\varepsilon^2$** (the factor $\varepsilon^{-2}$ on transverse diffusion;
$\varepsilon$ appears nowhere at odd power once $w=\varepsilon U\hat w$ is
absorbed).  Hence the solution map preserves evenness in $\varepsilon$.
**(iii) Data parity in $\varepsilon$.** (D2-weak): an
$\varepsilon$-independent $z$-even datum, evaluated at $z=\varepsilon L
z^*$, expands in even powers of $\varepsilon$ only (Lemma 5.2).  By (ii) the
expansion retains only even powers for all time; the odd coefficients,
solving linear homogeneous problems with zero data, vanish identically.
$\square$ **LEADING** (rigorous up to the uniqueness horizon; beyond, (T1)).

**Lemma 5.2 (data preparation).**  If $f_{\rm init}$ is $\varepsilon$-
independent and $z$-even ($C^2$), then $f_{\rm init}(x,y,\varepsilon
Lz^*)=f_{\rm init}(x,y,0)+\tfrac{(\varepsilon L z^*)^2}{2}\partial_z^2
f_{\rm init}(x,y,0)+O(\varepsilon^4)$: the leading datum is
$z^*$-independent and the first $z^*$-structure enters at $\varepsilon^2$.
*Strong form (class $\mathcal C$):* $\partial_zf_{\rm init}\equiv0$ gives
$f_0(0)=f_{\rm init}$, $f_n(0)=0$ for all $n\ge1$ — the initial fluctuation
vanishes **identically**, not merely at order. **EXACT.**

### 5.2 Regime (a): transverse-diffusion hierarchy and Fredholm closure

Assume $\Pi$-(ii)(a).  Orders are collected in the topologically correct
sequence — velocities and scalars first, then $w$, then **and only then**
the pressure.

**(a1) $\varepsilon^{-2},\varepsilon^{-1}$:** the only terms at these orders
are $(1/\mathrm{Re})\partial^2_{z^*}u_{0,1}$ (Table 4.3; no pressure — its
gradient is in-plane and $O(1)$).  Hence $\partial^2_{z^*}u_0=0$; with
parity, $u_0=u_0(x^*,y^*,t^*)$. Likewise $v_0$, $T_0$ (transverse
conduction, Neumann (P3)), $Y_{k,0}$ (P4), and the RANS scalars (P4).
$u_1\equiv0$ by Theorem 5.1.

**(a2) $w$-profile.** Continuity at $\varepsilon^0$ with $\rho_0,u_0,v_0$
$z^*$-independent gives $\partial_{z^*}\hat w_1$ $z^*$-independent;
integrating with parity ($w(0)=0$) and the trace Lemma 3.4:

$$
\boxed{\;w=-(\mathbf V\!\cdot\!\mathbf G)\,z+O(\varepsilon^3 U)\;}
\qquad
\partial_zw=-\mathbf V\!\cdot\!\mathbf G,\qquad
\partial_xw=-z\,\partial_x(\mathbf V\!\cdot\!\mathbf G)\ \ (z\text{-odd}),
$$

and the impermeability-compatibility yields the Q2D mass equation as the
solvability condition (first instance of the Fredholm pattern). **LEADING.**

**(a3) Lemma 5.5 ($(\tau_{33})_0$ is $z^*$-independent).**
$\tau_{33}=2\mu\partial_zw-\tfrac23\mu(\nabla_{2\mathrm D}\!\cdot\!\mathbf V
+\partial_zw)$; at leading order every factor — $\mu(T_0)$,
$\partial_zw|_0=-(\mathbf V_0\!\cdot\!\mathbf G)$,
$\nabla_{2\mathrm D}\!\cdot\!\mathbf V_0$ — is $z^*$-independent by (a1)–(a2).
$\square$ **LEADING.**

**(a4) Pressure.** $z$-momentum at $\varepsilon^0$:
$-\partial_{z^*}p_0+\mathrm{Re}^{-1}\partial_{z^*}(\tau_{33})_0=0$; the
second term vanishes by Lemma 5.5, so $p_0=p_0(x,y,t)$, and by the EOS
$\rho_0$ is $z^*$-independent.  (This replaces the historical "viscous terms
subleading" argument; Appendix C, item C-3.) **LEADING.**

**(a5) Second order and Fredholm.** In the wall-fitted coordinate
$\zeta=2z/h\in[-1,1]$ the $\varepsilon^0$ $x$-momentum reduces to

$$
\frac{4}{\mathrm{Re}\,h^2}\,\partial^2_\zeta u_2=G_0(x,y,t),\qquad
G_0=\rho_0\big[\partial_tu_0+u_0\partial_xu_0+v_0\partial_yu_0\big]
+\partial_xp_0-\frac{1}{\mathrm{Re}}\nabla^2_{2\mathrm D}u_0
-\frac{1}{3\mathrm{Re}}\partial_x\big[\nabla_{2\mathrm D}\!\cdot\!\mathbf
V_0-(\mathbf V_0\!\cdot\!\mathbf G)\big],
$$

with $G_0$ $\zeta$-independent, so $u_2=\tfrac{\mathrm{Re}h^2}{8}G_0\zeta^2+
\mathrm{const}$ and $\partial_\zeta u_2(\pm1)=\pm\tfrac{\mathrm{Re}h^2}{4}
G_0$.  The slip condition (§4.6) prescribes the **complete** Neumann datum

$$
\partial_\zeta u_2(+1)=\frac{h}{2\varepsilon^2}\Big\{
\frac{h}{2}\,\partial_x(\mathbf V_0\!\cdot\!\mathbf G)
+h_x\big[\partial_xu_0+(\mathbf V_0\!\cdot\!\mathbf G)\big]
+\frac{h_y}{2}\big(\partial_yu_0+\partial_xv_0\big)\Big\},
$$

— *all three groups are the same order $O(\varepsilon\,U/L)$ before the
$\varepsilon^{-2}$ rescaling; a historical simplification that kept only the
first group is superseded (Appendix C, item C-4)*.  Equating the ODE and BC
determinations (in fully non-dimensional variables the $\varepsilon^{-2}$
factor cancels identically against the $h=\varepsilon L\hat h$ scaling —
no residual $\varepsilon$; legacy gap "C.4" resolution, machine-checked):

$$
\mathrm{Re}\,G_0=\partial_x(\mathbf V_0\!\cdot\!\mathbf G)
-2G_x\big[\partial_xu_0+(\mathbf V_0\!\cdot\!\mathbf G)\big]
-G_y\big(\partial_yu_0+\partial_xv_0\big).
$$

**Theorem 5.6 (Fredholm $\equiv$ model).**  This solvability condition is
**identically** the Q2D $x$-momentum equation of §8 ($\Delta\equiv0$ by
direct expansion; SymPy twin F-1).  The same scheme closes $v_2$, $T_2$
(Neumann from (P3)), $Y_{k,2}$ (P4), $k_2,\omega_2$ (P4/N12), each
solvability condition reproducing the corresponding model equation. $\square$
**LEADING.**  Interpretation: the discarded order **exists and is bounded**
along solutions of the model — the admissibility input (L1) of Theorem 8.2.

### 5.3 Regime (b): inviscid-dominant — full derivation

Assume $\Pi$-(ii)(b): $\mathrm{Re}_{\rm eff}\,\varepsilon^2\ge\Lambda$, so
$1/\mathrm{Re}_{\rm eff}=o(\varepsilon^2)$: at orders $\varepsilon^0$ and
$\varepsilon^2$ the hierarchy is **Euler**; all viscous terms (and the
tangential-traction content of (P2)) enter below $\varepsilon^2$.

**(b0) Pressure.** The $z$-momentum dominance
$-\varepsilon^{-2}\partial_{z^*}p$ is regime-independent (Table 4.4: the
inertia is $\varepsilon^2$, the viscous entries vanish along (b)); hence
$\partial_{z^*}p_0=\partial_{z^*}p_1=0$ directly — in (b) the pressure
argument is *simpler* than in (a). **LEADING.**

**(b1) Leading datum.** By Lemma 5.2, $f_0(\cdot,0)$ is $z^*$-independent
(class $\mathcal C$: exactly; weak (D2): up to $O(\varepsilon^2)$ absorbed
in $f_2(0)$).

**(b2) Propagation theorem.**  *The $\varepsilon^0$ Euler system preserves
$z^*$-independence of $(\rho_0,u_0,v_0,T_0,Y_{k,0},p_0)$ with $\hat w_1$
linear in $z^*$.*
*Proof.* Candidate: let $(\rho_0,u_0,v_0,p_0,\dots)(x,y,t)$ solve the Q2D
system of §8 with data $f_0(\cdot,0)$, and set $\hat w_1=-(\mathbf
V_0\!\cdot\!\mathbf G)\,z^*L/U$ (trace-compatible by Lemma 3.4).
Substitution into the $\varepsilon^0$ Euler hierarchy: continuity is the Q2D
mass equation (the term $\rho_0\partial_{z^*}\hat w_1$ is $z^*$-independent);
in $x$-momentum the transport term $\hat w_1\partial_{z^*}u_0$ vanishes
because $u_0$ is $z^*$-independent, leaving the 2D Euler $x$-momentum —
satisfied; likewise energy and species.  Impermeability holds by
construction; **no other wall condition exists at this order** (the
tangential-traction condition is carried by terms of size
$1/\mathrm{Re}_{\rm eff}\prec\varepsilon^2$: the slip BC *degenerates*, as
it must for Euler).  Uniqueness of classical solutions of compressible Euler
with slip walls (Schochet 1986; Kato 1975/Majda 1984) identifies the
candidate with the leading order on the smooth horizon. $\square$
**LEADING** (constants $\varepsilon$- and $\mathrm{Re}$-uniform under (T5)).

**(b3) Second order: linear transport, uniform bounds.**  Collecting
$\varepsilon^2$ in the Euler hierarchy gives a **linear symmetric-hyperbolic
system** for $U_2=(\rho_2,u_2,v_2,T_2,Y_{k,2})$, transported by the leading
flow:

$$
\partial_t U_2+\mathbf V_0\!\cdot\!\nabla_{2\mathrm D}U_2
+\hat w_1\,\partial_{z^*}U_2
+\mathcal L(U_0)\,U_2=S_2(x,y,z^*,t),
$$

where $\mathcal L(U_0)$ is the (zeroth-order) linearisation of the Euler
fluxes about $U_0$ and the forcing $S_2$ collects (i) quadratic
interactions of the expansion, (ii) geometric chain-rule terms, (iii) the
non-hydrostatic pressure $p_2$, obtained from $z$-momentum at
$\varepsilon^2$: $\partial_{z^*}p_2=-\rho_0\,\mathrm D\hat w_1/\mathrm Dt$
— **quadratic in $z^*$** since $\hat w_1$ is linear.  Structure: every
$z^*$-dependence in coefficients and forcing enters through $\hat w_1\propto
z^*$ and $\partial_{z^*}$ (the operator $z^*\partial_{z^*}$ preserves
polynomial degree), so if $U_2(0)$ is polynomial in $z^*$ of degree $\le2$
(Lemma 5.2: degree exactly $\le2$; class $\mathcal C$: degree $0$), then
$U_2(t)$ **remains polynomial of degree $\le2$ in $z^*$** for all smooth
time.  Linear well-posedness gives, with $K=O(\lVert\nabla\mathbf
V_0\rVert_\infty+\lVert\mathcal L\rVert_\infty)$ from (T5),

$$
\lVert U_2(t)\rVert\ \le\ e^{Kt}\Big(\lVert U_2(0)\rVert+\int_0^t
e^{-Ks}\lVert S_2(s)\rVert\,\mathrm ds\Big),
$$

**independent of $\mathrm{Re}_{\rm eff}$** — the uniformity statement that
regime (a) could not provide (there, $u_2\propto\mathrm{Re}$; here, no
Reynolds enters).  This closes, at the level of a complete derivation, the
historical sketch for regime (b). $\square$ **LEADING.**

**(b4) Status of the viscous geometric sources in (b).**  They are *not*
produced by the (b)-hierarchy (size $1/\mathrm{Re}_{\rm eff}\prec
\varepsilon^2$); they are **retained by the truncation rule (M3)**, being
$\varepsilon^0\cdot g(\mathrm{Re})$ terms constructed in regime (a), and are
asymptotically inert in (b), becoming leading wherever
$\mathrm{Re}_{\rm eff}$ drops locally (turbulent zones).  At Euler order the
momentum wall term reduces to the pressure thrust $h_x\,p$ (Shapiro) —
consistent with §7.3.

**(b5) The relaxation dichotomy (honest, load-bearing for (L5)).**  In
regime (a), spurious initial $z$-structure decays by transverse diffusion on
$t_{\rm relax}\sim h^2/\nu=\varepsilon^2\mathrm{Re}\,(L/U)\ll L/U$: the
model self-heals.  In regime (b) **no mechanism removes $z$-structure**:
Euler transports it indefinitely.  The legitimacy of the truncation in the
RDE regime therefore rests **entirely** on the data hypothesis (D2) — which
is why class $\mathcal C$ (zero initial out-of-plane gradients) is the
certified envelope, and why $z$-structured injection (discrete jets across
the gap) is a *permanent*, not transient, violation (§10.3).

### 5.4 Lemma (kinematic compatibility — slope forced by the model's own identities)

Under (M1) (fields vary on $L$) the bulk continuity balance gives
$|\partial_zw|=O(U/L)$, $z$-independent at leading order; the exact identity
of Lemma 3.5 gives $\langle\partial_zw\rangle_z=-\mathbf
V\!\cdot\!\mathbf G$ with $|\mathbf V\!\cdot\!\mathbf G|\le U|\mathbf G|
\le U\sigma/(\varepsilon L c_0/1)$ — more sharply $|\mathbf
V\!\cdot\!\mathbf G|\le U\,\sigma/h_{\min}$.  Mutual consistency requires

$$
\frac{U\sigma}{\varepsilon L}\lesssim\frac UL
\quad\Longleftrightarrow\quad
\boxed{\;\sigma=O(\varepsilon)\;}.
$$

Conversely, if $\sigma/\varepsilon\to\infty$, no flow with field scale $L$
can satisfy both continuity and impermeability: either the field scale
collapses to $L_h\ll L$ (contradicting (M1)) or $z$-layers form
(contradicting (M1)/(T3)).  Hence the slope clause of (G1) is not an
optional regularity refinement: it is the **unique consistency condition
between continuity and impermeability** — the formal content of the
intuition "the model is one of small out-of-plane variations". $\square$
**LEADING.**

---

## 6. Depth-averaged strain, vorticity and stress

### 6.1 Strain and vorticity (from §5 profiles)

With $u,v$ $z$-independent and $w=-(\mathbf V\!\cdot\!\mathbf G)z$ at
leading order:

$$
S_{11}=\partial_xu,\quad S_{22}=\partial_yv,\quad
\boxed{S_{33}=\partial_zw=-\mathbf V\!\cdot\!\mathbf G},\quad
S_{12}=\tfrac12(\partial_yu+\partial_xv),
$$
$$
S_{13}=-\tfrac z2\,\partial_x(\mathbf V\!\cdot\!\mathbf G),\quad
S_{23}=-\tfrac z2\,\partial_y(\mathbf V\!\cdot\!\mathbf G)\quad
(z\text{-odd},\ \langle S_{i3}\rangle_z=0),
$$

$\mathrm{tr}\,S=\nabla_{2\mathrm D}\!\cdot\!\mathbf V-\mathbf
V\!\cdot\!\mathbf G$.  Vorticity: $\omega_3=\partial_xv-\partial_yu$
$z$-independent; $\omega_1,\omega_2$ $z$-odd, so
$\langle|\boldsymbol\omega|^2\rangle_z=\omega_3^2+O(\varepsilon^2)$ — no
leading-order area-variation correction to vorticity-based closures.
**LEADING.**  The SST invariant:

$$
\langle|S|^2\rangle_z=2\sum_{ij}S_{ij}^2
=\underbrace{2\big[(\partial_xu)^2+(\partial_yv)^2\big]
+(\partial_yu+\partial_xv)^2}_{|S|^2_{2\mathrm D}}
+\,2\,(\mathbf V\!\cdot\!\mathbf G)^2+O(\varepsilon^2),
$$

so the turbulence production carries the **same-order** (not subleading)
component $2\mu_t(\mathbf V\!\cdot\!\mathbf G)^2$.  (SymPy twin F-3.)

### 6.2 Definition A versus Definition B

Two inequivalent reductions of the stress are conceivable:
$\mathbf A:\ T_{ij}:=\langle\tau^{3\mathrm D}_{ij}\rangle_z$;
$\mathbf B:\ T_{ij}:=\tau_{ij}(\langle u\rangle,\langle v\rangle,\langle
w\rangle{=}0)$.  They differ only through $S_{33}$: under A,
$\langle S_{33}\rangle_z=[w]/h=-\mathbf V\!\cdot\!\mathbf G$ (Lemma 3.5,
trace-exact); under B, $S_{33}=0$.  The difference is $O(U/L)$ — **not**
subleading.

**Proposition 6.2.** If the 2D model is the $z$-integration of the 3D
conservation laws, **Definition A is mathematically forced**: the integrated
viscous fluxes contain $\langle\tau_{ij}\rangle_z$ and the wall jumps
$[\tau_{i3}]$, never $\tau_{ij}(\langle\cdot\rangle)$.  Declaring B
(equivalently $S_{33}=0$) produces four $O(1)$ inconsistencies: stress error
$\tfrac43\mu(\mathbf V\!\cdot\!\mathbf G)$ in $\tau_{33}$; loss of the mass
source $\rho(\mathbf V\!\cdot\!\mathbf G)$ and of exact $h$-weighted
conservation; momentum and energy conservation violated at
$O(\mathbf V\!\cdot\!\mathbf G)$. $\square$ **LEADING.**

### 6.3 Theorem 6.3 (depth-averaged Newton stress)

Under (G1),(P1),(P2),(D1),(M1), with the Reynolds–Favre factorisation
$\langle\mu X\rangle_z=\langle\mu\rangle_z\langle X\rangle_z+O(\varepsilon^4)$
(products of $O(\varepsilon^2)$ fluctuations):

$$
\begin{aligned}
T_{11}&=2\mu\,\partial_xu-\tfrac23\mu\big(\nabla_{2\mathrm
D}\!\cdot\!\mathbf V-\mathbf V\!\cdot\!\mathbf G\big), &
T_{22}&=2\mu\,\partial_yv-\tfrac23\mu\big(\nabla_{2\mathrm
D}\!\cdot\!\mathbf V-\mathbf V\!\cdot\!\mathbf G\big),\\
T_{12}&=\mu\,(\partial_yu+\partial_xv), &
T_{33}&=-\tfrac43\mu\,(\mathbf V\!\cdot\!\mathbf G)
-\tfrac23\mu\,\nabla_{2\mathrm D}\!\cdot\!\mathbf V,
\end{aligned}
$$

with $\mathrm{tr}\,T=0$ for **any** $\langle\mathrm{tr}S\rangle_z$
(Stokes: $2\mu+3\lambda=0$).  RANS (Boussinesq): same with
$\mu\to\mu_{\rm eff}$ and $-\tfrac23\rho k$ added to each diagonal; the $k$
contribution cancels in every difference $T_{ii}-T_{jj}$ and contributes
$+\tfrac23\rho k(\mathbf V\!\cdot\!\mathbf G)$ to the energy source through
$-\mathbf V\!\cdot\!(\boldsymbol\tau_{\rm tot}\!\cdot\!\mathbf G)$ —
counted **once** (no $P^*$ double count; SymPy twin F-4). $\square$
**LEADING** ($O(\varepsilon^4)$ closure error). (SymPy twin F-2.)

---

## 7. $z$-integration of the conservation laws

Each 2D equation is the master identity (3.6) applied to the corresponding
3D flux, with wall terms evaluated under (P2)–(P4).  All equations are
stated in the converged convention: thermodynamic pressure $p$ in the
convective fluxes, Anderson-coherent signs, $\boldsymbol\tau_{\rm tot}$ with
$\mu_{\rm eff}$ and $-\tfrac23\rho k$ diagonals.

### 7.1 Continuity — exact $h$-weighted conservation

$\mathbf F=\rho\mathbf V$, $\mathbf F\!\cdot\!\hat{\mathbf n}|_\pm=0$
exactly (P2); the Leibniz trace terms cancel **identically** against
$[\rho w]$ (Lemma 3.4):

$$
\boxed{\;\partial_t\rho+\partial_x(\rho u)+\partial_y(\rho v)
=\rho\,(\mathbf V\!\cdot\!\mathbf G)\;}
\qquad\Longleftrightarrow\qquad
\partial_t(h\rho)+\partial_x(h\rho u)+\partial_y(h\rho v)=0 .
$$

The $h$-weighted form, in the averaged-flux variables
$\langle\rho u\rangle$, is **EXACT** (no expansion used); the per-area form
with model variables carries the $O(\varepsilon^4)$ product-closure error.

### 7.2 The geometric source is an algebraic residue, not a wall flux

$\mathbf G\!\cdot\!\langle\mathbf F\rangle_z$ arises from converting the
globally conservative $h$-weighted form to per-area form — the 2D extension
of quasi-1D duct flow.  The walls are impermeable: nothing crosses them
except momentum (pressure + normal viscous traction).

### 7.3 $x$-momentum and the pressure cancellation

$\mathbf F^{\rho u}=\rho u\mathbf V+p\,\mathbf e_x-\boldsymbol\tau\!\cdot\!
\mathbf e_x$.  Wall traction (P2 + §4.6):
$\mathbf F^{\rho u}\!\cdot\!\hat{\mathbf n}|_\pm
=\tfrac{h_x}{2}(\tau_{33}-p)\big|_\pm+O(\sigma^3)$.  Master identity, trace
$\to$ average replacement ($O(\sigma\varepsilon^2)$, (D1)+(M1)), and the
**exact algebraic cancellation** of $p\,G_x$ between
$\mathbf G\!\cdot\!\langle\mathbf F\rangle$ and the wall term:

$$
\boxed{\;\partial_t(\rho u)+\partial_x(\rho u^2+p-\tau_{\mathrm{tot},11})
+\partial_y(\rho uv-\tau_{\mathrm{tot},12})
=\rho u(\mathbf V\!\cdot\!\mathbf G)
+G_x(\tau_{\mathrm{tot},33}-\tau_{\mathrm{tot},11})
-G_y\,\tau_{\mathrm{tot},12}\;}
$$

Inviscid limit: RHS $\to\rho u(\mathbf V\!\cdot\!\mathbf G)$ (Shapiro).
Viscous residue: $T_{33}-T_{11}=-2\mu_{\rm eff}\big(\mathbf
V\!\cdot\!\mathbf G+\partial_xu\big)$; the $-\tfrac23\rho k$ terms cancel in
the difference. **LEADING.** (SymPy twins F-1, harness §8.3.)

### 7.4 $y$-momentum

$$
\boxed{\;\partial_t(\rho v)+\partial_x(\rho uv-\tau_{\mathrm{tot},12})
+\partial_y(\rho v^2+p-\tau_{\mathrm{tot},22})
=\rho v(\mathbf V\!\cdot\!\mathbf G)-G_x\tau_{\mathrm{tot},12}
+G_y(\tau_{\mathrm{tot},33}-\tau_{\mathrm{tot},22})\;}
$$

### 7.5 Total energy

$E_{\rm res}=e+\tfrac12(u^2+v^2)$ (the $\tfrac12\langle w^2\rangle_z=
\tfrac1{24}h^2(\mathbf V\!\cdot\!\mathbf G)^2$ correction is
$O(\varepsilon^2)$, discarded — §9).  The wall flux vanishes **exactly**:
impermeability kills convection and pressure work; slip makes
$\boldsymbol\tau\!\cdot\!\hat{\mathbf n}\parallel\hat{\mathbf n}$ so
$\mathbf V\!\cdot\!(\boldsymbol\tau\!\cdot\!\hat{\mathbf n})=
\sigma_n(\mathbf V\!\cdot\!\hat{\mathbf n})=0$; (P3)–(P4) kill conduction
and diffusion.  Hence (with $\Phi_k$ the $k$-coupling, $\mathbf q_{\rm tot}
=-(\kappa_l+\mu_tc_p/\mathrm{Pr}_t)\nabla T$, $\mathbf j^{\rm eff}_k$
Hirschfelder-corrected):

$$
\boxed{\;\partial_t(\rho E_{\rm res})+\partial_xF^E_x+\partial_yF^E_y
=(\rho E_{\rm res}+p)(\mathbf V\!\cdot\!\mathbf G)
-\mathbf V\!\cdot\!(\boldsymbol\tau_{\rm tot}\!\cdot\!\mathbf G)
+\mathbf G\!\cdot\!\mathbf q_{\rm tot}
+\sum_k h_k(\mathbf G\!\cdot\!\mathbf j^{\rm eff}_k)+\Phi_k\;}
$$

The $h$-weighted energy identity in averaged fluxes is **EXACT** (wall flux
identically zero); the per-area closure carries $O(\varepsilon^3)$ (viscous
commutation) and $O(\varepsilon^4)$ (products) — this resolves the
historical "ESATTA vs leading-order" labelling ambiguity (Appendix C, item
C-6).

### 7.6 Species and turbulence scalars

$$
\boxed{\;\partial_t(\rho Y_k)+\nabla_{2\mathrm D}\!\cdot\!\big(\rho\mathbf
VY_k+\mathbf j^{\rm eff}_k\big)
=\rho Y_k(\mathbf V\!\cdot\!\mathbf G)+\mathbf G\!\cdot\!\mathbf
j^{\rm eff}_k+\omega_k\;}
$$

with $\sum_k\mathbf j^{\rm eff}_k=0$ exactly (Hirschfelder, any $D_k$);
for $\phi\in\{k,\omega,\tilde\nu\}$ with diffusivity $D_\phi$ and
$\partial_n\phi=0$ at the walls (P4):

$$
\boxed{\;\partial_t(\rho\phi)+\nabla_{2\mathrm D}\!\cdot\!\big(\rho\mathbf
V\phi-D_\phi\nabla\phi\big)
=\rho\phi(\mathbf V\!\cdot\!\mathbf G)-D_\phi(\mathbf
G\!\cdot\!\nabla\phi)+S^\phi\;}
$$

— the diffusive geometric source carries the **negative** sign, identical in
structure to $+\mathbf G\!\cdot\!\mathbf q$ with $\mathbf q=-\kappa\nabla T$
(SymPy twin: harness V3-equivalent).  Production:
$P_k^{\rm Q2D}=\mu_t\big[|S|^2_{2\mathrm D}+2(\mathbf V\!\cdot\!\mathbf
G)^2\big]$ — a leading component, not a correction (§6.1).

### 7.7 Chemistry closure with explicit Zel'dovich constant

**Lemma 7.7.**  Let $\omega_k\in C^2$ in $(T,Y,\rho)$ and let the fields
satisfy (M1) with $T=T_0+\varepsilon^2T_2+O(\varepsilon^4)$ etc.  Then

$$
\big\langle\omega_k(T,Y,\rho)\big\rangle_z-\omega_k\big(\langle
T\rangle,\langle Y\rangle,\langle\rho\rangle\big)
=\frac{\varepsilon^4}{2}\Big[\partial^2_{TT}\omega_k\,
\mathrm{Var}_{z^*}(T_2)+\cdots\Big]+O(\varepsilon^6),
$$

(first-order terms cancel identically between the two Taylor expansions;
only **second derivatives times fluctuation variances** survive).  For
Arrhenius $\omega\propto\exp(-E_a/RT)$,
$\partial^2_{TT}\omega=\omega\,\big[\beta^2/T_0^2\big]\,(1+O(1/\beta))$ with
$\beta=E_a/(RT_0)$, so the **relative** closure error is

$$
\frac{|\langle\omega\rangle-\omega(\langle\cdot\rangle)|}{|\omega|}
\;\le\;C\,\beta^2\,\varepsilon^4\,
\Big(\frac{\lVert T_2\rVert}{T_0}\Big)^2+\text{(lower order in }\beta) .
$$

The closure is subleading **iff $\beta^2\varepsilon^4\ll1$** — hypothesis
(M2-iv).  At $\varepsilon=0.3$: benign for $\beta\lesssim10$ (error
$\lesssim1\%$), marginal for $\beta\sim30$ (high-activation-energy
detonability limits).  The historical statement "$O(\varepsilon^4)$"
concealed the $\beta^2$ amplification (Appendix C, item C-7). $\square$
**LEADING.**

---

## 8. The closed model and the Truncation Admissibility Theorem

### 8.1 The closed Q2D system

Equations 7.1, 7.3–7.6 with the tensor of §6.3, $P^{\rm Q2D}_k$, the SST
cross-diffusion in averaged variables, and $\Phi_k=(-P_k+\beta^*\rho\omega
k)\,\mathbb 1_{k\text{-coupling}}$ constitute the model — identical to the
consolidated historical system ("§N16") with all errata absorbed.  Per
(M3), the system retains all $\varepsilon^0$ terms uniformly over regimes
(a)/(b).

### 8.2 Theorem (truncation admissibility)

**Discarding the second order is legitimate iff the following seven
conditions hold.**  Under them, on the smooth horizon,

$$
\big\lVert\langle f^\varepsilon\rangle_z-f^{\rm Q2D}\big\rVert
= O(\varepsilon^2)\qquad\text{(formal; functional version: §12)} .
$$

| # | Condition | Content | Status here |
|---|---|---|---|
| **L0** | Parity | (G1)+(D1) $\Rightarrow$ odd orders vanish (Thm 5.1): first discarded order is the **second** — without symmetry the model is only $O(\varepsilon)$ | proved (uniqueness-local + (T1)) |
| **L1** | Admissibility of the discarded order | $f_2$ exists and is bounded along model solutions: regime (a) by Fredholm (Thm 5.6); regime (b) by linear transport (§5.3-b3) | proved per regime |
| **L2** | Uniformity of constants | effective parameter $\delta=\varepsilon^2\cdot C$ with $C$ tracked: $C=O(1)$ in (b); $C=O(\mathrm{Re}h^2|G_0|)$ in (a); $C\propto\beta^2$ for chemistry; $C\propto(\sigma/\varepsilon)^2$ for geometry; transition zone excluded by (M2-ii) | proved under (T3)+(T5) |
| **L3** | Observable non-degeneracy | legitimacy is relative to the **question**: valid for observables $O(1)$ at leading order; for $O(\varepsilon^2)$ observables (gap-normal waves, $z$-profiles, transverse energy) the model returns **zero**, not a 2 % error | declared (§9.3) |
| **L4** | Time horizon | error grows as $\varepsilon^2e^{Kt}$ (Gronwall): valid for $T=O(L/U)$; multi-rotation validity requires the leading dynamics to be **attracting** (errors saturate) — stability postulate, validated a posteriori | postulate (§10.5) |
| **L5** | Prepared data | (D2): $z$-structure of data $O(\varepsilon^2)$ (class $\mathcal C$: zero).  Regime (a) self-heals in $t\sim h^2/\nu$; regime (b) has **no relaxation mechanism** — (D2) is the sole defence | hypothesis; dichotomy §5.3-b5 |
| **L6** | Regular perturbation | the discarded terms carry no principal derivatives and do not change the type of the system (hyperbolic–parabolic structure preserved): truncation is a **regular**, not singular, perturbation | verified by inspection of §9.1 list |
| **L7** | Complete bookkeeping | every discarded term is provably $O(\varepsilon^2)$ under (G1)+$\Pi$ — in particular the slip-projection terms are $O(\sigma^2)=O(\varepsilon^2)$ **only because of the slope clause** | proved (§4.6, §9.1) |

*Proof.* L0 = Theorem 5.1.  L1 = Theorem 5.6 / §5.3-b3.  L2: collect the
constants of §5.2(a5), §5.3(b3), Lemma 7.7, §4.6.  L3, L5, L6 by the cited
inspections.  L4: the energy/relative-entropy estimate of §12 gives the
$e^{Kt}$ envelope.  L7: tables §4.3–4.6 and the discarded list §9.1.
$\square$

**Negative counterpart.**  The truncation is **not** legitimate — and no
error-bar smaller than $O(1)$ may be claimed — whenever any condition
fails: $z$-structured data in regime (b) (L5); observables at
$O(\varepsilon^2)$ (L3); $\sigma\gg\varepsilon$ (L7/CE-1); zonal
$\mathrm{Re}_{\rm eff}\varepsilon^2=\Theta(1)$ (L2); secular horizons
without attraction (L4); broken $z$-symmetry (L0 $\to$ first-order model);
front interiors (M1 locally violated — captured, not resolved).

---

## 9. Loss budget — what the truncation discards, and its physics

### 9.1 Itemised discarded terms

| Discarded term | Order (two-parameter) | Physics lost | When it matters |
|---|---|---|---|
| Non-hydrostatic $p_2$ ($\partial_{z^*}p_2=-\rho_0\mathrm D\hat w_1/\mathrm Dt$) | $\varepsilon^2$ | **gap-normal pressure dynamics**: transverse acoustics, wall-to-wall shock reflection | RDE: the single most significant structural loss — gap-resolved wave systems are unrepresentable *by construction* |
| $\tfrac12\langle w^2\rangle_z=\tfrac1{24}h^2(\mathbf V\!\cdot\!\mathbf G)^2$ in $E$ | $\varepsilon^2$ | transverse kinetic energy | strong contractions |
| $\langle\tau_{13}w\rangle$ in $F^E$ | $\varepsilon^2/\mathrm{Re}$ | transverse shear work | negligible along $\Pi$ |
| Profile-defect correlations ($\langle\rho u^2\rangle-\langle\rho\rangle\langle u\rangle^2$, $\langle\mu S\rangle-\langle\mu\rangle\langle S\rangle$) | $\varepsilon^4$ | shape factors; **the plug profile kills them at leading order** — this is *why* free-slip yields a clean model (no $6/5$ factor); at next order they return as Taylor-dispersion-like transport | only beyond $O(\varepsilon^2)$ |
| Slip-projection corrections ($(1-h_x^2/4)$, $\tau_{23}h_xh_y/4$) | $\sigma^2$ relative | wall-traction misprojection, grows with slope | throats, sharp area changes; $=O(\varepsilon^2)$ only under (G1) |
| $S_{13},S_{23}$ content of $\langle|S|^2\rangle$ | $\varepsilon^2$ | transverse-shear turbulence production | strong $\nabla(\mathbf V\!\cdot\!\mathbf G)$ |
| Chemistry de-correlation | $\beta^2\varepsilon^4$ | $\langle\omega\rangle\ne\omega(\langle\cdot\rangle)$ | detonability margins (M2-iv) |

### 9.2 Excluded regimes (lost applicability, not lost terms)

Corrugated geometry ($\sigma\gg\varepsilon$: homogenisation); viscous
shallow water ($\mathrm{Re}\varepsilon^2=\Theta(1)$); low-Mach and
hypersonic paths; $\mathrm{Da}\to\infty$ front interiors; $z$-structured
forcing in regime (b).

### 9.3 Observable-relative validity (L3 made concrete)

Valid questions (leading-order observables): depth-averaged propagation
speed, pressure/thrust history, mean temperature and composition fields,
integrated performance.  Invalid questions (the model returns zero, not an
approximation): gap-normal wave timing, $z$-profiles of any field,
transverse turbulence anisotropy in $z$, wall-normal heat flux structure.
If $z$-symmetry is broken (curved mid-surface, one-sided injection), all
error statements degrade from $O(\varepsilon^2)$ to $O(\varepsilon)$ (L0).

---

## 10. Validity envelope (quantitative, class $\mathcal C$)

### 10.1 Geometry

$\varepsilon\in[0.05,0.30]$ recommended ($\varepsilon^2$: 0.25–9 % nominal
model error); slope $\sigma\le M\varepsilon$ with moderate $M$ —
**monitorable at runtime**: $\sigma_{\rm loc}=h\,|\mathbf G|$ is exactly
computable cell-wise from the stored geometry
(`Cross_Section` $\times$ `area_geom_term`); recommended warning threshold
$\sigma_{\rm loc}>0.3$.  (This re-instates, on a corrected basis, the
validity monitor historically removed with flawed reasoning — Appendix C,
item C-5.)

### 10.2 Reynolds zoning

$\delta_{\rm local}=\max\{\varepsilon^2,\ (\mathrm{Re}_{\rm eff}
\varepsilon^2)^{-1}\}$ per zone; regime (b) for laminar/moderately turbulent
RDE zones ($\mathrm{Re}_{\rm eff}\varepsilon^2\sim10$–$10^4$); local
transition pockets ($\mu_t/\mu_l\sim10^4$) contribute $\delta_{\rm
local}\lesssim10\%$, confined and declared.

### 10.3 Data (class $\mathcal C$)

Initial fields with **zero out-of-plane gradients** (strong (D2)): the
initial fluctuation vanishes identically; in regime (b) this is the *only*
mechanism guaranteeing (L5) for all time (§5.3-b5).  Any $z$-structured
in-plane forcing (discrete injector jets) is a **permanent** (L5) violation
whose induced error must be assessed separately (3D or experimental
anchoring).

### 10.4 Mach

Certified path: $\mathrm{Ma}=\Theta(1)$ (transonic to moderately
supersonic; CJ cells locally violate (M1) and are *captured*).
$\mathrm{Ma}<0.3$ and $\mathrm{Ma}>5$ sustained lie off $\Pi$-(iii):
order tables must be re-derived (low-Mach pressure scaling; compressibility
corrections for the closures).

### 10.5 Time horizon (honest restatement — supersedes the historical $T_{\rm loc}$ estimate)

The local Cauchy theory (Matsumura–Nishida / Kato–Majda–Schochet) guarantees
the parity and hierarchy constructions on $T_{\rm loc}\sim C_0/K_\infty$.
With detonation-strength gradients the Sobolev size is $K_\infty\sim10^3$
(in $T_{\rm ref}=L/U$ units), giving $T_{\rm loc}\sim10^{-3}\,T_{\rm ref}$ —
**microseconds, far below one RDE rotation**.  (The historical claim
"$T_{\rm loc}\sim10$ ms $\sim$ 10 rotations" contained a $10^3$ arithmetic
inconsistency and is withdrawn; Appendix C, item C-8.)  Consequently the
multi-rotation applicability of the model does **not** rest on Cauchy
theory.  It rests, declaredly, on:
(i) the URANS-mean character of the fields (mean gradients are
turbulence-limited, lowering $K_\infty$ by orders of magnitude for the
*mean* problem); (ii) the **attractor postulate** (L4): for stable
limit-cycle operation the leading dynamics contracts errors, so the
$\varepsilon^2 e^{Kt}$ envelope saturates instead of compounding;
(iii) a-posteriori validation (§14).  This is a postulate-backed,
validation-anchored statement — not a theorem — and is declared as such.

### 10.6 The idealisation gap (P2–P4 vs real walls)

All error statements above are **relative to the idealised free-slip 3D
reference**.  Real RDE gap walls are no-slip and conducting: boundary layers
of thickness $\delta_{\rm BL}$ remove momentum
($\sim\delta_{\rm BL}/h$ relative) and extract heat.  For the model to be
*physically* (not only asymptotically) accurate one needs
$\delta_{\rm BL}/h\lesssim\varepsilon^2$ — typically **not** satisfied in
small gaps; then the slip-idealisation error, not the $O(\varepsilon^2)$
truncation error, dominates the comparison with experiment.  Semi-empirical
wall closures (friction/heat source terms) are the declared remedy and are
**outside** this derivation. **DECLARED MODEL-FORM ERROR.**

---

## 11. Global conservation and well-posedness

Integrating the $h$-weighted form over $\Omega_{2\mathrm D}$ with
impermeable in-plane boundaries:

- **Mass:** $\dfrac{\mathrm d}{\mathrm dt}\displaystyle\int h\rho =0$ — **EXACT**.
- **Energy:** $\dfrac{\mathrm d}{\mathrm dt}\displaystyle\int h\rho E_{\rm res}=\int h\,\Phi_k$
  — exact in averaged fluxes (wall flux identically zero, §7.5); closure
  error $O(\varepsilon^3)$ in model variables.
- **Species:** $\dfrac{\mathrm d}{\mathrm dt}\displaystyle\int h\rho Y_k=\int h\,\omega_k$.
- **$x$-momentum:**
  $\dfrac{\mathrm d}{\mathrm dt}\displaystyle\int h\rho u
  =\int \big[h_x\,(p-\tau_{\mathrm{tot},33})-h_y\,\tau_{\mathrm{tot},12}\big]$
  — pressure + viscous thrust on the inclined walls (Shapiro extended to
  2D).  The $-h_y\tau_{12}$ term, omitted in a historical statement, is
  restored (Appendix C, item C-9); the $y$-analogue is
  $h_y(p-\tau_{\mathrm{tot},33})-h_x\tau_{\mathrm{tot},12}$.

The truncated system preserves the hyperbolic–parabolic principal structure
of 2D compressible (U)RANS — the geometric sources are zeroth-order — so the
standard local well-posedness theory applies unchanged (L6).

---

## 12. Convergence to the 3D reference: status

**CONJECTURE** (declared).  $\lVert\langle f^\varepsilon\rangle_z-
f^{\rm Q2D}\rVert_{L^2(\Omega_{2\mathrm D}\times[0,T])}=O(\varepsilon^2)$ as
$\varepsilon\to0$ along $\Pi$.  Supporting structure: (i) the constructive
hierarchy of §5 (orders exist, are bounded, and their solvability conditions
are the model — Theorems 5.1/5.6, §5.3); (ii) a modulated relative-entropy
programme (kinetic + Boltzmann-$\psi$ thermal + Bregman pressure terms,
Gronwall with well-prepared data $H(0)=O(\varepsilon^4)$), valid for
$\mu(T)$ where BD-entropy is not (the BD relation requires $\mu(\rho)$);
its residual bound $\lVert{\rm res}\rVert_{H^{-1}}\le C\varepsilon^2$ is
established **under (M1)+(T5), not unconditionally** — a scholarly outline,
suitable for a dedicated mathematical companion, not a theorem of this page.
References: Dafermos (1979; 2016 Ch. 5), Germain (2011), Feireisl–Jin–
Novotný (2012), Feireisl–Novotný (2012) ($\mu(T)$ heat-conducting case),
Saint-Raymond (2009), Brenier (2000).

---

## 13. Registry of legitimate objections (Q&A of record)

**Q1 — "The model's variables are averages; how can it use pointwise wall
information?"**  It doesn't: the wall traces it uses are determined by the
BCs (P2)–(P4) as functions of the averages plus $O(\varepsilon^2)$
(Lemma 3.4–3.5).  The closure problem is genuine (§3); BC-consistency is the
route chosen, and §2.6/CE-1 delimits where it works.

**Q2 — "Is $S_{33}=-\mathbf V\!\cdot\!\mathbf G$ an assumption?"**  No: it
is the depth average of an **exact kinematic identity** (Lemma 3.5) forced
by impermeability; Definition A is then mathematically required
(Prop. 6.2), and declaring $S_{33}=0$ produces four $O(1)$ inconsistencies.

**Q3 — "Why is the error $O(\varepsilon^2)$ and not $O(\varepsilon)$?"**
Because the first-order coefficient vanishes **by symmetry** (Theorem 5.1)
— a theorem under (G1)+(D1), not bookkeeping.  Break the symmetry and the
model is first-order (L0).

**Q4 — "Isn't gentle slope implicit in thinness?"**  No: CE-1 is thin and
steep, and the model fails on it.  The slope clause of (G1) is independent,
and is *forced* by the model's own kinematics (Lemma 5.4).  Historically it
was intended verbally, used everywhere, and never declared (App. C, C-1).

**Q5 — "Under which hypotheses is discarding the second order
legitimate?"**  Exactly the seven conditions L0–L7 of Theorem 8.2 — no
fewer.  Each has a failure mode listed in the negative counterpart.

**Q6 — "The rigorous Fredholm closure lives in regime (a), but RDE lives in
regime (b)."**  Correct historically; resolved here: §5.3 derives regime
(b) per-equation (Euler propagation + linear-transport second order with
Re-uniform bounds), with the slip BC correctly degenerating at Euler order.
The two regimes yield the **same** closed system; the transition zone is
excluded (M2-ii) and error-bounded zonally (§10.2).

**Q7 — "In regime (b) you retain viscous sources smaller than terms you
discard."**  Yes — deliberately, by the truncation rule (M3): truncation is
in $\varepsilon$ only, uniformly over the Reynolds regimes, because
$\mathrm{Re}_{\rm eff}$ varies by orders of magnitude within one domain.
The retained terms are inert where small and leading where
$\mathrm{Re}_{\rm eff}$ drops; discarding them would be a *zonal* error.

**Q8 — "How long in time is the model justified?"**  By Cauchy theory:
microseconds at detonation gradients (§10.5 — the historical 10 ms claim is
withdrawn).  Operationally: by the attractor postulate (L4) + a-posteriori
validation.  Declared as such.

**Q9 — "What about the real walls?"**  The $O(\varepsilon^2)$ claim is
relative to the idealised slip reference; the idealisation gap
($\delta_{\rm BL}/h$) is a separate declared model-form error (§10.6) —
for small gaps it, not the truncation, dominates comparison with
experiment.

**Q10 — "Stokes hypothesis?"**  $\lambda=-\tfrac23\mu$ is assumed (P1);
for polyatomic mixtures with significant bulk viscosity the diagonal stress
acquires $\mu_v(\nabla\!\cdot\!\mathbf V)\delta_{ij}$ — a model-form choice
inherited from the solver family, orthogonal to the depth-averaging.

**Q11 — "Are the chemistry averages closed?"**  Yes at
$O(\beta^2\varepsilon^4)$ (Lemma 7.7) — with the Zel'dovich amplification
now explicit and bounded by (M2-iv).

**Q12 — "Favre or Reynolds averages?"**  Indistinguishable at the model's
order: $\tilde f-\langle f\rangle_z=\langle\rho'f'\rangle/\langle\rho\rangle
=O(\varepsilon^4)$ (products of $O(\varepsilon^2)$ fluctuations,
Cauchy–Schwarz).  The FV code's conservative reconstruction (Favre) is
consistent with the Reynolds-averaged derivation without modification.

**Q13 — "Is the wall treatment the same physics as the CFD viscous-symmetry
BC, and is it the best one can do?"**  Yes, verifiably: the lateral
`BC_Symmetry_Visc` realises the mirror formulation (F3) discretely — ghost
$\mathbf V_2=\mathbf V_1-2(\mathbf V_1\!\cdot\!\hat{\mathbf n})\hat{\mathbf
n}$ gives $\partial_nV_t=0$, $\partial_tV_n=0$, hence a **purely normal**
face traction $\tau_{nn}\hat{\mathbf n}$ (with the Definition-A tensor:
$\mathbf V\!\cdot\!\mathbf G$ is wired into the face stress), zero mass /
energy ($\mathbf V_t\perp\tau_{nn}\hat{\mathbf n}$, algebraically) / scalar
fluxes, plus the RSM odd-$R_{12}$ block (in-plane mirror only; under the
$z$-mirror the transported $R_{11},R_{12},R_{22},R_{33}$ are all even).  The
area-variation sources are the depth-integrated transposition of the *same*
continuum BC (the F1$\Leftrightarrow$F3 equivalence), transmitting exactly
$p$ and $\tau_{33}$ ($=\tau_{nn}+O(\sigma^2)$) through the tilt — the part a
Definition-B model would lose at $O(1)$.  "Better", rigorously, means
exactly two things: (i) the slope-exact traction projection (recovers the
$O(\sigma^2)$ terms of §4.6 — relevant only at $\sigma_{\rm loc}$ peaks);
(ii) the wall acoustic response captured laterally by the mirror Riemann
solve, which out-of-plane is the non-hydrostatic $p_2$ — irreducibly
second-order physics, recoverable only by a dispersive second-order model
(Green–Naghdi analogue) or a $z$-resolved 3D solver.  At leading order along
$\Pi$ the present closure is unique (Prop. 6.2 + §5).

**Q14 — "A slip wall has no gradients, hence no viscous stress — so why do
the sources carry $\tau_{33}$?"**  Half of the gradients vanish, not all:
the mirror ghost is a **parity selector**.  Even quantities ($V_t$, $T$,
$Y_k$, $k$, $\omega$) have odd normal derivative $\Rightarrow$ zero at the
face: no friction, no heat, no diffusive exchange.  The **normal velocity is
odd**, so its normal derivative is even and survives: $V_n=0$ *at* the wall
but $\partial_nV_n\neq0$ (continuity: in-plane dilatation forces normal
compression), generating the nonzero viscous **normal** stress
$\tau_{nn}$ — discretely visible in `BC_Symmetry_Visc`, where
$\mathbf V_2-\mathbf V_1=-2(\mathbf V\!\cdot\!\hat{\mathbf n})\hat{\mathbf
n}\neq0$ and the residual receives `Flux(nu:nv)`$\,\neq0$ while
`Flux(np)`$=0$.  In the depth-averaged setting the surviving gradient is
exactly $\partial_zw=-\mathbf V\!\cdot\!\mathbf G$ (Lemma 3.5) and its
stress is exactly $\tau_{33}$; setting it to zero is Definition B with its
four $O(1)$ inconsistencies (§6.2).  Flat wall ($\nabla h=0$): the
surviving gradient and the transmission both vanish — slip on a flat wall
transmits nothing, and the sources correctly live only where the area
varies.  Ghost hierarchy by information content: **connection** ghosts =
real neighbour data (everything passes: full
$\boldsymbol\tau\!\cdot\!\hat{\mathbf n}$, heat, species); **mirror** ghosts
= the slip parity selector (only $p+\tau_{nn}$ passes); **zero-gradient
extrapolated** ghosts = transparent boundary (nothing physical imposed).

**Q15 — "Would a 3D simulation with a single profile-following cell across
the gap see the same thing?"**  Yes — structurally, not accidentally: the
FV cell average **is** the depth average, and the one-cell balance **is**
the master identity (§3.6) in discrete form.  Lateral face areas $\propto h$
produce the $h$-weighted in-plane fluxes; the tilted top/bot wall BC
produces the wall thrust $h_x(p-\tau_{nn})$ from face geometry alone; and
the mirror ghost across the **tilted** normal regenerates the surviving
gradient by itself:
$w^{\rm top}_{\rm ghost}=-2(\mathbf V\!\cdot\!\hat{\mathbf n})n_z=
-h(\mathbf V\!\cdot\!\mathbf G)$,
$w^{\rm bot}_{\rm ghost}=+h(\mathbf V\!\cdot\!\mathbf G)$, hence the
discrete $\partial_zw=-(\mathbf V\!\cdot\!\mathbf G)$ — Lemma 3.5
reconstructed from ghost geometry, with the same $\tau_{33}$ and Definition-A
tensor.  One-cell 3D is therefore the **conservative-form** discretisation
of the same model of which Q2D is the **source-form** discretisation: same
physics content, differences only at discretisation level ($h$-weighted
conservation machine-exact vs $O(\Delta x^2)$; Riemann wall pressure vs
$\langle p\rangle$, an $O(\varepsilon^2)$ transient).  Neither sees any
$z$-structure ($p_2$, gap acoustics): one cell carries the identical
$O(\varepsilon^2)$ truncation.  Practical corollary — the cleanest
validation design: (i) one-cell 3D vs Q2D = cross-verification of the two
discretisations (cheap, mutual); (ii) $z$-refinement in 3D
(1$\to$4$\to$16 cells across the gap) measures the **model error**
directly, and its scaling with $\varepsilon$ is the open Phase-C
$\varepsilon^2$ test.  Caveat: the 3D run must use the *viscous* symmetry
(mirror w.r.t. the true tilted normal); a purely inviscid slip wall drops
the $\tau_{nn}$ transmission.

---

## 14. Implementation mapping and verification status

| Formal result | Code site | Verification |
|---|---|---|
| $S_{33}=-\mathbf V\!\cdot\!\mathbf G$ (Lemma 3.5, §6.1) | `Strain_Tensor_Std` (optional `V_dot_G`) | symbolic + in-solver MMS (GOLD) |
| $T_{33}$, $\mathrm{Divervel}=\nabla_{2\mathrm D}\!\cdot\!\mathbf V-\mathbf V\!\cdot\!\mathbf G$ (Thm 6.3) | `Stress_Tensor_Std/_2eq` | symbolic (F-2) + MMS |
| $x/y$-momentum sources (§7.3–7.4) | `Source_Area_Visc_Blk` | GOLD (harness §8.3–8.4 + $\mu=20$ MMS order 2) |
| Energy source (§7.5) incl. $k$ single-count | `Source_Area_Visc_Blk` + `rhokdelta` | symbolic-exact (F-4, harness) |
| Species + Hirschfelder (§7.6) | `Source_Area_Visc_Blk` | symbolic-exact ($\Sigma_k j^{\rm eff}_k=0$, any $D_k$) |
| RANS $-D_\phi(\mathbf G\!\cdot\!\nabla\phi)$ (§7.6) | `*_Source_Diffusive_Q2D` per model | **panel-only — OPEN (Phase C MMS)** |
| Eulerian $\rho\phi(\mathbf V\!\cdot\!\mathbf G)$ | `Source_Area_Eul_Blk` | GOLD (nozzle Richardson $p=1.99$) |
| Face $\mathbf V\!\cdot\!\mathbf G$ | `V_dot_G_face_avg` $\to$ diffusive fluxes/BCs | structural + planar $G=0$ bit-exact |
| Slope monitor $\sigma_{\rm loc}=h|\mathbf G|$ (§10.1) | **to be implemented** (`Cross_Section`$\times$`area_geom_term`) | — |

Code-verification record: `docs/review/area_variation_vv_campaign.md`
(12/12 symbolic, planar-limit bit-exact, conservation $O(\Delta x^{3.2})$,
nozzle Richardson, in-solver viscous MMS order 2 at stressed source weight).
Physical validation against a 3D slip-wall reference ($\varepsilon^2$
scaling in $\varepsilon$) remains the open Phase-C item.

---

## Appendix A — Selected proofs

**A.1 (Lemma 3.2).** $\partial_x\!\int_{-h/2}^{h/2}f\,\mathrm dz
=\int_{-h/2}^{h/2}\partial_xf\,\mathrm dz+\tfrac{h_x}{2}(f|_++f|_-)$;
divide by $h$, rearrange, use $h_x/h=-G_x$. $\square$

**A.2 (pressure cancellation, §7.3).**  Wall term
$-\tfrac1h\sum_\pm\mathbf F\!\cdot\!\hat{\mathbf n}|_\pm
=-G_x(\tau_{33}-p)+O(\sigma^3,\sigma\varepsilon^2)$; geometric source
$\mathbf G\!\cdot\!\langle\mathbf F\rangle=\rho u(\mathbf V\!\cdot\!\mathbf
G)+pG_x-T_{11}G_x-T_{12}G_y$; the $\pm pG_x$ cancel **algebraically**,
leaving $\rho u(\mathbf V\!\cdot\!\mathbf G)+G_x(T_{33}-T_{11})-G_yT_{12}$.
$\square$ (Machine-checked.)

**A.3 (Theorem 5.6, $\Delta\equiv0$).**  Substitute
$\mathrm{Re}\,G_0$ from its definition into the solvability condition,
subtract the quasi-linear Q2D $x$-momentum (mass equation used to cancel
$\rho_0u_0(\mathbf V\!\cdot\!\mathbf G)$ on both sides); the residual

$$
\Delta=\partial^2_xu_0+\partial^2_{xy}v_0
-\partial_x\big[\nabla_{2\mathrm D}\!\cdot\!\mathbf V_0-(\mathbf
V_0\!\cdot\!\mathbf G)\big]-\partial_x(\mathbf V_0\!\cdot\!\mathbf G)
\equiv0 .
$$

$\square$ (SymPy twin F-1.)

**A.4 (degree preservation, §5.3-b3).**  With $\hat w_1=a(x,y,t)\,z^*$,
the only $z^*$-operators in the $\varepsilon^2$ system are $z^*\partial_{z^*}$
(degree-preserving) and multiplications by $z^*$-polynomials of degree
$\le2$ (forcing, $p_2$); the solution of a linear transport equation with
polynomial-in-$z^*$ data and forcing of degree $\le2$ remains polynomial of
degree $\le2$ (project on the monomial basis $1,z^*,z^{*2}$: the system
block-triangularises). $\square$ (SymPy twin F-6.)

## Appendix B — Machine verification inventory

All scripts under `tools/audit/`, runnable without the solver unless noted.

| ID | Script | Verifies |
|---|---|---|
| H-1 | `mms_area_variation.py` (LAYER 1, 12 checks) | master identity vs coded sources, all equations; Stokes; Boussinesq trace; $k$ single-count; Hirschfelder $\Sigma=0$ |
| H-2 | `mms_area_variation.py` (LAYER 2–3) | discrete order 2.00; soot block |
| F-1 | `foundation_checks.py` | Fredholm $\Delta\equiv0$ (Thm 5.6 / A.3) |
| F-2 | `foundation_checks.py` | code stress mapping $\equiv$ Thm 6.3; $\mathrm{tr}\,T=0$ |
| F-3 | `foundation_checks.py` | $\langle|S|^2\rangle=|S|^2_{2\mathrm D}+2(\mathbf V\!\cdot\!\mathbf G)^2$ |
| F-4 | `foundation_checks.py` | $k$ single-count (E4-bis): $\partial(\text{energy src})/\partial(\rho k)=\tfrac23\mathbf V\!\cdot\!\mathbf G$; momentum $k$-free |
| F-5 | `foundation_checks.py` | energy species term $\equiv\sum h_k(\mathbf G\!\cdot\!\mathbf j^{\rm eff}_k)$; $u_2$ zero-mean/BC |
| F-6 | `foundation_checks.py` | regime-(b) degree preservation (A.4) |
| F-7 | `foundation_checks.py` | slip truncation: discarded terms $O(\sigma^2)$ relative (§4.6) |
| F-8 | `foundation_checks.py` | chemistry constant $\propto\beta^2$ (Lemma 7.7, symbolic + Arrhenius numeric) |

## Appendix C — Concordance with the historical documents (errata absorbed)

The historical corpus (`teoria_finale_Q2D.md` v6.x, `teoria_Q2D_integrazione_v7.md`,
`guida_formulazione_rigorosa.md`, the seven JFM-style review/rebuttal rounds)
remains unmodified as audit trail.  This page absorbs:

| # | Historical defect | Resolution here |
|---|---|---|
| C-1 | Slope hypothesis $\sigma=O(\varepsilon)$ used everywhere, never declared; "da (H4): $h_x=O(\varepsilon)$" non sequitur; rebuttal/prompt statement "$|\mathbf G|L=O(\varepsilon)$" dimensionally inconsistent (true value $O(1)$) | (G1) slope clause + Lemma 2.1 + Lemma 5.4 + CE-1 |
| C-2 | 14-term enumeration of $\partial_j\tau_{1j}$ (sound) | inherited, §4.3 |
| C-3 | $z$-momentum: "$O(1/\mathrm{Re})$ subleading in $\varepsilon$" false at fixed $\mathrm{Re}$; forward-referenced ordering | Lemma 5.5 + reordered §5.2(a1–a4) |
| C-4 | English doc §6.1 Step 4: incomplete slip-BC constant for $u_2$ (slip-stress groups dropped); Step 2 parity$\Rightarrow$even-exponent non sequitur | full Neumann datum in §5.2(a5); Theorem 5.1(ii–iii) |
| C-5 | Runtime validity monitor removed partly on the C-1 erroneous bookkeeping | corrected monitor spec §10.1 |
| C-6 | "ESATTA" vs "(leading order)" labelling ambiguity for global conservation | exactness assigned to objects: averaged-flux identities EXACT; model-variable closures LEADING (§7.1, §7.5, §11) |
| C-7 | Chemistry closure "$O(\varepsilon^4)$" concealing the Zel'dovich amplification | Lemma 7.7, hypothesis (M2-iv) |
| C-8 | $T_{\rm loc}$ arithmetic ($K_\infty$ $10^3\to10^2$; $10^{-2}T_{\rm ref}\to$ "10 ms") and the consequent RDE-horizon claim | withdrawn and restated honestly, §10.5 |
| C-9 | Global $x$-momentum residual missing $-h_y\tau_{12}$ (errata E3, never applied inline) | §11 |
| C-10 | Energy-source sign errata (E1) and RANS diffusive sign errata (E2/E2-bis), $P^*$ double-count (E4-bis) | converged convention of §7 (machine-checked) |
| C-11 | Citation "Smith–Majda 1984" (unverifiable) supporting leading-order regularity | replaced by Kato (1975), Majda (1984, monograph), Schochet (1986) |
| C-12 | Regime (b) closure only sketched while RDE lives in (b) | full derivation §5.3 |
| C-13 | (H_b)/(G_b) labels used but never defined | (T5), §2.5 |

## Appendix D — Notation

$\varepsilon=h_{\max}/L$ (thinness); $\sigma=\lVert\nabla h\rVert_\infty$
(slope); $\hat h=h/(\varepsilon L)$; $\mathbf G=-\nabla h/h$;
$\mathbf V\!\cdot\!\mathbf G=uG_x+vG_y$; $L_h=\inf h/|\nabla h|$;
$\langle\cdot\rangle_z$ depth average; $[\cdot]$ wall jump; $f|_\pm$ wall
traces; $z^*=z/(\varepsilon L)$; $\zeta=2z/h$;
$\mu_{\rm eff}=\mu_l+\mu_t$; $\boldsymbol\tau_{\rm tot}$ Boussinesq-augmented
stress; $E_{\rm res}=e+\tfrac12(u^2+v^2)$; $\beta=E_a/(RT_0)$;
$\Pi$ the distinguished path (§2.4); regimes (a)/(b) per (M2-ii);
labels: **EXACT / LEADING / POSTULATE / CONJECTURE**.

## References

1. A. Matsumura, T. Nishida, *J. Math. Kyoto Univ.* **20**(1):67–104, 1980.
2. A. Valli, *Ann. Sc. Norm. Super. Pisa* (compressible NS, bounded domains), 1982.
3. T. Kato, "The Cauchy problem for quasi-linear symmetric hyperbolic
   systems", *Arch. Rational Mech. Anal.* **58**:181–205, 1975.
4. A. Majda, *Compressible Fluid Flow and Systems of Conservation Laws in
   Several Space Variables*, Springer Appl. Math. Sci. **53**, 1984.
5. S. Schochet, "The compressible Euler equations in a bounded domain:
   existence of solutions and the incompressible limit", *Comm. Math.
   Phys.* **104**:49–75, 1986.
6. J.-F. Gerbeau, B. Perthame, *DCDS-B* **1**(1):89–102, 2001.
7. D. Bresch, P. Noble, *Methods Appl. Anal.* **14**(2):87–118, 2007.
8. C. M. Dafermos, *Arch. Rational Mech. Anal.* **70**:167–179, 1979; and
   *Hyperbolic Conservation Laws in Continuum Physics*, 4th ed., Springer,
   2016, Ch. 5.
9. E. Feireisl, B. J. Jin, A. Novotný, *J. Math. Fluid Mech.*
   **14**:717–730, 2012; E. Feireisl, A. Novotný, *Arch. Rational Mech.
   Anal.* **204**:683–706, 2012.
10. J. D. Anderson, *Modern Compressible Flow*, 3rd ed., McGraw-Hill, 2003.
11. S. B. Pope, *Turbulent Flows*, CUP, 2000.
12. R. B. Bird, W. E. Stewart, E. N. Lightfoot, *Transport Phenomena*,
    2nd ed., Wiley, 2002.



