# V&V Coverage Matrix

This page is the maintained, honest map of what the Q2D test suite
actually covers: for every ctest case, which features it exercises and
against which class of reference it is gated; and for every solver
feature, which tests (if any) exercise it. Every cell is verifiable in
the repository (`test/CMakeLists.txt`, the per-case `check.py` /
`verify.py`, and the registry `docs/user/registry.md`).

**Reference-type hierarchy** (strongest first):

1. **Analytic exact** — closed-form solution evaluated in the check
   script itself (or from `test/common/pycheck/q2d_exact1d.py`).
2. **Published data** — experimental or independently computed results
   from the literature, committed as digitized data files.
3. **Code-to-code** — comparison against an independent solver
   (MOSE, SU2, Wind-US, CFL3D) or an independent model in the same run
   (RSM vs SST).
4. **Committed regression baseline** — a snapshot produced by Q2D
   itself with documented provenance; catches regressions, not model
   errors.
5. **Smoke** — existence/finiteness/positivity only.

## Per-test matrix

Tier `fast` runs in CI (`ctest -LE heavy`); `heavy` runs locally
(`ctest -L heavy`, 20 threads). NC = self-corrupting negative control.

| Test | Tier | Features exercised | Reference type | Gate quality |
|---|---|---|---|---|
| Sod79 | fast | 2D plane, SLAU, mc, MUSCL, RK3 time-accurate | analytic exact (Riemann) | probe 1% + whole-domain L1 + shock position |
| Ein91 | fast | 2D plane, HLLE, vanleer, MUSCL, RK3 time-accurate | analytic exact (Riemann) | probe 5% + L1 norms |
| Tor99 | fast | 2D plane, HLLE++, vanleer, MUSCL, RK3 time-accurate | analytic exact (Riemann) | probe 1% + L1 + shock position |
| Noh87 | fast | 2D plane, HLLE, first-order, RK3 time-accurate | analytic exact (Noh planar) | probe + L1 + both shock positions |
| ShuOsher-LIMO3 | fast | LIMO3 limiter, HLLC | committed limiter baselines (TV metric) | TV ordering + permanent anti-fallback tripwire |
| Sod79/Ein91/Tor99/Noh87-convergence | fast | grid-refinement ladder x1–x8 | analytic exact per cell | observed L1 order gate + monotonicity |
| noh-cylindrical | fast | **axisymmetric source terms with analytic-exact reference** (axis = BC 300, radius = j/y) | analytic exact (Noh cylindrical, m=1) | rel L1(rho) 8e-2 (measured 4.77e-2) + plateau + shock position + x-uniformity gates |
| chang-nozzle | fast | axisymmetric, HLLC, mc | published (Chang 1996, digitized) | Mach wall/axis + hard convergence sentinel |
| chang-nozzle-mg2 | fast | axisymmetric + multigrid 2-level | published + in-run restriction consistency | + coarse/fine gates |
| prandtl-meyer | fast | 2D plane, point probes | analytic exact (Prandtl-Meyer) | exit Mach 0.5% |
| prandtl-meyer-lineprobe | fast | line probes (registry-driven I/O) | self-consistency vs field + analytic level | bit-exact geometry + tolerance overlay |
| planar-nozzle | fast | **Q2D area-variation regime** | analytic (quasi-1D area law) | tiered Mach gates |
| planar-nozzle-mg2 | fast | area-variation + multigrid (Coarse_Cross_Section) | analytic + restriction consistency | 3-layer acceptance |
| borda-channel | fast | BC 420 Borda injector | analytic (Fievisohn 2016 eqs.) | 3 gates + branch-entered guard, NC |
| borda-channel-negctl | fast | negative-control harness | self (corrupted expectation) | must-fail verification |
| hypersonic-cylinder | fast | MUSCL-SD, HLLE, strong bow shock | analytic (normal shock) + Billig correlation | p0/T0 + carbuncle symmetry + standoff |
| ramp-channel | fast | MUSCL-SD shock detector, minmod | analytic (theta-beta-Mach) | hard ratio gates + SD localization, NC (upstream) |
| oblique-shock | fast | 2-block connectivity, MUSCL-SD, vanalbada | analytic (theta-beta-Mach) | ratio gates 1.5% + convergence |
| supersonic-forward-step | fast | time-accurate 2-block, transient | committed regression trace (provenance in reference/README.md) | 5% trace + shock position |
| rocket-nozzle-frozen | fast | axisymmetric, 13-species frozen chemistry, SLAU | committed reference field | Linf/L2 field gates + gas-law consistency |
| rocket-nozzle-equilibrium | fast | CEA equilibrium chemistry, RK2 | committed reference field | same |
| rocket-nozzle-finite-rate | fast | finite-rate chemistry, sdirk4b ODE, 19 species | committed reference field | same |
| SRM-rocket-motor | fast | BC 502 Vieille burn-rate, the products at total Taf | analytic (Vieille balance) | 15% (3.9%) + total temperature 5e-3 + convergence sentinel |
| SRM-rocket-motor-inflow | fast | BC 502 in mode 2 (an inflow) | analytic (Vieille balance) | as SRM-rocket-motor |
| SRM-grain-state | fast | BC 502: the gas let in over one step, both modes, on a low and a high face; the pressure of the grain face, the gas moving away from it; the refused lines | the face state and total enthalpy from the tables | 1e-9; nine refusals |
| rotating-channel | fast | rotating frame (rigid-body), BC 201 periodic | analytic (geostrophic balance) | dp/dy + drift gates |
| rotating-channel-rsm | fast | rotating frame + SSGLRR (exact Coriolis production G_ij, W_abs C5), time-accurate RK3 | homogeneous ODE (RK4 reference) | per-cell tracking 2 % + no-G_ij discriminance |
| multi-solver-103 | fast | BC 103 dispatcher | smoke (quiescent equilibrium) | explicit smoke-only scope |
| swbli | heavy | SA turbulence, shock/boundary-layer | code-to-code (MOSE) + published envelope (Schulein/SU2/Wind-US) | 0.1% Cf + separation envelope; under-converged = ctest Skipped |
| flat-plate-laminar | heavy | laminar Navier-Stokes | analytic (Blasius/Howarth table) | 0.15 profile gate + anti-vacuous guard |
| flat-plate-turbulent-coarse-kw | heavy | SST turbulence | published (NASA TMR, `cf-sst.dat` FUN3D zone) | primary physical-Cf 15% + shape 10% |
| ssglrr-chang-nozzle-turbulent | heavy | SSGLRR RSM x axisymmetric | code-to-code (RSM vs SST) + realizability | G1–G5 code verification |
| ssglrr-planar-nozzle-turbulent | fast | SSGLRR RSM x area-variation | code-to-code + realizability | G1–G5 code verification |

Unregistered auxiliary studies (not in ctest, run manually):
`flat-plate-turbulent/{coarse-SA,fine-SA}` (smoke checks; deferred for
the documented SA stiffness issue), `{coarse-rsm,fine-rsm}` (offline
NASA TMR SSG/LRR-w2012 comparison scripts with committed reference
data).

## Per-feature matrix

| Feature | Activation | Tests | Strongest reference | Gap |
|---|---|---|---|---|
| 2D plane regime | `[Q2D-Sources]` absent/`none` | 1D suite + most 2D | analytic exact | — |
| Axisymmetric regime | `active = axisymmetric` | noh-cylindrical, chang x2, rocket x3, SRM, ssglrr-chang | analytic exact (noh-cylindrical) | — |
| Q2D area-variation regime | `active = area_variation` | planar-nozzle x2, ssglrr-planar | analytic (area law) | — |
| Riemann: HLLC, HLLE, HLLE++, SLAU | `riemann-solver` | 27 tests (see per-test matrix) | analytic exact | — |
| Riemann: SLAU2, HLLC+, HLLEM, AUSM+, AUSM+-up, AUSM+-up2, exact | `riemann-solver` | **NONE** | — | untested registry values |
| Limiters: vanleer, mc, minmod, vanalbada, LIMO3 | `flux-limiter` | suite-wide + ShuOsher (LIMO3) | analytic/TV | — |
| Limiter: superbee | `flux-limiter` | **NONE** | — | untested |
| Space: MUSCL, MUSCL-SD, first-order | `space-reconstruction` | suite-wide | analytic | — |
| Time: RK2, RK3, time-accurate + steady | `time-scheme` | suite-wide | analytic | `euler` scheme untested |
| IRS (implicit residual smoothing) | `irs-beta > 0` auto-enables (`Assign_Setup.f90`), or explicit `irs = .true.` | most steady 2D cases (`irs-beta = 0.3-0.5`) | analytic (indirectly, via the gated cases) | no IRS-specific gate (exercised, not isolated) |
| Multigrid 2-level | `[Q2D-Multigrid] levels` | chang-mg2, planar-mg2 | in-run restriction consistency | 3+ levels untested |
| Turbulence: SA, SST, SSGLRR | `turbulence` | swbli, kw/rsm pairs, flat-plate | published (NASA TMR) / code-to-code | — |
| Turbulence: SA-R, SA-RC, SAcomp, Wilcox2006, SSGLRR-SD, -QCR2000, -blowcorr | model string (+suffix) | **NONE** | — | untested variants |
| Chemistry: frozen/equilibrium/finite-rate | `[Q2D-Physics] chemistry` | rocket-nozzle x3 | committed reference field | no analytic/published chemistry gate |
| ODE: sdirk4b | `ode-solver` | rocket-finite-rate | committed reference | H-radau5 (default) and ros4 never explicitly gated |
| Probes / line probes | `[Q2D-Probes]` / `[Q2D-LineProbes]` | prandtl-meyer (+lineprobe) | analytic + self-consistency | — |
| Rotating frame | `rotational-frame = rigid-body` | rotating-channel / rotating-channel-rsm | analytic (geostrophic) + homogeneous-ODE (SSGLRR G_ij, W_abs) | stationary-face path untested; SA-RC/SST-RC r-tilde frame term covered analytically only (the frame-invariance vortex analysis) |
| BC 300 (symmetry), 400 (extrapolation) | bc.txt | nearly all tests | mixed | — |
| BC 101 (connection) | bc.txt | oblique-shock, supersonic-forward-step, swbli | analytic | — |
| BC 103 (multi-solver), 201 (periodic), 420 (Borda), 502 (SRM Vieille) | bc.txt | multi-solver-103 / rotating-channel / borda-channel / SRM-rocket-motor, SRM-rocket-motor-inflow, SRM-grain-state | smoke / analytic | — |
| BC 401 (stagnation inlet) | bc.txt | chang x2, planar-nozzle family, rocket x3 | published/analytic | — |
| BC 405 (supersonic inlet) | bc.txt | hypersonic-cylinder, oblique-shock, prandtl-meyer x2, ramp-channel, forward-step, swbli | analytic | — |
| BC 301 (wall heat flux), 403 (mass-flux inlet), 407 (far-field) | bc.txt | flat-plate family only | published | heavy/unregistered tier only |
| BC 302 (wall temperature), 406 (back-pressure outflow) | bc.txt | swbli only | published envelope | heavy tier only |
| BC 102 (chimera), 200 (rotational periodic), 303 (melting wall), 304 (pyrolysis wall), 402 (time-varying p0), 404 (mass-flux + static T), 501 (manifold) | bc.txt | **NONE** | — | dispatched in `Mod_BC_Fluxes` but never exercised |
| Soot: LL91, LIN | `soot-generation` | **NONE** | — | untested |
| GSI (HDPE/PP/Carbon/HTPB) | BC-code driven (not registry-exposed) | **NONE** | — | untested + registry TODO |
| OpenMP do/task | `OMP-parallel` | all CI jobs (do default) | n/a | `task` strategy untested |
| MPI | compile-time `USE_MPI` | CI builds + runs the fast tier with the MPI binary (single rank) | n/a | multi-rank ctest missing |

## Known gaps (prioritized)

1. **Untested Riemann solvers** (SLAU2, HLLC+, HLLEM, AUSM+ family,
   exact): cheapest fix — parametrized Sod79 reruns per solver
   (seconds each).
2. **Untested limiter** superbee: same mechanism as ShuOsher/Sod.
3. **Turbulence variants without a gate** (SA-R, SA-RC, SAcomp,
   Wilcox2006, SSGLRR-SD, QCR2000, blowcorr): rotation/curvature
   corrections need a rotating/curved validation case; heavy-tier
   candidates.
4. **IRS is exercised but never isolated**: `irs-beta > 0`
   auto-enables smoothing (`Assign_Setup.f90`), so most steady cases
   run with IRS on — but no test compares IRS-on vs IRS-off
   convergence, so a smoothing regression only shows up indirectly.
5. **MPI multi-rank** never exercised by ctest (single-rank MPI binary
   only): needs a 2-block case under `mpirun -np 2`.
6. **Chemistry gates are regression-only**: no analytic/published
   anchor (e.g. an adiabatic-flame or shock-tube ignition case with an
   independent reference would upgrade the class).
7. **ODE solvers** H-radau5 (explicitly) and ros4: finite-rate rerun
   with each solver vs the same reference field.
8. **Soot models** (LL91/LIN) and **GSI walls** (301-304 fuels):
   no test at all; GSI additionally lacks registry exposure
   (`Register_Physics.f90` TODO).
9. **Multigrid 3+ levels** and **time-scheme `euler`**: single cheap
   smoke each.
10. **Never-exercised BC codes** 102 (chimera), 200 (rotational
    periodic), 303/304 (melting/pyrolysis walls), 402 (time-varying
    p0), 404 (mass-flux + static-T inlet), 501 (manifold): add minimal
    channel cases per BC family. Note also that the `BC_Fluxes`
    dispatcher has no `case default` — an unlisted code falls through
    silently (input validation in `IO_BC.f90` is the only guard).
