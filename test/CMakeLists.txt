# ---------------------------------------------------------------------------
#  CTest definitions for Q2D validation cases.
#
#  Tiering (drives CI):
#    fast   -- minutes-scale cases; run on every push/PR (GitHub CI runs
#              `ctest -LE heavy`).
#    heavy  -- RANS/viscous cases needing many threads and/or hours of wall
#              clock; excluded from hosted CI, run locally / on beefy nodes:
#              `ctest -L heavy` (see docs/development/testing.md).
#
#  Every test uses SKIP_RETURN_CODE 77: a check script that cannot run its
#  quantitative gate exits 77 and shows up as "Skipped" -- never as a silent
#  "Passed" (false green).
# ---------------------------------------------------------------------------

set(Q2D_SKIP_RC 77)

# Common environment for all solver tests;
# PYTHONIOENCODING=utf-8 is Q2D-extra for check.py Unicode safety in LC_ALL=C
# envs -- em-dash/Greek letters in print statements would fail otherwise).
set(_Q2D_TEST_ENV
    "OMP_NUM_THREADS=1"
    "KMP_STACKSIZE=100M"
    "PYTHONIOENCODING=utf-8"
)

# Helper: register one validation case.
#   q2d_add_case(<name>
#                DIRECTORY <dir relative to test/>
#                [COMMAND <full shell command>]      # default: solve + check.py
#                [TIMEOUT <sec>]                     # default 600
#                [LABELS <l1;l2;...>]
#                [ENVIRONMENT <e1;e2;...>]           # default ${_Q2D_TEST_ENV}
#                [PROCESSORS <n>])                    # default 1
function(q2d_add_case NAME)
  cmake_parse_arguments(ARG "" "DIRECTORY;COMMAND;TIMEOUT;PROCESSORS" "LABELS;ENVIRONMENT" ${ARGN})
  if(NOT ARG_COMMAND)
    set(ARG_COMMAND "./Q2D.sh solve && python3 -B check.py")
  endif()
  if(NOT ARG_TIMEOUT)
    set(ARG_TIMEOUT 600)
  endif()
  if(NOT ARG_ENVIRONMENT)
    set(ARG_ENVIRONMENT "${_Q2D_TEST_ENV}")
  endif()
  if(NOT ARG_PROCESSORS)
    set(ARG_PROCESSORS 1)
  endif()
  add_test(
    NAME    ${NAME}
    COMMAND bash -c "${ARG_COMMAND}"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/${ARG_DIRECTORY}"
  )
  set_tests_properties(${NAME} PROPERTIES
    LABELS           "${ARG_LABELS}"
    TIMEOUT          ${ARG_TIMEOUT}
    ENVIRONMENT      "${ARG_ENVIRONMENT}"
    PROCESSORS       ${ARG_PROCESSORS}
    SKIP_RETURN_CODE ${Q2D_SKIP_RC}
    # Tests sharing a case directory (e.g. borda-channel + its negative
    # control) write the same OUTPUT/; never run them concurrently.
    RESOURCE_LOCK    "${ARG_DIRECTORY}"
  )
endfunction()

# ---------------------------------------------------------------------------
# 1-D Riemann / canonical cases (analytic references)
# ---------------------------------------------------------------------------

# Sod 1978/79 shock tube
q2d_add_case(Sod79
  DIRECTORY 1D/Sod79
  TIMEOUT 300
  LABELS "1D;validation;fast")

# Einfeldt et al. 1991 strong double rarefaction (1-2-0-3)
q2d_add_case(Ein91
  DIRECTORY 1D/Ein91
  TIMEOUT 300
  LABELS "1D;validation;fast")

# Toro 1999 modified-Sod transonic rarefaction
q2d_add_case(Tor99
  DIRECTORY 1D/Tor99
  TIMEOUT 300
  LABELS "1D;validation;fast")

# Noh 1987 planar Cartesian symmetric collision
q2d_add_case(Noh87
  DIRECTORY 1D/Noh87
  TIMEOUT 300
  LABELS "1D;validation;fast")

# Shu-Osher shock + sine-wave interaction, LIMO3 limiter (TV-based gates with
# committed first-order/minmod/vanleer/MC baselines + negative control)
q2d_add_case(ShuOsher-LIMO3
  DIRECTORY 1D/ShuOsher-LIMO3
  TIMEOUT 300
  LABELS "1D;validation;fast")

# 1-D observed-order-of-convergence gates: x1/x2/x4/x8 refinement ladder
# (Q2D.sh test) + verify.py --gate (hard exit 0/1 on the observed L1 order).
q2d_add_case(Sod79-convergence
  DIRECTORY 1D/Sod79
  COMMAND "./Q2D.sh test && python3 -B verify.py --gate"
  TIMEOUT 900
  LABELS "1D;convergence;fast")

q2d_add_case(Ein91-convergence
  DIRECTORY 1D/Ein91
  COMMAND "./Q2D.sh test && python3 -B verify.py --gate"
  TIMEOUT 900
  LABELS "1D;convergence;fast")

q2d_add_case(Tor99-convergence
  DIRECTORY 1D/Tor99
  COMMAND "./Q2D.sh test && python3 -B verify.py --gate"
  TIMEOUT 900
  LABELS "1D;convergence;fast")

q2d_add_case(Noh87-convergence
  DIRECTORY 1D/Noh87
  COMMAND "./Q2D.sh test && python3 -B verify.py --gate"
  TIMEOUT 900
  LABELS "1D;convergence;fast")

# ---------------------------------------------------------------------------
# 2-D Euler cases
# ---------------------------------------------------------------------------

# Chang 1996 axisymmetric nozzle (Eul-Frozen); wall+axis Mach vs digitized
# published data.
q2d_add_case(chang-nozzle
  DIRECTORY 2D/euler/chang-nozzle
  LABELS "2D;validation;fast")

# Chang 1996 nozzle with 2-level multigrid: mirrors chang-nozzle config +
# [Q2D-Multigrid] levels=2 + coarse-grid bc2.txt generated via
# tools/gen_coarse_bc.py (factor-2 restriction). check.py additionally gates
# the coarse-level field against the fine restriction.
q2d_add_case(chang-nozzle-mg2
  DIRECTORY 2D/euler/chang-nozzle-mg2
  LABELS "2D;validation;multigrid;fast")

# Prandtl-Meyer supersonic expansion fan vs analytic isentropic result
q2d_add_case(prandtl-meyer
  DIRECTORY 2D/euler/prandtl-meyer
  LABELS "2D;validation;fast")

# Same expansion fan; adds a j-fixed line probe ([midrow]) and checks the
# registry-driven line-probe output (OUTPUT/midrow.tec) reproduces the
# converged field at the sampled grid line. Exercises Scan_Ini ->
# Register_LineProbes -> Setup_LineProbes -> writer end-to-end.
q2d_add_case(prandtl-meyer-lineprobe
  DIRECTORY 2D/euler/prandtl-meyer-lineprobe
  LABELS "2D;validation;lineprobe;fast")

# Planar C-D nozzle, Q2D area-variation regime vs quasi-1D isentropic area law
q2d_add_case(planar-nozzle
  DIRECTORY 2D/euler/planar-nozzle
  LABELS "2D;validation;fast")

# Same area-variation nozzle x 2-level multigrid ([Q2D-Multigrid] levels=2).
# Exercises Cross_Section restricted fine->coarse via R1 nodal sub-sampling
# (Coarse_Cross_Section); MGL=2 run must converge, match the quasi-1D area
# law, and the coarse level must reproduce the fine restriction.
q2d_add_case(planar-nozzle-mg2
  DIRECTORY 2D/euler/planar-nozzle-mg2
  LABELS "2D;validation;multigrid;fast")

# Borda sudden-expansion injector BC (BC_Nozzle_Borda_Flow, Fievisohn 2016
# Sec 3.3 eqs 3.26-3.31; sentinel BC_Mach==-20 on BC code 420, A1/A3 in the
# bc.txt pamb field). Straight plane channel, choked subsonic-station-3
# regime, const-gamma air. check.py validates the injected face state against
# the paper's equations (choked mass flux + isentropic choke state +
# stagnation-pressure-loss interpretation) with a Borda-branch-entered guard.
q2d_add_case(borda-channel
  DIRECTORY 2D/euler/borda-channel
  LABELS "2D;validation;fast")

# Negative control for borda-channel: with BORDA_NEG_CONTROL=1 the check
# corrupts its own expectation and must detect the mismatch (exits 0 iff the
# corrupted gate FAILS). Institutionalizes the false-green tripwire in CI.
q2d_add_case(borda-channel-negctl
  DIRECTORY 2D/euler/borda-channel
  COMMAND "./Q2D.sh solve && BORDA_NEG_CONTROL=1 python3 -B check.py"
  LABELS "2D;negative-control;fast")

# M=8.1 hypersonic cylinder bow shock (2D plane Euler)
q2d_add_case(hypersonic-cylinder
  DIRECTORY 2D/euler/hypersonic-cylinder
  LABELS "2D;validation;fast")

# Cylindrical Noh (1987) implosion, 2D-axi: the only axisymmetric case with
# an ANALYTIC exact solution -- verifies the axisymmetric geometric source
# terms (Lib_Source_Axisymmetry) against D=(g-1)/2, rho_post=((g+1)/(g-1))^2.
# TIMEOUT: measured 12 s single-thread (gfortran TESTING, 2026-07-15) incl.
# check.py; 300 s = the suite's minimum fast-tier timeout (>x20 headroom).
q2d_add_case(noh-cylindrical
  DIRECTORY 2D/euler/noh-cylindrical
  TIMEOUT 300
  LABELS "2D;validation;axisymmetric;fast")

# M=3 30-deg wedge oblique shock in a channel (2D plane Euler)
q2d_add_case(ramp-channel
  DIRECTORY 2D/euler/ramp-channel
  TIMEOUT 900
  LABELS "2D;validation;fast")

# M=4 wedge attached oblique shock, 2-block (theta-beta-Mach theory gates)
q2d_add_case(oblique-shock
  DIRECTORY 2D/euler/oblique-shock
  LABELS "2D;validation;fast")

# Woodward-Colella M=3 forward step, time-accurate, 2-block
q2d_add_case(supersonic-forward-step
  DIRECTORY 2D/euler/supersonic-forward-step
  TIMEOUT 1800
  LABELS "2D;validation;fast")

# Rocket nozzle, 2D-axi, multi-species (13) frozen chemistry
q2d_add_case(rocket-nozzle-frozen
  DIRECTORY 2D/euler/rocket-nozzle/frozen
  LABELS "2D;chemistry;validation;fast")

# Rocket nozzle, 2D-axi, multi-species (13) CEA equilibrium chemistry.
# TIMEOUT: measured 448 s single-thread (gfortran-12 TESTING, 2026-07-15) x4
# headroom for slower CI runners.
q2d_add_case(rocket-nozzle-equilibrium
  DIRECTORY 2D/euler/rocket-nozzle/equilibrium
  TIMEOUT 1800
  LABELS "2D;chemistry;validation;fast")

# Rocket nozzle, 2D-axi, multi-species (19) finite-rate chemistry (sdirk4b).
# TIMEOUT: measured 1300 s single-thread (gfortran-12 TESTING, 2026-07-15)
# x2.8 headroom for slower CI runners -- the longest fast-tier case.
q2d_add_case(rocket-nozzle-finite-rate
  DIRECTORY 2D/euler/rocket-nozzle/finite-rate
  TIMEOUT 3600
  LABELS "2D;chemistry;validation;fast")

# SRM rocket motor, BC 502 Vieille burn-rate: steady chamber pressure vs
# analytical Vieille balance.
q2d_add_case(SRM-rocket-motor
  DIRECTORY 2D/euler/SRM-rocket-motor
  TIMEOUT 900
  LABELS "2D;chemistry;validation;fast")

# Coriolis geostrophic balance (mirror of MOSE 3D coriolis-channel reduced to
# the 2D plane, axis=[0,0,1]); dp/dy = -2*omega*rho*U_inf verification.
# Periodic-x (BC 201) + walls-y (BC 300).
q2d_add_case(rotating-channel
  DIRECTORY 2D/euler/rotating-channel
  LABELS "2D;validation;fast")

# BC 103 multi-solver coupling dispatcher smoke test
q2d_add_case(multi-solver-103
  DIRECTORY 2D/euler/multi-solver-103
  TIMEOUT 300
  LABELS "2D;validation;fast")

# ---------------------------------------------------------------------------
# 2-D viscous / RANS cases (heavy tier: many threads, long wall clock; not run
# on hosted CI -- see docs/development/testing.md for the local workflow).
#
# History: B-11 3-layer fix (waldis MOSE-mirror gating + area_geom_term
# universal zero-fill + BC mass-flux aliasing mirror of MOSE f76c5ca);
# flat-plate-turbulent coarse-SA / fine-SA deferred post-v1.0 (B-12 SA
# stiffness at CFL=1.0+IRS=0.5 -- audit concluded config-specific stiffness,
# not a porting bug).
# ---------------------------------------------------------------------------

# SWBLI: Schulein Mach-5 shock-wave/boundary-layer interaction (SA model);
# wall tauX vs MOSE-SA reference + separation/reattachment locations vs
# Schulein experiment / SU2 / Wind-US envelope.
q2d_add_case(swbli
  DIRECTORY 2D/viscous/swbli
  COMMAND "./Q2D.sh -p 20 solve sa && python3 -B check.py"
  TIMEOUT 18000
  ENVIRONMENT "OMP_NUM_THREADS=20;PYTHONIOENCODING=utf-8"
  PROCESSORS 20
  LABELS "2D;viscous;RANS;validation;heavy")

# Laminar flat plate vs Blasius (Howarth 1938 table)
q2d_add_case(flat-plate-laminar
  DIRECTORY 2D/viscous/flat-plate-laminar
  COMMAND "./Q2D.sh -p 20 solve && python3 -B check.py"
  TIMEOUT 2400
  ENVIRONMENT "OMP_NUM_THREADS=20;PYTHONIOENCODING=utf-8"
  PROCESSORS 20
  LABELS "2D;viscous;validation;heavy")

# Turbulent flat plate, Menter SST (coarse-kw) vs NASA TMR CFL3D Cf
q2d_add_case(flat-plate-turbulent-coarse-kw
  DIRECTORY 2D/viscous/flat-plate-turbulent/coarse-kw
  COMMAND "./Q2D.sh -p 20 solve && python3 -B check.py"
  TIMEOUT 5400
  ENVIRONMENT "OMP_NUM_THREADS=20;PYTHONIOENCODING=utf-8"
  PROCESSORS 20
  LABELS "2D;viscous;RANS;validation;heavy")

# ---------------------------------------------------------------------------
# SSG/LRR-omega RSM per-regime code-verification pairs (kw sibling vs rsm; the
# only input diff is the turbulence keyword) sharing one family-root
# verify.py (gates G1-G5: convergence < 1e-8, realizability, RSM-vs-SST axis
# pressure, production-term regime reduction, anisotropy sentinel).
# ---------------------------------------------------------------------------

# 2-D-axi: chang-nozzle-turbulent SSGLRR V&V (kw + rsm)
q2d_add_case(ssglrr-chang-nozzle-turbulent
  DIRECTORY 2D/viscous/chang-nozzle-turbulent
  COMMAND "(cd kw && ./Q2D.sh -p 20 solve) && (cd rsm && ./Q2D.sh -p 20 solve) && python3 -B verify.py"
  TIMEOUT 1800
  ENVIRONMENT "OMP_NUM_THREADS=20;PYTHONIOENCODING=utf-8"
  PROCESSORS 20
  LABELS "2D;viscous;RANS;validation;ssglrr;heavy")

# Q2D area-variation: planar-nozzle-turbulent SSGLRR V&V (kw + rsm; Euler
# area-variation core flow; fast ~2 min total)
q2d_add_case(ssglrr-planar-nozzle-turbulent
  DIRECTORY 2D/euler/planar-nozzle-turbulent
  COMMAND "(cd kw && ./Q2D.sh solve) && (cd rsm && ./Q2D.sh solve) && python3 -B verify.py"
  TIMEOUT 900
  LABELS "2D;RANS;validation;ssglrr;fast")

# Frame-invariance gate for the SSGLRR rotating-frame closures (exact Coriolis
# production G_ij + W_abs C5): uniform relative flow + anisotropic R in a frame
# at Omega=10; every cell must track the homogeneous ODE (verify.py, RK4
# reference, frame-invariance proven by construction). The gate
# fails by ~3 orders of magnitude on a binary without G_ij (discriminant G3).
q2d_add_case(rotating-channel-rsm
  DIRECTORY 2D/euler/rotating-channel-rsm
  COMMAND "./Q2D.sh solve && python3 -B verify.py"
  TIMEOUT 900
  LABELS "2D;RANS;rotating-frame;validation;fast")
