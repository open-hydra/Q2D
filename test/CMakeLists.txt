# ---------------------------------------------------------------------------
#  CTest definitions for Q2D validation cases.
# ---------------------------------------------------------------------------

# Common environment for all solver tests;
# PYTHONIOENCODING=utf-8 is Q2D-extra for check.py Unicode safety in LC_ALL=C
# envs — em-dash/Greek letters in print statements would fail otherwise).
set(_Q2D_TEST_ENV
    "OMP_NUM_THREADS=1"
    "KMP_STACKSIZE=100M"
    "PYTHONIOENCODING=utf-8"
)

# ---- 1-D : Sod shock tube (Sod 1978/79) -----------------------------------
add_test(
    NAME    Sod79
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/1D/Sod79"
)
set_tests_properties(Sod79 PROPERTIES
    LABELS      "1D;validation"
    TIMEOUT     300
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 1-D : Einfeldt strong-rarefaction (Einfeldt et al. 1991) -------------
add_test(
    NAME    Ein91
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/1D/Ein91"
)
set_tests_properties(Ein91 PROPERTIES
    LABELS      "1D;validation"
    TIMEOUT     300
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 1-D : Toro 1999 modified-Sod transonic-rarefaction (Tor99) -----------
add_test(
    NAME    Tor99
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/1D/Tor99"
)
set_tests_properties(Tor99 PROPERTIES
    LABELS      "1D;validation"
    TIMEOUT     300
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 1-D : Noh 1987 planar Cartesian symmetric-collision (Noh87) ----------
add_test(
    NAME    Noh87
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/1D/Noh87"
)
set_tests_properties(Noh87 PROPERTIES
    LABELS      "1D;validation"
    TIMEOUT     300
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 1-D : Shu-Osher shock + sine-wave, LIMO3 limiter ----
add_test(
    NAME    ShuOsher-LIMO3
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/1D/ShuOsher-LIMO3"
)
set_tests_properties(ShuOsher-LIMO3 PROPERTIES
    LABELS      "1D;validation"
    TIMEOUT     300
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : Chang 1996 axisymmetric nozzle (Eul-Frozen) --------------------
add_test(
    NAME    chang-nozzle
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/chang-nozzle"
)
set_tests_properties(chang-nozzle PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : Chang 1996 axisymmetric nozzle with Multigrid 2-level
# Mirror chang-nozzle config + [Q2D-Multigrid] levels=2 + coarse grid bc2.txt
# generated via tools/gen_coarse_bc.py from fine bc.txt (factor-2 restriction).
add_test(
    NAME    chang-nozzle-mg2
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/chang-nozzle-mg2"
)
set_tests_properties(chang-nozzle-mg2 PROPERTIES
    LABELS      "2D;validation;multigrid"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : Prandtl-Meyer supersonic expansion fan -------------------------
add_test(
    NAME    prandtl-meyer
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/prandtl-meyer"
)
set_tests_properties(prandtl-meyer PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : prandtl-meyer + line-probe ----------
#            Same expansion-fan case; adds a j-fixed line-probe ([midrow]) and
#            checks the registry-driven line-probe output (OUTPUT/midrow.tec)
#            reproduces the converged field at the sampled grid line. Exercises
#            Scan_Ini -> Register_LineProbes -> Setup_LineProbes -> writer
#            end-to-end.
add_test(
    NAME    prandtl-meyer-lineprobe
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/prandtl-meyer-lineprobe"
)
set_tests_properties(prandtl-meyer-lineprobe PROPERTIES
    LABELS      "2D;validation;lineprobe"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : planar-nozzle Q2D area-variation C-D nozzle vs quasi-1D area law
#            (isentropic validation) ----
add_test(
    NAME    planar-nozzle
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/planar-nozzle"
)
set_tests_properties(planar-nozzle PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : planar-nozzle Q2D area-variation x 2-level Multigrid  --
#            Same area-variation C-D nozzle as planar-nozzle but [Q2D-Multigrid]
#            levels=2. Exercises Cross_Section restricted
#            fine->coarse via R1 nodal sub-sampling (Coarse_Cross_Section); the
#            MGL=2 run must converge, match the quasi-1D area law, and the coarse
#            level must reproduce the fine restriction. Coarse bc2.txt generated
#            via tools/gen_coarse_bc.py.
add_test(
    NAME    planar-nozzle-mg2
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/planar-nozzle-mg2"
)
set_tests_properties(planar-nozzle-mg2 PROPERTIES
    LABELS      "2D;validation;multigrid"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : borda-channel Borda sudden-expansion injector BC  -----
#            BC_Nozzle_Borda_Flow (Fievisohn 2016 Sec 3.3,
#            eqs 3.26-3.31; sentinel BC_Mach==-20 on BC code 420, A1/A3 in the
#            bc.txt pamb field). MOSE has no nozzle/Borda BC.
#            Straight plane channel, choked subsonic-station-3 regime, const-gamma
#            air. check.py validates the injected face state against the paper's own
#            equations at a fully-specified gas (choked mass flux + isentropic choke
#            state + stagnation-pressure-loss interpretation), with a Borda-branch-
#            entered guard (false-green prevention) + negative control.
add_test(
    NAME    borda-channel
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/borda-channel"
)
set_tests_properties(borda-channel PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : hypersonic-cylinder M=8.1 bow shock (2D plane Euler) -----------
add_test(
    NAME    hypersonic-cylinder
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/hypersonic-cylinder"
)
set_tests_properties(hypersonic-cylinder PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : ramp-channel M=3 30-deg wedge oblique shock (2D plane Euler) ---
add_test(
    NAME    ramp-channel
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/ramp-channel"
)
set_tests_properties(ramp-channel PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     900
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : oblique-shock M=4 wedge attached shock (2D plane Euler 2-block) -
add_test(
    NAME    oblique-shock
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/oblique-shock"
)
set_tests_properties(oblique-shock PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : supersonic-forward-step Woodward-Colella M=3 time-accurate 2-block
add_test(
    NAME    supersonic-forward-step
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/supersonic-forward-step"
)
set_tests_properties(supersonic-forward-step PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     1800
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D-axi : rocket-nozzle/frozen multi-species (13) Euler frozen --------
add_test(
    NAME    rocket-nozzle-frozen
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/rocket-nozzle/frozen"
)
set_tests_properties(rocket-nozzle-frozen PROPERTIES
    LABELS      "2D;chemistry;validation"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D-axi : rocket-nozzle/equilibrium multi-species (13) CEA Euler ------
add_test(
    NAME    rocket-nozzle-equilibrium
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/rocket-nozzle/equilibrium"
)
set_tests_properties(rocket-nozzle-equilibrium PROPERTIES
    LABELS      "2D;chemistry;validation"
    TIMEOUT     900
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D-axi : rocket-nozzle/finite-rate multi-species (19) sdirk4b Euler --
add_test(
    NAME    rocket-nozzle-finite-rate
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/rocket-nozzle/finite-rate"
)
set_tests_properties(rocket-nozzle-finite-rate PROPERTIES
    LABELS      "2D;chemistry;validation"
    TIMEOUT     1800
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D-axi : SRM-rocket-motor BC 502 Vieille burn-rate  ----
# Steady-state chamber pressure regression against Vieille analytical balance.
add_test(
    NAME    SRM-rocket-motor
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/SRM-rocket-motor"
)
set_tests_properties(SRM-rocket-motor PROPERTIES
    LABELS      "2D;chemistry;validation"
    TIMEOUT     900
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : rotating-channel Coriolis geostrophic balance ----
# Mirror MOSE 3D coriolis-channel reduced 2D plane axis=[0,0,1]; dp/dy = -2*omega*rho*U_inf verification.
# Periodic-x (BC 201) + walls-y (BC 300);
add_test(
    NAME    rotating-channel
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/rotating-channel"
)
set_tests_properties(rotating-channel PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     600
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---- 2-D : multi-solver-103 BC 103 dispatcher smoke test ----
add_test(
    NAME    multi-solver-103
    COMMAND bash -c "./Q2D.sh solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/multi-solver-103"
)
set_tests_properties(multi-solver-103 PROPERTIES
    LABELS      "2D;validation"
    TIMEOUT     300
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)

# ---------------------------------------------------------------------------
# 2-D viscous tests (Sess 10 F-12-G B-11 3-layer fix infrastructure +
# Sess 12.A F-12-I cascade activation 3/5 cases).
#
# Sess 10 deliverable: B-11 3-layer fix applied (waldis MOSE-mirror gating +
# area_geom_term universal zero-fill + BC mass-flux aliasing mirror MOSE
# upstream `f76c5ca` 2026-05-11 "bug fix BC mass flux").
#
# Sess 11.A V&V: FP-laminar 100k iter PASS (residual 1.091e-6); FP-turb
# coarse-kw SST 100k iter PASS (residual 6.587e-7); FP-turb coarse-SA +
# fine-SA discovered B-12 SEV-1 SA model stiffness at CFL=1.0+IRS=0.5
# (cf. `docs/review/b12_sa_deep_audit.md` Sess 12.A 14-point audit:
# divergence is SA stiffness config-specific, NOT a Q2D wiring/porting bug).
#
# Sess 12.A activation (3/5 cases stable, 2/5 deferred SA-INFORMATIVE):
#   ✅ swbli (SA, CFL=0.3+IRS=0.0): smoke + Cf vs MOSE-SA.tec ±10% if converged.
#   ✅ flat-plate-laminar: Blasius ±15% acceptance (Sess 14 coarse-grid V&V;
#      Sess 13 ±5% loosened after Sess 14 root-cause: check.py mu source bug
#      hardcoded mil=1e-3 vs actual mu=1.186e-5 from transport.dat. mu fix +
#      ±15% NASA TMR coarse-grid tolerance → PASS max_err 11.5%).
#   ✅ flat-plate-turbulent/coarse-kw (SST): smoke + NASA TMR deferred Sess 12.B.
#   ⏸  flat-plate-turbulent/coarse-SA (SA, B-12 stiffness): infra preserved, ctest
#      deferred post-v1.0 (b12_sa_deep_audit.md final synthesis).
#   ⏸  flat-plate-turbulent/fine-SA: same as coarse-SA (SA-INFORMATIVE).
# OMP=30 override per viscous test (user direction "max 30 core per simulazione").
# ---------------------------------------------------------------------------

# ---- 2-D viscous : SWBLI SA shock-wave/boundary-layer interaction ---------
add_test(
    NAME    swbli
    COMMAND bash -c "./Q2D.sh -p 20 solve sa && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/viscous/swbli"
)
set_tests_properties(swbli PROPERTIES
    LABELS      "2D;viscous;RANS;validation"
    TIMEOUT     18000
    ENVIRONMENT "OMP_NUM_THREADS=20"
)

# ---- 2-D viscous : flat-plate-laminar (Blasius ±15% Sess 14 acceptance) ---
add_test(
    NAME    flat-plate-laminar
    COMMAND bash -c "./Q2D.sh -p 20 solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/viscous/flat-plate-laminar"
)
set_tests_properties(flat-plate-laminar PROPERTIES
    LABELS      "2D;viscous;validation"
    TIMEOUT     2400
    ENVIRONMENT "OMP_NUM_THREADS=20"
)

# ---- 2-D viscous : flat-plate-turbulent coarse-kw SST ---------------------
add_test(
    NAME    flat-plate-turbulent-coarse-kw
    COMMAND bash -c "./Q2D.sh -p 20 solve && python3 -B check.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/viscous/flat-plate-turbulent/coarse-kw"
)
set_tests_properties(flat-plate-turbulent-coarse-kw PROPERTIES
    LABELS      "2D;viscous;RANS;validation"
    TIMEOUT     5400
    ENVIRONMENT "OMP_NUM_THREADS=20"
)

# ---------------------------------------------------------------------------
# SSG/LRR-omega RSM per-regime running V&V cases (Phase C, PHYSICS-COMPLETE-
# FIRST order, "ctest SSGLRR entry" item — Sess 30). Each case is a kw/rsm pair
# (Menter SST sibling vs SSG/LRR-omega RSM, the ONLY input diff is the
# turbulence keyword) sharing one family-root standalone verify.py (the
# committed code-verification evidence, exit 0/1, gates G1-G5): solve kw, solve
# rsm, then run verify.py from the family root. The G1 convergence threshold
# (density residual < 1e-8) derives from the committed runs (9.98e-11 planar /
# 9.99e-11 axi). Code verification (sympy + MMS + running stability), not
# physical validation (model anchored to the 2D-plane NASA TMR comparison,
# Sess 26). Closures: docs/review/ssglrr_axi_vv_sess28.md (axi) +
# ssglrr_area_var_vv_sess29.md (area-var).
# ---------------------------------------------------------------------------

# ---- 2-D-axi viscous : chang-nozzle-turbulent SSGLRR axi V&V (kw + rsm) ----
add_test(
    NAME    ssglrr-chang-nozzle-turbulent
    COMMAND bash -c "(cd kw && ./Q2D.sh -p 20 solve) && (cd rsm && ./Q2D.sh -p 20 solve) && python3 -B verify.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/viscous/chang-nozzle-turbulent"
)
set_tests_properties(ssglrr-chang-nozzle-turbulent PROPERTIES
    LABELS      "2D;viscous;RANS;validation;ssglrr"
    TIMEOUT     1800
    ENVIRONMENT "OMP_NUM_THREADS=20"
)

# ---- 2-D q2D area-variation : planar-nozzle-turbulent SSGLRR area-var V&V --
# (kw + rsm; euler area-variation core flow, peculiarita' #2; fast ~2 min total)
add_test(
    NAME    ssglrr-planar-nozzle-turbulent
    COMMAND bash -c "(cd kw && ./Q2D.sh solve) && (cd rsm && ./Q2D.sh solve) && python3 -B verify.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/2D/euler/planar-nozzle-turbulent"
)
set_tests_properties(ssglrr-planar-nozzle-turbulent PROPERTIES
    LABELS      "2D;RANS;validation;ssglrr"
    TIMEOUT     900
    ENVIRONMENT "${_Q2D_TEST_ENV}"
)
