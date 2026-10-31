#!/usr/bin/env bash
#===============================================================================
#  SRM-grain-state: the gas a grain (BC 502) lets in, over one step from rest,
#  in mode 1 (a wall that blows the products in) and in mode 2 (an inflow,
#  the inflow 404), on face 1 (a low face) and on face 4 (a high face, the
#  contour of the nozzle); the pressure of the grain face with the gas moving
#  away from it; and the property lines of a grain that Q2D refuses
#  (make_case.py, check.py).
#===============================================================================
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p INPUT
rm -f field-mode1.tec field-mode1-top.tec
fail=0

python3 -B make_case.py case || exit 1
if ./Q2D.sh solve > solve-mode1.log 2>&1 && [[ -f OUTPUT/field.tec ]]; then
  cp OUTPUT/field.tec field-mode1.tec
  echo "-- mode 1 (a wall that blows the products in)"
  python3 -B check.py OUTPUT/field.tec || fail=1
else
  echo "[FAIL] mode 1: the run did not complete (see solve-mode1.log)"; fail=1
fi

python3 -B make_case.py case inflow || exit 1
if ./Q2D.sh solve > solve-mode2.log 2>&1 && [[ -f OUTPUT/field.tec ]] && grep -a -q 'SRM inflow (as 404)' solve-mode2.log; then
  echo "-- mode 2 (an inflow, the inflow 404)"
  python3 -B check.py OUTPUT/field.tec field-mode1.tec || fail=1
else
  echo "[FAIL] mode 2: the run did not complete as an inflow (see solve-mode2.log)"; fail=1
fi

python3 -B make_case.py case top || exit 1
if ./Q2D.sh solve > solve-mode1-top.log 2>&1 && [[ -f OUTPUT/field.tec ]]; then
  cp OUTPUT/field.tec field-mode1-top.tec
  echo "-- mode 1 on face 4 (a high face: its normal points out of the domain)"
  python3 -B check.py --face 4 OUTPUT/field.tec || fail=1
else
  echo "[FAIL] mode 1 on face 4: the run did not complete (see solve-mode1-top.log)"; fail=1
fi

python3 -B make_case.py case inflow top || exit 1
if ./Q2D.sh solve > solve-mode2-top.log 2>&1 && [[ -f OUTPUT/field.tec ]] && grep -a -q 'SRM inflow (as 404)' solve-mode2-top.log; then
  echo "-- mode 2 on face 4"
  python3 -B check.py --face 4 OUTPUT/field.tec field-mode1-top.tec || fail=1
else
  echo "[FAIL] mode 2 on face 4: the run did not complete as an inflow (see solve-mode2-top.log)"; fail=1
fi

for mode in wall inflow; do   # the gas of the duct moving away from the grain of face 1: the pressure of its face
  python3 -B make_case.py case moving $([[ $mode == inflow ]] && echo inflow) || exit 1
  if ./Q2D.sh solve > "solve-moving-$mode.log" 2>&1 && [[ -f OUTPUT/field.tec ]]; then
    echo "-- the gas moving away from the grain, grain as $([[ $mode == wall ]] && echo 'a wall (mode 1)' || echo 'an inflow (mode 2)')"
    python3 -B check.py --moving OUTPUT/field.tec || fail=1
  else
    echo "[FAIL] the gas moving away from the grain, mode $mode: the run did not complete (see solve-moving-$mode.log)"; fail=1
  fi
done

echo "-- property lines of a grain that Q2D refuses"
bad() {   # bad <kind> <pattern>: the run stops at set-up with exit code 1, the message, no solution
  python3 -B make_case.py bad "$1" || exit 1
  rm -f OUTPUT/field.tec
  timeout 60 ./Q2D.sh solve > "solve-bad-$1.log" 2>&1
  local rc=$?
  if [[ $rc -eq 1 && ! -f OUTPUT/field.tec ]] && grep -a -q -E "$2" "solve-bad-$1.log"; then
    echo "  [ok  ] $1: stopped at set-up, exit code 1"
  else
    echo "  [FAIL] $1: exit code $rc, the line was not refused as it should (see solve-bad-$1.log)"; fail=1
  fi
}
h='bc.txt \(record 1, type 502\)'
bad mode3   "$h: mode = 3\.0000000000000000E\+00 is neither 1 \(a wall that blows the gas in\) nor 2 \(an inflow\)"
bad mode15  "$h: mode = 1\.5000000000000000E\+00 is neither 1"
bad pref0   "$h: p_ref = 0\.0000000000000000E\+00 is not positive"
bad mdotneg "$h: mdot_ref = -1\.0000000000000000E\+00 is negative"
bad inflow0 "$h: mdot_ref = 0 in mode 2: an inflow needs a mass flux"
bad taf0    "$h: T_af = 0\.0000000000000000E\+00 is not positive: no flame temperature \(the boundary-condition builder of ATLAS writes 0 for a grain without eq-CEA-file\)"
bad tafhigh "$h: T_af = 3\.5000000000000000E\+03 is above 3499 K, one kelvin below the last row of the tables of the gas"
bad nan     "$h: mdot_ref, n, p_ref, T_af and mode must be finite numbers"
bad short   "$h: the line is short or holds a value that is not a finite number"

if [[ $fail -ne 0 ]]; then
  echo "SRM-grain-state: FAILED"; exit 1
fi
echo "SRM-grain-state: all checks passed"
