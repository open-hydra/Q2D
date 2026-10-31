#!/usr/bin/env bash
#===============================================================================
#  SRM-rocket-motor in mode 2: the same chamber with its 40 grain records in
#  mode 2 (an inflow, the inflow 404 at the static temperature of total Taf
#  with the mass flux of the burn law), run in inflow/ with the deck, the mesh
#  and the gas of the case: the gates of check.py hold as in mode 1.
#===============================================================================
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
rm -rf inflow && mkdir -p inflow/INPUT
ln -s ../MESH inflow/MESH
ln -s ../reference inflow/reference
for f in ic.tec phase.txt thermo.dat; do ln -s "../../INPUT/$f" "inflow/INPUT/$f"; done
cp input.ini check.py Q2D.sh inflow/
sed -i 's#^MASTERDIR=../../../..#MASTERDIR=../../../../..#' inflow/Q2D.sh
sed -i 's#parents\[3\]#parents[4]#' inflow/check.py
# the fifth value of each 502 line, the mode, 1 -> 2
awk 'prev ~ / 502 *$/ {$5 = "2.000000E+00"; $0 = "   " $0} {print; prev = $0}' INPUT/bc.txt > inflow/INPUT/bc.txt
test "$(grep -c '2.000000E+00' inflow/INPUT/bc.txt)" -eq 40 || { echo "[FAIL] the 40 grain lines of mode 2 not written"; exit 1; }
cd inflow && ./Q2D.sh solve > solve.log 2>&1 && grep -a -q 'SRM inflow (as 404)' solve.log && python3 -B check.py
