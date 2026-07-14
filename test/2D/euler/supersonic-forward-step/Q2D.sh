#!/bin/bash -
# Q2D.sh — supersonic-forward-step (Woodward-Colella M=3) 2-block 2D plane
# Euler time-accurate.
function print_usage { echo "solve|kill, -p N parallel"; exit 1; }
MASTERDIR=../../../..
MASTER=$MASTERDIR/bin/Q2D
LOCAL=./bin/Q2D
BG=0; NTHREADS=1
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then shift; break; fi
  case $1 in -b|--background) BG=1; shift;; -p|--parallel) NTHREADS=$2; shift 2;; -h|--help) print_usage; shift;; *) break;; esac
done
[[ $# == 0 ]] && print_usage
if [[ $1 == solve ]]; then
  cd INPUT
  ln -sf ../../../../common/nondim-Air/phase.txt phase.txt
  ln -sf ../../../../common/nondim-Air/thermo.dat thermo.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp $MASTER $LOCAL; fi
  if [[ $BG == 0 ]]; then $LOCAL
  else $LOCAL 2>errors_file >logfile & echo $! > .ID
  fi
fi
if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
