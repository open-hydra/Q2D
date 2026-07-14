#!/bin/bash -
#===============================================================================
#
#          FILE: Q2D.sh
#         USAGE: run "./Q2D.sh [options] {solve|kill}" in the current shell
#   DESCRIPTION: Compile/run script for chang-nozzle-turbulent (axisymmetric
#                core-flow RANS, SSGLRR axi V&V family) (Q2D_new).
#
#===============================================================================

function print_usage {
  echo " solve  -->  run Q2D"
  echo " kill   -->  kill"
  echo " -p N   -->  parallel N threads"
  echo " -b     -->  background"
  exit 1
}

MASTERDIR=../../../../../
MASTER=$MASTERDIR/bin/Q2D
LOCAL=./bin/Q2D
BG=0; NTHREADS=1
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then shift; break; fi
  case $1 in
    -b|--background) BG=1; shift ;;
    -p|--parallel)   NTHREADS=$2; shift 2 ;;
    -h|--help)       print_usage; shift ;;
    *) break ;;
  esac
done
[[ $# == 0 ]] && print_usage

if [[ $1 == solve ]]; then
  cd INPUT
  ln -sf ../../../../../common/Air-viscous/phase.txt phase.txt
  ln -sf ../../../../../common/Air-viscous/thermo.dat thermo.dat
  ln -sf ../../../../../common/Air-viscous/transport.dat transport.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp -f $MASTER $LOCAL; fi
  if [[ $BG == 0 ]]; then $LOCAL
  else $LOCAL 2>errors_file >logfile & echo $! > .ID
  fi
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
