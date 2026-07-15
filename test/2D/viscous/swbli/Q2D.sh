#!/bin/bash -
#===============================================================================
#
#          FILE: Q2D.sh
#         USAGE: run "./Q2D.sh [options] {solve|kill}" in the current shell
#   DESCRIPTION: Compile/run script for swbli (Shock/Boundary Layer Interaction)
#                2D viscous RANS-SA test (Q2D_new). 3-block grid (850x70 + 850x54 +
#                850x43 cells, K=1 strip Cartesian, peculiarità #1 2D plane).
#                Acceptance: M=5 supersonic inflow + isothermal wall Tw=300K
#                + shock-induced boundary-layer separation (RANS-SA model).
#===============================================================================

function print_usage {
  echo " solve {sa|sst}  -->  run Q2D with chosen RANS variant"
  echo " kill            -->  kill"
  echo " -p N            -->  parallel N threads"
  echo " -b              -->  background"
  exit 1
}

MASTERDIR=../../../../
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
  ln -sf ../../../../common/Air-viscous/phase.txt phase.txt
  ln -sf ../../../../common/Air-viscous/thermo.dat thermo.dat
  ln -sf ../../../../common/Air-viscous/transport.dat transport.dat
  if [[ $2 == sa ]]; then
    ln -sf ic-sa.tec ic.tec
    ln -sf bc-sa.txt bc.txt
  elif [[ $2 == sst ]]; then
    ln -sf ic-sst.tec ic.tec
    ln -sf bc-sst.txt bc.txt
  else
    echo "You need to specify RANS model: sa, sst"
    exit
  fi
  cd ..
  mkdir -p OUTPUT bin
  rm -f OUTPUT/*
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp -f $MASTER $LOCAL; fi
  if [[ $BG == 0 ]]; then $LOCAL || exit 1
  else $LOCAL 2>errors_file >logfile & echo $! > .ID
  fi
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
