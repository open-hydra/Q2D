#!/bin/bash -
#===============================================================================
#
#          FILE: Q2D.sh
#
#         USAGE: run "./Q2D.sh [options] {solve|kill}" in the current shell
#
#   DESCRIPTION: Compile/run script for the borda-channel 2D steady-state test.
#                Exercises the Borda sudden-expansion injector inflow BC
#                (BC_Nozzle_Borda_Flow; Fievisohn 2016 Sec 3.3,
#                eqs 3.26-3.31; sentinel BC_Mach == -20 on BC code 420, A1/A3 in the
#                bc.txt pamb field). MOSE has no nozzle/Borda BC.
#                Straight plane channel, choked subsonic-station-3 regime; constant-
#                gamma air (common/Air table, gamma=1.4).
#===============================================================================

function print_usage {
  echo
  echo "Tasks"
  echo " solve                 -->     run Q2D on the default mesh"
  echo " kill                  -->     kill the process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

# Directories and files definition
MASTERDIR=../../../..
MASTER=$MASTERDIR/bin/Q2D
LOCAL=./bin/Q2D

# Default Options
BG=0
NTHREADS=1

# Parse command-line options
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    shift
    break
  fi
  case $1 in
    -b | --background) BG=1; shift ;;
    -p | --parallel)   NTHREADS=$2; shift 2 ;;
    -h | --help)       print_usage; shift ;;
    * )                break ;;
  esac
done
[[ $# == 0 ]] && print_usage

if [[ $1 == solve ]]; then
  cd INPUT
  ln -sf ../../../../common/Air/phase.txt phase.txt
  ln -sf ../../../../common/Air/thermo.dat thermo.dat
  cd ..
  mkdir -p OUTPUT bin
  rm -f OUTPUT/*
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi
  if [[ $BG == 0 ]]; then
    $LOCAL || exit 1
  else
    $LOCAL 2>errors_file >logfile &
    echo $! > .ID
  fi
fi

if [[ $1 == kill ]]; then
  kill -9 "$(cat .ID)"
fi
