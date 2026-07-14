#!/bin/bash -
#===============================================================================
#
#          FILE: Q2D.sh
#
#         USAGE: run "./Q2D.sh [options] {solve|kill}" in the current shell
#
#   DESCRIPTION: Compile/run script for ramp-channel 2D steady-state test
#                (Q2D_new).
#                2D plane Euler M=3 over a 30-deg wedge inside a channel
#                (peculiarita' #1: K=1 strip; NO axisymmetric source).
#===============================================================================

function print_usage {
  echo
  echo "Tasks"
  echo " solve                 -->     run Q2D on the default mesh"
  echo " test                  -->     run Q2D sweep across Riemann solvers (HLLE, HLLC, HLLC+, SLAU)"
  echo " kill                  -->     kill the process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

MASTERDIR=../../../..
MASTER=$MASTERDIR/bin/Q2D
LOCAL=./bin/Q2D

BG=0
NTHREADS=1

while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then shift; break; fi
  case $1 in
    -b | --background) BG=1; shift ;;
    -p | --parallel)   NTHREADS=$2; shift 2 ;;
    -h | --help)       print_usage; shift ;;
    * )                break ;;
  esac
done
[[ $# == 0 ]] && print_usage

DIR=$(pwd)

if [[ $1 == solve ]]; then
  cd INPUT
  ln -sf ../../../../common/Air/phase.txt phase.txt
  ln -sf ../../../../common/Air/thermo.dat thermo.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi
  if [[ $BG == 0 ]]; then
    $LOCAL
  else
    $LOCAL 2>errors_file >logfile &
    echo $! > .ID
  fi
fi

if [[ $1 == test ]]; then
  cd INPUT
  ln -sf ../../../../common/Air/phase.txt phase.txt
  ln -sf ../../../../common/Air/thermo.dat thermo.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi

  riemann_schemes=(HLLE HLLC HLLC+ SLAU)
  input_file="input.ini"
  original_riemann=$(grep -E '^riemann-solver\s*=' "$input_file" | head -n 1 || true)

  for scheme in "${riemann_schemes[@]}"; do
    echo "Running with riemann = $scheme"
    perl -i -pe "s/^riemann-solver\s*=\s*.*/riemann-solver = $scheme/" "$input_file"
    if [[ $BG == 0 ]]; then
      $LOCAL
    else
      $LOCAL 2>errors_file >logfile &
      echo $! > .ID
    fi
    mv OUTPUT/field.tec OUTPUT/field-$scheme.tec
    echo "Completed run for $scheme"
    echo "----------------------------------"
  done

  if [[ -n "$original_riemann" ]]; then
    perl -i -pe "s/^riemann-solver\s*=\s*.*/$original_riemann/" "$input_file"
  fi
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
