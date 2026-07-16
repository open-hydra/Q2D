#!/bin/bash -
#===============================================================================
#
#          FILE: Q2D.sh
#         USAGE: run "./Q2D.sh [options] {solve|kill}" in the current shell
#   DESCRIPTION: Compile/run script for flat-plate-laminar (Blasius boundary-layer)
#                2D viscous test.
#                Mach 0.2 subsonic laminar Navier-Stokes, single block 65x65 mesh
#                (single K strip, Cartesian). Acceptance:
#                Blasius self-similar profile u/U_inf vs eta = y*sqrt(U/(nu*x))
#                along the x>0 wall.
#
#===============================================================================

function print_usage {
  echo " solve  -->  run Q2D"
  echo " kill   -->  kill"
  echo " -p N   -->  parallel N threads"
  echo " -b     -->  background"
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
  # fp-laminar is an artificial low-Re Blasius validation: it uses a
  # case-local transport.dat (mu=1e-3 -> Re/L~80,600) so the laminar BL is
  # resolvable on the coarse 65x65 mesh. The shared common/Air-viscous
  # transport.dat (real-air mu=1.186e-5) would push Re/L to 6.8e6 and leave
  # the BL unresolved. The MOSE peer shares this same latent limitation.
  ln -sf transport-laminar.dat transport.dat
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
