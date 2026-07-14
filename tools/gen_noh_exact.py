#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_noh_exact.py — Noh 1987 planar Cartesian symmetric-collision analytical
exact reference.

Pre-shock state: rho_pre = rho0, u_pre = u0 * sign(x0 - x), p_pre ~ 0.
Post-shock (stagnation) state derived from Rankine-Hugoniot strong-shock jump:
  rho_post = rho_pre * (gamma + 1) / (gamma - 1)
  u_post   = 0
  U_s      = (gamma - 1) / 2 * |u0|        (shock speed magnitude, lab frame)
  p_post   = (gamma + 1) / 2 * rho_pre * u0**2
           (equivalent: 2/(gamma+1) * rho_pre * (|u0| + U_s)**2,
            via shock-frame momentum balance with |u0 - U_s_signed| = (gamma+1)/2 * |u0|.)

For gamma=1.4, |u0|=1: rho_post = 6, p_post = 1.2, U_s = 0.2.
For gamma=5/3, |u0|=1: rho_post = 4, p_post = 4/3, U_s = 1/3.

Output format: 5 cols (x, rho, u, p, sie), 1 row per cell, no header.
sie = p / (rho * (gamma-1)).

Usage:
  python3 tools/gen_noh_exact.py \\
      --rho0 1.0 --u0 1.0 --p0 1e-15 \\
      --gamma 1.4  --x0 0.5  --t 1.0 \\
      --xmin 0.0   --xmax 1.0 --nx 100 \\
      --out reference/exact.txt
"""

import argparse
import sys
from pathlib import Path


def noh_state(x, x0, t, rho0, u0, p0, gamma):
    """Noh planar symmetric exact state at (x, t)."""
    U_s = 0.5 * (gamma - 1.0) * u0  # = 0.2 for gamma=1.4, u0=1
    rho_post = rho0 * (gamma + 1.0) / (gamma - 1.0)
    p_post   = 0.5 * (gamma + 1.0) * rho0 * u0**2
    if abs(x - x0) < U_s * t:
        return rho_post, 0.0, p_post
    elif x < x0:
        return rho0, u0,  p0
    else:
        return rho0, -u0, p0


def main(argv):
    ap = argparse.ArgumentParser(prog='gen_noh_exact.py',
        description='Noh 1987 planar Cartesian symmetric-collision exact.')
    ap.add_argument('--rho0', type=float, default=1.0)
    ap.add_argument('--u0',   type=float, default=1.0, help='speed magnitude (each half flows toward x0)')
    ap.add_argument('--p0',   type=float, default=1.0e-15, help='near-vacuum pressure')
    ap.add_argument('--gamma',type=float, default=1.4)
    ap.add_argument('--x0',   type=float, default=0.5, help='symmetry axis')
    ap.add_argument('--t',    type=float, required=True)
    ap.add_argument('--xmin', type=float, default=0.0)
    ap.add_argument('--xmax', type=float, default=1.0)
    ap.add_argument('--nx',   type=int,   required=True)
    ap.add_argument('--out',  type=Path,  required=True)
    args = ap.parse_args(argv)

    U_s = 0.5 * (args.gamma - 1.0) * args.u0
    rho_post = args.rho0 * (args.gamma + 1.0) / (args.gamma - 1.0)
    p_post   = 0.5 * (args.gamma + 1.0) * args.rho0 * args.u0**2
    print(f'[gen_noh_exact] Shock speed |U_s| = {U_s:.6f}, '
          f'post-shock rho = {rho_post:.6f}, p = {p_post:.6f}')
    print(f'[gen_noh_exact] Shock positions at t={args.t}: '
          f'left = {args.x0 - U_s*args.t:.6f}, right = {args.x0 + U_s*args.t:.6f}')

    dx = (args.xmax - args.xmin) / args.nx
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w') as f:
        for i in range(args.nx):
            x = args.xmin + (i + 0.5) * dx
            rho, u, p = noh_state(x, args.x0, args.t,
                                  args.rho0, args.u0, args.p0, args.gamma)
            sie = p / (rho * (args.gamma - 1.0))
            f.write(f'  {x:14.6f}  {rho:14.6f}  {u:14.6f}  {p:14.6f}  {sie:14.6f}\n')
    print(f'[gen_noh_exact] wrote {args.nx} cells -> {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
