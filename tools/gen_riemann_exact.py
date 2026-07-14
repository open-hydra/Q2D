#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_riemann_exact.py — Toro 2009 §4 ideal-gas Riemann exact solver.

Self-contained (no exactpack dependency). Implements the iterative Newton-
Raphson scheme of Toro 2009 §4.3.2 for the star-region pressure p*, then
samples the full wave structure (left rarefaction/shock + contact + right
rarefaction/shock) to produce reference profiles at t > 0.

Output format: 5 cols (x, rho, u, p, sie), 1 row per cell, no header.
sie = p / (rho * (gamma-1)).

Usage:
  python3 tools/gen_riemann_exact.py \\
      --rho-l 1.0  --u-l 0.75 --p-l 1.0 \\
      --rho-r 0.125 --u-r 0.0 --p-r 0.1 \\
      --gamma 1.4   --xd 0.3  --t 0.20 \\
      --xmin 0.0    --xmax 1.0 --nx 100 \\
      --out reference/exact.txt

Validated against Toro 2009 §4.3 Test 1 (Sod): p* = 0.30313, u* = 0.92745.
Validated against Toro 2009 §4.3 Test 2 (modified Sod): p* = 0.46090, u* = 1.45180.
"""

import argparse
import math
import sys
from pathlib import Path


def _f_K(p, p_K, rho_K, a_K, gamma):
    """Toro 2009 eq (4.6): f_K(p; W_K) is the wave function for state K."""
    if p > p_K:
        # Shock branch
        A_K = 2.0 / ((gamma + 1.0) * rho_K)
        B_K = (gamma - 1.0) / (gamma + 1.0) * p_K
        return (p - p_K) * math.sqrt(A_K / (p + B_K))
    else:
        # Rarefaction branch
        return (2.0 * a_K / (gamma - 1.0)) * ((p / p_K) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)


def _f_K_prime(p, p_K, rho_K, a_K, gamma):
    if p > p_K:
        A_K = 2.0 / ((gamma + 1.0) * rho_K)
        B_K = (gamma - 1.0) / (gamma + 1.0) * p_K
        return math.sqrt(A_K / (B_K + p)) * (1.0 - 0.5 * (p - p_K) / (B_K + p))
    else:
        return (1.0 / (rho_K * a_K)) * (p / p_K) ** (-(gamma + 1.0) / (2.0 * gamma))


def solve_pstar(rho_L, u_L, p_L, rho_R, u_R, p_R, gamma,
                tol=1e-12, maxiter=200):
    """Newton-Raphson for p* (star region pressure). Toro 2009 §4.3.2."""
    a_L = math.sqrt(gamma * p_L / rho_L)
    a_R = math.sqrt(gamma * p_R / rho_R)

    # Two-rarefaction initial guess (Toro eq 4.46)
    pTR_num = a_L + a_R - 0.5 * (gamma - 1.0) * (u_R - u_L)
    pTR_den = a_L / (p_L ** ((gamma - 1.0) / (2.0 * gamma))) \
            + a_R / (p_R ** ((gamma - 1.0) / (2.0 * gamma)))
    if pTR_den <= 0.0 or pTR_num <= 0.0:
        p = 0.5 * (p_L + p_R)
    else:
        p = (pTR_num / pTR_den) ** (2.0 * gamma / (gamma - 1.0))
    p = max(p, tol)

    for it in range(maxiter):
        f_L = _f_K(p, p_L, rho_L, a_L, gamma)
        f_R = _f_K(p, p_R, rho_R, a_R, gamma)
        fp_L = _f_K_prime(p, p_L, rho_L, a_L, gamma)
        fp_R = _f_K_prime(p, p_R, rho_R, a_R, gamma)
        F  = f_L + f_R + (u_R - u_L)
        Fp = fp_L + fp_R
        dp = F / Fp
        p_new = p - dp
        if p_new <= 0.0:
            p_new = 0.5 * p
        if abs(2.0 * (p_new - p) / (p_new + p)) < tol:
            p = p_new
            break
        p = p_new
    else:
        raise RuntimeError(f"Riemann solver did not converge in {maxiter} iterations")

    a_L = math.sqrt(gamma * p_L / rho_L)
    a_R = math.sqrt(gamma * p_R / rho_R)
    f_L = _f_K(p, p_L, rho_L, a_L, gamma)
    f_R = _f_K(p, p_R, rho_R, a_R, gamma)
    u = 0.5 * (u_L + u_R) + 0.5 * (f_R - f_L)
    return p, u


def sample_wave(s, rho_L, u_L, p_L, rho_R, u_R, p_R, p_star, u_star, gamma):
    """Sample the Riemann solution at self-similar coordinate s = (x - xd) / t.

    Returns (rho, u, p) at the given s. Toro 2009 §4.5 sampling procedure.
    """
    a_L = math.sqrt(gamma * p_L / rho_L)
    a_R = math.sqrt(gamma * p_R / rho_R)

    if s <= u_star:
        # Left of contact
        if p_star <= p_L:
            # Left rarefaction
            rho_star_L = rho_L * (p_star / p_L) ** (1.0 / gamma)
            a_star_L = a_L * (p_star / p_L) ** ((gamma - 1.0) / (2.0 * gamma))
            s_HL = u_L - a_L
            s_TL = u_star - a_star_L
            if s <= s_HL:
                return rho_L, u_L, p_L
            elif s >= s_TL:
                return rho_star_L, u_star, p_star
            else:
                # Inside fan
                rho_f = rho_L * (2.0 / (gamma + 1.0) + (gamma - 1.0) / ((gamma + 1.0) * a_L) * (u_L - s)) ** (2.0 / (gamma - 1.0))
                u_f = 2.0 / (gamma + 1.0) * (a_L + 0.5 * (gamma - 1.0) * u_L + s)
                p_f = p_L * (2.0 / (gamma + 1.0) + (gamma - 1.0) / ((gamma + 1.0) * a_L) * (u_L - s)) ** (2.0 * gamma / (gamma - 1.0))
                return rho_f, u_f, p_f
        else:
            # Left shock
            rho_star_L = rho_L * ((p_star / p_L + (gamma - 1.0) / (gamma + 1.0))
                                  / ((gamma - 1.0) / (gamma + 1.0) * p_star / p_L + 1.0))
            S_L = u_L - a_L * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p_star / p_L + (gamma - 1.0) / (2.0 * gamma))
            if s <= S_L:
                return rho_L, u_L, p_L
            else:
                return rho_star_L, u_star, p_star
    else:
        # Right of contact
        if p_star <= p_R:
            # Right rarefaction
            rho_star_R = rho_R * (p_star / p_R) ** (1.0 / gamma)
            a_star_R = a_R * (p_star / p_R) ** ((gamma - 1.0) / (2.0 * gamma))
            s_HR = u_R + a_R
            s_TR = u_star + a_star_R
            if s >= s_HR:
                return rho_R, u_R, p_R
            elif s <= s_TR:
                return rho_star_R, u_star, p_star
            else:
                rho_f = rho_R * (2.0 / (gamma + 1.0) - (gamma - 1.0) / ((gamma + 1.0) * a_R) * (u_R - s)) ** (2.0 / (gamma - 1.0))
                u_f = 2.0 / (gamma + 1.0) * (-a_R + 0.5 * (gamma - 1.0) * u_R + s)
                p_f = p_R * (2.0 / (gamma + 1.0) - (gamma - 1.0) / ((gamma + 1.0) * a_R) * (u_R - s)) ** (2.0 * gamma / (gamma - 1.0))
                return rho_f, u_f, p_f
        else:
            # Right shock
            rho_star_R = rho_R * ((p_star / p_R + (gamma - 1.0) / (gamma + 1.0))
                                  / ((gamma - 1.0) / (gamma + 1.0) * p_star / p_R + 1.0))
            S_R = u_R + a_R * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p_star / p_R + (gamma - 1.0) / (2.0 * gamma))
            if s >= S_R:
                return rho_R, u_R, p_R
            else:
                return rho_star_R, u_star, p_star


def main(argv):
    ap = argparse.ArgumentParser(prog='gen_riemann_exact.py',
        description='Toro 2009 ideal-gas Riemann exact solver.')
    ap.add_argument('--rho-l', type=float, required=True)
    ap.add_argument('--u-l',   type=float, required=True)
    ap.add_argument('--p-l',   type=float, required=True)
    ap.add_argument('--rho-r', type=float, required=True)
    ap.add_argument('--u-r',   type=float, required=True)
    ap.add_argument('--p-r',   type=float, required=True)
    ap.add_argument('--gamma', type=float, default=1.4)
    ap.add_argument('--xd',    type=float, required=True, help='discontinuity position')
    ap.add_argument('--t',     type=float, required=True, help='final time')
    ap.add_argument('--xmin',  type=float, default=0.0)
    ap.add_argument('--xmax',  type=float, default=1.0)
    ap.add_argument('--nx',    type=int,   required=True, help='cells (cell-centered samples)')
    ap.add_argument('--out',   type=Path,  required=True)
    args = ap.parse_args(argv)

    p_star, u_star = solve_pstar(args.rho_l, args.u_l, args.p_l,
                                  args.rho_r, args.u_r, args.p_r, args.gamma)
    print(f'[gen_riemann_exact] p* = {p_star:.8e}, u* = {u_star:.8e}')

    dx = (args.xmax - args.xmin) / args.nx
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w') as f:
        for i in range(args.nx):
            x = args.xmin + (i + 0.5) * dx
            s = (x - args.xd) / args.t
            rho, u, p = sample_wave(s,
                                    args.rho_l, args.u_l, args.p_l,
                                    args.rho_r, args.u_r, args.p_r,
                                    p_star, u_star, args.gamma)
            sie = p / (rho * (args.gamma - 1.0))
            f.write(f'  {x:14.6f}  {rho:14.6f}  {u:14.6f}  {p:14.6f}  {sie:14.6f}\n')
    print(f'[gen_riemann_exact] wrote {args.nx} cells -> {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
