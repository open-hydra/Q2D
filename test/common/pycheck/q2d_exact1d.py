"""Analytic exact solutions for the 1D Euler verification ladders.

The committed ``reference/exact.txt`` tables (100 points, dx = 0.01) are fine
for the x4 regression gates in ``check.py``, but linear interpolation of a
table smears every discontinuity over one table spacing -- an O(1e-3) L1
error floor that dominates the x8 ladder level and destroys the observed
order of convergence. The ``verify.py --gate`` convergence gates therefore
evaluate the exact solution analytically, at every cell center of every
ladder level:

- :func:`riemann_exact`  -- exact ideal-gas Riemann solver (Toro 2009, ch. 4):
  Newton iteration for p*, then vectorized similarity sampling at (x - x0)/t.
  Covers Sod79, Ein91 (double rarefaction, near-vacuum star state) and Tor99.
- :func:`noh_exact` -- Noh (1987) uniform-implosion exact solution with a
  geometry parameter m (0 = planar symmetric collapse, 1 = cylindrical
  implosion onto the axis): infinite-strength shock at distance D t from
  the collapse centre with D = (g-1) u0 / 2, stagnation state
  rho = rho0 ((g+1)/(g-1))^(m+1), p = rho0 u0^2 (g+1)/2 ((g+1)/(g-1))^m,
  pre-shock geometric-convergence density rho0 (1 + u0 t / r)^m.
- :func:`noh_exact_planar` -- backward-compatible m=0 wrapper (planar Noh):
  bit-identical to the historical planar-only implementation.

numpy-only (python 3.6+); no matplotlib, no ExactPack.
"""
import numpy as np


def _fK(p, rhoK, pK, g):
    """Toro f_K(p) and its derivative for one side of the Riemann problem."""
    aK = np.sqrt(g * pK / rhoK)
    if p > pK:  # shock
        A = 2.0 / ((g + 1.0) * rhoK)
        B = (g - 1.0) / (g + 1.0) * pK
        sq = np.sqrt(A / (p + B))
        f = (p - pK) * sq
        df = sq * (1.0 - 0.5 * (p - pK) / (p + B))
    else:  # rarefaction
        f = 2.0 * aK / (g - 1.0) * ((p / pK) ** ((g - 1.0) / (2.0 * g)) - 1.0)
        df = (p / pK) ** (-(g + 1.0) / (2.0 * g)) / (rhoK * aK)
    return f, df


def _pstar(rhol, ul, pl, rhor, ur, pr, g, tol=1e-12, itmax=100):
    """Newton iteration for the star-region pressure (positivity-guarded)."""
    p = max(0.5 * (pl + pr), 1e-8)
    du = ur - ul
    for _ in range(itmax):
        fl, dfl = _fK(p, rhol, pl, g)
        fr, dfr = _fK(p, rhor, pr, g)
        dp = (fl + fr + du) / (dfl + dfr)
        p_new = p - dp
        if p_new <= 0.0:
            p_new = 0.5 * p  # bisect toward vacuum instead of overshooting
        if abs(p_new - p) < tol * max(p, 1e-30):
            return p_new
        p = p_new
    return p


def riemann_exact(x, t, x0, rhol, ul, pl, rhor, ur, pr, gamma=1.4):
    """Exact ideal-gas Riemann solution sampled at positions ``x``, time ``t``.

    Returns a dict with keys ``rho``, ``u``, ``p``, ``sie`` (numpy arrays)
    plus the star state ``pstar``/``ustar`` for diagnostics.
    """
    g = float(gamma)
    x = np.asarray(x, float)
    ps = _pstar(rhol, ul, pl, rhor, ur, pr, g)
    fl, _ = _fK(ps, rhol, pl, g)
    fr, _ = _fK(ps, rhor, pr, g)
    us = 0.5 * (ul + ur) + 0.5 * (fr - fl)

    al = np.sqrt(g * pl / rhol)
    ar = np.sqrt(g * pr / rhor)
    gp = (g + 1.0) / (2.0 * g)
    gm = (g - 1.0) / (2.0 * g)

    s = (x - x0) / t
    rho = np.empty_like(s)
    u = np.empty_like(s)
    p = np.empty_like(s)

    # ---- left of the contact -------------------------------------------
    left = s <= us
    if ps > pl:  # left shock
        SL = ul - al * np.sqrt(gp * ps / pl + gm)
        rsl = rhol * ((ps / pl + (g - 1.0) / (g + 1.0))
                      / ((g - 1.0) / (g + 1.0) * ps / pl + 1.0))
        pre = left & (s < SL)
        post = left & (s >= SL)
        rho[pre], u[pre], p[pre] = rhol, ul, pl
        rho[post], u[post], p[post] = rsl, us, ps
    else:  # left rarefaction
        asl = al * (ps / pl) ** ((g - 1.0) / (2.0 * g))
        head, tail = ul - al, us - asl
        pre = left & (s < head)
        fan = left & (s >= head) & (s < tail)
        star = left & (s >= tail)
        rho[pre], u[pre], p[pre] = rhol, ul, pl
        fac = 2.0 / (g + 1.0) + (g - 1.0) / ((g + 1.0) * al) * (ul - s[fan])
        rho[fan] = rhol * fac ** (2.0 / (g - 1.0))
        u[fan] = 2.0 / (g + 1.0) * (al + 0.5 * (g - 1.0) * ul + s[fan])
        p[fan] = pl * fac ** (2.0 * g / (g - 1.0))
        rho[star] = rhol * (ps / pl) ** (1.0 / g)
        u[star], p[star] = us, ps

    # ---- right of the contact ------------------------------------------
    right = ~left
    if ps > pr:  # right shock
        SR = ur + ar * np.sqrt(gp * ps / pr + gm)
        rsr = rhor * ((ps / pr + (g - 1.0) / (g + 1.0))
                      / ((g - 1.0) / (g + 1.0) * ps / pr + 1.0))
        post = right & (s <= SR)
        pre = right & (s > SR)
        rho[post], u[post], p[post] = rsr, us, ps
        rho[pre], u[pre], p[pre] = rhor, ur, pr
    else:  # right rarefaction
        asr = ar * (ps / pr) ** ((g - 1.0) / (2.0 * g))
        head, tail = ur + ar, us + asr
        star = right & (s <= tail)
        fan = right & (s > tail) & (s <= head)
        pre = right & (s > head)
        rho[star] = rhor * (ps / pr) ** (1.0 / g)
        u[star], p[star] = us, ps
        fac = 2.0 / (g + 1.0) - (g - 1.0) / ((g + 1.0) * ar) * (ur - s[fan])
        rho[fan] = rhor * fac ** (2.0 / (g - 1.0))
        u[fan] = 2.0 / (g + 1.0) * (-ar + 0.5 * (g - 1.0) * ur + s[fan])
        p[fan] = pr * fac ** (2.0 * g / (g - 1.0))
        rho[pre], u[pre], p[pre] = rhor, ur, pr

    sie = p / ((g - 1.0) * rho)
    return {"rho": rho, "u": u, "p": p, "sie": sie, "pstar": ps, "ustar": us}


def noh_exact(x, t, x0=0.5, rho0=1.0, u0=1.0, p0=0.0, gamma=1.4, m=0):
    """Noh (1987) uniform-implosion exact solution at time ``t``.

    Geometry parameter ``m``: 0 = planar symmetric collapse about ``x0``
    (two infinite-strength shocks), 1 = cylindrical implosion onto the
    axis (pass radii as ``x`` with ``x0 = 0``).  Standard strong-shock
    (p0 -> 0) Noh relations (Noh, JCP 72 (1987) 78, eq. 3.1):

    - shock speed          D  = (gamma - 1) u0 / 2,
    - post-shock state     rho = rho0 ((g+1)/(g-1))^(m+1),  u = 0,
                           p   = rho0 u0^2 (g+1)/2 ((g+1)/(g-1))^m,
    - pre-shock state      rho = rho0 (1 + u0 t / r)^m  (geometric
                           convergence factor; constant rho0 for m=0),
                           u = -u0 (toward the centre), p = p0.

    Cross-check (gamma = 5/3, rho0 = u0 = 1): D = 1/3 and, for m = 1,
    rho_post = 16, p_post = 16/3 (the classic cylindrical Noh values);
    for m = 0 the function is bit-identical to the historical
    ``noh_exact_planar`` (post-shock rho = 4, p = 4/3).
    """
    if m not in (0, 1):
        raise ValueError("noh_exact: geometry parameter m must be 0 or 1")
    g = float(gamma)
    x = np.asarray(x, float)
    D = 0.5 * (g - 1.0) * u0
    q = (g + 1.0) / (g - 1.0)          # strong-shock compression ratio
    # NOTE m=0 path: q**0 == 1.0 and X**0 == 1.0 exactly (also for inf/nan
    # operands at r=0), so every m=0 expression reduces bit-identically to
    # the historical planar implementation.
    rs = rho0 * (g + 1.0) / (g - 1.0) * q ** m
    ps = 0.5 * rho0 * u0 * u0 * (g + 1.0) * q ** m
    s = np.abs(x - x0)
    inside = s <= D * t
    with np.errstate(divide='ignore', invalid='ignore'):
        rho_pre = rho0 * (1.0 + u0 * t / s) ** m
    rho = np.where(inside, rs, rho_pre)
    u = np.where(inside, 0.0, -u0 * np.sign(x - x0))
    p = np.where(inside, ps, p0)
    sie = p / ((g - 1.0) * rho)
    return {"rho": rho, "u": u, "p": p, "sie": sie,
            "shock_left": x0 - D * t, "shock_right": x0 + D * t,
            "shock_speed": D}


def noh_exact_planar(x, t, x0=0.5, rho0=1.0, u0=1.0, p0=0.0, gamma=1.4):
    """Planar Noh symmetric-collapse exact solution at time ``t``.

    Backward-compatible wrapper of :func:`noh_exact` with m = 0 (two
    infinite-strength shocks moving out from ``x0`` at speed
    D = (gamma - 1) u0 / 2; stagnation state rho = rho0 (gamma+1)/(gamma-1),
    u = 0, p = rho0 u0^2 (gamma+1)/2).  Bit-identical to the historical
    planar-only implementation.
    """
    return noh_exact(x, t, x0=x0, rho0=rho0, u0=u0, p0=p0, gamma=gamma, m=0)


def observed_order(h, err):
    """Least-squares slope of log(err) vs log(h) (observed convergence order)."""
    lh = np.log(np.asarray(h, float))
    le = np.log(np.asarray(err, float))
    A = np.vstack([lh, np.ones_like(lh)]).T
    slope = np.linalg.lstsq(A, le, rcond=None)[0][0]
    return float(slope)
