#!/usr/bin/env python3
"""SRM-grain-state: the gas a grain lets in, from the changes of the cells next to it over one step from rest.

    check.py [--face 4 | --moving] <field.tec> [<field of mode 1>]

From rest the gas of the duct gives no flux through the faces between cells but the pressure, which the pressure of
the grain face, that of the cell, balances (the complete flux of the grain): over the step dt the cell next to the
grain gains dt (A/V) times the flux of the face, mass g = 1 kg/(m2 s), momentum g u along the inward normal of the
face, energy g h0. The ratios of the changes are those of the face, whatever dt and the cell: u = d(rho u_n) / d(rho),
h0 = d(rho E) / d(rho), u_n the velocity along the inward normal, the face of the nodes turned clockwise. The grain is
on face 1, the 40 cells i = 1 (a low face, its normal +x points into the domain), or with --face 4 on face 4, the 190
cells j = 40 along the contour of the nozzle (a high face, its normal points out of the domain; the inward normals from
45 degrees to one side of -y to 15 to the other).
  S0  every cell gains mass, d(rho) > 1e-9 rho (2e-4 rho measured; round-off is 1e-16 rho): the grain lets gas in
      (applied with the sign of the other side, it would take it out, or let none in)
  S1  the mass the cell gains is all CH4 (d(rho_CH4) / d(rho) = 1 to 1e-9), N2 unchanged
  S2  the gas constant of the face state, R = u p / (g T_w), with T_w the static temperature of total T_af of
      h(T_w) + u^2/2 = h(T_af) (u = g R T_w / p, the table of CH4 as Q2D interpolates it, bisection here), is that of
      CH4 (1e-9), not that of N2: the composition the grain imposes
  S3  h0 = h(T_af), the total enthalpy of the products (1e-9; with T_af as the static temperature it would be
      u^2/2 above it, 1.2e-5)
  S4  with <field of mode 1>: the changes of mode 2 (the inflow 404) are those of mode 1 (1e-10)
With --moving (face 1) the gas of the duct moves at U0 = 10 m/s along +x, away from the grain: the 38 cells next to it
away from the corners also lose the uniform gas through their other faces, whose fluxes are exact (a uniform state),
rho0 U0 per unit of dt A/V in all (the grain face is vertical), so that the N2 a cell loses gives dt A/V, and its
changes the flux of the grain face:
  S1m the mass flux of the grain, mdot_ref (p/p_ref)^n = 1 at the pressure of the cell (1e-9)
  S3m its total enthalpy, h(T_af) (1e-9)
  S5  the pressure on the grain face, that of the cell (1e-9): the complete flux of the grain. A reflective flux
      alongside the grain without its pressure, the former BC 502, puts there the pressure of the impermeable wall,
      which the gas moving away from it lowers: 2.4e-3 below that of the cell with SLAU
The tolerances follow the 15 digits of field.tec: the changes are about 4e-4 of the values, so that a ratio of changes
holds about 11 digits (1e-11 measured).
"""
import math
import sys

PB, T0, TAF, G = 1.0e5, 300.0, 3000.0, 1.0
U0 = 10.0                                   # the velocity of the gas of the moving case, along +x [m/s]
RU = 8314.51                                # the universal gas constant of the FLINT of Q2D (Runiv, Lib_ThermoTransport)
WM = (16.043, 28.014)                       # CH4, N2
R = tuple(RU / w for w in WM)


def table(path="thermo-ch4-n2.dat"):
    rows, zone = {}, None
    with open(path) as f:
        for ln in f:
            if ln.startswith("ZONE"):
                zone = ln.split('"')[1]
                rows[zone] = {}
                continue
            t = ln.split()
            if zone and len(t) == 5:
                rows[zone][int(float(t[0]))] = float(t[2])
    return rows


def h(rows, sp, T):
    i = int(T)
    return rows[sp][i] + (rows[sp][i + 1] - rows[sp][i]) * (T - i)


def read_field(path):
    with open(path) as f:
        lines = f.read().splitlines()
    names = lines[0].split("=", 1)[1].replace('"', " ").split()
    head = lines[1]
    ni = int(head.split("I=")[1].split(",")[0])
    nj = int(head.split("J=")[1].split(",")[0])
    nums = []
    for ln in lines[2:]:
        try:
            nums.extend(float(v) for v in ln.split())
        except ValueError:
            continue
    nn, nc, out, pos = ni * nj, (ni - 1) * (nj - 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 2 else nc
        out[name] = nums[pos:pos + cnt]
        pos += cnt
    return out, ni - 1, nj - 1


def grain_changes(path, rows, face):
    """Per cell next to the grain: d(rho), and per unit of it d(rho_CH4), d(rho_N2), d(rho u_n), d(rho E)."""
    f, nci, ncj = read_field(path)
    x, y, ni = f["x"], f["y"], nci + 1
    rho0 = (1.0e-20, PB / (R[1] * T0))
    e0 = rho0[1] * (h(rows, "N2", T0) - R[1] * T0) + rho0[0] * (h(rows, "CH4", T0) - R[0] * T0)
    if face == 1:   # the cell (1, j), its face between the nodes (1, j) and (1, j + 1)
        cells = [(0, j, (0, j), (0, j + 1)) for j in range(ncj)]
    else:           # the cell (i, ncj), its face between the nodes (i, ncj + 1) and (i + 1, ncj + 1)
        cells = [(i, ncj - 1, (i, ncj), (i + 1, ncj)) for i in range(nci)]
    out = []
    for i, j, a, b in cells:
        c = j * nci + i
        na, nb = a[1] * ni + a[0], b[1] * ni + b[0]
        tx, ty = x[nb] - x[na], y[nb] - y[na]
        nx, ny = ty / math.hypot(tx, ty), -tx / math.hypot(tx, ty)
        r1, r2 = f["rho(1)"][c], f["rho(2)"][c]
        u, v, p = f["u"][c], f["v"][c], f["p"][c]
        rho = r1 + r2
        T = p / (r1 * R[0] + r2 * R[1])
        e = r1 * (h(rows, "CH4", T) - R[0] * T) + r2 * (h(rows, "N2", T) - R[1] * T) + 0.5 * rho * (u * u + v * v)
        drho = rho - sum(rho0)
        if drho == 0.0:
            return None                               # no gas let in: the grain face was closed over the step
        out.append((drho, (r1 - rho0[0]) / drho, (r2 - rho0[1]) / drho, rho * (u * nx + v * ny) / drho,
                    (e - e0) / drho))
    return out


def moving_flux(path, rows, t_w):
    """The moving case: per cell next to the grain of face 1, the mass flux, total enthalpy and pressure of the grain
    face, from the changes of the cell and the exact flux of the uniform gas through its east face."""
    f, nci, ncj = read_field(path)
    rho0 = (1.0e-20, PB / (R[1] * T0))
    hn2 = h(rows, "N2", T0)
    e0 = rho0[1] * (hn2 - R[1] * T0) + rho0[0] * (h(rows, "CH4", T0) - R[0] * T0) + 0.5 * sum(rho0) * U0 * U0
    u_b = G * R[0] * t_w / PB                     # the velocity of the products on the face
    out = []
    # the cells (1, j) away from the corners, whose walls (faces 3 and 4, slopes 2e-7 and 2e-6) let none of the gas
    # through where the uniform flow would cross them (2e-7 and 2.5e-6 on the mass flux of those two cells)
    for j in range(1, ncj - 1):
        c = j * nci
        r1, r2 = f["rho(1)"][c], f["rho(2)"][c]
        u, v, p = f["u"][c], f["v"][c], f["p"][c]
        rho = r1 + r2
        T = p / (r1 * R[0] + r2 * R[1])
        e = r1 * (h(rows, "CH4", T) - R[0] * T) + r2 * (h(rows, "N2", T) - R[1] * T) + 0.5 * rho * (u * u + v * v)
        k = (rho0[1] - r2) / (rho0[1] * U0)       # dt A/V, from the N2 the cell loses through its east face
        if not k > 0.0:
            return None
        g = (r1 - rho0[0]) / k
        p_face = (rho * u - sum(rho0) * U0) / k - G * u_b + sum(rho0) * U0 * U0 + PB
        h_in = ((e - e0) / k + rho0[1] * U0 * (hn2 + 0.5 * U0 * U0)) / g
        out.append((g, p_face, h_in))
    return out


def main(argv):
    face, moving = 1, False
    while argv and argv[0].startswith("--"):
        if argv[0] == "--face":
            face, argv = int(argv[1]), argv[2:]
        elif argv[0] == "--moving":
            moving, argv = True, argv[1:]
        else:
            sys.exit(f"unknown option {argv[0]}")
    path, path1 = argv[0], (argv[1] if len(argv) > 1 else None)
    rows = table()
    h0 = h(rows, "CH4", TAF)
    lo, hi = TAF - 50.0, TAF
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if h(rows, "CH4", mid) + 0.5 * (G * R[0] * mid / PB) ** 2 < h0:
            lo = mid
        else:
            hi = mid
    t_w = 0.5 * (lo + hi)
    if moving:
        fl = moving_flux(path, rows, t_w)
        if fl is None:
            print(f"S1m-S5 FAIL: a cell next to the grain did not lose the gas of the duct ({path})")
            return 1
        s1 = max(abs(g / G - 1.0) for g, _, _ in fl)
        s3 = max(abs(hh / h0 - 1.0) for _, _, hh in fl)
        s5 = max(abs(p / PB - 1.0) for _, p, _ in fl)
        print(f"S1m the mass flux of the grain, the gas of the duct at {U0:g} m/s away from it ({len(fl)} cells): "
              f"{fl[0][0]:.12f} kg/(m2 s), largest difference {s1:.1e} {'PASS' if s1 < 1e-9 else 'FAIL'}")
        print(f"S3m its total enthalpy {fl[0][2]:.9e} J/kg, h(T_af) {h0:.9e}: largest difference {s3:.1e} "
              f"{'PASS' if s3 < 1e-9 else 'FAIL'}")
        print(f"S5 the pressure on the grain face {fl[0][1]:.6f} Pa, that of the cell {PB:g} Pa: largest difference "
              f"{s5:.1e} {'PASS' if s5 < 1e-9 else 'FAIL'}")
        return int(not (s1 < 1e-9 and s3 < 1e-9 and s5 < 1e-9))
    ch = grain_changes(path, rows, face)
    if ch is None:
        print(f"S0-S3 FAIL: no gas entered the cells next to the grain over the step ({path}: d(rho) = 0)")
        return 1
    s0 = min(d for d, _, _, _, _ in ch) / (PB / (R[1] * T0))     # d(rho) / rho
    s1 = max(max(abs(a - 1.0), abs(b)) for _, a, b, _, _ in ch)
    r_w = [u * PB / (G * t_w) for _, _, _, u, _ in ch]
    s2 = max(abs(r / R[0] - 1.0) for r in r_w)
    s3 = max(abs(e / h0 - 1.0) for _, _, _, _, e in ch)
    u = G * R[0] * t_w / PB
    static = 0.5 * u * u / h0
    rc = 0
    print(f"S0 the {len(ch)} cells next to the grain on face {face} gain mass: smallest d(rho)/rho {s0:.3e} "
          f"{'PASS' if s0 > 1e-9 else 'FAIL'}")
    print(f"S1 the mass gained is CH4: {s1:.1e} {'PASS' if s1 < 1e-9 else 'FAIL'}")
    print(f"S2 the face state at T_w = {t_w:.9f} K (total {TAF:g} K, u = {u:.4f} m/s): R = {r_w[0]:.9f} J/(kg K), "
          f"CH4 {R[0]:.9f}, N2 {R[1]:.9f}: largest difference {s2:.1e} {'PASS' if s2 < 1e-9 else 'FAIL'}")
    print(f"S3 the energy of the gas let in, h0 = {ch[0][4]:.9e} J/kg, h(T_af) = {h0:.9e}: largest difference {s3:.1e} "
          f"(T_af static: {static:.1e}) {'PASS' if s3 < 1e-9 else 'FAIL'}")
    rc |= not (s0 > 1e-9 and s1 < 1e-9 and s2 < 1e-9 and s3 < 1e-9)
    if path1:
        c1 = grain_changes(path1, rows, face)
        if c1 is None:
            print(f"S4 FAIL: no gas entered the cells next to the grain in {path1} (d(rho) = 0)")
            return 1
        s4 = max(abs(a - b) / max(abs(a), abs(b), 1e-300) for x, y in zip(ch, c1) for a, b in zip(x[3:], y[3:]))
        print(f"S4 the changes of mode 2 against those of mode 1: largest difference {s4:.1e} "
              f"{'PASS' if s4 < 1e-10 else 'FAIL'}")
        rc |= not s4 < 1e-10
    return int(rc)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
