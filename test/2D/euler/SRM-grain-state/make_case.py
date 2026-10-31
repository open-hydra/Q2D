#!/usr/bin/env python3
"""SRM-grain-state: the state of the gas a grain (BC 502) lets in, over one step from rest.

The mesh of SRM-rocket-motor (MESH, symlinked), planar (no axisymmetric source), at rest: N2 at T0 = 300 K and
PB = 1e5 Pa, the other species, CH4, 1e-20. Its 40 grain records (face 1, i = 1) blow in CH4 (phase-ch4-n2.txt,
thermo-ch4-n2.dat: the two species of MOSE's test/1D/diff1 up to 3500 K) with mdot_ref = 1 kg/(m2 s) at p_ref = PB,
so that the mass flux is 1, and the total temperature TAF; the other records are those of SRM-rocket-motor. With
top, the grain is the 190 records of face 4 (j = 41, the contour of the nozzle, a high face: its normal points out of
the domain), which are walls (300) in SRM-rocket-motor, and face 1 a wall (300). With moving, the gas of the duct
moves at U0 = 10 m/s along +x, away from the grain of face 1.

  python3 -B make_case.py case [inflow] [top] [moving]  INPUT/: ic.tec, bc.txt (mode 1; inflow: mode 2), phase.txt,
                                                        thermo.dat
  python3 -B make_case.py bad KIND                      the case of mode 1 with the 502 lines that BAD[KIND] gives
"""
import shutil
import sys

PB, T0, TAF = 1.0e5, 300.0, 3000.0
U0 = 10.0                                   # the velocity of the gas of the moving case, along +x [m/s]
RU = 8314.51                                # the universal gas constant of the FLINT of Q2D
WM = {"CH4": 16.043, "N2": 28.014}
SRM = "../SRM-rocket-motor/INPUT"
# the 502 lines refused (mdot_ref, n, p_ref, T_af, mode, the mass fractions of CH4 and N2)
BAD = {
    "mode3":   "1.0, 0.4, 1.0E+05, 3.0E+03, 3, 1.0, 0.0",
    "mode15":  "1.0, 0.4, 1.0E+05, 3.0E+03, 1.5, 1.0, 0.0",
    "pref0":   "1.0, 0.4, 0.0, 3.0E+03, 1, 1.0, 0.0",
    "mdotneg": "-1.0, 0.4, 1.0E+05, 3.0E+03, 1, 1.0, 0.0",
    "inflow0": "0.0, 0.4, 1.0E+05, 3.0E+03, 2, 1.0, 0.0",
    "taf0":    "1.0, 0.4, 1.0E+05, 0.0, 1, 1.0, 0.0",
    "tafhigh": "1.0, 0.4, 1.0E+05, 3.5E+03, 1, 1.0, 0.0",   # the last row of thermo-ch4-n2.dat: above 3499 K
    "nan":     "NaN, 0.4, 1.0E+05, 3.0E+03, 1, 1.0, 0.0",
    "short":   "1.0, 0.4, 1.0E+05, 3.0E+03, 1",             # no mass fractions
}


def header(ln, code):
    """The record line ln (block, i, j, face, type) with the type code."""
    return "".join(f"{v:8d}" for v in [int(t) for t in ln.split()[:4]] + [code]) + "\n"


def write_case(line502, top=False, u0=0.0):
    with open(f"{SRM}/ic.tec") as f:
        lines = f.read().splitlines(True)
    head = lines[1]
    ni = int(head.split("I=")[1].split(",")[0])
    nj = int(head.split("J=")[1].split(",")[0])
    nodes = lines[2:2 + 2 * ni * nj]
    nc = (ni - 1) * (nj - 1)
    rho_n2 = PB / (RU / WM["N2"] * T0)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y"  "rho1" "rho2" "u" "v" "p"\n')
        f.write(head.replace("[3-6]=CELLCENTERED", "[3-7]=CELLCENTERED"))
        f.write("".join(nodes))
        for value in (1.0e-20, rho_n2, u0, 0.0, PB):
            f.write(f" {value:.15E}\n" * nc)
    with open(f"{SRM}/bc.txt") as f:
        bc = f.read().splitlines(True)
    out, i = [], 0
    while i < len(bc):
        t = bc[i].split()
        if t[-1] == "502":                        # the grain of SRM-rocket-motor, face 1: a wall with top
            out += [header(bc[i], 300)] if top else [bc[i], line502 + "\n"]
            i += 2
            continue
        if top and t[-1] == "300" and t[3] == "4":
            out += [header(bc[i], 502), line502 + "\n"]
        else:
            out.append(bc[i])
        i += 1
    with open("INPUT/bc.txt", "w") as f:
        f.writelines(out)
    shutil.copy("phase-ch4-n2.txt", "INPUT/phase.txt")
    shutil.copy("thermo-ch4-n2.dat", "INPUT/thermo.dat")


def main():
    if sys.argv[1] == "case":
        mode = 2 if "inflow" in sys.argv[2:] else 1
        write_case(f"1.0, 0.4, {PB:.6E}, {TAF:.1E}, {mode}, 1.0, 0.0", top="top" in sys.argv[2:],
                   u0=U0 if "moving" in sys.argv[2:] else 0.0)
    elif sys.argv[1] == "bad":
        write_case(BAD[sys.argv[2]])
    else:
        sys.exit(f"unknown request {sys.argv[1]}")


if __name__ == "__main__":
    main()
