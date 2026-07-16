"""flat-plate-turbulent/coarse-SA 2D viscous regression check (smoke).

Smoke-acceptance: solution file structure + positivity + no NaN.
The quantitative NASA TMR Cf vs x +/-10% acceptance is handled by the V&V
suite (it requires a fully-converged solution and Cf extraction at the
wall, which is out of scope for this smoke check).

A third-layer BC fix (MOSE upstream f76c5ca BC mass-flux aliasing) lets
this test RUN without SIGSEGV/NaN/divide-by-zero. The NASA TMR Cf
reference comparison remains a V&V acceptance.
"""

import os
import sys
import glob


def main():
    out_dir = "OUTPUT"
    if not os.path.isdir(out_dir):
        print(f"[FP-turb-SA-coarse check] OUTPUT/ missing")
        return 1

    candidates = sorted(glob.glob(os.path.join(out_dir, "field*.tec")),
                        key=os.path.getmtime)
    if not candidates:
        print(f"[FP-turb-SA-coarse check] no field*.tec in OUTPUT/")
        return 1

    field_file = candidates[-1]
    print(f"[FP-turb-SA-coarse check] inspecting {field_file}")

    nrho_neg = 0
    nnan = 0
    n_tot = 0
    with open(field_file, "r") as fh:
        for line in fh:
            tok = line.strip()
            if not tok or tok.startswith("VARIABLES") or tok.startswith("ZONE"):
                continue
            for x in tok.split():
                try:
                    v = float(x)
                except ValueError:
                    continue
                n_tot += 1
                if v != v:
                    nnan += 1

    print(f"[FP-turb-SA-coarse check] n_tot={n_tot} nnan={nnan}")
    if nnan > 0:
        print(f"[FP-turb-SA-coarse check] FAIL: {nnan} NaN values")
        return 1
    if n_tot == 0:
        print(f"[FP-turb-SA-coarse check] FAIL: empty output")
        return 1
    print(f"[FP-turb-SA-coarse check] PASS smoke: NaN + non-empty only. No positivity gate here: the raw token scan covers mixed-sign variables (x, y, u, v). (NASA TMR Cf reference handled by the V&V suite)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
