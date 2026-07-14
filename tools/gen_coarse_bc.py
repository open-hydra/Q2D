#!/usr/bin/env python3
"""
Generate coarse-grid bc.txt for Q2D multigrid level N from fine bc.txt.

Q2D bc.txt format (cf. src/lib/io/IO_BC.f90:114-140 dispatch):
  Head record: 5 or 6 integers (block, i, j, [k,] face, BC_code).
  Property line(s) AFTER head depend on BC code:
    * BC 101, 103, 201, 301-309, 401-407, 420, 501-502: 1 property line
    * BC 102 (chimera): variable-length block (ci+cii lines)
    * BC 200, 300, 400: NO property line

Multigrid factor-rap restriction:
  i_coarse = (i_fine + rap - 1) // rap   (1-based: fine i={1,2}->coarse 1 with rap=2)
  j_coarse = (j_fine + rap - 1) // rap
  k_coarse = (k_fine + rap - 1) // rap
  Block, face, BC code preserved; property values preserved verbatim.

Deduplicate by (block, i, j, k, face, code) tuple after coarse-index remap.

Usage:
  python3 tools/gen_coarse_bc.py <fine_bc.txt> <coarse_bc.txt> [--rap 2]

Example:
  python3 tools/gen_coarse_bc.py \\
    test/2D/euler/chang-nozzle/INPUT/bc.txt \\
    test/2D/euler/chang-nozzle-mg2/INPUT/bc2.txt \\
    --rap 2
"""
import argparse
import sys

# BC codes with 1 property line after head (cf. IO_BC.f90:131)
SINGLE_PROP_CODES = (101, 103, 201,
                     301, 302, 303, 304, 305, 306, 307, 308, 309,
                     401, 402, 403, 404, 405, 406, 407,
                     420, 501, 502)
# BC codes with NO property line (cf. IO_BC.f90:138)
NO_PROP_CODES = (200, 300, 400)
# BC 102 chimera: variable property block (not handled here; rare in test cases)


def parse_bc_file(path):
    """Parse Q2D bc.txt into list of (head_tuple, property_lines) records.
    head_tuple = (block, i, j, k, face, code) — k defaulted to 1 if 5-field format.
    property_lines = list (0 or 1 entries) — verbatim strings."""
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    records = []
    n = len(lines)
    idx = 0
    while idx < n:
        head_line = lines[idx]
        parts = head_line.split()
        if len(parts) == 5:
            block, i, j, face, code = (int(p) for p in parts)
            k = 1
        elif len(parts) == 6:
            block, i, j, k, face, code = (int(p) for p in parts)
        else:
            sys.exit(f"Unexpected head format at line {idx+1}: '{head_line}'")
        head = (block, i, j, k, face, code)
        idx += 1
        props = []
        if code in SINGLE_PROP_CODES:
            if idx >= n:
                sys.exit(f"Missing property line for BC {code} at line {idx+1}")
            props.append(lines[idx])
            idx += 1
        elif code in NO_PROP_CODES:
            pass
        elif code == 102:
            sys.exit("BC 102 chimera variable-block not supported (extend script if needed)")
        else:
            sys.exit(f"Unknown BC code {code} at line {idx} — extend SINGLE_PROP_CODES/NO_PROP_CODES")
        records.append((head, tuple(props)))
    return records


def coarse_index(i_fine, rap):
    return (i_fine + rap - 1) // rap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="fine-grid bc.txt path")
    ap.add_argument("output", help="coarse-grid bcN.txt path")
    ap.add_argument("--rap", type=int, default=2,
                    help="restriction factor (default 2 = level 1->2)")
    args = ap.parse_args()

    fine = parse_bc_file(args.input)

    seen = {}
    order = []
    for head, props in fine:
        block, i, j, k, face, code = head
        i_c = coarse_index(i, args.rap)
        j_c = coarse_index(j, args.rap)
        k_c = coarse_index(k, args.rap)
        coarse_head = (block, i_c, j_c, k_c, face, code)
        if coarse_head not in seen:
            order.append(coarse_head)
            seen[coarse_head] = props
        else:
            if seen[coarse_head] != props:
                print(f"[WARN] coarse {coarse_head} props mismatch — keeping first.",
                      file=sys.stderr)

    with open(args.output, "w") as f:
        for head in order:
            block, i_c, j_c, k_c, face, code = head
            # Mirror fine bc.txt 6-field format (Q2D uses 5- or 6-field; default 6 with k explicit).
            f.write(f"  {block:>5d}  {i_c:>5d}  {j_c:>5d}  {k_c:>5d}  {face:>5d}  {code:>5d}\n")
            for p in seen[head]:
                f.write(f"{p}\n")

    print(f"Fine records: {len(fine)}")
    print(f"Coarse records (post-dedup): {len(order)}")
    code_counts = {}
    for head in order:
        code = head[5]
        code_counts[code] = code_counts.get(code, 0) + 1
    print(f"Coarse code distribution: {sorted(code_counts.items())}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
