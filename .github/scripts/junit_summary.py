#!/usr/bin/env python3
"""Render a CTest JUnit XML file as a GitHub Actions step-summary table."""
import sys
import xml.etree.ElementTree as ET


def main(path):
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as exc:
        print(f"_No JUnit report available ({exc})._")
        return 0

    suites = [root] if root.tag == "testsuite" else root.findall("testsuite")
    rows, n_pass = [], 0
    n_fail = n_skip = 0
    for suite in suites:
        for case in suite.findall("testcase"):
            name = case.get("name", "?")
            time = float(case.get("time", 0.0))
            if case.find("failure") is not None or case.find("error") is not None:
                status, n_fail = "❌ fail", n_fail + 1
            elif case.find("skipped") is not None:
                status, n_skip = "⏭️ skip", n_skip + 1
            else:
                status, n_pass = "✅ pass", n_pass + 1
            rows.append((name, status, time))

    total = len(rows)
    print("## CTest results")
    print()
    print(f"**{n_pass}/{total} passed** · {n_fail} failed · {n_skip} skipped")
    print()
    print("| Test | Result | Time (s) |")
    print("|---|---|---:|")
    for name, status, time in sorted(rows, key=lambda r: -r[2]):
        print(f"| {name} | {status} | {time:.1f} |")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "junit.xml"))
