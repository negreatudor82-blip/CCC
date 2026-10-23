#!/usr/bin/env python3
"""Build the single-file C++17 distribution from the modular optimizer."""
from pathlib import Path
import argparse


PREAMBLE = """// Molecule Factory, Level 3 -- standalone C++17 cost optimizer.
// Build: g++ -O2 -std=c++17 factory_optimized.cpp -o factory
// Run: ./factory INPUT [OUTPUT] [SECONDS=60] [SEED=987654321]
// A matching existing OUTPUT is validated and retained as an incumbent.
// Includes inventory search, integer neighborhoods, suffix repair and batch reuse.
// Every accepted plan is replayed from empty inventory. Large-instance results
// are best found; optimality is claimed only by completed exact inventory search.

"""


def build(directory: Path) -> str:
    seen = set()

    def expand(path: Path) -> str:
        if path.name in seen:
            return ""
        seen.add(path.name)
        output = []
        for line in path.read_text().splitlines():
            if line.startswith('#include "'):
                name = line.split('"')[1]
                child = directory / name
                if child.resolve().parent != directory.resolve():
                    raise ValueError(f"Include escapes source directory: {name}")
                output.append(f"// ---- {name} ----\n" + expand(child))
            elif line.strip() != "#pragma once":
                output.append(line)
        return "\n".join(output) + "\n"

    text = expand(directory / "factory.cpp")
    return PREAMBLE + text[text.index("#include <bits/stdc++.h>"):]


def main() -> int:
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=directory.parent / "factory_optimized.cpp")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build(directory)
    if args.check:
        if not args.output.exists() or args.output.read_text() != result:
            print(f"STALE: {args.output}; run optimizer/amalgamate.py")
            return 1
        print("Standalone source matches modular sources")
        return 0
    args.output.write_text(result)
    print(f"Wrote {args.output} ({len(result.encode())} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
