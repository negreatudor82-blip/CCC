#!/usr/bin/env python3
"""Merge independently validated candidates, retaining the cheapest case plans.

Nothing is written until every discovered candidate has been replayed and every
chosen case checked against its incumbent. This script has no Git side effects.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import tempfile

from validate_solutions import InvalidSolution, load_input, replay


def merge(inputs: Path, incumbent: Path, candidates: list[Path], destination: Path) -> list[dict]:
    prepared = []
    report = []
    for input_path in sorted(inputs.glob("in_level-3_*.txt")):
        problem = load_input(input_path)
        filename = input_path.name.replace("in_", "out_", 1)
        baseline = incumbent / filename
        baseline_cases = replay(problem, baseline)
        baseline_blocks = re.split(r"\n[ \t]*\n", baseline.read_text().strip())
        selected = [(case["cost"], block, str(baseline))
                    for case, block in zip(baseline_cases, baseline_blocks)]
        for directory in candidates:
            path = directory / filename
            if not path.exists():
                continue
            # Freeze one complete version if a refinement process is still
            # replacing its output atomically while the merger runs.
            with tempfile.TemporaryDirectory(prefix="factory-candidate-") as directory:
                snapshot = Path(directory) / filename
                snapshot.write_bytes(path.read_bytes())
                actual = replay(problem, snapshot)
                blocks = re.split(r"\n[ \t]*\n", snapshot.read_text().strip())
            for i, (case, block) in enumerate(zip(actual, blocks)):
                if case["cost"] < selected[i][0]:
                    selected[i] = case["cost"], block, str(path)
        text = "\n\n".join(block for _, block, _ in selected) + "\n"
        # Check the serialization, including paragraph boundaries, with the independent replay.
        with tempfile.TemporaryDirectory(prefix="factory-merge-") as directory:
            path = Path(directory) / filename
            path.write_text(text)
            checked = replay(problem, path)
        costs = [case["cost"] for case in checked]
        if any(a > b["cost"] for a, b in zip(costs, baseline_cases)):
            raise InvalidSolution(f"{filename}: merger attempted a cost regression")
        prepared.append((filename, text))
        report.append({"input": input_path.name, "output": filename,
                       "previous_cost": sum(x["cost"] for x in baseline_cases),
                       "cost": sum(costs), "case_costs": costs,
                       "improved_cases": sum(a < b["cost"] for a, b in zip(costs, baseline_cases)),
                       "sources": [source for _, _, source in selected]})
    if not prepared:
        raise InvalidSolution("No inputs found")
    destination.mkdir(parents=True, exist_ok=True)
    # Each complete output file is installed by an atomic rename in its own directory.
    for filename, text in prepared:
        descriptor, name = tempfile.mkstemp(prefix=".solution-", dir=destination)
        try:
            with os.fdopen(descriptor, "w") as stream:
                stream.write(text)
            os.replace(name, destination / filename)
        finally:
            if Path(name).exists():
                Path(name).unlink()
    return report


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=root / "inputs")
    parser.add_argument("--incumbent", type=Path, default=root / "solutions")
    parser.add_argument("--candidate", type=Path, action="append", default=[])
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    rows = merge(args.inputs, args.incumbent, args.candidate, args.destination)
    if args.report:
        args.report.write_text(json.dumps(rows, indent=2) + "\n")
    for row in rows:
        print(f"{row['output']}: {row['previous_cost']} -> {row['cost']}; "
              f"{row['improved_cases']} cases improved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
