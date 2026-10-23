#!/usr/bin/env python3
"""Independently replay Level 3 outputs, optionally requiring no cost regressions.

This validator deliberately does not import the optimizer's parser or tests.
Only the Python standard library is used. Run from any working directory:

    python3 optimizer/validate_solutions.py
    python3 optimizer/validate_solutions.py --outputs /tmp/candidates --compare solutions
    python3 optimizer/validate_solutions.py --metadata solutions/validation.json

The comparison is per target case, not just the sum for each input.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import sys


class InvalidSolution(ValueError):
    pass


@dataclass
class Input:
    prices: dict[str, int]
    reactions: dict[str, tuple[dict[str, int], dict[str, int]]]
    targets: list[tuple[int, str]]


def side(text: str) -> dict[str, int]:
    amounts: dict[str, int] = defaultdict(int)
    for term in text.split(","):
        fields = term.split()
        if len(fields) != 2 or not fields[0].isdigit() or int(fields[0]) <= 0:
            raise InvalidSolution(f"Invalid reaction term: {term!r}")
        amounts[fields[1]] += int(fields[0])
    return dict(amounts)


def load_input(path: Path) -> Input:
    lines = iter(x.strip() for x in path.read_text().splitlines() if x.strip())
    prices = {}
    reactions = {}
    try:
        for _ in range(int(next(lines))):
            amount, molecule = next(lines).split("->")
            value = int(amount)
            name = molecule.strip()
            if value < 0 or len(name.split()) != 1 or name in prices:
                raise InvalidSolution("Invalid or duplicate price entry")
            prices[name] = value
        for _ in range(int(next(lines))):
            text = next(lines)
            a, b = text.split("->")
            reactions[text] = side(a), side(b)
        targets = []
        for _ in range(int(next(lines))):
            amount, name = next(lines).split()
            if int(amount) <= 0:
                raise InvalidSolution("This reaction-paragraph format requires positive targets")
            targets.append((int(amount), name))
        if next(lines, None) is not None:
            raise InvalidSolution("Extra data after targets")
    except (StopIteration, ValueError) as exc:
        raise InvalidSolution(f"Cannot parse {path}: {exc}") from exc
    return Input(prices, reactions, targets)


def formula(name: str) -> Counter[str]:
    atoms: Counter[str] = Counter()
    end = 0
    for match in re.finditer(r"([A-Z][a-z]?)([0-9]*)", name):
        if match.start() != end:
            raise InvalidSolution(f"Unsupported molecular formula: {name}")
        atoms[match[1]] += int(match[2] or "1")
        end = match.end()
    if end != len(name):
        raise InvalidSolution(f"Unsupported molecular formula: {name}")
    return atoms


def atom_audit(problem: Input) -> dict:
    """Describe exact element conservation; it is not an extra legality rule."""
    formulas = {}
    elements = set()
    for ins, outs in problem.reactions.values():
        for name in ins.keys() | outs.keys():
            if name not in formulas:
                formulas[name] = formula(name)
                elements.update(formulas[name])
    conserved = set(elements)
    unbalanced = 0
    for ins, outs in problem.reactions.values():
        change: Counter[str] = Counter()
        for name, amount in ins.items():
            for element, count in formulas[name].items():
                change[element] -= amount * count
        for name, amount in outs.items():
            for element, count in formulas[name].items():
                change[element] += amount * count
        broken = {element for element, value in change.items() if value}
        conserved.difference_update(broken)
        unbalanced += bool(broken)
    return {"conserved_elements": sorted(conserved), "unbalanced_reactions": unbalanced}


def replay(problem: Input, output: Path) -> list[dict]:
    text = output.read_text().replace("\r\n", "\n").strip()
    blocks = re.split(r"\n[ \t]*\n", text) if text else []
    if len(blocks) != len(problem.targets):
        raise InvalidSolution(
            f"{output}: expected {len(problem.targets)} paragraphs, found {len(blocks)}"
        )
    results = []
    for case, (block, (required, target)) in enumerate(zip(blocks, problem.targets), 1):
        inventory: dict[str, int] = defaultdict(int)
        purchases: dict[str, int] = defaultdict(int)
        cost = 0
        actions = block.splitlines()
        for step, line in enumerate(actions, 1):
            reaction = problem.reactions.get(line.strip())
            if reaction is None:
                raise InvalidSolution(f"{output}, case {case}, step {step}: unknown reaction")
            ins, outs = reaction
            # Gross inputs are consumed first. No product is available early,
            # even when the same catalyst also appears on the output side.
            for name, amount in ins.items():
                missing = max(0, amount - inventory[name])
                if missing:
                    if name not in problem.prices:
                        raise InvalidSolution(
                            f"{output}, case {case}, step {step}: missing {missing} {name}"
                        )
                    purchases[name] += missing
                    cost += missing * problem.prices[name]
                    inventory[name] += missing
                inventory[name] -= amount
                if inventory[name] < 0:
                    raise InvalidSolution("Internal validator inventory error")
            for name, amount in outs.items():
                inventory[name] += amount
        if inventory[target] < required:
            raise InvalidSolution(
                f"{output}, case {case}: final {target}={inventory[target]}, required {required}"
            )
        # Python integers are unbounded; no fixed-width cost overflow is hidden.
        results.append({"case": case, "target": target, "required": required,
                        "produced": inventory[target], "cost": cost,
                        "reactions": len(actions), "purchases": dict(sorted(purchases.items()))})
    return results


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=root / "inputs")
    parser.add_argument("--outputs", type=Path, default=root / "solutions")
    parser.add_argument("--compare", type=Path, help="directory of previous valid outputs")
    parser.add_argument("--metadata", type=Path, help="require saved JSON cost totals to match")
    parser.add_argument("--report", type=Path, help="write full independent audit as JSON")
    args = parser.parse_args()
    inputs = sorted(args.inputs.glob("in_level-3_*.txt"))
    if not inputs:
        parser.error("No Level 3 inputs found")
    metadata = None
    if args.metadata:
        metadata = {row["output"]: row for row in json.loads(args.metadata.read_text())}
    report = []
    for input_path in inputs:
        problem = load_input(input_path)
        output = args.outputs / input_path.name.replace("in_", "out_", 1)
        cases = replay(problem, output)
        costs = [case["cost"] for case in cases]
        row = {"input": input_path.name, "output": output.name,
               "input_sha256": sha256(input_path), "output_sha256": sha256(output),
               "cost": sum(costs), "cases": cases, **atom_audit(problem)}
        if args.compare:
            previous = replay(problem, args.compare / output.name)
            old_costs = [case["cost"] for case in previous]
            for case, (current, old) in enumerate(zip(costs, old_costs), 1):
                if current > old:
                    raise InvalidSolution(f"{output}, case {case}: regression {old} -> {current}")
            row["previous_cost"] = sum(old_costs)
            row["improved_cases"] = sum(a < b for a, b in zip(costs, old_costs))
        if metadata is not None:
            saved = metadata.get(output.name)
            if saved is None or saved["case_costs"] != costs or saved["improved_cost"] != sum(costs):
                raise InvalidSolution(f"{output}: metadata costs do not match replay")
        report.append(row)
        print(f"VALID {output.name}: {len(cases)} cases, cost {sum(costs)}")
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"VALIDATED {sum(len(row['cases']) for row in report)} cases; "
          f"total cost {sum(row['cost'] for row in report)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (InvalidSolution, OSError, KeyError, TypeError) as error:
        print(f"INVALID: {error}", file=sys.stderr)
        raise SystemExit(1)
