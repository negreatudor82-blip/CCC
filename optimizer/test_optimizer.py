#!/usr/bin/env python3
"""Run independent cost/legality checks against a Level 3 solver.

Usage: python3 test_optimizer.py ./factory [--original /path/to/original]
Only the Python standard library is needed. The solver interface is:
    executable input_path output_path search_seconds
"""

from __future__ import annotations

import argparse
import heapq
import random
import re
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Problem:
    prices: dict[str, int]
    reactions: dict[str, tuple[dict[str, int], dict[str, int]]]
    targets: list[tuple[int, str]]


def parse_side(side: str) -> dict[str, int]:
    result: dict[str, int] = defaultdict(int)
    for part in side.split(","):
        fields = part.split()
        if len(fields) != 2 or int(fields[0]) <= 0:
            raise ValueError(f"Malformed reaction component: {part!r}")
        result[fields[1]] += int(fields[0])
    return dict(result)


def parse_problem(path: Path) -> Problem:
    lines = iter(line.strip() for line in path.read_text().splitlines() if line.strip())
    prices: dict[str, int] = {}
    for _ in range(int(next(lines))):
        cost, molecule = next(lines).split("->")
        prices[molecule.strip()] = int(cost)
    reactions = {}
    for _ in range(int(next(lines))):
        line = next(lines)
        left, right = line.split("->")
        reactions[line] = (parse_side(left), parse_side(right))
    targets = []
    for _ in range(int(next(lines))):
        quantity, molecule = next(lines).split()
        targets.append((int(quantity), molecule))
    if list(lines):
        raise ValueError("Unexpected extra input lines")
    return Problem(prices, reactions, targets)


def validate_output(problem: Problem, path: Path) -> list[int]:
    """Validate only emitted reactions; never trust solver-reported costs."""
    text = path.read_text().replace("\r\n", "\n").strip()
    paragraphs = re.split(r"\n[ \t]*\n", text) if text else []
    if len(paragraphs) != len(problem.targets):
        raise AssertionError(
            f"Expected {len(problem.targets)} reaction paragraphs separated by blank lines; "
            f"found {len(paragraphs)}"
        )
    costs = []
    for number, (paragraph, (needed, target)) in enumerate(
        zip(paragraphs, problem.targets), 1
    ):
        inventory: dict[str, int] = defaultdict(int)
        cost = 0
        for line in paragraph.splitlines():
            line = line.strip()
            if line not in problem.reactions:
                raise AssertionError(f"Case {number}: unknown reaction {line!r}")
            consumed, produced = problem.reactions[line]
            for molecule, quantity in consumed.items():
                shortfall = max(0, quantity - inventory[molecule])
                if shortfall:
                    if molecule not in problem.prices:
                        raise AssertionError(
                            f"Case {number}: unavailable, non-buyable input {molecule}"
                        )
                    cost += shortfall * problem.prices[molecule]
                    inventory[molecule] += shortfall
                inventory[molecule] -= quantity
            for molecule, quantity in produced.items():
                inventory[molecule] += quantity
        if inventory[target] < needed:
            raise AssertionError(
                f"Case {number}: produced {inventory[target]} {target}, needed {needed}"
            )
        costs.append(cost)
    return costs


def run_solver(executable: Path, path: Path, seconds: float = 0.3):
    with tempfile.TemporaryDirectory(prefix="factory-test-") as directory:
        output = Path(directory) / "solution.txt"
        result = subprocess.run(
            [str(executable), str(path), str(output), str(seconds)],
            capture_output=True,
            text=True,
            timeout=max(10, seconds + 5),
        )
        contents = output.read_text() if output.exists() else ""
        if result.returncode:
            return result, None, contents
        if not output.exists():
            raise AssertionError("Solver returned success without writing its output")
        return result, validate_output(parse_problem(path), output), contents


def random_batches(seed: int):
    """Single-product batch problems have a small, exact Dijkstra reference."""
    rng = random.Random(seed)
    target = rng.randint(6, 22)
    recipes = [(rng.randint(1, 12), rng.randint(1, 8)) for _ in range(4)]
    # A costs 1. Each pair is (input cost, target yield).
    lines = ["1", "1 -> A", "", str(len(recipes))]
    lines += [f"{cost} A -> {quantity} T" for cost, quantity in recipes]
    lines += ["", "1", f"{target} T", ""]
    distances = {0: 0}
    queue = [(0, 0)]
    while queue:
        cost, amount = heapq.heappop(queue)
        if distances[amount] != cost:
            continue
        if amount >= target:
            return "\n".join(lines), cost
        for input_cost, yield_ in recipes:
            next_amount = amount + yield_
            next_cost = cost + input_cost
            if next_cost < distances.get(next_amount, float("inf")):
                distances[next_amount] = next_cost
                heapq.heappush(queue, (next_cost, next_amount))
    raise AssertionError("Batch reference search unexpectedly exhausted")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("--original", type=Path)
    parser.add_argument("--random-cases", type=int, default=12)
    args = parser.parse_args()
    solver = args.solver.resolve()
    original = args.original.resolve() if args.original else None
    fixtures = Path(__file__).resolve().parent / "tests"
    tests = [
        ("example", [72, 64]),
        ("mixed_batches", [9]),
        ("seeded_growth", [14]),
        ("catalyst_return", [15]),
        ("byproducts", [9]),
        ("buy_remainder", [7]),
        ("indirect_coproduct_order", [2]),
        ("gross_zero_target", [2]),
        ("zero_cost", [0]),
        ("unreachable", None),
    ]
    failures = 0
    passed = 0

    def check(name, path, expected):
        nonlocal failures, passed
        try:
            result, costs, contents = run_solver(solver, path, 1 if name == "example" else 0.3)
            if expected is None:
                if result.returncode == 0:
                    raise AssertionError("Unreachable target must return a nonzero exit status")
                if contents.strip():
                    raise AssertionError("Unreachable target must not write a bogus solution")
                detail = "correctly rejected"
            else:
                if result.returncode:
                    raise AssertionError(f"Solver exited {result.returncode}: {result.stderr.strip()}")
                if costs != expected:
                    raise AssertionError(f"Costs {costs}; expected optimum {expected}")
                detail = f"cost {costs}"
            if original and expected is not None:
                try:
                    old_result, old_costs, _ = run_solver(original, path, 0.1)
                    detail += f"; original {old_costs if old_result.returncode == 0 else 'failed'}"
                except (AssertionError, subprocess.TimeoutExpired) as exc:
                    detail += f"; original invalid ({exc})"
            print(f"PASS {name}: {detail}")
            passed += 1
        except (AssertionError, ValueError, OSError, subprocess.TimeoutExpired) as exc:
            print(f"FAIL {name}: {exc}")
            failures += 1

    for name, expected in tests:
        check(name, fixtures / f"{name}.txt", expected)
    with tempfile.TemporaryDirectory(prefix="factory-random-") as directory:
        for seed in range(args.random_cases):
            text, optimum = random_batches(seed)
            path = Path(directory) / f"random_{seed}.txt"
            path.write_text(text)
            check(f"random_batches_{seed}", path, [optimum])
    print(f"{passed} passed, {failures} failed")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
