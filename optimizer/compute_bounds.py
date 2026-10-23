#!/usr/bin/env python3
"""Generate exact-checkable stoichiometric purchase-cost lower bounds.

Generation needs NumPy and SciPy. Verification by verify_bounds.py uses only the
Python standard library. The C++ optimizer has no dependency on this script.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from validate_solutions import InvalidSolution, load_input
from verify_bounds import verify_certificate


def generate(inputs: Path) -> dict:
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import lil_matrix

    certificate = {
        "schema_version": 1,
        "description": "Nonnegative stoichiometric dual potentials; exact rational "
                       "verification proves purchase-cost lower bounds, not attainability.",
        "inputs": [],
    }
    paths = sorted(inputs.glob("in_level-3_*.txt"))
    if not paths:
        raise InvalidSolution("No Level 3 inputs found")
    for path in paths:
        problem = load_input(path)
        molecules = set(problem.prices)
        for ins, outs in problem.reactions.values():
            molecules.update(ins)
            molecules.update(outs)
        molecules.update(target for _, target in problem.targets)
        names = sorted(molecules)
        index = {molecule: i for i, molecule in enumerate(names)}
        reactions = list(problem.reactions.values())
        buyable = list(problem.prices)
        nr = len(reactions)
        matrix = lil_matrix((len(names), nr + len(buyable)), dtype=float)
        for r, (ins, outs) in enumerate(reactions):
            for molecule, quantity in ins.items():
                matrix[index[molecule], r] += quantity
            for molecule, quantity in outs.items():
                matrix[index[molecule], r] -= quantity
        for j, molecule in enumerate(buyable):
            matrix[index[molecule], nr + j] = -1
        matrix = matrix.tocsr()
        cost = np.r_[np.zeros(nr), [problem.prices[molecule] for molecule in buyable]]
        entry = {"input": path.name,
                 "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "cases": []}
        for number, (required, target) in enumerate(problem.targets, 1):
            rhs = np.zeros(len(names))
            rhs[index[target]] = -required
            result = linprog(cost, A_ub=matrix, b_ub=rhs, bounds=(0, None), method="highs")
            if not result.success:
                raise InvalidSolution(f"{path.name}, case {number}: {result.message}")
            # Rational reconstruction is accepted only if the whole finished
            # certificate passes exact verification below. No epsilon proves it.
            potential = [Fraction(float(-value)).limit_denominator(100000)
                         for value in result.ineqlin.marginals]
            exact = required * potential[index[target]]
            encoded = {molecule: [value.numerator, value.denominator]
                       for molecule, value in zip(names, potential) if value}
            entry["cases"].append({"case": number, "target": target, "required": required,
                                   "lower_bound": math.ceil(exact), "potentials": encoded})
        certificate["inputs"].append(entry)
    verify_certificate(certificate, inputs)
    return certificate


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=root / "inputs")
    parser.add_argument("--output", type=Path, default=root / "solutions/lower_bounds.json")
    args = parser.parse_args()
    certificate = generate(args.inputs)
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    print(f"Saved exact-verified certificate to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
