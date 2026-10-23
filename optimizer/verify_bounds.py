#!/usr/bin/env python3
"""Verify purchase-cost lower bounds using exact rational arithmetic only.

No optimization package is needed for verification. The certificate binds each
input by SHA-256 and gives nonnegative molecular potentials y satisfying:

    sum(outputs * y) <= sum(inputs * y)       for every reaction
    y[molecule] <= purchase_price[molecule]  for every buyable molecule

Any valid sequence starts empty, so its purchase cost is at least
ceil(target_quantity * y[target]). This is a lower bound, not an executable plan
or a proof that the accompanying solution attains the minimum.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

from validate_solutions import InvalidSolution, load_input


def integer(value, description: str) -> int:
    if type(value) is not int:
        raise InvalidSolution(f"{description} must be an integer")
    if value.bit_length() > 4096:
        raise InvalidSolution(f"{description} is unreasonably large")
    return value


def rational(value, description: str) -> Fraction:
    if not isinstance(value, list) or len(value) != 2:
        raise InvalidSolution(f"{description} must be [numerator, denominator]")
    numerator = integer(value[0], description)
    denominator = integer(value[1], description)
    if denominator <= 0:
        raise InvalidSolution(f"{description} has a nonpositive denominator")
    result = Fraction(numerator, denominator)
    if result < 0:
        raise InvalidSolution(f"{description} must be nonnegative")
    return result


def verify_certificate(certificate: dict, inputs: Path) -> list[dict]:
    if certificate.get("schema_version") != 1:
        raise InvalidSolution("Unsupported lower-bound certificate version")
    entries = certificate.get("inputs")
    if not isinstance(entries, list) or not entries:
        raise InvalidSolution("Missing certificate inputs")
    expected_names = {p.name for p in inputs.glob("in_level-3_*.txt")}
    actual_names = [entry.get("input") for entry in entries]
    if len(set(actual_names)) != len(actual_names) or set(actual_names) != expected_names:
        raise InvalidSolution("Certificate input membership differs from the supplied inputs")
    report = []
    for entry in entries:
        name = entry["input"]
        if Path(name).name != name:
            raise InvalidSolution("Certificate input must be a filename")
        path = inputs / name
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if entry.get("input_sha256") != digest:
            raise InvalidSolution(f"{name}: input SHA-256 does not match certificate")
        problem = load_input(path)
        known = set(problem.prices)
        for ins, outs in problem.reactions.values():
            known.update(ins)
            known.update(outs)
        known.update(target for _, target in problem.targets)
        cases = entry.get("cases")
        if not isinstance(cases, list) or len(cases) != len(problem.targets):
            raise InvalidSolution(f"{name}: certificate case count mismatch")
        lower_bounds = []
        for number, (saved, (required, target)) in enumerate(zip(cases, problem.targets), 1):
            label = f"{name}, case {number}"
            if (saved.get("case") != number or saved.get("target") != target
                    or saved.get("required") != required):
                raise InvalidSolution(f"{label}: target binding mismatch")
            encoded = saved.get("potentials")
            if not isinstance(encoded, dict) or not set(encoded).issubset(known):
                raise InvalidSolution(f"{label}: unknown molecular potential")
            potential = {molecule: rational(value, f"{label}, {molecule}")
                         for molecule, value in encoded.items()}
            zero = Fraction(0)
            for molecule, price in problem.prices.items():
                if potential.get(molecule, zero) > price:
                    raise InvalidSolution(f"{label}: potential exceeds price of {molecule}")
            for reaction, (ins, outs) in problem.reactions.items():
                consumed = sum(quantity * potential.get(molecule, zero)
                               for molecule, quantity in ins.items())
                produced = sum(quantity * potential.get(molecule, zero)
                               for molecule, quantity in outs.items())
                if produced > consumed:
                    raise InvalidSolution(f"{label}: potential increases in {reaction}")
            exact = required * potential.get(target, zero)
            bound = math.ceil(exact)
            if integer(saved.get("lower_bound"), f"{label}, lower_bound") != bound:
                raise InvalidSolution(f"{label}: claimed bound differs from exact calculation")
            lower_bounds.append(bound)
        report.append({"input": name, "cases": len(cases), "lower_bound": sum(lower_bounds),
                       "case_lower_bounds": lower_bounds})
    return report


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", nargs="?", type=Path,
                        default=root / "solutions/lower_bounds.json")
    parser.add_argument("--inputs", type=Path, default=root / "inputs")
    args = parser.parse_args()
    report = verify_certificate(json.loads(args.certificate.read_text()), args.inputs)
    for row in report:
        print(f"CERTIFIED {row['input']}: {row['cases']} cases, lower bound {row['lower_bound']}")
    print(f"CERTIFIED {sum(row['cases'] for row in report)} cases; "
          f"total lower bound {sum(row['lower_bound'] for row in report)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (InvalidSolution, OSError, ValueError, KeyError, TypeError) as error:
        print(f"INVALID CERTIFICATE: {error}", file=sys.stderr)
        raise SystemExit(1)
