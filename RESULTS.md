# Verified results

All 67 target cases were independently replayed from empty inventories. Every reaction is allowed by its input, gross reactants are available or buyable, only missing buyable reactants are purchased, every requested quantity is reached, and paragraph separators are valid. No case became more expensive.

| Input | Cases | Original program | Previous GitHub outputs | New verified cost | Certified lower bound |
| --- | ---: | ---: | ---: | ---: | ---: |
| example | 2 | 136 | 136 | 136 | 116 |
| small | 10 | 34,038 | 28,576 | 26,684 | 20,926 |
| medium | 15 | 26,762 | 22,366 | 21,316 | 13,837 |
| large | 20 | 36,529 | 27,094 | 24,970 | 17,550 |
| enormous | 20 | 40,636 | 34,421 | 31,394 | 19,256 |

The four full inputs improved from **112,457 to 104,364** (7.20% less than the previous GitHub outputs, and 24.35% less than the original program). **64 of 65 cases strictly improved**; the other case and the example retained their previous costs.

The previous GitHub baseline is commit `b082cbcba50f221900feb81dc9d2169702600bb7`. The original program was evaluated with a 60-second budget. New results combine 240-second and 180-second C++ searches, exact per-case replay, cross-case production-batch reuse, and optional integer-planning refinement from actual prefix inventories. Final restricted refinement continued until a complete round produced no additional savings on each input.

These are **best verified solutions, not proven global minima**. Lower bounds come from input-SHA-bound molecular potentials, with every inequality verified using exact rational arithmetic. The lower bounds remain below achieved costs, so they do not certify global optimality.

Validation completed:

- The exact standalone C++17 source compiled with `-Wall -Wextra` and passed 40 regression tests.
- An AddressSanitizer and UndefinedBehaviorSanitizer build completed an actual-input refinement and resume check without diagnostics.
- Independent replay checked all 67 output cases and per-case nonregression against the previous commit.
- All 67 rational lower-bound certificates passed exact verification; tampered outputs and certificates were rejected in separate checks.
- Download archives are checked for ZIP integrity and byte equality with their committed source files.

Run `python3 optimizer/validate_solutions.py --metadata solutions/validation.json` and `python3 optimizer/verify_bounds.py` to verify the release. The C++ executable itself has no Python or external optimization-library dependency.
