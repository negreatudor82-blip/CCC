# Verified results

All 67 target cases were independently replayed from empty inventories. Every output uses allowed reactions, satisfies gross reactant availability and target quantities, and retains mandatory case separators. No case became more expensive.

| Input | Cases | Previous cost | New verified cost | Certified lower bound |
| --- | ---: | ---: | ---: | ---: |
| example | 2 | 136 | 136 | 116 |
| small | 10 | 26,684 | 26,660 | 20,926 |
| medium | 15 | 21,316 | 21,296 | 13,837 |
| large | 20 | 24,970 | 24,949 | 17,550 |
| enormous | 20 | 31,394 | 31,227 | 19,256 |

Total cost decreased from **104,500 to 104,268** (232 saved, 0.222%). **19 of 67 cases strictly improved**; all other cases retained their costs. The baseline is commit `c5a6fa17006a146d01cc523c917be0fc26e59f04`.

Candidates were generated with `optimizer/refine_with_milp.py`: an initial 1-second-per-solve pass at prefix fractions `.7,.35,0`, followed by restricted refinement with 0.4-second solve budgets, up to two rounds, and fractions `.95,.9,.85,.75,.65,.55,.5,.45,.3,.25,.15,.1,.05`. The restricted searches were stopped after retaining their latest complete candidates. `optimizer/merge_solutions.py` selected the cheapest independently verified plan per case. These are best verified solutions, not proven global minima.

The C++17 executable passed all 40 regression tests. Independent replay validated all 67 cases and per-case nonregression. All 67 unchanged rational lower-bound certificates passed exact verification.

Run `python3 optimizer/validate_solutions.py --metadata solutions/validation.json` and `python3 optimizer/verify_bounds.py` to verify the results.
