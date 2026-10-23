# Verified results

All outputs were independently replayed from empty inventories. Every reaction is from its corresponding input, nonbuyable reactants were available, purchases cover only missing buyable reactants, every target quantity was met, and blank paragraph separators were checked.

| Input | Cases | Original cost | Best verified cost |
| --- | ---: | ---: | ---: |
| example | 2 | 136 | 136 |
| small | 10 | 34,038 | 28,576 |
| medium | 15 | 26,762 | 22,366 |
| large | 20 | 36,529 | 27,094 |
| enormous | 20 | 40,636 | 34,421 |

The original program was run with a 60-second budget. Best results combine validated candidates from the revised C++ search, independent integer-optimization experiments, and per-case incumbent retention. These larger-input results are best found, not proven global minima. The standalone C++ program can resume the supplied output files and retain or improve them.

The final standalone C++17 file compiled with -Wall -Wextra and passed 40 regression tests. ASan/UBSan checks on an earlier revision passed 22 tests plus an actual input and resume check.
