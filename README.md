# Molecule Factory — Level 3

`factory_optimized.cpp` is a complete, standalone C++17 source file. It needs no
Python packages, external solver, or extra headers.

[Download all output files](all_output_files.zip), browse the
[individual solutions](solutions/), or open the [C++ source](factory_optimized.cpp).
The supplied inputs are in [inputs/](inputs/) and the original problem and code
are preserved in [reference/](reference/).

```sh
g++ -O2 -std=c++17 factory_optimized.cpp -o factory
./factory inputs/in_level-3_4-enormous.txt solutions/out_level-3_4-enormous.txt 120
```

Arguments are `INPUT [OUTPUT] [SECONDS=60] [SEED=987654321]`. The time budget is
for the complete input file. Omit OUTPUT to derive its name from INPUT. A matching
existing output is independently checked and used as a starting solution, so you
can rerun with more time or a different seed while retaining the best result for
each target:

```sh
./factory inputs/in_level-3_4-enormous.txt solutions/out_level-3_4-enormous.txt 300 2026
```

The output file contains only allowed reaction lines and the mandatory blank
line between target cases. Progress and measured costs go to the console.
Every case starts with an empty inventory. A reaction buys only missing buyable
inputs; all other inputs must already exist. The finished sequence is replayed
before writing. If any case has no valid solution, the program reports failure
and leaves the existing output untouched.

The optimizer retains the original recipe search and adds recipe selection by
quantity and inventory, mixed production and purchases, catalyst/seed handling,
integer reaction-count neighborhoods, and replacement of a sequence's remainder
from its actual intermediate inventory. It also reuses valid production batches
across cases, reorders reactions and removes redundant production. Each case
still starts with an empty inventory. Small networks use exact inventory search. Only
an exhausted exact search is labelled `optimal`; larger searches report
`best found`. More runtime can improve a result but does not prove a global
minimum.

The four supplied inputs contain 65 target cases. Their provided solutions were
also checked by an independent Python inventory simulator. Regression tests
cover mixed batch sizes, buying a remainder, catalyst recycling, seeded growth,
byproducts, zero costs, and unreachable targets. The original example retains
its optimal costs of 72 and 64.

The `optimizer/` directory contains modular development sources and independent
tests; it is not needed to compile `factory_optimized.cpp`.

To run the independent regression checks (Python standard library only):

```sh
python3 optimizer/test_optimizer.py ./factory --random-cases 30
```

See [RESULTS.md](RESULTS.md) for verified costs and their limitations.

The independent release checks also replay every saved solution, reject any
per-case cost regression, and verify rational lower-bound certificates with
exact arithmetic:

```sh
python3 optimizer/validate_solutions.py --metadata solutions/validation.json
python3 optimizer/verify_bounds.py
```

Lower bounds measure what further savings might remain possible; they do not
mean that the saved outputs achieve the global minimum. Certificate verification
uses only Python's standard library. Regenerating the certificates with
`optimizer/compute_bounds.py` additionally needs SciPy.

An optional offline refinement helper explores additional integer-planning
neighborhoods from intermediate inventories. It needs NumPy and SciPy; these
packages are not required by the C++ executable. Like the C++ program, it retains
each existing case unless full inventory replay finds a cheaper valid sequence:

```sh
python3 optimizer/refine_with_milp.py inputs/in_level-3_4-enormous.txt \
    solutions/out_level-3_4-enormous.txt refined-enormous.txt \
    --seconds-per-solve 5 --rounds 3 --restricted-only \
    --fractions .8,.6,.4,.2,0
```

Inspect and independently validate candidate outputs before replacing published
files. `optimizer/merge_solutions.py` can merge complete candidate directories,
retaining the cheapest verified result separately for each target case.

To rebuild the standalone C++ file after changing modular sources:

```sh
python3 optimizer/amalgamate.py
python3 optimizer/amalgamate.py --check
```
