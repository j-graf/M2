# Systematic SymmetricRings benchmarks

This directory is independent of the historical `../benchmarks.m2` notebook.
The notebook remains unchanged and is not loaded or interpreted by this system.

The catalog separates mathematical inputs from executable operations:

- `partitions.m2` defines stable named partitions and pair families.
- `cases.m2` combines those inputs with operation, ring, and timing-tier data.
- `operations.m2` is the only file that translates operation names into code.
- `runner.m2` executes or lists one-process benchmark cases.
- `run-benchmarks.sh` creates a fresh M2 process for each cold repetition.
- `validate.m2` checks identifiers, references, operations, tiers, and weights.
- `summarize-results.awk` calculates per-case minimum, median, and maximum CPU time.
- `compare-results.awk` compares current medians with the best valid and most
  recent accepted baselines.
- `accept-summary.awk` prepares reviewed baseline rows.
- `collect-system-info.sh` records hardware, OS, M2, build, and Git metadata.
- `make-report.sh` combines system metadata and comparisons into Markdown.
- `baselines.tsv` is reserved for reviewed, accepted historical medians.

## Running

From the nested Macaulay2 source tree, first install the current package:

```sh
CCACHE_DISABLE=1 cmake --build BUILD/build --target install-SymmetricRings -j2
```

List the complete catalog:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh --list
```

Validate the catalog without running its computations:

```sh
BUILD/build/M2 --script \
  Macaulay2/packages/SymmetricRings/extras/benchmarks/validate.m2
```

Run one case three times in fresh processes:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --case pleth-three-row2-combined
```

Run a family or timing tier:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --family SchurPlethysm --tier Medium --repetitions 3

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --tier Small --repetitions 3 --output /tmp/symmetricrings-small.tsv
```

The supported tiers are `Small`, `Medium`, `Large`, and `Stress`. The initial
catalog avoids putting correctness smoke tests into the timing data; package
tests remain the correctness authority. The runner nevertheless checks the
expected homogeneous result weight by default. Use `--no-verify` only for
diagnostic measurements.

## Result format

Results are append-only TSV records with the following columns:

```text
case_id family operation tier coefficient_ring repetition
cpu_seconds wall_seconds result_terms result_weight input_group
lambda mu probe lambda_weight lambda_length mu_weight mu_length
expected_weight pair_class
```

`input_group` ties combined computations to their split stages. For example,
the three records for a plethysm group time combined `@`, production in power
sums, and subsequent conversion to Schur separately.

Each completed run automatically writes three files:

```text
results-TIMESTAMP.tsv
results-TIMESTAMP-summary.tsv
results-TIMESTAMP-comparison.tsv
results-TIMESTAMP-system.tsv
results-TIMESTAMP-report.md
```

The baseline comparison is also printed in the terminal. Cases are classified
as `stable`, `improvement`, `regression`, or `new`. It reports changes from
both the best valid historical median and the most recent accepted median. The
default classification threshold is ten percent; set `THRESHOLD_PERCENT` to
override it. Use `--baselines FILE` to compare against a different baseline
table.

Summarize a result file:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/summarize-results.awk \
  results-20260712-183000.tsv > summary.tsv
```

Use the median CPU time as the primary comparison. Keep the individual runs,
and use at least five repetitions before classifying a difference below ten
percent as a regression. Warm-cache studies should be separate case families;
the standard runner intentionally provides cold M2 processes.

Cases remain `new` until a baseline is deliberately accepted. After reviewing
a generated summary, prepare rows for appending with:

```sh
ACCEPTED_DATE=2026-07-12 NOTES='initial systematic baseline' \
  Macaulay2/packages/SymmetricRings/extras/benchmarks/accept-summary.awk \
  results-TIMESTAMP-summary.tsv >> \
  Macaulay2/packages/SymmetricRings/extras/benchmarks/baselines.tsv
```

Running benchmarks never silently changes their baseline.

The system report records the computer model, CPU/chip, core counts, RAM,
operating system and kernel, architecture, M2 version and executable, CMake
build type, Git commit and branch, and whether the worktree is dirty. It does
not record hardware serial numbers, UUIDs, or provisioning identifiers.

## Extending the catalog

Prefer adding a named mathematical shape to `partitions.m2` and referring to
its stable ID from cases. Add raw partition lists directly to a case only for a
localized crossover grid. New operation families receive a symbolic operation
name in `cases.m2` and one implementation branch in `operations.m2`.

Do not add implementation-route names to case IDs. Cases describe mathematics;
the dispatcher remains free to change algorithms. Forced-route experiments
belong in separate diagnostic runs and should not be accepted as production
baselines.
