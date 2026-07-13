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
- `estimate-run.awk` estimates a selected run without executing its cases.
- `accept-summary.awk` prepares reviewed baseline rows.
- `collect-system-info.sh` records hardware, OS, M2, build, and Git metadata.
- `collect-run-conditions.sh` takes unprivileged pre-run and post-run resource
  snapshots.
- `summarize-run-conditions.awk` assesses snapshots and family calibration
  drift.
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
  --tier Small --repetitions 3 --output symmetricrings-small
```

Run a curated cross-family performance sample:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-fixed light

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-fixed thorough
```

Or sample independently within every eligible family:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-random standard

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-random standard --seed 42
```

The `light`, `standard`, and `thorough` names describe breadth, not timing
tiers. Fixed profiles use nested curated lists, so `light` is a subset of
`standard`, which is a subset of `thorough`. Random profiles select at most 1,
2, or 4 cases per eligible family, respectively. A random seed is generated
and recorded when omitted; supplying `--seed` reproduces the selection.
Family and tier filters may be combined with either varied mode.

Run only cases that have never produced a recorded raw result and have never
appeared in the selected baseline table:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --new --list

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --new --repetitions 3
```

`--new` may be combined with family, tier, or varied selection. It checks case
IDs in the table selected by `--baselines` and in only the most recently
modified `results/*/raw.tsv`. Older non-baseline runs do not affect selection.
The reference run is recorded in `conditions.tsv` and the report.

Estimate the wall time for any selection without running it:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-fixed standard --estimate

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --family SchurPlethysm --tier Medium --repetitions 5 --estimate
```

`--estimate` applies the same case, family, tier, varied, seed, and new-case
filters as a real run. For each case it prefers the median wall time from the
single most recent result run, then falls back to the most recent accepted
baseline median. The total also models one fresh M2 process per repetition and
two calibration processes per family. This is an estimate, not a benchmark;
it creates no result directory. The default process and calibration constants
may be overridden with `SYMRINGS_BENCH_PROCESS_OVERHEAD_SECONDS` and
`SYMRINGS_BENCH_CALIBRATION_SECONDS`.

The supported tiers are `Small`, `Medium`, `Large`, and `Stress`. The initial
catalog avoids putting correctness smoke tests into the timing data; package
tests remain the correctness authority. The runner nevertheless checks the
expected homogeneous result weight by default. Use `--no-verify` only for
diagnostic measurements.

Optional breadth cases live in `BasisConversion-extra`, `SchurProduct-extra`,
`HallLittlewoodProduct-extra`, and `InnerProduct-extra`. They cover structured
linear combinations, direct LR expansion, post-plethysm multiplication,
additional Hall–Littlewood product bases, and several inner-product dispatcher
routes. Run an extra family directly with `--family NAME`, or let a varied
profile sample it.

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

Each completed run creates one timestamped directory under `results/`:

```text
results/TIMESTAMP/
  raw.tsv
  summary.tsv
  comparison.tsv
  system.tsv
  conditions.tsv
  report.md
```

Pass `--output RUN-NAME` to choose the subdirectory name under `results/`
explicitly. Run names must be a single path component. The runner refuses to
reuse an existing directory, so separate invocations cannot accidentally mix
their raw repetitions or derived reports.

The baseline comparison is also printed in the terminal. Cases are classified
as `stable`, `improvement`, `regression`, or `new`. It reports changes from
both the best valid historical median and the most recent accepted median. The
The Markdown report starts with an overall summary, then groups baseline
comparisons into one table per benchmark family, followed by system information
and run conditions. Each family table is preceded by status counts and details
of improvements and regressions. New cases are counted but not enumerated in
prose. The overall and per-family summaries use a header row and one value row
for quick scanning.
The default classification threshold is ten percent; set
`THRESHOLD_PERCENT` to override it. Use `--baselines FILE` to compare against a
different baseline table.

The runner executes a fixed `conv-S-three-p` calibration probe in a fresh M2
process immediately before and after each selected family. Nothing monitors or
samples the machine concurrently with timed benchmark cases. The unprivileged
condition snapshots record power source, Low Power Mode, coarse macOS thermal
state when available, memory availability and pageouts, and load averages.
The report classifies run health as `clean`, `warning`, or `compromised` and
includes the meaning of all three classifications. Calibration drift of at
least 10% produces a warning and drift of at least 20% marks a run compromised;
the assessment also accounts for thermal state, Low Power Mode, memory
pressure, and pageouts. Pageout increases are informational below 64 MiB. They
produce a warning only at 64 MiB or more together with at least 1 MiB/s activity
or memory below 15%, and mark a run compromised only at 512 MiB or more together
with at least 5 MiB/s activity or memory below 5%.

Summarize a result file:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/summarize-results.awk \
  results/20260712-183000/raw.tsv > summary.tsv
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
  results/TIMESTAMP/summary.tsv >> \
  Macaulay2/packages/SymmetricRings/extras/benchmarks/baselines.tsv
```

Running benchmarks never silently changes their baseline.

Run artifacts are historical snapshots. After accepting a summary, do not
regenerate the originating run's comparison or report against the updated
baseline table. Its classifications should continue to describe what was known
when the run occurred.

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
