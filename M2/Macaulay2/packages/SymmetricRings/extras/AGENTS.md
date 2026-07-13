# SymmetricRings Extras: Agent Instructions

This directory contains development notes and benchmark tooling, not package
code loaded by `SymmetricRings.m2`. Keep `README.md` as the detailed runbook;
this file is the concise operational authority for work in this subtree.

## Build and working directory

Run commands from:

```sh
cd /Users/johngraf/M2Dev/Project-SymFcns/M2/M2
```

Install the package before measuring:

```sh
CCACHE_DISABLE=1 cmake --build BUILD/build \
  --target install-SymmetricRings -j2
```

After implementation changes, also run:

```sh
CCACHE_DISABLE=1 cmake --build BUILD/build \
  --target install-SymmetricRings check-SymmetricRings -j2
```

Use `BUILD/build/M2` unless `M2_BIN` deliberately selects another executable.

## Systematic benchmark suite

The suite is in `extras/benchmarks/`. Validate before cataloged runs:

```sh
BUILD/build/M2 --script \
  Macaulay2/packages/SymmetricRings/extras/benchmarks/validate.m2
```

Inspect and run selections with `run-benchmarks.sh`:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh --list
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --case pleth-three-row2-combined
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --family SchurPlethysm --tier Medium --repetitions 3
```

Supported timing tiers are `Small`, `Medium`, `Large`, and `Stress`.
Filters compose; use `--list` before a potentially long run.

Cross-family breadth profiles have fixed and random modes:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-fixed light
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-random standard --seed 42
```

Levels are `light`, `standard`, and `thorough`; they describe breadth, not
timing. Fixed profiles are nested and reproducible. Random profiles select at
most 1, 2, or 4 cases per eligible family; record or supply the seed.

Use `--new` for cases absent from both the selected baseline table and the
single most recently modified result run. Older non-baseline runs are ignored.

Use `--estimate` for a non-running wall-time estimate with any selection. It
prefers the most recent run, falls back to accepted baselines, models fresh M2
processes and family calibration, and creates no result directory.

Optional breadth cases belong in descriptive `Family-extra` families.

## Execution and result artifacts

Every timed repetition uses a fresh M2 process, a temporary `HOME`, and
sequential execution. Never parallelize cold timing processes.

`--output RUN-NAME` chooses a new single-component directory under
`benchmarks/results/`; never reuse an existing run directory. A run contains:

```text
raw.tsv  summary.tsv  comparison.tsv
system.tsv  conditions.tsv  report.md
```

Keep expected-weight verification enabled for production measurements;
`--no-verify` is diagnostic only.

Reports contain an overall summary, family comparison tables, a comparison of
CPU/RAM/OS with the selected accepted baselines, system data, and run
conditions. Put the baseline-system result beside run health in the overall
summary. Count new cases but do not enumerate them in prose. Never add serial
numbers, hardware UUIDs, or provisioning identifiers.

After a suite run, `report.md` is the only result file to share or link in
chat. Keep every TSV locally for reproducibility and diagnosis.

Condition collection must remain unprivileged. Run one fixed calibration probe
in a fresh process immediately before and after each family; run no monitor
concurrently with timed cases. Health is `clean`, `warning`, or `compromised`.
Review warnings and never accept a compromised run. Small pageout changes are
informational; use the magnitude, rate, and memory thresholds in README.

## Baseline policy

`benchmarks/baselines.tsv` contains reviewed medians. Acceptance is always
deliberate and never automatic. A baseline must:

- use automatic production routing, not forced or traced diagnostics;
- match the mathematical case and coefficient ring;
- use fresh sequential processes and at least three repetitions;
- come from a reviewed, non-compromised run;
- preserve every raw repetition;
- invalidate a correctness-broken baseline rather than deleting it.

Use median CPU time as primary. Compare against both the best valid and most
recent accepted medians. A fastest individual repetition is not a baseline.
For changes below about 10%, use at least five repetitions before classifying.

Accept a reviewed summary explicitly with `accept-summary.awk`; pass the run's
`system.tsv` before its `summary.tsv`, provide an acceptance date and useful
note, and append its output to `baselines.tsv`. This records the coarse system
configuration with each accepted timing.

Run directories are historical snapshots. After accepting a baseline, do not
regenerate that run's `comparison.tsv` or `report.md`; a case that was new when
measured must remain new in its originating report.

## Catalog maintenance

- Put reusable shapes and pair families in `benchmarks/partitions.m2`.
- Put declarative cases and varied profiles in `benchmarks/cases.m2`.
- Register operations once in `benchmarkKnownOperations` and implement them in
  `benchmarks/operations.m2`.
- Keep case IDs mathematical, not tied to the selected implementation route.
- Give combined and split stages the same `InputGroup`.
- Keep fixed breadth profiles nested and validate every referenced case ID.
- After edits, validate and run at least the affected family or tier.

## One-off benchmark rules

Use one-off timings for new expressions, selector boundaries, forced-route
comparisons, and component diagnosis before adding a catalog case.

Use a temporary `HOME`, no preload, and a fresh process per repetition:

```sh
env HOME=/private/tmp/symmetricrings-bench BUILD/build/M2 \
  --no-preload --silent --stop -q \
  -e 'needsPackage "SymmetricRings"; R=symmetricRing(frac(QQ[t])); time G=S_{5,3}@S_{2,1}; exit 0'
```

- Run cold repetitions sequentially.
- Keep setup and verification outside timing unless intentionally measured.
- Match rings, normalization, route settings, and timed scope exactly.
- Use at least three repetitions initially and five for sub-10% differences.
- Retain all CPU and wall times; compare medians, not fastest samples.
- Record the expression, ring, revision, routes, term counts, and min/median/max.
- Never mix traced, forced, exploratory, or invalid data with baselines.
- Take condition snapshots between groups, never monitor during short timings.

## Routing diagnostics

Tracing changes timing and remains diagnostic:

```sh
M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1
M2_SYMMETRIC_RINGS_TRACE_INNER_PRODUCT=1
```

Forced `p -> S` routes are:

```sh
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=abacus-rim-hooks
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=border-strips
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=via-complete
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=grouped-characters
```

Run each forced route in a separate cold process. Avoid grouped characters on
medium or large inputs without a timeout. Forced routes never become baselines.

Other diagnostic controls include:

```sh
M2_SYMMETRIC_RINGS_FORCE_HALL_LITTLEWOOD_PIPELINE=grouped
M2_SYMMETRIC_RINGS_FORCE_HALL_LITTLEWOOD_PIPELINE=fallback
M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_ROUTE=dual-basis-coefficient
```

## Isolating performance changes

- Compare `QQ` with `frac(QQ[t])` to isolate shadow-QQ lift and promotion.
- Split plethysm production from `toS`, and compare `f @ g` with the split path.
- Compare forced `p -> S` routes on the same guaranteed power-sum expression.
- Inspect input/output terms, weight, partition lengths, and coefficient shape.
- Test hooks, rectangles, repeated parts, balanced shapes, sparse, and dense data.
- With `debug needsPackage "SymmetricRings"`, time private QQ lift, computation,
  and promotion helpers separately; private timings never replace public ones.

When a one-off case becomes important, add it to the systematic catalog rather
than maintaining an independent timing history.
