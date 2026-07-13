# Agent Notes for SymmetricRings Extras

This directory contains development notes and benchmarking tools. These files
are useful context for active development, but they are not package code loaded
by `SymmetricRings.m2`.

## Benchmarking Runbook

Run commands from the nested Macaulay2 source tree:

```sh
cd /Users/johngraf/M2Dev/Project-SymFcns/M2/M2
```

Build and install the current package before measuring:

```sh
CCACHE_DISABLE=1 cmake --build BUILD/build \
  --target install-SymmetricRings -j2
```

After implementation work, run the complete package verification as well:

```sh
CCACHE_DISABLE=1 cmake --build BUILD/build \
  --target install-SymmetricRings check-SymmetricRings -j2
```

The local executable is normally `BUILD/build/M2`.

## Systematic benchmark suite

The recorded benchmark system lives in `extras/benchmarks/`. This file is the
single authority for running it, maintaining its catalog, accepting baselines,
and interpreting or sharing results.

Validate case identifiers, partition references, operation names, coefficient
rings, timing tiers, varied profiles, and homogeneous weights without running
the computations, or list the selected case identifiers:

```sh
BUILD/build/M2 --script \
  Macaulay2/packages/SymmetricRings/extras/benchmarks/validate.m2

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh --list
```

Run one case, a mathematical family, or a timing tier:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --case pleth-three-row2-combined

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --family SchurPlethysm --repetitions 3

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --tier Medium --repetitions 3
```

Run a broad but bounded cross-family check with either a stable curated set or
a per-family random sample:

```sh
Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-fixed light

Macaulay2/packages/SymmetricRings/extras/benchmarks/run-benchmarks.sh \
  --varied-random standard --seed 42
```

The breadth levels are `light`, `standard`, and `thorough`. Fixed profiles are
nested and always select the same cases. Random profiles select at most 1, 2,
or 4 cases per eligible family; omit `--seed` to generate and record one, or
provide it to reproduce the selection. These levels describe coverage breadth,
not the existing timing tiers.

Use `--new` alone or with these filters to run cases found in neither the
selected baseline table nor the single most recently modified result run.
Older non-baseline runs are ignored; the report records which recent run was
used for the comparison.

Use `--estimate` with any selection filters to estimate total wall time without
running cases or creating a result directory. It prefers the most recent run's
medians, falls back to accepted baselines, and includes modeled fresh-process
and family-calibration overhead.

Filters may be combined. Use `--list` with the filters first when checking the
scope of a potentially long run. The supported tiers are `Small`, `Medium`,
`Large`, and `Stress`.

Optional breadth cases use `Family-extra` names, including additional basis
conversions, expanded and post-plethysm products, Hall–Littlewood basis
products, and inner-product dispatcher shapes. They may be run directly by
family or sampled through a varied profile.

Every completed timed invocation creates a new directory under
`extras/benchmarks/results/`. Use `--output RUN-NAME` to choose its
single-component name; an existing directory is never reused. The directory
contains:

```text
raw.tsv
summary.tsv
comparison.tsv
system.tsv
conditions.tsv
report.md
```

`raw.tsv` retains every cold repetition. `summary.tsv` reports per-case timing
statistics. `comparison.tsv` compares medians with the best valid and most
recent accepted baselines. The Markdown report begins with an overall summary,
then gives one summarized comparison table per benchmark family, followed by
system information and run conditions. Family tables give status counts and
details for improvements and regressions. Keep new-case counts in summary
tables, but do not enumerate new cases in prose.

After a benchmark suite is run, `report.md` is the only result file that should
be shared or linked in chat. Keep the TSV artifacts in the run directory for
reproducibility and local diagnosis, but do not share them separately in the
conversation. Never include hardware serial numbers, UUIDs, or provisioning
identifiers in system reports.

The runner uses only unprivileged condition checks. It runs a fixed calibration
probe in a fresh M2 process immediately before and after each selected family;
nothing monitors the machine concurrently with timed cases. Reports classify
run health as `clean`, `warning`, or `compromised`. Review warnings before
accepting a baseline, and do not accept a compromised run.

Small pageout deltas are informational and must not trigger a warning merely
because the cumulative counter increased. Use the magnitude, rate, and
concurrent memory thresholds documented in `benchmarks/README.md`.

The worker performs an expected-weight check after timing. Keep verification
enabled for production measurements. `--no-verify` is diagnostic only.

### Baseline acceptance

`benchmarks/baselines.tsv` contains reviewed accepted medians. Baseline
acceptance is deliberate and separate from running the suite; never allow a
new run to redefine the comparison standard automatically. A baseline must:

- come from automatic production routing rather than a forced diagnostic;
- use the same mathematical case and coefficient ring;
- use fresh sequential processes and at least three repetitions;
- come from a reviewed, non-compromised run;
- preserve every raw result file;
- mark a correctness-invalidated baseline invalid rather than deleting it.

Use median CPU time as the primary baseline statistic. Compare changes with
both the best valid historical median and the most recent accepted median. A
fastest individual repetition is not a baseline. For differences below roughly
10%, use at least five repetitions before classifying a regression or
improvement.

After reviewing a summary, deliberately append accepted rows with:

```sh
ACCEPTED_DATE=2026-07-12 NOTES='initial systematic baseline' \
  Macaulay2/packages/SymmetricRings/extras/benchmarks/accept-summary.awk \
  Macaulay2/packages/SymmetricRings/extras/benchmarks/results/TIMESTAMP/summary.tsv \
  >> Macaulay2/packages/SymmetricRings/extras/benchmarks/baselines.tsv
```

Treat every completed run directory as a historical snapshot of the baseline
state at execution time. Accepting its summary must not be followed by
regenerating that run's `comparison.tsv` or `report.md`; cases that were new
when measured remain `new` in the originating report.

### Maintaining the catalog

- Add reusable mathematical shapes and pair families to `benchmarks/partitions.m2`.
- Add declarative cases to `benchmarks/cases.m2`.
- Add each symbolic operation to `benchmarkKnownOperations` and implement it
  once in `benchmarks/operations.m2`.
- Keep case IDs mathematical; do not encode the currently selected algorithm.
- Connect combined operations and split stages with the same `InputGroup`.
- Put optional breadth cases in a descriptive `Family-extra` family and add
  representatives to fixed varied profiles rather than making every profile
  exhaustive.
- Keep fixed varied profiles nested and validate every referenced case ID.
- Re-run `benchmarks/validate.m2`, then execute at least the affected tier or
  family.

## General rules for one-off benchmarks

One-off measurements remain useful for new expressions, selector boundaries,
forced-route comparisons, and component-level diagnosis before a case belongs
in the systematic catalog.

Use a temporary `HOME` so startup files and user cache state do not affect the
measurement:

```sh
env HOME=/private/tmp/symmetricrings-bench \
  BUILD/build/M2 --no-preload --silent --stop -q \
  -e 'needsPackage "SymmetricRings"; R=symmetricRing(frac(QQ[t])); time G=S_{5,3}@S_{2,1}; exit 0'
```

For cold timings:

- Use a fresh M2 process for every repetition.
- Run repetitions sequentially; never parallelize timed processes.
- Keep construction, conversion, and verification outside the timed expression
  unless they are intentionally part of the operation being measured.
- Keep coefficient rings, normalization options, and timed scope identical
  when comparing results.
- Use at least three repetitions for an initial comparison. For differences
  below roughly 10%, use at least five cold repetitions before calling the
  change a regression or improvement.
- Prefer median CPU time for comparisons, but retain every individual CPU and
  wall time. Thread and GC time are useful when M2 rebuilding or
  coefficient-ring transfer is suspected.
- Record the exact expression, coefficient ring, route settings, package/Git
  revision, result term count, and minimum/median/maximum timings in the
  discussion.
- Do not mix exploratory, traced, forced-route, or invalid results with
  production baselines.

Do not run a continuous system monitor alongside short timed expressions. For
a longer one-off study, take unprivileged condition snapshots and calibration
measurements between groups of timings, not during them.

## Conversion and routing diagnostics

Print the selected conversion pipeline, route, target, basis guarantees, term
count, and weight with:

```sh
M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1
```

Tracing adds overhead. Never record a traced run as the benchmark timing.

Force a `p -> S` route for comparative diagnosis with one of:

```sh
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=abacus-rim-hooks
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=border-strips
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=via-complete
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=grouped-characters
```

Run forced routes in separate cold processes. `via-complete` and
`grouped-characters` are retained for forced diagnostics rather than normal
production selection. Grouped characters can become extremely slow outside
its intended tiny-weight domain; use a timeout or avoid forcing it on medium
and large examples.

The Hall-Littlewood generator pipeline can be forced with:

```sh
M2_SYMMETRIC_RINGS_FORCE_HALL_LITTLEWOOD_PIPELINE=grouped
M2_SYMMETRIC_RINGS_FORCE_HALL_LITTLEWOOD_PIPELINE=fallback
```

The fallback setting is for correctness and algorithm comparisons, not
standard production timings.

Inner-product selection can be inspected or forced with:

```sh
M2_SYMMETRIC_RINGS_TRACE_INNER_PRODUCT=1
M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_ROUTE=dual-basis-coefficient
```

As with conversion tracing and forcing, keep these results diagnostic unless
the automatic production route is being timed without tracing.

## Isolating performance changes

Useful ways to localize a slowdown include:

- Compare `R = symmetricRing(QQ)` with `R = symmetricRing(frac(QQ[t]))`.
  Direct QQ bypasses the constant-QQ shadow lift and result promotion.
- Split `F = plethysm(f,g)` from `toS F` to distinguish plethysm production
  from the subsequent conversion.
- Compare `f @ g` with `toS plethysm(f,g)` to distinguish the combined entry
  from a later conversion. Account for a possible QQ-shadow lift of a large
  power-sum expansion.
- Convert a Schur product to power sums first, then compare forced `p -> S`
  routes on the same guaranteed power-sum expression.
- Inspect input and output term counts. Promotion and rebuilding costs often
  scale with the number of terms rather than weight alone.
- Compare equal weights with different partition lengths, repeated parts,
  hooks, rectangles, balanced shapes, and sparse versus dense expansions.

Package-private helpers can be inspected in an exploratory M2 process with:

```m2
debug needsPackage "SymmetricRings"
```

For example, `constantQQLiftElement`, `plethysmToBasisFast`, and
`constantQQPromoteElement` can be timed separately to distinguish QQ
computation from term-by-term promotion back to the original ring. These are
diagnostic helpers, not public API, and component timings should not replace
the public operation benchmark.

When a one-off expression becomes important for regression tracking, add a
mathematical case to the systematic catalog rather than maintaining an
independent timing history. Keep case identifiers mathematical and independent
of whichever implementation route the dispatcher currently selects.
