# SymmetricRings systematic benchmark suite

This directory contains the reproducible benchmark system for the
`SymmetricRings` package. It measures representative mathematical operations,
compares a run with the fastest qualifying records, records the machine and run
conditions, and writes a human-readable report.

The suite is designed so that a contributor can run it without elevated
privileges and without knowing the implementation details of every algorithm.

## Quick start

Run commands from this directory.  The runner's help and listing modes are the
authoritative source for the exact command-line spelling:

```sh
./run-benchmarks.sh --help
./run-benchmarks.sh --list
```

A normal selected run has the form:

```sh
./run-benchmarks.sh --family SchurPlethysm
```

Each run creates one timestamped directory under `results/`.  Its final
`report.md` contains the comparison intended for people to read.  When sharing
a run in chat, share **only `report.md`**, not all of the auxiliary files.

Before a long run, use estimate mode with the same selection options. It reads
the latest available timings and fastest records but does not execute tests
or create a results directory:

```sh
./run-benchmarks.sh --family SchurPlethysm --estimate
```

## How cases are organized

The suite is table-driven.  Partition data, case definitions, and operation
definitions are kept separately so that mathematical inputs can be reused.
Each case has a stable name, a family, and a tier.

Core families cover routinely monitored operations.  Additional breadth is
placed in families whose names end in `-extra`, including extra basis
conversions, Schur and Hall--Littlewood products, and inner products.  The
extra families are useful for audits and algorithm work without making every
ordinary run exhaustive.

Use listing mode to see the current cases and their metadata.  Do not copy a
case list into this README: the table files are the maintained source of truth.

## Selecting work

Selections can be made by case, family, or tier.  Multiple selectors can be
used to focus a run.  Representative cross-family checks are available at
three coverage levels:

```sh
./run-benchmarks.sh --varied-fixed light
./run-benchmarks.sh --varied-fixed standard
./run-benchmarks.sh --varied-fixed thorough
```

The fixed form always selects the same representative cases, making it useful
for repeated checks during development.

The random form samples within the selected families:

```sh
./run-benchmarks.sh --varied-random standard --seed 12345
```

Record the seed when sharing a random run.  The coverage levels select
increasing numbers of cases per family; they describe breadth, not a promise
about execution speed.

To execute only cases that have no fastest record and no result in the single
most recently modified run, use:

```sh
./run-benchmarks.sh --new
```

“New” deliberately checks only `records.tsv` and that one latest run. It does
not scan all historical result directories.

## What happens during a run

For each selected case and repetition, the suite launches a fresh Macaulay2
process.  This avoids accidental sharing of caches between unrelated cases and
makes individual failures easier to diagnose.

The suite records:

- raw per-repetition measurements;
- summarized median timings;
- comparison with the pre-run fastest records;
- system information;
- calibration probes before and after each family;
- a run-condition classification;
- a Markdown report.

After writing the report, the suite automatically updates `records.tsv` with
any lower median based on at least three repetitions.

No performance monitor runs during a timed test.  Calibration occurs only
before and after families, so the probes do not compete with benchmark work.
All system and memory information is collected without elevated privileges.

The time estimate uses the most recent run median when one is available,
otherwise the fastest record, plus the suite's startup and calibration
model.  It is an estimate rather than a scheduling guarantee.

## Results directories

A completed run creates a directory like:

```text
results/
  YYYYMMDD-HHMMSS/
    raw...
    summary...
    comparison...
    system...
    conditions...
    report.md
```

When `--output LABEL` is supplied, the directory is named
`YYYYMMDD-HHMMSS-LABEL`; the timestamp is never omitted.

Exact auxiliary filenames may evolve.  Treat `report.md` as the stable human
entry point and retain the entire directory locally as the audit record.

The report is ordered for review:

1. overall summary;
2. fastest-record comparison, broken into one table per family;
3. system information;
4. run conditions.

The overall and family summaries use compact two-row tables. They count
improvements, regressions, stable cases, and new cases relative to the fastest
record; they also count records set. Detailed family tables contain per-case
classifications. The report does not repeat improvements or regressions in
prose or list new cases in prose.

## Repetitions and timing

The reported comparison uses median CPU time. A record-qualifying run
uses at least three repetitions.  Very short or noisy comparisons should use
five or more.  Keep the number of repetitions the same when comparing closely
matched runs.

The benchmarked operation should contain the mathematical work named by the
case, not report formatting or printing.  A test must also verify enough of its
result to catch an incorrect shortcut; fast wrong answers are not performance
improvements.

Each repetition uses the same mathematical input and coefficient ring recorded
by the case.  Changing a partition, ring, or construction method creates a
different benchmark even if the displayed operation looks similar.

## Run conditions

Every report classifies the run:

| Classification | Meaning |
|---|---|
| **Clean** | Calibration and memory indicators show no material disturbance. |
| **Warning** | Some condition may have influenced timings; interpret small differences cautiously. |
| **Compromised** | Strong evidence of interference or instability; rerun before accepting conclusions. |

The calibration probes compare stable work before and after each family.  A
change of roughly 10 percent raises a warning; roughly 20 percent compromises
the run.

Pageouts are evaluated by transferred volume and rate rather than by a tiny
nonzero count.  Less than 64 MiB is informational.  Larger activity becomes a
warning or compromised condition only when paired with substantial transfer
rate or low available memory.  This avoids labeling harmless background
pageouts as a bad run.

These classifications are evidence, not automatic explanations.  A clean run
can still contain a real software regression, while a warning may affect only
one family.  Read the family calibration and system notes alongside the
timings.

For cleaner measurements:

- connect a laptop to power;
- close CPU- and memory-heavy applications;
- allow the machine to reach a stable temperature;
- avoid builds, indexing, backups, and large downloads;
- rerun important regressions independently.

## Fastest records

`records.tsv` contains one fastest qualifying median per case and coefficient
ring. A qualifying value is the median CPU time of at least three repetitions;
an unusually fast individual repetition never becomes a record.

For each run, `comparison.tsv` and `report.md` use the record table as it stood
before that run. This preserves a meaningful comparison when the current run
sets a new record. After the report is complete, the runner automatically
updates `records.tsv`, recording the median, repetition count, capture time,
and source run directory.

Records are the active development comparison. They measure distance from the
fastest qualifying run observed for the same case and coefficient ring.

The record classification uses the configured percentage threshold.
`record_status` is separate and exact: it says whether an eligible run sets,
ties, or misses the record. Thus a small record-setting change may correctly
have `record_status = new-record` while its thresholded record classification
is `stable`.

## Adding a benchmark case

To add a case:

1. Choose or add named partition data.
2. Choose an existing mathematical operation, or add an operation definition
   with a precise result check.
3. Give the case a stable descriptive name.
4. Assign the most specific family and an appropriate tier.
5. Use an `-extra` family when the case adds breadth but need not run in every
   routine check.
6. Validate the tables and list the case before executing it.
7. Run it with enough repetitions and inspect its report.
8. Run at least three repetitions so the case can establish a record.

Cases should collectively vary shapes, weights, support densities, coefficient
rings, and semantic constructions.  Prefer mathematically recognizable inputs
such as Schur products, Pieri products, plethysms, and round-trip basis
conversions alongside synthetic stress cases.

Do not make a benchmark depend on the current dispatcher's choice unless the
case is explicitly a selector test.  Most cases should state a mathematical
operation and allow production routing to evolve.

## Algorithm comparisons and dispatcher work

Forced routes and selector tracing are diagnostic tools.  They are useful for
answering two different questions:

- Which mathematical kernel is fastest for this fixed expanded input?
- Does automatic dispatch choose well for realistic inputs?

Keep these questions separate. First compare exact outputs of forced routes,
then time them, then run the ordinary unforced case. General catalog records
should normally reflect the ordinary selector.

For a selector change, include examples on both sides of every proposed
threshold and examples from several semantic sources.  Weight alone is often
insufficient: term count, support density, largest part, common parts,
coefficient ring, and combinatorial tags can all predict behavior.

## Files in this directory

The suite is divided into small components for maintainability:

- partition tables define reusable mathematical indices;
- case tables combine inputs, operations, families, and tiers;
- operation code constructs and verifies computations;
- the runner and shell wrapper select and execute cases;
- validation checks the tables before expensive work;
- summarization and comparison compute medians and record differences;
- `update-records.awk` updates fastest qualifying medians after reporting;
- estimation predicts duration without running cases;
- system, condition, and calibration code describe the environment;
- report code renders the final Markdown document;
- `records.tsv` stores automatically maintained fastest qualifying medians.

When changing one component, preserve the distinction between raw measurement,
statistical summary, comparison policy, and presentation.  In particular, do
not encode report wording in a mathematical operation or archive policy in a
case definition.

## Troubleshooting

If a case fails, run that single case with one repetition and inspect its raw
output.  Confirm that Macaulay2 loads the intended working-tree package rather
than another installed copy.

If timings unexpectedly regress:

1. check the run-condition classification;
2. compare family calibration probes;
3. confirm the coefficient ring and exact case definition;
4. use dispatcher tracing outside the timed region;
5. compare a forced old and new route on the same expanded input;
6. repeat the case in a fresh run;
7. inspect whether input construction, shadow-QQ conversion, or the target
   operation is actually responsible.

If a new case is not classified as new, remember that `--new` consults both
`records.tsv` and the single most recent raw run. A result in either is enough
to make the case non-new.

## Reporting a benchmark result

For routine collaboration, the generated `report.md` is the only run artifact
that should be shared in chat.  It already contains the overall and family
summaries, detailed comparisons, system information, and run conditions.

Keep the remaining files in the timestamped results directory.  They provide
the local evidence needed for deeper investigation without overwhelming the
initial review.
