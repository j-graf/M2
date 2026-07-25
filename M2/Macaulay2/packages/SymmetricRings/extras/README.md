# SymmetricRings extras: experiments and performance work

This directory contains development material for the `SymmetricRings`
package.  It is not part of the public package API.  Files here support
benchmarking, algorithm comparisons, profiling, and one-off mathematical
experiments.

The intended audience is package maintainers, including mathematicians who
want to test a proposed formula or compare algorithms before changing the
engine.

## Interactive conversion-plan comparisons

Load the package in development mode to use the private `toBasisBench`
function from an ordinary M2 session:

```m2
debug needsPackage "SymmetricRings";
R = symmetricRing QQ;
F = p_{8,4,2} + 2*p_{7,4,3} + p_{6,5,3};
report = toBasisBench(
    F, S,
    "Plans" => {
        "Automatic",
        "PowerSum->Schur:abacus-rim-hooks",
        "PowerSum->Schur:via-complete-basis"
        },
    "Repetitions" => 5,
    "Warmups" => 1,
    "Track" => true);
report#"Summary"
```

The function compares only built-in engine conversions. It verifies that the
requested plans agree, rotates their timed order, and returns all CPU and wall
measurements together with median summaries. Tracking is performed in separate
untimed executions. These are warm-session diagnostic measurements and never
update the systematic benchmark history. A forced-plan comparison requires a
product-free expansion in one source basis, so the named top-level plan cannot
accidentally affect a nested multiplication conversion.

Ordinary `toBasis` has no benchmark mode and ignores ambient conversion
forcing, conversion tracing, exhaustive-plan checking, and forced
multiplication settings. All conversion-plan controls above are explicit and
local to the `toBasisBench` call.

## What belongs here

Use `extras` for work such as:

- timing two mathematically equivalent algorithms;
- studying which partition statistics predict a faster route;
- checking an experimental conversion before integrating it;
- collecting profiling or dispatcher-trace information;
- reproducing a performance regression;
- systematic package benchmarks in the `benchmarks` subdirectory.

Correctness tests that should run with the package belong in the package test
suite, not only here.  Public documentation belongs in the package source.
Stable C++ algorithms belong in `e/symmetric-rings`.

## Two kinds of benchmark

There are two complementary workflows.

### The systematic suite

`benchmarks/` contains a table-driven suite with automatic fastest-median
records, machine and run-condition records, family summaries, calibrated
timing, fixed and random representative selections, and estimates of running
time.

See `benchmarks/README.md` for the commands and file format.

### One-off investigations

A small, focused script is often better while designing an algorithm or
investigating one expression.  Such a benchmark should still be disciplined:

- verify that every compared route gives the same result before timing it;
- warm up code paths whose first call includes setup or compilation;
- use several repetitions and report a median;
- time CPU work rather than including printing large expressions;
- keep construction of the test input outside the timed region when the
  construction is not under study;
- record the coefficient ring, bases, partitions, tags, and forced-route
  settings;
- include examples where the proposed method is expected to lose as well as
  win;
- avoid other substantial work on the machine while collecting results.

The retained cold-regression notebook is
`manual-regression-benchmarks.m2`. It should contain only manual comparisons
that are not represented by the systematic suite.

For very short operations, batch many identical or comparable operations so
that the total measurement is meaningful.  For long operations, fewer
repetitions may be appropriate, but always retain a correctness comparison.

## Choosing mathematical examples

Symmetric-function performance is rarely determined by weight alone.  A useful
study varies several of the following:

- weight and partition length;
- largest part and repeated parts;
- sparse versus dense linear combinations;
- homogeneous versus mixed-degree expressions;
- single basis elements, products, sums, and plethysms;
- inputs arising naturally from LR, Pieri, border-strip, or plethysm
  computations;
- coefficient rings, especially `QQ` versus rings that cannot use the
  constant-QQ shadow ring;
- cold versus already-constructed intermediate data.

For power-sum-to-Schur work, useful statistics include the number of terms,
support density within all partitions of the degree, largest cycle part,
common cycle parts, coefficient ring, and semantic tags.  Measure these
features directly rather than labeling an expression only as “large.”

Expressions that mathematicians actually form are especially valuable.  For
example:

```m2
S_lambda * S_mu
S_lambda * h_r
S_lambda * p_r
S_lambda @ S_mu
(S_lambda @ S_mu) * S_nu
toBasis(p, toBasis(S, S_lambda * S_mu))
```

These may exercise a different dispatcher region from a randomly generated
power-sum sum of the same degree.

## Comparing conversion algorithms

The engine provides forced conversion routes for diagnosis and benchmarking.
Use them to compare kernels on exactly the same expanded input.  Do not leave
a forced route enabled when measuring production behavior.

A sound conversion comparison has three layers:

1. establish exact equality among all applicable routes;
2. time the individual forced routes to understand the algorithmic tradeoff;
3. time automatic routing to decide whether the selector makes the right
   choice.

Trace output can confirm the selected pipeline and route.  Keep trace and
printing outside the timed interval.

The older border-strip and grouped-character `p -> S` routes are intentionally
retained for forced comparisons even though automatic dispatch normally
chooses abacus rim hooks or conversion via complete functions.  They are
valuable independent checks.

## Interpreting timing changes

Compare with the fastest qualifying median in the selected device's
`benchmarks/history/PROFILE/records.tsv`, not only the immediately preceding
run, and inspect whether the mathematical route, coefficient ring, or machine
conditions changed. The parallel `latest.tsv` is the most recent verified
median for each case on that device. This is especially important after
changing an implementation or selector.

An unchanged total time can conceal opposing changes in conversion and input
construction.  A regression may come from:

- selecting a different route;
- conversion to or from the constant-QQ shadow ring;
- recomputing metadata or partition statistics;
- loss of a factorized or tagged representation;
- a cold cache or one-time table construction;
- memory pressure, pageouts, or thermal behavior;
- timing a different portion of the computation.

Prefer explanations based on traces, profiles, and controlled comparisons.
Do not tune a general selector to a single named benchmark without checking the
larger mathematical pattern.

## From experiment to package change

When an experiment succeeds:

1. write down the mathematical formula and its preconditions;
2. identify the independent algorithm used to verify it;
3. determine whether it is a kernel, a conversion workflow, or only a new
   selector rule;
4. add exact package tests before relying on performance tests;
5. benchmark varied examples through both forced and automatic routing;
6. update the appropriate maintainer README and user documentation;
7. preserve the fallback and diagnostic routes unless there is strong reason
   to remove them.

The engine architecture and the detailed procedure for adding a `p -> S`
route are described in the
[engine maintainer guide](../../../e/symmetric-rings/README.md).

## Sharing results

For a systematic-suite run, share only its generated `report.md` in chat or a
review discussion.  The raw and intermediate files remain in the run folder so
that maintainers can inspect them locally when needed.

For a one-off benchmark, share enough information to reproduce the result:
the expression, ring, routes, repetitions, median times, correctness outcome,
and relevant system or run-condition notes.
