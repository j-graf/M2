-- One-process worker for the systematic SymmetricRings benchmark catalog.
-- run-benchmarks.sh starts this file in a fresh M2 process for every timed
-- repetition. It can also list filtered case identifiers.

benchmarkDirectory = currentFileDirectory
-- The catalog includes development-only strict-kernel operations.  Debug
-- loading makes those private benchmark boundaries visible without exporting
-- them from the user package.
debug needsPackage "SymmetricRings"
load(benchmarkDirectory | "benchmark-helpers.m2")
load(benchmarkDirectory | "partitions.m2")
load(benchmarkDirectory | "operations.m2")
load(benchmarkDirectory | "cases.m2")

benchmarkEnvironment = key -> (
    value := try getenv key else null;
    if value === "" then null else value
    )

benchmarkMode = benchmarkEnvironment "SYMRINGS_BENCH_MODE"
if benchmarkMode === null then benchmarkMode = "run"
benchmarkFamilyFilter = benchmarkEnvironment "SYMRINGS_BENCH_FAMILY"
benchmarkTierFilter = benchmarkEnvironment "SYMRINGS_BENCH_TIER"
benchmarkIdFilter = benchmarkEnvironment "SYMRINGS_BENCH_CASE"
benchmarkVariedFixedFilter = benchmarkEnvironment "SYMRINGS_BENCH_VARIED_FIXED"
benchmarkTraceWorkflows =
    benchmarkEnvironment "SYMRINGS_BENCH_TRACE_WORKFLOW" === "1"

if benchmarkMode == "list" then (
    scan(benchmarkSelectCases(benchmarkFamilyFilter,
                              benchmarkTierFilter,
                              benchmarkIdFilter,
                              benchmarkVariedFixedFilter),
         case -> print case#"ID");
    exit 0
    )

if benchmarkMode == "plan" then (
    scan(benchmarkSelectCases(benchmarkFamilyFilter,
                              benchmarkTierFilter,
                              benchmarkIdFilter,
                              benchmarkVariedFixedFilter),
         case -> print(case#"ID" | "|" | case#"Family" | "|" |
                       case#"CoefficientRing" | "|" | case#"Tier"));
    exit 0
    )

-- This fixed, moderately sized conversion is a machine-performance probe,
-- not a catalog result. The shell runner invokes it in fresh processes only
-- before and after complete benchmark families.
if benchmarkMode == "calibration" then (
    calibrationCase := benchmarkCase "conv-S-three-p";
    Rbenchmark = benchmarkMakeRing(calibrationCase#"CoefficientRing");
    calibrationCPUStart := cpuTime();
    calibrationTiming := elapsedTiming (benchmarkExecute calibrationCase);
    print("CALIBRATION|conv-S-three-p|" |
          toString(cpuTime() - calibrationCPUStart) | "|" |
          toString(calibrationTiming#0));
    exit 0
    )

if benchmarkIdFilter === null then
    error "set SYMRINGS_BENCH_CASE to one benchmark identifier"

case = benchmarkCase benchmarkIdFilter
Rbenchmark = benchmarkMakeRing(case#"CoefficientRing")

cpuStart = cpuTime()
wallTiming = elapsedTiming (benchmarkExecute case)
cpuSeconds = cpuTime() - cpuStart
wallSeconds = wallTiming#0
benchmarkResult = wallTiming#1

benchmarkVerify = benchmarkEnvironment "SYMRINGS_BENCH_VERIFY"
if benchmarkVerify === "1" then
    benchmarkVerifyResult(case, benchmarkResult)

resultTerms = if instance(benchmarkResult, SymmetricRingElement)
    then #rawTerms benchmarkResult else 0
resultWeight = if instance(benchmarkResult, SymmetricRingElement)
    then weight benchmarkResult else 0
repetition = benchmarkEnvironment "SYMRINGS_BENCH_REPETITION"
if repetition === null then repetition = "1"
lambdaValue = if case#?"Lambda" then benchmarkInputPartition(case, "Lambda") else {}
muValue = if case#?"Mu" then benchmarkInputPartition(case, "Mu") else {}
probeValue = if case#?"Probe" then case#"Probe" else {}
pairClass = if case#?"PairClass" then case#"PairClass" else ""

print("BENCH|" |
      case#"ID" | "|" |
      case#"Family" | "|" |
      case#"Operation" | "|" |
      case#"Tier" | "|" |
      case#"CoefficientRing" | "|" |
      repetition | "|" |
      toString cpuSeconds | "|" |
      toString wallSeconds | "|" |
      toString resultTerms | "|" |
      toString resultWeight | "|" |
      case#"InputGroup" | "|" |
      toString lambdaValue | "|" |
      toString muValue | "|" |
      toString probeValue | "|" |
      toString(sum lambdaValue) | "|" |
      toString(#lambdaValue) | "|" |
      toString(sum muValue) | "|" |
      toString(#muValue) | "|" |
      toString(benchmarkExpectedWeight case) | "|" |
      pairClass)

exit 0
