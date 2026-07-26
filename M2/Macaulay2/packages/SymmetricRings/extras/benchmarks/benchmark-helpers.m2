-- Development-only comparison helpers for SymmetricRings engine routes.
-- Load this file explicitly after:
--     debug needsPackage "SymmetricRings"
-- It is intentionally not loaded by the package.

-- ============================================================================
-- Basis-Conversion Benchmarks
-- ============================================================================

toBasisBenchOptionDefaults = hashTable {
    "Plans" => {"Automatic"},
    "Repetitions" => 3,
    "Warmups" => 1,
    "Track" => false,
    "Verify" => true
    }

-- Gives Automatic a stable report label and represents it by the empty raw
-- identifier expected by the development engine entry point.
toBasisBenchPlanLabel = plan -> (
    if class plan === String then plan
    else error "expected each toBasisBench plan to be a string"
    )

toBasisBenchPlanIdentifier = plan -> (
    if plan == "Automatic" then "" else toBasisBenchPlanLabel plan
    )

-- Computes the median of a nonempty list of real timing measurements.
toBasisBenchMedian = measurements -> (
    ordered := sort measurements;
    middle := (#ordered) // 2;
    if odd(#ordered) then ordered#middle
    else (ordered#(middle - 1) + ordered#middle) / 2
    )

-- Executes one automatic or forced engine conversion. This mirrors the
-- production working-ring decision, including the constant-QQ shadow, but
-- calls the separate development raw entry point rather than public toBasis.
toBasisBenchExecute = (f, B, plan, track) -> (
    R0 := ring f;
    dispatchData := conversionDispatchData(f, B);
    if not isNativeBasisConversionApplicable dispatchData then
        error "toBasisBench supports only built-in engine conversions";
    planIdentifier := toBasisBenchPlanIdentifier plan;
    executeInRing := (Fwork, Bwork) ->
        userSymmetricElement(
            ring Fwork,
            rawSymmetricRingsToBasisBench(
                raw Fwork,
                Bwork#"BasisId",
                planIdentifier,
                track));

    if dispatchData#"PreferConstantQQ" then (
        workingInputs := toConstantQQIfPossible({f}, true);
        if isConstantQQBasisConversionApplicable(R0, workingInputs, B) then (
            Rwork := ring workingInputs#0;
            Bwork := basis(Rwork, basisKey B);
            if track then
                stderr << "SymmetricRings toBasisBench: working-ring=QQ-shadow"
                    << endl;
            return returnFromConstantQQ(
                executeInRing(workingInputs#0, Bwork),
                R0);
            );
        );
    if track then
        stderr << "SymmetricRings toBasisBench: working-ring=original"
            << endl;
    executeInRing(f, B)
    )

-- A forced top-level plan is unambiguous only for a product-free expansion in
-- one source basis. Automatic-only timings may still exercise the full mixed
-- and product-aware workflow.
validateToBasisBenchForcedInput = f -> (
    sourceBasisIds := {};
    scan(rawTerms f, term -> (
            atoms := term#1;
            if #atoms > 1 then
                error "forced toBasisBench plans require product-free input";
            if #atoms == 1 then
                sourceBasisIds = append(
                    sourceBasisIds,
                    (atoms#0)#"BasisId");
            ));
    sourceBasisIds = unique sourceBasisIds;
    if #sourceBasisIds != 1 then
        error "forced toBasisBench plans require one source basis";
    )

-- Compares complete plans on one fixed mathematical input. Warmups and trace
-- runs are untimed; timed plan order rotates between repetitions. The result
-- is an in-session diagnostic record and never updates systematic benchmarks.
toBasisBench = args -> (
    L := argumentList args;
    if #L < 2 then
        error "expected a symmetric function and a target basis";
    f := L#0;
    target := L#1;
    if not instance(f, SymmetricRingElement) then
        error "expected a symmetric function";
    opts := parseStringOptions(
        toBasisBenchOptionDefaults,
        drop(L, 2),
        "toBasisBench");
    plans := opts#"Plans";
    if class plans === String then plans = {plans};
    if not instance(plans, List) or #plans == 0 then
        error "expected Plans to be a nonempty list";
    planLabels := apply(plans, toBasisBenchPlanLabel);
    if #(unique planLabels) != #planLabels then
        error "expected distinct toBasisBench plans";
    if any(planLabels, label -> label != "Automatic") then
        validateToBasisBenchForcedInput f;

    repetitions := opts#"Repetitions";
    warmups := opts#"Warmups";
    track := opts#"Track";
    verify := opts#"Verify";
    if class repetitions =!= ZZ or repetitions <= 0 then
        error "expected Repetitions to be a positive integer";
    if class warmups =!= ZZ or warmups < 0 then
        error "expected Warmups to be a nonnegative integer";
    if class track =!= Boolean then
        error "expected Track to be Boolean";
    if class verify =!= Boolean then
        error "expected Verify to be Boolean";

    R0 := ring f;
    B := targetBasisOnRing(R0, target);
    if not isNativeBasisConversionApplicable(
        conversionDispatchData(f, B)) then
        error "toBasisBench supports only built-in engine conversions";

    -- Warm every requested plan equally before collecting measurements.
    for warmup from 1 to warmups do
        scan(plans, plan -> toBasisBenchExecute(f, B, plan, false));

    timings := new MutableHashTable;
    results := new MutableHashTable;
    scan(planLabels, label -> timings#label = {});
    for repetition from 0 to repetitions - 1 do
        for offset from 0 to #plans - 1 do (
            position := (repetition + offset) % #plans;
            plan := plans#position;
            label := planLabels#position;
            cpuStart := cpuTime();
            timed := elapsedTiming (
                result := toBasisBenchExecute(f, B, plan, false));
            cpuSeconds := cpuTime() - cpuStart;
            wallSeconds := timed#0;
            results#label = result;
            timings#label = append(
                timings#label,
                hashTable {
                    "Repetition" => repetition + 1,
                    "CPUSeconds" => cpuSeconds,
                    "WallSeconds" => wallSeconds
                    });
            );

    verified := null;
    if verify then (
        reference := results#(planLabels#0);
        scan(drop(planLabels, 1), label ->
            if results#label != reference then
                error("toBasisBench plans disagree: ", planLabels#0,
                    " and ", label));
        verified = true;
        );

    -- Tracking is intentionally separated from every timed repetition.
    if track then scan(plans, plan -> (
            stderr << "SymmetricRings toBasisBench: plan="
                << toBasisBenchPlanLabel plan << endl;
            toBasisBenchExecute(f, B, plan, true);
            ));

    summaries := hashTable apply(planLabels, label -> (
            planTimings := timings#label;
            cpuMeasurements := apply(
                planTimings, item -> item#"CPUSeconds");
            wallMeasurements := apply(
                planTimings, item -> item#"WallSeconds");
            orderedCPU := sort cpuMeasurements;
            orderedWall := sort wallMeasurements;
            label => hashTable {
                "MinimumCPUSeconds" => first orderedCPU,
                "MedianCPUSeconds" =>
                    toBasisBenchMedian cpuMeasurements,
                "MaximumCPUSeconds" => last orderedCPU,
                "MinimumWallSeconds" => first orderedWall,
                "MedianWallSeconds" =>
                    toBasisBenchMedian wallMeasurements,
                "MaximumWallSeconds" => last orderedWall
                }
            ));
    immutableTimings := hashTable apply(
        planLabels, label -> label => timings#label);
    immutableResults := hashTable apply(
        planLabels, label -> label => results#label);
    hashTable {
        "Input" => f,
        "TargetBasis" => B,
        "Plans" => planLabels,
        "Repetitions" => repetitions,
        "Warmups" => warmups,
        "Tracked" => track,
        "Verified" => verified,
        "Result" =>
            if #planLabels == 1
            then results#(planLabels#0)
            else null,
        "Results" => immutableResults,
        "Timings" => immutableTimings,
        "Summary" => summaries
        }
    )

-- ============================================================================
-- Strict-Binary Multiplication Benchmarks
-- ============================================================================

multiplyToBasisBenchOptionDefaults = hashTable {
    "Kernels" => {"Automatic"},
    "Repetitions" => 3,
    "Warmups" => 1,
    "Track" => false,
    "Verify" => true
    }

multiplyToBasisBenchExecute = (f, g, B, kernel, track) -> (
    forcedKernel :=
        if kernel == "Automatic" or kernel == "PowerSumReference"
        then ""
        else kernel;
    usePowerSumReference := kernel == "PowerSumReference";
    userSymmetricElement(
        ring f,
        rawSymmetricRingsMultiplyToBasisBench(
            raw f,
            raw g,
            B#"BasisId",
            forcedKernel,
            usePowerSumReference,
            track))
    )

-- Executes the product-free bilinear workflow with a request-scoped strategy
-- and optional tracing. Four arguments retain the automatic strategy.
multiplyExpressionsToBasisBenchExecute = args -> (
    L := argumentList args;
    if #L != 4 and #L != 5 then
        error "expected two symmetric functions, a target basis, an optional strategy, and a trace flag";
    f := L#0;
    g := L#1;
    B := L#2;
    strategy := if #L == 4 then "Automatic" else L#3;
    track := if #L == 4 then L#3 else L#4;
    userSymmetricElement(
        ring f,
        rawSymmetricRingsMultiplyExpressionsToBasisBench(
            raw f,
            raw g,
            B#"BasisId",
            strategy,
            track))
    )

multiplyToBasisBench = args -> (
    L := argumentList args;
    if #L < 3 then
        error "expected two symmetric functions and a target basis";
    f := L#0;
    g := L#1;
    target := L#2;
    if not instance(f, SymmetricRingElement)
        or not instance(g, SymmetricRingElement) then
        error "expected two symmetric functions";
    if ring f =!= ring g then
        error "expected elements in the same symmetric ring";
    opts := parseStringOptions(
        multiplyToBasisBenchOptionDefaults,
        drop(L, 3),
        "multiplyToBasisBench");
    kernels := opts#"Kernels";
    if class kernels === String then kernels = {kernels};
    if not instance(kernels, List) or #kernels == 0
        or any(kernels, kernel -> class kernel =!= String) then
        error "expected Kernels to be a nonempty list of strings";
    if #(unique kernels) != #kernels then
        error "expected distinct multiplyToBasisBench kernels";

    repetitions := opts#"Repetitions";
    warmups := opts#"Warmups";
    track := opts#"Track";
    verify := opts#"Verify";
    if class repetitions =!= ZZ or repetitions <= 0 then
        error "expected Repetitions to be a positive integer";
    if class warmups =!= ZZ or warmups < 0 then
        error "expected Warmups to be a nonnegative integer";
    if class track =!= Boolean then
        error "expected Track to be Boolean";
    if class verify =!= Boolean then
        error "expected Verify to be Boolean";

    R := ring f;
    B := targetBasisOnRing(R, target);
    leftDispatchData := conversionDispatchData(f, B);
    rightDispatchData := conversionDispatchData(g, B);
    if not leftDispatchData#"SourceUsesOnlyEngineBases"
        or not rightDispatchData#"SourceUsesOnlyEngineBases"
        or not leftDispatchData#"TargetIsEngineReadable"
        or not rightDispatchData#"TargetIsEngineReadable" then
        error "multiplyToBasisBench supports only built-in engine bases";

    for warmup from 1 to warmups do
        scan(kernels, kernel ->
            multiplyToBasisBenchExecute(f, g, B, kernel, false));

    timings := new MutableHashTable;
    results := new MutableHashTable;
    scan(kernels, kernel -> timings#kernel = {});
    for repetition from 0 to repetitions - 1 do
        for offset from 0 to #kernels - 1 do (
            position := (repetition + offset) % #kernels;
            kernel := kernels#position;
            cpuStart := cpuTime();
            timed := elapsedTiming (
                result :=
                    multiplyToBasisBenchExecute(
                        f, g, B, kernel, false));
            results#kernel = result;
            timings#kernel = append(
                timings#kernel,
                hashTable {
                    "Repetition" => repetition + 1,
                    "CPUSeconds" => cpuTime() - cpuStart,
                    "WallSeconds" => timed#0
                    });
            );

    verified := null;
    if verify then (
        reference := results#(kernels#0);
        scan(drop(kernels, 1), kernel ->
            if results#kernel != reference then
                error("multiplyToBasisBench kernels disagree: ",
                    kernels#0, " and ", kernel));
        verified = true;
        );

    if track then scan(kernels, kernel -> (
            stderr << "SymmetricRings multiplyToBasisBench: kernel="
                << kernel << endl;
            multiplyToBasisBenchExecute(f, g, B, kernel, true);
            ));

    summaries := hashTable apply(kernels, kernel -> (
            kernelTimings := timings#kernel;
            cpuMeasurements := apply(
                kernelTimings, item -> item#"CPUSeconds");
            wallMeasurements := apply(
                kernelTimings, item -> item#"WallSeconds");
            kernel => hashTable {
                "MinimumCPUSeconds" => first sort cpuMeasurements,
                "MedianCPUSeconds" =>
                    toBasisBenchMedian cpuMeasurements,
                "MaximumCPUSeconds" => last sort cpuMeasurements,
                "MinimumWallSeconds" => first sort wallMeasurements,
                "MedianWallSeconds" =>
                    toBasisBenchMedian wallMeasurements,
                "MaximumWallSeconds" => last sort wallMeasurements
                }
            ));
    hashTable {
        "LeftInput" => f,
        "RightInput" => g,
        "TargetBasis" => B,
        "Kernels" => kernels,
        "Repetitions" => repetitions,
        "Warmups" => warmups,
        "Tracked" => track,
        "Verified" => verified,
        "Result" =>
            if #kernels == 1 then results#(kernels#0) else null,
        "Results" =>
            hashTable apply(kernels, kernel -> kernel => results#kernel),
        "Timings" =>
            hashTable apply(kernels, kernel -> kernel => timings#kernel),
        "Summary" => summaries
        }
    )
