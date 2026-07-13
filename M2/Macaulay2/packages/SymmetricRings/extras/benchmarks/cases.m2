-- Declarative benchmark case catalog. Load partitions.m2 first.

benchmarkCaseRecord = (id, family, operation, ringKey, tier, fields) ->
    hashTable(({
        ("ID", id),
        ("Family", family),
        ("Operation", operation),
        ("CoefficientRing", ringKey),
        ("Tier", tier)
        }) | fields)

benchmarkUnarySpecs = {
    {"h16-p", "Complete", {16}, "QQ", "Small"},
    {"h24-p", "Complete", {24}, "FracQQt", "Small"},
    {"h30-p", "Complete", "w30-row", "FracQQt", "Medium"},
    {"e16-p", "Elementary", {16}, "QQ", "Small"},
    {"e24-p", "Elementary", {24}, "FracQQt", "Small"},
    {"e30-p", "Elementary", "w30-row", "FracQQt", "Medium"},
    {"S-hook-p", "Schur", "w12-hook", "FracQQt", "Small"},
    {"S-rectangle-p", "Schur", "w12-rectangle", "QQ", "Small"},
    {"S-three-p", "Schur", "w28-three-generic", "FracQQt", "Medium"},
    {"Somega-three-p", "SchurOmega", {11,8,6}, "FracQQt", "Medium"},
    {"m-generic-p", "Monomial", {6,3,1}, "FracQQt", "Small"},
    {"ff-generic-p", "Forgotten", {6,3,1}, "FracQQt", "Small"},
    {"q-row-p", "HallLittlewoodQGenerator", {21}, "FracQQt", "Medium"},
    {"b-row-p", "HallLittlewoodBGenerator", {21}, "FracQQt", "Medium"},
    {"Q-three-p", "HallLittlewoodQ", {14,4,2}, "FracQQt", "Medium"},
    {"B-three-p", "HallLittlewoodB", {14,4,2}, "FracQQt", "Medium"},
    {"P-three-p", "HallLittlewoodP", {14,4,2}, "FracQQt", "Medium"},
    {"Pomega-three-p", "HallLittlewoodPOmega", {14,4,2}, "FracQQt", "Medium"}
    }

benchmarkUnaryCases = apply(benchmarkUnarySpecs, spec ->
    benchmarkCaseRecord(
        "conv-"|spec#0,
        "BasisConversion",
        "BasisToPowerSums",
        spec#3,
        spec#4,
        {("SourceBasis", spec#1), ("Lambda", spec#2), ("InputGroup", spec#0)}))

benchmarkPowerSumTargetSpecs = {
    {"p30-h", "Complete", {30}, "FracQQt", "Medium"},
    {"p30-e", "Elementary", {30}, "FracQQt", "Medium"},
    {"p-four-S", "Schur", "w24-four-generic", "FracQQt", "Small"},
    {"p-four-Somega", "SchurOmega", "w24-four-generic", "FracQQt", "Small"},
    {"p-ones-m", "Monomial", {1,1,1,1,1,1,1,1,1,1}, "FracQQt", "Small"},
    {"p-ones-ff", "Forgotten", {1,1,1,1,1,1,1,1,1,1}, "FracQQt", "Small"},
    {"p18-q", "HallLittlewoodQGenerator", {18}, "FracQQt", "Small"},
    {"p18-b", "HallLittlewoodBGenerator", {18}, "FracQQt", "Small"},
    {"p8-Q", "HallLittlewoodQ", {8}, "FracQQt", "Small"},
    {"p8-B", "HallLittlewoodB", {8}, "FracQQt", "Small"},
    {"p8-P", "HallLittlewoodP", {8}, "FracQQt", "Small"},
    {"p8-Pomega", "HallLittlewoodPOmega", {8}, "FracQQt", "Small"}
    }

benchmarkPowerSumTargetCases = apply(benchmarkPowerSumTargetSpecs, spec ->
    benchmarkCaseRecord(
        "conv-"|spec#0,
        "BasisConversion",
        "PowerSumsToBasis",
        spec#3,
        spec#4,
        {("TargetBasis", spec#1), ("Lambda", spec#2), ("InputGroup", spec#0)}))

benchmarkRoundTripSpecs = {
    {"S-hook", "Schur", "w12-hook", "QQ", "Small"},
    {"S-three", "Schur", "w20-three-generic", "FracQQt", "Medium"},
    {"h-row", "Complete", {24}, "FracQQt", "Large"},
    {"e-row", "Elementary", {24}, "FracQQt", "Large"},
    {"Q-three", "HallLittlewoodQ", {8,4,2}, "FracQQt", "Medium"},
    {"B-three", "HallLittlewoodB", {8,4,2}, "FracQQt", "Medium"},
    {"P-three", "HallLittlewoodP", {8,4,2}, "FracQQt", "Medium"},
    {"Pomega-three", "HallLittlewoodPOmega", {8,4,2}, "FracQQt", "Medium"}
    }

benchmarkRoundTripCases = apply(benchmarkRoundTripSpecs, spec ->
    benchmarkCaseRecord(
        "roundtrip-"|spec#0,
        "RoundTrip",
        "BasisRoundTrip",
        spec#3,
        spec#4,
        {("SourceBasis", spec#1), ("Lambda", spec#2), ("InputGroup", spec#0)}))

benchmarkProductOperations = {
    {"product", "SchurProduct"},
    {"to-p", "SchurProductToPowerSums"},
    {"p-to-S", "SchurProductRoundTrip"}
    }

benchmarkProductCases = flatten apply(benchmarkProductPairs, pairData ->
    apply(benchmarkProductOperations, opData ->
        benchmarkCaseRecord(
            pairData#"ID"|"-"|opData#0,
            "SchurProduct",
            opData#1,
            "FracQQt",
            pairData#"Tier",
            {("Lambda", pairData#"Left"), ("Mu", pairData#"Right"),
             ("InputGroup", pairData#"ID"), ("PairClass", pairData#"PairClass")})))

benchmarkPieriSpecs = {
    {"pieri-hook-r1", "w12-hook", 1, "Small"},
    {"pieri-three-r2", "w13-three-generic", 2, "Medium"},
    {"pieri-four-r3", "w16-four-generic", 3, "Medium"}
    }

benchmarkPieriCases = flatten apply(benchmarkPieriSpecs, spec ->
    apply({{"horizontal", "HorizontalPieri"},
           {"vertical", "VerticalPieri"},
           {"border", "BorderStripProduct"}}, opData ->
        benchmarkCaseRecord(
            spec#0|"-"|opData#0,
            "PieriProduct",
            opData#1,
            "FracQQt",
            spec#3,
            {("Lambda", spec#1), ("Degree", spec#2), ("InputGroup", spec#0)})))

benchmarkPlethysmOperations = {
    {"combined", "SchurPlethysm"},
    {"to-p", "SchurPlethysmToPowerSums"},
    {"split-to-S", "SchurPlethysmSplit"}
    }

benchmarkPlethysmCases = flatten apply(benchmarkPlethysmPairs, pairData ->
    apply(benchmarkPlethysmOperations, opData ->
        benchmarkCaseRecord(
            pairData#"ID"|"-"|opData#0,
            "SchurPlethysm",
            opData#1,
            "FracQQt",
            pairData#"Tier",
            {("Lambda", pairData#"Left"), ("Mu", pairData#"Right"),
             ("InputGroup", pairData#"ID"), ("PairClass", pairData#"PairClass")})))

benchmarkPlethysmCrossoverSpecs = {
    {"cross-L2-r2", {6,3}, {2}, "Small"},
    {"cross-L3-r2", {7,4,2}, {2}, "Medium"},
    {"cross-L4-r2", {5,4,3,1}, {2}, "Large"},
    {"cross-L2-r3", {4,2}, {3}, "Medium"},
    {"cross-L3-r3", {5,3,1}, {3}, "Large"},
    {"cross-L2-r4", {4,2}, {4}, "Medium"},
    {"cross-L3-r4", {3,2,1}, {4}, "Medium"},
    {"cross-L2-r6", {3,2}, {6}, "Large"},
    {"cross-L2-r8", {2,1}, {8}, "Medium"}
    }

benchmarkPlethysmCrossoverCases = flatten apply(benchmarkPlethysmCrossoverSpecs, spec ->
    apply({{"combined", "SchurPlethysm"}, {"split-to-S", "SchurPlethysmSplit"}}, opData ->
        benchmarkCaseRecord(
            spec#0|"-"|opData#0,
            "PlethysmCrossover",
            opData#1,
            "QQ",
            spec#3,
            {("Lambda", spec#1), ("Mu", spec#2), ("InputGroup", spec#0),
             ("PairClass", "determinant-dilation-crossover")})))

benchmarkHallProductSpecs = {
    {"hl-small", {4,2}, {3,1}, "Small"},
    {"hl-row-three", {6}, {4,2}, "Medium"},
    {"hl-three-three", {5,3,1}, {4,3,1}, "Large"}
    }

benchmarkHallProductCases = flatten apply(benchmarkHallProductSpecs, spec ->
    apply({{"ordinary", "HallLittlewoodProduct"},
           {"retained", "HallLittlewoodProductRetained"}}, opData ->
        benchmarkCaseRecord(
            spec#0|"-"|opData#0,
            "HallLittlewoodProduct",
            opData#1,
            "FracQQt",
            spec#3,
            {("Lambda", spec#1), ("Mu", spec#2), ("InputGroup", spec#0)})))

benchmarkInnerProductCases = {
    benchmarkCaseRecord("inner-hl-pleth-small", "InnerProduct", "HallLittlewoodPlethysmInnerProduct", "FracQQt", "Medium",
        {("Lambda", {1,1}), ("Mu", {6,2}), ("Probe", {10,4,2}), ("InputGroup", "inner-hl-small")}),
    benchmarkCaseRecord("inner-hl-pleth-large", "InnerProduct", "HallLittlewoodPlethysmInnerProduct", "FracQQt", "Large",
        {("Lambda", {1,1}), ("Mu", {8,2}), ("Probe", {14,4,2}), ("InputGroup", "inner-hl-large")}),
    benchmarkCaseRecord("inner-p-S-sparse", "InnerProduct", "PowerSumsSchurInnerProduct", "FracQQt", "Small",
        {("Lambda", {8,4,2}), ("Mu", {8,4,2}), ("InputGroup", "inner-p-S-sparse")}),
    benchmarkCaseRecord("inner-p-S-cycle-heavy", "InnerProduct", "PowerSumsSchurInnerProduct", "FracQQt", "Small",
        {("Lambda", {12,2}), ("Mu", {8,4,2}), ("InputGroup", "inner-p-S-cycle-heavy")})
    }

-- Optional cases that broaden mathematical structure without making every
-- routine run pay for the full cross-product of bases and input shapes.
benchmarkExtraCases = {
    benchmarkCaseRecord("extra-conv-S-to-h", "BasisConversion-extra", "DirectBasisConversion", "FracQQt", "Small",
        {("SourceBasis", "Schur"), ("TargetBasis", "Complete"), ("Lambda", {5,4,2,1}), ("InputGroup", "extra-conv-direct-classical")}),
    benchmarkCaseRecord("extra-conv-h-to-S", "BasisConversion-extra", "DirectBasisConversion", "FracQQt", "Small",
        {("SourceBasis", "Complete"), ("TargetBasis", "Schur"), ("Lambda", {5,4,2,1}), ("InputGroup", "extra-conv-direct-classical")}),
    benchmarkCaseRecord("extra-conv-p-combination-S", "BasisConversion-extra", "PowerSumCombinationToSchur", "FracQQt", "Medium",
        {("Lambda", {8,4,2}), ("Mu", {7,4,3}), ("Probe", {6,5,3}), ("InputGroup", "extra-conv-p-combination")}),
    benchmarkCaseRecord("extra-conv-p-parameter-S", "BasisConversion-extra", "ParameterPowerSumCombinationToSchur", "FracQQt", "Medium",
        {("Lambda", {8,4,2}), ("Mu", {7,4,3}), ("Probe", {6,5,3}), ("InputGroup", "extra-conv-p-combination")}),
    benchmarkCaseRecord("extra-conv-S-combination-p", "BasisConversion-extra", "SchurCombinationToPowerSums", "FracQQt", "Medium",
        {("Lambda", {8,4,2}), ("Mu", {7,4,3}), ("Probe", {6,5,3}), ("InputGroup", "extra-conv-S-combination")}),

    benchmarkCaseRecord("extra-prod-direct-lr", "SchurProduct-extra", "SchurProductExpanded", "FracQQt", "Medium",
        {("Lambda", {4,3,1}), ("Mu", {3,3,2}), ("InputGroup", "extra-prod-lr")}),
    benchmarkCaseRecord("extra-prod-multiply-to-S", "SchurProduct-extra", "SchurProductMultiplyToBasis", "FracQQt", "Medium",
        {("Lambda", {4,3,1}), ("Mu", {3,3,2}), ("InputGroup", "extra-prod-lr")}),
    benchmarkCaseRecord("extra-prod-post-plethysm", "SchurProduct-extra", "PlethysmSchurProductToSchur", "FracQQt", "Large",
        {("Lambda", {4,2}), ("Mu", {2}), ("Probe", {2,1}), ("InputGroup", "extra-prod-post-plethysm")}),

    benchmarkCaseRecord("extra-hl-P-ordinary", "HallLittlewoodProduct-extra", "HallLittlewoodBasisProduct", "FracQQt", "Medium",
        {("SourceBasis", "HallLittlewoodP"), ("Lambda", {4,2}), ("Mu", {3,2,1}), ("InputGroup", "extra-hl-P")}),
    benchmarkCaseRecord("extra-hl-P-retained", "HallLittlewoodProduct-extra", "HallLittlewoodBasisProductRetained", "FracQQt", "Medium",
        {("SourceBasis", "HallLittlewoodP"), ("Lambda", {4,2}), ("Mu", {3,2,1}), ("InputGroup", "extra-hl-P")}),
    benchmarkCaseRecord("extra-hl-B-ordinary", "HallLittlewoodProduct-extra", "HallLittlewoodBasisProduct", "FracQQt", "Medium",
        {("SourceBasis", "HallLittlewoodB"), ("Lambda", {4,2}), ("Mu", {3,2,1}), ("InputGroup", "extra-hl-B")}),
    benchmarkCaseRecord("extra-hl-B-retained", "HallLittlewoodProduct-extra", "HallLittlewoodBasisProductRetained", "FracQQt", "Medium",
        {("SourceBasis", "HallLittlewoodB"), ("Lambda", {4,2}), ("Mu", {3,2,1}), ("InputGroup", "extra-hl-B")}),
    benchmarkCaseRecord("extra-hl-Pomega-ordinary", "HallLittlewoodProduct-extra", "HallLittlewoodBasisProduct", "FracQQt", "Medium",
        {("SourceBasis", "HallLittlewoodPOmega"), ("Lambda", {4,2}), ("Mu", {3,2,1}), ("InputGroup", "extra-hl-Pomega")}),
    benchmarkCaseRecord("extra-hl-Pomega-retained", "HallLittlewoodProduct-extra", "HallLittlewoodBasisProductRetained", "FracQQt", "Medium",
        {("SourceBasis", "HallLittlewoodPOmega"), ("Lambda", {4,2}), ("Mu", {3,2,1}), ("InputGroup", "extra-hl-Pomega")}),

    benchmarkCaseRecord("extra-inner-S-h", "InnerProduct-extra", "SchurClassicalInnerProduct", "QQ", "Medium",
        {("LeftBasis", "Schur"), ("RightBasis", "Complete"), ("Lambda", {8,4,2}), ("Mu", {6,5,3}), ("InputGroup", "extra-inner-Kostka")}),
    benchmarkCaseRecord("extra-inner-S-e", "InnerProduct-extra", "SchurClassicalInnerProduct", "QQ", "Medium",
        {("LeftBasis", "Schur"), ("RightBasis", "Elementary"), ("Lambda", {8,4,2}), ("Mu", {6,5,3}), ("InputGroup", "extra-inner-Kostka")}),
    benchmarkCaseRecord("extra-inner-Somega-h", "InnerProduct-extra", "SchurClassicalInnerProduct", "QQ", "Medium",
        {("LeftBasis", "SchurOmega"), ("RightBasis", "Complete"), ("Lambda", {8,4,2}), ("Mu", {6,5,3}), ("InputGroup", "extra-inner-Kostka-omega")}),
    benchmarkCaseRecord("extra-inner-Somega-e", "InnerProduct-extra", "SchurClassicalInnerProduct", "QQ", "Medium",
        {("LeftBasis", "SchurOmega"), ("RightBasis", "Elementary"), ("Lambda", {8,4,2}), ("Mu", {6,5,3}), ("InputGroup", "extra-inner-Kostka-omega")}),
    benchmarkCaseRecord("extra-inner-Q-P-diagonal", "InnerProduct-extra", "HallLittlewoodDiagonalInnerProduct", "FracQQt", "Medium",
        {("LeftBasis", "HallLittlewoodQ"), ("RightBasis", "HallLittlewoodP"), ("Lambda", {6,2}), ("Mu", {5,3}), ("InputGroup", "extra-inner-HL-diagonal")}),
    benchmarkCaseRecord("extra-inner-B-Pomega-diagonal", "InnerProduct-extra", "HallLittlewoodDiagonalInnerProduct", "FracQQt", "Medium",
        {("LeftBasis", "HallLittlewoodB"), ("RightBasis", "HallLittlewoodPOmega"), ("Lambda", {6,2}), ("Mu", {5,3}), ("InputGroup", "extra-inner-HL-diagonal")}),
    benchmarkCaseRecord("extra-inner-q-m-diagonal", "InnerProduct-extra", "GeneratorDualInnerProduct", "FracQQt", "Small",
        {("LeftBasis", "HallLittlewoodQGenerator"), ("RightBasis", "Monomial"), ("Lambda", {6,2}), ("Mu", {5,3}), ("InputGroup", "extra-inner-generator-diagonal")}),
    benchmarkCaseRecord("extra-inner-b-ff-diagonal", "InnerProduct-extra", "GeneratorDualInnerProduct", "FracQQt", "Small",
        {("LeftBasis", "HallLittlewoodBGenerator"), ("RightBasis", "Forgotten"), ("Lambda", {6,2}), ("Mu", {5,3}), ("InputGroup", "extra-inner-generator-diagonal")}),
    benchmarkCaseRecord("extra-inner-ordinary-expansions", "InnerProduct-extra", "OrdinaryExpansionInnerProduct", "QQ", "Medium",
        {("Lambda", {7,3,2}), ("Mu", {6,4,2}), ("InputGroup", "extra-inner-fallback")}),
    benchmarkCaseRecord("extra-inner-targeted-Q-P", "InnerProduct-extra", "TargetedHallLittlewoodInnerProduct", "FracQQt", "Medium",
        {("SourceBasis", "HallLittlewoodQ"), ("ProbeBasis", "HallLittlewoodP"), ("Lambda", {6,3,1}), ("Mu", {5,4,1}), ("Probe", {5,3,2}), ("InputGroup", "extra-inner-targeted-HL")})
    }

benchmarkCases = benchmarkUnaryCases |
                 benchmarkPowerSumTargetCases |
                 benchmarkRoundTripCases |
                 benchmarkProductCases |
                 benchmarkPieriCases |
                 benchmarkPlethysmCases |
                 benchmarkPlethysmCrossoverCases |
                 benchmarkHallProductCases |
                 benchmarkInnerProductCases |
                 benchmarkExtraCases

benchmarkCaseById = hashTable apply(benchmarkCases, x -> (x#"ID", x))

benchmarkCase = id -> (
    if not benchmarkCaseById#?id then error("unknown benchmark case: ", id);
    benchmarkCaseById#id
    )

benchmarkVariedLight = {
    "conv-p-four-S", "roundtrip-S-hook", "prod-small-balanced-p-to-S",
    "pieri-three-r2-horizontal", "pleth-two-row2-combined",
    "cross-L2-r3-combined", "hl-small-retained", "inner-p-S-sparse",
    "extra-prod-direct-lr", "extra-inner-Q-P-diagonal"
    }

benchmarkVariedStandard = unique(benchmarkVariedLight | {
    "conv-S-three-p", "conv-p30-h", "roundtrip-Q-three",
    "prod-four-three-to-p", "prod-unequal-weight-p-to-S",
    "pieri-four-r3-border", "pleth-three-row2-split-to-S",
    "pleth-hook-row-combined", "cross-L3-r4-split-to-S",
    "hl-row-three-ordinary", "inner-hl-pleth-small",
    "extra-conv-p-combination-S", "extra-conv-p-parameter-S",
    "extra-prod-post-plethysm", "extra-hl-P-retained",
    "extra-inner-S-h", "extra-inner-q-m-diagonal"
    })

benchmarkVariedThorough = unique(benchmarkVariedStandard | {
    "conv-Somega-three-p", "conv-Q-three-p", "conv-p8-Pomega",
    "roundtrip-h-row", "roundtrip-Pomega-three",
    "prod-large-three-product", "prod-large-three-p-to-S",
    "pieri-four-r3-horizontal", "pieri-four-r3-vertical",
    "pleth-three-balanced-row2-combined", "pleth-four-row2-split-to-S",
    "pleth-two-multirow-to-p", "cross-L2-r6-combined",
    "cross-L4-r2-split-to-S", "hl-three-three-retained",
    "inner-hl-pleth-large", "inner-p-S-cycle-heavy"
    } | apply(benchmarkExtraCases, case -> case#"ID"))

benchmarkVariedProfiles = hashTable {
    "light" => benchmarkVariedLight,
    "standard" => benchmarkVariedStandard,
    "thorough" => benchmarkVariedThorough
    }

benchmarkVariedCaseIds = profileName -> (
    if not benchmarkVariedProfiles#?profileName then
        error("unknown varied profile: ", profileName,
              "; expected light, standard, or thorough");
    benchmarkVariedProfiles#profileName
    )

benchmarkSelectCases = (familyFilter, tierFilter, idFilter, variedFilter) -> (
    variedIds := if variedFilter === null then null
                 else benchmarkVariedCaseIds variedFilter;
    select(benchmarkCases, case ->
        (familyFilter === null or case#"Family" == familyFilter) and
        (tierFilter === null or case#"Tier" == tierFilter) and
        (idFilter === null or case#"ID" == idFilter) and
        (variedIds === null or member(case#"ID", variedIds)))
    )
