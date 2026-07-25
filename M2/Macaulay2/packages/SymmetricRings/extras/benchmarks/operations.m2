-- Operation implementations for systematic SymmetricRings benchmarks.
-- Case records contain symbolic operation names; this file is the only place
-- where those names are translated into executable mathematics.

benchmarkKnownOperations = {
    "BasisToPowerSums", "PowerSumsToBasis", "BasisRoundTrip",
    "SchurProduct", "SchurProductToPowerSums", "SchurProductRoundTrip",
    "HorizontalPieri", "VerticalPieri", "BorderStripProduct",
    "SchurPlethysm", "SchurPlethysmToPowerSums", "SchurPlethysmSplit",
    "HallLittlewoodProduct", "HallLittlewoodProductRetained",
    "HallLittlewoodPlethysmInnerProduct", "PowerSumsSchurInnerProduct",
    "DirectBasisConversion", "PowerSumCombinationToSchur",
    "ParameterPowerSumCombinationToSchur", "SchurCombinationToPowerSums",
    "SchurProductExpanded", "SchurProductMultiplyToBasis",
    "StrictBinaryMultiplication", "BilinearMultiplicationToBasis",
    "CompleteProductTermToBasis",
    "PlethysmSchurProductToSchur", "HallLittlewoodBasisProduct",
    "HallLittlewoodBasisProductRetained", "SchurClassicalInnerProduct",
    "HallLittlewoodDiagonalInnerProduct", "GeneratorDualInnerProduct",
    "OrdinaryExpansionInnerProduct", "TargetedHallLittlewoodInnerProduct"
    }

-- Keep case IDs mathematical while routing every catalog case through the
-- public conversion and multiplication workflows.
benchmarkToBasis = (F, B) -> toBasis(F, B)

benchmarkMultiplyToBasis = (F, G, B) -> multiplyToBasis(F, G, B)

benchmarkMakeRing = ringKey -> (
    if ringKey == "QQ" then symmetricRing(QQ)
    else if ringKey == "FracQQt" then symmetricRing(frac(QQ[t]))
    else if ringKey == "QQt" then symmetricRing(QQ[t])
    else error("unknown benchmark coefficient ring: ", ringKey)
    )

benchmarkBasisElement = (basisKeyString, lambda) -> (
    if basisKeyString == "Schur" then S_lambda
    else if basisKeyString == "SchurOmega" then Somega_lambda
    else if basisKeyString == "Complete" then h_lambda
    else if basisKeyString == "Elementary" then e_lambda
    else if basisKeyString == "PowerSum" then p_lambda
    else if basisKeyString == "Monomial" then m_lambda
    else if basisKeyString == "Forgotten" then ff_lambda
    else if basisKeyString == "HallLittlewoodQGenerator" then q_lambda
    else if basisKeyString == "HallLittlewoodBGenerator" then b_lambda
    else if basisKeyString == "HallLittlewoodQ" then Q_lambda
    else if basisKeyString == "HallLittlewoodB" then B_lambda
    else if basisKeyString == "HallLittlewoodP" then P_lambda
    else if basisKeyString == "HallLittlewoodPOmega" then Pomega_lambda
    else error("unknown benchmark basis: ", basisKeyString)
    )

benchmarkBasis = basisKeyString -> basis(Rbenchmark, basisKeyString)

benchmarkLinearCombination = (basisKeyString, partitions) ->
    sum apply(partitions,
        lambda -> benchmarkBasisElement(basisKeyString, lambda))

benchmarkProductFromFactors = (basisKeys, partitions) -> (
    if #basisKeys != #partitions then
        error "benchmark factor bases and partitions have different lengths";
    product apply(#basisKeys,
        i -> benchmarkBasisElement(basisKeys#i, partitions#i))
    )

benchmarkInputPartition = (case, key) -> (
    value := case#key;
    if instance(value, String) then (benchmarkPartition value)#"Partition" else value
    )

benchmarkExecute = case -> (
    op := case#"Operation";
    lambda := if case#?"Lambda" then benchmarkInputPartition(case, "Lambda") else {};
    mu := if case#?"Mu" then benchmarkInputPartition(case, "Mu") else {};
    if op == "BasisToPowerSums" then
        benchmarkToBasis(benchmarkBasisElement(case#"SourceBasis", lambda), p)
    else if op == "PowerSumsToBasis" then
        benchmarkToBasis(p_lambda, benchmarkBasis(case#"TargetBasis"))
    else if op == "BasisRoundTrip" then (
        sourceElement := benchmarkBasisElement(case#"SourceBasis", lambda);
        benchmarkToBasis(
            benchmarkToBasis(sourceElement, p),
            benchmarkBasis(case#"SourceBasis"))
        )
    else if op == "SchurProduct" then S_lambda*S_mu
    else if op == "SchurProductToPowerSums" then
        benchmarkToBasis(S_lambda*S_mu, p)
    else if op == "SchurProductRoundTrip" then
        benchmarkToBasis(benchmarkToBasis(S_lambda*S_mu, p), S)
    else if op == "HorizontalPieri" then
        benchmarkToBasis(S_lambda*h_(case#"Degree"), S)
    else if op == "VerticalPieri" then
        benchmarkToBasis(S_lambda*e_(case#"Degree"), S)
    else if op == "BorderStripProduct" then
        benchmarkToBasis(S_lambda*p_(case#"Degree"), S)
    else if op == "SchurPlethysm" then S_lambda@S_mu
    else if op == "SchurPlethysmToPowerSums" then plethysm(S_lambda, S_mu)
    else if op == "SchurPlethysmSplit" then
        benchmarkToBasis(plethysm(S_lambda, S_mu), S)
    else if op == "HallLittlewoodProduct" then
        benchmarkToBasis(Q_lambda*Q_mu, Q)
    else if op == "HallLittlewoodProductRetained" then
        benchmarkMultiplyToBasis(Q_lambda, Q_mu, Q)
    else if op == "HallLittlewoodPlethysmInnerProduct" then
        hallInnerProduct(plethysm(Q_lambda, Q_mu), P_(case#"Probe"))
    else if op == "PowerSumsSchurInnerProduct" then
        hallInnerProduct(p_lambda, S_mu)
    else if op == "DirectBasisConversion" then
        benchmarkToBasis(
            benchmarkBasisElement(case#"SourceBasis", lambda),
            benchmarkBasis(case#"TargetBasis"))
    else if op == "PowerSumCombinationToSchur" then
        benchmarkToBasis(p_lambda + 2*p_mu + p_(case#"Probe"), S)
    else if op == "ParameterPowerSumCombinationToSchur" then
        benchmarkToBasis(
            p_lambda + t*p_mu + (t+1)*p_(case#"Probe"), S)
    else if op == "SchurCombinationToPowerSums" then
        benchmarkToBasis(S_lambda + 2*S_mu + S_(case#"Probe"), p)
    else if op == "SchurProductExpanded" then
        benchmarkToBasis(S_lambda*S_mu, S)
    else if op == "SchurProductMultiplyToBasis" then
        benchmarkMultiplyToBasis(S_lambda, S_mu, S)
    else if op == "StrictBinaryMultiplication" then
        multiplyToBasisBenchExecute(
            benchmarkBasisElement(case#"LeftBasis", lambda),
            benchmarkBasisElement(case#"RightBasis", mu),
            benchmarkBasis(case#"TargetBasis"),
            case#"Kernel",
            benchmarkTraceWorkflows)
    else if op == "BilinearMultiplicationToBasis" then
        multiplyExpressionsToBasisBenchExecute(
            benchmarkLinearCombination(
                case#"LeftBasis", case#"LeftPartitions"),
            benchmarkLinearCombination(
                case#"RightBasis", case#"RightPartitions"),
            benchmarkBasis(case#"TargetBasis"),
            benchmarkTraceWorkflows)
    else if op == "CompleteProductTermToBasis" then
        toBasisBenchExecute(
            benchmarkProductFromFactors(
                case#"FactorBases", case#"FactorPartitions"),
            benchmarkBasis(case#"TargetBasis"),
            "Automatic",
            benchmarkTraceWorkflows)
    else if op == "PlethysmSchurProductToSchur" then
        benchmarkToBasis(
            plethysm(S_lambda, S_mu)*S_(case#"Probe"), S)
    else if op == "HallLittlewoodBasisProduct" then (
        Bsource := benchmarkBasis(case#"SourceBasis");
        benchmarkToBasis(
            benchmarkBasisElement(case#"SourceBasis", lambda) *
            benchmarkBasisElement(case#"SourceBasis", mu), Bsource)
        )
    else if op == "HallLittlewoodBasisProductRetained" then (
        Bretained := benchmarkBasis(case#"SourceBasis");
        benchmarkMultiplyToBasis(
            benchmarkBasisElement(case#"SourceBasis", lambda),
            benchmarkBasisElement(case#"SourceBasis", mu), Bretained)
        )
    else if op == "SchurClassicalInnerProduct" then
        hallInnerProduct(benchmarkBasisElement(case#"LeftBasis", lambda),
                         benchmarkBasisElement(case#"RightBasis", mu),
                         "InnerProduct" => "Ordinary")
    else if op == "HallLittlewoodDiagonalInnerProduct" then
        hallInnerProduct(
            benchmarkBasisElement(case#"LeftBasis", lambda) +
                2*benchmarkBasisElement(case#"LeftBasis", mu),
            3*benchmarkBasisElement(case#"RightBasis", lambda) +
                benchmarkBasisElement(case#"RightBasis", mu))
    else if op == "GeneratorDualInnerProduct" then
        hallInnerProduct(
            benchmarkBasisElement(case#"LeftBasis", lambda) +
                2*benchmarkBasisElement(case#"LeftBasis", mu),
            3*benchmarkBasisElement(case#"RightBasis", lambda) +
                benchmarkBasisElement(case#"RightBasis", mu))
    else if op == "OrdinaryExpansionInnerProduct" then
        hallInnerProduct(S_lambda + S_mu, h_lambda + 2*h_mu,
                         "InnerProduct" => "Ordinary")
    else if op == "TargetedHallLittlewoodInnerProduct" then
        hallInnerProduct(benchmarkToBasis(
                benchmarkBasisElement(case#"SourceBasis", lambda) +
                benchmarkBasisElement(case#"SourceBasis", mu), p),
            benchmarkBasisElement(case#"ProbeBasis", case#"Probe"))
    else error("unknown benchmark operation: ", op)
    )

benchmarkExpectedWeight = case -> (
    op := case#"Operation";
    leftWeight := if case#?"Lambda" then sum benchmarkInputPartition(case, "Lambda") else 0;
    rightWeight := if case#?"Mu" then sum benchmarkInputPartition(case, "Mu") else 0;
    if case#?"ExpectedWeight" then case#"ExpectedWeight"
    else if member(op, {"SchurPlethysm", "SchurPlethysmToPowerSums", "SchurPlethysmSplit"}) then
        leftWeight*rightWeight
    else if member(op, {"SchurProduct", "SchurProductToPowerSums", "SchurProductRoundTrip",
                        "HallLittlewoodProduct", "HallLittlewoodProductRetained"}) then
        leftWeight+rightWeight
    else if member(op, {"HorizontalPieri", "VerticalPieri", "BorderStripProduct"}) then
        leftWeight+case#"Degree"
    else if op == "PlethysmSchurProductToSchur" then
        leftWeight*rightWeight + sum(case#"Probe")
    else if member(op, {"SchurProductExpanded", "SchurProductMultiplyToBasis",
                        "HallLittlewoodBasisProduct", "HallLittlewoodBasisProductRetained",
                        "StrictBinaryMultiplication"}) then
        leftWeight+rightWeight
    else if member(op, {"HallLittlewoodPlethysmInnerProduct", "PowerSumsSchurInnerProduct",
                        "SchurClassicalInnerProduct", "HallLittlewoodDiagonalInnerProduct",
                        "GeneratorDualInnerProduct", "OrdinaryExpansionInnerProduct",
                        "TargetedHallLittlewoodInnerProduct"}) then 0
    else leftWeight
    )

benchmarkVerifyResult = (case, result) -> (
    expectedWeight := benchmarkExpectedWeight case;
    if expectedWeight > 0 and weight result != expectedWeight then
        error("benchmark result has unexpected weight for ", case#"ID");
    true
    )
