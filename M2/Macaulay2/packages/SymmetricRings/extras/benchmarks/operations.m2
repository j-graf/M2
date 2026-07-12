-- Operation implementations for systematic SymmetricRings benchmarks.
-- Case records contain symbolic operation names; this file is the only place
-- where those names are translated into executable mathematics.

benchmarkKnownOperations = {
    "BasisToPowerSums", "PowerSumsToBasis", "BasisRoundTrip",
    "SchurProduct", "SchurProductToPowerSums", "SchurProductRoundTrip",
    "HorizontalPieri", "VerticalPieri", "BorderStripProduct",
    "SchurPlethysm", "SchurPlethysmToPowerSums", "SchurPlethysmSplit",
    "HallLittlewoodProduct", "HallLittlewoodProductRetained",
    "HallLittlewoodPlethysmInnerProduct", "PowerSumsSchurInnerProduct"
    }

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

benchmarkInputPartition = (case, key) -> (
    value := case#key;
    if instance(value, String) then (benchmarkPartition value)#"Partition" else value
    )

benchmarkExecute = case -> (
    op := case#"Operation";
    lambda := if case#?"Lambda" then benchmarkInputPartition(case, "Lambda") else {};
    mu := if case#?"Mu" then benchmarkInputPartition(case, "Mu") else {};
    if op == "BasisToPowerSums" then
        toBasis(benchmarkBasisElement(case#"SourceBasis", lambda), p)
    else if op == "PowerSumsToBasis" then
        toBasis(p_lambda, benchmarkBasis(case#"TargetBasis"))
    else if op == "BasisRoundTrip" then (
        sourceElement := benchmarkBasisElement(case#"SourceBasis", lambda);
        toBasis(toBasis(sourceElement, p), benchmarkBasis(case#"SourceBasis"))
        )
    else if op == "SchurProduct" then S_lambda*S_mu
    else if op == "SchurProductToPowerSums" then toBasis(S_lambda*S_mu, p)
    else if op == "SchurProductRoundTrip" then toS toBasis(S_lambda*S_mu, p)
    else if op == "HorizontalPieri" then toS(S_lambda*h_(case#"Degree"))
    else if op == "VerticalPieri" then toS(S_lambda*e_(case#"Degree"))
    else if op == "BorderStripProduct" then toS(S_lambda*p_(case#"Degree"))
    else if op == "SchurPlethysm" then S_lambda@S_mu
    else if op == "SchurPlethysmToPowerSums" then plethysm(S_lambda, S_mu)
    else if op == "SchurPlethysmSplit" then toS plethysm(S_lambda, S_mu)
    else if op == "HallLittlewoodProduct" then
        toBasis(Q_lambda*Q_mu, Q)
    else if op == "HallLittlewoodProductRetained" then
        multiplyToBasis(Q_lambda, Q_mu, Q)
    else if op == "HallLittlewoodPlethysmInnerProduct" then
        hallInnerProduct(plethysm(Q_lambda, Q_mu), P_(case#"Probe"))
    else if op == "PowerSumsSchurInnerProduct" then
        hallInnerProduct(p_lambda, S_mu)
    else error("unknown benchmark operation: ", op)
    )

benchmarkExpectedWeight = case -> (
    op := case#"Operation";
    leftWeight := if case#?"Lambda" then sum benchmarkInputPartition(case, "Lambda") else 0;
    rightWeight := if case#?"Mu" then sum benchmarkInputPartition(case, "Mu") else 0;
    if member(op, {"SchurPlethysm", "SchurPlethysmToPowerSums", "SchurPlethysmSplit"}) then
        leftWeight*rightWeight
    else if member(op, {"SchurProduct", "SchurProductToPowerSums", "SchurProductRoundTrip",
                        "HallLittlewoodProduct", "HallLittlewoodProductRetained"}) then
        leftWeight+rightWeight
    else if member(op, {"HorizontalPieri", "VerticalPieri", "BorderStripProduct"}) then
        leftWeight+case#"Degree"
    else if member(op, {"HallLittlewoodPlethysmInnerProduct", "PowerSumsSchurInnerProduct"}) then 0
    else leftWeight
    )

benchmarkVerifyResult = (case, result) -> (
    expectedWeight := benchmarkExpectedWeight case;
    if expectedWeight > 0 and weight result != expectedWeight then
        error("benchmark result has unexpected weight for ", case#"ID");
    true
    )
