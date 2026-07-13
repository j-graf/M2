-- Fast structural validation of the systematic benchmark catalog.

benchmarkDirectory = currentFileDirectory
needsPackage "SymmetricRings"
load(benchmarkDirectory | "partitions.m2")
load(benchmarkDirectory | "operations.m2")
load(benchmarkDirectory | "cases.m2")

benchmarkCaseIds = apply(benchmarkCases, case -> case#"ID")
assert(#benchmarkCaseIds == #unique benchmarkCaseIds)
assert(all(benchmarkCases, case ->
    member(case#"Operation", benchmarkKnownOperations)))
assert(all(benchmarkCases, case ->
    member(case#"CoefficientRing", {"QQ", "FracQQt", "QQt"})))
assert(all(benchmarkCases, case ->
    member(case#"Tier", {"Small", "Medium", "Large", "Stress"})))
assert(all({"light", "standard", "thorough"}, profileName ->
    all(benchmarkVariedCaseIds profileName, caseId -> member(caseId, benchmarkCaseIds))))
assert(all(benchmarkVariedLight, caseId -> member(caseId, benchmarkVariedStandard)))
assert(all(benchmarkVariedStandard, caseId -> member(caseId, benchmarkVariedThorough)))

scan(benchmarkCases, case -> (
        if case#?"Lambda" then benchmarkInputPartition(case, "Lambda");
        if case#?"Mu" then benchmarkInputPartition(case, "Mu");
        benchmarkExpectedWeight case;
        ))

-- Inner-product probes must have the same total weight as their inputs.
scan(select(benchmarkCases,
        case -> case#"Operation" == "HallLittlewoodPlethysmInnerProduct"),
    case -> assert(sum(case#"Probe") ==
        (sum benchmarkInputPartition(case, "Lambda")) *
        (sum benchmarkInputPartition(case, "Mu"))))

print("validated " | toString(#benchmarkCases) | " benchmark cases")
exit 0
