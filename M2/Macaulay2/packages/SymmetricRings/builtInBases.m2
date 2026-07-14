-- ============================================================================
-- Built-In Pairing, Specialization, And Basis Metadata
-- ============================================================================

-- The built-in declarations are kept after transformed-basis support is loaded
-- so they can register standard specializations and known transformed outcomes.

-- Creates a basis-specialization map to a target basis
basisSpecializationMap = targetKey -> (R0, idx) -> (
    B := basis(R0, targetKey);
    if instance(idx, List) and #idx == 2 and instance(idx#0, List) and instance(idx#1, List) then
        makeSkewElement(B, idx#0, idx#1)
    else B_idx
    )

-- Specialization metadata for Hall-Littlewood bases at t=0
hallLittlewoodZeroSpecialization = targetKey -> {
    hashTable {
        "Parameter" => "HallLittlewoodParameter",
        "Value" => 0,
        "Map" => basisSpecializationMap targetKey
        }
    }

-- ============================================================================
-- Built-In Pairing Helpers
-- ============================================================================

-- Unit diagonal pairing used for dual bases like h/m and e/ff
unitPairing = (R0, idx) -> 1_(coefficientRing R0)

-- Integer factorial used in z_lambda
factorialZZ = n -> if n <= 1 then 1 else product toList(1..n)

-- Counts occurrences of a part in a partition-like index.
indexMultiplicity = (lambda, part) -> #select(lambda, i -> i == part)

-- Computes z_lambda for the ordinary power-sum inner product.
zValueIndex = lambda -> (
    parts := unique lambda;
    if #parts == 0 then 1
    else product(parts, part -> part^(indexMultiplicity(lambda, part)) * factorialZZ(indexMultiplicity(lambda, part)))
    )

-- Ordinary Hall inner product on the power-sum basis.
ordinaryPowerSumPairing = (R0, idx) -> promote(zValueIndex idx, coefficientRing R0)

-- Hall-Littlewood deformation of the power-sum inner product.
hallLittlewoodPowerSumPairing = (R0, idx) -> (
    A := coefficientRing R0;
    t0 := R0#"HallLittlewoodParameter";
    if t0 === null then return ordinaryPowerSumPairing(R0, idx);
    result := promote(zValueIndex idx, A);
    scan(idx, part -> result = result / (1_A - t0^part));
    result
    )

-- Builds metadata for a basis with unit pairing against a dual basis.
ordinaryDualData = dualKey -> hashTable {
    "DualBasis" => dualKey,
    "Pairing" => unitPairing,
    "EngineKind" => "Dual"
    }

-- Hall-Littlewood unit dual data has the same shape as ordinary dual data.
hallLittlewoodDualData = ordinaryDualData

-- ============================================================================
-- Built-In Basis Declarations
-- ============================================================================

-- Central registration symbols for built-in bases. The stable BasisKey is the
-- identity; these values provide the initial public notation.
builtinBasisSymbols = hashTable {
    "HallLittlewoodQ" => "Q",
    "HallLittlewoodB" => "B",
    "HallLittlewoodP" => "P",
    "HallLittlewoodPOmega" => "Pomega",
    "HallLittlewoodQGenerator" => "q",
    "HallLittlewoodBGenerator" => "b",
    "Schur" => "S",
    "SchurOmega" => "Somega",
    "Complete" => "h",
    "Elementary" => "e",
    "PowerSum" => "p",
    "Monomial" => "m",
    "Forgotten" => "ff"
    }

-- Central presentation priority for built-in bases. Basis declarations below
-- read from this table so the complete order is visible in one place.
builtinBasisDisplayOrders = hashTable {
    "HallLittlewoodQ" => 90,
    "HallLittlewoodB" => 89,
    "HallLittlewoodP" => 88,
    "HallLittlewoodPOmega" => 87,
    "HallLittlewoodQGenerator" => 80,
    "HallLittlewoodBGenerator" => 79,
    "Schur" => 70,
    "SchurOmega" => 69,
    "Complete" => 60,
    "Elementary" => 59,
    "PowerSum" => 50,
    "Monomial" => 40,
    "Forgotten" => 39
    }

-- Built-in basis registrations and their standard metadata.
-- DisplayOrder fixes presentation order; declaration order fixes the engine
-- basis ids assigned by makeBasis. Reordering or inserting built-ins here can
-- change id-based dispatch assumptions and tests that compare basis metadata.
p = makeBuiltinBasis(builtinBasisSymbols#"PowerSum", "BasisKey" => "PowerSum", "DisplayName" => "power sum basis", "DisplayOrder" => builtinBasisDisplayOrders#"PowerSum", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "Omega" => "PowerSum", "InnerProductData" => hashTable {
        "Ordinary" => hashTable {"DualBasis" => "PowerSum", "Pairing" => ordinaryPowerSumPairing, "EngineKind" => "PowerSum"},
        "HallLittlewood" => hashTable {"DualBasis" => "PowerSum", "Pairing" => hallLittlewoodPowerSumPairing, "EngineKind" => "PowerSum"}
        })
h = makeBuiltinBasis(builtinBasisSymbols#"Complete", "BasisKey" => "Complete", "DisplayName" => "complete homogeneous basis", "DisplayOrder" => builtinBasisDisplayOrders#"Complete", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "Elementary", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "Monomial"})
e = makeBuiltinBasis(builtinBasisSymbols#"Elementary", "BasisKey" => "Elementary", "DisplayName" => "elementary basis", "DisplayOrder" => builtinBasisDisplayOrders#"Elementary", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "Complete", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "Forgotten"})
m = makeBuiltinBasis(builtinBasisSymbols#"Monomial", "BasisKey" => "Monomial", "DisplayName" => "monomial basis", "DisplayOrder" => builtinBasisDisplayOrders#"Monomial", "Omega" => "Forgotten", "InnerProductData" => hashTable {
        "Ordinary" => ordinaryDualData "Complete",
        "HallLittlewood" => hallLittlewoodDualData "HallLittlewoodQGenerator"
        })
ff = makeBuiltinBasis(builtinBasisSymbols#"Forgotten", "BasisKey" => "Forgotten", "DisplayName" => "forgotten basis", "DisplayOrder" => builtinBasisDisplayOrders#"Forgotten", "Omega" => "Monomial", "InnerProductData" => hashTable {
        "Ordinary" => ordinaryDualData "Elementary",
        "HallLittlewood" => hallLittlewoodDualData "HallLittlewoodBGenerator"
        })
S = makeBuiltinBasis(builtinBasisSymbols#"Schur", "BasisKey" => "Schur", "DisplayName" => "Schur basis", "DisplayOrder" => builtinBasisDisplayOrders#"Schur", "CanBeSkew" => true, "Omega" => "SchurOmega", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "Schur"})
Somega = makeBuiltinBasis(builtinBasisSymbols#"SchurOmega", "BasisKey" => "SchurOmega", "DisplayName" => "Schur Omega basis", "DisplayOrder" => builtinBasisDisplayOrders#"SchurOmega", "CanBeSkew" => true, "Omega" => "Schur", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "SchurOmega"})
q = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodQGenerator", "BasisKey" => "HallLittlewoodQGenerator", "DisplayName" => "Hall-Littlewood q basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodQGenerator", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "HallLittlewoodBGenerator", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Complete", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "Monomial"})
b = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodBGenerator", "BasisKey" => "HallLittlewoodBGenerator", "DisplayName" => "Hall-Littlewood b basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodBGenerator", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "HallLittlewoodQGenerator", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Elementary", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "Forgotten"})
Q = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodQ", "BasisKey" => "HallLittlewoodQ", "DisplayName" => "Hall-Littlewood Q basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodQ", "CanBeSkew" => true, "Omega" => "HallLittlewoodB", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Schur", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "HallLittlewoodP"})
B = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodB", "BasisKey" => "HallLittlewoodB", "DisplayName" => "Hall-Littlewood B basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodB", "CanBeSkew" => true, "Omega" => "HallLittlewoodQ", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "SchurOmega", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "HallLittlewoodPOmega"})
P = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodP", "BasisKey" => "HallLittlewoodP", "DisplayName" => "Hall-Littlewood P basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodP", "CanBeSkew" => true, "Omega" => "HallLittlewoodPOmega", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Schur", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "HallLittlewoodQ"})
Pomega = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodPOmega", "BasisKey" => "HallLittlewoodPOmega", "DisplayName" => "Hall-Littlewood P Omega basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodPOmega", "CanBeSkew" => true, "Omega" => "HallLittlewoodP", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "SchurOmega", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "HallLittlewoodB"})

registerKnownTransformedBasisOutcome hashTable {
    "SourceBasis" => "Complete",
    "Alphabet" => "(1-t)*X",
    "SumOver" => "SameIndex",
    "TermTransformKind" => "Identity",
    "EquivalentBasis" => "HallLittlewoodQGenerator"
    }

registerKnownTransformedBasisOutcome hashTable {
    "SourceBasis" => "Elementary",
    "Alphabet" => "(1-t)*X",
    "SumOver" => "SameIndex",
    "TermTransformKind" => "Identity",
    "EquivalentBasis" => "HallLittlewoodBGenerator"
    }
