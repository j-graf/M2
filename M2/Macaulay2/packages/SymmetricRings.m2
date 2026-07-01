-- -*- coding: utf-8 -*-
newPackage(
    "SymmetricRings",
    Version => "0.1",
    Date => "June 30, 2026",
    Authors => {
        {Name => "John Graf"}
        },
    Headline => "formal mixed-basis symmetric function rings",
    Keywords => {"Representation Theory"},
    AuxiliaryFiles => true,
    DebuggingMode => false
    )

export {
    "SymmetricRing",
    "SymmetricRingElement",
    "SymmetricBasis",
    "symmetricRing",
    "registerBasis",
    "basisData",
    "bases",
    "symmetricEquals",
    "straighten",
    "toBasis",
    "plethysm",
    "omegaInvolution",
    "hallInnerProduct",
    "specializeParameters",
    "hJacobiTrudi",
    "eJacobiTrudi",
    "weight",
    "partitionWeight",
    "partitionLength",
    "Bases",
    "Parameters",
    "HallLittlewoodParameter",
    "MacdonaldParameters",
    "ParameterSpecialization",
    "DefaultSeriesVariables",
    "AvailableWhen",
    "DisplayName",
    "DisplayOrder",
    "CanBeSkew",
    "IndexNormalizer",
    "IndexValidator",
    "Constructor",
    "IsMultiplicativeIndex",
    "MultiplicativeIndex",
    "Straighten",
    "ToPowerSums",
    "FromPowerSums",
    "TriangularData",
    "Omega",
    "Specialization",
    "PromoteSpecializedRing",
    "InnerProductData",
    "PlethysmBehavior",
    "Display",
    "Documentation",
    "ZeroIndexIsOne",
    "ZeroOnNegative"
    }

exportMutable {
    "builtinSymmetricBases",
    "userDefinedSymmetricBases",
    "availableSymmetricBases"
    }

importFrom(Core, {
    "rawSymmetricRing",
    "rawSymmetricRingsSetHallLittlewoodParameter",
    "rawSymmetricRingsRememberBasis",
    "rawSymmetricRingsBasisElement",
    "rawSymmetricRingsSum",
    "rawSymmetricRingsProduct",
    "rawSymmetricRingsJacobiTrudi",
    "rawSymmetricRingsToBasis",
    "rawSymmetricRingsPlethysm",
    "rawSymmetricRingsPlethysmToBasis",
    "rawSymmetricRingsSingleBasisId",
    "rawSymmetricRingsOmega",
    "rawSymmetricRingsStraighten",
    "rawSymmetricRingsHallInnerProduct",
    "rawSymmetricRingsTermCount",
    "rawSymmetricRingsTermCoefficient",
    "rawSymmetricRingsTermMonomial",
    "rawSymmetricRingsElementToString",
    "rawSymmetricRingsElementWeight",
    "raw",
    "RawRing",
    "commonEngineRingInitializations"
    })

load "SymmetricRings/symmetricRingsCore.m2"

beginDocumentation()

load "SymmetricRings/documentation.m2"

load "SymmetricRings/tests.m2"

end
