-- -*- coding: utf-8 -*-
newPackage(
    "SymmetricRings",
    Version => "0.1",
    Date => "June 30, 2026",
    Authors => {
        {Name => "John Graf", Email => "jrgraf@udel.edu", HomePage => "https://j-graf.github.io/"}
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
    "SymmetricFunctionOperator",
    "RaisingOperator",
    "symmetricRing",
    "registerTransformedBasis",
    "registerSpecializedBasis",
    "applyOperator",
    "operatorData",
    "raisingOperator",
    "basisData",
    "bases",
    "aliases",
    "omegaPartners",
    "specializations",
    "innerProductPairings",
    "straighten",
    "toBasis",
    "toS",
    "toH",
    "toE",
    "toP",
    "toM",
    "toFF",
    "plethysm",
    "omegaInvolution",
    "hallInnerProduct",
    "specializeParameters",
    "hJacobiTrudi",
    "eJacobiTrudi",
    "weight",
    "rawTerms",
    "partitionWeight",
    "partitionLength"
    }

importFrom(Core, {
    "rawSymmetricRing",
    "rawSymmetricRingsSetHallLittlewoodParameter",
    "rawSymmetricRingsRememberBasis",
    "rawSymmetricRingsBasisElement",
    "rawSymmetricRingsSum",
    "rawSymmetricRingsPromoteCollected",
    "rawSymmetricRingsLiftCollected",
    "rawSymmetricRingsProduct",
    "rawSymmetricRingsJacobiTrudi",
    "rawSymmetricRingsToBasis",
    "rawSymmetricRingsProductToBasisDispatch",
    "rawSymmetricRingsPlethysm",
    "rawSymmetricRingsPlethysmToBasis",
    "rawSymmetricRingsSingleBasisId",
    "rawSymmetricRingsHasPlethysmProvenance",
    "rawSymmetricRingsCopyConversionMetadata",
    "rawSymmetricRingsOmega",
    "rawSymmetricRingsStraighten",
    "rawSymmetricRingsHallInnerProduct",
    "rawSymmetricRingsTermCount",
    "rawSymmetricRingsTermCoefficient",
    "rawSymmetricRingsTermMonomial",
    "rawSymmetricRingsElementToString",
    "rawSymmetricRingsElementToStringLimited",
    "rawSymmetricRingsElementWeight",
    "raw",
    "RawRing",
    "commonEngineRingInitializations"
    })

load "SymmetricRings/registeringBases.m2"
load "SymmetricRings/operators.m2"
load "SymmetricRings/symmetricRingsAndElements.m2"
load "SymmetricRings/computations.m2"

beginDocumentation()

load "SymmetricRings/documentation.m2"

load "SymmetricRings/tests.m2"

end
