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
    "registerTransformedBasis",
    "basisData",
    "bases",
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
    "rawSymmetricRingsProduct",
    "rawSymmetricRingsJacobiTrudi",
    "rawSymmetricRingsToBasis",
    "rawSymmetricRingsToSchurViaHRecursive",
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
    "rawSymmetricRingsElementToStringLimited",
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
