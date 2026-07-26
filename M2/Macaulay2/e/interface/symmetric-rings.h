// Copyright 2026

#ifndef M2_INTERFACE_SYMMETRIC_RINGS_H_
#define M2_INTERFACE_SYMMETRIC_RINGS_H_

#include "engine-includes.hpp"

#if defined(__cplusplus)
class Ring;
class RingElement;
#else
typedef struct Ring Ring;
typedef struct RingElement RingElement;
#endif

#if defined(__cplusplus)
extern "C" {
#endif

// ============================================================================
// Ring, Basis, And Arithmetic Interface
// ============================================================================

const Ring *rawSymmetricRing(const Ring *A);
bool rawSymmetricRingsSetHallLittlewoodParameter(const Ring *R,
                                                 const RingElement *t);
bool rawSymmetricRingsSetComputationLimits(
    const Ring *R,
    int maxWeight,
    int maxEnumeratedPartitions,
    int maxGeneratedTerms,
    int maxRecursiveStates,
    int maxCacheEntries,
    int maxCharacterCacheEntries,
    int maxDeterminantStates,
    int maxEstimatedMemoryMB);
bool rawSymmetricRingsRememberBasis(const Ring *R,
                                    int basisId,
                                    M2_string canonicalBasisKey,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative);
const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                int innerLength,
                                                M2_arrayint index);
const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements);
const RingElement *rawSymmetricRingsPromoteCollected(const Ring *R,
                                                     const RingElement *f);
const RingElement *rawSymmetricRingsLiftCollected(const Ring *R,
                                                  const RingElement *f);
const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements);
const RingElement *rawSymmetricRingsJacobiTrudi(const Ring *R,
                                                int basisId,
                                                M2_arrayint outer,
                                                M2_arrayint inner);

// ============================================================================
// Conversion, Product, And Plethysm Interface
// ============================================================================

const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int targetBasisId);
const RingElement *rawSymmetricRingsToBasisBench(
    const RingElement *f,
    int targetBasisId,
    M2_string conversionPlan,
    bool traceConversion);
const RingElement *rawSymmetricRingsMultiplyToBasis(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId);
const RingElement *rawSymmetricRingsMultiplyToBasisBench(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId,
    M2_string forcedKernel,
    bool usePowerSumReference,
    bool traceWorkflow);
const RingElement *rawSymmetricRingsMultiplyExpressionsToBasisBench(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId,
    M2_string strategy,
    bool traceWorkflow);
const RingElement *rawSymmetricRingsPlethysm(const RingElement *f,
                                             const RingElement *g);
const RingElement *rawSymmetricRingsPlethysmToBasis(const RingElement *f,
                                                    const RingElement *g,
                                                    int targetBasisId);

// ============================================================================
// Conversion Metadata Interface
// ============================================================================

int rawSymmetricRingsSingleBasisId(const RingElement *f);
bool rawSymmetricRingsHasPlethysmProvenance(const RingElement *f);
bool rawSymmetricRingsCopyConversionMetadata(const RingElement *source,
                                             const RingElement *target);

// ============================================================================
// General Expression Helpers
// ============================================================================

engine_RawRingElementArray rawSymmetricRingsHomogeneousComponents(
    const RingElement *f);
const RingElement *rawSymmetricRingsHomogeneousComponent(
    const RingElement *f,
    int weight);
M2_arrayint rawSymmetricRingsWeightSupport(const RingElement *f);
const RingElement *rawSymmetricRingsTruncateWeights(
    const RingElement *f,
    int minimumWeight,
    int maximumWeight);
const RingElement *rawSymmetricRingsNormalizeExpression(
    const RingElement *f,
    bool straightenIndices,
    bool expandSkewFactors,
    int productTargetBasisId);
const RingElement *rawSymmetricRingsExpandSkewFactors(
    const RingElement *f);
const RingElement *rawSymmetricRingsExpandProductsInBasis(
    const RingElement *f,
    int targetBasisId);
M2_arrayint rawSymmetricRingsExpressionShape(const RingElement *f);
M2_arrayint rawSymmetricRingsBasisSupport(const RingElement *f);
engine_RawRingElementArray rawSymmetricRingsBasisComponents(
    const RingElement *f);
bool rawSymmetricRingsIsBasisExpansion(
    const RingElement *f,
    int basisId);
bool rawSymmetricRingsIsLinearCombinationOfBasisElements(
    const RingElement *f);
const RingElement *rawSymmetricRingsCoefficientsInBasis(
    const RingElement *f,
    int basisId,
    bool convert);
engine_RawRingElementArray
rawSymmetricRingsHomogeneousBasisComponents(
    const RingElement *f);
engine_RawRingElementArray rawSymmetricRingsSinglePartitionIndexedTerms(
    const RingElement *f);

// ============================================================================
// Omega, Straightening, And Pairing Interface
// ============================================================================

const RingElement *rawSymmetricRingsOmega(const RingElement *f,
                                          M2_arrayint omegaMap,
                                          bool useSomega);
const RingElement *rawSymmetricRingsStraighten(const RingElement *f);
const RingElement *rawSymmetricRingsHallInnerProduct(const RingElement *f,
                                                    const RingElement *g,
                                                    int innerProductKind,
                                                    M2_arrayint innerProductMap);
const RingElement *rawSymmetricRingsBasisCoefficient(
    const RingElement *f,
    const RingElement *targetBasisElement);

// ============================================================================
// Term And Presentation Introspection
// ============================================================================

int rawSymmetricRingsTermCount(const RingElement *f);
const RingElement *rawSymmetricRingsTermCoefficient(const RingElement *f, int i);
M2_arrayint rawSymmetricRingsTermMonomial(const RingElement *f, int i);
M2_arrayint rawSymmetricRingsPresentationTermIndices(const RingElement *f,
                                                     int maxTerms);
M2_arrayint rawSymmetricRingsPresentationTermMonomial(const RingElement *f,
                                                      int i);
M2_string rawSymmetricRingsElementToString(const RingElement *f);
int rawSymmetricRingsElementWeight(const RingElement *f);

#if defined(__cplusplus)
}
#endif

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
