// Copyright 2026

#include "interface/symmetric-rings.h"

#include "symmetric-rings/raw-interface.hpp"

// ============================================================================
// Ring, Basis, And Arithmetic Interface
// ============================================================================

const Ring *rawSymmetricRing(const Ring *A)
{
  return symmetric_rings::rawSymmetricRing(A);
}

bool rawSymmetricRingsSetHallLittlewoodParameter(const Ring *R,
                                                 const RingElement *t)
{
  return symmetric_rings::rawSymmetricRingsSetHallLittlewoodParameter(R, t);
}

bool rawSymmetricRingsSetComputationLimits(
    const Ring *R,
    int maxWeight,
    int maxEnumeratedPartitions,
    int maxGeneratedTerms,
    int maxRecursiveStates,
    int maxCacheEntries,
    int maxCharacterCacheEntries,
    int maxDeterminantStates,
    int maxEstimatedMemoryMB)
{
  return symmetric_rings::rawSymmetricRingsSetComputationLimits(
      R,
      maxWeight,
      maxEnumeratedPartitions,
      maxGeneratedTerms,
      maxRecursiveStates,
      maxCacheEntries,
      maxCharacterCacheEntries,
      maxDeterminantStates,
      maxEstimatedMemoryMB);
}

bool rawSymmetricRingsRememberBasis(const Ring *R,
                                    int basisId,
                                    M2_string canonicalBasisKey,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative)
{
  return symmetric_rings::rawSymmetricRingsRememberBasis(R,
                                                         basisId,
                                                         canonicalBasisKey,
                                                         displaySymbol,
                                                         displayOrder,
                                                         isMultiplicative);
}

const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                int innerLength,
                                                M2_arrayint index)
{
  return symmetric_rings::rawSymmetricRingsBasisElement(R,
                                                        basisId,
                                                        innerLength,
                                                        index);
}

const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements)
{
  return symmetric_rings::rawSymmetricRingsSum(R, elements);
}

const RingElement *rawSymmetricRingsPromoteCollected(const Ring *R,
                                                     const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsPromoteCollected(R, f);
}

const RingElement *rawSymmetricRingsLiftCollected(const Ring *R,
                                                  const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsLiftCollected(R, f);
}

const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements)
{
  return symmetric_rings::rawSymmetricRingsProduct(R, elements);
}

const RingElement *rawSymmetricRingsJacobiTrudi(const Ring *R,
                                                int basisId,
                                                M2_arrayint outer,
                                                M2_arrayint inner)
{
  return symmetric_rings::rawSymmetricRingsJacobiTrudi(R,
                                                       basisId,
                                                       outer,
                                                       inner);
}

// ============================================================================
// Conversion, Product, And Plethysm Interface
// ============================================================================

const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int targetBasisId)
{
  return symmetric_rings::rawSymmetricRingsToBasis(f, targetBasisId);
}

const RingElement *rawSymmetricRingsToBasisBench(
    const RingElement *f,
    int targetBasisId,
    M2_string conversionPlan,
    bool traceConversion)
{
  return symmetric_rings::rawSymmetricRingsToBasisBench(
      f,
      targetBasisId,
      conversionPlan,
      traceConversion);
}

const RingElement *rawSymmetricRingsMultiplyToBasis(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId)
{
  return symmetric_rings::rawSymmetricRingsMultiplyToBasis(
      f,
      g,
      targetBasisId);
}

const RingElement *rawSymmetricRingsMultiplyToBasisBench(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId,
    M2_string forcedKernel,
    bool usePowerSumReference,
    bool traceWorkflow)
{
  return symmetric_rings::rawSymmetricRingsMultiplyToBasisBench(
      f,
      g,
      targetBasisId,
      forcedKernel,
      usePowerSumReference,
      traceWorkflow);
}

const RingElement *rawSymmetricRingsMultiplyExpressionsToBasisBench(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId,
    M2_string strategy,
    bool traceWorkflow)
{
  return symmetric_rings::
      rawSymmetricRingsMultiplyExpressionsToBasisBench(
          f, g, targetBasisId, strategy, traceWorkflow);
}

const RingElement *rawSymmetricRingsPlethysm(const RingElement *f,
                                             const RingElement *g)
{
  return symmetric_rings::rawSymmetricRingsPlethysm(f, g);
}

const RingElement *rawSymmetricRingsPlethysmToBasis(const RingElement *f,
                                                    const RingElement *g,
                                                    int targetBasisId)
{
  return symmetric_rings::rawSymmetricRingsPlethysmToBasis(
      f, g, targetBasisId);
}

// ============================================================================
// Conversion Metadata Interface
// ============================================================================

int rawSymmetricRingsSingleBasisId(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsSingleBasisId(f);
}

bool rawSymmetricRingsHasPlethysmProvenance(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsHasPlethysmProvenance(f);
}

bool rawSymmetricRingsCopyConversionMetadata(const RingElement *source,
                                             const RingElement *target)
{
  return symmetric_rings::rawSymmetricRingsCopyConversionMetadata(
      source, target);
}

// ============================================================================
// General Expression Helpers
// ============================================================================

engine_RawRingElementArray rawSymmetricRingsHomogeneousComponents(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsHomogeneousComponents(f);
}

const RingElement *rawSymmetricRingsHomogeneousComponent(
    const RingElement *f,
    int weight)
{
  return symmetric_rings::
      rawSymmetricRingsHomogeneousComponent(f, weight);
}

M2_arrayint rawSymmetricRingsWeightSupport(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsWeightSupport(f);
}

const RingElement *rawSymmetricRingsTruncateWeights(
    const RingElement *f,
    int minimumWeight,
    int maximumWeight)
{
  return symmetric_rings::rawSymmetricRingsTruncateWeights(
      f, minimumWeight, maximumWeight);
}

const RingElement *rawSymmetricRingsNormalizeExpression(
    const RingElement *f,
    bool straightenIndices,
    bool expandSkewFactors,
    int productTargetBasisId)
{
  return symmetric_rings::rawSymmetricRingsNormalizeExpression(
      f,
      straightenIndices,
      expandSkewFactors,
      productTargetBasisId);
}

const RingElement *rawSymmetricRingsExpandSkewFactors(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsExpandSkewFactors(f);
}

const RingElement *rawSymmetricRingsExpandProductsInBasis(
    const RingElement *f,
    int targetBasisId)
{
  return symmetric_rings::
      rawSymmetricRingsExpandProductsInBasis(
          f, targetBasisId);
}

M2_arrayint rawSymmetricRingsExpressionShape(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsExpressionShape(f);
}

M2_arrayint rawSymmetricRingsBasisSupport(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsBasisSupport(f);
}

engine_RawRingElementArray rawSymmetricRingsBasisComponents(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsBasisComponents(f);
}

bool rawSymmetricRingsIsBasisExpansion(
    const RingElement *f,
    int basisId)
{
  return symmetric_rings::
      rawSymmetricRingsIsBasisExpansion(
          f, basisId);
}

bool rawSymmetricRingsIsLinearCombinationOfBasisElements(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsIsLinearCombinationOfBasisElements(f);
}

const RingElement *rawSymmetricRingsCoefficientsInBasis(
    const RingElement *f,
    int basisId,
    bool convert)
{
  return symmetric_rings::
      rawSymmetricRingsCoefficientsInBasis(
          f, basisId, convert);
}

engine_RawRingElementArray
rawSymmetricRingsHomogeneousBasisComponents(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsHomogeneousBasisComponents(f);
}

engine_RawRingElementArray rawSymmetricRingsSinglePartitionIndexedTerms(
    const RingElement *f)
{
  return symmetric_rings::
      rawSymmetricRingsSinglePartitionIndexedTerms(f);
}

// ============================================================================
// Omega, Straightening, And Pairing Interface
// ============================================================================

const RingElement *rawSymmetricRingsOmega(const RingElement *f,
                                          M2_arrayint omegaMap,
                                          bool useSomega)
{
  return symmetric_rings::rawSymmetricRingsOmega(f, omegaMap, useSomega);
}

const RingElement *rawSymmetricRingsStraighten(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsStraighten(f);
}

const RingElement *rawSymmetricRingsHallInnerProduct(const RingElement *f,
                                                    const RingElement *g,
                                                    int innerProductKind,
                                                    M2_arrayint innerProductMap)
{
  return symmetric_rings::rawSymmetricRingsHallInnerProduct(f,
                                                            g,
                                                            innerProductKind,
                                                            innerProductMap);
}

const RingElement *rawSymmetricRingsBasisCoefficient(
    const RingElement *f,
    const RingElement *targetBasisElement)
{
  return symmetric_rings::rawSymmetricRingsBasisCoefficient(
      f, targetBasisElement);
}

// ============================================================================
// Term And Presentation Introspection
// ============================================================================

int rawSymmetricRingsTermCount(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsTermCount(f);
}

const RingElement *rawSymmetricRingsTermCoefficient(const RingElement *f, int i)
{
  return symmetric_rings::rawSymmetricRingsTermCoefficient(f, i);
}

M2_arrayint rawSymmetricRingsTermMonomial(const RingElement *f, int i)
{
  return symmetric_rings::rawSymmetricRingsTermMonomial(f, i);
}

M2_arrayint rawSymmetricRingsPresentationTermIndices(const RingElement *f,
                                                     int maxTerms)
{
  return symmetric_rings::rawSymmetricRingsPresentationTermIndices(f, maxTerms);
}

M2_arrayint rawSymmetricRingsPresentationTermMonomial(const RingElement *f,
                                                      int i)
{
  return symmetric_rings::rawSymmetricRingsPresentationTermMonomial(f, i);
}

M2_string rawSymmetricRingsElementToString(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsElementToString(f);
}

int rawSymmetricRingsElementWeight(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsElementWeight(f);
}

// Local Variables:
// indent-tabs-mode: nil
// End:
