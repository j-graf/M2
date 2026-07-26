// Copyright 2026

#include "symmetric-rings/raw-interface.hpp"

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"

#include <optional>
#include <string>

namespace symmetric_rings {

namespace {

Partition partitionFromM2Array(M2_arrayint a)
{
    Partition result;
    if (a == nullptr) return result;
    result.reserve(a->len);
    for (int i = 0; i < a->len; ++i) result.push_back(a->array[i]);
    return result;
}

engine_RawRingElementArray ringElementArray(
    const SymmetricEngineRing *S,
    const std::vector<ring_elem>& elements)
{
    engine_RawRingElementArray result =
        getmemarraytype(
            engine_RawRingElementArray,
            static_cast<int>(elements.size()));
    result->len = static_cast<int>(elements.size());
    for (int i = 0; i < result->len; ++i)
      result->array[i] =
          RingElement::make_raw(S, elements[i]);
    return result;
}

M2_arrayint integerArray(const std::vector<int>& values)
{
    M2_arrayint result =
        M2_makearrayint(static_cast<int>(values.size()));
    for (size_t i = 0; i < values.size(); ++i)
      result->array[i] = values[i];
    return result;
}

} // namespace

// ============================================================================
// Ring, Basis, And Arithmetic Interface
// ============================================================================

const Ring *rawSymmetricRing(const Ring *A)
{
  try
    {
      return SymmetricEngineRing::create(A);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

bool rawSymmetricRingsSetHallLittlewoodParameter(const Ring *R,
                                                 const RingElement *t)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return false;
      bool ok = S->setHallLittlewoodParameter(t);
      if (!ok) ERROR("expected a Hall-Littlewood parameter in the coefficient ring");
      return ok;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
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
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return false;
      if (maxWeight <= 0 || maxEnumeratedPartitions <= 0 ||
          maxGeneratedTerms <= 0 ||
          maxRecursiveStates <= 0 || maxCacheEntries <= 0 ||
          maxCharacterCacheEntries <= 0 ||
          maxDeterminantStates <= 0 || maxEstimatedMemoryMB <= 0)
        {
          ERROR("symmetric-ring computation limits must be positive integers");
          return false;
        }
      S->setComputationLimits(
          static_cast<size_t>(maxWeight),
          static_cast<size_t>(maxEnumeratedPartitions),
          static_cast<size_t>(maxGeneratedTerms),
          static_cast<size_t>(maxRecursiveStates),
          static_cast<size_t>(maxCacheEntries),
          static_cast<size_t>(maxCharacterCacheEntries),
          static_cast<size_t>(maxDeterminantStates),
          static_cast<size_t>(maxEstimatedMemoryMB));
      return true;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

bool rawSymmetricRingsRememberBasis(const Ring *R,
                                    int basisId,
                                    M2_string canonicalBasisKey,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return false;
      S->rememberBasisMetadata(basisId,
                               fromM2String(canonicalBasisKey),
                               fromM2String(displaySymbol),
                               displayOrder,
                               isMultiplicative);
      if (error()) return false;
      return true;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                int innerLength,
                                                M2_arrayint index)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->basisElement(basisId, innerLength, index);
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      return RingElement::make_raw(S, S->batchSum(elements));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsPromoteCollected(const Ring *R,
                                                     const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result;
      if (!S->promoteCollectedExpansion(f, result)) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsLiftCollected(const Ring *R,
                                                  const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result;
      if (!S->liftCollectedExpansion(f, result)) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      return RingElement::make_raw(S, S->batchProduct(elements));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsJacobiTrudi(const Ring *R,
                                                int basisId,
                                                M2_arrayint outer,
                                                M2_arrayint inner)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->jacobiTrudiBasis(
          basisId, partitionFromM2Array(outer), partitionFromM2Array(inner));
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

// ============================================================================
// Conversion, Product, And Plethysm Interface
// ============================================================================

const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int targetBasisId)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->toBasis(f->get_value(), targetBasisId);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsToBasisBench(
    const RingElement *f,
    int targetBasisId,
    M2_string conversionPlan,
    bool traceConversion)
{
  try
    {
      const std::string requestedPlan =
          fromM2String(conversionPlan);
      const std::optional<std::string> forcedPlan =
          requestedPlan.empty()
              ? std::nullopt
              : std::optional<std::string>(requestedPlan);
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->toBasisBench(
          f->get_value(),
          targetBasisId,
          forcedPlan,
          traceConversion);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsMultiplyToBasis(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->multiplyToBasis(
          f->get_value(),
          g->get_value(),
          targetBasisId);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsMultiplyToBasisBench(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId,
    M2_string forcedKernel,
    bool usePowerSumReference,
    bool traceWorkflow)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      const std::string identifier =
          fromM2String(forcedKernel);
      const std::optional<std::string> forced =
          identifier.empty()
              ? std::nullopt
              : std::optional<std::string>(identifier);
      ring_elem result = S->multiplyToBasisBench(
          f->get_value(),
          g->get_value(),
          targetBasisId,
          forced,
          usePowerSumReference,
          traceWorkflow);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsMultiplyExpressionsToBasisBench(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId,
    M2_string strategy,
    bool traceWorkflow)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result =
          S->multiplyExpressionsToBasisBench(
              f->get_value(),
              g->get_value(),
              targetBasisId,
              fromM2String(strategy),
              traceWorkflow);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsPlethysm(const RingElement *f,
                                             const RingElement *g)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->plethysm(f->get_value(), g->get_value());
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsPlethysmToBasis(const RingElement *f,
                                                    const RingElement *g,
                                                    int targetBasisId)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->plethysmToBasisDispatch(f->get_value(),
                                                     g->get_value(),
                                                     targetBasisId);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

// ============================================================================
// Expression-Facts Cache Interface
// ============================================================================

int rawSymmetricRingsSingleBasisId(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return -1;
      return S->uniformBasisId(f->get_value());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return -1;
    }
}

bool rawSymmetricRingsHasPlethysmProvenance(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return false;
      (void) S;
      return hasCombinatorialTag(
          polyValue(f->get_value())->combinatorialTags,
          CombinatorialTag::Plethysm);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
  }

bool rawSymmetricRingsCopyConversionMetadata(const RingElement *source,
                                             const RingElement *target)
{
  // Keep the established raw name for M2 compatibility. The transported
  // payload is now the shared expression-facts cache plus provenance tags.
  try
    {
      const auto *sourceRing = symmetricRingFromElement(source);
      if (error()) return false;
      const auto *targetRing = symmetricRingFromElement(target);
      if (error()) return false;
      auto cache = polyValue(source->get_value())->expressionFactsCache;
      if (cache && sourceRing != targetRing)
        {
          // Basis IDs are ring-local. The M2 fallback reconstructs atoms on
          // the target ring, so basis identity and the exact cached facts
          // must be recomputed together.
          cache->discardRingLocalBasisFacts();
        }
      mutablePolyValue(target->get_value())->expressionFactsCache =
          std::move(cache);
      mutablePolyValue(target->get_value())->combinatorialTags =
          polyValue(source->get_value())->combinatorialTags;
      return true;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

// ============================================================================
// General Expression Helpers
// ============================================================================

engine_RawRingElementArray rawSymmetricRingsHomogeneousComponents(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      const auto components =
          S->homogeneousComponents(f->get_value());
      std::vector<ring_elem> expressions;
      expressions.reserve(components.size());
      for (const auto& component : components)
        expressions.push_back(component.expression);
      return ringElementArray(S, expressions);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsHomogeneousComponent(
    const RingElement *f,
    int weight)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return RingElement::make_raw(
          S,
          S->homogeneousComponent(
              f->get_value(), weight));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsWeightSupport(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return integerArray(
          S->weightSupport(f->get_value()));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsTruncateWeights(
    const RingElement *f,
    int minimumWeight,
    int maximumWeight)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->truncateWeights(
          f->get_value(),
          minimumWeight,
          maximumWeight);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsNormalizeExpression(
    const RingElement *f,
    bool straightenIndices,
    bool expandSkewFactors,
    int productTargetBasisId)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      SymmetricEngineRing::ExpressionNormalizationOptions options;
      options.straightenIndices = straightenIndices;
      options.expandSkewFactors = expandSkewFactors;
      options.productTargetBasisId = productTargetBasisId;
      ring_elem result =
          S->normalizeExpressionWithOptions(
              f->get_value(), options);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsExpandSkewFactors(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result =
          S->expandSkewFactors(f->get_value());
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsExpandProductsInBasis(
    const RingElement *f,
    int targetBasisId)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->expandProductsInBasis(
          f->get_value(), targetBasisId);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsExpressionShape(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      const auto shape =
          S->expressionShape(f->get_value());
      std::vector<int> values{
          2,
          static_cast<int>(shape.termCount),
          static_cast<int>(shape.scalarTermCount),
          static_cast<int>(shape.singleFactorTermCount),
          static_cast<int>(shape.productTermCount),
          static_cast<int>(shape.maximumFactorsPerTerm),
          static_cast<int>(shape.skewFactorCount),
          static_cast<int>(shape.maximumPartitionLength),
          shape.normalized ? 1 : 0,
          shape.skewFree ? 1 : 0,
          shape.collected ? 1 : 0,
          shape.homogeneousWeight.value_or(-1),
          shape.pureBasis.value_or(-1),
          shape.expandedBasis.value_or(-1),
          static_cast<int>(shape.basisIds.size())};
      values.insert(
          values.end(),
          shape.basisIds.begin(),
          shape.basisIds.end());
      values.push_back(
          static_cast<int>(shape.weights.size()));
      values.insert(
          values.end(),
          shape.weights.begin(),
          shape.weights.end());
      values.push_back(shape.hasMetadata ? 1 : 0);
      values.push_back(
          shape.metadataFactsComplete ? 1 : 0);
      values.push_back(
          shape.metadataNormalized ? 1 : 0);
      values.push_back(
          shape.metadataSkewFree ? 1 : 0);
      values.push_back(
          shape.metadataCollected ? 1 : 0);
      return integerArray(values);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsBasisSupport(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return integerArray(
          S->basisSupport(f->get_value()));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

engine_RawRingElementArray rawSymmetricRingsBasisComponents(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      const auto components =
          S->basisComponents(f->get_value());
      if (error()) return nullptr;
      std::vector<ring_elem> expressions;
      expressions.reserve(components.size());
      for (const auto& component : components)
        expressions.push_back(component.expression);
      return ringElementArray(S, expressions);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

bool rawSymmetricRingsIsBasisExpansion(
    const RingElement *f,
    int basisId)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return false;
      return S->isBasisExpansion(
          f->get_value(), basisId);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

bool rawSymmetricRingsIsLinearCombinationOfBasisElements(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return false;
      return S->isLinearCombinationOfBasisElements(
          f->get_value());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

const RingElement *rawSymmetricRingsCoefficientsInBasis(
    const RingElement *f,
    int basisId,
    bool convert)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result =
          S->coefficientsInBasisExpression(
              f->get_value(), basisId, convert);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

engine_RawRingElementArray
rawSymmetricRingsHomogeneousBasisComponents(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      const auto components =
          S->homogeneousBasisComponents(
              f->get_value());
      if (error()) return nullptr;
      std::vector<ring_elem> expressions;
      expressions.reserve(components.size());
      for (const auto& component : components)
        expressions.push_back(component.expression);
      return ringElementArray(S, expressions);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

engine_RawRingElementArray rawSymmetricRingsSinglePartitionIndexedTerms(
    const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return ringElementArray(
          S, S->singlePartitionIndexedTerms(f->get_value()));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

// ============================================================================
// Omega, Straightening, And Pairing Interface
// ============================================================================

const RingElement *rawSymmetricRingsOmega(const RingElement *f,
                                          M2_arrayint omegaMap,
                                          bool useSomega)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->omegaInvolution(f->get_value(), omegaMap, useSomega);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsStraighten(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->straighten(f->get_value());
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsHallInnerProduct(const RingElement *f,
                                                    const RingElement *g,
                                                    int innerProductKind,
                                                    M2_arrayint innerProductMap)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->hallInnerProduct(f->get_value(),
                                             g->get_value(),
                                             innerProductKind,
                                             innerProductMap);
      if (error()) return nullptr;
      return RingElement::make_raw(S->getCoefficientRing(), result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsBasisCoefficient(
    const RingElement *f,
    const RingElement *targetBasisElement)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (targetBasisElement->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->basisCoefficient(
          f->get_value(), targetBasisElement->get_value());
      if (error()) return nullptr;
      return RingElement::make_raw(S->getCoefficientRing(), result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

// ============================================================================
// Term And Presentation Introspection
// ============================================================================

int rawSymmetricRingsTermCount(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return 0;
      (void) S;
      return static_cast<int>(polyValue(f->get_value())->terms.size());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return 0;
    }
}

const RingElement *rawSymmetricRingsTermCoefficient(const RingElement *f, int i)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      const auto *poly = polyValue(f->get_value());
      if (i < 0 || static_cast<size_t>(i) >= poly->terms.size())
        {
          ERROR("symmetric-ring term index out of range");
          return nullptr;
        }
      const Ring *A = S->getCoefficientRing();
      return RingElement::make_raw(A, A->copy(poly->terms[i].coeff));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsTermMonomial(const RingElement *f, int i)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      (void) S;
      const auto *poly = polyValue(f->get_value());
      if (i < 0 || static_cast<size_t>(i) >= poly->terms.size())
        {
          ERROR("symmetric-ring term index out of range");
          return nullptr;
        }
      const auto& data = poly->terms[i].monomial.data;
      M2_arrayint result = M2_makearrayint(static_cast<int>(data.size()));
      for (size_t j = 0; j < data.size(); ++j)
        result->array[j] = data[j];
      return result;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsPresentationTermIndices(const RingElement *f,
                                                     int maxTerms)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      std::vector<size_t> order = S->presentationTermOrder(f->get_value());
      // Presentation always sorts the complete view before applying its limit.
      size_t count = order.size();
      if (maxTerms >= 0)
        count = std::min(count, static_cast<size_t>(maxTerms));
      M2_arrayint result = M2_makearrayint(static_cast<int>(count));
      for (size_t i = 0; i < count; ++i)
        result->array[i] = static_cast<int>(order[i]);
      return result;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsPresentationTermMonomial(const RingElement *f,
                                                      int i)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      const auto *poly = polyValue(f->get_value());
      if (i < 0 || static_cast<size_t>(i) >= poly->terms.size())
        {
          ERROR("symmetric-ring term index out of range");
          return nullptr;
        }
      std::vector<int> data =
          S->presentationMonomialData(poly->terms[i].monomial);
      M2_arrayint result = M2_makearrayint(static_cast<int>(data.size()));
      for (size_t j = 0; j < data.size(); ++j) result->array[j] = data[j];
      return result;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_string rawSymmetricRingsElementToString(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return toM2String(S->elementString(f->get_value()));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

int rawSymmetricRingsElementWeight(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return -1;
      return S->elementWeight(f->get_value());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return -1;
    }
}


} // namespace symmetric_rings
