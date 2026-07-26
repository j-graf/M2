// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"

#include <algorithm>
#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Expression-Fact Inference
// ============================================================================
// Facts in this section are exact descriptions of canonical engine storage.
// Derived predicates live on ExpressionFacts so selectors do not maintain
// redundant flags that can disagree with the underlying counts.

SymmetricEngineRing::ExpressionFacts
SymmetricEngineRing::inferExpressionFacts(
    ring_elem f,
    std::vector<size_t> *termFactorCounts) const
{
    ExpressionFacts facts;
    const auto *poly = polyValue(f);
    facts.combinatorialTags = poly->combinatorialTags;
    facts.termCount = poly->terms.size();
    std::optional<int> firstWeight;
    bool multipleWeights = false;
    if (termFactorCounts != nullptr)
      termFactorCounts->clear();

    size_t termPosition = 0;
    for (const auto& term : poly->terms)
      {
        size_t factorCount = 0;
        long long termWeight = 0;
        size_t pos = 0;
        while (pos < term.monomial.data.size())
          {
            ++factorCount;
            int basisId = atomBasisIdAt(term.monomial, pos);
            if (std::find(
                    facts.factorBases.begin(),
                    facts.factorBases.end(),
                    basisId) == facts.factorBases.end())
              facts.factorBases.push_back(basisId);
            if (atomIsSkewAt(term.monomial, pos))
              {
                facts.skewFree = false;
                ++facts.skewFactorCount;
                Partition outer = basisElementOuterIndex(term.monomial, pos);
                Partition inner = basisElementInnerIndex(term.monomial, pos);
                facts.maximumPartitionLength =
                    std::max(facts.maximumPartitionLength, outer.size());
                termWeight +=
                    static_cast<long long>(partitionWeight(outer)) -
                    static_cast<long long>(partitionWeight(inner));
                facts.normalized =
                    facts.normalized && isPartitionIndex(outer) &&
                    isPartitionIndex(inner) && partitionContains(outer, inner);
              }
            else
              {
                Partition index = basisElementIndex(term.monomial, pos);
                facts.maximumPartitionLength =
                    std::max(facts.maximumPartitionLength, index.size());
                termWeight += partitionWeight(index);
                if (!isMultiplicativeBasis(basisId))
                  facts.normalized =
                      facts.normalized && isPartitionIndex(index);
              }
            pos += atomLengthAt(term.monomial, pos);
          }

        facts.maximumFactorsPerTerm =
            std::max(facts.maximumFactorsPerTerm, factorCount);
        if (termFactorCounts != nullptr)
          {
            if (!termFactorCounts->empty())
              termFactorCounts->push_back(factorCount);
            else if (factorCount > 1)
              {
                termFactorCounts->reserve(facts.termCount);
                for (size_t i = 0; i < termPosition; ++i)
                  termFactorCounts->push_back(
                      poly->terms[i].monomial.data.empty() ? 0 : 1);
                termFactorCounts->push_back(factorCount);
              }
          }
        if (factorCount == 0)
          ++facts.scalarTermCount;
        else if (factorCount == 1)
          ++facts.singleFactorTermCount;
        else
          ++facts.productTermCount;
        if (termWeight < std::numeric_limits<int>::min() ||
            termWeight > std::numeric_limits<int>::max())
          throw exc::engine_error(
              "symmetric-function weight exceeds the supported integer range");
        const int weight = static_cast<int>(termWeight);
        if (!firstWeight)
          firstWeight = weight;
        else if (*firstWeight != weight)
          multipleWeights = true;
        ++termPosition;
      }

    std::sort(facts.factorBases.begin(), facts.factorBases.end());
    if (facts.factorBases.size() == 1)
      facts.pureBasis = facts.factorBases.front();
    if (facts.productTermCount == 0 &&
        facts.singleFactorTermCount + facts.scalarTermCount == facts.termCount &&
        facts.factorBases.size() == 1)
      facts.expandedBasis = facts.factorBases.front();
    if (firstWeight && !multipleWeights)
      facts.homogeneousWeight = *firstWeight;

    if (facts.singleBasisElement())
      {
        const auto& term = poly->terms.front();
        const auto& monomial = term.monomial;
        int basisId = atomBasisIdAt(monomial, 0);
        facts.singleBasisElementId = basisId;
        facts.singleBasisElementCoefficientOne =
            coefficientRing->is_equal(
                term.coeff, coefficientRing->one());
        if (!atomIsSkewAt(monomial, 0))
          {
            Partition index = basisElementIndex(monomial, 0);
            facts.singleBasisElementIndex = index;
          }
      }
    return facts;
  }

SymmetricEngineRing::ExpressionFacts
SymmetricEngineRing::inferCanonicalExpansionFacts(
    ring_elem f,
    int basisId,
    const std::optional<int>& knownHomogeneousWeight) const
{
    // A conversion plan or strict binary-kernel contract promises one
    // normalized, non-skew target atom per nonscalar term. Derive the common
    // output facts from that guarantee and inspect only support that cannot be
    // predicted, rather than running the general multi-factor analyzer.
    ExpressionFacts facts;
    const auto *poly = polyValue(f);
    facts.combinatorialTags = poly->combinatorialTags;
    facts.termCount = poly->terms.size();
    std::optional<int> firstWeight;
    bool multipleWeights = false;

    for (const auto& term : poly->terms)
      {
        if (term.monomial.data.empty())
          {
            ++facts.scalarTermCount;
            if (!knownHomogeneousWeight)
              {
                if (!firstWeight)
                  firstWeight = 0;
                else if (*firstWeight != 0)
                  multipleWeights = true;
              }
            continue;
          }
        if (atomLengthAt(term.monomial, 0) !=
                term.monomial.data.size() ||
            atomBasisIdAt(term.monomial, 0) != basisId ||
            atomIsSkewAt(term.monomial, 0))
          {
            facts.normalized = false;
            facts.skewFree = false;
            ++facts.productTermCount;
            continue;
          }
        facts.maximumFactorsPerTerm = 1;
        ++facts.singleFactorTermCount;
        facts.maximumPartitionLength = std::max(
            facts.maximumPartitionLength,
            static_cast<size_t>(
                atomIndexLengthAt(term.monomial, 0)));
        if (!knownHomogeneousWeight)
          {
            const int weight = monomialWeight(term.monomial);
            if (!firstWeight)
              firstWeight = weight;
            else if (*firstWeight != weight)
              multipleWeights = true;
          }
      }
    if (facts.singleFactorTermCount != 0)
      {
        facts.factorBases = {basisId};
        facts.pureBasis = basisId;
        facts.expandedBasis = basisId;
      }
    if (knownHomogeneousWeight)
      facts.homogeneousWeight = *knownHomogeneousWeight;
    else if (firstWeight && !multipleWeights)
      facts.homogeneousWeight = *firstWeight;

    if (facts.singleBasisElement())
      {
        const auto& term = poly->terms.front();
        facts.singleBasisElementId = basisId;
        facts.singleBasisElementIndex =
            basisElementIndex(term.monomial, 0);
        facts.singleBasisElementCoefficientOne =
            coefficientRing->is_equal(
                term.coeff, coefficientRing->one());
      }
    return facts;
  }

SymmetricEngineRing::ExpressionFacts
SymmetricEngineRing::basisElementFactsFromAtom(
    const SymmetricMonomial& monomial,
    size_t pos,
    ring_elem coefficient) const
{
    ExpressionFacts facts;
    const int basisId = atomBasisIdAt(monomial, pos);
    facts.termCount = 1;
    facts.singleFactorTermCount = 1;
    facts.maximumFactorsPerTerm = 1;
    facts.factorBases = {basisId};
    facts.pureBasis = basisId;
    facts.expandedBasis = basisId;
    facts.singleBasisElementId = basisId;
    facts.singleBasisElementCoefficientOne =
        coefficientRing->is_equal(
            coefficient, coefficientRing->one());
    int weight = 0;

    if (atomIsSkewAt(monomial, pos))
      {
        facts.skewFree = false;
        facts.skewFactorCount = 1;
        Partition outer = basisElementOuterIndex(monomial, pos);
        Partition inner = basisElementInnerIndex(monomial, pos);
        facts.maximumPartitionLength = outer.size();
        facts.normalized =
            isPartitionIndex(outer) &&
            isPartitionIndex(inner) &&
            partitionContains(outer, inner);
        weight =
            partitionWeight(outer) - partitionWeight(inner);
      }
    else
      {
        Partition index = basisElementIndex(monomial, pos);
        facts.singleBasisElementIndex = index;
        facts.maximumPartitionLength = index.size();
        if (!isMultiplicativeBasis(basisId))
          facts.normalized = isPartitionIndex(index);
        weight = partitionWeight(index);
      }
    facts.homogeneousWeight = weight;
    return facts;
  }

std::optional<SymmetricEngineRing::ExpressionFacts>
SymmetricEngineRing::expressionFactsFromMetadata(ring_elem f) const
{
    const auto *poly = polyValue(f);
    if (!poly->conversionMetadata) return std::nullopt;
    const auto& metadata = *poly->conversionMetadata;
    if (!metadata.expressionFactsComplete ||
        !metadata.normalized || !metadata.skewFree || !metadata.collected ||
        !metadata.termCount)
      return std::nullopt;

    ExpressionFacts facts;
    facts.normalized = true;
    facts.skewFree = true;
    facts.collected = true;
    facts.combinatorialTags = poly->combinatorialTags;
    facts.termCount = *metadata.termCount;
    if (!metadata.scalarTermCount ||
        !metadata.singleFactorTermCount)
      return std::nullopt;
    facts.scalarTermCount = *metadata.scalarTermCount;
    facts.singleFactorTermCount = *metadata.singleFactorTermCount;
    facts.productTermCount = 0;
    facts.maximumFactorsPerTerm =
        facts.singleFactorTermCount == 0 ? 0 : 1;
    if (facts.termCount != poly->terms.size() ||
        facts.scalarTermCount + facts.singleFactorTermCount +
            facts.productTermCount != facts.termCount ||
        facts.productTermCount != 0)
      return std::nullopt;
    // A zero or scalar expansion has no basis support. A product-free mixed
    // expansion is also pipeline-ready, but has no single expanded basis.
    if (!metadata.factorBases)
      return std::nullopt;
    facts.factorBases = *metadata.factorBases;
    if ((facts.singleFactorTermCount == 0 &&
         (!facts.factorBases.empty() || metadata.expandedBasis)) ||
        (facts.singleFactorTermCount != 0 &&
         facts.factorBases.empty()))
      return std::nullopt;
    if (metadata.expandedBasis &&
        (facts.factorBases.size() != 1 ||
         facts.factorBases.front() !=
             *metadata.expandedBasis))
      return std::nullopt;
    if (!metadata.expandedBasis &&
        facts.factorBases.size() == 1)
      return std::nullopt;
    if (metadata.pureBasis &&
        (facts.factorBases.size() != 1 ||
         facts.factorBases.front() != *metadata.pureBasis))
      return std::nullopt;
    facts.pureBasis = metadata.pureBasis;
    facts.expandedBasis = metadata.expandedBasis;
    facts.skewFactorCount = 0;
    facts.homogeneousWeight = metadata.homogeneousWeight;
    facts.maximumPartitionLength =
        metadata.maximumPartitionLength.value_or(0);
    facts.singleBasisElementId =
        metadata.singleBasisElementId;
    facts.singleBasisElementIndex =
        metadata.singleBasisElementIndex;
    facts.singleBasisElementCoefficientOne =
        metadata.singleBasisElementCoefficientOne;
    const bool singleBasisElement =
        facts.termCount == 1 &&
        facts.scalarTermCount == 0 &&
        facts.singleFactorTermCount == 1;
    if (singleBasisElement !=
            static_cast<bool>(facts.singleBasisElementId) ||
        (facts.singleBasisElementIndex &&
         !facts.singleBasisElementId) ||
        (facts.singleBasisElementCoefficientOne &&
         !facts.singleBasisElementId) ||
        (facts.singleBasisElementId &&
         (!facts.expandedBasis ||
          *facts.singleBasisElementId != *facts.expandedBasis)))
      return std::nullopt;
    return facts;
  }

// ============================================================================
// Metadata Lifecycle
// ============================================================================
// Arithmetic may preserve inexpensive hints even when cancellation or product
// formation prevents it from proving the complete canonical core. All such
// construction and invalidation goes through these helpers so no operation
// invents its own interpretation of SymmetricConversionMetadata.

void SymmetricEngineRing::invalidateExactExpressionFacts(
    SymmetricConversionMetadataSlot& metadata) const
{
    if (!metadata) return;
    metadata->invalidateExactExpressionFacts();
  }

void SymmetricEngineRing::refreshSingleBasisElementCoefficientFact(
    SymmetricConversionMetadata& metadata,
    const SymmetricRingPoly *poly) const
{
    if (!metadata.expressionFactsComplete) return;
    const bool singleBasisElement =
        metadata.termCount.value_or(0) == 1 &&
        metadata.singleFactorTermCount.value_or(0) == 1 &&
        metadata.scalarTermCount.value_or(0) == 0 &&
        poly->terms.size() == 1 &&
        !poly->terms.front().monomial.data.empty();
    if (!singleBasisElement)
      {
        metadata.singleBasisElementCoefficientOne.reset();
        return;
      }
    metadata.singleBasisElementCoefficientOne =
        coefficientRing->is_equal(
            poly->terms.front().coeff, coefficientRing->one());
  }

SymmetricConversionMetadata SymmetricEngineRing::metadataAfterAddition(
    const SymmetricConversionMetadata& left,
    const SymmetricConversionMetadata& right,
    const SymmetricRingPoly *result) const
{
    SymmetricConversionMetadata metadata;
    if (left.pureBasis == right.pureBasis)
      metadata.pureBasis = left.pureBasis;
    if (left.expandedBasis == right.expandedBasis)
      metadata.expandedBasis = left.expandedBasis;
    if (left.homogeneousWeight == right.homogeneousWeight)
      metadata.homogeneousWeight = left.homogeneousWeight;
    // Cancellation can only decrease this maximum, so the retained value is a
    // conservative cost hint even though the exact core is invalidated below.
    if (left.maximumPartitionLength && right.maximumPartitionLength)
      metadata.maximumPartitionLength = std::max(
          *left.maximumPartitionLength, *right.maximumPartitionLength);
    if (left.factorBases && right.factorBases)
      {
        std::vector<int> factors = *left.factorBases;
        factors.insert(
            factors.end(), right.factorBases->begin(), right.factorBases->end());
        std::sort(factors.begin(), factors.end());
        factors.erase(std::unique(factors.begin(), factors.end()), factors.end());
        metadata.factorBases = std::move(factors);
      }
    metadata.termCount = result->terms.size();
    metadata.normalized = left.normalized && right.normalized;
    metadata.skewFree = left.skewFree && right.skewFree;
    metadata.collected = true;
    // Cancellation can change every exact support count. The retained fields
    // above are hints only and cannot justify the canonical metadata bypass.
    metadata.invalidateExactExpressionFacts();
    return metadata;
  }

SymmetricConversionMetadata SymmetricEngineRing::metadataAfterProduct(
    const SymmetricConversionMetadata& left,
    const SymmetricConversionMetadata& right,
    const SymmetricRingPoly *result) const
{
    SymmetricConversionMetadata metadata;
    if (left.pureBasis && left.pureBasis == right.pureBasis)
      metadata.pureBasis = left.pureBasis;
    if (left.homogeneousWeight && right.homogeneousWeight)
      metadata.homogeneousWeight =
          *left.homogeneousWeight + *right.homogeneousWeight;
    // Product indices concatenate or merge. Neither max(left, right) nor
    // their sum is exact for every basis/storage combination, so leave this
    // optional hint absent until the realized result is inspected.
    if (left.factorBases && right.factorBases)
      {
        std::vector<int> factors = *left.factorBases;
        factors.insert(
            factors.end(), right.factorBases->begin(), right.factorBases->end());
        std::sort(factors.begin(), factors.end());
        factors.erase(std::unique(factors.begin(), factors.end()), factors.end());
        metadata.factorBases = std::move(factors);
      }
    metadata.normalized = left.normalized && right.normalized;
    metadata.skewFree = left.skewFree && right.skewFree;
    if (metadata.pureBasis && isMultiplicativeBasis(*metadata.pureBasis) &&
        metadata.skewFree)
      metadata.expandedBasis = metadata.pureBasis;
    metadata.termCount = result->terms.size();
    metadata.collected = true;
    // General products may merge factors or cancel coefficients. Exact facts
    // are attached only by a workflow that has inspected the realized result.
    metadata.invalidateExactExpressionFacts();
    return metadata;
  }

// ============================================================================
// Metadata Attachment
// ============================================================================
// Facts populate the storage envelope directly. The complete marker is set
// only when the exact canonical bypass core is present.

void SymmetricEngineRing::attachExpressionFacts(
    ring_elem f,
    const ExpressionFacts& facts,
    int targetBasisId,
    CombinatorialTags combinatorialTags) const
{
    if (!facts.canonicalExpansionInBasis(targetBasisId))
      {
        ERROR("only exact canonical target-basis facts can be attached "
              "as public conversion metadata");
        return;
      }
    SymmetricConversionMetadata metadata;
    metadata.pureBasis = facts.pureBasis;
    metadata.expandedBasis = facts.expandedBasis;
    metadata.homogeneousWeight = facts.homogeneousWeight;
    metadata.termCount = facts.termCount;
    metadata.maximumPartitionLength = facts.maximumPartitionLength;
    metadata.factorBases = facts.factorBases;
    metadata.normalized = facts.normalized;
    metadata.skewFree = facts.skewFree;
    metadata.collected = facts.collected;
    metadata.expressionFactsComplete = true;
    metadata.scalarTermCount = facts.scalarTermCount;
    metadata.singleFactorTermCount = facts.singleFactorTermCount;
    metadata.singleBasisElementId = facts.singleBasisElementId;
    metadata.singleBasisElementIndex = facts.singleBasisElementIndex;
    metadata.singleBasisElementCoefficientOne =
        facts.singleBasisElementCoefficientOne;
    auto *poly = mutablePolyValue(f);
    poly->combinatorialTags = combinatorialTags;
    poly->conversionMetadata = std::move(metadata);
  }

// ============================================================================
// Basis-Element Reconstruction
// ============================================================================

ring_elem SymmetricEngineRing::expressionFromBasisElement(
    const SymmetricMonomial& monomial,
    size_t position) const
{
    VECTOR(SymmetricTerm) terms{
        {coefficientRing->one(),
         monomialFromKey(
             atomBlockAt(monomial, position))}};
    return fromTermVector(terms, true);
  }

// ============================================================================
// Basis-Expansion Inspection
// ============================================================================
// Stored expressions are recognized as ordinary basis expansions or scaled basis elements.

bool SymmetricEngineRing::singleScaledBasisElement(
    ring_elem f,
    int basisId,
    Partition& index,
    ring_elem& coefficient) const
{
    if (basisId < 0) return false;
    const auto *poly = polyValue(f);
    if (poly->terms.size() != 1) return false;
    if (!singleBasisIndexFromMonomial(
            poly->terms.front().monomial, basisId, index))
      return false;
    coefficient = poly->terms.front().coeff;
    return isPartitionIndex(index);
  }


CoeffMap SymmetricEngineRing::coefficientsInBasis(ring_elem f, int basisId) const
{
    CoeffMap result;
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (term.monomial.data.empty())
          index = Partition{};
        else if (!singleBasisIndexFromMonomial(term.monomial, basisId, index))
          {
            ERROR("expected an expression in a single symmetric-function basis");
            return CoeffMap{};
          }
        addCoeff(result, index, term.coeff);
      }
    return result;
  }

bool SymmetricEngineRing::coefficientsInBasisIfPossible(
    ring_elem f,
    int basisId,
    CoeffMap& result) const
{
    result.clear();
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (term.monomial.data.empty())
          index = Partition{};
        else if (!singleBasisIndexFromMonomial(term.monomial, basisId, index))
          {
            result.clear();
            return false;
          }
        addCoeff(result, index, term.coeff);
      }
    return true;
  }



// ============================================================================
// Shared Result Construction
// ============================================================================

ring_elem SymmetricEngineRing::expressionHelperFromTerms(
    VECTOR(SymmetricTerm)& terms,
    CombinatorialTags tags) const
{
    ring_elem result = fromTermVector(terms, false);
    mutablePolyValue(result)->combinatorialTags = tags;
    const ExpressionFacts facts = inferExpressionFacts(result);
    attachExpressionHelperFacts(result, facts);
    return result;
  }

void SymmetricEngineRing::attachExpressionHelperFacts(
    ring_elem expression,
    const ExpressionFacts& facts) const
{
    SymmetricConversionMetadata metadata;
    metadata.pureBasis = facts.pureBasis;
    if (facts.expandedBasis &&
        facts.canonicalExpansionInBasis(
            *facts.expandedBasis))
      metadata.expandedBasis = facts.expandedBasis;
    metadata.homogeneousWeight = facts.homogeneousWeight;
    metadata.termCount = facts.termCount;
    metadata.scalarTermCount = facts.scalarTermCount;
    metadata.singleFactorTermCount =
        facts.singleFactorTermCount;
    metadata.maximumPartitionLength =
        facts.maximumPartitionLength;
    metadata.factorBases = facts.factorBases;
    metadata.normalized = facts.normalized;
    metadata.skewFree = facts.skewFree;
    metadata.collected = facts.collected;
    // The existing exact-metadata bypass represents product-free canonical
    // storage. Product-bearing results retain valid step postconditions and
    // support hints, but consumers must still inspect their factor counts.
    metadata.expressionFactsComplete =
        facts.normalized &&
        facts.skewFree &&
        facts.collected &&
        facts.noProducts();
    if (metadata.expressionFactsComplete &&
        facts.singleBasisElement())
      {
        metadata.singleBasisElementId =
            facts.singleBasisElementId;
        metadata.singleBasisElementIndex =
            facts.singleBasisElementIndex;
        metadata.singleBasisElementCoefficientOne =
            facts.singleBasisElementCoefficientOne;
      }
    auto *poly = mutablePolyValue(expression);
    poly->combinatorialTags =
        facts.combinatorialTags;
    poly->conversionMetadata =
        std::move(metadata);
  }

// ============================================================================
// Homogeneous Decomposition
// ============================================================================

std::vector<SymmetricEngineRing::HomogeneousComponent>
SymmetricEngineRing::homogeneousComponents(
    ring_elem expression) const
{
    std::map<int, VECTOR(SymmetricTerm)> termsByWeight;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      termsByWeight[monomialWeight(term.monomial)].push_back(term);

    std::vector<HomogeneousComponent> result;
    result.reserve(termsByWeight.size());
    for (auto& item : termsByWeight)
      result.push_back({
          item.first,
          expressionHelperFromTerms(
              item.second, poly->combinatorialTags)});
    return result;
  }

ring_elem SymmetricEngineRing::homogeneousComponent(
    ring_elem expression,
    int weight) const
{
    VECTOR(SymmetricTerm) selected;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      if (monomialWeight(term.monomial) == weight)
        selected.push_back(term);
    return expressionHelperFromTerms(
        selected, poly->combinatorialTags);
  }

std::vector<int> SymmetricEngineRing::weightSupport(
    ring_elem expression) const
{
    std::vector<int> result;
    const auto components = homogeneousComponents(expression);
    result.reserve(components.size());
    for (const auto& component : components)
      result.push_back(component.weight);
    return result;
  }

ring_elem SymmetricEngineRing::truncateWeights(
    ring_elem expression,
    int minimumWeight,
    int maximumWeight) const
{
    if (minimumWeight > maximumWeight)
      {
        ERROR("expected the minimum weight to be at most the maximum weight");
        return zero();
      }
    VECTOR(SymmetricTerm) selected;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      {
        const int weight = monomialWeight(term.monomial);
        if (minimumWeight <= weight && weight <= maximumWeight)
          selected.push_back(term);
      }
    return expressionHelperFromTerms(
        selected, poly->combinatorialTags);
  }

// ============================================================================
// Normalization And Product Resolution
// ============================================================================

SymmetricEngineRing::ExpressionNormalizationOptions
SymmetricEngineRing::pipelinePreparationNormalizationOptions()
{
    ExpressionNormalizationOptions options;
    options.straightenIndices = true;
    options.expandSkewFactors = true;
    return options;
  }

ring_elem SymmetricEngineRing::normalizeExpressionWithOptions(
    ring_elem expression,
    const ExpressionNormalizationOptions& options) const
{
    return normalizeExpressionWithOptionsDetailed(
        expression, options).expression;
  }

SymmetricEngineRing::ExpressionNormalizationResult
SymmetricEngineRing::normalizeExpressionWithOptionsDetailed(
    ring_elem expression,
    const ExpressionNormalizationOptions& options) const
{
    ExpressionNormalizationResult result;
    const CombinatorialTags inputTags =
        polyValue(expression)->combinatorialTags;
    if (options.productTargetBasisId >= 0)
      {
        (void) requireBasis(options.productTargetBasisId);
        if (error())
          {
            result.expression = zero();
            return result;
          }
        result.expression = expandProductsInBasis(
            expression, options.productTargetBasisId);
        if (error())
          {
            result.expression = zero();
            return result;
          }
        result.resolvedProducts = true;
        auto exactFacts =
            expressionFactsFromMetadata(result.expression);
        if (exactFacts)
          result.facts = std::move(*exactFacts);
        else
          result.facts = inferExpressionFacts(
              result.expression,
              &result.factorsPerTerm);
        attachExpressionHelperFacts(
            result.expression, result.facts);
        return result;
      }

    auto exactFacts =
        expressionFactsFromMetadata(expression);
    if (exactFacts)
      {
        result.expression = copy(expression);
        result.facts = std::move(*exactFacts);
        result.usedCompleteMetadataBypass = true;
        result.bypassedStraighteningFromMetadata =
            options.straightenIndices;
        result.bypassedSkewExpansionFromMetadata =
            options.expandSkewFactors;
        attachExpressionHelperFacts(
            result.expression, result.facts);
        return result;
      }

    const auto *inputPoly = polyValue(expression);
    const bool metadataProvesNormalized =
        inputPoly->conversionMetadata &&
        inputPoly->conversionMetadata->normalized;
    const bool metadataProvesSkewFree =
        inputPoly->conversionMetadata &&
        inputPoly->conversionMetadata->skewFree;
    std::vector<size_t> inputFactorsPerTerm;
    ExpressionFacts inputFacts =
        inferExpressionFacts(
            expression, &inputFactorsPerTerm);
    const bool indicesAlreadyStraightened =
        metadataProvesNormalized ||
        inputFacts.normalized;
    const bool alreadySkewFree =
        metadataProvesSkewFree ||
        inputFacts.skewFree;
    result.bypassedStraighteningFromMetadata =
        options.straightenIndices &&
        metadataProvesNormalized;
    result.bypassedSkewExpansionFromMetadata =
        options.expandSkewFactors &&
        metadataProvesSkewFree;
    if ((!options.straightenIndices ||
         indicesAlreadyStraightened) &&
        (!options.expandSkewFactors ||
         alreadySkewFree))
      {
        result.expression = copy(expression);
        result.facts = std::move(inputFacts);
        result.factorsPerTerm =
            std::move(inputFactorsPerTerm);
        result.usedAlreadyPreparedBypass = true;
        attachExpressionHelperFacts(
            result.expression, result.facts);
        return result;
      }

    result.expression = copy(expression);
    if (options.straightenIndices)
      {
        if (!indicesAlreadyStraightened)
          {
            result.expression =
                straighten(result.expression);
            result.performedStraightening = true;
            if (error())
              {
                result.expression = zero();
                return result;
              }
          }
      }
    if (options.expandSkewFactors)
      {
        if (!alreadySkewFree)
          {
            result.expression =
                expandSkewFactors(result.expression);
            result.performedSkewExpansion = true;
            if (error())
              {
                result.expression = zero();
                return result;
              }
          }
      }
    result.facts = inferExpressionFacts(
        result.expression,
        &result.factorsPerTerm);
    result.facts.combinatorialTags = inputTags;
    mutablePolyValue(
        result.expression)->combinatorialTags =
            inputTags;
    if (options.straightenIndices &&
        !result.facts.normalized)
      {
        ERROR("normalization did not establish its straightened-index contract");
        result.expression = zero();
        return result;
      }
    if (options.expandSkewFactors &&
        !result.facts.skewFree)
      {
        ERROR("normalization did not establish its skew-free contract");
        result.expression = zero();
        return result;
      }
    if (!result.facts.collected)
      {
        ERROR("normalization did not establish its collected-term contract");
        result.expression = zero();
        return result;
      }
    attachExpressionHelperFacts(
        result.expression, result.facts);
    return result;
  }

ring_elem SymmetricEngineRing::expandSkewFactors(
    ring_elem expression) const
{
    VECTOR(SymmetricTerm) expandedTerms;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      {
        ring_elem expanded = fromCoeff(term.coeff);
        size_t position = 0;
        while (position < term.monomial.data.size())
          {
            ring_elem factor =
                atomIsSkewAt(term.monomial, position)
                ? expandSkewBasisElement(
                      term.monomial, position)
                : expressionFromBasisElement(
                      term.monomial, position);
            if (error()) return zero();
            expanded = mult(expanded, factor);
            if (error()) return zero();
            position += atomLengthAt(
                term.monomial, position);
          }
        const auto *expandedPoly = polyValue(expanded);
        expandedTerms.insert(
            expandedTerms.end(),
            expandedPoly->terms.begin(),
            expandedPoly->terms.end());
      }
    ring_elem result = expressionHelperFromTerms(
        expandedTerms, poly->combinatorialTags);
    const ExpressionFacts facts =
        inferExpressionFacts(result);
    if (!facts.skewFree)
      {
        ERROR("skew expansion did not establish its skew-free contract");
        return zero();
      }
    attachExpressionHelperFacts(result, facts);
    return result;
  }

ring_elem SymmetricEngineRing::expandProductsInBasis(
    ring_elem expression,
    int targetBasisId) const
{
    (void) requireBasis(targetBasisId);
    if (error()) return zero();
    const ExpressionNormalizationOptions options =
        pipelinePreparationNormalizationOptions();
    ExpressionNormalizationResult prepared =
        normalizeExpressionWithOptionsDetailed(
            expression, options);
    if (error()) return zero();
    if (prepared.facts.noProducts())
      return prepared.expression;

    const auto& terms =
        polyValue(prepared.expression)->terms;
    if (prepared.factorsPerTerm.size() !=
        terms.size())
      {
        ERROR("product resolution received an inconsistent factor-count profile");
        return zero();
      }

    VECTOR(SymmetricTerm) productFreeTerms;
    for (size_t termPosition = 0;
         termPosition < terms.size();
         ++termPosition)
      {
        const auto& term = terms[termPosition];
        if (prepared.factorsPerTerm[termPosition] <= 1)
          {
            productFreeTerms.push_back(term);
            continue;
          }
        ResolvedProductTerm resolved =
            multiplyTermToBasis(
                term, targetBasisId);
        if (error()) return zero();
        for (const auto& resolvedTerm :
             polyValue(resolved.expression)->terms)
          appendTermIfNonZero(
              productFreeTerms,
              resolvedTerm.coeff,
              resolvedTerm.monomial);
      }

    ring_elem result =
        fromTermVector(productFreeTerms, false);
    auto *resultPoly = mutablePolyValue(result);
    resultPoly->combinatorialTags =
        prepared.facts.combinatorialTags;
    ExpressionFacts facts =
        inferExpressionFacts(result);
    facts.combinatorialTags =
        prepared.facts.combinatorialTags;
    if (!facts.normalized ||
        !facts.skewFree ||
        !facts.collected ||
        !facts.noProducts())
      {
        ERROR("product resolution did not establish its canonical-term contract");
        return zero();
      }
    attachExpressionHelperFacts(result, facts);
    return result;
  }

// ============================================================================
// Stable Expression Inspection
// ============================================================================

SymmetricEngineRing::ExpressionShape
SymmetricEngineRing::expressionShape(
    ring_elem expression) const
{
    const ExpressionFacts facts =
        inferExpressionFacts(expression);
    ExpressionShape result;
    result.weights = weightSupport(expression);
    result.basisIds = facts.factorBases;
    result.termCount = facts.termCount;
    result.scalarTermCount = facts.scalarTermCount;
    result.singleFactorTermCount =
        facts.singleFactorTermCount;
    result.productTermCount = facts.productTermCount;
    result.maximumFactorsPerTerm =
        facts.maximumFactorsPerTerm;
    result.skewFactorCount = facts.skewFactorCount;
    result.maximumPartitionLength =
        facts.maximumPartitionLength;
    result.normalized = facts.normalized;
    result.skewFree = facts.skewFree;
    result.collected = facts.collected;
    result.homogeneousWeight =
        facts.homogeneousWeight;
    result.pureBasis = facts.pureBasis;
    result.expandedBasis = facts.expandedBasis;
    const auto *poly = polyValue(expression);
    if (poly->conversionMetadata)
      {
        const auto& metadata =
            *poly->conversionMetadata;
        result.hasMetadata = true;
        result.metadataFactsComplete =
            metadata.expressionFactsComplete;
        result.metadataNormalized =
            metadata.normalized;
        result.metadataSkewFree =
            metadata.skewFree;
        result.metadataCollected =
            metadata.collected;
      }
    return result;
  }

std::vector<int> SymmetricEngineRing::basisSupport(
    ring_elem expression) const
{
    return inferExpressionFacts(expression).factorBases;
  }

bool SymmetricEngineRing::isBasisExpansion(
    ring_elem expression,
    int basisId) const
{
    (void) requireBasis(basisId);
    if (error()) return false;
    return inferExpressionFacts(expression)
        .canonicalExpansionInBasis(basisId);
  }

bool SymmetricEngineRing::isLinearCombinationOfBasisElements(
    ring_elem expression) const
{
    if (expressionFactsFromMetadata(expression))
      return true;
    const ExpressionFacts facts =
        inferExpressionFacts(expression);
    return facts.normalized &&
           facts.skewFree &&
           facts.collected &&
           facts.noProducts();
  }

ring_elem SymmetricEngineRing::coefficientsInBasisExpression(
    ring_elem expression,
    int basisId,
    bool convert) const
{
    (void) requireBasis(basisId);
    if (error()) return zero();
    if (convert) return toBasis(expression, basisId);
    if (!isBasisExpansion(expression, basisId))
      {
        ERROR("expected a canonical expansion in the requested basis");
        return zero();
      }
    return copy(expression);
  }

// ============================================================================
// Basis And Term Decomposition
// ============================================================================

std::vector<SymmetricEngineRing::BasisComponent>
SymmetricEngineRing::basisComponents(
    ring_elem expression) const
{
    std::map<int, VECTOR(SymmetricTerm)> termsByBasis;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      {
        int basisId = 0;
        if (!term.monomial.data.empty())
          {
            if (atomLengthAt(term.monomial, 0) !=
                term.monomial.data.size())
              {
                ERROR("basisComponents requires a product-free expression");
                return {};
              }
            basisId = atomBasisIdAt(term.monomial, 0);
          }
        termsByBasis[basisId].push_back(term);
      }

    std::vector<BasisComponent> result;
    result.reserve(termsByBasis.size());
    for (auto& item : termsByBasis)
      result.push_back({
          item.first,
          expressionHelperFromTerms(
              item.second, poly->combinatorialTags)});
    return result;
  }

std::vector<SymmetricEngineRing::HomogeneousBasisComponent>
SymmetricEngineRing::homogeneousBasisComponents(
    ring_elem expression) const
{
    std::map<
        std::pair<int, int>,
        VECTOR(SymmetricTerm)> groupedTerms;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      {
        int basisId = 0;
        if (!term.monomial.data.empty())
          {
            if (atomLengthAt(term.monomial, 0) !=
                term.monomial.data.size())
              {
                ERROR(
                    "homogeneousBasisComponents requires "
                    "a product-free expression");
                return {};
              }
            basisId = atomBasisIdAt(term.monomial, 0);
          }
        groupedTerms[{
            monomialWeight(term.monomial),
            basisId}].push_back(term);
      }

    std::vector<HomogeneousBasisComponent> result;
    result.reserve(groupedTerms.size());
    for (auto& item : groupedTerms)
      result.push_back({
          item.first.first,
          item.first.second,
          expressionHelperFromTerms(
              item.second, poly->combinatorialTags)});
    return result;
  }

std::vector<ring_elem> SymmetricEngineRing::singlePartitionIndexedTerms(
    ring_elem expression) const
{
    if (!isLinearCombinationOfBasisElements(expression))
      {
        ERROR("singlePartitionIndexedTerms requires a straightened, skew-free, "
              "product-free expression");
        return {};
      }
    std::vector<ring_elem> result;
    const auto *poly = polyValue(expression);
    result.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        VECTOR(SymmetricTerm) oneTerm{term};
        result.push_back(
            expressionHelperFromTerms(
                oneTerm, poly->combinatorialTags));
      }
    return result;
  }

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
