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

ExpressionFacts
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
    if (facts.normalized &&
        facts.skewFree &&
        facts.collected &&
        facts.productTermCount == 0 &&
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

ExpressionFacts
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
    bool contractViolated = false;

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
            contractViolated = true;
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
    // This fast inspector relies on a canonical target-basis contract. If a
    // caller violates that contract, return the general exact facts rather
    // than manufacturing partially contradictory facts.
    if (contractViolated)
      return inferExpressionFacts(f);

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

ExpressionFacts
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
    if (facts.singleBasisElement())
      {
        facts.expandedBasis = basisId;
        facts.singleBasisElementId = basisId;
        facts.singleBasisElementIndex =
            basisElementIndex(monomial, pos);
        facts.singleBasisElementCoefficientOne =
            coefficientRing->is_equal(
                coefficient, coefficientRing->one());
      }
    return facts;
  }

std::optional<ExpressionFacts>
SymmetricEngineRing::canonicalExpressionFactsFromCache(
    ring_elem f) const
{
    const auto *poly = polyValue(f);
    if (!poly->expressionFactsCache)
      return std::nullopt;
    const auto& cache =
        *poly->expressionFactsCache;
    if (!cache.hasCompleteCanonicalFacts())
      return std::nullopt;

    ExpressionFacts facts = cache.facts;
    facts.combinatorialTags = poly->combinatorialTags;
    if (facts.termCount != poly->terms.size() ||
        facts.scalarTermCount + facts.singleFactorTermCount +
            facts.productTermCount != facts.termCount)
      return std::nullopt;
    if (!facts.normalized ||
        !facts.skewFree ||
        !facts.collected ||
        !facts.noProducts())
      return std::nullopt;
    const size_t expectedMaximumFactors =
        facts.singleFactorTermCount == 0 ? 0 : 1;
    if (facts.maximumFactorsPerTerm != expectedMaximumFactors ||
        facts.skewFactorCount != 0 ||
        !std::is_sorted(
            facts.factorBases.begin(), facts.factorBases.end()) ||
        std::adjacent_find(
            facts.factorBases.begin(), facts.factorBases.end()) !=
                facts.factorBases.end())
      return std::nullopt;
    const std::optional<int> expectedBasis =
        facts.factorBases.size() == 1
            ? std::optional<int>{facts.factorBases.front()}
            : std::nullopt;
    if ((facts.singleFactorTermCount == 0 &&
         !facts.factorBases.empty()) ||
        (facts.singleFactorTermCount != 0 &&
         facts.factorBases.empty()) ||
        facts.pureBasis != expectedBasis ||
        facts.expandedBasis != expectedBasis)
      return std::nullopt;
    const bool singleBasisElement =
        facts.termCount == 1 &&
        facts.scalarTermCount == 0 &&
        facts.singleFactorTermCount == 1;
    if (singleBasisElement !=
            static_cast<bool>(facts.singleBasisElementId) ||
        singleBasisElement !=
            static_cast<bool>(facts.singleBasisElementIndex) ||
        singleBasisElement !=
            static_cast<bool>(
                facts.singleBasisElementCoefficientOne) ||
        (facts.singleBasisElementId &&
         (!facts.expandedBasis ||
          *facts.singleBasisElementId != *facts.expandedBasis)))
      return std::nullopt;
    return facts;
  }

ExpressionFacts SymmetricEngineRing::exactExpressionFacts(
    ring_elem f) const
{
    auto cached =
        canonicalExpressionFactsFromCache(f);
    return cached ? std::move(*cached) : inferExpressionFacts(f);
  }

// ============================================================================
// Expression-Facts Cache Lifecycle
// ============================================================================
// Arithmetic may preserve inexpensive positive facts even when cancellation
// or product formation prevents it from proving the complete canonical
// contract. All such construction and invalidation goes through these helpers
// so no operation invents its own interpretation of the shared cache.

void SymmetricEngineRing::discardSupportDependentFacts(
    ExpressionFactsCacheSlot& cache) const
{
    if (!cache) return;
    cache->discardSupportDependentFacts();
  }

void SymmetricEngineRing::refreshSingleBasisElementCoefficientFact(
    ExpressionFactsCache& cache,
    const SymmetricRingPoly *poly) const
{
    if (!cache.hasCompleteCanonicalFacts()) return;
    const bool singleBasisElement =
        cache.facts.singleBasisElement() &&
        poly->terms.size() == 1 &&
        !poly->terms.front().monomial.data.empty();
    if (!singleBasisElement)
      {
        cache.facts.singleBasisElementCoefficientOne.reset();
        return;
      }
    cache.facts.singleBasisElementCoefficientOne =
        coefficientRing->is_equal(
            poly->terms.front().coeff, coefficientRing->one());
  }

ExpressionFactsCache SymmetricEngineRing::expressionFactsCacheAfterAddition(
    const ExpressionFactsCache& left,
    const ExpressionFactsCache& right) const
{
    ExpressionFactsCache cache;
    if (left.knows(ExpressionFactKnowledge::PureBasis) &&
        right.knows(ExpressionFactKnowledge::PureBasis) &&
        left.facts.pureBasis &&
        left.facts.pureBasis == right.facts.pureBasis)
      {
        cache.facts.pureBasis = left.facts.pureBasis;
        cache.remember(ExpressionFactKnowledge::PureBasis);
      }
    if (left.knows(ExpressionFactKnowledge::ExpandedBasis) &&
        right.knows(ExpressionFactKnowledge::ExpandedBasis) &&
        left.facts.expandedBasis &&
        left.facts.expandedBasis == right.facts.expandedBasis)
      {
        cache.facts.expandedBasis = left.facts.expandedBasis;
        cache.remember(ExpressionFactKnowledge::ExpandedBasis);
      }
    if (left.knows(ExpressionFactKnowledge::HomogeneousWeight) &&
        right.knows(ExpressionFactKnowledge::HomogeneousWeight) &&
        left.facts.homogeneousWeight &&
        left.facts.homogeneousWeight == right.facts.homogeneousWeight)
      {
        cache.facts.homogeneousWeight = left.facts.homogeneousWeight;
        cache.remember(ExpressionFactKnowledge::HomogeneousWeight);
      }
    if (left.knows(ExpressionFactKnowledge::Normalized) &&
        right.knows(ExpressionFactKnowledge::Normalized) &&
        left.facts.normalized &&
        right.facts.normalized)
      {
        cache.facts.normalized = true;
        cache.remember(ExpressionFactKnowledge::Normalized);
      }
    if (left.knows(ExpressionFactKnowledge::SkewFree) &&
        right.knows(ExpressionFactKnowledge::SkewFree) &&
        left.facts.skewFree &&
        right.facts.skewFree)
      {
        cache.facts.skewFree = true;
        cache.remember(ExpressionFactKnowledge::SkewFree);
      }
    cache.facts.collected = true;
    cache.remember(ExpressionFactKnowledge::Collected);
    // Cancellation can change exact support, counts, and negative structural
    // claims. Retain only positive properties shared by both summands.
    return cache;
  }

ExpressionFactsCache SymmetricEngineRing::expressionFactsCacheAfterProduct(
    const ExpressionFactsCache& left,
    const ExpressionFactsCache& right) const
{
    ExpressionFactsCache cache;
    if (left.knows(ExpressionFactKnowledge::PureBasis) &&
        right.knows(ExpressionFactKnowledge::PureBasis) &&
        left.facts.pureBasis &&
        left.facts.pureBasis == right.facts.pureBasis)
      {
        cache.facts.pureBasis = left.facts.pureBasis;
        cache.remember(ExpressionFactKnowledge::PureBasis);
      }
    if (left.knows(ExpressionFactKnowledge::HomogeneousWeight) &&
        right.knows(ExpressionFactKnowledge::HomogeneousWeight) &&
        left.facts.homogeneousWeight &&
        right.facts.homogeneousWeight)
      {
        cache.facts.homogeneousWeight =
            *left.facts.homogeneousWeight +
            *right.facts.homogeneousWeight;
        cache.remember(ExpressionFactKnowledge::HomogeneousWeight);
      }
    // Product indices concatenate or merge, so exact support and maximum
    // partition length remain unknown until the realized result is inspected.
    if (left.knows(ExpressionFactKnowledge::Normalized) &&
        right.knows(ExpressionFactKnowledge::Normalized) &&
        left.facts.normalized &&
        right.facts.normalized)
      {
        cache.facts.normalized = true;
        cache.remember(ExpressionFactKnowledge::Normalized);
      }
    if (left.knows(ExpressionFactKnowledge::SkewFree) &&
        right.knows(ExpressionFactKnowledge::SkewFree) &&
        left.facts.skewFree &&
        right.facts.skewFree)
      {
        cache.facts.skewFree = true;
        cache.remember(ExpressionFactKnowledge::SkewFree);
      }
    if (cache.knows(ExpressionFactKnowledge::PureBasis) &&
        cache.facts.pureBasis &&
        isMultiplicativeBasis(*cache.facts.pureBasis) &&
        cache.knows(ExpressionFactKnowledge::Normalized) &&
        cache.facts.normalized &&
        cache.knows(ExpressionFactKnowledge::SkewFree) &&
        cache.facts.skewFree)
      {
        cache.facts.expandedBasis = cache.facts.pureBasis;
        cache.remember(ExpressionFactKnowledge::ExpandedBasis);
      }
    cache.facts.collected = true;
    cache.remember(ExpressionFactKnowledge::Collected);
    // General products may merge factors or cancel coefficients. Exact facts
    // are attached only by a workflow that has inspected the realized result.
    return cache;
  }

// ============================================================================
// Expression-Facts Cache Attachment
// ============================================================================
// Facts populate the storage envelope directly. The complete canonical marker
// is set only when the full workflow-bypass contract is present.

void SymmetricEngineRing::attachCanonicalExpansionFacts(
    ring_elem f,
    const ExpressionFacts& facts,
    int targetBasisId,
    CombinatorialTags combinatorialTags) const
{
    if (!facts.canonicalExpansionInBasis(targetBasisId))
      {
        ERROR("only exact canonical target-basis facts can be attached "
              "to the expression-facts cache");
        return;
      }
    ExpressionFacts attachedFacts = facts;
    attachedFacts.combinatorialTags = combinatorialTags;
    mutablePolyValue(f)->combinatorialTags = combinatorialTags;
    attachInspectedExpressionFacts(f, attachedFacts);
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

ring_elem SymmetricEngineRing::expressionFromTermsWithFacts(
    VECTOR(SymmetricTerm)& terms,
    CombinatorialTags tags) const
{
    ring_elem result = fromTermVector(terms, false);
    mutablePolyValue(result)->combinatorialTags = tags;
    const ExpressionFacts facts = inferExpressionFacts(result);
    attachInspectedExpressionFacts(result, facts);
    return result;
  }

void SymmetricEngineRing::attachInspectedExpressionFacts(
    ring_elem expression,
    const ExpressionFacts& facts) const
{
    ExpressionFactsCache cache;
    cache.facts = facts;
    cache.remember(ExpressionFactKnowledge::Normalized);
    cache.remember(ExpressionFactKnowledge::SkewFree);
    cache.remember(ExpressionFactKnowledge::Collected);
    cache.remember(ExpressionFactKnowledge::PureBasis);
    cache.remember(ExpressionFactKnowledge::ExpandedBasis);
    cache.remember(ExpressionFactKnowledge::HomogeneousWeight);
    cache.remember(ExpressionFactKnowledge::MaximumPartitionLength);
    cache.remember(ExpressionFactKnowledge::FactorBases);
    // The complete canonical bypass represents product-free canonical storage.
    // Product-bearing results retain exact individual fields and step
    // postconditions, but consumers must still inspect their factor counts.
    if (facts.normalized &&
        facts.skewFree &&
        facts.collected &&
        facts.noProducts())
      cache.setCompleteCanonicalFacts(facts);
    auto *poly = mutablePolyValue(expression);
    poly->combinatorialTags =
        facts.combinatorialTags;
    poly->expressionFactsCache =
        std::move(cache);
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
          expressionFromTermsWithFacts(
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
    return expressionFromTermsWithFacts(
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
    return expressionFromTermsWithFacts(
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
            canonicalExpressionFactsFromCache(
                result.expression);
        if (exactFacts)
          result.facts = std::move(*exactFacts);
        else
          result.facts = inferExpressionFacts(
              result.expression,
              &result.factorsPerTerm);
        attachInspectedExpressionFacts(
            result.expression, result.facts);
        return result;
      }

    auto exactFacts =
        canonicalExpressionFactsFromCache(expression);
    if (exactFacts)
      {
        result.expression = copy(expression);
        result.facts = std::move(*exactFacts);
        result.usedCompleteCanonicalFactsBypass = true;
        result.bypassedStraighteningFromCache =
            options.straightenIndices;
        result.bypassedSkewExpansionFromCache =
            options.expandSkewFactors;
        attachInspectedExpressionFacts(
            result.expression, result.facts);
        return result;
      }

    const auto *inputPoly = polyValue(expression);
    const bool cacheProvesNormalized =
        inputPoly->expressionFactsCache &&
        inputPoly->expressionFactsCache->knows(
            ExpressionFactKnowledge::Normalized) &&
        inputPoly->expressionFactsCache->facts.normalized;
    const bool cacheProvesSkewFree =
        inputPoly->expressionFactsCache &&
        inputPoly->expressionFactsCache->knows(
            ExpressionFactKnowledge::SkewFree) &&
        inputPoly->expressionFactsCache->facts.skewFree;
    std::vector<size_t> inputFactorsPerTerm;
    ExpressionFacts inputFacts =
        inferExpressionFacts(
            expression, &inputFactorsPerTerm);
    const bool indicesAlreadyStraightened =
        cacheProvesNormalized ||
        inputFacts.normalized;
    const bool alreadySkewFree =
        cacheProvesSkewFree ||
        inputFacts.skewFree;
    result.bypassedStraighteningFromCache =
        options.straightenIndices &&
        cacheProvesNormalized;
    result.bypassedSkewExpansionFromCache =
        options.expandSkewFactors &&
        cacheProvesSkewFree;
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
        attachInspectedExpressionFacts(
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
    attachInspectedExpressionFacts(
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
    ring_elem result = expressionFromTermsWithFacts(
        expandedTerms, poly->combinatorialTags);
    const ExpressionFacts facts =
        inferExpressionFacts(result);
    if (!facts.skewFree)
      {
        ERROR("skew expansion did not establish its skew-free contract");
        return zero();
      }
    attachInspectedExpressionFacts(result, facts);
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
    attachInspectedExpressionFacts(result, facts);
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
        exactExpressionFacts(expression);
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
    if (poly->expressionFactsCache)
      {
        const auto& cache =
            *poly->expressionFactsCache;
        result.hasMetadata = true;
        result.metadataFactsComplete =
            cache.hasCompleteCanonicalFacts() &&
            canonicalExpressionFactsFromCache(
                expression).has_value();
        result.metadataNormalized =
            cache.knows(ExpressionFactKnowledge::Normalized) &&
            cache.facts.normalized;
        result.metadataSkewFree =
            cache.knows(ExpressionFactKnowledge::SkewFree) &&
            cache.facts.skewFree;
        result.metadataCollected =
            cache.knows(ExpressionFactKnowledge::Collected) &&
            cache.facts.collected;
      }
    return result;
  }

std::vector<int> SymmetricEngineRing::basisSupport(
    ring_elem expression) const
{
    return exactExpressionFacts(expression).factorBases;
  }

bool SymmetricEngineRing::isBasisExpansion(
    ring_elem expression,
    int basisId) const
{
    (void) requireBasis(basisId);
    if (error()) return false;
    return exactExpressionFacts(expression)
        .canonicalExpansionInBasis(basisId);
  }

bool SymmetricEngineRing::isLinearCombinationOfBasisElements(
    ring_elem expression) const
{
    const ExpressionFacts facts =
        exactExpressionFacts(expression);
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
          expressionFromTermsWithFacts(
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
          expressionFromTermsWithFacts(
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
            expressionFromTermsWithFacts(
                oneTerm, poly->combinatorialTags));
      }
    return result;
  }

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
