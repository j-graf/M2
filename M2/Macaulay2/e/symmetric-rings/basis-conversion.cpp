// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "basic-rings/aring-glue.hpp"
#include "error.h"
#include "symmetric-rings/basis-conversion-policy.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <utility>

namespace symmetric_rings {

// Contributor map:
//   1. exact expression facts and persistent metadata;
//   2. normalization and condition-piece contexts;
//   3. declarative conversion plans, contracts, selection, and execution;
//   4. canonical source-group conversion;
//   5. multiplication plans and the three multiplication workflows;
//   6. the sole public built-in conversion workflow, toBasis.
//
// The matching declaration fragment follows the same order. Mathematical
// formulas live in basis-conversion-kernels.cpp and product combinatorics live
// in basis-conversion-products.cpp; this file owns orchestration, not formulas.

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
            facts.factorKindMask |=
                uint32_t{1} << static_cast<size_t>(
                    basisKindForId(basisId));
            ++facts.factorCount;

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
    // A registered conversion or multiplication plan promises one normalized,
    // non-skew target atom per nonscalar term. Derive the common output facts
    // from that contract and inspect only support that cannot be predicted,
    // rather than running the general multi-factor analyzer.
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
        facts.factorKindMask =
            uint32_t{1} << static_cast<size_t>(
                basisKindForId(basisId));
        facts.factorCount = facts.singleFactorTermCount;
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
    facts.factorKindMask =
        uint32_t{1} << static_cast<size_t>(
            basisKindForId(basisId));
    facts.factorCount = 1;
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

void SymmetricEngineRing::enrichBasisConversionPlanSelectionFacts(
    ring_elem f,
    int targetBasisId,
    ExpressionFacts& facts) const
{
    if (facts.planSelectionTargetBasisId == targetBasisId)
      return;
    if (!facts.expandedBasis ||
        basisKindForId(*facts.expandedBasis) != BasisKind::PowerSum)
      return;
    BasisKind targetKind = basisKindForId(targetBasisId);
    const bool schurTarget =
        targetKind == BasisKind::Schur ||
        targetKind == BasisKind::SchurOmega;
    const bool hallLittlewoodTarget =
        isHallLittlewoodCapitalBasisKind(targetKind);
    if (!schurTarget && !hallLittlewoodTarget) return;
    facts.planSelectionTargetBasisId = targetBasisId;

    if (facts.singleBasisElementIndex &&
        facts.singleBasisElement())
      {
        // Strict basis-element entry points already decoded the one index.
        // Derive every relevant p profile from that cached shape instead of
        // reparsing the expression immediately afterward.
        const Partition& index = *facts.singleBasisElementIndex;
        if (hallLittlewoodTarget)
          {
            facts.allPowerSumTermsSingleCycles =
                index.size() == 1;
          }
        if (schurTarget)
          {
            facts.completeFriendlyPowerSumTermCount =
                completeFriendlyPowerSumIndexForConversionSelection(
                    index, partitionWeight(index))
                    ? 1
                    : 0;
            std::vector<int> commonParts(index.begin(), index.end());
            std::sort(commonParts.begin(), commonParts.end());
            commonParts.erase(
                std::unique(commonParts.begin(), commonParts.end()),
                commonParts.end());
            facts.commonPowerSumParts = std::move(commonParts);
          }
        return;
      }

    bool allSingleCycles = true;
    size_t nonscalarTerms = 0;
    size_t completeFriendlyTerms = 0;
    std::vector<int> commonParts;
    bool firstIndex = true;
    for (const auto& term : polyValue(f)->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            allSingleCycles = false;
            continue;
          }
        if (index.empty()) continue;
        ++nonscalarTerms;
        if (index.size() != 1)
          allSingleCycles = false;
        if (!schurTarget) continue;
        if (completeFriendlyPowerSumIndexForConversionSelection(
                index, partitionWeight(index)))
          ++completeFriendlyTerms;

        std::vector<int> distinctParts(index.begin(), index.end());
        std::sort(distinctParts.begin(), distinctParts.end());
        distinctParts.erase(
            std::unique(distinctParts.begin(), distinctParts.end()),
            distinctParts.end());
        if (firstIndex)
          {
            commonParts = std::move(distinctParts);
            firstIndex = false;
          }
        else
          {
            std::vector<int> intersection;
            std::set_intersection(
                commonParts.begin(),
                commonParts.end(),
                distinctParts.begin(),
                distinctParts.end(),
                std::back_inserter(intersection));
            commonParts = std::move(intersection);
          }
      }
    if (nonscalarTerms == 0) return;
    if (hallLittlewoodTarget)
      facts.allPowerSumTermsSingleCycles = allSingleCycles;
    if (schurTarget)
      {
        facts.completeFriendlyPowerSumTermCount =
            completeFriendlyTerms;
        facts.commonPowerSumParts = std::move(commonParts);
      }
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
        !metadata.singleFactorTermCount ||
        !metadata.productTermCount ||
        !metadata.maximumFactorsPerTerm)
      return std::nullopt;
    facts.scalarTermCount = *metadata.scalarTermCount;
    facts.singleFactorTermCount = *metadata.singleFactorTermCount;
    facts.productTermCount = *metadata.productTermCount;
    facts.maximumFactorsPerTerm = *metadata.maximumFactorsPerTerm;
    if (facts.termCount != poly->terms.size() ||
        facts.scalarTermCount + facts.singleFactorTermCount +
            facts.productTermCount != facts.termCount ||
        facts.productTermCount != 0 ||
        facts.maximumFactorsPerTerm !=
            (facts.singleFactorTermCount == 0 ? 0 : 1) ||
        metadata.skewFactorCount.value_or(0) != 0)
      return std::nullopt;
    // A zero or scalar expansion has no distinguished basis. Every
    // nonscalar canonical expansion must name its unique source basis.
    if ((facts.singleFactorTermCount != 0 && !metadata.expandedBasis) ||
        (facts.singleFactorTermCount == 0 && metadata.expandedBasis))
      return std::nullopt;
    facts.pureBasis = metadata.pureBasis
        ? metadata.pureBasis : metadata.expandedBasis;
    facts.expandedBasis = metadata.expandedBasis;
    if (metadata.factorBases)
      facts.factorBases = *metadata.factorBases;
    else if (metadata.expandedBasis)
      facts.factorBases = {*metadata.expandedBasis};
    if (metadata.expandedBasis &&
        (facts.factorBases.size() != 1 ||
         facts.factorBases.front() != *metadata.expandedBasis))
      return std::nullopt;
    if (metadata.expandedBasis &&
        facts.singleFactorTermCount != 0)
      {
        facts.factorKindMask =
            uint32_t{1} << static_cast<size_t>(
                basisKindForId(*metadata.expandedBasis));
        facts.factorCount = facts.singleFactorTermCount;
      }
    facts.skewFactorCount = metadata.skewFactorCount.value_or(0);
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
        metadata.productTermCount.value_or(0) == 0 &&
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
    metadata.productTermCount = facts.productTermCount;
    metadata.maximumFactorsPerTerm = facts.maximumFactorsPerTerm;
    metadata.skewFactorCount = facts.skewFactorCount;
    metadata.singleBasisElementId = facts.singleBasisElementId;
    metadata.singleBasisElementIndex = facts.singleBasisElementIndex;
    metadata.singleBasisElementCoefficientOne =
        facts.singleBasisElementCoefficientOne;
    auto *poly = mutablePolyValue(f);
    poly->combinatorialTags = combinatorialTags;
    poly->conversionMetadata = std::move(metadata);
  }

// ============================================================================
// Canonicalization And Skew Expansion
// ============================================================================
// Normalization straightens every index and replaces every skew atom by a
// mathematically equivalent non-skew expansion. Products may remain; resolving
// them belongs to multiplyTermToBasis, the next owning workflow stage.

ring_elem SymmetricEngineRing::expressionFromAtom(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    VECTOR(SymmetricTerm) terms{
        {coefficientRing->one(),
         monomialFromKey(atomBlockAt(monomial, pos))}};
    return fromTermVector(terms, true);
  }

ring_elem SymmetricEngineRing::skewBasisElementExpansion(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    if (!atomIsSkewAt(monomial, pos))
      {
        ERROR("the skew-expansion helper requires a skew basis element");
        return zero();
      }

    BasisKind kind = basisKindForId(atomBasisIdAt(monomial, pos));
    Partition outer = basisElementOuterIndex(monomial, pos);
    Partition inner = basisElementInnerIndex(monomial, pos);
    if (kind == BasisKind::Schur ||
        kind == BasisKind::SchurOmega)
      {
        // Jacobi--Trudi naturally returns products in h or e. Leaving those
        // products visible lets the ordinary product-resolution stage choose
        // the target-aware multiplication workflow.
        BasisKind generatorKind =
            kind == BasisKind::Schur
                ? BasisKind::Complete
                : BasisKind::Elementary;
        int generatorId = requiredBasisIdForKind(generatorKind);
        if (error()) return zero();
        return jacobiTrudi(outer, inner, generatorId);
      }
    if (isHallLittlewoodCapitalBasisKind(kind))
      return skewHallLittlewoodToPowerSums(outer, inner, kind);

    ERROR("the normalization stage cannot expand this skew basis kind");
    return zero();
  }

ring_elem SymmetricEngineRing::normalizeExpression(
    ring_elem f,
    ExpressionFacts& resultFacts,
    std::vector<size_t> *termFactorCounts) const
{
    // Canonical engine expressions are immutable.  Inspect once and retain the
    // original value when no straightening or skew expansion is required,
    // rather than copying, sorting, and then inspecting the same terms again.
    resultFacts = inferExpressionFacts(f, termFactorCounts);
    if (resultFacts.normalized &&
        resultFacts.skewFree &&
        resultFacts.collected)
      return f;

    VECTOR(SymmetricTerm) normalizedTerms;
    for (const auto& term : polyValue(f)->terms)
      {
        const ring_elem straightened =
            scaled(term.coeff, straightenMonomial(term.monomial));
        if (error()) return zero();
        for (const auto& straightenedTerm : polyValue(straightened)->terms)
          {
            bool hasSkewFactor = false;
            size_t position = 0;
            while (position < straightenedTerm.monomial.data.size())
              {
                hasSkewFactor =
                    hasSkewFactor ||
                    atomIsSkewAt(straightenedTerm.monomial, position);
                position +=
                    atomLengthAt(straightenedTerm.monomial, position);
              }
            if (!hasSkewFactor)
              {
                normalizedTerms.push_back(straightenedTerm);
                continue;
              }

            ring_elem expanded = fromCoeff(straightenedTerm.coeff);
            position = 0;
            while (position < straightenedTerm.monomial.data.size())
              {
                ring_elem factor;
                if (!atomIsSkewAt(straightenedTerm.monomial, position))
                  factor = expressionFromAtom(
                      straightenedTerm.monomial, position);
                else
                  factor = skewBasisElementExpansion(
                      straightenedTerm.monomial, position);
                if (error()) return zero();
                expanded = mult(expanded, factor);
                if (error()) return zero();
                position +=
                    atomLengthAt(straightenedTerm.monomial, position);
              }
            const auto *expandedPoly = polyValue(expanded);
            normalizedTerms.insert(
                normalizedTerms.end(),
                expandedPoly->terms.begin(),
                expandedPoly->terms.end());
          }
      }
    ring_elem result = fromTermVector(normalizedTerms, false);
    mutablePolyValue(result)->combinatorialTags =
        polyValue(f)->combinatorialTags;
    resultFacts = inferExpressionFacts(result, termFactorCounts);
    if (!resultFacts.normalized ||
        !resultFacts.skewFree ||
        !resultFacts.collected)
      {
        ERROR("the normalization stage did not establish its "
              "canonical-factor contract");
        return zero();
      }
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr)
      std::fprintf(
          stderr,
          "SymmetricRings conversion-stage: stage=normalize "
          "input-terms=%zu output-terms=%zu products=%zu\n",
          polyValue(f)->terms.size(),
          resultFacts.termCount,
          resultFacts.productTermCount);
    return result;
  }

// ============================================================================
// Expression-Condition Context Construction
// ============================================================================
// Conditions remain independent inspectable values. This adapter supplies
// their exact mathematical context from canonical symmetric-function storage.

std::vector<ExpressionConditionContext>
SymmetricEngineRing::buildExpressionConditionContexts(
    ring_elem expression,
    ExpressionPieceKind pieces,
    const ExpressionFacts& facts) const
{
    const auto& terms = polyValue(expression)->terms;
    auto commonContext = [&](ExpressionPieceKind kind) {
      ExpressionConditionContext context;
      context.pieceKind = kind;
      context.combinatorialTags = facts.combinatorialTags;
      context.coefficientRingIsQQ = coefficientRing == globalQQ;
      return context;
    };

    if (pieces == ExpressionPieceKind::WholeExpression)
      {
        ExpressionConditionContext context =
            commonContext(pieces);
        context.termPositions.resize(terms.size());
        std::iota(
            context.termPositions.begin(),
            context.termPositions.end(),
            size_t{0});
        context.weight = facts.homogeneousWeight;
        context.termCount = facts.termCount;
        context.singleBasisElement = facts.singleBasisElement();
        context.allPowerSumTermsSingleCycles =
            facts.allPowerSumTermsSingleCycles;
        context.completeFriendlyPowerSumTermCount =
            facts.completeFriendlyPowerSumTermCount;
        context.index = facts.singleBasisElementIndex;
        std::optional<int> firstWeight;
        bool multipleWeights = false;
        for (const auto& term : terms)
          {
            const int weight = monomialWeight(term.monomial);
            if (!firstWeight)
              firstWeight = weight;
            else if (*firstWeight != weight)
              {
                multipleWeights = true;
                break;
              }
          }
        context.expressionHasMultipleWeights = multipleWeights;
        return {std::move(context)};
      }

    if (pieces == ExpressionPieceKind::Terms)
      {
        std::vector<ExpressionConditionContext> contexts;
        contexts.reserve(terms.size());
        for (size_t position = 0; position < terms.size(); ++position)
          {
            const auto& term = terms[position];
            ExpressionConditionContext context =
                commonContext(pieces);
            context.termPositions = {position};
            context.weight = monomialWeight(term.monomial);
            context.termCount = 1;
            if (!term.monomial.data.empty() &&
                atomLengthAt(term.monomial, 0) ==
                    term.monomial.data.size() &&
                !atomIsSkewAt(term.monomial, 0))
              {
                Partition index =
                    basisElementIndex(term.monomial, 0);
                context.index = index;
                if (basisKindForId(
                        atomBasisIdAt(term.monomial, 0)) ==
                    BasisKind::PowerSum)
                  context.powerSumIndexCompleteFriendly =
                      completeFriendlyPowerSumIndexForConversionSelection(
                          index, partitionWeight(index));
              }
            contexts.push_back(std::move(context));
          }
        return contexts;
      }

    std::map<int, std::vector<size_t>> positionsByWeight;
    for (size_t position = 0; position < terms.size(); ++position)
      positionsByWeight[monomialWeight(
          terms[position].monomial)].push_back(position);
    std::vector<ExpressionConditionContext> contexts;
    contexts.reserve(positionsByWeight.size());
    for (auto& item : positionsByWeight)
      {
        ExpressionConditionContext context =
            commonContext(pieces);
        context.termPositions = std::move(item.second);
        context.weight = item.first;
        context.termCount = context.termPositions.size();
        context.possibleTermCount =
            partitionCountForConversionSelection(item.first);
        size_t completeFriendlyTerms = 0;
        std::vector<int> commonParts;
        bool firstPowerSumIndex = true;
        bool allTermsArePowerSums = true;
        for (size_t position : context.termPositions)
          {
            Partition index;
            if (!powerSumIndexFromMonomial(
                    terms[position].monomial, index))
              {
                allTermsArePowerSums = false;
                break;
              }
            if (index.empty()) continue;
            if (completeFriendlyPowerSumIndexForConversionSelection(
                    index, partitionWeight(index)))
              ++completeFriendlyTerms;
            std::vector<int> distinctParts(
                index.begin(), index.end());
            std::sort(
                distinctParts.begin(), distinctParts.end());
            distinctParts.erase(
                std::unique(
                    distinctParts.begin(), distinctParts.end()),
                distinctParts.end());
            if (firstPowerSumIndex)
              {
                commonParts = std::move(distinctParts);
                firstPowerSumIndex = false;
              }
            else
              {
                std::vector<int> intersection;
                std::set_intersection(
                    commonParts.begin(),
                    commonParts.end(),
                    distinctParts.begin(),
                    distinctParts.end(),
                    std::back_inserter(intersection));
                commonParts = std::move(intersection);
              }
          }
        if (allTermsArePowerSums)
          {
            context.completeFriendlyPowerSumTermCount =
                completeFriendlyTerms;
            context.commonPowerSumParts = std::move(commonParts);
          }
        contexts.push_back(std::move(context));
      }
    return contexts;
  }

// ============================================================================
// Declarative Basis-Conversion Plan Database
// ============================================================================
// This is the single contributor-facing catalog for atomic, delegated,
// composed, and piecewise conversion formulas. Every entry is consumed
// directly by the generic executor.

SymmetricEngineRing::ConversionFormula
SymmetricEngineRing::kernelFormula(BasisConversionKernel kernel)
{
    ConversionFormula result;
    result.kind = ConversionFormulaKind::AtomicKernel;
    result.kernel = kernel;
    return result;
  }

SymmetricEngineRing::ConversionFormula
SymmetricEngineRing::planFormula(
    BasisConversionPlanId childPlan)
{
    ConversionFormula result;
    result.kind = ConversionFormulaKind::PlanComposition;
    result.childPlans.push_back(std::move(childPlan));
    return result;
  }

SymmetricEngineRing::ConversionFormula
SymmetricEngineRing::compositionFormula(
    std::initializer_list<BasisConversionPlanId> childPlans)
{
    if (childPlans.size() < 2)
      throw exc::engine_error(
          "a basis-conversion plan composition requires at least two "
          "named child plans");
    ConversionFormula result;
    result.kind = ConversionFormulaKind::PlanComposition;
    result.childPlans.assign(childPlans.begin(), childPlans.end());
    return result;
  }

const std::vector<SymmetricEngineRing::BasisConversionPlanDefinition>&
SymmetricEngineRing::basisConversionPlanDatabase()
{
    static const std::vector<BasisConversionPlanDefinition> database = [] {
      std::vector<BasisConversionPlanDefinition> plans;
      auto addAtomic =
          [&](std::string_view id,
              BasisKind source,
              BasisKind target,
              BasisConversionKernel kernel,
              ExpressionCondition applicability = always()) {
            plans.push_back({
                {std::string(id)},
                source,
                target,
                std::move(applicability),
                ExpressionPieceKind::WholeExpression,
                {{otherwise(), kernelFormula(kernel)}}});
          };

      // Built-in atomic formulas to power sums.
      addAtomic(
          "Complete->PowerSum:classical-formula",
          BasisKind::Complete,
          BasisKind::PowerSum,
          BasisConversionKernel::CompleteToPowerSumsClassical);
      addAtomic(
          "Elementary->PowerSum:classical-formula",
          BasisKind::Elementary,
          BasisKind::PowerSum,
          BasisConversionKernel::ElementaryToPowerSumsClassical);
      addAtomic(
          "Schur->PowerSum:characters",
          BasisKind::Schur,
          BasisKind::PowerSum,
          BasisConversionKernel::SchurToPowerSumsCharacters);
      addAtomic(
          "SchurOmega->PowerSum:characters",
          BasisKind::SchurOmega,
          BasisKind::PowerSum,
          BasisConversionKernel::SchurOmegaToPowerSumsCharacters);
      addAtomic(
          "Monomial->PowerSum:transition",
          BasisKind::Monomial,
          BasisKind::PowerSum,
          BasisConversionKernel::MonomialToPowerSumsTransition);
      addAtomic(
          "Forgotten->PowerSum:transition",
          BasisKind::Forgotten,
          BasisKind::PowerSum,
          BasisConversionKernel::ForgottenToPowerSumsTransition);
      addAtomic(
          "HallLittlewoodQGenerator->PowerSum:classical-formula",
          BasisKind::HallLittlewoodQGenerator,
          BasisKind::PowerSum,
          BasisConversionKernel::
              HallLittlewoodQGeneratorToPowerSumsClassical);
      addAtomic(
          "HallLittlewoodBGenerator->PowerSum:classical-formula",
          BasisKind::HallLittlewoodBGenerator,
          BasisKind::PowerSum,
          BasisConversionKernel::
              HallLittlewoodBGeneratorToPowerSumsClassical);
      addAtomic(
          "HallLittlewoodQ->PowerSum:raising-operators",
          BasisKind::HallLittlewoodQ,
          BasisKind::PowerSum,
          BasisConversionKernel::HallLittlewoodQToPowerSumsRaising);
      addAtomic(
          "HallLittlewoodB->PowerSum:raising-operators",
          BasisKind::HallLittlewoodB,
          BasisKind::PowerSum,
          BasisConversionKernel::HallLittlewoodBToPowerSumsRaising);
      addAtomic(
          "HallLittlewoodP->PowerSum:capital-normalization",
          BasisKind::HallLittlewoodP,
          BasisKind::PowerSum,
          BasisConversionKernel::
              HallLittlewoodPToPowerSumsNormalization);
      addAtomic(
          "HallLittlewoodPOmega->PowerSum:capital-normalization",
          BasisKind::HallLittlewoodPOmega,
          BasisKind::PowerSum,
          BasisConversionKernel::
              HallLittlewoodPOmegaToPowerSumsNormalization);

      // Atomic formulas from power sums.
      addAtomic(
          "PowerSum->Complete:logarithm-formula",
          BasisKind::PowerSum,
          BasisKind::Complete,
          BasisConversionKernel::PowerSumsToCompleteLogarithm);
      addAtomic(
          "PowerSum->Elementary:logarithm-formula",
          BasisKind::PowerSum,
          BasisKind::Elementary,
          BasisConversionKernel::PowerSumsToElementaryLogarithm);
      addAtomic(
          "PowerSum->Schur:abacus-rim-hooks",
          BasisKind::PowerSum,
          BasisKind::Schur,
          BasisConversionKernel::PowerSumsToSchurAbacusRimHooks);
      addAtomic(
          "PowerSum->Schur:grouped-characters",
          BasisKind::PowerSum,
          BasisKind::Schur,
          BasisConversionKernel::PowerSumsToSchurCharacters);
      addAtomic(
          "PowerSum->Schur:border-strips",
          BasisKind::PowerSum,
          BasisKind::Schur,
          BasisConversionKernel::PowerSumsToSchurBorderStrips);
      addAtomic(
          "PowerSum->SchurOmega:abacus-rim-hooks",
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          BasisConversionKernel::
              PowerSumsToSchurOmegaAbacusRimHooks);
      addAtomic(
          "PowerSum->SchurOmega:grouped-characters",
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          BasisConversionKernel::PowerSumsToSchurOmegaCharacters);
      addAtomic(
          "PowerSum->SchurOmega:border-strips",
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          BasisConversionKernel::PowerSumsToSchurOmegaBorderStrips);
      addAtomic(
          "PowerSum->HallLittlewoodQGenerator:logarithm-formula",
          BasisKind::PowerSum,
          BasisKind::HallLittlewoodQGenerator,
          BasisConversionKernel::
              PowerSumsToHallLittlewoodQGeneratorLogarithm);
      addAtomic(
          "PowerSum->HallLittlewoodBGenerator:logarithm-formula",
          BasisKind::PowerSum,
          BasisKind::HallLittlewoodBGenerator,
          BasisConversionKernel::
              PowerSumsToHallLittlewoodBGeneratorLogarithm);

      auto addPowerSumsToHallLittlewood =
          [&](BasisKind target, std::string_view targetName) {
            const std::string prefix =
                "PowerSum->" + std::string(targetName);
            addAtomic(
                std::string(prefix + ":single-cycles-Green"),
                BasisKind::PowerSum,
                target,
                BasisConversionKernel::
                    PowerSumSingleCyclesToHallLittlewoodGreen,
                allPowerSumTermsAreSingleCycles());
            addAtomic(
                std::string(prefix + ":triangular-reduction"),
                BasisKind::PowerSum,
                target,
                BasisConversionKernel::
                    PowerSumsToHallLittlewoodTriangular);
            addAtomic(
                std::string(prefix + ":Green-duality"),
                BasisKind::PowerSum,
                target,
                BasisConversionKernel::
                    PowerSumIndexToHallLittlewoodGreenDuality,
                expressionIsSingleBasisElement());
          };
      addPowerSumsToHallLittlewood(
          BasisKind::HallLittlewoodQ, "HallLittlewoodQ");
      addPowerSumsToHallLittlewood(
          BasisKind::HallLittlewoodB, "HallLittlewoodB");
      addPowerSumsToHallLittlewood(
          BasisKind::HallLittlewoodP, "HallLittlewoodP");
      addPowerSumsToHallLittlewood(
          BasisKind::HallLittlewoodPOmega,
          "HallLittlewoodPOmega");
      addAtomic(
          "PowerSum->Monomial:transition",
          BasisKind::PowerSum,
          BasisKind::Monomial,
          BasisConversionKernel::PowerSumsToMonomialTransition);
      addAtomic(
          "PowerSum->Forgotten:transition",
          BasisKind::PowerSum,
          BasisKind::Forgotten,
          BasisConversionKernel::PowerSumsToForgottenTransition);

      // Direct normalization, involution, and generator transitions.
      addAtomic(
          "HallLittlewoodQ->HallLittlewoodP:diagonal-scaling",
          BasisKind::HallLittlewoodQ,
          BasisKind::HallLittlewoodP,
          BasisConversionKernel::HallLittlewoodNormalization);
      addAtomic(
          "HallLittlewoodP->HallLittlewoodQ:diagonal-scaling",
          BasisKind::HallLittlewoodP,
          BasisKind::HallLittlewoodQ,
          BasisConversionKernel::HallLittlewoodNormalization);
      addAtomic(
          "HallLittlewoodB->HallLittlewoodPOmega:diagonal-scaling",
          BasisKind::HallLittlewoodB,
          BasisKind::HallLittlewoodPOmega,
          BasisConversionKernel::HallLittlewoodNormalization);
      addAtomic(
          "HallLittlewoodPOmega->HallLittlewoodB:diagonal-scaling",
          BasisKind::HallLittlewoodPOmega,
          BasisKind::HallLittlewoodB,
          BasisConversionKernel::HallLittlewoodNormalization);
      addAtomic(
          "Schur->SchurOmega:partition-conjugation",
          BasisKind::Schur,
          BasisKind::SchurOmega,
          BasisConversionKernel::SchurOmegaConjugation);
      addAtomic(
          "SchurOmega->Schur:partition-conjugation",
          BasisKind::SchurOmega,
          BasisKind::Schur,
          BasisConversionKernel::SchurOmegaConjugation);
      addAtomic(
          "Schur->Complete:Jacobi-Trudi",
          BasisKind::Schur,
          BasisKind::Complete,
          BasisConversionKernel::SchurToCompleteJacobiTrudi);
      addAtomic(
          "SchurOmega->Elementary:Jacobi-Trudi",
          BasisKind::SchurOmega,
          BasisKind::Elementary,
          BasisConversionKernel::
              SchurOmegaToElementaryJacobiTrudi);
      addAtomic(
          "Complete->Schur:recursive-transition",
          BasisKind::Complete,
          BasisKind::Schur,
          BasisConversionKernel::CompleteToSchurRecursive);
      addAtomic(
          "HallLittlewoodQGenerator->HallLittlewoodQ:"
          "triangular-reduction",
          BasisKind::HallLittlewoodQGenerator,
          BasisKind::HallLittlewoodQ,
          BasisConversionKernel::
              HallLittlewoodGeneratorToCapitalTriangular);
      addAtomic(
          "HallLittlewoodQGenerator->HallLittlewoodP:"
          "triangular-reduction",
          BasisKind::HallLittlewoodQGenerator,
          BasisKind::HallLittlewoodP,
          BasisConversionKernel::
              HallLittlewoodGeneratorToCapitalTriangular);
      addAtomic(
          "HallLittlewoodBGenerator->HallLittlewoodB:"
          "triangular-reduction",
          BasisKind::HallLittlewoodBGenerator,
          BasisKind::HallLittlewoodB,
          BasisConversionKernel::
              HallLittlewoodGeneratorToCapitalTriangular);
      addAtomic(
          "HallLittlewoodBGenerator->HallLittlewoodPOmega:"
          "triangular-reduction",
          BasisKind::HallLittlewoodBGenerator,
          BasisKind::HallLittlewoodPOmega,
          BasisConversionKernel::
              HallLittlewoodGeneratorToCapitalTriangular);

      // The established p -> Schur policy applies independently to each
      // homogeneous component of a nonhomogeneous expression. Encoding that
      // decision tree as ordered cases makes the split an ordinary plan:
      // each case fixes its kernel, and execution performs no child choice.
      auto schurComponentPolicyCases = [&] {
            const auto completeFormula =
                compositionFormula({
                    {"PowerSum->Complete:logarithm-formula"},
                    {"Complete->Schur:recursive-transition"}});
            const auto abacusFormula =
                kernelFormula(
                    BasisConversionKernel::
                        PowerSumsToSchurAbacusRimHooks);
            const auto termHybridFormula =
                planFormula({
                    "PowerSum->Schur:"
                    "complete-friendly-hybrid-plan"});
            const uint32_t plethysmTag =
                combinatorialTagMask(
                    CombinatorialTag::Plethysm);
            const uint32_t littlewoodRichardsonTag =
                combinatorialTagMask(
                    CombinatorialTag::
                        LittlewoodRichardson);
            const uint32_t horizontalPieriTag =
                combinatorialTagMask(
                    CombinatorialTag::HorizontalPieri);
            const uint32_t verticalPieriTag =
                combinatorialTagMask(
                    CombinatorialTag::VerticalPieri);
            const uint32_t borderStripsTag =
                combinatorialTagMask(
                    CombinatorialTag::BorderStrips);

            const auto supportPreference =
                componentSupportSquareFavorsComplete();
            const auto ordinaryCompletePreference =
                !coefficientRingIsQQ() ||
                componentDensityAtLeast(1, 4) ||
                supportPreference;
            const auto plethysmComplete =
                hasCombinatorialTag(plethysmTag) &&
                (supportPreference ||
                 (componentWeightGreaterThan(7) &&
                  componentTermCountAtLeast(2)));
            const auto littlewoodRichardsonComplete =
                !hasCombinatorialTag(plethysmTag) &&
                hasCombinatorialTag(
                    littlewoodRichardsonTag) &&
                componentWeightGreaterThan(6) &&
                ordinaryCompletePreference;
            const auto pieriComplete =
                !hasCombinatorialTag(plethysmTag) &&
                !hasCombinatorialTag(
                    littlewoodRichardsonTag) &&
                (hasCombinatorialTag(horizontalPieriTag) ||
                 hasCombinatorialTag(verticalPieriTag)) &&
                componentWeightGreaterThan(7) &&
                ordinaryCompletePreference;
            const auto smallCommonCycle =
                (coefficientRingIsQQ() &&
                 componentHasCommonPowerSumPartAtMostPercent(36)) ||
                (!coefficientRingIsQQ() &&
                 componentHasCommonPowerSumPartAtMostPercent(41));
            const auto largeBorderStripWithoutCommonOne =
                combinatorialTagsEqual(borderStripsTag) &&
                componentWeightGreaterThan(15) &&
                componentTermCountAtLeast(8) &&
                !componentHasCommonPowerSumPartOne() &&
                smallCommonCycle;
            const auto borderStripWithCommonOne =
                combinatorialTagsEqual(borderStripsTag) &&
                componentWeightGreaterThan(7) &&
                componentTermCountAtLeast(8) &&
                componentHasCommonPowerSumPartOne() &&
                ordinaryCompletePreference;
            const auto ordinaryOrBorderStrip =
                combinatorialTagsEqual(0) ||
                combinatorialTagsEqual(borderStripsTag);
            const auto allCompleteFriendly =
                ordinaryOrBorderStrip &&
                componentWeightGreaterThan(13) &&
                componentTermCountAtLeast(2) &&
                componentAllPowerSumTermsCompleteFriendly();
            const auto hybrid =
                ordinaryOrBorderStrip &&
                componentWeightGreaterThan(13) &&
                componentTermCountAtLeast(8) &&
                componentCompleteFriendlyFractionAtLeast(
                    1, 4, 8) &&
                componentHasMixedCompleteFriendlyPowerSumTerms();

            return std::vector<BasisConversionPlanCase>{
                {plethysmComplete,
                 completeFormula},
                {littlewoodRichardsonComplete,
                 completeFormula},
                {pieriComplete,
                 completeFormula},
                {largeBorderStripWithoutCommonOne,
                 completeFormula},
                {borderStripWithCommonOne,
                 completeFormula},
                {allCompleteFriendly,
                 completeFormula},
                // Delegate the complete component to the fixed term-level
                // hybrid. That child partitions its realized component
                // without re-entering the picker.
                {hybrid, termHybridFormula},
                {otherwise(), abacusFormula}};
          };
      plans.push_back({
          {"PowerSum->Schur:default-policy"},
          BasisKind::PowerSum,
          BasisKind::Schur,
          always(),
          ExpressionPieceKind::HomogeneousComponents,
          schurComponentPolicyCases()});

      // Named compositions and hybrids use the same plan representation as
      // direct kernels. Child identifiers completely determine execution.
      plans.push_back({
          {"PowerSum->Schur:via-complete"},
          BasisKind::PowerSum,
          BasisKind::Schur,
          always(),
          ExpressionPieceKind::WholeExpression,
          {{otherwise(),
            compositionFormula({
                {"PowerSum->Complete:logarithm-formula"},
                {"Complete->Schur:recursive-transition"}})}}});
      plans.push_back({
          {"PowerSum->Schur:complete-friendly-hybrid-plan"},
          BasisKind::PowerSum,
          BasisKind::Schur,
          always(),
          ExpressionPieceKind::Terms,
          {{powerSumIndexIsCompleteFriendly(),
            compositionFormula({
                {"PowerSum->Complete:logarithm-formula"},
                {"Complete->Schur:recursive-transition"}})},
           {otherwise(),
            kernelFormula(
                BasisConversionKernel::
                    PowerSumsToSchurAbacusRimHooks)}}});
      plans.push_back({
          {"PowerSum->SchurOmega:default-policy"},
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          always(),
          ExpressionPieceKind::WholeExpression,
          {{otherwise(),
            compositionFormula({
                {"PowerSum->Schur:default-policy"},
                {"Schur->SchurOmega:partition-conjugation"}})}}});
      plans.push_back({
          {"PowerSum->SchurOmega:via-complete"},
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          always(),
          ExpressionPieceKind::WholeExpression,
          {{otherwise(),
            compositionFormula({
                {"PowerSum->Schur:via-complete"},
                {"Schur->SchurOmega:partition-conjugation"}})}}});
      plans.push_back({
          {"PowerSum->SchurOmega:complete-friendly-hybrid-plan"},
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          always(),
          ExpressionPieceKind::WholeExpression,
          {{otherwise(),
            compositionFormula({
                {"PowerSum->Schur:complete-friendly-hybrid-plan"},
                {"Schur->SchurOmega:partition-conjugation"}})}}});

      // Every broad built-in X -> Y fallback is a named top-level plan. The
      // intermediate power-sum plans below are fixed child identifiers, so
      // realizing their intermediate expressions never invokes the picker.
      struct BroadPlanFamily
      {
        BasisKind kind;
        std::string_view name;
        std::string_view toPowerSums;
        std::string_view fromPowerSums;
        // Dense expansions in this source family use the established
        // source -> p -> h -> Schur composition. A false value delegates the
        // realized power-sum expression to the support-sensitive Schur plan.
        bool schurViaComplete;
      };
      const std::vector<BroadPlanFamily> families{
          {BasisKind::Complete,
           "Complete",
           "Complete->PowerSum:classical-formula",
           "PowerSum->Complete:logarithm-formula",
           true},
          {BasisKind::Elementary,
           "Elementary",
           "Elementary->PowerSum:classical-formula",
           "PowerSum->Elementary:logarithm-formula",
           true},
          {BasisKind::Schur,
           "Schur",
           "Schur->PowerSum:characters",
           "PowerSum->Schur:default-policy",
           true},
          {BasisKind::SchurOmega,
           "SchurOmega",
           "SchurOmega->PowerSum:characters",
           "PowerSum->SchurOmega:default-policy",
           true},
          {BasisKind::Monomial,
           "Monomial",
           "Monomial->PowerSum:transition",
           "PowerSum->Monomial:transition",
           false},
          {BasisKind::Forgotten,
           "Forgotten",
           "Forgotten->PowerSum:transition",
           "PowerSum->Forgotten:transition",
           false},
          {BasisKind::HallLittlewoodQGenerator,
           "HallLittlewoodQGenerator",
           "HallLittlewoodQGenerator->PowerSum:classical-formula",
           "PowerSum->HallLittlewoodQGenerator:logarithm-formula",
           true},
          {BasisKind::HallLittlewoodBGenerator,
           "HallLittlewoodBGenerator",
           "HallLittlewoodBGenerator->PowerSum:classical-formula",
           "PowerSum->HallLittlewoodBGenerator:logarithm-formula",
           true},
          {BasisKind::HallLittlewoodQ,
           "HallLittlewoodQ",
           "HallLittlewoodQ->PowerSum:raising-operators",
           "PowerSum->HallLittlewoodQ:triangular-reduction",
           true},
          {BasisKind::HallLittlewoodB,
           "HallLittlewoodB",
           "HallLittlewoodB->PowerSum:raising-operators",
           "PowerSum->HallLittlewoodB:triangular-reduction",
           true},
          {BasisKind::HallLittlewoodP,
           "HallLittlewoodP",
           "HallLittlewoodP->PowerSum:capital-normalization",
           "PowerSum->HallLittlewoodP:triangular-reduction",
           true},
          {BasisKind::HallLittlewoodPOmega,
           "HallLittlewoodPOmega",
           "HallLittlewoodPOmega->PowerSum:capital-normalization",
           "PowerSum->HallLittlewoodPOmega:triangular-reduction",
           true}};
      for (const auto& source : families)
        for (const auto& target : families)
          {
            if (source.kind == target.kind) continue;
            ConversionFormula broadFormula;
            if (target.kind == BasisKind::Schur &&
                source.schurViaComplete)
              {
                // Dense built-in sources historically convert to Schur
                // through complete functions. Monomial and forgotten
                // expansions are usually sparse after conversion to power
                // sums, so they deliberately retain the support-sensitive
                // p -> Schur child plan below.
                broadFormula = compositionFormula({
                    {std::string(source.toPowerSums)},
                    {"PowerSum->Complete:logarithm-formula"},
                    {"Complete->Schur:recursive-transition"}});
              }
            else
              broadFormula = compositionFormula({
                  {std::string(source.toPowerSums)},
                  {std::string(target.fromPowerSums)}});
            plans.push_back({
                {std::string(source.name) + "->" +
                 std::string(target.name) +
                 ":via-PowerSum-default"},
                source.kind,
                target.kind,
                always(),
                ExpressionPieceKind::WholeExpression,
                {{otherwise(), std::move(broadFormula)}}});
          }
      return plans;
    }();
    return database;
  }

// ============================================================================
// Declarative Plan Contracts And Validation
// ============================================================================

const SymmetricEngineRing::BasisConversionPlanDefinition *
SymmetricEngineRing::basisConversionPlanDefinition(
    const BasisConversionPlanId& id)
{
    const auto& plansById = basisConversionPlansById();
    const auto found = plansById.find(id.value);
    return found == plansById.end() ? nullptr : found->second;
  }

const std::map<
    std::string,
    const SymmetricEngineRing::BasisConversionPlanDefinition *>&
SymmetricEngineRing::basisConversionPlansById()
{
    static const std::map<
        std::string,
        const BasisConversionPlanDefinition *> index = [] {
      std::map<
          std::string,
          const BasisConversionPlanDefinition *> result;
      for (const auto& plan : basisConversionPlanDatabase())
        if (plan.id.value.empty() ||
            !result.emplace(plan.id.value, &plan).second)
          throw exc::engine_error(
              "basis-conversion plan identifiers must be "
              "nonempty and unique");
      return result;
    }();
    return index;
  }

const std::map<
    std::pair<
        SymmetricEngineRing::BasisKind,
        SymmetricEngineRing::BasisKind>,
    std::vector<
        const SymmetricEngineRing::BasisConversionPlanDefinition *>>&
SymmetricEngineRing::basisConversionPlansByEndpoints()
{
    static const std::map<
        std::pair<BasisKind, BasisKind>,
        std::vector<const BasisConversionPlanDefinition *>>
        index = [] {
          std::map<
              std::pair<BasisKind, BasisKind>,
              std::vector<const BasisConversionPlanDefinition *>>
              result;
          for (const auto& plan : basisConversionPlanDatabase())
            result[{plan.sourceBasisKind, plan.targetBasisKind}].
                push_back(&plan);
          return result;
        }();
    return index;
  }

SymmetricEngineRing::BasisConversionKernelContract
SymmetricEngineRing::basisConversionKernelContract(
    BasisConversionKernel kernel,
    BasisKind source,
    BasisKind target)
{
    switch (kernel)
      {
        case BasisConversionKernel::Unavailable:
          return false;
        case BasisConversionKernel::CompleteToPowerSumsClassical:
          return source == BasisKind::Complete &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            ElementaryToPowerSumsClassical:
          return source == BasisKind::Elementary &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::SchurToPowerSumsCharacters:
          return source == BasisKind::Schur &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            SchurOmegaToPowerSumsCharacters:
          return source == BasisKind::SchurOmega &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            MonomialToPowerSumsTransition:
          return source == BasisKind::Monomial &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            ForgottenToPowerSumsTransition:
          return source == BasisKind::Forgotten &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            HallLittlewoodQGeneratorToPowerSumsClassical:
          return source ==
                     BasisKind::HallLittlewoodQGenerator &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            HallLittlewoodBGeneratorToPowerSumsClassical:
          return source ==
                     BasisKind::HallLittlewoodBGenerator &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            HallLittlewoodQToPowerSumsRaising:
          return source == BasisKind::HallLittlewoodQ &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            HallLittlewoodBToPowerSumsRaising:
          return source == BasisKind::HallLittlewoodB &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            HallLittlewoodPToPowerSumsNormalization:
          return source == BasisKind::HallLittlewoodP &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            HallLittlewoodPOmegaToPowerSumsNormalization:
          return source == BasisKind::HallLittlewoodPOmega &&
                 target == BasisKind::PowerSum;
        case BasisConversionKernel::
            PowerSumsToCompleteLogarithm:
          return source == BasisKind::PowerSum &&
                 target == BasisKind::Complete;
        case BasisConversionKernel::
            PowerSumsToElementaryLogarithm:
          return source == BasisKind::PowerSum &&
                 target == BasisKind::Elementary;
        case BasisConversionKernel::
            PowerSumsToSchurBorderStrips:
        case BasisConversionKernel::
            PowerSumsToSchurAbacusRimHooks:
        case BasisConversionKernel::
            PowerSumsToSchurCharacters:
          return source == BasisKind::PowerSum &&
                 target == BasisKind::Schur;
        case BasisConversionKernel::
            PowerSumsToSchurOmegaBorderStrips:
        case BasisConversionKernel::
            PowerSumsToSchurOmegaAbacusRimHooks:
        case BasisConversionKernel::
            PowerSumsToSchurOmegaCharacters:
          return source == BasisKind::PowerSum &&
                 target == BasisKind::SchurOmega;
        case BasisConversionKernel::
            PowerSumsToHallLittlewoodQGeneratorLogarithm:
          return source == BasisKind::PowerSum &&
                 target ==
                     BasisKind::HallLittlewoodQGenerator;
        case BasisConversionKernel::
            PowerSumsToHallLittlewoodBGeneratorLogarithm:
          return source == BasisKind::PowerSum &&
                 target ==
                     BasisKind::HallLittlewoodBGenerator;
        case BasisConversionKernel::
            PowerSumSingleCyclesToHallLittlewoodGreen:
        case BasisConversionKernel::
            PowerSumIndexToHallLittlewoodGreenDuality:
        case BasisConversionKernel::
            PowerSumsToHallLittlewoodTriangular:
          return source == BasisKind::PowerSum &&
                 (target == BasisKind::HallLittlewoodQ ||
                  target == BasisKind::HallLittlewoodB ||
                  target == BasisKind::HallLittlewoodP ||
                  target == BasisKind::HallLittlewoodPOmega);
        case BasisConversionKernel::
            PowerSumsToMonomialTransition:
          return source == BasisKind::PowerSum &&
                 target == BasisKind::Monomial;
        case BasisConversionKernel::
            PowerSumsToForgottenTransition:
          return source == BasisKind::PowerSum &&
                 target == BasisKind::Forgotten;
        case BasisConversionKernel::HallLittlewoodNormalization:
          return
              (source == BasisKind::HallLittlewoodQ &&
               target == BasisKind::HallLittlewoodP) ||
              (source == BasisKind::HallLittlewoodP &&
               target == BasisKind::HallLittlewoodQ) ||
              (source == BasisKind::HallLittlewoodB &&
               target == BasisKind::HallLittlewoodPOmega) ||
              (source == BasisKind::HallLittlewoodPOmega &&
               target == BasisKind::HallLittlewoodB);
        case BasisConversionKernel::SchurOmegaConjugation:
          return
              (source == BasisKind::Schur &&
               target == BasisKind::SchurOmega) ||
              (source == BasisKind::SchurOmega &&
               target == BasisKind::Schur);
        case BasisConversionKernel::SchurToCompleteJacobiTrudi:
          return source == BasisKind::Schur &&
                 target == BasisKind::Complete;
        case BasisConversionKernel::
            SchurOmegaToElementaryJacobiTrudi:
          return source == BasisKind::SchurOmega &&
                 target == BasisKind::Elementary;
        case BasisConversionKernel::CompleteToSchurRecursive:
          return source == BasisKind::Complete &&
                 target == BasisKind::Schur;
        case BasisConversionKernel::
            HallLittlewoodGeneratorToCapitalTriangular:
          return
              (source ==
                   BasisKind::HallLittlewoodQGenerator &&
               (target == BasisKind::HallLittlewoodQ ||
                target == BasisKind::HallLittlewoodP)) ||
              (source ==
                   BasisKind::HallLittlewoodBGenerator &&
               (target == BasisKind::HallLittlewoodB ||
                target == BasisKind::HallLittlewoodPOmega));
      }
    return false;
  }

void SymmetricEngineRing::validateBasisConversionPlanDatabase()
{
    static const bool validated = [] {
      const auto& database = basisConversionPlanDatabase();
      const auto& plansById = basisConversionPlansById();
      for (const auto& plan : database)
        {
          try
            {
              validateExpressionCondition(
                  plan.applicability,
                  ExpressionPieceKind::WholeExpression);
            }
          catch (const std::exception& exception)
            {
              throw exc::engine_error(
                  "invalid applicability condition for plan " +
                  plan.id.value + ": " + exception.what());
            }
          try
            {
              validateExpressionCondition(
                  plan.outputGuarantee,
                  ExpressionPieceKind::WholeExpression);
            }
          catch (const std::exception& exception)
            {
              throw exc::engine_error(
                  "invalid output guarantee for plan " +
                  plan.id.value + ": " + exception.what());
            }
          if (plan.cases.empty())
            throw exc::engine_error(
                "basis-conversion plan " + plan.id.value +
                " has no cases");
          for (size_t i = 0; i < plan.cases.size(); ++i)
            {
              const auto& item = plan.cases[i];
              const bool isOtherwise =
                  item.condition.kind ==
                  ExpressionConditionKind::Otherwise;
              if (isOtherwise !=
                  (i + 1 == plan.cases.size()))
                throw exc::engine_error(
                    "basis-conversion plan " + plan.id.value +
                    " must end in exactly one otherwise case");
              try
                {
                  validateExpressionCondition(
                      item.condition,
                      plan.pieces,
                      isOtherwise);
                }
              catch (const std::exception& exception)
                {
                  throw exc::engine_error(
                      "invalid case condition for plan " +
                      plan.id.value + ": " + exception.what());
                }
              if (item.formula.kind ==
                  ConversionFormulaKind::AtomicKernel)
                {
                  const auto kernelContract =
                      basisConversionKernelContract(
                          item.formula.kernel,
                          plan.sourceBasisKind,
                          plan.targetBasisKind);
                  if (!item.formula.childPlans.empty() ||
                      !kernelContract.supportsEndpoints ||
                      !expressionConditionImplies(
                          kernelContract.outputGuarantee,
                          plan.outputGuarantee))
                    throw exc::engine_error(
                        "basis-conversion plan " +
                        plan.id.value +
                        " has an invalid atomic-kernel contract");
                }
              else if (
                  item.formula.kernel !=
                      BasisConversionKernel::Unavailable ||
                  item.formula.childPlans.empty())
                throw exc::engine_error(
                    "basis-conversion plan " + plan.id.value +
                    " has an invalid named-plan formula");
              item.formula.resolvedChildPlans.clear();
            }
        }

      std::map<std::string, std::vector<std::string>>
          dependencies;
      for (const auto& plan : database)
        for (const auto& item : plan.cases)
          {
            if (item.formula.kind !=
                ConversionFormulaKind::PlanComposition)
              continue;
            BasisKind current = plan.sourceBasisKind;
            ExpressionCondition guaranteedCondition =
                plan.pieces ==
                        ExpressionPieceKind::WholeExpression
                    ? plan.applicability
                    : always();
            if (plan.pieces ==
                    ExpressionPieceKind::WholeExpression &&
                item.condition.kind !=
                    ExpressionConditionKind::Otherwise)
              guaranteedCondition =
                  guaranteedCondition && item.condition;
            for (const auto& childId :
                 item.formula.childPlans)
              {
                auto found = plansById.find(childId.value);
                if (found == plansById.end())
                  throw exc::engine_error(
                      "basis-conversion plan " +
                      plan.id.value +
                      " references unknown child plan " +
                      childId.value);
                const auto& child = *found->second;
                item.formula.resolvedChildPlans.push_back(
                    found->second);
                if (child.sourceBasisKind != current)
                  throw exc::engine_error(
                      "basis-conversion plan " +
                      plan.id.value +
                      " has incompatible child-plan endpoints");
                if (!expressionConditionImplies(
                        guaranteedCondition,
                        child.applicability))
                  throw exc::engine_error(
                      "basis-conversion composition " +
                      plan.id.value +
                      " references a child with an unproved "
                      "applicability requirement");
                current = child.targetBasisKind;
                guaranteedCondition =
                    child.outputGuarantee;
                dependencies[plan.id.value].push_back(
                    child.id.value);
              }
            if (current != plan.targetBasisKind)
              throw exc::engine_error(
                  "basis-conversion plan " + plan.id.value +
                  " has the wrong composition target");
            if (!expressionConditionImplies(
                    guaranteedCondition,
                    plan.outputGuarantee))
              throw exc::engine_error(
                  "basis-conversion composition " +
                  plan.id.value +
                  " does not prove its declared output guarantee");
          }

      enum class VisitState { Unseen, Active, Complete };
      std::map<std::string, VisitState> state;
      std::function<void(const std::string&)> visit =
          [&](const std::string& id) {
            if (state[id] == VisitState::Complete) return;
            if (state[id] == VisitState::Active)
              throw exc::engine_error(
                  "basis-conversion plan dependencies are cyclic at " +
                  id);
            state[id] = VisitState::Active;
            for (const auto& child : dependencies[id])
              visit(child);
            state[id] = VisitState::Complete;
          };
      for (const auto& plan : database)
        visit(plan.id.value);
      return true;
    }();
    (void) validated;
  }

const std::vector<SymmetricEngineRing::ResolvedBasisConversionPlan>&
SymmetricEngineRing::registeredBasisConversionPlans(
    int sourceBasisId,
    int targetBasisId) const
{
    static const std::vector<ResolvedBasisConversionPlan> empty;
    requireBasis(sourceBasisId);
    requireBasis(targetBasisId);
    if (error()) return empty;
    const std::pair<int, int> endpoint{
        sourceBasisId, targetBasisId};
    auto cached =
        resolvedBasisConversionPlanCache.find(endpoint);
    if (cached != resolvedBasisConversionPlanCache.end())
      return cached->second;
    validateBasisConversionPlanDatabase();
    const BasisKind source = basisKindForId(sourceBasisId);
    const BasisKind target = basisKindForId(targetBasisId);
    std::vector<ResolvedBasisConversionPlan> result;
    const auto& plansByEndpoints =
        basisConversionPlansByEndpoints();
    auto found = plansByEndpoints.find({source, target});
    if (found != plansByEndpoints.end())
      {
        result.reserve(found->second.size());
        for (const auto *plan : found->second)
          result.push_back(
              {plan, sourceBasisId, targetBasisId});
      }
    auto inserted =
        resolvedBasisConversionPlanCache.emplace(
            endpoint, std::move(result));
    return inserted.first->second;
  }

SymmetricEngineRing::ResolvedBasisConversionPlan
SymmetricEngineRing::resolveBasisConversionPlan(
    const BasisConversionPlanId& id) const
{
    validateBasisConversionPlanDatabase();
    const auto *definition =
        basisConversionPlanDefinition(id);
    if (definition == nullptr)
      {
        ERROR("unknown basis-conversion plan: ",
              id.value.c_str());
        return {};
      }
    return resolveBasisConversionPlanDefinition(definition);
  }

SymmetricEngineRing::ResolvedBasisConversionPlan
SymmetricEngineRing::resolveBasisConversionPlanDefinition(
    const BasisConversionPlanDefinition *definition) const
{
    if (definition == nullptr) return {};
    auto cached =
        resolvedBasisConversionEndpointCache.find(
            definition->id.value);
    if (cached !=
        resolvedBasisConversionEndpointCache.end())
      return {
          definition,
          cached->second.first,
          cached->second.second};
    const int sourceBasisId =
        requiredBasisIdForKind(
            definition->sourceBasisKind);
    const int targetBasisId =
        requiredBasisIdForKind(
            definition->targetBasisKind);
    if (error()) return {};
    resolvedBasisConversionEndpointCache.emplace(
        definition->id.value,
        std::make_pair(sourceBasisId, targetBasisId));
    return {definition, sourceBasisId, targetBasisId};
  }

// ============================================================================
// Complete-Plan Applicability And Performance Policy
// ============================================================================

bool SymmetricEngineRing::basisConversionPlanApplicable(
    const ResolvedBasisConversionPlan& plan,
    ring_elem expression,
    const ExpressionFacts& facts) const
{
    if (!plan.valid() ||
        !facts.canonicalExpansionInBasis(
            plan.sourceBasisId))
      return false;
    if (plan.definition->applicability.kind ==
        ExpressionConditionKind::Always)
      return true;
    const auto contexts =
        buildExpressionConditionContexts(
            expression,
            ExpressionPieceKind::WholeExpression,
            facts);
    return contexts.size() == 1 &&
           expressionConditionHolds(
               plan.definition->applicability,
               contexts.front());
  }

size_t SymmetricEngineRing::basisConversionPlanCost(
    const ResolvedBasisConversionPlan& plan,
    const ExpressionFacts& facts,
    CombinatorialTags combinatorialTags) const
{
    if (!plan.valid())
      return std::numeric_limits<size_t>::max();
    size_t terms = std::max<size_t>(1, facts.termCount);
    size_t weight = static_cast<size_t>(
        std::max(
            0, facts.homogeneousWeight.value_or(0)));
    size_t cost = 100 + terms + weight;
    bool hasComposition = false;
    for (const auto& item : plan.definition->cases)
      {
        hasComposition =
            hasComposition ||
            item.formula.kind ==
                ConversionFormulaKind::PlanComposition;
        if (item.formula.kind !=
            ConversionFormulaKind::AtomicKernel)
          continue;
        switch (item.formula.kernel)
          {
            case BasisConversionKernel::
                HallLittlewoodNormalization:
            case BasisConversionKernel::SchurOmegaConjugation:
            case BasisConversionKernel::CompleteToSchurRecursive:
            case BasisConversionKernel::
                HallLittlewoodGeneratorToCapitalTriangular:
              cost = std::min(cost, 20 + terms);
              break;
            case BasisConversionKernel::
                PowerSumsToSchurCharacters:
            case BasisConversionKernel::
                PowerSumsToSchurOmegaCharacters:
              cost += 4 * terms + weight;
              break;
            default:
              break;
          }
      }
    if (hasComposition) cost += 1000;
    if (hasCombinatorialTag(
            combinatorialTags,
            CombinatorialTag::Plethysm) &&
        cost > 10)
      cost -= 10;
    return cost;
  }

SymmetricEngineRing::ResolvedBasisConversionPlan
SymmetricEngineRing::selectBasisConversionPlan(
    ring_elem expression,
    int sourceBasisId,
    int targetBasisId,
    const ExpressionFacts& facts,
    CombinatorialTags combinatorialTags,
    const std::optional<std::string>& forcedIdentifier,
    ExpressionFacts *selectionFacts) const
{
    if (basisConversionPlanExecutionDepth != 0)
      {
        ERROR("basis-conversion plan execution attempted to invoke "
              "the performance picker");
        return {};
      }
    if (!facts.canonicalExpansionInBasis(sourceBasisId))
      {
        ERROR("the basis-conversion picker requires a canonical "
              "single-source-basis expansion");
        return {};
      }
    const auto& plans =
        registeredBasisConversionPlans(
            sourceBasisId, targetBasisId);
    const BasisKind sourceKind =
        basisKindForId(sourceBasisId);
    const BasisKind targetKind =
        basisKindForId(targetBasisId);
    std::optional<std::string> requested =
        forcedIdentifier;
    if (!requested)
      {
        const char *forced =
            std::getenv(
                "M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN");
        if (forced != nullptr)
          requested = forced;
      }
    ExpressionFacts selectionProfile = facts;
    selectionProfile.combinatorialTags = combinatorialTags;
    auto publishSelectionFacts = [&] {
      if (selectionFacts != nullptr)
        *selectionFacts = selectionProfile;
    };
    // Most endpoints have one unconditional mathematical plan. Once forcing
    // has been ruled out, return that stable entry without allocating an
    // enriched profile or running general ranking policy.
    if (!requested &&
        plans.size() == 1 &&
        plans.front().definition->applicability.kind ==
            ExpressionConditionKind::Always)
      {
        publishSelectionFacts();
        return plans.front();
      }

    if (requested)
      {
        for (const auto& plan : plans)
          if (plan.definition->id.value == *requested)
            {
              // Unconditional forced plans need only the canonical source
              // contract already checked above. Conditional plans request
              // their endpoint-specific profile before applicability is
              // evaluated.
              if (plan.definition->applicability.kind !=
                  ExpressionConditionKind::Always)
                enrichBasisConversionPlanSelectionFacts(
                    expression, targetBasisId, selectionProfile);
              publishSelectionFacts();
              if (basisConversionPlanApplicable(
                      plan, expression, selectionProfile))
                return plan;
              break;
            }
        ERROR("the requested conversion plan is unknown or "
              "inapplicable: ",
              requested->c_str());
        return {};
      }

    if (sourceKind == BasisKind::PowerSum &&
        (targetKind == BasisKind::Schur ||
         targetKind == BasisKind::SchurOmega))
      {
        // The component plan itself owns all p -> Schur policy facts. Select
        // it before whole-expression enrichment so its component partitioner
        // is the only code that scans those profiles.
        const std::string_view defaultId =
            targetKind == BasisKind::SchurOmega
                ? "PowerSum->SchurOmega:default-policy"
                : "PowerSum->Schur:default-policy";
        for (const auto& plan : plans)
          if (plan.definition->id.value == defaultId)
            {
              publishSelectionFacts();
              return plan;
            }
      }

    enrichBasisConversionPlanSelectionFacts(
        expression, targetBasisId, selectionProfile);
    publishSelectionFacts();
    auto findId =
        [&](const std::string& id)
            -> std::optional<ResolvedBasisConversionPlan> {
          for (const auto& plan : plans)
            if (plan.definition->id.value == id &&
                basisConversionPlanApplicable(
                    plan, expression, selectionProfile))
              return plan;
          return std::nullopt;
        };

    if (sourceKind == BasisKind::PowerSum &&
        isHallLittlewoodCapitalBasisKind(targetKind))
      {
        std::string targetName;
        switch (targetKind)
          {
            case BasisKind::HallLittlewoodQ:
              targetName = "HallLittlewoodQ";
              break;
            case BasisKind::HallLittlewoodB:
              targetName = "HallLittlewoodB";
              break;
            case BasisKind::HallLittlewoodP:
              targetName = "HallLittlewoodP";
              break;
            case BasisKind::HallLittlewoodPOmega:
              targetName = "HallLittlewoodPOmega";
              break;
            default:
              break;
          }
        const std::string suffix =
            selectionProfile.allPowerSumTermsSingleCycles.
                    value_or(false)
                ? ":single-cycles-Green"
                : (selectionProfile.singleTerm() &&
                           selectionProfile.singleBasisElement()
                       ? ":Green-duality"
                       : ":triangular-reduction");
        auto selected = findId(
            "PowerSum->" + targetName + suffix);
        if (selected) return *selected;
      }

    const ResolvedBasisConversionPlan *best = nullptr;
    size_t bestCost =
        std::numeric_limits<size_t>::max();
    for (const auto& plan : plans)
      {
        if (!basisConversionPlanApplicable(
                plan, expression, selectionProfile))
          continue;
        const size_t cost = basisConversionPlanCost(
            plan, selectionProfile, combinatorialTags);
        if (best == nullptr || cost < bestCost)
          {
            best = &plan;
            bestCost = cost;
          }
      }
    if (best != nullptr) return *best;
    ERROR("the basis-conversion registry has no applicable "
          "complete plan from ",
          basisKeyForId(sourceBasisId).c_str(),
          " to ",
          basisKeyForId(targetBasisId).c_str());
    return {};
  }

// ============================================================================
// Atomic Basis-Conversion Kernel Executor
// ============================================================================
// This executor contains no applicability or performance policy. The switch is
// grouped by mathematical source/target family so a contributor can locate the
// one kernel associated with a registry entry.

ring_elem SymmetricEngineRing::executeBasisConversionKernel(
    BasisConversionKernel kernel,
    int sourceBasisId,
    int targetBasisId,
    ring_elem expression,
    CombinatorialTags combinatorialTags,
    ExpressionFacts *resultFacts,
    const std::optional<int>& knownHomogeneousWeight) const
{
    auto powerSumsToSchurOmega =
        [&](BasisConversionKernel kernel) {
          ring_elem omegaInput = omegaPowerSums(expression);
          if (error()) return zero();
          int schurId = requiredBasisIdForKind(BasisKind::Schur);
          if (error()) return zero();
          std::string schurDisplay = displayForBasis(schurId);
          int schurOrder = basisOrderForId(schurId);
          ring_elem inSchur;
          switch (kernel)
            {
              case BasisConversionKernel::
                  PowerSumsToSchurOmegaBorderStrips:
                inSchur = powerSumsToSchurViaBorderStrips(
                    omegaInput, schurId, schurDisplay, schurOrder);
                break;
              case BasisConversionKernel::
                  PowerSumsToSchurOmegaAbacusRimHooks:
                inSchur = powerSumsToSchurViaAbacusRimHooks(
                    omegaInput, schurId, schurDisplay, schurOrder);
                break;
              case BasisConversionKernel::
                  PowerSumsToSchurOmegaCharacters:
                inSchur = powerSumsToSchurLikeViaCharacters(
                    omegaInput,
                    schurId,
                    schurOrder,
                    schurDisplay,
                    false);
                break;
              default:
                ERROR("expected a Schur Omega conversion kernel");
                return zero();
            }
          if (error()) return zero();
          return replaceSingleBasis(
              inSchur, schurId, targetBasisId);
        };

    const std::string targetDisplay =
        displayForBasis(targetBasisId);
    const int targetOrder =
        basisOrderForId(targetBasisId);
    ring_elem result;
    switch (kernel)
      {
        case BasisConversionKernel::Unavailable:
          ERROR("an unavailable basis-conversion plan was executed");
          return zero();

        // Canonical built-in expansions to power sums.
        case BasisConversionKernel::CompleteToPowerSumsClassical:
        case BasisConversionKernel::ElementaryToPowerSumsClassical:
        case BasisConversionKernel::SchurToPowerSumsCharacters:
        case BasisConversionKernel::SchurOmegaToPowerSumsCharacters:
        case BasisConversionKernel::MonomialToPowerSumsTransition:
        case BasisConversionKernel::ForgottenToPowerSumsTransition:
        case BasisConversionKernel::
            HallLittlewoodQGeneratorToPowerSumsClassical:
        case BasisConversionKernel::
            HallLittlewoodBGeneratorToPowerSumsClassical:
        case BasisConversionKernel::
            HallLittlewoodQToPowerSumsRaising:
        case BasisConversionKernel::
            HallLittlewoodBToPowerSumsRaising:
        case BasisConversionKernel::
            HallLittlewoodPToPowerSumsNormalization:
        case BasisConversionKernel::
            HallLittlewoodPOmegaToPowerSumsNormalization:
          result = canonicalExpressionToPowerSumsViaBasisFormulas(
              expression,
              basisKindForId(sourceBasisId));
          break;

        // Classical power-sum targets.
        case BasisConversionKernel::PowerSumsToCompleteLogarithm:
          result = powerSumsToCompleteViaLogarithmFormula(
              expression, targetBasisId, targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToElementaryLogarithm:
          result = powerSumsToElementaryViaLogarithmFormula(
              expression, targetBasisId, targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToSchurBorderStrips:
          result = powerSumsToSchurViaBorderStrips(
              expression,
              targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToSchurAbacusRimHooks:
          result = powerSumsToSchurViaAbacusRimHooks(
              expression,
              targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToSchurCharacters:
          result = powerSumsToSchurLikeViaCharacters(
              expression,
              targetBasisId,
              targetOrder,
              targetDisplay,
              false);
          break;

        // Schur Omega targets share the Schur kernels after omega on p.
        case BasisConversionKernel::
            PowerSumsToSchurOmegaBorderStrips:
        case BasisConversionKernel::
            PowerSumsToSchurOmegaAbacusRimHooks:
        case BasisConversionKernel::
            PowerSumsToSchurOmegaCharacters:
          result = powerSumsToSchurOmega(kernel);
          break;

        // Hall--Littlewood generator and capital targets.
        case BasisConversionKernel::
            PowerSumsToHallLittlewoodQGeneratorLogarithm:
        case BasisConversionKernel::
            PowerSumsToHallLittlewoodBGeneratorLogarithm:
          {
            CoeffMap generators =
                powerSumsToHallGeneratorMapViaLogarithmFormula(
                    expression,
                    kernel == BasisConversionKernel::
                        PowerSumsToHallLittlewoodBGeneratorLogarithm);
            result = coeffMapToElement(
                generators,
                targetBasisId,
                targetDisplay,
                targetOrder,
                true);
            break;
          }
        case BasisConversionKernel::
            PowerSumSingleCyclesToHallLittlewoodGreen:
          result =
              powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials(
                  expression,
                  targetBasisId,
                  targetDisplay,
                  targetOrder);
          break;
        case BasisConversionKernel::
            PowerSumIndexToHallLittlewoodGreenDuality:
          result =
              powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality(
                  expression,
                  targetBasisId,
                  targetDisplay,
                  targetOrder);
          break;
        case BasisConversionKernel::
            PowerSumsToHallLittlewoodTriangular:
          result = powerSumsToHallLittlewoodViaTriangularReduction(
              expression,
              targetBasisId,
              targetDisplay,
              targetOrder);
          break;

        // Monomial and forgotten transition matrices.
        case BasisConversionKernel::PowerSumsToMonomialTransition:
        case BasisConversionKernel::PowerSumsToForgottenTransition:
          result = powerSumsToTargetViaTermwiseConversion(
              expression,
              targetBasisId,
              targetDisplay,
              targetOrder,
              isMultiplicativeBasis(targetBasisId));
          break;

        // Direct normalization, involution, and triangular transitions.
        case BasisConversionKernel::HallLittlewoodNormalization:
          {
            BasisKind sourceKind =
                basisKindForId(sourceBasisId);
            result =
                hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
                    expression,
                    sourceBasisId,
                    targetBasisId,
                    targetDisplay,
                    targetOrder,
                    sourceKind == BasisKind::HallLittlewoodQ ||
                        sourceKind == BasisKind::HallLittlewoodB);
            break;
          }
        case BasisConversionKernel::SchurOmegaConjugation:
          result = schurOmegaConversionViaPartitionConjugation(
              expression,
              sourceBasisId,
              targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::SchurToCompleteJacobiTrudi:
        case BasisConversionKernel::
            SchurOmegaToElementaryJacobiTrudi:
          result =
              canonicalSchurLikeExpressionToGeneratorsViaJacobiTrudi(
                  expression,
                  sourceBasisId,
                  targetBasisId);
          break;
        case BasisConversionKernel::CompleteToSchurRecursive:
          result = completeToSchurViaRecursiveTransition(
              expression,
              sourceBasisId,
              displayForBasis(sourceBasisId),
              basisOrderForId(sourceBasisId),
              isMultiplicativeBasis(sourceBasisId),
              targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::
            HallLittlewoodGeneratorToCapitalTriangular:
          if (!tryExpressionToHallLittlewoodViaTriangularReduction(
                  expression,
                  targetBasisId,
                  targetDisplay,
                  targetOrder,
                  result))
            {
              ERROR("a selected Hall-Littlewood triangular kernel "
                    "violated its input contract");
              return zero();
            }
          break;
      }
    if (error()) return zero();

    // Intermediate stages need exact facts to select their successor. A final
    // transient result does not: the selected kernel's output contract already
    // guarantees a canonical target expansion, and the owning workflow will
    // infer once after collecting all groups.
    if (resultFacts != nullptr)
      {
        *resultFacts = inferCanonicalExpansionFacts(
            result,
            targetBasisId,
            knownHomogeneousWeight);
        resultFacts->combinatorialTags = combinatorialTags;
        if (!resultFacts->canonicalExpansionInBasis(targetBasisId))
          {
            ERROR("a basis-conversion plan violated its canonical output "
                  "contract");
            return zero();
          }
      }
    mutablePolyValue(result)->combinatorialTags = combinatorialTags;
    return result;
  }

// ============================================================================
// Generic Declarative Formula And Plan Executor
// ============================================================================
// The executor interprets the already-selected plan. Atomic formulas call one
// kernel; compositions resolve their fixed child identifiers recursively.
// Neither path invokes performance policy or substitutes another plan.

ring_elem SymmetricEngineRing::executeBasisConversionFormula(
    const ConversionFormula& formula,
    int sourceBasisId,
    int targetBasisId,
    ring_elem expression,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts& inputFacts,
    ExpressionFacts *resultFacts) const
{
    if (!inputFacts.canonicalExpansionInBasis(sourceBasisId))
      {
        ERROR("a basis-conversion formula received input outside "
              "its canonical source-basis contract");
        return zero();
      }
    if (formula.kind == ConversionFormulaKind::AtomicKernel)
      {
        const auto kernelContract =
            basisConversionKernelContract(
                formula.kernel,
                basisKindForId(sourceBasisId),
                basisKindForId(targetBasisId));
        if (!kernelContract.supportsEndpoints)
          {
            ERROR("an atomic basis-conversion formula has invalid "
                  "runtime endpoints");
            return zero();
          }
        return executeBasisConversionKernel(
            formula.kernel,
            sourceBasisId,
            targetBasisId,
            expression,
            combinatorialTags,
            resultFacts,
            inputFacts.homogeneousWeight);
      }

    ring_elem current = expression;
    ExpressionFacts currentFacts = inputFacts;
    int currentBasisId = sourceBasisId;
    if (formula.resolvedChildPlans.size() !=
        formula.childPlans.size())
      {
        ERROR("a basis-conversion composition reached execution "
              "before its child plans were resolved");
        return zero();
      }
    for (const auto *childDefinition :
         formula.resolvedChildPlans)
      {
        ResolvedBasisConversionPlan child =
            resolveBasisConversionPlanDefinition(
                childDefinition);
        if (error()) return zero();
        if (child.sourceBasisId != currentBasisId)
          {
            ERROR("a basis-conversion composition has incompatible "
                  "runtime child endpoints");
            return zero();
          }
        current = executeBasisConversionPlan(
            child,
            current,
            combinatorialTags,
            currentFacts,
            nullptr);
        if (error()) return zero();
        currentBasisId = child.targetBasisId;
        // A child contract proves canonical target output, but ExpressionFacts
        // represents exact support rather than an abstract guarantee. Inspect
        // the realized intermediate before passing it to the next child so no
        // condition or contract check can observe fabricated term counts.
        currentFacts = inferCanonicalExpansionFacts(
            current,
            currentBasisId,
            inputFacts.homogeneousWeight);
        currentFacts.combinatorialTags = combinatorialTags;
      }
    if (currentBasisId != targetBasisId ||
        !currentFacts.canonicalExpansionInBasis(targetBasisId))
      {
        ERROR("a basis-conversion composition violated its target "
              "contract");
        return zero();
      }
    if (resultFacts != nullptr)
      {
        *resultFacts = std::move(currentFacts);
        resultFacts->combinatorialTags =
            combinatorialTags;
      }
    return current;
  }

ring_elem SymmetricEngineRing::executeBasisConversionPlan(
    const ResolvedBasisConversionPlan& plan,
    ring_elem expression,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts& inputFacts,
    ExpressionFacts *resultFacts) const
{
    if (!plan.valid() ||
        !inputFacts.canonicalExpansionInBasis(
            plan.sourceBasisId))
      {
        ERROR("an invalid complete basis-conversion plan reached "
              "the executor");
        return zero();
      }
    struct ExecutionScope
    {
      size_t& depth;
      explicit ExecutionScope(size_t& value) : depth(value)
      {
        ++depth;
      }
      ~ExecutionScope() { --depth; }
    } executionScope(basisConversionPlanExecutionDepth);

    const bool wholeExpressionOtherwise =
        plan.definition->pieces ==
            ExpressionPieceKind::WholeExpression &&
        plan.definition->cases.size() == 1 &&
        plan.definition->cases.front().condition.kind ==
            ExpressionConditionKind::Otherwise;
    ExpressionFacts executionFacts = inputFacts;
    executionFacts.combinatorialTags = combinatorialTags;
    const bool needsSelectionProfiles =
        plan.definition->applicability.kind !=
            ExpressionConditionKind::Always ||
        (plan.definition->pieces ==
             ExpressionPieceKind::WholeExpression &&
         !wholeExpressionOtherwise);
    // Term and homogeneous-component contexts derive their own profiles from
    // the realized pieces. Only whole-expression conditions consume the
    // endpoint-specific profile stored on ExpressionFacts.
    if (needsSelectionProfiles)
      enrichBasisConversionPlanSelectionFacts(
          expression,
          plan.targetBasisId,
          executionFacts);
    if (!basisConversionPlanApplicable(
            plan, expression, executionFacts))
      {
        ERROR("a selected basis-conversion plan is inapplicable "
              "to its runtime input: ",
              plan.definition->id.value.c_str());
        return zero();
      }
    if (std::getenv(
            "M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr)
      std::fprintf(
          stderr,
          "SymmetricRings conversion-plan: source=%s target=%s "
          "plan=%s\n",
          basisKeyForId(plan.sourceBasisId).c_str(),
          basisKeyForId(plan.targetBasisId).c_str(),
          plan.definition->id.value.c_str());

    if (wholeExpressionOtherwise)
      {
        // Most registered plans are one whole-expression formula. Their
        // declarative meaning is exactly a direct call on the original input;
        // constructing a one-piece partition and then copying and recollecting
        // the complete expression would add work without changing semantics.
        const bool needDirectFacts =
            resultFacts != nullptr ||
            plan.definition->outputGuarantee.kind !=
                ExpressionConditionKind::Always;
        ExpressionFacts directFacts;
        ring_elem directResult =
            executeBasisConversionFormula(
                plan.definition->cases.front().formula,
                plan.sourceBasisId,
                plan.targetBasisId,
                expression,
                combinatorialTags,
                executionFacts,
                needDirectFacts ? &directFacts : nullptr);
        if (error()) return zero();
        if (needDirectFacts &&
            !directFacts.canonicalExpansionInBasis(
                plan.targetBasisId))
          {
            ERROR("a complete basis-conversion plan violated its "
                  "canonical output contract");
            return zero();
          }
        if (plan.definition->outputGuarantee.kind !=
            ExpressionConditionKind::Always)
          {
            try
              {
                const auto outputContexts =
                    buildExpressionConditionContexts(
                        directResult,
                        ExpressionPieceKind::WholeExpression,
                        directFacts);
                if (outputContexts.size() != 1 ||
                    !expressionConditionHolds(
                        plan.definition->outputGuarantee,
                        outputContexts.front()))
                  {
                    ERROR("a complete basis-conversion plan violated its "
                          "declared output guarantee");
                    return zero();
                  }
              }
            catch (const std::exception& exception)
              {
                ERROR("basis-conversion output-guarantee validation "
                      "failed: ",
                      exception.what());
                return zero();
              }
          }
        mutablePolyValue(directResult)->combinatorialTags =
            combinatorialTags;
        if (resultFacts != nullptr)
          *resultFacts = std::move(directFacts);
        return directResult;
      }

    std::vector<ExpressionConditionContext> pieces;
    ExpressionConditionPartition partition;
    try
      {
        pieces = buildExpressionConditionContexts(
            expression,
            plan.definition->pieces,
            executionFacts);
        std::vector<ExpressionCondition> conditions;
        conditions.reserve(
            plan.definition->cases.size());
        for (const auto& item : plan.definition->cases)
          conditions.push_back(item.condition);
        partition = partitionExpressionConditionContexts(
            pieces,
            conditions,
            polyValue(expression)->terms.size());
      }
    catch (const std::exception& exception)
      {
        ERROR("basis-conversion plan partitioning failed: ",
              exception.what());
        return zero();
      }

    const auto& sourceTerms =
        polyValue(expression)->terms;
    VECTOR(SymmetricTerm) resultTerms;
    for (size_t caseIndex = 0;
         caseIndex < plan.definition->cases.size();
         ++caseIndex)
      {
        const auto& termPositions =
            partition.termPositionsByCase[caseIndex];
        if (termPositions.empty()) continue;

        // A case formula receives the complete subexpression assigned to that
        // case, not one invocation per term or component. This preserves the
        // declarative semantics and lets fixed child plans inspect the actual
        // aggregate support produced for the case.
        VECTOR(SymmetricTerm) caseTerms;
        caseTerms.reserve(termPositions.size());
        for (size_t position : termPositions)
          {
            if (position >= sourceTerms.size())
              {
                ERROR("a basis-conversion piece contains an invalid "
                      "term position");
                return zero();
              }
            caseTerms.push_back(sourceTerms[position]);
          }
        ring_elem caseExpression =
            fromTermVector(caseTerms, true);
        ExpressionFacts caseFacts;
        if (plan.definition->pieces ==
            ExpressionPieceKind::WholeExpression)
          caseFacts = executionFacts;
        else
          {
            caseFacts = inferCanonicalExpansionFacts(
                caseExpression,
                plan.sourceBasisId);
            caseFacts.combinatorialTags =
                combinatorialTags;
          }
        ExpressionFacts convertedFacts;
        ring_elem converted =
            executeBasisConversionFormula(
                plan.definition->cases[caseIndex].formula,
                plan.sourceBasisId,
                plan.targetBasisId,
                caseExpression,
                combinatorialTags,
                caseFacts,
                &convertedFacts);
        if (error()) return zero();
        if (!convertedFacts.canonicalExpansionInBasis(
                plan.targetBasisId))
          {
            ERROR("a basis-conversion case violated its canonical "
                  "target contract");
            return zero();
          }
        for (const auto& term :
             polyValue(converted)->terms)
          appendTermIfNonZero(
              resultTerms, term.coeff, term.monomial);
      }

    ring_elem result =
        fromTermVector(resultTerms, false);
    ExpressionFacts finalFacts =
        inferCanonicalExpansionFacts(
            result,
            plan.targetBasisId,
            executionFacts.homogeneousWeight);
    finalFacts.combinatorialTags = combinatorialTags;
    if (!finalFacts.canonicalExpansionInBasis(
            plan.targetBasisId))
      {
        ERROR("a complete basis-conversion plan violated its "
              "canonical output contract");
        return zero();
      }
    if (plan.definition->outputGuarantee.kind !=
        ExpressionConditionKind::Always)
      {
        try
          {
            const auto outputContexts =
                buildExpressionConditionContexts(
                    result,
                    ExpressionPieceKind::WholeExpression,
                    finalFacts);
            if (outputContexts.size() != 1 ||
                !expressionConditionHolds(
                    plan.definition->outputGuarantee,
                    outputContexts.front()))
              {
                ERROR("a complete basis-conversion plan violated its "
                      "declared output guarantee");
                return zero();
              }
          }
        catch (const std::exception& exception)
          {
            ERROR("basis-conversion output-guarantee validation "
                  "failed: ",
                  exception.what());
            return zero();
          }
      }
    mutablePolyValue(result)->combinatorialTags =
        combinatorialTags;
    if (resultFacts != nullptr)
      *resultFacts = std::move(finalFacts);
    return result;
  }

bool SymmetricEngineRing::checkAllApplicableBasisConversionPlans(
    ring_elem expression,
    int sourceBasisId,
    int targetBasisId,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts& inputFacts,
    ring_elem expectedResult) const
{
    if (std::getenv(
            "M2_SYMMETRIC_RINGS_CHECK_ALL_CONVERSION_PLANS") ==
        nullptr)
      return true;

    ExpressionFacts facts = inputFacts;
    facts.combinatorialTags = combinatorialTags;
    enrichBasisConversionPlanSelectionFacts(
        expression, targetBasisId, facts);
    size_t checked = 0;
    for (const auto& plan :
         registeredBasisConversionPlans(
             sourceBasisId, targetBasisId))
      {
        if (!basisConversionPlanApplicable(
                plan, expression, facts))
          continue;
        ring_elem candidate = executeBasisConversionPlan(
            plan,
            expression,
            combinatorialTags,
            facts);
        if (error()) return false;
        ++checked;
        if (!is_equal(candidate, expectedResult))
          {
            ERROR("registered basis-conversion plan disagrees with "
                  "the selected plan: ",
                  plan.definition->id.value.c_str());
            return false;
          }
      }
    if (checked == 0)
      {
        ERROR("no applicable basis-conversion plan was available "
              "for exhaustive checking");
        return false;
      }
    return true;
  }

// ============================================================================
// Canonical Group Conversion Helper
// ============================================================================
// This helper owns source-basis grouping, selected-plan execution, and final
// collection for input that is already normalized, skew-free, and product-free.
// Basis-specific partitioning belongs to the selected declarative plan.

ring_elem SymmetricEngineRing::convertCanonicalExpressionToBasis(
    ring_elem f,
    int targetBasisId,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts *knownFacts) const
{
    const bool trace =
        std::getenv(
            "M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr;

    // Exact zero/scalar metadata satisfies the target contract without basis
    // grouping or plan selection.
    if (knownFacts != nullptr &&
        knownFacts->singleFactorTermCount == 0)
      {
        ring_elem result = copyPolyValue(polyValue(f));
        attachExpressionFacts(
            result,
            *knownFacts,
            targetBasisId,
            combinatorialTags);
        return result;
      }

    // A pure canonical expansion can enter complete-plan selection directly.
    // The selected declarative plan owns its basis-specific partitioning.
    if (knownFacts != nullptr && knownFacts->expandedBasis)
      {
        int sourceBasisId = *knownFacts->expandedBasis;
        ExpressionFacts resultFacts;
        ring_elem result;
        if (sourceBasisId == targetBasisId)
          {
            result = copyPolyValue(polyValue(f));
            resultFacts = *knownFacts;
          }
        else
          {
            ExpressionFacts selectedInputFacts = *knownFacts;
            ResolvedBasisConversionPlan selection =
                selectBasisConversionPlan(
                    f,
                    sourceBasisId,
                    targetBasisId,
                    *knownFacts,
                    combinatorialTags,
                    std::nullopt,
                    &selectedInputFacts);
            if (error()) return zero();
            result = executeBasisConversionPlan(
                selection,
                f,
                combinatorialTags,
                selectedInputFacts,
                &resultFacts);
            if (error()) return zero();
            if (!checkAllApplicableBasisConversionPlans(
                    f,
                    sourceBasisId,
                    targetBasisId,
                    combinatorialTags,
                    selectedInputFacts,
                    result))
              return zero();
          }
        if (!resultFacts.canonicalExpansionInBasis(targetBasisId))
          {
            ERROR("the canonical metadata bypass did not produce "
                  "a target-basis expansion");
            return zero();
          }
        attachExpressionFacts(
            result,
            resultFacts,
            targetBasisId,
            combinatorialTags);
        return result;
      }

    struct PreparedConversion
    {
      ring_elem expression;
      ExpressionFacts facts;
      std::optional<ResolvedBasisConversionPlan> selection;
    };

    std::vector<std::pair<int, VECTOR(SymmetricTerm)>> groups;
    VECTOR(SymmetricTerm) scalarTerms;
    for (const auto& term : polyValue(f)->terms)
      {
        if (term.monomial.data.empty())
          {
            scalarTerms.push_back(term);
            continue;
          }
        if (atomLengthAt(term.monomial, 0) != term.monomial.data.size())
          {
            ERROR("the conversion workflow expected product-free input");
            return zero();
          }
        const int basisId = atomBasisIdAt(term.monomial, 0);
        auto group = std::find_if(
            groups.begin(),
            groups.end(),
            [&](const auto& candidate) {
              return candidate.first == basisId;
            });
        if (group == groups.end())
          {
            groups.push_back({basisId, {}});
            group = std::prev(groups.end());
          }
        group->second.push_back(term);
      }

    std::vector<PreparedConversion> preparedConversions;
    auto prepareGroup = [&](int sourceBasisId,
                            ring_elem expression,
                            ExpressionFacts facts) {
      std::optional<ResolvedBasisConversionPlan> selection;
      if (sourceBasisId != targetBasisId)
        {
          selection = selectBasisConversionPlan(
              expression,
              sourceBasisId,
              targetBasisId,
              facts,
              combinatorialTags,
              std::nullopt,
              &facts);
          if (error()) return;
        }
      preparedConversions.push_back({
          expression, std::move(facts), std::move(selection)});
    };

    if (trace)
      std::fprintf(
          stderr,
          "SymmetricRings conversion-stage: stage=group-by-basis "
          "groups=%zu scalar-terms=%zu\n",
          groups.size(),
          scalarTerms.size());
    for (auto& group : groups)
      {
        ring_elem groupExpression = fromTermVector(group.second, true);
        bool exactKnownGroup =
            knownFacts != nullptr &&
            knownFacts->expandedBasis == group.first &&
            knownFacts->scalarTermCount == 0 &&
            groups.size() == 1;
        ExpressionFacts groupFacts = exactKnownGroup
            ? *knownFacts
            : inferExpressionFacts(groupExpression);
        prepareGroup(
            group.first, groupExpression, std::move(groupFacts));
        if (error()) return zero();
      }

    // Selection for every source group is complete before any group executes.
    // Each selected definition fixes every case, kernel, and named child.
    VECTOR(SymmetricTerm) resultTerms = std::move(scalarTerms);
    for (const auto& prepared : preparedConversions)
      {
        ring_elem converted = prepared.selection
            ? executeBasisConversionPlan(
                  *prepared.selection,
                  prepared.expression,
                  combinatorialTags,
                  prepared.facts)
            : prepared.expression;
        if (error()) return zero();
        if (prepared.selection &&
            !checkAllApplicableBasisConversionPlans(
                prepared.expression,
                prepared.selection->sourceBasisId,
                prepared.selection->targetBasisId,
                combinatorialTags,
                prepared.facts,
                converted))
          return zero();
        for (const auto& convertedTerm : polyValue(converted)->terms)
          appendTermIfNonZero(
              resultTerms,
              convertedTerm.coeff,
              convertedTerm.monomial);
      }

    // Source groups can overlap in the target. Merge them once, after every
    // independently selected complete plan has executed.
    ring_elem result = fromTermVector(resultTerms, false);
    ExpressionFacts resultFacts =
        inferCanonicalExpansionFacts(result, targetBasisId);
    if (!resultFacts.canonicalExpansionInBasis(targetBasisId))
      {
        ERROR("the canonical group-conversion helper did not produce "
              "a target-basis expansion");
        return zero();
      }
    mutablePolyValue(result)->combinatorialTags = combinatorialTags;
    attachExpressionFacts(
        result,
        resultFacts,
        targetBasisId,
        combinatorialTags);
    return result;
  }

// ============================================================================
// Multiplication Plan Catalog
// ============================================================================
// Plans are listed from most structured to broadest. The final power-sum plan
// is the independent correctness path for every built-in operand pair.
//
// Multiplication occupies the remainder of this file before the public
// toBasis definition because product resolution in toBasis delegates to these
// helpers. This dependency order keeps the owning workflow acyclic while the
// block headers keep conversion and multiplication topics distinct.

bool SymmetricEngineRing::isSchurCompatibleBasisKind(BasisKind kind)
{
    return kind == BasisKind::Schur ||
           kind == BasisKind::Complete ||
           kind == BasisKind::Elementary ||
           kind == BasisKind::PowerSum;
  }

bool SymmetricEngineRing::isMonomialLikeCompatibleBasisKind(
    BasisKind kind,
    BasisKind target)
{
    return kind == target ||
           kind == BasisKind::Complete ||
           kind == BasisKind::Elementary ||
           kind == BasisKind::PowerSum;
  }

std::string_view SymmetricEngineRing::multiplicationPlanIdentifier(
    MultiplicationPlanId id)
{
    switch (id)
      {
        case MultiplicationPlanId::Unavailable:
          return "unavailable";
        case MultiplicationPlanId::SchurLittlewoodRichardson:
          return "Schur-product:littlewood-richardson";
        case MultiplicationPlanId::SchurHorizontalPieri:
          return "Schur-product:horizontal-pieri";
        case MultiplicationPlanId::SchurVerticalPieri:
          return "Schur-product:vertical-pieri";
        case MultiplicationPlanId::SchurBorderStrips:
          return "Schur-product:border-strips";
        case MultiplicationPlanId::SchurCompatibleFactorRules:
          return "Schur-product:compatible-factor-rules";
        case MultiplicationPlanId::MonomialLikeExponentSplittings:
          return "monomial-like-product:exponent-splittings";
        case MultiplicationPlanId::HallLittlewoodGenerators:
          return "Hall-Littlewood-product:generators";
        case MultiplicationPlanId::MultiplicativeTarget:
          return "product:multiplicative-target";
        case MultiplicationPlanId::PowerSums:
          return "product:power-sums";
      }
    return "unavailable";
  }

const std::vector<SymmetricEngineRing::MultiplicationPlan>&
SymmetricEngineRing::registeredMultiplicationPlans(
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId) const
{
    const int leftBasisId = leftFacts.expandedBasis.value_or(-1);
    const int rightBasisId = rightFacts.expandedBasis.value_or(-1);
    const std::tuple<int, int, int> key{
        leftBasisId, rightBasisId, targetBasisId};
    auto found = multiplicationPlanRegistryCache.find(key);
    if (found != multiplicationPlanRegistryCache.end())
      return found->second;
    auto inserted = multiplicationPlanRegistryCache.emplace(
        key,
        buildMultiplicationPlans(
            leftFacts, rightFacts, targetBasisId));
    return inserted.first->second;
  }

std::vector<SymmetricEngineRing::MultiplicationPlan>
SymmetricEngineRing::buildMultiplicationPlans(
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId) const
{
    std::vector<MultiplicationPlan> plans;
    int powerSumBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return plans;
    BasisKind targetKind = basisKindForId(targetBasisId);
    std::optional<BasisKind> leftKind;
    std::optional<BasisKind> rightKind;
    if (leftFacts.expandedBasis)
      leftKind = basisKindForId(*leftFacts.expandedBasis);
    if (rightFacts.expandedBasis)
      rightKind = basisKindForId(*rightFacts.expandedBasis);

    if (targetKind == BasisKind::Schur && leftKind && rightKind)
      {
        MultiplicationPlanId planId =
            MultiplicationPlanId::SchurLittlewoodRichardson;
        int specializedKinds =
            static_cast<int>(*leftKind == BasisKind::Complete ||
                             *rightKind == BasisKind::Complete) +
            static_cast<int>(*leftKind == BasisKind::Elementary ||
                             *rightKind == BasisKind::Elementary) +
            static_cast<int>(*leftKind == BasisKind::PowerSum ||
                             *rightKind == BasisKind::PowerSum);
        if (specializedKinds > 1)
          planId =
              MultiplicationPlanId::SchurCompatibleFactorRules;
        else if (*leftKind == BasisKind::Complete ||
                 *rightKind == BasisKind::Complete)
          planId = MultiplicationPlanId::SchurHorizontalPieri;
        else if (*leftKind == BasisKind::Elementary ||
                 *rightKind == BasisKind::Elementary)
          planId = MultiplicationPlanId::SchurVerticalPieri;
        else if (*leftKind == BasisKind::PowerSum ||
                 *rightKind == BasisKind::PowerSum)
          planId = MultiplicationPlanId::SchurBorderStrips;
        plans.push_back({
            planId,
            MultiplicationKernel::SchurCompatibleFactors,
            -1,
            -1,
            targetBasisId,
            true});
      }

    if ((targetKind == BasisKind::Monomial ||
         targetKind == BasisKind::Forgotten))
      plans.push_back({
          MultiplicationPlanId::MonomialLikeExponentSplittings,
          MultiplicationKernel::MonomialLikeExpansion,
          -1,
          -1,
          targetBasisId,
          true});

    if (isHallLittlewoodCapitalBasisKind(targetKind))
      {
        BasisKind generatorKind =
            targetKind == BasisKind::HallLittlewoodQ ||
                    targetKind == BasisKind::HallLittlewoodP
                ? BasisKind::HallLittlewoodQGenerator
                : BasisKind::HallLittlewoodBGenerator;
        int generatorId = requiredBasisIdForKind(generatorKind);
        if (error()) return plans;
        plans.push_back({
            MultiplicationPlanId::HallLittlewoodGenerators,
            MultiplicationKernel::CanonicalBasisProduct,
            generatorId,
            generatorId,
            generatorId,
            true});
      }

    if (isMultiplicativeBasis(targetBasisId))
      plans.push_back({
          MultiplicationPlanId::MultiplicativeTarget,
          MultiplicationKernel::CanonicalBasisProduct,
          targetBasisId,
          targetBasisId,
          targetBasisId,
          true});
    plans.push_back({
        MultiplicationPlanId::PowerSums,
        MultiplicationKernel::CanonicalBasisProduct,
        powerSumBasisId,
        powerSumBasisId,
        powerSumBasisId,
        true});
    return plans;
  }

// ============================================================================
// Multiplication-Plan Applicability Contracts
// ============================================================================

bool SymmetricEngineRing::multiplicationPlanApplicable(
    const MultiplicationPlan& plan,
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId) const
{
    if (!leftFacts.singleBasisElementId ||
        !rightFacts.singleBasisElementId ||
        !leftFacts.singleBasisElementCoefficientOne.value_or(false) ||
        !rightFacts.singleBasisElementCoefficientOne.value_or(false) ||
        !leftFacts.canonicalExpansionInBasis(
            *leftFacts.singleBasisElementId) ||
        !rightFacts.canonicalExpansionInBasis(
            *rightFacts.singleBasisElementId))
      return false;
    const bool schurPlan =
        plan.id == MultiplicationPlanId::SchurLittlewoodRichardson ||
        plan.id == MultiplicationPlanId::SchurHorizontalPieri ||
        plan.id == MultiplicationPlanId::SchurVerticalPieri ||
        plan.id == MultiplicationPlanId::SchurBorderStrips ||
        plan.id == MultiplicationPlanId::SchurCompatibleFactorRules;
    if ((schurPlan &&
         plan.kernel != MultiplicationKernel::SchurCompatibleFactors) ||
        (plan.id ==
             MultiplicationPlanId::MonomialLikeExponentSplittings &&
         plan.kernel != MultiplicationKernel::MonomialLikeExpansion) ||
        ((plan.id == MultiplicationPlanId::HallLittlewoodGenerators ||
          plan.id == MultiplicationPlanId::MultiplicativeTarget ||
          plan.id == MultiplicationPlanId::PowerSums) &&
         plan.kernel != MultiplicationKernel::CanonicalBasisProduct) ||
        plan.id == MultiplicationPlanId::Unavailable)
      return false;
    if (plan.kernel == MultiplicationKernel::SchurCompatibleFactors)
      {
        if (!leftFacts.expandedBasis || !rightFacts.expandedBasis)
          return false;
        BasisKind leftKind = basisKindForId(*leftFacts.expandedBasis);
        BasisKind rightKind = basisKindForId(*rightFacts.expandedBasis);
        BasisKind targetKind = basisKindForId(targetBasisId);
        return plan.kernelOutputBasisId == targetBasisId &&
               targetKind == BasisKind::Schur &&
               isSchurCompatibleBasisKind(leftKind) &&
               isSchurCompatibleBasisKind(rightKind) &&
               leftFacts.skewFree && rightFacts.skewFree;
      }
    if (plan.kernel == MultiplicationKernel::MonomialLikeExpansion)
      {
        if (!leftFacts.expandedBasis || !rightFacts.expandedBasis)
          return false;
        BasisKind leftKind = basisKindForId(*leftFacts.expandedBasis);
        BasisKind rightKind = basisKindForId(*rightFacts.expandedBasis);
        BasisKind targetKind = basisKindForId(targetBasisId);
        return plan.kernelOutputBasisId == targetBasisId &&
               (targetKind == BasisKind::Monomial ||
                targetKind == BasisKind::Forgotten) &&
               isMonomialLikeCompatibleBasisKind(
                   leftKind, targetKind) &&
               isMonomialLikeCompatibleBasisKind(
                   rightKind, targetKind) &&
               leftFacts.skewFree && rightFacts.skewFree;
      }
    if (plan.kernel != MultiplicationKernel::CanonicalBasisProduct ||
        !plan.kernelOutputCanonical ||
        plan.leftKernelBasisId <= 0 ||
        plan.rightKernelBasisId <= 0)
      return false;
    switch (plan.id)
      {
        case MultiplicationPlanId::HallLittlewoodGenerators:
          {
            if (!leftFacts.expandedBasis || !rightFacts.expandedBasis)
              return false;
            const BasisKind targetKind = basisKindForId(targetBasisId);
            const BasisKind leftKind =
                basisKindForId(*leftFacts.expandedBasis);
            const BasisKind rightKind =
                basisKindForId(*rightFacts.expandedBasis);
            return isHallLittlewoodCapitalBasisKind(targetKind) &&
                   leftKind != BasisKind::Custom &&
                   rightKind != BasisKind::Custom;
          }
        case MultiplicationPlanId::MultiplicativeTarget:
          return plan.leftKernelBasisId == targetBasisId &&
                 plan.rightKernelBasisId == targetBasisId &&
                 plan.kernelOutputBasisId == targetBasisId &&
                 isMultiplicativeBasis(targetBasisId);
        case MultiplicationPlanId::PowerSums:
          {
            const int powerSumBasisId =
                requiredBasisIdForKind(BasisKind::PowerSum);
            return !error() &&
                   plan.leftKernelBasisId == powerSumBasisId &&
                   plan.rightKernelBasisId == powerSumBasisId &&
                   plan.kernelOutputBasisId == powerSumBasisId;
          }
        case MultiplicationPlanId::Unavailable:
        case MultiplicationPlanId::SchurLittlewoodRichardson:
        case MultiplicationPlanId::SchurHorizontalPieri:
        case MultiplicationPlanId::SchurVerticalPieri:
        case MultiplicationPlanId::SchurBorderStrips:
        case MultiplicationPlanId::SchurCompatibleFactorRules:
        case MultiplicationPlanId::MonomialLikeExponentSplittings:
          return false;
      }
    return false;
  }

// ============================================================================
// Multiplication-Plan Picker
// ============================================================================
// Operand conversions are selected here as complete plans together with the
// product kernel. The later post-kernel conversion is intentionally selected
// by the owning workflow only after the kernel's realized output has been
// normalized and inspected.

SymmetricEngineRing::MultiplicationPlanSelection
SymmetricEngineRing::selectMultiplicationPlan(
    ring_elem f,
    ring_elem g,
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId) const
{
    const auto& plans =
        registeredMultiplicationPlans(
            leftFacts, rightFacts, targetBasisId);
    auto completeSelection =
        [&](const MultiplicationPlan& plan) {
          MultiplicationPlanSelection selection{
              plan, std::nullopt, std::nullopt};
          if (plan.leftKernelBasisId > 0 &&
              *leftFacts.expandedBasis != plan.leftKernelBasisId)
            selection.leftConversion =
                selectBasisConversionPlan(
                    f,
                    *leftFacts.expandedBasis,
                    plan.leftKernelBasisId,
                    leftFacts,
                    0);
          if (error()) return selection;
          if (plan.rightKernelBasisId > 0 &&
              *rightFacts.expandedBasis != plan.rightKernelBasisId)
            selection.rightConversion =
                selectBasisConversionPlan(
                    g,
                    *rightFacts.expandedBasis,
                    plan.rightKernelBasisId,
                    rightFacts,
                    0);
          return selection;
        };
    const char *forced =
        std::getenv("M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN");
    if (forced != nullptr)
      {
        for (const auto& plan : plans)
          if (multiplicationPlanIdentifier(plan.id) == forced &&
              multiplicationPlanApplicable(
                  plan, leftFacts, rightFacts, targetBasisId))
            return completeSelection(plan);
        ERROR("the requested multiplication plan is unknown or inapplicable: ",
              forced);
        int powerSumId = requiredBasisIdForKind(BasisKind::PowerSum);
        return {{MultiplicationPlanId::Unavailable,
                 MultiplicationKernel::CanonicalBasisProduct,
                 powerSumId,
                 powerSumId,
                 powerSumId,
                 true},
                std::nullopt,
                std::nullopt};
      }

    for (const auto& plan : plans)
      if (multiplicationPlanApplicable(
              plan, leftFacts, rightFacts, targetBasisId))
        return completeSelection(plan);
    ERROR("the multiplication registry has no applicable plan");
    int powerSumId = requiredBasisIdForKind(BasisKind::PowerSum);
    return {{MultiplicationPlanId::Unavailable,
             MultiplicationKernel::CanonicalBasisProduct,
             powerSumId,
             powerSumId,
             powerSumId,
             true},
            std::nullopt,
            std::nullopt};
  }

// ============================================================================
// Policy-Free Multiplication Kernel Executor
// ============================================================================
// Try-style mathematical helpers revalidate the selected plan's structural
// preconditions. A false result here is a contract violation, not a signal to
// select a fallback plan.

ring_elem SymmetricEngineRing::executeMultiplicationKernel(
    const MultiplicationPlan& plan,
    ring_elem f,
    ring_elem g,
    int targetBasisId) const
{
    ring_elem kernelResult;
    if (plan.kernel == MultiplicationKernel::SchurCompatibleFactors)
      {
        bool completed = tryProductToSchurViaCompatibleFactors(
            f,
            g,
            targetBasisId,
            displayForBasis(targetBasisId),
            basisOrderForId(targetBasisId),
            kernelResult);
        if (error()) return zero();
        if (!completed)
          {
            ERROR("a selected multiplication plan violated its "
                  "applicability contract");
            return zero();
          }
      }
    else if (plan.kernel ==
             MultiplicationKernel::MonomialLikeExpansion)
      {
        bool completed = tryProductToMonomialLikeTarget(
            f,
            g,
            targetBasisId,
            displayForBasis(targetBasisId),
            basisOrderForId(targetBasisId),
            isMultiplicativeBasis(targetBasisId),
            kernelResult);
        if (error()) return zero();
        if (!completed)
          {
            ERROR("a selected multiplication kernel violated its "
                  "applicability contract");
            return zero();
          }
      }
    else
      kernelResult = mult(f, g);
    return error() ? zero() : kernelResult;
  }

// ============================================================================
// Binary Multiplication Workflow
// ============================================================================
// Stages are deliberately visible and linear: convert operands, execute the
// selected kernel, normalize its declared output, select and execute one
// complete post-kernel conversion plan, and validate the target contract.

ring_elem SymmetricEngineRing::runMultiplicationWorkflow(
    const MultiplicationPlanSelection& selection,
    ring_elem f,
    ring_elem g,
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId,
    bool attachResultFacts) const
{
    const bool trace =
        std::getenv(
            "M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr;
    const MultiplicationPlan& plan = selection.plan;
    CombinatorialTags tags = selectMultiplicationTags(f, g);
    ring_elem left = f;
    ring_elem right = g;
    if (selection.leftConversion)
      {
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-stage: "
              "stage=convert-left plan=%s\n",
              selection.leftConversion->
                  definition->id.value.c_str());
        left = executeBasisConversionPlan(
            *selection.leftConversion,
            f,
            polyValue(f)->combinatorialTags,
            leftFacts);
        if (error()) return zero();
      }
    if (selection.rightConversion)
      {
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-stage: "
              "stage=convert-right plan=%s\n",
              selection.rightConversion->
                  definition->id.value.c_str());
        right = executeBasisConversionPlan(
            *selection.rightConversion,
            g,
            polyValue(g)->combinatorialTags,
            rightFacts);
        if (error()) return zero();
      }

    if (trace)
      std::fprintf(
          stderr,
          "SymmetricRings multiplication-stage: stage=kernel "
          "plan=%s output-basis=%s canonical=%s\n",
          multiplicationPlanIdentifier(plan.id).data(),
          basisKeyForId(plan.kernelOutputBasisId).c_str(),
          plan.kernelOutputCanonical ? "yes" : "no");
    ring_elem kernelResult = executeMultiplicationKernel(
        plan, left, right, targetBasisId);
    if (error()) return zero();

    ExpressionFacts kernelFacts;
    ring_elem canonicalKernelResult;
    if (plan.kernelOutputCanonical)
      {
        canonicalKernelResult = kernelResult;
        // Exact kernel support is necessary only for a post-kernel picker or
        // for public metadata. A transient target-basis pair product can rely
        // on the registered kernel's canonical-output contract.
        if (plan.kernelOutputBasisId != targetBasisId ||
            attachResultFacts)
          {
            auto metadataFacts =
                expressionFactsFromMetadata(canonicalKernelResult);
            kernelFacts = metadataFacts
                ? *metadataFacts
                : inferCanonicalExpansionFacts(
                      canonicalKernelResult,
                      plan.kernelOutputBasisId);
          }
      }
    else
      canonicalKernelResult =
          normalizeExpression(kernelResult, kernelFacts);
    if (error()) return zero();
    if ((plan.kernelOutputBasisId != targetBasisId ||
         attachResultFacts) &&
        !kernelFacts.canonicalExpansionInBasis(
            plan.kernelOutputBasisId))
      {
        ERROR("a multiplication kernel violated its declared output contract");
        return zero();
      }

    mutablePolyValue(canonicalKernelResult)->combinatorialTags = tags;
    if (plan.kernelOutputBasisId == targetBasisId &&
        !attachResultFacts)
      return canonicalKernelResult;
    kernelFacts.combinatorialTags = tags;
    ring_elem result = canonicalKernelResult;
    ExpressionFacts resultFacts;
    if (plan.kernelOutputBasisId != targetBasisId)
      {
        ResolvedBasisConversionPlan postKernelSelection =
            selectBasisConversionPlan(
                canonicalKernelResult,
                plan.kernelOutputBasisId,
                targetBasisId,
                kernelFacts,
                tags,
                std::nullopt,
                &kernelFacts);
        if (error()) return zero();
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-stage: "
              "stage=post-kernel-conversion plan=%s\n",
              postKernelSelection.
                  definition->id.value.c_str());
        result = executeBasisConversionPlan(
            postKernelSelection,
            canonicalKernelResult,
            tags,
            kernelFacts,
            attachResultFacts ? &resultFacts : nullptr);
      }
    else if (attachResultFacts)
      resultFacts = std::move(kernelFacts);
    if (error()) return zero();

    if (!attachResultFacts) return result;
    if (!resultFacts.canonicalExpansionInBasis(targetBasisId))
      {
        ERROR("the binary multiplication workflow did not produce "
              "a canonical target-basis expansion");
        return zero();
      }
    attachExpressionFacts(
        result,
        resultFacts,
        targetBasisId,
        tags);
    return result;
  }

// ============================================================================
// Strict Basis-Element Multiplication Helper
// ============================================================================
// Coefficients are distributed by the outer wrapper. This helper therefore
// accepts exactly two coefficient-one, normalized, non-skew basis elements.

ring_elem SymmetricEngineRing::multiplyBasisElementsWithFacts(
    ring_elem f,
    ring_elem g,
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId,
    bool attachResultFacts) const
{
    MultiplicationPlanSelection selection =
        selectMultiplicationPlan(
            f,
            g,
            leftFacts, rightFacts, targetBasisId);
    if (error()) return zero();
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr)
      std::fprintf(
          stderr,
          "SymmetricRings multiplication-plan: target=%s plan=%s\n",
          displayForBasis(targetBasisId).c_str(),
          multiplicationPlanIdentifier(selection.plan.id).data());
    return runMultiplicationWorkflow(
        selection,
        f,
        g,
        leftFacts,
        rightFacts,
        targetBasisId,
        attachResultFacts);
}

// ============================================================================
// One-Term Product Resolver
// ============================================================================
// Canonical storage has already removed zero coefficients. This helper removes
// identity factors, applies the surviving coefficient once at the end, and
// uses the multiplicative-target or pairwise branch prescribed by the design.

SymmetricEngineRing::ResolvedProductTerm
SymmetricEngineRing::multiplyTermToBasis(
    const SymmetricTerm& term,
    int targetBasisId) const
{
    const bool trace =
        std::getenv(
            "M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr;
    struct ProductFactor
    {
      ring_elem expression;
      ExpressionFacts facts;
    };
    std::vector<ProductFactor> factors;
    size_t pos = 0;
    while (pos < term.monomial.data.size())
      {
        bool identityFactor = false;
        if (!atomIsSkewAt(term.monomial, pos))
          {
            // The input has already been normalized. Test the encoded index
            // in place instead of allocating a Partition merely to trim
            // trailing zeroes and ask whether it is empty.
            identityFactor = true;
            const size_t indexLength =
                static_cast<size_t>(
                    atomIndexLengthAt(term.monomial, pos));
            for (size_t i = 0; i < indexLength; ++i)
              if (term.monomial.data[
                      pos + atomHeaderSize + i] != 0)
                {
                  identityFactor = false;
                  break;
                }
          }
        if (!identityFactor)
          factors.push_back({
              expressionFromAtom(term.monomial, pos),
              basisElementFactsFromAtom(
                  term.monomial,
                  pos,
                  coefficientRing->one())});
        pos += atomLengthAt(term.monomial, pos);
      }

    // Factor conversions are internal stage values. Keep their exact facts in
    // the workflow instead of serializing a full metadata record onto each
    // temporary polynomial; toBasis attaches only the final public result.
    auto convertFactor =
        [&](const ProductFactor& factor) {
          if (!factor.facts.expandedBasis)
            {
              ERROR("a product factor was not a canonical basis element");
              return zero();
            }
          const int sourceBasisId = *factor.facts.expandedBasis;
          CombinatorialTags factorTags =
              polyValue(factor.expression)->combinatorialTags;
          if (sourceBasisId == targetBasisId)
            return factor.expression;
          ResolvedBasisConversionPlan selection =
              selectBasisConversionPlan(
                  factor.expression,
                  sourceBasisId,
                  targetBasisId,
                  factor.facts,
                  factorTags);
          if (error()) return zero();
          return executeBasisConversionPlan(
              selection,
              factor.expression,
              factorTags,
              factor.facts);
        };

    ring_elem result;
    if (factors.empty())
      {
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings product-resolution: "
              "branch=scalar target=%s\n",
              basisKeyForId(targetBasisId).c_str());
        result = one();
      }
    else if (factors.size() == 1)
      {
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings product-resolution: "
              "branch=single-factor target=%s\n",
              basisKeyForId(targetBasisId).c_str());
        result = convertFactor(factors.front());
      }
    else if (isMultiplicativeBasis(targetBasisId))
      {
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings product-resolution: "
              "branch=multiplicative-target factors=%zu target=%s\n",
              factors.size(),
              basisKeyForId(targetBasisId).c_str());
        std::vector<ring_elem> convertedFactors;
        convertedFactors.reserve(factors.size());
        for (const auto& factor : factors)
          {
            ring_elem converted = convertFactor(factor);
            if (error()) return {zero(), {}};
            convertedFactors.push_back(converted);
          }
        // A balanced fold limits the size disparity of intermediate products
        // and avoids multiplying every expansion through a steadily growing
        // accumulator.
        while (convertedFactors.size() > 1)
          {
            std::vector<ring_elem> nextLevel;
            nextLevel.reserve((convertedFactors.size() + 1) / 2);
            for (size_t i = 0; i < convertedFactors.size(); i += 2)
              {
                if (i + 1 == convertedFactors.size())
                  nextLevel.push_back(convertedFactors[i]);
                else
                  {
                    ring_elem product = mult(
                        convertedFactors[i], convertedFactors[i + 1]);
                    if (error()) return {zero(), {}};
                    nextLevel.push_back(product);
                  }
              }
            convertedFactors = std::move(nextLevel);
          }
        result = convertedFactors.front();
      }
    else
      {
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings product-resolution: "
              "branch=pairwise-plans factors=%zu target=%s\n",
              factors.size(),
              basisKeyForId(targetBasisId).c_str());
        result = multiplyBasisElementsWithFacts(
            factors[0].expression,
            factors[1].expression,
            factors[0].facts,
            factors[1].facts,
            targetBasisId,
            true);
        for (size_t i = 2; i < factors.size() && !error(); ++i)
          result = multiplyToBasis(
              result, factors[i].expression, targetBasisId);
      }
    if (error()) return {zero(), {}};
    // Coefficients are immutable ring elements. Avoid walking and copying a
    // potentially large expansion when the product term already has unit
    // coefficient.
    if (!coefficientRing->is_equal(
            term.coeff, coefficientRing->one()))
      result = scaled(term.coeff, result);
    if (error()) return {zero(), {}};

    auto metadataFacts = expressionFactsFromMetadata(result);
    ExpressionFacts resultFacts = metadataFacts
        ? *metadataFacts
        : inferCanonicalExpansionFacts(
              result, targetBasisId);
    if (!resultFacts.canonicalExpansionInBasis(targetBasisId))
      {
        ERROR("a resolved product term violated its canonical "
              "target-basis contract");
        return {zero(), {}};
      }
    resultFacts.combinatorialTags =
        polyValue(result)->combinatorialTags;
    if (!metadataFacts)
      attachExpressionFacts(
          result,
          resultFacts,
          targetBasisId,
          resultFacts.combinatorialTags);
    return {result, std::move(resultFacts)};
  }

// ============================================================================
// Distributive Multiplication Entry Workflow
// ============================================================================
// The public API accepts product-free linear combinations. It distributes
// coefficients and delegates each nonscalar pair to the strict helper above.

ring_elem SymmetricEngineRing::multiplyToBasis(
    ring_elem f,
    ring_elem g,
    int targetBasisId) const
{
    requireBasis(targetBasisId);
    if (error()) return zero();
    auto strictInputFacts =
        [&](ring_elem operand)
            -> std::optional<ExpressionFacts> {
          const auto *poly = polyValue(operand);
          if (poly->terms.size() != 1 ||
              poly->terms.front().monomial.data.empty() ||
              atomLengthAt(poly->terms.front().monomial, 0) !=
                  poly->terms.front().monomial.data.size())
            return std::nullopt;
          ExpressionFacts facts = basisElementFactsFromAtom(
              poly->terms.front().monomial,
              0,
              poly->terms.front().coeff);
          facts.combinatorialTags = poly->combinatorialTags;
          if (!facts.singleBasisElementId ||
              !facts.singleBasisElementCoefficientOne.value_or(false) ||
              !facts.canonicalExpansionInBasis(
                  *facts.singleBasisElementId))
            return std::nullopt;
          return facts;
        };
    auto strictLeftFacts = strictInputFacts(f);
    auto strictRightFacts = strictInputFacts(g);
    if (strictLeftFacts && strictRightFacts)
      {
        MultiplicationPlanSelection selection =
            selectMultiplicationPlan(
                f,
                g,
                *strictLeftFacts, *strictRightFacts, targetBasisId);
        if (error()) return zero();
        return runMultiplicationWorkflow(
            selection,
            f,
            g,
            *strictLeftFacts,
            *strictRightFacts,
            targetBasisId,
            true);
      }

    ExpressionFacts inferredLeftFacts;
    auto knownLeftFacts = expressionFactsFromMetadata(f);
    ring_elem left = f;
    if (!knownLeftFacts)
      left = normalizeExpression(f, inferredLeftFacts);
    if (error()) return zero();
    const ExpressionFacts& leftFacts = knownLeftFacts
        ? *knownLeftFacts
        : inferredLeftFacts;

    ExpressionFacts inferredRightFacts;
    auto knownRightFacts = expressionFactsFromMetadata(g);
    ring_elem right = g;
    if (!knownRightFacts)
      right = normalizeExpression(g, inferredRightFacts);
    if (error()) return zero();
    const ExpressionFacts& rightFacts = knownRightFacts
        ? *knownRightFacts
        : inferredRightFacts;
    if (!leftFacts.noProducts() || !rightFacts.noProducts())
      {
        ERROR("multiplyToBasis expects product-free operands");
        return zero();
      }

    CombinatorialTags tags = selectMultiplicationTags(left, right);
    if (leftFacts.singleBasisElement() &&
        rightFacts.singleBasisElement() &&
        leftFacts.singleBasisElementCoefficientOne.value_or(false) &&
        rightFacts.singleBasisElementCoefficientOne.value_or(false))
      {
        // The public entry point most often receives exactly the strict kernel
        // contract. Reuse the facts established above instead
        // of extracting two temporary atoms and inferring them again.
        MultiplicationPlanSelection selection =
            selectMultiplicationPlan(
                left,
                right,
                leftFacts, rightFacts, targetBasisId);
        if (error()) return zero();
        return runMultiplicationWorkflow(
            selection,
            left,
            right,
            leftFacts,
            rightFacts,
            targetBasisId,
            true);
      }

    // Materialize and inspect each canonical atom once. In the distributive
    // cross product a left atom is otherwise rebuilt and rescanned once per
    // right term (and conversely), even though its facts are immutable.
    struct PreparedOperandTerm
    {
      const SymmetricTerm *term;
      bool scalar;
      ring_elem factor;
      ExpressionFacts facts;
    };
    auto prepareOperand =
        [&](ring_elem operand) {
          std::vector<PreparedOperandTerm> prepared;
          prepared.reserve(polyValue(operand)->terms.size());
          for (const auto& term : polyValue(operand)->terms)
            {
              bool scalar = term.monomial.data.empty();
              ring_elem factor = scalar
                  ? one()
                  : expressionFromAtom(term.monomial, 0);
              ExpressionFacts facts;
              if (!scalar)
                facts = basisElementFactsFromAtom(
                    term.monomial, 0, coefficientRing->one());
              prepared.push_back(
                  {&term, scalar, factor, std::move(facts)});
            }
          return prepared;
        };
    std::vector<PreparedOperandTerm> preparedLeft =
        prepareOperand(left);
    std::vector<PreparedOperandTerm> preparedRight =
        prepareOperand(right);

    auto convertTransientFactor =
        [&](const PreparedOperandTerm& prepared) {
          if (!prepared.facts.expandedBasis)
            {
              ERROR("a prepared multiplication atom has no source basis");
              return zero();
            }
          const int sourceBasisId = *prepared.facts.expandedBasis;
          if (sourceBasisId == targetBasisId)
            return prepared.factor;
          ResolvedBasisConversionPlan conversion =
              selectBasisConversionPlan(
                  prepared.factor,
                  sourceBasisId,
                  targetBasisId,
                  prepared.facts,
                  tags);
          if (error()) return zero();
          return executeBasisConversionPlan(
              conversion,
              prepared.factor,
              tags,
              prepared.facts);
        };

    // Structured kernels need no operand conversion plans, and their selection
    // depends only on the two source bases and target. Reuse that completed
    // selection across equal-basis pairs in this distributive request.
    std::map<
        std::tuple<int, int, int>,
        MultiplicationPlanSelection>
        reusablePairSelections;
    VECTOR(SymmetricTerm) resultTerms;
    for (const auto& leftTerm : preparedLeft)
      for (const auto& rightTerm : preparedRight)
        {
          ring_elem coefficient =
              coefficientRing->mult(
                  leftTerm.term->coeff, rightTerm.term->coeff);
          if (coefficientRing->is_zero(coefficient)) continue;

          ring_elem product;
          if (leftTerm.scalar && rightTerm.scalar)
            product = one();
          else if (leftTerm.scalar || rightTerm.scalar)
            {
              product = convertTransientFactor(
                  leftTerm.scalar ? rightTerm : leftTerm);
            }
          else
            {
              const auto key = std::make_tuple(
                  *leftTerm.facts.expandedBasis,
                  *rightTerm.facts.expandedBasis,
                  targetBasisId);
              auto reusable = reusablePairSelections.find(key);
              if (reusable != reusablePairSelections.end())
                product = runMultiplicationWorkflow(
                    reusable->second,
                    leftTerm.factor,
                    rightTerm.factor,
                    leftTerm.facts,
                    rightTerm.facts,
                    targetBasisId,
                    false);
              else
                {
                  MultiplicationPlanSelection selection =
                      selectMultiplicationPlan(
                          leftTerm.factor,
                          rightTerm.factor,
                          leftTerm.facts,
                          rightTerm.facts,
                          targetBasisId);
                  if (error()) return zero();
                  if (!selection.leftConversion &&
                      !selection.rightConversion)
                    reusablePairSelections.emplace(key, selection);
                  product = runMultiplicationWorkflow(
                      selection,
                      leftTerm.factor,
                      rightTerm.factor,
                      leftTerm.facts,
                      rightTerm.facts,
                      targetBasisId,
                      false);
                }
            }
          if (error()) return zero();
          for (const auto& productTerm : polyValue(product)->terms)
            appendTermIfNonZero(
                resultTerms,
                coefficientRing->mult(
                    coefficient, productTerm.coeff),
                productTerm.monomial);
        }

    // Distribution can create many overlapping target terms. Collect once,
    // after every strict binary workflow has finished, instead of repeatedly
    // copying and merging a growing polynomial after each operand pair.
    ring_elem result = fromTermVector(resultTerms, false);
    ExpressionFacts resultFacts =
        inferCanonicalExpansionFacts(result, targetBasisId);
    attachExpressionFacts(
        result, resultFacts, targetBasisId, tags);
    return result;
  }

// ============================================================================
// Unified Basis-Conversion Entry Workflow
// ============================================================================
// This is the sole owning engine workflow for built-in conversion. Its stages
// are: exact-metadata bypass or normalization, product resolution, source-basis
// grouping, complete-plan selection and execution, and collection.

ring_elem SymmetricEngineRing::toBasis(
    ring_elem f,
    int targetBasisId) const
{
    const bool trace =
        std::getenv(
            "M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr;
    requireBasis(targetBasisId);
    if (error()) return zero();
    CombinatorialTags combinatorialTags =
        polyValue(f)->combinatorialTags;
    const auto *inputPoly = polyValue(f);
    if (inputPoly->terms.size() == 1 &&
        !inputPoly->terms.front().monomial.data.empty() &&
        atomLengthAt(inputPoly->terms.front().monomial, 0) ==
            inputPoly->terms.front().monomial.data.size())
      {
        ExpressionFacts strictFacts =
            basisElementFactsFromAtom(
                inputPoly->terms.front().monomial,
                0,
                inputPoly->terms.front().coeff);
        strictFacts.combinatorialTags = combinatorialTags;
        if (strictFacts.singleBasisElementId &&
            strictFacts.canonicalExpansionInBasis(
                *strictFacts.singleBasisElementId))
          return convertCanonicalExpressionToBasis(
              f,
              targetBasisId,
              combinatorialTags,
              &strictFacts);
      }
    auto knownFacts = expressionFactsFromMetadata(f);
    if (knownFacts)
      {
        if (trace)
          std::fprintf(stderr,
                       "SymmetricRings toBasis: "
                       "bypass=canonical-metadata target=%s\n",
                       displayForBasis(targetBasisId).c_str());
        return convertCanonicalExpressionToBasis(
            f, targetBasisId, combinatorialTags, &*knownFacts);
      }

    // Stage 1: establish normalized, non-skew factors.
    ExpressionFacts normalizedFacts;
    std::vector<size_t> normalizedFactorCounts;
    ring_elem normalized =
        normalizeExpression(
            f, normalizedFacts, &normalizedFactorCounts);
    if (error()) return zero();

    if (trace)
      {
        std::fprintf(stderr,
                     "SymmetricRings toBasis: target=%s terms=%zu "
                     "products=%zu bases=%zu\n",
                     displayForBasis(targetBasisId).c_str(),
                     normalizedFacts.termCount,
                     normalizedFacts.productTermCount,
                     normalizedFacts.factorBases.size());
      }

    // Normalization already produced the exact canonical facts needed by the
    // group helper.  Preserve that state instead of rebuilding and rescanning
    // the overwhelmingly common product-free input.
    if (normalizedFacts.noProducts())
      return convertCanonicalExpressionToBasis(
          normalized,
          targetBasisId,
          combinatorialTags,
          &normalizedFacts);

    // Stage 2: resolve every remaining multifactor term into the target.
    const auto& normalizedTerms = polyValue(normalized)->terms;
    if (normalizedFactorCounts.size() != normalizedTerms.size())
      {
        ERROR("normalized factor-count profile is inconsistent");
        return zero();
      }
    if (normalizedTerms.size() == 1 &&
        normalizedFactorCounts.front() > 1)
      {
        ResolvedProductTerm resolved =
            multiplyTermToBasis(
                normalizedTerms.front(),
                targetBasisId);
        if (error()) return zero();
        resolved.facts.combinatorialTags =
            combinatorialTags;
        // multiplyTermToBasis already attached the same exact structural
        // facts. Only the semantic tag belongs to the owning outer operation.
        mutablePolyValue(
            resolved.expression)->combinatorialTags =
                combinatorialTags;
        return resolved.expression;
      }

    VECTOR(SymmetricTerm) preparedTerms;
    size_t resolvedProductTerms = 0;
    bool passthroughTermsAlreadyInTarget = true;
    for (size_t termPosition = 0;
         termPosition < normalizedTerms.size();
         ++termPosition)
      {
        const auto& term = normalizedTerms[termPosition];
        const size_t factorCount =
            normalizedFactorCounts[termPosition];
        if (factorCount <= 1)
          {
            preparedTerms.push_back(term);
            if (factorCount == 1 &&
                atomBasisIdAt(term.monomial, 0) != targetBasisId)
              passthroughTermsAlreadyInTarget = false;
            continue;
          }
        ResolvedProductTerm resolved =
            multiplyTermToBasis(term, targetBasisId);
        if (error()) return zero();
        for (const auto& convertedTerm :
             polyValue(resolved.expression)->terms)
          appendTermIfNonZero(
              preparedTerms,
              convertedTerm.coeff,
              convertedTerm.monomial);
        ++resolvedProductTerms;
      }

    if (trace)
      std::fprintf(
          stderr,
          "SymmetricRings conversion-stage: "
          "stage=resolve-products resolved=%zu passthrough=%zu\n",
          resolvedProductTerms,
          preparedTerms.size());

    // Stages 3--5 are owned by the canonical group helper: group by source
    // basis/weight, select every plan, then execute and collect in the target.
    ring_elem prepared = fromTermVector(preparedTerms, false);
    mutablePolyValue(prepared)->combinatorialTags = combinatorialTags;
    if (passthroughTermsAlreadyInTarget)
      {
        // Every multifactor term was resolved into the target above, and all
        // scalar/single-factor passthrough terms were already target-native.
        // The combined expression is therefore target-closed; do not regroup
        // it merely to execute an identity plan.
        ExpressionFacts preparedFacts =
            inferCanonicalExpansionFacts(
                prepared, targetBasisId);
        if (!preparedFacts.canonicalExpansionInBasis(targetBasisId))
          {
            ERROR("resolved product terms did not satisfy the target "
                  "basis contract");
            return zero();
          }
        attachExpressionFacts(
            prepared,
            preparedFacts,
            targetBasisId,
            combinatorialTags);
        return prepared;
      }
    return convertCanonicalExpressionToBasis(
        prepared, targetBasisId, combinatorialTags);
  }

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
