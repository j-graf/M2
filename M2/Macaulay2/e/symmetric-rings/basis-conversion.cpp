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
//   2. normalization and condition-piece pieceFacts;
//   3. plan indexing, structural contracts, and generic execution;
//   4. canonical source-group conversion;
//   5. multiplication plans and the three multiplication workflows;
//   6. the sole public built-in conversion workflow, toBasis.
//
// The matching declaration fragment follows the same order. Mathematical
// formulas live in basis-conversion-kernels.cpp, complete conversion plans in
// basis-conversion-plans.cpp, and endpoint policy in
// basis-conversion-picker.cpp. Product combinatorics live in
// basis-conversion-products.cpp; this file owns shared orchestration.

namespace {

// Benchmark controls are scoped to one synchronous engine request. Ordinary
// toBasis installs an empty context, so ambient diagnostic variables cannot
// force plans, enable tracing, or alter nested multiplication.
struct BasisConversionDiagnosticContext
{
  std::optional<std::string> forcedPlanIdentifier;
  bool traceConversion = false;
};

thread_local const BasisConversionDiagnosticContext
    *activeBasisConversionDiagnostics = nullptr;

class ScopedBasisConversionDiagnostics
{
 public:
  ScopedBasisConversionDiagnostics(
      const BasisConversionDiagnosticContext& context,
      bool replaceExisting)
      : mPrevious(activeBasisConversionDiagnostics)
  {
    if (replaceExisting ||
        activeBasisConversionDiagnostics == nullptr)
      {
        activeBasisConversionDiagnostics = &context;
        mInstalled = true;
      }
  }

  ~ScopedBasisConversionDiagnostics()
  {
    if (mInstalled)
      activeBasisConversionDiagnostics = mPrevious;
  }

  ScopedBasisConversionDiagnostics(
      const ScopedBasisConversionDiagnostics&) = delete;
  ScopedBasisConversionDiagnostics& operator=(
      const ScopedBasisConversionDiagnostics&) = delete;

 private:
  const BasisConversionDiagnosticContext *mPrevious;
  bool mInstalled = false;
};

const std::optional<std::string>&
forcedBasisConversionPlanIdentifier()
{
  static const std::optional<std::string> noForcedPlan;
  return activeBasisConversionDiagnostics == nullptr
      ? noForcedPlan
      : activeBasisConversionDiagnostics->
            forcedPlanIdentifier;
}

} // namespace

bool SymmetricEngineRing::basisConversionContextActive() const
{
  return activeBasisConversionDiagnostics != nullptr;
}

bool SymmetricEngineRing::basisConversionTraceEnabled() const
{
  if (activeBasisConversionDiagnostics != nullptr)
    return activeBasisConversionDiagnostics->traceConversion;
  // Other development-only operations retain their existing outer-pipeline
  // trace. Any ordinary toBasis call they make installs the empty context
  // above and therefore cannot inherit this environment setting.
  return std::getenv(
      "M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr;
}

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
    // A conversion or multiplication plan promises one normalized,
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
            facts.mostlyShortCyclePowerSumTermCount =
                powerSumIndexHasMostlyShortCycles(
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
    size_t mostlyShortCycleTerms = 0;
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
        if (powerSumIndexHasMostlyShortCycles(
                index, partitionWeight(index)))
          ++mostlyShortCycleTerms;

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
        facts.mostlyShortCyclePowerSumTermCount =
            mostlyShortCycleTerms;
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
    if (basisConversionTraceEnabled())
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
// Expression-Piece Inspection
// ============================================================================
// Conditions remain independent inspectable values. Each piece records the
// exact mathematical facts available from canonical symmetric-function
// storage.

std::vector<ExpressionPieceFacts>
SymmetricEngineRing::inspectExpressionPieces(
    ring_elem expression,
    ExpressionPieceKind pieceKind,
    const ExpressionFacts& facts) const
{
    const auto& terms = polyValue(expression)->terms;
    auto commonPieceFacts = [&](ExpressionPieceKind kind) {
      ExpressionPieceFacts piece;
      piece.pieceKind = kind;
      piece.combinatorialTags = facts.combinatorialTags;
      piece.coefficientRingIsQQ = coefficientRing == globalQQ;
      return piece;
    };

    if (pieceKind == ExpressionPieceKind::WholeExpression)
      {
        ExpressionPieceFacts piece =
            commonPieceFacts(pieceKind);
        piece.termPositions.resize(terms.size());
        std::iota(
            piece.termPositions.begin(),
            piece.termPositions.end(),
            size_t{0});
        piece.weight = facts.homogeneousWeight;
        piece.termCount = facts.termCount;
        piece.singleBasisElement = facts.singleBasisElement();
        piece.allPowerSumTermsSingleCycles =
            facts.allPowerSumTermsSingleCycles;
        piece.mostlyShortCyclePowerSumTermCount =
            facts.mostlyShortCyclePowerSumTermCount;
        piece.index = facts.singleBasisElementIndex;
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
        piece.expressionHasMultipleWeights = multipleWeights;
        return {std::move(piece)};
      }

    if (pieceKind == ExpressionPieceKind::IndividualTerms)
      {
        std::vector<ExpressionPieceFacts> result;
        result.reserve(terms.size());
        for (size_t position = 0; position < terms.size(); ++position)
          {
            const auto& term = terms[position];
            ExpressionPieceFacts piece =
                commonPieceFacts(pieceKind);
            piece.termPositions = {position};
            piece.weight = monomialWeight(term.monomial);
            piece.termCount = 1;
            if (!term.monomial.data.empty() &&
                atomLengthAt(term.monomial, 0) ==
                    term.monomial.data.size() &&
                !atomIsSkewAt(term.monomial, 0))
              {
                Partition index =
                    basisElementIndex(term.monomial, 0);
                piece.index = index;
                if (basisKindForId(
                        atomBasisIdAt(term.monomial, 0)) ==
                    BasisKind::PowerSum)
                  piece.powerSumIndexHasMostlyShortCycles =
                      powerSumIndexHasMostlyShortCycles(
                          index, partitionWeight(index));
              }
            result.push_back(std::move(piece));
          }
        return result;
      }

    std::map<int, std::vector<size_t>> positionsByWeight;
    for (size_t position = 0; position < terms.size(); ++position)
      positionsByWeight[monomialWeight(
          terms[position].monomial)].push_back(position);
    std::vector<ExpressionPieceFacts> result;
    result.reserve(positionsByWeight.size());
    for (auto& item : positionsByWeight)
      {
        ExpressionPieceFacts piece =
            commonPieceFacts(pieceKind);
        piece.termPositions = std::move(item.second);
        piece.weight = item.first;
        piece.termCount = piece.termPositions.size();
        piece.possibleTermCount =
            partitionCountForConversionSelection(item.first);
        size_t mostlyShortCycleTerms = 0;
        std::vector<int> commonParts;
        bool firstPowerSumIndex = true;
        bool allTermsArePowerSums = true;
        for (size_t position : piece.termPositions)
          {
            Partition index;
            if (!powerSumIndexFromMonomial(
                    terms[position].monomial, index))
              {
                allTermsArePowerSums = false;
                break;
              }
            if (index.empty()) continue;
            if (powerSumIndexHasMostlyShortCycles(
                    index, partitionWeight(index)))
              ++mostlyShortCycleTerms;
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
            piece.mostlyShortCyclePowerSumTermCount =
                mostlyShortCycleTerms;
            piece.commonPowerSumParts = std::move(commonParts);
          }
        result.push_back(std::move(piece));
      }
    return result;
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
SymmetricEngineRing::basisConversionPlansBySourceAndTarget()
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

void SymmetricEngineRing::validateBasisConversionPlanDatabase()
{
    static const bool validated = [] {
      const auto& database = basisConversionPlanDatabase();
      const auto& plansById = basisConversionPlansById();
      const auto& plansBySourceAndTarget =
          basisConversionPlansBySourceAndTarget();

      // The concrete plans together with the one generic u -> p -> v plan
      // cover the complete directed graph of built-in bases. Identity arrows
      // remain a generic workflow bypass.
      const std::vector<BasisKind> builtInConversionBases{
          BasisKind::PowerSum,
          BasisKind::Complete,
          BasisKind::Elementary,
          BasisKind::Monomial,
          BasisKind::Forgotten,
          BasisKind::Schur,
          BasisKind::SchurOmega,
          BasisKind::HallLittlewoodQGenerator,
          BasisKind::HallLittlewoodBGenerator,
          BasisKind::HallLittlewoodQ,
          BasisKind::HallLittlewoodB,
          BasisKind::HallLittlewoodP,
          BasisKind::HallLittlewoodPOmega};
      for (BasisKind source : builtInConversionBases)
        for (BasisKind target : builtInConversionBases)
          if (source != target)
            {
              const bool hasSpecificPlan =
                  plansBySourceAndTarget.find({source, target}) !=
                  plansBySourceAndTarget.end();
              const bool hasGenericPlan =
                  source != BasisKind::PowerSum &&
                  target != BasisKind::PowerSum;
              if (!hasSpecificPlan && !hasGenericPlan)
                throw exc::engine_error(
                    "basis-conversion plans do not cover every "
                    "ordered pair of built-in bases");
            }

      if (plansById.find(
              genericPowerSumFallbackPlanIdentifier().value) !=
          plansById.end())
        throw exc::engine_error(
            "the generic power-sum composition identifier must not "
            "belong to an endpoint-specific plan");
      std::map<BasisKind, const PowerSumFallbackPlanPair *>
          powerSumFallbacks;
      for (const auto& fallback :
           powerSumFallbackPlanDatabase())
        {
          if (fallback.basisKind == BasisKind::Custom ||
              fallback.basisKind == BasisKind::PowerSum ||
              !powerSumFallbacks.emplace(
                   fallback.basisKind, &fallback).second)
            throw exc::engine_error(
                "power-sum fallback bases must be distinct built-in "
                "non-power-sum bases");
          auto toPowerSums =
              plansById.find(fallback.toPowerSums.value);
          auto fromPowerSums =
              plansById.find(fallback.fromPowerSums.value);
          if (toPowerSums == plansById.end() ||
              fromPowerSums == plansById.end() ||
              toPowerSums->second->sourceBasisKind !=
                  fallback.basisKind ||
              toPowerSums->second->targetBasisKind !=
                  BasisKind::PowerSum ||
              fromPowerSums->second->sourceBasisKind !=
                  BasisKind::PowerSum ||
              fromPowerSums->second->targetBasisKind !=
                  fallback.basisKind ||
              toPowerSums->second->applicability.kind !=
                  ExpressionConditionKind::Always ||
              fromPowerSums->second->applicability.kind !=
                  ExpressionConditionKind::Always)
            throw exc::engine_error(
                "a power-sum fallback must name unconditional plans "
                "with the stated endpoints");
        }
      if (powerSumFallbacks.size() + 1 !=
          builtInConversionBases.size())
        throw exc::engine_error(
            "the generic power-sum composition does not cover every "
            "non-power-sum built-in basis");

      for (const auto& plan : database)
        {
          if (plan.sourceBasisKind == BasisKind::Custom ||
              plan.targetBasisKind == BasisKind::Custom ||
              plan.sourceBasisKind == plan.targetBasisKind)
            throw exc::engine_error(
                "a basis-conversion plan must have two distinct "
                "built-in endpoints");
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
                      plan.pieceKind,
                      isOtherwise);
                }
              catch (const std::exception& exception)
                {
                  throw exc::engine_error(
                      "invalid case condition for plan " +
                      plan.id.value + ": " + exception.what());
                }
              if (const auto *kernelCase =
                      std::get_if<KernelPlan>(
                          &item.formula))
                {
                  if (kernelCase->name.empty() ||
                      kernelCase->kernel == nullptr ||
                      !expressionConditionImplies(
                          kernelCase->outputGuarantee,
                          plan.outputGuarantee))
                    throw exc::engine_error(
                        "basis-conversion plan " +
                        plan.id.value +
                        " has an invalid kernel formula");
                  try
                    {
                      validateExpressionCondition(
                          kernelCase->outputGuarantee,
                          ExpressionPieceKind::WholeExpression);
                    }
                  catch (const std::exception& exception)
                    {
                      throw exc::engine_error(
                          "invalid kernel output guarantee for plan " +
                          plan.id.value + ": " + exception.what());
                    }
                }
              else
                {
                  const auto& composition =
                      std::get<
                          CompositionPlan>(
                          item.formula);
                  if (composition.plans.empty())
                    throw exc::engine_error(
                        "basis-conversion plan " + plan.id.value +
                        " has an empty composition");
                  composition.planDefinitions.clear();
                }
            }
        }

      std::map<std::string, std::vector<std::string>>
          dependencies;
      for (const auto& plan : database)
        for (const auto& item : plan.cases)
          {
            if (std::holds_alternative<
                    KernelPlan>(
                    item.formula))
              continue;
            auto& composition =
                std::get<
                    CompositionPlan>(
                    item.formula);
            BasisKind current = plan.sourceBasisKind;
            ExpressionCondition guaranteedCondition =
                plan.pieceKind ==
                        ExpressionPieceKind::WholeExpression
                    ? plan.applicability
                    : always();
            if (plan.pieceKind ==
                    ExpressionPieceKind::WholeExpression &&
                item.condition.kind !=
                    ExpressionConditionKind::Otherwise)
              guaranteedCondition =
                  guaranteedCondition && item.condition;
            for (const auto& componentPlanId : composition.plans)
              {
                auto found =
                    plansById.find(componentPlanId.value);
                if (found == plansById.end())
                  throw exc::engine_error(
                      "basis-conversion plan " +
                      plan.id.value +
                      " references unknown component plan " +
                      componentPlanId.value);
                const auto& componentPlan = *found->second;
                composition.planDefinitions.push_back(
                    found->second);
                if (componentPlan.sourceBasisKind != current)
                  throw exc::engine_error(
                      "basis-conversion plan " +
                      plan.id.value +
                      " has incompatible composition endpoints");
                if (!expressionConditionImplies(
                        guaranteedCondition,
                        componentPlan.applicability))
                  throw exc::engine_error(
                      "basis-conversion composition " +
                      plan.id.value +
                      "references a plan with an unproved "
                      "applicability requirement");
                current = componentPlan.targetBasisKind;
                guaranteedCondition =
                    componentPlan.outputGuarantee;
                dependencies[plan.id.value].push_back(
                    componentPlan.id.value);
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

const std::vector<SymmetricEngineRing::RingBasisConversionPlan>&
SymmetricEngineRing::basisConversionPlansFor(
    int sourceBasisId,
    int targetBasisId) const
{
    static const std::vector<RingBasisConversionPlan> empty;
    requireBasis(sourceBasisId);
    requireBasis(targetBasisId);
    if (error()) return empty;
    const std::pair<int, int> endpoint{
        sourceBasisId, targetBasisId};
    auto cached =
        ringBasisConversionPlanCache.find(endpoint);
    if (cached != ringBasisConversionPlanCache.end())
      return cached->second;
    validateBasisConversionPlanDatabase();
    const BasisKind source = basisKindForId(sourceBasisId);
    const BasisKind target = basisKindForId(targetBasisId);
    std::vector<RingBasisConversionPlan> result;
    const auto& plansBySourceAndTarget =
        basisConversionPlansBySourceAndTarget();
    auto found = plansBySourceAndTarget.find({source, target});
    if (found != plansBySourceAndTarget.end())
      {
        result.reserve(found->second.size());
        for (const auto *plan : found->second)
          result.push_back(
              {plan, sourceBasisId, targetBasisId});
      }
    auto inserted =
        ringBasisConversionPlanCache.emplace(
            endpoint, std::move(result));
    return inserted.first->second;
  }

SymmetricEngineRing::RingBasisConversionPlan
SymmetricEngineRing::basisConversionPlanForRing(
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
    return basisConversionPlanForRing(definition);
  }

SymmetricEngineRing::RingBasisConversionPlan
SymmetricEngineRing::basisConversionPlanForRing(
    const BasisConversionPlanDefinition *definition) const
{
    if (definition == nullptr) return {};
    auto cached =
        basisConversionEndpointsForRingCache.find(
            definition->id.value);
    if (cached !=
        basisConversionEndpointsForRingCache.end())
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
    basisConversionEndpointsForRingCache.emplace(
        definition->id.value,
        std::make_pair(sourceBasisId, targetBasisId));
    return {definition, sourceBasisId, targetBasisId};
  }

// ============================================================================
// Generic Declarative Formula And Plan Executor
// ============================================================================
// The executor interprets the already-selected plan. Kernel formulas call one
// kernel; compositions identify their fixed component plans recursively.
// Neither path invokes performance policy or substitutes another plan.

ring_elem SymmetricEngineRing::executeBasisConversionPlanFormula(
    const BasisConversionPlanFormula& formula,
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
    if (const auto *kernelCase =
            std::get_if<KernelPlan>(
                &formula))
      {
        if (kernelCase->kernel == nullptr)
          {
            ERROR("a basis-conversion kernel formula has no callable");
            return zero();
          }
        const auto& sourceDescriptor =
            requireBasis(sourceBasisId);
        const auto& targetDescriptor =
            requireBasis(targetBasisId);
        if (error()) return zero();
        BasisConversionInput input{
            expression,
            {sourceDescriptor.kind,
             sourceBasisId,
             &sourceDescriptor.displaySymbol,
             sourceDescriptor.displayOrder},
            {targetDescriptor.kind,
             targetBasisId,
             &targetDescriptor.displaySymbol,
             targetDescriptor.displayOrder},
            inputFacts.homogeneousWeight};
        ring_elem result = (this->*kernelCase->kernel)(input);
        if (error()) return zero();

        const bool needFacts =
            resultFacts != nullptr ||
            kernelCase->outputGuarantee.kind !=
                ExpressionConditionKind::Always;
        ExpressionFacts kernelFacts;
        if (needFacts)
          {
            kernelFacts = inferCanonicalExpansionFacts(
                result,
                targetBasisId,
                inputFacts.homogeneousWeight);
            kernelFacts.combinatorialTags =
                combinatorialTags;
            if (!kernelFacts.canonicalExpansionInBasis(
                    targetBasisId))
              {
                ERROR("basis-conversion kernel ",
                      kernelCase->name.c_str(),
                      " violated its canonical target contract");
                return zero();
              }
          }
        if (kernelCase->outputGuarantee.kind !=
            ExpressionConditionKind::Always)
          {
            const auto pieceFacts =
                inspectExpressionPieces(
                    result,
                    ExpressionPieceKind::WholeExpression,
                    kernelFacts);
            if (pieceFacts.size() != 1 ||
                !expressionConditionHolds(
                    kernelCase->outputGuarantee,
                    pieceFacts.front()))
              {
                ERROR("basis-conversion kernel ",
                      kernelCase->name.c_str(),
                      " violated its declared output guarantee");
                return zero();
              }
          }
        mutablePolyValue(result)->combinatorialTags =
            combinatorialTags;
        if (resultFacts != nullptr)
          *resultFacts = std::move(kernelFacts);
        return result;
      }

    ring_elem current = expression;
    ExpressionFacts currentFacts = inputFacts;
    int currentBasisId = sourceBasisId;
    auto executeComponentPlan =
        [&](const BasisConversionPlanDefinition *componentDefinition)
            -> bool {
          if (componentDefinition == nullptr)
            {
              ERROR("a basis-conversion composition reached execution "
                    "before a component plan was identified");
              return false;
            }
          RingBasisConversionPlan componentPlan =
              basisConversionPlanForRing(
                  componentDefinition);
          if (error()) return false;
          if (componentPlan.sourceBasisId != currentBasisId)
            {
              ERROR("a basis-conversion composition has "
                    "incompatible runtime endpoints");
              return false;
            }
          current = executeBasisConversionPlan(
              componentPlan,
              current,
              combinatorialTags,
              currentFacts,
              nullptr);
          if (error()) return false;
          currentBasisId = componentPlan.targetBasisId;
          // Conditions in the next fixed component plan apply to the actual
          // realized intermediate, never to a predicted support profile.
          currentFacts = inferCanonicalExpansionFacts(
              current,
              currentBasisId,
              inputFacts.homogeneousWeight);
          currentFacts.combinatorialTags =
              combinatorialTags;
          return true;
        };

    const auto& composition =
        std::get<
            CompositionPlan>(
            formula);
    if (composition.planDefinitions.size() !=
        composition.plans.size())
      {
        ERROR("a basis-conversion composition reached execution "
              "before its named plans were identified");
        return zero();
      }
    for (const auto *componentDefinition :
         composition.planDefinitions)
      if (!executeComponentPlan(componentDefinition))
        return zero();
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
    const RingBasisConversionPlan& plan,
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
        plan.definition->pieceKind ==
            ExpressionPieceKind::WholeExpression &&
        plan.definition->cases.size() == 1 &&
        plan.definition->cases.front().condition.kind ==
            ExpressionConditionKind::Otherwise;
    ExpressionFacts executionFacts = inputFacts;
    executionFacts.combinatorialTags = combinatorialTags;
    const bool needsSelectionProfiles =
        plan.definition->applicability.kind !=
            ExpressionConditionKind::Always ||
        (plan.definition->pieceKind ==
             ExpressionPieceKind::WholeExpression &&
         !wholeExpressionOtherwise);
    // Term and homogeneous-component pieceFacts derive their own profiles from
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
    if (basisConversionTraceEnabled())
      std::fprintf(
          stderr,
          "SymmetricRings conversion-plan: source=%s target=%s "
          "plan=%s\n",
          basisKeyForId(plan.sourceBasisId).c_str(),
          basisKeyForId(plan.targetBasisId).c_str(),
          plan.definition->id.value.c_str());

    if (wholeExpressionOtherwise)
      {
        // Most available plans are one whole-expression formula. Their
        // declarative meaning is exactly a direct call on the original input;
        // constructing a one-piece partition and then copying and recollecting
        // the complete expression would add work without changing semantics.
        const bool needDirectFacts =
            resultFacts != nullptr ||
            plan.definition->outputGuarantee.kind !=
                ExpressionConditionKind::Always;
        ExpressionFacts directFacts;
        ring_elem directResult =
            executeBasisConversionPlanFormula(
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
                const auto outputPieceFacts =
                    inspectExpressionPieces(
                        directResult,
                        ExpressionPieceKind::WholeExpression,
                        directFacts);
                if (outputPieceFacts.size() != 1 ||
                    !expressionConditionHolds(
                        plan.definition->outputGuarantee,
                        outputPieceFacts.front()))
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

    std::vector<ExpressionPieceFacts> pieces;
    PlanCaseAssignments partition;
    try
      {
        pieces = inspectExpressionPieces(
            expression,
            plan.definition->pieceKind,
            executionFacts);
        std::vector<ExpressionCondition> conditions;
        conditions.reserve(
            plan.definition->cases.size());
        for (const auto& item : plan.definition->cases)
          conditions.push_back(item.condition);
        partition = assignExpressionPiecesToCases(
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
            partition.termPositionsForCase[caseIndex];
        if (termPositions.empty()) continue;

        // A case formula receives the complete subexpression assigned to that
        // case, not one invocation per term or component. This preserves the
        // declarative semantics and lets fixed component plans inspect the
        // actual aggregate support produced for the case.
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
        if (plan.definition->pieceKind ==
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
            executeBasisConversionPlanFormula(
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
            const auto outputPieceFacts =
                inspectExpressionPieces(
                    result,
                    ExpressionPieceKind::WholeExpression,
                    finalFacts);
            if (outputPieceFacts.size() != 1 ||
                !expressionConditionHolds(
                    plan.definition->outputGuarantee,
                    outputPieceFacts.front()))
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
    const bool trace = basisConversionTraceEnabled();

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
            RingBasisConversionPlan selection =
                selectBasisConversionPlan(
                    f,
                    sourceBasisId,
                    targetBasisId,
                    *knownFacts,
                    combinatorialTags,
                    forcedBasisConversionPlanIdentifier(),
                    &selectedInputFacts);
            if (error()) return zero();
            result = executeBasisConversionPlan(
                selection,
                f,
                combinatorialTags,
                selectedInputFacts,
                &resultFacts);
            if (error()) return zero();
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
      std::optional<RingBasisConversionPlan> selection;
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
      std::optional<RingBasisConversionPlan> selection;
      if (sourceBasisId != targetBasisId)
        {
          selection = selectBasisConversionPlan(
              expression,
              sourceBasisId,
              targetBasisId,
              facts,
              combinatorialTags,
              forcedBasisConversionPlanIdentifier(),
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
    // Each selected definition fixes every case, kernel, and component plan.
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
// Multiplication Plans
// ============================================================================
// Plans are listed from most structured to broadest. The final power-sum plan
// is the independent correctness path for every built-in operand pair.
//
// Multiplication occupies the remainder of this file before the public
// toBasis definition because product resolution in toBasis calls these
// helpers. This dependency order keeps the owning workflow acyclic while the
// block headers keep conversion and multiplication topics distinct.

bool SymmetricEngineRing::hasSchurFactorFormula(BasisKind kind)
{
    return kind == BasisKind::Schur ||
           kind == BasisKind::Complete ||
           kind == BasisKind::Elementary ||
           kind == BasisKind::PowerSum;
  }

bool SymmetricEngineRing::hasMonomialLikeFactorFormula(
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
          return "Schur-product:Littlewood-Richardson";
        case MultiplicationPlanId::SchurHorizontalPieri:
          return "Schur-product:horizontal-Pieri";
        case MultiplicationPlanId::SchurVerticalPieri:
          return "Schur-product:vertical-Pieri";
        case MultiplicationPlanId::SchurMurnaghanNakayama:
          return "Schur-product:Murnaghan-Nakayama";
        case MultiplicationPlanId::SchurMixedPieriAndMurnaghanNakayama:
          return "Schur-product:mixed-Pieri-and-Murnaghan-Nakayama";
        case MultiplicationPlanId::MonomialLikeExponentSplittings:
          return "monomial-like-product:exponent-splittings";
        case MultiplicationPlanId::ViaHallLittlewoodGenerators:
          return "Hall-Littlewood-product:via-generator-basis";
        case MultiplicationPlanId::InTargetBasis:
          return "product:in-target-basis";
        case MultiplicationPlanId::ViaPowerSums:
          return "product:via-power-sums";
      }
    return "unavailable";
  }

const std::vector<SymmetricEngineRing::MultiplicationPlan>&
SymmetricEngineRing::multiplicationPlansFor(
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId) const
{
    const int leftBasisId = leftFacts.expandedBasis.value_or(-1);
    const int rightBasisId = rightFacts.expandedBasis.value_or(-1);
    const std::tuple<int, int, int> key{
        leftBasisId, rightBasisId, targetBasisId};
    auto found = multiplicationPlansCache.find(key);
    if (found != multiplicationPlansCache.end())
      return found->second;
    auto inserted = multiplicationPlansCache.emplace(
        key,
        buildMultiplicationPlansFor(
            leftFacts, rightFacts, targetBasisId));
    return inserted.first->second;
  }

std::vector<SymmetricEngineRing::MultiplicationPlan>
SymmetricEngineRing::buildMultiplicationPlansFor(
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
        int specializedFactorTypes =
            static_cast<int>(*leftKind == BasisKind::Complete ||
                             *rightKind == BasisKind::Complete) +
            static_cast<int>(*leftKind == BasisKind::Elementary ||
                             *rightKind == BasisKind::Elementary) +
            static_cast<int>(*leftKind == BasisKind::PowerSum ||
                             *rightKind == BasisKind::PowerSum);
        if (specializedFactorTypes > 1)
          planId =
              MultiplicationPlanId::SchurMixedPieriAndMurnaghanNakayama;
        else if (*leftKind == BasisKind::Complete ||
                 *rightKind == BasisKind::Complete)
          planId = MultiplicationPlanId::SchurHorizontalPieri;
        else if (*leftKind == BasisKind::Elementary ||
                 *rightKind == BasisKind::Elementary)
          planId = MultiplicationPlanId::SchurVerticalPieri;
        else if (*leftKind == BasisKind::PowerSum ||
                 *rightKind == BasisKind::PowerSum)
          planId = MultiplicationPlanId::SchurMurnaghanNakayama;
        plans.push_back({
            planId,
            MultiplicationKernel::SchurFactorFormulas,
            -1,
            -1,
            targetBasisId,
            true});
      }

    if ((targetKind == BasisKind::Monomial ||
         targetKind == BasisKind::Forgotten))
      plans.push_back({
          MultiplicationPlanId::MonomialLikeExponentSplittings,
          MultiplicationKernel::MonomialExponentSplittings,
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
            MultiplicationPlanId::ViaHallLittlewoodGenerators,
            MultiplicationKernel::ProductInSingleBasis,
            generatorId,
            generatorId,
            generatorId,
            true});
      }

    if (isMultiplicativeBasis(targetBasisId))
      plans.push_back({
          MultiplicationPlanId::InTargetBasis,
          MultiplicationKernel::ProductInSingleBasis,
          targetBasisId,
          targetBasisId,
          targetBasisId,
          true});
    plans.push_back({
        MultiplicationPlanId::ViaPowerSums,
        MultiplicationKernel::ProductInSingleBasis,
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
        plan.id == MultiplicationPlanId::SchurMurnaghanNakayama ||
        plan.id == MultiplicationPlanId::SchurMixedPieriAndMurnaghanNakayama;
    if ((schurPlan &&
         plan.kernel != MultiplicationKernel::SchurFactorFormulas) ||
        (plan.id ==
             MultiplicationPlanId::MonomialLikeExponentSplittings &&
         plan.kernel != MultiplicationKernel::MonomialExponentSplittings) ||
        ((plan.id == MultiplicationPlanId::ViaHallLittlewoodGenerators ||
          plan.id == MultiplicationPlanId::InTargetBasis ||
          plan.id == MultiplicationPlanId::ViaPowerSums) &&
         plan.kernel != MultiplicationKernel::ProductInSingleBasis) ||
        plan.id == MultiplicationPlanId::Unavailable)
      return false;
    if (plan.kernel == MultiplicationKernel::SchurFactorFormulas)
      {
        if (!leftFacts.expandedBasis || !rightFacts.expandedBasis)
          return false;
        BasisKind leftKind = basisKindForId(*leftFacts.expandedBasis);
        BasisKind rightKind = basisKindForId(*rightFacts.expandedBasis);
        BasisKind targetKind = basisKindForId(targetBasisId);
        return plan.productBasisId == targetBasisId &&
               targetKind == BasisKind::Schur &&
               hasSchurFactorFormula(leftKind) &&
               hasSchurFactorFormula(rightKind) &&
               leftFacts.skewFree && rightFacts.skewFree;
      }
    if (plan.kernel == MultiplicationKernel::MonomialExponentSplittings)
      {
        if (!leftFacts.expandedBasis || !rightFacts.expandedBasis)
          return false;
        BasisKind leftKind = basisKindForId(*leftFacts.expandedBasis);
        BasisKind rightKind = basisKindForId(*rightFacts.expandedBasis);
        BasisKind targetKind = basisKindForId(targetBasisId);
        return plan.productBasisId == targetBasisId &&
               (targetKind == BasisKind::Monomial ||
                targetKind == BasisKind::Forgotten) &&
               hasMonomialLikeFactorFormula(
                   leftKind, targetKind) &&
               hasMonomialLikeFactorFormula(
                   rightKind, targetKind) &&
               leftFacts.skewFree && rightFacts.skewFree;
      }
    if (plan.kernel != MultiplicationKernel::ProductInSingleBasis ||
        !plan.productIsCanonical ||
        plan.leftOperandBasisId <= 0 ||
        plan.rightOperandBasisId <= 0)
      return false;
    switch (plan.id)
      {
        case MultiplicationPlanId::ViaHallLittlewoodGenerators:
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
        case MultiplicationPlanId::InTargetBasis:
          return plan.leftOperandBasisId == targetBasisId &&
                 plan.rightOperandBasisId == targetBasisId &&
                 plan.productBasisId == targetBasisId &&
                 isMultiplicativeBasis(targetBasisId);
        case MultiplicationPlanId::ViaPowerSums:
          {
            const int powerSumBasisId =
                requiredBasisIdForKind(BasisKind::PowerSum);
            return !error() &&
                   plan.leftOperandBasisId == powerSumBasisId &&
                   plan.rightOperandBasisId == powerSumBasisId &&
                   plan.productBasisId == powerSumBasisId;
          }
        case MultiplicationPlanId::Unavailable:
        case MultiplicationPlanId::SchurLittlewoodRichardson:
        case MultiplicationPlanId::SchurHorizontalPieri:
        case MultiplicationPlanId::SchurVerticalPieri:
        case MultiplicationPlanId::SchurMurnaghanNakayama:
        case MultiplicationPlanId::SchurMixedPieriAndMurnaghanNakayama:
        case MultiplicationPlanId::MonomialLikeExponentSplittings:
          return false;
      }
    return false;
  }

// ============================================================================
// Multiplication-Plan Picker
// ============================================================================
// Operand conversions are selected here as complete plans together with the
// product kernel. The final conversion is intentionally selected by the
// owning workflow only after the realized product has been normalized and
// inspected.

SymmetricEngineRing::MultiplicationPlanSelection
SymmetricEngineRing::selectMultiplicationPlan(
    ring_elem f,
    ring_elem g,
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId) const
{
    const auto& plans =
        multiplicationPlansFor(
            leftFacts, rightFacts, targetBasisId);
    auto withOperandConversions =
        [&](const MultiplicationPlan& plan) {
          MultiplicationPlanSelection selection{
              plan, std::nullopt, std::nullopt};
          if (plan.leftOperandBasisId > 0 &&
              *leftFacts.expandedBasis != plan.leftOperandBasisId)
            selection.leftConversion =
                selectBasisConversionPlan(
                    f,
                    *leftFacts.expandedBasis,
                    plan.leftOperandBasisId,
                    leftFacts,
                    0);
          if (error()) return selection;
          if (plan.rightOperandBasisId > 0 &&
              *rightFacts.expandedBasis != plan.rightOperandBasisId)
            selection.rightConversion =
                selectBasisConversionPlan(
                    g,
                    *rightFacts.expandedBasis,
                    plan.rightOperandBasisId,
                    rightFacts,
                    0);
          return selection;
        };
    const char *forced =
        basisConversionContextActive()
            ? nullptr
            : std::getenv(
                  "M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN");
    if (forced != nullptr)
      {
        for (const auto& plan : plans)
          if (multiplicationPlanIdentifier(plan.id) == forced &&
              multiplicationPlanApplicable(
                  plan, leftFacts, rightFacts, targetBasisId))
            return withOperandConversions(plan);
        ERROR("the requested multiplication plan is unknown or inapplicable: ",
              forced);
        int powerSumId = requiredBasisIdForKind(BasisKind::PowerSum);
        return {{MultiplicationPlanId::Unavailable,
                 MultiplicationKernel::ProductInSingleBasis,
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
        return withOperandConversions(plan);
    ERROR("no applicable multiplication plan is available");
    int powerSumId = requiredBasisIdForKind(BasisKind::PowerSum);
    return {{MultiplicationPlanId::Unavailable,
             MultiplicationKernel::ProductInSingleBasis,
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
    ring_elem product;
    if (plan.kernel == MultiplicationKernel::SchurFactorFormulas)
      {
        bool completed = tryProductToSchurViaCompatibleFactors(
            f,
            g,
            targetBasisId,
            displayForBasis(targetBasisId),
            basisOrderForId(targetBasisId),
            product);
        if (error()) return zero();
        if (!completed)
          {
            ERROR("a selected multiplication plan violated its "
                  "applicability contract");
            return zero();
          }
      }
    else if (plan.kernel ==
             MultiplicationKernel::MonomialExponentSplittings)
      {
        bool completed = tryProductToMonomialLikeTarget(
            f,
            g,
            targetBasisId,
            displayForBasis(targetBasisId),
            basisOrderForId(targetBasisId),
            isMultiplicativeBasis(targetBasisId),
            product);
        if (error()) return zero();
        if (!completed)
          {
            ERROR("a selected multiplication kernel violated its "
                  "applicability contract");
            return zero();
          }
      }
    else
      product = mult(f, g);
    return error() ? zero() : product;
  }

// ============================================================================
// Binary Multiplication Workflow
// ============================================================================
// Stages are deliberately visible and linear: convert operands, execute the
// selected kernel, normalize its product, select and execute one complete
// final conversion plan, and validate the target contract.

ring_elem SymmetricEngineRing::runMultiplicationWorkflow(
    const MultiplicationPlanSelection& selection,
    ring_elem f,
    ring_elem g,
    const ExpressionFacts& leftFacts,
    const ExpressionFacts& rightFacts,
    int targetBasisId,
    bool attachResultFacts) const
{
    const bool trace = basisConversionTraceEnabled();
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
          basisKeyForId(plan.productBasisId).c_str(),
          plan.productIsCanonical ? "yes" : "no");
    ring_elem product = executeMultiplicationKernel(
        plan, left, right, targetBasisId);
    if (error()) return zero();

    ExpressionFacts productFacts;
    ring_elem canonicalProduct;
    if (plan.productIsCanonical)
      {
        canonicalProduct = product;
        // Exact product support is necessary only for a final conversion or
        // for public metadata. A transient target-basis pair product can rely
        // on the kernel's canonical-output contract.
        if (plan.productBasisId != targetBasisId ||
            attachResultFacts)
          {
            auto metadataFacts =
                expressionFactsFromMetadata(canonicalProduct);
            productFacts = metadataFacts
                ? *metadataFacts
                : inferCanonicalExpansionFacts(
                      canonicalProduct,
                      plan.productBasisId);
          }
      }
    else
      canonicalProduct =
          normalizeExpression(product, productFacts);
    if (error()) return zero();
    if ((plan.productBasisId != targetBasisId ||
         attachResultFacts) &&
        !productFacts.canonicalExpansionInBasis(
            plan.productBasisId))
      {
        ERROR("a multiplication kernel violated its declared output contract");
        return zero();
      }

    mutablePolyValue(canonicalProduct)->combinatorialTags = tags;
    if (plan.productBasisId == targetBasisId &&
        !attachResultFacts)
      return canonicalProduct;
    productFacts.combinatorialTags = tags;
    ring_elem result = canonicalProduct;
    ExpressionFacts resultFacts;
    if (plan.productBasisId != targetBasisId)
      {
        RingBasisConversionPlan finalConversion =
            selectBasisConversionPlan(
                canonicalProduct,
                plan.productBasisId,
                targetBasisId,
                productFacts,
                tags,
                std::nullopt,
                &productFacts);
        if (error()) return zero();
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-stage: "
              "stage=convert-product plan=%s\n",
              finalConversion.
                  definition->id.value.c_str());
        result = executeBasisConversionPlan(
            finalConversion,
            canonicalProduct,
            tags,
            productFacts,
            attachResultFacts ? &resultFacts : nullptr);
      }
    else if (attachResultFacts)
      resultFacts = std::move(productFacts);
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
    if (basisConversionTraceEnabled())
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
// uses target-basis multiplication or the pairwise plan prescribed by the
// design.

SymmetricEngineRing::ResolvedProductTerm
SymmetricEngineRing::multiplyTermToBasis(
    const SymmetricTerm& term,
    int targetBasisId) const
{
    const bool trace = basisConversionTraceEnabled();
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
          RingBasisConversionPlan selection =
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
              "branch=target-basis-product factors=%zu target=%s\n",
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
// coefficients and sends each nonscalar pair to the strict helper above.

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
          RingBasisConversionPlan conversion =
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
// Basis-Conversion Input Preparation
// ============================================================================
// This stage establishes the normalized, skew-free factor expansion required
// before products can be resolved. A canonical element or exact metadata proves
// that the same mathematical preparation has already been completed.

SymmetricEngineRing::PreparedBasisConversionInput
SymmetricEngineRing::prepareBasisConversionInput(
    ring_elem expression,
    int targetBasisId) const
{
    const bool trace = basisConversionTraceEnabled();
    CombinatorialTags combinatorialTags =
        polyValue(expression)->combinatorialTags;
    const auto *inputPoly = polyValue(expression);
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
          return {
              expression,
              std::move(strictFacts),
              {},
              combinatorialTags};
      }
    auto knownFacts =
        expressionFactsFromMetadata(expression);
    if (knownFacts)
      {
        if (trace)
          std::fprintf(stderr,
                       "SymmetricRings toBasis: "
                       "bypass=canonical-metadata target=%s\n",
                       displayForBasis(targetBasisId).c_str());
        return {
            expression,
            std::move(*knownFacts),
            {},
            combinatorialTags};
      }

    ExpressionFacts normalizedFacts;
    std::vector<size_t> normalizedFactorCounts;
    ring_elem normalized =
        normalizeExpression(
            expression,
            normalizedFacts,
            &normalizedFactorCounts);
    if (error())
      return {zero(), {}, {}, combinatorialTags};

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
    return {
        normalized,
        std::move(normalizedFacts),
        std::move(normalizedFactorCounts),
        combinatorialTags};
  }

// ============================================================================
// Product Resolution For Basis Conversion
// ============================================================================
// Conversion plans act linearly on canonical source-basis expansions. This
// stage replaces every multifactor term by its canonical target-basis product
// while retaining scalar and single-factor terms. Its output is therefore
// product-free and ready for source-basis decomposition.

SymmetricEngineRing::ProductFreeBasisConversionInput
SymmetricEngineRing::resolveProductsForBasisConversion(
    PreparedBasisConversionInput prepared,
    int targetBasisId) const
{
    if (prepared.facts.noProducts())
      return {
          prepared.expression,
          std::move(prepared.facts),
          false};

    const auto& normalizedTerms =
        polyValue(prepared.expression)->terms;
    if (prepared.factorsPerTerm.size() !=
        normalizedTerms.size())
      {
        ERROR("normalized factor-count profile is inconsistent");
        return {zero(), std::nullopt, false};
      }
    if (normalizedTerms.size() == 1 &&
        prepared.factorsPerTerm.front() > 1)
      {
        ResolvedProductTerm resolved =
            multiplyTermToBasis(
                normalizedTerms.front(),
                targetBasisId);
        if (error())
          return {zero(), std::nullopt, false};
        resolved.facts.combinatorialTags =
            prepared.combinatorialTags;
        // multiplyTermToBasis already attached the same exact structural
        // facts. Only the semantic tag belongs to the owning outer operation.
        mutablePolyValue(
            resolved.expression)->combinatorialTags =
                prepared.combinatorialTags;
        return {
            resolved.expression,
            std::move(resolved.facts),
            true};
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
            prepared.factorsPerTerm[termPosition];
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
        if (error())
          return {zero(), std::nullopt, false};
        for (const auto& convertedTerm :
             polyValue(resolved.expression)->terms)
          appendTermIfNonZero(
              preparedTerms,
              convertedTerm.coeff,
              convertedTerm.monomial);
        ++resolvedProductTerms;
      }

    if (basisConversionTraceEnabled())
      std::fprintf(
          stderr,
          "SymmetricRings conversion-stage: "
          "stage=resolve-products resolved=%zu passthrough=%zu\n",
          resolvedProductTerms,
          preparedTerms.size());

    ring_elem productFree =
        fromTermVector(preparedTerms, false);
    mutablePolyValue(productFree)->combinatorialTags =
        prepared.combinatorialTags;
    if (passthroughTermsAlreadyInTarget)
      {
        // Every multifactor term was resolved into the target above, and all
        // scalar/single-factor passthrough terms were already target-native.
        // The combined expression is therefore target-closed; do not regroup
        // it merely to execute an identity plan.
        ExpressionFacts productFreeFacts =
            inferCanonicalExpansionFacts(
                productFree, targetBasisId);
        if (!productFreeFacts.canonicalExpansionInBasis(
                targetBasisId))
          {
            ERROR("resolved product terms did not satisfy the target "
                  "basis contract");
            return {zero(), std::nullopt, false};
          }
        attachExpressionFacts(
            productFree,
            productFreeFacts,
            targetBasisId,
            prepared.combinatorialTags);
        return {
            productFree,
            std::move(productFreeFacts),
            true};
      }
    return {productFree, std::nullopt, false};
  }

// ============================================================================
// Unified Basis-Conversion Entry Workflow
// ============================================================================
// This is the sole owning engine workflow for built-in conversion. Its body
// follows the mathematical decomposition: prepare canonical factors, resolve
// products, decompose by source basis, apply one complete plan to each source
// expansion, and collect the target expansion.

ring_elem SymmetricEngineRing::toBasis(
    ring_elem f,
    int targetBasisId) const
{
    const BasisConversionDiagnosticContext ordinaryConversion;
    const ScopedBasisConversionDiagnostics diagnosticScope(
        ordinaryConversion,
        false);
    requireBasis(targetBasisId);
    if (error()) return zero();

    PreparedBasisConversionInput prepared =
        prepareBasisConversionInput(f, targetBasisId);
    if (error()) return zero();
    const CombinatorialTags combinatorialTags =
        prepared.combinatorialTags;

    ProductFreeBasisConversionInput productFree =
        resolveProductsForBasisConversion(
            std::move(prepared), targetBasisId);
    if (error()) return zero();
    if (productFree.productResolutionCompletedInTarget)
      return productFree.expression;

    // Linearity now applies: the canonical helper decomposes the expression
    // by source basis, selects one complete source-to-target plan for each
    // summand, executes those plans, and collects their target expansions.
    const ExpressionFacts *exactFacts =
        productFree.exactFacts
            ? &*productFree.exactFacts
            : nullptr;
    return convertCanonicalExpressionToBasis(
        productFree.expression,
        targetBasisId,
        combinatorialTags,
        exactFacts);
  }

ring_elem SymmetricEngineRing::toBasisBench(
    ring_elem f,
    int targetBasisId,
    const std::optional<std::string>& forcedPlanIdentifier,
    bool traceConversion) const
{
    const BasisConversionDiagnosticContext benchmarkConversion{
        forcedPlanIdentifier,
        traceConversion};
    const ScopedBasisConversionDiagnostics diagnosticScope(
        benchmarkConversion,
        true);
    return toBasis(f, targetBasisId);
  }

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
