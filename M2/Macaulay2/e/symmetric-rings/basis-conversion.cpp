// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "basic-rings/aring-glue.hpp"
#include "error.h"
#include "symmetric-rings/basis-conversion-policy.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <queue>
#include <sstream>
#include <utility>

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

void SymmetricEngineRing::enrichConversionSelectionFacts(
    ring_elem f,
    int targetBasisId,
    ExpressionFacts& facts) const
{
    if (!facts.expandedBasis ||
        basisKindForId(*facts.expandedBasis) != BasisKind::PowerSum)
      return;
    BasisKind targetKind = basisKindForId(targetBasisId);
    const bool schurTarget =
        targetKind == BasisKind::Schur ||
        targetKind == BasisKind::SchurOmega;
    const bool hallLittlewoodTarget =
        targetKind == BasisKind::HallLittlewoodQ ||
        targetKind == BasisKind::HallLittlewoodB ||
        targetKind == BasisKind::HallLittlewoodP ||
        targetKind == BasisKind::HallLittlewoodPOmega;
    if (!schurTarget && !hallLittlewoodTarget) return;

    if (schurTarget && !facts.homogeneousWeight &&
        facts.weightTermPositions.empty())
      {
        size_t position = 0;
        for (const auto& term : polyValue(f)->terms)
          facts.weightTermPositions[
              monomialWeight(term.monomial)].push_back(position++);
      }
    if (schurTarget && facts.homogeneousWeight &&
        !facts.possibleTermCount)
      {
        const size_t possible =
            partitionCountForConversionSelection(*facts.homogeneousWeight);
        facts.possibleTermCount = possible;
        if (possible != 0)
          {
            facts.density =
                static_cast<double>(facts.termCount) /
                static_cast<double>(possible);
          }
      }

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
        !metadata.noProducts || !*metadata.noProducts ||
        !metadata.termCount || !metadata.singleTerm ||
        !metadata.singleBasisElement)
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
    // A zero or scalar expansion has no distinguished basis. Every
    // nonscalar canonical expansion must name its unique source basis.
    if (facts.singleFactorTermCount != 0 && !metadata.expandedBasis)
      return std::nullopt;
    facts.pureBasis = metadata.pureBasis
        ? metadata.pureBasis : metadata.expandedBasis;
    facts.expandedBasis = metadata.expandedBasis;
    if (metadata.factorBases)
      facts.factorBases = *metadata.factorBases;
    else if (metadata.expandedBasis)
      facts.factorBases = {*metadata.expandedBasis};
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
    facts.density = metadata.density;
    facts.singleBasisElementId =
        metadata.singleBasisElementId;
    facts.singleBasisElementIndex =
        metadata.singleBasisElementIndex;
    facts.singleBasisElementCoefficientOne =
        metadata.singleBasisElementCoefficientOne;
    return facts;
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
    (void) targetBasisId;
    SymmetricConversionMetadata metadata;
    metadata.pureBasis = facts.pureBasis;
    metadata.expandedBasis = facts.expandedBasis;
    metadata.homogeneousWeight = facts.homogeneousWeight;
    metadata.termCount = facts.termCount;
    metadata.maximumPartitionLength = facts.maximumPartitionLength;
    metadata.density = facts.density;
    metadata.factorBases = facts.factorBases;
    metadata.singleBasisElement = facts.singleBasisElement();
    metadata.singleTerm = facts.singleTerm();
    metadata.noProducts = facts.noProducts();
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
    if (kind == BasisKind::HallLittlewoodQ ||
        kind == BasisKind::HallLittlewoodB ||
        kind == BasisKind::HallLittlewoodP ||
        kind == BasisKind::HallLittlewoodPOmega)
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

SymmetricEngineRing::ExpressionFacts
SymmetricEngineRing::conservativeExpressionFactsForBasis(
    const ExpressionFacts& sourceFacts,
    int basisId) const
{
    // Plan selection for an unexecuted intermediate edge may rely only on
    // invariants guaranteed by every preceding plan. Shape- and support-based
    // facts intentionally remain unknown until that intermediate exists.
    ExpressionFacts facts;
    facts.normalized = true;
    facts.skewFree = true;
    facts.collected = true;
    facts.combinatorialTags = sourceFacts.combinatorialTags;
    facts.maximumFactorsPerTerm = 1;
    facts.factorBases = {basisId};
    facts.pureBasis = basisId;
    facts.expandedBasis = basisId;
    facts.homogeneousWeight = sourceFacts.homogeneousWeight;
    return facts;
  }

// ============================================================================
// Direct Basis-Conversion Plan Catalog
// ============================================================================
// This function is the single contributor-facing catalog. It declares which
// mathematical kernels connect one adjacent basis pair and the exact input
// shape each kernel accepts. It contains no cost or preference policy.

const std::vector<SymmetricEngineRing::BasisConversionPlan>&
SymmetricEngineRing::registeredBasisConversionPlans(
    int sourceBasisId,
    int targetBasisId) const
{
    const std::pair<int, int> key{sourceBasisId, targetBasisId};
    auto found = basisConversionPlanRegistryCache.find(key);
    if (found != basisConversionPlanRegistryCache.end())
      return found->second;
    auto inserted = basisConversionPlanRegistryCache.emplace(
        key,
        buildBasisConversionPlans(sourceBasisId, targetBasisId));
    return inserted.first->second;
  }

std::vector<SymmetricEngineRing::BasisConversionPlan>
SymmetricEngineRing::buildBasisConversionPlans(
    int sourceBasisId,
    int targetBasisId) const
{
    requireBasis(sourceBasisId);
    requireBasis(targetBasisId);
    if (error()) return {};

    std::vector<BasisConversionPlan> plans;
    auto addPlan =
        [&](std::string_view identifier,
            BasisConversionKernel kernel,
            BasisConversionInputShape inputShape =
                BasisConversionInputShape::CanonicalExpansion) {
          plans.push_back({
              identifier,
              sourceBasisId,
              targetBasisId,
              kernel,
              inputShape});
    };

    int powerSumBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return {};

    // Every built-in basis has one broad, policy-free formula to power sums.
    if (targetBasisId == powerSumBasisId &&
        sourceBasisId != powerSumBasisId &&
        basisKindForId(sourceBasisId) != BasisKind::Custom)
      {
        BasisConversionKernel kernel =
            BasisConversionKernel::Unavailable;
        std::string_view identifier;
        switch (basisKindForId(sourceBasisId))
          {
            case BasisKind::Complete:
              kernel =
                  BasisConversionKernel::CompleteToPowerSumsClassical;
              identifier = "Complete->PowerSum:classical-formula";
              break;
            case BasisKind::Elementary:
              kernel =
                  BasisConversionKernel::ElementaryToPowerSumsClassical;
              identifier = "Elementary->PowerSum:classical-formula";
              break;
            case BasisKind::Schur:
              kernel =
                  BasisConversionKernel::SchurToPowerSumsCharacters;
              identifier = "Schur->PowerSum:characters";
              break;
            case BasisKind::SchurOmega:
              kernel =
                  BasisConversionKernel::SchurOmegaToPowerSumsCharacters;
              identifier = "SchurOmega->PowerSum:characters";
              break;
            case BasisKind::Monomial:
              kernel =
                  BasisConversionKernel::MonomialToPowerSumsTransition;
              identifier = "Monomial->PowerSum:transition";
              break;
            case BasisKind::Forgotten:
              kernel =
                  BasisConversionKernel::ForgottenToPowerSumsTransition;
              identifier = "Forgotten->PowerSum:transition";
              break;
            case BasisKind::HallLittlewoodQGenerator:
              kernel = BasisConversionKernel::
                  HallLittlewoodQGeneratorToPowerSumsClassical;
              identifier =
                  "HallLittlewoodQGenerator->PowerSum:classical-formula";
              break;
            case BasisKind::HallLittlewoodBGenerator:
              kernel = BasisConversionKernel::
                  HallLittlewoodBGeneratorToPowerSumsClassical;
              identifier =
                  "HallLittlewoodBGenerator->PowerSum:classical-formula";
              break;
            case BasisKind::HallLittlewoodQ:
              kernel = BasisConversionKernel::
                  HallLittlewoodQToPowerSumsRaising;
              identifier =
                  "HallLittlewoodQ->PowerSum:raising-operators";
              break;
            case BasisKind::HallLittlewoodB:
              kernel = BasisConversionKernel::
                  HallLittlewoodBToPowerSumsRaising;
              identifier =
                  "HallLittlewoodB->PowerSum:raising-operators";
              break;
            case BasisKind::HallLittlewoodP:
              kernel = BasisConversionKernel::
                  HallLittlewoodPToPowerSumsNormalization;
              identifier =
                  "HallLittlewoodP->PowerSum:capital-normalization";
              break;
            case BasisKind::HallLittlewoodPOmega:
              kernel = BasisConversionKernel::
                  HallLittlewoodPOmegaToPowerSumsNormalization;
              identifier =
                  "HallLittlewoodPOmega->PowerSum:capital-normalization";
              break;
            case BasisKind::PowerSum:
            case BasisKind::Custom:
              return plans;
          }
        addPlan(identifier, kernel);
      }

    BasisKind sourceKind = basisKindForId(sourceBasisId);
    BasisKind targetKind = basisKindForId(targetBasisId);

    // Direct changes of normalization or involution avoid an intermediate
    // basis while preserving the canonical-expansion contract.
    if ((sourceKind == BasisKind::HallLittlewoodQ &&
         targetKind == BasisKind::HallLittlewoodP) ||
        (sourceKind == BasisKind::HallLittlewoodP &&
         targetKind == BasisKind::HallLittlewoodQ) ||
        (sourceKind == BasisKind::HallLittlewoodB &&
         targetKind == BasisKind::HallLittlewoodPOmega) ||
        (sourceKind == BasisKind::HallLittlewoodPOmega &&
         targetKind == BasisKind::HallLittlewoodB))
      addPlan(
          "Hall-Littlewood-capital-normalized:diagonal-scaling",
          BasisConversionKernel::HallLittlewoodNormalization);
    if ((sourceKind == BasisKind::Schur &&
         targetKind == BasisKind::SchurOmega) ||
        (sourceKind == BasisKind::SchurOmega &&
         targetKind == BasisKind::Schur))
      addPlan(
          "Schur-omega:partition-conjugation",
          BasisConversionKernel::SchurOmegaConjugation);
    if (sourceKind == BasisKind::Complete &&
        targetKind == BasisKind::Schur)
      addPlan(
          "h->S:recursive-transition",
          BasisConversionKernel::CompleteToSchurRecursive);
    bool hallLittlewoodGeneratorTransition =
        ((sourceKind == BasisKind::HallLittlewoodQGenerator) &&
         (targetKind == BasisKind::HallLittlewoodQ ||
          targetKind == BasisKind::HallLittlewoodP)) ||
        ((sourceKind == BasisKind::HallLittlewoodBGenerator) &&
         (targetKind == BasisKind::HallLittlewoodB ||
          targetKind == BasisKind::HallLittlewoodPOmega));
    if (hallLittlewoodGeneratorTransition)
      {
        // Plan identifiers are immutable registry labels. Spell out the four
        // edges so the registry can store non-owning string views safely.
        std::string_view identifier;
        if (sourceKind == BasisKind::HallLittlewoodQGenerator)
          identifier =
              targetKind == BasisKind::HallLittlewoodQ
                  ? "HallLittlewoodQGenerator->HallLittlewoodQ:"
                    "triangular-reduction"
                  : "HallLittlewoodQGenerator->HallLittlewoodP:"
                    "triangular-reduction";
        else
          identifier =
              targetKind == BasisKind::HallLittlewoodB
                  ? "HallLittlewoodBGenerator->HallLittlewoodB:"
                    "triangular-reduction"
                  : "HallLittlewoodBGenerator->HallLittlewoodPOmega:"
                    "triangular-reduction";
        addPlan(
            identifier,
            BasisConversionKernel::
                HallLittlewoodGeneratorToCapitalTriangular);
      }

    // All remaining target families start from power sums.
    if (sourceBasisId != powerSumBasisId ||
        targetKind == BasisKind::Custom)
      return plans;

    switch (targetKind)
      {
        case BasisKind::PowerSum:
          break;
        case BasisKind::Complete:
          addPlan("p->h:logarithm-formula",
                  BasisConversionKernel::
                      PowerSumsToCompleteLogarithm);
          break;
        case BasisKind::Elementary:
          addPlan("p->e:logarithm-formula",
                  BasisConversionKernel::
                      PowerSumsToElementaryLogarithm);
          break;
        case BasisKind::Schur:
          addPlan("p->S:abacus-rim-hooks",
                  BasisConversionKernel::
                      PowerSumsToSchurAbacusRimHooks);
          addPlan("p->S:abacus-rim-hooks+complete",
                  BasisConversionKernel::
                      PowerSumsToSchurAbacusAndComplete,
                  BasisConversionInputShape::
                      MixedCompleteFriendlyPowerSumExpansion);
          addPlan("p->h->S:recursive-transition",
                  BasisConversionKernel::
                      PowerSumsToSchurComplete);
          addPlan("p->S:grouped-characters",
                  BasisConversionKernel::
                      PowerSumsToSchurCharacters);
          addPlan("p->S:border-strips",
                  BasisConversionKernel::
                      PowerSumsToSchurBorderStrips);
          break;
        case BasisKind::SchurOmega:
          addPlan("p->omega(p)->S->Somega:abacus-rim-hooks",
                  BasisConversionKernel::
                      PowerSumsToSchurOmegaAbacusRimHooks);
          addPlan("p->omega(p)->S->Somega:abacus-rim-hooks+complete",
                  BasisConversionKernel::
                      PowerSumsToSchurOmegaAbacusAndComplete,
                  BasisConversionInputShape::
                      MixedCompleteFriendlyPowerSumExpansion);
          addPlan("p->omega(p)->h->S->Somega",
                  BasisConversionKernel::
                      PowerSumsToSchurOmegaComplete);
          addPlan("p->omega(p)->S->Somega:grouped-characters",
                  BasisConversionKernel::
                      PowerSumsToSchurOmegaCharacters);
          addPlan("p->omega(p)->S->Somega:border-strips",
                  BasisConversionKernel::
                      PowerSumsToSchurOmegaBorderStrips);
          break;
        case BasisKind::HallLittlewoodQGenerator:
          addPlan("p->q:logarithm-formula",
                  BasisConversionKernel::
                      PowerSumsToHallLittlewoodQGeneratorLogarithm);
          break;
        case BasisKind::HallLittlewoodBGenerator:
          addPlan("p->b:logarithm-formula",
                  BasisConversionKernel::
                      PowerSumsToHallLittlewoodBGeneratorLogarithm);
          break;
        case BasisKind::HallLittlewoodQ:
        case BasisKind::HallLittlewoodB:
        case BasisKind::HallLittlewoodP:
        case BasisKind::HallLittlewoodPOmega:
          addPlan(
              "single-cycle-p-terms->Hall-Littlewood:Green-polynomials",
              BasisConversionKernel::
                  PowerSumSingleCyclesToHallLittlewoodGreen,
              BasisConversionInputShape::
                  SingleCyclePowerSumExpansion);
          addPlan("p->Hall-Littlewood:triangular-reduction",
                  BasisConversionKernel::
                      PowerSumsToHallLittlewoodTriangular);
          addPlan(
              "p_mu->Hall-Littlewood:Green-polynomials-via-duality",
              BasisConversionKernel::
                  PowerSumIndexToHallLittlewoodGreenDuality,
              BasisConversionInputShape::SingleBasisElement);
          break;
        case BasisKind::Monomial:
          addPlan("p->m:transition",
                  BasisConversionKernel::
                      PowerSumsToMonomialTransition);
          break;
        case BasisKind::Forgotten:
          addPlan("p->ff:transition",
                  BasisConversionKernel::
                      PowerSumsToForgottenTransition);
          break;
        case BasisKind::Custom:
          return plans;
      }
    return plans;
  }

// ============================================================================
// Direct-Plan Applicability Contracts
// ============================================================================
// Applicability is mathematical, not a performance preference. Once this
// returns true, execution must either complete or report a mathematical/engine
// error; it must never select a replacement plan.

bool SymmetricEngineRing::basisConversionPlanApplicable(
    const BasisConversionPlan& plan,
    const ExpressionFacts& facts) const
{
    if (plan.kernel == BasisConversionKernel::Unavailable ||
        !facts.canonicalExpansionInBasis(plan.sourceBasisId))
      return false;

    switch (plan.inputShape)
      {
        case BasisConversionInputShape::CanonicalExpansion:
          return true;
        case BasisConversionInputShape::SingleBasisElement:
          return facts.singleTerm() && facts.singleBasisElement();
        case BasisConversionInputShape::
            SingleCyclePowerSumExpansion:
          return facts.allPowerSumTermsSingleCycles.value_or(false);
        case BasisConversionInputShape::
            MixedCompleteFriendlyPowerSumExpansion:
          {
            size_t completeFriendly =
                facts.completeFriendlyPowerSumTermCount.value_or(0);
            return facts.homogeneousWeight.has_value() &&
                   completeFriendly > 0 &&
                   completeFriendly < facts.termCount;
          }
      }
    return false;
  }

// ============================================================================
// Direct-Plan Cost Policy
// ============================================================================
// These structural estimates only rank mathematically applicable plans.
// Numerical tuning is intentionally deferred until the workflows are
// benchmarked as complete units.

size_t SymmetricEngineRing::basisConversionPlanCost(
    const BasisConversionPlan& plan,
    const ExpressionFacts& facts,
    CombinatorialTags combinatorialTags) const
{
    size_t terms = std::max<size_t>(1, facts.termCount);
    size_t weight = static_cast<size_t>(
        std::max(0, facts.homogeneousWeight.value_or(0)));
    size_t cost = 100 + terms + weight;
    switch (plan.kernel)
      {
        case BasisConversionKernel::Unavailable:
          return std::numeric_limits<size_t>::max();
        case BasisConversionKernel::HallLittlewoodNormalization:
        case BasisConversionKernel::SchurOmegaConjugation:
        case BasisConversionKernel::CompleteToSchurRecursive:
        case BasisConversionKernel::
            HallLittlewoodGeneratorToCapitalTriangular:
          return 20 + terms;
        case BasisConversionKernel::PowerSumsToSchurCharacters:
        case BasisConversionKernel::
            PowerSumsToSchurOmegaCharacters:
          cost += 4 * terms + weight;
          break;
        case BasisConversionKernel::
            PowerSumsToSchurAbacusAndComplete:
        case BasisConversionKernel::
            PowerSumsToSchurOmegaAbacusAndComplete:
          cost += 2 * terms;
          break;
        default:
          break;
      }
    if (hasCombinatorialTag(
            combinatorialTags, CombinatorialTag::Plethysm) &&
        cost > 10)
      cost -= 10;
    return cost;
  }

// ============================================================================
// Direct-Plan Picker
// ============================================================================
// Pair-local policy chooses one applicable plan. Composition policy remains in
// pickBasisConversionPlans so registry entries never choose intermediates.

SymmetricEngineRing::BasisConversionPlan
SymmetricEngineRing::selectBasisConversionPlan(
    int sourceBasisId,
    int targetBasisId,
    const ExpressionFacts& facts,
    CombinatorialTags combinatorialTags,
    const std::optional<std::string>& forcedIdentifier) const
{
    const auto& plans =
        registeredBasisConversionPlans(sourceBasisId, targetBasisId);

    if (forcedIdentifier)
      {
        for (const auto& plan : plans)
          if (plan.identifier == *forcedIdentifier &&
              basisConversionPlanApplicable(plan, facts))
            return plan;
        ERROR("the requested conversion plan is unknown or inapplicable: ",
              forcedIdentifier->c_str());
        return {"unavailable",
                sourceBasisId,
                targetBasisId,
                BasisConversionKernel::Unavailable,
                BasisConversionInputShape::CanonicalExpansion};
      }

    BasisKind targetKind = basisKindForId(targetBasisId);
    if (sourceBasisId == requiredBasisIdForKind(BasisKind::PowerSum) &&
        (targetKind == BasisKind::Schur ||
         targetKind == BasisKind::SchurOmega))
      {
        auto findKernel =
            [&](BasisConversionKernel kernel)
                -> std::optional<BasisConversionPlan> {
              for (const auto& plan : plans)
                if (plan.kernel == kernel &&
                    basisConversionPlanApplicable(plan, facts))
                  return plan;
              return std::nullopt;
            };
        const bool omega = targetKind == BasisKind::SchurOmega;
        const size_t termCount = facts.termCount;
        const size_t possibleTerms = facts.possibleTermCount.value_or(
            partitionCountForConversionSelection(
                facts.homogeneousWeight.value_or(0)));
        const long double termsSquared =
            static_cast<long double>(termCount) *
            static_cast<long double>(termCount);
        const bool supportSquareFavorsComplete =
            2.0L * termsSquared >=
                35.0L * static_cast<long double>(possibleTerms);
        const int weight = facts.homogeneousWeight.value_or(0);
        const double density = possibleTerms == 0
            ? 0.0
            : static_cast<double>(termCount) /
                  static_cast<double>(possibleTerms);
        bool plethysm =
            hasCombinatorialTag(combinatorialTags,
                                CombinatorialTag::Plethysm);
        const CombinatorialTags borderStripsTag =
            combinatorialTagMask(CombinatorialTag::BorderStrips);
        size_t completeFriendly =
            facts.completeFriendlyPowerSumTermCount.value_or(0);
        size_t minimumCompleteGroup =
            std::max<size_t>(8, (termCount + 3) / 4);

        BasisConversionKernel preferredKernel;
        // This decision order and every predicate mirror
        // selectPowerSumsToSchurRoute.  Keeping the policy equal is essential:
        // benchmarks should measure workflow overhead, not different kernels.
        bool completePreferred = false;
        if (plethysm)
          completePreferred =
              supportSquareFavorsComplete ||
              (weight >= 8 && termCount >= 2);
        else if (hasCombinatorialTag(
                     combinatorialTags,
                     CombinatorialTag::LittlewoodRichardson))
          completePreferred =
              weight >= 7 &&
              (coefficientRing != globalQQ || density >= 0.25 ||
               supportSquareFavorsComplete);
        else if (hasCombinatorialTag(
                     combinatorialTags,
                     CombinatorialTag::HorizontalPieri) ||
                 hasCombinatorialTag(
                     combinatorialTags,
                     CombinatorialTag::VerticalPieri))
          completePreferred =
              weight >= 8 &&
              (coefficientRing != globalQQ || density >= 0.25 ||
               supportSquareFavorsComplete);

        static const std::vector<int> noCommonParts;
        const auto& commonParts = facts.commonPowerSumParts
            ? *facts.commonPowerSumParts
            : noCommonParts;
        bool hasSmallCommonCycle = std::any_of(
            commonParts.begin(),
            commonParts.end(),
            [&](int part) {
              int crossoverPercent =
                  coefficientRing == globalQQ ? 36 : 41;
              return 100 * part <= crossoverPercent * weight;
            });
        bool hasCommonOne =
            std::find(commonParts.begin(), commonParts.end(), 1) !=
            commonParts.end();
        if (!completePreferred &&
            combinatorialTags == borderStripsTag &&
            weight >= 16 && termCount >= 8 &&
            !hasCommonOne && hasSmallCommonCycle)
          completePreferred = true;
        if (!completePreferred &&
            combinatorialTags == borderStripsTag &&
            weight >= 8 && termCount >= 8 && hasCommonOne &&
            (coefficientRing != globalQQ || density >= 0.25 ||
             supportSquareFavorsComplete))
          completePreferred = true;
        if (!completePreferred &&
            (combinatorialTags == 0 ||
             combinatorialTags == borderStripsTag) &&
            weight >= 14 && termCount >= 2 &&
            completeFriendly == termCount)
          completePreferred = true;

        bool hybridPreferred =
            !completePreferred &&
            (combinatorialTags == 0 ||
             combinatorialTags == borderStripsTag) &&
            weight >= 14 && termCount >= 8 &&
            completeFriendly >= minimumCompleteGroup &&
            completeFriendly < termCount;

        if (completePreferred)
          preferredKernel =
              omega
                  ? BasisConversionKernel::
                        PowerSumsToSchurOmegaComplete
                  : BasisConversionKernel::
                        PowerSumsToSchurComplete;
        else if (hybridPreferred)
          preferredKernel =
              omega
                  ? BasisConversionKernel::
                        PowerSumsToSchurOmegaAbacusAndComplete
                  : BasisConversionKernel::
                        PowerSumsToSchurAbacusAndComplete;
        else
          preferredKernel =
              omega
                  ? BasisConversionKernel::
                        PowerSumsToSchurOmegaAbacusRimHooks
                  : BasisConversionKernel::
                        PowerSumsToSchurAbacusRimHooks;
        auto preferred = findKernel(preferredKernel);
        if (preferred) return *preferred;
      }

    if (sourceBasisId == requiredBasisIdForKind(BasisKind::PowerSum) &&
        (targetKind == BasisKind::HallLittlewoodQ ||
         targetKind == BasisKind::HallLittlewoodB ||
         targetKind == BasisKind::HallLittlewoodP ||
         targetKind == BasisKind::HallLittlewoodPOmega))
      {
        BasisConversionKernel preferredKernel =
            facts.allPowerSumTermsSingleCycles.value_or(false)
                ? BasisConversionKernel::
                      PowerSumSingleCyclesToHallLittlewoodGreen
                : (facts.singleTerm() && facts.singleBasisElement()
                       ? BasisConversionKernel::
                             PowerSumIndexToHallLittlewoodGreenDuality
                       : BasisConversionKernel::
                             PowerSumsToHallLittlewoodTriangular);
        for (const auto& plan : plans)
          if (plan.kernel == preferredKernel &&
              basisConversionPlanApplicable(plan, facts))
            return plan;
      }

    const BasisConversionPlan *best = nullptr;
    size_t bestCost = std::numeric_limits<size_t>::max();
    for (const auto& plan : plans)
      {
        if (!basisConversionPlanApplicable(plan, facts))
          continue;
        const size_t cost = basisConversionPlanCost(
            plan, facts, combinatorialTags);
        if (best == nullptr || cost < bestCost)
          {
            best = &plan;
            bestCost = cost;
          }
      }
    if (best != nullptr) return *best;
    ERROR("the basis-conversion registry has no applicable plan from ",
          basisKeyForId(sourceBasisId).c_str(),
          " to ",
          basisKeyForId(targetBasisId).c_str());
    return {"unavailable",
            sourceBasisId,
            targetBasisId,
            BasisConversionKernel::Unavailable,
            BasisConversionInputShape::CanonicalExpansion};
  }

// ============================================================================
// Basis-Composition Picker
// ============================================================================
// The picker returns the entire ordered composition before execution. Exact
// forcing controls are contract-test bypasses; automatic mode searches the
// registered basis graph using guarantees for unexecuted intermediates. A
// data-dependent plan chosen from those guarantees is provisional and is
// reselected from exact facts when that intermediate has been realized.

SymmetricEngineRing::BasisConversionPlanSelection
SymmetricEngineRing::pickBasisConversionPlans(
    int sourceBasisId,
    int targetBasisId,
    const ExpressionFacts& facts,
    CombinatorialTags combinatorialTags) const
{
    BasisConversionPlanSelection selection;
    if (!facts.canonicalExpansionInBasis(sourceBasisId))
      {
        ERROR("the basis-composition picker requires a canonical "
              "single-source-basis expansion");
        return selection;
      }
    int powerSumBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return selection;

    auto basisIdFromKey = [&](const std::string& key) {
      for (const auto& item : basisDescriptors)
        if (item.second.canonicalKey == key) return item.first;
      return -1;
    };
    auto split = [](const std::string& value, char delimiter) {
      std::vector<std::string> pieces;
      std::stringstream stream(value);
      std::string piece;
      while (std::getline(stream, piece, delimiter))
        pieces.push_back(piece);
      return pieces;
    };
    auto hasApplicablePlan =
        [&](int from,
            int to,
            const ExpressionFacts& edgeFacts) {
          const auto& plans =
              registeredBasisConversionPlans(from, to);
          return std::any_of(
              plans.begin(),
              plans.end(),
              [&](const BasisConversionPlan& plan) {
                return basisConversionPlanApplicable(
                    plan, edgeFacts);
              });
        };

    // Stage 1: choose or validate the complete basis composition.
    std::vector<int> composition;
    const char *forcedComposition =
        std::getenv(
            "M2_SYMMETRIC_RINGS_FORCE_CONVERSION_COMPOSITION");
    if (forcedComposition != nullptr)
      {
        std::string requested(forcedComposition);
        if (requested == "direct")
          composition = sourceBasisId == targetBasisId
              ? std::vector<int>{sourceBasisId}
              : std::vector<int>{sourceBasisId, targetBasisId};
        else if (requested == "via-power-sums")
          composition =
              {sourceBasisId, powerSumBasisId, targetBasisId};
        else
          {
            for (const auto& key : split(requested, '>'))
              {
                std::string normalizedKey = key;
                if (!normalizedKey.empty() &&
                    normalizedKey.back() == '-')
                  normalizedKey.pop_back();
                int basisId = basisIdFromKey(normalizedKey);
                if (basisId < 0)
                  {
                    ERROR("the requested conversion composition "
                          "contains an unknown basis: ",
                          normalizedKey.c_str());
                    return selection;
                  }
                composition.push_back(basisId);
              }
          }
        if (composition.empty() ||
            composition.front() != sourceBasisId ||
            composition.back() != targetBasisId ||
            (composition.size() == 1 &&
             sourceBasisId != targetBasisId))
          {
            ERROR("the requested conversion composition has the "
                  "wrong endpoints");
            return selection;
          }
        for (size_t i = 1; i < composition.size(); ++i)
          if (composition[i - 1] == composition[i])
            {
              ERROR("a requested conversion composition contains "
                    "a repeated adjacent basis");
              return selection;
            }
      }
    else if (sourceBasisId == targetBasisId)
      composition = {sourceBasisId};
    else if (hasApplicablePlan(sourceBasisId, targetBasisId, facts))
      {
        // Avoid allocating and searching the complete basis graph when an
        // applicable direct kernel already supplies the preferred edge.
        composition = {sourceBasisId, targetBasisId};
      }
    else
      {
        const size_t infinity = std::numeric_limits<size_t>::max();
        std::map<int, size_t> distances;
        std::map<int, int> previous;
        using QueueEntry = std::pair<size_t, int>;
        std::priority_queue<
            QueueEntry,
            std::vector<QueueEntry>,
            std::greater<QueueEntry>> queue;
        for (const auto& item : basisDescriptors)
          distances[item.first] = infinity;
        distances[sourceBasisId] = 0;
        queue.push({0, sourceBasisId});

        while (!queue.empty())
          {
            auto current = queue.top();
            queue.pop();
            size_t distance = current.first;
            int from = current.second;
            if (distance != distances[from]) continue;
            if (from == targetBasisId) break;
            ExpressionFacts conservativeFacts;
            const ExpressionFacts *edgeFacts = &facts;
            if (from != sourceBasisId)
              {
                conservativeFacts =
                    conservativeExpressionFactsForBasis(
                        facts, from);
                edgeFacts = &conservativeFacts;
              }
            for (const auto& item : basisDescriptors)
              {
                int to = item.first;
                if (to == from ||
                    basisKindForId(to) == BasisKind::Custom)
                  continue;
                const BasisConversionPlan *best = nullptr;
                size_t bestPlanCost =
                    std::numeric_limits<size_t>::max();
                for (const auto& plan :
                     registeredBasisConversionPlans(from, to))
                  {
                    if (!basisConversionPlanApplicable(
                            plan, *edgeFacts))
                      continue;
                    const size_t planCost =
                        basisConversionPlanCost(
                            plan, *edgeFacts, combinatorialTags);
                    if (planCost < bestPlanCost)
                      {
                        best = &plan;
                        bestPlanCost = planCost;
                      }
                  }
                if (best == nullptr) continue;
                size_t edgeCost = 1000 + bestPlanCost;
                if (distance <= infinity - edgeCost &&
                    distance + edgeCost < distances[to])
                  {
                    distances[to] = distance + edgeCost;
                    previous[to] = from;
                    queue.push({distances[to], to});
                  }
              }
          }
        if (!distances.count(targetBasisId) ||
            distances[targetBasisId] == infinity)
          {
            ERROR("the basis-conversion registry has no composition "
                  "for the requested basis pair");
            return selection;
          }
        for (int current = targetBasisId;; current = previous[current])
          {
            composition.push_back(current);
            if (current == sourceBasisId) break;
          }
        std::reverse(composition.begin(), composition.end());
      }

    // Stage 2: validate optional ordered plan forcing for the chosen edges.
    std::vector<std::string> orderedForcedPlans;
    const char *forcedPlans =
        std::getenv(
            "M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLANS");
    if (forcedPlans != nullptr)
      orderedForcedPlans = split(forcedPlans, ',');
    if (!orderedForcedPlans.empty() &&
        orderedForcedPlans.size() + 1 != composition.size())
      {
        ERROR("the requested ordered conversion plans do not match "
              "the selected composition length");
        return selection;
      }

    const char *singleForcedPlan =
        std::getenv("M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN");
    bool singleForcedPlanUsed = singleForcedPlan == nullptr;
    // Stage 3: select the exact first edge and provisional later edges.
    selection.basisComposition = composition;
    const ExpressionFacts *edgeFacts = &facts;
    ExpressionFacts laterEdgeFacts;
    for (size_t i = 1; i < composition.size(); ++i)
      {
        std::optional<std::string> forcedIdentifier;
        if (!orderedForcedPlans.empty())
          forcedIdentifier = orderedForcedPlans[i - 1];
        else if (singleForcedPlan != nullptr)
          {
            for (const auto& plan :
                 registeredBasisConversionPlans(
                     composition[i - 1], composition[i]))
              if (plan.identifier == singleForcedPlan &&
                  basisConversionPlanApplicable(plan, *edgeFacts))
                {
                  forcedIdentifier = singleForcedPlan;
                  singleForcedPlanUsed = true;
                  break;
                }
          }
        BasisConversionPlan plan =
            selectBasisConversionPlan(
                composition[i - 1],
                composition[i],
                *edgeFacts,
                combinatorialTags,
                forcedIdentifier);
        if (error()) return selection;
        selection.plans.push_back(std::move(plan));
        selection.reselectFromExactIntermediate.push_back(
            i > 1 && !forcedIdentifier.has_value());
        laterEdgeFacts = conservativeExpressionFactsForBasis(
            facts, composition[i]);
        edgeFacts = &laterEdgeFacts;
      }
    if (!singleForcedPlanUsed)
      {
        ERROR("the requested conversion plan is not part of the "
              "selected composition: ",
              singleForcedPlan);
        return selection;
      }
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr)
      {
        std::string compositionTrace;
        for (size_t i = 0; i < composition.size(); ++i)
          {
            if (i != 0) compositionTrace += "->";
            compositionTrace += basisKeyForId(composition[i]);
          }
        std::fprintf(
            stderr,
            "SymmetricRings conversion-selection: composition=%s "
            "plans=%zu terms=%zu weight=%d density=%.6f\n",
            compositionTrace.c_str(),
            selection.plans.size(),
            facts.termCount,
            facts.homogeneousWeight.value_or(-1),
            facts.density.value_or(-1.0));
      }
    return selection;
  }

// ============================================================================
// Composite Schur Kernels
// ============================================================================
// These kernels combine the complete-function transition with the direct
// Schur kernels. Selection policy remains in the plan picker above.

ring_elem SymmetricEngineRing::powerSumsToSchurViaComplete(
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetOrder) const
{
    CoeffMap inComplete = powerSumsToCompleteMapViaLogarithmFormula(f);
    if (error()) return zero();
    CoeffMap result =
        completeToSchurCoefficientsViaRecursiveTransition(inComplete);
    if (error()) return zero();
    return coeffMapToElement(
        result, targetBasisId, targetDisplay, targetOrder, false);
}

ring_elem SymmetricEngineRing::powerSumsToSchurViaAbacusAndComplete(
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetOrder) const
{
    CoeffMap abacusResult;
    CoeffMap completeInput;
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during hybrid basis conversion");
            return zero();
          }
        const int weight = partitionWeight(index);
        if (completeFriendlyPowerSumIndexForConversionSelection(index, weight))
          addPowerSumIndexToCompleteMapViaLogarithmFormula(
              index, term.coeff, completeInput);
        else
          addPowerSumIndexToSchurMapViaAbacusRimHooks(
              index, term.coeff, abacusResult);
      }

    if (!completeInput.empty())
      {
        CoeffMap completeResult =
            completeToSchurCoefficientsViaRecursiveTransition(completeInput);
        if (error()) return zero();
        for (const auto& term : completeResult)
          addNormalizedCoeff(abacusResult, term.first, term.second);
      }
    return coeffMapToElement(
        abacusResult, targetBasisId, targetDisplay, targetOrder, false);
}

// ============================================================================
// Direct-Plan Executor
// ============================================================================
// This executor contains no applicability or performance policy. The switch is
// grouped by mathematical source/target family so a contributor can locate the
// one kernel associated with a registry entry.

ring_elem SymmetricEngineRing::executeBasisConversionPlan(
    const BasisConversionPlan& plan,
    ring_elem expression,
    CombinatorialTags combinatorialTags,
    ExpressionFacts *resultFacts,
    const std::optional<int>& knownHomogeneousWeight) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr)
      std::fprintf(stderr,
                   "SymmetricRings conversion-plan: source=%s target=%s "
                   "plan=%s\n",
                   displayForBasis(plan.sourceBasisId).c_str(),
                   displayForBasis(plan.targetBasisId).c_str(),
                   plan.identifier.data());

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
                  PowerSumsToSchurOmegaAbacusAndComplete:
                inSchur = powerSumsToSchurViaAbacusAndComplete(
                    omegaInput, schurId, schurDisplay, schurOrder);
                break;
              case BasisConversionKernel::
                  PowerSumsToSchurOmegaComplete:
                inSchur = powerSumsToSchurViaComplete(
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
              inSchur, schurId, plan.targetBasisId);
        };

    const std::string targetDisplay =
        displayForBasis(plan.targetBasisId);
    const int targetOrder =
        basisOrderForId(plan.targetBasisId);
    ring_elem result;
    switch (plan.kernel)
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
              basisKindForId(plan.sourceBasisId));
          break;

        // Classical power-sum targets.
        case BasisConversionKernel::PowerSumsToCompleteLogarithm:
          result = powerSumsToCompleteViaLogarithmFormula(
              expression, plan.targetBasisId, targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToElementaryLogarithm:
          result = powerSumsToElementaryViaLogarithmFormula(
              expression, plan.targetBasisId, targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToSchurBorderStrips:
          result = powerSumsToSchurViaBorderStrips(
              expression,
              plan.targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToSchurAbacusRimHooks:
          result = powerSumsToSchurViaAbacusRimHooks(
              expression,
              plan.targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToSchurAbacusAndComplete:
          result = powerSumsToSchurViaAbacusAndComplete(
              expression,
              plan.targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToSchurComplete:
          result = powerSumsToSchurViaComplete(
              expression,
              plan.targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::PowerSumsToSchurCharacters:
          result = powerSumsToSchurLikeViaCharacters(
              expression,
              plan.targetBasisId,
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
            PowerSumsToSchurOmegaAbacusAndComplete:
        case BasisConversionKernel::
            PowerSumsToSchurOmegaComplete:
        case BasisConversionKernel::
            PowerSumsToSchurOmegaCharacters:
          result = powerSumsToSchurOmega(plan.kernel);
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
                    plan.kernel == BasisConversionKernel::
                        PowerSumsToHallLittlewoodBGeneratorLogarithm);
            result = coeffMapToElement(
                generators,
                plan.targetBasisId,
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
                  plan.targetBasisId,
                  targetDisplay,
                  targetOrder);
          break;
        case BasisConversionKernel::
            PowerSumIndexToHallLittlewoodGreenDuality:
          result =
              powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality(
                  expression,
                  plan.targetBasisId,
                  targetDisplay,
                  targetOrder);
          break;
        case BasisConversionKernel::
            PowerSumsToHallLittlewoodTriangular:
          result = powerSumsToHallLittlewoodViaTriangularReduction(
              expression,
              plan.targetBasisId,
              targetDisplay,
              targetOrder);
          break;

        // Monomial and forgotten transition matrices.
        case BasisConversionKernel::PowerSumsToMonomialTransition:
        case BasisConversionKernel::PowerSumsToForgottenTransition:
          result = powerSumsToTargetViaTermwiseConversion(
              expression,
              plan.targetBasisId,
              targetDisplay,
              targetOrder,
              isMultiplicativeBasis(plan.targetBasisId));
          break;

        // Direct normalization, involution, and triangular transitions.
        case BasisConversionKernel::HallLittlewoodNormalization:
          {
            BasisKind sourceKind =
                basisKindForId(plan.sourceBasisId);
            result =
                hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
                    expression,
                    plan.sourceBasisId,
                    plan.targetBasisId,
                    targetDisplay,
                    targetOrder,
                    sourceKind == BasisKind::HallLittlewoodQ ||
                        sourceKind == BasisKind::HallLittlewoodB);
            break;
          }
        case BasisConversionKernel::SchurOmegaConjugation:
          result = schurOmegaConversionViaPartitionConjugation(
              expression,
              plan.sourceBasisId,
              plan.targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::CompleteToSchurRecursive:
          result = completeToSchurViaRecursiveTransition(
              expression,
              plan.sourceBasisId,
              displayForBasis(plan.sourceBasisId),
              basisOrderForId(plan.sourceBasisId),
              isMultiplicativeBasis(plan.sourceBasisId),
              plan.targetBasisId,
              targetDisplay,
              targetOrder);
          break;
        case BasisConversionKernel::
            HallLittlewoodGeneratorToCapitalTriangular:
          if (!tryExpressionToHallLittlewoodViaTriangularReduction(
                  expression,
                  plan.targetBasisId,
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
            plan.targetBasisId,
            knownHomogeneousWeight);
        resultFacts->combinatorialTags = combinatorialTags;
        if (!resultFacts->canonicalExpansionInBasis(plan.targetBasisId))
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
// Ordered-Composition Workflow
// ============================================================================
// The direct plan executor above is policy-free. This owning workflow performs
// the one policy action that cannot be correct before execution: selecting a
// data-dependent later-edge kernel from the exact intermediate support.

ring_elem SymmetricEngineRing::executeBasisConversionPlans(
    const BasisConversionPlanSelection& selection,
    ring_elem expression,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts& inputFacts,
    ExpressionFacts *resultFacts) const
{
    ring_elem current = expression;
    const ExpressionFacts *facts = &inputFacts;
    ExpressionFacts realizedFacts;
    if (selection.basisComposition.empty() ||
        selection.plans.size() + 1 !=
            selection.basisComposition.size() ||
        selection.reselectFromExactIntermediate.size() !=
            selection.plans.size() ||
        !facts->canonicalExpansionInBasis(
            selection.basisComposition.front()))
      {
        ERROR("an invalid basis-conversion selection reached the executor");
        return zero();
      }
    for (size_t i = 0; i < selection.plans.size(); ++i)
      {
        BasisConversionPlan plan = selection.plans[i];
        if (selection.reselectFromExactIntermediate[i])
          {
            enrichConversionSelectionFacts(
                current,
                selection.basisComposition[i + 1],
                realizedFacts);
            plan = selectBasisConversionPlan(
                selection.basisComposition[i],
                selection.basisComposition[i + 1],
                *facts,
                combinatorialTags);
            if (error()) return zero();
          }
        if (plan.sourceBasisId != selection.basisComposition[i] ||
            plan.targetBasisId != selection.basisComposition[i + 1])
          {
            ERROR("a basis-conversion plan does not match its selected "
                  "composition edge");
            return zero();
          }
        if (!facts->canonicalExpansionInBasis(plan.sourceBasisId))
          {
            ERROR("a basis-conversion plan received input outside its "
                  "canonical source-basis contract");
            return zero();
          }
        const bool needsRealizedFacts =
            i + 1 < selection.plans.size() || resultFacts != nullptr;
        current = executeBasisConversionPlan(
            plan,
            current,
            combinatorialTags,
            needsRealizedFacts ? &realizedFacts : nullptr,
            facts->homogeneousWeight);
        if (error()) return zero();
        if (needsRealizedFacts) facts = &realizedFacts;
      }
    if (resultFacts != nullptr)
      {
        if (selection.plans.empty())
          *resultFacts = inputFacts;
        else
          *resultFacts = std::move(realizedFacts);
      }
    return current;
  }

// ============================================================================
// Canonical Group Conversion Helper
// ============================================================================
// This helper owns grouping, weight-block selection, execution, and collection
// for input that is already normalized, skew-free, and product-free.

ring_elem SymmetricEngineRing::convertCanonicalExpressionToBasis(
    ring_elem f,
    int targetBasisId,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts *knownFacts) const
{
    const bool trace =
        std::getenv(
            "M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr;
    BasisKind targetKind = basisKindForId(targetBasisId);

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

    // A pure canonical expansion can enter the picker directly. The sole
    // exception is nonhomogeneous p-to-Schur conversion, whose weight blocks
    // deliberately receive independent selections below.
    if (knownFacts != nullptr && knownFacts->expandedBasis)
      {
        int sourceBasisId = *knownFacts->expandedBasis;
        bool needsPowerSumWeightBlocks =
            basisKindForId(sourceBasisId) == BasisKind::PowerSum &&
            (targetKind == BasisKind::Schur ||
             targetKind == BasisKind::SchurOmega) &&
            !knownFacts->homogeneousWeight;
        if (!needsPowerSumWeightBlocks)
          {
            // Most direct conversions consume only the exact canonical core.
            // Copy the comparatively rich facts object only for power-sum
            // policies that attach a transient selector profile.
            ExpressionFacts enrichedFacts;
            const ExpressionFacts *selectionFacts = knownFacts;
            if (basisKindForId(sourceBasisId) == BasisKind::PowerSum)
              {
                enrichedFacts = *knownFacts;
                enrichConversionSelectionFacts(
                    f, targetBasisId, enrichedFacts);
                selectionFacts = &enrichedFacts;
              }
            BasisConversionPlanSelection selection =
                pickBasisConversionPlans(
                    sourceBasisId,
                    targetBasisId,
                    *selectionFacts,
                    combinatorialTags);
            if (error()) return zero();
            ExpressionFacts resultFacts;
            ring_elem result = executeBasisConversionPlans(
                selection,
                f,
                combinatorialTags,
                *selectionFacts,
                &resultFacts);
            if (error()) return zero();
            if (selection.plans.empty())
              {
                result = copyPolyValue(polyValue(result));
                resultFacts = *knownFacts;
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
      }

    struct PreparedConversion
    {
      ring_elem expression;
      ExpressionFacts facts;
      BasisConversionPlanSelection selection;
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
      BasisConversionPlanSelection selection =
          pickBasisConversionPlans(
              sourceBasisId,
              targetBasisId,
              facts,
              combinatorialTags);
      if (error()) return;
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
        enrichConversionSelectionFacts(
            groupExpression, targetBasisId, groupFacts);
        bool splitPowerSumWeights =
            basisKindForId(group.first) == BasisKind::PowerSum &&
            (targetKind == BasisKind::Schur ||
             targetKind == BasisKind::SchurOmega) &&
            !groupFacts.homogeneousWeight;
        if (!splitPowerSumWeights)
          {
            prepareGroup(
                group.first, groupExpression, std::move(groupFacts));
            if (error()) return zero();
            continue;
          }

        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings conversion-stage: source=%s target=%s "
              "stage=weight-blocks blocks=%zu\n",
              displayForBasis(group.first).c_str(),
              displayForBasis(targetBasisId).c_str(),
              groupFacts.weightTermPositions.size());
        for (const auto& weightPositions :
             groupFacts.weightTermPositions)
          {
            VECTOR(SymmetricTerm) blockTerms;
            blockTerms.reserve(weightPositions.second.size());
            for (size_t position : weightPositions.second)
              {
                if (position >= group.second.size())
                  {
                    ERROR("cached weight-block metadata is invalid");
                    return zero();
                  }
                blockTerms.push_back(group.second[position]);
              }
            ring_elem blockExpression =
                fromTermVector(blockTerms, true);
            ExpressionFacts blockFacts =
                inferCanonicalExpansionFacts(
                    blockExpression,
                    group.first,
                    weightPositions.first);
            if (!blockFacts.homogeneousWeight)
              {
                ERROR("a conversion weight block was not homogeneous");
                return zero();
              }
            enrichConversionSelectionFacts(
                blockExpression, targetBasisId, blockFacts);
            prepareGroup(
                group.first,
                blockExpression,
                std::move(blockFacts));
            if (error()) return zero();
          }
      }

    // Selection for every source/weight group is complete before any plan
    // executes, so execution cannot introduce a hidden policy decision.
    VECTOR(SymmetricTerm) resultTerms = std::move(scalarTerms);
    for (const auto& prepared : preparedConversions)
      {
        ring_elem converted =
            executeBasisConversionPlans(
                prepared.selection,
                prepared.expression,
                combinatorialTags,
                prepared.facts);
        if (error()) return zero();
        for (const auto& convertedTerm : polyValue(converted)->terms)
          appendTermIfNonZero(
              resultTerms,
              convertedTerm.coeff,
              convertedTerm.monomial);
      }

    // Source and degree groups can overlap in the target. Merge them once,
    // after every independently selected composition has executed.
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

    auto isSchurCompatible = [](BasisKind kind) {
      return kind == BasisKind::Schur ||
             kind == BasisKind::Complete ||
             kind == BasisKind::Elementary ||
             kind == BasisKind::PowerSum;
    };
    if (targetKind == BasisKind::Schur && leftKind && rightKind &&
        isSchurCompatible(*leftKind) &&
        isSchurCompatible(*rightKind))
      {
        std::string_view identifier =
            "Schur-product:littlewood-richardson";
        int specializedKinds =
            static_cast<int>(*leftKind == BasisKind::Complete ||
                             *rightKind == BasisKind::Complete) +
            static_cast<int>(*leftKind == BasisKind::Elementary ||
                             *rightKind == BasisKind::Elementary) +
            static_cast<int>(*leftKind == BasisKind::PowerSum ||
                             *rightKind == BasisKind::PowerSum);
        if (specializedKinds > 1)
          identifier = "Schur-product:compatible-factor-rules";
        else if (*leftKind == BasisKind::Complete ||
                 *rightKind == BasisKind::Complete)
          identifier = "Schur-product:horizontal-pieri";
        else if (*leftKind == BasisKind::Elementary ||
                 *rightKind == BasisKind::Elementary)
          identifier = "Schur-product:vertical-pieri";
        else if (*leftKind == BasisKind::PowerSum ||
                 *rightKind == BasisKind::PowerSum)
          identifier = "Schur-product:border-strips";
        plans.push_back({
            identifier,
            MultiplicationKernel::SchurCompatibleFactors,
            -1,
            -1,
            targetBasisId,
            true});
      }

    auto isMonomialLikeCompatible =
        [&](BasisKind kind) {
          return kind == targetKind ||
                 kind == BasisKind::Complete ||
                 kind == BasisKind::Elementary ||
                 kind == BasisKind::PowerSum;
        };
    if ((targetKind == BasisKind::Monomial ||
         targetKind == BasisKind::Forgotten) &&
        leftKind && rightKind &&
        isMonomialLikeCompatible(*leftKind) &&
        isMonomialLikeCompatible(*rightKind))
      plans.push_back({
          "monomial-like-product:exponent-splittings",
          MultiplicationKernel::MonomialLikeExpansion,
          -1,
          -1,
          targetBasisId,
          true});

    bool hallTarget =
        targetKind == BasisKind::HallLittlewoodQ ||
        targetKind == BasisKind::HallLittlewoodB ||
        targetKind == BasisKind::HallLittlewoodP ||
        targetKind == BasisKind::HallLittlewoodPOmega;
    if (hallTarget && leftKind && rightKind &&
        *leftKind != BasisKind::Custom &&
        *rightKind != BasisKind::Custom)
      {
        BasisKind generatorKind =
            targetKind == BasisKind::HallLittlewoodQ ||
                    targetKind == BasisKind::HallLittlewoodP
                ? BasisKind::HallLittlewoodQGenerator
                : BasisKind::HallLittlewoodBGenerator;
        int generatorId = requiredBasisIdForKind(generatorKind);
        if (error()) return plans;
        plans.push_back({
            "Hall-Littlewood-product:generators",
            MultiplicationKernel::CanonicalBasisProduct,
            generatorId,
            generatorId,
            generatorId,
            true});
      }

    if (isMultiplicativeBasis(targetBasisId))
      plans.push_back({
          "product:multiplicative-target",
          MultiplicationKernel::CanonicalBasisProduct,
          targetBasisId,
          targetBasisId,
          targetBasisId,
          true});
    plans.push_back({
        "product:power-sums",
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
    if (plan.kernel == MultiplicationKernel::SchurCompatibleFactors)
      {
        if (!leftFacts.expandedBasis || !rightFacts.expandedBasis)
          return false;
        BasisKind leftKind = basisKindForId(*leftFacts.expandedBasis);
        BasisKind rightKind = basisKindForId(*rightFacts.expandedBasis);
        BasisKind targetKind = basisKindForId(targetBasisId);
        auto compatible = [](BasisKind kind) {
          return kind == BasisKind::Schur ||
                 kind == BasisKind::Complete ||
                 kind == BasisKind::Elementary ||
                 kind == BasisKind::PowerSum;
        };
        return plan.kernelOutputBasisId == targetBasisId &&
               targetKind == BasisKind::Schur &&
               compatible(leftKind) && compatible(rightKind) &&
               leftFacts.skewFree && rightFacts.skewFree;
      }
    if (plan.kernel == MultiplicationKernel::MonomialLikeExpansion)
      {
        if (!leftFacts.expandedBasis || !rightFacts.expandedBasis)
          return false;
        BasisKind leftKind = basisKindForId(*leftFacts.expandedBasis);
        BasisKind rightKind = basisKindForId(*rightFacts.expandedBasis);
        BasisKind targetKind = basisKindForId(targetBasisId);
        auto compatible = [&](BasisKind kind) {
          return kind == targetKind ||
                 kind == BasisKind::Complete ||
                 kind == BasisKind::Elementary ||
                 kind == BasisKind::PowerSum;
        };
        return plan.kernelOutputBasisId == targetBasisId &&
               (targetKind == BasisKind::Monomial ||
                targetKind == BasisKind::Forgotten) &&
               compatible(leftKind) && compatible(rightKind) &&
               leftFacts.skewFree && rightFacts.skewFree;
      }
    if (plan.kernel != MultiplicationKernel::CanonicalBasisProduct ||
        !plan.kernelOutputCanonical ||
        plan.leftKernelBasisId <= 0 ||
        plan.rightKernelBasisId <= 0)
      return false;
    return basisKindForId(plan.leftKernelBasisId) != BasisKind::Custom &&
           basisKindForId(plan.rightKernelBasisId) != BasisKind::Custom;
  }

// ============================================================================
// Multiplication-Plan Picker
// ============================================================================
// Operand basis compositions are selected here with the product kernel.
// Support-dependent conversion edges are finalized from exact facts at their
// execution stage boundaries by the owning composition workflow.

SymmetricEngineRing::MultiplicationPlanSelection
SymmetricEngineRing::selectMultiplicationPlan(
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
                pickBasisConversionPlans(
                    *leftFacts.expandedBasis,
                    plan.leftKernelBasisId,
                    leftFacts,
                    0);
          if (error()) return selection;
          if (plan.rightKernelBasisId > 0 &&
              *rightFacts.expandedBasis != plan.rightKernelBasisId)
            selection.rightConversion =
                pickBasisConversionPlans(
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
          if (plan.identifier == forced &&
              multiplicationPlanApplicable(
                  plan, leftFacts, rightFacts, targetBasisId))
            return completeSelection(plan);
        ERROR("the requested multiplication plan is unknown or inapplicable: ",
              forced);
        int powerSumId = requiredBasisIdForKind(BasisKind::PowerSum);
        return {{"unavailable",
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
    return {{"unavailable",
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
// The try-style shared helpers predate this workflow. A false result here is a
// contract violation, not a signal to select a fallback plan.

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
// selected kernel, normalize its declared output, select the post-kernel
// composition, execute it, and validate the public target contract.

ring_elem SymmetricEngineRing::executeMultiplicationWorkflow(
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
              "stage=convert-left plans=%zu\n",
              selection.leftConversion->plans.size());
        left = executeBasisConversionPlans(
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
              "stage=convert-right plans=%zu\n",
              selection.rightConversion->plans.size());
        right = executeBasisConversionPlans(
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
          plan.identifier.data(),
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
        enrichConversionSelectionFacts(
            canonicalKernelResult, targetBasisId, kernelFacts);
        BasisConversionPlanSelection postKernelSelection =
            pickBasisConversionPlans(
                plan.kernelOutputBasisId,
                targetBasisId,
                kernelFacts,
                tags);
        if (error()) return zero();
        if (trace)
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-stage: "
              "stage=post-kernel-conversion plans=%zu\n",
              postKernelSelection.plans.size());
        result = executeBasisConversionPlans(
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
            leftFacts, rightFacts, targetBasisId);
    if (error()) return zero();
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr)
      std::fprintf(
          stderr,
          "SymmetricRings multiplication-plan: target=%s plan=%s\n",
          displayForBasis(targetBasisId).c_str(),
          selection.plan.identifier.data());
    return executeMultiplicationWorkflow(
        selection,
        f,
        g,
        leftFacts,
        rightFacts,
        targetBasisId,
        attachResultFacts);
}

ring_elem SymmetricEngineRing::multiplyBasisElementsToBasis(
    ring_elem f,
    ring_elem g,
    int targetBasisId) const
{
    auto strictFacts =
        [&](ring_elem operand) {
          const auto *poly = polyValue(operand);
          if (poly->terms.size() != 1 ||
              poly->terms.front().monomial.data.empty() ||
              atomLengthAt(poly->terms.front().monomial, 0) !=
                  poly->terms.front().monomial.data.size())
            return ExpressionFacts{};
          return basisElementFactsFromAtom(
              poly->terms.front().monomial,
              0,
              poly->terms.front().coeff);
        };
    ExpressionFacts leftFacts = strictFacts(f);
    ExpressionFacts rightFacts = strictFacts(g);
    if (!leftFacts.singleBasisElement() ||
        !rightFacts.singleBasisElement() ||
        !leftFacts.singleBasisElementId ||
        !rightFacts.singleBasisElementId ||
        !leftFacts.singleBasisElementCoefficientOne.value_or(false) ||
        !rightFacts.singleBasisElementCoefficientOne.value_or(false) ||
        !leftFacts.canonicalExpansionInBasis(
            *leftFacts.singleBasisElementId) ||
        !rightFacts.canonicalExpansionInBasis(
            *rightFacts.singleBasisElementId))
      {
        ERROR("the binary multiplication workflow expects exactly two "
              "coefficient-one normalized non-skew basis elements");
        return zero();
      }

    return multiplyBasisElementsWithFacts(
        f,
        g,
        leftFacts,
        rightFacts,
        targetBasisId,
        true);
  }

// ============================================================================
// One-Term Product Resolver
// ============================================================================
// Canonical storage has already removed zero coefficients. This helper removes
// identity factors, applies the surviving coefficient once at the end, and
// uses the multiplicative-target or pairwise branch prescribed by the design.

ring_elem SymmetricEngineRing::multiplyTermToBasis(
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
          ExpressionFacts enrichedFacts;
          const ExpressionFacts *selectionFacts = &factor.facts;
          if (basisKindForId(sourceBasisId) == BasisKind::PowerSum)
            {
              enrichedFacts = factor.facts;
              enrichConversionSelectionFacts(
                  factor.expression, targetBasisId, enrichedFacts);
              selectionFacts = &enrichedFacts;
            }
          BasisConversionPlanSelection selection =
              pickBasisConversionPlans(
                  sourceBasisId,
                  targetBasisId,
                  *selectionFacts,
                  factorTags);
          if (error()) return zero();
          return executeBasisConversionPlans(
              selection,
              factor.expression,
              factorTags,
              *selectionFacts);
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
            if (error()) return zero();
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
                    if (error()) return zero();
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
            factors.size() > 2);
        for (size_t i = 2; i < factors.size() && !error(); ++i)
          result = multiplyToBasis(
              result, factors[i].expression, targetBasisId);
      }
    if (error()) return zero();
    // Coefficients are immutable ring elements. Avoid walking and copying a
    // potentially large expansion when the product term already has unit
    // coefficient.
    if (coefficientRing->is_equal(
            term.coeff, coefficientRing->one()))
      return result;
    return scaled(term.coeff, result);
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
                *strictLeftFacts, *strictRightFacts, targetBasisId);
        if (error()) return zero();
        return executeMultiplicationWorkflow(
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
                leftFacts, rightFacts, targetBasisId);
        if (error()) return zero();
        return executeMultiplicationWorkflow(
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
          ExpressionFacts enrichedFacts;
          const ExpressionFacts *selectionFacts =
              &prepared.facts;
          if (basisKindForId(sourceBasisId) == BasisKind::PowerSum)
            {
              enrichedFacts = prepared.facts;
              enrichConversionSelectionFacts(
                  prepared.factor, targetBasisId, enrichedFacts);
              selectionFacts = &enrichedFacts;
            }
          BasisConversionPlanSelection conversion =
              pickBasisConversionPlans(
                  sourceBasisId,
                  targetBasisId,
                  *selectionFacts,
                  tags);
          if (error()) return zero();
          return executeBasisConversionPlans(
              conversion,
              prepared.factor,
              tags,
              *selectionFacts);
        };

    // Structured kernels need no operand compositions, and their selection
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
                product = executeMultiplicationWorkflow(
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
                          leftTerm.facts,
                          rightTerm.facts,
                          targetBasisId);
                  if (error()) return zero();
                  if (!selection.leftConversion &&
                      !selection.rightConversion)
                    reusablePairSelections.emplace(key, selection);
                  product = executeMultiplicationWorkflow(
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
// grouping, composition selection, exact stage-edge selection and execution,
// and collection.

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
    VECTOR(SymmetricTerm) preparedTerms;
    size_t resolvedProductTerms = 0;
    bool passthroughTermsAlreadyInTarget = true;
    const auto& normalizedTerms = polyValue(normalized)->terms;
    if (normalizedFactorCounts.size() != normalizedTerms.size())
      {
        ERROR("normalized factor-count profile is inconsistent");
        return zero();
      }
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
        ring_elem converted =
            multiplyTermToBasis(term, targetBasisId);
        if (error()) return zero();
        for (const auto& convertedTerm : polyValue(converted)->terms)
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
        // it merely to execute an identity composition.
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
