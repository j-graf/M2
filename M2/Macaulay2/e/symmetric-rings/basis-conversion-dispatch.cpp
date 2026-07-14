// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "basic-rings/aring-glue.hpp"
#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"
#include "rings/ZZ.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace symmetric_rings {

namespace {
size_t partitionCount(int weight)
{
    if (weight < 0) return 0;
    std::vector<size_t> counts(static_cast<size_t>(weight) + 1, 0);
    counts[0] = 1;
    for (int part = 1; part <= weight; ++part)
      for (int total = part; total <= weight; ++total)
        {
          size_t addend = counts[static_cast<size_t>(total - part)];
          size_t& value = counts[static_cast<size_t>(total)];
          if (value > std::numeric_limits<size_t>::max() - addend)
            value = std::numeric_limits<size_t>::max();
          else
            value += addend;
        }
    return counts[static_cast<size_t>(weight)];
}

bool completeFriendlyPowerSumIndex(const Partition& index, int weight)
{
    if (weight <= 0 || index.size() < 4) return false;
    int smallCycleWeight = 0;
    for (int part : index)
      if (part <= 4) smallCycleWeight += part;
    // Repeated short cycles create broad rim-hook frontiers, but their
    // logarithmic h-expansions remain narrow.  Hybrid conversion has a fixed
    // split-and-merge cost, so only strongly short-cycle indices should enter
    // its complete subgroup; weaker two-thirds candidates were observed to
    // make scattered and middle-slice supports slower than whole-input
    // abacus conversion.
    return 4 * smallCycleWeight >= 3 * weight;
}

std::string combinatorialTagNames(CombinatorialTags tags)
{
    if (tags == 0) return "none";
    std::string result;
    auto append = [&](CombinatorialTag tag, const char *name) {
      if (!hasCombinatorialTag(tags, tag)) return;
      if (!result.empty()) result += ",";
      result += name;
    };
    append(CombinatorialTag::Plethysm, "plethysm");
    append(CombinatorialTag::LittlewoodRichardson,
           "littlewood-richardson");
    append(CombinatorialTag::HorizontalPieri, "horizontal-pieri");
    append(CombinatorialTag::VerticalPieri, "vertical-pieri");
    append(CombinatorialTag::BorderStrips, "border-strips");
    return result;
}
}

// ============================================================================
// Conversion Guarantees And Profiles
// ============================================================================
// Expression facts are inferred once and updated as pipeline stages change them.

SymmetricEngineRing::ConversionGuarantees
SymmetricEngineRing::inferConversionGuarantees(
        ring_elem f,
        int targetBasisId) const
{
    ConversionGuarantees guarantees;
    const auto *poly = polyValue(f);
    if (poly->conversionMetadata)
      {
        const auto& metadata = *poly->conversionMetadata;
        guarantees.pureBasis = metadata.pureBasis;
        guarantees.expandedBasis = metadata.expandedBasis;
        guarantees.homogeneousWeight = metadata.homogeneousWeight;
        guarantees.termCount = metadata.termCount;
        guarantees.maximumPartitionLength = metadata.maximumPartitionLength;
        guarantees.density = metadata.density;
        guarantees.factorBases = metadata.factorBases;
        if (metadata.singleBasisElement)
          guarantees.singleBasisElement = *metadata.singleBasisElement
              ? KnownState::True : KnownState::False;
        if (metadata.singleTerm)
          guarantees.singleTerm = *metadata.singleTerm
              ? KnownState::True : KnownState::False;
        if (metadata.noProducts)
          guarantees.noProducts = *metadata.noProducts
              ? KnownState::True : KnownState::False;
        if (metadata.normalized) guarantees.normalized = KnownState::True;
        if (metadata.skewFree) guarantees.skewFree = KnownState::True;
        if (metadata.collected) guarantees.collected = KnownState::True;
        return strengthenConversionGuarantees(std::move(guarantees), targetBasisId);
      }
    guarantees.termCount = poly->terms.size();

    bool noProducts = true;
    bool normalized = true;
    bool skewFree = true;
    int pureBasis = 0;
    bool mixedBases = false;
    size_t totalBasisElements = 0;
    std::optional<int> homogeneousWeight;
    bool weightsAgree = true;
    std::vector<int> factorBases;

    for (const auto& term : poly->terms)
      {
        size_t termBasisElements = 0;
        int termWeight = 0;
        size_t pos = 0;
        while (pos < term.monomial.data.size())
          {
            ++termBasisElements;
            ++totalBasisElements;
            int basisId = atomBasisIdAt(term.monomial, pos);
            if (std::find(factorBases.begin(), factorBases.end(), basisId) ==
                factorBases.end())
              factorBases.push_back(basisId);
            if (pureBasis == 0)
              pureBasis = basisId;
            else if (pureBasis != basisId)
              mixedBases = true;

            if (atomIsSkewAt(term.monomial, pos))
              {
                skewFree = false;
                Partition outer =
                    basisElementOuterIndex(term.monomial, pos);
                Partition inner =
                    basisElementInnerIndex(term.monomial, pos);
                termWeight += partitionWeight(outer) - partitionWeight(inner);
                normalized = normalized && isPartitionIndex(outer) &&
                             isPartitionIndex(inner);
              }
            else
              {
                Partition index = basisElementIndex(term.monomial, pos);
                termWeight += partitionWeight(index);
                if (!isMultiplicativeBasis(basisId))
                  normalized = normalized && isPartitionIndex(index);
              }
            pos += atomLengthAt(term.monomial, pos);
          }
        if (termBasisElements > 1) noProducts = false;
        if (!homogeneousWeight)
          homogeneousWeight = termWeight;
        else if (*homogeneousWeight != termWeight)
          weightsAgree = false;
      }

    guarantees.singleTerm = poly->terms.size() == 1
        ? KnownState::True
        : KnownState::False;
    guarantees.singleBasisElement =
        poly->terms.size() == 1 && totalBasisElements == 1
        ? KnownState::True
        : KnownState::False;
    guarantees.noProducts = noProducts ? KnownState::True : KnownState::False;
    guarantees.normalized = normalized ? KnownState::True : KnownState::False;
    guarantees.skewFree = skewFree ? KnownState::True : KnownState::False;
    guarantees.collected = KnownState::True;
    if (weightsAgree && homogeneousWeight)
      guarantees.homogeneousWeight = homogeneousWeight;
    std::sort(factorBases.begin(), factorBases.end());
    guarantees.factorBases = std::move(factorBases);

    if (!mixedBases && pureBasis != 0) guarantees.pureBasis = pureBasis;
    if (guarantees.pureBasis && noProducts && skewFree)
      guarantees.expandedBasis = guarantees.pureBasis;

    if (poly->terms.empty() || totalBasisElements == 0)
      guarantees.targetClosed = KnownState::True;
    else if (guarantees.expandedBasis == targetBasisId && normalized)
      guarantees.targetClosed = KnownState::True;
    else
      guarantees.targetClosed = KnownState::False;

    return strengthenConversionGuarantees(std::move(guarantees), targetBasisId);
  }

SymmetricEngineRing::ConversionGuarantees
SymmetricEngineRing::strengthenConversionGuarantees(
        ConversionGuarantees guarantees,
        int targetBasisId) const
{
    auto setIfUnknown = [](KnownState& state, KnownState value) {
      if (state == KnownState::Unknown) state = value;
    };

    if (guarantees.singleBasisElement == KnownState::True)
      {
        setIfUnknown(guarantees.singleTerm, KnownState::True);
        setIfUnknown(guarantees.noProducts, KnownState::True);
        if (!guarantees.termCount) guarantees.termCount = 1;
      }

    if (guarantees.expandedBasis)
      {
        if (!guarantees.pureBasis)
          guarantees.pureBasis = guarantees.expandedBasis;
        setIfUnknown(guarantees.noProducts, KnownState::True);
        setIfUnknown(guarantees.skewFree, KnownState::True);
      }

    if (guarantees.pureBasis && !guarantees.factorBases)
      guarantees.factorBases = std::vector<int>{*guarantees.pureBasis};

    if (guarantees.pureBasis &&
        isMultiplicativeBasis(*guarantees.pureBasis) &&
        guarantees.skewFree == KnownState::True)
      {
        if (!guarantees.expandedBasis)
          guarantees.expandedBasis = guarantees.pureBasis;
        setIfUnknown(guarantees.noProducts, KnownState::True);
      }

    if (guarantees.expandedBasis == targetBasisId &&
        guarantees.normalized == KnownState::True)
      setIfUnknown(guarantees.targetClosed, KnownState::True);

    return guarantees;
  }

SymmetricEngineRing::ConversionGuarantees
SymmetricEngineRing::ensureConversionProfile(
        ring_elem f,
        ConversionGuarantees guarantees,
        ConversionProfileFact fact) const
{
    if (fact == ConversionProfileFact::Density)
      {
        if (!guarantees.density && guarantees.homogeneousWeight &&
            guarantees.termCount)
          {
            size_t possible = partitionCount(*guarantees.homogeneousWeight);
            if (possible != 0)
              guarantees.density = static_cast<double>(*guarantees.termCount) /
                                   static_cast<double>(possible);
          }
        return guarantees;
      }

    if (guarantees.maximumPartitionLength) return guarantees;
    size_t maximum = 0;
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        size_t pos = 0;
        while (pos < term.monomial.data.size())
          {
            if (atomIsSkewAt(term.monomial, pos))
              maximum = std::max(
                  maximum,
                  static_cast<size_t>(std::max(
                      partitionLength(
                          basisElementOuterIndex(term.monomial, pos)),
                      partitionLength(
                          basisElementInnerIndex(term.monomial, pos)))));
            else
              maximum = std::max(
                  maximum,
                  static_cast<size_t>(partitionLength(
                      basisElementIndex(term.monomial, pos))));
            pos += atomLengthAt(term.monomial, pos);
          }
      }
    guarantees.maximumPartitionLength = maximum;
    return guarantees;
  }

SymmetricEngineRing::ConversionGuarantees
SymmetricEngineRing::guaranteesForFactorizedProduct(
        ring_elem f,
        ring_elem g,
        int targetBasisId) const
{
    ConversionGuarantees left = inferConversionGuarantees(f, targetBasisId);
    ConversionGuarantees right = inferConversionGuarantees(g, targetBasisId);
    ConversionGuarantees result;
    result.factorizedProduct = KnownState::True;

    if (left.pureBasis && right.pureBasis &&
        left.pureBasis == right.pureBasis)
      result.pureBasis = left.pureBasis;
    if (left.factorBases && right.factorBases)
      {
        std::vector<int> factors = *left.factorBases;
        factors.insert(factors.end(),
                       right.factorBases->begin(),
                       right.factorBases->end());
        std::sort(factors.begin(), factors.end());
        factors.erase(std::unique(factors.begin(), factors.end()), factors.end());
        result.factorBases = std::move(factors);
      }
    if (left.homogeneousWeight && right.homogeneousWeight)
      result.homogeneousWeight =
          *left.homogeneousWeight + *right.homogeneousWeight;
    if (left.termCount && right.termCount)
      result.termCount = *left.termCount * *right.termCount;
    if (left.maximumPartitionLength && right.maximumPartitionLength)
      result.maximumPartitionLength = std::max(
          *left.maximumPartitionLength, *right.maximumPartitionLength);

    auto conjunction = [](KnownState a, KnownState b) {
      if (a == KnownState::False || b == KnownState::False)
        return KnownState::False;
      if (a == KnownState::True && b == KnownState::True)
        return KnownState::True;
      return KnownState::Unknown;
    };
    result.singleTerm = conjunction(left.singleTerm, right.singleTerm);
    result.normalized = conjunction(left.normalized, right.normalized);
    result.skewFree = conjunction(left.skewFree, right.skewFree);
    result.collected = conjunction(left.collected, right.collected);
    result.targetClosed = KnownState::Unknown;
    return strengthenConversionGuarantees(std::move(result), targetBasisId);
  }

SymmetricEngineRing::ConversionGuarantees
SymmetricEngineRing::guaranteesAfterNormalization(
        ring_elem f,
        int targetBasisId) const
{
    ConversionGuarantees result = inferConversionGuarantees(f, targetBasisId);
    result.normalized = KnownState::True;
    result.collected = KnownState::True;
    result.factorizedProduct = KnownState::False;
    return strengthenConversionGuarantees(std::move(result), targetBasisId);
  }

SymmetricEngineRing::ConversionGuarantees
SymmetricEngineRing::guaranteesAfterProductExpansion(
        ring_elem f,
        int targetBasisId) const
{
    ConversionGuarantees result = inferConversionGuarantees(f, targetBasisId);
    result.factorizedProduct = KnownState::False;
    return strengthenConversionGuarantees(std::move(result), targetBasisId);
  }

SymmetricEngineRing::ConversionGuarantees
SymmetricEngineRing::guaranteesAfterSourceTargetConversion(
        const ConversionGuarantees& inputGuarantees,
        ring_elem result,
        int targetBasisId) const
{
    ConversionGuarantees guarantees;
    const auto *poly = polyValue(result);
    guarantees.pureBasis = targetBasisId;
    guarantees.expandedBasis = targetBasisId;
    guarantees.factorBases = std::vector<int>{targetBasisId};
    guarantees.homogeneousWeight = inputGuarantees.homogeneousWeight;
    guarantees.termCount = poly->terms.size();
    guarantees.singleTerm = poly->terms.size() == 1
        ? KnownState::True : KnownState::False;
    guarantees.singleBasisElement = poly->terms.size() == 1 &&
                            !poly->terms[0].monomial.data.empty()
        ? KnownState::True : KnownState::False;
    guarantees.noProducts = KnownState::True;
    guarantees.normalized = KnownState::True;
    guarantees.skewFree = KnownState::True;
    guarantees.collected = KnownState::True;
    guarantees.targetClosed = KnownState::True;
    guarantees.factorizedProduct = KnownState::False;
    return strengthenConversionGuarantees(std::move(guarantees), targetBasisId);
  }

SymmetricEngineRing::ConversionGuarantees
SymmetricEngineRing::guaranteesAfterAddition(
        const ConversionGuarantees& inputGuarantees,
        ring_elem result,
        int targetBasisId) const
{
    return guaranteesAfterSourceTargetConversion(
        inputGuarantees, result, targetBasisId);
  }

void SymmetricEngineRing::attachConversionGuarantees(
        ring_elem f,
        const ConversionGuarantees& guarantees,
        CombinatorialTags combinatorialTags) const
{
    SymmetricConversionMetadata metadata;
    metadata.pureBasis = guarantees.pureBasis;
    metadata.expandedBasis = guarantees.expandedBasis;
    metadata.homogeneousWeight = guarantees.homogeneousWeight;
    metadata.termCount = guarantees.termCount;
    metadata.maximumPartitionLength = guarantees.maximumPartitionLength;
    metadata.density = guarantees.density;
    metadata.factorBases = guarantees.factorBases;
    if (guarantees.singleBasisElement != KnownState::Unknown)
      metadata.singleBasisElement =
          guarantees.singleBasisElement == KnownState::True;
    if (guarantees.singleTerm != KnownState::Unknown)
      metadata.singleTerm = guarantees.singleTerm == KnownState::True;
    if (guarantees.noProducts != KnownState::Unknown)
      metadata.noProducts = guarantees.noProducts == KnownState::True;
    metadata.normalized = guarantees.normalized == KnownState::True;
    metadata.skewFree = guarantees.skewFree == KnownState::True;
    metadata.collected = guarantees.collected == KnownState::True;
    auto *poly = mutablePolyValue(f);
    poly->combinatorialTags = combinatorialTags;
    poly->conversionMetadata = std::move(metadata);
  }

// ============================================================================
// Top-Level Conversion Pipeline
// ============================================================================
// The selector chooses a workflow; execution delegates to one named pipeline.

SymmetricEngineRing::ConversionPipeline SymmetricEngineRing::selectConversionPipeline(
        const ConversionRequest& request,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        bool targetIsMultiplicative) const
{
    if (request.kind == ConversionRequestKind::FactorizedProduct)
      return ConversionPipeline::FactorizedProduct;
    if (request.kind == ConversionRequestKind::PostPlethysm)
      return ConversionPipeline::PostPlethysm;

    const ConversionInput& input = request.input;
    if (input.guarantees.expandedBasis == pBasisId)
      return ConversionPipeline::PowerSums;

    if (targetIsMultiplicative &&
        input.guarantees.normalized == KnownState::True &&
        input.guarantees.skewFree == KnownState::True)
      return ConversionPipeline::GroupedMultiplicativeTarget;

    if (canUseGroupedHallLittlewoodPipeline(input.guarantees,
                                            targetBasisId,
                                            targetDisplay))
      return ConversionPipeline::GroupedHallLittlewood;

    if (canUseWholeExpressionPipeline(input.guarantees,
                                      targetBasisId,
                                      targetDisplay,
                                      targetIsMultiplicative))
      return ConversionPipeline::WholeExpression;

    return ConversionPipeline::FallbackTerm;
  }

const char *SymmetricEngineRing::conversionPipelineName(
        ConversionPipeline pipeline) const
{
    switch (pipeline)
      {
        case ConversionPipeline::WholeExpression: return "whole-expression";
        case ConversionPipeline::GroupedMultiplicativeTarget:
          return "grouped-multiplicative-target";
        case ConversionPipeline::GroupedHallLittlewood:
          return "grouped-hall-littlewood";
        case ConversionPipeline::PowerSums: return "power-sums";
        case ConversionPipeline::PostPlethysmPowerSums:
          return "post-plethysm-power-sums";
        case ConversionPipeline::FallbackTerm: return "fallback-term";
        case ConversionPipeline::FactorizedProduct: return "factorized-product";
        case ConversionPipeline::PostPlethysm: return "post-plethysm";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceConversionSelection(
        ConversionPipeline pipeline,
        const ConversionRequest& request,
        int targetBasisId,
        const std::string& targetDisplay) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    const ConversionInput& input = request.input;

    auto knownInt = [](const std::optional<int>& value) {
      return value ? std::to_string(*value) : std::string("unknown");
    };
    auto knownSize = [](const std::optional<size_t>& value) {
      return value ? std::to_string(*value) : std::string("unknown");
    };
    auto knownDouble = [](const std::optional<double>& value) {
      return value ? std::to_string(*value) : std::string("unknown");
    };
    const char *route = "pipeline-defined";
    if ((pipeline == ConversionPipeline::PowerSums ||
         pipeline == ConversionPipeline::PostPlethysmPowerSums))
      {
        int pBasisId = basisIdForKind(BasisKind::PowerSum);
        if (pBasisId >= 0 && targetBasisId >= 0)
          route = powerSumsToTargetRouteName(
              selectPowerSumsToTargetRoute(
                  input, pBasisId, targetBasisId, targetDisplay));
      }

    std::fprintf(stderr,
                 "SymmetricRings conversion: pipeline=%s route=%s target=%s "
                 "pureBasis=%s expandedBasis=%s terms=%s weight=%s "
                 "maxLength=%s density=%s tags=%s\n",
                 conversionPipelineName(pipeline),
                 route,
                 targetDisplay.c_str(),
                 knownInt(input.guarantees.pureBasis).c_str(),
                 knownInt(input.guarantees.expandedBasis).c_str(),
                 knownSize(input.guarantees.termCount).c_str(),
                 knownInt(input.guarantees.homogeneousWeight).c_str(),
                 knownSize(input.guarantees.maximumPartitionLength).c_str(),
                 knownDouble(input.guarantees.density).c_str(),
                 combinatorialTagNames(input.combinatorialTags).c_str());
  }

ring_elem SymmetricEngineRing::executeConversionPipeline(
        ConversionPipeline pipeline,
        const ConversionRequest& request,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    const ConversionInput& input = request.input;
    switch (pipeline)
      {
        case ConversionPipeline::WholeExpression:
          return runWholeExpressionPipeline(input,
                                         pBasisId,
                                         pDisplay,
                                         pOrder,
                                         pIsMultiplicative,
                                         targetBasisId,
                                         targetDisplay,
                                         targetOrder,
                                         targetIsMultiplicative);
        case ConversionPipeline::GroupedMultiplicativeTarget:
          return runGroupedMultiplicativeTargetPipeline(input,
                                         pBasisId,
                                         pDisplay,
                                         pOrder,
                                         pIsMultiplicative,
                                         targetBasisId,
                                         targetDisplay,
                                         targetOrder,
                                         targetIsMultiplicative);
        case ConversionPipeline::GroupedHallLittlewood:
          return runGroupedHallLittlewoodPipeline(input,
                                         pBasisId,
                                         pDisplay,
                                         pOrder,
                                         pIsMultiplicative,
                                         targetBasisId,
                                         targetDisplay,
                                         targetOrder,
                                         targetIsMultiplicative);
        case ConversionPipeline::PowerSums:
          return runPowerSumsPipeline(input,
                                  pBasisId,
                                  pDisplay,
                                  pOrder,
                                  pIsMultiplicative,
                                  targetBasisId,
                                  targetDisplay,
                                  targetOrder,
                                  targetIsMultiplicative);
        case ConversionPipeline::PostPlethysmPowerSums:
          return runPostPlethysmPowerSumsPipeline(input,
                                  pBasisId,
                                  pDisplay,
                                  pOrder,
                                  pIsMultiplicative,
                                  targetBasisId,
                                  targetDisplay,
                                  targetOrder,
                                  targetIsMultiplicative);
        case ConversionPipeline::FallbackTerm:
          return runFallbackTermPipeline(input,
                                      pBasisId,
                                      pDisplay,
                                      pOrder,
                                      pIsMultiplicative,
                                      targetBasisId,
                                      targetDisplay,
                                      targetOrder,
                                      targetIsMultiplicative);
        case ConversionPipeline::FactorizedProduct:
          if (!request.leftOperand || !request.rightOperand)
            {
              ERROR("factorized-product conversion requires two operands");
              return zero();
            }
          return runFactorizedProductPipeline(input,
                                           *request.leftOperand,
                                           *request.rightOperand,
                                           pBasisId,
                                           pDisplay,
                                           pOrder,
                                           pIsMultiplicative,
                                           targetBasisId,
                                           targetDisplay,
                                           targetOrder,
                                           targetIsMultiplicative);
        case ConversionPipeline::PostPlethysm:
          if (!request.leftOperand || !request.rightOperand)
            {
              ERROR("post-plethysm conversion requires two operands");
              return zero();
            }
          return runPostPlethysmPipeline(input,
                                         *request.leftOperand,
                                         *request.rightOperand,
                                         pBasisId,
                                         pDisplay,
                                         pOrder,
                                         pIsMultiplicative,
                                         targetBasisId,
                                         targetDisplay,
                                         targetOrder,
                                         targetIsMultiplicative);
      }

    ERROR("unknown basis conversion pipeline");
    return zero();
  }

// ============================================================================
// Whole-Expression Pipeline
// ============================================================================
// Whole-expression routes avoid termwise work when guarantees permit it.

bool SymmetricEngineRing::canUseWholeExpressionPipeline(
        const ConversionGuarantees& guarantees,
        int targetBasisId,
        const std::string& targetDisplay,
        bool targetIsMultiplicative) const
{
    if (guarantees.targetClosed == KnownState::True) return true;

    if (guarantees.expandedBasis &&
        guarantees.normalized == KnownState::True &&
        guarantees.skewFree == KnownState::True)
      return true;

    if (!guarantees.factorBases) return false;

    auto allFactorKindsAre = [&](const std::vector<BasisKind>& kinds) {
      for (int basisId : *guarantees.factorBases)
        {
          BasisKind kind = basisKindForId(basisId);
          if (std::find(kinds.begin(), kinds.end(), kind) == kinds.end())
            return false;
        }
      return true;
    };

    bool ordinaryNormalized =
        guarantees.skewFree == KnownState::True &&
        guarantees.normalized == KnownState::True;

    const BasisKind targetKind = basisKindForId(targetBasisId);
    if (targetKind == BasisKind::Schur)
      {
        if (guarantees.pureBasis == targetBasisId) return true;
        return ordinaryNormalized &&
               allFactorKindsAre({BasisKind::Schur, BasisKind::Complete,
                                  BasisKind::Elementary, BasisKind::PowerSum});
      }

    if (targetKind == BasisKind::SchurOmega)
      return ordinaryNormalized && allFactorKindsAre({BasisKind::Elementary});

    if (targetKind == BasisKind::HallLittlewoodQ ||
        targetKind == BasisKind::HallLittlewoodP ||
        targetKind == BasisKind::HallLittlewoodB ||
        targetKind == BasisKind::HallLittlewoodPOmega)
      return false;

    if (targetKind != BasisKind::Complete &&
        targetKind != BasisKind::Elementary &&
        !targetIsMultiplicative)
      return false;

    if (!ordinaryNormalized) return false;
    return allFactorKindsAre({BasisKind::PowerSum, BasisKind::Complete,
                              BasisKind::Elementary,
                              BasisKind::HallLittlewoodQGenerator,
                              BasisKind::HallLittlewoodBGenerator,
                              BasisKind::Monomial,
                              BasisKind::Forgotten, BasisKind::Schur,
                              BasisKind::SchurOmega, BasisKind::HallLittlewoodQ,
                              BasisKind::HallLittlewoodB,
                              BasisKind::HallLittlewoodP,
                              BasisKind::HallLittlewoodPOmega});
  }

SymmetricEngineRing::WholeExpressionRoute
SymmetricEngineRing::selectWholeExpressionRoute(
        const ConversionInput& input,
        int targetBasisId,
        const std::string& targetDisplay) const
{
    if (input.guarantees.targetClosed == KnownState::True)
      return WholeExpressionRoute::AlreadyInTarget;
    if (input.guarantees.expandedBasis &&
        input.guarantees.normalized == KnownState::True &&
        input.guarantees.skewFree == KnownState::True)
      {
        const BasisKind sourceKind =
            basisKindForId(*input.guarantees.expandedBasis);
        const BasisKind targetKind = basisKindForId(targetBasisId);
        if ((sourceKind == BasisKind::HallLittlewoodQ && targetKind == BasisKind::HallLittlewoodP) ||
            (sourceKind == BasisKind::HallLittlewoodP && targetKind == BasisKind::HallLittlewoodQ) ||
            (sourceKind == BasisKind::HallLittlewoodB && targetKind == BasisKind::HallLittlewoodPOmega) ||
            (sourceKind == BasisKind::HallLittlewoodPOmega && targetKind == BasisKind::HallLittlewoodB))
          return WholeExpressionRoute::ViaHallLittlewoodNormalization;
        if ((sourceKind == BasisKind::Schur && targetKind == BasisKind::SchurOmega) ||
            (sourceKind == BasisKind::SchurOmega && targetKind == BasisKind::Schur))
          return WholeExpressionRoute::ViaSchurOmegaConjugation;
      }
    switch (basisKindForId(targetBasisId))
      {
      case BasisKind::Schur:
        {
        int completeId = requiredBasisIdForKind(BasisKind::Complete);
        if (error()) return WholeExpressionRoute::NoApplicableRoute;
        if (input.guarantees.pureBasis == completeId)
          return WholeExpressionRoute::ViaCompleteToSchurRecursiveTransition;
        return WholeExpressionRoute::ViaSchurCompatibleProducts;
        }
      case BasisKind::SchurOmega:
        return WholeExpressionRoute::ViaSchurTriangularReduction;
      case BasisKind::HallLittlewoodQ:
      case BasisKind::HallLittlewoodB:
      case BasisKind::HallLittlewoodP:
      case BasisKind::HallLittlewoodPOmega:
        return WholeExpressionRoute::ViaHallLittlewoodTriangularReduction;
      default:
        return WholeExpressionRoute::NoApplicableRoute;
      }
  }

const char *SymmetricEngineRing::wholeExpressionRouteName(
        WholeExpressionRoute route) const
{
    switch (route)
      {
        case WholeExpressionRoute::AlreadyInTarget:
          return "already-in-target";
        case WholeExpressionRoute::ViaCompleteToSchurRecursiveTransition:
          return "h->S:recursive-transition";
        case WholeExpressionRoute::ViaSchurCompatibleProducts:
          return "Schur-compatible-products";
        case WholeExpressionRoute::ViaSchurTriangularReduction:
          return "Schur-triangular-reduction";
        case WholeExpressionRoute::ViaHallLittlewoodNormalization:
          return "Hall-Littlewood-capital-normalized:diagonal-scaling";
        case WholeExpressionRoute::ViaSchurOmegaConjugation:
          return "Schur-omega:partition-conjugation";
        case WholeExpressionRoute::ViaHallLittlewoodTriangularReduction:
          return "Hall-Littlewood-triangular-reduction";
        case WholeExpressionRoute::NoApplicableRoute:
          return "source-target-fallback";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceWholeExpressionSelection(
        WholeExpressionRoute route,
        const ConversionInput& input,
        const std::string& targetDisplay) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::fprintf(stderr,
                 "SymmetricRings whole-expression: target=%s route=%s "
                 "terms=%zu\n",
                 targetDisplay.c_str(),
                 wholeExpressionRouteName(route),
                 input.guarantees.termCount.value_or(
                     polyValue(input.expression)->terms.size()));
  }

bool SymmetricEngineRing::executeWholeExpressionRoute(
        WholeExpressionRoute route,
        const ConversionInput& input,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        ring_elem& result) const
{
    switch (route)
      {
        case WholeExpressionRoute::AlreadyInTarget:
          result = copyPolyValue(polyValue(input.expression));
          return true;
        case WholeExpressionRoute::ViaCompleteToSchurRecursiveTransition:
          {
            int completeId = requiredBasisIdForKind(BasisKind::Complete);
            if (error()) return false;
            result = completeToSchurViaRecursiveTransition(
                input.expression,
                completeId,
                displayForBasis(completeId),
                basisOrderForId(completeId),
                isMultiplicativeBasis(completeId),
                targetBasisId,
                targetDisplay,
                targetOrder);
            return !error();
          }
        case WholeExpressionRoute::ViaSchurCompatibleProducts:
          return trySchurCompatibleExpressionToSchur(
              input.expression, targetBasisId, targetDisplay, targetOrder, result);
        case WholeExpressionRoute::ViaSchurTriangularReduction:
          return tryExpressionToSchurViaTriangularReduction(
              input.expression, targetBasisId, targetDisplay, targetOrder, result);
        case WholeExpressionRoute::ViaHallLittlewoodNormalization:
          {
            int sourceBasisId = *input.guarantees.expandedBasis;
            BasisKind sourceKind = basisKindForId(sourceBasisId);
            result = hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
                input.expression,
                sourceBasisId,
                targetBasisId,
                targetDisplay,
                targetOrder,
                sourceKind == BasisKind::HallLittlewoodQ ||
                    sourceKind == BasisKind::HallLittlewoodB);
            return !error();
          }
        case WholeExpressionRoute::ViaSchurOmegaConjugation:
          result = schurOmegaConversionViaPartitionConjugation(
              input.expression,
              *input.guarantees.expandedBasis,
              targetBasisId,
              targetDisplay,
              targetOrder);
          return !error();
        case WholeExpressionRoute::ViaHallLittlewoodTriangularReduction:
          return tryExpressionToHallLittlewoodViaTriangularReduction(
              input.expression, targetBasisId, targetDisplay, targetOrder, result);
        case WholeExpressionRoute::NoApplicableRoute:
          return false;
      }
    return false;
  }

ring_elem SymmetricEngineRing::runWholeExpressionPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    WholeExpressionRoute route = selectWholeExpressionRoute(
        input, targetBasisId, targetDisplay);
    if (error()) return zero();
    traceWholeExpressionSelection(route, input, targetDisplay);
    ring_elem result;
    if (executeWholeExpressionRoute(route,
                                     input,
                                     targetBasisId,
                                     targetDisplay,
                                     targetOrder,
                                     result))
      {
        ConversionGuarantees guarantees =
            guaranteesAfterSourceTargetConversion(
                input.guarantees, result, targetBasisId);
        attachConversionGuarantees(
            result, guarantees, input.combinatorialTags);
        return result;
      }
    if (error()) return zero();
    return sourceToTargetDispatch(input,
                                  pBasisId,
                                  pDisplay,
                                  pOrder,
                                  pIsMultiplicative,
                                  targetBasisId,
                                  targetDisplay,
                                  targetOrder,
                                  targetIsMultiplicative);
  }

// ============================================================================
// Power-Sums-To-Target Dispatch
// ============================================================================
// All p-to-target route choices are visible together in this dispatcher family.

SymmetricEngineRing::PowerSumsToTargetRoute
SymmetricEngineRing::selectPowerSumsToSchurRoute(
        const ConversionInput& input,
        bool omega) const
{
      const char *forcedRoute = std::getenv(
          "M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE");
      if (forcedRoute != nullptr)
        {
          std::string route(forcedRoute);
          if (route == "border-strips")
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips
                : PowerSumsToTargetRoute::ViaSchurBorderStrips;
          if (route == "abacus-rim-hooks" || route == "abacus")
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusRimHooks
                : PowerSumsToTargetRoute::ViaSchurAbacusRimHooks;
          if (route == "via-complete")
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurComplete
                : PowerSumsToTargetRoute::ViaSchurComplete;
          if (route == "grouped-characters")
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurCharacters
                : PowerSumsToTargetRoute::ViaSchurCharacters;
        }

      const ConversionGuarantees& guarantees = input.guarantees;
      if (!guarantees.homogeneousWeight)
        return PowerSumsToTargetRoute::ViaSchurDegreeBlocks;

      int weight = *guarantees.homogeneousWeight;
      size_t termCount = guarantees.termCount.value_or(
          polyValue(input.expression)->terms.size());
      size_t possibleTerms = partitionCount(weight);
      double density = possibleTerms == 0
          ? 0.0
          : static_cast<double>(termCount) /
            static_cast<double>(possibleTerms);
      long double termsSquared = static_cast<long double>(termCount) *
                                 static_cast<long double>(termCount);
      bool supportSquareFavorsComplete =
          2.0L * termsSquared >=
          35.0L * static_cast<long double>(possibleTerms);

      // Plethysm outputs have a different support profile from ordinary
      // power-sum expressions. Benchmarks place their crossover at
      // termCount/sqrt(partitionCount(weight)) approximately sqrt(35/2).
      const CombinatorialTags borderStripsTag =
          combinatorialTagMask(CombinatorialTag::BorderStrips);
      if (hasCombinatorialTag(input.combinatorialTags,
                              CombinatorialTag::Plethysm))
        {
          // Schur plethysm outputs acquire enough correlated short-cycle
          // support for the complete route well before the generic density
          // crossover.  Tiny outputs remain on abacus, where setup dominates.
          bool completeWins = supportSquareFavorsComplete ||
                              (weight >= 8 && termCount >= 2);
          return completeWins
              ? (omega
                  ? PowerSumsToTargetRoute::ViaOmegaThenSchurComplete
                  : PowerSumsToTargetRoute::ViaSchurComplete)
              : (omega
                  ? PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusRimHooks
                  : PowerSumsToTargetRoute::ViaSchurAbacusRimHooks);
        }

      if (hasCombinatorialTag(input.combinatorialTags,
                              CombinatorialTag::LittlewoodRichardson))
        {
          bool completeWins = weight >= 7 &&
              (coefficientRing != globalQQ || density >= 0.25 ||
               supportSquareFavorsComplete);
          if (completeWins)
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurComplete
                : PowerSumsToTargetRoute::ViaSchurComplete;
        }

      if (hasCombinatorialTag(input.combinatorialTags,
                              CombinatorialTag::HorizontalPieri) ||
          hasCombinatorialTag(input.combinatorialTags,
                              CombinatorialTag::VerticalPieri))
        {
          bool completeWins = weight >= 8 &&
              (coefficientRing != globalQQ || density >= 0.25 ||
               supportSquareFavorsComplete);
          if (completeWins)
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurComplete
                : PowerSumsToTargetRoute::ViaSchurComplete;
        }

      size_t completeFriendlyTerms = 0;
      std::vector<int> commonPowerSumParts;
      bool firstPowerSumIndex = true;
      const auto *poly = polyValue(input.expression);
      for (const auto& term : poly->terms)
        {
          Partition index;
          if (!powerSumIndexFromMonomial(term.monomial, index))
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusRimHooks
                : PowerSumsToTargetRoute::ViaSchurAbacusRimHooks;
          if (completeFriendlyPowerSumIndex(index, weight))
            ++completeFriendlyTerms;
          if (input.combinatorialTags == borderStripsTag)
            {
              if (firstPowerSumIndex)
                {
                  commonPowerSumParts.assign(index.begin(), index.end());
                  commonPowerSumParts.erase(
                      std::unique(commonPowerSumParts.begin(),
                                  commonPowerSumParts.end()),
                      commonPowerSumParts.end());
                  firstPowerSumIndex = false;
                }
              else
                commonPowerSumParts.erase(
                    std::remove_if(
                        commonPowerSumParts.begin(),
                        commonPowerSumParts.end(),
                        [&](int part) {
                          return std::find(index.begin(), index.end(), part) ==
                                 index.end();
                        }),
                    commonPowerSumParts.end());
            }
        }

      // A retained p_r factor appears as a common cycle in every index.
      // Complete conversion wins when that cycle is small relative to the
      // remaining Schur-shaped weight; a dominant retained cycle favors
      // termwise rim hooks instead.
      bool hasSmallCommonCycle = std::any_of(
          commonPowerSumParts.begin(),
          commonPowerSumParts.end(),
          [&](int part) {
            int crossoverPercent = coefficientRing == globalQQ ? 36 : 41;
            return 100 * part <= crossoverPercent * weight;
          });
      bool hasCommonOne = std::find(commonPowerSumParts.begin(),
                                    commonPowerSumParts.end(),
                                    1) != commonPowerSumParts.end();
      if (input.combinatorialTags == borderStripsTag && weight >= 16 &&
          termCount >= 8 && !hasCommonOne && hasSmallCommonCycle)
        return omega
            ? PowerSumsToTargetRoute::ViaOmegaThenSchurComplete
            : PowerSumsToTargetRoute::ViaSchurComplete;
      if (input.combinatorialTags == borderStripsTag && weight >= 8 &&
          termCount >= 8 && hasCommonOne &&
          (coefficientRing != globalQQ || density >= 0.25 ||
           supportSquareFavorsComplete))
        return omega
            ? PowerSumsToTargetRoute::ViaOmegaThenSchurComplete
            : PowerSumsToTargetRoute::ViaSchurComplete;

      if ((input.combinatorialTags == 0 ||
           input.combinatorialTags == borderStripsTag) &&
          weight >= 14)
        {
          if (termCount >= 2 && completeFriendlyTerms == termCount)
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurComplete
                : PowerSumsToTargetRoute::ViaSchurComplete;

          // Mixed supports should not force every term through one route.
          // A nontrivial short-cycle group is converted through h, while the
          // remaining long-cycle terms retain the lower-overhead abacus path.
          // A scattered minority does not amortize evaluating and merging two
          // different transition systems.  Require a substantial coherent
          // short-cycle subgroup; the stricter per-index predicate above
          // supplies the coherence test.
          size_t minimumCompleteGroup = std::max<size_t>(
              8, (termCount + 3) / 4);
          if (termCount >= 8 &&
              completeFriendlyTerms >= minimumCompleteGroup &&
              completeFriendlyTerms < termCount)
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusAndComplete
                : PowerSumsToTargetRoute::ViaSchurAbacusAndComplete;
        }

      return omega
          ? PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusRimHooks
          : PowerSumsToTargetRoute::ViaSchurAbacusRimHooks;
  }

SymmetricEngineRing::PowerSumsToTargetRoute
SymmetricEngineRing::selectPowerSumsToTargetRoute(
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay) const
{
    if (targetBasisId == pBasisId)
      return PowerSumsToTargetRoute::AlreadyInTarget;
    switch (basisKindForId(targetBasisId))
      {
        case BasisKind::Schur:
          return selectPowerSumsToSchurRoute(input, false);
        case BasisKind::SchurOmega:
          return selectPowerSumsToSchurRoute(input, true);
        case BasisKind::Complete:
          return PowerSumsToTargetRoute::ViaCompleteLogarithmFormula;
        case BasisKind::Elementary:
          return PowerSumsToTargetRoute::ViaElementaryLogarithmFormula;
        case BasisKind::HallLittlewoodQGenerator:
        case BasisKind::HallLittlewoodBGenerator:
          return PowerSumsToTargetRoute::ViaHallLittlewoodGeneratorLogarithmFormula;
        case BasisKind::HallLittlewoodQ:
        case BasisKind::HallLittlewoodP:
        case BasisKind::HallLittlewoodB:
        case BasisKind::HallLittlewoodPOmega:
          {
        const auto *poly = polyValue(input.expression);
        Partition singleIndex;
        bool hasSinglePowerSumIndex =
            poly->terms.size() == 1 &&
            powerSumIndexFromMonomial(
                poly->terms.front().monomial, singleIndex) &&
            !singleIndex.empty();
        const char *forcedRoute = std::getenv(
            "M2_SYMMETRIC_RINGS_FORCE_HALL_POWER_SUM_ROUTE");
        if (forcedRoute != nullptr)
          {
            std::string route(forcedRoute);
            if (route == "green-duality" && hasSinglePowerSumIndex)
              return PowerSumsToTargetRoute::ViaHallLittlewoodGreenPolynomialsAndDuality;
            if (route == "triangular")
              return PowerSumsToTargetRoute::ViaHallLittlewoodTriangularReduction;
          }
        bool allTermsAreSingleCycles = !poly->terms.empty();
        for (const auto& term : poly->terms)
          {
            Partition index;
            if (!powerSumIndexFromMonomial(term.monomial, index) ||
                index.size() > 1)
              {
                allTermsAreSingleCycles = false;
                break;
              }
          }
        if (allTermsAreSingleCycles)
          return PowerSumsToTargetRoute::ViaHallLittlewoodSingleCycleGreenPolynomials;
        if (hasSinglePowerSumIndex)
          return PowerSumsToTargetRoute::ViaHallLittlewoodGreenPolynomialsAndDuality;
            return PowerSumsToTargetRoute::ViaHallLittlewoodTriangularReduction;
          }
        case BasisKind::Monomial:
          return PowerSumsToTargetRoute::ViaMonomialTransition;
        case BasisKind::Forgotten:
          return PowerSumsToTargetRoute::ViaForgottenTransition;
        case BasisKind::PowerSum:
          return PowerSumsToTargetRoute::AlreadyInTarget;
        case BasisKind::Custom:
          return PowerSumsToTargetRoute::ViaTermwiseFallback;
      }
    return PowerSumsToTargetRoute::ViaTermwiseFallback;
  }

const char *SymmetricEngineRing::powerSumsToTargetRouteName(
    PowerSumsToTargetRoute route) const
{
    switch (route)
      {
        case PowerSumsToTargetRoute::AlreadyInTarget:
          return "p->p:identity";
        case PowerSumsToTargetRoute::ViaSchurDegreeBlocks:
          return "p->S:degree-block-dispatch";
        case PowerSumsToTargetRoute::ViaSchurBorderStrips:
          return "p->S:border-strips";
        case PowerSumsToTargetRoute::ViaSchurAbacusRimHooks:
          return "p->S:abacus-rim-hooks";
        case PowerSumsToTargetRoute::ViaSchurAbacusAndComplete:
          return "p->S:abacus-rim-hooks+complete";
        case PowerSumsToTargetRoute::ViaSchurComplete:
          return "p->h->S:recursive-transition";
        case PowerSumsToTargetRoute::ViaSchurCharacters:
          return "p->S:grouped-characters";
        case PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips:
          return "p->omega(p)->S->Somega:border-strips";
        case PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusRimHooks:
          return "p->omega(p)->S->Somega:abacus-rim-hooks";
        case PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusAndComplete:
          return "p->omega(p)->S->Somega:abacus-rim-hooks+complete";
        case PowerSumsToTargetRoute::ViaOmegaThenSchurComplete:
          return "p->omega(p)->h->S->Somega";
        case PowerSumsToTargetRoute::ViaOmegaThenSchurCharacters:
          return "p->omega(p)->S->Somega:grouped-characters";
        case PowerSumsToTargetRoute::ViaCompleteLogarithmFormula:
          return "p->h:logarithm-formula";
        case PowerSumsToTargetRoute::ViaElementaryLogarithmFormula:
          return "p->e:logarithm-formula";
        case PowerSumsToTargetRoute::ViaHallLittlewoodGeneratorLogarithmFormula:
          return "p->q/b:logarithm-formula";
        case PowerSumsToTargetRoute::ViaHallLittlewoodSingleCycleGreenPolynomials:
          return "single-cycle-p-terms->Hall-Littlewood:Green-polynomials";
        case PowerSumsToTargetRoute::ViaHallLittlewoodGreenPolynomialsAndDuality:
          return "p_mu->Hall-Littlewood:Green-polynomials-via-duality";
        case PowerSumsToTargetRoute::ViaHallLittlewoodTriangularReduction:
          return "p->Hall-Littlewood:triangular-reduction";
        case PowerSumsToTargetRoute::ViaMonomialTransition:
          return "p->m:transition";
        case PowerSumsToTargetRoute::ViaForgottenTransition:
          return "p->ff:transition";
        case PowerSumsToTargetRoute::ViaTermwiseFallback:
          return "p->target:termwise-fallback";
      }
    return "p->target:unknown";
  }

void SymmetricEngineRing::tracePowerSumsToTargetSelection(
        PowerSumsToTargetRoute route,
        const std::string& targetDisplay) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::fprintf(stderr,
                 "SymmetricRings power-sums-target: target=%s route=%s\n",
                 targetDisplay.c_str(),
                 powerSumsToTargetRouteName(route));
  }

ring_elem SymmetricEngineRing::executePowerSumsToTargetRoute(
        PowerSumsToTargetRoute route,
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder,
        bool targetIsMultiplicative) const
{
    auto viaSchurOmega = [&]() {
        ring_elem omegaInputExpression = omegaPowerSums(input.expression);
        if (error()) return zero();
        int schurId = requiredBasisIdForKind(BasisKind::Schur);
        if (error()) return zero();
        int schurOrder = basisOrderForId(schurId);
        std::string schurDisplay = displayForBasis(schurId);
        ring_elem inSchur;
        switch (route)
          {
            case PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips:
              inSchur = powerSumsToSchurViaBorderStrips(
                  omegaInputExpression, schurId, schurDisplay, schurOrder);
              break;
            case PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusRimHooks:
              inSchur = powerSumsToSchurViaAbacusRimHooks(
                  omegaInputExpression, schurId, schurDisplay, schurOrder);
              break;
            case PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusAndComplete:
              inSchur = powerSumsToSchurViaAbacusAndComplete(
                  omegaInputExpression, schurId, schurDisplay, schurOrder);
              break;
            case PowerSumsToTargetRoute::ViaOmegaThenSchurCharacters:
              inSchur = powerSumsToSchurLikeViaCharacters(
                  omegaInputExpression, schurId, schurOrder,
                  schurDisplay, false);
              break;
            case PowerSumsToTargetRoute::ViaOmegaThenSchurComplete:
              inSchur = powerSumsToSchurViaComplete(
                  omegaInputExpression, schurId, schurDisplay, schurOrder);
              break;
            default:
              ERROR("expected a Schur Omega power-sum route");
              return zero();
          }
        if (error()) return zero();
        return replaceSingleBasis(inSchur, schurId, targetBasisId);
      };

    switch (route)
      {
        case PowerSumsToTargetRoute::AlreadyInTarget:
          return copyPolyValue(polyValue(input.expression));
        case PowerSumsToTargetRoute::ViaSchurDegreeBlocks:
          return runPowerSumsToSchurDegreeBlockPipeline(
              input, pBasisId, targetBasisId, targetDisplay,
              targetDisplayOrder, targetIsMultiplicative);
        case PowerSumsToTargetRoute::ViaSchurBorderStrips:
          return powerSumsToSchurViaBorderStrips(
              input.expression, targetBasisId, targetDisplay, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaSchurAbacusRimHooks:
          return powerSumsToSchurViaAbacusRimHooks(
              input.expression, targetBasisId, targetDisplay, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaSchurAbacusAndComplete:
          return powerSumsToSchurViaAbacusAndComplete(
              input.expression, targetBasisId, targetDisplay, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaSchurCharacters:
          return powerSumsToSchurLikeViaCharacters(
              input.expression, targetBasisId, targetDisplayOrder,
              targetDisplay, false);
        case PowerSumsToTargetRoute::ViaSchurComplete:
          return powerSumsToSchurViaComplete(
              input.expression, targetBasisId, targetDisplay, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips:
        case PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusRimHooks:
        case PowerSumsToTargetRoute::ViaOmegaThenSchurAbacusAndComplete:
        case PowerSumsToTargetRoute::ViaOmegaThenSchurCharacters:
        case PowerSumsToTargetRoute::ViaOmegaThenSchurComplete:
          return viaSchurOmega();
        case PowerSumsToTargetRoute::ViaCompleteLogarithmFormula:
          return powerSumsToCompleteViaLogarithmFormula(
              input.expression, targetBasisId, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaElementaryLogarithmFormula:
          return powerSumsToElementaryViaLogarithmFormula(
              input.expression, targetBasisId, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaHallLittlewoodGeneratorLogarithmFormula:
          {
            CoeffMap generators = powerSumsToHallGeneratorMapViaLogarithmFormula(
                input.expression,
                basisKindForId(targetBasisId) ==
                    BasisKind::HallLittlewoodBGenerator);
            if (error()) return zero();
            return coeffMapToElement(generators,
                                     targetBasisId,
                                     targetDisplay,
                                     targetDisplayOrder,
                                     true);
          }
        case PowerSumsToTargetRoute::ViaHallLittlewoodSingleCycleGreenPolynomials:
          return powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials(
              input.expression, targetBasisId, targetDisplay, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaHallLittlewoodGreenPolynomialsAndDuality:
          return powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality(
              input.expression, targetBasisId, targetDisplay, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaHallLittlewoodTriangularReduction:
          return powerSumsToHallLittlewoodViaTriangularReduction(
              input.expression, targetBasisId, targetDisplay, targetDisplayOrder);
        case PowerSumsToTargetRoute::ViaMonomialTransition:
        case PowerSumsToTargetRoute::ViaForgottenTransition:
        case PowerSumsToTargetRoute::ViaTermwiseFallback:
          return powerSumsToTargetViaTermwiseConversion(
              input.expression, targetBasisId, targetDisplay,
              targetDisplayOrder, targetIsMultiplicative);
      }
    ERROR("unknown power-sums-to-target conversion route");
    return zero();
  }

ring_elem SymmetricEngineRing::powerSumsToTargetDispatch(
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder,
        bool targetIsMultiplicative) const
{
    PowerSumsToTargetRoute route = selectPowerSumsToTargetRoute(
        input, pBasisId, targetBasisId, targetDisplay);
    tracePowerSumsToTargetSelection(route, targetDisplay);
    return executePowerSumsToTargetRoute(route,
                                         input,
                                         pBasisId,
                                         targetBasisId,
                                         targetDisplay,
                                         targetDisplayOrder,
                                         targetIsMultiplicative);
  }

ring_elem SymmetricEngineRing::runPowerSumsToSchurDegreeBlockPipeline(
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder,
        bool targetIsMultiplicative) const
{
    GCMap<int, VECTOR(SymmetricTerm)> termsByDegree;
    const auto *poly = polyValue(input.expression);
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during degree-block conversion");
            return zero();
          }
        termsByDegree[partitionWeight(index)].push_back(
            {coefficientRing->copy(term.coeff), term.monomial});
      }

    ring_elem result = zero();
    for (auto& degreeTerms : termsByDegree)
      {
        ring_elem blockExpression = fromTermVector(degreeTerms.second, true);
        ConversionInput blockInput{
            blockExpression,
            inferConversionGuarantees(blockExpression, targetBasisId),
            input.combinatorialTags};
        blockInput.guarantees.pureBasis = pBasisId;
        blockInput.guarantees.expandedBasis = pBasisId;
        blockInput.guarantees = strengthenConversionGuarantees(
            std::move(blockInput.guarantees), targetBasisId);
        PowerSumsToTargetRoute blockRoute = selectPowerSumsToTargetRoute(
            blockInput, pBasisId, targetBasisId, targetDisplay);
        if (blockRoute == PowerSumsToTargetRoute::ViaSchurDegreeBlocks)
          {
            ERROR("power-sum degree block did not have a homogeneous weight");
            return zero();
          }
        tracePowerSumsToTargetSelection(blockRoute, targetDisplay);
        ring_elem converted = executePowerSumsToTargetRoute(
            blockRoute,
            blockInput,
            pBasisId,
            targetBasisId,
            targetDisplay,
            targetDisplayOrder,
            targetIsMultiplicative);
        if (error()) return zero();
        result = add(result, converted);
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumsToSchurViaComplete(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder) const
{
    CoeffMap inComplete = powerSumsToCompleteMapViaLogarithmFormula(f);
    if (error()) return zero();
    CoeffMap result = completeToSchurCoefficientsViaRecursiveTransition(
        inComplete);
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
        int weight = partitionWeight(index);
        if (completeFriendlyPowerSumIndex(index, weight))
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

ring_elem SymmetricEngineRing::runPowerSumsPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    return sourceToTargetDispatch(input,
                                  pBasisId,
                                  pDisplay,
                                  pOrder,
                                  pIsMultiplicative,
                                  targetBasisId,
                                  targetDisplay,
                                  targetOrder,
                                  targetIsMultiplicative);
  }

// ============================================================================
// General Source-To-Target Dispatch
// ============================================================================
// General conversion selects a direct route or composes source-to-p and p-to-target.

SymmetricEngineRing::SourceToTargetRoute
SymmetricEngineRing::selectSourceToTargetRoute(
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay) const
{
    if (input.guarantees.targetClosed == KnownState::True)
      return SourceToTargetRoute::AlreadyInTarget;

    if (input.guarantees.expandedBasis &&
        input.guarantees.normalized == KnownState::True &&
        input.guarantees.skewFree == KnownState::True)
      {
        const BasisKind sourceKind =
            basisKindForId(*input.guarantees.expandedBasis);
        const BasisKind targetKind = basisKindForId(targetBasisId);
        if ((sourceKind == BasisKind::HallLittlewoodQ && targetKind == BasisKind::HallLittlewoodP) ||
            (sourceKind == BasisKind::HallLittlewoodP && targetKind == BasisKind::HallLittlewoodQ) ||
            (sourceKind == BasisKind::HallLittlewoodB && targetKind == BasisKind::HallLittlewoodPOmega) ||
            (sourceKind == BasisKind::HallLittlewoodPOmega && targetKind == BasisKind::HallLittlewoodB))
          return SourceToTargetRoute::ViaHallLittlewoodNormalization;
        if ((sourceKind == BasisKind::Schur && targetKind == BasisKind::SchurOmega) ||
            (sourceKind == BasisKind::SchurOmega && targetKind == BasisKind::Schur))
          return SourceToTargetRoute::ViaSchurOmegaConjugation;
      }

    if (input.guarantees.expandedBasis == pBasisId)
      return SourceToTargetRoute::ViaSelectedPowerSumsToTargetRoute;

    if (basisKindForId(targetBasisId) == BasisKind::Schur)
      {
        int completeId = requiredBasisIdForKind(BasisKind::Complete);
        if (error()) return SourceToTargetRoute::ViaPowerSumsThenCompleteThenSchur;
        if (input.guarantees.pureBasis == completeId)
          return SourceToTargetRoute::ViaCompleteToSchurRecursiveTransition;
        // Monomial and forgotten transitions naturally produce relatively
        // sparse power-sum expansions.  Let the ordinary p-to-S selector use
        // that actual support instead of unconditionally forcing the
        // source->p->h->S composition used by denser transformed bases.
        std::optional<int> sourceBasis = input.guarantees.expandedBasis
            ? input.guarantees.expandedBasis
            : input.guarantees.pureBasis;
        if (sourceBasis)
          {
            BasisKind sourceKind = basisKindForId(*sourceBasis);
            if (sourceKind == BasisKind::Monomial ||
                sourceKind == BasisKind::Forgotten)
              return SourceToTargetRoute::ViaSourceToPowerSumsThenTarget;
          }
        return SourceToTargetRoute::ViaPowerSumsThenCompleteThenSchur;
      }

    return SourceToTargetRoute::ViaSourceToPowerSumsThenTarget;
  }

const char *SymmetricEngineRing::sourceToTargetRouteName(
    SourceToTargetRoute route) const
{
    switch (route)
      {
        case SourceToTargetRoute::AlreadyInTarget:
          return "already-in-target";
        case SourceToTargetRoute::ViaHallLittlewoodNormalization:
          return "Hall-Littlewood-capital-normalized:diagonal-scaling";
        case SourceToTargetRoute::ViaSchurOmegaConjugation:
          return "Schur-omega:partition-conjugation";
        case SourceToTargetRoute::ViaSelectedPowerSumsToTargetRoute:
          return "expanded-p->selected-p-target-route";
        case SourceToTargetRoute::ViaCompleteToSchurRecursiveTransition:
          return "h->S:recursive-transition";
        case SourceToTargetRoute::ViaPowerSumsThenCompleteThenSchur:
          return "source->p->h->S";
        case SourceToTargetRoute::ViaSourceToPowerSumsThenTarget:
          return "source->p->target";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceSourceToTargetSelection(
        SourceToTargetRoute route,
        const ConversionInput& input,
        const std::string& targetDisplay) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::string source = "mixed";
    if (input.guarantees.expandedBasis)
      source = displayForBasis(*input.guarantees.expandedBasis);
    else if (input.guarantees.pureBasis)
      source = displayForBasis(*input.guarantees.pureBasis);
    std::fprintf(stderr,
                 "SymmetricRings source-target: source=%s target=%s route=%s\n",
                 source.c_str(),
                 targetDisplay.c_str(),
                 sourceToTargetRouteName(route));
  }

ring_elem SymmetricEngineRing::executeSourceToTargetRoute(
        SourceToTargetRoute route,
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    auto finish = [&](ring_elem result) {
      if (!error())
        {
          ConversionGuarantees guarantees =
              guaranteesAfterSourceTargetConversion(
                  input.guarantees, result, targetBasisId);
          attachConversionGuarantees(
              result, guarantees, input.combinatorialTags);
        }
      return result;
    };

    switch (route)
      {
        case SourceToTargetRoute::AlreadyInTarget:
          return finish(copyPolyValue(polyValue(input.expression)));
        case SourceToTargetRoute::ViaHallLittlewoodNormalization:
          {
            int sourceBasisId = *input.guarantees.expandedBasis;
            BasisKind sourceKind = basisKindForId(sourceBasisId);
            bool capitalToNormalized =
                sourceKind == BasisKind::HallLittlewoodQ ||
                sourceKind == BasisKind::HallLittlewoodB;
            return finish(
                hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
                    input.expression, sourceBasisId, targetBasisId,
                    targetDisplay, targetOrder, capitalToNormalized));
          }
        case SourceToTargetRoute::ViaSchurOmegaConjugation:
          return finish(schurOmegaConversionViaPartitionConjugation(
              input.expression,
              *input.guarantees.expandedBasis,
              targetBasisId,
              targetDisplay,
              targetOrder));
        case SourceToTargetRoute::ViaSelectedPowerSumsToTargetRoute:
          return finish(powerSumsToTargetDispatch(
              input, pBasisId, targetBasisId, targetDisplay,
              targetOrder, targetIsMultiplicative));
        case SourceToTargetRoute::ViaCompleteToSchurRecursiveTransition:
          {
            int completeId = requiredBasisIdForKind(BasisKind::Complete);
            if (error()) return zero();
            return finish(completeToSchurViaRecursiveTransition(
                input.expression,
                completeId,
                displayForBasis(completeId),
                basisOrderForId(completeId),
                isMultiplicativeBasis(completeId),
                targetBasisId,
                targetDisplay,
                targetOrder));
          }
        case SourceToTargetRoute::ViaPowerSumsThenCompleteThenSchur:
        case SourceToTargetRoute::ViaSourceToPowerSumsThenTarget:
          break;
      }

    ring_elem inPowerSums =
        expressionToPowerSumsViaBasisElementRoutes(input.expression);
    if (error()) return zero();
    ConversionInput powerSumInput{
        inPowerSums,
        inferConversionGuarantees(inPowerSums, targetBasisId),
        input.combinatorialTags};
    powerSumInput.guarantees.pureBasis = pBasisId;
    powerSumInput.guarantees.expandedBasis = pBasisId;
    if (!powerSumInput.guarantees.homogeneousWeight &&
        input.guarantees.homogeneousWeight)
      powerSumInput.guarantees.homogeneousWeight =
          input.guarantees.homogeneousWeight;
    powerSumInput.guarantees = strengthenConversionGuarantees(
        std::move(powerSumInput.guarantees), targetBasisId);

    switch (route)
      {
        case SourceToTargetRoute::ViaPowerSumsThenCompleteThenSchur:
          {
            int completeId = requiredBasisIdForKind(BasisKind::Complete);
            if (error()) return zero();
            int completeOrder = basisOrderForId(completeId);
            bool completeIsMultiplicative = isMultiplicativeBasis(completeId);
            ring_elem inComplete = powerSumsToCompleteViaLogarithmFormula(
                inPowerSums, completeId, completeOrder);
            if (error()) return zero();
            return finish(completeToSchurViaRecursiveTransition(
                inComplete,
                completeId,
                displayForBasis(completeId),
                completeOrder,
                completeIsMultiplicative,
                targetBasisId,
                targetDisplay,
                targetOrder));
          }
        case SourceToTargetRoute::ViaSourceToPowerSumsThenTarget:
          return finish(powerSumsToTargetDispatch(
              powerSumInput, pBasisId, targetBasisId, targetDisplay,
              targetOrder, targetIsMultiplicative));
        case SourceToTargetRoute::AlreadyInTarget:
        case SourceToTargetRoute::ViaHallLittlewoodNormalization:
        case SourceToTargetRoute::ViaSchurOmegaConjugation:
        case SourceToTargetRoute::ViaSelectedPowerSumsToTargetRoute:
        case SourceToTargetRoute::ViaCompleteToSchurRecursiveTransition:
          ERROR("source-to-target route reached an invalid execution stage");
          return zero();
      }
    ERROR("unknown source-to-target route");
    return zero();
  }

ring_elem SymmetricEngineRing::sourceToTargetDispatch(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    requireBasis(pBasisId);
    requireBasis(targetBasisId);
    SourceToTargetRoute route = selectSourceToTargetRoute(
        input, pBasisId, targetBasisId, targetDisplay);
    if (error()) return zero();
    traceSourceToTargetSelection(route, input, targetDisplay);
    return executeSourceToTargetRoute(route,
                                      input,
                                      pBasisId,
                                      targetBasisId,
                                      targetDisplay,
                                      targetOrder,
                                      targetIsMultiplicative);
  }

// ============================================================================
// Grouped Conversion Pipelines
// ============================================================================
// Grouped pipelines preserve whole-expression structure for suitable targets.

bool SymmetricEngineRing::canUseGroupedHallLittlewoodPipeline(
        const ConversionGuarantees& guarantees,
        int targetBasisId,
        const std::string& targetDisplay) const
{
    const BasisKind targetKind = basisKindForId(targetBasisId);
    if (targetKind != BasisKind::HallLittlewoodQ &&
        targetKind != BasisKind::HallLittlewoodP &&
        targetKind != BasisKind::HallLittlewoodB &&
        targetKind != BasisKind::HallLittlewoodPOmega)
      return false;
    if (guarantees.normalized != KnownState::True ||
        guarantees.skewFree != KnownState::True || !guarantees.factorBases)
      return false;
    const BasisKind generatorKind =
        targetKind == BasisKind::HallLittlewoodQ ||
        targetKind == BasisKind::HallLittlewoodP
            ? BasisKind::HallLittlewoodQGenerator
            : BasisKind::HallLittlewoodBGenerator;
    for (int basisId : *guarantees.factorBases)
      if (basisKindForId(basisId) != generatorKind) return false;
    const char *forced = std::getenv(
        "M2_SYMMETRIC_RINGS_FORCE_HALL_LITTLEWOOD_PIPELINE");
    if (forced != nullptr)
      {
        std::string choice(forced);
        if (choice == "grouped") return true;
        if (choice == "fallback") return false;
      }
    return true;
  }

ring_elem SymmetricEngineRing::runGroupedMultiplicativeTargetPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    ring_elem result;
    if (tryExpressionToTarget(input.expression,
                              targetBasisId,
                              targetDisplay,
                              targetOrder,
                              targetIsMultiplicative,
                              result))
      {
        ConversionGuarantees guarantees =
            guaranteesAfterSourceTargetConversion(
                input.guarantees, result, targetBasisId);
        attachConversionGuarantees(
            result, guarantees, input.combinatorialTags);
        return result;
      }
    if (error()) return zero();
    return sourceToTargetDispatch(input,
                                  pBasisId,
                                  pDisplay,
                                  pOrder,
                                  pIsMultiplicative,
                                  targetBasisId,
                                  targetDisplay,
                                  targetOrder,
                                  targetIsMultiplicative);
  }

ring_elem SymmetricEngineRing::runGroupedHallLittlewoodPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    ring_elem result;
    if (tryExpressionToHallLittlewoodViaTriangularReduction(
            input.expression,
            targetBasisId,
            targetDisplay,
            targetOrder,
            result))
      {
        ConversionGuarantees guarantees =
            guaranteesAfterSourceTargetConversion(
                input.guarantees, result, targetBasisId);
        attachConversionGuarantees(
            result, guarantees, input.combinatorialTags);
        return result;
      }
    if (error()) return zero();
    return sourceToTargetDispatch(input,
                                  pBasisId,
                                  pDisplay,
                                  pOrder,
                                  pIsMultiplicative,
                                  targetBasisId,
                                  targetDisplay,
                                  targetOrder,
                                  targetIsMultiplicative);
  }

// ============================================================================
// Generic Basis-Element And Expression Dispatch
// ============================================================================
// Route selection and expression workflows used by the general fallback.

SymmetricEngineRing::BasisElementToPowerSumsRoute
SymmetricEngineRing::selectBasisElementToPowerSumsRoute(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    BasisKind kind = basisKindForId(atomBasisIdAt(monomial, pos));
    if (atomIsSkewAt(monomial, pos))
      {
        switch (kind)
          {
            case BasisKind::Schur:
              return BasisElementToPowerSumsRoute::ViaSkewSchurJacobiTrudiComplete;
            case BasisKind::SchurOmega:
              return BasisElementToPowerSumsRoute::ViaSkewSchurOmegaJacobiTrudiElementary;
            case BasisKind::HallLittlewoodQ:
            case BasisKind::HallLittlewoodB:
            case BasisKind::HallLittlewoodP:
            case BasisKind::HallLittlewoodPOmega:
              return BasisElementToPowerSumsRoute::ViaSkewHallLittlewood;
            case BasisKind::Custom:
            case BasisKind::PowerSum:
            case BasisKind::Complete:
            case BasisKind::Elementary:
            case BasisKind::Monomial:
            case BasisKind::Forgotten:
            case BasisKind::HallLittlewoodQGenerator:
            case BasisKind::HallLittlewoodBGenerator:
              return BasisElementToPowerSumsRoute::NoApplicableRoute;
          }
      }
    switch (kind)
      {
        case BasisKind::PowerSum:
          return BasisElementToPowerSumsRoute::AlreadyPowerSums;
        case BasisKind::Complete:
          return BasisElementToPowerSumsRoute::ViaCompleteClassicalFormula;
        case BasisKind::Elementary:
          return BasisElementToPowerSumsRoute::ViaElementaryClassicalFormula;
        case BasisKind::HallLittlewoodQGenerator:
        case BasisKind::HallLittlewoodBGenerator:
          return BasisElementToPowerSumsRoute::ViaHallLittlewoodGeneratorClassicalFormula;
        case BasisKind::Schur:
          return BasisElementToPowerSumsRoute::ViaSchurCharacters;
        case BasisKind::SchurOmega:
          return BasisElementToPowerSumsRoute::ViaSchurOmegaCharacters;
        case BasisKind::Monomial:
          return BasisElementToPowerSumsRoute::ViaMonomialTransition;
        case BasisKind::Forgotten:
          return BasisElementToPowerSumsRoute::ViaForgottenTransition;
        case BasisKind::HallLittlewoodQ:
        case BasisKind::HallLittlewoodB:
          return BasisElementToPowerSumsRoute::ViaHallLittlewoodRaisingOperators;
        case BasisKind::HallLittlewoodP:
        case BasisKind::HallLittlewoodPOmega:
          return BasisElementToPowerSumsRoute::ViaHallLittlewoodCapitalNormalization;
        case BasisKind::Custom:
          return BasisElementToPowerSumsRoute::NoApplicableRoute;
      }
    return BasisElementToPowerSumsRoute::NoApplicableRoute;
  }

const char *SymmetricEngineRing::basisElementToPowerSumsRouteName(
    BasisElementToPowerSumsRoute route) const
{
    switch (route)
      {
        case BasisElementToPowerSumsRoute::AlreadyPowerSums: return "p:identity";
        case BasisElementToPowerSumsRoute::ViaCompleteClassicalFormula: return "h->p:classical-formula";
        case BasisElementToPowerSumsRoute::ViaElementaryClassicalFormula: return "e->p:classical-formula";
        case BasisElementToPowerSumsRoute::ViaHallLittlewoodGeneratorClassicalFormula: return "q/b->p:classical-formula";
        case BasisElementToPowerSumsRoute::ViaSchurCharacters: return "S->p:characters";
        case BasisElementToPowerSumsRoute::ViaSchurOmegaCharacters: return "Somega->p:characters";
        case BasisElementToPowerSumsRoute::ViaMonomialTransition: return "m->p:transition";
        case BasisElementToPowerSumsRoute::ViaForgottenTransition: return "ff->p:transition";
        case BasisElementToPowerSumsRoute::ViaHallLittlewoodRaisingOperators: return "Q/B->p:raising-operators";
        case BasisElementToPowerSumsRoute::ViaHallLittlewoodCapitalNormalization: return "P/Pomega->Q/B->p";
        case BasisElementToPowerSumsRoute::ViaSkewSchurJacobiTrudiComplete: return "skew-S->h->p:Jacobi-Trudi";
        case BasisElementToPowerSumsRoute::ViaSkewSchurOmegaJacobiTrudiElementary: return "skew-Somega->e->p:Jacobi-Trudi";
        case BasisElementToPowerSumsRoute::ViaSkewHallLittlewood: return "skew-Hall-Littlewood->p";
        case BasisElementToPowerSumsRoute::NoApplicableRoute: return "not-applicable";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceBasisElementToPowerSumsSelection(
    BasisElementToPowerSumsRoute route,
    const std::string& sourceDisplay) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::fprintf(stderr,
                 "SymmetricRings basis-element-to-power-sums: source=%s route=%s\n",
                 sourceDisplay.c_str(),
                 basisElementToPowerSumsRouteName(route));
  }

ring_elem SymmetricEngineRing::executeBasisElementToPowerSumsRoute(
    BasisElementToPowerSumsRoute route,
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    std::string display = displayForBasis(atomBasisIdAt(monomial, pos));
    if (route == BasisElementToPowerSumsRoute::ViaSkewSchurJacobiTrudiComplete ||
        route == BasisElementToPowerSumsRoute::ViaSkewSchurOmegaJacobiTrudiElementary)
      {
        BasisKind generatorKind =
            route == BasisElementToPowerSumsRoute::ViaSkewSchurJacobiTrudiComplete
                ? BasisKind::Complete : BasisKind::Elementary;
        int generatorId = requiredBasisIdForKind(generatorKind);
        if (error()) return zero();
        std::string generatorDisplay = displayForBasis(generatorId);
        return expressionToPowerSumsViaBasisElementRoutes(
            jacobiTrudi(basisElementOuterIndex(monomial, pos),
                        basisElementInnerIndex(monomial, pos),
                        generatorId));
      }
    if (route == BasisElementToPowerSumsRoute::ViaSkewHallLittlewood)
      return skewHallLittlewoodToPowerSums(
                                           basisElementOuterIndex(monomial, pos),
                                           basisElementInnerIndex(monomial, pos),
                                           basisKindForId(
                                               atomBasisIdAt(monomial, pos)));

    Partition index = basisElementIndex(monomial, pos);
    if (route == BasisElementToPowerSumsRoute::AlreadyPowerSums)
      return basisElementFromIndex(registeredPowerSumBasisId(), index);
    if (route == BasisElementToPowerSumsRoute::ViaCompleteClassicalFormula ||
        route == BasisElementToPowerSumsRoute::ViaElementaryClassicalFormula)
      {
        ring_elem result = one();
        for (int part : index)
          {
            ring_elem factor =
                route == BasisElementToPowerSumsRoute::ViaCompleteClassicalFormula
                    ? completePartToPowerSumsViaClassicalFormula(part)
                    : elementaryPartToPowerSumsViaClassicalFormula(part);
            result = mult(result, factor);
          }
        return result;
      }
    if (route == BasisElementToPowerSumsRoute::ViaHallLittlewoodGeneratorClassicalFormula)
      {
        ring_elem result = one();
        for (int part : index)
          result = mult(result,
                        hallLittlewoodGeneratorPartToPowerSumsViaClassicalFormula(
                            part,
                            basisKindForId(atomBasisIdAt(monomial, pos)) ==
                                BasisKind::HallLittlewoodBGenerator));
        return result;
      }
    if (route == BasisElementToPowerSumsRoute::ViaSchurCharacters)
      return schurLikeToPowerSumsViaCharacters(index, false);
    if (route == BasisElementToPowerSumsRoute::ViaSchurOmegaCharacters)
      return schurLikeToPowerSumsViaCharacters(index, true);
    if (route == BasisElementToPowerSumsRoute::ViaMonomialTransition ||
        route == BasisElementToPowerSumsRoute::ViaForgottenTransition)
      return monomialToPowerSumsViaTransitionMatrix(
          index, route == BasisElementToPowerSumsRoute::ViaForgottenTransition);
    if (route == BasisElementToPowerSumsRoute::ViaHallLittlewoodRaisingOperators)
      return hallLittlewoodCapitalToPowerSumsViaRaisingOperators(
          index,
          basisKindForId(atomBasisIdAt(monomial, pos)) ==
              BasisKind::HallLittlewoodB);
    if (route == BasisElementToPowerSumsRoute::ViaHallLittlewoodCapitalNormalization)
      return hallLittlewoodNormalizedToPowerSumsViaCapitalNormalization(
          index,
          basisKindForId(atomBasisIdAt(monomial, pos)) ==
              BasisKind::HallLittlewoodPOmega);
    if (atomIsSkewAt(monomial, pos))
      ERROR("basis conversion for skew ", display.c_str(), " basis elements is not implemented yet");
    else
      ERROR("basis conversion to power sums is not implemented for basis ", display.c_str());
    return zero();
  }

ring_elem SymmetricEngineRing::basisElementToPowerSumsDispatch(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    BasisElementToPowerSumsRoute route =
        selectBasisElementToPowerSumsRoute(monomial, pos);
    std::string sourceDisplay = displayForBasis(atomBasisIdAt(monomial, pos));
    traceBasisElementToPowerSumsSelection(route, sourceDisplay);
    return executeBasisElementToPowerSumsRoute(route, monomial, pos);
  }

ring_elem SymmetricEngineRing::monomialToPowerSumsViaBasisElementRoutes(const SymmetricMonomial& monomial) const
{
    ring_elem result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        result = mult(result, basisElementToPowerSumsDispatch(monomial, pos));
        if (error()) return zero();
        pos += atomLengthAt(monomial, pos);
      }
    return result;
  }

ring_elem SymmetricEngineRing::expressionToPowerSumsViaBasisElementRoutes(ring_elem f) const
{
    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        ring_elem converted =
            monomialToPowerSumsViaBasisElementRoutes(term.monomial);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, converted));
      }
    return result;
  }

bool SymmetricEngineRing::tryBasisElementToTarget(const SymmetricMonomial& monomial,
                          size_t pos,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          bool targetIsMultiplicative,
                          ring_elem& result) const
{
    int basisId = atomBasisIdAt(monomial, pos);
    std::string display = displayForBasis(basisId);
    const BasisKind sourceKind = basisKindForId(basisId);
    const BasisKind targetKind = basisKindForId(targetBasisId);
    Partition index = basisElementIndex(monomial, pos);

    if (!atomIsSkewAt(monomial, pos) && basisId == targetBasisId)
      {
        result = basisElementFromIndex(targetBasisId, index);
        return true;
      }

    if (targetKind == BasisKind::Complete && sourceKind == BasisKind::Schur)
      {
        Partition outer = atomIsSkewAt(monomial, pos) ? basisElementOuterIndex(monomial, pos)
                                                      : index;
        Partition inner = atomIsSkewAt(monomial, pos) ? basisElementInnerIndex(monomial, pos)
                                                      : Partition{};
        result = jacobiTrudi(outer, inner, targetBasisId);
        return true;
      }

    if (targetKind == BasisKind::Elementary &&
        sourceKind == BasisKind::SchurOmega)
      {
        Partition outer = atomIsSkewAt(monomial, pos) ? basisElementOuterIndex(monomial, pos)
                                                      : index;
        Partition inner = atomIsSkewAt(monomial, pos) ? basisElementInnerIndex(monomial, pos)
                                                      : Partition{};
        result = jacobiTrudi(outer, inner, targetBasisId);
        return true;
      }

    if (atomIsSkewAt(monomial, pos)) return false;
    if (sourceKind == BasisKind::Custom) return false;
    ring_elem inPowerSums = basisElementToPowerSumsDispatch(monomial, pos);
    if (error()) return false;
    int pBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return false;
    ConversionInput input{
        inPowerSums,
        inferConversionGuarantees(inPowerSums, targetBasisId),
        0};
    input.guarantees.pureBasis = pBasisId;
    input.guarantees.expandedBasis = pBasisId;
    input.guarantees = strengthenConversionGuarantees(
        std::move(input.guarantees), targetBasisId);
    result = powerSumsToTargetDispatch(input,
                                       pBasisId,
                                       targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       targetIsMultiplicative);
    return !error();
  }

bool SymmetricEngineRing::tryMonomialToTarget(const SymmetricMonomial& monomial,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetDisplayOrder,
                              bool targetIsMultiplicative,
                              ring_elem& result) const
{
    result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        ring_elem factor;
        if (!tryBasisElementToTarget(monomial,
                                pos,
                                targetBasisId,
                                targetDisplay,
                                targetDisplayOrder,
                                targetIsMultiplicative,
                                factor))
          return false;
        result = mult(result, factor);
        pos += atomLengthAt(monomial, pos);
      }
    return true;
  }

SymmetricEngineRing::ExpressionToTargetRoute
SymmetricEngineRing::selectExpressionToTargetRoute(
    int targetBasisId,
    bool targetIsMultiplicative) const
{
    switch (basisKindForId(targetBasisId))
      {
      case BasisKind::Schur:
      case BasisKind::SchurOmega:
        return ExpressionToTargetRoute::ViaSchurTriangularReduction;
      case BasisKind::HallLittlewoodQ:
      case BasisKind::HallLittlewoodP:
      case BasisKind::HallLittlewoodB:
      case BasisKind::HallLittlewoodPOmega:
        return ExpressionToTargetRoute::ViaHallLittlewoodTriangularReduction;
      case BasisKind::Complete:
      case BasisKind::Elementary:
        return ExpressionToTargetRoute::ViaFactorwiseConversion;
      default:
        return targetIsMultiplicative
            ? ExpressionToTargetRoute::ViaFactorwiseConversion
            : ExpressionToTargetRoute::NoApplicableRoute;
      }
  }

const char *SymmetricEngineRing::expressionToTargetRouteName(
    ExpressionToTargetRoute route) const
{
    switch (route)
      {
        case ExpressionToTargetRoute::ViaSchurTriangularReduction:
          return "Schur-triangular-reduction";
        case ExpressionToTargetRoute::ViaHallLittlewoodTriangularReduction:
          return "Hall-Littlewood-triangular-reduction";
        case ExpressionToTargetRoute::ViaFactorwiseConversion:
          return "factorwise-conversion";
        case ExpressionToTargetRoute::NoApplicableRoute:
          return "not-applicable";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceExpressionToTargetSelection(
    ExpressionToTargetRoute route,
    const std::string& targetDisplay) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::fprintf(stderr,
                 "SymmetricRings expression-target: target=%s route=%s\n",
                 targetDisplay.c_str(),
                 expressionToTargetRouteName(route));
  }

bool SymmetricEngineRing::executeExpressionToTargetRoute(
    ExpressionToTargetRoute route,
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder,
    bool targetIsMultiplicative,
    ring_elem& result) const
{
    if (route == ExpressionToTargetRoute::ViaSchurTriangularReduction)
      return tryExpressionToSchurViaTriangularReduction(
          f, targetBasisId, targetDisplay, targetDisplayOrder, result);
    if (route == ExpressionToTargetRoute::ViaHallLittlewoodTriangularReduction)
      return tryExpressionToHallLittlewoodViaTriangularReduction(
          f, targetBasisId, targetDisplay, targetDisplayOrder, result);
    if (route != ExpressionToTargetRoute::ViaFactorwiseConversion)
      return false;
    result = zero();
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        ring_elem converted;
        if (!tryMonomialToTarget(term.monomial,
                                    targetBasisId,
                                    targetDisplay,
                                    targetDisplayOrder,
                                    targetIsMultiplicative,
                                    converted))
          return false;
        result = add(result, scaled(term.coeff, converted));
      }
    return true;
  }

bool SymmetricEngineRing::tryExpressionToTarget(ring_elem f,
                             int targetBasisId,
                             const std::string& targetDisplay,
                             int targetDisplayOrder,
                             bool targetIsMultiplicative,
                             ring_elem& result) const
{
    ExpressionToTargetRoute route = selectExpressionToTargetRoute(
        targetBasisId, targetIsMultiplicative);
    traceExpressionToTargetSelection(route, targetDisplay);
    return executeExpressionToTargetRoute(route,
                                           f,
                                           targetBasisId,
                                           targetDisplay,
                                           targetDisplayOrder,
                                           targetIsMultiplicative,
                                           result);
  }

// ============================================================================
// Fallback Term Pipeline
// ============================================================================
// The correctness baseline normalizes, classifies, expands products, and converts terms.

bool SymmetricEngineRing::normalizeTermForConversion(
    const SymmetricTerm& term,
    VECTOR(SymmetricTerm)& normalizedTerms) const
{
    normalizedTerms.clear();
    ring_elem straightened = straightenMonomial(term.monomial);
    if (error()) return false;

    const auto *poly = polyValue(straightened);
    normalizedTerms.reserve(poly->terms.size());
    for (const auto& straightenedTerm : poly->terms)
      {
        ring_elem coeff = coefficientRing->mult(term.coeff,
                                                straightenedTerm.coeff);
        if (coefficientRing->is_zero(coeff)) continue;
        normalizedTerms.push_back({coeff, straightenedTerm.monomial});
      }
    return true;
  }

SymmetricEngineRing::TermConversionClassification
SymmetricEngineRing::classifyNormalizedTerm(
        const SymmetricTerm& term,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder,
        bool targetIsMultiplicative) const
{
    int factorCount = 0;
    std::vector<BasisKind> factorKinds;
    bool hasSkewFactor = false;
    size_t pos = 0;
    while (pos < term.monomial.data.size())
      {
        ++factorCount;
        factorKinds.push_back(
            basisKindForId(atomBasisIdAt(term.monomial, pos)));
        hasSkewFactor = hasSkewFactor || atomIsSkewAt(term.monomial, pos);
        pos += atomLengthAt(term.monomial, pos);
      }

    return TermConversionClassification{
        factorCount == 0,
        factorCount == 1,
        factorCount > 1,
        targetIsMultiplicative,
        factorCount,
        singleBasisIdInMonomial(term.monomial),
        basisKindForId(targetBasisId),
        std::move(factorKinds),
        hasSkewFactor};
  }

SymmetricEngineRing::ProductExpansionMethod
SymmetricEngineRing::selectProductExpansionMethod(
    const TermConversionClassification& classification) const
{
    if (classification.targetKind == BasisKind::Schur)
      {
        bool compatible = true;
        bool hasSchur = false;
        bool hasComplete = false;
        bool hasElementary = false;
        bool hasPowerSum = false;
        for (BasisKind kind : classification.factorKinds)
          {
            hasSchur = hasSchur || kind == BasisKind::Schur;
            hasComplete = hasComplete || kind == BasisKind::Complete;
            hasElementary = hasElementary || kind == BasisKind::Elementary;
            hasPowerSum = hasPowerSum || kind == BasisKind::PowerSum;
            if (kind != BasisKind::Schur && kind != BasisKind::Complete &&
                kind != BasisKind::Elementary && kind != BasisKind::PowerSum)
              compatible = false;
          }
        if (classification.hasSkewFactor && !hasSchur) compatible = false;
        if (!compatible) return ProductExpansionMethod::NoApplicableMethod;

        int specializedKinds = static_cast<int>(hasComplete) +
                               static_cast<int>(hasElementary) +
                               static_cast<int>(hasPowerSum);
        if (classification.hasSkewFactor || specializedKinds > 1)
          return ProductExpansionMethod::ViaSchurCompatibleRules;
        if (hasPowerSum) return ProductExpansionMethod::ViaBorderStrips;
        if (hasElementary) return ProductExpansionMethod::ViaVerticalPieri;
        if (hasComplete) return ProductExpansionMethod::ViaHorizontalPieri;
        return ProductExpansionMethod::ViaLittlewoodRichardson;
      }

    if (classification.targetKind == BasisKind::Complete ||
        classification.targetKind == BasisKind::Elementary ||
        classification.targetIsMultiplicative)
      return ProductExpansionMethod::ViaFactorwiseConversion;

    return ProductExpansionMethod::NoApplicableMethod;
  }

const char *SymmetricEngineRing::productExpansionMethodName(
    ProductExpansionMethod method) const
{
    switch (method)
      {
        case ProductExpansionMethod::ViaLittlewoodRichardson:
          return "littlewood-richardson";
        case ProductExpansionMethod::ViaHorizontalPieri:
          return "horizontal-pieri";
        case ProductExpansionMethod::ViaVerticalPieri:
          return "vertical-pieri";
        case ProductExpansionMethod::ViaBorderStrips:
          return "border-strips";
        case ProductExpansionMethod::ViaSchurCompatibleRules:
          return "schur-compatible-rules";
        case ProductExpansionMethod::ViaFactorwiseConversion:
          return "factorwise-conversion";
        case ProductExpansionMethod::NoApplicableMethod:
          return "source-target-fallback";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceProductExpansionSelection(
        ProductExpansionMethod method,
        const TermConversionClassification& classification) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::ostringstream factors;
    for (size_t i = 0; i < classification.factorKinds.size(); ++i)
      {
        if (i != 0) factors << ",";
        factors << basisKindName(classification.factorKinds[i]);
      }
    std::fprintf(stderr,
                 "SymmetricRings product-expansion: factors=%s target=%s "
                 "method=%s\n",
                 factors.str().c_str(),
                 basisKindName(classification.targetKind),
                 productExpansionMethodName(method));
  }

bool SymmetricEngineRing::executeProductExpansionMethod(
                          ProductExpansionMethod method,
                          const SymmetricTerm& term,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          bool targetIsMultiplicative,
                          ConversionInput& result) const
{
    if (method == ProductExpansionMethod::NoApplicableMethod) return false;

    ring_elem expanded;
    if (method == ProductExpansionMethod::ViaLittlewoodRichardson)
      {
        if (!trySchurProductMonomialToSchurViaLittlewoodRichardson(term.monomial,
                                         targetBasisId,
                                         targetDisplay,
                                         targetDisplayOrder,
                                         expanded))
          return false;
      }
    else if (method == ProductExpansionMethod::ViaHorizontalPieri ||
             method == ProductExpansionMethod::ViaVerticalPieri ||
             method == ProductExpansionMethod::ViaBorderStrips ||
             method == ProductExpansionMethod::ViaSchurCompatibleRules)
      {
        if (!trySchurCompatibleMonomialToSchur(term.monomial,
                                               targetBasisId,
                                               targetDisplay,
                                               targetDisplayOrder,
                                               expanded))
          return false;
      }
    else if (method == ProductExpansionMethod::ViaFactorwiseConversion)
      {
        if (!tryMonomialToTarget(term.monomial,
                                targetBasisId,
                                targetDisplay,
                                targetDisplayOrder,
                                targetIsMultiplicative,
                                expanded))
          return false;
      }
    else
      return false;

    result.expression = scaled(term.coeff, expanded);
    result.guarantees = guaranteesAfterProductExpansion(
        result.expression, targetBasisId);
    return true;
  }

ring_elem SymmetricEngineRing::runFallbackTermPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    ring_elem f = input.expression;
    ring_elem result = zero();
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        VECTOR(SymmetricTerm) normalizedTerms;
        if (!normalizeTermForConversion(term, normalizedTerms)) return zero();

        for (const auto& normalizedTerm : normalizedTerms)
          {
            TermConversionClassification classification =
                classifyNormalizedTerm(normalizedTerm,
                                       targetBasisId,
                                       targetDisplay,
                                       targetOrder,
                                       targetIsMultiplicative);
            ProductExpansionMethod method =
                selectProductExpansionMethod(classification);
            traceProductExpansionSelection(method, classification);

            ConversionInput expanded;
            bool productExpanded = executeProductExpansionMethod(
                method,
                normalizedTerm,
                targetBasisId,
                targetDisplay,
                targetOrder,
                targetIsMultiplicative,
                expanded);
            if (error()) return zero();

            if (!productExpanded)
              {
                VECTOR(SymmetricTerm) oneTerm{normalizedTerm};
                expanded.expression = fromTermVector(oneTerm, true);
                expanded.guarantees = guaranteesAfterNormalization(
                    expanded.expression, targetBasisId);
              }

            ring_elem converted = sourceToTargetDispatch(
                expanded,
                pBasisId,
                pDisplay,
                pOrder,
                pIsMultiplicative,
                targetBasisId,
                targetDisplay,
                targetOrder,
                targetIsMultiplicative);
            if (error()) return zero();
            result = add(result, converted);
          }
      }
    ConversionGuarantees resultGuarantees =
        guaranteesAfterAddition(input.guarantees, result, targetBasisId);
    attachConversionGuarantees(
        result, resultGuarantees, input.combinatorialTags);
    return result;
  }

// ============================================================================
// Retained-Operand Pipelines
// ============================================================================
// These workflows preserve operands supplied by products or plethysm.

ring_elem SymmetricEngineRing::runFactorizedProductPipeline(
        const ConversionInput& requestInput,
        ring_elem f,
        ring_elem g,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    ConversionInput productInput = requestInput;

    ring_elem direct;
    if (tryProductToTarget(f,
                                 g,
                                 targetBasisId,
                                 targetDisplay,
                                 targetOrder,
                                 targetIsMultiplicative,
                                 direct))
      {
        ConversionGuarantees inputGuarantees =
            guaranteesForFactorizedProduct(f, g, targetBasisId);
        ConversionGuarantees outputGuarantees =
            guaranteesAfterSourceTargetConversion(
                inputGuarantees, direct, targetBasisId);
        attachConversionGuarantees(
            direct, outputGuarantees, productInput.combinatorialTags);
        return direct;
      }
    if (error()) return zero();

    ring_elem product = mult(f, g);
    if (error()) return zero();
    ConversionGuarantees productGuarantees =
        guaranteesForFactorizedProduct(f, g, targetBasisId);
    productGuarantees.factorizedProduct = KnownState::False;
    productGuarantees.termCount = polyValue(product)->terms.size();
    productGuarantees.singleTerm = polyValue(product)->terms.size() == 1
        ? KnownState::True : KnownState::False;
    productInput.expression = product;
    productInput.guarantees = strengthenConversionGuarantees(
        std::move(productGuarantees), targetBasisId);
    return toBasis(productInput,
                   pBasisId,
                   pDisplay,
                   pOrder,
                   pIsMultiplicative,
                   targetBasisId,
                   targetDisplay,
                   targetOrder,
                   targetIsMultiplicative);
  }

ring_elem SymmetricEngineRing::runPostPlethysmPipeline(
        const ConversionInput& requestInput,
        ring_elem f,
        ring_elem g,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    ring_elem specialized;
    if (trySchurPlethysmToSchurViaAdamsJacobiTrudi(f,
                                                   g,
                                                   targetBasisId,
                                                   targetDisplay,
                                                   targetOrder,
                                                   specialized))
      {
        ConversionGuarantees guarantees =
            guaranteesAfterSourceTargetConversion(
                requestInput.guarantees, specialized, targetBasisId);
        attachConversionGuarantees(
            specialized,
            guarantees,
            combinatorialTagMask(CombinatorialTag::Plethysm));
        return specialized;
      }
    if (error()) return zero();

    ring_elem result = plethysm(f, g);
    if (error()) return zero();
    ConversionInput input{
        result,
        inferConversionGuarantees(result, targetBasisId),
        combinatorialTagMask(CombinatorialTag::Plethysm)};
    return runPostPlethysmPowerSumsPipeline(input,
                                           pBasisId,
                                           pDisplay,
                                           pOrder,
                                           pIsMultiplicative,
                                           targetBasisId,
                                           targetDisplay,
                                           targetOrder,
                                           targetIsMultiplicative);
  }

ring_elem SymmetricEngineRing::runPostPlethysmPowerSumsPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const
{
    return sourceToTargetDispatch(input,
                                  pBasisId,
                                  pDisplay,
                                  pOrder,
                                  pIsMultiplicative,
                                  targetBasisId,
                                  targetDisplay,
                                  targetOrder,
                                  targetIsMultiplicative);
  }

// ============================================================================
// Basis-Coefficient Dispatch
// ============================================================================
// Targeted scalar transitions avoid constructing complete target-basis expansions.

SymmetricEngineRing::BasisCoefficientRoute
SymmetricEngineRing::selectBasisCoefficientRoute(
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay) const
{
    int sourceBasisId = singleBasisId(f);
    if (sourceBasisId == targetBasisId)
      {
        CoeffMap coefficients;
        if (coefficientsInBasisIfPossible(f, targetBasisId, coefficients))
          return BasisCoefficientRoute::AlreadyExpandedInTarget;
        return BasisCoefficientRoute::ViaFullBasisConversion;
      }
    int pBasisId = basisIdForKind(BasisKind::PowerSum);
    if (sourceBasisId != pBasisId)
      return BasisCoefficientRoute::ViaFullBasisConversion;
    switch (basisKindForId(targetBasisId))
      {
        case BasisKind::PowerSum:
          return BasisCoefficientRoute::ViaPowerSumLookup;
        case BasisKind::Schur:
          return BasisCoefficientRoute::ViaSchurCharacters;
        case BasisKind::SchurOmega:
          return BasisCoefficientRoute::ViaSchurOmegaCharacters;
        case BasisKind::Complete:
          return BasisCoefficientRoute::ViaCompleteLogarithmFormula;
        case BasisKind::Elementary:
          return BasisCoefficientRoute::ViaElementaryLogarithmFormula;
        case BasisKind::HallLittlewoodQGenerator:
        case BasisKind::HallLittlewoodBGenerator:
          return BasisCoefficientRoute::ViaHallLittlewoodGeneratorLogarithmFormula;
        case BasisKind::HallLittlewoodQ:
        case BasisKind::HallLittlewoodP:
        case BasisKind::HallLittlewoodB:
        case BasisKind::HallLittlewoodPOmega:
          return BasisCoefficientRoute::ViaHallLittlewoodCapitalGreenPolynomialDuality;
        case BasisKind::Monomial:
          return BasisCoefficientRoute::ViaMonomialTransition;
        case BasisKind::Forgotten:
          return BasisCoefficientRoute::ViaForgottenTransition;
        case BasisKind::Custom:
          return BasisCoefficientRoute::ViaFullBasisConversion;
      }
    return BasisCoefficientRoute::ViaFullBasisConversion;
  }

const char *SymmetricEngineRing::basisCoefficientRouteName(
    BasisCoefficientRoute route) const
{
    switch (route)
      {
        case BasisCoefficientRoute::AlreadyExpandedInTarget:
          return "already-expanded-in-target";
        case BasisCoefficientRoute::ViaPowerSumLookup:
          return "p:coefficient-lookup";
        case BasisCoefficientRoute::ViaSchurCharacters:
          return "p->S:characters";
        case BasisCoefficientRoute::ViaSchurOmegaCharacters:
          return "p->Somega:omega-characters";
        case BasisCoefficientRoute::ViaCompleteLogarithmFormula:
          return "p->h:logarithm-formula";
        case BasisCoefficientRoute::ViaElementaryLogarithmFormula:
          return "p->e:logarithm-formula";
        case BasisCoefficientRoute::ViaHallLittlewoodGeneratorLogarithmFormula:
          return "p->q/b:logarithm-formula";
        case BasisCoefficientRoute::ViaHallLittlewoodCapitalGreenPolynomialDuality:
          return "p->Q/P/B/Pomega:Green-polynomial-duality";
        case BasisCoefficientRoute::ViaMonomialTransition:
          return "p->m:transition";
        case BasisCoefficientRoute::ViaForgottenTransition:
          return "p->ff:transition";
        case BasisCoefficientRoute::ViaFullBasisConversion:
          return "full-basis-conversion";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceBasisCoefficientSelection(
    BasisCoefficientRoute route,
    const std::string& targetDisplay,
    const Partition& targetIndex) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::fprintf(stderr,
                 "SymmetricRings basis-coefficient: target=%s_%s route=%s\n",
                 targetDisplay.c_str(),
                 partitionKey(targetIndex).c_str(),
                 basisCoefficientRouteName(route));
  }

ring_elem SymmetricEngineRing::executeBasisCoefficientRoute(
    BasisCoefficientRoute route,
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetOrder,
    bool targetIsMultiplicative,
    const Partition& targetIndex) const
{
    auto coefficientFromMap = [&](const CoeffMap& coefficients) {
      auto found = coefficients.find(targetIndex);
      return found == coefficients.end() ? coefficientRing->zero()
                                         : found->second;
    };

    if (route == BasisCoefficientRoute::AlreadyExpandedInTarget)
      return coefficientFromMap(coefficientsInBasis(f, targetBasisId));

    int pBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return coefficientRing->zero();
    if (route == BasisCoefficientRoute::ViaFullBasisConversion)
      {
        ring_elem converted = toBasis(f, targetBasisId);
        if (error()) return coefficientRing->zero();
        return coefficientFromMap(coefficientsInBasis(converted, targetBasisId));
      }

    CoeffMap powerSumCoefficients = coefficientsInBasis(f, pBasisId);
    if (error()) return coefficientRing->zero();
    if (route == BasisCoefficientRoute::ViaPowerSumLookup)
      return coefficientFromMap(powerSumCoefficients);

    ring_elem result = coefficientRing->zero();
    const BasisKind targetKind = basisKindForId(targetBasisId);
    for (const auto& term : powerSumCoefficients)
      {
        ring_elem transition = coefficientRing->zero();
        if (route == BasisCoefficientRoute::ViaSchurCharacters ||
            route == BasisCoefficientRoute::ViaSchurOmegaCharacters)
          {
            int value = characterValue(targetIndex, term.first);
            if (route == BasisCoefficientRoute::ViaSchurOmegaCharacters &&
                (partitionWeight(term.first) - partitionLength(term.first)) % 2 != 0)
              value = -value;
            transition = coefficientRing->from_long(value);
          }
        else if (route == BasisCoefficientRoute::ViaCompleteLogarithmFormula)
          transition = coefficientFromMap(
              powerSumIndexToCompleteMapViaLogarithmFormula(term.first));
        else if (route == BasisCoefficientRoute::ViaElementaryLogarithmFormula)
          transition = coefficientFromMap(
              powerSumIndexToElementaryMapViaLogarithmFormula(term.first));
        else if (route ==
                 BasisCoefficientRoute::ViaHallLittlewoodGeneratorLogarithmFormula)
          transition = coefficientFromMap(
              powerSumIndexToHallGeneratorMapViaLogarithmFormula(
                  term.first,
                  targetKind == BasisKind::HallLittlewoodBGenerator));
        else if (route ==
                 BasisCoefficientRoute::ViaHallLittlewoodCapitalGreenPolynomialDuality)
          {
            bool omega = targetKind == BasisKind::HallLittlewoodB ||
                         targetKind == BasisKind::HallLittlewoodPOmega;
            transition =
                powerSumIndexToHallLittlewoodCapitalCoefficientViaGreenPolynomialDuality(
                    term.first, targetIndex);
            if (omega &&
                (partitionWeight(term.first) - partitionLength(term.first)) % 2 != 0)
              transition = coefficientRing->negate(transition);
            if (targetKind == BasisKind::HallLittlewoodP ||
                targetKind == BasisKind::HallLittlewoodPOmega)
              transition = coefficientRing->mult(
                  transition, hallLittlewoodCFactor(targetIndex));
          }
        else if (route == BasisCoefficientRoute::ViaMonomialTransition ||
                 route == BasisCoefficientRoute::ViaForgottenTransition)
          {
            long value = pToMonomialCoefficient(term.first, targetIndex);
            if (route == BasisCoefficientRoute::ViaForgottenTransition &&
                (partitionWeight(term.first) - partitionLength(term.first)) % 2 != 0)
              value = -value;
            transition = coefficientRing->from_long(value);
          }
        if (error()) return coefficientRing->zero();
        result = coefficientRing->add(
            result,
            coefficientRing->mult(term.second, transition));
      }
    return result;
  }

ring_elem SymmetricEngineRing::basisCoefficientDispatch(
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetOrder,
    bool targetIsMultiplicative,
    const Partition& targetIndex) const
{
    requireBasis(targetBasisId);
    BasisCoefficientRoute route =
        selectBasisCoefficientRoute(f, targetBasisId, targetDisplay);
    traceBasisCoefficientSelection(route, targetDisplay, targetIndex);
    return executeBasisCoefficientRoute(route,
                                        f,
                                        targetBasisId,
                                        targetDisplay,
                                        targetOrder,
                                        targetIsMultiplicative,
                                        normalizePartition(targetIndex));
  }

ring_elem SymmetricEngineRing::basisCoefficient(
    ring_elem f,
    ring_elem targetBasisElement) const
{
    const auto *target = polyValue(targetBasisElement);
    if (target->terms.size() != 1 ||
        !coefficientRing->is_equal(target->terms.front().coeff,
                                   coefficientRing->one()))
      {
        ERROR("expected one basis element with coefficient one");
        return coefficientRing->zero();
      }
    const SymmetricMonomial& monomial = target->terms.front().monomial;
    if (monomial.data.empty() || atomIsSkewAt(monomial, 0) ||
        atomLengthAt(monomial, 0) != monomial.data.size())
      {
        ERROR("expected one non-skew basis element");
        return coefficientRing->zero();
      }
    Partition targetIndex = basisElementIndex(monomial, 0);
    if (!isPartitionIndex(targetIndex))
      {
        ERROR("expected a partition-indexed basis element");
        return coefficientRing->zero();
      }
    int targetBasisId = atomBasisIdAt(monomial, 0);
    return basisCoefficientDispatch(f,
                                    targetBasisId,
                                    displayForBasis(targetBasisId),
                                    atomOrderAt(monomial, 0),
                                    isMultiplicativeBasis(targetBasisId),
                                    targetIndex);
  }

// ============================================================================
// Public Conversion Entry Points
// ============================================================================
// All external conversion requests enter the pipeline system here.

ring_elem SymmetricEngineRing::conversionRequestToBasisDispatch(
                    const ConversionRequest& suppliedRequest,
                    int pBasisId,
                    const std::string& pDisplay,
                    int pOrder,
                    bool pIsMultiplicative,
                    int targetBasisId,
                    const std::string& targetDisplay,
                    int targetOrder,
                    bool targetIsMultiplicative) const
{
    requireBasis(pBasisId);
    requireBasis(targetBasisId);

    ConversionRequest request = suppliedRequest;
    request.input.guarantees = strengthenConversionGuarantees(
        std::move(request.input.guarantees), targetBasisId);
    ConversionPipeline pipeline =
        selectConversionPipeline(request,
                                 pBasisId,
                                 targetBasisId,
                                 targetDisplay,
                                 targetIsMultiplicative);
    traceConversionSelection(pipeline, request, targetBasisId, targetDisplay);
    return executeConversionPipeline(pipeline,
                                 request,
                                 pBasisId,
                                 pDisplay,
                                 pOrder,
                                 pIsMultiplicative,
                                 targetBasisId,
                                 targetDisplay,
                                 targetOrder,
                                 targetIsMultiplicative);
  }

ring_elem SymmetricEngineRing::toBasis(
                    const ConversionInput& input,
                    int pBasisId,
                    const std::string& pDisplay,
                    int pOrder,
                    bool pIsMultiplicative,
                    int targetBasisId,
                    const std::string& targetDisplay,
                    int targetOrder,
                    bool targetIsMultiplicative) const
{
    ConversionRequest request{
        ConversionRequestKind::Expression,
        input,
        std::nullopt,
        std::nullopt};
    return conversionRequestToBasisDispatch(request,
                                            pBasisId,
                                            pDisplay,
                                            pOrder,
                                            pIsMultiplicative,
                                            targetBasisId,
                                            targetDisplay,
                                            targetOrder,
                                            targetIsMultiplicative);
  }

ring_elem SymmetricEngineRing::toBasis(ring_elem f, int targetBasisId) const
{
    int pBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    const auto& powerSums = requireBasis(pBasisId);
    const auto& target = requireBasis(targetBasisId);
    CombinatorialTags combinatorialTags = polyValue(f)->combinatorialTags;
    ConversionInput input{
        f,
        inferConversionGuarantees(f, targetBasisId),
        combinatorialTags};
    return toBasis(input,
                   pBasisId,
                   powerSums.displaySymbol,
                   powerSums.displayOrder,
                   powerSums.multiplicative,
                   targetBasisId,
                   target.displaySymbol,
                   target.displayOrder,
                   target.multiplicative);
  }

ring_elem SymmetricEngineRing::productToBasisDispatch(ring_elem f,
                    ring_elem g,
                    int targetBasisId) const
{
    int pBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    const auto& powerSums = requireBasis(pBasisId);
    const auto& target = requireBasis(targetBasisId);

    ConversionInput input{
        zero(),
        guaranteesForFactorizedProduct(f, g, targetBasisId),
        selectMultiplicationTags(f, g)};
    ConversionRequest request{
        ConversionRequestKind::FactorizedProduct,
        input,
        f,
        g};
    return conversionRequestToBasisDispatch(request,
                                            pBasisId,
                                            powerSums.displaySymbol,
                                            powerSums.displayOrder,
                                            powerSums.multiplicative,
                                            targetBasisId,
                                            target.displaySymbol,
                                            target.displayOrder,
                                            target.multiplicative);
  }

} // namespace symmetric_rings
