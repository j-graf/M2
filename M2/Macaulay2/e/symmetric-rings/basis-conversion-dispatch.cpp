// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"
#include "rings/ZZ.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace symmetric_rings {

namespace {
// Character-table conversion grows too quickly to use outside tiny weights.
constexpr int groupedPowerSumCharacterMaxWeight = 4;
constexpr size_t groupedPowerSumCharacterMinTermCount = 4;
constexpr int powerSumViaCompleteFirstWeight = 24;
constexpr size_t powerSumViaCompleteMinTermCountAtFirstWeight = 100;
constexpr size_t powerSumViaCompleteMinTermCountAboveFirstWeight = 150;
}

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
    if (input.origin == SymmetricConversionOrigin::Plethysm &&
        input.guarantees.expandedBasis == pBasisId)
      return ConversionPipeline::PostPlethysmPowerSum;
    if (input.guarantees.expandedBasis == pBasisId)
      return ConversionPipeline::PowerSum;

    if (canUseWholeExpressionPipeline(input.guarantees,
                                      targetBasisId,
                                      targetDisplay,
                                      targetIsMultiplicative))
      return ConversionPipeline::WholeExpression;

    return ConversionPipeline::FallbackTerm;
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
        case ConversionPipeline::PowerSum:
          return runPowerSumPipeline(input,
                                  pBasisId,
                                  pDisplay,
                                  pOrder,
                                  pIsMultiplicative,
                                  targetBasisId,
                                  targetDisplay,
                                  targetOrder,
                                  targetIsMultiplicative);
        case ConversionPipeline::PostPlethysmPowerSum:
          return runPostPlethysmPowerSumPipeline(input,
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
          return runProductToBasisPipeline(*request.leftOperand,
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

SymmetricEngineRing::PowerSumsToSchurMethod
SymmetricEngineRing::selectPowerSumsToSchurMethod(
        const ConversionGuarantees& guarantees) const
{
    const char *forcedRoute = std::getenv(
        "M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE");
    if (forcedRoute != nullptr)
      {
        std::string route(forcedRoute);
        if (route == "border-strips")
          return PowerSumsToSchurMethod::ViaBorderStrips;
        if (route == "via-complete")
          return PowerSumsToSchurMethod::ViaComplete;
        if (route == "grouped-characters")
          return PowerSumsToSchurMethod::ViaCharacters;
      }

    if (guarantees.homogeneousWeight && guarantees.termCount)
      {
        int weight = *guarantees.homogeneousWeight;
        size_t termCount = *guarantees.termCount;
        if (weight >= 0 && weight <= groupedPowerSumCharacterMaxWeight &&
            termCount >= groupedPowerSumCharacterMinTermCount)
          return PowerSumsToSchurMethod::ViaCharacters;
        bool denseAtFirstWeight =
            weight == powerSumViaCompleteFirstWeight &&
            termCount >= powerSumViaCompleteMinTermCountAtFirstWeight;
        bool denseAboveFirstWeight =
            weight > powerSumViaCompleteFirstWeight &&
            termCount >= powerSumViaCompleteMinTermCountAboveFirstWeight;
        if (denseAtFirstWeight || denseAboveFirstWeight)
          return PowerSumsToSchurMethod::ViaComplete;
      }

    return PowerSumsToSchurMethod::ViaBorderStrips;
  }

SymmetricEngineRing::SourceToTargetMethod
SymmetricEngineRing::selectSourceToTargetMethod(
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay) const
{
    if (input.guarantees.targetClosed == KnownState::True)
      return SourceToTargetMethod::AlreadyInTarget;

    if (input.guarantees.expandedBasis == pBasisId)
      return SourceToTargetMethod::ViaPowerSumKernels;

    if (targetDisplay == "S")
      {
        int completeId = requiredBasisIdForDisplay("h");
        if (error()) return SourceToTargetMethod::ViaComplete;
        if (input.guarantees.pureBasis == completeId)
          return SourceToTargetMethod::ViaCompleteRecursiveTransition;
        return SourceToTargetMethod::ViaComplete;
      }

    return SourceToTargetMethod::ViaPowerSums;
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
    rememberBasis(pBasisId, pDisplay, pOrder, pIsMultiplicative);
    rememberBasis(targetBasisId,
                  targetDisplay,
                  targetOrder,
                  targetIsMultiplicative);
    SourceToTargetMethod method = selectSourceToTargetMethod(
        input, pBasisId, targetBasisId, targetDisplay);
    if (error()) return zero();
    traceSourceToTargetSelection(method, input, targetDisplay);
    auto finish = [&](ring_elem result) {
      if (!error())
        {
          ConversionGuarantees guarantees =
              guaranteesAfterSourceTargetConversion(
                  input.guarantees, result, targetBasisId);
          attachConversionGuarantees(result, guarantees, input.origin);
        }
      return result;
    };

    if (method == SourceToTargetMethod::AlreadyInTarget)
      return finish(copyPolyValue(polyValue(input.expression)));

    if (method == SourceToTargetMethod::ViaPowerSumKernels)
      {
        if (targetDisplay == "S")
          return finish(powerSumsToSchurDispatch(
              input, targetBasisId, targetDisplay, targetOrder));
        return finish(powerSumsToTargetDispatch(input.expression,
                                                targetBasisId,
                                                targetDisplay,
                                                targetOrder,
                                                targetIsMultiplicative));
      }

    if (method == SourceToTargetMethod::ViaCompleteRecursiveTransition)
      {
        int completeId = requiredBasisIdForDisplay("h");
        if (error()) return zero();
        return finish(completeToSchurViaRecursiveTransition(
            input.expression,
            completeId,
            "h",
            basisOrderForId(completeId),
            isMultiplicativeBasis(completeId),
            targetBasisId,
            targetDisplay,
            targetOrder));
      }

    ring_elem inPowerSums = expressionToPowerSumsDispatch(input.expression);
    if (error()) return zero();
    ConversionInput powerSumInput{
        inPowerSums,
        inferConversionGuarantees(inPowerSums, targetBasisId),
        input.origin};
    powerSumInput.guarantees.pureBasis = pBasisId;
    powerSumInput.guarantees.expandedBasis = pBasisId;
    if (!powerSumInput.guarantees.homogeneousWeight &&
        input.guarantees.homogeneousWeight)
      powerSumInput.guarantees.homogeneousWeight =
          input.guarantees.homogeneousWeight;
    powerSumInput.guarantees = strengthenConversionGuarantees(
        std::move(powerSumInput.guarantees), targetBasisId);

    if (method == SourceToTargetMethod::ViaComplete)
      {
        int completeId = requiredBasisIdForDisplay("h");
        if (error()) return zero();
        int completeOrder = basisOrderForId(completeId);
        bool completeIsMultiplicative = isMultiplicativeBasis(completeId);
        ring_elem inComplete = powerSumsToCompleteViaNewtonRecurrence(
            inPowerSums, completeId, completeOrder);
        if (error()) return zero();
        return finish(completeToSchurViaRecursiveTransition(
            inComplete,
            completeId,
            "h",
            completeOrder,
            completeIsMultiplicative,
            targetBasisId,
            targetDisplay,
            targetOrder));
      }

    if (targetDisplay == "S")
      return finish(powerSumsToSchurDispatch(
          powerSumInput, targetBasisId, targetDisplay, targetOrder));
    return finish(powerSumsToTargetDispatch(inPowerSums,
                                            targetBasisId,
                                            targetDisplay,
                                            targetOrder,
                                            targetIsMultiplicative));
  }

const char *SymmetricEngineRing::powerSumsToSchurMethodName(
    PowerSumsToSchurMethod method) const
{
    switch (method)
      {
        case PowerSumsToSchurMethod::ViaBorderStrips:
          return "p->S:border-strips";
        case PowerSumsToSchurMethod::ViaComplete:
          return "p->h->S:recursive-transition";
        case PowerSumsToSchurMethod::ViaCharacters:
          return "p->S:grouped-characters";
      }
    return "p->S:unknown";
  }

const char *SymmetricEngineRing::sourceToTargetMethodName(
    SourceToTargetMethod method) const
{
    switch (method)
      {
        case SourceToTargetMethod::AlreadyInTarget:
          return "already-in-target";
        case SourceToTargetMethod::ViaPowerSumKernels:
          return "power-sum-kernels";
        case SourceToTargetMethod::ViaCompleteRecursiveTransition:
          return "h->S:recursive-transition";
        case SourceToTargetMethod::ViaComplete:
          return "source->p->h->S";
        case SourceToTargetMethod::ViaPowerSums:
          return "source->p->target";
      }
    return "unknown";
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

const char *SymmetricEngineRing::conversionPipelineName(
        ConversionPipeline pipeline) const
{
    switch (pipeline)
      {
        case ConversionPipeline::WholeExpression: return "whole-expression";
        case ConversionPipeline::PowerSum: return "power-sum";
        case ConversionPipeline::PostPlethysmPowerSum:
          return "post-plethysm-power-sum";
        case ConversionPipeline::FallbackTerm: return "fallback-term";
        case ConversionPipeline::FactorizedProduct: return "factorized-product";
        case ConversionPipeline::PostPlethysm: return "post-plethysm";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceConversionSelection(
        ConversionPipeline pipeline,
        const ConversionRequest& request,
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
    const char *route = "pipeline-defined";
    if ((pipeline == ConversionPipeline::PowerSum ||
         pipeline == ConversionPipeline::PostPlethysmPowerSum) &&
        targetDisplay == "S")
      route = powerSumsToSchurMethodName(
          selectPowerSumsToSchurMethod(input.guarantees));

    std::fprintf(stderr,
                 "SymmetricRings conversion: pipeline=%s route=%s target=%s "
                 "pureBasis=%s expandedBasis=%s terms=%s weight=%s maxLength=%s\n",
                 conversionPipelineName(pipeline),
                 route,
                 targetDisplay.c_str(),
                 knownInt(input.guarantees.pureBasis).c_str(),
                 knownInt(input.guarantees.expandedBasis).c_str(),
                 knownSize(input.guarantees.termCount).c_str(),
                 knownInt(input.guarantees.homogeneousWeight).c_str(),
                 knownSize(input.guarantees.maximumPartitionLength).c_str());
  }

void SymmetricEngineRing::traceSourceToTargetSelection(
        SourceToTargetMethod method,
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
                 "SymmetricRings source-target: source=%s target=%s method=%s\n",
                 source.c_str(),
                 targetDisplay.c_str(),
                 sourceToTargetMethodName(method));
  }

void SymmetricEngineRing::traceProductExpansionSelection(
        ProductExpansionMethod method,
        const TermConversionClassification& classification) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::ostringstream factors;
    for (size_t i = 0; i < classification.factorDisplays.size(); ++i)
      {
        if (i != 0) factors << ",";
        factors << classification.factorDisplays[i];
      }
    std::fprintf(stderr,
                 "SymmetricRings product-expansion: factors=%s target=%s "
                 "method=%s\n",
                 factors.str().c_str(),
                 classification.targetDisplay.c_str(),
                 productExpansionMethodName(method));
  }

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
        guarantees.factorBases = metadata.factorBases;
        if (metadata.singleAtom)
          guarantees.singleAtom = *metadata.singleAtom
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
    size_t totalAtoms = 0;
    std::optional<int> homogeneousWeight;
    bool weightsAgree = true;
    std::vector<int> factorBases;

    for (const auto& term : poly->terms)
      {
        size_t termAtoms = 0;
        int termWeight = 0;
        size_t pos = 0;
        while (pos < term.monomial.data.size())
          {
            ++termAtoms;
            ++totalAtoms;
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
                Partition outer = atomOuterIndex(term.monomial, pos);
                Partition inner = atomInnerIndex(term.monomial, pos);
                termWeight += partitionWeight(outer) - partitionWeight(inner);
                normalized = normalized && isPartitionIndex(outer) &&
                             isPartitionIndex(inner);
              }
            else
              {
                Partition index = atomIndex(term.monomial, pos);
                termWeight += partitionWeight(index);
                if (!isMultiplicativeBasis(basisId))
                  normalized = normalized && isPartitionIndex(index);
              }
            pos += atomLengthAt(term.monomial, pos);
          }
        if (termAtoms > 1) noProducts = false;
        if (!homogeneousWeight)
          homogeneousWeight = termWeight;
        else if (*homogeneousWeight != termWeight)
          weightsAgree = false;
      }

    guarantees.singleTerm = poly->terms.size() == 1
        ? KnownState::True
        : KnownState::False;
    guarantees.singleAtom = poly->terms.size() == 1 && totalAtoms == 1
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

    if (poly->terms.empty() || totalAtoms == 0)
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

    if (guarantees.singleAtom == KnownState::True)
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
    guarantees.singleAtom = poly->terms.size() == 1 &&
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
        SymmetricConversionOrigin origin) const
{
    SymmetricConversionMetadata metadata;
    metadata.pureBasis = guarantees.pureBasis;
    metadata.expandedBasis = guarantees.expandedBasis;
    metadata.homogeneousWeight = guarantees.homogeneousWeight;
    metadata.termCount = guarantees.termCount;
    metadata.maximumPartitionLength = guarantees.maximumPartitionLength;
    metadata.factorBases = guarantees.factorBases;
    if (guarantees.singleAtom != KnownState::Unknown)
      metadata.singleAtom = guarantees.singleAtom == KnownState::True;
    if (guarantees.singleTerm != KnownState::Unknown)
      metadata.singleTerm = guarantees.singleTerm == KnownState::True;
    if (guarantees.noProducts != KnownState::Unknown)
      metadata.noProducts = guarantees.noProducts == KnownState::True;
    metadata.normalized = guarantees.normalized == KnownState::True;
    metadata.skewFree = guarantees.skewFree == KnownState::True;
    metadata.collected = guarantees.collected == KnownState::True;
    metadata.origin = origin;
    mutablePolyValue(f)->conversionMetadata = std::move(metadata);
  }

bool SymmetricEngineRing::canUseWholeExpressionPipeline(
        const ConversionGuarantees& guarantees,
        int targetBasisId,
        const std::string& targetDisplay,
        bool targetIsMultiplicative) const
{
    if (guarantees.targetClosed == KnownState::True) return true;
    if (!guarantees.factorBases) return false;

    auto allFactorDisplaysAre = [&](const std::vector<std::string>& displays) {
      for (int basisId : *guarantees.factorBases)
        {
          std::string display = displayForBasis(basisId);
          if (std::find(displays.begin(), displays.end(), display) == displays.end())
            return false;
        }
      return true;
    };

    bool ordinaryNormalized =
        guarantees.skewFree == KnownState::True &&
        guarantees.normalized == KnownState::True;

    if (targetDisplay == "S")
      {
        if (guarantees.pureBasis == targetBasisId) return true;
        return ordinaryNormalized &&
               allFactorDisplaysAre({"S", "h", "e", "p"});
      }

    if (targetDisplay == "Somega")
      return ordinaryNormalized && allFactorDisplaysAre({"e"});

    if (targetDisplay == "Q" || targetDisplay == "P")
      return ordinaryNormalized && allFactorDisplaysAre({"q"});

    if (targetDisplay == "B" || targetDisplay == "R")
      return ordinaryNormalized && allFactorDisplaysAre({"b"});

    if (targetDisplay != "h" && targetDisplay != "e" &&
        !targetIsMultiplicative)
      return false;

    if (!ordinaryNormalized) return false;
    return allFactorDisplaysAre({"p", "h", "e", "q", "b", "m", "ff",
                                 "S", "Somega", "Q", "B", "P", "R"});
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
    if (targetDisplay == "S")
      {
        int hId = requiredBasisIdForDisplay("h");
        if (error()) return zero();
        if (input.guarantees.pureBasis == hId)
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

    if (input.guarantees.targetClosed == KnownState::True)
      return sourceToTargetDispatch(input,
                                    pBasisId,
                                    pDisplay,
                                    pOrder,
                                    pIsMultiplicative,
                                    targetBasisId,
                                    targetDisplay,
                                    targetOrder,
                                    targetIsMultiplicative);

    ring_elem result;
    if (tryWholeExpressionToTarget(input.expression,
                                         targetBasisId,
                                         targetDisplay,
                                         targetOrder,
                                         targetIsMultiplicative,
                                         result))
      {
        ConversionGuarantees guarantees =
            guaranteesAfterSourceTargetConversion(
                input.guarantees, result, targetBasisId);
        attachConversionGuarantees(result, guarantees, input.origin);
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

ring_elem SymmetricEngineRing::runPowerSumPipeline(
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

ring_elem SymmetricEngineRing::runPostPlethysmPowerSumPipeline(
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

ring_elem SymmetricEngineRing::powerSumsToSchurDispatch(
        const ConversionInput& input,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder) const
{
    PowerSumsToSchurMethod method =
        selectPowerSumsToSchurMethod(input.guarantees);
    if (method == PowerSumsToSchurMethod::ViaBorderStrips)
      return powerSumsToSchurViaBorderStrips(input.expression,
                                             targetBasisId,
                                             targetDisplay,
                                             targetOrder);
    if (method == PowerSumsToSchurMethod::ViaCharacters)
      return powerSumsToSchurViaCharacters(input.expression,
                                           targetBasisId,
                                           targetOrder);
    if (method == PowerSumsToSchurMethod::ViaComplete)
      return powerSumsToSchurViaComplete(input.expression,
                                         targetBasisId,
                                         targetDisplay,
                                         targetOrder);

    ERROR("unknown power-sums-to-Schur conversion method");
    return zero();
  }

ring_elem SymmetricEngineRing::powerSumsToSchurViaComplete(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder) const
{
    int hId = requiredBasisIdForDisplay("h");
    if (error()) return zero();
    int hOrder = basisOrderForId(hId);
    bool hIsMultiplicative = isMultiplicativeBasis(hId);
    ring_elem inComplete = powerSumsToCompleteViaNewtonRecurrence(f,
                                                                 hId,
                                                                 hOrder);
    if (error()) return zero();
    return completeToSchurViaRecursiveTransition(inComplete,
                                                 hId,
                                                 "h",
                                                 hOrder,
                                                 hIsMultiplicative,
                                                 targetBasisId,
                                                 targetDisplay,
                                                 targetOrder);
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
    attachConversionGuarantees(result, resultGuarantees, input.origin);
    return result;
  }

bool SymmetricEngineRing::tryWholeExpressionToTarget(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative,
        ring_elem& result) const
{
    if (targetDisplay == "S")
      return trySchurCompatibleExpressionToSchur(f,
                                targetBasisId,
                                targetDisplay,
                                targetOrder,
                                result);

    return tryExpressionToTarget(f,
                                 targetBasisId,
                                 targetDisplay,
                                 targetOrder,
                                 targetIsMultiplicative,
                                 result);
  }

ring_elem SymmetricEngineRing::runProductToBasisPipeline(
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
    ring_elem direct;
    if (tryProductToTargetDispatch(f,
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
        attachConversionGuarantees(direct, outputGuarantees);
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
    ConversionInput productInput{
        product,
        strengthenConversionGuarantees(
            std::move(productGuarantees), targetBasisId)};
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
            specialized, guarantees, SymmetricConversionOrigin::Plethysm);
        return specialized;
      }
    if (error()) return zero();

    ring_elem result = plethysm(
        f, g, pBasisId, pDisplay, pOrder, pIsMultiplicative);
    if (error()) return zero();
    ConversionInput input{
        result,
        inferConversionGuarantees(result, targetBasisId),
        SymmetricConversionOrigin::Plethysm};
    return runPostPlethysmPowerSumPipeline(input,
                                           pBasisId,
                                           pDisplay,
                                           pOrder,
                                           pIsMultiplicative,
                                           targetBasisId,
                                           targetDisplay,
                                           targetOrder,
                                           targetIsMultiplicative);
  }

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
    std::vector<std::string> factorDisplays;
    bool hasSkewFactor = false;
    size_t pos = 0;
    while (pos < term.monomial.data.size())
      {
        ++factorCount;
        factorDisplays.push_back(
            displayForBasis(atomBasisIdAt(term.monomial, pos)));
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
        targetDisplay,
        std::move(factorDisplays),
        hasSkewFactor};
  }

SymmetricEngineRing::ProductExpansionMethod
SymmetricEngineRing::selectProductExpansionMethod(
    const TermConversionClassification& classification) const
{
    if (classification.targetDisplay == "S")
      {
        bool compatible = true;
        bool hasSchur = false;
        bool hasComplete = false;
        bool hasElementary = false;
        bool hasPowerSum = false;
        for (const auto& display : classification.factorDisplays)
          {
            hasSchur = hasSchur || display == "S";
            hasComplete = hasComplete || display == "h";
            hasElementary = hasElementary || display == "e";
            hasPowerSum = hasPowerSum || display == "p";
            if (display != "S" && display != "h" &&
                display != "e" && display != "p")
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

    if (classification.targetDisplay == "h" ||
        classification.targetDisplay == "e" ||
        classification.targetIsMultiplicative)
      return ProductExpansionMethod::ViaFactorwiseConversion;

    return ProductExpansionMethod::NoApplicableMethod;
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
    rememberBasis(pBasisId, pDisplay, pOrder, pIsMultiplicative);
    rememberBasis(targetBasisId, targetDisplay, targetOrder, targetIsMultiplicative);

    ConversionRequest request = suppliedRequest;
    request.input.guarantees = strengthenConversionGuarantees(
        std::move(request.input.guarantees), targetBasisId);

    ConversionPipeline pipeline =
        selectConversionPipeline(request,
                                 pBasisId,
                                 targetBasisId,
                                 targetDisplay,
                                 targetIsMultiplicative);
    traceConversionSelection(pipeline, request, targetDisplay);
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

ring_elem SymmetricEngineRing::toBasis(ring_elem f,
                    int pBasisId,
                    const std::string& pDisplay,
                    int pOrder,
                    bool pIsMultiplicative,
                    int targetBasisId,
                    const std::string& targetDisplay,
                    int targetOrder,
                    bool targetIsMultiplicative) const
{
    rememberBasis(pBasisId, pDisplay, pOrder, pIsMultiplicative);
    rememberBasis(targetBasisId, targetDisplay, targetOrder, targetIsMultiplicative);
    SymmetricConversionOrigin origin = SymmetricConversionOrigin::Unknown;
    const auto& metadata = polyValue(f)->conversionMetadata;
    if (metadata) origin = metadata->origin;
    ConversionInput input{
        f,
        inferConversionGuarantees(f, targetBasisId),
        origin};
    return toBasis(input,
                   pBasisId,
                   pDisplay,
                   pOrder,
                   pIsMultiplicative,
                   targetBasisId,
                   targetDisplay,
                   targetOrder,
                   targetIsMultiplicative);
  }

ring_elem SymmetricEngineRing::productToBasisDispatch(ring_elem f,
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
    rememberBasis(pBasisId, pDisplay, pOrder, pIsMultiplicative);
    rememberBasis(targetBasisId, targetDisplay, targetOrder, targetIsMultiplicative);

    ConversionInput input{
        zero(),
        guaranteesForFactorizedProduct(f, g, targetBasisId)};
    ConversionRequest request{
        ConversionRequestKind::FactorizedProduct,
        input,
        f,
        g};
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

} // namespace symmetric_rings
