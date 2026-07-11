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
#include <limits>
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
      return ConversionPipeline::PostPlethysmPowerSums;
    if (input.guarantees.expandedBasis == pBasisId)
      return ConversionPipeline::PowerSums;

    if (targetIsMultiplicative &&
        input.guarantees.normalized == KnownState::True &&
        input.guarantees.skewFree == KnownState::True)
      return ConversionPipeline::GroupedMultiplicativeTarget;

    if (canUseGroupedHallLittlewoodPipeline(input.guarantees, targetDisplay))
      return ConversionPipeline::GroupedHallLittlewood;

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
          return runFactorizedProductPipeline(*request.leftOperand,
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

SymmetricEngineRing::PowerSumsToTargetRoute
SymmetricEngineRing::selectPowerSumsToTargetRoute(
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay) const
{
    auto schurRoute = [&](bool omega) {
      const char *forcedRoute = std::getenv(
          "M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE");
      if (forcedRoute != nullptr)
        {
          std::string route(forcedRoute);
          if (route == "border-strips")
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips
                : PowerSumsToTargetRoute::ViaSchurBorderStrips;
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
      if (guarantees.homogeneousWeight && guarantees.termCount)
        {
          int weight = *guarantees.homogeneousWeight;
          size_t termCount = *guarantees.termCount;
          if (weight >= 0 && weight <= groupedPowerSumCharacterMaxWeight &&
              termCount >= groupedPowerSumCharacterMinTermCount)
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurCharacters
                : PowerSumsToTargetRoute::ViaSchurCharacters;
          bool denseAtFirstWeight =
              weight == powerSumViaCompleteFirstWeight &&
              termCount >= powerSumViaCompleteMinTermCountAtFirstWeight;
          bool denseAboveFirstWeight =
              weight > powerSumViaCompleteFirstWeight &&
              termCount >= powerSumViaCompleteMinTermCountAboveFirstWeight;
          if (denseAtFirstWeight || denseAboveFirstWeight)
            return omega
                ? PowerSumsToTargetRoute::ViaOmegaThenSchurComplete
                : PowerSumsToTargetRoute::ViaSchurComplete;
        }
      return omega
          ? PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips
          : PowerSumsToTargetRoute::ViaSchurBorderStrips;
    };

    if (targetBasisId == pBasisId)
      return PowerSumsToTargetRoute::AlreadyInTarget;
    if (targetDisplay == "S")
      return schurRoute(false);
    if (targetDisplay == "Somega")
      return schurRoute(true);
    if (targetDisplay == "h")
      return PowerSumsToTargetRoute::ViaCompleteLogarithmFormula;
    if (targetDisplay == "e")
      return PowerSumsToTargetRoute::ViaElementaryLogarithmFormula;
    if (targetDisplay == "q" || targetDisplay == "b")
      return PowerSumsToTargetRoute::ViaHallLittlewoodGeneratorLogarithmFormula;
    if (targetDisplay == "Q" || targetDisplay == "P" ||
        targetDisplay == "B" || targetDisplay == "R")
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
    if (targetDisplay == "m")
      return PowerSumsToTargetRoute::ViaMonomialTransition;
    if (isForgottenDisplay(targetDisplay))
      return PowerSumsToTargetRoute::ViaForgottenTransition;
    return PowerSumsToTargetRoute::ViaTermwiseFallback;
  }

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
        const std::string sourceDisplay =
            displayForBasis(*input.guarantees.expandedBasis);
        if ((sourceDisplay == "Q" && targetDisplay == "P") ||
            (sourceDisplay == "P" && targetDisplay == "Q") ||
            (sourceDisplay == "B" && targetDisplay == "R") ||
            (sourceDisplay == "R" && targetDisplay == "B"))
          return SourceToTargetRoute::ViaHallLittlewoodNormalization;
        if ((sourceDisplay == "S" && targetDisplay == "Somega") ||
            (sourceDisplay == "Somega" && targetDisplay == "S"))
          return SourceToTargetRoute::ViaSchurOmegaConjugation;
      }

    if (input.guarantees.expandedBasis == pBasisId)
      return SourceToTargetRoute::ViaSelectedPowerSumsToTargetRoute;

    if (targetDisplay == "S")
      {
        int completeId = requiredBasisIdForDisplay("h");
        if (error()) return SourceToTargetRoute::ViaPowerSumsThenCompleteThenSchur;
        if (input.guarantees.pureBasis == completeId)
          return SourceToTargetRoute::ViaCompleteToSchurRecursiveTransition;
        return SourceToTargetRoute::ViaPowerSumsThenCompleteThenSchur;
      }

    return SourceToTargetRoute::ViaSourceToPowerSumsThenTarget;
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
          attachConversionGuarantees(result, guarantees, input.origin);
        }
      return result;
    };

    if (route == SourceToTargetRoute::AlreadyInTarget)
      return finish(copyPolyValue(polyValue(input.expression)));

    if (route == SourceToTargetRoute::ViaHallLittlewoodNormalization)
      {
        int sourceBasisId = *input.guarantees.expandedBasis;
        const std::string sourceDisplay = displayForBasis(sourceBasisId);
        bool capitalToNormalized = sourceDisplay == "Q" || sourceDisplay == "B";
        return finish(hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
            input.expression,
            sourceBasisId,
            targetBasisId,
            targetDisplay,
            targetOrder,
            capitalToNormalized));
      }

    if (route == SourceToTargetRoute::ViaSchurOmegaConjugation)
      return finish(schurOmegaConversionViaPartitionConjugation(
          input.expression,
          *input.guarantees.expandedBasis,
          targetBasisId,
          targetDisplay,
          targetOrder));

    if (route == SourceToTargetRoute::ViaSelectedPowerSumsToTargetRoute)
      {
        return finish(powerSumsToTargetDispatch(input,
                                                pBasisId,
                                                targetBasisId,
                                                targetDisplay,
                                                targetOrder,
                                                targetIsMultiplicative));
      }

    if (route == SourceToTargetRoute::ViaCompleteToSchurRecursiveTransition)
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

    ring_elem inPowerSums =
        expressionToPowerSumsViaBasisElementRoutes(input.expression);
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

    if (route == SourceToTargetRoute::ViaPowerSumsThenCompleteThenSchur)
      {
        int completeId = requiredBasisIdForDisplay("h");
        if (error()) return zero();
        int completeOrder = basisOrderForId(completeId);
        bool completeIsMultiplicative = isMultiplicativeBasis(completeId);
        ring_elem inComplete = powerSumsToCompleteViaLogarithmFormula(
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

    return finish(powerSumsToTargetDispatch(powerSumInput,
                                            pBasisId,
                                            targetBasisId,
                                            targetDisplay,
                                            targetOrder,
                                            targetIsMultiplicative));
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

const char *SymmetricEngineRing::powerSumsToTargetRouteName(
    PowerSumsToTargetRoute route) const
{
    switch (route)
      {
        case PowerSumsToTargetRoute::AlreadyInTarget:
          return "p->p:identity";
        case PowerSumsToTargetRoute::ViaSchurBorderStrips:
          return "p->S:border-strips";
        case PowerSumsToTargetRoute::ViaSchurComplete:
          return "p->h->S:recursive-transition";
        case PowerSumsToTargetRoute::ViaSchurCharacters:
          return "p->S:grouped-characters";
        case PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips:
          return "p->omega(p)->S->Somega:border-strips";
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
        int pBasisId = basisIdForDisplay("p");
        int targetBasisId = basisIdForDisplay(targetDisplay);
        if (pBasisId >= 0 && targetBasisId >= 0)
          route = powerSumsToTargetRouteName(
              selectPowerSumsToTargetRoute(
                  input, pBasisId, targetBasisId, targetDisplay));
      }

    std::fprintf(stderr,
                 "SymmetricRings conversion: pipeline=%s route=%s target=%s "
                 "pureBasis=%s expandedBasis=%s terms=%s weight=%s "
                 "maxLength=%s density=%s\n",
                 conversionPipelineName(pipeline),
                 route,
                 targetDisplay.c_str(),
                 knownInt(input.guarantees.pureBasis).c_str(),
                 knownInt(input.guarantees.expandedBasis).c_str(),
                 knownSize(input.guarantees.termCount).c_str(),
                 knownInt(input.guarantees.homogeneousWeight).c_str(),
                 knownSize(input.guarantees.maximumPartitionLength).c_str(),
                 knownDouble(input.guarantees.density).c_str());
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
        SymmetricConversionOrigin origin) const
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

    if (guarantees.expandedBasis &&
        guarantees.normalized == KnownState::True &&
        guarantees.skewFree == KnownState::True)
      return true;

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

    if (targetDisplay == "Q" || targetDisplay == "P" ||
        targetDisplay == "B" || targetDisplay == "R")
      return false;

    if (targetDisplay != "h" && targetDisplay != "e" &&
        !targetIsMultiplicative)
      return false;

    if (!ordinaryNormalized) return false;
    return allFactorDisplaysAre({"p", "h", "e", "q", "b", "m", "ff",
                                 "S", "Somega", "Q", "B", "P", "R"});
  }

bool SymmetricEngineRing::canUseGroupedHallLittlewoodPipeline(
        const ConversionGuarantees& guarantees,
        const std::string& targetDisplay) const
{
    if (targetDisplay != "Q" && targetDisplay != "P" &&
        targetDisplay != "B" && targetDisplay != "R")
      return false;
    if (guarantees.normalized != KnownState::True ||
        guarantees.skewFree != KnownState::True || !guarantees.factorBases)
      return false;
    const std::string generator =
        targetDisplay == "Q" || targetDisplay == "P" ? "q" : "b";
    for (int basisId : *guarantees.factorBases)
      if (displayForBasis(basisId) != generator) return false;
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
        const std::string sourceDisplay =
            displayForBasis(*input.guarantees.expandedBasis);
        if ((sourceDisplay == "Q" && targetDisplay == "P") ||
            (sourceDisplay == "P" && targetDisplay == "Q") ||
            (sourceDisplay == "B" && targetDisplay == "R") ||
            (sourceDisplay == "R" && targetDisplay == "B"))
          return WholeExpressionRoute::ViaHallLittlewoodNormalization;
        if ((sourceDisplay == "S" && targetDisplay == "Somega") ||
            (sourceDisplay == "Somega" && targetDisplay == "S"))
          return WholeExpressionRoute::ViaSchurOmegaConjugation;
      }
    if (targetDisplay == "S")
      {
        int completeId = requiredBasisIdForDisplay("h");
        if (error()) return WholeExpressionRoute::NoApplicableRoute;
        if (input.guarantees.pureBasis == completeId)
          return WholeExpressionRoute::ViaCompleteToSchurRecursiveTransition;
        return WholeExpressionRoute::ViaSchurCompatibleProducts;
      }
    if (targetDisplay == "Somega")
      return WholeExpressionRoute::ViaSchurTriangularReduction;
    if (targetDisplay == "Q" || targetDisplay == "B" ||
        targetDisplay == "P" || targetDisplay == "R")
      return WholeExpressionRoute::ViaHallLittlewoodTriangularReduction;
    (void) targetBasisId;
    return WholeExpressionRoute::NoApplicableRoute;
  }

bool SymmetricEngineRing::executeWholeExpressionRoute(
        WholeExpressionRoute route,
        const ConversionInput& input,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        ring_elem& result) const
{
    if (route == WholeExpressionRoute::AlreadyInTarget)
      {
        result = copyPolyValue(polyValue(input.expression));
        return true;
      }
    if (route == WholeExpressionRoute::ViaCompleteToSchurRecursiveTransition)
      {
        int completeId = requiredBasisIdForDisplay("h");
        if (error()) return false;
        result = completeToSchurViaRecursiveTransition(
            input.expression,
            completeId,
            "h",
            basisOrderForId(completeId),
            isMultiplicativeBasis(completeId),
            targetBasisId,
            targetDisplay,
            targetOrder);
        return !error();
      }
    if (route == WholeExpressionRoute::ViaSchurCompatibleProducts)
      return trySchurCompatibleExpressionToSchur(input.expression,
                                                  targetBasisId,
                                                  targetDisplay,
                                                  targetOrder,
                                                  result);
    if (route == WholeExpressionRoute::ViaSchurTriangularReduction)
      return tryExpressionToSchurViaTriangularReduction(input.expression,
                                                        targetBasisId,
                                                        targetDisplay,
                                                        targetOrder,
                                                        result);
    if (route == WholeExpressionRoute::ViaHallLittlewoodNormalization)
      {
        int sourceBasisId = *input.guarantees.expandedBasis;
        const std::string sourceDisplay = displayForBasis(sourceBasisId);
        result = hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
            input.expression,
            sourceBasisId,
            targetBasisId,
            targetDisplay,
            targetOrder,
            sourceDisplay == "Q" || sourceDisplay == "B");
        return !error();
      }
    if (route == WholeExpressionRoute::ViaSchurOmegaConjugation)
      {
        result = schurOmegaConversionViaPartitionConjugation(
            input.expression,
            *input.guarantees.expandedBasis,
            targetBasisId,
            targetDisplay,
            targetOrder);
        return !error();
      }
    if (route == WholeExpressionRoute::ViaHallLittlewoodTriangularReduction)
      return tryExpressionToHallLittlewoodViaTriangularReduction(
          input.expression,
          targetBasisId,
          targetDisplay,
          targetOrder,
          result);
    return false;
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

ring_elem SymmetricEngineRing::executePowerSumsToTargetRoute(
        PowerSumsToTargetRoute route,
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder,
        bool targetIsMultiplicative) const
{
    if (route == PowerSumsToTargetRoute::AlreadyInTarget)
      return copyPolyValue(polyValue(input.expression));
    if (route == PowerSumsToTargetRoute::ViaSchurBorderStrips)
      return powerSumsToSchurViaBorderStrips(input.expression,
                                             targetBasisId,
                                             targetDisplay,
                                             targetDisplayOrder);
    if (route == PowerSumsToTargetRoute::ViaSchurCharacters)
      return powerSumsToSchurLikeViaCharacters(input.expression,
                                               targetBasisId,
                                               targetDisplayOrder,
                                               targetDisplay,
                                               false);
    if (route == PowerSumsToTargetRoute::ViaSchurComplete)
      return powerSumsToSchurViaComplete(input.expression,
                                         targetBasisId,
                                         targetDisplay,
                                         targetDisplayOrder);
    if (route == PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips ||
        route == PowerSumsToTargetRoute::ViaOmegaThenSchurCharacters ||
        route == PowerSumsToTargetRoute::ViaOmegaThenSchurComplete)
      {
        ring_elem omegaInputExpression = omegaPowerSums(input.expression);
        if (error()) return zero();
        int schurId = requiredBasisIdForDisplay("S");
        if (error()) return zero();
        int schurOrder = basisOrderForId(schurId);
        ring_elem inSchur;
        if (route == PowerSumsToTargetRoute::ViaOmegaThenSchurBorderStrips)
          inSchur = powerSumsToSchurViaBorderStrips(
              omegaInputExpression, schurId, "S", schurOrder);
        else if (route == PowerSumsToTargetRoute::ViaOmegaThenSchurCharacters)
          inSchur = powerSumsToSchurLikeViaCharacters(
              omegaInputExpression, schurId, schurOrder, "S", false);
        else
          inSchur = powerSumsToSchurViaComplete(
              omegaInputExpression, schurId, "S", schurOrder);
        if (error()) return zero();
        return replaceSingleBasis(inSchur, schurId, targetDisplay);
      }
    if (route == PowerSumsToTargetRoute::ViaCompleteLogarithmFormula)
      return powerSumsToCompleteViaLogarithmFormula(
          input.expression, targetBasisId, targetDisplayOrder);
    if (route == PowerSumsToTargetRoute::ViaElementaryLogarithmFormula)
      return powerSumsToElementaryViaLogarithmFormula(
          input.expression, targetBasisId, targetDisplayOrder);
    if (route == PowerSumsToTargetRoute::ViaHallLittlewoodGeneratorLogarithmFormula)
      {
        CoeffMap generators = powerSumsToHallGeneratorMapViaLogarithmFormula(
            input.expression, targetDisplay == "b");
        if (error()) return zero();
        return coeffMapToElement(generators,
                                 targetBasisId,
                                 targetDisplay,
                                 targetDisplayOrder,
                                 true);
      }
    if (route == PowerSumsToTargetRoute::ViaHallLittlewoodSingleCycleGreenPolynomials)
      return powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials(
          input.expression,
          targetBasisId,
          targetDisplay,
          targetDisplayOrder);
    if (route == PowerSumsToTargetRoute::ViaHallLittlewoodGreenPolynomialsAndDuality)
      return powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality(
          input.expression,
          targetBasisId,
          targetDisplay,
          targetDisplayOrder);
    if (route == PowerSumsToTargetRoute::ViaHallLittlewoodTriangularReduction)
      return powerSumsToHallLittlewoodViaTriangularReduction(
          input.expression, targetBasisId, targetDisplay, targetDisplayOrder);
    if (route == PowerSumsToTargetRoute::ViaMonomialTransition ||
        route == PowerSumsToTargetRoute::ViaForgottenTransition ||
        route == PowerSumsToTargetRoute::ViaTermwiseFallback)
      return powerSumsToTargetViaTermwiseConversion(
          input.expression,
          targetBasisId,
          targetDisplay,
          targetDisplayOrder,
          targetIsMultiplicative);

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
    ring_elem inComplete = powerSumsToCompleteViaLogarithmFormula(f,
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

ring_elem SymmetricEngineRing::runFactorizedProductPipeline(
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
