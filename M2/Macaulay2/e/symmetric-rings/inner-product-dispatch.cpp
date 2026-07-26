// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Inner-Product Requests And Shared Expression Facts
// ============================================================================
// The request carries the same exact facts used by conversion, multiplication,
// normalization, and stored expression caches. Pairing data is the only
// inner-product-specific request state.


SymmetricEngineRing::InnerProductRequest
SymmetricEngineRing::buildInnerProductRequest(
    ring_elem f,
    ring_elem g,
    InnerProductKind kind,
    std::map<int, InnerProductTarget> pairings) const
{
    return InnerProductRequest{
        f,
        g,
        InnerProductContext{kind,
                            powerSumPairingKind(kind),
                            std::move(pairings)},
        exactExpressionFacts(f),
        exactExpressionFacts(g)};
  }

SymmetricEngineRing::InnerProductKind
SymmetricEngineRing::innerProductKindFromCode(int kindCode) const
{
    switch (kindCode)
      {
        case 0: return InnerProductKind::OrdinaryHall;
        case 1: return InnerProductKind::HallLittlewood;
        case 2: return InnerProductKind::SchurQ;
        case 3: return InnerProductKind::Macdonald;
      }
    ERROR("unknown inner-product kind");
    return InnerProductKind::OrdinaryHall;
  }

SymmetricEngineRing::PowerSumPairingKind
SymmetricEngineRing::powerSumPairingKind(InnerProductKind kind) const
{
    switch (kind)
      {
        case InnerProductKind::OrdinaryHall:
          return PowerSumPairingKind::OrdinaryHall;
        case InnerProductKind::HallLittlewood:
          return PowerSumPairingKind::HallLittlewood;
        case InnerProductKind::SchurQ:
          return PowerSumPairingKind::SchurQ;
        case InnerProductKind::Macdonald:
          return PowerSumPairingKind::Macdonald;
      }
    return PowerSumPairingKind::OrdinaryHall;
  }

const char *SymmetricEngineRing::innerProductKindName(
    InnerProductKind kind) const
{
    switch (kind)
      {
        case InnerProductKind::OrdinaryHall: return "ordinary";
        case InnerProductKind::HallLittlewood: return "Hall-Littlewood";
        case InnerProductKind::SchurQ: return "Schur-Q";
        case InnerProductKind::Macdonald: return "Macdonald";
      }
    return "unknown";
  }


// ============================================================================
// Candidate Cost And Selection
// ============================================================================
// Routes are compared using shared exact facts and transition-cache state.

SymmetricEngineRing::InnerProductCacheState
SymmetricEngineRing::innerProductTransitionCacheState(
    const InnerProductCandidate& candidate,
    const InnerProductRequest& request) const
{
    ring_elem expansion = candidate.orientation == InnerProductOrientation::Original
        ? request.left : request.right;
    const ExpressionFacts& expansionFacts =
        candidate.orientation == InnerProductOrientation::Original
            ? request.leftFacts : request.rightFacts;
    const ExpressionFacts& probeFacts =
        candidate.orientation == InnerProductOrientation::Original
            ? request.rightFacts : request.leftFacts;
    if (candidate.route == InnerProductRoute::ViaSchurCompleteKostkaNumbers ||
        candidate.route ==
            InnerProductRoute::ViaSchurElementaryConjugateKostkaNumbers ||
        candidate.route ==
            InnerProductRoute::ViaSchurOmegaCompleteConjugateKostkaNumbers ||
        candidate.route == InnerProductRoute::ViaSchurOmegaElementaryKostkaNumbers)
      {
        Partition shape = *expansionFacts.singleBasisElementIndex;
        if (candidate.route ==
                InnerProductRoute::ViaSchurElementaryConjugateKostkaNumbers ||
            candidate.route ==
                InnerProductRoute::ViaSchurOmegaCompleteConjugateKostkaNumbers)
          shape = conjugatePartition(shape);
        return kostkaNumberCache.find(
                   std::make_pair(normalizePartition(shape),
                                  normalizePartition(
                                      *probeFacts.singleBasisElementIndex))) !=
                kostkaNumberCache.end()
            ? InnerProductCacheState::Present : InnerProductCacheState::Absent;
      }

    int pId = basisIdForKind(BasisKind::PowerSum);
    if (candidate.route != InnerProductRoute::ViaDualBasisCoefficient ||
        expansionFacts.expandedBasis != pId ||
        !probeFacts.singleBasisElementId)
      return InnerProductCacheState::Unknown;
    std::string targetDisplay = displayForBasis(candidate.coefficientBasisId);
    BasisKind targetKind = basisKindForId(candidate.coefficientBasisId);
    if (targetKind != BasisKind::HallLittlewoodQ &&
        targetKind != BasisKind::HallLittlewoodP &&
        targetKind != BasisKind::HallLittlewoodB &&
        targetKind != BasisKind::HallLittlewoodPOmega)
      return InnerProductCacheState::Unknown;

    CoeffMap powerSums;
    if (!coefficientsInBasisIfPossible(expansion, pId, powerSums))
      return InnerProductCacheState::Unknown;
    Partition target =
        normalizePartition(*probeFacts.singleBasisElementIndex);
    for (const auto& term : powerSums)
      {
        Partition cycleType = normalizePartition(term.first);
        if (hallLittlewoodPowerSumToCapitalColumnCache.find(cycleType) !=
            hallLittlewoodPowerSumToCapitalColumnCache.end())
          continue;
        std::string key = partitionKey(target) + "|" + partitionKey(cycleType);
        if (hallLittlewoodPowerSumToCapitalCoefficientCache.find(key) ==
            hallLittlewoodPowerSumToCapitalCoefficientCache.end())
          return InnerProductCacheState::Absent;
      }
    return InnerProductCacheState::Present;
  }

size_t SymmetricEngineRing::estimateInnerProductCandidateCost(
    InnerProductCandidate& candidate,
    const InnerProductRequest& request) const
{
    candidate.transitionCached =
        innerProductTransitionCacheState(candidate, request);
    const ExpressionFacts& expansionFacts =
        candidate.orientation == InnerProductOrientation::Original
            ? request.leftFacts : request.rightFacts;
    const ExpressionFacts& probeFacts =
        candidate.orientation == InnerProductOrientation::Original
            ? request.rightFacts : request.leftFacts;
    size_t expansionTerms = expansionFacts.termCount;
    size_t probeTerms = probeFacts.termCount;
    size_t weight = static_cast<size_t>(std::max(
        1, expansionFacts.homogeneousWeight.value_or(
               probeFacts.homogeneousWeight.value_or(1))));

    switch (candidate.route)
      {
        case InnerProductRoute::ViaRegisteredDiagonalPairing:
        case InnerProductRoute::ViaPowerSumDiagonalPairing:
          return std::min(expansionTerms, probeTerms);
        case InnerProductRoute::ViaSchurCompleteKostkaNumbers:
        case InnerProductRoute::ViaSchurElementaryConjugateKostkaNumbers:
        case InnerProductRoute::ViaSchurOmegaCompleteConjugateKostkaNumbers:
        case InnerProductRoute::ViaSchurOmegaElementaryKostkaNumbers:
          if (candidate.transitionCached == InnerProductCacheState::Present) return 1;
          return weight * std::max(
              expansionFacts.maximumPartitionLength,
              probeFacts.maximumPartitionLength);
        case InnerProductRoute::ViaDualBasisCoefficient:
          if (expansionFacts.expandedBasis == candidate.coefficientBasisId)
            return expansionTerms;
          if (candidate.transitionCached == InnerProductCacheState::Present)
            return expansionTerms + 1;
          if (expansionFacts.expandedBasis ==
              basisIdForKind(BasisKind::PowerSum))
            // A single dual coefficient evaluates one targeted transition
            // functional on each p-term.  It does not construct the complete
            // transition column, so weight is not an appropriate multiplier.
            return expansionTerms *
                (1 + probeFacts.maximumPartitionLength);
          return 100 + expansionTerms * weight;
        case InnerProductRoute::ViaWeightedSchurCharacters:
          return expansionTerms * probeTerms *
              std::max<size_t>(
                  1, expansionFacts.maximumPartitionLength);
        case InnerProductRoute::ViaConvertBothToPowerSums:
          return 1000 + weight * (expansionTerms + probeTerms);
      }
    return std::numeric_limits<size_t>::max();
  }

SymmetricEngineRing::InnerProductCandidate
SymmetricEngineRing::selectInnerProductCandidate(
    std::vector<InnerProductCandidate> candidates,
    const InnerProductRequest& request) const
{
    for (auto& candidate : candidates)
      candidate.estimatedCost = estimateInnerProductCandidateCost(
          candidate, request);
    const char *forced = std::getenv("M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_ROUTE");
    if (forced != nullptr)
      {
        candidates.erase(
            std::remove_if(
                candidates.begin(), candidates.end(),
                [&](const InnerProductCandidate& candidate) {
                  return std::string(forced) !=
                      innerProductRouteName(candidate.route);
                }),
            candidates.end());
        if (candidates.empty())
          {
            ERROR("forced inner-product route is not applicable");
            return InnerProductCandidate{
                InnerProductPipeline::FallbackPowerSums,
                InnerProductRoute::ViaConvertBothToPowerSums,
                InnerProductOrientation::Original};
          }
      }
    return *std::min_element(
        candidates.begin(), candidates.end(),
        [](const InnerProductCandidate& a, const InnerProductCandidate& b) {
          return std::tie(a.estimatedCost,
                          a.route,
                          a.coefficientBasisId,
                          a.orientation) <
                 std::tie(b.estimatedCost,
                          b.route,
                          b.coefficientBasisId,
                          b.orientation);
        });
  }


const char *SymmetricEngineRing::innerProductRouteName(
    InnerProductRoute route) const
{
    switch (route)
      {
        case InnerProductRoute::ViaRegisteredDiagonalPairing:
          return "registered-diagonal-pairing";
        case InnerProductRoute::ViaDualBasisCoefficient:
          return "dual-basis-coefficient";
        case InnerProductRoute::ViaSchurCompleteKostkaNumbers:
          return "Schur-complete:Kostka-numbers";
        case InnerProductRoute::ViaSchurElementaryConjugateKostkaNumbers:
          return "Schur-elementary:conjugate-Kostka-numbers";
        case InnerProductRoute::ViaSchurOmegaCompleteConjugateKostkaNumbers:
          return "Schur Omega-complete:conjugate-Kostka-numbers";
        case InnerProductRoute::ViaSchurOmegaElementaryKostkaNumbers:
          return "Schur Omega-elementary:Kostka-numbers";
        case InnerProductRoute::ViaWeightedSchurCharacters:
          return "power-sums-Schur:weighted-characters";
        case InnerProductRoute::ViaPowerSumDiagonalPairing:
          return "power-sums:diagonal-pairing";
        case InnerProductRoute::ViaConvertBothToPowerSums:
          return "convert-both-to-power-sums";
      }
    return "unknown";
  }


void SymmetricEngineRing::traceInnerProductRouteSelection(
    const InnerProductCandidate& candidate) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_INNER_PRODUCT") == nullptr) return;
    std::fprintf(
        stderr,
        "SymmetricRings inner-product: route=%s orientation=%s estimatedCost=%zu cache=%s\n",
        innerProductRouteName(candidate.route),
        candidate.orientation == InnerProductOrientation::Original
            ? "original" : "swapped",
        candidate.estimatedCost,
        candidate.transitionCached == InnerProductCacheState::Present
            ? "warm"
            : candidate.transitionCached == InnerProductCacheState::Absent ? "cold" : "unknown");
  }


ring_elem SymmetricEngineRing::executeInnerProductRoute(
    const InnerProductCandidate& candidate,
    const InnerProductRequest& request) const
{
    ring_elem expansion = candidate.orientation == InnerProductOrientation::Original
        ? request.left : request.right;
    ring_elem probe = candidate.orientation == InnerProductOrientation::Original
        ? request.right : request.left;

    switch (candidate.route)
      {
        case InnerProductRoute::ViaRegisteredDiagonalPairing:
        case InnerProductRoute::ViaPowerSumDiagonalPairing:
          {
            CoeffMap sourceCoefficients =
                coefficientsInBasis(expansion, candidate.pairingSourceBasisId);
            if (error()) return coefficientRing->zero();
            CoeffMap dualCoefficients =
                coefficientsInBasis(probe, candidate.pairingDualBasisId);
            if (error()) return coefficientRing->zero();
            return candidate.pairingKind == InnerProductPairingKind::PowerSum
                ? powerSumPairing(sourceCoefficients,
                                  dualCoefficients,
                                  request.context)
                : coefficientPairing(sourceCoefficients, dualCoefficients);
          }
        case InnerProductRoute::ViaDualBasisCoefficient:
          {
            Partition index;
            ring_elem probeCoefficient;
            if (!singleScaledBasisElement(probe,
                                          candidate.pairingDualBasisId,
                                          index,
                                          probeCoefficient))
              {
                ERROR("dual-basis coefficient route requires one basis element");
                return coefficientRing->zero();
              }
            int basisId = candidate.coefficientBasisId;
            ring_elem coefficient = basisCoefficientDispatch(
                expansion,
                basisId,
                index);
            if (error()) return coefficientRing->zero();
            return coefficientRing->mult(probeCoefficient, coefficient);
          }
        case InnerProductRoute::ViaSchurCompleteKostkaNumbers:
          return schurCompleteInnerProductViaKostkaNumbers(expansion, probe);
        case InnerProductRoute::ViaSchurElementaryConjugateKostkaNumbers:
          return schurElementaryInnerProductViaConjugateKostkaNumbers(
              expansion, probe);
        case InnerProductRoute::ViaSchurOmegaCompleteConjugateKostkaNumbers:
          return schurOmegaCompleteInnerProductViaConjugateKostkaNumbers(
              expansion, probe);
        case InnerProductRoute::ViaSchurOmegaElementaryKostkaNumbers:
          return schurOmegaElementaryInnerProductViaKostkaNumbers(
              expansion, probe);
        case InnerProductRoute::ViaWeightedSchurCharacters:
          return powerSumsSchurInnerProductViaWeightedCharacters(
              expansion, probe, request.context);
        case InnerProductRoute::ViaConvertBothToPowerSums:
          return runFallbackPowerSumsPipeline(request);
      }
    ERROR("unknown inner product route");
    return coefficientRing->zero();
  }


// ============================================================================
// Top-Level Inner-Product Pipeline
// ============================================================================
// The selector chooses one workflow; execution applies one named pipeline.

SymmetricEngineRing::InnerProductPipeline
SymmetricEngineRing::selectInnerProductPipeline(
    const InnerProductRequest& request) const
{
    const char *forced =
        std::getenv("M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_PIPELINE");
    if (forced != nullptr)
      {
        std::string value(forced);
        if (value == "fallback-power-sums" || value == "FallbackPowerSums" ||
            value == "fallback")
          return InnerProductPipeline::FallbackPowerSums;
        ERROR("only the fallback-power-sums inner-product pipeline can be forced");
        return InnerProductPipeline::FallbackPowerSums;
      }
    if (request.leftFacts.expandedBasis &&
        request.rightFacts.expandedBasis)
      for (const auto& pairing : request.context.pairings)
        if ((*request.leftFacts.expandedBasis == pairing.first &&
             *request.rightFacts.expandedBasis ==
                 pairing.second.dualBasisId) ||
            (*request.rightFacts.expandedBasis == pairing.first &&
             *request.leftFacts.expandedBasis ==
                 pairing.second.dualBasisId))
          return InnerProductPipeline::DiagonalBasis;
    if (request.leftFacts.singleBasisElement() ||
        request.rightFacts.singleBasisElement())
      return InnerProductPipeline::SingleBasisElement;
    int pId = basisIdForKind(BasisKind::PowerSum);
    if (request.leftFacts.expandedBasis == pId ||
        request.rightFacts.expandedBasis == pId)
      return InnerProductPipeline::PowerSumsStructured;
    return InnerProductPipeline::FallbackPowerSums;
  }


const char *SymmetricEngineRing::innerProductPipelineName(
    InnerProductPipeline pipeline) const
{
    switch (pipeline)
      {
        case InnerProductPipeline::DiagonalBasis: return "diagonal-basis";
        case InnerProductPipeline::SingleBasisElement: return "single-basis-element";
        case InnerProductPipeline::PowerSumsStructured: return "power-sums-structured";
        case InnerProductPipeline::FallbackPowerSums: return "fallback-power-sums";
      }
    return "unknown";
  }


void SymmetricEngineRing::traceInnerProductPipelineSelection(
    InnerProductPipeline pipeline,
    const InnerProductRequest& request) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_INNER_PRODUCT") == nullptr) return;
    auto basisName = [&](const ExpressionFacts& facts) {
      return facts.expandedBasis
          ? displayForBasis(*facts.expandedBasis)
          : facts.pureBasis ? displayForBasis(*facts.pureBasis)
                            : std::string("unknown");
    };
    auto optionalValue = [](const auto& value) {
      return value ? std::to_string(*value) : std::string("unknown");
    };
    std::fprintf(
        stderr,
        "SymmetricRings inner-product: pipeline=%s context=%s leftBasis=%s rightBasis=%s leftWeight=%s rightWeight=%s leftTerms=%s rightTerms=%s\n",
        innerProductPipelineName(pipeline),
        innerProductKindName(request.context.kind),
        basisName(request.leftFacts).c_str(),
        basisName(request.rightFacts).c_str(),
        optionalValue(request.leftFacts.homogeneousWeight).c_str(),
        optionalValue(request.rightFacts.homogeneousWeight).c_str(),
        std::to_string(request.leftFacts.termCount).c_str(),
        std::to_string(request.rightFacts.termCount).c_str());
  }


ring_elem SymmetricEngineRing::executeInnerProductPipeline(
    InnerProductPipeline pipeline,
    const InnerProductRequest& request) const
{
    switch (pipeline)
      {
        case InnerProductPipeline::DiagonalBasis:
          return runDiagonalBasisPipeline(request);
        case InnerProductPipeline::SingleBasisElement:
          return runSingleBasisElementPipeline(request);
        case InnerProductPipeline::PowerSumsStructured:
          return runPowerSumsStructuredPipeline(request);
        case InnerProductPipeline::FallbackPowerSums:
          {
            InnerProductCandidate candidate{
                InnerProductPipeline::FallbackPowerSums,
                InnerProductRoute::ViaConvertBothToPowerSums,
                InnerProductOrientation::Original};
            candidate.estimatedCost =
                estimateInnerProductCandidateCost(candidate, request);
            traceInnerProductRouteSelection(candidate);
            return runFallbackPowerSumsPipeline(request);
          }
      }
    ERROR("unknown inner-product pipeline");
    return coefficientRing->zero();
  }


ring_elem SymmetricEngineRing::hallInnerProductDispatch(
    const InnerProductRequest& request) const
{
    if (request.leftFacts.homogeneousWeight &&
        request.rightFacts.homogeneousWeight &&
        *request.leftFacts.homogeneousWeight !=
            *request.rightFacts.homogeneousWeight)
      return coefficientRing->zero();
    InnerProductPipeline pipeline = selectInnerProductPipeline(request);
    traceInnerProductPipelineSelection(pipeline, request);
    return executeInnerProductPipeline(pipeline, request);
  }


// ============================================================================
// Diagonal-Basis Pipeline
// ============================================================================
// Registered dual or diagonal expansions pair by intersecting coefficient maps.

SymmetricEngineRing::InnerProductCandidate
SymmetricEngineRing::selectDiagonalBasisRoute(
    const InnerProductRequest& request) const
{
    InnerProductCandidate fallback{
        InnerProductPipeline::FallbackPowerSums,
        InnerProductRoute::ViaConvertBothToPowerSums,
        InnerProductOrientation::Original};
    fallback.estimatedCost = std::numeric_limits<size_t>::max();
    if (!request.leftFacts.expandedBasis ||
        !request.rightFacts.expandedBasis)
      return fallback;

    std::vector<InnerProductCandidate> candidates;
    for (const auto& pairing : request.context.pairings)
      {
        if (*request.leftFacts.expandedBasis == pairing.first &&
            *request.rightFacts.expandedBasis == pairing.second.dualBasisId)
          {
            InnerProductCandidate candidate{
              InnerProductPipeline::DiagonalBasis,
              pairing.second.kind == InnerProductPairingKind::PowerSum
                  ? InnerProductRoute::ViaPowerSumDiagonalPairing
                  : InnerProductRoute::ViaRegisteredDiagonalPairing,
              InnerProductOrientation::Original};
            candidate.pairingSourceBasisId = pairing.first;
            candidate.pairingDualBasisId = pairing.second.dualBasisId;
            candidate.pairingKind = pairing.second.kind;
            candidates.push_back(candidate);
          }
        if (*request.rightFacts.expandedBasis == pairing.first &&
            *request.leftFacts.expandedBasis == pairing.second.dualBasisId)
          {
            InnerProductCandidate candidate{
              InnerProductPipeline::DiagonalBasis,
              pairing.second.kind == InnerProductPairingKind::PowerSum
                  ? InnerProductRoute::ViaPowerSumDiagonalPairing
                  : InnerProductRoute::ViaRegisteredDiagonalPairing,
              InnerProductOrientation::Swapped};
            candidate.pairingSourceBasisId = pairing.first;
            candidate.pairingDualBasisId = pairing.second.dualBasisId;
            candidate.pairingKind = pairing.second.kind;
            candidates.push_back(candidate);
          }
      }
    if (candidates.empty()) return fallback;
    return selectInnerProductCandidate(std::move(candidates), request);
  }


ring_elem SymmetricEngineRing::runDiagonalBasisPipeline(
    const InnerProductRequest& request) const
{
    InnerProductCandidate candidate = selectDiagonalBasisRoute(request);
    traceInnerProductRouteSelection(candidate);
    return executeInnerProductRoute(candidate, request);
  }


// ============================================================================
// Single-Basis-Element Pipeline
// ============================================================================
// A scalar probe enables targeted coefficients or direct scalar formulas.

SymmetricEngineRing::InnerProductCandidate
SymmetricEngineRing::selectSingleBasisElementRoute(
    const InnerProductRequest& request) const
{
    std::vector<InnerProductCandidate> candidates;
    int schurId = basisIdForKind(BasisKind::Schur);
    int completeId = basisIdForKind(BasisKind::Complete);
    int elementaryId = basisIdForKind(BasisKind::Elementary);
    int schurOmegaId = basisIdForKind(BasisKind::SchurOmega);
    int pId = basisIdForKind(BasisKind::PowerSum);

    if (request.context.kind == InnerProductKind::OrdinaryHall &&
        request.leftFacts.singleBasisElement() &&
        request.rightFacts.singleBasisElement())
      {
        auto addKostkaCandidate = [&](int leftBasisId,
                                      int rightBasisId,
                                      InnerProductRoute route,
                                      InnerProductOrientation orientation) {
          const ExpressionFacts& left =
              orientation == InnerProductOrientation::Original
                  ? request.leftFacts : request.rightFacts;
          const ExpressionFacts& right =
              orientation == InnerProductOrientation::Original
                  ? request.rightFacts : request.leftFacts;
          if (left.singleBasisElementId == leftBasisId &&
              right.singleBasisElementId == rightBasisId)
            candidates.push_back(InnerProductCandidate{
                InnerProductPipeline::SingleBasisElement,
                route,
                orientation});
        };
        for (InnerProductOrientation orientation :
             {InnerProductOrientation::Original,
              InnerProductOrientation::Swapped})
          {
            addKostkaCandidate(schurId,
                               completeId,
                               InnerProductRoute::ViaSchurCompleteKostkaNumbers,
                               orientation);
            addKostkaCandidate(
                schurId,
                elementaryId,
                InnerProductRoute::ViaSchurElementaryConjugateKostkaNumbers,
                orientation);
            addKostkaCandidate(
                schurOmegaId,
                completeId,
                InnerProductRoute::ViaSchurOmegaCompleteConjugateKostkaNumbers,
                orientation);
            addKostkaCandidate(
                schurOmegaId,
                elementaryId,
                InnerProductRoute::ViaSchurOmegaElementaryKostkaNumbers,
                orientation);
          }
      }

    auto addDualCoefficientCandidates = [&](const ExpressionFacts& probe,
                                            InnerProductOrientation orientation) {
      if (!probe.singleBasisElement() ||
          !probe.singleBasisElementId)
        return;
      for (const auto& pairing : request.context.pairings)
        {
          if (pairing.second.kind != InnerProductPairingKind::Dual ||
              pairing.second.dualBasisId !=
                  *probe.singleBasisElementId)
            continue;
          InnerProductCandidate candidate{
              InnerProductPipeline::SingleBasisElement,
              InnerProductRoute::ViaDualBasisCoefficient,
              orientation};
          candidate.coefficientBasisId = pairing.first;
          candidate.pairingSourceBasisId = pairing.first;
          candidate.pairingDualBasisId = pairing.second.dualBasisId;
          candidate.pairingKind = pairing.second.kind;
          candidates.push_back(candidate);
        }
    };
    addDualCoefficientCandidates(request.rightFacts,
                                 InnerProductOrientation::Original);
    addDualCoefficientCandidates(request.leftFacts,
                                 InnerProductOrientation::Swapped);

    if (request.leftFacts.expandedBasis == pId &&
        request.rightFacts.singleBasisElementId == schurId)
      candidates.push_back(InnerProductCandidate{
          InnerProductPipeline::SingleBasisElement,
          InnerProductRoute::ViaWeightedSchurCharacters,
          InnerProductOrientation::Original});
    if (request.rightFacts.expandedBasis == pId &&
        request.leftFacts.singleBasisElementId == schurId)
      candidates.push_back(InnerProductCandidate{
          InnerProductPipeline::SingleBasisElement,
          InnerProductRoute::ViaWeightedSchurCharacters,
          InnerProductOrientation::Swapped});

    candidates.push_back(InnerProductCandidate{
        InnerProductPipeline::FallbackPowerSums,
        InnerProductRoute::ViaConvertBothToPowerSums,
        InnerProductOrientation::Original});
    return selectInnerProductCandidate(std::move(candidates), request);
  }


ring_elem SymmetricEngineRing::runSingleBasisElementPipeline(
    const InnerProductRequest& request) const
{
    InnerProductCandidate candidate = selectSingleBasisElementRoute(request);
    traceInnerProductRouteSelection(candidate);
    return executeInnerProductRoute(candidate, request);
  }


// ============================================================================
// Structured Power-Sum Pipeline
// ============================================================================
// Power-sum coordinates use targeted functionals before the general fallback.

SymmetricEngineRing::InnerProductCandidate
SymmetricEngineRing::selectPowerSumsStructuredRoute(
    const InnerProductRequest& request) const
{
    int pId = basisIdForKind(BasisKind::PowerSum);
    int schurId = basisIdForKind(BasisKind::Schur);
    std::vector<InnerProductCandidate> candidates;
    if (request.leftFacts.expandedBasis == pId &&
        request.rightFacts.expandedBasis == schurId)
      candidates.push_back(InnerProductCandidate{
          InnerProductPipeline::PowerSumsStructured,
          InnerProductRoute::ViaWeightedSchurCharacters,
          InnerProductOrientation::Original});
    if (request.rightFacts.expandedBasis == pId &&
        request.leftFacts.expandedBasis == schurId)
      candidates.push_back(InnerProductCandidate{
          InnerProductPipeline::PowerSumsStructured,
          InnerProductRoute::ViaWeightedSchurCharacters,
          InnerProductOrientation::Swapped});
    candidates.push_back(InnerProductCandidate{
        InnerProductPipeline::FallbackPowerSums,
        InnerProductRoute::ViaConvertBothToPowerSums,
        InnerProductOrientation::Original});
    return selectInnerProductCandidate(std::move(candidates), request);
  }


ring_elem SymmetricEngineRing::runPowerSumsStructuredPipeline(
    const InnerProductRequest& request) const
{
    InnerProductCandidate candidate = selectPowerSumsStructuredRoute(request);
    traceInnerProductRouteSelection(candidate);
    return executeInnerProductRoute(candidate, request);
  }


// ============================================================================
// Fallback Power-Sum Pipeline
// ============================================================================
// The correctness baseline converts both arguments to power sums and pairs diagonally.

ring_elem SymmetricEngineRing::runFallbackPowerSumsPipeline(
    const InnerProductRequest& request) const
{
    int pId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return coefficientRing->zero();

    ring_elem fP = toBasis(request.left, pId);
    if (error()) return coefficientRing->zero();
    ring_elem gP = toBasis(request.right, pId);
    if (error()) return coefficientRing->zero();

    CoeffMap fCoeffs = coefficientsInBasis(fP, pId);
    if (error()) return coefficientRing->zero();
    CoeffMap gCoeffs = coefficientsInBasis(gP, pId);
    if (error()) return coefficientRing->zero();

    return powerSumPairing(fCoeffs, gCoeffs, request.context);
  }


// ============================================================================
// Public Inner-Product Entry Points
// ============================================================================
// External metadata is decoded and all requests enter the shared dispatcher.

std::map<int, InnerProductTarget>
SymmetricEngineRing::innerProductTargetMap(M2_arrayint innerProductMap) const
{
    std::map<int, InnerProductTarget> result;
    if (innerProductMap == nullptr) return result;
    if (innerProductMap->len % 3 != 0)
      {
        ERROR("invalid inner product metadata map");
        return result;
      }
    for (int i = 0; i < innerProductMap->len; i += 3)
      {
        int kindCode = innerProductMap->array[i + 2];
        if (kindCode != static_cast<int>(InnerProductPairingKind::Dual) &&
            kindCode != static_cast<int>(InnerProductPairingKind::PowerSum))
          {
            ERROR("invalid inner product pairing kind");
            return {};
          }
        result[innerProductMap->array[i]] = InnerProductTarget{
            innerProductMap->array[i + 1],
            static_cast<InnerProductPairingKind>(kindCode)};
      }
    return result;
  }

ring_elem SymmetricEngineRing::hallInnerProductElements(
    ring_elem f,
    ring_elem g,
    InnerProductKind kind,
    const std::map<int, InnerProductTarget>& pairings) const
{
    return hallInnerProductDispatch(
        buildInnerProductRequest(f, g, kind, pairings));
  }

ring_elem SymmetricEngineRing::hallInnerProduct(
    ring_elem f,
    ring_elem g,
    int kindCode,
    M2_arrayint innerProductMap) const
{
    InnerProductKind kind = innerProductKindFromCode(kindCode);
    if (error()) return coefficientRing->zero();
    std::map<int, InnerProductTarget> pairings =
        innerProductTargetMap(innerProductMap);
    if (error()) return coefficientRing->zero();
    return hallInnerProductDispatch(
        buildInnerProductRequest(
            f, g, kind, std::move(pairings)));
  }

} // namespace symmetric_rings
