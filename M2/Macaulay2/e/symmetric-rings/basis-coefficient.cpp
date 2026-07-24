// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"
#include "rings/ZZ.hpp"

#include <cstdio>
#include <cstdlib>

namespace symmetric_rings {

// ============================================================================
// Basis-Coefficient Selection
// ============================================================================
// Targeted scalar transitions avoid constructing complete target expansions.
// When no scalar kernel applies, the broad route uses the complete
// conversion to the shared basis-conversion workflow.

SymmetricEngineRing::BasisCoefficientRoute
SymmetricEngineRing::selectBasisCoefficientRoute(
    ring_elem f,
    int targetBasisId) const
{
    int sourceBasisId = singleBasisId(f);
    if (sourceBasisId == targetBasisId)
      {
        CoeffMap coefficients;
        if (coefficientsInBasisIfPossible(f, targetBasisId, coefficients))
          return BasisCoefficientRoute::AlreadyExpandedInTarget;
        return BasisCoefficientRoute::ViaFullBasisConversion;
      }
    int powerSumBasisId = basisIdForKind(BasisKind::PowerSum);
    if (sourceBasisId != powerSumBasisId)
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
          return BasisCoefficientRoute::
              ViaHallLittlewoodGeneratorLogarithmFormula;
        case BasisKind::HallLittlewoodQ:
        case BasisKind::HallLittlewoodP:
        case BasisKind::HallLittlewoodB:
        case BasisKind::HallLittlewoodPOmega:
          return BasisCoefficientRoute::
              ViaHallLittlewoodCapitalGreenPolynomialDuality;
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
        case BasisCoefficientRoute::
            ViaHallLittlewoodGeneratorLogarithmFormula:
          return "p->q/b:logarithm-formula";
        case BasisCoefficientRoute::
            ViaHallLittlewoodCapitalGreenPolynomialDuality:
          return "p->Q/P/B/Pomega:Green-polynomial-duality";
        case BasisCoefficientRoute::ViaMonomialTransition:
          return "p->m:transition-matrix";
        case BasisCoefficientRoute::ViaForgottenTransition:
          return "p->ff:transition-matrix";
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
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr)
      return;
    std::fprintf(
        stderr,
        "SymmetricRings basis-coefficient: target=%s_%s route=%s\n",
        targetDisplay.c_str(),
        partitionKey(targetIndex).c_str(),
        basisCoefficientRouteName(route));
  }

// ============================================================================
// Basis-Coefficient Execution
// ============================================================================

ring_elem SymmetricEngineRing::executeBasisCoefficientRoute(
    BasisCoefficientRoute route,
    ring_elem f,
    int targetBasisId,
    const Partition& targetIndex) const
{
    auto coefficientFromMap = [&](const CoeffMap& coefficients) {
      auto found = coefficients.find(targetIndex);
      return found == coefficients.end() ? coefficientRing->zero()
                                         : found->second;
    };

    if (route == BasisCoefficientRoute::AlreadyExpandedInTarget)
      return coefficientFromMap(coefficientsInBasis(f, targetBasisId));

    int powerSumBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return coefficientRing->zero();
    if (route == BasisCoefficientRoute::ViaFullBasisConversion)
      {
        ring_elem converted = toBasis(f, targetBasisId);
        if (error()) return coefficientRing->zero();
        return coefficientFromMap(
            coefficientsInBasis(converted, targetBasisId));
      }

    CoeffMap powerSumCoefficients =
        coefficientsInBasis(f, powerSumBasisId);
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
            mpz_class value =
                characterValueWithinLimits(targetIndex, term.first);
            if (route == BasisCoefficientRoute::ViaSchurOmegaCharacters &&
                (partitionWeight(term.first) -
                 partitionLength(term.first)) %
                        2 !=
                    0)
              value = -value;
            transition = coefficientRing->from_int(value.get_mpz_t());
          }
        else if (route ==
                 BasisCoefficientRoute::ViaCompleteLogarithmFormula)
          transition = coefficientFromMap(
              powerSumIndexToCompleteMapViaLogarithmFormula(term.first));
        else if (route ==
                 BasisCoefficientRoute::ViaElementaryLogarithmFormula)
          transition = coefficientFromMap(
              powerSumIndexToElementaryMapViaLogarithmFormula(term.first));
        else if (
            route ==
            BasisCoefficientRoute::
                ViaHallLittlewoodGeneratorLogarithmFormula)
          transition = coefficientFromMap(
              powerSumIndexToHallGeneratorMapViaLogarithmFormula(
                  term.first,
                  targetKind == BasisKind::HallLittlewoodBGenerator));
        else if (
            route ==
            BasisCoefficientRoute::
                ViaHallLittlewoodCapitalGreenPolynomialDuality)
          {
            bool omega = targetKind == BasisKind::HallLittlewoodB ||
                         targetKind == BasisKind::HallLittlewoodPOmega;
            transition =
                powerSumIndexToHallLittlewoodCapitalCoefficientViaGreenPolynomialDuality(
                    term.first, targetIndex);
            if (omega &&
                (partitionWeight(term.first) -
                 partitionLength(term.first)) %
                        2 !=
                    0)
              transition = coefficientRing->negate(transition);
            if (targetKind == BasisKind::HallLittlewoodP ||
                targetKind == BasisKind::HallLittlewoodPOmega)
              transition = coefficientRing->mult(
                  transition, hallLittlewoodCFactor(targetIndex));
          }
        else if (
            route == BasisCoefficientRoute::ViaMonomialTransition ||
            route == BasisCoefficientRoute::ViaForgottenTransition)
          {
            mpz_class value = pToMonomialCoefficientWithinLimits(
                term.first, targetIndex);
            if (route ==
                    BasisCoefficientRoute::ViaForgottenTransition &&
                (partitionWeight(term.first) -
                 partitionLength(term.first)) %
                        2 !=
                    0)
              value = -value;
            transition = coefficientRing->from_int(value.get_mpz_t());
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
    const Partition& targetIndex) const
{
    const auto& target = requireBasis(targetBasisId);
    if (error()) return coefficientRing->zero();
    BasisCoefficientRoute route =
        selectBasisCoefficientRoute(f, targetBasisId);
    Partition normalizedTargetIndex = normalizePartition(targetIndex);
    traceBasisCoefficientSelection(
        route, target.displaySymbol, normalizedTargetIndex);
    return executeBasisCoefficientRoute(
        route, f, targetBasisId, normalizedTargetIndex);
  }

// ============================================================================
// Public Basis-Coefficient Entry Point
// ============================================================================

ring_elem SymmetricEngineRing::basisCoefficient(
    ring_elem f,
    ring_elem targetBasisElement) const
{
    const auto *target = polyValue(targetBasisElement);
    if (target->terms.size() != 1 ||
        !coefficientRing->is_equal(
            target->terms.front().coeff, coefficientRing->one()))
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
    return basisCoefficientDispatch(
        f, atomBasisIdAt(monomial, 0), targetIndex);
  }

} // namespace symmetric_rings
