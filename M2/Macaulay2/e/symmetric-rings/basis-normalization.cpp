// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"

#include <cstdio>
#include <map>
#include <set>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Declarative Basis-Normalization Rules
// ============================================================================

const std::vector<SymmetricEngineRing::BasisNormalizationRule>&
SymmetricEngineRing::basisNormalizationRuleDatabase()
{
    static const std::vector<BasisNormalizationRule> rules{
        {BasisKind::Schur,
         &SymmetricEngineRing::straightenSchurFamilyBasisElement,
         &SymmetricEngineRing::expandSkewSchurFamilyBasisElement},
        {BasisKind::SchurOmega,
         &SymmetricEngineRing::straightenSchurFamilyBasisElement,
         &SymmetricEngineRing::expandSkewSchurFamilyBasisElement},
        {BasisKind::HallLittlewoodQ,
         &SymmetricEngineRing::
             straightenHallLittlewoodCapitalBasisElement,
         &SymmetricEngineRing::
             expandSkewHallLittlewoodBasisElement},
        {BasisKind::HallLittlewoodB,
         &SymmetricEngineRing::
             straightenHallLittlewoodCapitalBasisElement,
         &SymmetricEngineRing::
             expandSkewHallLittlewoodBasisElement},
        {BasisKind::HallLittlewoodP,
         &SymmetricEngineRing::
             straightenHallLittlewoodNormalizedBasisElement,
         &SymmetricEngineRing::
             expandSkewHallLittlewoodBasisElement},
        {BasisKind::HallLittlewoodPOmega,
         &SymmetricEngineRing::
             straightenHallLittlewoodNormalizedBasisElement,
         &SymmetricEngineRing::
             expandSkewHallLittlewoodBasisElement}};
    return rules;
  }

const std::map<
    SymmetricEngineRing::BasisKind,
    const SymmetricEngineRing::BasisNormalizationRule *>&
SymmetricEngineRing::basisNormalizationRulesByKind()
{
    static const std::map<BasisKind, const BasisNormalizationRule *> index =
        [] {
          std::map<BasisKind, const BasisNormalizationRule *> result;
          for (const auto& rule :
               basisNormalizationRuleDatabase())
            result.emplace(rule.basisKind, &rule);
          return result;
        }();
    return index;
  }

void SymmetricEngineRing::validateBasisNormalizationRuleDatabase()
{
    static const bool validated = [] {
      std::set<BasisKind> kinds;
      for (const auto& rule :
           basisNormalizationRuleDatabase())
        {
          if (!kinds.insert(rule.basisKind).second)
            throw exc::engine_error(
                "duplicate basis-normalization rule");
          if (rule.straighteningFormula == nullptr &&
              rule.skewExpansionFormula == nullptr)
            throw exc::engine_error(
                "a basis-normalization rule has no formula");
        }
      if (basisNormalizationRulesByKind().size() !=
          basisNormalizationRuleDatabase().size())
        throw exc::engine_error(
            "the basis-normalization rule index is incomplete");
      return true;
    }();
    (void) validated;
  }

const SymmetricEngineRing::BasisNormalizationRule *
SymmetricEngineRing::basisNormalizationRuleFor(
    BasisKind basisKind) const
{
    validateBasisNormalizationRuleDatabase();
    const auto& rules = basisNormalizationRulesByKind();
    auto found = rules.find(basisKind);
    return found == rules.end() ? nullptr : found->second;
  }

// ============================================================================
// Mathematical Normalization Formulas
// ============================================================================

ring_elem SymmetricEngineRing::straightenSchurFamilyBasisElement(
    const SymmetricMonomial& monomial,
    size_t position) const
{
    return straightenSchurBasisElement(
        basisElementIndex(monomial, position),
        atomBasisIdAt(monomial, position));
  }

ring_elem
SymmetricEngineRing::straightenHallLittlewoodCapitalBasisElement(
    const SymmetricMonomial& monomial,
    size_t position) const
{
    return straightenHallCapitalBasisElement(
        basisElementIndex(monomial, position),
        atomBasisIdAt(monomial, position));
  }

ring_elem
SymmetricEngineRing::straightenHallLittlewoodNormalizedBasisElement(
    const SymmetricMonomial& monomial,
    size_t position) const
{
    const int normalizedBasisId =
        atomBasisIdAt(monomial, position);
    const BasisKind normalizedKind =
        basisKindForId(normalizedBasisId);
    const Partition index =
        basisElementIndex(monomial, position);
    if (isPartitionIndex(index))
      return basisElementFromIndex(
          normalizedBasisId, index);

    const BasisKind capitalKind =
        normalizedKind == BasisKind::HallLittlewoodP
            ? BasisKind::HallLittlewoodQ
            : BasisKind::HallLittlewoodB;
    const int capitalBasisId =
        requiredBasisIdForKind(capitalKind);
    if (error()) return zero();
    ring_elem straightCapital =
        straightenHallCapitalBasisElement(
            index, capitalBasisId);
    if (error()) return zero();
    ring_elem scaledCapital = scaled(
        coefficientQuotient(
            coefficientRing->one(),
            hallLittlewoodCFactor(index)),
        straightCapital);
    if (error()) return zero();
    const BasisConversionPlanId normalizationPlan{
        normalizedKind == BasisKind::HallLittlewoodP
            ? "HallLittlewoodQ->HallLittlewoodP:diagonal-scaling"
            : "HallLittlewoodB->HallLittlewoodPOmega:"
              "diagonal-scaling"};
    return executeNamedBasisConversionPlan(
        normalizationPlan, scaledCapital);
  }

ring_elem SymmetricEngineRing::expandSkewSchurFamilyBasisElement(
    const SymmetricMonomial& monomial,
    size_t position) const
{
    const BasisKind basisKind =
        basisKindForId(atomBasisIdAt(monomial, position));
    const BasisKind generatorKind =
        basisKind == BasisKind::Schur
            ? BasisKind::Complete
            : BasisKind::Elementary;
    const int generatorBasisId =
        requiredBasisIdForKind(generatorKind);
    if (error()) return zero();
    // Jacobi--Trudi naturally returns products. Product resolution is the
    // following conversion stage, so normalization deliberately keeps them.
    return jacobiTrudi(
        basisElementOuterIndex(monomial, position),
        basisElementInnerIndex(monomial, position),
        generatorBasisId);
  }

ring_elem
SymmetricEngineRing::expandSkewHallLittlewoodBasisElement(
    const SymmetricMonomial& monomial,
    size_t position) const
{
    return skewHallLittlewoodToPowerSums(
        basisElementOuterIndex(monomial, position),
        basisElementInnerIndex(monomial, position),
        basisKindForId(atomBasisIdAt(monomial, position)));
  }

// ============================================================================
// Normalization-Step Execution
// ============================================================================

ring_elem SymmetricEngineRing::straightenBasisElement(
    const SymmetricMonomial& monomial,
    size_t position) const
{
    if (atomIsSkewAt(monomial, position))
      return expressionFromBasisElement(
          monomial, position);
    const int basisId =
        atomBasisIdAt(monomial, position);
    const BasisNormalizationRule *rule =
        basisNormalizationRuleFor(
            basisKindForId(basisId));
    if (rule == nullptr ||
        rule->straighteningFormula == nullptr)
      return basisElementFromIndex(
          basisId,
          basisElementIndex(monomial, position));
    return (this->*rule->straighteningFormula)(
        monomial, position);
  }

ring_elem SymmetricEngineRing::straightenMonomial(
    const SymmetricMonomial& monomial) const
{
    ring_elem result = one();
    size_t position = 0;
    while (position < monomial.data.size())
      {
        ring_elem factor =
            straightenBasisElement(
                monomial, position);
        if (error()) return zero();
        result = mult(result, factor);
        position +=
            atomLengthAt(monomial, position);
      }
    return result;
  }

ring_elem SymmetricEngineRing::straightenElement(
    ring_elem expression) const
{
    ring_elem result = zero();
    for (const auto& term :
         polyValue(expression)->terms)
      {
        ring_elem straightened =
            straightenMonomial(term.monomial);
        if (error()) return zero();
        result = add(
            result,
            scaled(term.coeff, straightened));
      }
    return result;
  }

ring_elem SymmetricEngineRing::expandSkewBasisElement(
    const SymmetricMonomial& monomial,
    size_t position) const
{
    if (!atomIsSkewAt(monomial, position))
      {
        ERROR("skew expansion requires a skew basis element");
        return zero();
      }
    const BasisNormalizationRule *rule =
        basisNormalizationRuleFor(
            basisKindForId(
                atomBasisIdAt(monomial, position)));
    if (rule == nullptr ||
        rule->skewExpansionFormula == nullptr)
      {
        ERROR("no skew-expansion formula is registered for this basis");
        return zero();
      }
    return (this->*rule->skewExpansionFormula)(
        monomial, position);
  }

ring_elem SymmetricEngineRing::straighten(
    ring_elem expression) const
{
    return straightenElement(expression);
  }

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
