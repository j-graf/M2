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
// Generic Normalization Workflow
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

ring_elem SymmetricEngineRing::normalizeExpression(
    ring_elem expression,
    ExpressionFacts& resultFacts,
    std::vector<size_t> *termFactorCounts) const
{
    resultFacts =
        inferExpressionFacts(
            expression, termFactorCounts);
    if (resultFacts.normalized &&
        resultFacts.skewFree &&
        resultFacts.collected)
      return expression;

    VECTOR(SymmetricTerm) normalizedTerms;
    for (const auto& term :
         polyValue(expression)->terms)
      {
        const ring_elem straightened =
            scaled(
                term.coeff,
                straightenMonomial(term.monomial));
        if (error()) return zero();
        for (const auto& straightenedTerm :
             polyValue(straightened)->terms)
          {
            bool hasSkewFactor = false;
            size_t position = 0;
            while (position <
                   straightenedTerm.monomial.data.size())
              {
                hasSkewFactor =
                    hasSkewFactor ||
                    atomIsSkewAt(
                        straightenedTerm.monomial,
                        position);
                position += atomLengthAt(
                    straightenedTerm.monomial,
                    position);
              }
            if (!hasSkewFactor)
              {
                normalizedTerms.push_back(
                    straightenedTerm);
                continue;
              }

            ring_elem expanded =
                fromCoeff(straightenedTerm.coeff);
            position = 0;
            while (position <
                   straightenedTerm.monomial.data.size())
              {
                ring_elem factor =
                    atomIsSkewAt(
                        straightenedTerm.monomial,
                        position)
                        ? expandSkewBasisElement(
                              straightenedTerm.monomial,
                              position)
                        : expressionFromBasisElement(
                              straightenedTerm.monomial,
                              position);
                if (error()) return zero();
                expanded = mult(expanded, factor);
                if (error()) return zero();
                position += atomLengthAt(
                    straightenedTerm.monomial,
                    position);
              }
            const auto *expandedPoly =
                polyValue(expanded);
            normalizedTerms.insert(
                normalizedTerms.end(),
                expandedPoly->terms.begin(),
                expandedPoly->terms.end());
          }
      }
    ring_elem result =
        fromTermVector(normalizedTerms, false);
    mutablePolyValue(result)->combinatorialTags =
        polyValue(expression)->combinatorialTags;
    resultFacts =
        inferExpressionFacts(
            result, termFactorCounts);
    if (!resultFacts.normalized ||
        !resultFacts.skewFree ||
        !resultFacts.collected)
      {
        ERROR("basis normalization did not establish its "
              "normalized, skew-free contract");
        return zero();
      }
    if (basisConversionTraceEnabled())
      std::fprintf(
          stderr,
          "SymmetricRings conversion-stage: stage=normalize "
          "input-terms=%zu output-terms=%zu products=%zu\n",
          polyValue(expression)->terms.size(),
          resultFacts.termCount,
          resultFacts.productTermCount);
    return result;
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
