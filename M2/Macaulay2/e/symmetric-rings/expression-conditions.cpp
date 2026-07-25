// Copyright 2026

#include "symmetric-rings/expression-conditions.hpp"

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace symmetric_rings {

// ============================================================================
// Condition Construction
// ============================================================================

namespace {

ExpressionCondition integerCondition(
    ExpressionConditionKind kind, int value)
{
  ExpressionCondition result;
  result.kind = kind;
  result.firstInteger = value;
  return result;
}

ExpressionCondition rationalCondition(
    ExpressionConditionKind kind, size_t numerator, size_t denominator)
{
  if (denominator == 0)
    throw std::invalid_argument(
        "an expression-condition denominator cannot be zero");
  if (numerator >
          static_cast<size_t>(std::numeric_limits<int>::max()) ||
      denominator >
          static_cast<size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument(
        "an expression-condition ratio exceeds the supported integer range");
  ExpressionCondition result;
  result.kind = kind;
  result.firstInteger = static_cast<int>(numerator);
  result.secondInteger = static_cast<int>(denominator);
  return result;
}

ExpressionCondition fractionWithMinimumCondition(
    ExpressionConditionKind kind,
    size_t numerator,
    size_t denominator,
    size_t minimumCount)
{
  ExpressionCondition result =
      rationalCondition(kind, numerator, denominator);
  if (minimumCount >
      static_cast<size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument(
        "an expression-condition minimum exceeds the supported "
        "integer range");
  result.thirdInteger = static_cast<int>(minimumCount);
  return result;
}

ExpressionCondition tagCondition(
    ExpressionConditionKind kind, uint32_t tagMask)
{
  ExpressionCondition result;
  result.kind = kind;
  result.tagMask = tagMask;
  return result;
}

ExpressionCondition logicalCondition(
    ExpressionConditionKind kind,
    std::vector<ExpressionCondition> operands)
{
  ExpressionCondition result;
  result.kind = kind;
  result.operands = std::move(operands);
  return result;
}

template<typename T>
const T& requiredValue(
    const std::optional<T>& value,
    const char *property,
    const ExpressionCondition& condition)
{
  if (!value)
    throw std::runtime_error(
        std::string("expression condition ") +
        expressionConditionToString(condition) +
        " requires unavailable " + property);
  return *value;
}

void requirePieceKind(
    const ExpressionPieceFacts& piece,
    ExpressionPieceKind expected,
    const ExpressionCondition& condition)
{
  if (piece.pieceKind != expected)
    throw std::runtime_error(
        std::string("expression condition ") +
        expressionConditionToString(condition) +
        " is not meaningful for " +
        expressionPieceKindName(piece.pieceKind));
}

bool containsOtherwise(const ExpressionCondition& condition)
{
  if (condition.kind == ExpressionConditionKind::Otherwise)
    return true;
  return std::any_of(
      condition.operands.begin(),
      condition.operands.end(),
      [](const ExpressionCondition& operand) {
        return containsOtherwise(operand);
      });
}

std::string logicalString(
    const ExpressionCondition& condition, const char *operation)
{
  std::ostringstream out;
  out << "(";
  for (size_t i = 0; i < condition.operands.size(); ++i)
    {
      if (i != 0) out << " " << operation << " ";
      out << expressionConditionToString(condition.operands[i]);
    }
  out << ")";
  return out.str();
}

} // namespace

ExpressionCondition always()
{
  ExpressionCondition result;
  result.kind = ExpressionConditionKind::Always;
  return result;
}

ExpressionCondition otherwise()
{
  ExpressionCondition result;
  result.kind = ExpressionConditionKind::Otherwise;
  return result;
}

ExpressionCondition termWeightAtMost(int weight)
{
  return integerCondition(
      ExpressionConditionKind::TermWeightAtMost, weight);
}

ExpressionCondition termWeightGreaterThan(int weight)
{
  return integerCondition(
      ExpressionConditionKind::TermWeightGreaterThan, weight);
}

ExpressionCondition componentWeightAtMost(int weight)
{
  return integerCondition(
      ExpressionConditionKind::ComponentWeightAtMost, weight);
}

ExpressionCondition componentWeightGreaterThan(int weight)
{
  return integerCondition(
      ExpressionConditionKind::ComponentWeightGreaterThan, weight);
}

ExpressionCondition componentTermCountAtLeast(size_t count)
{
  return rationalCondition(
      ExpressionConditionKind::ComponentTermCountAtLeast, count, 1);
}

ExpressionCondition componentDensityAtLeast(
    size_t numerator, size_t denominator)
{
  return rationalCondition(
      ExpressionConditionKind::ComponentDensityAtLeast,
      numerator,
      denominator);
}

ExpressionCondition componentSupportSquareRatioAtLeast(
    size_t numerator,
    size_t denominator)
{
  return rationalCondition(
      ExpressionConditionKind::
          ComponentSupportSquareRatioAtLeast,
      numerator,
      denominator);
}

ExpressionCondition componentAllPowerSumIndicesHaveMostlyShortCycles()
{
  return integerCondition(
      ExpressionConditionKind::
          ComponentAllPowerSumIndicesHaveMostlyShortCycles,
      0);
}

ExpressionCondition componentHasBothMostlyShortCycleAndOtherTerms()
{
  return integerCondition(
      ExpressionConditionKind::
          ComponentHasBothMostlyShortCycleAndOtherTerms,
      0);
}

ExpressionCondition componentMostlyShortCycleFractionAtLeast(
    size_t numerator,
    size_t denominator,
    size_t minimumCount)
{
  return fractionWithMinimumCondition(
      ExpressionConditionKind::
          ComponentMostlyShortCycleFractionAtLeast,
      numerator,
      denominator,
      minimumCount);
}

ExpressionCondition componentHasCommonPowerSumPartOne()
{
  return integerCondition(
      ExpressionConditionKind::
          ComponentHasCommonPowerSumPartOne,
      0);
}

ExpressionCondition componentHasCommonPowerSumPartAtMostPercentOfWeight(
    int percent)
{
  if (percent < 0)
    throw std::invalid_argument(
        "a common-part percentage cannot be negative");
  return integerCondition(
      ExpressionConditionKind::
          ComponentHasCommonPowerSumPartAtMostPercentOfWeight,
      percent);
}

ExpressionCondition indexIsHook()
{
  return integerCondition(ExpressionConditionKind::IndexIsHook, 0);
}

ExpressionCondition indexIsRectangle()
{
  return integerCondition(ExpressionConditionKind::IndexIsRectangle, 0);
}

ExpressionCondition indexIsSelfConjugate()
{
  return integerCondition(
      ExpressionConditionKind::IndexIsSelfConjugate, 0);
}

ExpressionCondition expressionHasMultipleWeights()
{
  return integerCondition(
      ExpressionConditionKind::ExpressionHasMultipleWeights, 0);
}

ExpressionCondition expressionIsSingleBasisElement()
{
  return integerCondition(
      ExpressionConditionKind::ExpressionIsSingleBasisElement, 0);
}

ExpressionCondition allPowerSumTermsAreSingleCycles()
{
  return integerCondition(
      ExpressionConditionKind::AllPowerSumTermsAreSingleCycles, 0);
}

ExpressionCondition mixedShortCyclePowerSumExpansion()
{
  return integerCondition(
      ExpressionConditionKind::MixedShortCyclePowerSumExpansion, 0);
}

ExpressionCondition powerSumIndexHasMostlyShortCycles()
{
  return integerCondition(
      ExpressionConditionKind::PowerSumIndexHasMostlyShortCycles, 0);
}

ExpressionCondition combinatorialTagsEqual(uint32_t tags)
{
  return tagCondition(
      ExpressionConditionKind::CombinatorialTagsEqual, tags);
}

ExpressionCondition hasCombinatorialTag(uint32_t tag)
{
  return tagCondition(
      ExpressionConditionKind::HasCombinatorialTag, tag);
}

ExpressionCondition coefficientRingIsQQ()
{
  return integerCondition(
      ExpressionConditionKind::CoefficientRingIsQQ, 0);
}

ExpressionCondition operator&&(
    ExpressionCondition left, ExpressionCondition right)
{
  return logicalCondition(
      ExpressionConditionKind::And,
      {std::move(left), std::move(right)});
}

ExpressionCondition operator||(
    ExpressionCondition left, ExpressionCondition right)
{
  return logicalCondition(
      ExpressionConditionKind::Or,
      {std::move(left), std::move(right)});
}

ExpressionCondition operator!(ExpressionCondition condition)
{
  return logicalCondition(
      ExpressionConditionKind::Not, {std::move(condition)});
}

// ============================================================================
// Condition Evaluation
// ============================================================================

bool expressionConditionHolds(
    const ExpressionCondition& condition,
    const ExpressionPieceFacts& piece)
{
  switch (condition.kind)
    {
      case ExpressionConditionKind::Always:
      case ExpressionConditionKind::Otherwise:
        return true;
      case ExpressionConditionKind::TermWeightAtMost:
        requirePieceKind(
            piece, ExpressionPieceKind::IndividualTerms, condition);
        return requiredValue(
                   piece.weight, "term weight", condition) <=
               condition.firstInteger;
      case ExpressionConditionKind::TermWeightGreaterThan:
        requirePieceKind(
            piece, ExpressionPieceKind::IndividualTerms, condition);
        return requiredValue(
                   piece.weight, "term weight", condition) >
               condition.firstInteger;
      case ExpressionConditionKind::ComponentWeightAtMost:
        requirePieceKind(
            piece,
            ExpressionPieceKind::HomogeneousComponents,
            condition);
        return requiredValue(
                   piece.weight, "component weight", condition) <=
               condition.firstInteger;
      case ExpressionConditionKind::ComponentWeightGreaterThan:
        requirePieceKind(
            piece,
            ExpressionPieceKind::HomogeneousComponents,
            condition);
        return requiredValue(
                   piece.weight, "component weight", condition) >
               condition.firstInteger;
      case ExpressionConditionKind::ComponentTermCountAtLeast:
        requirePieceKind(
            piece,
            ExpressionPieceKind::HomogeneousComponents,
            condition);
        return requiredValue(
                   piece.termCount, "component term count", condition) >=
               static_cast<size_t>(condition.firstInteger);
      case ExpressionConditionKind::ComponentDensityAtLeast:
        {
          requirePieceKind(
              piece,
              ExpressionPieceKind::HomogeneousComponents,
              condition);
          const size_t terms = requiredValue(
              piece.termCount, "component term count", condition);
          const size_t possible = requiredValue(
              piece.possibleTermCount,
              "component possible-term count",
              condition);
          return static_cast<long double>(terms) *
                     condition.secondInteger >=
                 static_cast<long double>(possible) *
                     condition.firstInteger;
        }
      case ExpressionConditionKind::
          ComponentSupportSquareRatioAtLeast:
        {
          requirePieceKind(
              piece,
              ExpressionPieceKind::HomogeneousComponents,
              condition);
          const size_t terms = requiredValue(
              piece.termCount, "component term count", condition);
          const size_t possible = requiredValue(
              piece.possibleTermCount,
              "component possible-term count",
              condition);
          const long double termsSquared =
              static_cast<long double>(terms) *
              static_cast<long double>(terms);
          return termsSquared * condition.secondInteger >=
                 static_cast<long double>(possible) *
                     condition.firstInteger;
        }
      case ExpressionConditionKind::
          ComponentAllPowerSumIndicesHaveMostlyShortCycles:
        {
          requirePieceKind(
              piece,
              ExpressionPieceKind::HomogeneousComponents,
              condition);
          const size_t shortCycleTerms = requiredValue(
              piece.mostlyShortCyclePowerSumTermCount,
              "mostly-short-cycle power-sum count",
              condition);
          return shortCycleTerms == requiredValue(
              piece.termCount, "component term count", condition);
        }
      case ExpressionConditionKind::
          ComponentHasBothMostlyShortCycleAndOtherTerms:
        {
          requirePieceKind(
              piece,
              ExpressionPieceKind::HomogeneousComponents,
              condition);
          const size_t shortCycleTerms = requiredValue(
              piece.mostlyShortCyclePowerSumTermCount,
              "mostly-short-cycle power-sum count",
              condition);
          const size_t terms = requiredValue(
              piece.termCount, "component term count", condition);
          return shortCycleTerms > 0 && shortCycleTerms < terms;
        }
      case ExpressionConditionKind::
          ComponentMostlyShortCycleFractionAtLeast:
        {
          requirePieceKind(
              piece,
              ExpressionPieceKind::HomogeneousComponents,
              condition);
          const size_t shortCycleTerms = requiredValue(
              piece.mostlyShortCyclePowerSumTermCount,
              "mostly-short-cycle power-sum count",
              condition);
          const size_t terms = requiredValue(
              piece.termCount, "component term count", condition);
          return shortCycleTerms >=
                     static_cast<size_t>(condition.thirdInteger) &&
                 static_cast<long double>(shortCycleTerms) *
                         condition.secondInteger >=
                     static_cast<long double>(terms) *
                         condition.firstInteger;
        }
      case ExpressionConditionKind::
          ComponentHasCommonPowerSumPartOne:
        {
          requirePieceKind(
              piece,
              ExpressionPieceKind::HomogeneousComponents,
              condition);
          const auto& parts = requiredValue(
              piece.commonPowerSumParts,
              "common power-sum parts",
              condition);
          return std::find(parts.begin(), parts.end(), 1) !=
                 parts.end();
        }
      case ExpressionConditionKind::
          ComponentHasCommonPowerSumPartAtMostPercentOfWeight:
        {
          requirePieceKind(
              piece,
              ExpressionPieceKind::HomogeneousComponents,
              condition);
          const auto& parts = requiredValue(
              piece.commonPowerSumParts,
              "common power-sum parts",
              condition);
          const int weight = requiredValue(
              piece.weight, "component weight", condition);
          return std::any_of(
              parts.begin(), parts.end(), [&](int part) {
                return 100 * static_cast<long long>(part) <=
                       static_cast<long long>(
                           condition.firstInteger) *
                           weight;
              });
        }
      case ExpressionConditionKind::IndexIsHook:
        requirePieceKind(
            piece, ExpressionPieceKind::IndividualTerms, condition);
        return isHookPartition(requiredValue(
            piece.index, "term index", condition));
      case ExpressionConditionKind::IndexIsRectangle:
        requirePieceKind(
            piece, ExpressionPieceKind::IndividualTerms, condition);
        return isRectanglePartition(requiredValue(
            piece.index, "term index", condition));
      case ExpressionConditionKind::IndexIsSelfConjugate:
        requirePieceKind(
            piece, ExpressionPieceKind::IndividualTerms, condition);
        return isSelfConjugatePartition(requiredValue(
            piece.index, "term index", condition));
      case ExpressionConditionKind::ExpressionHasMultipleWeights:
        requirePieceKind(
            piece,
            ExpressionPieceKind::WholeExpression,
            condition);
        return requiredValue(
            piece.expressionHasMultipleWeights,
            "multiple-weight fact",
            condition);
      case ExpressionConditionKind::ExpressionIsSingleBasisElement:
        requirePieceKind(
            piece,
            ExpressionPieceKind::WholeExpression,
            condition);
        return requiredValue(
            piece.singleBasisElement,
            "single-basis-element fact",
            condition);
      case ExpressionConditionKind::AllPowerSumTermsAreSingleCycles:
        requirePieceKind(
            piece,
            ExpressionPieceKind::WholeExpression,
            condition);
        return requiredValue(
            piece.allPowerSumTermsSingleCycles,
            "single-cycle power-sum fact",
            condition);
      case ExpressionConditionKind::
          MixedShortCyclePowerSumExpansion:
        {
          requirePieceKind(
              piece,
              ExpressionPieceKind::WholeExpression,
              condition);
          const size_t shortCycleTerms = requiredValue(
              piece.mostlyShortCyclePowerSumTermCount,
              "mostly-short-cycle power-sum count",
              condition);
          const size_t terms = requiredValue(
              piece.termCount, "expression term count", condition);
          return shortCycleTerms > 0 && shortCycleTerms < terms;
        }
      case ExpressionConditionKind::PowerSumIndexHasMostlyShortCycles:
        requirePieceKind(
            piece, ExpressionPieceKind::IndividualTerms, condition);
        return requiredValue(
            piece.powerSumIndexHasMostlyShortCycles,
            "mostly-short-cycle power-sum index fact",
            condition);
      case ExpressionConditionKind::CombinatorialTagsEqual:
        return requiredValue(
                   piece.combinatorialTags,
                   "combinatorial tags",
                   condition) ==
               condition.tagMask;
      case ExpressionConditionKind::HasCombinatorialTag:
        return (requiredValue(
                    piece.combinatorialTags,
                    "combinatorial tags",
                    condition) &
                condition.tagMask) != 0;
      case ExpressionConditionKind::CoefficientRingIsQQ:
        return requiredValue(
            piece.coefficientRingIsQQ,
            "coefficient-ring kind",
            condition);
      case ExpressionConditionKind::And:
        if (condition.operands.size() != 2)
          throw std::runtime_error(
              "an And expression condition requires two operands");
        return expressionConditionHolds(
                   condition.operands[0], piece) &&
               expressionConditionHolds(
                   condition.operands[1], piece);
      case ExpressionConditionKind::Or:
        if (condition.operands.size() != 2)
          throw std::runtime_error(
              "an Or expression condition requires two operands");
        return expressionConditionHolds(
                   condition.operands[0], piece) ||
               expressionConditionHolds(
                   condition.operands[1], piece);
      case ExpressionConditionKind::Not:
        if (condition.operands.size() != 1)
          throw std::runtime_error(
              "a Not expression condition requires one operand");
        return !expressionConditionHolds(
            condition.operands[0], piece);
    }
  throw std::runtime_error("unknown expression condition");
}

// ============================================================================
// Condition Contract Validation
// ============================================================================

void validateExpressionCondition(
    const ExpressionCondition& condition,
    ExpressionPieceKind pieceKind,
    bool allowOtherwise)
{
  auto requireKind =
      [&](ExpressionPieceKind expected) {
        if (pieceKind != expected)
          throw std::invalid_argument(
              "expression condition " +
              expressionConditionToString(condition) +
              " is not meaningful for " +
              expressionPieceKindName(pieceKind));
      };

  switch (condition.kind)
    {
      case ExpressionConditionKind::Always:
        if (!condition.operands.empty())
          throw std::invalid_argument(
              "always() cannot have operands");
        return;
      case ExpressionConditionKind::Otherwise:
        if (!allowOtherwise)
          throw std::invalid_argument(
              "otherwise() is allowed only as a top-level plan case");
        if (!condition.operands.empty())
          throw std::invalid_argument(
              "otherwise() cannot have operands");
        return;
      case ExpressionConditionKind::TermWeightAtMost:
      case ExpressionConditionKind::TermWeightGreaterThan:
      case ExpressionConditionKind::IndexIsHook:
      case ExpressionConditionKind::IndexIsRectangle:
      case ExpressionConditionKind::IndexIsSelfConjugate:
      case ExpressionConditionKind::PowerSumIndexHasMostlyShortCycles:
        requireKind(ExpressionPieceKind::IndividualTerms);
        break;
      case ExpressionConditionKind::ComponentWeightAtMost:
      case ExpressionConditionKind::ComponentWeightGreaterThan:
      case ExpressionConditionKind::ComponentTermCountAtLeast:
      case ExpressionConditionKind::ComponentDensityAtLeast:
      case ExpressionConditionKind::
          ComponentSupportSquareRatioAtLeast:
      case ExpressionConditionKind::
          ComponentAllPowerSumIndicesHaveMostlyShortCycles:
      case ExpressionConditionKind::
          ComponentHasBothMostlyShortCycleAndOtherTerms:
      case ExpressionConditionKind::
          ComponentMostlyShortCycleFractionAtLeast:
      case ExpressionConditionKind::
          ComponentHasCommonPowerSumPartOne:
      case ExpressionConditionKind::
          ComponentHasCommonPowerSumPartAtMostPercentOfWeight:
        requireKind(
            ExpressionPieceKind::HomogeneousComponents);
        break;
      case ExpressionConditionKind::ExpressionHasMultipleWeights:
      case ExpressionConditionKind::ExpressionIsSingleBasisElement:
      case ExpressionConditionKind::
          AllPowerSumTermsAreSingleCycles:
      case ExpressionConditionKind::
          MixedShortCyclePowerSumExpansion:
        requireKind(ExpressionPieceKind::WholeExpression);
        break;
      case ExpressionConditionKind::CombinatorialTagsEqual:
      case ExpressionConditionKind::HasCombinatorialTag:
      case ExpressionConditionKind::CoefficientRingIsQQ:
        break;
      case ExpressionConditionKind::And:
      case ExpressionConditionKind::Or:
        if (condition.operands.size() != 2)
          throw std::invalid_argument(
              "And and Or expression conditions require two operands");
        for (const auto& operand : condition.operands)
          validateExpressionCondition(
              operand, pieceKind, false);
        return;
      case ExpressionConditionKind::Not:
        if (condition.operands.size() != 1)
          throw std::invalid_argument(
              "a Not expression condition requires one operand");
        validateExpressionCondition(
            condition.operands.front(), pieceKind, false);
        return;
    }
  if (!condition.operands.empty())
    throw std::invalid_argument(
        "a leaf expression condition cannot have operands");
}

namespace {

bool sameExpressionCondition(
    const ExpressionCondition& left,
    const ExpressionCondition& right)
{
  if (left.kind != right.kind ||
      left.firstInteger != right.firstInteger ||
      left.secondInteger != right.secondInteger ||
      left.thirdInteger != right.thirdInteger ||
      left.tagMask != right.tagMask ||
      left.operands.size() != right.operands.size())
    return false;
  for (size_t i = 0; i < left.operands.size(); ++i)
    if (!sameExpressionCondition(
            left.operands[i], right.operands[i]))
      return false;
  return true;
}

bool leafConditionImplies(
    const ExpressionCondition& guarantee,
    const ExpressionCondition& requirement)
{
  if (guarantee.kind == requirement.kind)
    {
      switch (guarantee.kind)
        {
          case ExpressionConditionKind::TermWeightAtMost:
          case ExpressionConditionKind::ComponentWeightAtMost:
            return guarantee.firstInteger <=
                   requirement.firstInteger;
          case ExpressionConditionKind::TermWeightGreaterThan:
          case ExpressionConditionKind::ComponentWeightGreaterThan:
          case ExpressionConditionKind::ComponentTermCountAtLeast:
            return guarantee.firstInteger >=
                   requirement.firstInteger;
          case ExpressionConditionKind::ComponentDensityAtLeast:
            return static_cast<int64_t>(
                       guarantee.firstInteger) *
                       requirement.secondInteger >=
                   static_cast<int64_t>(
                       requirement.firstInteger) *
                       guarantee.secondInteger;
          default:
            break;
        }
    }
  if (guarantee.kind ==
          ExpressionConditionKind::CombinatorialTagsEqual &&
      requirement.kind ==
          ExpressionConditionKind::HasCombinatorialTag)
    return (guarantee.tagMask & requirement.tagMask) != 0;
  return false;
}

} // namespace

bool expressionConditionImplies(
    const ExpressionCondition& guarantee,
    const ExpressionCondition& requirement)
{
  if (requirement.kind == ExpressionConditionKind::Always ||
      sameExpressionCondition(guarantee, requirement))
    return true;

  if (requirement.kind == ExpressionConditionKind::And)
    return std::all_of(
        requirement.operands.begin(),
        requirement.operands.end(),
        [&](const ExpressionCondition& operand) {
          return expressionConditionImplies(
              guarantee, operand);
        });
  // A disjunction guarantees a requirement only when every branch does. A
  // conjunction guarantees anything already implied by either operand.
  if (guarantee.kind == ExpressionConditionKind::Or)
    return std::all_of(
        guarantee.operands.begin(),
        guarantee.operands.end(),
        [&](const ExpressionCondition& operand) {
          return expressionConditionImplies(
              operand, requirement);
        });
  if (guarantee.kind == ExpressionConditionKind::And)
    return std::any_of(
        guarantee.operands.begin(),
        guarantee.operands.end(),
        [&](const ExpressionCondition& operand) {
          return expressionConditionImplies(
              operand, requirement);
        });
  if (requirement.kind == ExpressionConditionKind::Or)
    return std::any_of(
        requirement.operands.begin(),
        requirement.operands.end(),
        [&](const ExpressionCondition& operand) {
          return expressionConditionImplies(
              guarantee, operand);
        });

  return leafConditionImplies(guarantee, requirement);
}

ExpressionFactRequirements expressionFactRequirements(
    const ExpressionCondition& condition)
{
  ExpressionFactRequirements result =
      noExpressionFactRequirements;
  switch (condition.kind)
    {
      case ExpressionConditionKind::AllPowerSumTermsAreSingleCycles:
        result |= requirePowerSumSingleCycleFacts;
        break;
      case ExpressionConditionKind::
          ComponentAllPowerSumIndicesHaveMostlyShortCycles:
      case ExpressionConditionKind::
          ComponentHasBothMostlyShortCycleAndOtherTerms:
      case ExpressionConditionKind::
          ComponentMostlyShortCycleFractionAtLeast:
      case ExpressionConditionKind::MixedShortCyclePowerSumExpansion:
      case ExpressionConditionKind::PowerSumIndexHasMostlyShortCycles:
        result |= requirePowerSumShortCycleFacts;
        break;
      case ExpressionConditionKind::
          ComponentHasCommonPowerSumPartOne:
      case ExpressionConditionKind::
          ComponentHasCommonPowerSumPartAtMostPercentOfWeight:
        result |= requirePowerSumCommonPartFacts;
        break;
      default:
        break;
    }
  for (const auto& operand : condition.operands)
    result |= expressionFactRequirements(operand);
  return result;
}

ExpressionFactRequirements expressionFactRequirements(
    const std::vector<ExpressionCondition>& conditions)
{
  ExpressionFactRequirements result =
      noExpressionFactRequirements;
  for (const auto& condition : conditions)
    result |= expressionFactRequirements(condition);
  return result;
}

// ============================================================================
// Ordered Piece Partitioning
// ============================================================================

PlanCaseAssignments assignExpressionPiecesToCases(
    const std::vector<ExpressionPieceFacts>& pieces,
    const std::vector<ExpressionCondition>& orderedConditions,
    size_t expectedTermCount)
{
  if (orderedConditions.empty())
    throw std::invalid_argument(
        "an expression partition requires at least one condition");
  for (size_t i = 0; i < orderedConditions.size(); ++i)
    {
      const bool isDefault =
          orderedConditions[i].kind ==
          ExpressionConditionKind::Otherwise;
      if (isDefault != (i + 1 == orderedConditions.size()))
        throw std::invalid_argument(
            "otherwise() must appear exactly once as the final condition");
      if (!isDefault && containsOtherwise(orderedConditions[i]))
        throw std::invalid_argument(
            "otherwise() cannot appear inside a logical condition");
    }

  PlanCaseAssignments result;
  result.termPositionsForCase.resize(orderedConditions.size());
  std::vector<bool> assigned(pieces.size(), false);
  std::vector<bool> coveredTerms(expectedTermCount, false);
  size_t assignedPieces = 0;
  size_t coveredTermCount = 0;

  for (size_t caseIndex = 0;
       caseIndex < orderedConditions.size();
       ++caseIndex)
    for (size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex)
      {
        if (assigned[pieceIndex]) continue;
        if (!expressionConditionHolds(
                orderedConditions[caseIndex], pieces[pieceIndex]))
          continue;
        assigned[pieceIndex] = true;
        ++assignedPieces;
        for (size_t termPosition : pieces[pieceIndex].termPositions)
          {
            if (termPosition >= expectedTermCount)
              throw std::runtime_error(
                  "an expression piece contains an invalid term position");
            if (coveredTerms[termPosition])
              throw std::runtime_error(
                  "expression pieces overlap in their term positions");
            coveredTerms[termPosition] = true;
            ++coveredTermCount;
            result.termPositionsForCase[caseIndex].push_back(termPosition);
          }
      }

  if (assignedPieces != pieces.size() ||
      coveredTermCount != expectedTermCount)
    throw std::runtime_error(
        "ordered expression conditions do not cover the complete input");
  return result;
}

// ============================================================================
// Condition Diagnostics
// ============================================================================

std::string expressionConditionToString(
    const ExpressionCondition& condition)
{
  switch (condition.kind)
    {
      case ExpressionConditionKind::Always:
        return "always";
      case ExpressionConditionKind::Otherwise:
        return "otherwise";
      case ExpressionConditionKind::TermWeightAtMost:
        return "term-weight <= " +
               std::to_string(condition.firstInteger);
      case ExpressionConditionKind::TermWeightGreaterThan:
        return "term-weight > " +
               std::to_string(condition.firstInteger);
      case ExpressionConditionKind::ComponentWeightAtMost:
        return "component-weight <= " +
               std::to_string(condition.firstInteger);
      case ExpressionConditionKind::ComponentWeightGreaterThan:
        return "component-weight > " +
               std::to_string(condition.firstInteger);
      case ExpressionConditionKind::ComponentTermCountAtLeast:
        return "component-terms >= " +
               std::to_string(condition.firstInteger);
      case ExpressionConditionKind::ComponentDensityAtLeast:
        return "component-density >= " +
               std::to_string(condition.firstInteger) + "/" +
               std::to_string(condition.secondInteger);
      case ExpressionConditionKind::
          ComponentSupportSquareRatioAtLeast:
        return "component-terms^2 / component-possible-terms >= " +
               std::to_string(condition.firstInteger) + "/" +
               std::to_string(condition.secondInteger);
      case ExpressionConditionKind::
          ComponentAllPowerSumIndicesHaveMostlyShortCycles:
        return "all-component-power-sum-indices-have-mostly-short-cycles";
      case ExpressionConditionKind::
          ComponentHasBothMostlyShortCycleAndOtherTerms:
        return "component-has-both-mostly-short-cycle-and-other-terms";
      case ExpressionConditionKind::
          ComponentMostlyShortCycleFractionAtLeast:
        return "component-mostly-short-cycle-terms >= max(" +
               std::to_string(condition.thirdInteger) + ", " +
               std::to_string(condition.firstInteger) + "/" +
               std::to_string(condition.secondInteger) +
               " * component-terms)";
      case ExpressionConditionKind::
          ComponentHasCommonPowerSumPartOne:
        return "component-common-power-sum-parts-contain-1";
      case ExpressionConditionKind::
          ComponentHasCommonPowerSumPartAtMostPercentOfWeight:
        return "component-has-common-power-sum-part <= " +
               std::to_string(condition.firstInteger) +
               "% of component-weight";
      case ExpressionConditionKind::IndexIsHook:
        return "index-is-hook";
      case ExpressionConditionKind::IndexIsRectangle:
        return "index-is-rectangle";
      case ExpressionConditionKind::IndexIsSelfConjugate:
        return "index-is-self-conjugate";
      case ExpressionConditionKind::ExpressionHasMultipleWeights:
        return "expression-has-multiple-weights";
      case ExpressionConditionKind::ExpressionIsSingleBasisElement:
        return "expression-is-single-basis-element";
      case ExpressionConditionKind::AllPowerSumTermsAreSingleCycles:
        return "all-power-sum-terms-are-single-cycles";
      case ExpressionConditionKind::
          MixedShortCyclePowerSumExpansion:
        return "mixed-short-cycle-power-sum-expansion";
      case ExpressionConditionKind::PowerSumIndexHasMostlyShortCycles:
        return "power-sum-index-has-mostly-short-cycles";
      case ExpressionConditionKind::CombinatorialTagsEqual:
        return "combinatorial-tags == " +
               std::to_string(condition.tagMask);
      case ExpressionConditionKind::HasCombinatorialTag:
        return "has-combinatorial-tag(" +
               std::to_string(condition.tagMask) + ")";
      case ExpressionConditionKind::CoefficientRingIsQQ:
        return "coefficient-ring-is-QQ";
      case ExpressionConditionKind::And:
        return logicalString(condition, "and");
      case ExpressionConditionKind::Or:
        return logicalString(condition, "or");
      case ExpressionConditionKind::Not:
        if (condition.operands.size() != 1) return "not(<invalid>)";
        return "not(" +
               expressionConditionToString(condition.operands[0]) +
               ")";
    }
  return "<unknown-condition>";
}

const char *expressionPieceKindName(ExpressionPieceKind kind)
{
  switch (kind)
    {
      case ExpressionPieceKind::WholeExpression:
        return "whole-expression pieces";
      case ExpressionPieceKind::IndividualTerms:
        return "term pieces";
      case ExpressionPieceKind::HomogeneousComponents:
        return "homogeneous-component pieces";
    }
  return "unknown pieces";
}

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
