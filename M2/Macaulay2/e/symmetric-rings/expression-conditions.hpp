// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_EXPRESSION_CONDITIONS_HPP_
#define M2_SYMMETRIC_RINGS_EXPRESSION_CONDITIONS_HPP_

#include "symmetric-rings/partitions.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Expression Pieces And Inspectable Conditions
// ============================================================================

enum class ExpressionPieceKind
{
  WholeExpression,
  IndividualTerms,
  HomogeneousComponents
};

enum class ExpressionConditionKind
{
  Always,
  Otherwise,
  TermWeightAtMost,
  TermWeightGreaterThan,
  ComponentWeightAtMost,
  ComponentWeightGreaterThan,
  ComponentTermCountAtLeast,
  ComponentDensityAtLeast,
  ComponentSupportSquareRatioAtLeast,
  ComponentAllPowerSumIndicesHaveMostlyShortCycles,
  ComponentHasBothMostlyShortCycleAndOtherTerms,
  ComponentMostlyShortCycleFractionAtLeast,
  ComponentHasCommonPowerSumPartOne,
  ComponentHasCommonPowerSumPartAtMostPercentOfWeight,
  IndexIsHook,
  IndexIsRectangle,
  IndexIsSelfConjugate,
  ExpressionHasMultipleWeights,
  ExpressionIsSingleBasisElement,
  AllPowerSumTermsAreSingleCycles,
  MixedShortCyclePowerSumExpansion,
  PowerSumIndexHasMostlyShortCycles,
  CombinatorialTagsEqual,
  HasCombinatorialTag,
  CoefficientRingIsQQ,
  And,
  Or,
  Not
};

struct ExpressionCondition
{
  ExpressionConditionKind kind = ExpressionConditionKind::Always;
  int firstInteger = 0;
  int secondInteger = 1;
  int thirdInteger = 0;
  uint32_t tagMask = 0;
  std::vector<ExpressionCondition> operands;
};

struct ExpressionPieceFacts
{
  ExpressionPieceKind pieceKind = ExpressionPieceKind::WholeExpression;
  std::vector<size_t> termPositions;
  std::optional<int> weight;
  std::optional<size_t> termCount;
  std::optional<size_t> possibleTermCount;
  std::optional<Partition> index;
  std::optional<bool> singleBasisElement;
  std::optional<bool> allPowerSumTermsSingleCycles;
  std::optional<size_t> mostlyShortCyclePowerSumTermCount;
  std::optional<std::vector<int>> commonPowerSumParts;
  std::optional<bool> powerSumIndexHasMostlyShortCycles;
  std::optional<bool> expressionHasMultipleWeights;
  std::optional<uint32_t> combinatorialTags;
  std::optional<bool> coefficientRingIsQQ;
};

struct PlanCaseAssignments
{
  // Each entry contains the source term positions assigned to the
  // corresponding ordered condition.
  std::vector<std::vector<size_t>> termPositionsForCase;
};

// ============================================================================
// Condition Construction
// ============================================================================

ExpressionCondition always();
ExpressionCondition otherwise();
ExpressionCondition termWeightAtMost(int weight);
ExpressionCondition termWeightGreaterThan(int weight);
ExpressionCondition componentWeightAtMost(int weight);
ExpressionCondition componentWeightGreaterThan(int weight);
ExpressionCondition componentTermCountAtLeast(size_t count);
ExpressionCondition componentDensityAtLeast(
    size_t numerator, size_t denominator);
ExpressionCondition componentSupportSquareRatioAtLeast(
    size_t numerator, size_t denominator);
ExpressionCondition componentAllPowerSumIndicesHaveMostlyShortCycles();
ExpressionCondition componentHasBothMostlyShortCycleAndOtherTerms();
ExpressionCondition componentMostlyShortCycleFractionAtLeast(
    size_t numerator, size_t denominator, size_t minimumCount);
ExpressionCondition componentHasCommonPowerSumPartOne();
ExpressionCondition componentHasCommonPowerSumPartAtMostPercentOfWeight(
    int percent);
ExpressionCondition indexIsHook();
ExpressionCondition indexIsRectangle();
ExpressionCondition indexIsSelfConjugate();
ExpressionCondition expressionHasMultipleWeights();
ExpressionCondition expressionIsSingleBasisElement();
ExpressionCondition allPowerSumTermsAreSingleCycles();
ExpressionCondition mixedShortCyclePowerSumExpansion();
ExpressionCondition powerSumIndexHasMostlyShortCycles();
ExpressionCondition combinatorialTagsEqual(uint32_t tags);
ExpressionCondition hasCombinatorialTag(uint32_t tag);
ExpressionCondition coefficientRingIsQQ();

ExpressionCondition operator&&(
    ExpressionCondition left, ExpressionCondition right);
ExpressionCondition operator||(
    ExpressionCondition left, ExpressionCondition right);
ExpressionCondition operator!(ExpressionCondition condition);

// ============================================================================
// Condition Evaluation, Partitioning, And Diagnostics
// ============================================================================

bool expressionConditionHolds(
    const ExpressionCondition& condition,
    const ExpressionPieceFacts& piece);

void validateExpressionCondition(
    const ExpressionCondition& condition,
    ExpressionPieceKind pieceKind,
    bool allowOtherwise = false);

// Returns true only when the inspectable condition language can prove the
// implication. False means "not proved", not necessarily mathematically false.
bool expressionConditionImplies(
    const ExpressionCondition& guarantee,
    const ExpressionCondition& requirement);

PlanCaseAssignments assignExpressionPiecesToCases(
    const std::vector<ExpressionPieceFacts>& pieces,
    const std::vector<ExpressionCondition>& orderedConditions,
    size_t expectedTermCount);

std::string expressionConditionToString(
    const ExpressionCondition& condition);
const char *expressionPieceKindName(ExpressionPieceKind kind);

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
