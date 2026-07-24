// Copyright 2026

#include "symmetric-rings/expression-conditions.hpp"

#include <gtest/gtest.h>

using namespace symmetric_rings;

TEST(ExpressionConditions, LogicalValuesPreserveShortCircuitEvaluation)
{
  ExpressionPieceFacts term;
  term.pieceKind = ExpressionPieceKind::IndividualTerms;
  term.weight = 12;

  EXPECT_TRUE(expressionConditionHolds(
      termWeightGreaterThan(10) ||
          indexIsHook(),
      term));
  EXPECT_FALSE(expressionConditionHolds(
      termWeightAtMost(10) &&
          indexIsHook(),
      term));
  EXPECT_EQ(
      "(term-weight > 10 and not(index-is-self-conjugate))",
      expressionConditionToString(
          termWeightGreaterThan(10) &&
          !indexIsSelfConjugate()));
}

TEST(ExpressionConditions, PartitionCasesAreOrderedAndComplete)
{
  ExpressionPieceFacts first;
  first.pieceKind = ExpressionPieceKind::IndividualTerms;
  first.termPositions = {0};
  first.weight = 4;
  first.index = Partition{3, 1};

  ExpressionPieceFacts second;
  second.pieceKind = ExpressionPieceKind::IndividualTerms;
  second.termPositions = {1};
  second.weight = 4;
  second.index = Partition{2, 2};

  ExpressionPieceFacts third;
  third.pieceKind = ExpressionPieceKind::IndividualTerms;
  third.termPositions = {2};
  third.weight = 5;
  third.index = Partition{2, 2, 1};

  const auto assignments = assignExpressionPiecesToCases(
      {first, second, third},
      {indexIsHook(), indexIsRectangle(), otherwise()},
      3);

  EXPECT_EQ(
      (std::vector<size_t>{0}),
      assignments.termPositionsForCase[0]);
  EXPECT_EQ(
      (std::vector<size_t>{1}),
      assignments.termPositionsForCase[1]);
  EXPECT_EQ(
      (std::vector<size_t>{2}),
      assignments.termPositionsForCase[2]);
}

TEST(ExpressionConditions, PartitionRejectsInvalidCaseDefinitions)
{
  ExpressionPieceFacts term;
  term.pieceKind = ExpressionPieceKind::IndividualTerms;
  term.termPositions = {0};
  term.weight = 1;

  EXPECT_THROW(
      assignExpressionPiecesToCases(
          {term}, {otherwise(), termWeightAtMost(2)}, 1),
      std::invalid_argument);
  EXPECT_THROW(
      assignExpressionPiecesToCases(
          {term}, {termWeightAtMost(2)}, 1),
      std::invalid_argument);
  EXPECT_THROW(
      assignExpressionPiecesToCases(
          {term},
          {termWeightAtMost(2) || otherwise(), otherwise()},
          1),
      std::invalid_argument);
}

TEST(ExpressionConditions, StaticValidationChecksPieceKindsAndArity)
{
  EXPECT_NO_THROW(validateExpressionCondition(
      indexIsHook() && termWeightGreaterThan(3),
      ExpressionPieceKind::IndividualTerms));
  EXPECT_THROW(
      validateExpressionCondition(
          indexIsHook(),
          ExpressionPieceKind::HomogeneousComponents),
      std::invalid_argument);
  EXPECT_THROW(
      validateExpressionCondition(
          otherwise(),
          ExpressionPieceKind::IndividualTerms),
      std::invalid_argument);

  ExpressionCondition invalidAnd;
  invalidAnd.kind = ExpressionConditionKind::And;
  invalidAnd.operands = {indexIsHook()};
  EXPECT_THROW(
      validateExpressionCondition(
          invalidAnd,
          ExpressionPieceKind::IndividualTerms),
      std::invalid_argument);
}

TEST(ExpressionConditions, ContractImplicationIsConservative)
{
  EXPECT_TRUE(expressionConditionImplies(
      termWeightAtMost(5),
      termWeightAtMost(10)));
  EXPECT_FALSE(expressionConditionImplies(
      termWeightAtMost(10),
      termWeightAtMost(5)));
  EXPECT_TRUE(expressionConditionImplies(
      componentDensityAtLeast(1, 2),
      componentDensityAtLeast(1, 4)));
  EXPECT_TRUE(expressionConditionImplies(
      coefficientRingIsQQ() &&
          expressionIsSingleBasisElement(),
      expressionIsSingleBasisElement()));
  EXPECT_TRUE(expressionConditionImplies(
      coefficientRingIsQQ() &&
          expressionIsSingleBasisElement(),
      expressionIsSingleBasisElement() &&
          coefficientRingIsQQ()));
  EXPECT_TRUE(expressionConditionImplies(
      expressionIsSingleBasisElement(),
      expressionIsSingleBasisElement() ||
          expressionHasMultipleWeights()));
  EXPECT_FALSE(expressionConditionImplies(
      expressionHasMultipleWeights(),
      expressionIsSingleBasisElement()));
}

TEST(ExpressionConditions, ComponentDensityUsesExactCrossMultiplication)
{
  ExpressionPieceFacts component;
  component.pieceKind =
      ExpressionPieceKind::HomogeneousComponents;
  component.termCount = 3;
  component.possibleTermCount = 12;

  EXPECT_TRUE(expressionConditionHolds(
      componentDensityAtLeast(1, 4), component));
  EXPECT_FALSE(expressionConditionHolds(
      componentDensityAtLeast(1, 3), component));
}

TEST(ExpressionConditions, ComponentProfilesAndFullTagMasksArePreserved)
{
  ExpressionPieceFacts component;
  component.pieceKind =
      ExpressionPieceKind::HomogeneousComponents;
  component.termCount = 40;
  component.mostlyShortCyclePowerSumTermCount = 10;
  component.combinatorialTags = uint32_t{1} << 31;

  EXPECT_TRUE(expressionConditionHolds(
      componentMostlyShortCycleFractionAtLeast(1, 4, 8),
      component));
  EXPECT_TRUE(expressionConditionHolds(
      combinatorialTagsEqual(uint32_t{1} << 31),
      component));
  EXPECT_TRUE(expressionConditionHolds(
      hasCombinatorialTag(uint32_t{1} << 31),
      component));
}

TEST(ExpressionConditions, PartitionShapePrimitivesAreInspectable)
{
  ExpressionPieceFacts term;
  term.pieceKind = ExpressionPieceKind::IndividualTerms;

  term.index = Partition{5, 1, 1};
  EXPECT_TRUE(expressionConditionHolds(indexIsHook(), term));
  EXPECT_FALSE(expressionConditionHolds(indexIsRectangle(), term));

  term.index = Partition{3, 3};
  EXPECT_TRUE(expressionConditionHolds(indexIsRectangle(), term));

  term.index = Partition{3, 2, 1};
  EXPECT_TRUE(
      expressionConditionHolds(indexIsSelfConjugate(), term));

  term.index = Partition{3, 3, 0};
  EXPECT_TRUE(expressionConditionHolds(indexIsRectangle(), term));

  term.index = Partition{3, 2, 1, 0};
  EXPECT_TRUE(
      expressionConditionHolds(indexIsSelfConjugate(), term));

  term.index = Partition{0};
  EXPECT_FALSE(expressionConditionHolds(indexIsHook(), term));
}

// Local Variables:
// indent-tabs-mode: nil
// End:
