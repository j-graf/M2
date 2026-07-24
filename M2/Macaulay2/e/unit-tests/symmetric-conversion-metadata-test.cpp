// Copyright 2026

#include "symmetric-rings/storage.hpp"

#include <gtest/gtest.h>

using namespace symmetric_rings;

TEST(SymmetricConversionMetadata, ExactInvalidationPreservesIndependentHints)
{
  SymmetricConversionMetadata metadata;
  metadata.expressionFactsComplete = true;
  metadata.pureBasis = 7;
  metadata.expandedBasis = 7;
  metadata.homogeneousWeight = 6;
  metadata.termCount = 2;
  metadata.singleBasisElementId = 7;
  metadata.singleBasisElementIndex = Partition{4, 2};
  metadata.singleBasisElementCoefficientOne = true;

  metadata.invalidateExactExpressionFacts();

  EXPECT_FALSE(metadata.expressionFactsComplete);
  EXPECT_EQ(7, metadata.pureBasis);
  EXPECT_EQ(7, metadata.expandedBasis);
  EXPECT_EQ(6, metadata.homogeneousWeight);
  EXPECT_EQ(2U, metadata.termCount);
  EXPECT_EQ(7, metadata.singleBasisElementId);
  EXPECT_EQ((Partition{4, 2}), metadata.singleBasisElementIndex);
  EXPECT_FALSE(metadata.singleBasisElementCoefficientOne);
}

TEST(SymmetricConversionMetadata, CrossRingInvalidationClearsBasisIdentity)
{
  SymmetricConversionMetadata metadata;
  metadata.expressionFactsComplete = true;
  metadata.pureBasis = 7;
  metadata.expandedBasis = 7;
  metadata.factorBases = std::vector<int>{7};
  metadata.homogeneousWeight = 6;
  metadata.termCount = 1;
  metadata.singleBasisElementId = 7;
  metadata.singleBasisElementIndex = Partition{4, 2};
  metadata.singleBasisElementCoefficientOne = true;

  metadata.discardRingLocalBasisFacts();

  EXPECT_FALSE(metadata.expressionFactsComplete);
  EXPECT_FALSE(metadata.pureBasis);
  EXPECT_FALSE(metadata.expandedBasis);
  EXPECT_FALSE(metadata.factorBases);
  EXPECT_FALSE(metadata.singleBasisElementId);
  EXPECT_FALSE(metadata.singleBasisElementIndex);
  EXPECT_FALSE(metadata.singleBasisElementCoefficientOne);
  EXPECT_EQ(6, metadata.homogeneousWeight);
  EXPECT_EQ(1U, metadata.termCount);
}

// Local Variables:
// indent-tabs-mode: nil
// End:
