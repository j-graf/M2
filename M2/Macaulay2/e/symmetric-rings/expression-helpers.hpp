// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_EXPRESSION_HELPERS_HPP_
#define M2_SYMMETRIC_RINGS_EXPRESSION_HELPERS_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// General Expression Helpers
// ============================================================================
// These helpers expose stable mathematical decompositions and normalization
// operations. They do not contain operation-specific selection policy.

 public:
  struct HomogeneousComponent
  {
    int weight;
    ring_elem expression;
  };

  struct BasisComponent
  {
    int basisId;
    ring_elem expression;
  };

  struct HomogeneousBasisComponent
  {
    int weight;
    int basisId;
    ring_elem expression;
  };

  struct ExpressionShape
  {
    std::vector<int> weights;
    std::vector<int> basisIds;
    size_t termCount = 0;
    size_t scalarTermCount = 0;
    size_t singleFactorTermCount = 0;
    size_t productTermCount = 0;
    size_t maximumFactorsPerTerm = 0;
    size_t skewFactorCount = 0;
    size_t maximumPartitionLength = 0;
    bool normalized = true;
    bool skewFree = true;
    bool collected = true;
    std::optional<int> homogeneousWeight;
    std::optional<int> pureBasis;
    std::optional<int> expandedBasis;
    bool hasMetadata = false;
    bool metadataFactsComplete = false;
    bool metadataNormalized = false;
    bool metadataSkewFree = false;
    bool metadataCollected = false;
  };

  struct ExpressionNormalizationOptions
  {
    bool straightenIndices = true;
    bool expandSkewFactors = false;
    int productTargetBasisId = -1;
  };

  struct ExpressionNormalizationResult
  {
    // This is the complete preparation payload required by the conversion
    // and multiplication boundaries. An empty factor-count vector proves
    // that no product term needs a per-term profile.
    ring_elem expression;
    ExpressionFacts facts;
    std::vector<size_t> factorsPerTerm;
    bool usedCompleteMetadataBypass = false;
    bool usedAlreadyPreparedBypass = false;
    bool bypassedStraighteningFromMetadata = false;
    bool bypassedSkewExpansionFromMetadata = false;
    bool performedStraightening = false;
    bool performedSkewExpansion = false;
    bool resolvedProducts = false;
  };

  static ExpressionNormalizationOptions
  pipelinePreparationNormalizationOptions();
  std::vector<HomogeneousComponent> homogeneousComponents(
      ring_elem expression) const;
  ring_elem homogeneousComponent(ring_elem expression, int weight) const;
  std::vector<int> weightSupport(ring_elem expression) const;
  ring_elem truncateWeights(
      ring_elem expression,
      int minimumWeight,
      int maximumWeight) const;
  ring_elem normalizeExpressionWithOptions(
      ring_elem expression,
      const ExpressionNormalizationOptions& options) const;
  ExpressionNormalizationResult normalizeExpressionWithOptionsDetailed(
      ring_elem expression,
      const ExpressionNormalizationOptions& options) const;
  ring_elem expandSkewFactors(ring_elem expression) const;
  ring_elem expandProductsInBasis(
      ring_elem expression,
      int targetBasisId) const;
  ExpressionShape expressionShape(ring_elem expression) const;
  std::vector<int> basisSupport(ring_elem expression) const;
  std::vector<BasisComponent> basisComponents(
      ring_elem expression) const;
  bool isBasisExpansion(
      ring_elem expression,
      int basisId) const;
  bool isLinearCombinationOfBasisElements(
      ring_elem expression) const;
  ring_elem coefficientsInBasisExpression(
      ring_elem expression,
      int basisId,
      bool convert) const;
  std::vector<HomogeneousBasisComponent> homogeneousBasisComponents(
      ring_elem expression) const;
  std::vector<ring_elem> singlePartitionIndexedTerms(
      ring_elem expression) const;

 private:
  ring_elem expressionHelperFromTerms(
      VECTOR(SymmetricTerm)& terms,
      CombinatorialTags tags) const;
  void attachExpressionHelperFacts(
      ring_elem expression,
      const ExpressionFacts& facts) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
