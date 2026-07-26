// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_EXPRESSION_HELPERS_HPP_
#define M2_SYMMETRIC_RINGS_EXPRESSION_HELPERS_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// General Expression Helpers
// ============================================================================
// These helpers expose stable mathematical decompositions and normalization
// operations. They do not contain operation-specific selection policy.

  struct ExpressionFacts
  {
    // Every ExpressionFacts value is an exact description of one realized
    // expression. Canonical-output contracts are plan/kernel properties and
    // must never be represented by invented support counts in this type.
    bool normalized = true;
    bool skewFree = true;
    bool collected = true;
    CombinatorialTags combinatorialTags = 0;
    size_t termCount = 0;
    size_t scalarTermCount = 0;
    size_t singleFactorTermCount = 0;
    size_t productTermCount = 0;
    size_t maximumFactorsPerTerm = 0;
    std::vector<int> factorBases;
    size_t skewFactorCount = 0;
    std::optional<int> pureBasis;
    std::optional<int> expandedBasis;
    std::optional<int> homogeneousWeight;
    size_t maximumPartitionLength = 0;
    std::optional<int> singleBasisElementId;
    std::optional<Partition> singleBasisElementIndex;
    std::optional<bool> singleBasisElementCoefficientOne;

    bool noProducts() const { return productTermCount == 0; }
    bool singleTerm() const { return termCount == 1; }
    bool mixedBasis() const { return factorBases.size() > 1; }
    bool provenanceUnknown() const { return combinatorialTags == 0; }
    bool provenanceMixed() const
    {
      return combinatorialTags != 0 &&
             (combinatorialTags & (combinatorialTags - 1)) != 0;
    }
    bool singleBasisElement() const
    {
      return termCount == 1 &&
             scalarTermCount == 0 &&
             singleFactorTermCount == 1 &&
             productTermCount == 0 &&
             normalized &&
             skewFree &&
             collected;
    }
    bool canonicalExpansionInBasis(int basisId) const
    {
      if (!normalized || !skewFree || !collected || !noProducts())
        return false;
      // Zero and scalar expressions are canonical in every basis.
      return singleFactorTermCount == 0 ||
             (expandedBasis && *expandedBasis == basisId);
    }
  };

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
  ExpressionFacts inferExpressionFacts(
      ring_elem f,
      std::vector<size_t> *termFactorCounts = nullptr) const;
  ExpressionFacts inferCanonicalExpansionFacts(
      ring_elem f,
      int basisId,
      const std::optional<int>& knownHomogeneousWeight =
          std::nullopt) const;
  ExpressionFacts basisElementFactsFromAtom(
      const SymmetricMonomial& monomial,
      size_t pos,
      ring_elem coefficient) const;
  std::optional<ExpressionFacts> expressionFactsFromMetadata(
      ring_elem f) const;
  void invalidateExactExpressionFacts(
      SymmetricConversionMetadataSlot& metadata) const;
  void refreshSingleBasisElementCoefficientFact(
      SymmetricConversionMetadata& metadata,
      const SymmetricRingPoly *poly) const;
  SymmetricConversionMetadata metadataAfterAddition(
      const SymmetricConversionMetadata& left,
      const SymmetricConversionMetadata& right,
      const SymmetricRingPoly *result) const;
  SymmetricConversionMetadata metadataAfterProduct(
      const SymmetricConversionMetadata& left,
      const SymmetricConversionMetadata& right,
      const SymmetricRingPoly *result) const;
  void attachExpressionFacts(
      ring_elem f,
      const ExpressionFacts& facts,
      int targetBasisId,
      CombinatorialTags combinatorialTags) const;
  bool singleScaledBasisElement(
      ring_elem f,
      int basisId,
      Partition& index,
      ring_elem& coefficient) const;
  CoeffMap coefficientsInBasis(ring_elem f, int basisId) const;
  bool coefficientsInBasisIfPossible(
      ring_elem f,
      int basisId,
      CoeffMap& result) const;
  ring_elem expressionFromBasisElement(
      const SymmetricMonomial& monomial,
      size_t position) const;
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
