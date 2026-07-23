// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_CONVERSION_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_CONVERSION_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Expression Facts
// ============================================================================
// Conversion and multiplication share one exact facts representation.

  struct ExpressionFacts
  {
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
    uint32_t factorKindMask = 0;
    size_t factorCount = 0;
    size_t skewFactorCount = 0;
    std::optional<int> pureBasis;
    std::optional<int> expandedBasis;
    std::map<int, std::vector<size_t>> weightTermPositions;
    std::optional<int> homogeneousWeight;
    size_t maximumPartitionLength = 0;
    std::optional<size_t> possibleTermCount;
    std::optional<double> density;
    std::optional<bool> allPowerSumTermsSingleCycles;
    std::optional<size_t> completeFriendlyPowerSumTermCount;
    std::optional<std::vector<int>> commonPowerSumParts;
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
    bool hasFactorKind(BasisKind kind) const
    {
      return (factorKindMask &
              (uint32_t{1} << static_cast<size_t>(kind))) != 0;
    }
    bool allFactorsSchurCompatible() const
    {
      uint32_t compatible = 0;
      for (BasisKind kind : {BasisKind::Schur,
                             BasisKind::Complete,
                             BasisKind::Elementary,
                             BasisKind::PowerSum})
        compatible |=
            uint32_t{1} << static_cast<size_t>(kind);
      return factorCount != 0 &&
             (factorKindMask & ~compatible) == 0;
    }
    bool allFactorsHallLittlewoodGenerators() const
    {
      uint32_t generators = 0;
      for (BasisKind kind : {BasisKind::HallLittlewoodQGenerator,
                             BasisKind::HallLittlewoodBGenerator})
        generators |=
            uint32_t{1} << static_cast<size_t>(kind);
      return factorCount != 0 &&
             (factorKindMask & ~generators) == 0;
    }
    bool singleBasisElement() const
    {
      return termCount == 1 && singleFactorTermCount == 1 && normalized;
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

// ============================================================================
// Expression Preparation And Metadata
// ============================================================================
// Preparation establishes the canonical, skew-free term representation
// consumed by both the conversion registry and multiplication workflows.

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
  void enrichConversionSelectionFacts(
      ring_elem f,
      int targetBasisId,
      ExpressionFacts& facts) const;
  std::optional<ExpressionFacts> expressionFactsFromMetadata(
      ring_elem f) const;
  void attachExpressionFacts(
      ring_elem f,
      const ExpressionFacts& facts,
      int targetBasisId,
      CombinatorialTags combinatorialTags) const;
  ring_elem expressionFromAtom(
      const SymmetricMonomial& monomial,
      size_t pos) const;
  ring_elem skewBasisElementExpansion(
      const SymmetricMonomial& monomial,
      size_t pos) const;
  ring_elem normalizeExpression(
      ring_elem f,
      ExpressionFacts& resultFacts,
      std::vector<size_t> *termFactorCounts = nullptr) const;
  ExpressionFacts conservativeExpressionFactsForBasis(
      const ExpressionFacts& sourceFacts,
      int basisId) const;

// ============================================================================
// Basis-Conversion Plan Contracts
// ============================================================================
// A plan is one complete conversion between adjacent bases.  The registry is a
// policy-free catalog, the picker chooses a composition, and the executor
// performs the already-selected plans without choosing a replacement.

  enum class BasisConversionKernel
  {
    // Error sentinel returned only after selection has reported an error.
    // It is never registered as a mathematical plan.
    Unavailable,
    CompleteToPowerSumsClassical,
    ElementaryToPowerSumsClassical,
    SchurToPowerSumsCharacters,
    SchurOmegaToPowerSumsCharacters,
    MonomialToPowerSumsTransition,
    ForgottenToPowerSumsTransition,
    HallLittlewoodQGeneratorToPowerSumsClassical,
    HallLittlewoodBGeneratorToPowerSumsClassical,
    HallLittlewoodQToPowerSumsRaising,
    HallLittlewoodBToPowerSumsRaising,
    HallLittlewoodPToPowerSumsNormalization,
    HallLittlewoodPOmegaToPowerSumsNormalization,
    PowerSumsToCompleteLogarithm,
    PowerSumsToElementaryLogarithm,
    PowerSumsToSchurBorderStrips,
    PowerSumsToSchurAbacusRimHooks,
    PowerSumsToSchurAbacusAndComplete,
    PowerSumsToSchurComplete,
    PowerSumsToSchurCharacters,
    PowerSumsToSchurOmegaBorderStrips,
    PowerSumsToSchurOmegaAbacusRimHooks,
    PowerSumsToSchurOmegaAbacusAndComplete,
    PowerSumsToSchurOmegaComplete,
    PowerSumsToSchurOmegaCharacters,
    PowerSumsToHallLittlewoodQGeneratorLogarithm,
    PowerSumsToHallLittlewoodBGeneratorLogarithm,
    PowerSumSingleCyclesToHallLittlewoodGreen,
    PowerSumIndexToHallLittlewoodGreenDuality,
    PowerSumsToHallLittlewoodTriangular,
    PowerSumsToMonomialTransition,
    PowerSumsToForgottenTransition,
    HallLittlewoodNormalization,
    SchurOmegaConjugation,
    CompleteToSchurRecursive,
    HallLittlewoodGeneratorToCapitalTriangular
  };

  enum class BasisConversionInputShape
  {
    CanonicalExpansion,
    SingleBasisElement,
    SingleCyclePowerSumExpansion,
    MixedCompleteFriendlyPowerSumExpansion
  };

  struct BasisConversionPlan
  {
    std::string_view identifier;
    int sourceBasisId;
    int targetBasisId;
    BasisConversionKernel kernel;
    BasisConversionInputShape inputShape;
  };

  struct BasisConversionPlanSelection
  {
    std::vector<int> basisComposition;
    std::vector<BasisConversionPlan> plans;
    // A later edge may have been provisionally selected from guarantees
    // because its exact intermediate expression did not yet exist. Such an
    // edge is reselected from exact facts at its execution stage boundary.
    // Explicitly forced plans are never reselected.
    std::vector<bool> reselectFromExactIntermediate;
  };

// ============================================================================
// Basis-Conversion Plan Registry
// ============================================================================
// Registration declares mathematical availability without selection policy.

  const std::vector<BasisConversionPlan>&
  registeredBasisConversionPlans(
      int sourceBasisId,
      int targetBasisId) const;
  std::vector<BasisConversionPlan> buildBasisConversionPlans(
      int sourceBasisId,
      int targetBasisId) const;

// ============================================================================
// Basis-Conversion Applicability And Selection
// ============================================================================
// Applicability is mathematical; costs and pickers own performance policy.

  bool basisConversionPlanApplicable(
      const BasisConversionPlan& plan,
      const ExpressionFacts& facts) const;
  size_t basisConversionPlanCost(
      const BasisConversionPlan& plan,
      const ExpressionFacts& facts,
      CombinatorialTags combinatorialTags) const;
  BasisConversionPlan selectBasisConversionPlan(
      int sourceBasisId,
      int targetBasisId,
      const ExpressionFacts& facts,
      CombinatorialTags combinatorialTags,
      const std::optional<std::string>& forcedIdentifier =
          std::nullopt) const;
  BasisConversionPlanSelection pickBasisConversionPlans(
      int sourceBasisId,
      int targetBasisId,
      const ExpressionFacts& facts,
      CombinatorialTags combinatorialTags) const;

// ============================================================================
// Composite Schur Kernels
// ============================================================================

  ring_elem powerSumsToSchurViaComplete(
      ring_elem f,
      int targetBasisId,
      const std::string& targetDisplay,
      int targetOrder) const;
  ring_elem powerSumsToSchurViaAbacusAndComplete(
      ring_elem f,
      int targetBasisId,
      const std::string& targetDisplay,
      int targetOrder) const;

// ============================================================================
// Basis-Conversion Execution
// ============================================================================
// Executors run selected plans and never choose a replacement.

  ring_elem executeBasisConversionPlan(
      const BasisConversionPlan& plan,
      ring_elem expression,
      CombinatorialTags combinatorialTags,
      ExpressionFacts *resultFacts,
      const std::optional<int>& knownHomogeneousWeight =
          std::nullopt) const;
  ring_elem executeBasisConversionPlans(
      const BasisConversionPlanSelection& selection,
      ring_elem expression,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts& inputFacts,
      ExpressionFacts *resultFacts = nullptr) const;
  ring_elem convertCanonicalExpressionToBasis(
      ring_elem f,
      int targetBasisId,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts *knownFacts = nullptr) const;
  mutable std::map<
      std::pair<int, int>,
      std::vector<BasisConversionPlan>>
      basisConversionPlanRegistryCache;

// ============================================================================
// Multiplication Plan Contracts
// ============================================================================
// A multiplication plan declares every operand conversion and the basis and
// canonical-form guarantee of its kernel output. Basis compositions are fixed
// before execution. Data-dependent kernels within those compositions are
// selected at their stage boundaries from exact realized facts, matching the
// owning workflow's staged-selection contract.

  enum class MultiplicationKernel
  {
    SchurCompatibleFactors,
    MonomialLikeExpansion,
    CanonicalBasisProduct
  };

  struct MultiplicationPlan
  {
    std::string_view identifier;
    MultiplicationKernel kernel;
    // A positive operand basis requests an explicit pre-kernel conversion.
    // Structured kernels use -1 because they consume the original compatible
    // basis kinds directly.
    int leftKernelBasisId;
    int rightKernelBasisId;
    // Every kernel declares both its mathematical output basis and whether
    // normalization can be bypassed for that output.
    int kernelOutputBasisId;
    bool kernelOutputCanonical;
  };

  struct MultiplicationPlanSelection
  {
    MultiplicationPlan plan;
    std::optional<BasisConversionPlanSelection> leftConversion;
    std::optional<BasisConversionPlanSelection> rightConversion;
  };

// ============================================================================
// Multiplication Registry And Workflow Helpers
// ============================================================================

  const std::vector<MultiplicationPlan>&
  registeredMultiplicationPlans(
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId) const;
  std::vector<MultiplicationPlan> buildMultiplicationPlans(
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId) const;
  bool multiplicationPlanApplicable(
      const MultiplicationPlan& plan,
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId) const;
  MultiplicationPlanSelection selectMultiplicationPlan(
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId) const;
  ring_elem executeMultiplicationKernel(
      const MultiplicationPlan& plan,
      ring_elem f,
      ring_elem g,
      int targetBasisId) const;
  ring_elem executeMultiplicationWorkflow(
      const MultiplicationPlanSelection& selection,
      ring_elem f,
      ring_elem g,
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId,
      bool attachResultFacts) const;
  ring_elem multiplyBasisElementsWithFacts(
      ring_elem f,
      ring_elem g,
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId,
      bool attachResultFacts) const;
  ring_elem multiplyBasisElementsToBasis(
      ring_elem f,
      ring_elem g,
      int targetBasisId) const;
  ring_elem multiplyTermToBasis(
      const SymmetricTerm& term,
      int targetBasisId) const;
  mutable std::map<
      std::tuple<int, int, int>,
      std::vector<MultiplicationPlan>>
      multiplicationPlanRegistryCache;

// ============================================================================
// Public Conversion And Multiplication Entry Points
// ============================================================================
// Definitions follow dependency order: multiplication helpers precede
// toBasis because product resolution delegates to them.

 public:
  ring_elem multiplyToBasis(
      ring_elem f,
      ring_elem g,
      int targetBasisId) const;
  ring_elem toBasis(ring_elem f, int targetBasisId) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
