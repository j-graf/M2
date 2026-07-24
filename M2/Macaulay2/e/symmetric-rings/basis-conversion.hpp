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
    uint32_t factorKindMask = 0;
    size_t factorCount = 0;
    size_t skewFactorCount = 0;
    std::optional<int> pureBasis;
    std::optional<int> expandedBasis;
    std::optional<int> homogeneousWeight;
    size_t maximumPartitionLength = 0;
    std::optional<bool> allPowerSumTermsSingleCycles;
    std::optional<size_t> completeFriendlyPowerSumTermCount;
    std::optional<std::vector<int>> commonPowerSumParts;
    // Selection profiles are expression- and target-specific. This marker
    // lets the selected executor reuse them without rescanning the same input.
    std::optional<int> planSelectionTargetBasisId;
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
  void enrichBasisConversionPlanSelectionFacts(
      ring_elem f,
      int targetBasisId,
      ExpressionFacts& facts) const;
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
  std::vector<ExpressionConditionContext>
  buildExpressionConditionContexts(
      ring_elem expression,
      ExpressionPieceKind pieces,
      const ExpressionFacts& facts) const;

// ============================================================================
// Basis-Conversion Plan Contracts
// ============================================================================
// A plan is one complete source-to-target conversion. The registry is a
// policy-free catalog, the picker chooses one top-level plan, and the generic
// executor follows that plan's fixed formulas without choosing a replacement.

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
    PowerSumsToSchurCharacters,
    PowerSumsToSchurOmegaBorderStrips,
    PowerSumsToSchurOmegaAbacusRimHooks,
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
    SchurToCompleteJacobiTrudi,
    SchurOmegaToElementaryJacobiTrudi,
    CompleteToSchurRecursive,
    HallLittlewoodGeneratorToCapitalTriangular
  };

  // The registry stores complete source-to-target plans. A case formula is an
  // atomic kernel or a fixed sequence of named plans. A one-plan sequence is
  // delegation; a longer sequence is mathematical composition.
  struct BasisConversionPlanId
  {
    std::string value;

    bool operator==(const BasisConversionPlanId& other) const
    {
      return value == other.value;
    }
  };

  enum class ConversionFormulaKind
  {
    AtomicKernel,
    PlanComposition
  };

  struct BasisConversionPlanDefinition;

  struct ConversionFormula
  {
    ConversionFormulaKind kind =
        ConversionFormulaKind::AtomicKernel;
    BasisConversionKernel kernel =
        BasisConversionKernel::Unavailable;
    std::vector<BasisConversionPlanId> childPlans;
    // Database validation resolves stable identifiers once. Runtime execution
    // then follows definition pointers and resolves only ring-local basis IDs.
    mutable std::vector<const BasisConversionPlanDefinition *>
        resolvedChildPlans;
  };

  struct BasisConversionPlanCase
  {
    ExpressionCondition condition;
    ConversionFormula formula;
  };

  struct BasisConversionPlanDefinition
  {
    BasisConversionPlanId id;
    BasisKind sourceBasisKind;
    BasisKind targetBasisKind;
    ExpressionCondition applicability;
    ExpressionPieceKind pieces;
    std::vector<BasisConversionPlanCase> cases;
    // Every plan already guarantees a canonical target expansion. This
    // condition records only a stronger mathematical postcondition, when one
    // is needed to prove a later child's applicability. always() means that
    // the plan promises no additional shape/profile condition.
    ExpressionCondition outputGuarantee = always();
  };

  struct ResolvedBasisConversionPlan
  {
    const BasisConversionPlanDefinition *definition = nullptr;
    int sourceBasisId = -1;
    int targetBasisId = -1;

    bool valid() const { return definition != nullptr; }
  };

// ============================================================================
// Basis-Conversion Plan Registry
// ============================================================================
// Registration declares mathematical availability without selection policy.

  static ConversionFormula kernelFormula(
      BasisConversionKernel kernel);
  static ConversionFormula planFormula(
      BasisConversionPlanId childPlan);
  static ConversionFormula compositionFormula(
      std::initializer_list<BasisConversionPlanId> childPlans);
  static const std::vector<BasisConversionPlanDefinition>&
  basisConversionPlanDatabase();
  static const std::map<
      std::string,
      const BasisConversionPlanDefinition *>&
  basisConversionPlansById();
  static const std::map<
      std::pair<BasisKind, BasisKind>,
      std::vector<const BasisConversionPlanDefinition *>>&
  basisConversionPlansByEndpoints();
  static const BasisConversionPlanDefinition *
  basisConversionPlanDefinition(
      const BasisConversionPlanId& id);
  struct BasisConversionKernelContract
  {
    BasisConversionKernelContract() = default;
    BasisConversionKernelContract(bool supports)
        : supportsEndpoints(supports)
    {
    }
    bool supportsEndpoints = false;
    // Conversion kernels always return canonical target expansions. This
    // condition records only a stronger mathematical postcondition available
    // to a following child plan.
    ExpressionCondition outputGuarantee = always();
  };
  static BasisConversionKernelContract basisConversionKernelContract(
      BasisConversionKernel kernel,
      BasisKind source,
      BasisKind target);
  static void validateBasisConversionPlanDatabase();
  const std::vector<ResolvedBasisConversionPlan>&
  registeredBasisConversionPlans(
      int sourceBasisId,
      int targetBasisId) const;
  ResolvedBasisConversionPlan resolveBasisConversionPlan(
      const BasisConversionPlanId& id) const;
  ResolvedBasisConversionPlan resolveBasisConversionPlanDefinition(
      const BasisConversionPlanDefinition *definition) const;
  mutable std::map<
      std::pair<int, int>,
      std::vector<ResolvedBasisConversionPlan>>
      resolvedBasisConversionPlanCache;

// ============================================================================
// Basis-Conversion Applicability And Selection
// ============================================================================
// Applicability is mathematical; costs and pickers own performance policy.

  bool basisConversionPlanApplicable(
      const ResolvedBasisConversionPlan& plan,
      ring_elem expression,
      const ExpressionFacts& facts) const;
  size_t basisConversionPlanCost(
      const ResolvedBasisConversionPlan& plan,
      const ExpressionFacts& facts,
      CombinatorialTags combinatorialTags) const;
  ResolvedBasisConversionPlan selectBasisConversionPlan(
      ring_elem expression,
      int sourceBasisId,
      int targetBasisId,
      const ExpressionFacts& facts,
      CombinatorialTags combinatorialTags,
      const std::optional<std::string>& forcedIdentifier =
          std::nullopt,
      ExpressionFacts *selectionFacts = nullptr) const;

// ============================================================================
// Basis-Conversion Execution
// ============================================================================
// Executors run selected plans and never choose a replacement.

  ring_elem executeBasisConversionKernel(
      BasisConversionKernel kernel,
      int sourceBasisId,
      int targetBasisId,
      ring_elem expression,
      CombinatorialTags combinatorialTags,
      ExpressionFacts *resultFacts,
      const std::optional<int>& knownHomogeneousWeight =
          std::nullopt) const;
  ring_elem executeBasisConversionFormula(
      const ConversionFormula& formula,
      int sourceBasisId,
      int targetBasisId,
      ring_elem expression,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts& inputFacts,
      ExpressionFacts *resultFacts = nullptr) const;
  ring_elem executeBasisConversionPlan(
      const ResolvedBasisConversionPlan& plan,
      ring_elem expression,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts& inputFacts,
      ExpressionFacts *resultFacts = nullptr) const;
  bool checkAllApplicableBasisConversionPlans(
      ring_elem expression,
      int sourceBasisId,
      int targetBasisId,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts& inputFacts,
      ring_elem expectedResult) const;
  ring_elem convertCanonicalExpressionToBasis(
      ring_elem f,
      int targetBasisId,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts *knownFacts = nullptr) const;

// ============================================================================
// Multiplication Plan Contracts
// ============================================================================
// A multiplication plan declares every operand conversion and the basis and
// canonical-form guarantee of its kernel output. Operand conversions are fixed
// before kernel execution; the owning workflow selects a support-dependent
// post-kernel conversion only after the declared output has been realized.

  enum class MultiplicationKernel
  {
    SchurCompatibleFactors,
    MonomialLikeExpansion,
    CanonicalBasisProduct
  };

  enum class MultiplicationPlanId
  {
    Unavailable,
    SchurLittlewoodRichardson,
    SchurHorizontalPieri,
    SchurVerticalPieri,
    SchurBorderStrips,
    SchurCompatibleFactorRules,
    MonomialLikeExponentSplittings,
    HallLittlewoodGenerators,
    MultiplicativeTarget,
    PowerSums
  };

  struct MultiplicationPlan
  {
    MultiplicationPlanId id;
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
    std::optional<ResolvedBasisConversionPlan> leftConversion;
    std::optional<ResolvedBasisConversionPlan> rightConversion;
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
  static bool isSchurCompatibleBasisKind(BasisKind kind);
  static bool isMonomialLikeCompatibleBasisKind(
      BasisKind kind, BasisKind target);
  static std::string_view multiplicationPlanIdentifier(
      MultiplicationPlanId id);
  bool multiplicationPlanApplicable(
      const MultiplicationPlan& plan,
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId) const;
  MultiplicationPlanSelection selectMultiplicationPlan(
      ring_elem f,
      ring_elem g,
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId) const;
  ring_elem executeMultiplicationKernel(
      const MultiplicationPlan& plan,
      ring_elem f,
      ring_elem g,
      int targetBasisId) const;
  ring_elem runMultiplicationWorkflow(
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
  struct ResolvedProductTerm
  {
    ring_elem expression;
    ExpressionFacts facts;
  };
  ResolvedProductTerm multiplyTermToBasis(
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
