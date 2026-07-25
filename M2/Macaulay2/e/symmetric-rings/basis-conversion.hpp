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
    std::optional<size_t> mostlyShortCyclePowerSumTermCount;
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
// consumed by both the conversion plan database and multiplication workflows.

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
  std::vector<ExpressionPieceFacts>
  inspectExpressionPieces(
      ring_elem expression,
      ExpressionPieceKind pieceKind,
      const ExpressionFacts& facts) const;

// ============================================================================
// Basis-Conversion Plan Contracts
// ============================================================================
// A plan is one complete source-to-target conversion. The plan database is
// free of top-level selection, the picker chooses one complete plan, and the
// generic executor follows that plan's fixed piecewise formulas without
// choosing a replacement.

  // The database stores complete source-to-target plans. A case uses either
  // one kernel or a nonempty fixed composition of named plans.
  struct BasisConversionPlanId
  {
    std::string value;

    bool operator==(const BasisConversionPlanId& other) const
    {
      return value == other.value;
    }
  };

  struct BasisConversionPlanDefinition;

  struct KernelPlan
  {
    std::string name;
    BasisConversionKernel kernel = nullptr;
    ExpressionCondition outputGuarantee = always();
  };

  struct CompositionPlan
  {
    std::vector<BasisConversionPlanId> plans;
    mutable std::vector<const BasisConversionPlanDefinition *>
        planDefinitions;
  };

  using BasisConversionPlanFormula =
      std::variant<
          KernelPlan,
          CompositionPlan>;

  struct BasisConversionPlanCase
  {
    ExpressionCondition condition;
    BasisConversionPlanFormula formula;
  };

  struct BasisConversionPlanDefinition
  {
    BasisConversionPlanId id;
    BasisKind sourceBasisKind;
    BasisKind targetBasisKind;
    ExpressionCondition applicability;
    ExpressionPieceKind pieceKind;
    std::vector<BasisConversionPlanCase> cases;
    // Every plan already guarantees a canonical target expansion. This
    // condition records only a stronger mathematical postcondition, when one
    // is needed to prove a later component plan's applicability. always()
    // means that the plan promises no additional shape/profile condition.
    ExpressionCondition outputGuarantee = always();
  };

  struct RingBasisConversionPlan
  {
    const BasisConversionPlanDefinition *definition = nullptr;
    int sourceBasisId = -1;
    int targetBasisId = -1;
    // Generic compositions are instantiated as ordinary immutable plan
    // definitions. Shared ownership keeps that definition stable while the
    // ring-specific plan is copied through selection and execution.
    std::shared_ptr<const BasisConversionPlanDefinition>
        instantiatedDefinition;

    RingBasisConversionPlan() = default;
    RingBasisConversionPlan(
        const BasisConversionPlanDefinition *definitionValue,
        int sourceBasisIdValue,
        int targetBasisIdValue,
        std::shared_ptr<const BasisConversionPlanDefinition>
            instantiatedDefinitionValue = {})
        : definition(definitionValue),
          sourceBasisId(sourceBasisIdValue),
          targetBasisId(targetBasisIdValue),
          instantiatedDefinition(
              std::move(instantiatedDefinitionValue))
    {
    }

    bool valid() const { return definition != nullptr; }
  };

  struct PowerSumFallbackPlanPair
  {
    BasisKind basisKind;
    std::string_view basisName;
    BasisConversionPlanId toPowerSums;
    BasisConversionPlanId fromPowerSums;
  };

// ============================================================================
// Basis-Conversion Plan Database
// ============================================================================
// Plan declarations record mathematical availability and any fixed piecewise
// formula policy. They do not choose among complete plans with the same
// endpoints.

  static BasisConversionPlanFormula useKernel(
      std::string name,
      BasisConversionKernel kernel,
      ExpressionCondition outputGuarantee = always());
  static BasisConversionPlanFormula composePlans(
      std::initializer_list<BasisConversionPlanId> plans);
  static const std::vector<BasisConversionPlanDefinition>&
  basisConversionPlanDatabase();
  static const std::vector<PowerSumFallbackPlanPair>&
  powerSumFallbackPlanDatabase();
  static const BasisConversionPlanId&
  genericPowerSumFallbackPlanIdentifier();
  static std::vector<BasisConversionPlanCase>
  powerSumsToSchurComponentFormulaCases();
  static const std::map<
      std::string,
      const BasisConversionPlanDefinition *>&
  basisConversionPlansById();
  static const std::map<
      std::pair<BasisKind, BasisKind>,
      std::vector<const BasisConversionPlanDefinition *>>&
  basisConversionPlansBySourceAndTarget();
  static const BasisConversionPlanDefinition *
  basisConversionPlanDefinition(
      const BasisConversionPlanId& id);
  static void validateBasisConversionPlanDatabase();
  const std::vector<RingBasisConversionPlan>&
  basisConversionPlansFor(
      int sourceBasisId,
      int targetBasisId) const;
  RingBasisConversionPlan basisConversionPlanForRing(
      const BasisConversionPlanId& id) const;
  RingBasisConversionPlan basisConversionPlanForRing(
      const BasisConversionPlanDefinition *definition) const;
  RingBasisConversionPlan basisConversionViaPowerSumsPlan(
      int sourceBasisId,
      int targetBasisId) const;
  mutable std::map<
      std::pair<int, int>,
      std::vector<RingBasisConversionPlan>>
      ringBasisConversionPlanCache;
  mutable std::map<
      std::pair<int, int>,
      RingBasisConversionPlan>
      powerSumFallbackPlanCache;

// ============================================================================
// Basis-Conversion Applicability And Selection
// ============================================================================
// Applicability is mathematical. Fixed plan cases may choose formulas for
// individual pieces by cost, while pickers own performance choices among
// complete plans.

  struct BasisConversionPreference
  {
    ExpressionCondition condition;
    BasisConversionPlanId plan;
    mutable const BasisConversionPlanDefinition *planDefinition =
        nullptr;
  };

  struct BasisConversionPicker
  {
    BasisKind sourceBasisKind;
    BasisKind targetBasisKind;
    std::vector<BasisConversionPreference> preferences;
    // Mathematically valid alternatives that remain available for forcing,
    // comparison, and testing but are not selected automatically.
    std::vector<BasisConversionPlanId> alternativePlans;
  };

  static const std::vector<BasisConversionPicker>&
  basisConversionPickerDatabase();
  static const std::map<
      std::pair<BasisKind, BasisKind>,
      const BasisConversionPicker *>&
  basisConversionPickersBySourceAndTarget();
  static const BasisConversionPicker *
  basisConversionPickerFor(
      BasisKind source,
      BasisKind target);
  static void validateBasisConversionPickerDatabase();
  bool basisConversionPlanApplicable(
      const RingBasisConversionPlan& plan,
      ring_elem expression,
      const ExpressionFacts& facts) const;
  RingBasisConversionPlan selectBasisConversionPlan(
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

  ring_elem executeBasisConversionPlanFormula(
      const BasisConversionPlanFormula& formula,
      int sourceBasisId,
      int targetBasisId,
      ring_elem expression,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts& inputFacts,
      ExpressionFacts *resultFacts = nullptr) const;
  ring_elem executeBasisConversionPlan(
      const RingBasisConversionPlan& plan,
      ring_elem expression,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts& inputFacts,
      ExpressionFacts *resultFacts = nullptr) const;
  // The development benchmark installs an explicit diagnostic context.
  // Ordinary conversion installs a non-diagnostic context so ambient
  // environment variables cannot alter its selection or timing.
  bool basisConversionContextActive() const;
  bool basisConversionTraceEnabled() const;
  ring_elem convertCanonicalExpressionToBasis(
      ring_elem f,
      int targetBasisId,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts *knownFacts = nullptr) const;

// ============================================================================
// Multiplication Plan Contracts
// ============================================================================
// A multiplication plan declares every operand conversion, the multiplication
// kernel, and the basis and canonical-form guarantee of the product. Operand
// conversions are fixed before multiplication; the owning workflow selects a
// support-dependent final conversion only after the product has been realized.

  enum class MultiplicationKernel
  {
    SchurFactorFormulas,
    MonomialExponentSplittings,
    ProductInSingleBasis
  };

  enum class MultiplicationPlanId
  {
    Unavailable,
    SchurLittlewoodRichardson,
    SchurHorizontalPieri,
    SchurVerticalPieri,
    SchurMurnaghanNakayama,
    SchurMixedPieriAndMurnaghanNakayama,
    MonomialLikeExponentSplittings,
    ViaHallLittlewoodGenerators,
    InTargetBasis,
    ViaPowerSums
  };

  struct MultiplicationPlan
  {
    MultiplicationPlanId id;
    MultiplicationKernel kernel;
    // A positive operand basis requests an explicit conversion before
    // multiplication. Structured kernels use -1 because they consume compatible
    // basis kinds directly.
    int leftOperandBasisId;
    int rightOperandBasisId;
    // Every plan declares the product basis and whether normalization can be
    // bypassed for that product.
    int productBasisId;
    bool productIsCanonical;
  };

  struct MultiplicationPlanSelection
  {
    MultiplicationPlan plan;
    std::optional<RingBasisConversionPlan> leftConversion;
    std::optional<RingBasisConversionPlan> rightConversion;
  };

// ============================================================================
// Multiplication Plans And Workflow Helpers
// ============================================================================

  const std::vector<MultiplicationPlan>&
  multiplicationPlansFor(
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId) const;
  std::vector<MultiplicationPlan> buildMultiplicationPlansFor(
      const ExpressionFacts& leftFacts,
      const ExpressionFacts& rightFacts,
      int targetBasisId) const;
  static bool hasSchurFactorFormula(BasisKind kind);
  static bool hasMonomialLikeFactorFormula(
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
      multiplicationPlansCache;

// ============================================================================
// Basis-Conversion Workflow Stages
// ============================================================================
// Preparation produces normalized, skew-free factors and exact facts. Product
// resolution then produces the product-free expression on which linear
// source-basis decomposition is mathematically valid.

  struct PreparedBasisConversionInput
  {
    ring_elem expression;
    ExpressionFacts facts;
    std::vector<size_t> factorsPerTerm;
    CombinatorialTags combinatorialTags = 0;
  };

  struct ProductFreeBasisConversionInput
  {
    ring_elem expression;
    std::optional<ExpressionFacts> exactFacts;
    // True means product resolution itself completed the requested conversion:
    // every surviving term is already a canonical target-basis term.
    bool productResolutionCompletedInTarget = false;
  };

  PreparedBasisConversionInput prepareBasisConversionInput(
      ring_elem expression,
      int targetBasisId) const;
  ProductFreeBasisConversionInput resolveProductsForBasisConversion(
      PreparedBasisConversionInput prepared,
      int targetBasisId) const;

// ============================================================================
// Public Conversion And Multiplication Entry Points
// ============================================================================
// Definitions follow dependency order: multiplication helpers precede
// toBasis because product resolution calls them.

 public:
  ring_elem multiplyToBasis(
      ring_elem f,
      ring_elem g,
      int targetBasisId) const;
  ring_elem toBasis(ring_elem f, int targetBasisId) const;
  // Development-only raw entry used by private M2 benchmark tooling.
  ring_elem toBasisBench(
      ring_elem f,
      int targetBasisId,
      const std::optional<std::string>& forcedPlanIdentifier,
      bool traceConversion) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
