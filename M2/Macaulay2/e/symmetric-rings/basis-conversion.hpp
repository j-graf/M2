// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_CONVERSION_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_CONVERSION_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Expression Preparation And Metadata
// ============================================================================
// Preparation establishes the canonical, skew-free term representation
// consumed by both the conversion plan database and multiplication workflows.

 private:
  struct PowerSumSupportFacts
  {
    std::optional<bool> allTermsAreSingleCycles;
    std::optional<size_t> mostlyShortCycleTermCount;
    std::optional<std::vector<int>> commonParts;
  };
  PowerSumSupportFacts inspectPowerSumSupportFacts(
      ring_elem expression,
      const std::vector<size_t>& termPositions,
      ExpressionFactRequirements requirements) const;
  std::vector<ExpressionPieceFacts>
  inspectExpressionPieces(
      ring_elem expression,
      ExpressionPieceKind pieceKind,
      const ExpressionFacts& facts,
      ExpressionFactRequirements requirements =
          noExpressionFactRequirements) const;

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
    std::string identifier;
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

  struct ResolvedBasisConversionPlan
  {
    const BasisConversionPlanDefinition *definition = nullptr;
    int sourceBasisId = -1;
    int targetBasisId = -1;
    // Generic compositions are instantiated as ordinary immutable plan
    // definitions. Shared ownership keeps that definition stable while the
    // ring-specific plan is copied through selection and execution.
    std::shared_ptr<const BasisConversionPlanDefinition>
        instantiatedDefinition;

    ResolvedBasisConversionPlan() = default;
    ResolvedBasisConversionPlan(
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
      std::string identifier,
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
  const std::vector<ResolvedBasisConversionPlan>&
  basisConversionPlansFor(
      int sourceBasisId,
      int targetBasisId) const;
  ResolvedBasisConversionPlan basisConversionPlanForRing(
      const BasisConversionPlanId& id) const;
  ResolvedBasisConversionPlan basisConversionPlanForRing(
      const BasisConversionPlanDefinition *definition) const;
  ResolvedBasisConversionPlan basisConversionViaPowerSumsPlan(
      int sourceBasisId,
      int targetBasisId) const;
  mutable std::map<
      std::pair<int, int>,
      std::vector<ResolvedBasisConversionPlan>>
      ringBasisConversionPlanCache;
  mutable std::map<
      std::pair<int, int>,
      ResolvedBasisConversionPlan>
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
      const ResolvedBasisConversionPlan& plan,
      ring_elem expression,
      const ExpressionFacts& facts) const;
  ResolvedBasisConversionPlan selectBasisConversionPlan(
      ring_elem expression,
      int sourceBasisId,
      int targetBasisId,
      const ExpressionFacts& facts,
      CombinatorialTags combinatorialTags,
      const std::optional<std::string>& forcedIdentifier =
          std::nullopt) const;

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
      const ResolvedBasisConversionPlan& plan,
      ring_elem expression,
      CombinatorialTags combinatorialTags,
      const ExpressionFacts& inputFacts,
      ExpressionFacts *resultFacts = nullptr) const;
  ring_elem executeNamedBasisConversionPlan(
      const BasisConversionPlanId& planIdentifier,
      ring_elem expression,
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
// Public Conversion Entry Points
// ============================================================================
// Shared normalization and term-local product resolution are declared in
// expression-helpers.hpp and consumed by toBasis before plan selection.

 public:
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
