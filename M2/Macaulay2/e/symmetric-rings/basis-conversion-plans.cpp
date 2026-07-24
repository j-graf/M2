// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include <utility>

namespace symmetric_rings {

// ============================================================================
// Complete Basis-Conversion Plans
// ============================================================================
// This is the single contributor-facing list of mathematically complete
// formulas. A plan uses one kernel, uses one named plan, composes named
// plans, or assigns fixed formulas to ordered expression pieces.

// ============================================================================
// Plan-Formula Constructors
// ============================================================================

SymmetricEngineRing::BasisConversionPlanFormula
SymmetricEngineRing::useKernel(
    std::string name,
    BasisConversionKernel kernel,
    ExpressionCondition outputGuarantee)
{
    return KernelPlanFormula{
        std::move(name),
        kernel,
        std::move(outputGuarantee)};
  }

SymmetricEngineRing::BasisConversionPlanFormula
SymmetricEngineRing::usePlan(
    BasisConversionPlanId plan)
{
    return NamedPlanFormula{
        std::move(plan)};
  }

SymmetricEngineRing::BasisConversionPlanFormula
SymmetricEngineRing::composePlans(
    std::initializer_list<BasisConversionPlanId> plans)
{
    if (plans.size() < 2)
      throw exc::engine_error(
          "a basis-conversion plan composition requires at least two "
          "named child plans");
    ComposedPlansFormula result;
    result.plans.assign(
        plans.begin(), plans.end());
    return result;
  }

// ============================================================================
// Canonical Power-Sum Fallback Plans
// ============================================================================
// The generic fallback P_{X,Y}^{(p)} needs one mathematically complete
// X -> p child and one p -> Y child. This table designates those children;
// it does not choose among them at runtime.

const std::vector<SymmetricEngineRing::PowerSumFallbackPlanPair>&
SymmetricEngineRing::powerSumFallbackPlanDatabase()
{
    static const std::vector<PowerSumFallbackPlanPair> fallbacks{
        {BasisKind::Complete,
         "Complete",
         {"Complete->PowerSum:Newton-identities"},
         {"PowerSum->Complete:logarithm-formula"}},
        {BasisKind::Elementary,
         "Elementary",
         {"Elementary->PowerSum:Newton-identities"},
         {"PowerSum->Elementary:logarithm-formula"}},
        {BasisKind::Schur,
         "Schur",
         {"Schur->PowerSum:Frobenius-character-formula"},
         {"PowerSum->Schur:homogeneous-component-formulas"}},
        {BasisKind::SchurOmega,
         "SchurOmega",
         {"SchurOmega->PowerSum:Frobenius-character-formula"},
         {"PowerSum->SchurOmega:homogeneous-component-formulas"}},
        {BasisKind::Monomial,
         "Monomial",
         {"Monomial->PowerSum:transition-matrix"},
         {"PowerSum->Monomial:transition-matrix"}},
        {BasisKind::Forgotten,
         "Forgotten",
         {"Forgotten->PowerSum:transition-matrix"},
         {"PowerSum->Forgotten:transition-matrix"}},
        {BasisKind::HallLittlewoodQGenerator,
         "HallLittlewoodQGenerator",
         {"HallLittlewoodQGenerator->PowerSum:generating-function"},
         {"PowerSum->HallLittlewoodQGenerator:logarithm-formula"}},
        {BasisKind::HallLittlewoodBGenerator,
         "HallLittlewoodBGenerator",
         {"HallLittlewoodBGenerator->PowerSum:generating-function"},
         {"PowerSum->HallLittlewoodBGenerator:logarithm-formula"}},
        {BasisKind::HallLittlewoodQ,
         "HallLittlewoodQ",
         {"HallLittlewoodQ->PowerSum:raising-operators"},
         {"PowerSum->HallLittlewoodQ:triangular-reduction"}},
        {BasisKind::HallLittlewoodB,
         "HallLittlewoodB",
         {"HallLittlewoodB->PowerSum:raising-operators"},
         {"PowerSum->HallLittlewoodB:triangular-reduction"}},
        {BasisKind::HallLittlewoodP,
         "HallLittlewoodP",
         {"HallLittlewoodP->PowerSum:via-Q-normalization"},
         {"PowerSum->HallLittlewoodP:triangular-reduction"}},
        {BasisKind::HallLittlewoodPOmega,
         "HallLittlewoodPOmega",
         {"HallLittlewoodPOmega->PowerSum:via-B-normalization"},
         {"PowerSum->HallLittlewoodPOmega:triangular-reduction"}}};
    return fallbacks;
  }

const SymmetricEngineRing::BasisConversionPlanId&
SymmetricEngineRing::genericPowerSumFallbackPlanIdentifier()
{
    static const BasisConversionPlanId id{
        "X->Y:via-power-sums"};
    return id;
  }

const std::vector<SymmetricEngineRing::BasisConversionPlanDefinition>&
SymmetricEngineRing::basisConversionPlanDatabase()
{
    static const std::vector<BasisConversionPlanDefinition> database = [] {
      std::vector<BasisConversionPlanDefinition> plans;
      auto addConversionPlan =
          [&](std::string_view id,
              BasisKind source,
              BasisKind target,
              BasisConversionKernel kernel,
              ExpressionCondition applicability = always()) {
            plans.push_back({
                {std::string(id)},
                source,
                target,
                std::move(applicability),
                ExpressionPieceKind::WholeExpression,
                {{otherwise(),
                  useKernel(
                      std::string(id), kernel)}}});
          };

      // ======================================================================
      // Conversion Plans To Power Sums
      // ======================================================================

      addConversionPlan(
          "Complete->PowerSum:Newton-identities",
          BasisKind::Complete,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              completeToPowerSumsViaNewtonIdentities);
      addConversionPlan(
          "Elementary->PowerSum:Newton-identities",
          BasisKind::Elementary,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              elementaryToPowerSumsViaNewtonIdentities);
      addConversionPlan(
          "Schur->PowerSum:Frobenius-character-formula",
          BasisKind::Schur,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              schurToPowerSumsViaFrobeniusCharacterFormula);
      addConversionPlan(
          "SchurOmega->PowerSum:Frobenius-character-formula",
          BasisKind::SchurOmega,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              schurOmegaToPowerSumsViaFrobeniusCharacterFormula);
      addConversionPlan(
          "Monomial->PowerSum:transition-matrix",
          BasisKind::Monomial,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              monomialToPowerSumsViaTransitionMatrix);
      addConversionPlan(
          "Forgotten->PowerSum:transition-matrix",
          BasisKind::Forgotten,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              forgottenToPowerSumsViaTransitionMatrix);
      addConversionPlan(
          "HallLittlewoodQGenerator->PowerSum:generating-function",
          BasisKind::HallLittlewoodQGenerator,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              hallLittlewoodQGeneratorsToPowerSumsViaGeneratingFunction);
      addConversionPlan(
          "HallLittlewoodBGenerator->PowerSum:generating-function",
          BasisKind::HallLittlewoodBGenerator,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              hallLittlewoodBGeneratorsToPowerSumsViaGeneratingFunction);
      addConversionPlan(
          "HallLittlewoodQ->PowerSum:raising-operators",
          BasisKind::HallLittlewoodQ,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              hallLittlewoodQToPowerSumsViaRaisingOperators);
      addConversionPlan(
          "HallLittlewoodB->PowerSum:raising-operators",
          BasisKind::HallLittlewoodB,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              hallLittlewoodBToPowerSumsViaRaisingOperators);
      addConversionPlan(
          "HallLittlewoodP->PowerSum:via-Q-normalization",
          BasisKind::HallLittlewoodP,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              hallLittlewoodPToPowerSumsViaQNormalization);
      addConversionPlan(
          "HallLittlewoodPOmega->PowerSum:via-B-normalization",
          BasisKind::HallLittlewoodPOmega,
          BasisKind::PowerSum,
          &SymmetricEngineRing::
              hallLittlewoodPOmegaToPowerSumsViaBNormalization);

      // ======================================================================
      // Conversion Plans From Power Sums
      // ======================================================================

      addConversionPlan(
          "PowerSum->Complete:logarithm-formula",
          BasisKind::PowerSum,
          BasisKind::Complete,
          &SymmetricEngineRing::
              powerSumsToCompleteViaLogarithmFormula);
      addConversionPlan(
          "PowerSum->Elementary:logarithm-formula",
          BasisKind::PowerSum,
          BasisKind::Elementary,
          &SymmetricEngineRing::
              powerSumsToElementaryViaLogarithmFormula);
      addConversionPlan(
          "PowerSum->Schur:abacus-rim-hooks",
          BasisKind::PowerSum,
          BasisKind::Schur,
          &SymmetricEngineRing::
              powerSumsToSchurViaAbacusRimHooks);
      addConversionPlan(
          "PowerSum->Schur:Frobenius-character-formula",
          BasisKind::PowerSum,
          BasisKind::Schur,
          &SymmetricEngineRing::
              powerSumsToSchurViaFrobeniusCharacterFormula);
      addConversionPlan(
          "PowerSum->Schur:Murnaghan-Nakayama",
          BasisKind::PowerSum,
          BasisKind::Schur,
          &SymmetricEngineRing::
              powerSumsToSchurViaMurnaghanNakayama);
      addConversionPlan(
          "PowerSum->SchurOmega:abacus-rim-hooks",
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          &SymmetricEngineRing::
              powerSumsToSchurOmegaViaAbacusRimHooks);
      addConversionPlan(
          "PowerSum->SchurOmega:Frobenius-character-formula",
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          &SymmetricEngineRing::
              powerSumsToSchurOmegaViaFrobeniusCharacterFormula);
      addConversionPlan(
          "PowerSum->SchurOmega:Murnaghan-Nakayama",
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          &SymmetricEngineRing::
              powerSumsToSchurOmegaViaMurnaghanNakayama);
      addConversionPlan(
          "PowerSum->HallLittlewoodQGenerator:logarithm-formula",
          BasisKind::PowerSum,
          BasisKind::HallLittlewoodQGenerator,
          &SymmetricEngineRing::
              powerSumsToHallLittlewoodGeneratorsViaLogarithmFormula);
      addConversionPlan(
          "PowerSum->HallLittlewoodBGenerator:logarithm-formula",
          BasisKind::PowerSum,
          BasisKind::HallLittlewoodBGenerator,
          &SymmetricEngineRing::
              powerSumsToHallLittlewoodGeneratorsViaLogarithmFormula);

      auto addPowerSumsToHallLittlewoodPlans =
          [&](BasisKind target, std::string_view targetName) {
            const std::string prefix =
                "PowerSum->" + std::string(targetName);
            addConversionPlan(
                std::string(prefix + ":single-cycles-via-Green-polynomials"),
                BasisKind::PowerSum,
                target,
                &SymmetricEngineRing::
                    powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials,
                allPowerSumTermsAreSingleCycles());
            addConversionPlan(
                std::string(prefix + ":triangular-reduction"),
                BasisKind::PowerSum,
                target,
                &SymmetricEngineRing::
                    powerSumsToHallLittlewoodViaTriangularReduction);
            addConversionPlan(
                std::string(prefix + ":Green-polynomials-and-duality"),
                BasisKind::PowerSum,
                target,
                &SymmetricEngineRing::
                    powerSumBasisElementToHallLittlewoodViaGreenPolynomialsAndDuality,
                expressionIsSingleBasisElement());
          };
      addPowerSumsToHallLittlewoodPlans(
          BasisKind::HallLittlewoodQ, "HallLittlewoodQ");
      addPowerSumsToHallLittlewoodPlans(
          BasisKind::HallLittlewoodB, "HallLittlewoodB");
      addPowerSumsToHallLittlewoodPlans(
          BasisKind::HallLittlewoodP, "HallLittlewoodP");
      addPowerSumsToHallLittlewoodPlans(
          BasisKind::HallLittlewoodPOmega,
          "HallLittlewoodPOmega");
      addConversionPlan(
          "PowerSum->Monomial:transition-matrix",
          BasisKind::PowerSum,
          BasisKind::Monomial,
          &SymmetricEngineRing::
              powerSumsToMonomialOrForgottenViaTransitionMatrix);
      addConversionPlan(
          "PowerSum->Forgotten:transition-matrix",
          BasisKind::PowerSum,
          BasisKind::Forgotten,
          &SymmetricEngineRing::
              powerSumsToMonomialOrForgottenViaTransitionMatrix);

      // ======================================================================
      // Direct Basis-Family Transitions
      // ======================================================================

      addConversionPlan(
          "HallLittlewoodQ->HallLittlewoodP:diagonal-scaling",
          BasisKind::HallLittlewoodQ,
          BasisKind::HallLittlewoodP,
          &SymmetricEngineRing::
              hallLittlewoodPairedBasesViaDiagonalScaling);
      addConversionPlan(
          "HallLittlewoodP->HallLittlewoodQ:diagonal-scaling",
          BasisKind::HallLittlewoodP,
          BasisKind::HallLittlewoodQ,
          &SymmetricEngineRing::
              hallLittlewoodPairedBasesViaDiagonalScaling);
      addConversionPlan(
          "HallLittlewoodB->HallLittlewoodPOmega:diagonal-scaling",
          BasisKind::HallLittlewoodB,
          BasisKind::HallLittlewoodPOmega,
          &SymmetricEngineRing::
              hallLittlewoodPairedBasesViaDiagonalScaling);
      addConversionPlan(
          "HallLittlewoodPOmega->HallLittlewoodB:diagonal-scaling",
          BasisKind::HallLittlewoodPOmega,
          BasisKind::HallLittlewoodB,
          &SymmetricEngineRing::
              hallLittlewoodPairedBasesViaDiagonalScaling);
      addConversionPlan(
          "Schur->SchurOmega:partition-conjugation",
          BasisKind::Schur,
          BasisKind::SchurOmega,
          &SymmetricEngineRing::
              schurAndSchurOmegaViaPartitionConjugation);
      addConversionPlan(
          "SchurOmega->Schur:partition-conjugation",
          BasisKind::SchurOmega,
          BasisKind::Schur,
          &SymmetricEngineRing::
              schurAndSchurOmegaViaPartitionConjugation);
      addConversionPlan(
          "Schur->Complete:Jacobi-Trudi",
          BasisKind::Schur,
          BasisKind::Complete,
          &SymmetricEngineRing::
              schurAndSchurOmegaToGeneratorsViaJacobiTrudi);
      addConversionPlan(
          "SchurOmega->Elementary:Jacobi-Trudi",
          BasisKind::SchurOmega,
          BasisKind::Elementary,
          &SymmetricEngineRing::
              schurAndSchurOmegaToGeneratorsViaJacobiTrudi);
      addConversionPlan(
          "Complete->Schur:horizontal-Pieri",
          BasisKind::Complete,
          BasisKind::Schur,
          &SymmetricEngineRing::
              completeToSchurViaHorizontalPieri);
      addConversionPlan(
          "HallLittlewoodQGenerator->HallLittlewoodQ:"
          "triangular-reduction",
          BasisKind::HallLittlewoodQGenerator,
          BasisKind::HallLittlewoodQ,
          &SymmetricEngineRing::
              hallLittlewoodGeneratorsToCapitalBasesViaTriangularReduction);
      addConversionPlan(
          "HallLittlewoodQGenerator->HallLittlewoodP:"
          "triangular-reduction",
          BasisKind::HallLittlewoodQGenerator,
          BasisKind::HallLittlewoodP,
          &SymmetricEngineRing::
              hallLittlewoodGeneratorsToCapitalBasesViaTriangularReduction);
      addConversionPlan(
          "HallLittlewoodBGenerator->HallLittlewoodB:"
          "triangular-reduction",
          BasisKind::HallLittlewoodBGenerator,
          BasisKind::HallLittlewoodB,
          &SymmetricEngineRing::
              hallLittlewoodGeneratorsToCapitalBasesViaTriangularReduction);
      addConversionPlan(
          "HallLittlewoodBGenerator->HallLittlewoodPOmega:"
          "triangular-reduction",
          BasisKind::HallLittlewoodBGenerator,
          BasisKind::HallLittlewoodPOmega,
          &SymmetricEngineRing::
              hallLittlewoodGeneratorsToCapitalBasesViaTriangularReduction);

      // ======================================================================
      // Power Sums To Schur Piecewise Plans
      // ======================================================================
      // The complete formula applies independently to each homogeneous
      // component of a nonhomogeneous power-sum expansion. Its ordered cases
      // fix one formula for every component, so execution makes no additional
      // plan choice.
      auto powerSumsToSchurHomogeneousCases = [&] {
            const auto viaCompleteBasis =
                composePlans({
                    {"PowerSum->Complete:logarithm-formula"},
                    {"Complete->Schur:horizontal-Pieri"}});
            const auto viaAbacusRimHooks =
                useKernel(
                    "PowerSum->Schur:abacus-rim-hooks",
                    &SymmetricEngineRing::
                        powerSumsToSchurViaAbacusRimHooks);
            const auto shortCycleHybrid =
                usePlan({
                    "PowerSum->Schur:"
                    "short-cycle-hybrid"});
            const uint32_t plethysmTag =
                combinatorialTagMask(
                    CombinatorialTag::Plethysm);
            const uint32_t littlewoodRichardsonTag =
                combinatorialTagMask(
                    CombinatorialTag::
                        LittlewoodRichardson);
            const uint32_t horizontalPieriTag =
                combinatorialTagMask(
                    CombinatorialTag::HorizontalPieri);
            const uint32_t verticalPieriTag =
                combinatorialTagMask(
                    CombinatorialTag::VerticalPieri);
            const uint32_t borderStripsTag =
                combinatorialTagMask(
                    CombinatorialTag::BorderStrips);

            const auto supportIsDenseForCompleteConversion =
                componentSupportSquareRatioAtLeast(35, 2);
            const auto useCompleteForOrdinaryComponent =
                !coefficientRingIsQQ() ||
                componentDensityAtLeast(1, 4) ||
                supportIsDenseForCompleteConversion;
            const auto useCompleteForPlethysm =
                hasCombinatorialTag(plethysmTag) &&
                (supportIsDenseForCompleteConversion ||
                 (componentWeightGreaterThan(7) &&
                  componentTermCountAtLeast(2)));
            const auto useCompleteForLittlewoodRichardson =
                !hasCombinatorialTag(plethysmTag) &&
                hasCombinatorialTag(
                    littlewoodRichardsonTag) &&
                componentWeightGreaterThan(6) &&
                useCompleteForOrdinaryComponent;
            const auto useCompleteForPieri =
                !hasCombinatorialTag(plethysmTag) &&
                !hasCombinatorialTag(
                    littlewoodRichardsonTag) &&
                (hasCombinatorialTag(horizontalPieriTag) ||
                 hasCombinatorialTag(verticalPieriTag)) &&
                componentWeightGreaterThan(7) &&
                useCompleteForOrdinaryComponent;
            const auto hasSmallCommonCycle =
                (coefficientRingIsQQ() &&
                 componentHasCommonPowerSumPartAtMostPercentOfWeight(36)) ||
                (!coefficientRingIsQQ() &&
                 componentHasCommonPowerSumPartAtMostPercentOfWeight(41));
            const auto useCompleteForLargeBorderStrip =
                combinatorialTagsEqual(borderStripsTag) &&
                componentWeightGreaterThan(15) &&
                componentTermCountAtLeast(8) &&
                !componentHasCommonPowerSumPartOne() &&
                hasSmallCommonCycle;
            const auto useCompleteForBorderStripWithCommonOne =
                combinatorialTagsEqual(borderStripsTag) &&
                componentWeightGreaterThan(7) &&
                componentTermCountAtLeast(8) &&
                componentHasCommonPowerSumPartOne() &&
                useCompleteForOrdinaryComponent;
            const auto ordinaryOrBorderStrip =
                combinatorialTagsEqual(0) ||
                combinatorialTagsEqual(borderStripsTag);
            const auto useCompleteForMostlyShortCycleTerms =
                ordinaryOrBorderStrip &&
                componentWeightGreaterThan(13) &&
                componentTermCountAtLeast(2) &&
                componentAllPowerSumIndicesHaveMostlyShortCycles();
            const auto useTermwiseHybrid =
                ordinaryOrBorderStrip &&
                componentWeightGreaterThan(13) &&
                componentTermCountAtLeast(8) &&
                componentMostlyShortCycleFractionAtLeast(
                    1, 4, 8) &&
                componentHasBothMostlyShortCycleAndOtherTerms();

            return std::vector<BasisConversionPlanCase>{
                {useCompleteForPlethysm,
                 viaCompleteBasis},
                {useCompleteForLittlewoodRichardson,
                 viaCompleteBasis},
                {useCompleteForPieri,
                 viaCompleteBasis},
                {useCompleteForLargeBorderStrip,
                 viaCompleteBasis},
                {useCompleteForBorderStripWithCommonOne,
                 viaCompleteBasis},
                {useCompleteForMostlyShortCycleTerms,
                 viaCompleteBasis},
                // Use the fixed term-level hybrid for this component. That
                // named plan partitions its realized input without
                // re-entering the picker.
                {useTermwiseHybrid, shortCycleHybrid},
                {otherwise(), viaAbacusRimHooks}};
          };
      plans.push_back({
          {"PowerSum->Schur:homogeneous-component-formulas"},
          BasisKind::PowerSum,
          BasisKind::Schur,
          always(),
          ExpressionPieceKind::HomogeneousComponents,
          powerSumsToSchurHomogeneousCases()});

      // ======================================================================
      // Named Schur Compositions And Hybrids
      // ======================================================================
      // Named compositions and hybrids use the same plan representation as
      // direct kernels. Child identifiers completely determine execution.
      plans.push_back({
          {"PowerSum->Schur:via-complete-basis"},
          BasisKind::PowerSum,
          BasisKind::Schur,
          always(),
          ExpressionPieceKind::WholeExpression,
          {{otherwise(),
            composePlans({
                {"PowerSum->Complete:logarithm-formula"},
                {"Complete->Schur:horizontal-Pieri"}})}}});
      plans.push_back({
          {"PowerSum->Schur:short-cycle-hybrid"},
          BasisKind::PowerSum,
          BasisKind::Schur,
          always(),
          ExpressionPieceKind::IndividualTerms,
          {{powerSumIndexHasMostlyShortCycles(),
            composePlans({
                {"PowerSum->Complete:logarithm-formula"},
                {"Complete->Schur:horizontal-Pieri"}})},
           {otherwise(),
            useKernel(
                "PowerSum->Schur:abacus-rim-hooks",
                &SymmetricEngineRing::
                    powerSumsToSchurViaAbacusRimHooks)}}});
      plans.push_back({
          {"PowerSum->SchurOmega:homogeneous-component-formulas"},
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          always(),
          ExpressionPieceKind::WholeExpression,
          {{otherwise(),
            composePlans({
                {"PowerSum->Schur:homogeneous-component-formulas"},
                {"Schur->SchurOmega:partition-conjugation"}})}}});
      plans.push_back({
          {"PowerSum->SchurOmega:via-complete-basis"},
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          always(),
          ExpressionPieceKind::WholeExpression,
          {{otherwise(),
            composePlans({
                {"PowerSum->Schur:via-complete-basis"},
                {"Schur->SchurOmega:partition-conjugation"}})}}});
      plans.push_back({
          {"PowerSum->SchurOmega:short-cycle-hybrid"},
          BasisKind::PowerSum,
          BasisKind::SchurOmega,
          always(),
          ExpressionPieceKind::WholeExpression,
          {{otherwise(),
            composePlans({
                {"PowerSum->Schur:short-cycle-hybrid"},
                {"Schur->SchurOmega:partition-conjugation"}})}}});

      return plans;
    }();
    return database;
  }

// ============================================================================
// Generic Power-Sum Composition Plan
// ============================================================================
// For non-power-sum endpoints this instantiates the single mathematical plan
//
//             P_{X,Y}^{(p)} = P_{p,Y} o P_{X,p}.
//
// The result is an ordinary immutable plan definition with both named plans
// identified. The generic executor therefore treats it exactly like a
// composition written explicitly in the database.

SymmetricEngineRing::RingBasisConversionPlan
SymmetricEngineRing::basisConversionViaPowerSumsPlan(
    int sourceBasisId,
    int targetBasisId) const
{
    requireBasis(sourceBasisId);
    requireBasis(targetBasisId);
    if (error()) return {};
    const BasisKind sourceKind =
        basisKindForId(sourceBasisId);
    const BasisKind targetKind =
        basisKindForId(targetBasisId);
    if (sourceKind == targetKind ||
        sourceKind == BasisKind::Custom ||
        targetKind == BasisKind::Custom ||
        sourceKind == BasisKind::PowerSum ||
        targetKind == BasisKind::PowerSum)
      return {};

    const std::pair<int, int> endpoint{
        sourceBasisId, targetBasisId};
    auto cached =
        powerSumFallbackPlanCache.find(
            endpoint);
    if (cached !=
        powerSumFallbackPlanCache.end())
      return cached->second;

    validateBasisConversionPlanDatabase();
    const PowerSumFallbackPlanPair *sourceFallback = nullptr;
    const PowerSumFallbackPlanPair *targetFallback = nullptr;
    for (const auto& fallback :
         powerSumFallbackPlanDatabase())
      {
        if (fallback.basisKind == sourceKind)
          sourceFallback = &fallback;
        if (fallback.basisKind == targetKind)
          targetFallback = &fallback;
      }
    if (sourceFallback == nullptr ||
        targetFallback == nullptr)
      {
        ERROR("the generic power-sum composition has no fallback "
              "for its endpoints");
        return {};
      }

    const auto *toPowerSums =
        basisConversionPlanDefinition(
            sourceFallback->toPowerSums);
    const auto *fromPowerSums =
        basisConversionPlanDefinition(
            targetFallback->fromPowerSums);
    if (toPowerSums == nullptr ||
        fromPowerSums == nullptr)
      {
        ERROR("the generic power-sum composition has missing "
              "named plans");
        return {};
      }
    BasisConversionPlanFormula formula =
        composePlans({
            sourceFallback->toPowerSums,
            targetFallback->fromPowerSums});
    auto& composition =
        std::get<ComposedPlansFormula>(
            formula);
    composition.planDefinitions = {
        toPowerSums, fromPowerSums};
    auto instantiatedDefinition =
        std::make_shared<BasisConversionPlanDefinition>(
            BasisConversionPlanDefinition{
                genericPowerSumFallbackPlanIdentifier(),
                sourceKind,
                targetKind,
                always(),
                ExpressionPieceKind::WholeExpression,
                {{otherwise(), std::move(formula)}}});
    RingBasisConversionPlan result{
        instantiatedDefinition.get(),
        sourceBasisId,
        targetBasisId,
        std::move(instantiatedDefinition)};
    auto inserted =
        powerSumFallbackPlanCache.emplace(
            endpoint, std::move(result));
    return inserted.first->second;
  }


} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
