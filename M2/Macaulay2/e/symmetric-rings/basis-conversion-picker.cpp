// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include <utility>

namespace symmetric_rings {

// ============================================================================
// Performance Choices Among Complete Plans
// ============================================================================
// Each source-to-target entry is an ordered condition-to-plan list. Other
// valid complete plans are listed explicitly as non-automatic alternatives.
//
// If an endpoint has no specific automatic plan, selection uses the single
// parameterized u -> p -> v plan. Its two named plans are designated in the
// plan file and identified before execution. Endpoint pickers below therefore
// describe only genuine choices among specific mathematical formulas.

const std::vector<
    SymmetricEngineRing::BasisConversionPicker>&
SymmetricEngineRing::basisConversionPickerDatabase()
{
    static const std::vector<BasisConversionPicker>
        database = [] {
          std::vector<BasisConversionPicker> pickers;
          auto addPreferredConversionPlan =
              [&](BasisKind source,
                  BasisKind target,
                  std::string preferredPlan) {
                BasisConversionPicker picker{
                    source,
                    target,
                    {{otherwise(), {std::move(preferredPlan)}}},
                    {genericPowerSumFallbackPlanIdentifier()}};
                pickers.push_back(std::move(picker));
              };

          // ==================================================================
          // Power Sums To Schur And Schur Omega
          // ==================================================================
          // The homogeneous-component plan is the automatic choice. The other
          // complete formulas remain available by name for mathematical
          // comparison, tests, and benchmarks.
          pickers.push_back({
              BasisKind::PowerSum,
              BasisKind::Schur,
              {{otherwise(),
                {"PowerSum->Schur:homogeneous-component-formulas"}}},
              {{"PowerSum->Schur:abacus-rim-hooks"},
               {"PowerSum->Schur:Frobenius-character-formula"},
               {"PowerSum->Schur:Murnaghan-Nakayama"},
               {"PowerSum->Schur:via-complete-basis"},
               {"PowerSum->Schur:short-cycle-hybrid"}}});
          pickers.push_back({
              BasisKind::PowerSum,
              BasisKind::SchurOmega,
              {{otherwise(),
                {"PowerSum->SchurOmega:homogeneous-component-formulas"}}},
              {{"PowerSum->SchurOmega:abacus-rim-hooks"},
               {"PowerSum->SchurOmega:Frobenius-character-formula"},
               {"PowerSum->SchurOmega:Murnaghan-Nakayama"},
               {"PowerSum->SchurOmega:via-complete-basis"},
               {"PowerSum->SchurOmega:"
                "short-cycle-hybrid"}}});

          // ==================================================================
          // Power Sums To Hall--Littlewood Capital Bases
          // ==================================================================
          // These targets use the same ordered mathematical cases with
          // target-specific complete plan identifiers.
          auto addPowerSumsToHallLittlewoodPicker =
              [&](BasisKind target,
                  std::string_view targetName) {
                const std::string prefix =
                    "PowerSum->" +
                    std::string(targetName);
                pickers.push_back({
                    BasisKind::PowerSum,
                    target,
                    {{allPowerSumTermsAreSingleCycles(),
                      {prefix + ":single-cycles-via-Green-polynomials"}},
                     {expressionIsSingleBasisElement(),
                      {prefix + ":Green-polynomials-and-duality"}},
                     {otherwise(),
                      {prefix + ":triangular-reduction"}}},
                    {}});
              };
          addPowerSumsToHallLittlewoodPicker(
              BasisKind::HallLittlewoodQ,
              "HallLittlewoodQ");
          addPowerSumsToHallLittlewoodPicker(
              BasisKind::HallLittlewoodB,
              "HallLittlewoodB");
          addPowerSumsToHallLittlewoodPicker(
              BasisKind::HallLittlewoodP,
              "HallLittlewoodP");
          addPowerSumsToHallLittlewoodPicker(
              BasisKind::HallLittlewoodPOmega,
              "HallLittlewoodPOmega");

          // ==================================================================
          // Direct Family Transitions Versus Broad Fallbacks
          // ==================================================================
          // Direct involutions, normalizations, Jacobi--Trudi formulas, and
          // triangular reductions always outrank their broad p fallbacks.
          addPreferredConversionPlan(
              BasisKind::Complete,
              BasisKind::Schur,
              "Complete->Schur:horizontal-Pieri");
          addPreferredConversionPlan(
              BasisKind::Schur,
              BasisKind::Complete,
              "Schur->Complete:Jacobi-Trudi");
          addPreferredConversionPlan(
              BasisKind::SchurOmega,
              BasisKind::Elementary,
              "SchurOmega->Elementary:Jacobi-Trudi");
          addPreferredConversionPlan(
              BasisKind::Schur,
              BasisKind::SchurOmega,
              "Schur->SchurOmega:partition-conjugation");
          addPreferredConversionPlan(
              BasisKind::SchurOmega,
              BasisKind::Schur,
              "SchurOmega->Schur:partition-conjugation");
          addPreferredConversionPlan(
              BasisKind::HallLittlewoodQ,
              BasisKind::HallLittlewoodP,
              "HallLittlewoodQ->HallLittlewoodP:diagonal-scaling");
          addPreferredConversionPlan(
              BasisKind::HallLittlewoodP,
              BasisKind::HallLittlewoodQ,
              "HallLittlewoodP->HallLittlewoodQ:diagonal-scaling");
          addPreferredConversionPlan(
              BasisKind::HallLittlewoodB,
              BasisKind::HallLittlewoodPOmega,
              "HallLittlewoodB->HallLittlewoodPOmega:"
              "diagonal-scaling");
          addPreferredConversionPlan(
              BasisKind::HallLittlewoodPOmega,
              BasisKind::HallLittlewoodB,
              "HallLittlewoodPOmega->HallLittlewoodB:"
              "diagonal-scaling");
          addPreferredConversionPlan(
              BasisKind::HallLittlewoodQGenerator,
              BasisKind::HallLittlewoodQ,
              "HallLittlewoodQGenerator->HallLittlewoodQ:"
              "triangular-reduction");
          addPreferredConversionPlan(
              BasisKind::HallLittlewoodQGenerator,
              BasisKind::HallLittlewoodP,
              "HallLittlewoodQGenerator->HallLittlewoodP:"
              "triangular-reduction");
          addPreferredConversionPlan(
              BasisKind::HallLittlewoodBGenerator,
              BasisKind::HallLittlewoodB,
              "HallLittlewoodBGenerator->HallLittlewoodB:"
              "triangular-reduction");
          addPreferredConversionPlan(
              BasisKind::HallLittlewoodBGenerator,
              BasisKind::HallLittlewoodPOmega,
              "HallLittlewoodBGenerator->HallLittlewoodPOmega:"
              "triangular-reduction");
          return pickers;
        }();
    return database;
  }

// ============================================================================
// Picker Indexing And Structural Validation
// ============================================================================

const std::map<
    std::pair<
        SymmetricEngineRing::BasisKind,
        SymmetricEngineRing::BasisKind>,
    const SymmetricEngineRing::BasisConversionPicker *>&
SymmetricEngineRing::basisConversionPickersBySourceAndTarget()
{
    static const std::map<
        std::pair<BasisKind, BasisKind>,
        const BasisConversionPicker *> index = [] {
          std::map<
              std::pair<BasisKind, BasisKind>,
              const BasisConversionPicker *> result;
          for (const auto& picker :
               basisConversionPickerDatabase())
            if (!result.emplace(
                     std::make_pair(
                         picker.sourceBasisKind,
                         picker.targetBasisKind),
                     &picker).second)
              throw exc::engine_error(
                  "basis-conversion picker endpoints must be unique");
          return result;
        }();
    return index;
  }

const SymmetricEngineRing::BasisConversionPicker *
SymmetricEngineRing::basisConversionPickerFor(
    BasisKind source,
    BasisKind target)
{
    const auto& index =
        basisConversionPickersBySourceAndTarget();
    auto found = index.find({source, target});
    return found == index.end() ? nullptr : found->second;
  }

void SymmetricEngineRing::validateBasisConversionPickerDatabase()
{
    static const bool validated = [] {
      validateBasisConversionPlanDatabase();
      const auto& plansById =
          basisConversionPlansById();
      const auto& plansBySourceAndTarget =
          basisConversionPlansBySourceAndTarget();
      const auto& pickersBySourceAndTarget =
          basisConversionPickersBySourceAndTarget();
      auto powerSumFallbackAvailable =
          [](BasisKind source, BasisKind target) {
            return source != target &&
                   source != BasisKind::Custom &&
                   target != BasisKind::Custom &&
                   source != BasisKind::PowerSum &&
                   target != BasisKind::PowerSum;
          };

      for (const auto& picker :
           basisConversionPickerDatabase())
        {
          if (picker.preferences.empty())
            throw exc::engine_error(
                "a basis-conversion picker has no preferences");
          std::map<std::string, bool> accountedFor;
          for (size_t i = 0; i < picker.preferences.size(); ++i)
            {
              const auto& item = picker.preferences[i];
              const bool isOtherwise =
                  item.condition.kind ==
                  ExpressionConditionKind::Otherwise;
              if (isOtherwise !=
                  (i + 1 == picker.preferences.size()))
                throw exc::engine_error(
                    "a basis-conversion picker must end in exactly "
                    "one otherwise case");
              validateExpressionCondition(
                  item.condition,
                  ExpressionPieceKind::WholeExpression,
                  isOtherwise);
              if (!accountedFor.emplace(
                       item.plan.value, true).second)
                throw exc::engine_error(
                    "a basis-conversion picker references one plan "
                    "more than once");
              if (item.plan.value ==
                  genericPowerSumFallbackPlanIdentifier().
                      value)
                {
                  if (!powerSumFallbackAvailable(
                          picker.sourceBasisKind,
                          picker.targetBasisKind))
                    throw exc::engine_error(
                        "a basis-conversion picker references the "
                        "generic power-sum plan at incompatible endpoints");
                  item.planDefinition = nullptr;
                }
              else
                {
                  auto found =
                      plansById.find(item.plan.value);
                  if (found == plansById.end() ||
                      found->second->sourceBasisKind !=
                          picker.sourceBasisKind ||
                      found->second->targetBasisKind !=
                          picker.targetBasisKind)
                    throw exc::engine_error(
                        "a basis-conversion picker references a plan "
                        "with missing or incompatible endpoints");
                  item.planDefinition = found->second;
                  if (isOtherwise &&
                      found->second->applicability.kind !=
                          ExpressionConditionKind::Always)
                    throw exc::engine_error(
                        "a basis-conversion picker's otherwise case "
                        "must name an unconditional broad plan");
                }
            }
          for (const auto& id :
               picker.alternativePlans)
            {
              auto found = plansById.find(id.value);
              const bool usesPowerSumFallback =
                  id.value ==
                  genericPowerSumFallbackPlanIdentifier().
                      value;
              const bool validPowerSumFallback =
                  usesPowerSumFallback &&
                  powerSumFallbackAvailable(
                      picker.sourceBasisKind,
                      picker.targetBasisKind);
              const bool validSpecific =
                  !usesPowerSumFallback &&
                  found != plansById.end() &&
                  found->second->sourceBasisKind ==
                      picker.sourceBasisKind &&
                  found->second->targetBasisKind ==
                      picker.targetBasisKind;
              if ((!validPowerSumFallback && !validSpecific) ||
                  !accountedFor.emplace(
                      id.value, true).second)
                throw exc::engine_error(
                    "invalid alternative basis-conversion plan");
            }
          const auto endpoint =
              std::make_pair(
                  picker.sourceBasisKind,
                  picker.targetBasisKind);
          auto endpointPlans =
              plansBySourceAndTarget.find(endpoint);
          const size_t specificPlanCount =
              endpointPlans == plansBySourceAndTarget.end()
                  ? 0
                  : endpointPlans->second.size();
          const size_t completePlanCount =
              specificPlanCount +
              (powerSumFallbackAvailable(
                   picker.sourceBasisKind,
                   picker.targetBasisKind)
                   ? 1
                   : 0);
          if (accountedFor.size() != completePlanCount)
            throw exc::engine_error(
                "a basis-conversion picker does not account for "
                "every plan at its endpoints");
        }

      for (const auto& endpoint : plansBySourceAndTarget)
        {
          const bool needsPicker =
              endpoint.second.size() +
                  (powerSumFallbackAvailable(
                       endpoint.first.first,
                       endpoint.first.second)
                       ? 1
                       : 0) >
              1;
          const bool hasPicker =
              pickersBySourceAndTarget.find(endpoint.first) !=
              pickersBySourceAndTarget.end();
          if (needsPicker != hasPicker)
            throw exc::engine_error(
                needsPicker
                    ? "multiple basis-conversion plans require an "
                      "explicit picker"
                    : "a sole basis-conversion plan must not have a "
                      "redundant picker");
        }
      return true;
    }();
    (void) validated;
  }

// ============================================================================
// Complete-Plan Applicability And Performance Policy
// ============================================================================

bool SymmetricEngineRing::basisConversionPlanApplicable(
    const ResolvedBasisConversionPlan& plan,
    ring_elem expression,
    const ExpressionFacts& facts) const
{
    if (!plan.valid() ||
        !facts.canonicalExpansionInBasis(
            plan.sourceBasisId))
      return false;
    if (plan.definition->applicability.kind ==
        ExpressionConditionKind::Always)
      return true;
    const auto pieceFacts =
        inspectExpressionPieces(
            expression,
            ExpressionPieceKind::WholeExpression,
            facts,
            expressionFactRequirements(
                plan.definition->applicability));
    return pieceFacts.size() == 1 &&
           expressionConditionHolds(
               plan.definition->applicability,
               pieceFacts.front());
  }

SymmetricEngineRing::ResolvedBasisConversionPlan
SymmetricEngineRing::selectBasisConversionPlan(
    ring_elem expression,
    int sourceBasisId,
    int targetBasisId,
    const ExpressionFacts& facts,
    CombinatorialTags combinatorialTags,
    const std::optional<std::string>& forcedIdentifier) const
{
    if (basisConversionPlanExecutionDepth != 0)
      {
        ERROR("basis-conversion plan execution attempted to invoke "
              "the performance picker");
        return {};
      }
    if (!facts.canonicalExpansionInBasis(sourceBasisId))
      {
        ERROR("the basis-conversion picker requires a canonical "
              "single-source-basis expansion");
        return {};
      }
    validateBasisConversionPickerDatabase();
    const auto& plans =
        basisConversionPlansFor(
            sourceBasisId, targetBasisId);
    ResolvedBasisConversionPlan powerSumFallback =
        basisConversionViaPowerSumsPlan(
            sourceBasisId, targetBasisId);
    if (error()) return {};
    const BasisKind sourceKind =
        basisKindForId(sourceBasisId);
    const BasisKind targetKind =
        basisKindForId(targetBasisId);
    const std::optional<std::string>& requested =
        forcedIdentifier;
    ExpressionFacts selectionProfile = facts;
    selectionProfile.combinatorialTags = combinatorialTags;

    if (!requested &&
        plans.empty() &&
        powerSumFallback.valid())
      {
        return powerSumFallback;
      }
    // Most remaining endpoints have one unconditional specific mathematical
    // plan. Return it without allocating a profile or running general
    // conditional policy.
    if (!requested &&
        plans.size() == 1 &&
        !powerSumFallback.valid() &&
        plans.front().definition->
                applicability.kind ==
            ExpressionConditionKind::Always)
      {
        return plans.front();
      }

    if (requested)
      {
        if (powerSumFallback.valid() &&
            *requested ==
                genericPowerSumFallbackPlanIdentifier().
                    value)
          {
            return powerSumFallback;
          }
        for (const auto& plan : plans)
          if (plan.definition->id.value == *requested)
            {
              if (basisConversionPlanApplicable(
                      plan, expression, selectionProfile))
                return plan;
              break;
            }
        ERROR("the requested conversion plan is unknown or "
              "inapplicable: ",
              requested->c_str());
        return {};
      }

    const auto *picker =
        basisConversionPickerFor(
            sourceKind, targetKind);
    if (picker == nullptr)
      {
        ERROR("multiple basis-conversion plans reached selection "
              "without a declarative endpoint picker");
        return {};
      }

    bool needsProfile = false;
    for (const auto& item : picker->preferences)
      if (item.condition.kind !=
          ExpressionConditionKind::Otherwise)
        {
          needsProfile = true;
          break;
        }
    std::vector<ExpressionPieceFacts> pieceFacts;
    if (needsProfile)
      {
        std::vector<ExpressionCondition> preferenceConditions;
        preferenceConditions.reserve(
            picker->preferences.size());
        for (const auto& item : picker->preferences)
          preferenceConditions.push_back(item.condition);
        pieceFacts = inspectExpressionPieces(
            expression,
            ExpressionPieceKind::WholeExpression,
            selectionProfile,
            expressionFactRequirements(
                preferenceConditions));
        if (pieceFacts.size() != 1)
          {
            ERROR("a basis-conversion picker could not construct "
                  "one whole-expression facts record");
            return {};
          }
      }

    for (const auto& item : picker->preferences)
      {
        const bool preferred =
            item.condition.kind ==
                ExpressionConditionKind::Otherwise ||
            (needsProfile &&
             expressionConditionHolds(
                 item.condition, pieceFacts.front()));
        if (!preferred)
          continue;
        const bool usesPowerSumFallback =
            item.plan.value ==
            genericPowerSumFallbackPlanIdentifier().
                value;
        if (usesPowerSumFallback)
          {
            if (powerSumFallback.valid())
              return powerSumFallback;
            ERROR("a basis-conversion picker selected the power-sum "
                  "fallback at incompatible endpoints");
            return {};
          }
        if (item.planDefinition == nullptr)
          {
            ERROR("a basis-conversion picker reached selection "
                  "before its plan definitions were identified");
            return {};
          }
        ResolvedBasisConversionPlan selected{
            item.planDefinition,
            sourceBasisId,
            targetBasisId};
        if (basisConversionPlanApplicable(
                selected,
                expression,
                selectionProfile))
          return selected;
      }

    ERROR("the basis-conversion picker has no applicable complete "
          "plan from ",
          basisKeyForId(sourceBasisId).c_str(),
          " to ",
          basisKeyForId(targetBasisId).c_str());
    return {};
  }


} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
