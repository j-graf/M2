// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "basic-rings/aring-glue.hpp"
#include "error.h"
#include "symmetric-rings/basis-conversion-policy.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <utility>

namespace symmetric_rings {

// Contributor map:
//   1. shared exact expression facts and their stored cache;
//   2. condition-driven expression-piece facts;
//   3. plan indexing, structural contracts, and generic execution;
//   4. canonical source-group conversion;
//   5. the sole public built-in conversion workflow, toBasis.
//
// The matching declaration fragment follows the same order. Mathematical
// formulas live in basis-conversion-kernels.cpp, normalization rules in
// basis-normalization.cpp, complete conversion plans in
// basis-conversion-plans.cpp, and endpoint policy in
// basis-conversion-picker.cpp. Binary kernels, their commutative picker, and
// multiplication workflows live in the corresponding multiplication files.

namespace {

// Benchmark controls are scoped to one synchronous engine request. Ordinary
// toBasis installs an empty context, so ambient diagnostic variables cannot
// force plans, enable tracing, or alter nested multiplication.
struct BasisConversionDiagnosticContext
{
  std::optional<std::string> forcedPlanIdentifier;
  bool traceConversion = false;
};

thread_local const BasisConversionDiagnosticContext
    *activeBasisConversionDiagnostics = nullptr;

class ScopedBasisConversionDiagnostics
{
 public:
  ScopedBasisConversionDiagnostics(
      const BasisConversionDiagnosticContext& context,
      bool replaceExisting)
      : mPrevious(activeBasisConversionDiagnostics)
  {
    if (replaceExisting ||
        activeBasisConversionDiagnostics == nullptr)
      {
        activeBasisConversionDiagnostics = &context;
        mInstalled = true;
      }
  }

  ~ScopedBasisConversionDiagnostics()
  {
    if (mInstalled)
      activeBasisConversionDiagnostics = mPrevious;
  }

  ScopedBasisConversionDiagnostics(
      const ScopedBasisConversionDiagnostics&) = delete;
  ScopedBasisConversionDiagnostics& operator=(
      const ScopedBasisConversionDiagnostics&) = delete;

 private:
  const BasisConversionDiagnosticContext *mPrevious;
  bool mInstalled = false;
};

const std::optional<std::string>&
forcedBasisConversionPlanIdentifier()
{
  static const std::optional<std::string> noForcedPlan;
  return activeBasisConversionDiagnostics == nullptr
      ? noForcedPlan
      : activeBasisConversionDiagnostics->
            forcedPlanIdentifier;
}

} // namespace

bool SymmetricEngineRing::basisConversionContextActive() const
{
  return activeBasisConversionDiagnostics != nullptr;
}

bool SymmetricEngineRing::basisConversionTraceEnabled() const
{
  if (activeBasisConversionDiagnostics != nullptr)
    return activeBasisConversionDiagnostics->traceConversion;
  // Other development-only operations retain their existing outer-pipeline
  // trace. Any ordinary toBasis call they make installs the empty context
  // above and therefore cannot inherit this environment setting.
  return std::getenv(
      "M2_SYMMETRIC_RINGS_TRACE_CONVERSION") != nullptr;
}

// ============================================================================
// Conversion-Specific Expression Facts
// ============================================================================

SymmetricEngineRing::PowerSumSupportFacts
SymmetricEngineRing::inspectPowerSumSupportFacts(
    ring_elem expression,
    const std::vector<size_t>& termPositions,
    ExpressionFactRequirements requirements) const
{
    PowerSumSupportFacts result;
    if (requirements == noExpressionFactRequirements)
      return result;
    const auto& terms = polyValue(expression)->terms;
    bool allSingleCycles = true;
    size_t nonscalarTerms = 0;
    size_t mostlyShortCycleTerms = 0;
    std::vector<int> commonParts;
    bool firstIndex = true;
    bool allTermsArePowerSums = true;
    for (size_t position : termPositions)
      {
        if (position >= terms.size())
          {
            allTermsArePowerSums = false;
            break;
          }
        Partition index;
        if (!powerSumIndexFromMonomial(
                terms[position].monomial, index))
          {
            allTermsArePowerSums = false;
            break;
          }
        if (index.empty()) continue;
        ++nonscalarTerms;
        if ((requirements &
             requirePowerSumSingleCycleFacts) != 0 &&
            index.size() != 1)
          allSingleCycles = false;
        if ((requirements &
             requirePowerSumShortCycleFacts) != 0 &&
            powerSumIndexHasMostlyShortCycles(
                index, partitionWeight(index)))
          ++mostlyShortCycleTerms;
        if ((requirements &
             requirePowerSumCommonPartFacts) != 0)
          {
            std::vector<int> distinctParts(
                index.begin(), index.end());
            std::sort(
                distinctParts.begin(),
                distinctParts.end());
            distinctParts.erase(
                std::unique(
                    distinctParts.begin(),
                    distinctParts.end()),
                distinctParts.end());
            if (firstIndex)
              {
                commonParts =
                    std::move(distinctParts);
                firstIndex = false;
              }
            else
              {
                std::vector<int> intersection;
                std::set_intersection(
                    commonParts.begin(),
                    commonParts.end(),
                    distinctParts.begin(),
                    distinctParts.end(),
                    std::back_inserter(intersection));
                commonParts =
                    std::move(intersection);
              }
          }
      }
    if (!allTermsArePowerSums ||
        nonscalarTerms == 0)
      return result;
    if ((requirements &
         requirePowerSumSingleCycleFacts) != 0)
      result.allTermsAreSingleCycles =
          allSingleCycles;
    if ((requirements &
         requirePowerSumShortCycleFacts) != 0)
      result.mostlyShortCycleTermCount =
          mostlyShortCycleTerms;
    if ((requirements &
         requirePowerSumCommonPartFacts) != 0)
      result.commonParts = std::move(commonParts);
    return result;
  }

// ============================================================================
// Expression-Piece Inspection
// ============================================================================
// Conditions remain independent inspectable values. Each piece records the
// exact mathematical facts available from canonical symmetric-function
// storage.

std::vector<ExpressionPieceFacts>
SymmetricEngineRing::inspectExpressionPieces(
    ring_elem expression,
    ExpressionPieceKind pieceKind,
    const ExpressionFacts& facts,
    ExpressionFactRequirements requirements) const
{
    const auto& terms = polyValue(expression)->terms;
    auto commonPieceFacts = [&](ExpressionPieceKind kind) {
      ExpressionPieceFacts piece;
      piece.pieceKind = kind;
      piece.combinatorialTags = facts.combinatorialTags;
      piece.coefficientRingIsQQ = coefficientRing == globalQQ;
      return piece;
    };

    if (pieceKind == ExpressionPieceKind::WholeExpression)
      {
        ExpressionPieceFacts piece =
            commonPieceFacts(pieceKind);
        piece.termPositions.resize(terms.size());
        std::iota(
            piece.termPositions.begin(),
            piece.termPositions.end(),
            size_t{0});
        piece.weight = facts.homogeneousWeight;
        piece.termCount = facts.termCount;
        piece.singleBasisElement = facts.singleBasisElement();
        const PowerSumSupportFacts support =
            inspectPowerSumSupportFacts(
                expression,
                piece.termPositions,
                requirements);
        piece.allPowerSumTermsSingleCycles =
            support.allTermsAreSingleCycles;
        piece.mostlyShortCyclePowerSumTermCount =
            support.mostlyShortCycleTermCount;
        piece.commonPowerSumParts =
            support.commonParts;
        piece.index = facts.singleBasisElementIndex;
        std::optional<int> firstWeight;
        bool multipleWeights = false;
        for (const auto& term : terms)
          {
            const int weight = monomialWeight(term.monomial);
            if (!firstWeight)
              firstWeight = weight;
            else if (*firstWeight != weight)
              {
                multipleWeights = true;
                break;
              }
          }
        piece.expressionHasMultipleWeights = multipleWeights;
        return {std::move(piece)};
      }

    if (pieceKind == ExpressionPieceKind::IndividualTerms)
      {
        std::vector<ExpressionPieceFacts> result;
        result.reserve(terms.size());
        for (size_t position = 0; position < terms.size(); ++position)
          {
            const auto& term = terms[position];
            ExpressionPieceFacts piece =
                commonPieceFacts(pieceKind);
            piece.termPositions = {position};
            piece.weight = monomialWeight(term.monomial);
            piece.termCount = 1;
            if (!term.monomial.data.empty() &&
                atomLengthAt(term.monomial, 0) ==
                    term.monomial.data.size() &&
                !atomIsSkewAt(term.monomial, 0))
              {
                Partition index =
                    basisElementIndex(term.monomial, 0);
                piece.index = index;
                const PowerSumSupportFacts support =
                    inspectPowerSumSupportFacts(
                        expression,
                        piece.termPositions,
                        requirements);
                if (support.mostlyShortCycleTermCount)
                  piece.powerSumIndexHasMostlyShortCycles =
                      *support.mostlyShortCycleTermCount == 1;
              }
            result.push_back(std::move(piece));
          }
        return result;
      }

    std::map<int, std::vector<size_t>> positionsByWeight;
    for (size_t position = 0; position < terms.size(); ++position)
      positionsByWeight[monomialWeight(
          terms[position].monomial)].push_back(position);
    std::vector<ExpressionPieceFacts> result;
    result.reserve(positionsByWeight.size());
    for (auto& item : positionsByWeight)
      {
        ExpressionPieceFacts piece =
            commonPieceFacts(pieceKind);
        piece.termPositions = std::move(item.second);
        piece.weight = item.first;
        piece.termCount = piece.termPositions.size();
        piece.possibleTermCount =
            partitionCountForConversionSelection(item.first);
        const PowerSumSupportFacts support =
            inspectPowerSumSupportFacts(
                expression,
                piece.termPositions,
                requirements);
        piece.allPowerSumTermsSingleCycles =
            support.allTermsAreSingleCycles;
        piece.mostlyShortCyclePowerSumTermCount =
            support.mostlyShortCycleTermCount;
        piece.commonPowerSumParts =
            support.commonParts;
        result.push_back(std::move(piece));
      }
    return result;
  }

// ============================================================================
// Declarative Plan Contracts And Validation
// ============================================================================

const SymmetricEngineRing::BasisConversionPlanDefinition *
SymmetricEngineRing::basisConversionPlanDefinition(
    const BasisConversionPlanId& id)
{
    const auto& plansById = basisConversionPlansById();
    const auto found = plansById.find(id.value);
    return found == plansById.end() ? nullptr : found->second;
  }

const std::map<
    std::string,
    const SymmetricEngineRing::BasisConversionPlanDefinition *>&
SymmetricEngineRing::basisConversionPlansById()
{
    static const std::map<
        std::string,
        const BasisConversionPlanDefinition *> index = [] {
      std::map<
          std::string,
          const BasisConversionPlanDefinition *> result;
      for (const auto& plan : basisConversionPlanDatabase())
        if (plan.id.value.empty() ||
            !result.emplace(plan.id.value, &plan).second)
          throw exc::engine_error(
              "basis-conversion plan identifiers must be "
              "nonempty and unique");
      return result;
    }();
    return index;
  }

const std::map<
    std::pair<
        SymmetricEngineRing::BasisKind,
        SymmetricEngineRing::BasisKind>,
    std::vector<
        const SymmetricEngineRing::BasisConversionPlanDefinition *>>&
SymmetricEngineRing::basisConversionPlansBySourceAndTarget()
{
    static const std::map<
        std::pair<BasisKind, BasisKind>,
        std::vector<const BasisConversionPlanDefinition *>>
        index = [] {
          std::map<
              std::pair<BasisKind, BasisKind>,
              std::vector<const BasisConversionPlanDefinition *>>
              result;
          for (const auto& plan : basisConversionPlanDatabase())
            result[{plan.sourceBasisKind, plan.targetBasisKind}].
                push_back(&plan);
          return result;
        }();
    return index;
  }

void SymmetricEngineRing::validateBasisConversionPlanDatabase()
{
    static const bool validated = [] {
      const auto& database = basisConversionPlanDatabase();
      const auto& plansById = basisConversionPlansById();
      const auto& plansBySourceAndTarget =
          basisConversionPlansBySourceAndTarget();

      // The concrete plans together with the one generic u -> p -> v plan
      // cover the complete directed graph of built-in bases. Identity arrows
      // remain a generic workflow bypass.
      // The fallback descriptors are the registry of conversion-supported
      // built-in bases. Deriving the coverage universe from them prevents a
      // new basis from being listed once for fallback and again in validator
      // infrastructure.
      std::vector<BasisKind> builtInConversionBases{
          BasisKind::PowerSum};
      for (const auto& fallback :
           powerSumFallbackPlanDatabase())
        builtInConversionBases.push_back(
            fallback.basisKind);
      std::sort(
          builtInConversionBases.begin(),
          builtInConversionBases.end(),
          [](BasisKind left, BasisKind right) {
            return static_cast<int>(left) <
                   static_cast<int>(right);
          });
      builtInConversionBases.erase(
          std::unique(
              builtInConversionBases.begin(),
              builtInConversionBases.end()),
          builtInConversionBases.end());
      for (BasisKind source : builtInConversionBases)
        for (BasisKind target : builtInConversionBases)
          if (source != target)
            {
              const bool hasSpecificPlan =
                  plansBySourceAndTarget.find({source, target}) !=
                  plansBySourceAndTarget.end();
              const bool hasGenericPlan =
                  source != BasisKind::PowerSum &&
                  target != BasisKind::PowerSum;
              if (!hasSpecificPlan && !hasGenericPlan)
                throw exc::engine_error(
                    "basis-conversion plans do not cover every "
                    "ordered pair of built-in bases");
            }

      if (plansById.find(
              genericPowerSumFallbackPlanIdentifier().value) !=
          plansById.end())
        throw exc::engine_error(
            "the generic power-sum composition identifier must not "
            "belong to an endpoint-specific plan");
      std::map<BasisKind, const PowerSumFallbackPlanPair *>
          powerSumFallbacks;
      for (const auto& fallback :
           powerSumFallbackPlanDatabase())
        {
          if (fallback.basisKind == BasisKind::Custom ||
              fallback.basisKind == BasisKind::PowerSum ||
              !powerSumFallbacks.emplace(
                   fallback.basisKind, &fallback).second)
            throw exc::engine_error(
                "power-sum fallback bases must be distinct built-in "
                "non-power-sum bases");
          auto toPowerSums =
              plansById.find(fallback.toPowerSums.value);
          auto fromPowerSums =
              plansById.find(fallback.fromPowerSums.value);
          if (toPowerSums == plansById.end() ||
              fromPowerSums == plansById.end() ||
              toPowerSums->second->sourceBasisKind !=
                  fallback.basisKind ||
              toPowerSums->second->targetBasisKind !=
                  BasisKind::PowerSum ||
              fromPowerSums->second->sourceBasisKind !=
                  BasisKind::PowerSum ||
              fromPowerSums->second->targetBasisKind !=
                  fallback.basisKind ||
              toPowerSums->second->applicability.kind !=
                  ExpressionConditionKind::Always ||
              fromPowerSums->second->applicability.kind !=
                  ExpressionConditionKind::Always)
            throw exc::engine_error(
                "a power-sum fallback must name unconditional plans "
                "with the stated endpoints");
        }
      if (powerSumFallbacks.size() + 1 !=
          builtInConversionBases.size())
        throw exc::engine_error(
            "the generic power-sum composition does not cover every "
            "non-power-sum built-in basis");

      for (const auto& plan : database)
        {
          if (plan.sourceBasisKind == BasisKind::Custom ||
              plan.targetBasisKind == BasisKind::Custom ||
              plan.sourceBasisKind == plan.targetBasisKind)
            throw exc::engine_error(
                "a basis-conversion plan must have two distinct "
                "built-in endpoints");
          try
            {
              validateExpressionCondition(
                  plan.applicability,
                  ExpressionPieceKind::WholeExpression);
            }
          catch (const std::exception& exception)
            {
              throw exc::engine_error(
                  "invalid applicability condition for plan " +
                  plan.id.value + ": " + exception.what());
            }
          try
            {
              validateExpressionCondition(
                  plan.outputGuarantee,
                  ExpressionPieceKind::WholeExpression);
            }
          catch (const std::exception& exception)
            {
              throw exc::engine_error(
                  "invalid output guarantee for plan " +
                  plan.id.value + ": " + exception.what());
            }
          if (plan.cases.empty())
            throw exc::engine_error(
                "basis-conversion plan " + plan.id.value +
                " has no cases");
          for (size_t i = 0; i < plan.cases.size(); ++i)
            {
              const auto& item = plan.cases[i];
              const bool isOtherwise =
                  item.condition.kind ==
                  ExpressionConditionKind::Otherwise;
              if (isOtherwise !=
                  (i + 1 == plan.cases.size()))
                throw exc::engine_error(
                    "basis-conversion plan " + plan.id.value +
                    " must end in exactly one otherwise case");
              try
                {
                  validateExpressionCondition(
                      item.condition,
                      plan.pieceKind,
                      isOtherwise);
                }
              catch (const std::exception& exception)
                {
                  throw exc::engine_error(
                      "invalid case condition for plan " +
                      plan.id.value + ": " + exception.what());
                }
              if (const auto *kernelCase =
                      std::get_if<KernelPlan>(
                          &item.formula))
                {
                  if (kernelCase->identifier.empty() ||
                      kernelCase->kernel == nullptr ||
                      !expressionConditionImplies(
                          kernelCase->outputGuarantee,
                          plan.outputGuarantee))
                    throw exc::engine_error(
                        "basis-conversion plan " +
                        plan.id.value +
                        " has an invalid kernel formula");
                  try
                    {
                      validateExpressionCondition(
                          kernelCase->outputGuarantee,
                          ExpressionPieceKind::WholeExpression);
                    }
                  catch (const std::exception& exception)
                    {
                      throw exc::engine_error(
                          "invalid kernel output guarantee for plan " +
                          plan.id.value + ": " + exception.what());
                    }
                }
              else
                {
                  const auto& composition =
                      std::get<
                          CompositionPlan>(
                          item.formula);
                  if (composition.plans.empty())
                    throw exc::engine_error(
                        "basis-conversion plan " + plan.id.value +
                        " has an empty composition");
                  composition.planDefinitions.clear();
                }
            }
        }

      std::map<std::string, std::vector<std::string>>
          dependencies;
      for (const auto& plan : database)
        for (const auto& item : plan.cases)
          {
            if (std::holds_alternative<
                    KernelPlan>(
                    item.formula))
              continue;
            auto& composition =
                std::get<
                    CompositionPlan>(
                    item.formula);
            BasisKind current = plan.sourceBasisKind;
            ExpressionCondition guaranteedCondition =
                plan.pieceKind ==
                        ExpressionPieceKind::WholeExpression
                    ? plan.applicability
                    : always();
            if (plan.pieceKind ==
                    ExpressionPieceKind::WholeExpression &&
                item.condition.kind !=
                    ExpressionConditionKind::Otherwise)
              guaranteedCondition =
                  guaranteedCondition && item.condition;
            for (const auto& componentPlanId : composition.plans)
              {
                auto found =
                    plansById.find(componentPlanId.value);
                if (found == plansById.end())
                  throw exc::engine_error(
                      "basis-conversion plan " +
                      plan.id.value +
                      " references unknown component plan " +
                      componentPlanId.value);
                const auto& componentPlan = *found->second;
                composition.planDefinitions.push_back(
                    found->second);
                if (componentPlan.sourceBasisKind != current)
                  throw exc::engine_error(
                      "basis-conversion plan " +
                      plan.id.value +
                      " has incompatible composition endpoints");
                if (!expressionConditionImplies(
                        guaranteedCondition,
                        componentPlan.applicability))
                  throw exc::engine_error(
                      "basis-conversion composition " +
                      plan.id.value +
                      "references a plan with an unproved "
                      "applicability requirement");
                current = componentPlan.targetBasisKind;
                guaranteedCondition =
                    componentPlan.outputGuarantee;
                dependencies[plan.id.value].push_back(
                    componentPlan.id.value);
              }
            if (current != plan.targetBasisKind)
              throw exc::engine_error(
                  "basis-conversion plan " + plan.id.value +
                  " has the wrong composition target");
            if (!expressionConditionImplies(
                    guaranteedCondition,
                    plan.outputGuarantee))
              throw exc::engine_error(
                  "basis-conversion composition " +
                  plan.id.value +
                  " does not prove its declared output guarantee");
          }

      enum class VisitState { Unseen, Active, Complete };
      std::map<std::string, VisitState> state;
      std::function<void(const std::string&)> visit =
          [&](const std::string& id) {
            if (state[id] == VisitState::Complete) return;
            if (state[id] == VisitState::Active)
              throw exc::engine_error(
                  "basis-conversion plan dependencies are cyclic at " +
                  id);
            state[id] = VisitState::Active;
            for (const auto& dependency : dependencies[id])
              visit(dependency);
            state[id] = VisitState::Complete;
          };
      for (const auto& plan : database)
        visit(plan.id.value);
      return true;
    }();
    (void) validated;
  }

const std::vector<SymmetricEngineRing::ResolvedBasisConversionPlan>&
SymmetricEngineRing::basisConversionPlansFor(
    int sourceBasisId,
    int targetBasisId) const
{
    static const std::vector<ResolvedBasisConversionPlan> empty;
    requireBasis(sourceBasisId);
    requireBasis(targetBasisId);
    if (error()) return empty;
    const std::pair<int, int> endpoint{
        sourceBasisId, targetBasisId};
    auto cached =
        ringBasisConversionPlanCache.find(endpoint);
    if (cached != ringBasisConversionPlanCache.end())
      return cached->second;
    validateBasisConversionPlanDatabase();
    const BasisKind source = basisKindForId(sourceBasisId);
    const BasisKind target = basisKindForId(targetBasisId);
    std::vector<ResolvedBasisConversionPlan> result;
    const auto& plansBySourceAndTarget =
        basisConversionPlansBySourceAndTarget();
    auto found = plansBySourceAndTarget.find({source, target});
    if (found != plansBySourceAndTarget.end())
      {
        result.reserve(found->second.size());
        for (const auto *plan : found->second)
          result.push_back(
              {plan, sourceBasisId, targetBasisId});
      }
    auto inserted =
        ringBasisConversionPlanCache.emplace(
            endpoint, std::move(result));
    return inserted.first->second;
  }

SymmetricEngineRing::ResolvedBasisConversionPlan
SymmetricEngineRing::basisConversionPlanForRing(
    const BasisConversionPlanId& id) const
{
    validateBasisConversionPlanDatabase();
    const auto *definition =
        basisConversionPlanDefinition(id);
    if (definition == nullptr)
      {
        ERROR("unknown basis-conversion plan: ",
              id.value.c_str());
        return {};
      }
    return basisConversionPlanForRing(definition);
  }

SymmetricEngineRing::ResolvedBasisConversionPlan
SymmetricEngineRing::basisConversionPlanForRing(
    const BasisConversionPlanDefinition *definition) const
{
    if (definition == nullptr) return {};
    auto cached =
        basisConversionEndpointsForRingCache.find(
            definition->id.value);
    if (cached !=
        basisConversionEndpointsForRingCache.end())
      return {
          definition,
          cached->second.first,
          cached->second.second};
    const int sourceBasisId =
        requiredBasisIdForKind(
            definition->sourceBasisKind);
    const int targetBasisId =
        requiredBasisIdForKind(
            definition->targetBasisKind);
    if (error()) return {};
    basisConversionEndpointsForRingCache.emplace(
        definition->id.value,
        std::make_pair(sourceBasisId, targetBasisId));
    return {definition, sourceBasisId, targetBasisId};
  }

// ============================================================================
// Generic Declarative Formula And Plan Executor
// ============================================================================
// The executor interprets the already-selected plan. Kernel formulas call one
// kernel; compositions identify their fixed component plans recursively.
// Neither path invokes performance policy or substitutes another plan.

ring_elem SymmetricEngineRing::executeBasisConversionPlanFormula(
    const BasisConversionPlanFormula& formula,
    int sourceBasisId,
    int targetBasisId,
    ring_elem expression,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts& inputFacts,
    ExpressionFacts *resultFacts) const
{
    if (!inputFacts.canonicalExpansionInBasis(sourceBasisId))
      {
        ERROR("a basis-conversion formula received input outside "
              "its canonical source-basis contract");
        return zero();
      }
    if (const auto *kernelCase =
            std::get_if<KernelPlan>(
                &formula))
      {
        if (kernelCase->kernel == nullptr)
          {
            ERROR("a basis-conversion kernel formula has no callable");
            return zero();
          }
        const auto& sourceDescriptor =
            requireBasis(sourceBasisId);
        const auto& targetDescriptor =
            requireBasis(targetBasisId);
        if (error()) return zero();
        BasisConversionInput input{
            expression,
            {sourceDescriptor.kind,
             sourceBasisId,
             &sourceDescriptor.displaySymbol,
             sourceDescriptor.displayOrder},
            {targetDescriptor.kind,
             targetBasisId,
             &targetDescriptor.displaySymbol,
             targetDescriptor.displayOrder},
            inputFacts.homogeneousWeight};
        ring_elem result = (this->*kernelCase->kernel)(input);
        if (error()) return zero();

        const bool needFacts =
            resultFacts != nullptr ||
            kernelCase->outputGuarantee.kind !=
                ExpressionConditionKind::Always;
        ExpressionFacts kernelFacts;
        if (needFacts)
          {
            kernelFacts = inferCanonicalExpansionFacts(
                result,
                targetBasisId,
                inputFacts.homogeneousWeight);
            kernelFacts.combinatorialTags =
                combinatorialTags;
            if (!kernelFacts.canonicalExpansionInBasis(
                    targetBasisId))
              {
                ERROR("basis-conversion kernel ",
                      kernelCase->identifier.c_str(),
                      " violated its canonical target contract");
                return zero();
              }
          }
        if (kernelCase->outputGuarantee.kind !=
            ExpressionConditionKind::Always)
          {
            const auto pieceFacts =
                inspectExpressionPieces(
                    result,
                    ExpressionPieceKind::WholeExpression,
                    kernelFacts,
                    expressionFactRequirements(
                        kernelCase->outputGuarantee));
            if (pieceFacts.size() != 1 ||
                !expressionConditionHolds(
                    kernelCase->outputGuarantee,
                    pieceFacts.front()))
              {
                ERROR("basis-conversion kernel ",
                      kernelCase->identifier.c_str(),
                      " violated its declared output guarantee");
                return zero();
              }
          }
        mutablePolyValue(result)->combinatorialTags =
            combinatorialTags;
        if (resultFacts != nullptr)
          *resultFacts = std::move(kernelFacts);
        return result;
      }

    ring_elem current = expression;
    ExpressionFacts currentFacts = inputFacts;
    int currentBasisId = sourceBasisId;
    auto executeComponentPlan =
        [&](const BasisConversionPlanDefinition *componentDefinition)
            -> bool {
          if (componentDefinition == nullptr)
            {
              ERROR("a basis-conversion composition reached execution "
                    "before a component plan was identified");
              return false;
            }
          ResolvedBasisConversionPlan componentPlan =
              basisConversionPlanForRing(
                  componentDefinition);
          if (error()) return false;
          if (componentPlan.sourceBasisId != currentBasisId)
            {
              ERROR("a basis-conversion composition has "
                    "incompatible runtime endpoints");
              return false;
            }
          current = executeBasisConversionPlan(
              componentPlan,
              current,
              combinatorialTags,
              currentFacts,
              nullptr);
          if (error()) return false;
          currentBasisId = componentPlan.targetBasisId;
          // Conditions in the next fixed component plan apply to the actual
          // realized intermediate, never to a predicted support profile.
          currentFacts = inferCanonicalExpansionFacts(
              current,
              currentBasisId,
              inputFacts.homogeneousWeight);
          currentFacts.combinatorialTags =
              combinatorialTags;
          return true;
        };

    const auto& composition =
        std::get<
            CompositionPlan>(
            formula);
    if (composition.planDefinitions.size() !=
        composition.plans.size())
      {
        ERROR("a basis-conversion composition reached execution "
              "before its named plans were identified");
        return zero();
      }
    for (const auto *componentDefinition :
         composition.planDefinitions)
      if (!executeComponentPlan(componentDefinition))
        return zero();
    if (currentBasisId != targetBasisId ||
        !currentFacts.canonicalExpansionInBasis(targetBasisId))
      {
        ERROR("a basis-conversion composition violated its target "
              "contract");
        return zero();
      }
    if (resultFacts != nullptr)
      {
        *resultFacts = std::move(currentFacts);
        resultFacts->combinatorialTags =
            combinatorialTags;
      }
    return current;
  }

ring_elem SymmetricEngineRing::executeBasisConversionPlan(
    const ResolvedBasisConversionPlan& plan,
    ring_elem expression,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts& inputFacts,
    ExpressionFacts *resultFacts) const
{
    if (!plan.valid() ||
        !inputFacts.canonicalExpansionInBasis(
            plan.sourceBasisId))
      {
        ERROR("an invalid complete basis-conversion plan reached "
              "the executor");
        return zero();
      }
    struct ExecutionScope
    {
      size_t& depth;
      explicit ExecutionScope(size_t& value) : depth(value)
      {
        ++depth;
      }
      ~ExecutionScope() { --depth; }
    } executionScope(basisConversionPlanExecutionDepth);

    const bool wholeExpressionOtherwise =
        plan.definition->pieceKind ==
            ExpressionPieceKind::WholeExpression &&
        plan.definition->cases.size() == 1 &&
        plan.definition->cases.front().condition.kind ==
            ExpressionConditionKind::Otherwise;
    ExpressionFacts executionFacts = inputFacts;
    executionFacts.combinatorialTags = combinatorialTags;
    if (!basisConversionPlanApplicable(
            plan, expression, executionFacts))
      {
        ERROR("a selected basis-conversion plan is inapplicable "
              "to its runtime input: ",
              plan.definition->id.value.c_str());
        return zero();
      }
    if (basisConversionTraceEnabled())
      std::fprintf(
          stderr,
          "SymmetricRings conversion-plan: source=%s target=%s "
          "plan=%s\n",
          basisKeyForId(plan.sourceBasisId).c_str(),
          basisKeyForId(plan.targetBasisId).c_str(),
          plan.definition->id.value.c_str());

    if (wholeExpressionOtherwise)
      {
        // Most available plans are one whole-expression formula. Their
        // declarative meaning is exactly a direct call on the original input;
        // constructing a one-piece partition and then copying and recollecting
        // the complete expression would add work without changing semantics.
        const bool needDirectFacts =
            resultFacts != nullptr ||
            plan.definition->outputGuarantee.kind !=
                ExpressionConditionKind::Always;
        ExpressionFacts directFacts;
        ring_elem directResult =
            executeBasisConversionPlanFormula(
                plan.definition->cases.front().formula,
                plan.sourceBasisId,
                plan.targetBasisId,
                expression,
                combinatorialTags,
                executionFacts,
                needDirectFacts ? &directFacts : nullptr);
        if (error()) return zero();
        if (needDirectFacts &&
            !directFacts.canonicalExpansionInBasis(
                plan.targetBasisId))
          {
            ERROR("a complete basis-conversion plan violated its "
                  "canonical output contract");
            return zero();
          }
        if (plan.definition->outputGuarantee.kind !=
            ExpressionConditionKind::Always)
          {
            try
              {
                const auto outputPieceFacts =
                    inspectExpressionPieces(
                        directResult,
                        ExpressionPieceKind::WholeExpression,
                        directFacts,
                        expressionFactRequirements(
                            plan.definition->
                                outputGuarantee));
                if (outputPieceFacts.size() != 1 ||
                    !expressionConditionHolds(
                        plan.definition->outputGuarantee,
                        outputPieceFacts.front()))
                  {
                    ERROR("a complete basis-conversion plan violated its "
                          "declared output guarantee");
                    return zero();
                  }
              }
            catch (const std::exception& exception)
              {
                ERROR("basis-conversion output-guarantee validation "
                      "failed: ",
                      exception.what());
                return zero();
              }
          }
        mutablePolyValue(directResult)->combinatorialTags =
            combinatorialTags;
        if (resultFacts != nullptr)
          *resultFacts = std::move(directFacts);
        return directResult;
      }

    std::vector<ExpressionPieceFacts> pieces;
    PlanCaseAssignments partition;
    try
      {
        std::vector<ExpressionCondition> conditions;
        conditions.reserve(
            plan.definition->cases.size());
        for (const auto& item : plan.definition->cases)
          conditions.push_back(item.condition);
        pieces = inspectExpressionPieces(
            expression,
            plan.definition->pieceKind,
            executionFacts,
            expressionFactRequirements(conditions));
        partition = assignExpressionPiecesToCases(
            pieces,
            conditions,
            polyValue(expression)->terms.size());
      }
    catch (const std::exception& exception)
      {
        ERROR("basis-conversion plan partitioning failed: ",
              exception.what());
        return zero();
      }

    const auto& sourceTerms =
        polyValue(expression)->terms;
    VECTOR(SymmetricTerm) resultTerms;
    for (size_t caseIndex = 0;
         caseIndex < plan.definition->cases.size();
         ++caseIndex)
      {
        const auto& termPositions =
            partition.termPositionsForCase[caseIndex];
        if (termPositions.empty()) continue;

        // A case formula receives the complete subexpression assigned to that
        // case, not one invocation per term or component. This preserves the
        // declarative semantics and lets fixed component plans inspect the
        // actual aggregate support produced for the case.
        VECTOR(SymmetricTerm) caseTerms;
        caseTerms.reserve(termPositions.size());
        for (size_t position : termPositions)
          {
            if (position >= sourceTerms.size())
              {
                ERROR("a basis-conversion piece contains an invalid "
                      "term position");
                return zero();
              }
            caseTerms.push_back(sourceTerms[position]);
          }
        ring_elem caseExpression =
            fromTermVector(caseTerms, true);
        ExpressionFacts caseFacts;
        if (plan.definition->pieceKind ==
            ExpressionPieceKind::WholeExpression)
          caseFacts = executionFacts;
        else
          {
            caseFacts = inferCanonicalExpansionFacts(
                caseExpression,
                plan.sourceBasisId);
            caseFacts.combinatorialTags =
                combinatorialTags;
          }
        ExpressionFacts convertedFacts;
        ring_elem converted =
            executeBasisConversionPlanFormula(
                plan.definition->cases[caseIndex].formula,
                plan.sourceBasisId,
                plan.targetBasisId,
                caseExpression,
                combinatorialTags,
                caseFacts,
                &convertedFacts);
        if (error()) return zero();
        if (!convertedFacts.canonicalExpansionInBasis(
                plan.targetBasisId))
          {
            ERROR("a basis-conversion case violated its canonical "
                  "target contract");
            return zero();
          }
        for (const auto& term :
             polyValue(converted)->terms)
          appendTermIfNonZero(
              resultTerms, term.coeff, term.monomial);
      }

    ring_elem result =
        fromTermVector(resultTerms, false);
    ExpressionFacts finalFacts =
        inferCanonicalExpansionFacts(
            result,
            plan.targetBasisId,
            executionFacts.homogeneousWeight);
    finalFacts.combinatorialTags = combinatorialTags;
    if (!finalFacts.canonicalExpansionInBasis(
            plan.targetBasisId))
      {
        ERROR("a complete basis-conversion plan violated its "
              "canonical output contract");
        return zero();
      }
    if (plan.definition->outputGuarantee.kind !=
        ExpressionConditionKind::Always)
      {
        try
          {
            const auto outputPieceFacts =
                inspectExpressionPieces(
                    result,
                    ExpressionPieceKind::WholeExpression,
                    finalFacts,
                    expressionFactRequirements(
                        plan.definition->
                            outputGuarantee));
            if (outputPieceFacts.size() != 1 ||
                !expressionConditionHolds(
                    plan.definition->outputGuarantee,
                    outputPieceFacts.front()))
              {
                ERROR("a complete basis-conversion plan violated its "
                      "declared output guarantee");
                return zero();
              }
          }
        catch (const std::exception& exception)
          {
            ERROR("basis-conversion output-guarantee validation "
                  "failed: ",
                  exception.what());
            return zero();
          }
      }
    mutablePolyValue(result)->combinatorialTags =
        combinatorialTags;
    if (resultFacts != nullptr)
      *resultFacts = std::move(finalFacts);
    return result;
  }

ring_elem SymmetricEngineRing::executeNamedBasisConversionPlan(
    const BasisConversionPlanId& planIdentifier,
    ring_elem expression,
    ExpressionFacts *resultFacts) const
{
    const ResolvedBasisConversionPlan plan =
        basisConversionPlanForRing(planIdentifier);
    if (error() || !plan.valid()) return zero();
    ExpressionFacts inputFacts =
        inferCanonicalExpansionFacts(
            expression, plan.sourceBasisId);
    inputFacts.combinatorialTags =
        polyValue(expression)->combinatorialTags;
    if (!inputFacts.canonicalExpansionInBasis(
            plan.sourceBasisId))
      {
        ERROR("named basis-conversion plan ",
              planIdentifier.value.c_str(),
              " received a noncanonical source expansion");
        return zero();
      }
    return executeBasisConversionPlan(
        plan,
        expression,
        inputFacts.combinatorialTags,
        inputFacts,
        resultFacts);
  }

// ============================================================================
// Canonical Group Conversion Helper
// ============================================================================
// This helper owns source-basis grouping, selected-plan execution, and final
// collection for input that is already normalized, skew-free, and product-free.
// Basis-specific partitioning belongs to the selected declarative plan.

ring_elem SymmetricEngineRing::convertCanonicalExpressionToBasis(
    ring_elem f,
    int targetBasisId,
    CombinatorialTags combinatorialTags,
    const ExpressionFacts *knownFacts) const
{
    if (binaryMultiplicationKernelExecutionDepth != 0)
      {
        ERROR("a binary multiplication kernel attempted to invoke "
              "basis conversion");
        return zero();
      }
    const bool trace = basisConversionTraceEnabled();

    // Exact zero/scalar facts satisfy the target contract without basis
    // grouping or plan selection.
    if (knownFacts != nullptr &&
        knownFacts->singleFactorTermCount == 0)
      {
        ring_elem result = copyPolyValue(polyValue(f));
        attachCanonicalExpansionFacts(
            result,
            *knownFacts,
            targetBasisId,
            combinatorialTags);
        return result;
      }

    // A pure canonical expansion can enter complete-plan selection directly.
    // The selected declarative plan owns its basis-specific partitioning.
    if (knownFacts != nullptr && knownFacts->expandedBasis)
      {
        int sourceBasisId = *knownFacts->expandedBasis;
        ExpressionFacts resultFacts;
        ring_elem result;
        if (sourceBasisId == targetBasisId)
          {
            result = copyPolyValue(polyValue(f));
            resultFacts = *knownFacts;
          }
        else
          {
            ExpressionFacts selectedInputFacts = *knownFacts;
            ResolvedBasisConversionPlan selection =
                selectBasisConversionPlan(
                    f,
                    sourceBasisId,
                    targetBasisId,
                    *knownFacts,
                    combinatorialTags,
                    forcedBasisConversionPlanIdentifier());
            if (error()) return zero();
            result = executeBasisConversionPlan(
                selection,
                f,
                combinatorialTags,
                selectedInputFacts,
                &resultFacts);
            if (error()) return zero();
          }
        if (!resultFacts.canonicalExpansionInBasis(targetBasisId))
          {
            ERROR("the canonical fact-cache bypass did not produce "
                  "a target-basis expansion");
            return zero();
          }
        attachCanonicalExpansionFacts(
            result,
            resultFacts,
            targetBasisId,
            combinatorialTags);
        return result;
      }

    struct PreparedConversion
    {
      ring_elem expression;
      ExpressionFacts facts;
      std::optional<ResolvedBasisConversionPlan> selection;
    };

    std::vector<std::pair<int, VECTOR(SymmetricTerm)>> groups;
    VECTOR(SymmetricTerm) scalarTerms;
    for (const auto& term : polyValue(f)->terms)
      {
        if (term.monomial.data.empty())
          {
            scalarTerms.push_back(term);
            continue;
          }
        if (atomLengthAt(term.monomial, 0) != term.monomial.data.size())
          {
            ERROR("the conversion workflow expected product-free input");
            return zero();
          }
        const int basisId = atomBasisIdAt(term.monomial, 0);
        auto group = std::find_if(
            groups.begin(),
            groups.end(),
            [&](const auto& candidate) {
              return candidate.first == basisId;
            });
        if (group == groups.end())
          {
            groups.push_back({basisId, {}});
            group = std::prev(groups.end());
          }
        group->second.push_back(term);
      }

    std::vector<PreparedConversion> preparedConversions;
    auto prepareGroup = [&](int sourceBasisId,
                            ring_elem expression,
                            ExpressionFacts facts) {
      std::optional<ResolvedBasisConversionPlan> selection;
      if (sourceBasisId != targetBasisId)
        {
          selection = selectBasisConversionPlan(
              expression,
              sourceBasisId,
              targetBasisId,
              facts,
              combinatorialTags,
              forcedBasisConversionPlanIdentifier());
          if (error()) return;
        }
      preparedConversions.push_back({
          expression, std::move(facts), std::move(selection)});
    };

    if (trace)
      std::fprintf(
          stderr,
          "SymmetricRings conversion-stage: stage=group-by-basis "
          "groups=%zu scalar-terms=%zu\n",
          groups.size(),
          scalarTerms.size());
    for (auto& group : groups)
      {
        ring_elem groupExpression = fromTermVector(group.second, true);
        bool exactKnownGroup =
            knownFacts != nullptr &&
            knownFacts->expandedBasis == group.first &&
            knownFacts->scalarTermCount == 0 &&
            groups.size() == 1;
        ExpressionFacts groupFacts = exactKnownGroup
            ? *knownFacts
            : inferExpressionFacts(groupExpression);
        prepareGroup(
            group.first, groupExpression, std::move(groupFacts));
        if (error()) return zero();
      }

    // Selection for every source group is complete before any group executes.
    // Each selected definition fixes every case, kernel, and component plan.
    VECTOR(SymmetricTerm) resultTerms = std::move(scalarTerms);
    for (const auto& prepared : preparedConversions)
      {
        ring_elem converted = prepared.selection
            ? executeBasisConversionPlan(
                  *prepared.selection,
                  prepared.expression,
                  combinatorialTags,
                  prepared.facts)
            : prepared.expression;
        if (error()) return zero();
        for (const auto& convertedTerm : polyValue(converted)->terms)
          appendTermIfNonZero(
              resultTerms,
              convertedTerm.coeff,
              convertedTerm.monomial);
      }

    // Source groups can overlap in the target. Merge them once, after every
    // independently selected complete plan has executed.
    ring_elem result = fromTermVector(resultTerms, false);
    ExpressionFacts resultFacts =
        inferCanonicalExpansionFacts(result, targetBasisId);
    if (!resultFacts.canonicalExpansionInBasis(targetBasisId))
      {
        ERROR("the canonical group-conversion helper did not produce "
              "a target-basis expansion");
        return zero();
      }
    mutablePolyValue(result)->combinatorialTags = combinatorialTags;
    attachCanonicalExpansionFacts(
        result,
        resultFacts,
        targetBasisId,
        combinatorialTags);
    return result;
  }

// ============================================================================
// Unified Basis-Conversion Entry Workflow
// ============================================================================
// This is the sole owning engine workflow for built-in conversion. Its body
// follows the mathematical decomposition: prepare canonical factors, resolve
// products, decompose by source basis, apply one complete plan to each source
// expansion, and collect the target expansion.

ring_elem SymmetricEngineRing::toBasis(
    ring_elem f,
    int targetBasisId) const
{
    const BasisConversionDiagnosticContext ordinaryConversion;
    const ScopedBasisConversionDiagnostics diagnosticScope(
        ordinaryConversion,
        false);

    ExpressionNormalizationOptions options =
        pipelinePreparationNormalizationOptions();
    options.productTargetBasisId = targetBasisId;
    ExpressionNormalizationResult normalized =
        normalizeExpressionWithOptionsDetailed(f, options);
    if (error()) return zero();
    const CombinatorialTags combinatorialTags =
        normalized.facts.combinatorialTags;
    if (basisConversionTraceEnabled())
      {
        if (normalized.usedCompleteCanonicalFactsBypass)
          std::fprintf(
              stderr,
              "SymmetricRings toBasis: "
              "bypass=canonical-facts-cache target=%s\n",
              displayForBasis(targetBasisId).c_str());
        else
          std::fprintf(
              stderr,
              "SymmetricRings toBasis: target=%s terms=%zu bases=%zu "
              "straightened=%s skew-expanded=%s "
              "product-terms-resolved=%zu\n",
              displayForBasis(targetBasisId).c_str(),
              normalized.facts.termCount,
              normalized.facts.factorBases.size(),
              normalized.performedStraightening ? "yes" : "no",
              normalized.performedSkewExpansion ? "yes" : "no",
              normalized.resolvedProductTermCount);
      }

    if (normalized.resolvedProducts &&
        normalized.facts.canonicalExpansionInBasis(
            targetBasisId))
      {
        attachCanonicalExpansionFacts(
            normalized.expression,
            normalized.facts,
            targetBasisId,
            combinatorialTags);
        return normalized.expression;
      }

    // Linearity now applies: the canonical helper decomposes the expression
    // by source basis, selects one complete source-to-target plan for each
    // summand, executes those plans, and collects their target expansions.
    return convertCanonicalExpressionToBasis(
        normalized.expression,
        targetBasisId,
        combinatorialTags,
        &normalized.facts);
  }

ring_elem SymmetricEngineRing::toBasisBench(
    ring_elem f,
    int targetBasisId,
    const std::optional<std::string>& forcedPlanIdentifier,
    bool traceConversion) const
{
    const BasisConversionDiagnosticContext benchmarkConversion{
        forcedPlanIdentifier,
        traceConversion};
    const ScopedBasisConversionDiagnostics diagnosticScope(
        benchmarkConversion,
        true);
    return toBasis(f, targetBasisId);
  }

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
