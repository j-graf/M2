// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"

#include <algorithm>
#include <cstdio>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace symmetric_rings {

// This file is the single inventory and performance policy for strict binary
// multiplication. Endpoints are commutative unordered pairs; no expression or
// multifactor strategy is selected here.

// ============================================================================
// Commutative Endpoints And Strict-Pair Facts
// ============================================================================

SymmetricEngineRing::UnorderedBasisPair
SymmetricEngineRing::UnorderedBasisPair::from(
    BasisKind left, BasisKind right)
{
    if (static_cast<int>(right) < static_cast<int>(left))
      std::swap(left, right);
    return {left, right};
  }

int SymmetricEngineRing::BinaryMultiplicationFacts::combinedWeight() const
{
    return partitionWeight(left->index) + partitionWeight(right->index);
  }

const Partition *
SymmetricEngineRing::BinaryMultiplicationFacts::uniqueIndexInBasis(
    BasisKind basis) const
{
    const bool leftMatches = left->basis.kind == basis;
    const bool rightMatches = right->basis.kind == basis;
    if (leftMatches == rightMatches) return nullptr;
    return leftMatches ? &left->index : &right->index;
  }

SymmetricEngineRing::BinaryMultiplicationFacts
SymmetricEngineRing::binaryMultiplicationFacts(
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right,
    BasisKind targetBasis) const
{
    return {
        &left,
        &right,
        UnorderedBasisPair::from(
            left.basis.kind, right.basis.kind),
        targetBasis};
  }

// ============================================================================
// Inspectable Condition Language
// ============================================================================

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::binaryAlways()
{
    return {};
  }

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::binaryOtherwise()
{
    BinaryMultiplicationCondition result;
    result.kind = BinaryMultiplicationConditionKind::Otherwise;
    return result;
  }

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::completeIndexIsOneRow()
{
    BinaryMultiplicationCondition result;
    result.kind =
        BinaryMultiplicationConditionKind::
            BasisIndexHasAtMostOnePart;
    result.namedBasis = BasisKind::Complete;
    return result;
  }

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::elementaryIndexIsOneColumn()
{
    BinaryMultiplicationCondition result;
    result.kind =
        BinaryMultiplicationConditionKind::
            BasisIndexHasAtMostOnePart;
    result.namedBasis = BasisKind::Elementary;
    return result;
  }

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::powerSumIndexIsSingleCycle()
{
    BinaryMultiplicationCondition result;
    result.kind =
        BinaryMultiplicationConditionKind::
            BasisIndexHasAtMostOnePart;
    result.namedBasis = BasisKind::PowerSum;
    return result;
  }

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::combinedWeightAtMost(int weight)
{
    BinaryMultiplicationCondition result;
    result.kind =
        BinaryMultiplicationConditionKind::CombinedWeightAtMost;
    result.integer = weight;
    return result;
  }

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::binaryAnd(
    BinaryMultiplicationCondition left,
    BinaryMultiplicationCondition right)
{
    BinaryMultiplicationCondition result;
    result.kind = BinaryMultiplicationConditionKind::And;
    result.operands = {std::move(left), std::move(right)};
    return result;
  }

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::binaryOr(
    BinaryMultiplicationCondition left,
    BinaryMultiplicationCondition right)
{
    BinaryMultiplicationCondition result;
    result.kind = BinaryMultiplicationConditionKind::Or;
    result.operands = {std::move(left), std::move(right)};
    return result;
  }

SymmetricEngineRing::BinaryMultiplicationCondition
SymmetricEngineRing::binaryNot(
    BinaryMultiplicationCondition condition)
{
    BinaryMultiplicationCondition result;
    result.kind = BinaryMultiplicationConditionKind::Not;
    result.operands = {std::move(condition)};
    return result;
  }

bool SymmetricEngineRing::binaryMultiplicationConditionHolds(
    const BinaryMultiplicationCondition& condition,
    const BinaryMultiplicationFacts& facts)
{
    switch (condition.kind)
      {
        case BinaryMultiplicationConditionKind::Always:
          return true;
        case BinaryMultiplicationConditionKind::Otherwise:
          throw std::logic_error(
              "otherwise is an ordered binary-picker marker, "
              "not an evaluable condition");
        case BinaryMultiplicationConditionKind::
            BasisIndexHasAtMostOnePart:
          {
            const Partition *index =
                facts.uniqueIndexInBasis(condition.namedBasis);
            if (index == nullptr)
              throw std::logic_error(
                  "a basis-named binary condition does not "
                  "identify exactly one factor");
            return index->size() <= 1;
          }
        case BinaryMultiplicationConditionKind::
            CombinedWeightAtMost:
          return facts.combinedWeight() <= condition.integer;
        case BinaryMultiplicationConditionKind::And:
          return std::all_of(
              condition.operands.begin(),
              condition.operands.end(),
              [&](const BinaryMultiplicationCondition& operand) {
                return binaryMultiplicationConditionHolds(
                    operand, facts);
              });
        case BinaryMultiplicationConditionKind::Or:
          return std::any_of(
              condition.operands.begin(),
              condition.operands.end(),
              [&](const BinaryMultiplicationCondition& operand) {
                return binaryMultiplicationConditionHolds(
                    operand, facts);
              });
        case BinaryMultiplicationConditionKind::Not:
          return !binaryMultiplicationConditionHolds(
              condition.operands.front(), facts);
      }
    return false;
  }

std::string
SymmetricEngineRing::binaryMultiplicationConditionToString(
    const BinaryMultiplicationCondition& condition)
{
    switch (condition.kind)
      {
        case BinaryMultiplicationConditionKind::Always:
          return "always";
        case BinaryMultiplicationConditionKind::Otherwise:
          return "otherwise";
        case BinaryMultiplicationConditionKind::
            BasisIndexHasAtMostOnePart:
          if (condition.namedBasis == BasisKind::Complete)
            return "the complete index is one row";
          if (condition.namedBasis == BasisKind::Elementary)
            return "the elementary index is one column";
          if (condition.namedBasis == BasisKind::PowerSum)
            return "the power-sum index is one cycle";
          return "the named-basis index has at most one part";
        case BinaryMultiplicationConditionKind::
            CombinedWeightAtMost:
          return "combined weight at most " +
                 std::to_string(condition.integer);
        case BinaryMultiplicationConditionKind::And:
        case BinaryMultiplicationConditionKind::Or:
          {
            const char *operation =
                condition.kind ==
                    BinaryMultiplicationConditionKind::And
                    ? " and "
                    : " or ";
            std::ostringstream out;
            out << "(";
            for (size_t i = 0; i < condition.operands.size(); ++i)
              {
                if (i != 0) out << operation;
                out << binaryMultiplicationConditionToString(
                    condition.operands[i]);
              }
            out << ")";
            return out.str();
          }
        case BinaryMultiplicationConditionKind::Not:
          return "not (" +
                 binaryMultiplicationConditionToString(
                     condition.operands.front()) +
                 ")";
      }
    return "unknown";
  }

void SymmetricEngineRing::validateBinaryMultiplicationCondition(
    const BinaryMultiplicationCondition& condition,
    const UnorderedBasisPair& endpoint,
    bool allowOtherwise)
{
    if (condition.kind ==
            BinaryMultiplicationConditionKind::Otherwise &&
        !allowOtherwise)
      throw std::logic_error(
          "otherwise is valid only as the final binary-picker "
          "preference");

    if (condition.kind ==
        BinaryMultiplicationConditionKind::
            BasisIndexHasAtMostOnePart)
      {
        const int occurrences =
            static_cast<int>(
                endpoint.first == condition.namedBasis) +
            static_cast<int>(
                endpoint.second == condition.namedBasis);
        if (occurrences != 1)
          throw std::logic_error(
              "a basis-named binary condition must identify "
              "exactly one endpoint factor");
      }

    const bool logical =
        condition.kind == BinaryMultiplicationConditionKind::And ||
        condition.kind == BinaryMultiplicationConditionKind::Or ||
        condition.kind == BinaryMultiplicationConditionKind::Not;
    if (logical)
      {
        const size_t expected =
            condition.kind == BinaryMultiplicationConditionKind::Not
                ? 1
                : 2;
        if (condition.operands.size() != expected)
          throw std::logic_error(
              "a binary logical condition has invalid arity");
      }
    else if (!condition.operands.empty())
      throw std::logic_error(
          "a primitive binary condition cannot have operands");

    for (const auto& operand : condition.operands)
      validateBinaryMultiplicationCondition(
          operand, endpoint, false);
  }

// ============================================================================
// Mathematical Kernel Inventory
// ============================================================================

const std::vector<
    SymmetricEngineRing::BinaryMultiplicationKernelDefinition>&
SymmetricEngineRing::binaryMultiplicationKernelDatabase()
{
    static const std::vector<
        BinaryMultiplicationKernelDefinition> definitions = [] {
      using K = BasisKind;
      std::vector<BinaryMultiplicationKernelDefinition> result;
      auto addMultiplicationKernel =
          [&](std::string identifier,
              K first,
              K second,
              K target,
              BinaryMultiplicationCondition applicableWhen,
              BinaryMultiplicationKernel kernel) {
            result.push_back({
                std::move(identifier),
                first,
                second,
                target,
                std::move(applicableWhen),
                kernel});
          };

      // Schur products.
      addMultiplicationKernel("Schur*Schur->Schur:"
          "Littlewood-Richardson-tableaux",
          K::Schur,
          K::Schur,
          K::Schur,
          binaryAlways(),
          &SymmetricEngineRing::
              schurTimesSchurViaTableauLittlewoodRichardson);
      addMultiplicationKernel("Schur*Schur->Schur:"
          "Littlewood-Richardson-coefficients",
          K::Schur,
          K::Schur,
          K::Schur,
          binaryAlways(),
          &SymmetricEngineRing::
              schurTimesSchurViaCoefficientLittlewoodRichardson);

      addMultiplicationKernel("Schur*Complete->Schur:horizontal-Pieri",
          K::Schur,
          K::Complete,
          K::Schur,
          completeIndexIsOneRow(),
          &SymmetricEngineRing::
              schurTimesCompleteViaHorizontalPieri);
      addMultiplicationKernel("Schur*Complete->Schur:"
          "repeated-horizontal-Pieri",
          K::Schur,
          K::Complete,
          K::Schur,
          binaryAlways(),
          &SymmetricEngineRing::
              schurTimesCompleteViaRepeatedHorizontalPieri);

      addMultiplicationKernel("Schur*Elementary->Schur:vertical-Pieri",
          K::Schur,
          K::Elementary,
          K::Schur,
          elementaryIndexIsOneColumn(),
          &SymmetricEngineRing::
              schurTimesElementaryViaVerticalPieri);
      addMultiplicationKernel("Schur*Elementary->Schur:"
          "repeated-vertical-Pieri",
          K::Schur,
          K::Elementary,
          K::Schur,
          binaryAlways(),
          &SymmetricEngineRing::
              schurTimesElementaryViaRepeatedVerticalPieri);

      addMultiplicationKernel("Schur*PowerSum->Schur:Murnaghan-Nakayama",
          K::Schur,
          K::PowerSum,
          K::Schur,
          powerSumIndexIsSingleCycle(),
          &SymmetricEngineRing::
              schurTimesPowerSumViaMurnaghanNakayama);
      addMultiplicationKernel("Schur*PowerSum->Schur:"
          "repeated-Murnaghan-Nakayama",
          K::Schur,
          K::PowerSum,
          K::Schur,
          binaryAlways(),
          &SymmetricEngineRing::
              schurTimesPowerSumsViaRepeatedMurnaghanNakayama);

      auto addSchurCompatible =
          [&](const char *identifier, K first, K second) {
            addMultiplicationKernel(identifier,
                first,
                second,
                K::Schur,
                binaryAlways(),
                &SymmetricEngineRing::
                    schurCompatibleBasisTermsViaPieriAndMurnaghanNakayama);
          };
      addSchurCompatible(
          "Complete*Complete->Schur:repeated-horizontal-Pieri",
          K::Complete,
          K::Complete);
      addSchurCompatible(
          "Complete*Elementary->Schur:repeated-Pieri",
          K::Complete,
          K::Elementary);
      addSchurCompatible(
          "Complete*PowerSum->Schur:"
          "Pieri-and-Murnaghan-Nakayama",
          K::Complete,
          K::PowerSum);
      addSchurCompatible(
          "Elementary*Elementary->Schur:repeated-vertical-Pieri",
          K::Elementary,
          K::Elementary);
      addSchurCompatible(
          "Elementary*PowerSum->Schur:"
          "Pieri-and-Murnaghan-Nakayama",
          K::Elementary,
          K::PowerSum);
      addSchurCompatible(
          "PowerSum*PowerSum->Schur:"
          "repeated-Murnaghan-Nakayama",
          K::PowerSum,
          K::PowerSum);

      // Hall--Littlewood capital products.  Each endpoint uses the same
      // raising-operator/generator-product/triangular-reduction formula
      // specialized to its declared normalization and omega family.
      auto addHallLittlewoodCapital =
          [&](const char *identifier, K basis) {
            addMultiplicationKernel(identifier,
                basis,
                basis,
                basis,
                binaryAlways(),
                &SymmetricEngineRing::
                    hallLittlewoodCapitalTermsViaGeneratorTriangularFormula);
          };
      addHallLittlewoodCapital(
          "HallLittlewoodQ*HallLittlewoodQ->HallLittlewoodQ:"
          "generator-triangular-formula",
          K::HallLittlewoodQ);
      addHallLittlewoodCapital(
          "HallLittlewoodB*HallLittlewoodB->HallLittlewoodB:"
          "generator-triangular-formula",
          K::HallLittlewoodB);
      addHallLittlewoodCapital(
          "HallLittlewoodP*HallLittlewoodP->HallLittlewoodP:"
          "generator-triangular-formula",
          K::HallLittlewoodP);
      addHallLittlewoodCapital(
          "HallLittlewoodPOmega*HallLittlewoodPOmega"
          "->HallLittlewoodPOmega:generator-triangular-formula",
          K::HallLittlewoodPOmega);

      // Monomial and forgotten products. The same exponent-splitting
      // formula family is declared explicitly for each mathematical endpoint.
      auto addMonomialOrForgottenKernel =
          [&](const char *identifier,
              K first,
              K second,
              K target) {
            addMultiplicationKernel(identifier,
                first,
                second,
                target,
                binaryAlways(),
                &SymmetricEngineRing::
                    monomialAndForgottenTermsViaExponentSplittings);
          };
      addMonomialOrForgottenKernel(
          "Monomial*Monomial->Monomial:exponent-splittings",
          K::Monomial, K::Monomial, K::Monomial);
      addMonomialOrForgottenKernel(
          "Monomial*Complete->Monomial:exponent-splittings",
          K::Monomial, K::Complete, K::Monomial);
      addMonomialOrForgottenKernel(
          "Monomial*Elementary->Monomial:exponent-splittings",
          K::Monomial, K::Elementary, K::Monomial);
      addMonomialOrForgottenKernel(
          "Monomial*PowerSum->Monomial:exponent-splittings",
          K::Monomial, K::PowerSum, K::Monomial);
      addMonomialOrForgottenKernel(
          "Complete*Complete->Monomial:exponent-splittings",
          K::Complete, K::Complete, K::Monomial);
      addMonomialOrForgottenKernel(
          "Complete*Elementary->Monomial:exponent-splittings",
          K::Complete, K::Elementary, K::Monomial);
      addMonomialOrForgottenKernel(
          "Complete*PowerSum->Monomial:exponent-splittings",
          K::Complete, K::PowerSum, K::Monomial);
      addMonomialOrForgottenKernel(
          "Elementary*Elementary->Monomial:exponent-splittings",
          K::Elementary, K::Elementary, K::Monomial);
      addMonomialOrForgottenKernel(
          "Elementary*PowerSum->Monomial:exponent-splittings",
          K::Elementary, K::PowerSum, K::Monomial);
      addMonomialOrForgottenKernel(
          "PowerSum*PowerSum->Monomial:exponent-splittings",
          K::PowerSum, K::PowerSum, K::Monomial);

      addMonomialOrForgottenKernel(
          "Forgotten*Forgotten->Forgotten:exponent-splittings",
          K::Forgotten, K::Forgotten, K::Forgotten);
      addMonomialOrForgottenKernel(
          "Forgotten*Complete->Forgotten:exponent-splittings",
          K::Forgotten, K::Complete, K::Forgotten);
      addMonomialOrForgottenKernel(
          "Forgotten*Elementary->Forgotten:exponent-splittings",
          K::Forgotten, K::Elementary, K::Forgotten);
      addMonomialOrForgottenKernel(
          "Forgotten*PowerSum->Forgotten:exponent-splittings",
          K::Forgotten, K::PowerSum, K::Forgotten);
      addMonomialOrForgottenKernel(
          "Complete*Complete->Forgotten:exponent-splittings",
          K::Complete, K::Complete, K::Forgotten);
      addMonomialOrForgottenKernel(
          "Complete*Elementary->Forgotten:exponent-splittings",
          K::Complete, K::Elementary, K::Forgotten);
      addMonomialOrForgottenKernel(
          "Complete*PowerSum->Forgotten:exponent-splittings",
          K::Complete, K::PowerSum, K::Forgotten);
      addMonomialOrForgottenKernel(
          "Elementary*Elementary->Forgotten:exponent-splittings",
          K::Elementary, K::Elementary, K::Forgotten);
      addMonomialOrForgottenKernel(
          "Elementary*PowerSum->Forgotten:exponent-splittings",
          K::Elementary, K::PowerSum, K::Forgotten);
      addMonomialOrForgottenKernel(
          "PowerSum*PowerSum->Forgotten:exponent-splittings",
          K::PowerSum, K::PowerSum, K::Forgotten);
      return result;
    }();
    return definitions;
  }

// ============================================================================
// Ordered Performance Policy
// ============================================================================

const std::vector<SymmetricEngineRing::BinaryMultiplicationPicker>&
SymmetricEngineRing::binaryMultiplicationPickerDatabase()
{
    static const std::vector<BinaryMultiplicationPicker> pickers = [] {
      using K = BasisKind;
      std::vector<BinaryMultiplicationPicker> result;
      auto addMultiplicationPicker =
          [&](K first,
              K second,
              K target,
              std::vector<std::string> automatic,
              std::vector<std::string> alternatives = {}) {
            result.push_back({
                UnorderedBasisPair::from(first, second),
                target,
                {{binaryOtherwise(), std::move(automatic), {}}},
                std::move(alternatives),
                {}});
          };

      addMultiplicationPicker(K::Schur,
          K::Schur,
          K::Schur,
          {"Schur*Schur->Schur:"
           "Littlewood-Richardson-tableaux"},
          {"Schur*Schur->Schur:"
           "Littlewood-Richardson-coefficients"});
      addMultiplicationPicker(K::Schur,
          K::Complete,
          K::Schur,
          {"Schur*Complete->Schur:horizontal-Pieri",
           "Schur*Complete->Schur:"
           "repeated-horizontal-Pieri"});
      addMultiplicationPicker(K::Schur,
          K::Elementary,
          K::Schur,
          {"Schur*Elementary->Schur:vertical-Pieri",
           "Schur*Elementary->Schur:"
           "repeated-vertical-Pieri"});
      addMultiplicationPicker(K::Schur,
          K::PowerSum,
          K::Schur,
          {"Schur*PowerSum->Schur:Murnaghan-Nakayama",
           "Schur*PowerSum->Schur:"
           "repeated-Murnaghan-Nakayama"});

      addMultiplicationPicker(K::Complete,
          K::Complete,
          K::Schur,
          {"Complete*Complete->Schur:"
           "repeated-horizontal-Pieri"});
      addMultiplicationPicker(K::Complete,
          K::Elementary,
          K::Schur,
          {"Complete*Elementary->Schur:repeated-Pieri"});
      addMultiplicationPicker(K::Complete,
          K::PowerSum,
          K::Schur,
          {"Complete*PowerSum->Schur:"
           "Pieri-and-Murnaghan-Nakayama"});
      addMultiplicationPicker(K::Elementary,
          K::Elementary,
          K::Schur,
          {"Elementary*Elementary->Schur:"
           "repeated-vertical-Pieri"});
      addMultiplicationPicker(K::Elementary,
          K::PowerSum,
          K::Schur,
          {"Elementary*PowerSum->Schur:"
           "Pieri-and-Murnaghan-Nakayama"});
      addMultiplicationPicker(K::PowerSum,
          K::PowerSum,
          K::Schur,
          {"PowerSum*PowerSum->Schur:"
           "repeated-Murnaghan-Nakayama"});

      // Hall--Littlewood capital targets.
      addMultiplicationPicker(K::HallLittlewoodQ,
          K::HallLittlewoodQ,
          K::HallLittlewoodQ,
          {"HallLittlewoodQ*HallLittlewoodQ->HallLittlewoodQ:"
           "generator-triangular-formula"});
      addMultiplicationPicker(K::HallLittlewoodB,
          K::HallLittlewoodB,
          K::HallLittlewoodB,
          {"HallLittlewoodB*HallLittlewoodB->HallLittlewoodB:"
           "generator-triangular-formula"});
      addMultiplicationPicker(K::HallLittlewoodP,
          K::HallLittlewoodP,
          K::HallLittlewoodP,
          {"HallLittlewoodP*HallLittlewoodP->HallLittlewoodP:"
           "generator-triangular-formula"});
      addMultiplicationPicker(K::HallLittlewoodPOmega,
          K::HallLittlewoodPOmega,
          K::HallLittlewoodPOmega,
          {"HallLittlewoodPOmega*HallLittlewoodPOmega"
           "->HallLittlewoodPOmega:generator-triangular-formula"});

      // Monomial target.
      addMultiplicationPicker(K::Monomial, K::Monomial, K::Monomial,
          {"Monomial*Monomial->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::Monomial, K::Complete, K::Monomial,
          {"Monomial*Complete->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::Monomial, K::Elementary, K::Monomial,
          {"Monomial*Elementary->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::Monomial, K::PowerSum, K::Monomial,
          {"Monomial*PowerSum->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::Complete, K::Complete, K::Monomial,
          {"Complete*Complete->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::Complete, K::Elementary, K::Monomial,
          {"Complete*Elementary->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::Complete, K::PowerSum, K::Monomial,
          {"Complete*PowerSum->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::Elementary, K::Elementary, K::Monomial,
          {"Elementary*Elementary->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::Elementary, K::PowerSum, K::Monomial,
          {"Elementary*PowerSum->Monomial:exponent-splittings"});
      addMultiplicationPicker(K::PowerSum, K::PowerSum, K::Monomial,
          {"PowerSum*PowerSum->Monomial:exponent-splittings"});

      // Forgotten target.
      addMultiplicationPicker(K::Forgotten, K::Forgotten, K::Forgotten,
          {"Forgotten*Forgotten->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::Forgotten, K::Complete, K::Forgotten,
          {"Forgotten*Complete->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::Forgotten, K::Elementary, K::Forgotten,
          {"Forgotten*Elementary->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::Forgotten, K::PowerSum, K::Forgotten,
          {"Forgotten*PowerSum->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::Complete, K::Complete, K::Forgotten,
          {"Complete*Complete->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::Complete, K::Elementary, K::Forgotten,
          {"Complete*Elementary->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::Complete, K::PowerSum, K::Forgotten,
          {"Complete*PowerSum->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::Elementary, K::Elementary, K::Forgotten,
          {"Elementary*Elementary->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::Elementary, K::PowerSum, K::Forgotten,
          {"Elementary*PowerSum->Forgotten:exponent-splittings"});
      addMultiplicationPicker(K::PowerSum, K::PowerSum, K::Forgotten,
          {"PowerSum*PowerSum->Forgotten:exponent-splittings"});
      return result;
    }();
    return pickers;
  }

const std::map<
    std::string,
    const SymmetricEngineRing::BinaryMultiplicationKernelDefinition *>&
SymmetricEngineRing::binaryMultiplicationKernelsByIdentifier()
{
    static const std::map<
        std::string,
        const BinaryMultiplicationKernelDefinition *> index = [] {
      std::map<
          std::string,
          const BinaryMultiplicationKernelDefinition *> result;
      for (const auto& definition :
           binaryMultiplicationKernelDatabase())
        result.emplace(definition.identifier, &definition);
      return result;
    }();
    return index;
  }

const std::map<
    std::pair<SymmetricEngineRing::UnorderedBasisPair,
              SymmetricEngineRing::BasisKind>,
    const SymmetricEngineRing::BinaryMultiplicationPicker *>&
SymmetricEngineRing::binaryMultiplicationPickersByEndpoint()
{
    static const std::map<
        std::pair<UnorderedBasisPair, BasisKind>,
        const BinaryMultiplicationPicker *> index = [] {
      std::map<
          std::pair<UnorderedBasisPair, BasisKind>,
          const BinaryMultiplicationPicker *> result;
      for (const auto& picker :
           binaryMultiplicationPickerDatabase())
        result.emplace(
            std::make_pair(
                picker.factorBases, picker.targetBasis),
            &picker);
      return result;
    }();
    return index;
  }

// ============================================================================
// Database Validation
// ============================================================================

void SymmetricEngineRing::validateBinaryMultiplicationDatabase()
{
    static const bool validated = [] {
      const auto& definitions =
          binaryMultiplicationKernelDatabase();
      const auto& byIdentifier =
          binaryMultiplicationKernelsByIdentifier();
      if (definitions.size() != byIdentifier.size())
        throw std::logic_error(
            "binary multiplication has duplicate kernel identifiers");

      for (const auto& definition : definitions)
        {
          if (definition.identifier.empty() ||
              definition.kernel == nullptr ||
              definition.firstArgumentBasis == BasisKind::Custom ||
              definition.secondArgumentBasis == BasisKind::Custom ||
              definition.targetBasis == BasisKind::Custom)
            throw std::logic_error(
                "a binary multiplication kernel definition is incomplete");
          validateBinaryMultiplicationCondition(
              definition.applicableWhen,
              UnorderedBasisPair::from(
                  definition.firstArgumentBasis,
                  definition.secondArgumentBasis));
        }

      const auto& pickers =
          binaryMultiplicationPickerDatabase();
      const auto& byEndpoint =
          binaryMultiplicationPickersByEndpoint();
      if (pickers.size() != byEndpoint.size())
        throw std::logic_error(
            "binary multiplication has duplicate or reversed "
            "picker endpoints");

      std::set<std::string> reachable;
      for (const auto& picker : pickers)
        {
          if (picker.factorBases.first == BasisKind::Custom ||
              picker.factorBases.second == BasisKind::Custom ||
              picker.targetBasis == BasisKind::Custom)
            throw std::logic_error(
                "a binary multiplication picker has an invalid endpoint");
          if (picker.preferences.empty() &&
              picker.alternativeKernels.empty())
            throw std::logic_error(
                "a binary multiplication picker has no kernels");
          if (!picker.preferences.empty() &&
              picker.preferences.back().condition.kind !=
                  BinaryMultiplicationConditionKind::Otherwise)
            throw std::logic_error(
                "an automatic binary picker must end in otherwise");

          bool sawOtherwise = false;
          std::set<std::string> automaticHere;
          for (size_t i = 0; i < picker.preferences.size(); ++i)
            {
              const auto& preference =
                  picker.preferences[i];
              const bool isOtherwise =
                  preference.condition.kind ==
                  BinaryMultiplicationConditionKind::Otherwise;
              if (isOtherwise && i + 1 != picker.preferences.size())
                throw std::logic_error(
                    "otherwise must be the final binary preference");
              if (sawOtherwise)
                throw std::logic_error(
                    "a binary picker has preferences after otherwise");
              sawOtherwise = sawOtherwise || isOtherwise;
              validateBinaryMultiplicationCondition(
                  preference.condition,
                  picker.factorBases,
                  isOtherwise);
              if (preference.kernelsInOrder.empty())
                throw std::logic_error(
                    "a binary preference has no ordered kernels");
              for (const auto& identifier :
                   preference.kernelsInOrder)
                {
                  auto found = byIdentifier.find(identifier);
                  if (found == byIdentifier.end())
                    throw std::logic_error(
                        "a binary picker references an unknown kernel");
                  const auto *definition = found->second;
                  if (UnorderedBasisPair::from(
                          definition->firstArgumentBasis,
                          definition->secondArgumentBasis) !=
                          picker.factorBases ||
                      definition->targetBasis != picker.targetBasis)
                    throw std::logic_error(
                        "a binary picker references a kernel "
                        "with different endpoints");
                  if (!automaticHere.insert(identifier).second)
                    throw std::logic_error(
                        "a binary picker repeats an automatic kernel");
                  preference.kernelDefinitions.push_back(definition);
                  reachable.insert(identifier);
                }
            }

          std::set<std::string> alternativesHere;
          for (const auto& identifier :
               picker.alternativeKernels)
            {
              auto found = byIdentifier.find(identifier);
              if (found == byIdentifier.end())
                throw std::logic_error(
                    "a binary picker references an unknown alternative");
              const auto *definition = found->second;
              if (UnorderedBasisPair::from(
                      definition->firstArgumentBasis,
                      definition->secondArgumentBasis) !=
                      picker.factorBases ||
                  definition->targetBasis != picker.targetBasis)
                throw std::logic_error(
                    "a binary alternative has different endpoints");
              if (automaticHere.count(identifier) != 0 ||
                  !alternativesHere.insert(identifier).second)
                throw std::logic_error(
                    "a binary alternative is duplicated or automatic");
              picker.alternativeKernelDefinitions.push_back(
                  definition);
              reachable.insert(identifier);
            }
        }

      if (reachable.size() != definitions.size())
        throw std::logic_error(
            "a binary multiplication kernel is orphaned");

      return true;
    }();
    (void) validated;
  }

// ============================================================================
// Selection, Availability, And Execution
// ============================================================================

std::optional<
    SymmetricEngineRing::SelectedBinaryMultiplicationKernel>
SymmetricEngineRing::selectBinaryMultiplicationKernel(
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right,
    int targetBasisId,
    const BinaryMultiplicationRequest& request,
    bool traceSelection) const
{
    validateBinaryMultiplicationDatabase();
    if (binaryMultiplicationKernelExecutionDepth != 0)
      {
        ERROR("a binary multiplication kernel attempted to invoke "
              "the binary picker");
        return std::nullopt;
      }
    const auto& targetDescriptor = requireBasis(targetBasisId);
    if (error()) return std::nullopt;
    const RingBasis target{
        targetDescriptor.kind,
        targetBasisId,
        &targetDescriptor.displaySymbol,
        targetDescriptor.displayOrder};
    const BinaryMultiplicationFacts facts =
        binaryMultiplicationFacts(
            left, right, targetDescriptor.kind);
    const auto endpoint =
        std::make_pair(facts.factorBases, targetDescriptor.kind);

    auto orient =
        [&](const BinaryMultiplicationKernelDefinition *definition)
            -> std::optional<SelectedBinaryMultiplicationKernel> {
          bool reversed = false;
          if (left.basis.kind == definition->firstArgumentBasis &&
              right.basis.kind == definition->secondArgumentBasis)
            reversed = false;
          else if (
              left.basis.kind == definition->secondArgumentBasis &&
              right.basis.kind == definition->firstArgumentBasis)
            reversed = true;
          else
            return std::nullopt;

          BinaryMultiplicationFacts orientedFacts =
              reversed
                  ? binaryMultiplicationFacts(
                        right,
                        left,
                        targetDescriptor.kind)
                  : facts;
          if (!binaryMultiplicationConditionHolds(
                  definition->applicableWhen,
                  orientedFacts))
            return std::nullopt;
          return SelectedBinaryMultiplicationKernel{
              definition, target, reversed};
        };

    if (request.forcedKernel)
      {
        auto found =
            binaryMultiplicationKernelsByIdentifier().find(
                *request.forcedKernel);
        if (found == binaryMultiplicationKernelsByIdentifier().end() ||
            UnorderedBasisPair::from(
                found->second->firstArgumentBasis,
                found->second->secondArgumentBasis) !=
                facts.factorBases ||
            found->second->targetBasis != targetDescriptor.kind)
          {
            ERROR("the requested binary multiplication kernel "
                  "has incompatible endpoints: ",
                  request.forcedKernel->c_str());
            return std::nullopt;
          }
        auto selected = orient(found->second);
        if (!selected)
          ERROR("the requested binary multiplication kernel "
                "is mathematically inapplicable: ",
                request.forcedKernel->c_str());
        return selected;
      }

    auto picker =
        binaryMultiplicationPickersByEndpoint().find(endpoint);
    if (picker ==
        binaryMultiplicationPickersByEndpoint().end())
      return std::nullopt;

    for (const auto& preference :
         picker->second->preferences)
      {
        const bool applies =
            preference.condition.kind ==
                BinaryMultiplicationConditionKind::Otherwise ||
            binaryMultiplicationConditionHolds(
                preference.condition, facts);
        if (!applies) continue;
        for (const auto *definition :
             preference.kernelDefinitions)
          if (auto selected = orient(definition))
            {
              if (traceSelection &&
                  (request.traceWorkflow ||
                   basisConversionTraceEnabled()))
                std::fprintf(
                    stderr,
                    "SymmetricRings binary-kernel: endpoint={%s,%s}->%s "
                    "kernel=%s reversed=%s\n",
                    basisKindName(facts.factorBases.first),
                    basisKindName(facts.factorBases.second),
                    basisKindName(targetDescriptor.kind),
                    definition->identifier.c_str(),
                    selected->callableTakesReversedArguments
                        ? "yes"
                        : "no");
              return selected;
            }
      }
    return std::nullopt;
  }

bool
SymmetricEngineRing::automaticBinaryMultiplicationKernelAvailable(
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right,
    int targetBasisId) const
{
    return selectBinaryMultiplicationKernel(
               left,
               right,
               targetBasisId,
               BinaryMultiplicationRequest{},
               false)
        .has_value();
  }

} // namespace symmetric_rings
