// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>

namespace symmetric_rings {

// This file owns explicit mathematical closure declarations for complete
// product terms. It does not execute multiplication or choose binary kernels.

// ============================================================================
// Mathematical Family Database
// ============================================================================

const std::vector<
    SymmetricEngineRing::TargetClosedMultiplicationFamily>&
SymmetricEngineRing::targetClosedMultiplicationFamilyDatabase()
{
    using K = BasisKind;
    static const std::vector<TargetClosedMultiplicationFamily> families{
        {
            "Schur-compatible",
            K::Schur,
            {K::Schur, K::Elementary, K::Complete, K::PowerSum}
        },
        {
            "HallLittlewoodQ-same-basis",
            K::HallLittlewoodQ,
            {K::HallLittlewoodQ}
        },
        {
            "HallLittlewoodB-same-basis",
            K::HallLittlewoodB,
            {K::HallLittlewoodB}
        },
        {
            "HallLittlewoodP-same-basis",
            K::HallLittlewoodP,
            {K::HallLittlewoodP}
        },
        {
            "HallLittlewoodPOmega-same-basis",
            K::HallLittlewoodPOmega,
            {K::HallLittlewoodPOmega}
        }};
    return families;
  }

// ============================================================================
// Closure-Contract Validation
// ============================================================================

void
SymmetricEngineRing::validateTargetClosedMultiplicationFamilyDatabase()
{
    static const bool validated = [] {
      validateBinaryMultiplicationDatabase();
      const auto& binaryPickers =
          binaryMultiplicationPickersByEndpoint();
      std::set<std::string> identifiers;
      std::set<BasisKind> targets;

      for (const auto& family :
           targetClosedMultiplicationFamilyDatabase())
        {
          if (family.identifier.empty() ||
              family.targetBasis == BasisKind::Custom ||
              family.factorBases.empty())
            throw std::logic_error(
                "a target-closed multiplication family is incomplete");
          if (!identifiers.insert(family.identifier).second)
            throw std::logic_error(
                "target-closed multiplication families have duplicate "
                "identifiers");
          if (!targets.insert(family.targetBasis).second)
            throw std::logic_error(
                "a multiplication target has more than one target-closed "
                "family");

          std::set<BasisKind> factorBases;
          for (BasisKind factorBasis : family.factorBases)
            {
              if (factorBasis == BasisKind::Custom ||
                  !factorBases.insert(factorBasis).second)
                throw std::logic_error(
                    "a target-closed multiplication family has an invalid "
                    "or repeated factor basis");

              const auto endpoint = std::make_pair(
                  UnorderedBasisPair::from(
                      family.targetBasis, factorBasis),
                  family.targetBasis);
              auto picker = binaryPickers.find(endpoint);
              if (picker == binaryPickers.end())
                throw std::logic_error(
                    "a target-closed multiplication family lacks a required "
                    "binary endpoint");

              // A fold cannot recover through the broad fallback after it
              // starts. The final otherwise preference must therefore contain
              // a kernel whose own applicability is unconditional.
              const auto& preferences =
                  picker->second->preferences;
              if (preferences.empty() ||
                  preferences.back().condition.kind !=
                      BinaryMultiplicationConditionKind::Otherwise)
                throw std::logic_error(
                    "a target-closed multiplication family lacks a total "
                    "automatic binary preference");
              const bool hasUnconditionalKernel =
                  std::any_of(
                      preferences.back().
                          kernelDefinitions.begin(),
                      preferences.back().
                          kernelDefinitions.end(),
                      [](const BinaryMultiplicationKernelDefinition
                             *definition) {
                        return definition->applicableWhen.kind ==
                            BinaryMultiplicationConditionKind::Always;
                      });
              if (!hasUnconditionalKernel)
                throw std::logic_error(
                    "a target-closed multiplication family lacks a total "
                    "automatic binary kernel");
            }

          if (factorBases.count(family.targetBasis) == 0)
            throw std::logic_error(
                "a target-closed multiplication family must contain its "
                "target basis");
        }
      return true;
    }();
    (void) validated;
  }

// ============================================================================
// Complete-Factor-List Selection
// ============================================================================

const SymmetricEngineRing::TargetClosedMultiplicationFamily *
SymmetricEngineRing::selectTargetClosedMultiplicationFamily(
    const std::vector<int>& factorBasisIds,
    int targetBasisId) const
{
    validateTargetClosedMultiplicationFamilyDatabase();
    const BasisKind targetBasis =
        basisKindForId(targetBasisId);
    for (const auto& family :
         targetClosedMultiplicationFamilyDatabase())
      {
        if (family.targetBasis != targetBasis)
          continue;
        const bool everyFactorBelongs =
            std::all_of(
                factorBasisIds.begin(),
                factorBasisIds.end(),
                [&](int factorBasisId) {
                  const BasisKind factorBasis =
                      basisKindForId(factorBasisId);
                  return std::find(
                             family.factorBases.begin(),
                             family.factorBases.end(),
                             factorBasis) !=
                         family.factorBases.end();
                });
        return everyFactorBelongs ? &family : nullptr;
      }
    return nullptr;
  }

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
