// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"
#include "rings/ZZ.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Direct Power-Sum And Schur Kernels
// ============================================================================

ring_elem SymmetricEngineRing::omegaPowerSums(ring_elem f) const
{
    const auto *poly = polyValue(f);
    VECTOR(SymmetricTerm) terms;
    terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during omega");
            return zero();
          }
        long sign = ((partitionWeight(index) - partitionLength(index)) % 2 == 0) ? 1 : -1;
        terms.push_back({coefficientRing->mult(coefficientRing->from_long(sign),
                                               term.coeff),
                         term.monomial});
      }
    return fromTermVector(terms, true);
  }

ring_elem SymmetricEngineRing::schurOmegaBasisElementAsSchur(
    const Partition& alpha) const
{
    auto straightened = straightenSchurIndex(alpha);
    if (straightened.first == 0) return zero();
    int schurId = requiredBasisIdForKind(BasisKind::Schur);
    if (error()) return zero();
    ring_elem term = basisElementFromIndex(
        schurId, conjugatePartition(straightened.second));
    if (straightened.first < 0) term = negate(term);
    return term;
  }

// ============================================================================
// Metadata-Driven Basis And Monomial Mapping
// ============================================================================

ring_elem SymmetricEngineRing::omegaBasisElementDirect(
                                               const SymmetricMonomial& monomial,
                            size_t pos,
                            const std::map<int, OmegaTarget>& omegaTargets,
                            bool useSomega) const
{
    int basisId = atomBasisIdAt(monomial, pos);
    BasisKind kind = basisKindForId(basisId);
    if (basisId == registeredPowerSumBasisId())
      {
        Partition index = basisElementIndex(monomial, pos);
        long sign = ((partitionWeight(index) - partitionLength(index)) % 2 == 0) ? 1 : -1;
        ring_elem term =
            basisElementFromIndex(registeredPowerSumBasisId(), index);
        return scaled(coefficientRing->from_long(sign), term);
      }

    if (!useSomega && kind == BasisKind::Schur)
      {
        int schurId = requiredBasisIdForKind(BasisKind::Schur);
        if (error()) return zero();
        if (atomIsSkewAt(monomial, pos))
          return basisElementFromSkewIndex(
              schurId,
              conjugatePartition(basisElementOuterIndex(monomial, pos)),
              conjugatePartition(basisElementInnerIndex(monomial, pos)));
        return schurOmegaBasisElementAsSchur(
            basisElementIndex(monomial, pos));
      }

    if (!useSomega && kind == BasisKind::SchurOmega &&
        !atomIsSkewAt(monomial, pos))
      return straightenSchurBasisElement(
          basisElementIndex(monomial, pos),
          requiredBasisIdForKind(BasisKind::Schur));

    auto target = omegaTargets.find(basisId);
    if (target != omegaTargets.end())
      {
        Partition payload = basisElementIndex(monomial, pos);
        const auto& targetBasis = requireBasis(target->second.basisId);
        auto *poly = new SymmetricRingPoly;
        SymmetricMonomial targetMonomial;
        appendAtomBlock(targetMonomial,
                        makeAtomBlock(targetBasis.displayOrder,
                                      target->second.basisId,
                                      atomInnerLengthAt(monomial, pos),
                                      payload));
        poly->terms.push_back({coefficientRing->one(),
                               canonicalMonomial(targetMonomial)});
        return makePolyValue(poly);
      }

    // Metadata mappings above are exact omega formulas. When no such formula
    // exists, make the basis element explicit and let the shared conversion
    // registry own the route to power sums.
    ring_elem basisElement = expressionFromAtom(monomial, pos);
    int powerSumBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return zero();
    ring_elem inPowerSums = toBasis(basisElement, powerSumBasisId);
    if (error()) return zero();
    return omegaPowerSums(inPowerSums);
  }

ring_elem SymmetricEngineRing::omegaMonomial(const SymmetricMonomial& monomial,
                          const std::map<int, OmegaTarget>& omegaTargets,
                          bool useSomega) const
{
    ring_elem result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        ring_elem factor =
            omegaBasisElementDirect(monomial, pos, omegaTargets, useSomega);
        if (error()) return zero();
        result = mult(result, factor);
        pos += atomLengthAt(monomial, pos);
      }
    return result;
  }

std::map<int, OmegaTarget> SymmetricEngineRing::omegaTargetMap(M2_arrayint omegaMap) const
{
    std::map<int, OmegaTarget> result;
    if (omegaMap == nullptr) return result;
    if (omegaMap->len % 2 != 0)
      {
        ERROR("invalid omega metadata map");
        return result;
      }
    for (int i = 0; i < omegaMap->len; i += 2)
      result[omegaMap->array[i]] =
          OmegaTarget{omegaMap->array[i + 1]};
    return result;
  }

// ============================================================================
// Public Omega Entry Point
// ============================================================================

ring_elem SymmetricEngineRing::omegaInvolution(ring_elem f, M2_arrayint omegaMap, bool useSomega) const
{
    std::map<int, OmegaTarget> targets = omegaTargetMap(omegaMap);
    if (error()) return zero();
    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        ring_elem converted = omegaMonomial(term.monomial, targets, useSomega);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, converted));
      }
    return result;
  }

} // namespace symmetric_rings
