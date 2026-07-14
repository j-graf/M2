// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"

namespace symmetric_rings {

// ============================================================================
// Basis-Expansion Inspection
// ============================================================================
// Stored expressions are recognized as ordinary basis expansions or scaled basis elements.

bool SymmetricEngineRing::singleScaledBasisElement(
    ring_elem f,
    int basisId,
    Partition& index,
    ring_elem& coefficient) const
{
    if (basisId < 0) return false;
    const auto *poly = polyValue(f);
    if (poly->terms.size() != 1) return false;
    if (!singleBasisIndexFromMonomial(
            poly->terms.front().monomial, basisId, index))
      return false;
    coefficient = poly->terms.front().coeff;
    return isPartitionIndex(index);
  }


CoeffMap SymmetricEngineRing::coefficientsInBasis(ring_elem f, int basisId) const
{
    CoeffMap result;
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (term.monomial.data.empty())
          index = Partition{};
        else if (!singleBasisIndexFromMonomial(term.monomial, basisId, index))
          {
            ERROR("expected an expression in a single symmetric-function basis");
            return CoeffMap{};
          }
        addCoeff(result, index, term.coeff);
      }
    return result;
  }

bool SymmetricEngineRing::coefficientsInBasisIfPossible(
    ring_elem f,
    int basisId,
    CoeffMap& result) const
{
    result.clear();
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (term.monomial.data.empty())
          index = Partition{};
        else if (!singleBasisIndexFromMonomial(term.monomial, basisId, index))
          {
            result.clear();
            return false;
          }
        addCoeff(result, index, term.coeff);
      }
    return true;
  }



} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
