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
// Adams-Operation Plethysm
// ============================================================================

ring_elem SymmetricEngineRing::powerSumsViaAdamsOperation(ring_elem f, int multiplier) const
{
    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during plethysm");
            return zero();
          }
        Partition transformed;
        transformed.reserve(index.size());
        for (int part : index)
          if (part * multiplier > 0) transformed.push_back(part * multiplier);
        ring_elem termElement = powerSumElementFromIndex(transformed);
        result = add(result, scaled(term.coeff, termElement));
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumPlethysmViaAdamsOperations(ring_elem fPowerSums, ring_elem gPowerSums) const
{
    const auto *poly = polyValue(fPowerSums);
    GCMap<int, ring_elem> adamsCache;
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during plethysm");
            return zero();
          }
        ring_elem substituted = one();
        for (int part : index)
          {
            auto cached = adamsCache.find(part);
            if (cached == adamsCache.end())
              {
                ring_elem factor = powerSumsViaAdamsOperation(gPowerSums, part);
                if (error()) return zero();
                cached = adamsCache.emplace(part, factor).first;
              }
            ring_elem factor = copyPolyValue(polyValue(cached->second));
            substituted = mult(substituted, factor);
          }
        result = add(result, scaled(term.coeff, substituted));
      }
    return result;
  }

// ============================================================================
// Specialized Schur Plethysm
// ============================================================================

bool SymmetricEngineRing::singleSchurPartition(ring_elem f, int schurId, Partition& lambda) const
{
    const auto *poly = polyValue(f);
    if (poly->terms.size() != 1) return false;
    if (!coefficientRing->is_equal(poly->terms[0].coeff, coefficientRing->one()))
      return false;
    if (!singleBasisIndexFromMonomial(poly->terms[0].monomial, schurId, lambda))
      return false;
    return isPartitionIndex(lambda);
  }

// Specialized Schur-plethysm path used for one-row inner Schur functions.
std::string SymmetricEngineRing::completePlethysmCacheKey(
                                       const std::string& algorithm,
                                       int n,
                                       int schurId,
                                       const Partition& inner) const
{
    return algorithm + "|" + std::to_string(n) + "|" +
           std::to_string(schurId) + "|" + partitionKey(inner);
  }

ring_elem SymmetricEngineRing::completePlethysmViaAdamsRecurrence(int n,
                                  const Partition& inner,
                                  ring_elem innerPowerSums,
                                  int schurId,
                                  const std::string& schurDisplay,
                                  int schurOrder) const
{
    if (n < 0) return zero();
    if (n == 0) return one();

    std::string key = completePlethysmCacheKey("adams", n, schurId, inner);
    auto cached = schurCompletePlethysmCache.find(key);
    if (cached != schurCompletePlethysmCache.end())
      return copyPolyValue(polyValue(cached->second));

    ring_elem total = zero();
    for (int i = 1; i <= n; ++i)
      {
        ring_elem pIAtInner = powerSumsViaAdamsOperation(innerPowerSums, i);
        if (error()) return zero();
        // In this specialized route the retained inner Schur function is a
        // short row, so every Adams image has only partitionCount(|inner|)
        // power-sum terms even though its total degree is i*|inner|.  Applying
        // rim hooks to that sparse support avoids constructing and scanning
        // character data for the much larger dilated degree.
        ring_elem pIAtInnerSchur = powerSumsToSchurViaAbacusRimHooks(
            pIAtInner, schurId, schurDisplay, schurOrder);
        if (error()) return zero();
        ring_elem rest = completePlethysmViaAdamsRecurrence(n - i,
                                                      inner,
                                                      innerPowerSums,
                                                      schurId,
                                                      schurDisplay,
                                                      schurOrder);
        if (error()) return zero();
        ring_elem product =
            multiplySchurExpansionsViaLittlewoodRichardson(pIAtInnerSchur, rest, schurId, schurDisplay, schurOrder);
        if (error()) return zero();
        total = add(total, product);
      }

    ring_elem reciprocal = rationalCoefficient(1, n);
    if (error()) return zero();
    ring_elem result = scaled(reciprocal, total);
    schurCompletePlethysmCache[key] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::schurPlethysmToSchurViaAdamsJacobiTrudi(const Partition& outer,
                                     const Partition& inner,
                                     int schurId,
                                     const std::string& schurDisplay,
                                     int schurOrder) const
{
    size_t n = outer.size();
    if (n == 0) return one();
    if (n >= 8 * sizeof(size_t))
      {
        ERROR("Jacobi-Trudi determinant is too large");
        return zero();
      }

    ring_elem innerPowerSums = schurLikeToPowerSumsViaCharacters(inner, false);
    if (error()) return zero();

    RingElemMatrix matrix(n, RingElemVector(n));
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        {
          int degree = outer[i] - static_cast<int>(i) + static_cast<int>(j);
          matrix[i][j] = completePlethysmViaAdamsRecurrence(degree,
                                                      inner,
                                                      innerPowerSums,
                                                      schurId,
                                                      schurDisplay,
                                                      schurOrder);
          if (error()) return zero();
        }

    size_t limit = static_cast<size_t>(1) << n;
    RingElemVector dp;
    dp.reserve(limit);
    for (size_t i = 0; i < limit; ++i) dp.push_back(zero());
    dp[0] = one();

    for (size_t mask = 0; mask < limit; ++mask)
      {
        if (is_zero(dp[mask])) continue;
        int row = popcountMask(mask);
        if (row >= static_cast<int>(n)) continue;
        for (size_t col = 0; col < n; ++col)
          {
            size_t bit = static_cast<size_t>(1) << col;
            if ((mask & bit) != 0) continue;
            if (is_zero(matrix[row][col])) continue;
            ring_elem term = multiplySchurExpansionsViaLittlewoodRichardson(dp[mask],
                                                   matrix[row][col],
                                                   schurId,
                                                   schurDisplay,
                                                   schurOrder);
            if (error()) return zero();
            if (selectedGreaterThan(mask, col, n) % 2 == 1) term = negate(term);
            size_t next = mask | bit;
            dp[next] = add(dp[next], term);
          }
      }

    return dp[limit - 1];
  }

bool SymmetricEngineRing::trySchurPlethysmToSchurViaAdamsJacobiTrudi(ring_elem f,
                            ring_elem g,
                            int targetBasisId,
                            const std::string& targetDisplay,
                            int targetOrder,
                            ring_elem& result) const
{
    if (!hasBasisKind(targetBasisId, BasisKind::Schur)) return false;
    Partition outer;
    Partition inner;
    if (!singleSchurPartition(f, targetBasisId, outer)) return false;
    if (!singleSchurPartition(g, targetBasisId, inner)) return false;
    if (inner.size() != 1) return false;
    // The guaranteed post-plethysm power-sum pipeline is faster for one-row
    // outer input.
    if (outer.size() == 1) return false;
    // Adams/Jacobi-Trudi grows both with the determinant dimension and with
    // the Adams dilation supplied by the one-row inner shape. Crossover
    // sweeps give a conservative stable region: two-row outers through inner
    // size four, and three-row outers for inner size two. Outside that region
    // materialization followed by the shared p -> S dispatcher avoids sharp
    // high-inner and four-row regressions.
    int innerPart = inner.front();
    if (outer.size() == 2 && innerPart > 4) return false;
    if (outer.size() == 3 && innerPart > 2) return false;
    if (outer.size() > 3) return false;
    result = schurPlethysmToSchurViaAdamsJacobiTrudi(outer,
                                      inner,
                                      targetBasisId,
                                      targetDisplay,
                                      targetOrder);
    return !error();
  }

// ============================================================================
// Public Plethysm Entry Points
// ============================================================================

ring_elem SymmetricEngineRing::plethysm(ring_elem f, ring_elem g) const
{
    int pBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    ring_elem fPowerSums = toBasis(f, pBasisId);
    if (error()) return zero();
    ring_elem gPowerSums = toBasis(g, pBasisId);
    if (error()) return zero();
    ring_elem result = powerSumPlethysmViaAdamsOperations(fPowerSums, gPowerSums);
    if (error()) return zero();

    SymmetricConversionMetadata metadata;
    metadata.pureBasis = pBasisId;
    metadata.expandedBasis = pBasisId;
    metadata.termCount = polyValue(result)->terms.size();
    metadata.factorBases = std::vector<int>{pBasisId};
    metadata.singleTerm = polyValue(result)->terms.size() == 1;
    metadata.noProducts = true;
    metadata.singleBasisElement = polyValue(result)->terms.size() == 1 &&
                          !polyValue(result)->terms[0].monomial.data.empty();
    metadata.normalized = true;
    metadata.skewFree = true;
    metadata.collected = true;
    int outerWeight = elementWeight(f);
    int innerWeight = elementWeight(g);
    if (outerWeight >= 0 && innerWeight >= 0)
      metadata.homogeneousWeight = outerWeight * innerWeight;
    auto *resultPoly = mutablePolyValue(result);
    resultPoly->combinatorialTags =
        combinatorialTagMask(CombinatorialTag::Plethysm);
    resultPoly->conversionMetadata = metadata;
    return result;
  }

ring_elem SymmetricEngineRing::plethysmToBasisDispatch(ring_elem f,
                            ring_elem g,
                            int targetBasisId) const
{
    int pBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    const auto& powerSums = requireBasis(pBasisId);
    const auto& target = requireBasis(targetBasisId);

    ConversionGuarantees guarantees;
    ConversionGuarantees outerGuarantees =
        inferConversionGuarantees(f, targetBasisId);
    ConversionGuarantees innerGuarantees =
        inferConversionGuarantees(g, targetBasisId);
    if (outerGuarantees.homogeneousWeight &&
        innerGuarantees.homogeneousWeight)
      guarantees.homogeneousWeight =
          *outerGuarantees.homogeneousWeight *
          *innerGuarantees.homogeneousWeight;
    ConversionInput input{
        zero(),
        std::move(guarantees),
        combinatorialTagMask(CombinatorialTag::Plethysm)};
    ConversionRequest request{
        ConversionRequestKind::PostPlethysm,
        input,
        f,
        g};
    return conversionRequestToBasisDispatch(request,
                                            pBasisId,
                                            powerSums.displaySymbol,
                                            powerSums.displayOrder,
                                            powerSums.multiplicative,
                                            targetBasisId,
                                            target.displaySymbol,
                                            target.displayOrder,
                                            target.multiplicative);
  }
} // namespace symmetric_rings
