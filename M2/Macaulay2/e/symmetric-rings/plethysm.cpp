// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"
#include "rings/ZZ.hpp"

#include <algorithm>
#include <cstdio>
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

bool SymmetricEngineRing::
schurPlethysmToSchurViaAdamsJacobiTrudiApplicable(
    ring_elem f,
    ring_elem g,
    int targetBasisId,
    Partition& outer,
    Partition& inner) const
{
    if (!hasBasisKind(targetBasisId, BasisKind::Schur)) return false;
    if (!singleSchurPartition(f, targetBasisId, outer)) return false;
    if (!singleSchurPartition(g, targetBasisId, inner)) return false;
    if (inner.size() != 1) return false;
    // The materialized power-sum route is faster for one-row outer input.
    if (outer.size() == 1) return false;

    // Adams/Jacobi-Trudi grows both with the determinant dimension and with
    // the Adams dilation supplied by the one-row inner shape. Crossover
    // sweeps give a conservative stable region: two-row outers through inner
    // size four, and three-row outers for inner size two.
    int innerPart = inner.front();
    if (outer.size() == 2 && innerPart > 4) return false;
    if (outer.size() == 3 && innerPart > 2) return false;
    return outer.size() <= 3;
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
        // The realized Adams image is an ordinary canonical power-sum
        // expansion. Preserve its plethysm provenance and let the shared
        // p -> Schur picker choose and execute the complete conversion plan.
        mutablePolyValue(pIAtInner)->combinatorialTags =
            combinatorialTagMask(CombinatorialTag::Plethysm);
        ring_elem pIAtInnerSchur = toBasis(pIAtInner, schurId);
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
    requireCacheEntryCapacity("Schur plethysm cache");
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

    ring_elem innerSchur = basisElementFromIndex(schurId, inner);
    mutablePolyValue(innerSchur)->combinatorialTags =
        combinatorialTagMask(CombinatorialTag::Plethysm);
    const int powerSumId =
        requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return zero();
    ring_elem innerPowerSums =
        toBasis(innerSchur, powerSumId);
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
    if (!determinantStatesWithinLimit(
            limit, "plethysm Jacobi-Trudi determinant"))
      return zero();
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

// ============================================================================
// Plethysm-To-Basis Selection And Execution
// ============================================================================
// A combined request makes exactly one complete decision before performing
// algebra. The broad route materializes the canonical power-sum plethysm and
// sends its only conversion decision to the shared conversion plan database.

SymmetricEngineRing::PlethysmToBasisRoute
SymmetricEngineRing::selectPlethysmToBasisRoute(
    ring_elem f,
    ring_elem g,
    int targetBasisId) const
{
    Partition outer;
    Partition inner;
    if (schurPlethysmToSchurViaAdamsJacobiTrudiApplicable(
            f, g, targetBasisId, outer, inner))
      return PlethysmToBasisRoute::ViaSchurAdamsJacobiTrudi;
    return PlethysmToBasisRoute::ViaPowerSumsThenBasisConversion;
  }

const char *SymmetricEngineRing::plethysmToBasisRouteName(
    PlethysmToBasisRoute route) const
{
    switch (route)
      {
        case PlethysmToBasisRoute::ViaSchurAdamsJacobiTrudi:
          return "S->S:Adams-Jacobi-Trudi";
        case PlethysmToBasisRoute::ViaPowerSumsThenBasisConversion:
          return "p-materialization->target:default-conversion";
      }
    return "unknown";
  }

void SymmetricEngineRing::tracePlethysmToBasisSelection(
    PlethysmToBasisRoute route,
    int targetBasisId) const
{
    if (!basisConversionTraceEnabled())
      return;
    std::fprintf(
        stderr,
        "SymmetricRings plethysm-to-basis: target=%s route=%s\n",
        basisKeyForId(targetBasisId).c_str(),
        plethysmToBasisRouteName(route));
  }

ring_elem SymmetricEngineRing::executePlethysmToBasisRoute(
    PlethysmToBasisRoute route,
    ring_elem f,
    ring_elem g,
    int targetBasisId) const
{
    if (route == PlethysmToBasisRoute::ViaSchurAdamsJacobiTrudi)
      {
        Partition outer;
        Partition inner;
        if (!schurPlethysmToSchurViaAdamsJacobiTrudiApplicable(
                f, g, targetBasisId, outer, inner))
          {
            ERROR("a selected plethysm route violated its applicability contract");
            return zero();
          }
        const auto& target = requireBasis(targetBasisId);
        ring_elem result = schurPlethysmToSchurViaAdamsJacobiTrudi(
            outer,
            inner,
            targetBasisId,
            target.displaySymbol,
            target.displayOrder);
        if (error()) return zero();
        ExpressionFacts facts = inferCanonicalExpansionFacts(
            result,
            targetBasisId,
            partitionWeight(outer) * partitionWeight(inner));
        if (error()) return zero();
        attachExpressionFacts(
            result,
            facts,
            targetBasisId,
            combinatorialTagMask(CombinatorialTag::Plethysm));
        return result;
      }

    ring_elem powerSums = plethysm(f, g);
    if (error()) return zero();
    return toBasis(powerSums, targetBasisId);
  }

// ============================================================================
// Public Plethysm Entry Points
// ============================================================================

void SymmetricEngineRing::requirePlethysmWithinWeightLimit(
    ring_elem f, ring_elem g) const
{
    auto maximumAbsoluteWeight = [](const SymmetricRingPoly *poly) {
      long long maximum = 0;
      for (const auto& term : poly->terms)
        {
          const long long weight = monomialWeight(term.monomial);
          const long long magnitude = weight < 0 ? -weight : weight;
          maximum = std::max(maximum, magnitude);
        }
      return maximum;
    };
    const long long outerWeight = maximumAbsoluteWeight(polyValue(f));
    const long long innerWeight = maximumAbsoluteWeight(polyValue(g));
    requireWeightWithinLimit(
        outerWeight * innerWeight, "symmetric-function plethysm");
  }

ring_elem SymmetricEngineRing::plethysm(ring_elem f, ring_elem g) const
{
    requirePlethysmWithinWeightLimit(f, g);
    int pBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    ring_elem fPowerSums = toBasis(f, pBasisId);
    if (error()) return zero();
    ring_elem gPowerSums = toBasis(g, pBasisId);
    if (error()) return zero();
    ring_elem result = powerSumPlethysmViaAdamsOperations(fPowerSums, gPowerSums);
    if (error()) return zero();

    std::optional<int> resultWeight;
    int outerWeight = elementWeight(f);
    int innerWeight = elementWeight(g);
    if (outerWeight >= 0 && innerWeight >= 0)
      resultWeight = outerWeight * innerWeight;
    ExpressionFacts resultFacts = inferCanonicalExpansionFacts(
        result, pBasisId, resultWeight);
    if (error()) return zero();
    const CombinatorialTags tags =
        combinatorialTagMask(CombinatorialTag::Plethysm);
    auto *resultPoly = mutablePolyValue(result);
    resultPoly->combinatorialTags = tags;
    resultFacts.combinatorialTags = tags;
    attachExpressionFacts(result, resultFacts, pBasisId, tags);
    return result;
  }

ring_elem SymmetricEngineRing::plethysmToBasisDispatch(ring_elem f,
                            ring_elem g,
                            int targetBasisId) const
{
    requirePlethysmWithinWeightLimit(f, g);
    requireBasis(targetBasisId);
    if (error()) return zero();
    PlethysmToBasisRoute route =
        selectPlethysmToBasisRoute(f, g, targetBasisId);
    tracePlethysmToBasisSelection(route, targetBasisId);
    return executePlethysmToBasisRoute(route, f, g, targetBasisId);
  }
} // namespace symmetric_rings
