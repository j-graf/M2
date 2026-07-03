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

ring_elem SymmetricEngineRing::adamsPowerSums(ring_elem f, int multiplier) const
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

ring_elem SymmetricEngineRing::plethysmPowerSums(ring_elem fPowerSums, ring_elem gPowerSums) const
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
                ring_elem factor = adamsPowerSums(gPowerSums, part);
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

int SymmetricEngineRing::singleBasisIdInMonomial(const SymmetricMonomial& monomial) const
{
    int result = 0;
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        int basisId = atomBasisIdAt(monomial, pos);
        if (result == 0)
          result = basisId;
        else if (result != basisId)
          return -1;
        pos += atomLengthAt(monomial, pos);
      }
    return result;
  }

int SymmetricEngineRing::singleBasisId(ring_elem f) const
{
    const auto *poly = polyValue(f);
    int result = 0;
    for (const auto& term : poly->terms)
      {
        int termBasisId = singleBasisIdInMonomial(term.monomial);
        if (termBasisId < 0) return -1;
        if (termBasisId == 0) continue;
        if (result == 0)
          result = termBasisId;
        else if (result != termBasisId)
          return -1;
      }
    return result;
  }

bool SymmetricEngineRing::powerSumIndexFromMonomial(const SymmetricMonomial& monomial,
                                 Partition& index) const
{
    index.clear();
    if (monomial.data.empty()) return true;
    if (!isSinglePowerSumBlock(monomial)) return false;
    index = atomIndex(monomial, 0);
    return true;
  }

ring_elem SymmetricEngineRing::powerSumMonomialToTarget(const Partition& index,
                                     const std::string& targetDisplay,
                                     int targetBasisId,
                                     int targetDisplayOrder,
                                     bool targetIsMultiplicative) const
{
    if (targetDisplay == "p")
      return basisElementFromIndex(targetBasisId,
                                   targetDisplay,
                                   targetDisplayOrder,
                                   targetIsMultiplicative,
                                   index);

    if (targetDisplay == "h")
      {
        ring_elem result = one();
        for (int part : index)
          result = mult(result, powerSumPartToComplete(part, targetBasisId, targetDisplayOrder));
        return result;
      }

    if (targetDisplay == "e")
      {
        ring_elem result = one();
        for (int part : index)
          result = mult(result, powerSumPartToElementary(part, targetBasisId, targetDisplayOrder));
        return result;
      }

    if (targetDisplay == "q" || targetDisplay == "b")
      {
        ring_elem result = one();
        bool omega = targetDisplay == "b";
        for (int part : index)
          result = mult(result,
                        powerSumPartToHallGenerator(part,
                                                    targetBasisId,
                                                    targetDisplay,
                                                    targetDisplayOrder,
                                                    omega));
        return result;
      }

    if (targetDisplay == "m" || isForgottenDisplay(targetDisplay))
      return powerSumIndexToMonomialTarget(index,
                                           targetBasisId,
                                           targetDisplay,
                                           targetDisplayOrder,
                                           isForgottenDisplay(targetDisplay));

    if (targetDisplay == "S")
      return powerSumToSchur(index, targetBasisId, targetDisplayOrder);

    if (targetDisplay == "Somega")
      {
        long sign = ((partitionWeight(index) - partitionLength(index)) % 2 == 0) ? 1 : -1;
        return powerSumToSchurLike(index,
                                   targetBasisId,
                                   targetDisplayOrder,
                                   "Somega",
                                   sign);
      }

    ERROR("basis conversion from power sums is not implemented for basis ", targetDisplay.c_str());
    return zero();
  }

ring_elem SymmetricEngineRing::powerSumsToTarget(ring_elem f,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetDisplayOrder,
                              bool targetIsMultiplicative) const
{
    if (targetDisplay == "Q" || targetDisplay == "B" ||
        targetDisplay == "P" || targetDisplay == "R")
      return powerSumsToHallCapitalTarget(f,
                                          targetBasisId,
                                          targetDisplay,
                                          targetDisplayOrder);
    if (targetDisplay == "S")
      return powerSumsToSchurLike(f,
                                  targetBasisId,
                                  targetDisplayOrder,
                                  targetDisplay,
                                  false);
    if (targetDisplay == "Somega")
      return powerSumsToSchurLike(f,
                                  targetBasisId,
                                  targetDisplayOrder,
                                  targetDisplay,
                                  true);

    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during basis conversion");
            return zero();
          }
        ring_elem converted = powerSumMonomialToTarget(index,
                                                       targetDisplay,
                                                       targetBasisId,
                                                       targetDisplayOrder,
                                                       targetIsMultiplicative);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, converted));
      }
    return result;
  }

bool SymmetricEngineRing::atomToDirectTarget(const SymmetricMonomial& monomial,
                          size_t pos,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          bool targetIsMultiplicative,
                          ring_elem& result) const
{
    int basisId = atomBasisIdAt(monomial, pos);
    std::string display = displayForBasis(basisId);
    Partition index = atomIndex(monomial, pos);

    if (!atomIsSkewAt(monomial, pos) && display == targetDisplay)
      {
        result = basisElementFromIndex(targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       targetIsMultiplicative,
                                       index);
        return true;
      }

    if (targetDisplay == "h" && display == "S")
      {
        Partition outer = atomIsSkewAt(monomial, pos) ? atomOuterIndex(monomial, pos)
                                                      : index;
        Partition inner = atomIsSkewAt(monomial, pos) ? atomInnerIndex(monomial, pos)
                                                      : Partition{};
        result = jacobiTrudi(outer,
                             inner,
                             targetBasisId,
                             targetDisplay,
                             targetDisplayOrder,
                             targetIsMultiplicative);
        return true;
      }

    if (targetDisplay == "e" && display == "Somega")
      {
        Partition outer = atomIsSkewAt(monomial, pos) ? atomOuterIndex(monomial, pos)
                                                      : index;
        Partition inner = atomIsSkewAt(monomial, pos) ? atomInnerIndex(monomial, pos)
                                                      : Partition{};
        result = jacobiTrudi(outer,
                             inner,
                             targetBasisId,
                             targetDisplay,
                             targetDisplayOrder,
                             targetIsMultiplicative);
        return true;
      }

    if (atomIsSkewAt(monomial, pos)) return false;
    if (display != "p" && display != "h" && display != "e" &&
        display != "q" && display != "b" && display != "m" &&
        !isForgottenDisplay(display) && display != "S" && display != "Somega" &&
        display != "Q" && display != "B" && display != "P" && display != "R")
      return false;
    ring_elem inPowerSums = atomToPowerSums(monomial, pos);
    if (error()) return false;
    result = powerSumsToTarget(inPowerSums,
                               targetBasisId,
                               targetDisplay,
                               targetDisplayOrder,
                               targetIsMultiplicative);
    return !error();
  }

bool SymmetricEngineRing::monomialToDirectTarget(const SymmetricMonomial& monomial,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetDisplayOrder,
                              bool targetIsMultiplicative,
                              ring_elem& result) const
{
    result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        ring_elem factor;
        if (!atomToDirectTarget(monomial,
                                pos,
                                targetBasisId,
                                targetDisplay,
                                targetDisplayOrder,
                                targetIsMultiplicative,
                                factor))
          return false;
        result = mult(result, factor);
        pos += atomLengthAt(monomial, pos);
      }
    return true;
  }

bool SymmetricEngineRing::elementToDirectTarget(ring_elem f,
                             int targetBasisId,
                             const std::string& targetDisplay,
                             int targetDisplayOrder,
                             bool targetIsMultiplicative,
                             ring_elem& result) const
{
    if (directTriangularSchurConversion(f,
                                        targetBasisId,
                                        targetDisplay,
                                        targetDisplayOrder,
                                        result))
      return true;
    if (error()) return false;
    if (directTriangularHallConversion(f,
                                       targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       result))
      return true;
    if (error()) return false;

    if (targetDisplay != "h" && targetDisplay != "e" && !targetIsMultiplicative)
      return false;
    result = zero();
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        ring_elem converted;
        if (!monomialToDirectTarget(term.monomial,
                                    targetBasisId,
                                    targetDisplay,
                                    targetDisplayOrder,
                                    targetIsMultiplicative,
                                    converted))
          return false;
        result = add(result, scaled(term.coeff, converted));
      }
    return true;
  }

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

std::string SymmetricEngineRing::schurCompletePlethysmKey(int n,
                                       int schurId,
                                       const Partition& inner) const
{
    return std::to_string(n) + "|" + std::to_string(schurId) + "|" +
           partitionKey(inner);
  }

ring_elem SymmetricEngineRing::schurCompletePlethysm(int n,
                                  const Partition& inner,
                                  int powerSumBasisId0,
                                  const std::string& powerSumDisplay,
                                  int powerSumOrder,
                                  bool powerSumIsMultiplicative,
                                  int schurId,
                                  const std::string& schurDisplay,
                                  int schurOrder) const
{
    if (n < 0) return zero();
    if (n == 0) return one();

    std::string key = schurCompletePlethysmKey(n, schurId, inner);
    auto cached = schurCompletePlethysmCache.find(key);
    if (cached != schurCompletePlethysmCache.end())
      return copyPolyValue(polyValue(cached->second));

    ring_elem innerSchur =
        basisElementFromIndex(schurId, schurDisplay, schurOrder, false, inner);
    ring_elem total = zero();
    for (int i = 1; i <= n; ++i)
      {
        ring_elem pI = basisElementFromIndex(powerSumBasisId0,
                                             powerSumDisplay,
                                             powerSumOrder,
                                             powerSumIsMultiplicative,
                                             Partition{i});
        ring_elem pIAtInner = plethysm(pI,
                                       innerSchur,
                                       powerSumBasisId0,
                                       powerSumDisplay,
                                       powerSumOrder,
                                       powerSumIsMultiplicative);
        if (error()) return zero();
        ring_elem pIAtInnerSchur = toBasis(pIAtInner,
                                           powerSumBasisId0,
                                           powerSumDisplay,
                                           powerSumOrder,
                                           powerSumIsMultiplicative,
                                           schurId,
                                           schurDisplay,
                                           schurOrder,
                                           false);
        if (error()) return zero();
        ring_elem rest = schurCompletePlethysm(n - i,
                                               inner,
                                               powerSumBasisId0,
                                               powerSumDisplay,
                                               powerSumOrder,
                                               powerSumIsMultiplicative,
                                               schurId,
                                               schurDisplay,
                                               schurOrder);
        if (error()) return zero();
        ring_elem product =
            multiplySchurElements(pIAtInnerSchur, rest, schurId, schurDisplay, schurOrder);
        if (error()) return zero();
        total = add(total, product);
      }

    ring_elem reciprocal = rationalCoefficient(1, n);
    if (error()) return zero();
    ring_elem result = scaled(reciprocal, total);
    schurCompletePlethysmCache[key] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::schurPlethysmJacobiTrudi(const Partition& outer,
                                     const Partition& inner,
                                     int powerSumBasisId0,
                                     const std::string& powerSumDisplay,
                                     int powerSumOrder,
                                     bool powerSumIsMultiplicative,
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

    RingElemMatrix matrix(n, RingElemVector(n));
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        {
          int degree = outer[i] - static_cast<int>(i) + static_cast<int>(j);
          matrix[i][j] = schurCompletePlethysm(degree,
                                               inner,
                                               powerSumBasisId0,
                                               powerSumDisplay,
                                               powerSumOrder,
                                               powerSumIsMultiplicative,
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
            ring_elem term = multiplySchurElements(dp[mask],
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

bool SymmetricEngineRing::schurPlethysmToSchur(ring_elem f,
                            ring_elem g,
                            int powerSumBasisId0,
                            const std::string& powerSumDisplay,
                            int powerSumOrder,
                            bool powerSumIsMultiplicative,
                            int targetBasisId,
                            const std::string& targetDisplay,
                            int targetOrder,
                            ring_elem& result) const
{
    if (targetDisplay != "S") return false;
    Partition outer;
    Partition inner;
    if (!singleSchurPartition(f, targetBasisId, outer)) return false;
    if (!singleSchurPartition(g, targetBasisId, inner)) return false;
    result = schurPlethysmJacobiTrudi(outer,
                                      inner,
                                      powerSumBasisId0,
                                      powerSumDisplay,
                                      powerSumOrder,
                                      powerSumIsMultiplicative,
                                      targetBasisId,
                                      targetDisplay,
                                      targetOrder);
    return !error();
  }

ring_elem SymmetricEngineRing::plethysm(ring_elem f,
                     ring_elem g,
                     int pBasisId,
                     const std::string& pDisplay,
                     int pOrder,
                     bool pIsMultiplicative) const
{
    rememberBasis(pBasisId, pDisplay, pOrder, pIsMultiplicative);
    ring_elem fPowerSums = toBasis(f,
                                   pBasisId,
                                   pDisplay,
                                   pOrder,
                                   pIsMultiplicative,
                                   pBasisId,
                                   pDisplay,
                                   pOrder,
                                   pIsMultiplicative);
    if (error()) return zero();
    ring_elem gPowerSums = toBasis(g,
                                   pBasisId,
                                   pDisplay,
                                   pOrder,
                                   pIsMultiplicative,
                                   pBasisId,
                                   pDisplay,
                                   pOrder,
                                   pIsMultiplicative);
    if (error()) return zero();
    return plethysmPowerSums(fPowerSums, gPowerSums);
  }

ring_elem SymmetricEngineRing::plethysmToBasis(ring_elem f,
                            ring_elem g,
                            int pBasisId,
                            const std::string& pDisplay,
                            int pOrder,
                            bool pIsMultiplicative,
                            int targetBasisId,
                            const std::string& targetDisplay,
                            int targetOrder,
                            bool targetIsMultiplicative) const
{
    rememberBasis(pBasisId, pDisplay, pOrder, pIsMultiplicative);
    rememberBasis(targetBasisId, targetDisplay, targetOrder, targetIsMultiplicative);

    ring_elem specialized;
    if (schurPlethysmToSchur(f,
                             g,
                             pBasisId,
                             pDisplay,
                             pOrder,
                             pIsMultiplicative,
                             targetBasisId,
                             targetDisplay,
                             targetOrder,
                             specialized))
      return specialized;
    if (error()) return zero();

    ring_elem result = plethysm(f, g, pBasisId, pDisplay, pOrder, pIsMultiplicative);
    if (error()) return zero();
    return toBasis(result,
                   pBasisId,
                   pDisplay,
                   pOrder,
                   pIsMultiplicative,
                   targetBasisId,
                   targetDisplay,
                   targetOrder,
                   targetIsMultiplicative);
  }

} // namespace symmetric_rings
