// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"
#include "rings/ZZ.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace symmetric_rings {

long SymmetricEngineRing::littlewoodRichardsonCoefficientViaTableaux(const Partition& lambda,
                     const Partition& content,
                     const Partition& nu) const
{
    if (!partitionContains(nu, lambda)) return 0;
    if (partitionWeight(nu) != partitionWeight(lambda) + partitionWeight(content))
      return 0;
    if (content.empty()) return trimTrailingZerosPartition(nu) == lambda ? 1 : 0;

    std::vector<std::pair<int, int>> cells;
    for (size_t r = 0; r < nu.size(); ++r)
      {
        int inner = partitionPart(lambda, r);
        for (int c = nu[r]; c > inner; --c)
          cells.push_back({static_cast<int>(r), c});
      }

    std::vector<std::vector<int>> tableau(nu.size());
    for (size_t r = 0; r < nu.size(); ++r)
      tableau[r].assign(nu[r] + 2, 0);

    int alphabet = static_cast<int>(content.size());
    std::vector<int> remaining(content.begin(), content.end());
    std::vector<int> used(alphabet, 0);
    long count = 0;

    std::function<void(size_t)> fill = [&](size_t k) {
      if (k == cells.size())
        {
          ++count;
          return;
        }

      int row = cells[k].first;
      int col = cells[k].second;
      int right = (col + 1 <= partitionPart(nu, row) &&
                   col + 1 > partitionPart(lambda, row))
          ? tableau[row][col + 1]
          : 0;
      int above = (row > 0 &&
                   col <= partitionPart(nu, static_cast<size_t>(row - 1)) &&
                   col > partitionPart(lambda, static_cast<size_t>(row - 1)))
          ? tableau[row - 1][col]
          : 0;

      for (int value = 1; value <= alphabet; ++value)
        {
          int idx = value - 1;
          if (remaining[idx] == 0) continue;
          if (right != 0 && value > right) continue;
          if (above != 0 && above >= value) continue;

          --remaining[idx];
          ++used[idx];
          bool lattice = true;
          for (int i = 0; i + 1 < alphabet; ++i)
            if (used[i] < used[i + 1])
              {
                lattice = false;
                break;
              }
          if (lattice)
            {
              tableau[row][col] = value;
              fill(k + 1);
              tableau[row][col] = 0;
            }
          --used[idx];
          ++remaining[idx];
        }
    };

    fill(0);
    return count;
  }

void SymmetricEngineRing::partitionsContainingRec(const Partition& lambda,
                               int addedWeight,
                               size_t row,
                               int previousPart,
                               Partition& current,
                               std::vector<Partition>& result) const
{
    if (addedWeight == 0 && row >= lambda.size())
      {
        result.push_back(trimTrailingZerosPartition(current));
        return;
      }

    int base = partitionPart(lambda, row);
    if (base == 0 && addedWeight == 0)
      {
        result.push_back(trimTrailingZerosPartition(current));
        return;
      }
    if (base == 0 && previousPart == 0 && addedWeight > 0) return;

    int upper = std::min(previousPart, base + addedWeight);
    for (int part = upper; part >= base; --part)
      {
        int extra = part - base;
        if (extra > addedWeight) continue;
        current.push_back(part);
        partitionsContainingRec(lambda,
                                addedWeight - extra,
                                row + 1,
                                part,
                                current,
                                result);
        current.pop_back();
      }
  }

std::vector<Partition> SymmetricEngineRing::partitionsContaining(const Partition& lambda,
                                              int addedWeight) const
{
    std::vector<Partition> result;
    Partition current;
    int firstBound = lambda.empty() ? addedWeight : lambda.front() + addedWeight;
    partitionsContainingRec(lambda, addedWeight, 0, firstBound, current, result);
    return result;
  }

const std::vector<LRProductTerm>& SymmetricEngineRing::littlewoodRichardsonProductViaCoefficientEnumeration(const Partition& a,
                                              const Partition& b) const
{
    Partition lambda = normalizePartition(a);
    Partition mu = normalizePartition(b);
    std::string key = littlewoodRichardsonProductCacheKey(lambda, mu);
    auto cached = littlewoodRichardsonCoefficientProductCache.find(key);
    if (cached != littlewoodRichardsonCoefficientProductCache.end()) return cached->second;

    std::vector<LRProductTerm> result;
    if (lambda.empty())
      {
        result.push_back({mu, 1});
      }
    else if (mu.empty())
      {
        result.push_back({lambda, 1});
      }
    else
      {
        if (partitionWeight(lambda) < partitionWeight(mu)) std::swap(lambda, mu);
        for (const auto& nu : partitionsContaining(lambda, partitionWeight(mu)))
          {
            long c = littlewoodRichardsonCoefficientViaTableaux(lambda, mu, nu);
            if (c != 0) result.push_back({nu, c});
          }
      }

    auto inserted = littlewoodRichardsonCoefficientProductCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

const std::vector<LRProductTerm>& SymmetricEngineRing::littlewoodRichardsonProductViaTableauEnumeration(
    const Partition& a,
    const Partition& b) const
{
    Partition first = normalizePartition(a);
    Partition second = normalizePartition(b);
    if (lexLessPartition(second, first)) std::swap(first, second);

    std::pair<Partition, Partition> key{first, second};
    auto cached = littlewoodRichardsonTableauProductCache.find(key);
    if (cached != littlewoodRichardsonTableauProductCache.end()) return cached->second;

    std::vector<LRProductTerm> result;
    if (first.empty())
      {
        result.push_back({second, 1});
      }
    else if (second.empty())
      {
        result.push_back({first, 1});
      }
    else
      {
        int firstRows = static_cast<int>(first.size());
        int secondRows = static_cast<int>(second.size());
        int maxRows = firstRows + secondRows;
        int offset = second.front();

        std::vector<int> outer(maxRows + 1, 0);
        std::vector<int> inner(maxRows + 1, 0);
        for (int i = 1; i <= firstRows; ++i)
          {
            outer[i] = offset + first[static_cast<size_t>(i - 1)];
            inner[i] = offset;
          }
        for (int i = firstRows + 1; i <= maxRows; ++i)
          outer[i] = second[static_cast<size_t>(i - firstRows - 1)];

        int finalWeight = 0;
        for (int i = 1; i <= maxRows; ++i) finalWeight += outer[i] - inner[i];

        std::vector<int> cellRow(finalWeight + 1, 0);
        std::vector<int> cellCol(finalWeight + 1, 0);
        int nextCell = 1;
        for (int row = 1; row <= maxRows; ++row)
          for (int col = outer[row]; col > inner[row]; --col)
            {
              cellRow[nextCell] = row;
              cellCol[nextCell] = col;
              ++nextCell;
            }

        std::vector<int> counts(maxRows + 1, 0);
        std::vector<int> valueAtCell(finalWeight + 1, 0);
        std::vector<int> countAtCell(finalWeight + 1, 0);
        std::vector<Partition> emitted;

        std::function<void(int)> fill = [&](int current) {
          if (current == finalWeight)
            {
              Partition nu;
              nu.reserve(static_cast<size_t>(maxRows));
              for (int i = 1; i <= maxRows; ++i)
                if (counts[i] > 0) nu.push_back(counts[i]);
              emitted.push_back(nu);
              return;
            }

          int k = current + 1;
          int row = cellRow[k];
          int col = cellCol[k];

          int hi;
          if (col == outer[row])
            {
              hi = maxRows;
              for (int i = 1; i <= maxRows; ++i)
                if (counts[i] == 0)
                  {
                    hi = i;
                    break;
                  }
            }
          else
            hi = valueAtCell[k - 1];

          int lo = 1;
          if (row > 1 && col > inner[row - 1])
            {
              int above = k - outer[row] + inner[row - 1];
              int aboveValue = valueAtCell[above];
              int aboveCount = countAtCell[above];
              lo = hi + 1;
              for (int i = aboveValue + 1; i <= hi; ++i)
                if (counts[i] < aboveCount)
                  {
                    lo = i;
                    break;
                  }
            }

          int thisCount = 32000;
          for (int value = lo; value <= hi; ++value)
            {
              int previousCount = thisCount;
              thisCount = counts[value];
              if (previousCount <= thisCount) continue;
              ++counts[value];
              valueAtCell[k] = value;
              countAtCell[k] = counts[value];
              fill(k);
              --counts[value];
            }
        };

        fill(0);
        std::sort(emitted.begin(), emitted.end());
        for (size_t i = 0; i < emitted.size();)
          {
            size_t j = i + 1;
            while (j < emitted.size() && emitted[j] == emitted[i]) ++j;
            result.push_back({emitted[i], static_cast<long>(j - i)});
            i = j;
          }
      }

    auto inserted = littlewoodRichardsonTableauProductCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

const std::vector<LRProductTerm>& SymmetricEngineRing::skewSchurToSchurViaLittlewoodRichardson(
    const Partition& outer0,
    const Partition& inner0) const
{
    Partition outer = normalizePartition(outer0);
    Partition inner = normalizePartition(inner0);
    std::pair<Partition, Partition> key{outer, inner};
    auto cached = skewSchurToSchurViaLittlewoodRichardsonCache.find(key);
    if (cached != skewSchurToSchurViaLittlewoodRichardsonCache.end()) return cached->second;

    std::vector<LRProductTerm> result;
    if (!partitionContains(outer, inner))
      {
        auto inserted = skewSchurToSchurViaLittlewoodRichardsonCache.emplace(key, std::move(result));
        return inserted.first->second;
      }

    int weightDifference = partitionWeight(outer) - partitionWeight(inner);
    if (weightDifference == 0)
      {
        result.push_back({Partition{}, 1});
      }
    else if (inner.empty())
      {
        result.push_back({outer, 1});
      }
    else
      {
        for (const auto& alpha : partitionsOf(weightDifference))
          {
            long coefficient = littlewoodRichardsonCoefficientViaTableaux(inner, alpha, outer);
            if (coefficient != 0) result.push_back({alpha, coefficient});
          }
      }

    auto inserted = skewSchurToSchurViaLittlewoodRichardsonCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

bool SymmetricEngineRing::addedBorderStripCell(const Partition& lambda,
                                               const Partition& nu,
                                               int row,
                                               int col) const
{
    if (row < 0 || static_cast<size_t>(row) >= nu.size()) return false;
    if (col < 1 || col > nu[static_cast<size_t>(row)]) return false;
    return col > partitionPart(lambda, static_cast<size_t>(row));
  }

bool SymmetricEngineRing::addedBorderStripConnected(const Partition& lambda,
                                                    const Partition& nu) const
{
    int total = 0;
    int startRow = -1;
    int startCol = -1;
    for (size_t r = 0; r < nu.size(); ++r)
      for (int c = partitionPart(lambda, r) + 1; c <= nu[r]; ++c)
        {
          ++total;
          if (startRow < 0)
            {
              startRow = static_cast<int>(r);
              startCol = c;
            }
        }
    if (total == 0) return false;

    std::vector<std::pair<int, int>> stack{{startRow, startCol}};
    std::map<std::pair<int, int>, bool> seen;
    seen[stack.back()] = true;
    int visited = 0;
    while (!stack.empty())
      {
        auto cell = stack.back();
        stack.pop_back();
        ++visited;
        const int dr[4] = {1, -1, 0, 0};
        const int dc[4] = {0, 0, 1, -1};
        for (int i = 0; i < 4; ++i)
          {
            std::pair<int, int> next{cell.first + dr[i], cell.second + dc[i]};
            if (!seen[next] &&
                addedBorderStripCell(lambda, nu, next.first, next.second))
              {
                seen[next] = true;
                stack.push_back(next);
              }
          }
      }
    return visited == total;
  }

bool SymmetricEngineRing::addedBorderStripHasNoTwoByTwo(
    const Partition& lambda,
    const Partition& nu) const
{
    for (size_t r = 0; r + 1 < nu.size(); ++r)
      for (int c = 1; c <= std::max(partitionPart(nu, r),
                                    partitionPart(nu, r + 1)); ++c)
        if (addedBorderStripCell(lambda, nu, static_cast<int>(r), c) &&
            addedBorderStripCell(lambda, nu, static_cast<int>(r) + 1, c) &&
            addedBorderStripCell(lambda, nu, static_cast<int>(r), c + 1) &&
            addedBorderStripCell(lambda, nu, static_cast<int>(r) + 1, c + 1))
          return false;
    return true;
  }

const std::vector<LRProductTerm>& SymmetricEngineRing::schurTimesPowerSumViaBorderStrips(
    const Partition& lambda0,
    int part) const
{
    Partition lambda = normalizePartition(lambda0);
    std::pair<Partition, int> key{lambda, part};
    auto cached = schurTimesPowerSumViaBorderStripsCache.find(key);
    if (cached != schurTimesPowerSumViaBorderStripsCache.end()) return cached->second;

    std::vector<LRProductTerm> result;
    if (part < 0)
      {
        auto inserted = schurTimesPowerSumViaBorderStripsCache.emplace(key, std::move(result));
        return inserted.first->second;
      }
    if (part == 0)
      {
        result.push_back({lambda, 1});
        auto inserted = schurTimesPowerSumViaBorderStripsCache.emplace(key, std::move(result));
        return inserted.first->second;
      }

    for (const auto& nu : partitionsContaining(lambda, part))
      {
        if (!addedBorderStripConnected(lambda, nu)) continue;
        if (!addedBorderStripHasNoTwoByTwo(lambda, nu)) continue;
        int rows = 0;
        for (size_t r = 0; r < nu.size(); ++r)
          if (partitionPart(nu, r) > partitionPart(lambda, r)) ++rows;
        long sign = ((rows - 1) % 2 == 0) ? 1 : -1;
        result.push_back({nu, sign});
      }

    auto inserted = schurTimesPowerSumViaBorderStripsCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

long SymmetricEngineRing::monomialProductCoefficientViaExponentSplittings(
    const Partition& lambda0,
    const Partition& mu0,
    const Partition& nu0) const
{
    Partition lambda = normalizePartition(lambda0);
    Partition mu = normalizePartition(mu0);
    Partition nu = normalizePartition(nu0);
    if (partitionWeight(lambda) + partitionWeight(mu) != partitionWeight(nu))
      return 0;

    std::map<int, int> lambdaCounts;
    std::map<int, int> muCounts;
    for (int part : lambda) ++lambdaCounts[part];
    for (int part : mu) ++muCounts[part];

    auto allUsed = [](const std::map<int, int>& counts) {
      for (const auto& item : counts)
        if (item.second != 0) return false;
      return true;
    };

    long total = 0;
    std::function<void(size_t)> split = [&](size_t pos) {
      if (pos == nu.size())
        {
          if (allUsed(lambdaCounts) && allUsed(muCounts)) ++total;
          return;
        }

      int part = nu[pos];
      for (int left = 0; left <= part; ++left)
        {
          int right = part - left;
          if (left > 0)
            {
              auto found = lambdaCounts.find(left);
              if (found == lambdaCounts.end() || found->second == 0) continue;
              --found->second;
            }
          if (right > 0)
            {
              auto found = muCounts.find(right);
              if (found == muCounts.end() || found->second == 0)
                {
                  if (left > 0) ++lambdaCounts[left];
                  continue;
                }
              --found->second;
            }

          split(pos + 1);

          if (right > 0) ++muCounts[right];
          if (left > 0) ++lambdaCounts[left];
        }
    };

    split(0);
    return total;
  }

const std::vector<LRProductTerm>& SymmetricEngineRing::monomialProductViaExponentSplittings(
    const Partition& a,
    const Partition& b) const
{
    Partition first = normalizePartition(a);
    Partition second = normalizePartition(b);
    if (lexLessPartition(second, first)) std::swap(first, second);
    std::pair<Partition, Partition> key{first, second};
    auto cached = monomialProductViaExponentSplittingsCache.find(key);
    if (cached != monomialProductViaExponentSplittingsCache.end()) return cached->second;

    std::vector<LRProductTerm> result;
    if (first.empty())
      {
        result.push_back({second, 1});
      }
    else if (second.empty())
      {
        result.push_back({first, 1});
      }
    else
      {
        int totalWeight = partitionWeight(first) + partitionWeight(second);
        for (const auto& nu : partitionsOf(totalWeight))
          {
            if (partitionLength(nu) > partitionLength(first) + partitionLength(second))
              continue;
            long coefficient = monomialProductCoefficientViaExponentSplittings(first, second, nu);
            if (coefficient != 0) result.push_back({nu, coefficient});
          }
      }

    auto inserted = monomialProductViaExponentSplittingsCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

std::vector<Partition> SymmetricEngineRing::schurTimesCompleteViaHorizontalPieri(
    const Partition& lambda,
    int row) const
{
    std::vector<Partition> result;
    if (row < 0) return result;
    if (row == 0)
      {
        result.push_back(normalizePartition(lambda));
        return result;
      }
    Partition base = normalizePartition(lambda);
    for (const auto& nu : partitionsContaining(base, row))
      {
        bool horizontal = true;
        for (size_t i = 0; i + 1 < nu.size(); ++i)
          if (partitionPart(nu, i + 1) > partitionPart(base, i))
            {
              horizontal = false;
              break;
            }
        if (horizontal) result.push_back(nu);
      }
    return result;
  }

std::vector<Partition> SymmetricEngineRing::schurTimesElementaryViaVerticalPieri(
    const Partition& lambda,
    int col) const
{
    std::vector<Partition> result;
    if (col < 0) return result;
    Partition conjugate = conjugatePartition(lambda);
    for (const auto& nu : schurTimesCompleteViaHorizontalPieri(conjugate, col))
      result.push_back(conjugatePartition(nu));
    return result;
  }

ring_elem SymmetricEngineRing::multiplySchurExpansionsViaLittlewoodRichardson(ring_elem f,
                                  ring_elem g,
                                  int schurId,
                                  const std::string& schurDisplay,
                                  int schurOrder) const
{
    CoeffMap left = coefficientsInBasis(f, schurId);
    if (error()) return zero();
    CoeffMap right = coefficientsInBasis(g, schurId);
    if (error()) return zero();

    CoeffMap result;
    for (const auto& a : left)
      for (const auto& b : right)
        {
          ring_elem baseCoeff = coefficientRing->mult(a.second, b.second);
          if (coefficientRing->is_zero(baseCoeff)) continue;
          for (const auto& product : littlewoodRichardsonProductViaCoefficientEnumeration(a.first, b.first))
            {
              ring_elem coeff = product.coefficient == 1
                  ? baseCoeff
                  : coefficientRing->mult(coefficientRing->from_long(product.coefficient),
                                          baseCoeff);
              addCoeff(result, product.nu, coeff);
            }
        }
    return coeffMapToElement(result, schurId, schurDisplay, schurOrder, false);
  }

bool SymmetricEngineRing::tryAtomToSchurFactors(
    const SymmetricMonomial& monomial,
    size_t pos,
    int targetBasisId,
    std::vector<Partition>& factors) const
{
    if (atomIsSkewAt(monomial, pos)) return false;

    int basisId = atomBasisIdAt(monomial, pos);
    Partition index = atomIndex(monomial, pos);
    if (basisId == targetBasisId)
      {
        if (!isPartitionIndex(index)) return false;
        factors.push_back(trimTrailingZerosPartition(index));
        return true;
      }

    std::string display = displayForBasis(basisId);
    if (display != "h" && display != "e") return false;

    for (int part : index)
      {
        if (part < 0) return false;
        if (part == 0) continue;
        if (display == "h")
          factors.push_back(Partition{part});
        else
          factors.push_back(Partition(static_cast<size_t>(part), 1));
      }
    return true;
  }

bool SymmetricEngineRing::trySchurProductMonomialToSchurViaLittlewoodRichardson(
    const SymmetricMonomial& monomial,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder,
    ring_elem& result) const
{
    if (targetDisplay != "S") return false;
    if (monomial.data.empty())
      {
        result = one();
        return true;
      }

    CoeffMap current = oneCoeffMap();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        std::vector<Partition> factors;
        if (!tryAtomToSchurFactors(monomial, pos, targetBasisId, factors))
          return false;

        for (const auto& index : factors)
          {
            CoeffMap next;
            for (const auto& term : current)
              for (const auto& product : littlewoodRichardsonProductViaCoefficientEnumeration(term.first, index))
                {
                  ring_elem coeff = product.coefficient == 1
                      ? term.second
                      : coefficientRing->mult(coefficientRing->from_long(product.coefficient),
                                              term.second);
                  addCoeff(next, product.nu, coeff);
                }
            current = next;
          }
        pos += atomLengthAt(monomial, pos);
      }

    result = coeffMapToElement(current,
                               targetBasisId,
                               targetDisplay,
                               targetDisplayOrder,
                               false);
    return true;
  }

bool SymmetricEngineRing::trySchurCompatibleMonomialToSchur(
    const SymmetricMonomial& monomial,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder,
    ring_elem& result) const
{
    if (targetDisplay != "S") return false;
    std::vector<SchurCompatibleFactor> factors;
    if (!trySchurCompatibleFactorsFromMonomial(monomial, targetBasisId, factors))
      return false;
    return schurCompatibleFactorsToSchurDispatch(std::move(factors),
                                                targetBasisId,
                                                targetDisplay,
                                                targetDisplayOrder,
                                                result);
  }

bool SymmetricEngineRing::trySchurCompatibleFactorsFromMonomial(
    const SymmetricMonomial& monomial,
    int targetBasisId,
    std::vector<SchurCompatibleFactor>& factors) const
{
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        int basisId = atomBasisIdAt(monomial, pos);
        if (atomIsSkewAt(monomial, pos))
          {
            if (basisId != targetBasisId) return false;
            CoeffMap expansion;
            int weight = partitionWeight(atomOuterIndex(monomial, pos)) -
                         partitionWeight(atomInnerIndex(monomial, pos));
            for (const auto& item : skewSchurToSchurViaLittlewoodRichardson(atomOuterIndex(monomial, pos),
                                                       atomInnerIndex(monomial, pos)))
              addCoeff(expansion, item.nu, cachedInteger(item.coefficient));
            factors.push_back({SchurCompatibleFactor::SchurExpansion,
                               Partition{},
                               expansion,
                               weight});
            pos += atomLengthAt(monomial, pos);
            continue;
          }

        Partition index = atomIndex(monomial, pos);
        if (basisId == targetBasisId)
          {
            auto straightened = straightenSchurIndex(index);
            if (straightened.first == 1)
              {
                factors.push_back({SchurCompatibleFactor::General,
                                   straightened.second,
                                   CoeffMap{},
                                   partitionWeight(straightened.second)});
              }
            else
              {
                CoeffMap expansion;
                if (straightened.first != 0)
                  addCoeff(expansion,
                           straightened.second,
                           cachedInteger(straightened.first));
                factors.push_back({SchurCompatibleFactor::SchurExpansion,
                                   Partition{},
                                   expansion,
                                   partitionWeight(index)});
              }
          }
        else
          {
            std::string display = displayForBasis(basisId);
            if (display != "h" && display != "e" && display != "p") return false;
            for (int part : index)
              {
                if (part < 0) return false;
                if (part == 0) continue;
                SchurCompatibleFactor::Kind kind = SchurCompatibleFactor::PowerSum;
                if (display == "h")
                  kind = SchurCompatibleFactor::Horizontal;
                else if (display == "e")
                  kind = SchurCompatibleFactor::Vertical;
                factors.push_back({kind, Partition{part}, CoeffMap{}, part});
              }
          }
        pos += atomLengthAt(monomial, pos);
      }
    return true;
  }

SymmetricEngineRing::SchurFactorMethod
SymmetricEngineRing::selectSchurFactorMethod(
    const SchurCompatibleFactor& factor) const
{
    switch (factor.kind)
      {
        case SchurCompatibleFactor::Horizontal:
          return SchurFactorMethod::ViaHorizontalPieri;
        case SchurCompatibleFactor::Vertical:
          return SchurFactorMethod::ViaVerticalPieri;
        case SchurCompatibleFactor::PowerSum:
          return SchurFactorMethod::ViaBorderStrips;
        case SchurCompatibleFactor::SchurExpansion:
          return SchurFactorMethod::ViaLittlewoodRichardsonExpansion;
        case SchurCompatibleFactor::General:
          return SchurFactorMethod::ViaLittlewoodRichardson;
      }
    return SchurFactorMethod::ViaLittlewoodRichardson;
  }

const char *SymmetricEngineRing::schurFactorMethodName(
    SchurFactorMethod method) const
{
    switch (method)
      {
        case SchurFactorMethod::AlreadySchur:
          return "already-Schur";
        case SchurFactorMethod::ViaLittlewoodRichardson:
          return "littlewood-richardson";
        case SchurFactorMethod::ViaHorizontalPieri:
          return "horizontal-pieri";
        case SchurFactorMethod::ViaVerticalPieri:
          return "vertical-pieri";
        case SchurFactorMethod::ViaBorderStrips:
          return "border-strips";
        case SchurFactorMethod::ViaLittlewoodRichardsonExpansion:
          return "littlewood-richardson-expansion";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceSchurFactorMethod(
    SchurFactorMethod method,
    const SchurCompatibleFactor& factor) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::fprintf(stderr,
                 "SymmetricRings Schur factor: weight=%d method=%s\n",
                 factor.weight,
                 schurFactorMethodName(method));
  }

bool SymmetricEngineRing::schurCompatibleFactorsToSchurDispatch(
    std::vector<SchurCompatibleFactor> factors,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder,
    ring_elem& result) const
{
    if (targetDisplay != "S") return false;
    if (factors.empty())
      {
        result = one();
        return true;
      }

    if (factors.size() == 1 &&
        factors.front().kind == SchurCompatibleFactor::General)
      {
        traceSchurFactorMethod(SchurFactorMethod::AlreadySchur,
                               factors.front());
        result = basisElementFromIndex(targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       false,
                                       factors.front().index);
        return true;
      }

    if (factors.size() == 2 &&
        factors[0].kind == SchurCompatibleFactor::General &&
        factors[1].kind == SchurCompatibleFactor::General)
      {
        traceSchurFactorMethod(SchurFactorMethod::ViaLittlewoodRichardson,
                               factors[1]);
        const auto& product = littlewoodRichardsonProductViaTableauEnumeration(factors[0].index, factors[1].index);
        VECTOR(SymmetricTerm) terms;
        terms.reserve(product.size());
        rememberBasis(targetBasisId, targetDisplay, targetDisplayOrder, false);
        for (const auto& item : product)
          {
            SymmetricMonomial termMonomial;
            appendAtomBlock(termMonomial,
                            makeAtomBlock(targetDisplayOrder,
                                          targetBasisId,
                                          0,
                                          item.nu));
            terms.push_back({cachedInteger(item.coefficient),
                             canonicalMonomial(termMonomial)});
          }
        result = fromTermVector(terms, false);
        return true;
      }

    std::stable_sort(factors.begin(),
                     factors.end(),
                     [](const SchurCompatibleFactor& a,
                        const SchurCompatibleFactor& b) {
                       if (a.weight != b.weight) return a.weight < b.weight;
                       int aRank =
                           (a.kind == SchurCompatibleFactor::General ||
                            a.kind == SchurCompatibleFactor::SchurExpansion) ? 0 :
                           a.kind == SchurCompatibleFactor::PowerSum ? 2 : 1;
                       int bRank =
                           (b.kind == SchurCompatibleFactor::General ||
                            b.kind == SchurCompatibleFactor::SchurExpansion) ? 0 :
                           b.kind == SchurCompatibleFactor::PowerSum ? 2 : 1;
                       if (aRank != bRank) return aRank < bRank;
                       if (a.kind != b.kind) return a.kind < b.kind;
                       return lexLessPartition(a.index, b.index);
                     });

    CoeffMap current = oneCoeffMap();
    for (const auto& factor : factors)
      {
        SchurFactorMethod method = selectSchurFactorMethod(factor);
        traceSchurFactorMethod(method, factor);
        CoeffMap next;
        for (const auto& term : current)
          {
            if (method == SchurFactorMethod::ViaHorizontalPieri)
              {
                int row = factor.index.empty() ? 0 : factor.index.front();
                for (const auto& nu : schurTimesCompleteViaHorizontalPieri(term.first, row))
                  addCoeff(next, nu, term.second);
              }
            else if (method == SchurFactorMethod::ViaVerticalPieri)
              {
                int col = factor.index.empty() ? 0 : factor.index.front();
                for (const auto& nu : schurTimesElementaryViaVerticalPieri(term.first, col))
                  addCoeff(next, nu, term.second);
              }
            else if (method == SchurFactorMethod::ViaBorderStrips)
              {
                int part = factor.index.empty() ? 0 : factor.index.front();
                for (const auto& product : schurTimesPowerSumViaBorderStrips(term.first, part))
                  {
                    ring_elem coeff = product.coefficient == 1
                        ? term.second
                        : coefficientRing->mult(cachedInteger(product.coefficient),
                                                term.second);
                    addCoeff(next, product.nu, coeff);
                  }
              }
            else if (method == SchurFactorMethod::ViaLittlewoodRichardsonExpansion)
              {
                for (const auto& expansionTerm : factor.expansion)
                  {
                    ring_elem expansionCoeff = coefficientRing->mult(term.second,
                                                                     expansionTerm.second);
                    if (coefficientRing->is_zero(expansionCoeff)) continue;
                    for (const auto& product : littlewoodRichardsonProductViaTableauEnumeration(term.first,
                                                               expansionTerm.first))
                      {
                        ring_elem coeff = product.coefficient == 1
                            ? expansionCoeff
                            : coefficientRing->mult(cachedInteger(product.coefficient),
                                                    expansionCoeff);
                        addCoeff(next, product.nu, coeff);
                      }
                  }
              }
            else
              {
                for (const auto& product : littlewoodRichardsonProductViaTableauEnumeration(term.first, factor.index))
                  {
                    ring_elem coeff = product.coefficient == 1
                        ? term.second
                        : coefficientRing->mult(cachedInteger(product.coefficient),
                                                term.second);
                    addCoeff(next, product.nu, coeff);
                  }
              }
          }
        current = next;
      }

    result = coeffMapToElement(current,
                               targetBasisId,
                               targetDisplay,
                               targetDisplayOrder,
                               false);
    return true;
  }

bool SymmetricEngineRing::trySchurCompatibleExpressionToSchur(ring_elem f,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          ring_elem& result) const
{
    if (targetDisplay != "S") return false;
    const auto *poly = polyValue(f);
    VECTOR(SymmetricTerm) terms;
    for (const auto& term : poly->terms)
      {
        ring_elem converted;
        if (!trySchurCompatibleMonomialToSchur(term.monomial,
                                             targetBasisId,
                                             targetDisplay,
                                             targetDisplayOrder,
                                             converted))
          return false;
        const auto *convertedPoly = polyValue(converted);
        terms.reserve(terms.size() + convertedPoly->terms.size());
        for (const auto& convertedTerm : convertedPoly->terms)
          terms.push_back({coefficientRing->mult(term.coeff, convertedTerm.coeff),
                           convertedTerm.monomial});
      }
    result = fromTermVector(terms, false);
    return true;
  }

ring_elem SymmetricEngineRing::powerSumsToSchurViaBorderStrips(
                          ring_elem f,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder) const
{
    const auto *poly = polyValue(f);
    VECTOR(SymmetricTerm) terms;
    for (const auto& term : poly->terms)
      {
        std::vector<SchurCompatibleFactor> factors;
        if (!term.monomial.data.empty())
          {
            Partition index = atomIndex(term.monomial, 0);
            factors.reserve(index.size());
            for (int part : index)
              if (part > 0)
                factors.push_back({SchurCompatibleFactor::PowerSum,
                                   Partition{part},
                                   CoeffMap{},
                                   part});
          }

        ring_elem converted;
        if (!schurCompatibleFactorsToSchurDispatch(std::move(factors),
                                                    targetBasisId,
                                                    targetDisplay,
                                                    targetDisplayOrder,
                                                    converted))
          return zero();
        const auto *convertedPoly = polyValue(converted);
        terms.reserve(terms.size() + convertedPoly->terms.size());
        for (const auto& convertedTerm : convertedPoly->terms)
          terms.push_back({coefficientRing->mult(term.coeff,
                                                 convertedTerm.coeff),
                           convertedTerm.monomial});
      }
    return fromTermVector(terms, false);
  }

bool SymmetricEngineRing::tryProductToSchurViaCompatibleFactors(ring_elem f,
                          ring_elem g,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          ring_elem& result) const
{
    if (targetDisplay != "S") return false;
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    VECTOR(SymmetricTerm) terms;
    for (const auto& leftTerm : left->terms)
      for (const auto& rightTerm : right->terms)
        {
          ring_elem baseCoeff = coefficientRing->mult(leftTerm.coeff,
                                                      rightTerm.coeff);
          if (coefficientRing->is_zero(baseCoeff)) continue;

          SymmetricMonomial productMonomial =
              multiplyMonomials(leftTerm.monomial, rightTerm.monomial);
          ring_elem converted;
          if (!trySchurCompatibleMonomialToSchur(productMonomial,
                                               targetBasisId,
                                               targetDisplay,
                                               targetDisplayOrder,
                                               converted))
            return false;

          const auto *convertedPoly = polyValue(converted);
          terms.reserve(terms.size() + convertedPoly->terms.size());
          for (const auto& convertedTerm : convertedPoly->terms)
            terms.push_back({coefficientRing->mult(baseCoeff,
                                                   convertedTerm.coeff),
                             convertedTerm.monomial});
        }

    result = fromTermVector(terms, false);
    return true;
  }

bool SymmetricEngineRing::tryMonomialLikeAtomToCoeffMap(
    const SymmetricMonomial& monomial,
    size_t pos,
    const std::string& targetDisplay,
    CoeffMap& result) const
{
    if (atomIsSkewAt(monomial, pos)) return false;
    bool forgottenTarget = isForgottenDisplay(targetDisplay);
    if (targetDisplay != "m" && !forgottenTarget) return false;

    std::string display = displayForBasis(atomBasisIdAt(monomial, pos));
    Partition index = atomIndex(monomial, pos);

    auto partMap = [&](const std::string& kind, int part, long sign) {
      CoeffMap map;
      if (part < 0) return map;
      if (part == 0)
        {
          addCoeff(map, Partition{}, cachedInteger(sign));
          return map;
        }
      if (kind == "h")
        {
          for (const auto& lambda : partitionsOf(part))
            addCoeff(map, lambda, cachedInteger(sign));
        }
      else if (kind == "e")
        {
          addCoeff(map, Partition(static_cast<size_t>(part), 1),
                   cachedInteger(sign));
        }
      else if (kind == "p")
        {
          addCoeff(map, Partition{part}, cachedInteger(sign));
        }
      return map;
    };

    if (!forgottenTarget && display == "m")
      {
        addCoeff(result, index, coefficientRing->one());
        return true;
      }
    if (forgottenTarget && isForgottenDisplay(display))
      {
        addCoeff(result, index, coefficientRing->one());
        return true;
      }

    std::string mappedDisplay = display;
    bool powerSumOmegaSign = false;
    if (forgottenTarget)
      {
        if (display == "h")
          mappedDisplay = "e";
        else if (display == "e")
          mappedDisplay = "h";
        else if (display == "p")
          {
            mappedDisplay = "p";
            powerSumOmegaSign = true;
          }
        else
          return false;
      }
    else if (display != "h" && display != "e" && display != "p")
      return false;

    result = oneCoeffMap();
    for (int part : index)
      {
        if (part < 0) return false;
        long sign = 1;
        if (powerSumOmegaSign && part % 2 == 0) sign = -1;
        CoeffMap factor = partMap(mappedDisplay, part, sign);
        result = multiplyMonomialCoeffMaps(result, factor);
      }
    return true;
  }

bool SymmetricEngineRing::tryMonomialLikeMonomialToTarget(
    const SymmetricMonomial& monomial,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder,
    bool targetIsMultiplicative,
    ring_elem& result) const
{
    if (targetDisplay != "m" && !isForgottenDisplay(targetDisplay))
      return false;

    CoeffMap current = oneCoeffMap();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        CoeffMap factor;
        if (!tryMonomialLikeAtomToCoeffMap(monomial, pos, targetDisplay, factor))
          return false;
        current = multiplyMonomialCoeffMaps(current, factor);
        pos += atomLengthAt(monomial, pos);
      }

    result = coeffMapToElement(current,
                               targetBasisId,
                               targetDisplay,
                               targetDisplayOrder,
                               targetIsMultiplicative);
    return true;
  }

bool SymmetricEngineRing::tryProductToMonomialLikeTarget(
    ring_elem f,
    ring_elem g,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder,
    bool targetIsMultiplicative,
    ring_elem& result) const
{
    if (targetDisplay != "m" && !isForgottenDisplay(targetDisplay))
      return false;

    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    VECTOR(SymmetricTerm) terms;
    for (const auto& leftTerm : left->terms)
      for (const auto& rightTerm : right->terms)
        {
          ring_elem baseCoeff = coefficientRing->mult(leftTerm.coeff,
                                                      rightTerm.coeff);
          if (coefficientRing->is_zero(baseCoeff)) continue;

          SymmetricMonomial productMonomial =
              multiplyMonomials(leftTerm.monomial, rightTerm.monomial);
          ring_elem converted;
          if (!tryMonomialLikeMonomialToTarget(productMonomial,
                                            targetBasisId,
                                            targetDisplay,
                                            targetDisplayOrder,
                                            targetIsMultiplicative,
                                            converted))
            return false;

          const auto *convertedPoly = polyValue(converted);
          terms.reserve(terms.size() + convertedPoly->terms.size());
          for (const auto& convertedTerm : convertedPoly->terms)
            terms.push_back({coefficientRing->mult(baseCoeff,
                                                   convertedTerm.coeff),
                             convertedTerm.monomial});
        }

    result = fromTermVector(terms, false);
    return true;
  }

bool SymmetricEngineRing::tryProductToTargetDispatch(ring_elem f,
                          ring_elem g,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          bool targetIsMultiplicative,
                          ring_elem& result) const
{
    rememberBasis(targetBasisId,
                  targetDisplay,
                  targetDisplayOrder,
                  targetIsMultiplicative);

    if (targetDisplay == "S")
      return tryProductToSchurViaCompatibleFactors(f,
                                g,
                                targetBasisId,
                                targetDisplay,
                                targetDisplayOrder,
                                result);

    if (targetDisplay == "m" || isForgottenDisplay(targetDisplay))
      return tryProductToMonomialLikeTarget(f,
                                       g,
                                       targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       targetIsMultiplicative,
                                       result);

    int leftBasis = singleBasisId(f);
    int rightBasis = singleBasisId(g);
    if (targetIsMultiplicative && leftBasis == targetBasisId)
      {
        ring_elem convertedRight;
        if (tryExpressionToTarget(g,
                                  targetBasisId,
                                  targetDisplay,
                                  targetDisplayOrder,
                                  targetIsMultiplicative,
                                  convertedRight))
          {
            result = mult(f, convertedRight);
            return !error();
          }
        if (error()) return false;
      }
    if (targetIsMultiplicative && rightBasis == targetBasisId)
      {
        ring_elem convertedLeft;
        if (tryExpressionToTarget(f,
                                  targetBasisId,
                                  targetDisplay,
                                  targetDisplayOrder,
                                  targetIsMultiplicative,
                                  convertedLeft))
          {
            result = mult(convertedLeft, g);
            return !error();
          }
        if (error()) return false;
      }

    bool leftNative = leftBasis == 0 || leftBasis == targetBasisId;
    bool rightNative = rightBasis == 0 || rightBasis == targetBasisId;
    if (leftNative && rightNative)
      {
        result = mult(f, g);
        return !error();
      }

    return false;
  }

} // namespace symmetric_rings
