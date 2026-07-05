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
std::string SymmetricEngineRing::jacobiTrudiCacheKey(int basisId,
                                  const Partition& outer,
                                  const Partition& inner) const
{
    return std::to_string(basisId) + "|" + partitionKey(outer) + "/" +
           partitionKey(inner);
  }

int SymmetricEngineRing::popcountMask(size_t mask) const
{
    int result = 0;
    while (mask != 0)
      {
        result += static_cast<int>(mask & 1);
        mask >>= 1;
      }
    return result;
  }

int SymmetricEngineRing::selectedGreaterThan(size_t mask, size_t col, size_t n) const
{
    int result = 0;
    for (size_t j = col + 1; j < n; ++j)
      if ((mask & (static_cast<size_t>(1) << j)) != 0) ++result;
    return result;
  }

ring_elem SymmetricEngineRing::jacobiTrudi(const Partition& outer,
                        const Partition& inner,
                        int basisId,
                        const std::string& display,
                        int order,
                        bool isMultiplicative) const
{
    rememberBasis(basisId, display, order, isMultiplicative);
    auto& cache = display == "e" ? eJacobiTrudiCache : hJacobiTrudiCache;
    std::string cacheKey = jacobiTrudiCacheKey(basisId, outer, inner);
    auto cached = cache.find(cacheKey);
    if (cached != cache.end()) return copyPolyValue(polyValue(cached->second));

    size_t n = std::max(outer.size(), inner.size());
    if (n == 0) return one();
    if (n >= 8 * sizeof(size_t))
      {
        ERROR("Jacobi-Trudi determinant is too large");
        return zero();
      }

    RingElemMatrix matrix(n, RingElemVector(n));
    for (size_t i = 0; i < n; ++i)
      {
        int lambdaI = i < outer.size() ? outer[i] : 0;
        for (size_t j = 0; j < n; ++j)
          {
            int muJ = j < inner.size() ? inner[j] : 0;
            int degree = lambdaI - muJ - static_cast<int>(i) + static_cast<int>(j);
            matrix[i][j] = basisPartElement(basisId,
                                            display,
                                            order,
                                            isMultiplicative,
                                            degree);
          }
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
            ring_elem term = mult(dp[mask], matrix[row][col]);
            if (selectedGreaterThan(mask, col, n) % 2 == 1) term = negate(term);
            size_t next = mask | bit;
            dp[next] = add(dp[next], term);
          }
      }

    ring_elem result = dp[limit - 1];
    cache[cacheKey] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::scaled(ring_elem coeff, ring_elem f) const
{
    return makePolyValue(multByCoefficient(coeff, polyValue(f)));
  }

ring_elem SymmetricEngineRing::coefficientQuotient(ring_elem numerator, ring_elem denominator) const
{
    ring_elem quotient = coefficientRing->divide(numerator, denominator);
    ring_elem check = coefficientRing->mult(quotient, denominator);
    if (coefficientRing->is_equal(check, numerator)) return quotient;
    ERROR("coefficient division failed during basis conversion; use a coefficient ring where the required denominators are invertible, for example frac(QQ[t]) instead of QQ[t]");
    return coefficientRing->zero();
  }

ring_elem SymmetricEngineRing::hallLittlewoodCFactor(const Partition& lambda) const
{
    std::map<int, int> multiplicities;
    for (int part : normalizePartition(lambda)) multiplicities[part]++;
    ring_elem result = coefficientRing->one();
    for (const auto& item : multiplicities)
      for (int j = 1; j <= item.second; ++j)
        {
          ring_elem tPower = coefficientRing->power(hallLittlewoodParameter, j);
          result = coefficientRing->mult(
              result,
              coefficientRing->subtract(coefficientRing->one(), tPower));
        }
    return result;
  }

void SymmetricEngineRing::addCoeff(CoeffMap& target, const Partition& index, ring_elem coeff) const
{
    if (coefficientRing->is_zero(coeff)) return;
    Partition key = normalizePartition(index);
    auto existing = target.find(key);
    if (existing == target.end())
      {
        target[key] = coeff;
        return;
      }
    ring_elem sum = coefficientRing->add(existing->second, coeff);
    if (coefficientRing->is_zero(sum))
      target.erase(existing);
    else
      existing->second = sum;
  }

CoeffMap SymmetricEngineRing::scaledCoeffMap(ring_elem coeff, const CoeffMap& source) const
{
    CoeffMap result;
    if (coefficientRing->is_zero(coeff)) return result;
    for (const auto& item : source)
      addCoeff(result, item.first, coefficientRing->mult(coeff, item.second));
    return result;
  }

CoeffMap SymmetricEngineRing::addCoeffMaps(const CoeffMap& a, const CoeffMap& b) const
{
    CoeffMap result = a;
    for (const auto& item : b) addCoeff(result, item.first, item.second);
    return result;
  }

CoeffMap SymmetricEngineRing::multiplyCoeffMaps(const CoeffMap& a, const CoeffMap& b) const
{
    CoeffMap result;
    for (const auto& left : a)
      for (const auto& right : b)
        {
          Partition index = left.first;
          index.insert(index.end(), right.first.begin(), right.first.end());
          ring_elem coeff = coefficientRing->mult(left.second, right.second);
          addCoeff(result, index, coeff);
        }
    return result;
  }

CoeffMap SymmetricEngineRing::multiplyMonomialCoeffMaps(const CoeffMap& a,
                                                        const CoeffMap& b) const
{
    CoeffMap result;
    for (const auto& left : a)
      for (const auto& right : b)
        {
          ring_elem baseCoeff = coefficientRing->mult(left.second, right.second);
          if (coefficientRing->is_zero(baseCoeff)) continue;
          for (const auto& product : monomialProduct(left.first, right.first))
            {
              ring_elem coeff = product.coefficient == 1
                  ? baseCoeff
                  : coefficientRing->mult(cachedInteger(product.coefficient),
                                          baseCoeff);
              addCoeff(result, product.nu, coeff);
            }
        }
    return result;
  }

CoeffMap SymmetricEngineRing::oneCoeffMap() const
{
    return CoeffMap{{Partition{}, coefficientRing->one()}};
  }

bool SymmetricEngineRing::isPartitionIndex(const Partition& p) const
{
    return trimTrailingZerosPartition(p) == normalizePartition(p);
  }

int SymmetricEngineRing::partitionPart(const Partition& p, size_t i) const
{
    return i < p.size() ? p[i] : 0;
  }

bool SymmetricEngineRing::partitionContains(const Partition& outer, const Partition& inner) const
{
    for (size_t i = 0; i < inner.size(); ++i)
      if (partitionPart(outer, i) < inner[i]) return false;
    return true;
  }

std::string SymmetricEngineRing::lrProductKey(const Partition& lambda, const Partition& mu) const
{
    if (lexLessPartition(mu, lambda))
      return partitionKey(mu) + "*" + partitionKey(lambda);
    return partitionKey(lambda) + "*" + partitionKey(mu);
  }

ring_elem SymmetricEngineRing::cachedInteger(long n) const
{
    auto found = smallIntegerCoeffCache.find(n);
    if (found != smallIntegerCoeffCache.end()) return found->second;
    ring_elem value = coefficientRing->from_long(n);
    smallIntegerCoeffCache[n] = value;
    return value;
  }

long SymmetricEngineRing::lrCoefficient(const Partition& lambda,
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

const std::vector<LRProductTerm>& SymmetricEngineRing::lrProduct(const Partition& a,
                                              const Partition& b) const
{
    Partition lambda = normalizePartition(a);
    Partition mu = normalizePartition(b);
    std::string key = lrProductKey(lambda, mu);
    auto cached = lrProductCache.find(key);
    if (cached != lrProductCache.end()) return cached->second;

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
            long c = lrCoefficient(lambda, mu, nu);
            if (c != 0) result.push_back({nu, c});
          }
      }

    auto inserted = lrProductCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

const std::vector<LRProductTerm>& SymmetricEngineRing::directLRProduct(
    const Partition& a,
    const Partition& b) const
{
    Partition first = normalizePartition(a);
    Partition second = normalizePartition(b);
    if (lexLessPartition(second, first)) std::swap(first, second);

    std::pair<Partition, Partition> key{first, second};
    auto cached = directLRProductCache.find(key);
    if (cached != directLRProductCache.end()) return cached->second;

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

    auto inserted = directLRProductCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

const std::vector<LRProductTerm>& SymmetricEngineRing::skewSchurExpansion(
    const Partition& outer0,
    const Partition& inner0) const
{
    Partition outer = normalizePartition(outer0);
    Partition inner = normalizePartition(inner0);
    std::pair<Partition, Partition> key{outer, inner};
    auto cached = skewSchurExpansionCache.find(key);
    if (cached != skewSchurExpansionCache.end()) return cached->second;

    std::vector<LRProductTerm> result;
    if (!partitionContains(outer, inner))
      {
        auto inserted = skewSchurExpansionCache.emplace(key, std::move(result));
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
            long coefficient = lrCoefficient(inner, alpha, outer);
            if (coefficient != 0) result.push_back({alpha, coefficient});
          }
      }

    auto inserted = skewSchurExpansionCache.emplace(key, std::move(result));
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

const std::vector<LRProductTerm>& SymmetricEngineRing::powerSumSchurProduct(
    const Partition& lambda0,
    int part) const
{
    Partition lambda = normalizePartition(lambda0);
    std::pair<Partition, int> key{lambda, part};
    auto cached = powerSumSchurProductCache.find(key);
    if (cached != powerSumSchurProductCache.end()) return cached->second;

    std::vector<LRProductTerm> result;
    if (part < 0)
      {
        auto inserted = powerSumSchurProductCache.emplace(key, std::move(result));
        return inserted.first->second;
      }
    if (part == 0)
      {
        result.push_back({lambda, 1});
        auto inserted = powerSumSchurProductCache.emplace(key, std::move(result));
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

    auto inserted = powerSumSchurProductCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

long SymmetricEngineRing::monomialProductCoefficient(
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

const std::vector<LRProductTerm>& SymmetricEngineRing::monomialProduct(
    const Partition& a,
    const Partition& b) const
{
    Partition first = normalizePartition(a);
    Partition second = normalizePartition(b);
    if (lexLessPartition(second, first)) std::swap(first, second);
    std::pair<Partition, Partition> key{first, second};
    auto cached = monomialProductCache.find(key);
    if (cached != monomialProductCache.end()) return cached->second;

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
            long coefficient = monomialProductCoefficient(first, second, nu);
            if (coefficient != 0) result.push_back({nu, coefficient});
          }
      }

    auto inserted = monomialProductCache.emplace(key, std::move(result));
    return inserted.first->second;
  }

std::vector<Partition> SymmetricEngineRing::horizontalStripProducts(
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

std::vector<Partition> SymmetricEngineRing::verticalStripProducts(
    const Partition& lambda,
    int col) const
{
    std::vector<Partition> result;
    if (col < 0) return result;
    Partition conjugate = conjugatePartition(lambda);
    for (const auto& nu : horizontalStripProducts(conjugate, col))
      result.push_back(conjugatePartition(nu));
    return result;
  }

ring_elem SymmetricEngineRing::multiplySchurElements(ring_elem f,
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
          for (const auto& product : lrProduct(a.first, b.first))
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

std::string SymmetricEngineRing::powerSumsToSchurRecipeKey(int degree,
                                        const std::vector<Partition>& inputPartitions,
                                        bool omegaStyle) const
{
    std::ostringstream out;
    out << degree << "|" << (omegaStyle ? 1 : 0);
    for (const auto& mu : inputPartitions)
      out << ";" << partitionKey(mu);
    return out.str();
  }

const std::vector<SchurConversionRecipeEntry>&
  SymmetricEngineRing::powerSumsToSchurRecipe(int degree,
                         const std::vector<Partition>& inputPartitions,
                         bool omegaStyle) const
{
    std::string key =
        powerSumsToSchurRecipeKey(degree, inputPartitions, omegaStyle);
    auto cached = powerSumsToSchurRecipeCache.find(key);
    if (cached != powerSumsToSchurRecipeCache.end()) return cached->second;

    const CharacterTable& table = characterTable(degree);
    std::vector<SchurConversionRecipeEntry> recipe;
    for (size_t row = 0; row < table.partitions.size(); ++row)
      {
        SchurConversionRecipeEntry entry;
        entry.lambda = table.partitions[row];
        for (size_t input = 0; input < inputPartitions.size(); ++input)
          {
            const Partition& mu = inputPartitions[input];
            auto muCol = table.partitionRows.find(mu);
            int chi = muCol == table.partitionRows.end()
                ? characterValue(table.partitions[row], mu)
                : characterTableValue(table, row, muCol->second);
            if (chi == 0) continue;
            if (omegaStyle && ((degree - partitionLength(mu)) % 2 != 0))
              chi = -chi;
            entry.contributions.push_back({input, chi});
          }
        if (!entry.contributions.empty()) recipe.push_back(std::move(entry));
      }

    auto inserted = powerSumsToSchurRecipeCache.emplace(key, std::move(recipe));
    return inserted.first->second;
  }
} // namespace symmetric_rings
