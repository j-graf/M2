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

Partition SymmetricEngineRing::leadingPartition(const CoeffMap& H) const
{
    if (H.empty()) return Partition{};
    Partition result = H.begin()->first;
    for (const auto& item : H)
      if (lexLessPartition(item.first, result)) result = item.first;
    return result;
  }

ring_elem SymmetricEngineRing::hPartToPowerSums(int n) const
{
    if (n == 0) return one();
    auto cached = hToPowerSumCache.find(n);
    if (cached != hToPowerSumCache.end()) return copyPolyValue(polyValue(cached->second));
    ring_elem result = zero();
    for (const auto& mu : partitionsOf(n))
      {
        ring_elem coeff = rationalCoefficient(1, zValue(mu));
        if (error()) return zero();
        ring_elem term = basisElementFromIndex(powerSumBasisId, "p", 10, true, mu);
        result = add(result, scaled(coeff, term));
      }
    hToPowerSumCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::ePartToPowerSums(int n) const
{
    if (n == 0) return one();
    auto cached = eToPowerSumCache.find(n);
    if (cached != eToPowerSumCache.end()) return copyPolyValue(polyValue(cached->second));
    ring_elem result = zero();
    for (const auto& mu : partitionsOf(n))
      {
        long sign = ((n - static_cast<int>(mu.size())) % 2 == 0) ? 1 : -1;
        ring_elem coeff = rationalCoefficient(sign, zValue(mu));
        if (error()) return zero();
        ring_elem term = basisElementFromIndex(powerSumBasisId, "p", 10, true, mu);
        result = add(result, scaled(coeff, term));
      }
    eToPowerSumCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::hallLittlewoodPartToPowerSums(int n, bool omega) const
{
    if (n == 0) return one();
    auto& cache = omega ? bToPowerSumCache : qToPowerSumCache;
    auto cached = cache.find(n);
    if (cached != cache.end()) return copyPolyValue(polyValue(cached->second));
    ring_elem result = zero();
    for (const auto& mu : partitionsOf(n))
      {
        long sign = 1;
        if (omega && ((n - static_cast<int>(mu.size())) % 2 == 1)) sign = -1;
        ring_elem rational = rationalCoefficient(sign, zValue(mu));
        if (error()) return zero();
        ring_elem coeff = coefficientRing->mult(rational, hallLittlewoodFactor(mu));
        ring_elem term = basisElementFromIndex(powerSumBasisId, "p", 10, true, mu);
        result = add(result, scaled(coeff, term));
      }
    cache[n] = result;
    return copyPolyValue(polyValue(result));
  }

CoeffMap SymmetricEngineRing::raisingExpansion(const Partition& lambda) const
{
    Partition trimmed = lambda;
    while (!trimmed.empty() && trimmed.back() == 0) trimmed.pop_back();
    CoeffMap current{{trimmed, coefficientRing->one()}};
    size_t ell = trimmed.size();
    for (size_t i = 0; i + 1 < ell; ++i)
      for (size_t j = i + 1; j < ell; ++j)
        {
          CoeffMap next;
          for (const auto& item : current)
            {
              const Partition& comp = item.first;
              ring_elem coeff = item.second;
              int maxRaise = std::max(comp[j], 0);
              for (int k = 0; k <= maxRaise; ++k)
                {
                  Partition newComp = comp;
                  if (k > 0)
                    {
                      newComp[i] += k;
                      newComp[j] -= k;
                    }
                  ring_elem factor = coefficientRing->one();
                  if (k > 0)
                    {
                      ring_elem tMinusOne =
                          coefficientRing->subtract(hallLittlewoodParameter,
                                                    coefficientRing->one());
                      factor = coefficientRing->mult(
                          tMinusOne,
                          coefficientRing->power(hallLittlewoodParameter, k - 1));
                    }
                  ring_elem contribution = coefficientRing->mult(coeff, factor);
                  auto existing = next.find(newComp);
                  if (existing == next.end())
                    next[newComp] = contribution;
                  else
                    {
                      ring_elem sum = coefficientRing->add(existing->second, contribution);
                      if (coefficientRing->is_zero(sum))
                        next.erase(existing);
                      else
                        existing->second = sum;
                    }
                }
            }
          current = next;
        }
    return current;
  }

CoeffMap SymmetricEngineRing::raisingGeneratorMap(const Partition& lambda) const
{
    CoeffMap result;
    for (const auto& item : raisingExpansion(lambda))
      addCoeff(result, item.first, item.second);
    return result;
  }

ring_elem SymmetricEngineRing::hallCapitalToPowerSums(const Partition& lambda, bool omega) const
{
    ring_elem result = zero();
    for (const auto& item : raisingExpansion(lambda))
      {
        bool invalid = false;
        for (int part : item.first)
          if (part < 0) invalid = true;
        if (invalid) continue;

        ring_elem term = one();
        for (int part : item.first)
          term = mult(term, hallLittlewoodPartToPowerSums(part, omega));
        result = add(result, scaled(item.second, term));
      }
    return result;
  }

ring_elem SymmetricEngineRing::hallPToPowerSums(const Partition& lambda, bool omega) const
{
    ring_elem numerator = hallCapitalToPowerSums(lambda, omega);
    return scaled(coefficientQuotient(coefficientRing->one(),
                                      hallLittlewoodCFactor(lambda)),
                  numerator);
  }

ring_elem SymmetricEngineRing::schurToPowerSums(const Partition& lambda) const
{
    int n = 0;
    for (int part : lambda) n += part;
    ring_elem result = zero();
    const CharacterTable& table = characterTable(n);
    auto lambdaRow = table.partitionRows.find(lambda);
    if (lambdaRow == table.partitionRows.end())
      {
        for (const auto& mu : table.partitions)
          {
            int chi = characterValue(lambda, mu);
            if (chi == 0) continue;
            ring_elem coeff = rationalCoefficient(chi, zValue(mu));
            if (error()) return zero();
            ring_elem term = basisElementFromIndex(powerSumBasisId, "p", 10, true, mu);
            result = add(result, scaled(coeff, term));
          }
      }
    else
      {
        size_t row = lambdaRow->second;
        for (size_t col = 0; col < table.partitions.size(); ++col)
          {
            int chi = characterTableValue(table, row, col);
            if (chi == 0) continue;
            ring_elem coeff = rationalCoefficient(chi, table.zValues[col]);
            if (error()) return zero();
            ring_elem term = basisElementFromIndex(powerSumBasisId,
                                                   "p",
                                                   10,
                                                   true,
                                                   table.partitions[col]);
            result = add(result, scaled(coeff, term));
          }
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumPartToComplete(int n, int hId, int hOrder) const
{
    if (n == 0) return one();
    auto cached = powerSumToCompleteCache.find(n);
    if (cached != powerSumToCompleteCache.end()) return copyPolyValue(polyValue(cached->second));
    ring_elem result =
        scaled(coefficientRing->from_long(n),
               basisElementFromIndex(hId, "h", hOrder, true, Partition{n}));
    for (int i = 1; i < n; ++i)
      {
        ring_elem pI = powerSumPartToComplete(i, hId, hOrder);
        ring_elem hRest = basisElementFromIndex(hId, "h", hOrder, true, Partition{n - i});
        result = subtract(result, mult(pI, hRest));
      }
    powerSumToCompleteCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::powerSumPartToElementary(int n, int eId, int eOrder) const
{
    if (n == 0) return one();
    auto cached = powerSumToElementaryCache.find(n);
    if (cached != powerSumToElementaryCache.end()) return copyPolyValue(polyValue(cached->second));
    ring_elem result =
        scaled(coefficientRing->from_long(n),
               basisElementFromIndex(eId, "e", eOrder, true, Partition{n}));
    for (int i = 1; i < n; ++i)
      {
        long sign = (i % 2 == 1) ? 1 : -1;
        ring_elem pI = powerSumPartToElementary(i, eId, eOrder);
        ring_elem eRest = basisElementFromIndex(eId, "e", eOrder, true, Partition{n - i});
        result = subtract(result, scaled(coefficientRing->from_long(sign), mult(pI, eRest)));
      }
    if (n % 2 == 0) result = negate(result);
    powerSumToElementaryCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::powerSumPartToHallGenerator(int n,
                                        int generatorId,
                                        const std::string& display,
                                        int generatorOrder,
                                        bool omega) const
{
    if (n == 0) return one();
    auto& cache = omega ? powerSumToBGeneratorCache : powerSumToQGeneratorCache;
    auto cached = cache.find(n);
    if (cached != cache.end()) return copyPolyValue(polyValue(cached->second));

    ring_elem tPower = coefficientRing->power(hallLittlewoodParameter, n);
    ring_elem factor = coefficientRing->subtract(coefficientRing->one(), tPower);
    if (omega && n % 2 == 0) factor = coefficientRing->negate(factor);
    ring_elem leadingCoeff = coefficientQuotient(coefficientRing->from_long(n), factor);
    if (error()) return zero();
    ring_elem result =
        scaled(leadingCoeff,
               basisElementFromIndex(generatorId,
                                     display,
                                     generatorOrder,
                                     true,
                                     Partition{n}));

    for (int i = 1; i < n; ++i)
      {
        ring_elem numerator =
            coefficientRing->subtract(coefficientRing->one(),
                                      coefficientRing->power(hallLittlewoodParameter, i));
        if (omega && i % 2 == 0) numerator = coefficientRing->negate(numerator);
        ring_elem coeff = coefficientRing->negate(coefficientQuotient(numerator, factor));
        if (error()) return zero();
        ring_elem rest = basisElementFromIndex(generatorId,
                                               display,
                                               generatorOrder,
                                               true,
                                               Partition{n - i});
        ring_elem previous = powerSumPartToHallGenerator(i,
                                                         generatorId,
                                                         display,
                                                         generatorOrder,
                                                         omega);
        if (error()) return zero();
        result = add(result, scaled(coeff, mult(rest, previous)));
      }

    cache[n] = result;
    return copyPolyValue(polyValue(result));
  }

CoeffMap SymmetricEngineRing::powerSumPartToHallGeneratorMap(int n, bool omega) const
{
    if (n == 0) return oneCoeffMap();
    ring_elem tPower = coefficientRing->power(hallLittlewoodParameter, n);
    ring_elem factor = coefficientRing->subtract(coefficientRing->one(), tPower);
    if (omega && n % 2 == 0) factor = coefficientRing->negate(factor);

    CoeffMap result;
    ring_elem leadingCoeff = coefficientQuotient(coefficientRing->from_long(n), factor);
    if (error()) return CoeffMap{};
    addCoeff(result, Partition{n}, leadingCoeff);

    for (int i = 1; i < n; ++i)
      {
        ring_elem numerator =
            coefficientRing->subtract(coefficientRing->one(),
                                      coefficientRing->power(hallLittlewoodParameter, i));
        if (omega && i % 2 == 0) numerator = coefficientRing->negate(numerator);
        ring_elem coeff = coefficientRing->negate(coefficientQuotient(numerator, factor));
        if (error()) return CoeffMap{};
        CoeffMap rest{{Partition{n - i}, coefficientRing->one()}};
        CoeffMap previous = powerSumPartToHallGeneratorMap(i, omega);
        if (error()) return CoeffMap{};
        CoeffMap product = multiplyCoeffMaps(rest, previous);
        result = addCoeffMaps(result, scaledCoeffMap(coeff, product));
      }

    return result;
  }

CoeffMap SymmetricEngineRing::powerSumIndexToHallGeneratorMap(const Partition& index, bool omega) const
{
    CoeffMap result = oneCoeffMap();
    for (int part : index)
      result = multiplyCoeffMaps(result, powerSumPartToHallGeneratorMap(part, omega));
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumsToHallGeneratorMap(ring_elem f, bool omega) const
{
    const auto *poly = polyValue(f);
    CoeffMap result;
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during basis conversion");
            return CoeffMap{};
          }
        CoeffMap converted = powerSumIndexToHallGeneratorMap(index, omega);
        result = addCoeffMaps(result, scaledCoeffMap(term.coeff, converted));
      }
    return result;
  }

CoeffMap SymmetricEngineRing::triangularReduceHallCapital(const CoeffMap& generatorMap,
                                       bool omega) const
{
    (void)omega;
    CoeffMap current = generatorMap;
    CoeffMap result;
    while (!current.empty())
      {
        Partition lambda = leadingPartition(current);
        CoeffMap expansion = raisingGeneratorMap(lambda);
        auto lead = expansion.find(lambda);
        if (lead == expansion.end() || coefficientRing->is_zero(lead->second))
          {
            ERROR("triangular expansion has zero leading coefficient");
            return CoeffMap{};
          }
        ring_elem c = coefficientQuotient(current[lambda], lead->second);
        addCoeff(result, lambda, c);
        current = addCoeffMaps(current,
                               scaledCoeffMap(coefficientRing->negate(c), expansion));
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumToSchurLike(const Partition& mu,
                                int schurId,
                                int schurOrder,
                                const std::string& display,
                                long sign) const
{
    int n = 0;
    for (int part : mu) n += part;
    ring_elem result = zero();
    const CharacterTable& table = characterTable(n);
    auto muCol = table.partitionRows.find(mu);
    if (muCol == table.partitionRows.end())
      {
        for (const auto& lambda : table.partitions)
          {
            int chi = characterValue(lambda, mu);
            if (chi == 0) continue;
            ring_elem term = basisElementFromIndex(schurId, display, schurOrder, false, lambda);
            result = add(result, scaled(coefficientRing->from_long(sign * chi), term));
          }
      }
    else
      {
        size_t col = muCol->second;
        for (size_t row = 0; row < table.partitions.size(); ++row)
          {
            int chi = characterTableValue(table, row, col);
            if (chi == 0) continue;
            ring_elem term = basisElementFromIndex(schurId,
                                                   display,
                                                   schurOrder,
                                                   false,
                                                   table.partitions[row]);
            result = add(result, scaled(coefficientRing->from_long(sign * chi), term));
          }
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumToSchur(const Partition& mu, int schurId, int schurOrder) const
{
    return powerSumToSchurLike(mu, schurId, schurOrder, "S", 1);
  }

ring_elem SymmetricEngineRing::powerSumsToSchurLike(ring_elem f,
                                 int schurId,
                                 int schurOrder,
                                 const std::string& display,
                                 bool omegaStyle) const
{
    const auto *poly = polyValue(f);

    GCMap<int, CoeffMap> powerSumCoeffsByDegree;
    for (const auto& term : poly->terms)
      {
        Partition mu;
        if (!powerSumIndexFromMonomial(term.monomial, mu))
          {
            ERROR("expected a pure power-sum expression during basis conversion");
            return zero();
          }
        addCoeff(powerSumCoeffsByDegree[partitionWeight(mu)], mu, term.coeff);
      }

    CoeffMap schurCoeffs;
    for (const auto& degreeData : powerSumCoeffsByDegree)
      {
        int degree = degreeData.first;
        const CoeffMap& powerSumCoeffs = degreeData.second;
        std::vector<Partition> inputPartitions;
        RingElemVector inputCoeffs;
        inputPartitions.reserve(powerSumCoeffs.size());
        inputCoeffs.reserve(powerSumCoeffs.size());
        for (const auto& muCoeff : powerSumCoeffs)
          {
            inputPartitions.push_back(muCoeff.first);
            inputCoeffs.push_back(muCoeff.second);
          }

        const std::vector<SchurConversionRecipeEntry>& recipe =
            powerSumsToSchurRecipe(degree, inputPartitions, omegaStyle);
        for (const auto& entry : recipe)
          {
            ring_elem coeff = coefficientRing->zero();
            for (const auto& contribution : entry.contributions)
              {
                size_t input = contribution.first;
                int chi = contribution.second;
                ring_elem termCoeff =
                    coefficientRing->mult(coefficientRing->from_long(chi),
                                          inputCoeffs[input]);
                coeff = coefficientRing->add(coeff, termCoeff);
              }
            addCoeff(schurCoeffs, entry.lambda, coeff);
          }
      }

    return coeffMapToElement(schurCoeffs, schurId, display, schurOrder, false);
  }

RingElemVector SymmetricEngineRing::solveSquareSystem(RingElemMatrix M,
                                   RingElemVector v) const
{
    size_t n = v.size();
    for (size_t col = 0; col < n; ++col)
      {
        size_t pivot = n;
        for (size_t r = col; r < n; ++r)
          if (!coefficientRing->is_zero(M[r][col]))
            {
              pivot = r;
              break;
            }
        if (pivot == n)
          {
            ERROR("basis conversion matrix is singular");
            return {};
          }
        if (pivot != col)
          {
            std::swap(M[pivot], M[col]);
            std::swap(v[pivot], v[col]);
          }
        ring_elem pivotValue = M[col][col];
        for (size_t c = col; c < n; ++c)
          M[col][c] = coefficientQuotient(M[col][c], pivotValue);
        v[col] = coefficientQuotient(v[col], pivotValue);
        for (size_t r = 0; r < n; ++r)
          if (r != col && !coefficientRing->is_zero(M[r][col]))
            {
              ring_elem factor = M[r][col];
              for (size_t c = col; c < n; ++c)
                M[r][c] = coefficientRing->subtract(
                    M[r][c],
                    coefficientRing->mult(factor, M[col][c]));
              v[r] = coefficientRing->subtract(v[r],
                                               coefficientRing->mult(factor, v[col]));
            }
      }
    return v;
  }

ring_elem SymmetricEngineRing::monomialBasisToPowerSums(const Partition& lambda, bool forgotten) const
{
    Partition key = normalizePartition(lambda);
    int d = partitionWeight(key);
    std::string cacheKey = partitionKey(key);
    auto& degreeCache = monomialToPowerSumCache[d];
    auto cached = degreeCache.find(cacheKey);
    ring_elem monomialResult;
    if (cached != degreeCache.end())
      {
        monomialResult = copyPolyValue(polyValue(cached->second));
      }
    else
      {
        std::vector<Partition> parts = partitionsOf(d);
        size_t n = parts.size();
        RingElemMatrix M(n, RingElemVector(n, coefficientRing->zero()));
        RingElemVector v(n, coefficientRing->zero());
        for (size_t row = 0; row < n; ++row)
          {
            if (parts[row] == key) v[row] = coefficientRing->one();
            for (size_t col = 0; col < n; ++col)
              M[row][col] =
                  coefficientRing->from_long(pToMonomialCoefficient(parts[col], parts[row]));
          }
        RingElemVector coeffs = solveSquareSystem(M, v);
        if (error()) return zero();
        ring_elem result = zero();
        for (size_t i = 0; i < n; ++i)
          if (!coefficientRing->is_zero(coeffs[i]))
            {
              ring_elem term = basisElementFromIndex(powerSumBasisId,
                                                     "p",
                                                     10,
                                                     true,
                                                     parts[i]);
              result = add(result, scaled(coeffs[i], term));
            }
        degreeCache[cacheKey] = result;
        monomialResult = copyPolyValue(polyValue(result));
      }
    return forgotten ? omegaPowerSums(monomialResult) : monomialResult;
  }

ring_elem SymmetricEngineRing::powerSumIndexToMonomialTarget(const Partition& lambda,
                                          int targetBasisId,
                                          const std::string& targetDisplay,
                                          int targetDisplayOrder,
                                          bool forgotten) const
{
    int d = partitionWeight(lambda);
    ring_elem result = zero();
    long sign = ((d - partitionLength(lambda)) % 2 == 0) ? 1 : -1;
    for (const auto& mu : partitionsOf(d))
      {
        long count = pToMonomialCoefficient(lambda, mu);
        if (count == 0) continue;
        ring_elem coeff = coefficientRing->from_long(forgotten ? sign * count : count);
        ring_elem term = basisElementFromIndex(targetBasisId,
                                               targetDisplay,
                                               targetDisplayOrder,
                                               false,
                                               mu);
        result = add(result, scaled(coeff, term));
      }
    return result;
  }

ring_elem SymmetricEngineRing::coeffMapToElement(const CoeffMap& H,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetDisplayOrder,
                              bool targetIsMultiplicative) const
{
    VECTOR(SymmetricTerm) terms;
    terms.reserve(H.size());
    for (const auto& item : H)
      {
        if (coefficientRing->is_zero(item.second)) continue;
        rememberBasis(targetBasisId,
                      targetDisplay,
                      targetDisplayOrder,
                      targetIsMultiplicative);
        SymmetricMonomial monomial;
        appendAtomBlock(monomial,
                        makeAtomBlock(targetDisplayOrder,
                                      targetBasisId,
                                      0,
                                      item.first));
        terms.push_back({item.second, canonicalMonomial(monomial)});
      }
    return fromTermVector(terms, false);
  }

int SymmetricEngineRing::requiredBasisIdForDisplay(const std::string& display) const
{
    int id = basisIdForDisplay(display);
    if (id < 0) ERROR("basis metadata for ", display.c_str(), " is not available");
    return id;
  }

ring_elem SymmetricEngineRing::basisElementForDisplay(const std::string& display,
                                   const Partition& index) const
{
    int id = requiredBasisIdForDisplay(display);
    if (error()) return zero();
    return basisElementFromIndex(id,
                                 display,
                                 basisOrderForId(id),
                                 isMultiplicativeBasis(id),
                                 index);
  }

ring_elem SymmetricEngineRing::replaceSingleBasis(ring_elem f,
                               int sourceBasisId,
                               const std::string& targetDisplay) const
{
    int targetId = requiredBasisIdForDisplay(targetDisplay);
    if (error()) return zero();
    CoeffMap coeffs = coefficientsInBasis(f, sourceBasisId);
    if (error()) return zero();
    return coeffMapToElement(coeffs,
                             targetId,
                             targetDisplay,
                             basisOrderForId(targetId),
                             isMultiplicativeBasis(targetId));
  }

ring_elem SymmetricEngineRing::powerSumsToHallCapitalTarget(ring_elem f,
                                         int targetBasisId,
                                         const std::string& targetDisplay,
                                         int targetDisplayOrder) const
{
    bool omega = targetDisplay == "B" || targetDisplay == "R";
    bool normalized = targetDisplay == "P" || targetDisplay == "R";
    CoeffMap generators = powerSumsToHallGeneratorMap(f, omega);
    if (error()) return zero();
    CoeffMap capitals = triangularReduceHallCapital(generators, omega);
    if (error()) return zero();
    if (normalized)
      {
        CoeffMap adjusted;
        for (const auto& item : capitals)
          addCoeff(adjusted,
                   item.first,
                   coefficientRing->mult(item.second, hallLittlewoodCFactor(item.first)));
        capitals = adjusted;
      }
    return coeffMapToElement(capitals,
                             targetBasisId,
                             targetDisplay,
                             targetDisplayOrder,
                             false);
  }

CoeffMap SymmetricEngineRing::schurGeneratorMap(const Partition& lambda,
                             bool omegaStyle,
                             int generatorId,
                             const std::string& generatorDisplay) const
{
    ring_elem expansion = jacobiTrudi(lambda,
                                      Partition{},
                                      generatorId,
                                      generatorDisplay,
                                      basisOrderForId(generatorId),
                                      isMultiplicativeBasis(generatorId));
    if (error()) return CoeffMap{};
    return coefficientsInBasis(expansion, generatorId);
  }

CoeffMap SymmetricEngineRing::triangularReduceSchur(const CoeffMap& generatorMap,
                                 bool omegaStyle,
                                 int generatorId,
                                 const std::string& generatorDisplay) const
{
    CoeffMap current = generatorMap;
    CoeffMap result;
    while (!current.empty())
      {
        Partition lambda = leadingPartition(current);
        CoeffMap expansion =
            schurGeneratorMap(lambda, omegaStyle, generatorId, generatorDisplay);
        if (error()) return CoeffMap{};
        auto lead = expansion.find(lambda);
        if (lead == expansion.end() || coefficientRing->is_zero(lead->second))
          {
            ERROR("triangular expansion has zero leading coefficient");
            return CoeffMap{};
          }
        ring_elem c = coefficientQuotient(current[lambda], lead->second);
        addCoeff(result, lambda, c);
        current = addCoeffMaps(current,
                               scaledCoeffMap(coefficientRing->negate(c), expansion));
      }
    return result;
  }

bool SymmetricEngineRing::directTriangularSchurConversion(ring_elem f,
                                       int targetBasisId,
                                       const std::string& targetDisplay,
                                       int targetDisplayOrder,
                                       ring_elem& result) const
{
    std::string generatorDisplay;
    if (targetDisplay == "S")
      generatorDisplay = "h";
    else if (targetDisplay == "Somega")
      generatorDisplay = "e";
    else
      return false;

    int generatorId = requiredBasisIdForDisplay(generatorDisplay);
    if (error()) return false;
    CoeffMap generatorCoeffs;
    if (!coefficientsInBasisIfPossible(f, generatorId, generatorCoeffs))
      return false;
    CoeffMap targetCoeffs = triangularReduceSchur(generatorCoeffs,
                                                  targetDisplay == "Somega",
                                                  generatorId,
                                                  generatorDisplay);
    if (error()) return false;
    result = coeffMapToElement(targetCoeffs,
                               targetBasisId,
                               targetDisplay,
                               targetDisplayOrder,
                               false);
    return true;
  }

bool SymmetricEngineRing::directTriangularHallConversion(ring_elem f,
                                      int targetBasisId,
                                      const std::string& targetDisplay,
                                      int targetDisplayOrder,
                                      ring_elem& result) const
{
    std::string generatorDisplay;
    bool omega = false;
    bool normalized = false;
    if (targetDisplay == "Q" || targetDisplay == "P")
      {
        generatorDisplay = "q";
        normalized = targetDisplay == "P";
      }
    else if (targetDisplay == "B" || targetDisplay == "R")
      {
        generatorDisplay = "b";
        omega = true;
        normalized = targetDisplay == "R";
      }
    else
      return false;

    int generatorId = requiredBasisIdForDisplay(generatorDisplay);
    if (error()) return false;
    CoeffMap generatorCoeffs;
    if (!coefficientsInBasisIfPossible(f, generatorId, generatorCoeffs))
      return false;
    CoeffMap capitals = triangularReduceHallCapital(generatorCoeffs, omega);
    if (error()) return false;
    if (normalized)
      {
        CoeffMap adjusted;
        for (const auto& item : capitals)
          addCoeff(adjusted,
                   item.first,
                   coefficientRing->mult(item.second, hallLittlewoodCFactor(item.first)));
        capitals = adjusted;
      }
    result = coeffMapToElement(capitals,
                               targetBasisId,
                               targetDisplay,
                               targetDisplayOrder,
                               false);
    return true;
  }

Partition SymmetricEngineRing::atomIndex(const SymmetricMonomial& monomial, size_t pos) const
{
    Partition result;
    int n = atomIndexLengthAt(monomial, pos);
    result.reserve(n);
    for (int i = 0; i < n; ++i) result.push_back(monomial.data[pos + atomHeaderSize + i]);
    return result;
  }

Partition SymmetricEngineRing::atomOuterIndex(const SymmetricMonomial& monomial, size_t pos) const
{
    Partition result;
    int n = atomOuterLengthAt(monomial, pos);
    result.reserve(n);
    for (int i = 0; i < n; ++i) result.push_back(monomial.data[pos + atomHeaderSize + i]);
    return result;
  }

Partition SymmetricEngineRing::atomInnerIndex(const SymmetricMonomial& monomial, size_t pos) const
{
    Partition result;
    int outerLength = atomOuterLengthAt(monomial, pos);
    int innerLength = atomInnerLengthAt(monomial, pos);
    result.reserve(innerLength);
    for (int i = 0; i < innerLength; ++i)
      result.push_back(monomial.data[pos + atomHeaderSize + outerLength + i]);
    return result;
  }

ring_elem SymmetricEngineRing::atomToPowerSums(const SymmetricMonomial& monomial, size_t pos) const
{
    int basisId = atomBasisIdAt(monomial, pos);
    std::string display = displayForBasis(basisId);
    if (atomIsSkewAt(monomial, pos))
      {
        Partition outer = atomOuterIndex(monomial, pos);
        Partition inner = atomInnerIndex(monomial, pos);
        if (display == "S")
          {
            int hId = requiredBasisIdForDisplay("h");
            if (error()) return zero();
            return elementToPowerSums(jacobiTrudi(outer,
                                                  inner,
                                                  hId,
                                                  "h",
                                                  basisOrderForId(hId),
                                                  isMultiplicativeBasis(hId)));
          }
        if (display == "Somega")
          {
            int eId = requiredBasisIdForDisplay("e");
            if (error()) return zero();
            return elementToPowerSums(jacobiTrudi(outer,
                                                  inner,
                                                  eId,
                                                  "e",
                                                  basisOrderForId(eId),
                                                  isMultiplicativeBasis(eId)));
          }
        if (display == "Q" || display == "B" || display == "P" || display == "R")
          return skewHallLittlewoodToPowerSums(outer, inner, display);
        ERROR("basis conversion for skew ", display.c_str(), " atoms is not implemented yet");
        return zero();
      }
    Partition index = atomIndex(monomial, pos);

    if (display == "p")
      return basisElementFromIndex(powerSumBasisId, "p", 10, true, index);

    if (display == "h" || display == "e")
      {
        ring_elem result = one();
        for (int part : index)
          {
            ring_elem factor = (display == "h") ? hPartToPowerSums(part)
                                                : ePartToPowerSums(part);
            result = mult(result, factor);
          }
        return result;
      }

    if (display == "q" || display == "b")
      {
        ring_elem result = one();
        for (int part : index)
          result = mult(result, hallLittlewoodPartToPowerSums(part, display == "b"));
        return result;
      }

    if (display == "S")
      return schurToPowerSums(index);

    if (display == "Somega")
      return omegaPowerSums(schurToPowerSums(index));

    if (display == "m" || isForgottenDisplay(display))
      return monomialBasisToPowerSums(index, isForgottenDisplay(display));

    if (display == "Q" || display == "B")
      return hallCapitalToPowerSums(index, display == "B");

    if (display == "P" || display == "R")
      return hallPToPowerSums(index, display == "R");

    ERROR("basis conversion to power sums is not implemented for basis ", display.c_str());
    return zero();
  }

ring_elem SymmetricEngineRing::monomialToPowerSums(const SymmetricMonomial& monomial) const
{
    ring_elem result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        result = mult(result, atomToPowerSums(monomial, pos));
        if (error()) return zero();
        pos += atomLengthAt(monomial, pos);
      }
    return result;
  }

ring_elem SymmetricEngineRing::elementToPowerSums(ring_elem f) const
{
    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        ring_elem converted = monomialToPowerSums(term.monomial);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, converted));
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumElementFromIndex(const Partition& index) const
{
    Partition normalized = normalizePartition(index);
    if (normalized.empty()) return one();
    return basisElementFromIndex(powerSumBasisId, "p", 10, true, normalized);
  }

CoeffMap SymmetricEngineRing::multiplySchurCoeffMapByRow(const CoeffMap& source, int row) const
{
    CoeffMap result;
    if (row == 0) return source;
    if (row < 0)
      {
        ERROR("expected nonnegative h index during recursive h-to-Schur conversion");
        return result;
      }
    Partition rowPartition{row};
    for (const auto& term : source)
      for (const auto& product : lrProduct(term.first, rowPartition))
        {
          ring_elem coeff = product.coefficient == 1
              ? term.second
              : coefficientRing->mult(coefficientRing->from_long(product.coefficient),
                                      term.second);
          addCoeff(result, product.nu, coeff);
        }
    return result;
  }

CoeffMap SymmetricEngineRing::hToSchurRecTransCoeffs(const CoeffMap& hCoeffs) const
{
    int lead = 0;
    for (const auto& item : hCoeffs)
      if (!item.first.empty())
        lead = std::max(lead, item.first.front());

    if (lead == 0)
      {
        CoeffMap result;
        for (const auto& item : hCoeffs)
          if (trimTrailingZerosPartition(item.first).empty())
            addCoeff(result, Partition{}, item.second);
          else
            {
              ERROR("invalid h-basis monomial during recursive h-to-Schur conversion");
              return CoeffMap{};
            }
        return result;
      }

    std::map<int, CoeffMap, std::greater<int>> grouped;
    int maxExponent = 0;
    for (const auto& item : hCoeffs)
      {
        Partition rest;
        int exponent = 0;
        for (int part : item.first)
          {
            if (part == lead)
              ++exponent;
            else
              rest.push_back(part);
          }
        maxExponent = std::max(maxExponent, exponent);
        addCoeff(grouped[exponent], rest, item.second);
      }

    CoeffMap result;
    for (int exponent = maxExponent; exponent >= 0; --exponent)
      {
        result = multiplySchurCoeffMapByRow(result, lead);
        if (error()) return CoeffMap{};
        auto found = grouped.find(exponent);
        if (found != grouped.end())
          result = addCoeffMaps(result, hToSchurRecTransCoeffs(found->second));
        if (error()) return CoeffMap{};
      }
    return result;
  }

ring_elem SymmetricEngineRing::hToSchurViaRecTrans(ring_elem f,
                                int hBasisId,
                                const std::string& hDisplay,
                                int hOrder,
                                bool hIsMultiplicative,
                                int schurId,
                                const std::string& schurDisplay,
                                int schurOrder) const
{
    rememberBasis(hBasisId, hDisplay, hOrder, hIsMultiplicative);
    rememberBasis(schurId, schurDisplay, schurOrder, false);

    CoeffMap hCoeffs = coefficientsInBasis(f, hBasisId);
    if (error()) return zero();
    CoeffMap result = hToSchurRecTransCoeffs(hCoeffs);
    if (error()) return zero();
    return coeffMapToElement(result, schurId, schurDisplay, schurOrder, false);
  }

ring_elem SymmetricEngineRing::jacobiTrudiBasis(int basisId,
                             const std::string& display,
                             int order,
                             bool isMultiplicative,
                             const Partition& outer,
                             const Partition& inner) const
{
    return jacobiTrudi(outer, inner, basisId, display, order, isMultiplicative);
  }

ring_elem SymmetricEngineRing::toBasis(ring_elem f,
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
    ring_elem direct;
    if (elementToDirectTarget(f,
                              targetBasisId,
                              targetDisplay,
                              targetOrder,
                              targetIsMultiplicative,
                              direct))
      return direct;
    if (error()) return zero();
    ring_elem inPowerSums = elementToPowerSums(f);
    if (error()) return zero();
    return powerSumsToTarget(inPowerSums,
                             targetBasisId,
                             targetDisplay,
                             targetOrder,
                             targetIsMultiplicative);
  }

} // namespace symmetric_rings
