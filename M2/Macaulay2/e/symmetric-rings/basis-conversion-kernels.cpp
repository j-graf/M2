// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"
#include "rings/ZZ.hpp"

#include <gmpxx.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace symmetric_rings {

namespace {

using IntegerPolynomial = std::vector<mpz_class>;
using IntegerPolynomialMap = std::map<Partition, IntegerPolynomial>;

void trimIntegerPolynomial(IntegerPolynomial& polynomial)
{
    while (!polynomial.empty() && polynomial.back() == 0)
      polynomial.pop_back();
}

void addRaisingFactorProduct(IntegerPolynomial& target,
                             const IntegerPolynomial& source,
                             int raise)
{
    if (source.empty()) return;
    size_t needed = source.size() + static_cast<size_t>(raise);
    if (target.size() < needed) target.resize(needed);
    for (size_t degree = 0; degree < source.size(); ++degree)
      {
        if (raise == 0)
          target[degree] += source[degree];
        else
          {
            target[degree + static_cast<size_t>(raise)] += source[degree];
            target[degree + static_cast<size_t>(raise - 1)] -= source[degree];
          }
      }
    trimIntegerPolynomial(target);
}

mpz_class labeledCycleAssignmentsToCompleteFactors(
    const Partition& cycleType,
    const Partition& factorDegrees)
{
    using State = std::pair<size_t, Partition>;
    std::map<State, mpz_class> memo;
    std::function<mpz_class(size_t, Partition&)> count =
        [&](size_t position, Partition& remaining) -> mpz_class {
          if (position == cycleType.size())
            {
              for (int value : remaining)
                if (value != 0) return 0;
              return 1;
            }
          State state{position, remaining};
          auto cached = memo.find(state);
          if (cached != memo.end()) return cached->second;
          mpz_class result = 0;
          int cycle = cycleType[position];
          for (size_t factor = 0; factor < remaining.size(); ++factor)
            if (remaining[factor] >= cycle)
              {
                remaining[factor] -= cycle;
                result += count(position + 1, remaining);
                remaining[factor] += cycle;
              }
          memo.emplace(std::move(state), result);
          return result;
        };
    Partition remaining = factorDegrees;
    return count(0, remaining);
}

} // namespace

// ============================================================================
// Shared Conversion Utilities
// ============================================================================
// Coefficient maps, basis metadata, and stored basis-element inspection.

bool SymmetricEngineRing::singleBasisIndexFromMonomial(
    const SymmetricMonomial& monomial,
    int basisId,
    Partition& index) const
{
    index.clear();
    if (monomial.data.empty()) return false;
    if (atomIsSkewAt(monomial, 0)) return false;
    if (atomBasisIdAt(monomial, 0) != basisId) return false;
    if (atomLengthAt(monomial, 0) != monomial.data.size()) return false;
    index = basisElementIndex(monomial, 0);
    return true;
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

void SymmetricEngineRing::addScaledCoeffMap(CoeffMap& target,
                                             ring_elem coeff,
                                             const CoeffMap& source) const
{
    if (coefficientRing->is_zero(coeff)) return;
    for (const auto& item : source)
      addCoeff(target,
               item.first,
               coefficientRing->mult(coeff, item.second));
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

std::string SymmetricEngineRing::littlewoodRichardsonProductCacheKey(const Partition& lambda, const Partition& mu) const
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

Partition SymmetricEngineRing::leadingPartition(const CoeffMap& H) const
{
    if (H.empty()) return Partition{};
    return H.begin()->first;
  }

ring_elem SymmetricEngineRing::coeffMapToElement(const CoeffMap& H,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetDisplayOrder,
                              bool targetIsMultiplicative) const
{
    rememberBasis(targetBasisId,
                  targetDisplay,
                  targetDisplayOrder,
                  targetIsMultiplicative);
    VECTOR(SymmetricTerm) terms;
    terms.reserve(H.size());
    for (const auto& item : H)
      {
        if (coefficientRing->is_zero(item.second)) continue;
        SymmetricMonomial monomial;
        if (!item.first.empty())
          appendAtomBlock(monomial,
                          makeAtomBlock(targetDisplayOrder,
                                        targetBasisId,
                                        0,
                                        item.first));
        terms.push_back({item.second, std::move(monomial)});
      }
    return fromTermVector(terms, false);
  }

ring_elem SymmetricEngineRing::basisElementForKind(BasisKind kind,
                                   const Partition& index) const
{
    int id = requiredBasisIdForKind(kind);
    if (error()) return zero();
    return basisElementFromIndex(id,
                                 displayForBasis(id),
                                 basisOrderForId(id),
                                 isMultiplicativeBasis(id),
                                 index);
  }

ring_elem SymmetricEngineRing::replaceSingleBasis(ring_elem f,
                               int sourceBasisId,
                               int targetBasisId) const
{
    CoeffMap coeffs = coefficientsInBasis(f, sourceBasisId);
    if (error()) return zero();
    return coeffMapToElement(coeffs,
                             targetBasisId,
                             displayForBasis(targetBasisId),
                             basisOrderForId(targetBasisId),
                             isMultiplicativeBasis(targetBasisId));
  }

Partition SymmetricEngineRing::basisElementIndex(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    Partition result;
    int n = atomIndexLengthAt(monomial, pos);
    result.reserve(n);
    for (int i = 0; i < n; ++i) result.push_back(monomial.data[pos + atomHeaderSize + i]);
    return result;
  }

Partition SymmetricEngineRing::basisElementOuterIndex(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    Partition result;
    int n = atomOuterLengthAt(monomial, pos);
    result.reserve(n);
    for (int i = 0; i < n; ++i) result.push_back(monomial.data[pos + atomHeaderSize + i]);
    return result;
  }

Partition SymmetricEngineRing::basisElementInnerIndex(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    Partition result;
    int outerLength = atomOuterLengthAt(monomial, pos);
    int innerLength = atomInnerLengthAt(monomial, pos);
    result.reserve(innerLength);
    for (int i = 0; i < innerLength; ++i)
      result.push_back(monomial.data[pos + atomHeaderSize + outerLength + i]);
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
    index = basisElementIndex(monomial, 0);
    return true;
  }

ring_elem SymmetricEngineRing::powerSumElementFromIndex(const Partition& index) const
{
    Partition normalized = normalizePartition(index);
    if (normalized.empty()) return one();
    return basisElementFromIndex(powerSumBasisId,
                                 displayForBasis(powerSumBasisId),
                                 basisOrderForId(powerSumBasisId),
                                 true,
                                 normalized);
  }

// ============================================================================
// Jacobi-Trudi And Determinant Utilities
// ============================================================================
// Determinantal constructions shared by Schur-style conversion and straightening.

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
    auto& cache = basisKindForId(basisId) == BasisKind::Elementary
                      ? eJacobiTrudiCache : hJacobiTrudiCache;
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

ring_elem SymmetricEngineRing::jacobiTrudiBasis(int basisId,
                             const std::string& display,
                             int order,
                             bool isMultiplicative,
                             const Partition& outer,
                             const Partition& inner) const
{
    return jacobiTrudi(outer, inner, basisId, display, order, isMultiplicative);
  }

// ============================================================================
// Complete And Elementary Bases
// ============================================================================
// Classical h/e-to-p formulas, logarithmic inverse formulas, and h-to-S transition.

ring_elem SymmetricEngineRing::completePartToPowerSumsViaClassicalFormula(int n) const
{
    if (n == 0) return one();
    auto cached = completeToPowerSumsCache.find(n);
    if (cached != completeToPowerSumsCache.end()) return copyPolyValue(polyValue(cached->second));
    CoeffMap coefficients;
    for (const auto& mu : partitionsOf(n))
      {
        ring_elem coeff = rationalCoefficient(1, zValue(mu));
        if (error()) return zero();
        addCoeff(coefficients, mu, coeff);
      }
    ring_elem result = coeffMapToElement(
        coefficients,
        powerSumBasisId,
        displayForBasis(powerSumBasisId),
        basisOrderForId(powerSumBasisId),
        true);
    completeToPowerSumsCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::elementaryPartToPowerSumsViaClassicalFormula(int n) const
{
    if (n == 0) return one();
    auto cached = elementaryToPowerSumsCache.find(n);
    if (cached != elementaryToPowerSumsCache.end()) return copyPolyValue(polyValue(cached->second));
    CoeffMap coefficients;
    for (const auto& mu : partitionsOf(n))
      {
        long sign = ((n - static_cast<int>(mu.size())) % 2 == 0) ? 1 : -1;
        ring_elem coeff = rationalCoefficient(sign, zValue(mu));
        if (error()) return zero();
        addCoeff(coefficients, mu, coeff);
      }
    ring_elem result = coeffMapToElement(
        coefficients,
        powerSumBasisId,
        displayForBasis(powerSumBasisId),
        basisOrderForId(powerSumBasisId),
        true);
    elementaryToPowerSumsCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::powerSumLogarithmCoefficient(
    const Partition& lambda) const
{
    mpq_t value;
    mpq_init(value);
    mpz_fac_ui(mpq_numref(value),
               lambda.empty() ? 0 : lambda.size() - 1);
    if (lambda.size() % 2 == 0) mpz_neg(mpq_numref(value), mpq_numref(value));

    mpz_set_ui(mpq_denref(value), 1);
    std::map<int, int> multiplicities;
    for (int part : lambda) ++multiplicities[part];
    mpz_t factorial;
    mpz_init(factorial);
    for (const auto& item : multiplicities)
      {
        mpz_fac_ui(factorial, item.second);
        mpz_mul(mpq_denref(value), mpq_denref(value), factorial);
      }
    mpz_clear(factorial);
    mpq_canonicalize(value);

    ring_elem result;
    if (!coefficientRing->from_rational(value, result))
      {
        ERROR("coefficient division failed during logarithm basis conversion");
        result = coefficientRing->zero();
      }
    mpq_clear(value);
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumPartToGeneratorMapViaLogarithmFormula(
    int n,
    ring_elem common) const
{
    CoeffMap result;
    for (const auto& lambda : partitionsOf(n))
      {
        ring_elem logarithmCoeff = powerSumLogarithmCoefficient(lambda);
        if (error()) return CoeffMap{};
        addCoeff(result,
                 lambda,
                 coefficientRing->mult(common, logarithmCoeff));
      }
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumPartToCompleteMapViaLogarithmFormula(int n) const
{
    if (n == 0) return oneCoeffMap();
    auto cached = powerSumToCompleteMapCache.find(n);
    if (cached != powerSumToCompleteMapCache.end()) return cached->second;
    CoeffMap result = powerSumPartToGeneratorMapViaLogarithmFormula(
        n, cachedInteger(n));
    powerSumToCompleteMapCache[n] = result;
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumPartToElementaryMapViaLogarithmFormula(int n) const
{
    if (n == 0) return oneCoeffMap();
    auto cached = powerSumToElementaryMapCache.find(n);
    if (cached != powerSumToElementaryMapCache.end()) return cached->second;
    ring_elem common = cachedInteger(n);
    if (n % 2 == 0) common = coefficientRing->negate(common);
    CoeffMap result = powerSumPartToGeneratorMapViaLogarithmFormula(n, common);
    powerSumToElementaryMapCache[n] = result;
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumIndexToCompleteMapViaLogarithmFormula(
    const Partition& index) const
{
    CoeffMap result = oneCoeffMap();
    for (int part : index)
      result = multiplyCoeffMaps(
          result, powerSumPartToCompleteMapViaLogarithmFormula(part));
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumIndexToElementaryMapViaLogarithmFormula(
    const Partition& index) const
{
    CoeffMap result = oneCoeffMap();
    for (int part : index)
      result = multiplyCoeffMaps(
          result, powerSumPartToElementaryMapViaLogarithmFormula(part));
    return result;
  }

ring_elem SymmetricEngineRing::powerSumsToCompleteViaLogarithmFormula(
    ring_elem f, int completeId, int completeOrder) const
{
    CoeffMap result;
    for (const auto& term : polyValue(f)->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during basis conversion");
            return zero();
          }
        addScaledCoeffMap(
            result,
            term.coeff,
            powerSumIndexToCompleteMapViaLogarithmFormula(index));
      }
    return coeffMapToElement(
        result, completeId, displayForBasis(completeId), completeOrder, true);
  }

ring_elem SymmetricEngineRing::powerSumsToElementaryViaLogarithmFormula(
    ring_elem f, int elementaryId, int elementaryOrder) const
{
    CoeffMap result;
    for (const auto& term : polyValue(f)->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during basis conversion");
            return zero();
          }
        addScaledCoeffMap(
            result,
            term.coeff,
            powerSumIndexToElementaryMapViaLogarithmFormula(index));
      }
    return coeffMapToElement(result,
                             elementaryId,
                             displayForBasis(elementaryId),
                             elementaryOrder,
                             true);
  }

CoeffMap SymmetricEngineRing::multiplySchurExpansionViaRowPieri(const CoeffMap& source, int row) const
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
      for (const auto& product : littlewoodRichardsonProductViaCoefficientEnumeration(term.first, rowPartition))
        {
          ring_elem coeff = product.coefficient == 1
              ? term.second
              : coefficientRing->mult(coefficientRing->from_long(product.coefficient),
                                      term.second);
          addCoeff(result, product.nu, coeff);
        }
    return result;
  }

CoeffMap SymmetricEngineRing::completeToSchurCoefficientsViaRecursiveTransition(const CoeffMap& hCoeffs) const
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

    GCMap<int, CoeffMap, std::greater<int>> grouped;
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
        result = multiplySchurExpansionViaRowPieri(result, lead);
        if (error()) return CoeffMap{};
        auto found = grouped.find(exponent);
        if (found != grouped.end())
          result = addCoeffMaps(result, completeToSchurCoefficientsViaRecursiveTransition(found->second));
        if (error()) return CoeffMap{};
      }
    return result;
  }

ring_elem SymmetricEngineRing::completeToSchurViaRecursiveTransition(ring_elem f,
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
    CoeffMap result = completeToSchurCoefficientsViaRecursiveTransition(hCoeffs);
    if (error()) return zero();
    return coeffMapToElement(result, schurId, schurDisplay, schurOrder, false);
  }

// ============================================================================
// Schur And Omega-Schur Bases
// ============================================================================
// Character, omega, and triangular conversions for S and Somega.

ring_elem SymmetricEngineRing::schurOmegaConversionViaPartitionConjugation(
    ring_elem f,
    int sourceBasisId,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetOrder) const
{
    CoeffMap source = coefficientsInBasis(f, sourceBasisId);
    if (error()) return zero();
    CoeffMap result;
    for (const auto& item : source)
      addCoeff(result, conjugatePartition(item.first), item.second);
    return coeffMapToElement(
        result, targetBasisId, targetDisplay, targetOrder, false);
  }

ring_elem SymmetricEngineRing::schurLikeToPowerSumsViaCharacters(
    const Partition& lambda,
    bool omegaStyle) const
{
    int n = 0;
    for (int part : lambda) n += part;
    CoeffMap coefficients;
    const CharacterTable& table = characterTable(n);
    auto lambdaRow = table.partitionRows.find(lambda);
    if (lambdaRow == table.partitionRows.end())
      {
        for (const auto& mu : table.partitions)
          {
            int chi = characterValue(lambda, mu);
            if (chi == 0) continue;
            if (omegaStyle && ((n - static_cast<int>(mu.size())) % 2 == 1))
              chi = -chi;
            ring_elem coeff = rationalCoefficient(chi, zValue(mu));
            if (error()) return zero();
            addCoeff(coefficients, mu, coeff);
          }
      }
    else
      {
        size_t row = lambdaRow->second;
        for (size_t col = 0; col < table.partitions.size(); ++col)
          {
            int chi = characterTableValue(table, row, col);
            if (chi == 0) continue;
            if (omegaStyle &&
                ((n - static_cast<int>(table.partitions[col].size())) % 2 == 1))
              chi = -chi;
            ring_elem coeff = rationalCoefficient(chi, table.zValues[col]);
            if (error()) return zero();
            addCoeff(coefficients, table.partitions[col], coeff);
          }
      }
    return coeffMapToElement(coefficients,
                             powerSumBasisId,
                             displayForBasis(powerSumBasisId),
                             basisOrderForId(powerSumBasisId),
                             true);
  }

ring_elem SymmetricEngineRing::powerSumIndexToSchurLikeViaCharacters(const Partition& mu,
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

ring_elem SymmetricEngineRing::powerSumIndexToSchurViaCharacters(const Partition& mu, int schurId, int schurOrder) const
{
    return powerSumIndexToSchurLikeViaCharacters(
        mu, schurId, schurOrder, displayForBasis(schurId), 1);
  }

ring_elem SymmetricEngineRing::powerSumsToSchurLikeViaCharacters(ring_elem f,
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

ring_elem SymmetricEngineRing::powerSumsToSchurViaCharacters(
    ring_elem f, int schurId, int schurOrder) const
{
    return powerSumsToSchurLikeViaCharacters(
        f, schurId, schurOrder, displayForBasis(schurId), false);
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
        addScaledCoeffMap(current, coefficientRing->negate(c), expansion);
      }
    return result;
  }

bool SymmetricEngineRing::tryExpressionToSchurViaTriangularReduction(ring_elem f,
                                       int targetBasisId,
                                       const std::string& targetDisplay,
                                       int targetDisplayOrder,
                                       ring_elem& result) const
{
    const BasisKind targetKind = basisKindForId(targetBasisId);
    BasisKind generatorKind;
    switch (targetKind)
      {
      case BasisKind::Schur:
        generatorKind = BasisKind::Complete;
        break;
      case BasisKind::SchurOmega:
        generatorKind = BasisKind::Elementary;
        break;
      default:
        return false;
      }
    int generatorId = requiredBasisIdForKind(generatorKind);
    if (error()) return false;
    const std::string generatorDisplay = displayForBasis(generatorId);
    CoeffMap generatorCoeffs;
    if (!coefficientsInBasisIfPossible(f, generatorId, generatorCoeffs))
      return false;
    CoeffMap targetCoeffs = triangularReduceSchur(generatorCoeffs,
                                                  targetKind == BasisKind::SchurOmega,
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

// ============================================================================
// Monomial And Forgotten Bases
// ============================================================================
// Transition matrices and coefficient-map products for m and ff.

CoeffMap SymmetricEngineRing::multiplyMonomialCoeffMaps(const CoeffMap& a,
                                                        const CoeffMap& b) const
{
    CoeffMap result;
    for (const auto& left : a)
      for (const auto& right : b)
        {
          ring_elem baseCoeff = coefficientRing->mult(left.second, right.second);
          if (coefficientRing->is_zero(baseCoeff)) continue;
          for (const auto& product : monomialProductViaExponentSplittings(left.first, right.first))
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

ring_elem SymmetricEngineRing::monomialToPowerSumsViaTransitionMatrix(const Partition& lambda, bool forgotten) const
{
    Partition key = normalizePartition(lambda);
    int d = partitionWeight(key);
    std::string cacheKey = partitionKey(key);
    auto& degreeCache = forgotten ? forgottenToPowerSumCache[d]
                                  : monomialToPowerSumCache[d];
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
        auto keyPosition = std::find(parts.begin(), parts.end(), key);
        if (keyPosition == parts.end())
          {
            ERROR("invalid monomial partition during basis conversion");
            return zero();
          }
        size_t keyRow = static_cast<size_t>(keyPosition - parts.begin());
        RingElemVector coeffs(n, coefficientRing->zero());
        // Coarsenings precede refinements, making the p-to-m matrix triangular.
        for (size_t offset = 0; offset <= keyRow; ++offset)
          {
            size_t row = keyRow - offset;
            ring_elem rhs = row == keyRow
                ? coefficientRing->one()
                : coefficientRing->zero();
            for (size_t col = row + 1; col <= keyRow; ++col)
              {
                long transition = pToMonomialCoefficient(parts[col], parts[row]);
                if (transition == 0 || coefficientRing->is_zero(coeffs[col])) continue;
                rhs = coefficientRing->subtract(
                    rhs,
                    coefficientRing->mult(cachedInteger(transition), coeffs[col]));
              }
            long diagonal = pToMonomialCoefficient(parts[row], parts[row]);
            coeffs[row] = coefficientQuotient(rhs, cachedInteger(diagonal));
            if (error()) return zero();
          }
        CoeffMap coefficients;
        for (size_t i = 0; i <= keyRow; ++i)
          if (!coefficientRing->is_zero(coeffs[i]))
            {
              ring_elem coeff = coeffs[i];
              if (forgotten && ((d - partitionLength(parts[i])) % 2 == 1))
                coeff = coefficientRing->negate(coeff);
              addCoeff(coefficients, parts[i], coeff);
            }
        ring_elem result = coeffMapToElement(
            coefficients,
            powerSumBasisId,
            displayForBasis(powerSumBasisId),
            basisOrderForId(powerSumBasisId),
            true);
        degreeCache[cacheKey] = result;
        monomialResult = copyPolyValue(polyValue(result));
      }
    return monomialResult;
  }

ring_elem SymmetricEngineRing::powerSumIndexToMonomialViaTransitionMatrix(const Partition& lambda,
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

// ============================================================================
// Hall-Littlewood Generator Bases
// ============================================================================
// Classical and logarithmic conversions for the multiplicative q and b bases.

ring_elem SymmetricEngineRing::hallLittlewoodGeneratorPartToPowerSumsViaClassicalFormula(int n, bool omega) const
{
    if (n == 0) return one();
    auto& cache = omega ? hallLittlewoodBGeneratorToPowerSumsCache : hallLittlewoodQGeneratorToPowerSumsCache;
    auto cached = cache.find(n);
    if (cached != cache.end()) return copyPolyValue(polyValue(cached->second));
    CoeffMap coefficients;
    for (const auto& mu : partitionsOf(n))
      {
        long sign = 1;
        if (omega && ((n - static_cast<int>(mu.size())) % 2 == 1)) sign = -1;
        ring_elem rational = rationalCoefficient(sign, zValue(mu));
        if (error()) return zero();
        ring_elem coeff = coefficientRing->mult(rational, hallLittlewoodFactor(mu));
        addCoeff(coefficients, mu, coeff);
      }
    ring_elem result = coeffMapToElement(
        coefficients,
        powerSumBasisId,
        displayForBasis(powerSumBasisId),
        basisOrderForId(powerSumBasisId),
        true);
    cache[n] = result;
    return copyPolyValue(polyValue(result));
  }

CoeffMap SymmetricEngineRing::hallLittlewoodGeneratorPartToPowerSumsQuotientMapViaClassicalFormula(
    int n,
    bool omega) const
{
    if (n == 0) return oneCoeffMap();
    auto& cache = omega ? hallLittlewoodBGeneratorToPowerSumsQuotientMapCache
                        : hallLittlewoodQGeneratorToPowerSumsQuotientMapCache;
    auto cached = cache.find(n);
    if (cached != cache.end()) return cached->second;

    CoeffMap result;
    for (const auto& mu : partitionsOf(n))
      {
        long sign = 1;
        if (omega && ((n - static_cast<int>(mu.size())) % 2 == 1)) sign = -1;
        ring_elem rational = rationalCoefficient(sign, zValue(mu));
        if (error()) return CoeffMap{};
        addCoeff(result, mu, rational);
      }
    cache[n] = result;
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumPartToHallGeneratorMapViaLogarithmFormula(
    int n,
    bool omega) const
{
    if (n == 0) return oneCoeffMap();
    auto& cache = omega ? powerSumToBGeneratorMapCache
                        : powerSumToQGeneratorMapCache;
    auto cached = cache.find(n);
    if (cached != cache.end()) return cached->second;

    ring_elem factor = hallLittlewoodFactor(Partition{n});
    if (omega && n % 2 == 0) factor = coefficientRing->negate(factor);
    ring_elem common = coefficientQuotient(cachedInteger(n), factor);
    if (error()) return CoeffMap{};
    CoeffMap result =
        powerSumPartToGeneratorMapViaLogarithmFormula(n, common);

    cache[n] = result;
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumIndexToHallGeneratorMapViaLogarithmFormula(
    const Partition& index,
    bool omega) const
{
    CoeffMap result = oneCoeffMap();
    for (int part : index)
      result = multiplyCoeffMaps(
          result,
          powerSumPartToHallGeneratorMapViaLogarithmFormula(part, omega));
    return result;
  }

CoeffMap SymmetricEngineRing::powerSumsToHallGeneratorMapViaLogarithmFormula(
    ring_elem f,
    bool omega) const
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
        CoeffMap converted =
            powerSumIndexToHallGeneratorMapViaLogarithmFormula(index, omega);
        addScaledCoeffMap(result, term.coeff, converted);
      }
    return result;
  }

// ============================================================================
// Hall-Littlewood Capital And Normalized Bases
// ============================================================================
// Raising operators, Green polynomials, normalization, skew functions, and triangular reduction.

ring_elem SymmetricEngineRing::hallLittlewoodCFactor(const Partition& lambda) const
{
    Partition normalized = normalizePartition(lambda);
    auto cached = hallLittlewoodCFactorCache.find(normalized);
    if (cached != hallLittlewoodCFactorCache.end()) return cached->second;
    std::map<int, int> multiplicities;
    for (int part : normalized) multiplicities[part]++;
    ring_elem result = coefficientRing->one();
    for (const auto& item : multiplicities)
      for (int j = 1; j <= item.second; ++j)
        {
          ring_elem tPower = coefficientRing->power(hallLittlewoodParameter, j);
          result = coefficientRing->mult(
              result,
              coefficientRing->subtract(coefficientRing->one(), tPower));
        }
    hallLittlewoodCFactorCache[normalized] = result;
    return result;
  }

CoeffMap SymmetricEngineRing::hallLittlewoodSingleCycleGreenMap(
    int n,
    bool normalized) const
{
    auto& cache = normalized
        ? hallLittlewoodSingleCycleNormalizedGreenMapCache
        : hallLittlewoodSingleCycleCapitalGreenMapCache;
    auto cached = cache.find(n);
    if (cached != cache.end()) return cached->second;

    CoeffMap result;
    for (const Partition& lambda : partitionsOf(n))
      {
        int length = static_cast<int>(lambda.size());
        int nStatistic = 0;
        for (int i = 1; i < length; ++i) nStatistic += i * lambda[i];
        int exponent = nStatistic - length * (length - 1) / 2;
        ring_elem green = coefficientRing->power(
            hallLittlewoodParameter, exponent);
        for (int i = 1; i < length; ++i)
          green = coefficientRing->mult(
              green,
              coefficientRing->subtract(
                  coefficientRing->power(hallLittlewoodParameter, i),
                  coefficientRing->one()));

        ring_elem coefficient = normalized
            ? green
            : coefficientQuotient(green, hallLittlewoodCFactor(lambda));
        if (error()) return CoeffMap{};
        addCoeff(result, lambda, coefficient);
      }
    cache[n] = result;
    return result;
  }

ring_elem
SymmetricEngineRing::hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
    ring_elem f,
    int sourceBasisId,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetOrder,
    bool capitalToNormalized) const
{
    CoeffMap source = coefficientsInBasis(f, sourceBasisId);
    if (error()) return zero();
    CoeffMap result;
    for (const auto& item : source)
      {
        ring_elem factor = hallLittlewoodCFactor(item.first);
        ring_elem coefficient = capitalToNormalized
            ? coefficientRing->mult(item.second, factor)
            : coefficientQuotient(item.second, factor);
        if (error()) return zero();
        addCoeff(result, item.first, coefficient);
      }
    return coeffMapToElement(
        result, targetBasisId, targetDisplay, targetOrder, false);
  }

CoeffMap SymmetricEngineRing::raisingExpansion(const Partition& lambda) const
{
    Partition trimmed = lambda;
    while (!trimmed.empty() && trimmed.back() == 0) trimmed.pop_back();
    IntegerPolynomialMap current{{trimmed, IntegerPolynomial{mpz_class(1)}}};
    size_t ell = trimmed.size();
    // Replenish each donor coordinate before later operators drain it.
    for (size_t offset = 1; offset < ell; ++offset)
      {
        size_t i = ell - 1 - offset;
        for (size_t j = i + 1; j < ell; ++j)
        {
          IntegerPolynomialMap next;
          for (const auto& item : current)
            {
              const Partition& comp = item.first;
              int maxRaise = std::max(comp[j], 0);
              for (int k = 0; k <= maxRaise; ++k)
                {
                  Partition newComp = comp;
                  if (k > 0)
                    {
                      newComp[i] += k;
                      newComp[j] -= k;
                    }
                  auto inserted = next.emplace(newComp, IntegerPolynomial{});
                  addRaisingFactorProduct(inserted.first->second, item.second, k);
                  if (inserted.first->second.empty()) next.erase(inserted.first);
                }
            }
          current = next;
        }
      }

    IntegerPolynomialMap normalized;
    for (const auto& item : current)
      {
        Partition index = normalizePartition(item.first);
        auto inserted = normalized.emplace(index, IntegerPolynomial{});
        addRaisingFactorProduct(inserted.first->second, item.second, 0);
        if (inserted.first->second.empty()) normalized.erase(inserted.first);
      }

    CoeffMap result;
    for (const auto& item : normalized)
      {
        ring_elem coefficient = coefficientRing->zero();
        for (auto degree = item.second.rbegin(); degree != item.second.rend(); ++degree)
          {
            coefficient = coefficientRing->mult(
                coefficient, hallLittlewoodParameter);
            if (*degree != 0)
              coefficient = coefficientRing->add(
                  coefficient,
                  coefficientRing->from_int(degree->get_mpz_t()));
          }
        addCoeff(result, item.first, coefficient);
      }
    return result;
  }

CoeffMap SymmetricEngineRing::raisingGeneratorMap(const Partition& lambda) const
{
    auto cached = hallLittlewoodRaisingGeneratorMapCache.find(lambda);
    if (cached != hallLittlewoodRaisingGeneratorMapCache.end())
      return cached->second;
    CoeffMap result;
    for (const auto& item : raisingExpansion(lambda))
      addCoeff(result, item.first, item.second);
    hallLittlewoodRaisingGeneratorMapCache[lambda] = result;
    return result;
  }

ring_elem SymmetricEngineRing::hallLittlewoodCapitalToPowerSumsViaRaisingOperators(const Partition& lambda, bool omega) const
{
    CoeffMap result;
    for (const auto& item : raisingGeneratorMap(lambda))
      {
        CoeffMap term = oneCoeffMap();
        for (int part : item.first)
          term = multiplyCoeffMaps(
              term,
              hallLittlewoodGeneratorPartToPowerSumsQuotientMapViaClassicalFormula(part,
                                                                                    omega));
        addScaledCoeffMap(result, item.second, term);
      }
    // Every assembly of p_mu has the same product of (1-t^part) factors.
    for (auto& item : result)
      item.second = coefficientRing->mult(item.second,
                                          hallLittlewoodFactor(item.first));
    return coeffMapToElement(result,
                             powerSumBasisId,
                             displayForBasis(powerSumBasisId),
                             basisOrderForId(powerSumBasisId),
                             true);
  }

ring_elem SymmetricEngineRing::hallLittlewoodNormalizedToPowerSumsViaCapitalNormalization(const Partition& lambda, bool omega) const
{
    ring_elem numerator = hallLittlewoodCapitalToPowerSumsViaRaisingOperators(lambda, omega);
    return scaled(coefficientQuotient(coefficientRing->one(),
                                      hallLittlewoodCFactor(lambda)),
                  numerator);
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
        if (lead == expansion.end() ||
            !coefficientRing->is_equal(lead->second, coefficientRing->one()))
          {
            ERROR("Hall-Littlewood transition is not unitriangular");
            return CoeffMap{};
          }
        ring_elem c = current[lambda];
        addCoeff(result, lambda, c);
        current.erase(lambda);
        ring_elem negative = coefficientRing->negate(c);
        for (const auto& item : expansion)
          if (item.first != lambda)
            addCoeff(current,
                     item.first,
                     coefficientRing->mult(negative, item.second));
      }
    return result;
  }

ring_elem SymmetricEngineRing::skewQOrBFunction(
    const Partition& lambda,
    const Partition& mu,
    bool omega) const
{
    int d = partitionWeight(lambda) - partitionWeight(mu);
    if (d < 0) return zero();

    ring_elem lambdaTerm = basisElementForKind(
        omega ? BasisKind::HallLittlewoodB : BasisKind::HallLittlewoodQ,
        lambda);
    ring_elem muTerm = basisElementForKind(
        omega ? BasisKind::HallLittlewoodPOmega : BasisKind::HallLittlewoodP,
        mu);
    if (error()) return zero();

    int targetId = requiredBasisIdForKind(
        omega ? BasisKind::Forgotten : BasisKind::Monomial);
    std::string targetDisplay = displayForBasis(targetId);
    if (error()) return zero();

    ring_elem result = zero();
    for (const auto& nu : partitionsOf(d))
      {
        ring_elem generator = basisElementForKind(
            omega ? BasisKind::HallLittlewoodBGenerator
                  : BasisKind::HallLittlewoodQGenerator,
            nu);
        ring_elem test = mult(muTerm, generator);
        ring_elem coeff = hallInnerProductElements(
            lambdaTerm, test, InnerProductKind::HallLittlewood);
        if (error()) return zero();
        if (coefficientRing->is_zero(coeff)) continue;
        ring_elem term = basisElementFromIndex(targetId,
                                               targetDisplay,
                                               basisOrderForId(targetId),
                                               isMultiplicativeBasis(targetId),
                                               nu);
        result = add(result, scaled(coeff, term));
      }
    return result;
  }

ring_elem SymmetricEngineRing::skewPOrPOmegaToPowerSums(
    const Partition& lambda,
    const Partition& mu,
    bool omega) const
{
    ring_elem skewCapital = skewQOrBFunction(lambda, mu, omega);
    if (error()) return zero();
    int pId = requiredBasisIdForKind(BasisKind::PowerSum);
    const BasisKind capitalKind = omega ? BasisKind::HallLittlewoodB
                                        : BasisKind::HallLittlewoodQ;
    int capitalId = requiredBasisIdForKind(capitalKind);
    if (error()) return zero();
    const std::string pDisplay = displayForBasis(pId);
    const std::string capitalDisplay = displayForBasis(capitalId);

    ring_elem inCapital = toBasis(skewCapital,
                                  pId,
                                  pDisplay,
                                  basisOrderForId(pId),
                                  isMultiplicativeBasis(pId),
                                  capitalId,
                                  capitalDisplay,
                                  basisOrderForId(capitalId),
                                  isMultiplicativeBasis(capitalId));
    if (error()) return zero();
    const BasisKind normalizedKind = omega ? BasisKind::HallLittlewoodPOmega
                                           : BasisKind::HallLittlewoodP;
    const int normalizedId = requiredBasisIdForKind(normalizedKind);
    if (error()) return zero();
    ring_elem normalized = replaceSingleBasis(inCapital, capitalId, normalizedId);
    if (error()) return zero();
    return expressionToPowerSumsViaBasisElementRoutes(normalized);
  }

ring_elem SymmetricEngineRing::skewHallLittlewoodToPowerSums(
    const Partition& lambda,
    const Partition& mu,
    BasisKind basisKind) const
{
    switch (basisKind)
      {
      case BasisKind::HallLittlewoodQ:
      case BasisKind::HallLittlewoodB:
        return expressionToPowerSumsViaBasisElementRoutes(
            skewQOrBFunction(lambda, mu,
                             basisKind == BasisKind::HallLittlewoodB));
      case BasisKind::HallLittlewoodP:
      case BasisKind::HallLittlewoodPOmega:
        return skewPOrPOmegaToPowerSums(
            lambda, mu, basisKind == BasisKind::HallLittlewoodPOmega);
      default:
        ERROR("expected a skew Hall-Littlewood basis element");
        return zero();
      }
  }

ring_elem SymmetricEngineRing::powerSumsToHallLittlewoodCapitalViaTriangularReduction(ring_elem f,
                                         int targetBasisId,
                                         const std::string& targetDisplay,
                                         int targetDisplayOrder) const
{
    const BasisKind targetKind = basisKindForId(targetBasisId);
    bool omega = targetKind == BasisKind::HallLittlewoodB ||
                 targetKind == BasisKind::HallLittlewoodPOmega;
    bool normalized = targetKind == BasisKind::HallLittlewoodP ||
                      targetKind == BasisKind::HallLittlewoodPOmega;
    const auto *poly = polyValue(f);
    CoeffMap capitals;
    Partition singleIndex;
    bool singlePowerSumIndex =
        poly->terms.size() == 1 &&
        powerSumIndexFromMonomial(poly->terms.front().monomial, singleIndex) &&
        !singleIndex.empty();
    if (singlePowerSumIndex)
      {
        auto cached = hallLittlewoodPowerSumToCapitalColumnCache.find(singleIndex);
        if (cached == hallLittlewoodPowerSumToCapitalColumnCache.end())
          {
            CoeffMap generators =
                powerSumIndexToHallGeneratorMapViaLogarithmFormula(
                    singleIndex, false);
            if (error()) return zero();
            CoeffMap column = triangularReduceHallCapital(generators, false);
            if (error()) return zero();
            cached = hallLittlewoodPowerSumToCapitalColumnCache.emplace(
                singleIndex, std::move(column)).first;
          }
        ring_elem coefficient = poly->terms.front().coeff;
        if (omega &&
            (partitionWeight(singleIndex) - partitionLength(singleIndex)) % 2 != 0)
          coefficient = coefficientRing->negate(coefficient);
        addScaledCoeffMap(capitals, coefficient, cached->second);
      }
    else
      {
        VECTOR(SymmetricTerm) uncachedTerms;
        for (const auto& term : poly->terms)
          {
            Partition index;
            if (!powerSumIndexFromMonomial(term.monomial, index))
              {
                ERROR("expected a pure power-sum expression during basis conversion");
                return zero();
              }
            if (index.empty())
              {
                addCoeff(capitals, Partition{}, term.coeff);
                continue;
              }
            auto cached = hallLittlewoodPowerSumToCapitalColumnCache.find(index);
            if (cached == hallLittlewoodPowerSumToCapitalColumnCache.end())
              {
                uncachedTerms.push_back(term);
                continue;
              }
            ring_elem coefficient = term.coeff;
            if (omega &&
                (partitionWeight(index) - partitionLength(index)) % 2 != 0)
              coefficient = coefficientRing->negate(coefficient);
            addScaledCoeffMap(capitals, coefficient, cached->second);
          }
        if (!uncachedTerms.empty())
          {
            ring_elem uncachedExpression = fromTermVector(uncachedTerms, true);
            CoeffMap generators =
                powerSumsToHallGeneratorMapViaLogarithmFormula(
                    uncachedExpression, omega);
            if (error()) return zero();
            CoeffMap uncachedCapitals =
                triangularReduceHallCapital(generators, omega);
            if (error()) return zero();
            capitals = addCoeffMaps(capitals, uncachedCapitals);
          }
      }
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

ring_elem
SymmetricEngineRing::powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials(
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder) const
{
    const auto *poly = polyValue(f);
    const BasisKind targetKind = basisKindForId(targetBasisId);
    bool omega = targetKind == BasisKind::HallLittlewoodB ||
                 targetKind == BasisKind::HallLittlewoodPOmega;
    bool normalized = targetKind == BasisKind::HallLittlewoodP ||
                      targetKind == BasisKind::HallLittlewoodPOmega;
    CoeffMap result;
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index) ||
            index.size() > 1)
          {
            ERROR("single-cycle Green-polynomial conversion requires only p_n terms");
            return zero();
          }
        if (index.empty())
          {
            addCoeff(result, Partition{}, term.coeff);
            continue;
          }
        int n = index.front();
        ring_elem inputCoefficient = term.coeff;
        if (omega && n % 2 == 0)
          inputCoefficient = coefficientRing->negate(inputCoefficient);
        addScaledCoeffMap(
            result,
            inputCoefficient,
            hallLittlewoodSingleCycleGreenMap(n, normalized));
        if (error()) return zero();
      }
    return coeffMapToElement(result,
                             targetBasisId,
                             targetDisplay,
                             targetDisplayOrder,
                             false);
  }

ring_elem
SymmetricEngineRing::powerSumIndexToHallLittlewoodCapitalCoefficientViaGreenPolynomialDuality(
    const Partition& cycleType,
    const Partition& targetIndex) const
{
    Partition mu = normalizePartition(cycleType);
    Partition lambda = normalizePartition(targetIndex);
    if (partitionWeight(mu) != partitionWeight(lambda))
      return coefficientRing->zero();

    auto column = hallLittlewoodPowerSumToCapitalColumnCache.find(mu);
    if (column != hallLittlewoodPowerSumToCapitalColumnCache.end())
      {
        auto coefficient = column->second.find(lambda);
        return coefficient == column->second.end()
            ? coefficientRing->zero()
            : coefficient->second;
      }

    std::string cacheKey = partitionKey(lambda) + "|" + partitionKey(mu);
    auto cached = hallLittlewoodPowerSumToCapitalCoefficientCache.find(cacheKey);
    if (cached != hallLittlewoodPowerSumToCapitalCoefficientCache.end())
      return cached->second;

    ring_elem green = coefficientRing->zero();
    for (const auto& raisingTerm : raisingGeneratorMap(lambda))
      {
        mpz_class multiplicity = labeledCycleAssignmentsToCompleteFactors(
            mu, raisingTerm.first);
        if (multiplicity == 0) continue;
        green = coefficientRing->add(
            green,
            coefficientRing->mult(
                raisingTerm.second,
                coefficientRing->from_int(multiplicity.get_mpz_t())));
      }
    ring_elem result = coefficientRing->is_zero(green)
        ? coefficientRing->zero()
        : coefficientQuotient(green, hallLittlewoodCFactor(lambda));
    if (error()) return coefficientRing->zero();
    hallLittlewoodPowerSumToCapitalCoefficientCache[cacheKey] = result;
    return result;
  }

ring_elem
SymmetricEngineRing::powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality(
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder) const
{
    const auto *poly = polyValue(f);
    Partition cycleType;
    if (poly->terms.size() != 1 ||
        !powerSumIndexFromMonomial(poly->terms.front().monomial, cycleType) ||
        cycleType.empty())
      {
        ERROR("Green-polynomial duality conversion requires one p_mu term");
        return zero();
      }

    auto cached = hallLittlewoodPowerSumToCapitalColumnCache.find(cycleType);
    if (cached == hallLittlewoodPowerSumToCapitalColumnCache.end())
      {
        CoeffMap capitalColumn;
        for (const Partition& lambda : partitionsOf(partitionWeight(cycleType)))
          {
            ring_elem coefficient =
                powerSumIndexToHallLittlewoodCapitalCoefficientViaGreenPolynomialDuality(
                    cycleType, lambda);
            if (error()) return zero();
            addCoeff(capitalColumn, lambda, coefficient);
          }
        cached = hallLittlewoodPowerSumToCapitalColumnCache.emplace(
            cycleType, std::move(capitalColumn)).first;
      }

    const BasisKind targetKind = basisKindForId(targetBasisId);
    bool omega = targetKind == BasisKind::HallLittlewoodB ||
                 targetKind == BasisKind::HallLittlewoodPOmega;
    bool normalized = targetKind == BasisKind::HallLittlewoodP ||
                      targetKind == BasisKind::HallLittlewoodPOmega;
    ring_elem inputCoefficient = poly->terms.front().coeff;
    if (omega &&
        (partitionWeight(cycleType) - partitionLength(cycleType)) % 2 != 0)
      inputCoefficient = coefficientRing->negate(inputCoefficient);
    CoeffMap result;
    for (const auto& item : cached->second)
      {
        ring_elem coefficient = item.second;
        if (normalized)
          coefficient = coefficientRing->mult(
              coefficient, hallLittlewoodCFactor(item.first));
        addCoeff(result,
                 item.first,
                 coefficientRing->mult(inputCoefficient, coefficient));
      }
    return coeffMapToElement(result,
                             targetBasisId,
                             targetDisplay,
                             targetDisplayOrder,
                             false);
  }

ring_elem SymmetricEngineRing::powerSumsToHallLittlewoodViaTriangularReduction(
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder) const
{
    return powerSumsToHallLittlewoodCapitalViaTriangularReduction(f,
                                        targetBasisId,
                                        targetDisplay,
                                        targetDisplayOrder);
  }

bool SymmetricEngineRing::tryExpressionToHallLittlewoodViaTriangularReduction(ring_elem f,
                                      int targetBasisId,
                                      const std::string& targetDisplay,
                                      int targetDisplayOrder,
                                      ring_elem& result) const
{
    BasisKind generatorKind;
    bool omega = false;
    bool normalized = false;
    switch (basisKindForId(targetBasisId))
      {
      case BasisKind::HallLittlewoodQ:
        generatorKind = BasisKind::HallLittlewoodQGenerator;
        break;
      case BasisKind::HallLittlewoodP:
        generatorKind = BasisKind::HallLittlewoodQGenerator;
        normalized = true;
        break;
      case BasisKind::HallLittlewoodB:
        generatorKind = BasisKind::HallLittlewoodBGenerator;
        omega = true;
        break;
      case BasisKind::HallLittlewoodPOmega:
        generatorKind = BasisKind::HallLittlewoodBGenerator;
        omega = true;
        normalized = true;
        break;
      default:
        return false;
      }

    int generatorId = requiredBasisIdForKind(generatorKind);
    if (error()) return false;
    const std::string generatorDisplay = displayForBasis(generatorId);
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

// ============================================================================
// Generic Basis-Element Conversion
// ============================================================================
// Per-basis-element conversion feeds monomial, expression, and termwise fallback routes.

SymmetricEngineRing::BasisElementToPowerSumsRoute
SymmetricEngineRing::selectBasisElementToPowerSumsRoute(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    BasisKind kind = basisKindForId(atomBasisIdAt(monomial, pos));
    if (atomIsSkewAt(monomial, pos))
      {
        switch (kind)
          {
            case BasisKind::Schur:
              return BasisElementToPowerSumsRoute::ViaSkewSchurJacobiTrudiComplete;
            case BasisKind::SchurOmega:
              return BasisElementToPowerSumsRoute::ViaSkewSchurOmegaJacobiTrudiElementary;
            case BasisKind::HallLittlewoodQ:
            case BasisKind::HallLittlewoodB:
            case BasisKind::HallLittlewoodP:
            case BasisKind::HallLittlewoodPOmega:
              return BasisElementToPowerSumsRoute::ViaSkewHallLittlewood;
            case BasisKind::Custom:
            case BasisKind::PowerSum:
            case BasisKind::Complete:
            case BasisKind::Elementary:
            case BasisKind::Monomial:
            case BasisKind::Forgotten:
            case BasisKind::HallLittlewoodQGenerator:
            case BasisKind::HallLittlewoodBGenerator:
              return BasisElementToPowerSumsRoute::NoApplicableRoute;
          }
      }
    switch (kind)
      {
        case BasisKind::PowerSum:
          return BasisElementToPowerSumsRoute::AlreadyPowerSums;
        case BasisKind::Complete:
          return BasisElementToPowerSumsRoute::ViaCompleteClassicalFormula;
        case BasisKind::Elementary:
          return BasisElementToPowerSumsRoute::ViaElementaryClassicalFormula;
        case BasisKind::HallLittlewoodQGenerator:
        case BasisKind::HallLittlewoodBGenerator:
          return BasisElementToPowerSumsRoute::ViaHallLittlewoodGeneratorClassicalFormula;
        case BasisKind::Schur:
          return BasisElementToPowerSumsRoute::ViaSchurCharacters;
        case BasisKind::SchurOmega:
          return BasisElementToPowerSumsRoute::ViaSchurOmegaCharacters;
        case BasisKind::Monomial:
          return BasisElementToPowerSumsRoute::ViaMonomialTransition;
        case BasisKind::Forgotten:
          return BasisElementToPowerSumsRoute::ViaForgottenTransition;
        case BasisKind::HallLittlewoodQ:
        case BasisKind::HallLittlewoodB:
          return BasisElementToPowerSumsRoute::ViaHallLittlewoodRaisingOperators;
        case BasisKind::HallLittlewoodP:
        case BasisKind::HallLittlewoodPOmega:
          return BasisElementToPowerSumsRoute::ViaHallLittlewoodCapitalNormalization;
        case BasisKind::Custom:
          return BasisElementToPowerSumsRoute::NoApplicableRoute;
      }
    return BasisElementToPowerSumsRoute::NoApplicableRoute;
  }

const char *SymmetricEngineRing::basisElementToPowerSumsRouteName(
    BasisElementToPowerSumsRoute route) const
{
    switch (route)
      {
        case BasisElementToPowerSumsRoute::AlreadyPowerSums: return "p:identity";
        case BasisElementToPowerSumsRoute::ViaCompleteClassicalFormula: return "h->p:classical-formula";
        case BasisElementToPowerSumsRoute::ViaElementaryClassicalFormula: return "e->p:classical-formula";
        case BasisElementToPowerSumsRoute::ViaHallLittlewoodGeneratorClassicalFormula: return "q/b->p:classical-formula";
        case BasisElementToPowerSumsRoute::ViaSchurCharacters: return "S->p:characters";
        case BasisElementToPowerSumsRoute::ViaSchurOmegaCharacters: return "Somega->p:characters";
        case BasisElementToPowerSumsRoute::ViaMonomialTransition: return "m->p:transition";
        case BasisElementToPowerSumsRoute::ViaForgottenTransition: return "ff->p:transition";
        case BasisElementToPowerSumsRoute::ViaHallLittlewoodRaisingOperators: return "Q/B->p:raising-operators";
        case BasisElementToPowerSumsRoute::ViaHallLittlewoodCapitalNormalization: return "P/Pomega->Q/B->p";
        case BasisElementToPowerSumsRoute::ViaSkewSchurJacobiTrudiComplete: return "skew-S->h->p:Jacobi-Trudi";
        case BasisElementToPowerSumsRoute::ViaSkewSchurOmegaJacobiTrudiElementary: return "skew-Somega->e->p:Jacobi-Trudi";
        case BasisElementToPowerSumsRoute::ViaSkewHallLittlewood: return "skew-Hall-Littlewood->p";
        case BasisElementToPowerSumsRoute::NoApplicableRoute: return "not-applicable";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceBasisElementToPowerSumsSelection(
    BasisElementToPowerSumsRoute route,
    const std::string& sourceDisplay) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::fprintf(stderr,
                 "SymmetricRings basis-element-to-power-sums: source=%s route=%s\n",
                 sourceDisplay.c_str(),
                 basisElementToPowerSumsRouteName(route));
  }

ring_elem SymmetricEngineRing::executeBasisElementToPowerSumsRoute(
    BasisElementToPowerSumsRoute route,
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    std::string display = displayForBasis(atomBasisIdAt(monomial, pos));
    if (route == BasisElementToPowerSumsRoute::ViaSkewSchurJacobiTrudiComplete ||
        route == BasisElementToPowerSumsRoute::ViaSkewSchurOmegaJacobiTrudiElementary)
      {
        BasisKind generatorKind =
            route == BasisElementToPowerSumsRoute::ViaSkewSchurJacobiTrudiComplete
                ? BasisKind::Complete : BasisKind::Elementary;
        int generatorId = requiredBasisIdForKind(generatorKind);
        if (error()) return zero();
        std::string generatorDisplay = displayForBasis(generatorId);
        return expressionToPowerSumsViaBasisElementRoutes(
            jacobiTrudi(basisElementOuterIndex(monomial, pos),
                        basisElementInnerIndex(monomial, pos),
                        generatorId,
                        generatorDisplay,
                        basisOrderForId(generatorId),
                        isMultiplicativeBasis(generatorId)));
      }
    if (route == BasisElementToPowerSumsRoute::ViaSkewHallLittlewood)
      return skewHallLittlewoodToPowerSums(
                                           basisElementOuterIndex(monomial, pos),
                                           basisElementInnerIndex(monomial, pos),
                                           basisKindForId(
                                               atomBasisIdAt(monomial, pos)));

    Partition index = basisElementIndex(monomial, pos);
    if (route == BasisElementToPowerSumsRoute::AlreadyPowerSums)
      return basisElementFromIndex(powerSumBasisId,
                                   displayForBasis(powerSumBasisId),
                                   basisOrderForId(powerSumBasisId),
                                   true,
                                   index);
    if (route == BasisElementToPowerSumsRoute::ViaCompleteClassicalFormula ||
        route == BasisElementToPowerSumsRoute::ViaElementaryClassicalFormula)
      {
        ring_elem result = one();
        for (int part : index)
          {
            ring_elem factor =
                route == BasisElementToPowerSumsRoute::ViaCompleteClassicalFormula
                    ? completePartToPowerSumsViaClassicalFormula(part)
                    : elementaryPartToPowerSumsViaClassicalFormula(part);
            result = mult(result, factor);
          }
        return result;
      }
    if (route == BasisElementToPowerSumsRoute::ViaHallLittlewoodGeneratorClassicalFormula)
      {
        ring_elem result = one();
        for (int part : index)
          result = mult(result,
                        hallLittlewoodGeneratorPartToPowerSumsViaClassicalFormula(
                            part,
                            basisKindForId(atomBasisIdAt(monomial, pos)) ==
                                BasisKind::HallLittlewoodBGenerator));
        return result;
      }
    if (route == BasisElementToPowerSumsRoute::ViaSchurCharacters)
      return schurLikeToPowerSumsViaCharacters(index, false);
    if (route == BasisElementToPowerSumsRoute::ViaSchurOmegaCharacters)
      return schurLikeToPowerSumsViaCharacters(index, true);
    if (route == BasisElementToPowerSumsRoute::ViaMonomialTransition ||
        route == BasisElementToPowerSumsRoute::ViaForgottenTransition)
      return monomialToPowerSumsViaTransitionMatrix(
          index, route == BasisElementToPowerSumsRoute::ViaForgottenTransition);
    if (route == BasisElementToPowerSumsRoute::ViaHallLittlewoodRaisingOperators)
      return hallLittlewoodCapitalToPowerSumsViaRaisingOperators(
          index,
          basisKindForId(atomBasisIdAt(monomial, pos)) ==
              BasisKind::HallLittlewoodB);
    if (route == BasisElementToPowerSumsRoute::ViaHallLittlewoodCapitalNormalization)
      return hallLittlewoodNormalizedToPowerSumsViaCapitalNormalization(
          index,
          basisKindForId(atomBasisIdAt(monomial, pos)) ==
              BasisKind::HallLittlewoodPOmega);
    if (atomIsSkewAt(monomial, pos))
      ERROR("basis conversion for skew ", display.c_str(), " basis elements is not implemented yet");
    else
      ERROR("basis conversion to power sums is not implemented for basis ", display.c_str());
    return zero();
  }

ring_elem SymmetricEngineRing::basisElementToPowerSumsDispatch(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    BasisElementToPowerSumsRoute route =
        selectBasisElementToPowerSumsRoute(monomial, pos);
    std::string sourceDisplay = displayForBasis(atomBasisIdAt(monomial, pos));
    traceBasisElementToPowerSumsSelection(route, sourceDisplay);
    return executeBasisElementToPowerSumsRoute(route, monomial, pos);
  }

ring_elem SymmetricEngineRing::monomialToPowerSumsViaBasisElementRoutes(const SymmetricMonomial& monomial) const
{
    ring_elem result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        result = mult(result, basisElementToPowerSumsDispatch(monomial, pos));
        if (error()) return zero();
        pos += atomLengthAt(monomial, pos);
      }
    return result;
  }

ring_elem SymmetricEngineRing::expressionToPowerSumsViaBasisElementRoutes(ring_elem f) const
{
    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        ring_elem converted =
            monomialToPowerSumsViaBasisElementRoutes(term.monomial);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, converted));
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumIndexToTargetViaTermwiseKernel(const Partition& index,
                                     const std::string& targetDisplay,
                                     int targetBasisId,
                                     int targetDisplayOrder,
                                     bool targetIsMultiplicative) const
{
    BasisKind targetKind = basisKindForId(targetBasisId);
    switch (targetKind)
      {
        case BasisKind::PowerSum:
          return basisElementFromIndex(targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       targetIsMultiplicative,
                                       index);
        case BasisKind::Complete:
          return coeffMapToElement(
              powerSumIndexToCompleteMapViaLogarithmFormula(index),
              targetBasisId, targetDisplay, targetDisplayOrder, true);
        case BasisKind::Elementary:
          return coeffMapToElement(
              powerSumIndexToElementaryMapViaLogarithmFormula(index),
              targetBasisId, targetDisplay, targetDisplayOrder, true);
        case BasisKind::HallLittlewoodQGenerator:
        case BasisKind::HallLittlewoodBGenerator:
          return coeffMapToElement(
              powerSumIndexToHallGeneratorMapViaLogarithmFormula(
                  index,
                  targetKind == BasisKind::HallLittlewoodBGenerator),
              targetBasisId, targetDisplay, targetDisplayOrder, true);
        case BasisKind::Monomial:
        case BasisKind::Forgotten:
          return powerSumIndexToMonomialViaTransitionMatrix(
              index, targetBasisId, targetDisplay, targetDisplayOrder,
              targetKind == BasisKind::Forgotten);
        case BasisKind::Schur:
          return powerSumIndexToSchurLikeViaCharacters(
              index, targetBasisId, targetDisplayOrder, targetDisplay, 1);
        case BasisKind::SchurOmega:
          {
            long sign =
                ((partitionWeight(index) - partitionLength(index)) % 2 == 0)
                ? 1 : -1;
            return powerSumIndexToSchurLikeViaCharacters(
                index, targetBasisId, targetDisplayOrder, targetDisplay, sign);
          }
        case BasisKind::Custom:
        case BasisKind::HallLittlewoodQ:
        case BasisKind::HallLittlewoodB:
        case BasisKind::HallLittlewoodP:
        case BasisKind::HallLittlewoodPOmega:
          ERROR("basis conversion from power sums is not implemented for basis ",
                targetDisplay.c_str());
          return zero();
      }
    ERROR("unknown target basis kind");
    return zero();
  }

ring_elem SymmetricEngineRing::powerSumsToTargetViaTermwiseConversion(ring_elem f,
                                      int targetBasisId,
                                      const std::string& targetDisplay,
                                      int targetDisplayOrder,
                                      bool targetIsMultiplicative) const
{
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
        ring_elem converted = powerSumIndexToTargetViaTermwiseKernel(index,
                                                       targetDisplay,
                                                       targetBasisId,
                                                       targetDisplayOrder,
                                                       targetIsMultiplicative);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, converted));
      }
    return result;
  }

bool SymmetricEngineRing::tryBasisElementToTarget(const SymmetricMonomial& monomial,
                          size_t pos,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          bool targetIsMultiplicative,
                          ring_elem& result) const
{
    int basisId = atomBasisIdAt(monomial, pos);
    std::string display = displayForBasis(basisId);
    const BasisKind sourceKind = basisKindForId(basisId);
    const BasisKind targetKind = basisKindForId(targetBasisId);
    Partition index = basisElementIndex(monomial, pos);

    if (!atomIsSkewAt(monomial, pos) && basisId == targetBasisId)
      {
        result = basisElementFromIndex(targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       targetIsMultiplicative,
                                       index);
        return true;
      }

    if (targetKind == BasisKind::Complete && sourceKind == BasisKind::Schur)
      {
        Partition outer = atomIsSkewAt(monomial, pos) ? basisElementOuterIndex(monomial, pos)
                                                      : index;
        Partition inner = atomIsSkewAt(monomial, pos) ? basisElementInnerIndex(monomial, pos)
                                                      : Partition{};
        result = jacobiTrudi(outer,
                             inner,
                             targetBasisId,
                             targetDisplay,
                             targetDisplayOrder,
                             targetIsMultiplicative);
        return true;
      }

    if (targetKind == BasisKind::Elementary &&
        sourceKind == BasisKind::SchurOmega)
      {
        Partition outer = atomIsSkewAt(monomial, pos) ? basisElementOuterIndex(monomial, pos)
                                                      : index;
        Partition inner = atomIsSkewAt(monomial, pos) ? basisElementInnerIndex(monomial, pos)
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
    if (sourceKind == BasisKind::Custom) return false;
    ring_elem inPowerSums = basisElementToPowerSumsDispatch(monomial, pos);
    if (error()) return false;
    int pBasisId = requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return false;
    ConversionInput input{
        inPowerSums,
        inferConversionGuarantees(inPowerSums, targetBasisId),
        0};
    input.guarantees.pureBasis = pBasisId;
    input.guarantees.expandedBasis = pBasisId;
    input.guarantees = strengthenConversionGuarantees(
        std::move(input.guarantees), targetBasisId);
    result = powerSumsToTargetDispatch(input,
                                       pBasisId,
                                       targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       targetIsMultiplicative);
    return !error();
  }

bool SymmetricEngineRing::tryMonomialToTarget(const SymmetricMonomial& monomial,
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
        if (!tryBasisElementToTarget(monomial,
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

SymmetricEngineRing::ExpressionToTargetRoute
SymmetricEngineRing::selectExpressionToTargetRoute(
    int targetBasisId,
    bool targetIsMultiplicative) const
{
    switch (basisKindForId(targetBasisId))
      {
      case BasisKind::Schur:
      case BasisKind::SchurOmega:
        return ExpressionToTargetRoute::ViaSchurTriangularReduction;
      case BasisKind::HallLittlewoodQ:
      case BasisKind::HallLittlewoodP:
      case BasisKind::HallLittlewoodB:
      case BasisKind::HallLittlewoodPOmega:
        return ExpressionToTargetRoute::ViaHallLittlewoodTriangularReduction;
      case BasisKind::Complete:
      case BasisKind::Elementary:
        return ExpressionToTargetRoute::ViaFactorwiseConversion;
      default:
        return targetIsMultiplicative
            ? ExpressionToTargetRoute::ViaFactorwiseConversion
            : ExpressionToTargetRoute::NoApplicableRoute;
      }
  }

const char *SymmetricEngineRing::expressionToTargetRouteName(
    ExpressionToTargetRoute route) const
{
    switch (route)
      {
        case ExpressionToTargetRoute::ViaSchurTriangularReduction:
          return "Schur-triangular-reduction";
        case ExpressionToTargetRoute::ViaHallLittlewoodTriangularReduction:
          return "Hall-Littlewood-triangular-reduction";
        case ExpressionToTargetRoute::ViaFactorwiseConversion:
          return "factorwise-conversion";
        case ExpressionToTargetRoute::NoApplicableRoute:
          return "not-applicable";
      }
    return "unknown";
  }

void SymmetricEngineRing::traceExpressionToTargetSelection(
    ExpressionToTargetRoute route,
    const std::string& targetDisplay) const
{
    if (std::getenv("M2_SYMMETRIC_RINGS_TRACE_CONVERSION") == nullptr) return;
    std::fprintf(stderr,
                 "SymmetricRings expression-target: target=%s route=%s\n",
                 targetDisplay.c_str(),
                 expressionToTargetRouteName(route));
  }

bool SymmetricEngineRing::executeExpressionToTargetRoute(
    ExpressionToTargetRoute route,
    ring_elem f,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder,
    bool targetIsMultiplicative,
    ring_elem& result) const
{
    if (route == ExpressionToTargetRoute::ViaSchurTriangularReduction)
      return tryExpressionToSchurViaTriangularReduction(
          f, targetBasisId, targetDisplay, targetDisplayOrder, result);
    if (route == ExpressionToTargetRoute::ViaHallLittlewoodTriangularReduction)
      return tryExpressionToHallLittlewoodViaTriangularReduction(
          f, targetBasisId, targetDisplay, targetDisplayOrder, result);
    if (route != ExpressionToTargetRoute::ViaFactorwiseConversion)
      return false;
    result = zero();
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        ring_elem converted;
        if (!tryMonomialToTarget(term.monomial,
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

bool SymmetricEngineRing::tryExpressionToTarget(ring_elem f,
                             int targetBasisId,
                             const std::string& targetDisplay,
                             int targetDisplayOrder,
                             bool targetIsMultiplicative,
                             ring_elem& result) const
{
    ExpressionToTargetRoute route = selectExpressionToTargetRoute(
        targetBasisId, targetIsMultiplicative);
    traceExpressionToTargetSelection(route, targetDisplay);
    return executeExpressionToTargetRoute(route,
                                           f,
                                           targetBasisId,
                                           targetDisplay,
                                           targetDisplayOrder,
                                           targetIsMultiplicative,
                                           result);
  }

// ============================================================================
// Straightening
// ============================================================================
// Composition-indexed basis elements are rewritten into canonical basis expansions.

Partition SymmetricEngineRing::replaceAdjacentPair(const Partition& alpha,
                                size_t pos,
                                int first,
                                int second) const
{
    Partition result = alpha;
    result[pos] = first;
    result[pos + 1] = second;
    return result;
  }

ring_elem SymmetricEngineRing::straightenSchurBasisElement(const Partition& alpha,
                                int basisId) const
{
    auto straightened = straightenSchurIndex(alpha);
    if (straightened.first == 0) return zero();
    ring_elem term = basisElementFromIndex(basisId,
                                           displayForBasis(basisId),
                                           basisOrderForId(basisId),
                                           isMultiplicativeBasis(basisId),
                                           straightened.second);
    if (straightened.first < 0) term = negate(term);
    return term;
  }

ring_elem SymmetricEngineRing::straightenHallCapitalBasisElement(const Partition& alpha,
                                      int basisId) const
{
    Partition trimmed = trimTrailingZerosPartition(alpha);
    if (trimmed.empty()) return one();
    size_t bad = trimmed.size();
    for (size_t i = 0; i + 1 < trimmed.size(); ++i)
      if (trimmed[i] < trimmed[i + 1])
        {
          bad = i;
          break;
        }
    if (bad == trimmed.size())
      return basisElementFromIndex(basisId,
                                   displayForBasis(basisId),
                                   basisOrderForId(basisId),
                                   isMultiplicativeBasis(basisId),
                                   trimmed);

    int s = trimmed[bad];
    int r = trimmed[bad + 1];
    int diff = r - s;
    int top = diff / 2;
    ring_elem result =
        scaled(hallLittlewoodParameter,
               straightenHallCapitalBasisElement(replaceAdjacentPair(trimmed, bad, r, s),
                                          basisId));
    for (int i = 1; i <= top; ++i)
      {
        ring_elem coeff;
        if (diff % 2 == 0 && i == top)
          coeff = coefficientRing->subtract(
              coefficientRing->power(hallLittlewoodParameter, i),
              coefficientRing->power(hallLittlewoodParameter, i - 1));
        else
          coeff = coefficientRing->subtract(
              coefficientRing->power(hallLittlewoodParameter, i + 1),
              coefficientRing->power(hallLittlewoodParameter, i - 1));
        result = add(result,
                     scaled(coeff,
                            straightenHallCapitalBasisElement(
                                replaceAdjacentPair(trimmed, bad, r - i, s + i),
                                basisId)));
      }
    return result;
  }

ring_elem SymmetricEngineRing::straightenBasisElement(const SymmetricMonomial& monomial, size_t pos) const
{
    const int basisId = atomBasisIdAt(monomial, pos);
    const BasisKind basisKind = basisKindForId(basisId);
    std::string display = displayForBasis(basisId);
    if (atomIsSkewAt(monomial, pos))
      {
        auto *poly = new SymmetricRingPoly;
        poly->terms.push_back({coefficientRing->one(), monomialFromKey(atomBlockAt(monomial, pos))});
        return makePolyValue(poly);
      }
    Partition index = basisElementIndex(monomial, pos);
    switch (basisKind)
      {
      case BasisKind::Schur:
      case BasisKind::SchurOmega:
        return straightenSchurBasisElement(index, basisId);
      case BasisKind::HallLittlewoodQ:
      case BasisKind::HallLittlewoodB:
        return straightenHallCapitalBasisElement(index, basisId);
      case BasisKind::HallLittlewoodP:
      case BasisKind::HallLittlewoodPOmega:
        {
        if (isPartitionIndex(index))
          return basisElementFromIndex(basisId,
                                       display,
                                       atomOrderAt(monomial, pos),
                                       false,
                                       index);
        const BasisKind capitalKind = basisKind == BasisKind::HallLittlewoodP
                                          ? BasisKind::HallLittlewoodQ
                                          : BasisKind::HallLittlewoodB;
        const int capitalId = requiredBasisIdForKind(capitalKind);
        if (error()) return zero();
        const std::string capitalDisplay = displayForBasis(capitalId);
        ring_elem straightCapital =
            straightenHallCapitalBasisElement(index, capitalId);
        if (error()) return zero();
        ring_elem normalizedCapital = scaled(
            coefficientQuotient(coefficientRing->one(),
                                hallLittlewoodCFactor(index)),
            straightCapital);
        if (error()) return zero();
        return powerSumsToHallLittlewoodCapitalViaTriangularReduction(
            expressionToPowerSumsViaBasisElementRoutes(normalizedCapital),
            basisId,
            display,
            basisOrderForId(basisId));
        }
      default:
        break;
      }
    return basisElementFromIndex(basisId,
                                 display,
                                 atomOrderAt(monomial, pos),
                                 isMultiplicativeBasis(basisId),
                                 index);
  }

ring_elem SymmetricEngineRing::straightenMonomial(const SymmetricMonomial& monomial) const
{
    ring_elem result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        ring_elem factor = straightenBasisElement(monomial, pos);
        if (error()) return zero();
        result = mult(result, factor);
        pos += atomLengthAt(monomial, pos);
      }
    return result;
  }

ring_elem SymmetricEngineRing::straightenElement(ring_elem f) const
{
    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        ring_elem straightened = straightenMonomial(term.monomial);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, straightened));
      }
    return result;
  }

ring_elem SymmetricEngineRing::straighten(ring_elem f) const
{
    return straightenElement(f);
  }

} // namespace symmetric_rings
