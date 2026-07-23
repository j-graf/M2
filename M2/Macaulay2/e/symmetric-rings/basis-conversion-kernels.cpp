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
#include <iterator>
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
    const Partition& factorDegrees,
    size_t maxStates,
    bool& limitExceeded)
{
    using State = std::pair<size_t, Partition>;
    std::map<State, mpz_class> memo;
    std::function<mpz_class(size_t, Partition&)> count =
        [&](size_t position, Partition& remaining) -> mpz_class {
          if (limitExceeded) return 0;
          if (position == cycleType.size())
            {
              for (int value : remaining)
                if (value != 0) return 0;
              return 1;
            }
          State state{position, remaining};
          auto cached = memo.find(state);
          if (cached != memo.end()) return cached->second;
          if (memo.size() >= maxStates)
            {
              limitExceeded = true;
              return 0;
            }
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
    limitExceeded = false;
    return count(0, remaining);
}

} // namespace

// ============================================================================
// Shared Conversion Utilities
// ============================================================================
// Coefficient maps, basis metadata, stored basis-element inspection, and the
// policy-free assembler for canonical built-in expansions to power sums.

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
    ERROR(
        "coefficient division failed during basis conversion; use a coefficient "
        "ring where the required denominators are invertible, for example "
        "frac(QQ[t]) instead of QQ[t]");
    return coefficientRing->zero();
  }

void SymmetricEngineRing::addNormalizedCoeff(CoeffMap& target,
                                             Partition index,
                                             ring_elem coeff) const
{
    if (coefficientRing->is_zero(coeff)) return;
    auto existing = target.find(index);
    if (existing == target.end())
      {
        if (target.size() >= computationLimits.maxGeneratedTerms)
          throw exc::engine_error(
              "generated-term limit exceeded during basis conversion; increase "
              "MaxGeneratedTerms in the symmetricRing ComputationLimits option");
        target.emplace(std::move(index), coeff);
        return;
      }
    ring_elem sum = coefficientRing->add(existing->second, coeff);
    if (coefficientRing->is_zero(sum))
      target.erase(existing);
    else
      existing->second = sum;
  }

void SymmetricEngineRing::addCoeff(CoeffMap& target,
                                   const Partition& index,
                                   ring_elem coeff) const
{
    addNormalizedCoeff(target, normalizePartition(index), coeff);
  }

void SymmetricEngineRing::addScaledCoeffMap(CoeffMap& target,
                                             ring_elem coeff,
                                             const CoeffMap& source) const
{
    if (coefficientRing->is_zero(coeff)) return;
    for (const auto& item : source)
      addNormalizedCoeff(target,
                         item.first,
                         coefficientRing->mult(coeff, item.second));
  }

CoeffMap SymmetricEngineRing::addCoeffMaps(const CoeffMap& a, const CoeffMap& b) const
{
    CoeffMap result = a;
    for (const auto& item : b)
      addNormalizedCoeff(result, item.first, item.second);
    return result;
  }

CoeffMap SymmetricEngineRing::multiplyCoeffMaps(const CoeffMap& a, const CoeffMap& b) const
{
    CoeffMap result;
    for (const auto& left : a)
      for (const auto& right : b)
        {
          Partition index;
          index.reserve(left.first.size() + right.first.size());
          std::merge(left.first.begin(),
                     left.first.end(),
                     right.first.begin(),
                     right.first.end(),
                     std::back_inserter(index),
                     std::greater<int>());
          ring_elem coeff = coefficientRing->mult(left.second, right.second);
          addNormalizedCoeff(result, std::move(index), coeff);
        }
    return result;
  }

CoeffMap SymmetricEngineRing::oneCoeffMap() const
{
    return CoeffMap{{Partition{}, coefficientRing->one()}};
  }

ring_elem SymmetricEngineRing::cachedInteger(long n) const
{
    auto found = smallIntegerCoeffCache.find(n);
    if (found != smallIntegerCoeffCache.end()) return found->second;
    ring_elem value = coefficientRing->from_long(n);
    requireCacheEntryCapacity("small-integer cache");
    smallIntegerCoeffCache[n] = value;
    return value;
  }

ring_elem SymmetricEngineRing::cachedInteger(const mpz_class& n) const
{
    if (mpz_fits_slong_p(n.get_mpz_t()))
      return cachedInteger(n.get_si());
    return coefficientRing->from_int(n.get_mpz_t());
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
    requireBasis(targetBasisId);
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
    return basisElementFromIndex(id, index);
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
    return basisElementFromIndex(registeredPowerSumBasisId(), normalized);
  }

ring_elem
SymmetricEngineRing::canonicalExpressionToPowerSumsViaBasisFormulas(
    ring_elem f,
    BasisKind sourceKind) const
{
    // This is deterministic mathematical dispatch by basis family, not
    // performance policy. Callers must already guarantee one canonical source
    // basis element per nonscalar term.
    ring_elem result = zero();
    for (const auto& term : polyValue(f)->terms)
      {
        if (term.monomial.data.empty())
          {
            result = add(result, fromCoeff(term.coeff));
            continue;
          }
        if (atomLengthAt(term.monomial, 0) != term.monomial.data.size() ||
            atomIsSkewAt(term.monomial, 0) ||
            basisKindForId(atomBasisIdAt(term.monomial, 0)) != sourceKind)
          {
            ERROR("a source-to-power-sums basis formula received "
                  "noncanonical input");
            return zero();
          }

        Partition index = basisElementIndex(term.monomial, 0);
        ring_elem convertedElement;
        switch (sourceKind)
          {
            case BasisKind::PowerSum:
              convertedElement = basisElementFromIndex(
                  atomBasisIdAt(term.monomial, 0), index);
              break;
            case BasisKind::Complete:
            case BasisKind::Elementary:
              convertedElement = one();
              for (int part : index)
                {
                  ring_elem factor =
                      sourceKind == BasisKind::Complete
                          ? completePartToPowerSumsViaClassicalFormula(part)
                          : elementaryPartToPowerSumsViaClassicalFormula(part);
                  convertedElement = mult(convertedElement, factor);
                }
              break;
            case BasisKind::HallLittlewoodQGenerator:
            case BasisKind::HallLittlewoodBGenerator:
              convertedElement = one();
              for (int part : index)
                convertedElement = mult(
                    convertedElement,
                    hallLittlewoodGeneratorPartToPowerSumsViaClassicalFormula(
                        part,
                        sourceKind ==
                            BasisKind::HallLittlewoodBGenerator));
              break;
            case BasisKind::Schur:
              convertedElement =
                  schurLikeToPowerSumsViaCharacters(index, false);
              break;
            case BasisKind::SchurOmega:
              convertedElement =
                  schurLikeToPowerSumsViaCharacters(index, true);
              break;
            case BasisKind::Monomial:
            case BasisKind::Forgotten:
              convertedElement =
                  monomialToPowerSumsViaTransitionMatrix(
                      index, sourceKind == BasisKind::Forgotten);
              break;
            case BasisKind::HallLittlewoodQ:
            case BasisKind::HallLittlewoodB:
              convertedElement =
                  hallLittlewoodCapitalToPowerSumsViaRaisingOperators(
                      index, sourceKind == BasisKind::HallLittlewoodB);
              break;
            case BasisKind::HallLittlewoodP:
            case BasisKind::HallLittlewoodPOmega:
              convertedElement =
                  hallLittlewoodNormalizedToPowerSumsViaCapitalNormalization(
                      index,
                      sourceKind == BasisKind::HallLittlewoodPOmega);
              break;
            case BasisKind::Custom:
              ERROR("a source-to-power-sums basis formula does not support "
                    "custom bases");
              return zero();
          }
        if (error()) return zero();
        result = add(
            result, scaled(term.coeff, convertedElement));
        if (error()) return zero();
      }
    return result;
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
                        int basisId) const
{
    requireBasis(basisId);
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
            matrix[i][j] = basisPartElement(basisId, degree);
          }
      }

    size_t limit = static_cast<size_t>(1) << n;
    if (!determinantStatesWithinLimit(limit, "Jacobi-Trudi determinant"))
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
            ring_elem term = mult(dp[mask], matrix[row][col]);
            if (selectedGreaterThan(mask, col, n) % 2 == 1) term = negate(term);
            size_t next = mask | bit;
            dp[next] = add(dp[next], term);
          }
      }

    ring_elem result = dp[limit - 1];
    requireCacheEntryCapacity("Jacobi-Trudi cache");
    cache[cacheKey] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::jacobiTrudiBasis(int basisId,
                             const Partition& outer,
                             const Partition& inner) const
{
    return jacobiTrudi(outer, inner, basisId);
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
    for (const auto& mu : partitionsOfWithinLimits(n, "power-sum conversion"))
      {
        ring_elem coeff = rationalCoefficient(1, zValue(mu));
        if (error()) return zero();
        addCoeff(coefficients, mu, coeff);
      }
    ring_elem result = coeffMapToElement(
        coefficients,
        registeredPowerSumBasisId(),
        displayForBasis(registeredPowerSumBasisId()),
        basisOrderForId(registeredPowerSumBasisId()),
        true);
    requireCacheEntryCapacity("complete-to-power-sum cache");
    completeToPowerSumsCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

ring_elem SymmetricEngineRing::elementaryPartToPowerSumsViaClassicalFormula(int n) const
{
    if (n == 0) return one();
    auto cached = elementaryToPowerSumsCache.find(n);
    if (cached != elementaryToPowerSumsCache.end()) return copyPolyValue(polyValue(cached->second));
    CoeffMap coefficients;
    for (const auto& mu : partitionsOfWithinLimits(n, "power-sum conversion"))
      {
        long sign = ((n - static_cast<int>(mu.size())) % 2 == 0) ? 1 : -1;
        ring_elem coeff = rationalCoefficient(sign, zValue(mu));
        if (error()) return zero();
        addCoeff(coefficients, mu, coeff);
      }
    ring_elem result = coeffMapToElement(
        coefficients,
        registeredPowerSumBasisId(),
        displayForBasis(registeredPowerSumBasisId()),
        basisOrderForId(registeredPowerSumBasisId()),
        true);
    requireCacheEntryCapacity("elementary-to-power-sum cache");
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
    for (const auto& lambda : partitionsOfWithinLimits(n, "monomial conversion"))
      {
        ring_elem logarithmCoeff = powerSumLogarithmCoefficient(lambda);
        if (error()) return CoeffMap{};
        addCoeff(result,
                 lambda,
                 coefficientRing->mult(common, logarithmCoeff));
      }
    return result;
  }

CoeffMap
SymmetricEngineRing::powerSumPartToIntegralGeneratorMapViaLogarithmFormula(
    int n,
    int commonSign) const
{
    CoeffMap result;
    mpz_t coefficient;
    mpz_t divisor;
    mpz_init(coefficient);
    mpz_init(divisor);
    for (const auto& lambda : partitionsOfWithinLimits(n, "forgotten conversion"))
      {
        // [h_lambda] p_n = (-1)^(length(lambda)-1)
        //     n (length(lambda)-1)! / product_i multiplicity_i(lambda)!.
        // This is always integral, so avoid constructing a rational only to
        // multiply it by n and cancel its denominator in the coefficient ring.
        mpz_fac_ui(coefficient, lambda.empty() ? 0 : lambda.size() - 1);
        mpz_mul_ui(coefficient, coefficient, static_cast<unsigned long>(n));
        if (((lambda.size() - 1) % 2 == 1) != (commonSign < 0))
          mpz_neg(coefficient, coefficient);

        for (size_t first = 0; first < lambda.size();)
          {
            size_t last = first + 1;
            while (last < lambda.size() && lambda[last] == lambda[first]) ++last;
            mpz_fac_ui(divisor, last - first);
            mpz_divexact(coefficient, coefficient, divisor);
            first = last;
          }
        addCoeff(result, lambda, coefficientRing->from_int(coefficient));
      }
    mpz_clear(divisor);
    mpz_clear(coefficient);
    return result;
  }

const CoeffMap&
SymmetricEngineRing::powerSumPartToCompleteMapViaLogarithmFormula(int n) const
{
    auto cached = powerSumToCompleteMapCache.find(n);
    if (cached != powerSumToCompleteMapCache.end()) return cached->second;
    CoeffMap result = n == 0
        ? oneCoeffMap()
        : powerSumPartToIntegralGeneratorMapViaLogarithmFormula(n, 1);
    requireCacheEntryCapacity("power-sum-to-complete cache");
    auto inserted = powerSumToCompleteMapCache.emplace(n, std::move(result));
    return inserted.first->second;
}

const CoeffMap&
SymmetricEngineRing::powerSumPartToElementaryMapViaLogarithmFormula(int n) const
{
    auto cached = powerSumToElementaryMapCache.find(n);
    if (cached != powerSumToElementaryMapCache.end()) return cached->second;
    CoeffMap result = n == 0
        ? oneCoeffMap()
        : powerSumPartToIntegralGeneratorMapViaLogarithmFormula(
              n, n % 2 == 0 ? -1 : 1);
    requireCacheEntryCapacity("power-sum-to-elementary cache");
    auto inserted = powerSumToElementaryMapCache.emplace(n, std::move(result));
    return inserted.first->second;
}

CoeffMap SymmetricEngineRing::powerSumIndexToCompleteMapViaLogarithmFormula(
    const Partition& index) const
{
    CoeffMap result;
    bool firstFactor = true;
    for (auto part = index.rbegin(); part != index.rend(); ++part)
      {
        const CoeffMap& factor =
            powerSumPartToCompleteMapViaLogarithmFormula(*part);
        if (firstFactor)
          {
            result = factor;
            firstFactor = false;
          }
        else
          result = multiplyCoeffMaps(result, factor);
      }
    return firstFactor ? oneCoeffMap() : result;
}

CoeffMap SymmetricEngineRing::powerSumIndexToElementaryMapViaLogarithmFormula(
    const Partition& index) const
{
    CoeffMap result;
    bool firstFactor = true;
    for (auto part = index.rbegin(); part != index.rend(); ++part)
      {
        const CoeffMap& factor =
            powerSumPartToElementaryMapViaLogarithmFormula(*part);
        if (firstFactor)
          {
            result = factor;
            firstFactor = false;
          }
        else
          result = multiplyCoeffMaps(result, factor);
      }
    return firstFactor ? oneCoeffMap() : result;
}

void SymmetricEngineRing::addPowerSumIndexToCompleteMapViaLogarithmFormula(
    const Partition& index,
    ring_elem coefficient,
    CoeffMap& result) const
{
    addScaledCoeffMap(
        result, coefficient,
        powerSumIndexToCompleteMapViaLogarithmFormula(index));
}

CoeffMap SymmetricEngineRing::powerSumsToCompleteMapViaLogarithmFormula(
    ring_elem f) const
{
    CoeffMap result;
    for (const auto& term : polyValue(f)->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during basis conversion");
            return CoeffMap{};
          }
        addPowerSumIndexToCompleteMapViaLogarithmFormula(
            index, term.coeff, result);
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumsToCompleteViaLogarithmFormula(
    ring_elem f, int completeId, int completeOrder) const
{
    CoeffMap result = powerSumsToCompleteMapViaLogarithmFormula(f);
    if (error()) return zero();
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

CoeffMap SymmetricEngineRing::multiplySchurExpansionViaRowPieri(
    const CoeffMap& source,
    int row) const
{
    CoeffMap result;
    if (row == 0) return source;
    if (row < 0)
      {
        ERROR("expected nonnegative h index during recursive h-to-Schur conversion");
        return result;
      }
    for (const auto& term : source)
      for (const auto& nu : schurTimesCompleteViaHorizontalPieri(term.first, row))
        addNormalizedCoeff(result, nu, term.second);
    return result;
  }

CoeffMap SymmetricEngineRing::completeToSchurCoefficientsViaRecursiveTransition(
    const CoeffMap& hCoeffs) const
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
          {
            CoeffMap summand =
                completeToSchurCoefficientsViaRecursiveTransition(found->second);
            for (const auto& item : summand)
              addNormalizedCoeff(result, item.first, item.second);
          }
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
    requireBasis(hBasisId);
    requireBasis(schurId);

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
            mpz_class chi = characterValueWithinLimits(lambda, mu);
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
            mpz_class chi = characterTableValue(table, row, col);
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
                             registeredPowerSumBasisId(),
                             displayForBasis(registeredPowerSumBasisId()),
                             basisOrderForId(registeredPowerSumBasisId()),
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
            mpz_class chi = characterValueWithinLimits(lambda, mu);
            if (chi == 0) continue;
            ring_elem term = basisElementFromIndex(schurId, lambda);
            mpz_class signedChi = sign * chi;
            result = add(
                result,
                scaled(
                    coefficientRing->from_int(signedChi.get_mpz_t()), term));
          }
      }
    else
      {
        size_t col = muCol->second;
        for (size_t row = 0; row < table.partitions.size(); ++row)
          {
            mpz_class chi = characterTableValue(table, row, col);
            if (chi == 0) continue;
            ring_elem term =
                basisElementFromIndex(schurId, table.partitions[row]);
            mpz_class signedChi = sign * chi;
            result = add(
                result,
                scaled(
                    coefficientRing->from_int(signedChi.get_mpz_t()), term));
          }
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumIndexToSchurViaCharacters(
    const Partition& mu,
    int schurId,
    int schurOrder) const
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
                const mpz_class& chi = contribution.second;
                ring_elem termCoeff =
                    coefficientRing->mult(
                        coefficientRing->from_int(chi.get_mpz_t()),
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
    ring_elem expansion = jacobiTrudi(lambda, Partition{}, generatorId);
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
            mpz_class chi = muCol == table.partitionRows.end()
                ? characterValueWithinLimits(table.partitions[row], mu)
                : characterTableValue(table, row, muCol->second);
            if (chi == 0) continue;
            if (omegaStyle && ((degree - partitionLength(mu)) % 2 != 0))
              chi = -chi;
            entry.contributions.push_back({input, chi});
          }
        if (!entry.contributions.empty()) recipe.push_back(std::move(entry));
      }

    requireCacheEntryCapacity("power-sum-to-Schur recipe cache");
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
              addCoeff(result, product.partition, coeff);
            }
        }
    return result;
  }

ring_elem SymmetricEngineRing::monomialToPowerSumsViaTransitionMatrix(
    const Partition& lambda,
    bool forgotten) const
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
        std::vector<Partition> parts =
            partitionsOfWithinLimits(d, "Hall-Littlewood conversion");
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
                mpz_class transition =
                    pToMonomialCoefficientWithinLimits(
                        parts[col], parts[row]);
                if (transition == 0 || coefficientRing->is_zero(coeffs[col])) continue;
                rhs = coefficientRing->subtract(
                    rhs,
                    coefficientRing->mult(
                        coefficientRing->from_int(transition.get_mpz_t()),
                        coeffs[col]));
              }
            mpz_class diagonal =
                pToMonomialCoefficientWithinLimits(
                    parts[row], parts[row]);
            coeffs[row] = coefficientQuotient(
                rhs, coefficientRing->from_int(diagonal.get_mpz_t()));
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
            registeredPowerSumBasisId(),
            displayForBasis(registeredPowerSumBasisId()),
            basisOrderForId(registeredPowerSumBasisId()),
            true);
        requireCacheEntryCapacity("monomial transition cache");
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
    for (const auto& mu : partitionsOfWithinLimits(d, "Hall-Littlewood conversion"))
      {
        mpz_class count =
            pToMonomialCoefficientWithinLimits(lambda, mu);
        if (count == 0) continue;
        if (forgotten && sign < 0) count = -count;
        ring_elem coeff = coefficientRing->from_int(count.get_mpz_t());
        ring_elem term = basisElementFromIndex(targetBasisId, mu);
        result = add(result, scaled(coeff, term));
      }
    return result;
  }

// ============================================================================
// Hall-Littlewood Generator Bases
// ============================================================================
// Classical and logarithmic conversions for the multiplicative q and b bases.

ring_elem
SymmetricEngineRing::hallLittlewoodGeneratorPartToPowerSumsViaClassicalFormula(
    int n,
    bool omega) const
{
    if (n == 0) return one();
    auto& cache = omega ? hallLittlewoodBGeneratorToPowerSumsCache
                        : hallLittlewoodQGeneratorToPowerSumsCache;
    auto cached = cache.find(n);
    if (cached != cache.end()) return copyPolyValue(polyValue(cached->second));
    CoeffMap coefficients;
    for (const auto& mu : partitionsOfWithinLimits(n, "Hall-Littlewood conversion"))
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
        registeredPowerSumBasisId(),
        displayForBasis(registeredPowerSumBasisId()),
        basisOrderForId(registeredPowerSumBasisId()),
        true);
    requireCacheEntryCapacity("Hall-Littlewood generator cache");
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
    for (const auto& mu : partitionsOfWithinLimits(n, "Hall-Littlewood conversion"))
      {
        long sign = 1;
        if (omega && ((n - static_cast<int>(mu.size())) % 2 == 1)) sign = -1;
        ring_elem rational = rationalCoefficient(sign, zValue(mu));
        if (error()) return CoeffMap{};
        addCoeff(result, mu, rational);
      }
    requireCacheEntryCapacity("Hall-Littlewood generator quotient cache");
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

    requireCacheEntryCapacity("Hall-Littlewood logarithm cache");
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
    requireCacheEntryCapacity("Hall-Littlewood C-factor cache");
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
    for (const Partition& lambda : partitionsOfWithinLimits(n, "Hall-Littlewood conversion"))
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
    requireCacheEntryCapacity("Hall-Littlewood Green cache");
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
        if (std::any_of(item.first.begin(), item.first.end(),
                        [](int part) { return part < 0; }))
          continue;
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
    requireCacheEntryCapacity("Hall-Littlewood raising cache");
    hallLittlewoodRaisingGeneratorMapCache[lambda] = result;
    return result;
  }

ring_elem
SymmetricEngineRing::hallLittlewoodCapitalToPowerSumsViaRaisingOperators(
    const Partition& lambda,
    bool omega) const
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
                             registeredPowerSumBasisId(),
                             displayForBasis(registeredPowerSumBasisId()),
                             basisOrderForId(registeredPowerSumBasisId()),
                             true);
  }

ring_elem
SymmetricEngineRing::hallLittlewoodNormalizedToPowerSumsViaCapitalNormalization(
    const Partition& lambda,
    bool omega) const
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
    for (const auto& nu : partitionsOfWithinLimits(d, "Hall-Littlewood conversion"))
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
        ring_elem term = basisElementFromIndex(targetId, nu);
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
    const BasisKind coordinateKind =
        omega ? BasisKind::Forgotten : BasisKind::Monomial;
    ring_elem inPowerSums =
        canonicalExpressionToPowerSumsViaBasisFormulas(
            skewCapital, coordinateKind);
    if (error()) return zero();
    const BasisKind capitalKind = omega ? BasisKind::HallLittlewoodB
                                        : BasisKind::HallLittlewoodQ;
    int capitalId = requiredBasisIdForKind(capitalKind);
    if (error()) return zero();
    ring_elem inCapital =
        powerSumsToHallLittlewoodViaTriangularReduction(
            inPowerSums,
            capitalId,
            displayForBasis(capitalId),
            basisOrderForId(capitalId));
    if (error()) return zero();
    const BasisKind normalizedKind = omega ? BasisKind::HallLittlewoodPOmega
                                           : BasisKind::HallLittlewoodP;
    const int normalizedId = requiredBasisIdForKind(normalizedKind);
    if (error()) return zero();
    ring_elem normalized = replaceSingleBasis(inCapital, capitalId, normalizedId);
    if (error()) return zero();
    return canonicalExpressionToPowerSumsViaBasisFormulas(
        normalized, normalizedKind);
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
        {
          bool omega = basisKind == BasisKind::HallLittlewoodB;
          ring_elem skewCapital =
              skewQOrBFunction(lambda, mu, omega);
          if (error()) return zero();
          return canonicalExpressionToPowerSumsViaBasisFormulas(
              skewCapital,
              omega ? BasisKind::Forgotten : BasisKind::Monomial);
        }
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
            requireCacheEntryCapacity("Hall-Littlewood column cache");
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
        bool recursiveLimitExceeded = false;
        mpz_class multiplicity = labeledCycleAssignmentsToCompleteFactors(
            mu,
            raisingTerm.first,
            computationLimits.maxRecursiveStates,
            recursiveLimitExceeded);
        if (recursiveLimitExceeded)
          throw exc::engine_error(
              "Hall-Littlewood assignment recursion exceeds MaxRecursiveStates "
              "in the symmetricRing ComputationLimits option");
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
    requireCacheEntryCapacity("Hall-Littlewood coefficient cache");
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
        for (const Partition& lambda :
             partitionsOfWithinLimits(
                 partitionWeight(cycleType), "power-sum to Schur conversion"))
          {
            ring_elem coefficient =
                powerSumIndexToHallLittlewoodCapitalCoefficientViaGreenPolynomialDuality(
                    cycleType, lambda);
            if (error()) return zero();
            addCoeff(capitalColumn, lambda, coefficient);
          }
        requireCacheEntryCapacity("Hall-Littlewood column cache");
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
// Termwise Conversion Kernels
// ============================================================================
// Per-index formulas used by the general conversion fallback.

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
          return basisElementFromIndex(targetBasisId, index);
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
    ring_elem term = basisElementFromIndex(basisId, straightened.second);
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
      {
        if (trimmed.back() < 0) return zero();
        return basisElementFromIndex(basisId, trimmed);
      }

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

ring_elem SymmetricEngineRing::straightenBasisElement(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    const int basisId = atomBasisIdAt(monomial, pos);
    const BasisKind basisKind = basisKindForId(basisId);
    std::string display = displayForBasis(basisId);
    if (atomIsSkewAt(monomial, pos))
      {
        auto *poly = new SymmetricRingPoly;
        poly->terms.push_back(
            {coefficientRing->one(),
             monomialFromKey(atomBlockAt(monomial, pos))});
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
          return basisElementFromIndex(basisId, index);
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
        ring_elem capitalPowerSums =
            canonicalExpressionToPowerSumsViaBasisFormulas(normalizedCapital,
                                                           capitalKind);
        if (error()) return zero();
        return powerSumsToHallLittlewoodCapitalViaTriangularReduction(
            capitalPowerSums,
            basisId,
            display,
            basisOrderForId(basisId));
        }
      default:
        break;
      }
    return basisElementFromIndex(basisId, index);
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
