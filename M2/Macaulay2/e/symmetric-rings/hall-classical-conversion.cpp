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
} // namespace symmetric_rings
