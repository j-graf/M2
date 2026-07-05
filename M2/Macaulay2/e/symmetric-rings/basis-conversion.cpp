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

bool SymmetricEngineRing::atomToSchurFactors(
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

bool SymmetricEngineRing::schurProductMonomialToSchur(
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
        if (!atomToSchurFactors(monomial, pos, targetBasisId, factors))
          return false;

        for (const auto& index : factors)
          {
            CoeffMap next;
            for (const auto& term : current)
              for (const auto& product : lrProduct(term.first, index))
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

bool SymmetricEngineRing::fastSchurProductMonomialToSchur(
    const SymmetricMonomial& monomial,
    int targetBasisId,
    const std::string& targetDisplay,
    int targetDisplayOrder,
    ring_elem& result) const
{
    if (targetDisplay != "S") return false;
    std::vector<SchurCompatibleFactor> factors;
    if (!schurCompatibleFactorsFromMonomial(monomial, targetBasisId, factors))
      return false;
    return multiplySchurCompatibleFactorsToSchur(std::move(factors),
                                                targetBasisId,
                                                targetDisplay,
                                                targetDisplayOrder,
                                                result);
  }

bool SymmetricEngineRing::schurCompatibleFactorsFromMonomial(
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
            for (const auto& item : skewSchurExpansion(atomOuterIndex(monomial, pos),
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

bool SymmetricEngineRing::multiplySchurCompatibleFactorsToSchur(
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
        const auto& product = directLRProduct(factors[0].index, factors[1].index);
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
        CoeffMap next;
        for (const auto& term : current)
          {
            if (factor.kind == SchurCompatibleFactor::Horizontal)
              {
                int row = factor.index.empty() ? 0 : factor.index.front();
                for (const auto& nu : horizontalStripProducts(term.first, row))
                  addCoeff(next, nu, term.second);
              }
            else if (factor.kind == SchurCompatibleFactor::Vertical)
              {
                int col = factor.index.empty() ? 0 : factor.index.front();
                for (const auto& nu : verticalStripProducts(term.first, col))
                  addCoeff(next, nu, term.second);
              }
            else if (factor.kind == SchurCompatibleFactor::PowerSum)
              {
                int part = factor.index.empty() ? 0 : factor.index.front();
                for (const auto& product : powerSumSchurProduct(term.first, part))
                  {
                    ring_elem coeff = product.coefficient == 1
                        ? term.second
                        : coefficientRing->mult(cachedInteger(product.coefficient),
                                                term.second);
                    addCoeff(next, product.nu, coeff);
                  }
              }
            else if (factor.kind == SchurCompatibleFactor::SchurExpansion)
              {
                for (const auto& expansionTerm : factor.expansion)
                  {
                    ring_elem expansionCoeff = coefficientRing->mult(term.second,
                                                                     expansionTerm.second);
                    if (coefficientRing->is_zero(expansionCoeff)) continue;
                    for (const auto& product : directLRProduct(term.first,
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
                for (const auto& product : directLRProduct(term.first, factor.index))
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

bool SymmetricEngineRing::fastElementToSchur(ring_elem f,
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
        if (!fastSchurProductMonomialToSchur(term.monomial,
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

bool SymmetricEngineRing::fastProductToSchur(ring_elem f,
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
          if (!fastSchurProductMonomialToSchur(productMonomial,
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

bool SymmetricEngineRing::monomialLikeAtomToCoeffMap(
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

bool SymmetricEngineRing::monomialLikeMonomialToTarget(
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
        if (!monomialLikeAtomToCoeffMap(monomial, pos, targetDisplay, factor))
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

bool SymmetricEngineRing::fastProductToMonomialLike(
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
          if (!monomialLikeMonomialToTarget(productMonomial,
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

bool SymmetricEngineRing::fastProductToTarget(ring_elem f,
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
      return fastProductToSchur(f,
                                g,
                                targetBasisId,
                                targetDisplay,
                                targetDisplayOrder,
                                result);

    if (targetDisplay == "m" || isForgottenDisplay(targetDisplay))
      return fastProductToMonomialLike(f,
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
        if (elementToDirectTarget(g,
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
        if (elementToDirectTarget(f,
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

bool SymmetricEngineRing::termToDirectTarget(const SymmetricTerm& term,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          bool targetIsMultiplicative,
                          ring_elem& result) const
{
    if (schurProductMonomialToSchur(term.monomial,
                                    targetBasisId,
                                    targetDisplay,
                                    targetDisplayOrder,
                                    result))
      {
        result = scaled(term.coeff, result);
        return true;
      }

    if (targetDisplay != "h" && targetDisplay != "e" && !targetIsMultiplicative)
      return false;

    ring_elem converted;
    if (!monomialToDirectTarget(term.monomial,
                                targetBasisId,
                                targetDisplay,
                                targetDisplayOrder,
                                targetIsMultiplicative,
                                converted))
      return false;
    result = scaled(term.coeff, converted);
    return true;
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

    if (targetDisplay == "S")
      {
        int hId = requiredBasisIdForDisplay("h");
        if (error()) return zero();
        int hOrder = basisOrderForId(hId);
        bool hIsMultiplicative = isMultiplicativeBasis(hId);

        if (singleBasisId(f) == hId)
          return hToSchurViaRecTrans(f,
                                     hId,
                                     "h",
                                     hOrder,
                                     hIsMultiplicative,
                                     targetBasisId,
                                     targetDisplay,
                                     targetOrder);

        const auto *poly = polyValue(f);
        ring_elem result = zero();
        VECTOR(SymmetricTerm) remainderTerms;
        remainderTerms.reserve(poly->terms.size());
        bool convertedAny = false;

        for (const auto& term : poly->terms)
          {
            ring_elem converted;
            if (termToDirectTarget(term,
                                   targetBasisId,
                                   targetDisplay,
                                   targetOrder,
                                   targetIsMultiplicative,
                                   converted))
              {
                result = add(result, converted);
                convertedAny = true;
              }
            else
              {
                if (error()) return zero();
                remainderTerms.push_back(term);
              }
          }

        if (remainderTerms.empty()) return result;
        if (convertedAny)
          {
            ring_elem remainder = fromTermVector(remainderTerms, false);
            ring_elem convertedRemainder = toBasis(remainder,
                                                   pBasisId,
                                                   pDisplay,
                                                   pOrder,
                                                   pIsMultiplicative,
                                                   targetBasisId,
                                                   targetDisplay,
                                                   targetOrder,
                                                   targetIsMultiplicative);
            if (error()) return zero();
            return add(result, convertedRemainder);
          }

        ring_elem inComplete = toBasis(f,
                                       pBasisId,
                                       pDisplay,
                                       pOrder,
                                       pIsMultiplicative,
                                       hId,
                                       "h",
                                       hOrder,
                                       hIsMultiplicative);
        if (error()) return zero();
        return hToSchurViaRecTrans(inComplete,
                                   hId,
                                   "h",
                                   hOrder,
                                   hIsMultiplicative,
                                   targetBasisId,
                                   targetDisplay,
                                   targetOrder);
      }

    ring_elem direct;
    if (elementToDirectTarget(f,
                              targetBasisId,
                              targetDisplay,
                              targetOrder,
                              targetIsMultiplicative,
                              direct))
      return direct;
    if (error()) return zero();

    const auto *poly = polyValue(f);
    ring_elem result = zero();
    VECTOR(SymmetricTerm) remainderTerms;
    remainderTerms.reserve(poly->terms.size());
    bool convertedAny = false;

    for (const auto& term : poly->terms)
      {
        ring_elem converted;
        if (termToDirectTarget(term,
                               targetBasisId,
                               targetDisplay,
                               targetOrder,
                               targetIsMultiplicative,
                               converted))
          {
            result = add(result, converted);
            convertedAny = true;
          }
        else
          {
            if (error()) return zero();
            remainderTerms.push_back(term);
          }
      }

    if (remainderTerms.empty()) return result;
    if (!convertedAny)
      {
        ring_elem inPowerSums = elementToPowerSums(f);
        if (error()) return zero();
        return powerSumsToTarget(inPowerSums,
                                 targetBasisId,
                                 targetDisplay,
                                 targetOrder,
                                 targetIsMultiplicative);
      }

    ring_elem remainder = fromTermVector(remainderTerms, false);
    ring_elem convertedRemainder = toBasis(remainder,
                                           pBasisId,
                                           pDisplay,
                                           pOrder,
                                           pIsMultiplicative,
                                           targetBasisId,
                                           targetDisplay,
                                           targetOrder,
                                           targetIsMultiplicative);
    if (error()) return zero();
    return add(result, convertedRemainder);
  }

ring_elem SymmetricEngineRing::toSchurFast(ring_elem f,
                    int pBasisId,
                    const std::string& pDisplay,
                    int pOrder,
                    bool pIsMultiplicative,
                    int schurBasisId,
                    const std::string& schurDisplay,
                    int schurOrder) const
{
    rememberBasis(pBasisId, pDisplay, pOrder, pIsMultiplicative);
    rememberBasis(schurBasisId, schurDisplay, schurOrder, false);

    ring_elem direct;
    if (fastElementToSchur(f, schurBasisId, schurDisplay, schurOrder, direct))
      return direct;
    if (error()) return zero();

    return toBasis(f,
                   pBasisId,
                   pDisplay,
                   pOrder,
                   pIsMultiplicative,
                   schurBasisId,
                   schurDisplay,
                   schurOrder,
                   false);
  }

ring_elem SymmetricEngineRing::multiplyToSchurFast(ring_elem f,
                    ring_elem g,
                    int pBasisId,
                    const std::string& pDisplay,
                    int pOrder,
                    bool pIsMultiplicative,
                    int schurBasisId,
                    const std::string& schurDisplay,
                    int schurOrder) const
{
    return multiplyToBasisFast(f,
                               g,
                               pBasisId,
                               pDisplay,
                               pOrder,
                               pIsMultiplicative,
                               schurBasisId,
                               schurDisplay,
                               schurOrder,
                               false);
  }

ring_elem SymmetricEngineRing::multiplyToBasisFast(ring_elem f,
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

    ring_elem direct;
    if (fastProductToTarget(f,
                            g,
                            targetBasisId,
                            targetDisplay,
                            targetOrder,
                            targetIsMultiplicative,
                            direct))
      return direct;
    if (error()) return zero();

    ring_elem product = mult(f, g);
    if (error()) return zero();
    return toBasis(product,
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
