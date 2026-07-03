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

ring_elem SymmetricEngineRing::omegaPowerSums(ring_elem f) const
{
    const auto *poly = polyValue(f);
    VECTOR(SymmetricTerm) terms;
    terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during omega");
            return zero();
          }
        long sign = ((partitionWeight(index) - partitionLength(index)) % 2 == 0) ? 1 : -1;
        terms.push_back({coefficientRing->mult(coefficientRing->from_long(sign),
                                               term.coeff),
                         term.monomial});
      }
    return fromTermVector(terms, true);
  }

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

ring_elem SymmetricEngineRing::straightenSchurAtom(const Partition& alpha,
                                const std::string& display) const
{
    auto straightened = straightenSchurIndex(alpha);
    if (straightened.first == 0) return zero();
    int id = requiredBasisIdForDisplay(display);
    if (error()) return zero();
    ring_elem term = basisElementFromIndex(id,
                                           display,
                                           basisOrderForId(id),
                                           isMultiplicativeBasis(id),
                                           straightened.second);
    if (straightened.first < 0) term = negate(term);
    return term;
  }

ring_elem SymmetricEngineRing::omegaSchurAtomAsSchur(const Partition& alpha) const
{
    auto straightened = straightenSchurIndex(alpha);
    if (straightened.first == 0) return zero();
    int schurId = requiredBasisIdForDisplay("S");
    if (error()) return zero();
    ring_elem term = basisElementFromIndex(schurId,
                                           "S",
                                           basisOrderForId(schurId),
                                           isMultiplicativeBasis(schurId),
                                           conjugatePartition(straightened.second));
    if (straightened.first < 0) term = negate(term);
    return term;
  }

ring_elem SymmetricEngineRing::straightenHallCapitalAtom(const Partition& alpha,
                                      const std::string& display) const
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
      return basisElementForDisplay(display, trimmed);

    int s = trimmed[bad];
    int r = trimmed[bad + 1];
    int diff = r - s;
    int top = diff / 2;
    ring_elem result =
        scaled(hallLittlewoodParameter,
               straightenHallCapitalAtom(replaceAdjacentPair(trimmed, bad, r, s),
                                          display));
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
                            straightenHallCapitalAtom(
                                replaceAdjacentPair(trimmed, bad, r - i, s + i),
                                display)));
      }
    return result;
  }

ring_elem SymmetricEngineRing::straightenAtom(const SymmetricMonomial& monomial, size_t pos) const
{
    std::string display = displayForBasis(atomBasisIdAt(monomial, pos));
    if (atomIsSkewAt(monomial, pos))
      {
        auto *poly = new SymmetricRingPoly;
        poly->terms.push_back({coefficientRing->one(), monomialFromKey(atomBlockAt(monomial, pos))});
        return makePolyValue(poly);
      }
    Partition index = atomIndex(monomial, pos);
    if (display == "S" || display == "Somega")
      return straightenSchurAtom(index, display);
    if (display == "Q" || display == "B")
      return straightenHallCapitalAtom(index, display);
    if (display == "P" || display == "R")
      {
        std::string capitalDisplay = display == "P" ? "Q" : "B";
        ring_elem straightCapital = straightenHallCapitalAtom(index, capitalDisplay);
        if (error()) return zero();
        return powerSumsToHallCapitalTarget(elementToPowerSums(straightCapital),
                                            requiredBasisIdForDisplay(display),
                                            display,
                                            basisOrderForId(requiredBasisIdForDisplay(display)));
      }
    return basisElementFromIndex(atomBasisIdAt(monomial, pos),
                                 display,
                                 atomOrderAt(monomial, pos),
                                 isMultiplicativeBasis(atomBasisIdAt(monomial, pos)),
                                 index);
  }

ring_elem SymmetricEngineRing::straightenMonomial(const SymmetricMonomial& monomial) const
{
    ring_elem result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        ring_elem factor = straightenAtom(monomial, pos);
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

ring_elem SymmetricEngineRing::omegaDirectAtom(const SymmetricMonomial& monomial,
                            size_t pos,
                            const std::map<int, OmegaTarget>& omegaTargets,
                            bool useSomega) const
{
    int basisId = atomBasisIdAt(monomial, pos);
    std::string display = displayForBasis(basisId);
    if (basisId == powerSumBasisId)
      {
        Partition index = atomIndex(monomial, pos);
        long sign = ((partitionWeight(index) - partitionLength(index)) % 2 == 0) ? 1 : -1;
        ring_elem term = basisElementFromIndex(powerSumBasisId,
                                               "p",
                                               basisOrderForId(powerSumBasisId),
                                               true,
                                               index);
        return scaled(coefficientRing->from_long(sign), term);
      }

    if (!useSomega && display == "S")
      {
        int schurId = requiredBasisIdForDisplay("S");
        if (error()) return zero();
        if (atomIsSkewAt(monomial, pos))
          return basisElementFromSkewIndex(schurId,
                                           "S",
                                           basisOrderForId(schurId),
                                           false,
                                           conjugatePartition(atomOuterIndex(monomial, pos)),
                                           conjugatePartition(atomInnerIndex(monomial, pos)));
        return omegaSchurAtomAsSchur(atomIndex(monomial, pos));
      }

    if (!useSomega && display == "Somega" && !atomIsSkewAt(monomial, pos))
      return straightenSchurAtom(atomIndex(monomial, pos), "S");

    auto target = omegaTargets.find(basisId);
    if (target != omegaTargets.end())
      {
        Partition payload = atomIndex(monomial, pos);
        std::string targetDisplay = displayForBasis(target->second.basisId);
        rememberBasis(target->second.basisId,
                      targetDisplay,
                      target->second.order,
                      target->second.isMultiplicative);
        auto *poly = new SymmetricRingPoly;
        SymmetricMonomial targetMonomial;
        appendAtomBlock(targetMonomial,
                        makeAtomBlock(target->second.order,
                                      target->second.basisId,
                                      atomInnerLengthAt(monomial, pos),
                                      payload));
        poly->terms.push_back({coefficientRing->one(),
                               canonicalMonomial(targetMonomial)});
        return makePolyValue(poly);
      }

    ring_elem inPowerSums = atomToPowerSums(monomial, pos);
    if (error()) return zero();
    return omegaPowerSums(inPowerSums);
  }

ring_elem SymmetricEngineRing::omegaMonomial(const SymmetricMonomial& monomial,
                          const std::map<int, OmegaTarget>& omegaTargets,
                          bool useSomega) const
{
    ring_elem result = one();
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        ring_elem factor = omegaDirectAtom(monomial, pos, omegaTargets, useSomega);
        if (error()) return zero();
        result = mult(result, factor);
        pos += atomLengthAt(monomial, pos);
      }
    return result;
  }

std::map<int, OmegaTarget> SymmetricEngineRing::omegaTargetMap(M2_arrayint omegaMap) const
{
    std::map<int, OmegaTarget> result;
    if (omegaMap == nullptr) return result;
    if (omegaMap->len % 4 != 0)
      {
        ERROR("invalid omega metadata map");
        return result;
      }
    for (int i = 0; i < omegaMap->len; i += 4)
      result[omegaMap->array[i]] =
          OmegaTarget{omegaMap->array[i + 1],
                      omegaMap->array[i + 2],
                      omegaMap->array[i + 3] != 0};
    return result;
  }

ring_elem SymmetricEngineRing::omegaInvolution(ring_elem f, M2_arrayint omegaMap, bool useSomega) const
{
    std::map<int, OmegaTarget> targets = omegaTargetMap(omegaMap);
    if (error()) return zero();
    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        ring_elem converted = omegaMonomial(term.monomial, targets, useSomega);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, converted));
      }
    return result;
  }

ring_elem SymmetricEngineRing::straighten(ring_elem f) const
{
    return straightenElement(f);
  }

} // namespace symmetric_rings
