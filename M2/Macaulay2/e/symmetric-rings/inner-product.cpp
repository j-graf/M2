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

bool SymmetricEngineRing::singleBasisIndexFromMonomial(const SymmetricMonomial& monomial,
                                    int basisId,
                                    Partition& index) const
{
    index.clear();
    if (monomial.data.empty()) return false;
    if (atomIsSkewAt(monomial, 0)) return false;
    if (atomBasisIdAt(monomial, 0) != basisId) return false;
    if (atomLengthAt(monomial, 0) != monomial.data.size()) return false;
    index = atomIndex(monomial, 0);
    return true;
  }

CoeffMap SymmetricEngineRing::coefficientsInBasis(ring_elem f, int basisId) const
{
    CoeffMap result;
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (term.monomial.data.empty())
          {
            index = Partition{};
          }
        else if (!singleBasisIndexFromMonomial(term.monomial, basisId, index))
          {
            ERROR("expected an expression in a single symmetric-function basis");
            return CoeffMap{};
          }
        addCoeff(result, index, term.coeff);
      }
    return result;
  }

bool SymmetricEngineRing::coefficientsInBasisIfPossible(ring_elem f,
                                      int basisId,
                                      CoeffMap& result) const
{
    result.clear();
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (term.monomial.data.empty())
          {
            index = Partition{};
          }
        else if (!singleBasisIndexFromMonomial(term.monomial, basisId, index))
          {
            result.clear();
            return false;
          }
        addCoeff(result, index, term.coeff);
      }
    return true;
  }

ring_elem SymmetricEngineRing::coefficientPairing(const CoeffMap& fCoeffs,
                               const CoeffMap& gCoeffs) const
{
    ring_elem result = coefficientRing->zero();
    for (const auto& item : fCoeffs)
      {
        auto it = gCoeffs.find(item.first);
        if (it != gCoeffs.end())
          result = coefficientRing->add(
              result,
              coefficientRing->mult(item.second, it->second));
      }
    return result;
  }

ring_elem SymmetricEngineRing::powerSumInnerProductFactor(const Partition& lambda) const
{
    ring_elem numerator = rationalCoefficient(zValue(lambda), 1);
    if (error()) return coefficientRing->zero();
    return coefficientQuotient(numerator, hallLittlewoodFactor(lambda));
  }

ring_elem SymmetricEngineRing::powerSumPairing(const CoeffMap& fCoeffs,
                            const CoeffMap& gCoeffs) const
{
    ring_elem result = coefficientRing->zero();
    for (const auto& item : fCoeffs)
      {
        auto it = gCoeffs.find(item.first);
        if (it == gCoeffs.end()) continue;
        ring_elem factor = powerSumInnerProductFactor(item.first);
        if (error()) return coefficientRing->zero();
        result = coefficientRing->add(
            result,
            coefficientRing->mult(coefficientRing->mult(item.second, it->second),
                                  factor));
      }
    return result;
  }

std::map<int, InnerProductTarget> SymmetricEngineRing::innerProductTargetMap(M2_arrayint innerProductMap) const
{
    std::map<int, InnerProductTarget> result;
    if (innerProductMap == nullptr) return result;
    if (innerProductMap->len % 3 != 0)
      {
        ERROR("invalid inner product metadata map");
        return result;
      }
    for (int i = 0; i < innerProductMap->len; i += 3)
      result[innerProductMap->array[i]] =
          InnerProductTarget{innerProductMap->array[i + 1],
                             innerProductMap->array[i + 2]};
    return result;
  }

bool SymmetricEngineRing::directHallInnerProductFromMetadata(
      ring_elem f,
      ring_elem g,
      const std::map<int, InnerProductTarget>& metadata,
      ring_elem& result) const
{
    for (const auto& item : metadata)
      {
        CoeffMap fCoeffs;
        CoeffMap gCoeffs;
        if (!coefficientsInBasisIfPossible(f, item.first, fCoeffs)) continue;
        if (!coefficientsInBasisIfPossible(g, item.second.dualBasisId, gCoeffs)) continue;
        if (item.second.kind == 1)
          result = coefficientPairing(fCoeffs, gCoeffs);
        else if (item.second.kind == 2)
          result = powerSumPairing(fCoeffs, gCoeffs);
        else
          continue;
        return !error();
      }
    return false;
  }

bool SymmetricEngineRing::directPowerSumSchurInnerProduct(ring_elem f,
                                       ring_elem g,
                                       ring_elem& result) const
{
    int pId = requiredBasisIdForDisplay("p");
    int schurId = requiredBasisIdForDisplay("S");
    if (error()) return false;

    CoeffMap pCoeffs;
    CoeffMap schurCoeffs;
    if (!coefficientsInBasisIfPossible(f, pId, pCoeffs) ||
        !coefficientsInBasisIfPossible(g, schurId, schurCoeffs))
      {
        pCoeffs.clear();
        schurCoeffs.clear();
        if (!coefficientsInBasisIfPossible(g, pId, pCoeffs) ||
            !coefficientsInBasisIfPossible(f, schurId, schurCoeffs))
          return false;
      }

    result = coefficientRing->zero();
    for (const auto& pTerm : pCoeffs)
      {
        int degree = partitionWeight(pTerm.first);
        Partition mu = normalizePartition(pTerm.first);
        for (const auto& schurTerm : schurCoeffs)
          {
            if (partitionWeight(schurTerm.first) != degree) continue;
            if (!isPartitionIndex(schurTerm.first)) return false;
            int chi = characterValue(schurTerm.first, mu);
            if (chi == 0) continue;
            ring_elem contribution = coefficientRing->mult(pTerm.second,
                                                           schurTerm.second);
            if (chi != 1)
              contribution =
                  coefficientRing->mult(coefficientRing->from_long(chi),
                                        contribution);
            result = coefficientRing->add(result, contribution);
          }
      }
    return true;
  }

ring_elem SymmetricEngineRing::hallInnerProductElements(ring_elem f,
                                     ring_elem g,
                                     const std::map<int, InnerProductTarget>& metadata) const
{
    ring_elem direct;
    if (directHallInnerProductFromMetadata(f, g, metadata, direct))
      return direct;
    if (error()) return coefficientRing->zero();
    if (directPowerSumSchurInnerProduct(f, g, direct))
      return direct;
    if (error()) return coefficientRing->zero();

    int pId = requiredBasisIdForDisplay("p");
    if (error()) return coefficientRing->zero();

    ring_elem fP = toBasis(f,
                           pId,
                           "p",
                           basisOrderForId(pId),
                           isMultiplicativeBasis(pId),
                           pId,
                           "p",
                           basisOrderForId(pId),
                           isMultiplicativeBasis(pId));
    if (error()) return coefficientRing->zero();
    ring_elem gP = toBasis(g,
                           pId,
                           "p",
                           basisOrderForId(pId),
                           isMultiplicativeBasis(pId),
                           pId,
                           "p",
                           basisOrderForId(pId),
                           isMultiplicativeBasis(pId));
    if (error()) return coefficientRing->zero();

    CoeffMap fCoeffs = coefficientsInBasis(fP, pId);
    if (error()) return coefficientRing->zero();
    CoeffMap gCoeffs = coefficientsInBasis(gP, pId);
    if (error()) return coefficientRing->zero();

    return powerSumPairing(fCoeffs, gCoeffs);
  }

ring_elem SymmetricEngineRing::skewQOrBFunction(const Partition& lambda,
                             const Partition& mu,
                             bool omega) const
{
    int d = partitionWeight(lambda) - partitionWeight(mu);
    if (d < 0) return zero();

    ring_elem lambdaTerm = basisElementForDisplay(omega ? "B" : "Q", lambda);
    ring_elem muTerm = basisElementForDisplay(omega ? "R" : "P", mu);
    if (error()) return zero();

    int targetId = requiredBasisIdForDisplay(omega ? "ff" : "m");
    std::string targetDisplay = omega ? "ff" : "m";
    if (error()) return zero();

    ring_elem result = zero();
    for (const auto& nu : partitionsOf(d))
      {
        ring_elem generator = basisElementForDisplay(omega ? "b" : "q", nu);
        ring_elem test = mult(muTerm, generator);
        ring_elem coeff = hallInnerProductElements(lambdaTerm, test);
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

ring_elem SymmetricEngineRing::skewPOrRToPowerSums(const Partition& lambda,
                                const Partition& mu,
                                bool omega) const
{
    ring_elem skewCapital = skewQOrBFunction(lambda, mu, omega);
    if (error()) return zero();
    int pId = requiredBasisIdForDisplay("p");
    int capitalId = requiredBasisIdForDisplay(omega ? "B" : "Q");
    if (error()) return zero();

    ring_elem inCapital = toBasis(skewCapital,
                                  pId,
                                  "p",
                                  basisOrderForId(pId),
                                  isMultiplicativeBasis(pId),
                                  capitalId,
                                  omega ? "B" : "Q",
                                  basisOrderForId(capitalId),
                                  isMultiplicativeBasis(capitalId));
    if (error()) return zero();
    ring_elem normalized = replaceSingleBasis(inCapital,
                                              capitalId,
                                              omega ? "R" : "P");
    if (error()) return zero();
    return elementToPowerSums(normalized);
  }

ring_elem SymmetricEngineRing::skewHallLittlewoodToPowerSums(const Partition& lambda,
                                          const Partition& mu,
                                          const std::string& display) const
{
    if (display == "Q" || display == "B")
      return elementToPowerSums(skewQOrBFunction(lambda, mu, display == "B"));
    if (display == "P" || display == "R")
      return skewPOrRToPowerSums(lambda, mu, display == "R");
    ERROR("expected a skew Hall-Littlewood basis atom");
    return zero();
  }

ring_elem SymmetricEngineRing::hallInnerProduct(ring_elem f,
                             ring_elem g,
                             M2_arrayint innerProductMap) const
{
    std::map<int, InnerProductTarget> metadata =
        innerProductTargetMap(innerProductMap);
    if (error()) return coefficientRing->zero();
    return hallInnerProductElements(f, g, metadata);
  }

} // namespace symmetric_rings
