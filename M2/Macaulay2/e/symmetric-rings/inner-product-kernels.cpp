// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "rings/ZZ.hpp"

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Diagonal Pairing Kernels
// ============================================================================
// Coefficient intersections implement registered dual and power-sum diagonal pairings.

ring_elem SymmetricEngineRing::coefficientPairing(
    const CoeffMap& fCoeffs,
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

ring_elem SymmetricEngineRing::powerSumDiagonalFactor(
    const Partition& lambda,
    const InnerProductContext& context) const
{
    ring_elem numerator = rationalCoefficient(zValue(lambda), 1);
    if (error()) return coefficientRing->zero();
    switch (context.powerSumPairing)
      {
        case PowerSumPairingKind::OrdinaryHall:
          return numerator;
        case PowerSumPairingKind::HallLittlewood:
          return coefficientQuotient(numerator, hallLittlewoodFactor(lambda));
        case PowerSumPairingKind::SchurQ:
          ERROR("Schur-Q power-sum pairing is not implemented");
          return coefficientRing->zero();
        case PowerSumPairingKind::Macdonald:
          ERROR("Macdonald power-sum pairing is not implemented");
          return coefficientRing->zero();
      }
    ERROR("unknown power-sum pairing kind");
    return coefficientRing->zero();
  }

ring_elem SymmetricEngineRing::powerSumPairing(
    const CoeffMap& fCoeffs,
    const CoeffMap& gCoeffs,
    const InnerProductContext& context) const
{
    ring_elem result = coefficientRing->zero();
    for (const auto& item : fCoeffs)
      {
        auto it = gCoeffs.find(item.first);
        if (it == gCoeffs.end()) continue;
        ring_elem factor = powerSumDiagonalFactor(item.first, context);
        if (error()) return coefficientRing->zero();
        result = coefficientRing->add(
            result,
            coefficientRing->mult(coefficientRing->mult(item.second, it->second),
                                  factor));
      }
    return result;
  }


// ============================================================================
// Kostka-Number Pairings
// ============================================================================
// Ordinary Schur and Schur Omega pairings use Kostka and conjugate Kostka numbers.

mpz_class SymmetricEngineRing::kostkaNumberViaSemistandardTableaux(
    const Partition& shape,
    const Partition& content) const
{
    Partition lambda = normalizePartition(shape);
    Partition mu = normalizePartition(content);
    requirePartitionWithinWeightLimit(lambda, "Kostka shape");
    requirePartitionWithinWeightLimit(mu, "Kostka content");
    if (partitionWeight(lambda) != partitionWeight(mu)) return 0;
    int lambdaSum = 0;
    int muSum = 0;
    size_t dominanceLength = std::max(lambda.size(), mu.size());
    for (size_t i = 0; i < dominanceLength; ++i)
      {
        lambdaSum += partitionPart(lambda, i);
        muSum += partitionPart(mu, i);
        if (lambdaSum < muSum) return 0;
      }

    auto cacheKey = std::make_pair(lambda, mu);
    auto cached = kostkaNumberCache.find(cacheKey);
    if (cached != kostkaNumberCache.end()) return cached->second;

    const size_t cellCount = static_cast<size_t>(partitionWeight(lambda));
    if (!estimatedMemoryWithinLimit(
            cellCount,
            sizeof(std::pair<size_t, int>) + sizeof(int),
            "Kostka tableau"))
      return 0;
    std::vector<std::pair<size_t, int>> cells;
    cells.reserve(cellCount);
    for (size_t row = 0; row < lambda.size(); ++row)
      for (int col = 0; col < lambda[row]; ++col)
        cells.push_back({row, col});
    std::vector<std::vector<int>> tableau(lambda.size());
    for (size_t row = 0; row < lambda.size(); ++row)
      tableau[row].assign(lambda[row], 0);
    std::vector<int> remaining(mu.begin(), mu.end());
    mpz_class result = 0;
    size_t states = 0;
    std::function<void(size_t)> fill = [&](size_t cell) {
      if (!consumeRecursiveState(states, "Kostka tableau enumeration"))
        return;
      if (cell == cells.size())
        {
          ++result;
          return;
        }
      size_t row = cells[cell].first;
      int col = cells[cell].second;
      int left = col > 0 ? tableau[row][col - 1] : 0;
      int above = row > 0 && col < static_cast<int>(tableau[row - 1].size())
          ? tableau[row - 1][col] : 0;
      for (size_t value = 1; value <= remaining.size(); ++value)
        {
          if (remaining[value - 1] == 0) continue;
          if (left != 0 && static_cast<int>(value) < left) continue;
          if (above != 0 && static_cast<int>(value) <= above) continue;
          --remaining[value - 1];
          tableau[row][col] = static_cast<int>(value);
          fill(cell + 1);
          tableau[row][col] = 0;
          ++remaining[value - 1];
        }
    };
    fill(0);
    requireCacheEntryCapacity("Kostka cache");
    kostkaNumberCache[cacheKey] = result;
    return result;
  }

ring_elem SymmetricEngineRing::kostkaInnerProductForBasisElements(
    ring_elem schurTypeBasisElement,
    int schurTypeBasisId,
    ring_elem multiplicativeBasisElement,
    int multiplicativeBasisId,
    bool conjugateShape) const
{
    Partition shape;
    Partition content;
    ring_elem schurCoefficient;
    ring_elem multiplicativeCoefficient;
    if (!singleScaledBasisElement(
            schurTypeBasisElement,
            schurTypeBasisId,
            shape,
            schurCoefficient) ||
        !singleScaledBasisElement(
            multiplicativeBasisElement,
            multiplicativeBasisId,
            content,
            multiplicativeCoefficient))
      {
        ERROR("Kostka inner-product route requires two scaled basis elements");
        return coefficientRing->zero();
      }
    if (conjugateShape) shape = conjugatePartition(shape);
    mpz_class kostka =
        kostkaNumberViaSemistandardTableaux(shape, content);
    return coefficientRing->mult(
        coefficientRing->from_int(kostka.get_mpz_t()),
        coefficientRing->mult(schurCoefficient,
                              multiplicativeCoefficient));
  }

ring_elem SymmetricEngineRing::schurCompleteInnerProductViaKostkaNumbers(
    ring_elem schurBasisElement,
    ring_elem completeBasisElement) const
{
    return kostkaInnerProductForBasisElements(schurBasisElement,
                              requiredBasisIdForKind(BasisKind::Schur),
                              completeBasisElement,
                              requiredBasisIdForKind(BasisKind::Complete),
                              false);
  }

ring_elem
SymmetricEngineRing::schurElementaryInnerProductViaConjugateKostkaNumbers(
    ring_elem schurBasisElement,
    ring_elem elementaryBasisElement) const
{
    return kostkaInnerProductForBasisElements(schurBasisElement,
                              requiredBasisIdForKind(BasisKind::Schur),
                              elementaryBasisElement,
                              requiredBasisIdForKind(BasisKind::Elementary),
                              true);
  }

ring_elem
SymmetricEngineRing::schurOmegaCompleteInnerProductViaConjugateKostkaNumbers(
    ring_elem schurOmegaBasisElement,
    ring_elem completeBasisElement) const
{
    return kostkaInnerProductForBasisElements(schurOmegaBasisElement,
                              requiredBasisIdForKind(BasisKind::SchurOmega),
                              completeBasisElement,
                              requiredBasisIdForKind(BasisKind::Complete),
                              true);
  }

ring_elem
SymmetricEngineRing::schurOmegaElementaryInnerProductViaKostkaNumbers(
    ring_elem schurOmegaBasisElement,
    ring_elem elementaryBasisElement) const
{
    return kostkaInnerProductForBasisElements(schurOmegaBasisElement,
                              requiredBasisIdForKind(BasisKind::SchurOmega),
                              elementaryBasisElement,
                              requiredBasisIdForKind(BasisKind::Elementary),
                              false);
  }

// ============================================================================
// Power-Sum And Schur Pairings
// ============================================================================
// Weighted characters evaluate ordinary or Hall-Littlewood power-sum pairings with Schur expansions.

ring_elem
SymmetricEngineRing::powerSumsSchurInnerProductViaWeightedCharacters(
    ring_elem f,
    ring_elem g,
    const InnerProductContext& context) const
{
    int pId = requiredBasisIdForKind(BasisKind::PowerSum);
    int schurId = requiredBasisIdForKind(BasisKind::Schur);
    if (error()) return coefficientRing->zero();

    CoeffMap pCoeffs;
    CoeffMap schurCoeffs;
    if (!coefficientsInBasisIfPossible(f, pId, pCoeffs) ||
        !coefficientsInBasisIfPossible(g, schurId, schurCoeffs))
      {
        ERROR("weighted-character route requires power sums and a Schur expansion");
        return coefficientRing->zero();
      }

    ring_elem result = coefficientRing->zero();
    for (const auto& pTerm : pCoeffs)
      {
        int degree = partitionWeight(pTerm.first);
        Partition mu = normalizePartition(pTerm.first);
        for (const auto& schurTerm : schurCoeffs)
          {
            if (partitionWeight(schurTerm.first) != degree) continue;
            if (!isPartitionIndex(schurTerm.first))
              {
                ERROR("weighted-character route requires partition-indexed Schur elements");
                return coefficientRing->zero();
              }
            mpz_class chi =
                characterValueWithinLimits(schurTerm.first, mu);
            if (chi == 0) continue;
            ring_elem diagonal = powerSumDiagonalFactor(mu, context);
            if (error()) return coefficientRing->zero();
            ring_elem characterWeight = coefficientQuotient(
                diagonal, rationalCoefficient(zValue(mu), 1));
            ring_elem contribution = coefficientRing->mult(
                coefficientRing->mult(pTerm.second, schurTerm.second),
                characterWeight);
            if (error()) return coefficientRing->zero();
            if (chi != 1)
              contribution = coefficientRing->mult(
                  coefficientRing->from_int(chi.get_mpz_t()), contribution);
            result = coefficientRing->add(result, contribution);
          }
      }
    return result;
  }

} // namespace symmetric_rings
