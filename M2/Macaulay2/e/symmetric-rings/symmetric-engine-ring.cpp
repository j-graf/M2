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
#include <utility>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Basis Metadata
// ============================================================================

bool SymmetricEngineRing::isMultiplicativeBasis(int basisId) const
{
    return requireBasis(basisId).multiplicative;
  }

SymmetricEngineRing::BasisKind SymmetricEngineRing::basisKindForId(int basisId) const
{
    auto found = basisDescriptors.find(basisId);
    return found == basisDescriptors.end() ? BasisKind::Custom
                                           : found->second.kind;
  }

SymmetricEngineRing::BasisKind
SymmetricEngineRing::basisKindFromCanonicalKey(const std::string& key) const
{
    static const std::map<std::string, BasisKind> kinds = {
        {"PowerSum", BasisKind::PowerSum},
        {"Complete", BasisKind::Complete},
        {"Elementary", BasisKind::Elementary},
        {"Monomial", BasisKind::Monomial},
        {"Forgotten", BasisKind::Forgotten},
        {"Schur", BasisKind::Schur},
        {"SchurOmega", BasisKind::SchurOmega},
        {"HallLittlewoodQGenerator", BasisKind::HallLittlewoodQGenerator},
        {"HallLittlewoodBGenerator", BasisKind::HallLittlewoodBGenerator},
        {"HallLittlewoodQ", BasisKind::HallLittlewoodQ},
        {"HallLittlewoodB", BasisKind::HallLittlewoodB},
        {"HallLittlewoodP", BasisKind::HallLittlewoodP},
        {"HallLittlewoodPOmega", BasisKind::HallLittlewoodPOmega}};
    auto found = kinds.find(key);
    return found == kinds.end() ? BasisKind::Custom : found->second;
  }

bool SymmetricEngineRing::hasBasisKind(int basisId, BasisKind kind) const
{
    return basisKindForId(basisId) == kind;
  }

int SymmetricEngineRing::basisIdForKind(BasisKind kind) const
{
    auto found = basisIdsByKind.find(kind);
    return found == basisIdsByKind.end() ? -1 : found->second;
  }

int SymmetricEngineRing::registeredPowerSumBasisId() const
{
    return basisIdForKind(BasisKind::PowerSum);
  }

int SymmetricEngineRing::requiredBasisIdForKind(BasisKind kind) const
{
    int id = basisIdForKind(kind);
    if (id < 0)
      ERROR("basis metadata for ", basisKindName(kind), " is not available");
    return id;
  }

const char *SymmetricEngineRing::basisKindName(BasisKind kind) const
{
    switch (kind)
      {
        case BasisKind::PowerSum: return "PowerSum";
        case BasisKind::Complete: return "Complete";
        case BasisKind::Elementary: return "Elementary";
        case BasisKind::Monomial: return "Monomial";
        case BasisKind::Forgotten: return "Forgotten";
        case BasisKind::Schur: return "Schur";
        case BasisKind::SchurOmega: return "SchurOmega";
        case BasisKind::HallLittlewoodQGenerator: return "HallLittlewoodQGenerator";
        case BasisKind::HallLittlewoodBGenerator: return "HallLittlewoodBGenerator";
        case BasisKind::HallLittlewoodQ: return "HallLittlewoodQ";
        case BasisKind::HallLittlewoodB: return "HallLittlewoodB";
        case BasisKind::HallLittlewoodP: return "HallLittlewoodP";
        case BasisKind::HallLittlewoodPOmega: return "HallLittlewoodPOmega";
        case BasisKind::Custom: return "Custom";
      }
    return "Custom";
  }

const SymmetricEngineRing::BasisDescriptor&
SymmetricEngineRing::requireBasis(int basisId) const
{
    auto found = basisDescriptors.find(basisId);
    if (found != basisDescriptors.end()) return found->second;
    ERROR("metadata for symmetric-function basis id ", basisId,
          " is not available");
    static const BasisDescriptor missing{"", "", 0, false, BasisKind::Custom};
    return missing;
  }

void SymmetricEngineRing::rememberBasesFrom(const SymmetricEngineRing *R) const
{
    bool changed = false;
    for (const auto& item : R->basisDescriptors)
      if (basisDescriptors.find(item.first) == basisDescriptors.end())
        {
          basisDescriptors[item.first] = item.second;
          changed = true;
        }
    for (const auto& item : R->basisIdsByKind) basisIdsByKind[item.first] = item.second;
    if (changed)
      {
        basisConversionPlanRegistryCache.clear();
        multiplicationPlanRegistryCache.clear();
      }
  }

std::string SymmetricEngineRing::displayForBasis(int basisId) const
{
    return requireBasis(basisId).displaySymbol;
  }

std::string SymmetricEngineRing::basisKeyForId(int basisId) const
{
    return requireBasis(basisId).canonicalKey;
  }

int SymmetricEngineRing::basisOrderForId(int basisId) const
{
    return requireBasis(basisId).displayOrder;
  }

// ============================================================================
// Shared Character, Coefficient, And Hall-Littlewood State
// ============================================================================

const CharacterTable& SymmetricEngineRing::characterTable(int degree) const
{
    auto cached = characterTableCache.find(degree);
    if (cached != characterTableCache.end()) return cached->second;

    CharacterTable table;
    table.partitions = partitionsOfWithinLimits(degree, "character table");
    if (error())
      {
        static const CharacterTable empty;
        return empty;
      }
    table.zValues.reserve(table.partitions.size());
    for (size_t i = 0; i < table.partitions.size(); ++i)
      {
        table.partitionRows[table.partitions[i]] = i;
        table.zValues.push_back(zValue(table.partitions[i]));
      }

    requireCacheEntryCapacity("character-table cache");
    auto inserted = characterTableCache.emplace(degree, std::move(table));
    return inserted.first->second;
  }

mpz_class SymmetricEngineRing::characterTableValue(
    const CharacterTable& table, size_t row, size_t col) const
{
    const auto key = std::make_pair(row, col);
    auto cached = table.values.find(key);
    if (cached != table.values.end()) return cached->second;
    if (characterCacheEntryCount >= computationLimits.maxCharacterCacheEntries)
      throw exc::engine_error(
          "character cache entry limit exceeded; increase MaxCharacterCacheEntries in the symmetricRing ComputationLimits option");
    bool recursiveLimitExceeded = false;
    mpz_class value = characterValueWithLimit(
        table.partitions[row],
        table.partitions[col],
        computationLimits.maxRecursiveStates,
        recursiveLimitExceeded);
    if (recursiveLimitExceeded)
      throw exc::engine_error(
          "character recursion state limit exceeded; increase MaxRecursiveStates in the symmetricRing ComputationLimits option");
    requireCacheEntryCapacity("character-value cache");
    table.values.emplace(key, value);
    ++characterCacheEntryCount;
    return value;
  }

mpz_class SymmetricEngineRing::characterValueWithinLimits(
    const Partition& lambda, const Partition& mu) const
{
    bool recursiveLimitExceeded = false;
    mpz_class value = characterValueWithLimit(
        lambda,
        mu,
        computationLimits.maxRecursiveStates,
        recursiveLimitExceeded);
    if (recursiveLimitExceeded)
      throw exc::engine_error(
          "character recursion state limit exceeded; increase MaxRecursiveStates in the symmetricRing ComputationLimits option");
    return value;
  }

mpz_class SymmetricEngineRing::pToMonomialCoefficientWithinLimits(
    const Partition& lambda, const Partition& mu) const
{
    bool recursiveLimitExceeded = false;
    mpz_class value = pToMonomialCoefficientWithLimit(
        lambda,
        mu,
        computationLimits.maxRecursiveStates,
        recursiveLimitExceeded);
    if (recursiveLimitExceeded)
      throw exc::engine_error(
          "power-sum-to-monomial recursion state limit exceeded; increase MaxRecursiveStates in the symmetricRing ComputationLimits option");
    return value;
  }

bool SymmetricEngineRing::estimatedMemoryWithinLimit(
    size_t count, size_t bytesPerItem, const char *operation) const
{
    if (bytesPerItem != 0 &&
        count > computationLimits.maxEstimatedMemoryBytes / bytesPerItem)
      throw exc::engine_error(
          std::string(operation) +
          " exceeds the estimated memory limit; increase MaxEstimatedMemoryMB in the symmetricRing ComputationLimits option");
    return true;
  }

std::vector<Partition> SymmetricEngineRing::partitionsOfWithinLimits(
    int degree, const char *operation) const
{
    const size_t count =
        partitionCountUpToLimit(degree, computationLimits.maxEnumeratedPartitions);
    if (count > computationLimits.maxEnumeratedPartitions)
      throw exc::engine_error(
          std::string(operation) +
          " exceeds the partition-enumeration limit; increase MaxEnumeratedPartitions in the symmetricRing ComputationLimits option");
    // Vector nodes, integer payloads, and allocator overhead vary by platform.
    // Sixty-four bytes per partition is a deliberately conservative preflight.
    if (!estimatedMemoryWithinLimit(count, 64, operation)) return {};
    return partitionsOf(degree);
  }

bool SymmetricEngineRing::determinantStatesWithinLimit(
    size_t states, const char *operation) const
{
    if (states > computationLimits.maxDeterminantStates)
      throw exc::engine_error(
          std::string(operation) +
          " exceeds the determinant-state limit; increase MaxDeterminantStates in the symmetricRing ComputationLimits option");
    return estimatedMemoryWithinLimit(states, 2 * sizeof(ring_elem), operation);
  }

bool SymmetricEngineRing::consumeRecursiveState(
    size_t& states, const char *operation) const
{
    if (states >= computationLimits.maxRecursiveStates)
      throw exc::engine_error(
          std::string(operation) +
          " exceeds the recursion-state limit; increase MaxRecursiveStates in the symmetricRing ComputationLimits option");
    ++states;
    return true;
  }

void SymmetricEngineRing::requireCacheEntryCapacity(
    const char *operation) const
{
    if (computationCacheEntryCount >= computationLimits.maxCacheEntries)
      throw exc::engine_error(
          std::string(operation) +
          " exceeds MaxCacheEntries in the symmetricRing ComputationLimits option");
    ++computationCacheEntryCount;
  }

void SymmetricEngineRing::requireWeightWithinLimit(
    long long weight, const char *operation) const
{
    const unsigned long long magnitude =
        weight < 0
            ? static_cast<unsigned long long>(-(weight + 1)) + 1
            : static_cast<unsigned long long>(weight);
    if (magnitude > computationLimits.maxWeight)
      throw exc::engine_error(
          std::string(operation) +
          " exceeds MaxWeight in the symmetricRing ComputationLimits option");
  }

void SymmetricEngineRing::requirePartitionWithinWeightLimit(
    const Partition& index, const char *operation) const
{
    if (index.size() > computationLimits.maxWeight)
      throw exc::engine_error(
          std::string(operation) +
          " exceeds MaxWeight in the symmetricRing ComputationLimits option");
    unsigned long long complexity = 0;
    for (int part : index)
      {
        const long long widePart = part;
        const unsigned long long magnitude =
            widePart < 0
                ? static_cast<unsigned long long>(-(widePart + 1)) + 1
                : static_cast<unsigned long long>(widePart);
        if (magnitude > computationLimits.maxWeight - complexity)
          throw exc::engine_error(
              std::string(operation) +
              " exceeds MaxWeight in the symmetricRing ComputationLimits option");
        complexity += magnitude;
      }
  }

void SymmetricEngineRing::requireMonomialWithinWeightLimit(
    const SymmetricMonomial& monomial, const char *operation) const
{
    requireWeightWithinLimit(monomialWeight(monomial), operation);
  }

ring_elem SymmetricEngineRing::rationalCoefficient(long numerator, long denominator) const
{
    return rationalCoefficient(
        mpz_class(numerator), mpz_class(denominator));
  }

ring_elem SymmetricEngineRing::rationalCoefficient(
    const mpz_class& numerator,
    const mpz_class& denominator) const
{
    mpq_t q;
    mpq_init(q);
    mpz_set(mpq_numref(q), numerator.get_mpz_t());
    mpz_set(mpq_denref(q), denominator.get_mpz_t());
    mpq_canonicalize(q);
    ring_elem result;
    if (!coefficientRing->from_rational(q, result))
      {
        ERROR("coefficient division failed during basis conversion; use a coefficient ring where the required denominators are invertible, for example frac(QQ[t]) instead of QQ[t]");
        result = coefficientRing->zero();
      }
    mpq_clear(q);
    return result;
  }

void SymmetricEngineRing::clearHallLittlewoodCaches() const
{
    hallLittlewoodQGeneratorToPowerSumsCache.clear();
    hallLittlewoodBGeneratorToPowerSumsCache.clear();
    hallLittlewoodQGeneratorToPowerSumsQuotientMapCache.clear();
    hallLittlewoodBGeneratorToPowerSumsQuotientMapCache.clear();
    hallLittlewoodPartFactorCache.clear();
    hallLittlewoodCFactorCache.clear();
    hallLittlewoodSingleCycleCapitalGreenMapCache.clear();
    hallLittlewoodSingleCycleNormalizedGreenMapCache.clear();
    hallLittlewoodPowerSumToCapitalColumnCache.clear();
    hallLittlewoodPowerSumToCapitalCoefficientCache.clear();
    hallLittlewoodRaisingGeneratorMapCache.clear();
    powerSumToQGeneratorMapCache.clear();
    powerSumToBGeneratorMapCache.clear();
  }

ring_elem SymmetricEngineRing::hallLittlewoodFactor(const Partition& mu) const
{
    ring_elem result = coefficientRing->one();
    for (int part : mu)
      {
        auto cached = hallLittlewoodPartFactorCache.find(part);
        if (cached == hallLittlewoodPartFactorCache.end())
          {
            ring_elem tPower = coefficientRing->power(hallLittlewoodParameter, part);
            ring_elem oneMinus = coefficientRing->subtract(
                coefficientRing->one(), tPower);
            requireCacheEntryCapacity("Hall-Littlewood factor cache");
            cached = hallLittlewoodPartFactorCache.emplace(part, oneMinus).first;
          }
        result = coefficientRing->mult(result, cached->second);
      }
    return result;
  }

// ============================================================================
// Basis-Element Construction
// ============================================================================

ring_elem SymmetricEngineRing::basisElementFromIndex(
    int basisId,
    const Partition& index) const
{
    requirePartitionWithinWeightLimit(index, "basis element");
    const auto& basis = requireBasis(basisId);
    auto result = new SymmetricRingPoly;
    SymmetricMonomial monomial;
    appendAtomBlock(
        monomial, makeAtomBlock(basis.displayOrder, basisId, 0, index));
    result->terms.push_back({coefficientRing->one(), canonicalMonomial(monomial)});
    return makePolyValue(result);
  }

ring_elem SymmetricEngineRing::basisElementFromSkewIndex(
    int basisId,
    const Partition& outer,
    const Partition& inner) const
{
    requirePartitionWithinWeightLimit(outer, "skew outer index");
    requirePartitionWithinWeightLimit(inner, "skew inner index");
    requireWeightWithinLimit(
        static_cast<long long>(partitionWeight(outer)) -
            static_cast<long long>(partitionWeight(inner)),
        "skew basis element");
    const auto& basis = requireBasis(basisId);
    Partition payload = outer;
    payload.insert(payload.end(), inner.begin(), inner.end());
    auto result = new SymmetricRingPoly;
    SymmetricMonomial monomial;
    appendAtomBlock(monomial,
                    makeAtomBlock(basis.displayOrder,
                                  basisId,
                                  static_cast<int>(inner.size()),
                                  payload));
    result->terms.push_back({coefficientRing->one(), canonicalMonomial(monomial)});
    return makePolyValue(result);
  }

ring_elem SymmetricEngineRing::basisPartElement(int basisId, int n) const
{
    if (n < 0) return zero();
    if (n == 0) return one();
    return basisElementFromIndex(basisId, Partition{n});
  }

// ============================================================================
// Ring Lifecycle And Public Metadata
// ============================================================================

SymmetricEngineRing::SymmetricEngineRing(const Ring *A)
      : coefficientRing(A), hallLittlewoodParameter(A->from_long(0))
{}

SymmetricEngineRing *SymmetricEngineRing::create(const Ring *A)
{
    auto result = new SymmetricEngineRing(A);
    result->initialize_ring(A->characteristic());
    result->zeroV = result->from_long(0);
    result->oneV = result->from_long(1);
    result->minus_oneV = result->from_long(-1);
    return result;
  }

const Ring *SymmetricEngineRing::getCoefficientRing() const
{ return coefficientRing; }

void SymmetricEngineRing::rememberBasisMetadata(int basisId,
                             const std::string& canonicalBasisKey,
                             const std::string& display,
                             int order,
                             bool isMultiplicative) const
{
    BasisKind kind = basisKindFromCanonicalKey(canonicalBasisKey);
    BasisDescriptor descriptor{
        canonicalBasisKey, display, order, isMultiplicative, kind};
    auto found = basisDescriptors.find(basisId);
    if (found != basisDescriptors.end())
      {
        const auto& existing = found->second;
        if (existing.canonicalKey != descriptor.canonicalKey ||
            existing.displaySymbol != descriptor.displaySymbol ||
            existing.displayOrder != descriptor.displayOrder ||
            existing.multiplicative != descriptor.multiplicative ||
            existing.kind != descriptor.kind)
          {
            ERROR("conflicting metadata for symmetric-function basis id ",
                  basisId);
            return;
          }
      }
    else
      {
        basisDescriptors.emplace(basisId, std::move(descriptor));
        basisConversionPlanRegistryCache.clear();
        multiplicationPlanRegistryCache.clear();
      }
    if (kind != BasisKind::Custom) basisIdsByKind[kind] = basisId;
  }

bool SymmetricEngineRing::setHallLittlewoodParameter(const RingElement *t) const
{
    ring_elem promoted;
    if (t->get_ring() == coefficientRing)
      promoted = t->get_value();
    else if (!coefficientRing->promote(t->get_ring(), t->get_value(), promoted))
      return false;
    hallLittlewoodParameter = promoted;
    clearHallLittlewoodCaches();
    return true;
  }

void SymmetricEngineRing::setComputationLimits(
    size_t maxWeight,
    size_t maxEnumeratedPartitions,
    size_t maxGeneratedTerms,
    size_t maxRecursiveStates,
    size_t maxCacheEntries,
    size_t maxCharacterCacheEntries,
    size_t maxDeterminantStates,
    size_t maxEstimatedMemoryMB) const
{
    computationLimits.maxWeight = maxWeight;
    computationLimits.maxEnumeratedPartitions = maxEnumeratedPartitions;
    computationLimits.maxGeneratedTerms = maxGeneratedTerms;
    computationLimits.maxRecursiveStates = maxRecursiveStates;
    computationLimits.maxCacheEntries = maxCacheEntries;
    computationLimits.maxCharacterCacheEntries = maxCharacterCacheEntries;
    computationLimits.maxDeterminantStates = maxDeterminantStates;
    if (maxEstimatedMemoryMB >
        std::numeric_limits<size_t>::max() / (1024ULL * 1024ULL))
      computationLimits.maxEstimatedMemoryBytes =
          std::numeric_limits<size_t>::max();
    else
      computationLimits.maxEstimatedMemoryBytes =
          maxEstimatedMemoryMB * 1024ULL * 1024ULL;
  }

int SymmetricEngineRing::uniformBasisId(ring_elem f) const
{
    return singleBasisId(f);
  }

const SymmetricEngineRing *symmetricRingFromElement(const RingElement *f)
{
  auto result = dynamic_cast<const SymmetricEngineRing *>(f->get_ring());
  if (result == nullptr) ERROR("expected a SymmetricRings engine element");
  return result;
}

const SymmetricEngineRing *symmetricRingFromRing(const Ring *R)
{
  auto result = dynamic_cast<const SymmetricEngineRing *>(R);
  if (result == nullptr) ERROR("expected a SymmetricRings engine ring");
  return result;
}

} // namespace symmetric_rings
