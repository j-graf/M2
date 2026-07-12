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

bool SymmetricEngineRing::isMultiplicativeBasis(int basisId) const
{
    auto it = multiplicativeBases.find(basisId);
    return it != multiplicativeBases.end() && it->second;
  }

bool SymmetricEngineRing::isPowerSumBasis(int basisId) const
{
    return powerSumBasisId >= 0 && basisId == powerSumBasisId;
  }

SymmetricEngineRing::BasisKind SymmetricEngineRing::basisKindForId(int basisId) const
{
    auto found = basisKinds.find(basisId);
    return found == basisKinds.end() ? BasisKind::Custom : found->second;
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

void SymmetricEngineRing::rememberBasis(int basisId,
                     const std::string& display,
                     int order,
                     bool isMultiplicative) const
{
    if (!display.empty()) basisDisplays[basisId] = display;
    basisOrders[basisId] = order;
    multiplicativeBases[basisId] = isMultiplicative;
  }

void SymmetricEngineRing::rememberBasesFrom(const SymmetricEngineRing *R) const
{
    for (const auto& item : R->basisDisplays) basisDisplays[item.first] = item.second;
    for (const auto& item : R->basisKeys) basisKeys[item.first] = item.second;
    for (const auto& item : R->basisOrders) basisOrders[item.first] = item.second;
    for (const auto& item : R->multiplicativeBases)
      multiplicativeBases[item.first] = item.second;
    for (const auto& item : R->basisKinds) basisKinds[item.first] = item.second;
    for (const auto& item : R->basisIdsByKind) basisIdsByKind[item.first] = item.second;
    if (powerSumBasisId < 0) powerSumBasisId = R->powerSumBasisId;
  }

std::string SymmetricEngineRing::displayForBasis(int basisId) const
{
    auto it = basisDisplays.find(basisId);
    if (it != basisDisplays.end()) return it->second;
    return "basis" + std::to_string(basisId);
  }

std::string SymmetricEngineRing::basisKeyForId(int basisId) const
{
    auto it = basisKeys.find(basisId);
    if (it != basisKeys.end()) return it->second;
    ERROR("stable key metadata for symmetric-function basis id ", basisId,
          " is not available");
    return "";
  }

int SymmetricEngineRing::basisOrderForId(int basisId) const
{
    auto it = basisOrders.find(basisId);
    if (it != basisOrders.end()) return it->second;
    ERROR("display metadata for symmetric-function basis id ", basisId,
          " is not available");
    return 0;
  }

const CharacterTable& SymmetricEngineRing::characterTable(int degree) const
{
    auto cached = characterTableCache.find(degree);
    if (cached != characterTableCache.end()) return cached->second;

    CharacterTable table;
    table.partitions = partitionsOf(degree);
    table.zValues.reserve(table.partitions.size());
    for (size_t i = 0; i < table.partitions.size(); ++i)
      {
        table.partitionRows[table.partitions[i]] = i;
        table.zValues.push_back(zValue(table.partitions[i]));
      }

    table.values.resize(
        table.partitions.size(),
        std::vector<int>(table.partitions.size(), unknownCharacterValue));

    auto inserted = characterTableCache.emplace(degree, std::move(table));
    return inserted.first->second;
  }

int SymmetricEngineRing::characterTableValue(const CharacterTable& table, size_t row, size_t col) const
{
    int& value = table.values[row][col];
    if (value == unknownCharacterValue)
      value = characterValue(table.partitions[row], table.partitions[col]);
    return value;
  }

ring_elem SymmetricEngineRing::rationalCoefficient(long numerator, long denominator) const
{
    mpq_t q;
    mpq_init(q);
    mpq_set_si(q, numerator, denominator);
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
            cached = hallLittlewoodPartFactorCache.emplace(part, oneMinus).first;
          }
        result = coefficientRing->mult(result, cached->second);
      }
    return result;
  }

ring_elem SymmetricEngineRing::basisElementFromIndex(int basisId,
                                  const std::string& display,
                                  int order,
                                  bool isMultiplicative,
                                  const Partition& index) const
{
    rememberBasis(basisId, display, order, isMultiplicative);
    auto result = new SymmetricRingPoly;
    SymmetricMonomial monomial;
    appendAtomBlock(monomial, makeAtomBlock(order, basisId, 0, index));
    result->terms.push_back({coefficientRing->one(), canonicalMonomial(monomial)});
    return makePolyValue(result);
  }

ring_elem SymmetricEngineRing::basisElementFromSkewIndex(int basisId,
                                      const std::string& display,
                                      int order,
                                      bool isMultiplicative,
                                      const Partition& outer,
                                      const Partition& inner) const
{
    rememberBasis(basisId, display, order, isMultiplicative);
    Partition payload = outer;
    payload.insert(payload.end(), inner.begin(), inner.end());
    auto result = new SymmetricRingPoly;
    SymmetricMonomial monomial;
    appendAtomBlock(monomial,
                    makeAtomBlock(order,
                                  basisId,
                                  static_cast<int>(inner.size()),
                                  payload));
    result->terms.push_back({coefficientRing->one(), canonicalMonomial(monomial)});
    return makePolyValue(result);
  }

ring_elem SymmetricEngineRing::basisPartElement(int basisId,
                             const std::string& display,
                             int order,
                             bool isMultiplicative,
                             int n) const
{
    if (n < 0) return zero();
    if (n == 0) return one();
    return basisElementFromIndex(basisId, display, order, isMultiplicative, Partition{n});
  }

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
    rememberBasis(basisId, display, order, isMultiplicative);
    basisKeys[basisId] = canonicalBasisKey;
    BasisKind kind = basisKindFromCanonicalKey(canonicalBasisKey);
    basisKinds[basisId] = kind;
    if (kind != BasisKind::Custom) basisIdsByKind[kind] = basisId;
    if (kind == BasisKind::PowerSum) powerSumBasisId = basisId;
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
