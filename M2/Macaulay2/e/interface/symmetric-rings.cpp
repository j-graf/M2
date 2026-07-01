// Copyright 2026

#include "interface/symmetric-rings.h"

#include "buffer.hpp"
#include "error.h"
#include "exceptions.hpp"
#include "newdelete.hpp"
#include "ring-elements/ring-element.hpp"
#include "ringmap.hpp"
#include "rings/ZZ.hpp"
#include "rings/ring.hpp"
#include "rings/ringelem.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr size_t atomHeaderSize = 4;

struct SymmetricMonomial
{
  // Concatenated atom blocks:
  // [displayOrder, basisId, outerLength, innerLength, payload_1, ..., payload_n, ...]
  // For skew atoms, payload is outer shape followed by inner shape.
  VECTOR(int) data;
};

struct SymmetricTerm
{
  ring_elem coeff;
  SymmetricMonomial monomial;
};

class SymmetricRingPoly : public our_new_delete
{
 public:
  VECTOR(SymmetricTerm) terms;
};

struct BasisIndexKey
{
  int basisId;
  std::vector<int> index;

  bool operator<(const BasisIndexKey& other) const
  {
    if (basisId != other.basisId) return basisId < other.basisId;
    return index < other.index;
  }
};

using Partition = std::vector<int>;

template <typename K, typename V, typename Compare = std::less<K>>
using GCMap = std::map<K, V, Compare, gc_allocator<std::pair<const K, V>>>;

using CoeffMap = GCMap<Partition, ring_elem>;
using RingElemVector = VECTOR(ring_elem);
using RingElemMatrix = VECTOR(RingElemVector);

constexpr int unknownCharacterValue = std::numeric_limits<int>::min();

struct CharacterTable
{
  std::vector<Partition> partitions;
  std::vector<long> zValues;
  std::map<Partition, size_t> partitionRows;
  mutable std::vector<std::vector<int>> values;
};

struct SchurConversionRecipeEntry
{
  Partition lambda;
  std::vector<std::pair<size_t, int>> contributions;
};

struct LRProductTerm
{
  Partition nu;
  long coefficient;
};

struct OmegaTarget
{
  int basisId;
  int order;
  bool isMultiplicative;
};

struct InnerProductTarget
{
  int dualBasisId;
  int kind;
};

std::string partitionKey(const Partition& p)
{
  std::ostringstream out;
  for (size_t i = 0; i < p.size(); ++i)
    {
      if (i > 0) out << ",";
      out << p[i];
    }
  return out.str();
}

Partition normalizePartition(const Partition& p)
{
  Partition result;
  for (int part : p)
    if (part > 0) result.push_back(part);
  std::sort(result.begin(), result.end(), std::greater<int>());
  return result;
}

int partitionWeight(const Partition& p)
{
  int result = 0;
  for (int part : p) result += part;
  return result;
}

int partitionLength(const Partition& p)
{
  int result = 0;
  for (int part : p)
    if (part > 0) ++result;
  return result;
}

bool lexLessPartition(const Partition& a, const Partition& b)
{
  size_t n = std::max(a.size(), b.size());
  for (size_t i = 0; i < n; ++i)
    {
      int aa = i < a.size() ? a[i] : 0;
      int bb = i < b.size() ? b[i] : 0;
      if (aa < bb) return true;
      if (aa > bb) return false;
    }
  return false;
}

Partition trimTrailingZerosPartition(const Partition& p)
{
  Partition result = p;
  while (!result.empty() && result.back() == 0) result.pop_back();
  return result;
}

Partition conjugatePartition(const Partition& p)
{
  Partition lambda = normalizePartition(p);
  Partition result;
  if (lambda.empty()) return result;
  int maxPart = lambda.front();
  result.reserve(maxPart);
  for (int col = 1; col <= maxPart; ++col)
    {
      int count = 0;
      for (int part : lambda)
        if (part >= col) ++count;
      result.push_back(count);
    }
  return trimTrailingZerosPartition(result);
}

std::pair<int, Partition> straightenSchurIndex(const Partition& alpha)
{
  Partition trimmed = trimTrailingZerosPartition(alpha);
  size_t ell = trimmed.size();
  if (ell == 0) return {1, Partition{}};

  std::vector<int> shifted;
  shifted.reserve(ell);
  for (size_t i = 0; i < ell; ++i)
    shifted.push_back(trimmed[i] + static_cast<int>(ell - 1 - i));

  std::vector<int> sortedShifted = shifted;
  std::sort(sortedShifted.begin(), sortedShifted.end());
  for (size_t i = 1; i < sortedShifted.size(); ++i)
    if (sortedShifted[i] == sortedShifted[i - 1]) return {0, Partition{}};

  int inversions = 0;
  for (size_t i = 0; i + 1 < ell; ++i)
    for (size_t j = i + 1; j < ell; ++j)
      if (shifted[i] < shifted[j]) ++inversions;

  std::sort(sortedShifted.begin(), sortedShifted.end(), std::greater<int>());
  Partition beta;
  beta.reserve(ell);
  for (size_t i = 0; i < ell; ++i)
    beta.push_back(sortedShifted[i] - static_cast<int>(ell - 1 - i));
  return {(inversions % 2 == 0) ? 1 : -1, trimTrailingZerosPartition(beta)};
}

void partitionsRec(int n, int maxPart, Partition& current, std::vector<Partition>& result)
{
  if (n == 0)
    {
      result.push_back(current);
      return;
    }
  for (int part = std::min(n, maxPart); part >= 1; --part)
    {
      current.push_back(part);
      partitionsRec(n - part, part, current, result);
      current.pop_back();
    }
}

std::vector<Partition> partitionsOf(int n)
{
  std::vector<Partition> result;
  Partition current;
  partitionsRec(n, n, current, result);
  return result;
}

long assignmentCountRec(const Partition& parts, size_t pos, std::vector<int>& targets)
{
  if (pos == parts.size())
    {
      for (int target : targets)
        if (target != 0) return 0;
      return 1;
    }
  long total = 0;
  int part = parts[pos];
  for (size_t i = 0; i < targets.size(); ++i)
    if (targets[i] >= part)
      {
        targets[i] -= part;
        total += assignmentCountRec(parts, pos + 1, targets);
        targets[i] += part;
      }
  return total;
}

long pToMonomialCoefficient(const Partition& lambda, const Partition& mu)
{
  Partition normalizedLambda = normalizePartition(lambda);
  Partition normalizedMu = normalizePartition(mu);
  if (partitionWeight(normalizedLambda) != partitionWeight(normalizedMu)) return 0;
  std::vector<int> targets = normalizedMu;
  return assignmentCountRec(normalizedLambda, 0, targets);
}

long zValue(const Partition& lambda)
{
  std::map<int, int> multiplicities;
  for (int part : lambda) multiplicities[part]++;
  long result = 1;
  for (const auto& item : multiplicities)
    {
      for (int i = 0; i < item.second; ++i) result *= item.first;
      for (int i = 2; i <= item.second; ++i) result *= i;
    }
  return result;
}

bool isPartitionAfterRemoval(const Partition& lambda, const std::vector<int>& removed)
{
  int previous = -1;
  for (size_t i = 0; i < lambda.size(); ++i)
    {
      int remaining = lambda[i] - removed[i];
      if (remaining < 0) return false;
      if (previous >= 0 && remaining > previous) return false;
      previous = remaining;
    }
  return true;
}

Partition remainingPartition(const Partition& lambda, const std::vector<int>& removed)
{
  Partition result;
  for (size_t i = 0; i < lambda.size(); ++i)
    {
      int remaining = lambda[i] - removed[i];
      if (remaining > 0) result.push_back(remaining);
    }
  return result;
}

Partition partitionFromM2Array(M2_arrayint a)
{
  Partition result;
  if (a == nullptr) return result;
  result.reserve(a->len);
  for (int i = 0; i < a->len; ++i) result.push_back(a->array[i]);
  return result;
}

bool removedCell(const Partition& lambda,
                 const std::vector<int>& removed,
                 int row,
                 int col)
{
  if (row < 0 || static_cast<size_t>(row) >= lambda.size()) return false;
  if (col < 1 || col > lambda[row]) return false;
  return col > lambda[row] - removed[row];
}

bool isConnectedSkew(const Partition& lambda, const std::vector<int>& removed)
{
  int total = 0;
  int startRow = -1;
  int startCol = -1;
  for (size_t r = 0; r < lambda.size(); ++r)
    for (int c = lambda[r] - removed[r] + 1; c <= lambda[r]; ++c)
      {
        ++total;
        if (startRow < 0)
          {
            startRow = static_cast<int>(r);
            startCol = c;
          }
      }
  if (total == 0) return false;

  std::vector<std::pair<int, int>> stack{{startRow, startCol}};
  std::map<std::pair<int, int>, bool> seen;
  seen[stack.back()] = true;
  int visited = 0;
  while (!stack.empty())
    {
      auto cell = stack.back();
      stack.pop_back();
      ++visited;
      const int dr[4] = {1, -1, 0, 0};
      const int dc[4] = {0, 0, 1, -1};
      for (int i = 0; i < 4; ++i)
        {
          std::pair<int, int> next{cell.first + dr[i], cell.second + dc[i]};
          if (!seen[next] && removedCell(lambda, removed, next.first, next.second))
            {
              seen[next] = true;
              stack.push_back(next);
            }
        }
    }
  return visited == total;
}

bool hasNoTwoByTwo(const Partition& lambda, const std::vector<int>& removed)
{
  for (size_t r = 0; r + 1 < lambda.size(); ++r)
    for (int c = 1; c <= std::max(lambda[r], lambda[r + 1]); ++c)
      if (removedCell(lambda, removed, static_cast<int>(r), c) &&
          removedCell(lambda, removed, static_cast<int>(r) + 1, c) &&
          removedCell(lambda, removed, static_cast<int>(r), c + 1) &&
          removedCell(lambda, removed, static_cast<int>(r) + 1, c + 1))
        return false;
  return true;
}

struct RimHookRemoval
{
  Partition remaining;
  int height;
};

void rimHookRec(const Partition& lambda,
                int row,
                int remainingSize,
                std::vector<int>& removed,
                std::vector<RimHookRemoval>& result)
{
  if (row == static_cast<int>(lambda.size()))
    {
      if (remainingSize != 0) return;
      if (!isPartitionAfterRemoval(lambda, removed)) return;
      if (!isConnectedSkew(lambda, removed)) return;
      if (!hasNoTwoByTwo(lambda, removed)) return;
      int rows = 0;
      for (int count : removed)
        if (count > 0) ++rows;
      result.push_back({remainingPartition(lambda, removed), rows - 1});
      return;
    }
  for (int count = 0; count <= std::min(lambda[row], remainingSize); ++count)
    {
      removed[row] = count;
      rimHookRec(lambda, row + 1, remainingSize - count, removed, result);
    }
  removed[row] = 0;
}

std::vector<RimHookRemoval> rimHookRemovals(const Partition& lambda, int size)
{
  std::vector<RimHookRemoval> result;
  std::vector<int> removed(lambda.size(), 0);
  rimHookRec(lambda, 0, size, removed, result);
  return result;
}

int characterValueMemo(const Partition& lambda,
                       const Partition& mu,
                       std::map<std::string, int>& memo)
{
  if (mu.empty()) return lambda.empty() ? 1 : 0;
  std::string key = partitionKey(lambda) + "|" + partitionKey(mu);
  auto it = memo.find(key);
  if (it != memo.end()) return it->second;

  int total = 0;
  int hookSize = mu.front();
  Partition rest(mu.begin() + 1, mu.end());
  for (const auto& removal : rimHookRemovals(lambda, hookSize))
    {
      int sign = (removal.height % 2 == 0) ? 1 : -1;
      total += sign * characterValueMemo(removal.remaining, rest, memo);
    }
  memo[key] = total;
  return total;
}

int characterValue(const Partition& lambda, const Partition& mu)
{
  static std::map<std::string, int> memo;
  return characterValueMemo(lambda, mu, memo);
}

std::string fromM2String(M2_string s)
{
  if (s == nullptr) return "";
  return std::string(reinterpret_cast<const char *>(s->array), s->len);
}

M2_string toM2String(const std::string& s)
{
  return M2_tostring(s.c_str());
}

std::string join(const std::vector<std::string>& parts,
                 const std::string& delimiter)
{
  std::ostringstream out;
  for (size_t i = 0; i < parts.size(); ++i)
    {
      if (i > 0) out << delimiter;
      out << parts[i];
    }
  return out.str();
}

const SymmetricRingPoly *polyValue(ring_elem f)
{
  return reinterpret_cast<const SymmetricRingPoly *>(f.get_Poly());
}

ring_elem makePolyValue(SymmetricRingPoly *f)
{
  return ring_elem(reinterpret_cast<const void *>(f));
}

template <typename T>
void hashCombine(size_t& seed, const T& val)
{
  seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

std::string coeffToString(const Ring *R, ring_elem c)
{
  buffer o;
  R->elem_text_out(o, c);
  return std::string(o.str());
}

size_t atomLengthAt(const SymmetricMonomial& monomial, size_t pos)
{
  if (pos + 3 >= monomial.data.size()) return monomial.data.size() - pos;
  return atomHeaderSize + monomial.data[pos + 2] + monomial.data[pos + 3];
}

int atomOrderAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos];
}

int atomBasisIdAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos + 1];
}

bool atomIsSkewAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos + 3] > 0;
}

int atomOuterLengthAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos + 2];
}

int atomInnerLengthAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos + 3];
}

int atomIndexLengthAt(const SymmetricMonomial& monomial, size_t pos)
{
  return atomOuterLengthAt(monomial, pos) + atomInnerLengthAt(monomial, pos);
}

bool blockLess(const std::vector<int>& a, const std::vector<int>& b)
{
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

bool monomialLess(const SymmetricMonomial& a, const SymmetricMonomial& b)
{
  return std::lexicographical_compare(
      a.data.begin(), a.data.end(), b.data.begin(), b.data.end());
}

int compareMonomials(const SymmetricMonomial& a, const SymmetricMonomial& b)
{
  if (monomialLess(a, b)) return LT;
  if (monomialLess(b, a)) return GT;
  return EQ;
}

std::vector<int> atomBlockAt(const SymmetricMonomial& monomial, size_t pos)
{
  size_t len = atomLengthAt(monomial, pos);
  return std::vector<int>(monomial.data.begin() + pos,
                          monomial.data.begin() + pos + len);
}

void appendAtomBlock(SymmetricMonomial& monomial, const std::vector<int>& block)
{
  monomial.data.insert(monomial.data.end(), block.begin(), block.end());
}

std::vector<int> makeAtomBlock(int displayOrder,
                               int basisId,
                               int innerLength,
                               M2_arrayint index)
{
  std::vector<int> block;
  int payloadLength = index == nullptr ? 0 : index->len;
  block.push_back(displayOrder);
  block.push_back(basisId);
  block.push_back(payloadLength - innerLength);
  block.push_back(innerLength);
  if (index != nullptr)
    for (int i = 0; i < index->len; ++i) block.push_back(index->array[i]);
  return block;
}

std::vector<int> makeAtomBlock(int displayOrder,
                               int basisId,
                               int innerLength,
                               const std::vector<int>& index)
{
  std::vector<int> block;
  block.reserve(atomHeaderSize + index.size());
  block.push_back(displayOrder);
  block.push_back(basisId);
  block.push_back(static_cast<int>(index.size()) - innerLength);
  block.push_back(innerLength);
  block.insert(block.end(), index.begin(), index.end());
  return block;
}

std::vector<int> monomialKey(const SymmetricMonomial& monomial)
{
  return std::vector<int>(monomial.data.begin(), monomial.data.end());
}

SymmetricMonomial monomialFromKey(const std::vector<int>& key)
{
  SymmetricMonomial result;
  result.data.insert(result.data.end(), key.begin(), key.end());
  return result;
}

int monomialWeight(const SymmetricMonomial& monomial)
{
  int result = 0;
  size_t pos = 0;
  while (pos < monomial.data.size())
    {
      int n = atomIndexLengthAt(monomial, pos);
      if (atomIsSkewAt(monomial, pos))
        {
          int outerLength = atomOuterLengthAt(monomial, pos);
          int innerLength = atomInnerLengthAt(monomial, pos);
          if (outerLength + innerLength == n)
            {
              for (int i = 0; i < outerLength; ++i)
                result += monomial.data[pos + atomHeaderSize + i];
              for (int i = outerLength; i < n; ++i)
                result -= monomial.data[pos + atomHeaderSize + i];
            }
        }
      else
        {
          for (int i = 0; i < n; ++i) result += monomial.data[pos + atomHeaderSize + i];
        }
      pos += atomHeaderSize + n;
    }
  return result;
}

BasisIndexKey basisIndexKey(const SymmetricMonomial& monomial, size_t pos)
{
  BasisIndexKey key;
  key.basisId = atomBasisIdAt(monomial, pos);
  int n = atomIndexLengthAt(monomial, pos);
  key.index.reserve(n);
  for (int i = 0; i < n; ++i) key.index.push_back(monomial.data[pos + atomHeaderSize + i]);
  return key;
}

class SymmetricEngineRing : public Ring
{
 private:
  const Ring *coefficientRing;
  mutable std::map<int, std::string> basisDisplays;
  mutable std::map<int, int> basisOrders;
  mutable std::map<int, bool> multiplicativeBases;
  mutable int powerSumBasisId = -1;
  mutable GCMap<int, ring_elem> hToPowerSumCache;
  mutable GCMap<int, ring_elem> eToPowerSumCache;
  mutable GCMap<int, ring_elem> qToPowerSumCache;
  mutable GCMap<int, ring_elem> bToPowerSumCache;
  mutable std::map<int, CharacterTable> characterTableCache;
  mutable GCMap<int, ring_elem> powerSumToCompleteCache;
  mutable GCMap<int, ring_elem> powerSumToElementaryCache;
  mutable GCMap<int, ring_elem> powerSumToQGeneratorCache;
  mutable GCMap<int, ring_elem> powerSumToBGeneratorCache;
  mutable GCMap<int, GCMap<std::string, ring_elem>> monomialToPowerSumCache;
  mutable std::map<std::string, std::vector<SchurConversionRecipeEntry>>
      powerSumsToSchurRecipeCache;
  mutable std::map<std::string, std::vector<LRProductTerm>> lrProductCache;
  mutable GCMap<std::string, ring_elem> schurCompletePlethysmCache;
  mutable GCMap<std::string, ring_elem> hJacobiTrudiCache;
  mutable GCMap<std::string, ring_elem> eJacobiTrudiCache;
  mutable ring_elem hallLittlewoodParameter;

  bool isMultiplicativeBasis(int basisId) const
  {
    auto it = multiplicativeBases.find(basisId);
    return it != multiplicativeBases.end() && it->second;
  }

  bool isPowerSumBasis(int basisId) const
  {
    return powerSumBasisId >= 0 && basisId == powerSumBasisId;
  }

  bool isForgottenDisplay(const std::string& display) const
  {
    return display == "f" || display == "ff";
  }

  void rememberBasis(int basisId,
                     const std::string& display,
                     int order,
                     bool isMultiplicative) const
  {
    if (!display.empty()) basisDisplays[basisId] = display;
    basisOrders[basisId] = order;
    multiplicativeBases[basisId] = isMultiplicative;
    if (display == "p") powerSumBasisId = basisId;
  }

  void rememberBasesFrom(const SymmetricEngineRing *R) const
  {
    for (const auto& item : R->basisDisplays) basisDisplays[item.first] = item.second;
    for (const auto& item : R->basisOrders) basisOrders[item.first] = item.second;
    for (const auto& item : R->multiplicativeBases)
      multiplicativeBases[item.first] = item.second;
    if (powerSumBasisId < 0) powerSumBasisId = R->powerSumBasisId;
  }

  std::string displayForBasis(int basisId) const
  {
    auto it = basisDisplays.find(basisId);
    if (it != basisDisplays.end()) return it->second;
    return "basis" + std::to_string(basisId);
  }

  int basisIdForDisplay(const std::string& display) const
  {
    for (const auto& item : basisDisplays)
      if (item.second == display) return item.first;
    return -1;
  }

  int basisOrderForId(int basisId) const
  {
    auto it = basisOrders.find(basisId);
    if (it != basisOrders.end()) return it->second;
    if (basisId == basisIdForDisplay("p")) return 10;
    if (basisId == basisIdForDisplay("h")) return 20;
    if (basisId == basisIdForDisplay("e")) return 30;
    if (basisId == basisIdForDisplay("S")) return 60;
    return 100;
  }

  const CharacterTable& characterTable(int degree) const
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

  int characterTableValue(const CharacterTable& table, size_t row, size_t col) const
  {
    int& value = table.values[row][col];
    if (value == unknownCharacterValue)
      value = characterValue(table.partitions[row], table.partitions[col]);
    return value;
  }

  ring_elem rationalCoefficient(long numerator, long denominator) const
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

  void clearHallLittlewoodCaches() const
  {
    qToPowerSumCache.clear();
    bToPowerSumCache.clear();
    powerSumToQGeneratorCache.clear();
    powerSumToBGeneratorCache.clear();
  }

  ring_elem hallLittlewoodFactor(const Partition& mu) const
  {
    ring_elem result = coefficientRing->one();
    for (int part : mu)
      {
        ring_elem tPower = coefficientRing->power(hallLittlewoodParameter, part);
        ring_elem oneMinus = coefficientRing->subtract(coefficientRing->one(), tPower);
        result = coefficientRing->mult(result, oneMinus);
      }
    return result;
  }

  ring_elem basisElementFromIndex(int basisId,
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

  ring_elem basisElementFromSkewIndex(int basisId,
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

  ring_elem basisPartElement(int basisId,
                             const std::string& display,
                             int order,
                             bool isMultiplicative,
                             int n) const
  {
    if (n < 0) return zero();
    if (n == 0) return one();
    return basisElementFromIndex(basisId, display, order, isMultiplicative, Partition{n});
  }

  std::string jacobiTrudiCacheKey(int basisId,
                                  const Partition& outer,
                                  const Partition& inner) const
  {
    return std::to_string(basisId) + "|" + partitionKey(outer) + "/" +
           partitionKey(inner);
  }

  int popcountMask(size_t mask) const
  {
    int result = 0;
    while (mask != 0)
      {
        result += static_cast<int>(mask & 1);
        mask >>= 1;
      }
    return result;
  }

  int selectedGreaterThan(size_t mask, size_t col, size_t n) const
  {
    int result = 0;
    for (size_t j = col + 1; j < n; ++j)
      if ((mask & (static_cast<size_t>(1) << j)) != 0) ++result;
    return result;
  }

  ring_elem jacobiTrudi(const Partition& outer,
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

  ring_elem scaled(ring_elem coeff, ring_elem f) const
  {
    return makePolyValue(multByCoefficient(coeff, polyValue(f)));
  }

  ring_elem coefficientQuotient(ring_elem numerator, ring_elem denominator) const
  {
    ring_elem quotient = coefficientRing->divide(numerator, denominator);
    ring_elem check = coefficientRing->mult(quotient, denominator);
    if (coefficientRing->is_equal(check, numerator)) return quotient;
    ERROR("coefficient division failed during basis conversion; use a coefficient ring where the required denominators are invertible, for example frac(QQ[t]) instead of QQ[t]");
    return coefficientRing->zero();
  }

  ring_elem hallLittlewoodCFactor(const Partition& lambda) const
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

  void addCoeff(CoeffMap& target, const Partition& index, ring_elem coeff) const
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

  CoeffMap scaledCoeffMap(ring_elem coeff, const CoeffMap& source) const
  {
    CoeffMap result;
    if (coefficientRing->is_zero(coeff)) return result;
    for (const auto& item : source)
      addCoeff(result, item.first, coefficientRing->mult(coeff, item.second));
    return result;
  }

  CoeffMap addCoeffMaps(const CoeffMap& a, const CoeffMap& b) const
  {
    CoeffMap result = a;
    for (const auto& item : b) addCoeff(result, item.first, item.second);
    return result;
  }

  CoeffMap multiplyCoeffMaps(const CoeffMap& a, const CoeffMap& b) const
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

  CoeffMap oneCoeffMap() const
  {
    return CoeffMap{{Partition{}, coefficientRing->one()}};
  }

  bool isPartitionIndex(const Partition& p) const
  {
    return trimTrailingZerosPartition(p) == normalizePartition(p);
  }

  int partitionPart(const Partition& p, size_t i) const
  {
    return i < p.size() ? p[i] : 0;
  }

  bool partitionContains(const Partition& outer, const Partition& inner) const
  {
    for (size_t i = 0; i < inner.size(); ++i)
      if (partitionPart(outer, i) < inner[i]) return false;
    return true;
  }

  std::string lrProductKey(const Partition& lambda, const Partition& mu) const
  {
    if (lexLessPartition(mu, lambda))
      return partitionKey(mu) + "*" + partitionKey(lambda);
    return partitionKey(lambda) + "*" + partitionKey(mu);
  }

  long lrCoefficient(const Partition& lambda,
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

  void partitionsContainingRec(const Partition& lambda,
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

  std::vector<Partition> partitionsContaining(const Partition& lambda,
                                              int addedWeight) const
  {
    std::vector<Partition> result;
    Partition current;
    int firstBound = lambda.empty() ? addedWeight : lambda.front() + addedWeight;
    partitionsContainingRec(lambda, addedWeight, 0, firstBound, current, result);
    return result;
  }

  const std::vector<LRProductTerm>& lrProduct(const Partition& a,
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

  ring_elem multiplySchurElements(ring_elem f,
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

  std::string powerSumsToSchurRecipeKey(int degree,
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
  powerSumsToSchurRecipe(int degree,
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

  Partition leadingPartition(const CoeffMap& H) const
  {
    if (H.empty()) return Partition{};
    Partition result = H.begin()->first;
    for (const auto& item : H)
      if (lexLessPartition(item.first, result)) result = item.first;
    return result;
  }

  ring_elem hPartToPowerSums(int n) const
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

  ring_elem ePartToPowerSums(int n) const
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

  ring_elem hallLittlewoodPartToPowerSums(int n, bool omega) const
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

  CoeffMap raisingExpansion(const Partition& lambda) const
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

  CoeffMap raisingGeneratorMap(const Partition& lambda) const
  {
    CoeffMap result;
    for (const auto& item : raisingExpansion(lambda))
      addCoeff(result, item.first, item.second);
    return result;
  }

  ring_elem hallCapitalToPowerSums(const Partition& lambda, bool omega) const
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

  ring_elem hallPToPowerSums(const Partition& lambda, bool omega) const
  {
    ring_elem numerator = hallCapitalToPowerSums(lambda, omega);
    return scaled(coefficientQuotient(coefficientRing->one(),
                                      hallLittlewoodCFactor(lambda)),
                  numerator);
  }

  ring_elem schurToPowerSums(const Partition& lambda) const
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

  ring_elem powerSumPartToComplete(int n, int hId, int hOrder) const
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

  ring_elem powerSumPartToElementary(int n, int eId, int eOrder) const
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

  ring_elem powerSumPartToHallGenerator(int n,
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

  CoeffMap powerSumPartToHallGeneratorMap(int n, bool omega) const
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

  CoeffMap powerSumIndexToHallGeneratorMap(const Partition& index, bool omega) const
  {
    CoeffMap result = oneCoeffMap();
    for (int part : index)
      result = multiplyCoeffMaps(result, powerSumPartToHallGeneratorMap(part, omega));
    return result;
  }

  CoeffMap powerSumsToHallGeneratorMap(ring_elem f, bool omega) const
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

  CoeffMap triangularReduceHallCapital(const CoeffMap& generatorMap,
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

  ring_elem powerSumToSchurLike(const Partition& mu,
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

  ring_elem powerSumToSchur(const Partition& mu, int schurId, int schurOrder) const
  {
    return powerSumToSchurLike(mu, schurId, schurOrder, "S", 1);
  }

  ring_elem powerSumsToSchurLike(ring_elem f,
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

  RingElemVector solveSquareSystem(RingElemMatrix M,
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

  ring_elem omegaPowerSums(ring_elem f) const
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

  ring_elem monomialBasisToPowerSums(const Partition& lambda, bool forgotten) const
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

  ring_elem powerSumIndexToMonomialTarget(const Partition& lambda,
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

  ring_elem coeffMapToElement(const CoeffMap& H,
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

  Partition replaceAdjacentPair(const Partition& alpha,
                                size_t pos,
                                int first,
                                int second) const
  {
    Partition result = alpha;
    result[pos] = first;
    result[pos + 1] = second;
    return result;
  }

  ring_elem straightenSchurAtom(const Partition& alpha,
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

  ring_elem omegaSchurAtomAsSchur(const Partition& alpha) const
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

  ring_elem straightenHallCapitalAtom(const Partition& alpha,
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

  ring_elem straightenAtom(const SymmetricMonomial& monomial, size_t pos) const
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

  ring_elem straightenMonomial(const SymmetricMonomial& monomial) const
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

  ring_elem straightenElement(ring_elem f) const
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

  int requiredBasisIdForDisplay(const std::string& display) const
  {
    int id = basisIdForDisplay(display);
    if (id < 0) ERROR("basis metadata for ", display.c_str(), " is not available");
    return id;
  }

  bool singleBasisIndexFromMonomial(const SymmetricMonomial& monomial,
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

  CoeffMap coefficientsInBasis(ring_elem f, int basisId) const
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

  bool coefficientsInBasisIfPossible(ring_elem f,
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

  ring_elem coefficientPairing(const CoeffMap& fCoeffs,
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

  ring_elem powerSumInnerProductFactor(const Partition& lambda) const
  {
    ring_elem numerator = rationalCoefficient(zValue(lambda), 1);
    if (error()) return coefficientRing->zero();
    return coefficientQuotient(numerator, hallLittlewoodFactor(lambda));
  }

  ring_elem powerSumPairing(const CoeffMap& fCoeffs,
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

  std::map<int, InnerProductTarget> innerProductTargetMap(M2_arrayint innerProductMap) const
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

  bool directHallInnerProductFromMetadata(
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

  bool directPowerSumSchurInnerProduct(ring_elem f,
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

  ring_elem basisElementForDisplay(const std::string& display,
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

  ring_elem hallInnerProductElements(ring_elem f,
                                     ring_elem g,
                                     const std::map<int, InnerProductTarget>& metadata = {}) const
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

  ring_elem skewQOrBFunction(const Partition& lambda,
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

  ring_elem replaceSingleBasis(ring_elem f,
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

  ring_elem skewPOrRToPowerSums(const Partition& lambda,
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

  ring_elem skewHallLittlewoodToPowerSums(const Partition& lambda,
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

  ring_elem powerSumsToHallCapitalTarget(ring_elem f,
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

  CoeffMap schurGeneratorMap(const Partition& lambda,
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

  CoeffMap triangularReduceSchur(const CoeffMap& generatorMap,
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

  bool directTriangularSchurConversion(ring_elem f,
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

  bool directTriangularHallConversion(ring_elem f,
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

  Partition atomIndex(const SymmetricMonomial& monomial, size_t pos) const
  {
    Partition result;
    int n = atomIndexLengthAt(monomial, pos);
    result.reserve(n);
    for (int i = 0; i < n; ++i) result.push_back(monomial.data[pos + atomHeaderSize + i]);
    return result;
  }

  Partition atomOuterIndex(const SymmetricMonomial& monomial, size_t pos) const
  {
    Partition result;
    int n = atomOuterLengthAt(monomial, pos);
    result.reserve(n);
    for (int i = 0; i < n; ++i) result.push_back(monomial.data[pos + atomHeaderSize + i]);
    return result;
  }

  Partition atomInnerIndex(const SymmetricMonomial& monomial, size_t pos) const
  {
    Partition result;
    int outerLength = atomOuterLengthAt(monomial, pos);
    int innerLength = atomInnerLengthAt(monomial, pos);
    result.reserve(innerLength);
    for (int i = 0; i < innerLength; ++i)
      result.push_back(monomial.data[pos + atomHeaderSize + outerLength + i]);
    return result;
  }

  ring_elem atomToPowerSums(const SymmetricMonomial& monomial, size_t pos) const
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

  ring_elem monomialToPowerSums(const SymmetricMonomial& monomial) const
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

  ring_elem elementToPowerSums(ring_elem f) const
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

  ring_elem powerSumElementFromIndex(const Partition& index) const
  {
    Partition normalized = normalizePartition(index);
    if (normalized.empty()) return one();
    return basisElementFromIndex(powerSumBasisId, "p", 10, true, normalized);
  }

  ring_elem adamsPowerSums(ring_elem f, int multiplier) const
  {
    const auto *poly = polyValue(f);
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during plethysm");
            return zero();
          }
        Partition transformed;
        transformed.reserve(index.size());
        for (int part : index)
          if (part * multiplier > 0) transformed.push_back(part * multiplier);
        ring_elem termElement = powerSumElementFromIndex(transformed);
        result = add(result, scaled(term.coeff, termElement));
      }
    return result;
  }

  ring_elem plethysmPowerSums(ring_elem fPowerSums, ring_elem gPowerSums) const
  {
    const auto *poly = polyValue(fPowerSums);
    GCMap<int, ring_elem> adamsCache;
    ring_elem result = zero();
    for (const auto& term : poly->terms)
      {
        Partition index;
        if (!powerSumIndexFromMonomial(term.monomial, index))
          {
            ERROR("expected a pure power-sum expression during plethysm");
            return zero();
          }
        ring_elem substituted = one();
        for (int part : index)
          {
            auto cached = adamsCache.find(part);
            if (cached == adamsCache.end())
              {
                ring_elem factor = adamsPowerSums(gPowerSums, part);
                if (error()) return zero();
                cached = adamsCache.emplace(part, factor).first;
              }
            ring_elem factor = copyPolyValue(polyValue(cached->second));
            substituted = mult(substituted, factor);
          }
        result = add(result, scaled(term.coeff, substituted));
      }
    return result;
  }

  int singleBasisIdInMonomial(const SymmetricMonomial& monomial) const
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

  int singleBasisId(ring_elem f) const
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

  bool powerSumIndexFromMonomial(const SymmetricMonomial& monomial,
                                 Partition& index) const
  {
    index.clear();
    if (monomial.data.empty()) return true;
    if (!isSinglePowerSumBlock(monomial)) return false;
    index = atomIndex(monomial, 0);
    return true;
  }

  ring_elem powerSumMonomialToTarget(const Partition& index,
                                     const std::string& targetDisplay,
                                     int targetBasisId,
                                     int targetDisplayOrder,
                                     bool targetIsMultiplicative) const
  {
    if (targetDisplay == "p")
      return basisElementFromIndex(targetBasisId,
                                   targetDisplay,
                                   targetDisplayOrder,
                                   targetIsMultiplicative,
                                   index);

    if (targetDisplay == "h")
      {
        ring_elem result = one();
        for (int part : index)
          result = mult(result, powerSumPartToComplete(part, targetBasisId, targetDisplayOrder));
        return result;
      }

    if (targetDisplay == "e")
      {
        ring_elem result = one();
        for (int part : index)
          result = mult(result, powerSumPartToElementary(part, targetBasisId, targetDisplayOrder));
        return result;
      }

    if (targetDisplay == "q" || targetDisplay == "b")
      {
        ring_elem result = one();
        bool omega = targetDisplay == "b";
        for (int part : index)
          result = mult(result,
                        powerSumPartToHallGenerator(part,
                                                    targetBasisId,
                                                    targetDisplay,
                                                    targetDisplayOrder,
                                                    omega));
        return result;
      }

    if (targetDisplay == "m" || isForgottenDisplay(targetDisplay))
      return powerSumIndexToMonomialTarget(index,
                                           targetBasisId,
                                           targetDisplay,
                                           targetDisplayOrder,
                                           isForgottenDisplay(targetDisplay));

    if (targetDisplay == "S")
      return powerSumToSchur(index, targetBasisId, targetDisplayOrder);

    if (targetDisplay == "Somega")
      {
        long sign = ((partitionWeight(index) - partitionLength(index)) % 2 == 0) ? 1 : -1;
        return powerSumToSchurLike(index,
                                   targetBasisId,
                                   targetDisplayOrder,
                                   "Somega",
                                   sign);
      }

    ERROR("basis conversion from power sums is not implemented for basis ", targetDisplay.c_str());
    return zero();
  }

  ring_elem powerSumsToTarget(ring_elem f,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetDisplayOrder,
                              bool targetIsMultiplicative) const
  {
    if (targetDisplay == "Q" || targetDisplay == "B" ||
        targetDisplay == "P" || targetDisplay == "R")
      return powerSumsToHallCapitalTarget(f,
                                          targetBasisId,
                                          targetDisplay,
                                          targetDisplayOrder);
    if (targetDisplay == "S")
      return powerSumsToSchurLike(f,
                                  targetBasisId,
                                  targetDisplayOrder,
                                  targetDisplay,
                                  false);
    if (targetDisplay == "Somega")
      return powerSumsToSchurLike(f,
                                  targetBasisId,
                                  targetDisplayOrder,
                                  targetDisplay,
                                  true);

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
        ring_elem converted = powerSumMonomialToTarget(index,
                                                       targetDisplay,
                                                       targetBasisId,
                                                       targetDisplayOrder,
                                                       targetIsMultiplicative);
        if (error()) return zero();
        result = add(result, scaled(term.coeff, converted));
      }
    return result;
  }

  bool atomToDirectTarget(const SymmetricMonomial& monomial,
                          size_t pos,
                          int targetBasisId,
                          const std::string& targetDisplay,
                          int targetDisplayOrder,
                          bool targetIsMultiplicative,
                          ring_elem& result) const
  {
    int basisId = atomBasisIdAt(monomial, pos);
    std::string display = displayForBasis(basisId);
    Partition index = atomIndex(monomial, pos);

    if (!atomIsSkewAt(monomial, pos) && display == targetDisplay)
      {
        result = basisElementFromIndex(targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       targetIsMultiplicative,
                                       index);
        return true;
      }

    if (targetDisplay == "h" && display == "S")
      {
        Partition outer = atomIsSkewAt(monomial, pos) ? atomOuterIndex(monomial, pos)
                                                      : index;
        Partition inner = atomIsSkewAt(monomial, pos) ? atomInnerIndex(monomial, pos)
                                                      : Partition{};
        result = jacobiTrudi(outer,
                             inner,
                             targetBasisId,
                             targetDisplay,
                             targetDisplayOrder,
                             targetIsMultiplicative);
        return true;
      }

    if (targetDisplay == "e" && display == "Somega")
      {
        Partition outer = atomIsSkewAt(monomial, pos) ? atomOuterIndex(monomial, pos)
                                                      : index;
        Partition inner = atomIsSkewAt(monomial, pos) ? atomInnerIndex(monomial, pos)
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
    if (display != "p" && display != "h" && display != "e" &&
        display != "q" && display != "b" && display != "m" &&
        !isForgottenDisplay(display) && display != "S" && display != "Somega" &&
        display != "Q" && display != "B" && display != "P" && display != "R")
      return false;
    ring_elem inPowerSums = atomToPowerSums(monomial, pos);
    if (error()) return false;
    result = powerSumsToTarget(inPowerSums,
                               targetBasisId,
                               targetDisplay,
                               targetDisplayOrder,
                               targetIsMultiplicative);
    return !error();
  }

  bool monomialToDirectTarget(const SymmetricMonomial& monomial,
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
        if (!atomToDirectTarget(monomial,
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

  bool elementToDirectTarget(ring_elem f,
                             int targetBasisId,
                             const std::string& targetDisplay,
                             int targetDisplayOrder,
                             bool targetIsMultiplicative,
                             ring_elem& result) const
  {
    if (directTriangularSchurConversion(f,
                                        targetBasisId,
                                        targetDisplay,
                                        targetDisplayOrder,
                                        result))
      return true;
    if (error()) return false;
    if (directTriangularHallConversion(f,
                                       targetBasisId,
                                       targetDisplay,
                                       targetDisplayOrder,
                                       result))
      return true;
    if (error()) return false;

    if (targetDisplay != "h" && targetDisplay != "e" && !targetIsMultiplicative)
      return false;
    result = zero();
    const auto *poly = polyValue(f);
    for (const auto& term : poly->terms)
      {
        ring_elem converted;
        if (!monomialToDirectTarget(term.monomial,
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

  ring_elem omegaDirectAtom(const SymmetricMonomial& monomial,
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

  ring_elem omegaMonomial(const SymmetricMonomial& monomial,
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

  std::map<int, OmegaTarget> omegaTargetMap(M2_arrayint omegaMap) const
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

  SymmetricMonomial canonicalMonomial(const SymmetricMonomial& monomial) const
  {
    std::vector<std::vector<int>> blocks;
    std::map<std::pair<int, int>, std::vector<int>> compressed;

    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        int order = atomOrderAt(monomial, pos);
        int basisId = atomBasisIdAt(monomial, pos);
        int n = atomIndexLengthAt(monomial, pos);
        if (!atomIsSkewAt(monomial, pos) && isMultiplicativeBasis(basisId))
          {
            auto& index = compressed[std::make_pair(order, basisId)];
            for (int i = 0; i < n; ++i)
              index.push_back(monomial.data[pos + atomHeaderSize + i]);
          }
        else
          {
            blocks.push_back(atomBlockAt(monomial, pos));
          }
        pos += atomHeaderSize + n;
      }

    for (auto& item : compressed)
      {
        auto& index = item.second;
        std::sort(index.begin(), index.end(), std::greater<int>());
        blocks.push_back(makeAtomBlock(item.first.first, item.first.second, 0, index));
      }

    std::sort(blocks.begin(), blocks.end(), blockLess);
    SymmetricMonomial result;
    size_t totalLength = 0;
    for (const auto& block : blocks) totalLength += block.size();
    result.data.reserve(totalLength);
    for (const auto& block : blocks) appendAtomBlock(result, block);
    return result;
  }

  bool isSinglePowerSumBlock(const SymmetricMonomial& monomial) const
  {
    if (monomial.data.empty()) return false;
    if (atomIsSkewAt(monomial, 0)) return false;
    if (atomBasisIdAt(monomial, 0) != powerSumBasisId) return false;
    return atomLengthAt(monomial, 0) == monomial.data.size();
  }

  SymmetricMonomial multiplyPowerSumMonomials(const SymmetricMonomial& a,
                                             const SymmetricMonomial& b) const
  {
    std::vector<int> index;
    int order = atomOrderAt(a, 0);
    int nA = atomIndexLengthAt(a, 0);
    int nB = atomIndexLengthAt(b, 0);
    index.reserve(nA + nB);
    for (int i = 0; i < nA; ++i) index.push_back(a.data[atomHeaderSize + i]);
    for (int i = 0; i < nB; ++i) index.push_back(b.data[atomHeaderSize + i]);
    std::sort(index.begin(), index.end(), std::greater<int>());

    SymmetricMonomial result;
    appendAtomBlock(result, makeAtomBlock(order, powerSumBasisId, 0, index));
    return result;
  }

  SymmetricMonomial multiplyMonomials(const SymmetricMonomial& a,
                                      const SymmetricMonomial& b) const
  {
    if (a.data.empty()) return b;
    if (b.data.empty()) return a;
    if (isPowerSumBasis(powerSumBasisId) &&
        isSinglePowerSumBlock(a) &&
        isSinglePowerSumBlock(b))
      return multiplyPowerSumMonomials(a, b);

    SymmetricMonomial result;
    result.data.reserve(a.data.size() + b.data.size());
    result.data.insert(result.data.end(), a.data.begin(), a.data.end());
    result.data.insert(result.data.end(), b.data.begin(), b.data.end());
    return canonicalMonomial(result);
  }

  ring_elem copyPolyValue(const SymmetricRingPoly *poly) const
  {
    auto result = new SymmetricRingPoly;
    result->terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      result->terms.push_back({coefficientRing->copy(term.coeff), term.monomial});
    return makePolyValue(result);
  }

  void appendTermIfNonZero(VECTOR(SymmetricTerm)& terms,
                           ring_elem coeff,
                           const SymmetricMonomial& monomial) const
  {
    if (!coefficientRing->is_zero(coeff)) terms.push_back({coeff, monomial});
  }

  ring_elem fromTermVector(VECTOR(SymmetricTerm)& terms, bool isSorted) const
  {
    if (!isSorted)
      std::sort(terms.begin(), terms.end(),
                [](const SymmetricTerm& a, const SymmetricTerm& b) {
                  return compareMonomials(a.monomial, b.monomial) == LT;
                });

    auto result = new SymmetricRingPoly;
    result->terms.reserve(terms.size());
    size_t i = 0;
    while (i < terms.size())
      {
        ring_elem coeff = terms[i].coeff;
        size_t j = i + 1;
        while (j < terms.size() &&
               compareMonomials(terms[i].monomial, terms[j].monomial) == EQ)
          {
            coeff = coefficientRing->add(coeff, terms[j].coeff);
            ++j;
          }
        appendTermIfNonZero(result->terms, coeff, terms[i].monomial);
        i = j;
      }
    return makePolyValue(result);
  }

  ring_elem concatenateTerms(const SymmetricRingPoly *first,
                             const SymmetricRingPoly *second) const
  {
    auto result = new SymmetricRingPoly;
    result->terms.reserve(first->terms.size() + second->terms.size());
    result->terms.insert(result->terms.end(), first->terms.begin(), first->terms.end());
    result->terms.insert(result->terms.end(), second->terms.begin(), second->terms.end());
    return makePolyValue(result);
  }

  void addToAccumulator(GCMap<std::vector<int>, ring_elem>& accumulator,
                        const SymmetricMonomial& monomial,
                        ring_elem coeff) const
  {
    if (coefficientRing->is_zero(coeff)) return;
    auto key = monomialKey(monomial);
    auto existing = accumulator.find(key);
    if (existing == accumulator.end())
      {
        accumulator[key] = coeff;
        return;
      }
    ring_elem sum = coefficientRing->add(existing->second, coeff);
    if (coefficientRing->is_zero(sum))
      accumulator.erase(existing);
    else
      existing->second = sum;
  }

  ring_elem fromAccumulator(
      const GCMap<std::vector<int>, ring_elem>& accumulator) const
  {
    auto result = new SymmetricRingPoly;
    result->terms.reserve(accumulator.size());
    for (const auto& term : accumulator)
      if (!coefficientRing->is_zero(term.second))
        result->terms.push_back({term.second, monomialFromKey(term.first)});
    return makePolyValue(result);
  }

  bool promoteInputElement(const RingElement *input, ring_elem &result) const
  {
    const Ring *inputRing = input->get_ring();
    if (inputRing == this)
      {
        result = input->get_value();
        return true;
      }
    return promote(inputRing, input->get_value(), result);
  }

 public:
  explicit SymmetricEngineRing(const Ring *A)
      : coefficientRing(A), hallLittlewoodParameter(A->from_long(0))
  {}

  static SymmetricEngineRing *create(const Ring *A)
  {
    auto result = new SymmetricEngineRing(A);
    result->initialize_ring(A->characteristic());
    result->zeroV = result->from_long(0);
    result->oneV = result->from_long(1);
    result->minus_oneV = result->from_long(-1);
    return result;
  }

  const Ring *getCoefficientRing() const { return coefficientRing; }

  void rememberBasisMetadata(int basisId,
                             const std::string& display,
                             int order,
                             bool isMultiplicative) const
  {
    rememberBasis(basisId, display, order, isMultiplicative);
  }

  bool setHallLittlewoodParameter(const RingElement *t) const
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

  ring_elem fromCoeff(ring_elem coeff) const
  {
    auto result = new SymmetricRingPoly;
    if (!coefficientRing->is_zero(coeff)) result->terms.push_back({coeff, {}});
    return makePolyValue(result);
  }

  ring_elem basisElement(int basisId,
                         const std::string& display,
                         int order,
                         bool isMultiplicative,
                         int innerLength,
                         M2_arrayint index) const
  {
    int payloadLength = index == nullptr ? 0 : index->len;
    if (innerLength < 0 || innerLength > payloadLength)
      {
        ERROR("invalid skew inner shape length");
        return zero();
      }
    rememberBasis(basisId, display, order, isMultiplicative);
    auto result = new SymmetricRingPoly;
    SymmetricMonomial monomial;
    appendAtomBlock(monomial, makeAtomBlock(order, basisId, innerLength, index));
    result->terms.push_back({coefficientRing->one(), canonicalMonomial(monomial)});
    return makePolyValue(result);
  }

  bool getScalar(const SymmetricRingPoly *f, ring_elem &result) const
  {
    if (f->terms.empty())
      {
        result = coefficientRing->zero();
        return true;
      }
    if (f->terms.size() != 1) return false;
    if (!f->terms[0].monomial.data.empty()) return false;
    result = f->terms[0].coeff;
    return true;
  }

  bool hasPowerSumConversionHook(const SymmetricMonomial& monomial, size_t pos) const
  {
    auto key = basisIndexKey(monomial, pos);
    return isPowerSumBasis(key.basisId);
  }

  std::string displayIndex(const SymmetricMonomial& monomial, size_t pos) const
  {
    int n = atomIndexLengthAt(monomial, pos);
    if (n <= 0) return "{}";
    if (atomIsSkewAt(monomial, pos))
      {
        int outerLength = atomOuterLengthAt(monomial, pos);
        int innerLength = atomInnerLengthAt(monomial, pos);
        if (outerLength + innerLength != n) return "{}/{}";
        std::vector<std::string> lambdaParts;
        std::vector<std::string> muParts;
        lambdaParts.reserve(outerLength);
        muParts.reserve(innerLength);
        for (int i = 0; i < outerLength; ++i)
          lambdaParts.push_back(std::to_string(monomial.data[pos + atomHeaderSize + i]));
        for (int i = outerLength; i < n; ++i)
          muParts.push_back(std::to_string(monomial.data[pos + atomHeaderSize + i]));
        return "{{" + join(lambdaParts, ",") + "}/{" + join(muParts, ",") + "}}";
      }
    if (n == 1) return std::to_string(monomial.data[pos + atomHeaderSize]);
    std::vector<std::string> parts;
    parts.reserve(n);
    for (int i = 0; i < n; ++i)
      parts.push_back(std::to_string(monomial.data[pos + atomHeaderSize + i]));
    return "{" + join(parts, ",") + "}";
  }

  std::string displayAtom(const SymmetricMonomial& monomial, size_t pos) const
  {
    return displayForBasis(atomBasisIdAt(monomial, pos)) + "_" +
           displayIndex(monomial, pos);
  }

  std::string displayMonomial(const SymmetricMonomial& monomial) const
  {
    if (monomial.data.empty()) return "1";
    std::vector<std::string> factors;
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        factors.push_back(displayAtom(monomial, pos));
        pos += atomLengthAt(monomial, pos);
      }
    return join(factors, "*");
  }

  std::string elementString(ring_elem f) const
  {
    const auto *poly = polyValue(f);
    if (poly->terms.empty()) return "0";
    std::vector<std::string> pieces;
    pieces.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        if (term.monomial.data.empty())
          {
            pieces.push_back(coeffToString(coefficientRing, term.coeff));
          }
        else if (coefficientRing->is_equal(term.coeff, coefficientRing->one()))
          {
            pieces.push_back(displayMonomial(term.monomial));
          }
        else if (coefficientRing->is_equal(term.coeff, coefficientRing->minus_one()))
          {
            pieces.push_back("-" + displayMonomial(term.monomial));
          }
        else
          {
            pieces.push_back("(" + coeffToString(coefficientRing, term.coeff) +
                             ")*" + displayMonomial(term.monomial));
          }
      }
    return join(pieces, " + ");
  }

  int elementWeight(ring_elem f) const
  {
    const auto *poly = polyValue(f);
    bool found = false;
    int result = 0;
    for (const auto& term : poly->terms)
      {
        int wt = monomialWeight(term.monomial);
        if (!found)
          {
            result = wt;
            found = true;
          }
        else if (result != wt)
          return -1;
      }
    return found ? result : 0;
  }

  virtual unsigned int computeHashValue(const ring_elem a) const
  {
    size_t seed = 7919237;
    const auto *poly = polyValue(a);
    for (const auto& term : poly->terms)
      {
        hashCombine(seed, coefficientRing->computeHashValue(term.coeff));
        for (int value : term.monomial.data) hashCombine(seed, value);
      }
    return static_cast<unsigned int>(seed);
  }

  virtual void text_out(buffer &o) const
  {
    o << "SymmetricRing(";
    coefficientRing->text_out(o);
    o << ")";
  }

  virtual void elem_text_out(buffer &o,
                             const ring_elem f,
                             bool p_one = true,
                             bool p_plus = false,
                             bool p_parens = false) const
  {
    (void)p_one;
    if (p_plus && !is_zero(f)) o << '+';
    bool parens = p_parens && polyValue(f)->terms.size() > 1;
    if (parens) o << '(';
    o << elementString(f);
    if (parens) o << ')';
  }

  virtual ring_elem from_long(long n) const
  {
    return fromCoeff(coefficientRing->from_long(n));
  }

  virtual ring_elem from_int(mpz_srcptr n) const
  {
    return fromCoeff(coefficientRing->from_int(n));
  }

  virtual bool from_rational(mpq_srcptr q, ring_elem &result) const
  {
    ring_elem coeff;
    if (!coefficientRing->from_rational(q, coeff)) return false;
    result = fromCoeff(coeff);
    return true;
  }

  virtual bool promote(const Ring *Rf, const ring_elem f, ring_elem &result) const
  {
    if (Rf == this)
      {
        result = f;
        return true;
      }
    if (Rf == coefficientRing)
      {
        result = fromCoeff(f);
        return true;
      }
    if (Rf == globalZZ)
      {
        ring_elem coeff;
        if (!coefficientRing->promote(globalZZ, f, coeff)) return false;
        result = fromCoeff(coeff);
        return true;
      }
    const auto *Sf = dynamic_cast<const SymmetricEngineRing *>(Rf);
    if (Sf != nullptr)
      {
        rememberBasesFrom(Sf);
        GCMap<std::vector<int>, ring_elem> accumulator;
        const auto *poly = polyValue(f);
        for (const auto& term : poly->terms)
          {
            ring_elem coeff;
            if (!coefficientRing->promote(
                    Sf->getCoefficientRing(), term.coeff, coeff))
              return false;
            addToAccumulator(accumulator, term.monomial, coeff);
          }
        result = fromAccumulator(accumulator);
        return true;
      }
    return false;
  }

  virtual bool lift(const Ring *Rg, const ring_elem f, ring_elem &result) const
  {
    const auto *poly = polyValue(f);
    if (Rg == coefficientRing) return getScalar(poly, result);
    if (Rg == globalZZ)
      {
        ring_elem coeff;
        if (!getScalar(poly, coeff)) return false;
        return coefficientRing->lift(globalZZ, coeff, result);
      }
    const auto *Sg = dynamic_cast<const SymmetricEngineRing *>(Rg);
    if (Sg != nullptr)
      {
        Sg->rememberBasesFrom(this);
        GCMap<std::vector<int>, ring_elem> accumulator;
        for (const auto& term : poly->terms)
          {
            ring_elem coeff;
            if (!coefficientRing->lift(
                    Sg->getCoefficientRing(), term.coeff, coeff))
              return false;
            Sg->addToAccumulator(accumulator, term.monomial, coeff);
          }
        result = Sg->fromAccumulator(accumulator);
        return true;
      }
    return false;
  }

  virtual bool is_unit(const ring_elem f) const
  {
    ring_elem coeff;
    return getScalar(polyValue(f), coeff) && coefficientRing->is_unit(coeff);
  }

  virtual bool is_zero(const ring_elem f) const
  {
    return polyValue(f)->terms.empty();
  }

  virtual bool is_equal(const ring_elem f, const ring_elem g) const
  {
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    if (left->terms.size() != right->terms.size()) return false;
    for (size_t i = 0; i < left->terms.size(); ++i)
      {
        if (compareMonomials(left->terms[i].monomial, right->terms[i].monomial) != EQ)
          return false;
        if (!coefficientRing->is_equal(left->terms[i].coeff, right->terms[i].coeff))
          return false;
      }
    return true;
  }

  virtual int compare_elems(const ring_elem f, const ring_elem g) const
  {
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    size_t n = std::min(left->terms.size(), right->terms.size());
    for (size_t i = 0; i < n; ++i)
      {
        int monomialCmp = compareMonomials(left->terms[i].monomial, right->terms[i].monomial);
        if (monomialCmp != EQ) return monomialCmp;
        int coeffCmp =
            coefficientRing->compare_elems(left->terms[i].coeff, right->terms[i].coeff);
        if (coeffCmp != EQ) return coeffCmp;
      }
    if (left->terms.size() < right->terms.size()) return LT;
    if (left->terms.size() > right->terms.size()) return GT;
    return EQ;
  }

  virtual ring_elem copy(const ring_elem f) const
  {
    return copyPolyValue(polyValue(f));
  }

  virtual void remove(ring_elem &f) const { (void)f; }

  virtual ring_elem negate(const ring_elem f) const
  {
    auto result = new SymmetricRingPoly;
    const auto *poly = polyValue(f);
    result->terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      result->terms.push_back({coefficientRing->negate(term.coeff), term.monomial});
    return makePolyValue(result);
  }

  virtual ring_elem add(const ring_elem f, const ring_elem g) const
  {
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    if (left->terms.empty()) return copyPolyValue(right);
    if (right->terms.empty()) return copyPolyValue(left);

    if (compareMonomials(left->terms.back().monomial,
                         right->terms.front().monomial) == LT)
      return concatenateTerms(left, right);
    if (compareMonomials(right->terms.back().monomial,
                         left->terms.front().monomial) == LT)
      return concatenateTerms(right, left);

    auto result = new SymmetricRingPoly;
    result->terms.reserve(left->terms.size() + right->terms.size());
    size_t i = 0;
    size_t j = 0;
    while (i < left->terms.size() && j < right->terms.size())
      {
        int cmp = compareMonomials(left->terms[i].monomial,
                                   right->terms[j].monomial);
        if (cmp == LT)
          result->terms.push_back(left->terms[i++]);
        else if (cmp == GT)
          result->terms.push_back(right->terms[j++]);
        else
          {
            ring_elem coeff = coefficientRing->add(left->terms[i].coeff,
                                                   right->terms[j].coeff);
            appendTermIfNonZero(result->terms, coeff, left->terms[i].monomial);
            ++i;
            ++j;
          }
      }
    result->terms.insert(result->terms.end(), left->terms.begin() + i, left->terms.end());
    result->terms.insert(result->terms.end(), right->terms.begin() + j, right->terms.end());
    return makePolyValue(result);
  }

  virtual ring_elem subtract(const ring_elem f, const ring_elem g) const
  {
    return add(f, negate(g));
  }

  SymmetricRingPoly *multByCoefficient(ring_elem coeff,
                                       const SymmetricRingPoly *poly) const
  {
    auto result = new SymmetricRingPoly;
    if (coefficientRing->is_zero(coeff)) return result;
    result->terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        ring_elem c = coefficientRing->mult(coeff, term.coeff);
        if (!coefficientRing->is_zero(c)) result->terms.push_back({c, term.monomial});
      }
    return result;
  }

  virtual ring_elem mult(const ring_elem f, const ring_elem g) const
  {
    ring_elem scalar;
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    if (getScalar(left, scalar)) return makePolyValue(multByCoefficient(scalar, right));
    if (getScalar(right, scalar)) return makePolyValue(multByCoefficient(scalar, left));
    if (left->terms.size() == 1 && right->terms.size() == 1)
      {
        ring_elem coeff = coefficientRing->mult(left->terms[0].coeff,
                                                right->terms[0].coeff);
        auto result = new SymmetricRingPoly;
        if (!coefficientRing->is_zero(coeff))
          result->terms.push_back({coeff,
                                   multiplyMonomials(left->terms[0].monomial,
                                                     right->terms[0].monomial)});
        return makePolyValue(result);
      }
    VECTOR(SymmetricTerm) products;
    products.reserve(left->terms.size() * right->terms.size());
    for (const auto& lt : left->terms)
      for (const auto& rt : right->terms)
        {
          ring_elem coeff = coefficientRing->mult(lt.coeff, rt.coeff);
          if (!coefficientRing->is_zero(coeff))
            products.push_back(
                {coeff, multiplyMonomials(lt.monomial, rt.monomial)});
        }
    return fromTermVector(products, false);
  }

  ring_elem batchSum(engine_RawRingElementArray elements) const
  {
    if (elements == nullptr || elements->len == 0) return zero();

    size_t termCount = 0;
    for (int i = 0; i < elements->len; ++i)
      {
        ring_elem promoted;
        if (!promoteInputElement(elements->array[i], promoted)) return zero();
        termCount += polyValue(promoted)->terms.size();
      }

    VECTOR(SymmetricTerm) terms;
    terms.reserve(termCount);
    bool sorted = true;
    bool havePrevious = false;
    SymmetricMonomial previous;

    for (int i = 0; i < elements->len; ++i)
      {
        ring_elem promoted;
        if (!promoteInputElement(elements->array[i], promoted)) return zero();
        const auto *poly = polyValue(promoted);
        for (const auto& term : poly->terms)
          {
            if (havePrevious && compareMonomials(previous, term.monomial) == GT)
              sorted = false;
            previous = term.monomial;
            havePrevious = true;
            terms.push_back(term);
          }
      }
    return fromTermVector(terms, sorted);
  }

  ring_elem batchProduct(engine_RawRingElementArray elements) const
  {
    if (elements == nullptr || elements->len == 0) return one();

    ring_elem coeff = coefficientRing->one();
    SymmetricMonomial monomial;
    bool singleTermProduct = true;

    for (int i = 0; i < elements->len; ++i)
      {
        ring_elem promoted;
        if (!promoteInputElement(elements->array[i], promoted)) return zero();
        const auto *poly = polyValue(promoted);
        if (poly->terms.empty()) return zero();
        if (poly->terms.size() != 1)
          {
            singleTermProduct = false;
            break;
          }
        coeff = coefficientRing->mult(coeff, poly->terms[0].coeff);
        if (coefficientRing->is_zero(coeff)) return zero();
        monomial = multiplyMonomials(monomial, poly->terms[0].monomial);
      }

    if (singleTermProduct)
      {
        auto result = new SymmetricRingPoly;
        result->terms.push_back({coeff, monomial});
        return makePolyValue(result);
      }

    ring_elem result = one();
    for (int i = 0; i < elements->len; ++i)
      {
        ring_elem promoted;
        if (!promoteInputElement(elements->array[i], promoted)) return zero();
        result = mult(result, promoted);
      }
    return result;
  }

  bool singleSchurPartition(ring_elem f, int schurId, Partition& lambda) const
  {
    const auto *poly = polyValue(f);
    if (poly->terms.size() != 1) return false;
    if (!coefficientRing->is_equal(poly->terms[0].coeff, coefficientRing->one()))
      return false;
    if (!singleBasisIndexFromMonomial(poly->terms[0].monomial, schurId, lambda))
      return false;
    return isPartitionIndex(lambda);
  }

  std::string schurCompletePlethysmKey(int n,
                                       int schurId,
                                       const Partition& inner) const
  {
    return std::to_string(n) + "|" + std::to_string(schurId) + "|" +
           partitionKey(inner);
  }

  ring_elem schurCompletePlethysm(int n,
                                  const Partition& inner,
                                  int powerSumBasisId0,
                                  const std::string& powerSumDisplay,
                                  int powerSumOrder,
                                  bool powerSumIsMultiplicative,
                                  int schurId,
                                  const std::string& schurDisplay,
                                  int schurOrder) const
  {
    if (n < 0) return zero();
    if (n == 0) return one();

    std::string key = schurCompletePlethysmKey(n, schurId, inner);
    auto cached = schurCompletePlethysmCache.find(key);
    if (cached != schurCompletePlethysmCache.end())
      return copyPolyValue(polyValue(cached->second));

    ring_elem innerSchur =
        basisElementFromIndex(schurId, schurDisplay, schurOrder, false, inner);
    ring_elem total = zero();
    for (int i = 1; i <= n; ++i)
      {
        ring_elem pI = basisElementFromIndex(powerSumBasisId0,
                                             powerSumDisplay,
                                             powerSumOrder,
                                             powerSumIsMultiplicative,
                                             Partition{i});
        ring_elem pIAtInner = plethysm(pI,
                                       innerSchur,
                                       powerSumBasisId0,
                                       powerSumDisplay,
                                       powerSumOrder,
                                       powerSumIsMultiplicative);
        if (error()) return zero();
        ring_elem pIAtInnerSchur = toBasis(pIAtInner,
                                           powerSumBasisId0,
                                           powerSumDisplay,
                                           powerSumOrder,
                                           powerSumIsMultiplicative,
                                           schurId,
                                           schurDisplay,
                                           schurOrder,
                                           false);
        if (error()) return zero();
        ring_elem rest = schurCompletePlethysm(n - i,
                                               inner,
                                               powerSumBasisId0,
                                               powerSumDisplay,
                                               powerSumOrder,
                                               powerSumIsMultiplicative,
                                               schurId,
                                               schurDisplay,
                                               schurOrder);
        if (error()) return zero();
        ring_elem product =
            multiplySchurElements(pIAtInnerSchur, rest, schurId, schurDisplay, schurOrder);
        if (error()) return zero();
        total = add(total, product);
      }

    ring_elem reciprocal = rationalCoefficient(1, n);
    if (error()) return zero();
    ring_elem result = scaled(reciprocal, total);
    schurCompletePlethysmCache[key] = result;
    return copyPolyValue(polyValue(result));
  }

  ring_elem schurPlethysmJacobiTrudi(const Partition& outer,
                                     const Partition& inner,
                                     int powerSumBasisId0,
                                     const std::string& powerSumDisplay,
                                     int powerSumOrder,
                                     bool powerSumIsMultiplicative,
                                     int schurId,
                                     const std::string& schurDisplay,
                                     int schurOrder) const
  {
    size_t n = outer.size();
    if (n == 0) return one();
    if (n >= 8 * sizeof(size_t))
      {
        ERROR("Jacobi-Trudi determinant is too large");
        return zero();
      }

    RingElemMatrix matrix(n, RingElemVector(n));
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        {
          int degree = outer[i] - static_cast<int>(i) + static_cast<int>(j);
          matrix[i][j] = schurCompletePlethysm(degree,
                                               inner,
                                               powerSumBasisId0,
                                               powerSumDisplay,
                                               powerSumOrder,
                                               powerSumIsMultiplicative,
                                               schurId,
                                               schurDisplay,
                                               schurOrder);
          if (error()) return zero();
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
            ring_elem term = multiplySchurElements(dp[mask],
                                                   matrix[row][col],
                                                   schurId,
                                                   schurDisplay,
                                                   schurOrder);
            if (error()) return zero();
            if (selectedGreaterThan(mask, col, n) % 2 == 1) term = negate(term);
            size_t next = mask | bit;
            dp[next] = add(dp[next], term);
          }
      }

    return dp[limit - 1];
  }

  bool schurPlethysmToSchur(ring_elem f,
                            ring_elem g,
                            int powerSumBasisId0,
                            const std::string& powerSumDisplay,
                            int powerSumOrder,
                            bool powerSumIsMultiplicative,
                            int targetBasisId,
                            const std::string& targetDisplay,
                            int targetOrder,
                            ring_elem& result) const
  {
    if (targetDisplay != "S") return false;
    Partition outer;
    Partition inner;
    if (!singleSchurPartition(f, targetBasisId, outer)) return false;
    if (!singleSchurPartition(g, targetBasisId, inner)) return false;
    result = schurPlethysmJacobiTrudi(outer,
                                      inner,
                                      powerSumBasisId0,
                                      powerSumDisplay,
                                      powerSumOrder,
                                      powerSumIsMultiplicative,
                                      targetBasisId,
                                      targetDisplay,
                                      targetOrder);
    return !error();
  }

  ring_elem jacobiTrudiBasis(int basisId,
                             const std::string& display,
                             int order,
                             bool isMultiplicative,
                             const Partition& outer,
                             const Partition& inner) const
  {
    return jacobiTrudi(outer, inner, basisId, display, order, isMultiplicative);
  }

  ring_elem toBasis(ring_elem f,
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

  ring_elem plethysm(ring_elem f,
                     ring_elem g,
                     int pBasisId,
                     const std::string& pDisplay,
                     int pOrder,
                     bool pIsMultiplicative) const
  {
    rememberBasis(pBasisId, pDisplay, pOrder, pIsMultiplicative);
    ring_elem fPowerSums = toBasis(f,
                                   pBasisId,
                                   pDisplay,
                                   pOrder,
                                   pIsMultiplicative,
                                   pBasisId,
                                   pDisplay,
                                   pOrder,
                                   pIsMultiplicative);
    if (error()) return zero();
    ring_elem gPowerSums = toBasis(g,
                                   pBasisId,
                                   pDisplay,
                                   pOrder,
                                   pIsMultiplicative,
                                   pBasisId,
                                   pDisplay,
                                   pOrder,
                                   pIsMultiplicative);
    if (error()) return zero();
    return plethysmPowerSums(fPowerSums, gPowerSums);
  }

  ring_elem plethysmToBasis(ring_elem f,
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

    ring_elem specialized;
    if (schurPlethysmToSchur(f,
                             g,
                             pBasisId,
                             pDisplay,
                             pOrder,
                             pIsMultiplicative,
                             targetBasisId,
                             targetDisplay,
                             targetOrder,
                             specialized))
      return specialized;
    if (error()) return zero();

    ring_elem result = plethysm(f, g, pBasisId, pDisplay, pOrder, pIsMultiplicative);
    if (error()) return zero();
    return toBasis(result,
                   pBasisId,
                   pDisplay,
                   pOrder,
                   pIsMultiplicative,
                   targetBasisId,
                   targetDisplay,
                   targetOrder,
                   targetIsMultiplicative);
  }

  int uniformBasisId(ring_elem f) const
  {
    return singleBasisId(f);
  }

  ring_elem omegaInvolution(ring_elem f, M2_arrayint omegaMap, bool useSomega) const
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

  ring_elem straighten(ring_elem f) const
  {
    return straightenElement(f);
  }

  ring_elem hallInnerProduct(ring_elem f,
                             ring_elem g,
                             M2_arrayint innerProductMap) const
  {
    std::map<int, InnerProductTarget> metadata =
        innerProductTargetMap(innerProductMap);
    if (error()) return coefficientRing->zero();
    return hallInnerProductElements(f, g, metadata);
  }

  virtual ring_elem invert(const ring_elem f) const
  {
    (void)f;
    return zero();
  }

  virtual ring_elem divide(const ring_elem f, const ring_elem g) const
  {
    (void)f;
    (void)g;
    return zero();
  }

  virtual void syzygy(const ring_elem a,
                      const ring_elem b,
                      ring_elem &x,
                      ring_elem &y) const
  {
    (void)a;
    (void)b;
    x = zero();
    y = zero();
  }

  virtual ring_elem eval(const RingMap *map,
                         const ring_elem f,
                         int first_var) const
  {
    (void)f;
    (void)first_var;
    return map->get_ring()->zero();
  }
};

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

} // namespace

const Ring *rawSymmetricRing(const Ring *A)
{
  try
    {
      return SymmetricEngineRing::create(A);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

bool rawSymmetricRingsSetHallLittlewoodParameter(const Ring *R,
                                                 const RingElement *t)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return false;
      bool ok = S->setHallLittlewoodParameter(t);
      if (!ok) ERROR("expected a Hall-Littlewood parameter in the coefficient ring");
      return ok;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

bool rawSymmetricRingsRememberBasis(const Ring *R,
                                    int basisId,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return false;
      S->rememberBasisMetadata(basisId,
                               fromM2String(displaySymbol),
                               displayOrder,
                               isMultiplicative);
      return true;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                bool isMultiplicative,
                                                int innerLength,
                                                M2_arrayint index)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->basisElement(basisId,
                                         fromM2String(displaySymbol),
                                         displayOrder,
                                         isMultiplicative,
                                         innerLength,
                                         index);
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      return RingElement::make_raw(S, S->batchSum(elements));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      return RingElement::make_raw(S, S->batchProduct(elements));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsJacobiTrudi(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                bool isMultiplicative,
                                                M2_arrayint outer,
                                                M2_arrayint inner)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->jacobiTrudiBasis(basisId,
                                             fromM2String(displaySymbol),
                                             displayOrder,
                                             isMultiplicative,
                                             partitionFromM2Array(outer),
                                             partitionFromM2Array(inner));
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int powerSumBasisId,
                                            M2_string powerSumDisplaySymbol,
                                            int powerSumDisplayOrder,
                                            bool powerSumIsMultiplicative,
                                            int targetBasisId,
                                            M2_string targetDisplaySymbol,
                                            int targetDisplayOrder,
                                            bool targetIsMultiplicative)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->toBasis(f->get_value(),
                                    powerSumBasisId,
                                    fromM2String(powerSumDisplaySymbol),
                                    powerSumDisplayOrder,
                                    powerSumIsMultiplicative,
                                    targetBasisId,
                                    fromM2String(targetDisplaySymbol),
                                    targetDisplayOrder,
                                    targetIsMultiplicative);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsPlethysm(const RingElement *f,
                                             const RingElement *g,
                                             int powerSumBasisId,
                                             M2_string powerSumDisplaySymbol,
                                             int powerSumDisplayOrder,
                                             bool powerSumIsMultiplicative)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->plethysm(f->get_value(),
                                     g->get_value(),
                                     powerSumBasisId,
                                     fromM2String(powerSumDisplaySymbol),
                                     powerSumDisplayOrder,
                                     powerSumIsMultiplicative);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsPlethysmToBasis(const RingElement *f,
                                                    const RingElement *g,
                                                    int powerSumBasisId,
                                                    M2_string powerSumDisplaySymbol,
                                                    int powerSumDisplayOrder,
                                                    bool powerSumIsMultiplicative,
                                                    int targetBasisId,
                                                    M2_string targetDisplaySymbol,
                                                    int targetDisplayOrder,
                                                    bool targetIsMultiplicative)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->plethysmToBasis(f->get_value(),
                                            g->get_value(),
                                            powerSumBasisId,
                                            fromM2String(powerSumDisplaySymbol),
                                            powerSumDisplayOrder,
                                            powerSumIsMultiplicative,
                                            targetBasisId,
                                            fromM2String(targetDisplaySymbol),
                                            targetDisplayOrder,
                                            targetIsMultiplicative);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

int rawSymmetricRingsSingleBasisId(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return -1;
      return S->uniformBasisId(f->get_value());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return -1;
    }
}

const RingElement *rawSymmetricRingsOmega(const RingElement *f,
                                          M2_arrayint omegaMap,
                                          bool useSomega)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->omegaInvolution(f->get_value(), omegaMap, useSomega);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsStraighten(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->straighten(f->get_value());
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsHallInnerProduct(const RingElement *f,
                                                    const RingElement *g,
                                                    M2_arrayint innerProductMap)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->hallInnerProduct(f->get_value(),
                                             g->get_value(),
                                             innerProductMap);
      if (error()) return nullptr;
      return RingElement::make_raw(S->getCoefficientRing(), result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

int rawSymmetricRingsTermCount(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return 0;
      (void) S;
      return static_cast<int>(polyValue(f->get_value())->terms.size());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return 0;
    }
}

const RingElement *rawSymmetricRingsTermCoefficient(const RingElement *f, int i)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      const auto *poly = polyValue(f->get_value());
      if (i < 0 || static_cast<size_t>(i) >= poly->terms.size())
        {
          ERROR("symmetric-ring term index out of range");
          return nullptr;
        }
      const Ring *A = S->getCoefficientRing();
      return RingElement::make_raw(A, A->copy(poly->terms[i].coeff));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsTermMonomial(const RingElement *f, int i)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      (void) S;
      const auto *poly = polyValue(f->get_value());
      if (i < 0 || static_cast<size_t>(i) >= poly->terms.size())
        {
          ERROR("symmetric-ring term index out of range");
          return nullptr;
        }
      const auto& data = poly->terms[i].monomial.data;
      M2_arrayint result = M2_makearrayint(static_cast<int>(data.size()));
      for (size_t j = 0; j < data.size(); ++j)
        result->array[j] = data[j];
      return result;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_string rawSymmetricRingsElementToString(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return toM2String(S->elementString(f->get_value()));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

int rawSymmetricRingsElementWeight(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return -1;
      return S->elementWeight(f->get_value());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return -1;
    }
}
