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
  std::map<std::string, int> memo;
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
  mutable std::map<int, bool> multiplicativeBases;
  mutable int powerSumBasisId = -1;
  mutable std::map<int, ring_elem> hToPowerSumCache;
  mutable std::map<int, ring_elem> eToPowerSumCache;
  mutable std::map<std::string, ring_elem> schurToPowerSumCache;
  mutable std::map<int, ring_elem> powerSumToCompleteCache;
  mutable std::map<int, ring_elem> powerSumToElementaryCache;
  mutable std::map<std::string, ring_elem> powerSumToSchurCache;

  bool isMultiplicativeBasis(int basisId) const
  {
    auto it = multiplicativeBases.find(basisId);
    return it != multiplicativeBases.end() && it->second;
  }

  bool isPowerSumBasis(int basisId) const
  {
    return powerSumBasisId >= 0 && basisId == powerSumBasisId;
  }

  void rememberBasis(int basisId,
                     const std::string& display,
                     bool isMultiplicative) const
  {
    if (!display.empty()) basisDisplays[basisId] = display;
    multiplicativeBases[basisId] = isMultiplicative;
    if (display == "p") powerSumBasisId = basisId;
  }

  void rememberBasesFrom(const SymmetricEngineRing *R) const
  {
    for (const auto& item : R->basisDisplays) basisDisplays[item.first] = item.second;
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
    for (const auto& item : basisDisplays)
      {
        (void)item;
      }
    if (basisId == basisIdForDisplay("p")) return 10;
    if (basisId == basisIdForDisplay("h")) return 20;
    if (basisId == basisIdForDisplay("e")) return 30;
    if (basisId == basisIdForDisplay("S")) return 60;
    return 100;
  }

  ring_elem rationalCoefficient(long numerator, long denominator) const
  {
    mpq_t q;
    mpq_init(q);
    mpq_set_si(q, numerator, denominator);
    mpq_canonicalize(q);
    ring_elem result;
    if (!coefficientRing->from_rational(q, result)) result = coefficientRing->zero();
    mpq_clear(q);
    return result;
  }

  ring_elem basisElementFromIndex(int basisId,
                                  const std::string& display,
                                  int order,
                                  bool isMultiplicative,
                                  const Partition& index) const
  {
    rememberBasis(basisId, display, isMultiplicative);
    auto result = new SymmetricRingPoly;
    SymmetricMonomial monomial;
    appendAtomBlock(monomial, makeAtomBlock(order, basisId, 0, index));
    result->terms.push_back({coefficientRing->one(), canonicalMonomial(monomial)});
    return makePolyValue(result);
  }

  ring_elem scaled(ring_elem coeff, ring_elem f) const
  {
    return makePolyValue(multByCoefficient(coeff, polyValue(f)));
  }

  ring_elem hPartToPowerSums(int n) const
  {
    auto cached = hToPowerSumCache.find(n);
    if (cached != hToPowerSumCache.end()) return copyPolyValue(polyValue(cached->second));
    ring_elem result = zero();
    for (const auto& mu : partitionsOf(n))
      {
        ring_elem coeff = rationalCoefficient(1, zValue(mu));
        ring_elem term = basisElementFromIndex(powerSumBasisId, "p", 10, true, mu);
        result = add(result, scaled(coeff, term));
      }
    hToPowerSumCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

  ring_elem ePartToPowerSums(int n) const
  {
    auto cached = eToPowerSumCache.find(n);
    if (cached != eToPowerSumCache.end()) return copyPolyValue(polyValue(cached->second));
    ring_elem result = zero();
    for (const auto& mu : partitionsOf(n))
      {
        long sign = ((n - static_cast<int>(mu.size())) % 2 == 0) ? 1 : -1;
        ring_elem coeff = rationalCoefficient(sign, zValue(mu));
        ring_elem term = basisElementFromIndex(powerSumBasisId, "p", 10, true, mu);
        result = add(result, scaled(coeff, term));
      }
    eToPowerSumCache[n] = result;
    return copyPolyValue(polyValue(result));
  }

  ring_elem schurToPowerSums(const Partition& lambda) const
  {
    std::string cacheKey = partitionKey(lambda);
    auto cached = schurToPowerSumCache.find(cacheKey);
    if (cached != schurToPowerSumCache.end()) return copyPolyValue(polyValue(cached->second));
    int n = 0;
    for (int part : lambda) n += part;
    ring_elem result = zero();
    for (const auto& mu : partitionsOf(n))
      {
        int chi = characterValue(lambda, mu);
        if (chi == 0) continue;
        ring_elem coeff = rationalCoefficient(chi, zValue(mu));
        ring_elem term = basisElementFromIndex(powerSumBasisId, "p", 10, true, mu);
        result = add(result, scaled(coeff, term));
      }
    schurToPowerSumCache[cacheKey] = result;
    return copyPolyValue(polyValue(result));
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

  ring_elem powerSumToSchur(const Partition& mu, int schurId, int schurOrder) const
  {
    std::string cacheKey = partitionKey(mu);
    auto cached = powerSumToSchurCache.find(cacheKey);
    if (cached != powerSumToSchurCache.end()) return copyPolyValue(polyValue(cached->second));
    int n = 0;
    for (int part : mu) n += part;
    ring_elem result = zero();
    for (const auto& lambda : partitionsOf(n))
      {
        int chi = characterValue(lambda, mu);
        if (chi == 0) continue;
        ring_elem term = basisElementFromIndex(schurId, "S", schurOrder, false, lambda);
        result = add(result, scaled(coefficientRing->from_long(chi), term));
      }
    powerSumToSchurCache[cacheKey] = result;
    return copyPolyValue(polyValue(result));
  }

  Partition atomIndex(const SymmetricMonomial& monomial, size_t pos) const
  {
    Partition result;
    int n = atomIndexLengthAt(monomial, pos);
    result.reserve(n);
    for (int i = 0; i < n; ++i) result.push_back(monomial.data[pos + atomHeaderSize + i]);
    return result;
  }

  ring_elem atomToPowerSums(const SymmetricMonomial& monomial, size_t pos) const
  {
    int basisId = atomBasisIdAt(monomial, pos);
    std::string display = displayForBasis(basisId);
    if (atomIsSkewAt(monomial, pos))
      {
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

    if (display == "S")
      return schurToPowerSums(index);

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

    if (targetDisplay == "S")
      return powerSumToSchur(index, targetBasisId, targetDisplayOrder);

    ERROR("basis conversion from power sums is not implemented for basis ", targetDisplay.c_str());
    return zero();
  }

  ring_elem powerSumsToTarget(ring_elem f,
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
    result->terms.insert(result->terms.end(), poly->terms.begin(), poly->terms.end());
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

  void addToAccumulator(std::map<std::vector<int>, ring_elem>& accumulator,
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
      const std::map<std::vector<int>, ring_elem>& accumulator) const
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
  explicit SymmetricEngineRing(const Ring *A) : coefficientRing(A) {}

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
    rememberBasis(basisId, display, isMultiplicative);
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
        std::map<std::vector<int>, ring_elem> accumulator;
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
        std::map<std::vector<int>, ring_elem> accumulator;
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
    if (coefficientRing->is_equal(coeff, coefficientRing->one()))
      {
        auto result = new SymmetricRingPoly;
        result->terms.reserve(poly->terms.size());
        result->terms.insert(result->terms.end(), poly->terms.begin(), poly->terms.end());
        return result;
      }
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
    rememberBasis(pBasisId, pDisplay, pIsMultiplicative);
    rememberBasis(targetBasisId, targetDisplay, targetIsMultiplicative);
    ring_elem inPowerSums = elementToPowerSums(f);
    if (error()) return zero();
    return powerSumsToTarget(inPowerSums,
                             targetBasisId,
                             targetDisplay,
                             targetOrder,
                             targetIsMultiplicative);
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
