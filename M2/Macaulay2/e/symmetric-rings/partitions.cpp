// Copyright 2026

#include "symmetric-rings/partitions.hpp"

#include "error.h"
#include "exceptions.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <sstream>
#include <variant>

namespace symmetric_rings {

namespace {

// Most combinatorial coefficients at practical weights fit in a machine
// integer. Keep those recursive states allocation-free and promote exactly
// when a checked multiply or add overflows. Public results remain mpz_class,
// so this changes representation cost without changing mathematical range.
class CheckedExactInteger
{
 public:
  CheckedExactInteger(long value = 0) : value_(value) {}

  void addMultiple(long multiplier, const CheckedExactInteger& other)
  {
    if (std::holds_alternative<long>(value_) &&
        std::holds_alternative<long>(other.value_))
      {
        long product = 0;
        long sum = 0;
        if (!__builtin_mul_overflow(
                multiplier,
                std::get<long>(other.value_),
                &product) &&
            !__builtin_add_overflow(
                std::get<long>(value_), product, &sum))
          {
            value_ = sum;
            return;
          }
      }
    mpz_class sum = asMpz();
    mpz_class contribution = other.asMpz();
    contribution *= multiplier;
    sum += contribution;
    value_ = std::move(sum);
  }

  mpz_class asMpz() const
  {
    return std::holds_alternative<long>(value_)
        ? mpz_class(std::get<long>(value_))
        : std::get<mpz_class>(value_);
  }

 private:
  std::variant<long, mpz_class> value_;
};

} // namespace

// ============================================================================
// Basic Partition Operations
// ============================================================================

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
  long long result = 0;
  for (int part : p)
    {
      result += part;
      if (result < std::numeric_limits<int>::min() ||
          result > std::numeric_limits<int>::max())
        throw exc::engine_error(
            "partition weight exceeds the supported integer range");
    }
  return static_cast<int>(result);
}

int partitionLength(const Partition& p)
{
  int result = 0;
  for (int part : p)
    if (part > 0) ++result;
  return result;
}

bool hasNegativeTailWeight(const Partition& index)
{
  long long tailWeight = 0;
  for (auto part = index.rbegin(); part != index.rend(); ++part)
    {
      tailWeight += *part;
      if (tailWeight < 0) return true;
    }
  return false;
}

bool isPartitionIndex(const Partition& p)
{
  return trimTrailingZerosPartition(p) == normalizePartition(p);
}

bool isHookPartition(const Partition& p)
{
  if (!isPartitionIndex(p)) return false;
  const Partition index = trimTrailingZerosPartition(p);
  if (index.empty()) return false;
  for (size_t i = 1; i < index.size(); ++i)
    if (index[i] > 1) return false;
  return true;
}

bool isRectanglePartition(const Partition& p)
{
  if (!isPartitionIndex(p)) return false;
  const Partition index = trimTrailingZerosPartition(p);
  if (index.empty()) return false;
  return std::all_of(
      index.begin(),
      index.end(),
      [&](int part) { return part == index.front(); });
}

bool isSelfConjugatePartition(const Partition& p)
{
  if (!isPartitionIndex(p)) return false;
  const Partition index = trimTrailingZerosPartition(p);
  return conjugatePartition(index) == index;
}

int partitionPart(const Partition& p, size_t i)
{
  return i < p.size() ? p[i] : 0;
}

bool partitionContains(const Partition& outer, const Partition& inner)
{
  for (size_t i = 0; i < inner.size(); ++i)
    if (partitionPart(outer, i) < inner[i]) return false;
  return true;
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

// ============================================================================
// Enumeration And Power-Sum Statistics
// ============================================================================

size_t partitionCountUpToLimit(int n, size_t limit)
{
  if (n < 0) return 0;
  std::vector<mpz_class> counts;
  counts.push_back(1);
  const mpz_class threshold(limit);
  for (int degree = 1; degree <= n; ++degree)
    {
      mpz_class count = 0;
      for (int k = 1;; ++k)
        {
          const int pentagonal1 = k * (3 * k - 1) / 2;
          if (pentagonal1 > degree) break;
          const int sign = (k % 2 == 1) ? 1 : -1;
          count += sign * counts[degree - pentagonal1];
          const int pentagonal2 = k * (3 * k + 1) / 2;
          if (pentagonal2 <= degree)
            count += sign * counts[degree - pentagonal2];
        }
      if (count > threshold)
        return limit == std::numeric_limits<size_t>::max() ? limit : limit + 1;
      counts.push_back(count);
    }
  return counts[n].get_ui();
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

CheckedExactInteger assignmentCountRec(
    const Partition& parts,
    size_t pos,
    const Partition& targets,
    std::map<
        std::pair<size_t, Partition>,
        CheckedExactInteger>& memo,
    size_t maxStates,
    size_t& states,
    bool& limitExceeded)
{
  if (limitExceeded) return CheckedExactInteger{};
  if (states >= maxStates)
    {
      limitExceeded = true;
      return CheckedExactInteger{};
    }
  ++states;
  if (pos == parts.size())
    return CheckedExactInteger(targets.empty() ? 1 : 0);
  auto key = std::make_pair(pos, targets);
  auto cached = memo.find(key);
  if (cached != memo.end()) return cached->second;

  CheckedExactInteger total;
  int part = parts[pos];
  size_t i = 0;
  while (i < targets.size())
    {
      size_t j = i + 1;
      while (j < targets.size() && targets[j] == targets[i]) ++j;
      if (targets[i] >= part)
      {
        Partition next = targets;
        next[i] -= part;
        std::sort(next.begin(), next.end(), std::greater<int>());
        while (!next.empty() && next.back() == 0) next.pop_back();
        CheckedExactInteger child =
            assignmentCountRec(
                parts,
                pos + 1,
                next,
                memo,
                maxStates,
                states,
                limitExceeded);
        if (limitExceeded) return CheckedExactInteger{};
        total.addMultiple(
            static_cast<long>(j - i), child);
      }
      i = j;
    }
  memo[key] = total;
  return total;
}

mpz_class pToMonomialCoefficientWithLimit(
    const Partition& lambda,
    const Partition& mu,
    size_t maxStates,
    bool& limitExceeded)
{
  limitExceeded = false;
  Partition normalizedLambda = normalizePartition(lambda);
  Partition normalizedMu = normalizePartition(mu);
  if (partitionWeight(normalizedLambda) != partitionWeight(normalizedMu)) return 0;
  std::map<
      std::pair<size_t, Partition>,
      CheckedExactInteger> memo;
  size_t states = 0;
  return assignmentCountRec(
             normalizedLambda,
             0,
             normalizedMu,
             memo,
             maxStates,
             states,
             limitExceeded)
      .asMpz();
}

mpz_class zValue(const Partition& lambda)
{
  std::map<int, int> multiplicities;
  for (int part : lambda) multiplicities[part]++;
  mpz_class result = 1;
  for (const auto& item : multiplicities)
    {
      for (int i = 0; i < item.second; ++i) result *= item.first;
      for (int i = 2; i <= item.second; ++i) result *= i;
    }
  return result;
}

// ============================================================================
// Rim Hooks And Character Recursion
// ============================================================================
// Rim-hook removal validates connectivity and the absence of a two-by-two block.

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
                std::vector<RimHookRemoval>& result,
                size_t maxStates,
                size_t& states,
                bool& limitExceeded)
{
  if (limitExceeded) return;
  if (states >= maxStates)
    {
      limitExceeded = true;
      return;
    }
  ++states;
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
      rimHookRec(
          lambda,
          row + 1,
          remainingSize - count,
          removed,
          result,
          maxStates,
          states,
          limitExceeded);
      if (limitExceeded) return;
    }
  removed[row] = 0;
}

std::vector<RimHookRemoval> rimHookRemovals(
    const Partition& lambda,
    int size,
    size_t maxStates,
    size_t& states,
    bool& limitExceeded)
{
  std::vector<RimHookRemoval> result;
  std::vector<int> removed(lambda.size(), 0);
  rimHookRec(
      lambda,
      0,
      size,
      removed,
      result,
      maxStates,
      states,
      limitExceeded);
  return result;
}

CheckedExactInteger characterValueMemo(
    const Partition& lambda,
    const Partition& mu,
    std::map<std::string, CheckedExactInteger>& memo,
    size_t maxStates,
    size_t& states,
    bool& limitExceeded)
{
  if (limitExceeded) return CheckedExactInteger{};
  if (states >= maxStates)
    {
      limitExceeded = true;
      return CheckedExactInteger{};
    }
  ++states;
  if (mu.empty())
    return CheckedExactInteger(lambda.empty() ? 1 : 0);
  std::string key = partitionKey(lambda) + "|" + partitionKey(mu);
  auto it = memo.find(key);
  if (it != memo.end()) return it->second;
  CheckedExactInteger total;
  int hookSize = mu.front();
  Partition rest(mu.begin() + 1, mu.end());
  for (const auto& removal :
       rimHookRemovals(
           lambda,
           hookSize,
           maxStates,
           states,
           limitExceeded))
    {
      int sign = (removal.height % 2 == 0) ? 1 : -1;
      CheckedExactInteger child =
          characterValueMemo(
              removal.remaining,
              rest,
              memo,
              maxStates,
              states,
              limitExceeded);
      if (limitExceeded) return CheckedExactInteger{};
      total.addMultiple(sign, child);
    }
  // A row may reuse more distinct subproblems in aggregate than the
  // per-character state limit. Once the shared memo reaches that limit,
  // continue without adding entries; the state counter above remains the
  // authoritative limit for each requested character value.
  if (memo.size() < maxStates) memo[key] = total;
  return total;
}

mpz_class characterValueWithLimit(const Partition& lambda,
                                  const Partition& mu,
                                  size_t maxStates,
                                  bool& limitExceeded)
{
  std::map<std::string, CheckedExactInteger> memo;
  size_t states = 0;
  limitExceeded = false;
  return characterValueMemo(
             lambda, mu, memo, maxStates, states, limitExceeded)
      .asMpz();
}

std::vector<mpz_class> characterRowWithLimit(
    const Partition& lambda,
    const std::vector<Partition>& cycleTypes,
    size_t maxStatesPerValue,
    bool& limitExceeded)
{
  // Recursive character subproblems are independent of the coefficient ring.
  // Share them across the complete row while retaining the configured
  // per-character recursion bound used by characterValueWithLimit.
  std::map<std::string, CheckedExactInteger> memo;
  std::vector<mpz_class> result;
  result.reserve(cycleTypes.size());
  limitExceeded = false;
  for (const auto& mu : cycleTypes)
    {
      size_t states = 0;
      bool valueLimitExceeded = false;
      result.push_back(
          characterValueMemo(
              lambda,
              mu,
              memo,
              maxStatesPerValue,
              states,
              valueLimitExceeded)
              .asMpz());
      if (valueLimitExceeded)
        {
          limitExceeded = true;
          return {};
        }
    }
  return result;
}


} // namespace symmetric_rings
