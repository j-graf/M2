// Copyright 2026

#include "symmetric-rings/partitions.hpp"

#include <algorithm>
#include <functional>
#include <sstream>

namespace symmetric_rings {

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

long assignmentCountRec(
    const Partition& parts,
    size_t pos,
    const Partition& targets,
    std::map<std::pair<size_t, Partition>, long>& memo)
{
  if (pos == parts.size()) return targets.empty() ? 1 : 0;
  auto key = std::make_pair(pos, targets);
  auto cached = memo.find(key);
  if (cached != memo.end()) return cached->second;

  long total = 0;
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
        total += static_cast<long>(j - i) *
                 assignmentCountRec(parts, pos + 1, next, memo);
      }
      i = j;
    }
  memo[key] = total;
  return total;
}

long pToMonomialCoefficient(const Partition& lambda, const Partition& mu)
{
  Partition normalizedLambda = normalizePartition(lambda);
  Partition normalizedMu = normalizePartition(mu);
  if (partitionWeight(normalizedLambda) != partitionWeight(normalizedMu)) return 0;
  std::map<std::pair<size_t, Partition>, long> memo;
  return assignmentCountRec(normalizedLambda, 0, normalizedMu, memo);
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


} // namespace symmetric_rings
