// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_PARTITIONS_HPP_
#define M2_SYMMETRIC_RINGS_PARTITIONS_HPP_

#include <gmpxx.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace symmetric_rings {

using Partition = std::vector<int>;

struct CharacterTable
{
  std::vector<Partition> partitions;
  std::vector<mpz_class> zValues;
  std::map<Partition, size_t> partitionRows;
  // Rows are materialized lazily. Conversion kernels normally consume a
  // complete row, so row storage avoids an ordered-map lookup and allocation
  // for every character-table cell.
  mutable std::map<size_t, std::vector<mpz_class>> rows;
};

std::string partitionKey(const Partition& p);
Partition normalizePartition(const Partition& p);
int partitionWeight(const Partition& p);
int partitionLength(const Partition& p);
bool hasNegativeTailWeight(const Partition& index);
bool isPartitionIndex(const Partition& p);
bool isHookPartition(const Partition& p);
bool isRectanglePartition(const Partition& p);
bool isSelfConjugatePartition(const Partition& p);
int partitionPart(const Partition& p, size_t i);
bool partitionContains(const Partition& outer, const Partition& inner);
bool lexLessPartition(const Partition& a, const Partition& b);
Partition trimTrailingZerosPartition(const Partition& p);
Partition conjugatePartition(const Partition& p);
std::pair<int, Partition> straightenSchurIndex(const Partition& alpha);
// Returns limit + 1 when p(n) exceeds limit.  This counts without
// materializing the partitions and stops as soon as the limit is crossed.
size_t partitionCountUpToLimit(int n, size_t limit);
std::vector<Partition> partitionsOf(int n);
mpz_class pToMonomialCoefficientWithLimit(const Partition& lambda,
                                          const Partition& mu,
                                          size_t maxStates,
                                          bool& limitExceeded);
mpz_class zValue(const Partition& lambda);
mpz_class characterValueWithLimit(const Partition& lambda,
                                  const Partition& mu,
                                  size_t maxStates,
                                  bool& limitExceeded);
std::vector<mpz_class> characterRowWithLimit(
    const Partition& lambda,
    const std::vector<Partition>& cycleTypes,
    size_t maxStatesPerValue,
    bool& limitExceeded);

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
