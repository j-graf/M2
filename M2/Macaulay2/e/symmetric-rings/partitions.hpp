// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_PARTITIONS_HPP_
#define M2_SYMMETRIC_RINGS_PARTITIONS_HPP_

#include "interface/symmetric-rings.h"

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace symmetric_rings {

using Partition = std::vector<int>;

constexpr int unknownCharacterValue = std::numeric_limits<int>::min();

struct CharacterTable
{
  std::vector<Partition> partitions;
  std::vector<long> zValues;
  std::map<Partition, size_t> partitionRows;
  mutable std::vector<std::vector<int>> values;
};

std::string partitionKey(const Partition& p);
Partition normalizePartition(const Partition& p);
int partitionWeight(const Partition& p);
int partitionLength(const Partition& p);
bool lexLessPartition(const Partition& a, const Partition& b);
Partition trimTrailingZerosPartition(const Partition& p);
Partition conjugatePartition(const Partition& p);
std::pair<int, Partition> straightenSchurIndex(const Partition& alpha);
std::vector<Partition> partitionsOf(int n);
long pToMonomialCoefficient(const Partition& lambda, const Partition& mu);
long zValue(const Partition& lambda);
Partition partitionFromM2Array(M2_arrayint a);
int characterValue(const Partition& lambda, const Partition& mu);

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
