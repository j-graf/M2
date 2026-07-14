// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_OPERATION_RECORDS_HPP_
#define M2_SYMMETRIC_RINGS_OPERATION_RECORDS_HPP_

#include "symmetric-rings/partitions.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace symmetric_rings {

// Compact records shared between the central engine state and one operation
// module. They contain no storage behavior or route-selection policy.

struct SchurConversionRecipeEntry
{
  Partition lambda;
  std::vector<std::pair<size_t, int>> contributions;
};

struct PartitionCoefficientTerm
{
  Partition partition;
  long coefficient;
};

struct OmegaTarget
{
  int basisId;
};

enum class InnerProductPairingKind
{
  Dual = 1,
  PowerSum = 2
};

struct InnerProductTarget
{
  int dualBasisId;
  InnerProductPairingKind kind;
};

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
