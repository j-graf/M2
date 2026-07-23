// Copyright 2026

#include "symmetric-rings/basis-conversion-policy.hpp"

#include <limits>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Shared Conversion-Selection Facts
// ============================================================================

size_t partitionCountForConversionSelection(int weight)
{
  if (weight < 0) return 0;
  std::vector<size_t> counts(static_cast<size_t>(weight) + 1, 0);
  counts[0] = 1;
  for (int part = 1; part <= weight; ++part)
    for (int total = part; total <= weight; ++total)
      {
        const size_t addend =
            counts[static_cast<size_t>(total - part)];
        size_t& value = counts[static_cast<size_t>(total)];
        if (value > std::numeric_limits<size_t>::max() - addend)
          value = std::numeric_limits<size_t>::max();
        else
          value += addend;
      }
  return counts[static_cast<size_t>(weight)];
}

bool completeFriendlyPowerSumIndexForConversionSelection(
    const Partition& index,
    int weight)
{
  if (weight <= 0 || index.size() < 4) return false;
  int smallCycleWeight = 0;
  for (int part : index)
    if (part <= 4) smallCycleWeight += part;
  // Repeated short cycles create broad rim-hook frontiers, but their
  // logarithmic complete expansions remain narrow. Hybrid conversion has a
  // fixed split-and-merge cost, so only strongly short-cycle indices should
  // enter its complete-friendly subgroup.
  return 4 * smallCycleWeight >= 3 * weight;
}

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
