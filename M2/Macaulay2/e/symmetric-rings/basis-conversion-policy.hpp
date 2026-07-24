// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_CONVERSION_POLICY_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_CONVERSION_POLICY_HPP_

#include "symmetric-rings/partitions.hpp"

#include <cstddef>

namespace symmetric_rings {

// ============================================================================
// Conversion-Selection Facts
// ============================================================================
// These helpers provide performance evidence only; they never change a plan's
// mathematical applicability.

size_t partitionCountForConversionSelection(int weight);

bool powerSumIndexHasMostlyShortCycles(
    const Partition& index,
    int weight);

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
