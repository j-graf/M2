// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_PLETHYSM_HPP_
#define M2_SYMMETRIC_RINGS_PLETHYSM_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Adams-Operation Plethysm
// ============================================================================

  ring_elem powerSumsViaAdamsOperation(ring_elem f, int multiplier) const;
  ring_elem powerSumPlethysmViaAdamsOperations(ring_elem fPowerSums, ring_elem gPowerSums) const;

// ============================================================================
// Specialized Schur Plethysm
// ============================================================================

  bool singleSchurPartition(ring_elem f, int schurId, Partition& lambda) const;
  std::string completePlethysmCacheKey(const std::string& algorithm,
                                         int n,
                                         int schurId,
                                         const Partition& inner) const;
  ring_elem completePlethysmViaAdamsRecurrence(int n,
                                    const Partition& inner,
                                    ring_elem innerPowerSums,
                                    int schurId,
                                    const std::string& schurDisplay,
                                    int schurOrder) const;
  ring_elem schurPlethysmToSchurViaAdamsJacobiTrudi(const Partition& outer,
                                       const Partition& inner,
                                       int schurId,
                                       const std::string& schurDisplay,
                                       int schurOrder) const;
  bool trySchurPlethysmToSchurViaAdamsJacobiTrudi(ring_elem f,
                              ring_elem g,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetOrder,
                              ring_elem& result) const;

// ============================================================================
// Public Plethysm Entry Points
// ============================================================================

 public:
  ring_elem plethysm(ring_elem f, ring_elem g) const;
  ring_elem plethysmToBasisDispatch(ring_elem f,
                              ring_elem g,
                              int targetBasisId) const;

 private:

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
