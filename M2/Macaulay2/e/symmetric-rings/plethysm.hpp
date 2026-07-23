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

  enum class PlethysmToBasisRoute
  {
    ViaSchurAdamsJacobiTrudi,
    ViaPowerSumsThenBasisConversion
  };

  bool singleSchurPartition(ring_elem f, int schurId, Partition& lambda) const;
  bool schurPlethysmToSchurViaAdamsJacobiTrudiApplicable(
      ring_elem f,
      ring_elem g,
      int targetBasisId,
      Partition& outer,
      Partition& inner) const;
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
// ============================================================================
// Plethysm-To-Basis Selection And Execution
// ============================================================================

  PlethysmToBasisRoute selectPlethysmToBasisRoute(
      ring_elem f,
      ring_elem g,
      int targetBasisId) const;
  const char *plethysmToBasisRouteName(PlethysmToBasisRoute route) const;
  void tracePlethysmToBasisSelection(
      PlethysmToBasisRoute route,
      int targetBasisId) const;
  ring_elem executePlethysmToBasisRoute(
      PlethysmToBasisRoute route,
      ring_elem f,
      ring_elem g,
      int targetBasisId) const;
  void requirePlethysmWithinWeightLimit(ring_elem f, ring_elem g) const;

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
