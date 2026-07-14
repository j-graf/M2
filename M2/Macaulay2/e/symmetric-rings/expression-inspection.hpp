// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_EXPRESSION_INSPECTION_HPP_
#define M2_SYMMETRIC_RINGS_EXPRESSION_INSPECTION_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  // ============================================================================
  // Expression Inspection
  // ============================================================================
  // These shape probes do not perform conversion or pairing.

  bool singleScaledBasisElement(ring_elem f,
                                int basisId,
                                Partition& index,
                                ring_elem& coefficient) const;
  CoeffMap coefficientsInBasis(ring_elem f, int basisId) const;
  bool coefficientsInBasisIfPossible(ring_elem f,
                                      int basisId,
                                      CoeffMap& result) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
