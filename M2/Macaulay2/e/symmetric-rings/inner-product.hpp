// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_INNER_PRODUCT_HPP_
#define M2_SYMMETRIC_RINGS_INNER_PRODUCT_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  CoeffMap coefficientsInBasis(ring_elem f, int basisId) const;
  bool coefficientsInBasisIfPossible(ring_elem f,
                                        int basisId,
                                        CoeffMap& result) const;
  ring_elem coefficientPairing(const CoeffMap& fCoeffs,
                                 const CoeffMap& gCoeffs) const;
  ring_elem powerSumInnerProductFactor(const Partition& lambda) const;
  ring_elem powerSumPairing(const CoeffMap& fCoeffs,
                              const CoeffMap& gCoeffs) const;
  std::map<int, InnerProductTarget> innerProductTargetMap(M2_arrayint innerProductMap) const;
  bool directHallInnerProductFromMetadata(
        ring_elem f,
        ring_elem g,
        const std::map<int, InnerProductTarget>& metadata,
        ring_elem& result) const;
  bool directPowerSumSchurInnerProduct(ring_elem f,
                                         ring_elem g,
                                         ring_elem& result) const;
  ring_elem hallInnerProductElements(ring_elem f,
                                       ring_elem g,
                                       const std::map<int, InnerProductTarget>& metadata = {}) const;
 public:
  ring_elem hallInnerProduct(ring_elem f,
                               ring_elem g,
                               M2_arrayint innerProductMap) const;

 private:

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
