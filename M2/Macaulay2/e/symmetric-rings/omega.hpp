// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_OMEGA_HPP_
#define M2_SYMMETRIC_RINGS_OMEGA_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  ring_elem omegaPowerSums(ring_elem f) const;
  ring_elem omegaSchurBasisElementAsSchur(const Partition& alpha) const;
  ring_elem omegaBasisElementDirect(const SymmetricMonomial& monomial,
                              size_t pos,
                              const std::map<int, OmegaTarget>& omegaTargets,
                              bool useSomega) const;
  ring_elem omegaMonomial(const SymmetricMonomial& monomial,
                            const std::map<int, OmegaTarget>& omegaTargets,
                            bool useSomega) const;
  std::map<int, OmegaTarget> omegaTargetMap(M2_arrayint omegaMap) const;

 public:
  ring_elem omegaInvolution(ring_elem f, M2_arrayint omegaMap, bool useSomega) const;

 private:

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
