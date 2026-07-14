// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_PRESENTATION_HPP_
#define M2_SYMMETRIC_RINGS_PRESENTATION_HPP_

// Public declaration fragment included inside SymmetricEngineRing.

  // ============================================================================
  // Presentation And Display
  // ============================================================================

  std::string displayIndex(const SymmetricMonomial& monomial, size_t pos) const;
  std::string displayBasisElement(const SymmetricMonomial& monomial, size_t pos) const;
  std::string displayMonomial(const SymmetricMonomial& monomial) const;
  int compareBasisIdsForPresentation(int aBasis, int bBasis) const;
  std::vector<size_t> presentationAtomPositions(const SymmetricMonomial& monomial) const;
  std::vector<int> presentationMonomialData(const SymmetricMonomial& monomial) const;
  std::vector<size_t> presentationTermOrder(ring_elem f) const;
  std::string elementString(ring_elem f, int maxTerms = -1) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
