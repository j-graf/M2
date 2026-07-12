// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_INNER_PRODUCT_KERNELS_HPP_
#define M2_SYMMETRIC_RINGS_INNER_PRODUCT_KERNELS_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  // ============================================================================
  // Basis-Expansion Inspection
  // ============================================================================

  bool singleScaledBasisElement(ring_elem f,
                                int basisId,
                                Partition& index,
                                ring_elem& coefficient) const;
  CoeffMap coefficientsInBasis(ring_elem f, int basisId) const;
  bool coefficientsInBasisIfPossible(ring_elem f,
                                      int basisId,
                                      CoeffMap& result) const;

  // ============================================================================
  // Diagonal Pairing Kernels
  // ============================================================================

  ring_elem coefficientPairing(const CoeffMap& fCoeffs,
                               const CoeffMap& gCoeffs) const;
  ring_elem powerSumDiagonalFactor(
      const Partition& lambda,
      const InnerProductContext& context) const;
  ring_elem powerSumPairing(const CoeffMap& fCoeffs,
                            const CoeffMap& gCoeffs,
                            const InnerProductContext& context) const;

  // ============================================================================
  // Kostka-Number Pairings
  // ============================================================================

  long kostkaNumberViaSemistandardTableaux(
      const Partition& shape,
      const Partition& content) const;
  ring_elem kostkaInnerProductForBasisElements(
      ring_elem schurTypeBasisElement,
      int schurTypeBasisId,
      ring_elem multiplicativeBasisElement,
      int multiplicativeBasisId,
      bool conjugateShape) const;
  ring_elem schurCompleteInnerProductViaKostkaNumbers(
      ring_elem schurBasisElement,
      ring_elem completeBasisElement) const;
  ring_elem schurElementaryInnerProductViaConjugateKostkaNumbers(
      ring_elem schurBasisElement,
      ring_elem elementaryBasisElement) const;
  ring_elem schurOmegaCompleteInnerProductViaConjugateKostkaNumbers(
      ring_elem schurOmegaBasisElement,
      ring_elem completeBasisElement) const;
  ring_elem schurOmegaElementaryInnerProductViaKostkaNumbers(
      ring_elem schurOmegaBasisElement,
      ring_elem elementaryBasisElement) const;

  // ============================================================================
  // Power-Sum And Schur Pairings
  // ============================================================================

  ring_elem powerSumsSchurInnerProductViaWeightedCharacters(
      ring_elem powerSums,
      ring_elem schurExpansion,
      const InnerProductContext& context) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
