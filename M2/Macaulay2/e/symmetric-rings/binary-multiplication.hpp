// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BINARY_MULTIPLICATION_HPP_
#define M2_SYMMETRIC_RINGS_BINARY_MULTIPLICATION_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Strict Binary Multiplication
// ============================================================================

  CanonicalBasisTermView canonicalBasisTermView(
      ring_elem expression,
      const ExpressionFacts& facts) const;
  ring_elem executeBinaryMultiplicationKernel(
      const SelectedBinaryMultiplicationKernel& selection,
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right) const;
  ring_elem multiplyCanonicalBasisTerms(
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right,
      int targetBasisId) const;
  ring_elem multiplyCanonicalBasisTerms(
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right,
      int targetBasisId,
      const BinaryMultiplicationRequest& request) const;
  ring_elem multiplyBasisTermsViaMultiplicativeTarget(
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right,
      int targetBasisId) const;
  ring_elem multiplyBasisTermsViaPowerSums(
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right,
      int targetBasisId) const;

  void validateBinaryMultiplicationRequest(
      const BinaryMultiplicationRequest& request) const;
  void traceBinaryMultiplicationWorkflow(
      const char *workflow,
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right,
      int targetBasisId,
      const BinaryMultiplicationRequest& request) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
