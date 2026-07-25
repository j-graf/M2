// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_MULTIPLICATION_FOLDS_HPP_
#define M2_SYMMETRIC_RINGS_MULTIPLICATION_FOLDS_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Target-Closed Multiplication Families
// ============================================================================
// A family (w, A) states that an accumulator in w can be multiplied directly
// by a factor from every basis u in A and remain canonically in w. The
// multifactor workflow uses this declared mathematical closure; it never
// infers a family from the binary-kernel inventory.

  struct TargetClosedMultiplicationFamily
  {
    std::string identifier;
    BasisKind targetBasis = BasisKind::Custom;
    std::vector<BasisKind> factorBases;
  };

  static const std::vector<TargetClosedMultiplicationFamily>&
  targetClosedMultiplicationFamilyDatabase();
  static void validateTargetClosedMultiplicationFamilyDatabase();

  const TargetClosedMultiplicationFamily *
  selectTargetClosedMultiplicationFamily(
      const std::vector<int>& factorBasisIds,
      int targetBasisId) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
