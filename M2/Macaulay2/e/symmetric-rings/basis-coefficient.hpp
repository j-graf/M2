// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_COEFFICIENT_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_COEFFICIENT_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Basis-Coefficient Route Contract
// ============================================================================

  enum class BasisCoefficientRoute
  {
    AlreadyExpandedInTarget,
    ViaPowerSumLookup,
    ViaSchurCharacters,
    ViaSchurOmegaCharacters,
    ViaCompleteLogarithmFormula,
    ViaElementaryLogarithmFormula,
    ViaHallLittlewoodGeneratorLogarithmFormula,
    ViaHallLittlewoodCapitalGreenPolynomialDuality,
    ViaMonomialTransition,
    ViaForgottenTransition,
    ViaFullBasisConversion
  };

  BasisCoefficientRoute selectBasisCoefficientRoute(
      ring_elem f,
      int targetBasisId) const;
  const char *basisCoefficientRouteName(BasisCoefficientRoute route) const;
  void traceBasisCoefficientSelection(
      BasisCoefficientRoute route,
      const std::string& targetDisplay,
      const Partition& targetIndex) const;
  ring_elem executeBasisCoefficientRoute(
      BasisCoefficientRoute route,
      ring_elem f,
      int targetBasisId,
      const Partition& targetIndex) const;
  ring_elem basisCoefficientDispatch(
      ring_elem f,
      int targetBasisId,
      const Partition& targetIndex) const;

// ============================================================================
// Public Basis-Coefficient Entry Point
// ============================================================================

 public:
  ring_elem basisCoefficient(
      ring_elem f,
      ring_elem targetBasisElement) const;

 private:

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
