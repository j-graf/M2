// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_NORMALIZATION_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_NORMALIZATION_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Basis-Normalization Formulas
// ============================================================================
// A normalization formula acts on one stored basis element. Straightening
// formulas return a canonical expansion in the same basis. Skew-expansion
// formulas return a normalized, skew-free expression, with products allowed.

 private:
  using BasisNormalizationFormula = ring_elem (
      SymmetricEngineRing::*)(
          const SymmetricMonomial&, size_t) const;

  struct BasisNormalizationRule
  {
    BasisKind basisKind;
    BasisNormalizationFormula straighteningFormula = nullptr;
    BasisNormalizationFormula skewExpansionFormula = nullptr;
  };

  static const std::vector<BasisNormalizationRule>&
  basisNormalizationRuleDatabase();
  static const std::map<BasisKind, const BasisNormalizationRule *>&
  basisNormalizationRulesByKind();
  static void validateBasisNormalizationRuleDatabase();
  const BasisNormalizationRule *basisNormalizationRuleFor(
      BasisKind basisKind) const;

  ring_elem straightenSchurFamilyBasisElement(
      const SymmetricMonomial& monomial,
      size_t position) const;
  ring_elem straightenHallLittlewoodCapitalBasisElement(
      const SymmetricMonomial& monomial,
      size_t position) const;
  ring_elem straightenHallLittlewoodNormalizedBasisElement(
      const SymmetricMonomial& monomial,
      size_t position) const;
  ring_elem expandSkewSchurFamilyBasisElement(
      const SymmetricMonomial& monomial,
      size_t position) const;
  ring_elem expandSkewHallLittlewoodBasisElement(
      const SymmetricMonomial& monomial,
      size_t position) const;

// ============================================================================
// Generic Normalization Workflow
// ============================================================================

  ring_elem straightenBasisElement(
      const SymmetricMonomial& monomial,
      size_t position) const;
  ring_elem straightenMonomial(
      const SymmetricMonomial& monomial) const;
  ring_elem straightenElement(ring_elem expression) const;
  ring_elem expandSkewBasisElement(
      const SymmetricMonomial& monomial,
      size_t position) const;
  ring_elem normalizeExpression(
      ring_elem expression,
      ExpressionFacts& resultFacts,
      std::vector<size_t> *termFactorCounts = nullptr) const;

 public:
  ring_elem straighten(ring_elem expression) const;

 private:

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
