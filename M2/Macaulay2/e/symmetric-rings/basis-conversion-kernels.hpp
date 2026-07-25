// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_CONVERSION_KERNELS_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_CONVERSION_KERNELS_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Conversion-Kernel Interface
// ============================================================================
// Every plan-callable kernel receives the same mathematical input. Lower
// helpers remain free to use signatures natural to their formulas.

  struct BasisConversionInput
  {
    ring_elem expansion;
    RingBasis source;
    RingBasis target;
    std::optional<int> homogeneousWeight;
  };

  using BasisConversionKernel =
      ring_elem (SymmetricEngineRing::*)(
          const BasisConversionInput&) const;

// ============================================================================
// Plan-Callable Conversion Kernels
// ============================================================================
// These kernels are the only formula entry points named by conversion plans.
// Their names state the mathematical source, target, and formula; family
// helpers later in this file contain no plan selection.

  // Classical and Hall--Littlewood bases to power sums.
  ring_elem completeToPowerSumsViaNewtonIdentities(
      const BasisConversionInput& input) const;
  ring_elem elementaryToPowerSumsViaNewtonIdentities(
      const BasisConversionInput& input) const;
  ring_elem schurToPowerSumsViaFrobeniusCharacterFormula(
      const BasisConversionInput& input) const;
  ring_elem schurOmegaToPowerSumsViaFrobeniusCharacterFormula(
      const BasisConversionInput& input) const;
  ring_elem monomialToPowerSumsViaTransitionMatrix(
      const BasisConversionInput& input) const;
  ring_elem forgottenToPowerSumsViaTransitionMatrix(
      const BasisConversionInput& input) const;
  ring_elem hallLittlewoodQGeneratorsToPowerSumsViaGeneratingFunction(
      const BasisConversionInput& input) const;
  ring_elem hallLittlewoodBGeneratorsToPowerSumsViaGeneratingFunction(
      const BasisConversionInput& input) const;
  ring_elem hallLittlewoodQToPowerSumsViaRaisingOperators(
      const BasisConversionInput& input) const;
  ring_elem hallLittlewoodBToPowerSumsViaRaisingOperators(
      const BasisConversionInput& input) const;
  ring_elem hallLittlewoodPToPowerSumsViaQNormalization(
      const BasisConversionInput& input) const;
  ring_elem hallLittlewoodPOmegaToPowerSumsViaBNormalization(
      const BasisConversionInput& input) const;

  // Power sums to classical bases.
  ring_elem powerSumsToCompleteViaLogarithmFormula(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToElementaryViaLogarithmFormula(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToSchurViaMurnaghanNakayama(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToSchurViaAbacusRimHooks(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToSchurViaFrobeniusCharacterFormula(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToSchurOmegaViaMurnaghanNakayama(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToSchurOmegaViaAbacusRimHooks(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToSchurOmegaViaFrobeniusCharacterFormula(
      const BasisConversionInput& input) const;

  // Power sums to Hall--Littlewood, monomial, and forgotten bases.
  ring_elem powerSumsToHallLittlewoodGeneratorsViaLogarithmFormula(
      const BasisConversionInput& input) const;
  ring_elem powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials(
      const BasisConversionInput& input) const;
  ring_elem powerSumBasisElementToHallLittlewoodViaGreenPolynomialsAndDuality(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToHallLittlewoodViaTriangularReduction(
      const BasisConversionInput& input) const;
  ring_elem powerSumsToMonomialOrForgottenViaTransitionMatrix(
      const BasisConversionInput& input) const;

  // Direct involutions, normalizations, and triangular transitions.
  ring_elem hallLittlewoodPairedBasesViaDiagonalScaling(
      const BasisConversionInput& input) const;
  ring_elem schurAndSchurOmegaViaPartitionConjugation(
      const BasisConversionInput& input) const;
  ring_elem schurAndSchurOmegaToGeneratorsViaJacobiTrudi(
      const BasisConversionInput& input) const;
  ring_elem completeToSchurViaHorizontalPieri(
      const BasisConversionInput& input) const;
  ring_elem hallLittlewoodGeneratorsToCapitalBasesViaTriangularReduction(
      const BasisConversionInput& input) const;

// ============================================================================
// Shared Conversion Utilities
// ============================================================================
// Coefficient maps, basis metadata, stored basis-element inspection, and
// linear extension of basis-element formulas.

  ring_elem scaled(ring_elem coeff, ring_elem f) const;
  ring_elem coefficientQuotient(ring_elem numerator, ring_elem denominator) const;
  void addNormalizedCoeff(CoeffMap& target,
                          Partition index,
                          ring_elem coeff) const;
  void addCoeff(CoeffMap& target, const Partition& index, ring_elem coeff) const;
  void addScaledCoeffMap(CoeffMap& target,
                           ring_elem coeff,
                           const CoeffMap& source) const;
  CoeffMap addCoeffMaps(const CoeffMap& a, const CoeffMap& b) const;
  CoeffMap multiplyCoeffMaps(const CoeffMap& a, const CoeffMap& b) const;
  CoeffMap oneCoeffMap() const;
  ring_elem cachedInteger(long n) const;
  ring_elem cachedInteger(const mpz_class& n) const;
  Partition leadingPartition(const CoeffMap& H) const;
  ring_elem coeffMapToElement(const CoeffMap& H,
                                  int targetBasisId,
                                  const std::string& targetDisplay,
                                  int targetDisplayOrder,
                                  bool targetIsMultiplicative) const;
  ring_elem basisElementForKind(BasisKind kind,
                                       const Partition& index) const;
  ring_elem replaceSingleBasis(ring_elem f,
                                 int sourceBasisId,
                                 int targetBasisId) const;
  Partition basisElementIndex(
        const SymmetricMonomial& monomial,
        size_t pos) const;
  Partition basisElementOuterIndex(
        const SymmetricMonomial& monomial,
        size_t pos) const;
  Partition basisElementInnerIndex(
        const SymmetricMonomial& monomial,
        size_t pos) const;
  int singleBasisIdInMonomial(const SymmetricMonomial& monomial) const;
  int singleBasisId(ring_elem f) const;
  bool powerSumIndexFromMonomial(const SymmetricMonomial& monomial,
                                     Partition& index) const;
  ring_elem powerSumElementFromIndex(const Partition& index) const;
  ring_elem linearlyExtendToPowerSums(
        ring_elem expansion,
        BasisKind sourceKind,
        const std::function<ring_elem(const Partition&)>&
            basisElementFormula) const;
  bool singleBasisIndexFromMonomial(const SymmetricMonomial& monomial,
                                        int basisId,
                                        Partition& index) const;

// ============================================================================
// Jacobi-Trudi And Determinant Utilities
// ============================================================================
// Determinantal constructions shared by Schur-style conversion and straightening.

  std::string jacobiTrudiCacheKey(int basisId,
                                      const Partition& outer,
                                      const Partition& inner) const;
  int popcountMask(size_t mask) const;
  int selectedGreaterThan(size_t mask, size_t col, size_t n) const;
  ring_elem jacobiTrudi(const Partition& outer,
                        const Partition& inner,
                        int basisId) const;
  ring_elem canonicalSchurLikeExpressionToGeneratorsViaJacobiTrudi(
      ring_elem f,
      int sourceBasisId,
      int targetBasisId) const;
 public:
  ring_elem jacobiTrudiBasis(int basisId,
                                 const Partition& outer,
                                 const Partition& inner) const;

// ============================================================================
// Complete And Elementary Bases
// ============================================================================
// Classical h/e-to-p formulas, logarithmic inverse formulas, and h-to-S transition.

 private:
  ring_elem completePartToPowerSumsViaNewtonIdentities(int n) const;
  ring_elem elementaryPartToPowerSumsViaNewtonIdentities(int n) const;
  ring_elem powerSumLogarithmCoefficient(const Partition& lambda) const;
  CoeffMap powerSumPartToIntegralGeneratorMapViaLogarithmFormula(
        int n,
        int commonSign) const;
  CoeffMap powerSumPartToGeneratorMapViaLogarithmFormula(
        int n,
        ring_elem common) const;
  const CoeffMap& powerSumPartToCompleteMapViaLogarithmFormula(int n) const;
  const CoeffMap& powerSumPartToElementaryMapViaLogarithmFormula(int n) const;
  CoeffMap powerSumIndexToCompleteMapViaLogarithmFormula(
        const Partition& index) const;
  CoeffMap powerSumIndexToElementaryMapViaLogarithmFormula(
        const Partition& index) const;
  CoeffMap powerSumsToCompleteMapViaLogarithmFormula(ring_elem f) const;
  void addPowerSumIndexToCompleteMapViaLogarithmFormula(
        const Partition& index,
        ring_elem coefficient,
        CoeffMap& result) const;
  ring_elem powerSumsToCompleteViaLogarithmFormula(ring_elem f,
                                                      int completeId,
                                                      int completeOrder) const;
  ring_elem powerSumsToElementaryViaLogarithmFormula(ring_elem f,
                                                        int elementaryId,
                                                        int elementaryOrder) const;
 public:
  CoeffMap multiplySchurExpansionViaRowPieri(const CoeffMap& source, int row) const;
  CoeffMap completeToSchurCoefficientsViaHorizontalPieri(const CoeffMap& hCoeffs) const;
  ring_elem completeToSchurViaHorizontalPieri(ring_elem f,
                                    int hBasisId,
                                    const std::string& hDisplay,
                                    int hOrder,
                                    bool hIsMultiplicative,
                                    int schurId,
                                    const std::string& schurDisplay,
                                    int schurOrder) const;

// ============================================================================
// Schur And Omega-Schur Bases
// ============================================================================
// Character, omega, and triangular conversions for S and Somega.

 private:
  ring_elem schurOmegaConversionViaPartitionConjugation(
        ring_elem f,
        int sourceBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder) const;
  ring_elem schurLikeToPowerSumsViaFrobeniusCharacterFormula(
      const Partition& lambda,
      bool omegaStyle) const;
  ring_elem powerSumIndexToSchurLikeViaFrobeniusCharacterFormula(
      const Partition& mu,
      int schurId,
      int schurOrder,
      const std::string& display,
      long sign) const;
  ring_elem powerSumIndexToSchurViaFrobeniusCharacterFormula(
      const Partition& mu,
      int schurId,
      int schurOrder) const;
  ring_elem powerSumsToSchurLikeViaFrobeniusCharacterFormula(
      ring_elem f,
      int schurId,
      int schurOrder,
      const std::string& display,
      bool omegaStyle) const;
  ring_elem powerSumsToSchurViaFrobeniusCharacterFormula(
      ring_elem f,
      int schurId,
      int schurOrder) const;
  std::string powerSumsToSchurRecipeKey(int degree,
                                            const std::vector<Partition>& inputPartitions,
                                            bool omegaStyle) const;
  const std::vector<SchurConversionRecipeEntry>&
      powerSumsToSchurRecipe(int degree,
                             const std::vector<Partition>& inputPartitions,
                             bool omegaStyle) const;

// ============================================================================
// Monomial And Forgotten Bases
// ============================================================================
// Transition matrices and coefficient-map products for m and ff.

  CoeffMap multiplyMonomialCoeffMaps(const CoeffMap& a, const CoeffMap& b) const;
  ring_elem monomialToPowerSumsViaTransitionMatrix(const Partition& lambda, bool forgotten) const;
  ring_elem powerSumIndexToMonomialViaTransitionMatrix(const Partition& lambda,
                                              int targetBasisId,
                                              const std::string& targetDisplay,
                                              int targetDisplayOrder,
                                              bool forgotten) const;

// ============================================================================
// Hall-Littlewood Generator Bases
// ============================================================================
// Generating-function and logarithmic conversions for the multiplicative
// q and b bases.

  ring_elem hallLittlewoodGeneratorPartToPowerSumsViaGeneratingFunction(
      int n,
      bool omega) const;
  CoeffMap hallLittlewoodGeneratorPartToPowerSumsQuotientMapViaGeneratingFunction(
      int n,
      bool omega) const;
  CoeffMap powerSumPartToHallGeneratorMapViaLogarithmFormula(int n,
                                                               bool omega) const;
  CoeffMap powerSumIndexToHallGeneratorMapViaLogarithmFormula(
        const Partition& index,
        bool omega) const;
  CoeffMap powerSumsToHallGeneratorMapViaLogarithmFormula(ring_elem f,
                                                             bool omega) const;

// ============================================================================
// Hall-Littlewood Capital And Normalized Bases
// ============================================================================
// Raising operators, Green polynomials, normalization, skew functions, and triangular reduction.

  ring_elem hallLittlewoodCFactor(const Partition& lambda) const;
  CoeffMap hallLittlewoodSingleCycleGreenMap(int n, bool normalized) const;
  ring_elem hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
        ring_elem f,
        int sourceBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool capitalToNormalized) const;
  CoeffMap raisingExpansion(const Partition& lambda) const;
  CoeffMap raisingGeneratorMap(const Partition& lambda) const;
  ring_elem hallLittlewoodCapitalToPowerSumsViaRaisingOperators(
      const Partition& lambda,
      bool omega) const;
  ring_elem hallLittlewoodNormalizedToPowerSumsViaPairedBasisNormalization(
      const Partition& lambda,
      bool omega) const;
  CoeffMap triangularReduceHallCapital(const CoeffMap& generatorMap,
                                           bool omega) const;
  ring_elem skewQOrBFunction(const Partition& lambda,
                               const Partition& mu,
                               bool omega) const;
  ring_elem skewPOrPOmegaToPowerSums(const Partition& lambda,
                                  const Partition& mu,
                                  bool omega) const;
  ring_elem skewHallLittlewoodToPowerSums(const Partition& lambda,
                                            const Partition& mu,
                                            BasisKind basisKind) const;
  ring_elem powerSumsToHallLittlewoodCapitalViaTriangularReduction(ring_elem f,
                                           int targetBasisId,
                                           const std::string& targetDisplay,
                                           int targetDisplayOrder) const;
  ring_elem powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder) const;
  ring_elem powerSumIndexToHallLittlewoodCapitalCoefficientViaGreenPolynomialDuality(
        const Partition& cycleType,
        const Partition& targetIndex) const;
  ring_elem powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder) const;
  ring_elem powerSumsToHallLittlewoodViaTriangularReduction(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder) const;
  ring_elem hallLittlewoodGeneratorsToCapitalViaTriangularReduction(
      ring_elem f,
      int targetBasisId,
      const std::string& targetDisplay,
      int targetDisplayOrder) const;

// Shared recurrence formulas used by the declarative normalization rules.
  Partition replaceAdjacentPair(const Partition& alpha,
                                    size_t pos,
                                    int first,
                                    int second) const;
  ring_elem straightenSchurBasisElement(const Partition& alpha,
                                    int basisId) const;
  ring_elem straightenHallCapitalBasisElement(const Partition& alpha,
                                          int basisId) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
