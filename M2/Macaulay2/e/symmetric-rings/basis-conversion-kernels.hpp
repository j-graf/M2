// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_CONVERSION_KERNELS_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_CONVERSION_KERNELS_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  std::string jacobiTrudiCacheKey(int basisId,
                                    const Partition& outer,
                                    const Partition& inner) const;
  int popcountMask(size_t mask) const;
  int selectedGreaterThan(size_t mask, size_t col, size_t n) const;
  ring_elem jacobiTrudi(const Partition& outer,
                          const Partition& inner,
                          int basisId,
                          const std::string& display,
                          int order,
                          bool isMultiplicative) const;
  ring_elem scaled(ring_elem coeff, ring_elem f) const;
  ring_elem coefficientQuotient(ring_elem numerator, ring_elem denominator) const;
  ring_elem hallLittlewoodCFactor(const Partition& lambda) const;
  CoeffMap hallLittlewoodSingleCycleGreenMap(int n, bool normalized) const;
  ring_elem hallLittlewoodCapitalNormalizedConversionViaDiagonalScaling(
      ring_elem f,
      int sourceBasisId,
      int targetBasisId,
      const std::string& targetDisplay,
      int targetOrder,
      bool capitalToNormalized) const;
  ring_elem schurOmegaConversionViaPartitionConjugation(
      ring_elem f,
      int sourceBasisId,
      int targetBasisId,
      const std::string& targetDisplay,
      int targetOrder) const;
  void addCoeff(CoeffMap& target, const Partition& index, ring_elem coeff) const;
  void addScaledCoeffMap(CoeffMap& target,
                         ring_elem coeff,
                         const CoeffMap& source) const;
  CoeffMap addCoeffMaps(const CoeffMap& a, const CoeffMap& b) const;
  CoeffMap multiplyCoeffMaps(const CoeffMap& a, const CoeffMap& b) const;
  CoeffMap multiplyMonomialCoeffMaps(const CoeffMap& a, const CoeffMap& b) const;
  CoeffMap oneCoeffMap() const;
  bool isPartitionIndex(const Partition& p) const;
  int partitionPart(const Partition& p, size_t i) const;
  bool partitionContains(const Partition& outer, const Partition& inner) const;
  std::string powerSumsToSchurRecipeKey(int degree,
                                          const std::vector<Partition>& inputPartitions,
                                          bool omegaStyle) const;
  const std::vector<SchurConversionRecipeEntry>&
    powerSumsToSchurRecipe(int degree,
                           const std::vector<Partition>& inputPartitions,
                           bool omegaStyle) const;
  Partition leadingPartition(const CoeffMap& H) const;
  ring_elem completePartToPowerSumsViaClassicalFormula(int n) const;
  ring_elem elementaryPartToPowerSumsViaClassicalFormula(int n) const;
  ring_elem hallLittlewoodGeneratorPartToPowerSumsViaClassicalFormula(int n, bool omega) const;
  CoeffMap hallLittlewoodGeneratorPartToPowerSumsQuotientMapViaClassicalFormula(
      int n,
      bool omega) const;
  CoeffMap raisingExpansion(const Partition& lambda) const;
  CoeffMap raisingGeneratorMap(const Partition& lambda) const;
  ring_elem hallLittlewoodCapitalToPowerSumsViaRaisingOperators(const Partition& lambda, bool omega) const;
  ring_elem hallLittlewoodNormalizedToPowerSumsViaCapitalNormalization(const Partition& lambda, bool omega) const;
  ring_elem schurLikeToPowerSumsViaCharacters(const Partition& lambda,
                                               bool omegaStyle) const;
  CoeffMap powerSumPartToGeneratorMapViaLogarithmFormula(
      int n,
      ring_elem common) const;
  ring_elem powerSumLogarithmCoefficient(const Partition& lambda) const;
  CoeffMap powerSumPartToCompleteMapViaLogarithmFormula(int n) const;
  CoeffMap powerSumPartToElementaryMapViaLogarithmFormula(int n) const;
  CoeffMap powerSumIndexToCompleteMapViaLogarithmFormula(
      const Partition& index) const;
  CoeffMap powerSumIndexToElementaryMapViaLogarithmFormula(
      const Partition& index) const;
  CoeffMap powerSumPartToHallGeneratorMapViaLogarithmFormula(int n,
                                                             bool omega) const;
  CoeffMap powerSumIndexToHallGeneratorMapViaLogarithmFormula(
      const Partition& index,
      bool omega) const;
  CoeffMap powerSumsToHallGeneratorMapViaLogarithmFormula(ring_elem f,
                                                           bool omega) const;
  CoeffMap triangularReduceHallCapital(const CoeffMap& generatorMap,
                                         bool omega) const;
  ring_elem powerSumIndexToSchurLikeViaCharacters(const Partition& mu,
                                  int schurId,
                                  int schurOrder,
                                  const std::string& display,
                                  long sign) const;
  ring_elem powerSumIndexToSchurViaCharacters(const Partition& mu, int schurId, int schurOrder) const;
  ring_elem powerSumsToSchurLikeViaCharacters(ring_elem f,
                                   int schurId,
                                   int schurOrder,
                                   const std::string& display,
                                   bool omegaStyle) const;
  ring_elem powerSumsToSchurViaCharacters(ring_elem f,
                                           int schurId,
                                           int schurOrder) const;
  ring_elem monomialToPowerSumsViaTransitionMatrix(const Partition& lambda, bool forgotten) const;
  ring_elem powerSumIndexToMonomialViaTransitionMatrix(const Partition& lambda,
                                            int targetBasisId,
                                            const std::string& targetDisplay,
                                            int targetDisplayOrder,
                                            bool forgotten) const;
  ring_elem coeffMapToElement(const CoeffMap& H,
                                int targetBasisId,
                                const std::string& targetDisplay,
                                int targetDisplayOrder,
                                bool targetIsMultiplicative) const;
  Partition replaceAdjacentPair(const Partition& alpha,
                                  size_t pos,
                                  int first,
                                  int second) const;
  ring_elem straightenSchurAtom(const Partition& alpha,
                                  const std::string& display) const;
  ring_elem straightenHallCapitalAtom(const Partition& alpha,
                                        const std::string& display) const;
  ring_elem straightenAtom(const SymmetricMonomial& monomial, size_t pos) const;
  ring_elem straightenMonomial(const SymmetricMonomial& monomial) const;
  ring_elem straightenElement(ring_elem f) const;
  int requiredBasisIdForDisplay(const std::string& display) const;
  bool singleBasisIndexFromMonomial(const SymmetricMonomial& monomial,
                                      int basisId,
                                      Partition& index) const;
  ring_elem basisElementForDisplay(const std::string& display,
                                     const Partition& index) const;
  ring_elem powerSumsToHallLittlewoodCapitalViaTriangularReduction(ring_elem f,
                                         int targetBasisId,
                                         const std::string& targetDisplay,
                                         int targetDisplayOrder) const;
  ring_elem powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials(
      ring_elem f,
      int targetBasisId,
      const std::string& targetDisplay,
      int targetDisplayOrder) const;
  ring_elem powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality(
      ring_elem f,
      int targetBasisId,
      const std::string& targetDisplay,
      int targetDisplayOrder) const;
  CoeffMap schurGeneratorMap(const Partition& lambda,
                               bool omegaStyle,
                               int generatorId,
                               const std::string& generatorDisplay) const;
  CoeffMap triangularReduceSchur(const CoeffMap& generatorMap,
                                   bool omegaStyle,
                                   int generatorId,
                                   const std::string& generatorDisplay) const;
  bool tryExpressionToSchurViaTriangularReduction(ring_elem f,
                                         int targetBasisId,
                                         const std::string& targetDisplay,
                                         int targetDisplayOrder,
                                         ring_elem& result) const;
  bool tryExpressionToHallLittlewoodViaTriangularReduction(ring_elem f,
                                        int targetBasisId,
                                        const std::string& targetDisplay,
                                        int targetDisplayOrder,
                                        ring_elem& result) const;
  Partition atomIndex(const SymmetricMonomial& monomial, size_t pos) const;
  Partition atomOuterIndex(const SymmetricMonomial& monomial, size_t pos) const;
  Partition atomInnerIndex(const SymmetricMonomial& monomial, size_t pos) const;
  ring_elem atomToPowerSumsDispatch(const SymmetricMonomial& monomial, size_t pos) const;
  ring_elem monomialToPowerSumsDispatch(const SymmetricMonomial& monomial) const;
  ring_elem expressionToPowerSumsDispatch(ring_elem f) const;
  ring_elem skewQOrBFunction(const Partition& lambda,
                             const Partition& mu,
                             bool omega) const;
  ring_elem replaceSingleBasis(ring_elem f,
                               int sourceBasisId,
                               const std::string& targetDisplay) const;
  ring_elem skewPOrRToPowerSums(const Partition& lambda,
                                const Partition& mu,
                                bool omega) const;
  ring_elem skewHallLittlewoodToPowerSums(const Partition& lambda,
                                          const Partition& mu,
                                          const std::string& display) const;
  ring_elem powerSumElementFromIndex(const Partition& index) const;
  int singleBasisIdInMonomial(const SymmetricMonomial& monomial) const;
  int singleBasisId(ring_elem f) const;
  bool powerSumIndexFromMonomial(const SymmetricMonomial& monomial,
                                   Partition& index) const;
  ring_elem powerSumIndexToTargetDispatch(const Partition& index,
                                       const std::string& targetDisplay,
                                       int targetBasisId,
                                       int targetDisplayOrder,
                                       bool targetIsMultiplicative) const;
  ring_elem powerSumsToTargetViaTermwiseConversion(ring_elem f,
                                        int targetBasisId,
                                        const std::string& targetDisplay,
                                        int targetDisplayOrder,
                                        bool targetIsMultiplicative) const;
  ring_elem powerSumsToCompleteViaLogarithmFormula(ring_elem f,
                                                    int completeId,
                                                    int completeOrder) const;
  ring_elem powerSumsToElementaryViaLogarithmFormula(ring_elem f,
                                                      int elementaryId,
                                                      int elementaryOrder) const;
  ring_elem powerSumsToHallLittlewoodViaTriangularReduction(
      ring_elem f,
      int targetBasisId,
      const std::string& targetDisplay,
      int targetDisplayOrder) const;
  bool tryAtomToTarget(const SymmetricMonomial& monomial,
                            size_t pos,
                            int targetBasisId,
                            const std::string& targetDisplay,
                            int targetDisplayOrder,
                            bool targetIsMultiplicative,
                            ring_elem& result) const;
  bool tryMonomialToTarget(const SymmetricMonomial& monomial,
                                int targetBasisId,
                                const std::string& targetDisplay,
                                int targetDisplayOrder,
                                bool targetIsMultiplicative,
                                ring_elem& result) const;
  ExpressionToTargetMethod selectExpressionToTargetMethod(
                               const std::string& targetDisplay,
                               bool targetIsMultiplicative) const;
  bool executeExpressionToTargetMethod(
                               ExpressionToTargetMethod method,
                               ring_elem f,
                               int targetBasisId,
                               const std::string& targetDisplay,
                               int targetDisplayOrder,
                               bool targetIsMultiplicative,
                               ring_elem& result) const;
  const char *expressionToTargetMethodName(
                               ExpressionToTargetMethod method) const;
  void traceExpressionToTargetSelection(
                               ExpressionToTargetMethod method,
                               const std::string& targetDisplay) const;
  bool tryExpressionToTarget(ring_elem f,
                               int targetBasisId,
                               const std::string& targetDisplay,
                               int targetDisplayOrder,
                               bool targetIsMultiplicative,
                               ring_elem& result) const;

 public:
  CoeffMap multiplySchurExpansionViaRowPieri(const CoeffMap& source, int row) const;
  CoeffMap completeToSchurCoefficientsViaRecursiveTransition(const CoeffMap& hCoeffs) const;
  ring_elem completeToSchurViaRecursiveTransition(ring_elem f,
                                  int hBasisId,
                                  const std::string& hDisplay,
                                  int hOrder,
                                  bool hIsMultiplicative,
                                  int schurId,
                                  const std::string& schurDisplay,
                                  int schurOrder) const;
  ring_elem jacobiTrudiBasis(int basisId,
                               const std::string& display,
                               int order,
                               bool isMultiplicative,
                               const Partition& outer,
                               const Partition& inner) const;
  ring_elem straighten(ring_elem f) const;

 private:

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
