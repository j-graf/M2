// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_CONVERSION_PRODUCTS_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_CONVERSION_PRODUCTS_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  std::string littlewoodRichardsonProductCacheKey(const Partition& lambda, const Partition& mu) const;
  ring_elem cachedInteger(long n) const;
  long littlewoodRichardsonCoefficientViaTableaux(const Partition& lambda,
                       const Partition& content,
                       const Partition& nu) const;
  void partitionsContainingRec(const Partition& lambda,
                                 int addedWeight,
                                 size_t row,
                                 int previousPart,
                                 Partition& current,
                                 std::vector<Partition>& result) const;
  std::vector<Partition> partitionsContaining(const Partition& lambda,
                                                int addedWeight) const;
  const std::vector<LRProductTerm>& littlewoodRichardsonProductViaCoefficientEnumeration(const Partition& a,
                                                const Partition& b) const;
  const std::vector<LRProductTerm>& littlewoodRichardsonProductViaTableauEnumeration(const Partition& a,
                                                     const Partition& b) const;
  const std::vector<LRProductTerm>& skewSchurToSchurViaLittlewoodRichardson(const Partition& outer,
                                                        const Partition& inner) const;
  bool addedBorderStripCell(const Partition& lambda,
                              const Partition& nu,
                              int row,
                              int col) const;
  bool addedBorderStripConnected(const Partition& lambda,
                                   const Partition& nu) const;
  bool addedBorderStripHasNoTwoByTwo(const Partition& lambda,
                                       const Partition& nu) const;
  const std::vector<LRProductTerm>& schurTimesPowerSumViaBorderStrips(const Partition& lambda,
                                                          int part) const;
  long monomialProductCoefficientViaExponentSplittings(const Partition& lambda,
                                    const Partition& mu,
                                    const Partition& nu) const;
  const std::vector<LRProductTerm>& monomialProductViaExponentSplittings(const Partition& a,
                                                     const Partition& b) const;
  std::vector<Partition> schurTimesCompleteViaHorizontalPieri(const Partition& lambda,
                                                  int row) const;
  std::vector<Partition> schurTimesElementaryViaVerticalPieri(const Partition& lambda,
                                                int col) const;
  ring_elem multiplySchurExpansionsViaLittlewoodRichardson(ring_elem f,
                                    ring_elem g,
                                    int schurId,
                                    const std::string& schurDisplay,
                                    int schurOrder) const;
  bool tryAtomToSchurFactors(const SymmetricMonomial& monomial,
                            size_t pos,
                            int targetBasisId,
                            std::vector<Partition>& factors) const;
  bool trySchurCompatibleFactorsFromMonomial(
        const SymmetricMonomial& monomial,
        int targetBasisId,
        std::vector<SchurCompatibleFactor>& factors) const;
  SchurFactorMethod selectSchurFactorMethod(
        const SchurCompatibleFactor& factor) const;
  const char *schurFactorMethodName(SchurFactorMethod method) const;
  void traceSchurFactorMethod(
        SchurFactorMethod method,
        const SchurCompatibleFactor& factor) const;
  bool schurCompatibleFactorsToSchurDispatch(
        std::vector<SchurCompatibleFactor> factors,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder,
        ring_elem& result) const;
  bool trySchurProductMonomialToSchurViaLittlewoodRichardson(const SymmetricMonomial& monomial,
                                     int targetBasisId,
                                     const std::string& targetDisplay,
                                     int targetDisplayOrder,
                                     ring_elem& result) const;
  bool trySchurCompatibleMonomialToSchur(const SymmetricMonomial& monomial,
                                         int targetBasisId,
                                         const std::string& targetDisplay,
                                         int targetDisplayOrder,
                                         ring_elem& result) const;
  bool trySchurCompatibleExpressionToSchur(ring_elem f,
                            int targetBasisId,
                            const std::string& targetDisplay,
                            int targetDisplayOrder,
                            ring_elem& result) const;
  ring_elem powerSumsToSchurViaBorderStrips(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder) const;
  bool tryProductToSchurViaCompatibleFactors(ring_elem f,
                            ring_elem g,
                            int targetBasisId,
                            const std::string& targetDisplay,
                            int targetDisplayOrder,
                            ring_elem& result) const;
  bool tryMonomialLikeAtomToCoeffMap(const SymmetricMonomial& monomial,
                                    size_t pos,
                                    const std::string& targetDisplay,
                                    CoeffMap& result) const;
  bool tryMonomialLikeMonomialToTarget(const SymmetricMonomial& monomial,
                                      int targetBasisId,
                                      const std::string& targetDisplay,
                                      int targetDisplayOrder,
                                      bool targetIsMultiplicative,
                                      ring_elem& result) const;
  bool tryProductToMonomialLikeTarget(ring_elem f,
                                   ring_elem g,
                                   int targetBasisId,
                                   const std::string& targetDisplay,
                                   int targetDisplayOrder,
                                   bool targetIsMultiplicative,
                                   ring_elem& result) const;
  bool tryProductToHallLittlewoodViaGenerators(
                                   ring_elem f,
                                   ring_elem g,
                                   int targetBasisId,
                                   const std::string& targetDisplay,
                                   int targetDisplayOrder,
                                   ring_elem& result) const;
  ProductToTargetRoute selectProductToTargetRoute(
                                   ring_elem f,
                                   ring_elem g,
                                   int targetBasisId,
                                   const std::string& targetDisplay,
                                   bool targetIsMultiplicative) const;
  bool executeProductToTargetRoute(
                                   ProductToTargetRoute route,
                                   ring_elem f,
                                   ring_elem g,
                                   int targetBasisId,
                                   const std::string& targetDisplay,
                                   int targetDisplayOrder,
                                   bool targetIsMultiplicative,
                                   ring_elem& result) const;
  const char *productToTargetRouteName(ProductToTargetRoute route) const;
  void traceProductToTargetSelection(
                                   ProductToTargetRoute route,
                                   const std::string& targetDisplay) const;
  bool tryProductToTargetDispatch(ring_elem f,
                             ring_elem g,
                             int targetBasisId,
                             const std::string& targetDisplay,
                             int targetDisplayOrder,
                             bool targetIsMultiplicative,
                             ring_elem& result) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
