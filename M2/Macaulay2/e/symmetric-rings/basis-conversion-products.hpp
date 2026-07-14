// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_CONVERSION_PRODUCTS_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_CONVERSION_PRODUCTS_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  struct SchurCompatibleFactor
  {
    enum Kind
    {
      General,
      Horizontal,
      Vertical,
      PowerSum,
      PowerSumAbacus,
      SchurExpansion
    };
    Kind kind;
    Partition index;
    CoeffMap expansion;
    int weight;
  };

  enum class SchurFactorMethod
  {
    AlreadySchur,
    ViaLittlewoodRichardson,
    ViaHorizontalPieri,
    ViaVerticalPieri,
    ViaBorderStrips,
    ViaAbacusRimHooks,
    ViaLittlewoodRichardsonExpansion
  };

  enum class ProductToTargetRoute
  {
    ViaSchurCompatibleFactors,
    ViaMonomialLikeExpansion,
    ViaHallLittlewoodGenerators,
    ViaConvertRightFactor,
    ViaConvertLeftFactor,
    AlreadyInTarget,
    NoApplicableRoute
  };

// ============================================================================
// Littlewood-Richardson And Skew Schur Rules
// ============================================================================
// Tableau and coefficient enumeration for Schur products and skew expansion.

  std::string littlewoodRichardsonProductCacheKey(
      const Partition& lambda,
      const Partition& mu) const;

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
  const std::vector<PartitionCoefficientTerm>& littlewoodRichardsonProductViaCoefficientEnumeration(const Partition& a,
                                                  const Partition& b) const;
  const std::vector<PartitionCoefficientTerm>& littlewoodRichardsonProductViaTableauEnumeration(const Partition& a,
                                                       const Partition& b) const;
  const std::vector<PartitionCoefficientTerm>& skewSchurToSchurViaLittlewoodRichardson(const Partition& outer,
                                                          const Partition& inner) const;

// ============================================================================
// Pieri Rules
// ============================================================================
// Horizontal and vertical strips implement multiplication by h_n and e_n.

  std::vector<Partition> schurTimesCompleteViaHorizontalPieri(const Partition& lambda,
                                                    int row) const;
  std::vector<Partition> schurTimesElementaryViaVerticalPieri(const Partition& lambda,
                                                  int col) const;

// ============================================================================
// Border Strips And Murnaghan-Nakayama
// ============================================================================
// Border-strip validation and Schur multiplication by power sums.

  bool addedBorderStripCell(const Partition& lambda,
                                const Partition& nu,
                                int row,
                                int col) const;
  bool addedBorderStripConnected(const Partition& lambda,
                                     const Partition& nu) const;
  bool addedBorderStripHasNoTwoByTwo(const Partition& lambda,
                                         const Partition& nu) const;
  const std::vector<PartitionCoefficientTerm>& schurTimesPowerSumViaBorderStrips(const Partition& lambda,
                                                            int part) const;
  ring_elem powerSumsToSchurViaBorderStrips(
          ring_elem f,
          int targetBasisId,
          const std::string& targetDisplay,
          int targetDisplayOrder) const;
  const std::vector<PartitionCoefficientTerm>& schurTimesPowerSumViaAbacusRimHooks(
          const Partition& lambda,
          int part) const;
  void addPowerSumIndexToSchurMapViaAbacusRimHooks(
          const Partition& index,
          ring_elem coefficient,
          CoeffMap& result) const;
  ring_elem powerSumsToSchurViaAbacusRimHooks(
          ring_elem f,
          int targetBasisId,
          const std::string& targetDisplay,
          int targetDisplayOrder) const;

// ============================================================================
// Schur Product Planning And Execution
// ============================================================================
// Factor classification selects LR, Pieri, border-strip, or converted-factor methods.

  ring_elem multiplySchurExpansionsViaLittlewoodRichardson(ring_elem f,
                                      ring_elem g,
                                      int schurId,
                                      const std::string& schurDisplay,
                                      int schurOrder) const;
  bool tryBasisElementToSchurFactors(const SymmetricMonomial& monomial,
                              size_t pos,
                              int targetBasisId,
                              std::vector<Partition>& factors) const;
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
  bool trySchurCompatibleExpressionToSchur(ring_elem f,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetDisplayOrder,
                              ring_elem& result) const;
  bool tryProductToSchurViaCompatibleFactors(ring_elem f,
                              ring_elem g,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetDisplayOrder,
                              ring_elem& result) const;

// ============================================================================
// Monomial And Forgotten Products
// ============================================================================
// Exponent splittings and retained-product routes for m and ff.

  long monomialProductCoefficientViaExponentSplittings(const Partition& lambda,
                                      const Partition& mu,
                                      const Partition& nu) const;
  const std::vector<PartitionCoefficientTerm>& monomialProductViaExponentSplittings(const Partition& a,
                                                       const Partition& b) const;
  bool tryMonomialLikeBasisElementToCoeffMap(const SymmetricMonomial& monomial,
                                      size_t pos,
                                      int targetBasisId,
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

// ============================================================================
// Hall-Littlewood Products
// ============================================================================
// Hall-Littlewood products are evaluated through multiplicative generator bases.

  bool tryProductToHallLittlewoodViaGenerators(
                                     ring_elem f,
                                     ring_elem g,
                                     int targetBasisId,
                                     const std::string& targetDisplay,
                                     int targetDisplayOrder,
                                     ring_elem& result) const;

// ============================================================================
// Product-To-Target Dispatch
// ============================================================================
// The product dispatcher chooses and executes one visible retained-operand route.

  ProductToTargetRoute selectProductToTargetRoute(
                                     ring_elem f,
                                     ring_elem g,
                                     int targetBasisId,
                                     const std::string& targetDisplay,
                                     bool targetIsMultiplicative) const;
  const char *productToTargetRouteName(ProductToTargetRoute route) const;
  void traceProductToTargetSelection(
                                     ProductToTargetRoute route,
                                     const std::string& targetDisplay) const;
  bool executeProductToTargetRoute(
                                     ProductToTargetRoute route,
                                     ring_elem f,
                                     ring_elem g,
                                     int targetBasisId,
                                     const std::string& targetDisplay,
                                     int targetDisplayOrder,
                                     bool targetIsMultiplicative,
                                     ring_elem& result) const;
  bool tryProductToTarget(ring_elem f,
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
