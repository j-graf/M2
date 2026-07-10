// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_BASIS_CONVERSION_DISPATCH_HPP_
#define M2_SYMMETRIC_RINGS_BASIS_CONVERSION_DISPATCH_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  enum class KnownState { Unknown, True, False };
  enum class ConversionPipeline
  {
    WholeExpression,
    PowerSum,
    PostPlethysmPowerSum,
    FallbackTerm,
    FactorizedProduct,
    PostPlethysm
  };
  enum class ConversionRequestKind
  {
    Expression,
    FactorizedProduct,
    PostPlethysm
  };
  enum class PowerSumsToSchurMethod
  {
    ViaBorderStrips,
    ViaComplete,
    ViaCharacters
  };
  enum class SourceToTargetMethod
  {
    AlreadyInTarget,
    ViaPowerSumKernels,
    ViaCompleteRecursiveTransition,
    ViaComplete,
    ViaPowerSums
  };
  enum class ProductExpansionMethod
  {
    ViaLittlewoodRichardson,
    ViaHorizontalPieri,
    ViaVerticalPieri,
    ViaBorderStrips,
    ViaSchurCompatibleRules,
    ViaFactorwiseConversion,
    NoApplicableMethod
  };

  struct ConversionGuarantees
  {
    std::optional<int> pureBasis;
    std::optional<int> expandedBasis;
    std::optional<int> homogeneousWeight;
    std::optional<size_t> termCount;
    std::optional<size_t> maximumPartitionLength;
    std::optional<std::vector<int>> factorBases;

    KnownState singleAtom = KnownState::Unknown;
    KnownState singleTerm = KnownState::Unknown;
    KnownState noProducts = KnownState::Unknown;
    KnownState normalized = KnownState::Unknown;
    KnownState skewFree = KnownState::Unknown;
    KnownState collected = KnownState::Unknown;
    KnownState targetClosed = KnownState::Unknown;
    KnownState factorizedProduct = KnownState::Unknown;
  };

  struct ConversionInput
  {
    ring_elem expression;
    ConversionGuarantees guarantees;
    SymmetricConversionOrigin origin = SymmetricConversionOrigin::Unknown;
  };

  struct ConversionRequest
  {
    ConversionRequestKind kind;
    ConversionInput input;
    std::optional<ring_elem> leftOperand;
    std::optional<ring_elem> rightOperand;
  };

  struct TermConversionClassification
  {
    bool isScalar;
    bool isSingleFactor;
    bool isProduct;
    bool targetIsMultiplicative;
    int factorCount;
    int singleBasisId;
    std::string targetDisplay;
    std::vector<std::string> factorDisplays;
    bool hasSkewFactor;
  };

  ConversionPipeline selectConversionPipeline(
        const ConversionRequest& request,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay,
        bool targetIsMultiplicative) const;
  ring_elem executeConversionPipeline(
        ConversionPipeline pipeline,
        const ConversionRequest& request,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  ring_elem runWholeExpressionPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  ring_elem runPowerSumPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  ring_elem runPostPlethysmPowerSumPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  ring_elem runFallbackTermPipeline(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  ring_elem runPostPlethysmPipeline(
        const ConversionInput& requestInput,
        ring_elem f,
        ring_elem g,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  ConversionGuarantees inferConversionGuarantees(
        ring_elem f,
        int targetBasisId) const;
  ConversionGuarantees strengthenConversionGuarantees(
        ConversionGuarantees guarantees,
        int targetBasisId) const;
  ConversionGuarantees guaranteesForFactorizedProduct(
        ring_elem f,
        ring_elem g,
        int targetBasisId) const;
  ConversionGuarantees guaranteesAfterNormalization(
        ring_elem f,
        int targetBasisId) const;
  ConversionGuarantees guaranteesAfterProductExpansion(
        ring_elem f,
        int targetBasisId) const;
  ConversionGuarantees guaranteesAfterSourceTargetConversion(
        const ConversionGuarantees& inputGuarantees,
        ring_elem result,
        int targetBasisId) const;
  ConversionGuarantees guaranteesAfterAddition(
        const ConversionGuarantees& inputGuarantees,
        ring_elem result,
        int targetBasisId) const;
  void attachConversionGuarantees(
        ring_elem f,
        const ConversionGuarantees& guarantees,
        SymmetricConversionOrigin origin = SymmetricConversionOrigin::Unknown) const;
  bool canUseWholeExpressionPipeline(
        const ConversionGuarantees& guarantees,
        int targetBasisId,
        const std::string& targetDisplay,
        bool targetIsMultiplicative) const;
  PowerSumsToSchurMethod selectPowerSumsToSchurMethod(
        const ConversionGuarantees& guarantees) const;
  SourceToTargetMethod selectSourceToTargetMethod(
        const ConversionInput& input,
        int pBasisId,
        int targetBasisId,
        const std::string& targetDisplay) const;
  ring_elem sourceToTargetDispatch(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  ring_elem powerSumsToSchurDispatch(
        const ConversionInput& input,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder) const;
  ring_elem powerSumsToSchurViaComplete(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder) const;
  const char *powerSumsToSchurMethodName(
        PowerSumsToSchurMethod method) const;
  const char *sourceToTargetMethodName(SourceToTargetMethod method) const;
  const char *productExpansionMethodName(ProductExpansionMethod method) const;
  const char *conversionPipelineName(ConversionPipeline pipeline) const;
  void traceConversionSelection(
        ConversionPipeline pipeline,
        const ConversionRequest& request,
        const std::string& targetDisplay) const;
  void traceSourceToTargetSelection(
        SourceToTargetMethod method,
        const ConversionInput& input,
        const std::string& targetDisplay) const;
  void traceProductExpansionSelection(
        ProductExpansionMethod method,
        const TermConversionClassification& classification) const;
  bool tryWholeExpressionToTarget(
        ring_elem f,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative,
        ring_elem& result) const;
  ring_elem runProductToBasisPipeline(
        ring_elem f,
        ring_elem g,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  bool normalizeTermForConversion(const SymmetricTerm& term,
                                    VECTOR(SymmetricTerm)& normalizedTerms) const;
  TermConversionClassification classifyNormalizedTerm(
        const SymmetricTerm& term,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder,
        bool targetIsMultiplicative) const;
  ProductExpansionMethod selectProductExpansionMethod(
        const TermConversionClassification& classification) const;
  bool executeProductExpansionMethod(
        ProductExpansionMethod method,
        const SymmetricTerm& term,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetDisplayOrder,
        bool targetIsMultiplicative,
        ConversionInput& result) const;
  ring_elem conversionRequestToBasisDispatch(
        const ConversionRequest& request,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;
  ring_elem toBasis(
        const ConversionInput& input,
        int pBasisId,
        const std::string& pDisplay,
        int pOrder,
        bool pIsMultiplicative,
        int targetBasisId,
        const std::string& targetDisplay,
        int targetOrder,
        bool targetIsMultiplicative) const;

 public:
  ring_elem toBasis(ring_elem f,
                      int pBasisId,
                      const std::string& pDisplay,
                      int pOrder,
                      bool pIsMultiplicative,
                      int targetBasisId,
                      const std::string& targetDisplay,
                      int targetOrder,
                      bool targetIsMultiplicative) const;
  ring_elem productToBasisDispatch(ring_elem f,
                                   ring_elem g,
                                   int pBasisId,
                                   const std::string& pDisplay,
                                   int pOrder,
                                   bool pIsMultiplicative,
                                   int targetBasisId,
                                   const std::string& targetDisplay,
                                   int targetOrder,
                                   bool targetIsMultiplicative) const;

 private:

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
