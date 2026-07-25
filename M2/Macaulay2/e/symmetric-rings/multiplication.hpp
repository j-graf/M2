// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_MULTIPLICATION_HPP_
#define M2_SYMMETRIC_RINGS_MULTIPLICATION_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Linear Extension And Multifactor Resolution
// ============================================================================

  struct PreparedMultiplicationTerm
  {
    const SymmetricTerm *term = nullptr;
    bool scalar = false;
    ring_elem factor;
    ExpressionFacts facts;
  };

  struct PreparedMultiplicationOperands
  {
    ring_elem left;
    ring_elem right;
    ExpressionFacts leftFacts;
    ExpressionFacts rightFacts;
    std::vector<PreparedMultiplicationTerm> leftTerms;
    std::vector<PreparedMultiplicationTerm> rightTerms;
    CombinatorialTags combinatorialTags = 0;
  };

  struct ProductFactor
  {
    ring_elem expression;
    ExpressionFacts facts;
  };

  struct ResolvedProductTerm
  {
    ring_elem expression;
    ExpressionFacts facts;
  };

  enum class BilinearMultiplicationStrategy
  {
    Automatic,
    MultiplicativeTarget,
    KernelDistribution,
    PowerSumFallback
  };

  ExpressionFacts finalizeCanonicalMultiplicationResult(
      ring_elem result,
      int targetBasisId,
      CombinatorialTags combinatorialTags) const;
  ring_elem multiplyCanonicalExpansionsBalanced(
      std::vector<ring_elem> factors) const;

  PreparedMultiplicationOperands prepareMultiplicationOperands(
      ring_elem left,
      ring_elem right) const;
  ring_elem convertPreparedMultiplicationTerm(
      const PreparedMultiplicationTerm& term,
      int targetBasisId,
      CombinatorialTags tags) const;
  ring_elem convertCanonicalFactorToBasis(
      ring_elem factor,
      const ExpressionFacts& facts,
      int targetBasisId,
      CombinatorialTags tags) const;
  ring_elem multiplyExpressionsViaMultiplicativeTarget(
      const PreparedMultiplicationOperands& operands,
      int targetBasisId) const;
  ring_elem multiplyExpressionsViaPowerSums(
      const PreparedMultiplicationOperands& operands,
      int targetBasisId) const;
  bool everyNonscalarPairHasAutomaticKernel(
      const PreparedMultiplicationOperands& operands,
      int targetBasisId) const;
  ring_elem distributeMultiplicationOverTermPairs(
      const PreparedMultiplicationOperands& operands,
      int targetBasisId) const;
  ring_elem multiplyToBasisWithStrategy(
      ring_elem left,
      ring_elem right,
      int targetBasisId,
      BilinearMultiplicationStrategy strategy) const;

  std::vector<ProductFactor> productFactors(
      const SymmetricTerm& term) const;
  std::vector<int> basisIdsOfProductFactors(
      const std::vector<ProductFactor>& factors) const;
  ring_elem multiplyTermViaMultiplicativeTarget(
      const std::vector<ProductFactor>& factors,
      int targetBasisId) const;
  ring_elem multiplyTermViaTargetClosedFamily(
      const TargetClosedMultiplicationFamily& family,
      const std::vector<ProductFactor>& factors,
      int targetBasisId) const;
  ring_elem multiplyTermViaPowerSums(
      const std::vector<ProductFactor>& factors,
      int targetBasisId) const;
  ResolvedProductTerm multiplyTermToBasis(
      const SymmetricTerm& term,
      int targetBasisId) const;
  bool multiplicationWorkflowTraceEnabled() const;

 public:
  ring_elem multiplyToBasis(
      ring_elem f,
      ring_elem g,
      int targetBasisId) const;
  ring_elem multiplyToBasisBench(
      ring_elem f,
      ring_elem g,
      int targetBasisId,
      const std::optional<std::string>& forcedKernel,
      bool usePowerSumReference,
      bool traceWorkflow) const;
  ring_elem multiplyExpressionsToBasisBench(
      ring_elem f,
      ring_elem g,
      int targetBasisId,
      const std::string& strategy,
      bool traceWorkflow) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
