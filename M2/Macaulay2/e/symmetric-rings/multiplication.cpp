// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"

#include <algorithm>
#include <cstdio>
#include <utility>

namespace symmetric_rings {

// This file owns bilinear extension, complete product-term resolution, and the
// generic multifactor strategies. It does not declare kernels, multiplication
// families, or binary-picker records.

namespace {

// Development expression benchmarks install trace state for one synchronous
// call. Ordinary multiplication has no mutable or process-global forcing
// state and retains the existing outer-pipeline environment trace.
thread_local const bool *activeMultiplicationWorkflowTrace = nullptr;

class ScopedMultiplicationWorkflowTrace
{
 public:
  explicit ScopedMultiplicationWorkflowTrace(const bool& enabled)
      : mPrevious(activeMultiplicationWorkflowTrace)
  {
    activeMultiplicationWorkflowTrace = &enabled;
  }

  ~ScopedMultiplicationWorkflowTrace()
  {
    activeMultiplicationWorkflowTrace = mPrevious;
  }

 private:
  const bool *mPrevious;
};

} // namespace

bool SymmetricEngineRing::multiplicationWorkflowTraceEnabled() const
{
    return activeMultiplicationWorkflowTrace != nullptr
        ? *activeMultiplicationWorkflowTrace
        : basisConversionTraceEnabled();
  }

// ============================================================================
// Policy-Free Result And Product Helpers
// ============================================================================

ExpressionFacts
SymmetricEngineRing::finalizeCanonicalMultiplicationResult(
    ring_elem result,
    int targetBasisId,
    CombinatorialTags combinatorialTags) const
{
    ExpressionFacts facts =
        inferCanonicalExpansionFacts(
            result, targetBasisId);
    if (!facts.canonicalExpansionInBasis(
            targetBasisId))
      {
        ERROR("multiplication violated its canonical target-basis "
              "contract");
        return facts;
      }
    attachCanonicalExpansionFacts(
        result,
        facts,
        targetBasisId,
        combinatorialTags);
    return facts;
  }

ring_elem
SymmetricEngineRing::multiplyCanonicalExpansionsBalanced(
    std::vector<ring_elem> factors) const
{
    if (factors.empty())
      {
        ERROR("balanced multiplication requires at least one factor");
        return zero();
      }
    while (factors.size() > 1)
      {
        std::vector<ring_elem> next;
        next.reserve((factors.size() + 1) / 2);
        for (size_t position = 0;
             position < factors.size();
             position += 2)
          if (position + 1 == factors.size())
            next.push_back(factors[position]);
          else
            {
              ring_elem product = mult(
                  factors[position],
                  factors[position + 1]);
              if (error()) return zero();
              next.push_back(product);
            }
        factors = std::move(next);
      }
    return factors.front();
  }

// ============================================================================
// Product-Free Operand Preparation
// ============================================================================

SymmetricEngineRing::PreparedMultiplicationOperands
SymmetricEngineRing::prepareMultiplicationOperands(
    ring_elem leftInput,
    ring_elem rightInput) const
{
    PreparedMultiplicationOperands prepared;

    auto prepareExpression =
        [&](ring_elem input,
            ring_elem& expression,
            ExpressionFacts& facts,
            std::vector<PreparedMultiplicationTerm>& terms) {
          const ExpressionNormalizationResult normalized =
              normalizeExpressionWithOptionsDetailed(
                  input,
                  pipelinePreparationNormalizationOptions());
          expression = normalized.expression;
          facts = normalized.facts;
          if (error()) return;
          if (!facts.noProducts())
            {
              ERROR("multiplyToBasis expects product-free operands");
              return;
            }

          terms.reserve(polyValue(expression)->terms.size());
          for (const auto& term :
               polyValue(expression)->terms)
            {
              const bool scalar =
                  term.monomial.data.empty();
              if (!scalar &&
                  atomLengthAt(term.monomial, 0) !=
                      term.monomial.data.size())
                {
                  ERROR("multiplyToBasis encountered a product term "
                        "after product-free preparation");
                  return;
                }
              ring_elem factor =
                  scalar
                      ? one()
                      : expressionFromBasisElement(
                            term.monomial, 0);
              ExpressionFacts factorFacts;
              if (!scalar)
                factorFacts = basisElementFactsFromAtom(
                    term.monomial,
                    0,
                    coefficientRing->one());
              terms.push_back({
                  &term,
                  scalar,
                  factor,
                  std::move(factorFacts)});
            }
        };

    prepareExpression(
        leftInput,
        prepared.left,
        prepared.leftFacts,
        prepared.leftTerms);
    if (error()) return prepared;
    prepareExpression(
        rightInput,
        prepared.right,
        prepared.rightFacts,
        prepared.rightTerms);
    if (error()) return prepared;
    prepared.combinatorialTags =
        selectMultiplicationTags(
            prepared.left, prepared.right);
    return prepared;
  }

ring_elem SymmetricEngineRing::convertCanonicalFactorToBasis(
    ring_elem factor,
    const ExpressionFacts& facts,
    int targetBasisId,
    CombinatorialTags tags) const
{
    if (!facts.expandedBasis)
      {
        ERROR("a multiplication factor has no canonical source basis");
        return zero();
      }
    if (*facts.expandedBasis == targetBasisId)
      return factor;
    return convertCanonicalExpressionToBasis(
        factor,
        targetBasisId,
        tags,
        &facts);
  }

ring_elem SymmetricEngineRing::convertPreparedMultiplicationTerm(
    const PreparedMultiplicationTerm& term,
    int targetBasisId,
    CombinatorialTags tags) const
{
    if (term.scalar) return one();
    return convertCanonicalFactorToBasis(
        term.factor, term.facts, targetBasisId, tags);
  }

// ============================================================================
// Complete-Input Structural Strategies
// ============================================================================

ring_elem
SymmetricEngineRing::multiplyExpressionsViaMultiplicativeTarget(
    const PreparedMultiplicationOperands& operands,
    int targetBasisId) const
{
    ring_elem left =
        convertCanonicalExpressionToBasis(
            operands.left,
            targetBasisId,
            operands.combinatorialTags,
            &operands.leftFacts);
    if (error()) return zero();
    ring_elem right =
        convertCanonicalExpressionToBasis(
            operands.right,
            targetBasisId,
            operands.combinatorialTags,
            &operands.rightFacts);
    if (error()) return zero();
    ring_elem result = mult(left, right);
    if (error()) return zero();
    finalizeCanonicalMultiplicationResult(
        result,
        targetBasisId,
        operands.combinatorialTags);
    return result;
  }

ring_elem SymmetricEngineRing::multiplyExpressionsViaPowerSums(
    const PreparedMultiplicationOperands& operands,
    int targetBasisId) const
{
    const int powerSumBasisId =
        requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return zero();
    ring_elem left =
        convertCanonicalExpressionToBasis(
            operands.left,
            powerSumBasisId,
            operands.combinatorialTags,
            &operands.leftFacts);
    if (error()) return zero();
    ring_elem right =
        convertCanonicalExpressionToBasis(
            operands.right,
            powerSumBasisId,
            operands.combinatorialTags,
            &operands.rightFacts);
    if (error()) return zero();
    ring_elem product = mult(left, right);
    if (error()) return zero();
    ExpressionFacts powerSumFacts =
        inferCanonicalExpansionFacts(
            product, powerSumBasisId);
    if (powerSumBasisId == targetBasisId)
      {
        finalizeCanonicalMultiplicationResult(
            product,
            targetBasisId,
            operands.combinatorialTags);
        return product;
      }
    return convertCanonicalExpressionToBasis(
        product,
        targetBasisId,
        operands.combinatorialTags,
        &powerSumFacts);
  }

bool
SymmetricEngineRing::everyNonscalarPairHasAutomaticKernel(
    const PreparedMultiplicationOperands& operands,
    int targetBasisId) const
{
    for (const auto& left : operands.leftTerms)
      for (const auto& right : operands.rightTerms)
        {
          if (left.scalar || right.scalar) continue;
          const auto leftView =
              canonicalBasisTermView(
                  left.factor, left.facts);
          if (error()) return false;
          const auto rightView =
              canonicalBasisTermView(
                  right.factor, right.facts);
          if (error()) return false;
          if (!automaticBinaryMultiplicationKernelAvailable(
                  leftView, rightView, targetBasisId))
            return false;
        }
    return true;
  }

ring_elem SymmetricEngineRing::distributeMultiplicationOverTermPairs(
    const PreparedMultiplicationOperands& operands,
    int targetBasisId) const
{
    VECTOR(SymmetricTerm) resultTerms;
    for (const auto& left : operands.leftTerms)
      for (const auto& right : operands.rightTerms)
        {
          ring_elem coefficient =
              coefficientRing->mult(
                  left.term->coeff, right.term->coeff);
          if (coefficientRing->is_zero(coefficient))
            continue;

          ring_elem product;
          if (left.scalar && right.scalar)
            product = one();
          else if (left.scalar || right.scalar)
            product = convertPreparedMultiplicationTerm(
                left.scalar ? right : left,
                targetBasisId,
                operands.combinatorialTags);
          else
            {
              const auto leftView =
                  canonicalBasisTermView(
                      left.factor, left.facts);
              if (error()) return zero();
              const auto rightView =
                  canonicalBasisTermView(
                      right.factor, right.facts);
              if (error()) return zero();
              product = multiplyCanonicalBasisTerms(
                  leftView,
                  rightView,
                  targetBasisId);
            }
          if (error()) return zero();

          for (const auto& term :
               polyValue(product)->terms)
            appendTermIfNonZero(
                resultTerms,
                coefficientRing->mult(
                    coefficient, term.coeff),
                term.monomial);
        }

    ring_elem result =
        fromTermVector(resultTerms, false);
    finalizeCanonicalMultiplicationResult(
        result,
        targetBasisId,
        operands.combinatorialTags);
    return result;
  }

// ============================================================================
// Public Bilinear Workflow
// ============================================================================

ring_elem SymmetricEngineRing::multiplyToBasis(
    ring_elem left,
    ring_elem right,
    int targetBasisId) const
{
    return multiplyToBasisWithStrategy(
        left,
        right,
        targetBasisId,
        BilinearMultiplicationStrategy::Automatic);
  }

ring_elem SymmetricEngineRing::multiplyToBasisWithStrategy(
    ring_elem left,
    ring_elem right,
    int targetBasisId,
    BilinearMultiplicationStrategy strategy) const
{
    if (binaryMultiplicationKernelExecutionDepth != 0)
      {
        ERROR("a binary multiplication kernel attempted to call "
              "the bilinear multiplication workflow");
        return zero();
      }
    requireBasis(targetBasisId);
    if (error()) return zero();
    PreparedMultiplicationOperands operands =
        prepareMultiplicationOperands(left, right);
    if (error()) return zero();

    if (operands.leftTerms.empty() ||
        operands.rightTerms.empty())
      {
        ring_elem result = zero();
        finalizeCanonicalMultiplicationResult(
            result,
            targetBasisId,
            operands.combinatorialTags);
        return result;
      }

    const bool leftIsScalar =
        operands.leftFacts.singleFactorTermCount == 0;
    const bool rightIsScalar =
        operands.rightFacts.singleFactorTermCount == 0;
    if (leftIsScalar || rightIsScalar)
      {
        if (leftIsScalar && rightIsScalar)
          {
            ring_elem result =
                mult(operands.left, operands.right);
            finalizeCanonicalMultiplicationResult(
                result,
                targetBasisId,
                operands.combinatorialTags);
            return result;
          }
        const ring_elem scalar =
            leftIsScalar
                ? operands.leftTerms.front().term->coeff
                : operands.rightTerms.front().term->coeff;
        ring_elem expression =
            leftIsScalar ? operands.right : operands.left;
        const ExpressionFacts& expressionFacts =
            leftIsScalar
                ? operands.rightFacts
                : operands.leftFacts;
        ring_elem converted =
            convertCanonicalExpressionToBasis(
                expression,
                targetBasisId,
                operands.combinatorialTags,
                &expressionFacts);
        if (error()) return zero();
        ring_elem result = scaled(scalar, converted);
        finalizeCanonicalMultiplicationResult(
            result,
            targetBasisId,
            operands.combinatorialTags);
        return result;
      }

    const bool multiplicativeTarget =
        isMultiplicativeBasis(targetBasisId);
    if (strategy ==
            BilinearMultiplicationStrategy::
                MultiplicativeTarget &&
        !multiplicativeTarget)
      {
        ERROR("the forced multiplicative-target strategy requires "
              "a multiplicative target basis");
        return zero();
      }
    if (strategy ==
            BilinearMultiplicationStrategy::
                KernelDistribution &&
        multiplicativeTarget)
      {
        ERROR("the kernel-distribution benchmark strategy is not "
              "defined for a multiplicative target");
        return zero();
      }
    if (strategy ==
        BilinearMultiplicationStrategy::PowerSumFallback)
      {
        if (multiplicationWorkflowTraceEnabled())
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-expression: "
              "strategy=forced-complete-input-power-sums "
              "strict-binary=no target=%s\n",
              basisKeyForId(targetBasisId).c_str());
        return multiplyExpressionsViaPowerSums(
            operands, targetBasisId);
      }

    if (multiplicativeTarget &&
        strategy !=
            BilinearMultiplicationStrategy::
                KernelDistribution)
      {
        if (multiplicationWorkflowTraceEnabled())
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-expression: "
              "strategy=multiplicative-target "
              "strict-binary=no target=%s\n",
              basisKeyForId(targetBasisId).c_str());
        return multiplyExpressionsViaMultiplicativeTarget(
            operands, targetBasisId);
      }

    const bool everyPairHasKernel =
        everyNonscalarPairHasAutomaticKernel(
            operands, targetBasisId);
    if (error()) return zero();
    if (strategy ==
            BilinearMultiplicationStrategy::
                KernelDistribution &&
        !everyPairHasKernel)
      {
        ERROR("the forced kernel-distribution strategy requires "
              "an automatic binary kernel for every nonscalar term pair");
        return zero();
      }
    if (!everyPairHasKernel)
      {
        if (multiplicationWorkflowTraceEnabled())
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-expression: "
              "strategy=complete-input-power-sums "
              "strict-binary=no target=%s\n",
              basisKeyForId(targetBasisId).c_str());
        return multiplyExpressionsViaPowerSums(
            operands, targetBasisId);
      }

    if (multiplicationWorkflowTraceEnabled())
      std::fprintf(
          stderr,
          "SymmetricRings multiplication-expression: "
          "strategy=direct-kernel-distribution "
          "strict-binary=yes target=%s\n",
          basisKeyForId(targetBasisId).c_str());
    return distributeMultiplicationOverTermPairs(
        operands, targetBasisId);
  }

// ============================================================================
// Development Strict-Binary Entry Point
// ============================================================================

ring_elem SymmetricEngineRing::multiplyToBasisBench(
    ring_elem left,
    ring_elem right,
    int targetBasisId,
    const std::optional<std::string>& forcedKernel,
    bool usePowerSumReference,
    bool traceWorkflow) const
{
    PreparedMultiplicationOperands operands =
        prepareMultiplicationOperands(left, right);
    if (error()) return zero();
    if (!operands.leftFacts.singleBasisElement() ||
        !operands.rightFacts.singleBasisElement() ||
        !operands.leftFacts.
            singleBasisElementCoefficientOne.value_or(false) ||
        !operands.rightFacts.
            singleBasisElementCoefficientOne.value_or(false))
      {
        ERROR("multiplyToBasisBench requires two coefficient-one "
              "canonical basis terms");
        return zero();
      }
    const auto leftView =
        canonicalBasisTermView(
            operands.left, operands.leftFacts);
    if (error()) return zero();
    const auto rightView =
        canonicalBasisTermView(
            operands.right, operands.rightFacts);
    if (error()) return zero();
    return multiplyCanonicalBasisTerms(
        leftView,
        rightView,
        targetBasisId,
        {forcedKernel,
         usePowerSumReference,
         traceWorkflow});
  }

ring_elem SymmetricEngineRing::multiplyExpressionsToBasisBench(
    ring_elem left,
    ring_elem right,
    int targetBasisId,
    const std::string& strategy,
    bool traceWorkflow) const
{
    BilinearMultiplicationStrategy requestedStrategy =
        BilinearMultiplicationStrategy::Automatic;
    if (strategy == "Automatic")
      requestedStrategy =
          BilinearMultiplicationStrategy::Automatic;
    else if (strategy == "MultiplicativeTarget")
      requestedStrategy =
          BilinearMultiplicationStrategy::
              MultiplicativeTarget;
    else if (strategy == "KernelDistribution")
      requestedStrategy =
          BilinearMultiplicationStrategy::
              KernelDistribution;
    else if (strategy == "PowerSumFallback")
      requestedStrategy =
          BilinearMultiplicationStrategy::
              PowerSumFallback;
    else
      {
        ERROR("unknown bilinear multiplication benchmark strategy");
        return zero();
      }
    const ScopedMultiplicationWorkflowTrace traceScope(
        traceWorkflow);
    return multiplyToBasisWithStrategy(
        left,
        right,
        targetBasisId,
        requestedStrategy);
  }

// ============================================================================
// Complete Product-Term Preparation
// ============================================================================

std::vector<SymmetricEngineRing::ProductFactor>
SymmetricEngineRing::productFactors(
    const SymmetricTerm& term) const
{
    std::vector<ProductFactor> factors;
    size_t pos = 0;
    while (pos < term.monomial.data.size())
      {
        bool identity = false;
        if (!atomIsSkewAt(term.monomial, pos))
          {
            identity = true;
            const size_t indexLength =
                static_cast<size_t>(
                    atomIndexLengthAt(
                        term.monomial, pos));
            for (size_t i = 0; i < indexLength; ++i)
              if (term.monomial.data[
                      pos + atomHeaderSize + i] != 0)
                {
                  identity = false;
                  break;
                }
          }
        if (!identity)
          factors.push_back({
              expressionFromBasisElement(
                  term.monomial, pos),
              basisElementFactsFromAtom(
                  term.monomial,
                  pos,
                  coefficientRing->one())});
        pos += atomLengthAt(term.monomial, pos);
      }
    return factors;
  }

std::vector<int>
SymmetricEngineRing::basisIdsOfProductFactors(
    const std::vector<ProductFactor>& factors) const
{
    std::vector<int> basisIds;
    basisIds.reserve(factors.size());
    for (const auto& factor : factors)
      {
        if (!factor.facts.expandedBasis)
          {
            ERROR("a product factor has no canonical source basis");
            return {};
          }
        basisIds.push_back(*factor.facts.expandedBasis);
      }
    return basisIds;
  }

// ============================================================================
// Complete Product-Term Strategies
// ============================================================================

ring_elem
SymmetricEngineRing::multiplyTermViaMultiplicativeTarget(
    const std::vector<ProductFactor>& factors,
    int targetBasisId) const
{
    std::vector<ring_elem> level;
    level.reserve(factors.size());
    for (const auto& factor : factors)
      {
        ring_elem converted =
            convertCanonicalFactorToBasis(
                factor.expression,
                factor.facts,
                targetBasisId,
                polyValue(factor.expression)->
                    combinatorialTags);
        if (error()) return zero();
        level.push_back(converted);
      }

    return multiplyCanonicalExpansionsBalanced(
        std::move(level));
  }

ring_elem SymmetricEngineRing::multiplyTermViaTargetClosedFamily(
    const TargetClosedMultiplicationFamily& family,
    const std::vector<ProductFactor>& factors,
    int targetBasisId) const
{
    if (family.targetBasis != basisKindForId(targetBasisId))
      {
        ERROR("a target-closed multiplication family received the wrong "
              "target basis");
        return zero();
      }

    size_t seed = 0;
    // Every target-closed family uses the same fixed association: begin with
    // the first target-native factor when one exists, avoiding an unnecessary
    // conversion and keeping the first accumulator small. Otherwise convert
    // the first stored factor to the target. No factor ordering is searched.
    auto native = std::find_if(
        factors.begin(),
        factors.end(),
        [&](const ProductFactor& factor) {
          return factor.facts.expandedBasis &&
                 *factor.facts.expandedBasis == targetBasisId;
        });
    if (native != factors.end())
      seed = static_cast<size_t>(
          std::distance(factors.begin(), native));

    ring_elem result =
        convertCanonicalFactorToBasis(
            factors[seed].expression,
            factors[seed].facts,
            targetBasisId,
            polyValue(factors[seed].expression)->
                combinatorialTags);
    if (error()) return zero();
    for (size_t i = 0; i < factors.size(); ++i)
      {
        if (i == seed) continue;
        result = multiplyToBasis(
            result,
            factors[i].expression,
            targetBasisId);
        if (error()) return zero();
      }
    return result;
  }

ring_elem SymmetricEngineRing::multiplyTermViaPowerSums(
    const std::vector<ProductFactor>& factors,
    int targetBasisId) const
{
    const int powerSumBasisId =
        requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return zero();
    std::vector<ring_elem> powerSumFactors;
    powerSumFactors.reserve(factors.size());
    for (const auto& factor : factors)
      {
        ring_elem converted =
            convertCanonicalFactorToBasis(
                factor.expression,
                factor.facts,
                powerSumBasisId,
                polyValue(factor.expression)->
                    combinatorialTags);
        if (error()) return zero();
        powerSumFactors.push_back(converted);
      }

    ring_elem powerSumProduct =
        multiplyCanonicalExpansionsBalanced(
            std::move(powerSumFactors));
    if (error()) return zero();
    ExpressionFacts powerSumFacts =
        inferCanonicalExpansionFacts(
            powerSumProduct, powerSumBasisId);
    return convertCanonicalExpressionToBasis(
        powerSumProduct,
        targetBasisId,
        polyValue(powerSumProduct)->
            combinatorialTags,
        &powerSumFacts);
  }

// ============================================================================
// Complete Product-Term Workflow
// ============================================================================

SymmetricEngineRing::ResolvedProductTerm
SymmetricEngineRing::multiplyTermToBasis(
    const SymmetricTerm& term,
    int targetBasisId) const
{
    if (binaryMultiplicationKernelExecutionDepth != 0)
      {
        ERROR("a binary multiplication kernel attempted to call "
              "the multifactor multiplication workflow");
        return {zero(), {}};
      }
    std::vector<ProductFactor> factors =
        productFactors(term);
    const TargetClosedMultiplicationFamily *closedFamily =
        nullptr;
    if (factors.size() >= 2 &&
        !isMultiplicativeBasis(targetBasisId))
      {
        const std::vector<int> factorBasisIds =
            basisIdsOfProductFactors(factors);
        if (error()) return {zero(), {}};
        closedFamily =
            selectTargetClosedMultiplicationFamily(
                factorBasisIds, targetBasisId);
      }
    ring_elem result;

    if (factors.empty())
      {
        if (multiplicationWorkflowTraceEnabled())
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-term: "
              "strategy=scalar strict-binary=no target=%s\n",
              basisKeyForId(targetBasisId).c_str());
        result = one();
      }
    else if (factors.size() == 1)
      {
        if (multiplicationWorkflowTraceEnabled())
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-term: "
              "strategy=single-factor "
              "strict-binary=no target=%s\n",
              basisKeyForId(targetBasisId).c_str());
        result = convertCanonicalFactorToBasis(
            factors.front().expression,
            factors.front().facts,
            targetBasisId,
            polyValue(factors.front().expression)->
                combinatorialTags);
      }
    else if (isMultiplicativeBasis(targetBasisId))
      {
        if (multiplicationWorkflowTraceEnabled())
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-term: "
              "strategy=multiplicative-target-balanced "
              "strict-binary=no factors=%zu target=%s\n",
              factors.size(),
              basisKeyForId(targetBasisId).c_str());
        result = multiplyTermViaMultiplicativeTarget(
            factors, targetBasisId);
      }
    else if (closedFamily != nullptr)
      {
        if (multiplicationWorkflowTraceEnabled())
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-term: "
              "strategy=target-closed-fold "
              "strict-binary=yes factors=%zu family=%s target=%s\n",
              factors.size(),
              closedFamily->identifier.c_str(),
              basisKeyForId(targetBasisId).c_str());
        result = multiplyTermViaTargetClosedFamily(
            *closedFamily, factors, targetBasisId);
      }
    else
      {
        if (multiplicationWorkflowTraceEnabled())
          std::fprintf(
              stderr,
              "SymmetricRings multiplication-term: "
              "strategy=complete-term-power-sums "
              "strict-binary=no factors=%zu target=%s\n",
              factors.size(),
              basisKeyForId(targetBasisId).c_str());
        result = multiplyTermViaPowerSums(
            factors, targetBasisId);
      }
    if (error()) return {zero(), {}};

    if (!coefficientRing->is_equal(
            term.coeff, coefficientRing->one()))
      result = scaled(term.coeff, result);
    if (error()) return {zero(), {}};

    ExpressionFacts facts =
        finalizeCanonicalMultiplicationResult(
        result,
        targetBasisId,
        polyValue(result)->combinatorialTags);
    return {result, std::move(facts)};
  }

} // namespace symmetric_rings
