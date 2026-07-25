// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"

#include <cstdio>

namespace symmetric_rings {

// Strict binary multiplication accepts exactly two coefficient-one canonical
// basis terms. It owns no bilinear distribution and no multifactor policy.

// ============================================================================
// Strict Input Boundary
// ============================================================================

SymmetricEngineRing::CanonicalBasisTermView
SymmetricEngineRing::canonicalBasisTermView(
    ring_elem expression,
    const ExpressionFacts& facts) const
{
    if (!facts.singleBasisElement() ||
        !facts.singleBasisElementId ||
        !facts.singleBasisElementIndex ||
        !facts.singleBasisElementCoefficientOne.value_or(false) ||
        !facts.canonicalExpansionInBasis(
            *facts.singleBasisElementId))
      {
        ERROR("strict binary multiplication requires one "
              "coefficient-one canonical basis term");
        return {};
      }
    const auto& descriptor =
        requireBasis(*facts.singleBasisElementId);
    if (error()) return {};
    return {
        expression,
        {descriptor.kind,
         *facts.singleBasisElementId,
         &descriptor.displaySymbol,
         descriptor.displayOrder},
        *facts.singleBasisElementIndex};
  }

ring_elem SymmetricEngineRing::executeBinaryMultiplicationKernel(
    const SelectedBinaryMultiplicationKernel& selection,
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right) const
{
    if (!selection.valid() ||
        selection.definition->kernel == nullptr)
      {
        ERROR("cannot execute an empty binary multiplication selection");
        return zero();
      }
    const CanonicalBasisTermView& first =
        selection.callableTakesReversedArguments
            ? right : left;
    const CanonicalBasisTermView& second =
        selection.callableTakesReversedArguments
            ? left : right;
    const BinaryMultiplicationInput input{
        first, second, selection.target};

    struct KernelExecutionGuard
    {
      size_t& depth;
      explicit KernelExecutionGuard(size_t& value)
          : depth(value)
      {
        ++depth;
      }
      ~KernelExecutionGuard() { --depth; }
    } guard(binaryMultiplicationKernelExecutionDepth);
    ring_elem result =
        (this->*selection.definition->kernel)(input);
    if (error()) return zero();
    finalizeCanonicalMultiplicationResult(
        result,
        selection.target.id,
        selectMultiplicationTags(
            left.expression, right.expression));
    return error() ? zero() : result;
  }

void SymmetricEngineRing::validateBinaryMultiplicationRequest(
    const BinaryMultiplicationRequest& request) const
{
    if (request.forcedKernel &&
        request.usePowerSumReference)
      ERROR("a binary multiplication request cannot both force "
            "a kernel and request the power-sum reference");
  }

void SymmetricEngineRing::traceBinaryMultiplicationWorkflow(
    const char *workflow,
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right,
    int targetBasisId,
    const BinaryMultiplicationRequest& request) const
{
    if (!request.traceWorkflow &&
        !basisConversionTraceEnabled())
      return;
    std::fprintf(
        stderr,
        "SymmetricRings binary-multiplication: "
        "endpoint={%s,%s}->%s workflow=%s\n",
        basisKindName(
            UnorderedBasisPair::from(
                left.basis.kind,
                right.basis.kind).first),
        basisKindName(
            UnorderedBasisPair::from(
                left.basis.kind,
                right.basis.kind).second),
        basisKeyForId(targetBasisId).c_str(),
        workflow);
  }

// ============================================================================
// Fixed Structural Workflows
// ============================================================================

ring_elem
SymmetricEngineRing::multiplyBasisTermsViaMultiplicativeTarget(
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right,
    int targetBasisId) const
{
    ExpressionFacts leftFacts =
        inferCanonicalExpansionFacts(
            left.expression, left.basis.id);
    ExpressionFacts rightFacts =
        inferCanonicalExpansionFacts(
            right.expression, right.basis.id);
    ring_elem convertedLeft =
        convertCanonicalExpressionToBasis(
            left.expression,
            targetBasisId,
            polyValue(left.expression)->combinatorialTags,
            &leftFacts);
    if (error()) return zero();
    ring_elem convertedRight =
        convertCanonicalExpressionToBasis(
            right.expression,
            targetBasisId,
            polyValue(right.expression)->combinatorialTags,
            &rightFacts);
    if (error()) return zero();

    ring_elem result = mult(convertedLeft, convertedRight);
    if (error()) return zero();
    const CombinatorialTags tags =
        selectMultiplicationTags(
            left.expression, right.expression);
    finalizeCanonicalMultiplicationResult(
        result, targetBasisId, tags);
    return error() ? zero() : result;
  }

ring_elem SymmetricEngineRing::multiplyBasisTermsViaPowerSums(
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right,
    int targetBasisId) const
{
    const int powerSumBasisId =
        requiredBasisIdForKind(BasisKind::PowerSum);
    if (error()) return zero();

    ExpressionFacts leftFacts =
        inferCanonicalExpansionFacts(
            left.expression, left.basis.id);
    ExpressionFacts rightFacts =
        inferCanonicalExpansionFacts(
            right.expression, right.basis.id);
    ring_elem leftPowerSums =
        convertCanonicalExpressionToBasis(
            left.expression,
            powerSumBasisId,
            polyValue(left.expression)->combinatorialTags,
            &leftFacts);
    if (error()) return zero();
    ring_elem rightPowerSums =
        convertCanonicalExpressionToBasis(
            right.expression,
            powerSumBasisId,
            polyValue(right.expression)->combinatorialTags,
            &rightFacts);
    if (error()) return zero();

    ring_elem powerSumProduct =
        mult(leftPowerSums, rightPowerSums);
    if (error()) return zero();
    const CombinatorialTags tags =
        selectMultiplicationTags(
            left.expression, right.expression);
    ExpressionFacts productFacts =
        inferCanonicalExpansionFacts(
            powerSumProduct, powerSumBasisId);
    if (powerSumBasisId == targetBasisId)
      {
        finalizeCanonicalMultiplicationResult(
            powerSumProduct,
            targetBasisId,
            tags);
        return error() ? zero() : powerSumProduct;
      }
    return convertCanonicalExpressionToBasis(
        powerSumProduct,
        targetBasisId,
        tags,
        &productFacts);
  }

// ============================================================================
// Strict Binary Workflow
// ============================================================================

ring_elem SymmetricEngineRing::multiplyCanonicalBasisTerms(
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right,
    int targetBasisId) const
{
    return multiplyCanonicalBasisTerms(
        left,
        right,
        targetBasisId,
        BinaryMultiplicationRequest{});
  }

ring_elem SymmetricEngineRing::multiplyCanonicalBasisTerms(
    const CanonicalBasisTermView& left,
    const CanonicalBasisTermView& right,
    int targetBasisId,
    const BinaryMultiplicationRequest& request) const
{
    if (binaryMultiplicationKernelExecutionDepth != 0)
      {
        ERROR("a binary multiplication kernel attempted to reenter "
              "strict binary multiplication");
        return zero();
      }
    validateBinaryMultiplicationRequest(request);
    if (error()) return zero();
    requireBasis(targetBasisId);
    if (error()) return zero();

    if (request.forcedKernel)
      {
        auto selection =
            selectBinaryMultiplicationKernel(
                left, right, targetBasisId, request, true);
        if (error() || !selection) return zero();
        traceBinaryMultiplicationWorkflow(
            "forced-kernel",
            left,
            right,
            targetBasisId,
            request);
        ring_elem result =
            executeBinaryMultiplicationKernel(
                *selection, left, right);
        return error() ? zero() : result;
      }

    if (request.usePowerSumReference)
      {
        traceBinaryMultiplicationWorkflow(
            "power-sum-reference",
            left,
            right,
            targetBasisId,
            request);
        return multiplyBasisTermsViaPowerSums(
            left, right, targetBasisId);
      }

    if (isMultiplicativeBasis(targetBasisId))
      {
        traceBinaryMultiplicationWorkflow(
            "multiplicative-target",
            left,
            right,
            targetBasisId,
            request);
        return multiplyBasisTermsViaMultiplicativeTarget(
            left, right, targetBasisId);
      }

    if (auto selection =
            selectBinaryMultiplicationKernel(
                left, right, targetBasisId, request, true))
      {
        traceBinaryMultiplicationWorkflow(
            "direct-kernel",
            left,
            right,
            targetBasisId,
            request);
        ring_elem result =
            executeBinaryMultiplicationKernel(
                *selection, left, right);
        return error() ? zero() : result;
      }
    if (error()) return zero();

    traceBinaryMultiplicationWorkflow(
        "power-sum-fallback",
        left,
        right,
        targetBasisId,
        request);
    return multiplyBasisTermsViaPowerSums(
        left, right, targetBasisId);
  }

} // namespace symmetric_rings
