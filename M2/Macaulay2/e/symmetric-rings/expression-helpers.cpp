// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"

#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Shared Result Construction
// ============================================================================

ring_elem SymmetricEngineRing::expressionHelperFromTerms(
    VECTOR(SymmetricTerm)& terms,
    CombinatorialTags tags) const
{
    ring_elem result = fromTermVector(terms, false);
    mutablePolyValue(result)->combinatorialTags = tags;
    const ExpressionFacts facts = inferExpressionFacts(result);
    attachExpressionHelperFacts(result, facts);
    return result;
  }

void SymmetricEngineRing::attachExpressionHelperFacts(
    ring_elem expression,
    const ExpressionFacts& facts) const
{
    SymmetricConversionMetadata metadata;
    metadata.pureBasis = facts.pureBasis;
    if (facts.expandedBasis &&
        facts.canonicalExpansionInBasis(
            *facts.expandedBasis))
      metadata.expandedBasis = facts.expandedBasis;
    metadata.homogeneousWeight = facts.homogeneousWeight;
    metadata.termCount = facts.termCount;
    metadata.scalarTermCount = facts.scalarTermCount;
    metadata.singleFactorTermCount =
        facts.singleFactorTermCount;
    metadata.maximumPartitionLength =
        facts.maximumPartitionLength;
    metadata.factorBases = facts.factorBases;
    metadata.normalized = facts.normalized;
    metadata.skewFree = facts.skewFree;
    metadata.collected = facts.collected;
    // The existing exact-metadata bypass represents product-free canonical
    // storage. Product-bearing results retain valid step postconditions and
    // support hints, but consumers must still inspect their factor counts.
    metadata.expressionFactsComplete =
        facts.normalized &&
        facts.skewFree &&
        facts.collected &&
        facts.noProducts();
    if (metadata.expressionFactsComplete &&
        facts.singleBasisElement())
      {
        metadata.singleBasisElementId =
            facts.singleBasisElementId;
        metadata.singleBasisElementIndex =
            facts.singleBasisElementIndex;
        metadata.singleBasisElementCoefficientOne =
            facts.singleBasisElementCoefficientOne;
      }
    auto *poly = mutablePolyValue(expression);
    poly->combinatorialTags =
        facts.combinatorialTags;
    poly->conversionMetadata =
        std::move(metadata);
  }

// ============================================================================
// Homogeneous Decomposition
// ============================================================================

std::vector<SymmetricEngineRing::HomogeneousComponent>
SymmetricEngineRing::homogeneousComponents(
    ring_elem expression) const
{
    std::map<int, VECTOR(SymmetricTerm)> termsByWeight;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      termsByWeight[monomialWeight(term.monomial)].push_back(term);

    std::vector<HomogeneousComponent> result;
    result.reserve(termsByWeight.size());
    for (auto& item : termsByWeight)
      result.push_back({
          item.first,
          expressionHelperFromTerms(
              item.second, poly->combinatorialTags)});
    return result;
  }

ring_elem SymmetricEngineRing::homogeneousComponent(
    ring_elem expression,
    int weight) const
{
    VECTOR(SymmetricTerm) selected;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      if (monomialWeight(term.monomial) == weight)
        selected.push_back(term);
    return expressionHelperFromTerms(
        selected, poly->combinatorialTags);
  }

std::vector<int> SymmetricEngineRing::weightSupport(
    ring_elem expression) const
{
    std::vector<int> result;
    const auto components = homogeneousComponents(expression);
    result.reserve(components.size());
    for (const auto& component : components)
      result.push_back(component.weight);
    return result;
  }

ring_elem SymmetricEngineRing::truncateWeights(
    ring_elem expression,
    int minimumWeight,
    int maximumWeight) const
{
    if (minimumWeight > maximumWeight)
      {
        ERROR("expected the minimum weight to be at most the maximum weight");
        return zero();
      }
    VECTOR(SymmetricTerm) selected;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      {
        const int weight = monomialWeight(term.monomial);
        if (minimumWeight <= weight && weight <= maximumWeight)
          selected.push_back(term);
      }
    return expressionHelperFromTerms(
        selected, poly->combinatorialTags);
  }

// ============================================================================
// Normalization And Product Resolution
// ============================================================================

SymmetricEngineRing::ExpressionNormalizationOptions
SymmetricEngineRing::pipelinePreparationNormalizationOptions()
{
    ExpressionNormalizationOptions options;
    options.straightenIndices = true;
    options.expandSkewFactors = true;
    return options;
  }

ring_elem SymmetricEngineRing::normalizeExpressionWithOptions(
    ring_elem expression,
    const ExpressionNormalizationOptions& options) const
{
    return normalizeExpressionWithOptionsDetailed(
        expression, options).expression;
  }

SymmetricEngineRing::ExpressionNormalizationResult
SymmetricEngineRing::normalizeExpressionWithOptionsDetailed(
    ring_elem expression,
    const ExpressionNormalizationOptions& options) const
{
    ExpressionNormalizationResult result;
    const CombinatorialTags inputTags =
        polyValue(expression)->combinatorialTags;
    if (options.productTargetBasisId >= 0)
      {
        (void) requireBasis(options.productTargetBasisId);
        if (error())
          {
            result.expression = zero();
            return result;
          }
        result.expression = expandProductsInBasis(
            expression, options.productTargetBasisId);
        if (error())
          {
            result.expression = zero();
            return result;
          }
        result.resolvedProducts = true;
        auto exactFacts =
            expressionFactsFromMetadata(result.expression);
        if (exactFacts)
          result.facts = std::move(*exactFacts);
        else
          result.facts = inferExpressionFacts(
              result.expression,
              &result.factorsPerTerm);
        attachExpressionHelperFacts(
            result.expression, result.facts);
        return result;
      }

    auto exactFacts =
        expressionFactsFromMetadata(expression);
    if (exactFacts)
      {
        result.expression = copy(expression);
        result.facts = std::move(*exactFacts);
        result.usedCompleteMetadataBypass = true;
        result.bypassedStraighteningFromMetadata =
            options.straightenIndices;
        result.bypassedSkewExpansionFromMetadata =
            options.expandSkewFactors;
        attachExpressionHelperFacts(
            result.expression, result.facts);
        return result;
      }

    const auto *inputPoly = polyValue(expression);
    const bool metadataProvesNormalized =
        inputPoly->conversionMetadata &&
        inputPoly->conversionMetadata->normalized;
    const bool metadataProvesSkewFree =
        inputPoly->conversionMetadata &&
        inputPoly->conversionMetadata->skewFree;
    std::vector<size_t> inputFactorsPerTerm;
    ExpressionFacts inputFacts =
        inferExpressionFacts(
            expression, &inputFactorsPerTerm);
    const bool indicesAlreadyStraightened =
        metadataProvesNormalized ||
        inputFacts.normalized;
    const bool alreadySkewFree =
        metadataProvesSkewFree ||
        inputFacts.skewFree;
    result.bypassedStraighteningFromMetadata =
        options.straightenIndices &&
        metadataProvesNormalized;
    result.bypassedSkewExpansionFromMetadata =
        options.expandSkewFactors &&
        metadataProvesSkewFree;
    if ((!options.straightenIndices ||
         indicesAlreadyStraightened) &&
        (!options.expandSkewFactors ||
         alreadySkewFree))
      {
        result.expression = copy(expression);
        result.facts = std::move(inputFacts);
        result.factorsPerTerm =
            std::move(inputFactorsPerTerm);
        result.usedAlreadyPreparedBypass = true;
        attachExpressionHelperFacts(
            result.expression, result.facts);
        return result;
      }

    result.expression = copy(expression);
    if (options.straightenIndices)
      {
        if (!indicesAlreadyStraightened)
          {
            result.expression =
                straighten(result.expression);
            result.performedStraightening = true;
            if (error())
              {
                result.expression = zero();
                return result;
              }
          }
      }
    if (options.expandSkewFactors)
      {
        if (!alreadySkewFree)
          {
            result.expression =
                expandSkewFactors(result.expression);
            result.performedSkewExpansion = true;
            if (error())
              {
                result.expression = zero();
                return result;
              }
          }
      }
    result.facts = inferExpressionFacts(
        result.expression,
        &result.factorsPerTerm);
    result.facts.combinatorialTags = inputTags;
    mutablePolyValue(
        result.expression)->combinatorialTags =
            inputTags;
    if (options.straightenIndices &&
        !result.facts.normalized)
      {
        ERROR("normalization did not establish its straightened-index contract");
        result.expression = zero();
        return result;
      }
    if (options.expandSkewFactors &&
        !result.facts.skewFree)
      {
        ERROR("normalization did not establish its skew-free contract");
        result.expression = zero();
        return result;
      }
    if (!result.facts.collected)
      {
        ERROR("normalization did not establish its collected-term contract");
        result.expression = zero();
        return result;
      }
    attachExpressionHelperFacts(
        result.expression, result.facts);
    return result;
  }

ring_elem SymmetricEngineRing::expandSkewFactors(
    ring_elem expression) const
{
    VECTOR(SymmetricTerm) expandedTerms;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      {
        ring_elem expanded = fromCoeff(term.coeff);
        size_t position = 0;
        while (position < term.monomial.data.size())
          {
            ring_elem factor =
                atomIsSkewAt(term.monomial, position)
                ? expandSkewBasisElement(
                      term.monomial, position)
                : expressionFromBasisElement(
                      term.monomial, position);
            if (error()) return zero();
            expanded = mult(expanded, factor);
            if (error()) return zero();
            position += atomLengthAt(
                term.monomial, position);
          }
        const auto *expandedPoly = polyValue(expanded);
        expandedTerms.insert(
            expandedTerms.end(),
            expandedPoly->terms.begin(),
            expandedPoly->terms.end());
      }
    ring_elem result = expressionHelperFromTerms(
        expandedTerms, poly->combinatorialTags);
    const ExpressionFacts facts =
        inferExpressionFacts(result);
    if (!facts.skewFree)
      {
        ERROR("skew expansion did not establish its skew-free contract");
        return zero();
      }
    attachExpressionHelperFacts(result, facts);
    return result;
  }

ring_elem SymmetricEngineRing::expandProductsInBasis(
    ring_elem expression,
    int targetBasisId) const
{
    (void) requireBasis(targetBasisId);
    if (error()) return zero();
    const ExpressionNormalizationOptions options =
        pipelinePreparationNormalizationOptions();
    ExpressionNormalizationResult prepared =
        normalizeExpressionWithOptionsDetailed(
            expression, options);
    if (error()) return zero();
    if (prepared.facts.noProducts())
      return prepared.expression;

    const auto& terms =
        polyValue(prepared.expression)->terms;
    if (prepared.factorsPerTerm.size() !=
        terms.size())
      {
        ERROR("product resolution received an inconsistent factor-count profile");
        return zero();
      }

    VECTOR(SymmetricTerm) productFreeTerms;
    for (size_t termPosition = 0;
         termPosition < terms.size();
         ++termPosition)
      {
        const auto& term = terms[termPosition];
        if (prepared.factorsPerTerm[termPosition] <= 1)
          {
            productFreeTerms.push_back(term);
            continue;
          }
        ResolvedProductTerm resolved =
            multiplyTermToBasis(
                term, targetBasisId);
        if (error()) return zero();
        for (const auto& resolvedTerm :
             polyValue(resolved.expression)->terms)
          appendTermIfNonZero(
              productFreeTerms,
              resolvedTerm.coeff,
              resolvedTerm.monomial);
      }

    ring_elem result =
        fromTermVector(productFreeTerms, false);
    auto *resultPoly = mutablePolyValue(result);
    resultPoly->combinatorialTags =
        prepared.facts.combinatorialTags;
    ExpressionFacts facts =
        inferExpressionFacts(result);
    facts.combinatorialTags =
        prepared.facts.combinatorialTags;
    if (!facts.normalized ||
        !facts.skewFree ||
        !facts.collected ||
        !facts.noProducts())
      {
        ERROR("product resolution did not establish its canonical-term contract");
        return zero();
      }
    attachExpressionHelperFacts(result, facts);
    return result;
  }

// ============================================================================
// Stable Expression Inspection
// ============================================================================

SymmetricEngineRing::ExpressionShape
SymmetricEngineRing::expressionShape(
    ring_elem expression) const
{
    const ExpressionFacts facts =
        inferExpressionFacts(expression);
    ExpressionShape result;
    result.weights = weightSupport(expression);
    result.basisIds = facts.factorBases;
    result.termCount = facts.termCount;
    result.scalarTermCount = facts.scalarTermCount;
    result.singleFactorTermCount =
        facts.singleFactorTermCount;
    result.productTermCount = facts.productTermCount;
    result.maximumFactorsPerTerm =
        facts.maximumFactorsPerTerm;
    result.skewFactorCount = facts.skewFactorCount;
    result.maximumPartitionLength =
        facts.maximumPartitionLength;
    result.normalized = facts.normalized;
    result.skewFree = facts.skewFree;
    result.collected = facts.collected;
    result.homogeneousWeight =
        facts.homogeneousWeight;
    result.pureBasis = facts.pureBasis;
    result.expandedBasis = facts.expandedBasis;
    const auto *poly = polyValue(expression);
    if (poly->conversionMetadata)
      {
        const auto& metadata =
            *poly->conversionMetadata;
        result.hasMetadata = true;
        result.metadataFactsComplete =
            metadata.expressionFactsComplete;
        result.metadataNormalized =
            metadata.normalized;
        result.metadataSkewFree =
            metadata.skewFree;
        result.metadataCollected =
            metadata.collected;
      }
    return result;
  }

std::vector<int> SymmetricEngineRing::basisSupport(
    ring_elem expression) const
{
    return inferExpressionFacts(expression).factorBases;
  }

bool SymmetricEngineRing::isBasisExpansion(
    ring_elem expression,
    int basisId) const
{
    (void) requireBasis(basisId);
    if (error()) return false;
    return inferExpressionFacts(expression)
        .canonicalExpansionInBasis(basisId);
  }

bool SymmetricEngineRing::isLinearCombinationOfBasisElements(
    ring_elem expression) const
{
    if (expressionFactsFromMetadata(expression))
      return true;
    const ExpressionFacts facts =
        inferExpressionFacts(expression);
    return facts.normalized &&
           facts.skewFree &&
           facts.collected &&
           facts.noProducts();
  }

ring_elem SymmetricEngineRing::coefficientsInBasisExpression(
    ring_elem expression,
    int basisId,
    bool convert) const
{
    (void) requireBasis(basisId);
    if (error()) return zero();
    if (convert) return toBasis(expression, basisId);
    if (!isBasisExpansion(expression, basisId))
      {
        ERROR("expected a canonical expansion in the requested basis");
        return zero();
      }
    return copy(expression);
  }

// ============================================================================
// Basis And Term Decomposition
// ============================================================================

std::vector<SymmetricEngineRing::BasisComponent>
SymmetricEngineRing::basisComponents(
    ring_elem expression) const
{
    std::map<int, VECTOR(SymmetricTerm)> termsByBasis;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      {
        int basisId = 0;
        if (!term.monomial.data.empty())
          {
            if (atomLengthAt(term.monomial, 0) !=
                term.monomial.data.size())
              {
                ERROR("basisComponents requires a product-free expression");
                return {};
              }
            basisId = atomBasisIdAt(term.monomial, 0);
          }
        termsByBasis[basisId].push_back(term);
      }

    std::vector<BasisComponent> result;
    result.reserve(termsByBasis.size());
    for (auto& item : termsByBasis)
      result.push_back({
          item.first,
          expressionHelperFromTerms(
              item.second, poly->combinatorialTags)});
    return result;
  }

std::vector<SymmetricEngineRing::HomogeneousBasisComponent>
SymmetricEngineRing::homogeneousBasisComponents(
    ring_elem expression) const
{
    std::map<
        std::pair<int, int>,
        VECTOR(SymmetricTerm)> groupedTerms;
    const auto *poly = polyValue(expression);
    for (const auto& term : poly->terms)
      {
        int basisId = 0;
        if (!term.monomial.data.empty())
          {
            if (atomLengthAt(term.monomial, 0) !=
                term.monomial.data.size())
              {
                ERROR(
                    "homogeneousBasisComponents requires "
                    "a product-free expression");
                return {};
              }
            basisId = atomBasisIdAt(term.monomial, 0);
          }
        groupedTerms[{
            monomialWeight(term.monomial),
            basisId}].push_back(term);
      }

    std::vector<HomogeneousBasisComponent> result;
    result.reserve(groupedTerms.size());
    for (auto& item : groupedTerms)
      result.push_back({
          item.first.first,
          item.first.second,
          expressionHelperFromTerms(
              item.second, poly->combinatorialTags)});
    return result;
  }

std::vector<ring_elem> SymmetricEngineRing::singlePartitionIndexedTerms(
    ring_elem expression) const
{
    if (!isLinearCombinationOfBasisElements(expression))
      {
        ERROR("singlePartitionIndexedTerms requires a straightened, skew-free, "
              "product-free expression");
        return {};
      }
    std::vector<ring_elem> result;
    const auto *poly = polyValue(expression);
    result.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        VECTOR(SymmetricTerm) oneTerm{term};
        result.push_back(
            expressionHelperFromTerms(
                oneTerm, poly->combinatorialTags));
      }
    return result;
  }

} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
