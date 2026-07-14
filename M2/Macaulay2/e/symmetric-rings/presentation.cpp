// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Basis-Element And Monomial Presentation
// ============================================================================
// Presentation order is independent of canonical storage order. These helpers
// expose the stable user-facing order without rewriting stored terms.

std::string SymmetricEngineRing::displayIndex(const SymmetricMonomial& monomial, size_t pos) const
{
    int n = atomIndexLengthAt(monomial, pos);
    if (n <= 0) return "{}";
    if (atomIsSkewAt(monomial, pos))
      {
        int outerLength = atomOuterLengthAt(monomial, pos);
        int innerLength = atomInnerLengthAt(monomial, pos);
        if (outerLength + innerLength != n) return "{}/{}";
        std::vector<std::string> lambdaParts;
        std::vector<std::string> muParts;
        lambdaParts.reserve(outerLength);
        muParts.reserve(innerLength);
        for (int i = 0; i < outerLength; ++i)
          lambdaParts.push_back(std::to_string(monomial.data[pos + atomHeaderSize + i]));
        for (int i = outerLength; i < n; ++i)
          muParts.push_back(std::to_string(monomial.data[pos + atomHeaderSize + i]));
        return "{{" + join(lambdaParts, ",") + "}/{" + join(muParts, ",") + "}}";
      }
    if (n == 1) return std::to_string(monomial.data[pos + atomHeaderSize]);
    std::vector<std::string> parts;
    parts.reserve(n);
    for (int i = 0; i < n; ++i)
      parts.push_back(std::to_string(monomial.data[pos + atomHeaderSize + i]));
    return "{" + join(parts, ",") + "}";
  }

std::string SymmetricEngineRing::displayBasisElement(
    const SymmetricMonomial& monomial,
    size_t pos) const
{
    return displayForBasis(atomBasisIdAt(monomial, pos)) + "_" +
           displayIndex(monomial, pos);
  }

std::string SymmetricEngineRing::displayMonomial(const SymmetricMonomial& monomial) const
{
    if (monomial.data.empty()) return "1";
    std::vector<std::string> factors;
    for (size_t pos : presentationAtomPositions(monomial))
        factors.push_back(displayBasisElement(monomial, pos));
    return join(factors, "*");
  }

namespace {

int compareIndexSegmentsDescending(const SymmetricMonomial& a,
                                   size_t aStart,
                                   int aLength,
                                   const SymmetricMonomial& b,
                                   size_t bStart,
                                   int bLength)
{
    int common = std::min(aLength, bLength);
    for (int i = 0; i < common; ++i)
      {
        int av = a.data[aStart + static_cast<size_t>(i)];
        int bv = b.data[bStart + static_cast<size_t>(i)];
        if (av != bv) return av > bv ? LT : GT;
      }
    if (aLength != bLength) return aLength > bLength ? LT : GT;
    return EQ;
}

struct PresentationTermKey
{
  size_t storageIndex;
  int weight;
  std::vector<size_t> atomPositions;
  std::vector<int> basisSignature;
};

} // namespace

int SymmetricEngineRing::compareBasisIdsForPresentation(int aBasis,
                                                        int bBasis) const
{
    int aOrder = basisOrderForId(aBasis);
    int bOrder = basisOrderForId(bBasis);
    if (aOrder != bOrder) return aOrder > bOrder ? LT : GT;
    std::string aKey = basisKeyForId(aBasis);
    std::string bKey = basisKeyForId(bBasis);
    if (aKey != bKey) return aKey < bKey ? LT : GT;
    if (aBasis != bBasis) return aBasis < bBasis ? LT : GT;
    return EQ;
  }

std::vector<size_t> SymmetricEngineRing::presentationAtomPositions(
    const SymmetricMonomial& monomial) const
{
    std::vector<size_t> positions;
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        positions.push_back(pos);
        pos += atomLengthAt(monomial, pos);
      }

    auto compareAtoms = [&](size_t aPos, size_t bPos) {
      int aBasis = atomBasisIdAt(monomial, aPos);
      int bBasis = atomBasisIdAt(monomial, bPos);
      int basisCmp = compareBasisIdsForPresentation(aBasis, bBasis);
      if (basisCmp != EQ) return basisCmp == LT;

      int outerCmp = compareIndexSegmentsDescending(
          monomial,
          aPos + atomHeaderSize,
          atomOuterLengthAt(monomial, aPos),
          monomial,
          bPos + atomHeaderSize,
          atomOuterLengthAt(monomial, bPos));
      if (outerCmp != EQ) return outerCmp == LT;
      int innerCmp = compareIndexSegmentsDescending(
          monomial,
          aPos + atomHeaderSize + atomOuterLengthAt(monomial, aPos),
          atomInnerLengthAt(monomial, aPos),
          monomial,
          bPos + atomHeaderSize + atomOuterLengthAt(monomial, bPos),
          atomInnerLengthAt(monomial, bPos));
      return innerCmp == LT;
    };
    std::sort(positions.begin(), positions.end(), compareAtoms);
    return positions;
  }

std::vector<int> SymmetricEngineRing::presentationMonomialData(
    const SymmetricMonomial& monomial) const
{
    std::vector<int> result;
    result.reserve(monomial.data.size());
    for (size_t pos : presentationAtomPositions(monomial))
      {
        size_t length = atomLengthAt(monomial, pos);
        result.insert(result.end(),
                      monomial.data.begin() + pos,
                      monomial.data.begin() + pos + length);
      }
    return result;
  }

std::vector<size_t> SymmetricEngineRing::presentationTermOrder(ring_elem f) const
{
    const auto *poly = polyValue(f);
    std::vector<PresentationTermKey> keys;
    keys.reserve(poly->terms.size());
    const bool homogeneous = poly->conversionMetadata &&
                             poly->conversionMetadata->homogeneousWeight;
    const int homogeneousWeight = homogeneous
        ? *poly->conversionMetadata->homogeneousWeight
        : 0;

    for (size_t i = 0; i < poly->terms.size(); ++i)
      {
        const auto& monomial = poly->terms[i].monomial;
        PresentationTermKey key;
        key.storageIndex = i;
        key.weight = homogeneous ? homogeneousWeight : monomialWeight(monomial);
        key.atomPositions = presentationAtomPositions(monomial);
        for (size_t atomPos : key.atomPositions)
          {
            int basisId = atomBasisIdAt(monomial, atomPos);
            if (std::find(key.basisSignature.begin(), key.basisSignature.end(), basisId) ==
                key.basisSignature.end())
              key.basisSignature.push_back(basisId);
          }
        std::sort(key.basisSignature.begin(), key.basisSignature.end(),
                  [&](int aBasis, int bBasis) {
                    return compareBasisIdsForPresentation(aBasis, bBasis) == LT;
                  });
        keys.push_back(std::move(key));
      }

    auto compareAtoms = [&](const SymmetricMonomial& a,
                            size_t aPos,
                            const SymmetricMonomial& b,
                            size_t bPos) {
      int basisCmp = compareBasisIdsForPresentation(
          atomBasisIdAt(a, aPos), atomBasisIdAt(b, bPos));
      if (basisCmp != EQ) return basisCmp;
      int outerCmp = compareIndexSegmentsDescending(
          a, aPos + atomHeaderSize, atomOuterLengthAt(a, aPos),
          b, bPos + atomHeaderSize, atomOuterLengthAt(b, bPos));
      if (outerCmp != EQ) return outerCmp;
      return compareIndexSegmentsDescending(
          a,
          aPos + atomHeaderSize + atomOuterLengthAt(a, aPos),
          atomInnerLengthAt(a, aPos),
          b,
          bPos + atomHeaderSize + atomOuterLengthAt(b, bPos),
          atomInnerLengthAt(b, bPos));
    };

    auto presentationLess = [&](const PresentationTermKey& a,
                                const PresentationTermKey& b) {
      if (a.weight != b.weight) return a.weight > b.weight;
      bool aScalar = a.atomPositions.empty();
      bool bScalar = b.atomPositions.empty();
      if (aScalar != bScalar) return !aScalar;
      if (!aScalar)
        {
          int leadingBasisCmp = compareBasisIdsForPresentation(
              a.basisSignature.front(), b.basisSignature.front());
          if (leadingBasisCmp != EQ) return leadingBasisCmp == LT;
        }
      if (a.basisSignature.size() != b.basisSignature.size())
        return a.basisSignature.size() < b.basisSignature.size();
      for (size_t i = 0; i < a.basisSignature.size(); ++i)
        {
          int cmp = compareBasisIdsForPresentation(
              a.basisSignature[i], b.basisSignature[i]);
          if (cmp != EQ) return cmp == LT;
        }

      const auto& aMonomial = poly->terms[a.storageIndex].monomial;
      const auto& bMonomial = poly->terms[b.storageIndex].monomial;
      size_t common = std::min(a.atomPositions.size(), b.atomPositions.size());
      for (size_t i = 0; i < common; ++i)
        {
          int cmp = compareAtoms(aMonomial, a.atomPositions[i],
                                 bMonomial, b.atomPositions[i]);
          if (cmp != EQ) return cmp == LT;
        }
      if (a.atomPositions.size() != b.atomPositions.size())
        return a.atomPositions.size() < b.atomPositions.size();
      int storageCmp = compareMonomials(aMonomial, bMonomial);
      if (storageCmp != EQ) return storageCmp == LT;
      return a.storageIndex < b.storageIndex;
    };
    std::sort(keys.begin(), keys.end(), presentationLess);

    std::vector<size_t> result;
    result.reserve(keys.size());
    for (const auto& key : keys) result.push_back(key.storageIndex);
    return result;
  }

std::string SymmetricEngineRing::elementString(ring_elem f, int maxTerms) const
{
    const auto *poly = polyValue(f);
    if (poly->terms.empty()) return "0";
    std::vector<size_t> order = presentationTermOrder(f);
    size_t displayCount = order.size();
    if (maxTerms >= 0)
      displayCount = std::min(displayCount, static_cast<size_t>(maxTerms));
    std::vector<std::string> pieces;
    pieces.reserve(displayCount + 1);
    for (size_t i = 0; i < displayCount; ++i)
      {
        const auto& term = poly->terms[order[i]];
        if (term.monomial.data.empty())
          {
            pieces.push_back(coeffToString(coefficientRing, term.coeff));
          }
        else if (coefficientRing->is_equal(term.coeff, coefficientRing->one()))
          {
            pieces.push_back(displayMonomial(term.monomial));
          }
        else if (coefficientRing->is_equal(term.coeff, coefficientRing->minus_one()))
          {
            pieces.push_back("-" + displayMonomial(term.monomial));
          }
        else
          {
            pieces.push_back("(" + coeffToString(coefficientRing, term.coeff) +
                             ")*" + displayMonomial(term.monomial));
          }
      }
    if (displayCount < poly->terms.size())
      pieces.push_back(std::to_string(poly->terms.size() - displayCount) + " terms");
    return join(pieces, " + ");
  }


} // namespace symmetric_rings

// Local Variables:
// indent-tabs-mode: nil
// End:
