// Copyright 2026

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "basic-rings/aring-glue.hpp"
#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"
#include "rings/ZZ.hpp"
#include "rings/frac.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace symmetric_rings {

SymmetricMonomial SymmetricEngineRing::canonicalMonomial(const SymmetricMonomial& monomial) const
{
    std::vector<std::vector<int>> blocks;
    std::map<std::pair<int, int>, std::vector<int>> compressed;

    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        int order = atomOrderAt(monomial, pos);
        int basisId = atomBasisIdAt(monomial, pos);
        int n = atomIndexLengthAt(monomial, pos);
        if (!atomIsSkewAt(monomial, pos) && isMultiplicativeBasis(basisId))
          {
            auto& index = compressed[std::make_pair(order, basisId)];
            for (int i = 0; i < n; ++i)
              index.push_back(monomial.data[pos + atomHeaderSize + i]);
          }
        else
          {
            blocks.push_back(atomBlockAt(monomial, pos));
          }
        pos += atomHeaderSize + n;
      }

    for (auto& item : compressed)
      {
        auto& index = item.second;
        std::sort(index.begin(), index.end(), std::greater<int>());
        blocks.push_back(makeAtomBlock(item.first.first, item.first.second, 0, index));
      }

    std::sort(blocks.begin(), blocks.end(), blockLess);
    SymmetricMonomial result;
    size_t totalLength = 0;
    for (const auto& block : blocks) totalLength += block.size();
    result.data.reserve(totalLength);
    for (const auto& block : blocks) appendAtomBlock(result, block);
    return result;
  }

bool SymmetricEngineRing::isSinglePowerSumBlock(const SymmetricMonomial& monomial) const
{
    if (monomial.data.empty()) return false;
    if (atomIsSkewAt(monomial, 0)) return false;
    if (atomBasisIdAt(monomial, 0) != powerSumBasisId) return false;
    return atomLengthAt(monomial, 0) == monomial.data.size();
  }

SymmetricMonomial SymmetricEngineRing::multiplyPowerSumMonomials(const SymmetricMonomial& a,
                                             const SymmetricMonomial& b) const
{
    std::vector<int> index;
    int order = atomOrderAt(a, 0);
    int nA = atomIndexLengthAt(a, 0);
    int nB = atomIndexLengthAt(b, 0);
    index.reserve(nA + nB);
    for (int i = 0; i < nA; ++i) index.push_back(a.data[atomHeaderSize + i]);
    for (int i = 0; i < nB; ++i) index.push_back(b.data[atomHeaderSize + i]);
    std::sort(index.begin(), index.end(), std::greater<int>());

    SymmetricMonomial result;
    appendAtomBlock(result, makeAtomBlock(order, powerSumBasisId, 0, index));
    return result;
  }

SymmetricMonomial SymmetricEngineRing::multiplyMonomials(const SymmetricMonomial& a,
                                      const SymmetricMonomial& b) const
{
    if (a.data.empty()) return b;
    if (b.data.empty()) return a;
    if (isPowerSumBasis(powerSumBasisId) &&
        isSinglePowerSumBlock(a) &&
        isSinglePowerSumBlock(b))
      return multiplyPowerSumMonomials(a, b);

    SymmetricMonomial result;
    result.data.reserve(a.data.size() + b.data.size());
    result.data.insert(result.data.end(), a.data.begin(), a.data.end());
    result.data.insert(result.data.end(), b.data.begin(), b.data.end());
    return canonicalMonomial(result);
  }

ring_elem SymmetricEngineRing::copyPolyValue(const SymmetricRingPoly *poly) const
{
    auto result = new SymmetricRingPoly;
    result->conversionMetadata = poly->conversionMetadata;
    result->terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      result->terms.push_back({coefficientRing->copy(term.coeff), term.monomial});
    return makePolyValue(result);
  }

void SymmetricEngineRing::appendTermIfNonZero(VECTOR(SymmetricTerm)& terms,
                           ring_elem coeff,
                           const SymmetricMonomial& monomial) const
{
    if (!coefficientRing->is_zero(coeff)) terms.push_back({coeff, monomial});
  }

ring_elem SymmetricEngineRing::fromTermVector(VECTOR(SymmetricTerm)& terms, bool isSorted) const
{
    if (!isSorted)
      std::sort(terms.begin(), terms.end(),
                [](const SymmetricTerm& a, const SymmetricTerm& b) {
                  return compareMonomials(a.monomial, b.monomial) == LT;
                });

    auto result = new SymmetricRingPoly;
    result->terms.reserve(terms.size());
    size_t i = 0;
    while (i < terms.size())
      {
        ring_elem coeff = terms[i].coeff;
        size_t j = i + 1;
        while (j < terms.size() &&
               compareMonomials(terms[i].monomial, terms[j].monomial) == EQ)
          {
            coeff = coefficientRing->add(coeff, terms[j].coeff);
            ++j;
          }
        appendTermIfNonZero(result->terms, coeff, terms[i].monomial);
        i = j;
      }
    return makePolyValue(result);
  }

ring_elem SymmetricEngineRing::concatenateTerms(const SymmetricRingPoly *first,
                             const SymmetricRingPoly *second) const
{
    auto result = new SymmetricRingPoly;
    result->terms.reserve(first->terms.size() + second->terms.size());
    result->terms.insert(result->terms.end(), first->terms.begin(), first->terms.end());
    result->terms.insert(result->terms.end(), second->terms.begin(), second->terms.end());
    return makePolyValue(result);
  }

void SymmetricEngineRing::addToAccumulator(GCMap<std::vector<int>, ring_elem>& accumulator,
                        const SymmetricMonomial& monomial,
                        ring_elem coeff) const
{
    if (coefficientRing->is_zero(coeff)) return;
    auto key = monomialKey(monomial);
    auto existing = accumulator.find(key);
    if (existing == accumulator.end())
      {
        accumulator[key] = coeff;
        return;
      }
    ring_elem sum = coefficientRing->add(existing->second, coeff);
    if (coefficientRing->is_zero(sum))
      accumulator.erase(existing);
    else
      existing->second = sum;
  }

ring_elem SymmetricEngineRing::fromAccumulator(
      const GCMap<std::vector<int>, ring_elem>& accumulator) const
{
    auto result = new SymmetricRingPoly;
    result->terms.reserve(accumulator.size());
    for (const auto& term : accumulator)
      if (!coefficientRing->is_zero(term.second))
        result->terms.push_back({term.second, monomialFromKey(term.first)});
    return makePolyValue(result);
  }

bool SymmetricEngineRing::promoteInputElement(const RingElement *input, ring_elem &result) const
{
    const Ring *inputRing = input->get_ring();
    if (inputRing == this)
      {
        result = input->get_value();
        return true;
      }
    return promote(inputRing, input->get_value(), result);
  }

ring_elem SymmetricEngineRing::fromCoeff(ring_elem coeff) const
{
    auto result = new SymmetricRingPoly;
    if (!coefficientRing->is_zero(coeff)) result->terms.push_back({coeff, {}});
    return makePolyValue(result);
  }

ring_elem SymmetricEngineRing::basisElement(int basisId,
                         const std::string& display,
                         int order,
                         bool isMultiplicative,
                         int innerLength,
                         M2_arrayint index) const
{
    int payloadLength = index == nullptr ? 0 : index->len;
    if (innerLength < 0 || innerLength > payloadLength)
      {
        ERROR("invalid skew inner shape length");
        return zero();
      }
    rememberBasis(basisId, display, order, isMultiplicative);
    auto result = new SymmetricRingPoly;
    SymmetricMonomial monomial;
    appendAtomBlock(monomial, makeAtomBlock(order, basisId, innerLength, index));
    result->terms.push_back({coefficientRing->one(), canonicalMonomial(monomial)});
    return makePolyValue(result);
  }

bool SymmetricEngineRing::getScalar(const SymmetricRingPoly *f, ring_elem &result) const
{
    if (f->terms.empty())
      {
        result = coefficientRing->zero();
        return true;
      }
    if (f->terms.size() != 1) return false;
    if (!f->terms[0].monomial.data.empty()) return false;
    result = f->terms[0].coeff;
    return true;
  }

bool SymmetricEngineRing::hasPowerSumConversionHook(const SymmetricMonomial& monomial, size_t pos) const
{
    auto key = basisIndexKey(monomial, pos);
    return isPowerSumBasis(key.basisId);
  }

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

std::string SymmetricEngineRing::displayAtom(const SymmetricMonomial& monomial, size_t pos) const
{
    return displayForBasis(atomBasisIdAt(monomial, pos)) + "_" +
           displayIndex(monomial, pos);
  }

std::string SymmetricEngineRing::displayMonomial(const SymmetricMonomial& monomial) const
{
    if (monomial.data.empty()) return "1";
    std::vector<std::string> factors;
    size_t pos = 0;
    while (pos < monomial.data.size())
      {
        factors.push_back(displayAtom(monomial, pos));
        pos += atomLengthAt(monomial, pos);
      }
    return join(factors, "*");
  }

std::string SymmetricEngineRing::elementString(ring_elem f, int maxTerms) const
{
    const auto *poly = polyValue(f);
    if (poly->terms.empty()) return "0";
    size_t displayCount = poly->terms.size();
    if (maxTerms >= 0)
      displayCount = std::min(displayCount, static_cast<size_t>(maxTerms));
    std::vector<std::string> pieces;
    pieces.reserve(displayCount + 1);
    for (size_t i = 0; i < displayCount; ++i)
      {
        const auto& term = poly->terms[i];
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

int SymmetricEngineRing::elementWeight(ring_elem f) const
{
    const auto *poly = polyValue(f);
    bool found = false;
    int result = 0;
    for (const auto& term : poly->terms)
      {
        int wt = monomialWeight(term.monomial);
        if (!found)
          {
            result = wt;
            found = true;
          }
        else if (result != wt)
          return -1;
      }
    return found ? result : 0;
  }

unsigned int SymmetricEngineRing::computeHashValue(const ring_elem a) const
{
    size_t seed = 7919237;
    const auto *poly = polyValue(a);
    for (const auto& term : poly->terms)
      {
        hashCombine(seed, coefficientRing->computeHashValue(term.coeff));
        for (int value : term.monomial.data) hashCombine(seed, value);
      }
    return static_cast<unsigned int>(seed);
  }

void SymmetricEngineRing::text_out(buffer &o) const
{
    o << "SymmetricRing(";
    coefficientRing->text_out(o);
    o << ")";
  }

void SymmetricEngineRing::elem_text_out(buffer &o,
                             const ring_elem f,
                             bool p_one,
                             bool p_plus,
                             bool p_parens) const
{
    (void)p_one;
    if (p_plus && !is_zero(f)) o << '+';
    bool parens = p_parens && polyValue(f)->terms.size() > 1;
    if (parens) o << '(';
    o << elementString(f);
    if (parens) o << ')';
  }

ring_elem SymmetricEngineRing::from_long(long n) const
{
    return fromCoeff(coefficientRing->from_long(n));
  }

ring_elem SymmetricEngineRing::from_int(mpz_srcptr n) const
{
    return fromCoeff(coefficientRing->from_int(n));
  }

bool SymmetricEngineRing::from_rational(mpq_srcptr q, ring_elem &result) const
{
    ring_elem coeff;
    if (!coefficientRing->from_rational(q, coeff)) return false;
    result = fromCoeff(coeff);
    return true;
  }

bool SymmetricEngineRing::promote(const Ring *Rf, const ring_elem f, ring_elem &result) const
{
    if (Rf == this)
      {
        result = f;
        return true;
      }
    if (Rf == coefficientRing)
      {
        result = fromCoeff(f);
        return true;
      }
    if (Rf == globalZZ)
      {
        ring_elem coeff;
        if (!coefficientRing->promote(globalZZ, f, coeff)) return false;
        result = fromCoeff(coeff);
        return true;
      }
    const auto *Sf = dynamic_cast<const SymmetricEngineRing *>(Rf);
    if (Sf != nullptr)
      {
        rememberBasesFrom(Sf);
        GCMap<std::vector<int>, ring_elem> accumulator;
        const auto *poly = polyValue(f);
        for (const auto& term : poly->terms)
          {
            ring_elem coeff;
            if (!coefficientRing->promote(
                    Sf->getCoefficientRing(), term.coeff, coeff))
              return false;
            addToAccumulator(accumulator, term.monomial, coeff);
          }
        result = fromAccumulator(accumulator);
        return true;
      }
    return false;
  }

bool SymmetricEngineRing::lift(const Ring *Rg, const ring_elem f, ring_elem &result) const
{
    const auto *poly = polyValue(f);
    if (Rg == coefficientRing) return getScalar(poly, result);
    if (Rg == globalZZ)
      {
        ring_elem coeff;
        if (!getScalar(poly, coeff)) return false;
        return coefficientRing->lift(globalZZ, coeff, result);
      }
    const auto *Sg = dynamic_cast<const SymmetricEngineRing *>(Rg);
    if (Sg != nullptr)
      {
        Sg->rememberBasesFrom(this);
        GCMap<std::vector<int>, ring_elem> accumulator;
        for (const auto& term : poly->terms)
          {
            ring_elem coeff;
            if (!coefficientRing->lift(
                    Sg->getCoefficientRing(), term.coeff, coeff))
              return false;
            Sg->addToAccumulator(accumulator, term.monomial, coeff);
          }
        result = Sg->fromAccumulator(accumulator);
        return true;
      }
    return false;
  }

bool SymmetricEngineRing::is_unit(const ring_elem f) const
{
    ring_elem coeff;
    return getScalar(polyValue(f), coeff) && coefficientRing->is_unit(coeff);
  }

bool SymmetricEngineRing::is_zero(const ring_elem f) const
{
    return polyValue(f)->terms.empty();
  }

bool SymmetricEngineRing::is_equal(const ring_elem f, const ring_elem g) const
{
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    if (left->terms.size() != right->terms.size()) return false;
    for (size_t i = 0; i < left->terms.size(); ++i)
      {
        if (compareMonomials(left->terms[i].monomial, right->terms[i].monomial) != EQ)
          return false;
        if (!coefficientRing->is_equal(left->terms[i].coeff, right->terms[i].coeff))
          return false;
      }
    return true;
  }

int SymmetricEngineRing::compare_elems(const ring_elem f, const ring_elem g) const
{
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    size_t n = std::min(left->terms.size(), right->terms.size());
    for (size_t i = 0; i < n; ++i)
      {
        int monomialCmp = compareMonomials(left->terms[i].monomial, right->terms[i].monomial);
        if (monomialCmp != EQ) return monomialCmp;
        int coeffCmp =
            coefficientRing->compare_elems(left->terms[i].coeff, right->terms[i].coeff);
        if (coeffCmp != EQ) return coeffCmp;
      }
    if (left->terms.size() < right->terms.size()) return LT;
    if (left->terms.size() > right->terms.size()) return GT;
    return EQ;
  }

ring_elem SymmetricEngineRing::copy(const ring_elem f) const
{
    return copyPolyValue(polyValue(f));
  }

void SymmetricEngineRing::remove(ring_elem &f) const
{ (void)f; }

ring_elem SymmetricEngineRing::negate(const ring_elem f) const
{
    auto result = new SymmetricRingPoly;
    const auto *poly = polyValue(f);
    result->conversionMetadata = poly->conversionMetadata;
    result->terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      result->terms.push_back({coefficientRing->negate(term.coeff), term.monomial});
    return makePolyValue(result);
  }

ring_elem SymmetricEngineRing::add(const ring_elem f, const ring_elem g) const
{
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    auto preserveMetadata = [&](ring_elem value) {
      if (!left->conversionMetadata || !right->conversionMetadata) return value;
      const auto& leftMetadata = *left->conversionMetadata;
      const auto& rightMetadata = *right->conversionMetadata;
      SymmetricConversionMetadata metadata;
      if (leftMetadata.pureBasis == rightMetadata.pureBasis)
        metadata.pureBasis = leftMetadata.pureBasis;
      if (leftMetadata.expandedBasis == rightMetadata.expandedBasis)
        metadata.expandedBasis = leftMetadata.expandedBasis;
      if (leftMetadata.homogeneousWeight == rightMetadata.homogeneousWeight)
        metadata.homogeneousWeight = leftMetadata.homogeneousWeight;
      if (leftMetadata.maximumPartitionLength &&
          rightMetadata.maximumPartitionLength)
        metadata.maximumPartitionLength = std::max(
            *leftMetadata.maximumPartitionLength,
            *rightMetadata.maximumPartitionLength);
      if (leftMetadata.factorBases && rightMetadata.factorBases)
        {
          std::vector<int> factors = *leftMetadata.factorBases;
          factors.insert(factors.end(),
                         rightMetadata.factorBases->begin(),
                         rightMetadata.factorBases->end());
          std::sort(factors.begin(), factors.end());
          factors.erase(std::unique(factors.begin(), factors.end()), factors.end());
          metadata.factorBases = std::move(factors);
        }
      const auto *resultPoly = polyValue(value);
      metadata.termCount = resultPoly->terms.size();
      metadata.singleTerm = resultPoly->terms.size() == 1;
      metadata.singleAtom = resultPoly->terms.size() == 1 &&
                            !resultPoly->terms[0].monomial.data.empty();
      if (leftMetadata.noProducts && rightMetadata.noProducts)
        metadata.noProducts = *leftMetadata.noProducts &&
                              *rightMetadata.noProducts;
      metadata.normalized = leftMetadata.normalized && rightMetadata.normalized;
      metadata.skewFree = leftMetadata.skewFree && rightMetadata.skewFree;
      metadata.collected = true;
      mutablePolyValue(value)->conversionMetadata = std::move(metadata);
      return value;
    };
    if (left->terms.empty()) return copyPolyValue(right);
    if (right->terms.empty()) return copyPolyValue(left);

    if (compareMonomials(left->terms.back().monomial,
                         right->terms.front().monomial) == LT)
      return preserveMetadata(concatenateTerms(left, right));
    if (compareMonomials(right->terms.back().monomial,
                         left->terms.front().monomial) == LT)
      return preserveMetadata(concatenateTerms(right, left));

    auto result = new SymmetricRingPoly;
    result->terms.reserve(left->terms.size() + right->terms.size());
    size_t i = 0;
    size_t j = 0;
    while (i < left->terms.size() && j < right->terms.size())
      {
        int cmp = compareMonomials(left->terms[i].monomial,
                                   right->terms[j].monomial);
        if (cmp == LT)
          result->terms.push_back(left->terms[i++]);
        else if (cmp == GT)
          result->terms.push_back(right->terms[j++]);
        else
          {
            ring_elem coeff = coefficientRing->add(left->terms[i].coeff,
                                                   right->terms[j].coeff);
            appendTermIfNonZero(result->terms, coeff, left->terms[i].monomial);
            ++i;
            ++j;
          }
      }
    result->terms.insert(result->terms.end(), left->terms.begin() + i, left->terms.end());
    result->terms.insert(result->terms.end(), right->terms.begin() + j, right->terms.end());
    return preserveMetadata(makePolyValue(result));
  }

ring_elem SymmetricEngineRing::subtract(const ring_elem f, const ring_elem g) const
{
    return add(f, negate(g));
  }

SymmetricRingPoly *SymmetricEngineRing::multByCoefficient(ring_elem coeff,
                                       const SymmetricRingPoly *poly) const
{
    auto result = new SymmetricRingPoly;
    if (coefficientRing->is_zero(coeff)) return result;
    result->conversionMetadata = poly->conversionMetadata;
    result->terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        ring_elem c = coefficientRing->mult(coeff, term.coeff);
        if (!coefficientRing->is_zero(c)) result->terms.push_back({c, term.monomial});
      }
    return result;
  }

ring_elem SymmetricEngineRing::mult(const ring_elem f, const ring_elem g) const
{
    ring_elem scalar;
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    auto preserveProductMetadata = [&](ring_elem value) {
      if (!left->conversionMetadata || !right->conversionMetadata) return value;
      const auto& leftMetadata = *left->conversionMetadata;
      const auto& rightMetadata = *right->conversionMetadata;
      SymmetricConversionMetadata metadata;
      if (leftMetadata.pureBasis &&
          leftMetadata.pureBasis == rightMetadata.pureBasis)
        metadata.pureBasis = leftMetadata.pureBasis;
      if (leftMetadata.homogeneousWeight && rightMetadata.homogeneousWeight)
        metadata.homogeneousWeight = *leftMetadata.homogeneousWeight +
                                     *rightMetadata.homogeneousWeight;
      if (leftMetadata.maximumPartitionLength &&
          rightMetadata.maximumPartitionLength)
        metadata.maximumPartitionLength = std::max(
            *leftMetadata.maximumPartitionLength,
            *rightMetadata.maximumPartitionLength);
      if (leftMetadata.factorBases && rightMetadata.factorBases)
        {
          std::vector<int> factors = *leftMetadata.factorBases;
          factors.insert(factors.end(),
                         rightMetadata.factorBases->begin(),
                         rightMetadata.factorBases->end());
          std::sort(factors.begin(), factors.end());
          factors.erase(std::unique(factors.begin(), factors.end()), factors.end());
          metadata.factorBases = std::move(factors);
        }
      metadata.normalized = leftMetadata.normalized && rightMetadata.normalized;
      metadata.skewFree = leftMetadata.skewFree && rightMetadata.skewFree;
      if (metadata.pureBasis && isMultiplicativeBasis(*metadata.pureBasis) &&
          metadata.skewFree)
        {
          metadata.expandedBasis = metadata.pureBasis;
          metadata.noProducts = true;
        }
      const auto *resultPoly = polyValue(value);
      metadata.termCount = resultPoly->terms.size();
      metadata.singleTerm = resultPoly->terms.size() == 1;
      metadata.singleAtom = resultPoly->terms.size() == 1 &&
                            !resultPoly->terms[0].monomial.data.empty();
      metadata.collected = true;
      mutablePolyValue(value)->conversionMetadata = std::move(metadata);
      return value;
    };
    if (getScalar(left, scalar)) return makePolyValue(multByCoefficient(scalar, right));
    if (getScalar(right, scalar)) return makePolyValue(multByCoefficient(scalar, left));
    if (left->terms.size() == 1 && right->terms.size() == 1)
      {
        ring_elem coeff = coefficientRing->mult(left->terms[0].coeff,
                                                right->terms[0].coeff);
        auto result = new SymmetricRingPoly;
        if (!coefficientRing->is_zero(coeff))
          result->terms.push_back({coeff,
                                   multiplyMonomials(left->terms[0].monomial,
                                                     right->terms[0].monomial)});
        return preserveProductMetadata(makePolyValue(result));
      }
    VECTOR(SymmetricTerm) products;
    products.reserve(left->terms.size() * right->terms.size());
    for (const auto& lt : left->terms)
      for (const auto& rt : right->terms)
        {
          ring_elem coeff = coefficientRing->mult(lt.coeff, rt.coeff);
          if (!coefficientRing->is_zero(coeff))
            products.push_back(
                {coeff, multiplyMonomials(lt.monomial, rt.monomial)});
        }
    return preserveProductMetadata(fromTermVector(products, false));
  }

bool SymmetricEngineRing::promoteCollectedExpansion(
    const RingElement *input,
    ring_elem& result) const
{
    const auto *sourceRing =
        dynamic_cast<const SymmetricEngineRing *>(input->get_ring());
    if (sourceRing == nullptr) return false;
    if (sourceRing == this)
      {
        result = copyPolyValue(polyValue(input->get_value()));
        return true;
      }

    rememberBasesFrom(sourceRing);
    const auto *source = polyValue(input->get_value());
    auto target = new SymmetricRingPoly;
    // Coefficient promotion cannot disturb an already collected monomial order.
    target->terms.reserve(source->terms.size());
    target->conversionMetadata = source->conversionMetadata;
    for (const auto& term : source->terms)
      {
        ring_elem coeff;
        const Ring *sourceCoefficients = sourceRing->getCoefficientRing();
        bool promoted = sourceCoefficients == globalQQ
            ? coefficientRing->from_rational(term.coeff.get_mpq(), coeff)
            : coefficientRing->promote(sourceCoefficients, term.coeff, coeff);
        if (!promoted) return false;
        if (!coefficientRing->is_zero(coeff))
          target->terms.push_back({coeff, term.monomial});
      }
    result = makePolyValue(target);
    return true;
  }

bool SymmetricEngineRing::liftCollectedExpansion(
    const RingElement *input,
    ring_elem& result) const
{
    const auto *sourceRing =
        dynamic_cast<const SymmetricEngineRing *>(input->get_ring());
    if (sourceRing == nullptr) return false;
    if (sourceRing == this)
      {
        result = copyPolyValue(polyValue(input->get_value()));
        return true;
      }

    rememberBasesFrom(sourceRing);
    const auto *source = polyValue(input->get_value());
    const Ring *sourceCoefficients = sourceRing->getCoefficientRing();
    auto target = new SymmetricRingPoly;
    target->terms.reserve(source->terms.size());
    target->conversionMetadata = source->conversionMetadata;
    for (const auto& term : source->terms)
      {
        ring_elem coeff;
        bool lifted = sourceCoefficients == coefficientRing
            ? (coeff = coefficientRing->copy(term.coeff), true)
            : sourceCoefficients->lift(coefficientRing, term.coeff, coeff);
        if (!lifted && coefficientRing == globalQQ)
          {
            const auto *fractionSource =
                dynamic_cast<const FractionField *>(sourceCoefficients);
            if (fractionSource != nullptr)
              {
                const Ring *baseRing = fractionSource->get_ring();
                ring_elem numerator = fractionSource->numerator(term.coeff);
                ring_elem denominator = fractionSource->denominator(term.coeff);
                ring_elem numeratorQQ;
                ring_elem denominatorQQ;
                if (baseRing->lift(globalQQ, numerator, numeratorQQ) &&
                    baseRing->lift(globalQQ, denominator, denominatorQQ))
                  {
                    coeff = globalQQ->divide(numeratorQQ, denominatorQQ);
                    lifted = true;
                  }
              }
          }
        if (!lifted) return false;
        if (!coefficientRing->is_zero(coeff))
          target->terms.push_back({coeff, term.monomial});
      }
    result = makePolyValue(target);
    return true;
  }

ring_elem SymmetricEngineRing::batchSum(engine_RawRingElementArray elements) const
{
    if (elements == nullptr || elements->len == 0) return zero();

    size_t termCount = 0;
    for (int i = 0; i < elements->len; ++i)
      {
        ring_elem promoted;
        if (!promoteInputElement(elements->array[i], promoted)) return zero();
        termCount += polyValue(promoted)->terms.size();
      }

    VECTOR(SymmetricTerm) terms;
    terms.reserve(termCount);
    bool sorted = true;
    bool havePrevious = false;
    SymmetricMonomial previous;

    for (int i = 0; i < elements->len; ++i)
      {
        ring_elem promoted;
        if (!promoteInputElement(elements->array[i], promoted)) return zero();
        const auto *poly = polyValue(promoted);
        for (const auto& term : poly->terms)
          {
            if (havePrevious && compareMonomials(previous, term.monomial) == GT)
              sorted = false;
            previous = term.monomial;
            havePrevious = true;
            terms.push_back(term);
          }
      }
    return fromTermVector(terms, sorted);
  }

ring_elem SymmetricEngineRing::batchProduct(engine_RawRingElementArray elements) const
{
    if (elements == nullptr || elements->len == 0) return one();

    ring_elem coeff = coefficientRing->one();
    SymmetricMonomial monomial;
    bool singleTermProduct = true;

    for (int i = 0; i < elements->len; ++i)
      {
        ring_elem promoted;
        if (!promoteInputElement(elements->array[i], promoted)) return zero();
        const auto *poly = polyValue(promoted);
        if (poly->terms.empty()) return zero();
        if (poly->terms.size() != 1)
          {
            singleTermProduct = false;
            break;
          }
        coeff = coefficientRing->mult(coeff, poly->terms[0].coeff);
        if (coefficientRing->is_zero(coeff)) return zero();
        monomial = multiplyMonomials(monomial, poly->terms[0].monomial);
      }

    if (singleTermProduct)
      {
        auto result = new SymmetricRingPoly;
        result->terms.push_back({coeff, monomial});
        return makePolyValue(result);
      }

    ring_elem result = one();
    for (int i = 0; i < elements->len; ++i)
      {
        ring_elem promoted;
        if (!promoteInputElement(elements->array[i], promoted)) return zero();
        result = mult(result, promoted);
      }
    return result;
  }

ring_elem SymmetricEngineRing::invert(const ring_elem f) const
{
    (void)f;
    return zero();
  }

ring_elem SymmetricEngineRing::divide(const ring_elem f, const ring_elem g) const
{
    (void)f;
    (void)g;
    return zero();
  }

void SymmetricEngineRing::syzygy(const ring_elem a,
                      const ring_elem b,
                      ring_elem &x,
                      ring_elem &y) const
{
    (void)a;
    (void)b;
    x = zero();
    y = zero();
  }

ring_elem SymmetricEngineRing::eval(const RingMap *map,
                         const ring_elem f,
                         int first_var) const
{
    (void)f;
    (void)first_var;
    return map->get_ring()->zero();
  }

} // namespace symmetric_rings
