// Copyright 2026

#include "interface/symmetric-rings.h"

#include "buffer.hpp"
#include "error.h"
#include "exceptions.hpp"
#include "newdelete.hpp"
#include "ring-elements/ring-element.hpp"
#include "ringmap.hpp"
#include "rings/ZZ.hpp"
#include "rings/ring.hpp"
#include "rings/ringelem.hpp"

#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct SymmetricAtom
{
  // [displayOrder, basisId, indexLength, index_1, ..., index_n]
  VECTOR(int) data;
};

struct SymmetricMonomial
{
  VECTOR(SymmetricAtom) atoms;
};

struct SymmetricTerm
{
  ring_elem coeff;
  SymmetricMonomial monomial;
};

class SymmetricRingPoly : public our_new_delete
{
 public:
  VECTOR(SymmetricTerm) terms;
};

std::string fromM2String(M2_string s)
{
  if (s == nullptr) return "";
  return std::string(reinterpret_cast<const char *>(s->array), s->len);
}

M2_string toM2String(const std::string& s)
{
  return M2_tostring(s.c_str());
}

std::string join(const std::vector<std::string>& parts,
                 const std::string& delimiter)
{
  std::ostringstream out;
  for (size_t i = 0; i < parts.size(); ++i)
    {
      if (i > 0) out << delimiter;
      out << parts[i];
    }
  return out.str();
}

const SymmetricRingPoly *polyValue(ring_elem f)
{
  return reinterpret_cast<const SymmetricRingPoly *>(f.get_Poly());
}

ring_elem makePolyValue(SymmetricRingPoly *f)
{
  return ring_elem(reinterpret_cast<const void *>(f));
}

template <typename T>
void hashCombine(size_t& seed, const T& val)
{
  seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

std::string coeffToString(const Ring *R, ring_elem c)
{
  buffer o;
  R->elem_text_out(o, c);
  return std::string(o.str());
}

SymmetricAtom makeAtom(int displayOrder, int basisId, M2_arrayint index)
{
  SymmetricAtom atom;
  atom.data.push_back(displayOrder);
  atom.data.push_back(basisId);
  atom.data.push_back(index == nullptr ? 0 : index->len);
  if (index != nullptr)
    for (int i = 0; i < index->len; ++i) atom.data.push_back(index->array[i]);
  return atom;
}

bool atomLess(const SymmetricAtom& a, const SymmetricAtom& b)
{
  return std::lexicographical_compare(
      a.data.begin(), a.data.end(), b.data.begin(), b.data.end());
}

SymmetricMonomial canonicalMonomial(SymmetricMonomial monomial)
{
  std::sort(monomial.atoms.begin(), monomial.atoms.end(), atomLess);
  return monomial;
}

SymmetricMonomial mergeMonomials(const SymmetricMonomial& a,
                                 const SymmetricMonomial& b)
{
  SymmetricMonomial result;
  result.atoms.reserve(a.atoms.size() + b.atoms.size());
  size_t i = 0;
  size_t j = 0;
  while (i < a.atoms.size() && j < b.atoms.size())
    {
      if (atomLess(b.atoms[j], a.atoms[i]))
        result.atoms.push_back(b.atoms[j++]);
      else
        result.atoms.push_back(a.atoms[i++]);
    }
  result.atoms.insert(result.atoms.end(), a.atoms.begin() + i, a.atoms.end());
  result.atoms.insert(result.atoms.end(), b.atoms.begin() + j, b.atoms.end());
  return result;
}

std::vector<int> monomialKey(const SymmetricMonomial& monomial)
{
  std::vector<int> result;
  for (const auto& atom : monomial.atoms)
    {
      result.insert(result.end(), atom.data.begin(), atom.data.end());
    }
  return result;
}

SymmetricMonomial monomialFromKey(const std::vector<int>& key)
{
  SymmetricMonomial result;
  size_t pos = 0;
  while (pos < key.size())
    {
      if (pos + 2 >= key.size()) break;
      int blockLen = 3 + key[pos + 2];
      SymmetricAtom atom;
      for (int j = 0; j < blockLen && pos < key.size(); ++j)
        atom.data.push_back(key[pos++]);
      result.atoms.push_back(atom);
    }
  return result;
}

SymmetricMonomial multiplyMonomials(const SymmetricMonomial& a,
                                    const SymmetricMonomial& b)
{
  if (a.atoms.empty()) return b;
  if (b.atoms.empty()) return a;
  if (!atomLess(b.atoms.front(), a.atoms.back()))
    {
      SymmetricMonomial result;
      result.atoms.reserve(a.atoms.size() + b.atoms.size());
      result.atoms.insert(result.atoms.end(), a.atoms.begin(), a.atoms.end());
      result.atoms.insert(result.atoms.end(), b.atoms.begin(), b.atoms.end());
      return result;
    }
  if (atomLess(b.atoms.back(), a.atoms.front()))
    {
      SymmetricMonomial result;
      result.atoms.reserve(a.atoms.size() + b.atoms.size());
      result.atoms.insert(result.atoms.end(), b.atoms.begin(), b.atoms.end());
      result.atoms.insert(result.atoms.end(), a.atoms.begin(), a.atoms.end());
      return result;
    }
  return mergeMonomials(a, b);
}

int compareMonomials(const SymmetricMonomial& a, const SymmetricMonomial& b)
{
  size_t atomA = 0;
  size_t atomB = 0;
  size_t posA = 0;
  size_t posB = 0;
  while (atomA < a.atoms.size() && atomB < b.atoms.size())
    {
      const auto& dataA = a.atoms[atomA].data;
      const auto& dataB = b.atoms[atomB].data;
      while (posA < dataA.size() && posB < dataB.size())
        {
          if (dataA[posA] < dataB[posB]) return LT;
          if (dataB[posB] < dataA[posA]) return GT;
          ++posA;
          ++posB;
        }
      if (posA == dataA.size())
        {
          ++atomA;
          posA = 0;
        }
      if (posB == dataB.size())
        {
          ++atomB;
          posB = 0;
        }
    }
  if (atomA < a.atoms.size()) return GT;
  if (atomB < b.atoms.size()) return LT;
  return EQ;
}

int atomIndexLength(const SymmetricAtom& atom)
{
  return atom.data.size() >= 3 ? atom.data[2] : 0;
}

int atomBasisId(const SymmetricAtom& atom)
{
  return atom.data.size() >= 2 ? atom.data[1] : -1;
}

int atomIndexValue(const SymmetricAtom& atom, int i)
{
  return atom.data[3 + i];
}

std::string displayIndex(const SymmetricAtom& atom)
{
  int n = atomIndexLength(atom);
  if (n <= 0) return "{}";
  if (n == 1) return std::to_string(atomIndexValue(atom, 0));
  std::vector<std::string> parts;
  for (int i = 0; i < n; ++i) parts.push_back(std::to_string(atomIndexValue(atom, i)));
  return "{" + join(parts, ",") + "}";
}

int atomWeight(const SymmetricAtom& atom)
{
  int result = 0;
  for (int i = 0; i < atomIndexLength(atom); ++i) result += atomIndexValue(atom, i);
  return result;
}

class SymmetricEngineRing : public Ring
{
 private:
  const Ring *coefficientRing;
  mutable std::map<int, std::string> basisDisplays;

  void rememberBasis(int basisId, const std::string& display) const
  {
    if (!display.empty()) basisDisplays[basisId] = display;
  }

  void rememberBasesFrom(const SymmetricEngineRing *R) const
  {
    for (const auto& item : R->basisDisplays) basisDisplays[item.first] = item.second;
  }

  std::string displayForBasis(int basisId) const
  {
    auto it = basisDisplays.find(basisId);
    if (it != basisDisplays.end()) return it->second;
    return "basis" + std::to_string(basisId);
  }

  ring_elem copyPolyValue(const SymmetricRingPoly *poly) const
  {
    auto result = new SymmetricRingPoly;
    result->terms.reserve(poly->terms.size());
    result->terms.insert(result->terms.end(), poly->terms.begin(), poly->terms.end());
    return makePolyValue(result);
  }

  void appendTermIfNonZero(VECTOR(SymmetricTerm)& terms,
                           ring_elem coeff,
                           const SymmetricMonomial& monomial) const
  {
    if (!coefficientRing->is_zero(coeff)) terms.push_back({coeff, monomial});
  }

  ring_elem fromTermVector(VECTOR(SymmetricTerm)& terms, bool isSorted) const
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

  ring_elem concatenateTerms(const SymmetricRingPoly *first,
                             const SymmetricRingPoly *second) const
  {
    auto result = new SymmetricRingPoly;
    result->terms.reserve(first->terms.size() + second->terms.size());
    result->terms.insert(result->terms.end(), first->terms.begin(), first->terms.end());
    result->terms.insert(result->terms.end(), second->terms.begin(), second->terms.end());
    return makePolyValue(result);
  }

  void addToAccumulator(std::map<std::vector<int>, ring_elem>& accumulator,
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

  ring_elem fromAccumulator(
      const std::map<std::vector<int>, ring_elem>& accumulator) const
  {
    auto result = new SymmetricRingPoly;
    result->terms.reserve(accumulator.size());
    for (const auto& term : accumulator)
      if (!coefficientRing->is_zero(term.second))
        result->terms.push_back({term.second, monomialFromKey(term.first)});
    return makePolyValue(result);
  }

  bool promoteInputElement(const RingElement *input, ring_elem &result) const
  {
    const Ring *inputRing = input->get_ring();
    if (inputRing == this)
      {
        result = input->get_value();
        return true;
      }
    return promote(inputRing, input->get_value(), result);
  }

 public:
  explicit SymmetricEngineRing(const Ring *A) : coefficientRing(A) {}

  static SymmetricEngineRing *create(const Ring *A)
  {
    auto result = new SymmetricEngineRing(A);
    result->initialize_ring(A->characteristic());
    result->zeroV = result->from_long(0);
    result->oneV = result->from_long(1);
    result->minus_oneV = result->from_long(-1);
    return result;
  }

  const Ring *getCoefficientRing() const { return coefficientRing; }

  ring_elem fromCoeff(ring_elem coeff) const
  {
    auto result = new SymmetricRingPoly;
    if (!coefficientRing->is_zero(coeff)) result->terms.push_back({coeff, {}});
    return makePolyValue(result);
  }

  ring_elem basisElement(int basisId,
                         const std::string& display,
                         int order,
                         M2_arrayint index) const
  {
    rememberBasis(basisId, display);
    auto result = new SymmetricRingPoly;
    SymmetricMonomial monomial;
    monomial.atoms.push_back(makeAtom(order, basisId, index));
    result->terms.push_back({coefficientRing->one(), canonicalMonomial(monomial)});
    return makePolyValue(result);
  }

  bool getScalar(const SymmetricRingPoly *f, ring_elem &result) const
  {
    if (f->terms.empty())
      {
        result = coefficientRing->zero();
        return true;
      }
    if (f->terms.size() != 1) return false;
    if (!f->terms[0].monomial.atoms.empty()) return false;
    result = f->terms[0].coeff;
    return true;
  }

  std::string displayAtom(const SymmetricAtom& atom) const
  {
    return displayForBasis(atomBasisId(atom)) + "_" + displayIndex(atom);
  }

  std::string displayMonomial(const SymmetricMonomial& monomial) const
  {
    if (monomial.atoms.empty()) return "1";
    std::vector<std::string> factors;
    factors.reserve(monomial.atoms.size());
    for (const auto& atom : monomial.atoms) factors.push_back(displayAtom(atom));
    return join(factors, "*");
  }

  std::string elementString(ring_elem f) const
  {
    const auto *poly = polyValue(f);
    if (poly->terms.empty()) return "0";
    std::vector<std::string> pieces;
    pieces.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        if (term.monomial.atoms.empty())
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
    return join(pieces, " + ");
  }

  int elementWeight(ring_elem f) const
  {
    const auto *poly = polyValue(f);
    bool found = false;
    int result = 0;
    for (const auto& term : poly->terms)
      {
        int wt = 0;
        for (const auto& atom : term.monomial.atoms) wt += atomWeight(atom);
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

  virtual unsigned int computeHashValue(const ring_elem a) const
  {
    size_t seed = 7919237;
    const auto *poly = polyValue(a);
    for (const auto& term : poly->terms)
      {
        hashCombine(seed, coefficientRing->computeHashValue(term.coeff));
        for (int value : monomialKey(term.monomial)) hashCombine(seed, value);
      }
    return static_cast<unsigned int>(seed);
  }

  virtual void text_out(buffer &o) const
  {
    o << "SymmetricRing(";
    coefficientRing->text_out(o);
    o << ")";
  }

  virtual void elem_text_out(buffer &o,
                             const ring_elem f,
                             bool p_one = true,
                             bool p_plus = false,
                             bool p_parens = false) const
  {
    (void)p_one;
    if (p_plus && !is_zero(f)) o << '+';
    bool parens = p_parens && polyValue(f)->terms.size() > 1;
    if (parens) o << '(';
    o << elementString(f);
    if (parens) o << ')';
  }

  virtual ring_elem from_long(long n) const
  {
    return fromCoeff(coefficientRing->from_long(n));
  }

  virtual ring_elem from_int(mpz_srcptr n) const
  {
    return fromCoeff(coefficientRing->from_int(n));
  }

  virtual bool from_rational(mpq_srcptr q, ring_elem &result) const
  {
    ring_elem coeff;
    if (!coefficientRing->from_rational(q, coeff)) return false;
    result = fromCoeff(coeff);
    return true;
  }

  virtual bool promote(const Ring *Rf, const ring_elem f, ring_elem &result) const
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
        std::map<std::vector<int>, ring_elem> accumulator;
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

  virtual bool lift(const Ring *Rg, const ring_elem f, ring_elem &result) const
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
        std::map<std::vector<int>, ring_elem> accumulator;
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

  virtual bool is_unit(const ring_elem f) const
  {
    ring_elem coeff;
    return getScalar(polyValue(f), coeff) && coefficientRing->is_unit(coeff);
  }

  virtual bool is_zero(const ring_elem f) const
  {
    return polyValue(f)->terms.empty();
  }

  virtual bool is_equal(const ring_elem f, const ring_elem g) const
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

  virtual int compare_elems(const ring_elem f, const ring_elem g) const
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

  virtual ring_elem copy(const ring_elem f) const
  {
    return copyPolyValue(polyValue(f));
  }

  virtual void remove(ring_elem &f) const { (void)f; }

  virtual ring_elem negate(const ring_elem f) const
  {
    auto result = new SymmetricRingPoly;
    const auto *poly = polyValue(f);
    result->terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      result->terms.push_back({coefficientRing->negate(term.coeff), term.monomial});
    return makePolyValue(result);
  }

  virtual ring_elem add(const ring_elem f, const ring_elem g) const
  {
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
    if (left->terms.empty()) return copyPolyValue(right);
    if (right->terms.empty()) return copyPolyValue(left);

    if (compareMonomials(left->terms.back().monomial,
                         right->terms.front().monomial) == LT)
      return concatenateTerms(left, right);
    if (compareMonomials(right->terms.back().monomial,
                         left->terms.front().monomial) == LT)
      return concatenateTerms(right, left);

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
    return makePolyValue(result);
  }

  virtual ring_elem subtract(const ring_elem f, const ring_elem g) const
  {
    return add(f, negate(g));
  }

  SymmetricRingPoly *multByCoefficient(ring_elem coeff,
                                       const SymmetricRingPoly *poly) const
  {
    if (coefficientRing->is_equal(coeff, coefficientRing->one()))
      {
        auto result = new SymmetricRingPoly;
        result->terms.reserve(poly->terms.size());
        result->terms.insert(result->terms.end(), poly->terms.begin(), poly->terms.end());
        return result;
      }
    auto result = new SymmetricRingPoly;
    if (coefficientRing->is_zero(coeff)) return result;
    result->terms.reserve(poly->terms.size());
    for (const auto& term : poly->terms)
      {
        ring_elem c = coefficientRing->mult(coeff, term.coeff);
        if (!coefficientRing->is_zero(c)) result->terms.push_back({c, term.monomial});
      }
    return result;
  }

  virtual ring_elem mult(const ring_elem f, const ring_elem g) const
  {
    ring_elem scalar;
    const auto *left = polyValue(f);
    const auto *right = polyValue(g);
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
        return makePolyValue(result);
      }
    std::map<std::vector<int>, ring_elem> accumulator;
    for (const auto& lt : left->terms)
      for (const auto& rt : right->terms)
        {
          ring_elem coeff = coefficientRing->mult(lt.coeff, rt.coeff);
          auto monomial = multiplyMonomials(lt.monomial, rt.monomial);
          addToAccumulator(accumulator, monomial, coeff);
        }
    return fromAccumulator(accumulator);
  }

  ring_elem batchSum(engine_RawRingElementArray elements) const
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

  ring_elem batchProduct(engine_RawRingElementArray elements) const
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

  virtual ring_elem invert(const ring_elem f) const
  {
    (void)f;
    return zero();
  }

  virtual ring_elem divide(const ring_elem f, const ring_elem g) const
  {
    (void)f;
    (void)g;
    return zero();
  }

  virtual void syzygy(const ring_elem a,
                      const ring_elem b,
                      ring_elem &x,
                      ring_elem &y) const
  {
    (void)a;
    (void)b;
    x = zero();
    y = zero();
  }

  virtual ring_elem eval(const RingMap *map,
                         const ring_elem f,
                         int first_var) const
  {
    (void)f;
    (void)first_var;
    return map->get_ring()->zero();
  }
};

const SymmetricEngineRing *symmetricRingFromElement(const RingElement *f)
{
  auto result = dynamic_cast<const SymmetricEngineRing *>(f->get_ring());
  if (result == nullptr) ERROR("expected a SymmetricRings engine element");
  return result;
}

const SymmetricEngineRing *symmetricRingFromRing(const Ring *R)
{
  auto result = dynamic_cast<const SymmetricEngineRing *>(R);
  if (result == nullptr) ERROR("expected a SymmetricRings engine ring");
  return result;
}

} // namespace

const Ring *rawSymmetricRing(const Ring *A)
{
  try
    {
      return SymmetricEngineRing::create(A);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                M2_arrayint index)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->basisElement(basisId,
                                         fromM2String(displaySymbol),
                                         displayOrder,
                                         index);
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      return RingElement::make_raw(S, S->batchSum(elements));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      return RingElement::make_raw(S, S->batchProduct(elements));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_string rawSymmetricRingsElementToString(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return toM2String(S->elementString(f->get_value()));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

int rawSymmetricRingsElementWeight(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return -1;
      return S->elementWeight(f->get_value());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return -1;
    }
}
