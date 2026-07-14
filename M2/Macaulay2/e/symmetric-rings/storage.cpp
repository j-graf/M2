// Copyright 2026

#include "symmetric-rings/storage.hpp"

#include "buffer.hpp"

#include <algorithm>
#include <sstream>

namespace symmetric_rings {

// ============================================================================
// M2 And Engine Value Wrappers
// ============================================================================

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

SymmetricRingPoly *mutablePolyValue(ring_elem f)
{
  return const_cast<SymmetricRingPoly *>(polyValue(f));
}

ring_elem makePolyValue(SymmetricRingPoly *f)
{
  return ring_elem(reinterpret_cast<const void *>(f));
}

// ============================================================================
// Atom Access And Canonical Ordering
// ============================================================================

std::string coeffToString(const Ring *R, ring_elem c)
{
  buffer o;
  R->elem_text_out(o, c);
  return std::string(o.str());
}

size_t atomLengthAt(const SymmetricMonomial& monomial, size_t pos)
{
  if (pos + 3 >= monomial.data.size()) return monomial.data.size() - pos;
  return atomHeaderSize + monomial.data[pos + 2] + monomial.data[pos + 3];
}

int atomOrderAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos];
}

int atomBasisIdAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos + 1];
}

bool atomIsSkewAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos + 3] > 0;
}

int atomOuterLengthAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos + 2];
}

int atomInnerLengthAt(const SymmetricMonomial& monomial, size_t pos)
{
  return monomial.data[pos + 3];
}

int atomIndexLengthAt(const SymmetricMonomial& monomial, size_t pos)
{
  return atomOuterLengthAt(monomial, pos) + atomInnerLengthAt(monomial, pos);
}

bool blockLess(const std::vector<int>& a, const std::vector<int>& b)
{
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

bool monomialLess(const SymmetricMonomial& a, const SymmetricMonomial& b)
{
  return std::lexicographical_compare(
      a.data.begin(), a.data.end(), b.data.begin(), b.data.end());
}

int compareMonomials(const SymmetricMonomial& a, const SymmetricMonomial& b)
{
  if (monomialLess(a, b)) return LT;
  if (monomialLess(b, a)) return GT;
  return EQ;
}

// ============================================================================
// Atom Construction And Monomial Keys
// ============================================================================

std::vector<int> atomBlockAt(const SymmetricMonomial& monomial, size_t pos)
{
  size_t len = atomLengthAt(monomial, pos);
  return std::vector<int>(monomial.data.begin() + pos,
                          monomial.data.begin() + pos + len);
}

void appendAtomBlock(SymmetricMonomial& monomial, const std::vector<int>& block)
{
  monomial.data.insert(monomial.data.end(), block.begin(), block.end());
}

std::vector<int> makeAtomBlock(int displayOrder,
                               int basisId,
                               int innerLength,
                               M2_arrayint index)
{
  std::vector<int> block;
  int payloadLength = index == nullptr ? 0 : index->len;
  block.push_back(displayOrder);
  block.push_back(basisId);
  block.push_back(payloadLength - innerLength);
  block.push_back(innerLength);
  if (index != nullptr)
    for (int i = 0; i < index->len; ++i) block.push_back(index->array[i]);
  return block;
}

std::vector<int> makeAtomBlock(int displayOrder,
                               int basisId,
                               int innerLength,
                               const std::vector<int>& index)
{
  std::vector<int> block;
  block.reserve(atomHeaderSize + index.size());
  block.push_back(displayOrder);
  block.push_back(basisId);
  block.push_back(static_cast<int>(index.size()) - innerLength);
  block.push_back(innerLength);
  block.insert(block.end(), index.begin(), index.end());
  return block;
}

std::vector<int> monomialKey(const SymmetricMonomial& monomial)
{
  return std::vector<int>(monomial.data.begin(), monomial.data.end());
}

SymmetricMonomial monomialFromKey(const std::vector<int>& key)
{
  SymmetricMonomial result;
  result.data.insert(result.data.end(), key.begin(), key.end());
  return result;
}

int monomialWeight(const SymmetricMonomial& monomial)
{
  int result = 0;
  size_t pos = 0;
  while (pos < monomial.data.size())
    {
      int n = atomIndexLengthAt(monomial, pos);
      if (atomIsSkewAt(monomial, pos))
        {
          int outerLength = atomOuterLengthAt(monomial, pos);
          int innerLength = atomInnerLengthAt(monomial, pos);
          if (outerLength + innerLength == n)
            {
              for (int i = 0; i < outerLength; ++i)
                result += monomial.data[pos + atomHeaderSize + i];
              for (int i = outerLength; i < n; ++i)
                result -= monomial.data[pos + atomHeaderSize + i];
            }
        }
      else
        {
          for (int i = 0; i < n; ++i) result += monomial.data[pos + atomHeaderSize + i];
        }
      pos += atomHeaderSize + n;
    }
  return result;
}

BasisIndexKey basisIndexKey(const SymmetricMonomial& monomial, size_t pos)
{
  BasisIndexKey key;
  key.basisId = atomBasisIdAt(monomial, pos);
  int n = atomIndexLengthAt(monomial, pos);
  key.index.reserve(n);
  for (int i = 0; i < n; ++i) key.index.push_back(monomial.data[pos + atomHeaderSize + i]);
  return key;
}


} // namespace symmetric_rings
