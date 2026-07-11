// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_STORAGE_HPP_
#define M2_SYMMETRIC_RINGS_STORAGE_HPP_

#include "symmetric-rings/partitions.hpp"

#include "newdelete.hpp"
#include "rings/ring.hpp"
#include "rings/ringelem.hpp"

#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace symmetric_rings {

// ============================================================================
// Storage Types And Ordering
// ============================================================================

constexpr size_t atomHeaderSize = 4;

struct SymmetricMonomial
{
  // Concatenated atom blocks:
  // [displayOrder, basisId, outerLength, innerLength, payload_1, ..., payload_n, ...]
  // For skew atoms, payload is outer shape followed by inner shape.
  VECTOR(int) data;
};

struct SymmetricTerm
{
  ring_elem coeff;
  SymmetricMonomial monomial;
};

enum class SymmetricConversionOrigin
{
  Unknown,
  Plethysm
};

struct SymmetricConversionMetadata
{
  std::optional<int> pureBasis;
  std::optional<int> expandedBasis;
  std::optional<int> homogeneousWeight;
  std::optional<size_t> termCount;
  std::optional<size_t> maximumPartitionLength;
  std::optional<double> density;
  std::optional<std::vector<int>> factorBases;
  std::optional<bool> singleBasisElement;
  std::optional<bool> singleTerm;
  std::optional<bool> noProducts;
  bool normalized = false;
  bool skewFree = false;
  bool collected = false;
  SymmetricConversionOrigin origin = SymmetricConversionOrigin::Unknown;
};

class SymmetricRingPoly : public our_new_delete
{
 public:
  VECTOR(SymmetricTerm) terms;
  std::optional<SymmetricConversionMetadata> conversionMetadata;
};

struct BasisIndexKey
{
  int basisId;
  std::vector<int> index;

  bool operator<(const BasisIndexKey& other) const
  {
    if (basisId != other.basisId) return basisId < other.basisId;
    return index < other.index;
  }
};

template <typename K, typename V, typename Compare = std::less<K>>
using GCMap = std::map<K, V, Compare, gc_allocator<std::pair<const K, V>>>;

using CoeffMap = GCMap<Partition, ring_elem>;
using RingElemVector = VECTOR(ring_elem);
using RingElemMatrix = VECTOR(RingElemVector);

struct SchurConversionRecipeEntry
{
  Partition lambda;
  std::vector<std::pair<size_t, int>> contributions;
};

struct LRProductTerm
{
  Partition nu;
  long coefficient;
};

struct OmegaTarget
{
  int basisId;
  int order;
  bool isMultiplicative;
};

struct InnerProductTarget
{
  int dualBasisId;
  int kind;
};

std::string fromM2String(M2_string s);
M2_string toM2String(const std::string& s);
std::string join(const std::vector<std::string>& parts,
                 const std::string& delimiter);
const SymmetricRingPoly *polyValue(ring_elem f);
SymmetricRingPoly *mutablePolyValue(ring_elem f);
ring_elem makePolyValue(SymmetricRingPoly *f);

template <typename T>
void hashCombine(size_t& seed, const T& val)
{
  seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

std::string coeffToString(const Ring *R, ring_elem c);
size_t atomLengthAt(const SymmetricMonomial& monomial, size_t pos);
int atomOrderAt(const SymmetricMonomial& monomial, size_t pos);
int atomBasisIdAt(const SymmetricMonomial& monomial, size_t pos);
bool atomIsSkewAt(const SymmetricMonomial& monomial, size_t pos);
int atomOuterLengthAt(const SymmetricMonomial& monomial, size_t pos);
int atomInnerLengthAt(const SymmetricMonomial& monomial, size_t pos);
int atomIndexLengthAt(const SymmetricMonomial& monomial, size_t pos);
bool blockLess(const std::vector<int>& a, const std::vector<int>& b);
bool monomialLess(const SymmetricMonomial& a, const SymmetricMonomial& b);
int compareMonomials(const SymmetricMonomial& a, const SymmetricMonomial& b);
std::vector<int> atomBlockAt(const SymmetricMonomial& monomial, size_t pos);
void appendAtomBlock(SymmetricMonomial& monomial, const std::vector<int>& block);
std::vector<int> makeAtomBlock(int displayOrder,
                               int basisId,
                               int innerLength,
                               M2_arrayint index);
std::vector<int> makeAtomBlock(int displayOrder,
                               int basisId,
                               int innerLength,
                               const std::vector<int>& index);
std::vector<int> monomialKey(const SymmetricMonomial& monomial);
SymmetricMonomial monomialFromKey(const std::vector<int>& key);
int monomialWeight(const SymmetricMonomial& monomial);
BasisIndexKey basisIndexKey(const SymmetricMonomial& monomial, size_t pos);

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
