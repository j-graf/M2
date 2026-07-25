// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_STORAGE_HPP_
#define M2_SYMMETRIC_RINGS_STORAGE_HPP_

#include "symmetric-rings/partitions.hpp"

#include "engine-includes.hpp"
#include "newdelete.hpp"
#include "rings/ring.hpp"
#include "rings/ringelem.hpp"

#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <utility>
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

// ============================================================================
// Semantic Tags And Conversion Metadata
// ============================================================================

enum class CombinatorialTag : uint32_t
{
  Plethysm = 1u << 0,
  LittlewoodRichardson = 1u << 1,
  HorizontalPieri = 1u << 2,
  VerticalPieri = 1u << 3,
  BorderStrips = 1u << 4
};

using CombinatorialTags = uint32_t;

constexpr CombinatorialTags combinatorialTagMask(CombinatorialTag tag)
{
  return static_cast<CombinatorialTags>(tag);
}

constexpr bool hasCombinatorialTag(CombinatorialTags tags,
                                   CombinatorialTag tag)
{
  return (tags & combinatorialTagMask(tag)) != 0;
}

class SymmetricConversionMetadata : public our_gc_cleanup
{
 public:
  SymmetricConversionMetadata() = default;

  // A metadata object owns ordinary STL containers, so it is finalized
  // independently of its containing polynomial.  Copy and move construction
  // must run the default our_gc_cleanup constructor to register that new
  // allocation with the collector before assigning the stored facts.
  SymmetricConversionMetadata(const SymmetricConversionMetadata& other)
      : our_gc_cleanup()
  {
    *this = other;
  }
  SymmetricConversionMetadata(SymmetricConversionMetadata&& other)
      : our_gc_cleanup()
  {
    *this = std::move(other);
  }
  SymmetricConversionMetadata& operator=(
      const SymmetricConversionMetadata&) = default;
  SymmetricConversionMetadata& operator=(
      SymmetricConversionMetadata&&) = default;

  // Exact canonical facts form one indivisible contract. Arithmetic may keep
  // individually valid hints after invalidating that contract, but no caller
  // may use the remaining fields to justify a canonical-workflow bypass.
  void invalidateExactExpressionFacts()
  {
    expressionFactsComplete = false;
    singleBasisElementCoefficientOne.reset();
  }

  // Numeric basis IDs describe one SymmetricEngineRing. Cross-ring fallback
  // reconstruction may preserve shapes and weights without preserving that
  // ID map, so discard the complete basis-identity cluster together.
  void discardRingLocalBasisFacts()
  {
    invalidateExactExpressionFacts();
    pureBasis.reset();
    expandedBasis.reset();
    factorBases.reset();
    singleBasisElementId.reset();
    singleBasisElementIndex.reset();
  }

  // True only when the exact canonical core needed for a workflow bypass was
  // attached together. Expensive selector-only profiles remain optional and
  // are enriched lazily from the expression when a policy consumes them.
  // Arithmetic that can change the core leaves this false. Predicates such as
  // single-term, single-basis-element, and product-free are deliberately not
  // stored: consumers derive them from the counts below.
  bool expressionFactsComplete = false;
  std::optional<int> pureBasis;
  std::optional<int> expandedBasis;
  std::optional<int> homogeneousWeight;
  std::optional<size_t> termCount;
  std::optional<size_t> scalarTermCount;
  std::optional<size_t> singleFactorTermCount;
  std::optional<size_t> maximumPartitionLength;
  std::optional<std::vector<int>> factorBases;
  std::optional<int> singleBasisElementId;
  std::optional<Partition> singleBasisElementIndex;
  std::optional<bool> singleBasisElementCoefficientOne;
  bool normalized = false;
  bool skewFree = false;
  bool collected = false;
};

// Polynomial objects are the engine's pervasive value representation. Keep
// their uncommon, comparatively large conversion profile out of line and
// share it across ordinary copies. A mutable access detaches first, preserving
// the previous value semantics without copying profiles that are only read.
class SymmetricConversionMetadataSlot
{
 public:
  SymmetricConversionMetadataSlot() = default;
  SymmetricConversionMetadataSlot(const SymmetricConversionMetadata& value)
      : mValue(new SymmetricConversionMetadata(value))
  {
  }
  SymmetricConversionMetadataSlot(SymmetricConversionMetadata&& value)
      : mValue(new SymmetricConversionMetadata(std::move(value)))
  {
  }
  SymmetricConversionMetadataSlot(
      const SymmetricConversionMetadataSlot& other)
      : mValue(other.mValue)
  {
  }
  SymmetricConversionMetadataSlot(
      SymmetricConversionMetadataSlot&& other) noexcept
      : mValue(other.mValue)
  {
    other.mValue = nullptr;
  }

  SymmetricConversionMetadataSlot& operator=(
      const SymmetricConversionMetadataSlot& other)
  {
    mValue = other.mValue;
    return *this;
  }
  SymmetricConversionMetadataSlot& operator=(
      SymmetricConversionMetadataSlot&& other) noexcept
  {
    mValue = other.mValue;
    other.mValue = nullptr;
    return *this;
  }
  SymmetricConversionMetadataSlot& operator=(
      const SymmetricConversionMetadata& value)
  {
    mValue = new SymmetricConversionMetadata(value);
    return *this;
  }
  SymmetricConversionMetadataSlot& operator=(
      SymmetricConversionMetadata&& value)
  {
    mValue = new SymmetricConversionMetadata(std::move(value));
    return *this;
  }

  explicit operator bool() const { return mValue != nullptr; }
  SymmetricConversionMetadata& operator*()
  {
    detach();
    return *mValue;
  }
  const SymmetricConversionMetadata& operator*() const { return *mValue; }
  SymmetricConversionMetadata *operator->()
  {
    detach();
    return mValue;
  }
  const SymmetricConversionMetadata *operator->() const { return mValue; }

 private:
  void detach()
  {
    if (mValue == nullptr) return;
    mValue = new SymmetricConversionMetadata(*mValue);
  }

  SymmetricConversionMetadata *mValue = nullptr;
};

class SymmetricRingPoly : public our_new_delete
{
 public:
  VECTOR(SymmetricTerm) terms;
  // Semantic operation tags are independent of conversion-profile metadata:
  // ordinary arithmetic can create them even when no conversion has run.
  CombinatorialTags combinatorialTags = 0;
  SymmetricConversionMetadataSlot conversionMetadata;
};

// ============================================================================
// Shared Keys And Coefficient Containers
// ============================================================================

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

// ============================================================================
// Storage Access And Atom Encoding
// ============================================================================

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
