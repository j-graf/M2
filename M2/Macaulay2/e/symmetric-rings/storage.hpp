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
// Semantic Tags And Expression-Facts Cache
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

// Exact mathematical and structural facts about one realized expression.
// Operation-specific selectors consume this shared vocabulary rather than
// maintaining parallel profiles with subtly different derived predicates.
struct ExpressionFacts
{
  bool normalized = true;
  bool skewFree = true;
  bool collected = true;
  CombinatorialTags combinatorialTags = 0;
  size_t termCount = 0;
  size_t scalarTermCount = 0;
  size_t singleFactorTermCount = 0;
  size_t productTermCount = 0;
  size_t maximumFactorsPerTerm = 0;
  std::vector<int> factorBases;
  size_t skewFactorCount = 0;
  std::optional<int> pureBasis;
  std::optional<int> expandedBasis;
  std::optional<int> homogeneousWeight;
  size_t maximumPartitionLength = 0;
  std::optional<int> singleBasisElementId;
  std::optional<Partition> singleBasisElementIndex;
  std::optional<bool> singleBasisElementCoefficientOne;

  bool noProducts() const { return productTermCount == 0; }
  bool singleTerm() const { return termCount == 1; }
  bool mixedBasis() const { return factorBases.size() > 1; }
  bool provenanceUnknown() const { return combinatorialTags == 0; }
  bool provenanceMixed() const
  {
    return combinatorialTags != 0 &&
           (combinatorialTags & (combinatorialTags - 1)) != 0;
  }
  bool singleBasisElement() const
  {
    return termCount == 1 &&
           scalarTermCount == 0 &&
           singleFactorTermCount == 1 &&
           productTermCount == 0 &&
           normalized &&
           skewFree &&
           collected;
  }
  bool canonicalExpansionInBasis(int basisId) const
  {
    if (!normalized || !skewFree || !collected || !noProducts())
      return false;
    // Zero and scalar expressions are canonical in every basis.
    return singleFactorTermCount == 0 ||
           (expandedBasis && *expandedBasis == basisId);
  }
};

// A stored cache may contain either complete canonical ExpressionFacts or a
// small set of facts that arithmetic can preserve without rescanning its
// result. Only fields covered by this mask, or by CompleteCanonical as a
// whole-record contract, are authoritative.
enum class ExpressionFactKnowledge : uint32_t
{
  Normalized = 1u << 0,
  SkewFree = 1u << 1,
  Collected = 1u << 2,
  PureBasis = 1u << 3,
  ExpandedBasis = 1u << 4,
  HomogeneousWeight = 1u << 5,
  MaximumPartitionLength = 1u << 6,
  FactorBases = 1u << 7,
  CompleteCanonical = 1u << 8
};

constexpr uint32_t expressionFactKnowledgeMask(
    ExpressionFactKnowledge knowledge)
{
  return static_cast<uint32_t>(knowledge);
}

class ExpressionFactsCache : public our_gc_cleanup
{
 public:
  ExpressionFactsCache() = default;

  // A cache owns ordinary STL containers, so it is finalized
  // independently of its containing polynomial.  Copy and move construction
  // must run the default our_gc_cleanup constructor to register that new
  // allocation with the collector before assigning the stored facts.
  ExpressionFactsCache(const ExpressionFactsCache& other)
      : our_gc_cleanup()
  {
    *this = other;
  }
  ExpressionFactsCache(ExpressionFactsCache&& other)
      : our_gc_cleanup()
  {
    *this = std::move(other);
  }
  ExpressionFactsCache& operator=(
      const ExpressionFactsCache&) = default;
  ExpressionFactsCache& operator=(
      ExpressionFactsCache&&) = default;

  bool knows(ExpressionFactKnowledge fact) const
  {
    return (knowledge & expressionFactKnowledgeMask(fact)) != 0;
  }

  void remember(ExpressionFactKnowledge fact)
  {
    knowledge |= expressionFactKnowledgeMask(fact);
  }

  void forget(ExpressionFactKnowledge fact)
  {
    knowledge &= ~expressionFactKnowledgeMask(fact);
  }

  bool hasCompleteCanonicalFacts() const
  {
    return knows(
        ExpressionFactKnowledge::CompleteCanonical);
  }

  void setCompleteCanonicalFacts(ExpressionFacts value)
  {
    facts = std::move(value);
    knowledge =
        expressionFactKnowledgeMask(ExpressionFactKnowledge::Normalized) |
        expressionFactKnowledgeMask(ExpressionFactKnowledge::SkewFree) |
        expressionFactKnowledgeMask(ExpressionFactKnowledge::Collected) |
        expressionFactKnowledgeMask(ExpressionFactKnowledge::PureBasis) |
        expressionFactKnowledgeMask(ExpressionFactKnowledge::ExpandedBasis) |
        expressionFactKnowledgeMask(
            ExpressionFactKnowledge::HomogeneousWeight) |
        expressionFactKnowledgeMask(
            ExpressionFactKnowledge::MaximumPartitionLength) |
        expressionFactKnowledgeMask(ExpressionFactKnowledge::FactorBases) |
        expressionFactKnowledgeMask(
            ExpressionFactKnowledge::CompleteCanonical);
  }

  // Coefficient operations can remove terms without changing the surviving
  // monomials. Preserve structural postconditions, but discard every fact
  // whose exact value depends on which terms remain.
  void discardSupportDependentFacts()
  {
    forget(ExpressionFactKnowledge::CompleteCanonical);
    if (knows(ExpressionFactKnowledge::Normalized) &&
        !facts.normalized)
      {
        forget(ExpressionFactKnowledge::Normalized);
        facts.normalized = true;
      }
    if (knows(ExpressionFactKnowledge::SkewFree) &&
        !facts.skewFree)
      {
        forget(ExpressionFactKnowledge::SkewFree);
        facts.skewFree = true;
      }
    if (knows(ExpressionFactKnowledge::HomogeneousWeight) &&
        !facts.homogeneousWeight)
      {
        forget(ExpressionFactKnowledge::HomogeneousWeight);
        facts.homogeneousWeight.reset();
      }
    forget(ExpressionFactKnowledge::PureBasis);
    forget(ExpressionFactKnowledge::ExpandedBasis);
    forget(ExpressionFactKnowledge::MaximumPartitionLength);
    forget(ExpressionFactKnowledge::FactorBases);
    facts.pureBasis.reset();
    facts.expandedBasis.reset();
    facts.maximumPartitionLength = 0;
    facts.factorBases.clear();
    facts.singleBasisElementId.reset();
    facts.singleBasisElementIndex.reset();
    facts.singleBasisElementCoefficientOne.reset();
  }

  // Numeric basis IDs describe one SymmetricEngineRing. Cross-ring fallback
  // reconstruction may preserve shapes and weights without preserving that
  // ID map, so discard the complete basis-identity cluster together.
  void discardRingLocalBasisFacts()
  {
    forget(ExpressionFactKnowledge::CompleteCanonical);
    forget(ExpressionFactKnowledge::PureBasis);
    forget(ExpressionFactKnowledge::ExpandedBasis);
    forget(ExpressionFactKnowledge::FactorBases);
    facts.pureBasis.reset();
    facts.expandedBasis.reset();
    facts.factorBases.clear();
    facts.singleBasisElementId.reset();
    facts.singleBasisElementIndex.reset();
    facts.singleBasisElementCoefficientOne.reset();
  }

  ExpressionFacts facts;

 private:
  uint32_t knowledge = 0;
};

// Polynomial objects are the engine's pervasive value representation. Keep
// their uncommon, comparatively large expression-facts cache out of line and
// share it across ordinary copies. A mutable access detaches first, preserving
// the previous value semantics without copying caches that are only read.
class ExpressionFactsCacheSlot
{
 public:
  ExpressionFactsCacheSlot() = default;
  ExpressionFactsCacheSlot(const ExpressionFactsCache& value)
      : mValue(new ExpressionFactsCache(value))
  {
  }
  ExpressionFactsCacheSlot(ExpressionFactsCache&& value)
      : mValue(new ExpressionFactsCache(std::move(value)))
  {
  }
  ExpressionFactsCacheSlot(
      const ExpressionFactsCacheSlot& other)
      : mValue(other.mValue)
  {
  }
  ExpressionFactsCacheSlot(
      ExpressionFactsCacheSlot&& other) noexcept
      : mValue(other.mValue)
  {
    other.mValue = nullptr;
  }

  ExpressionFactsCacheSlot& operator=(
      const ExpressionFactsCacheSlot& other)
  {
    mValue = other.mValue;
    return *this;
  }
  ExpressionFactsCacheSlot& operator=(
      ExpressionFactsCacheSlot&& other) noexcept
  {
    mValue = other.mValue;
    other.mValue = nullptr;
    return *this;
  }
  ExpressionFactsCacheSlot& operator=(
      const ExpressionFactsCache& value)
  {
    mValue = new ExpressionFactsCache(value);
    return *this;
  }
  ExpressionFactsCacheSlot& operator=(
      ExpressionFactsCache&& value)
  {
    mValue = new ExpressionFactsCache(std::move(value));
    return *this;
  }

  explicit operator bool() const { return mValue != nullptr; }
  ExpressionFactsCache& operator*()
  {
    detach();
    return *mValue;
  }
  const ExpressionFactsCache& operator*() const { return *mValue; }
  ExpressionFactsCache *operator->()
  {
    detach();
    return mValue;
  }
  const ExpressionFactsCache *operator->() const { return mValue; }

 private:
  void detach()
  {
    if (mValue == nullptr) return;
    mValue = new ExpressionFactsCache(*mValue);
  }

  ExpressionFactsCache *mValue = nullptr;
};

class SymmetricRingPoly : public our_new_delete
{
 public:
  VECTOR(SymmetricTerm) terms;
  // Semantic operation tags are independent of structural fact caching:
  // ordinary arithmetic can create them even when no conversion has run.
  CombinatorialTags combinatorialTags = 0;
  ExpressionFactsCacheSlot expressionFactsCache;
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
