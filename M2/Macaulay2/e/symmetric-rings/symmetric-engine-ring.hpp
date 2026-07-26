// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_SYMMETRIC_ENGINE_RING_HPP_
#define M2_SYMMETRIC_RINGS_SYMMETRIC_ENGINE_RING_HPP_

#include "symmetric-rings/partitions.hpp"
#include "symmetric-rings/expression-conditions.hpp"
#include "symmetric-rings/operation-records.hpp"
#include "symmetric-rings/storage.hpp"

#include "buffer.hpp"
#include "ring-elements/ring-element.hpp"
#include "ringmap.hpp"
#include "rings/ring.hpp"
#include "rings/ringelem.hpp"

#include <array>
#include <functional>
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace symmetric_rings {

class SymmetricEngineRing : public Ring
{

 private:
  struct ComputationLimits
  {
    size_t maxWeight = 200;
    size_t maxEnumeratedPartitions = 250000;
    size_t maxGeneratedTerms = 250000;
    size_t maxRecursiveStates = 5000000;
    size_t maxCacheEntries = 250000;
    size_t maxCharacterCacheEntries = 2000000;
    size_t maxDeterminantStates = 262144;
    size_t maxEstimatedMemoryBytes = 512ULL * 1024ULL * 1024ULL;
  };

  enum class BasisKind
  {
    Custom = 0,
    PowerSum,
    Complete,
    Elementary,
    Monomial,
    Forgotten,
    Schur,
    SchurOmega,
    HallLittlewoodQGenerator,
    HallLittlewoodBGenerator,
    HallLittlewoodQ,
    HallLittlewoodB,
    HallLittlewoodP,
    HallLittlewoodPOmega
  };

  struct BasisDescriptor
  {
    std::string canonicalKey;
    std::string displaySymbol;
    int displayOrder;
    bool multiplicative;
    BasisKind kind;
  };

  // Ring-local realization of a mathematical basis.  Conversion and
  // multiplication share this neutral descriptor; neither pipeline owns a
  // second endpoint representation.
  struct RingBasis
  {
    BasisKind kind = BasisKind::Custom;
    int id = -1;
    const std::string *display = nullptr;
    int order = 0;

    const std::string& displayName() const
    {
      return *display;
    }
  };

  const Ring *coefficientRing;
  mutable ComputationLimits computationLimits;
  mutable std::map<int, BasisDescriptor> basisDescriptors;
  mutable std::map<BasisKind, int> basisIdsByKind;
  // Stable declarative plan IDs are process-wide, while their numeric basis
  // endpoints are ring-local. Cache that final resolution per ring.
  mutable std::map<std::string, std::pair<int, int>>
      basisConversionEndpointsForRingCache;

  // Classical and Hall-Littlewood conversion state.
  mutable GCMap<int, ring_elem> completeToPowerSumsCache;
  mutable GCMap<int, ring_elem> elementaryToPowerSumsCache;
  mutable GCMap<int, ring_elem> hallLittlewoodQGeneratorToPowerSumsCache;
  mutable GCMap<int, ring_elem> hallLittlewoodBGeneratorToPowerSumsCache;
  mutable GCMap<int, CoeffMap> hallLittlewoodQGeneratorToPowerSumsQuotientMapCache;
  mutable GCMap<int, CoeffMap> hallLittlewoodBGeneratorToPowerSumsQuotientMapCache;
  mutable GCMap<int, ring_elem> hallLittlewoodPartFactorCache;
  mutable GCMap<Partition, ring_elem> hallLittlewoodCFactorCache;
  mutable GCMap<int, CoeffMap> hallLittlewoodSingleCycleCapitalGreenMapCache;
  mutable GCMap<int, CoeffMap> hallLittlewoodSingleCycleNormalizedGreenMapCache;
  mutable GCMap<Partition, CoeffMap> hallLittlewoodPowerSumToCapitalColumnCache;
  mutable GCMap<std::string, ring_elem>
      hallLittlewoodPowerSumToCapitalCoefficientCache;
  mutable GCMap<Partition, CoeffMap> hallLittlewoodRaisingGeneratorMapCache;
  mutable std::map<int, CharacterTable> characterTableCache;
  mutable size_t characterCacheEntryCount = 0;
  mutable size_t computationCacheEntryCount = 0;
  mutable GCMap<int, CoeffMap> powerSumToCompleteMapCache;
  mutable GCMap<int, CoeffMap> powerSumToElementaryMapCache;
  mutable GCMap<int, CoeffMap> powerSumToQGeneratorMapCache;
  mutable GCMap<int, CoeffMap> powerSumToBGeneratorMapCache;
  mutable GCMap<int, GCMap<std::string, ring_elem>> monomialToPowerSumCache;
  mutable GCMap<int, GCMap<std::string, ring_elem>> forgottenToPowerSumCache;
  // A runtime invariant guard: executing a fixed conversion plan may recurse
  // into its component plans, but it must never invoke performance selection.
  mutable size_t basisConversionPlanExecutionDepth = 0;
  // A strict binary kernel is a terminal combinatorial formula. It must never
  // reenter binary selection or an outer multiplication workflow.
  mutable size_t binaryMultiplicationKernelExecutionDepth = 0;

  // Schur conversion and product state.
  mutable std::map<std::string, std::vector<SchurConversionRecipeEntry>>
      powerSumsToSchurRecipeCache;
  mutable std::map<std::string, std::vector<PartitionCoefficientTerm>> littlewoodRichardsonCoefficientProductCache;
  mutable std::map<std::pair<Partition, Partition>, std::vector<PartitionCoefficientTerm>>
      littlewoodRichardsonTableauProductCache;
  mutable std::map<std::pair<Partition, Partition>, std::vector<PartitionCoefficientTerm>>
      skewSchurToSchurViaLittlewoodRichardsonCache;
  mutable std::map<std::pair<Partition, Partition>, mpz_class> kostkaNumberCache;
  mutable std::map<std::pair<Partition, int>, std::vector<PartitionCoefficientTerm>>
      schurTimesPowerSumViaBorderStripsCache;
  mutable std::map<std::pair<Partition, int>, std::vector<PartitionCoefficientTerm>>
      schurTimesPowerSumViaAbacusRimHooksCache;
  mutable std::map<std::pair<Partition, Partition>, std::vector<PartitionCoefficientTerm>>
      monomialProductViaExponentSplittingsCache;

  // Shared scalar, plethysm, and determinant state.
  mutable GCMap<long, ring_elem> smallIntegerCoeffCache;
  // Cache for Schur plethysm via Adams operations and Jacobi-Trudi.
  mutable GCMap<std::string, ring_elem> schurCompletePlethysmCache;
  mutable GCMap<std::string, ring_elem> hJacobiTrudiCache;
  mutable GCMap<std::string, ring_elem> eJacobiTrudiCache;
  mutable ring_elem hallLittlewoodParameter;

  // ============================================================================
  // Basis Metadata
  // ============================================================================

  bool isMultiplicativeBasis(int basisId) const;
  BasisKind basisKindForId(int basisId) const;
  BasisKind basisKindFromCanonicalKey(const std::string& key) const;
  static bool isHallLittlewoodCapitalBasisKind(BasisKind kind);
  bool hasBasisKind(int basisId, BasisKind kind) const;
  int basisIdForKind(BasisKind kind) const;
  int registeredPowerSumBasisId() const;
  int requiredBasisIdForKind(BasisKind kind) const;
  const char *basisKindName(BasisKind kind) const;
  const BasisDescriptor& requireBasis(int basisId) const;
  void rememberBasesFrom(const SymmetricEngineRing *R) const;
  std::string displayForBasis(int basisId) const;
  std::string basisKeyForId(int basisId) const;
  int basisOrderForId(int basisId) const;
  const CharacterTable& characterTable(int degree) const;
  const std::vector<mpz_class>& characterTableRow(
      const CharacterTable& table, size_t row) const;
  mpz_class characterTableValue(
      const CharacterTable& table, size_t row, size_t col) const;
  mpz_class characterValueWithinLimits(
      const Partition& lambda, const Partition& mu) const;
  mpz_class pToMonomialCoefficientWithinLimits(
      const Partition& lambda, const Partition& mu) const;
  std::vector<Partition> partitionsOfWithinLimits(
      int degree, const char *operation) const;
  bool estimatedMemoryWithinLimit(
      size_t count, size_t bytesPerItem, const char *operation) const;
  bool determinantStatesWithinLimit(size_t states, const char *operation) const;
  [[noreturn]] void recursiveStateLimitExceeded(
      const char *operation) const;
  bool consumeRecursiveState(size_t& states, const char *operation) const
  {
    // This check occurs at every node in several combinatorial recursions.
    // Keep the successful path inline and leave diagnostic construction on
    // the cold, out-of-line failure path.
    if (states < computationLimits.maxRecursiveStates)
      {
        ++states;
        return true;
      }
    recursiveStateLimitExceeded(operation);
  }
  void requireCacheEntryCapacity(const char *operation) const;
  void requireCacheEntryCapacity(size_t count, const char *operation) const;
  void requirePartitionWithinWeightLimit(
      const Partition& index, const char *operation) const;
  void requireWeightWithinLimit(long long weight, const char *operation) const;
  void requireMonomialWithinWeightLimit(
      const SymmetricMonomial& monomial, const char *operation) const;
  ring_elem rationalCoefficient(long numerator, long denominator) const;
  ring_elem rationalCoefficient(
      const mpz_class& numerator,
      const mpz_class& denominator) const;
  void clearHallLittlewoodCaches() const;
  ring_elem hallLittlewoodFactor(const Partition& mu) const;
  ring_elem basisElementFromIndex(int basisId,
                                    const Partition& index) const;
  ring_elem basisElementFromSkewIndex(int basisId,
                                        const Partition& outer,
                                        const Partition& inner) const;
  ring_elem basisPartElement(int basisId, int n) const;

  // ============================================================================
  // Topic-Specific Engine Declarations
  // ============================================================================

#include "symmetric-rings/basis-conversion-kernels.hpp"
#include "symmetric-rings/multiplication-kernels.hpp"
#include "symmetric-rings/expression-helpers.hpp"
#include "symmetric-rings/basis-conversion.hpp"
#include "symmetric-rings/basis-normalization.hpp"
#include "symmetric-rings/multiplication-picker.hpp"
#include "symmetric-rings/multiplication-folds.hpp"
#include "symmetric-rings/binary-multiplication.hpp"
#include "symmetric-rings/multiplication.hpp"
#include "symmetric-rings/basis-coefficient.hpp"
#include "symmetric-rings/omega.hpp"
#include "symmetric-rings/plethysm.hpp"
#include "symmetric-rings/inner-product-dispatch.hpp"
#include "symmetric-rings/inner-product-kernels.hpp"

  // ============================================================================
  // Polynomial Arithmetic
  // ============================================================================

  SymmetricMonomial canonicalMonomial(const SymmetricMonomial& monomial) const;
  bool isSinglePowerSumBlock(const SymmetricMonomial& monomial) const;
  SymmetricMonomial multiplyPowerSumMonomials(const SymmetricMonomial& a,
                                               const SymmetricMonomial& b) const;
  SymmetricMonomial multiplyMonomials(const SymmetricMonomial& a,
                                        const SymmetricMonomial& b) const;
  ring_elem copyPolyValue(const SymmetricRingPoly *poly) const;
  void appendTermIfNonZero(VECTOR(SymmetricTerm)& terms,
                             ring_elem coeff,
                             const SymmetricMonomial& monomial) const;
  ring_elem fromTermVector(VECTOR(SymmetricTerm)& terms, bool isSorted) const;
  ring_elem concatenateTerms(const SymmetricRingPoly *first,
                               const SymmetricRingPoly *second) const;
  void addToAccumulator(GCMap<std::vector<int>, ring_elem>& accumulator,
                          const SymmetricMonomial& monomial,
                          ring_elem coeff) const;
  ring_elem fromAccumulator(
        const GCMap<std::vector<int>, ring_elem>& accumulator) const;
  bool promoteInputElement(const RingElement *input, ring_elem &result) const;
  CombinatorialTags selectMultiplicationTags(ring_elem f, ring_elem g) const;

 public:
  // ============================================================================
  // Public Engine Ring Methods
  // ============================================================================

  explicit SymmetricEngineRing(const Ring *A);
  static SymmetricEngineRing *create(const Ring *A);
  const Ring *getCoefficientRing() const;
  void rememberBasisMetadata(int basisId,
                               const std::string& canonicalBasisKey,
                               const std::string& display,
                               int order,
                               bool isMultiplicative) const;
  bool setHallLittlewoodParameter(const RingElement *t) const;
  void setComputationLimits(size_t maxWeight,
                            size_t maxEnumeratedPartitions,
                            size_t maxGeneratedTerms,
                            size_t maxRecursiveStates,
                            size_t maxCacheEntries,
                            size_t maxCharacterCacheEntries,
                            size_t maxDeterminantStates,
                            size_t maxEstimatedMemoryMB) const;
  ring_elem fromCoeff(ring_elem coeff) const;
  ring_elem basisElement(int basisId,
                           int innerLength,
                           M2_arrayint index) const;
  bool getScalar(const SymmetricRingPoly *f, ring_elem &result) const;
#include "symmetric-rings/presentation.hpp"
  int elementWeight(ring_elem f) const;
  virtual unsigned int computeHashValue(const ring_elem a) const;
  virtual void text_out(buffer &o) const;
  virtual void elem_text_out(buffer &o,
                               const ring_elem f,
                               bool p_one = true,
                               bool p_plus = false,
                               bool p_parens = false) const;
  virtual ring_elem from_long(long n) const;
  virtual ring_elem from_int(mpz_srcptr n) const;
  virtual bool from_rational(mpq_srcptr q, ring_elem &result) const;
  virtual bool promote(const Ring *Rf, const ring_elem f, ring_elem &result) const;
  virtual bool lift(const Ring *Rg, const ring_elem f, ring_elem &result) const;
  virtual bool is_unit(const ring_elem f) const;
  virtual bool is_zero(const ring_elem f) const;
  virtual bool is_equal(const ring_elem f, const ring_elem g) const;
  virtual int compare_elems(const ring_elem f, const ring_elem g) const;
  virtual ring_elem copy(const ring_elem f) const;
  virtual void remove(ring_elem &f) const;
  virtual ring_elem negate(const ring_elem f) const;
  virtual ring_elem add(const ring_elem f, const ring_elem g) const;
  virtual ring_elem subtract(const ring_elem f, const ring_elem g) const;
  SymmetricRingPoly *multByCoefficient(ring_elem coeff,
                                         const SymmetricRingPoly *poly) const;
  virtual ring_elem mult(const ring_elem f, const ring_elem g) const;
  bool promoteCollectedExpansion(const RingElement *input,
                                   ring_elem& result) const;
  bool liftCollectedExpansion(const RingElement *input,
                              ring_elem& result) const;
  ring_elem batchSum(engine_RawRingElementArray elements) const;
  ring_elem batchProduct(engine_RawRingElementArray elements) const;
  int uniformBasisId(ring_elem f) const;
  virtual ring_elem invert(const ring_elem f) const;
  virtual ring_elem divide(const ring_elem f, const ring_elem g) const;
  virtual void syzygy(const ring_elem a,
                        const ring_elem b,
                        ring_elem &x,
                        ring_elem &y) const;
  virtual ring_elem eval(const RingMap *map,
                           const ring_elem f,
                           int first_var) const;
};

const SymmetricEngineRing *symmetricRingFromElement(const RingElement *f);
const SymmetricEngineRing *symmetricRingFromRing(const Ring *R);

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
