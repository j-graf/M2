// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_SYMMETRIC_ENGINE_RING_HPP_
#define M2_SYMMETRIC_RINGS_SYMMETRIC_ENGINE_RING_HPP_

#include "symmetric-rings/partitions.hpp"
#include "symmetric-rings/storage.hpp"

#include "buffer.hpp"
#include "ring-elements/ring-element.hpp"
#include "ringmap.hpp"
#include "rings/ring.hpp"
#include "rings/ringelem.hpp"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace symmetric_rings {

class SymmetricEngineRing : public Ring
{

 private:
  struct SchurCompatibleFactor
  {
    enum Kind { General, Horizontal, Vertical, PowerSum, SchurExpansion };
    Kind kind;
    Partition index;
    CoeffMap expansion;
    int weight;
  };
  enum class SchurFactorMethod
  {
    AlreadySchur,
    ViaLittlewoodRichardson,
    ViaHorizontalPieri,
    ViaVerticalPieri,
    ViaBorderStrips,
    ViaLittlewoodRichardsonExpansion
  };
  enum class ProductToTargetRoute
  {
    ViaSchurCompatibleFactors,
    ViaMonomialLikeExpansion,
    ViaHallLittlewoodGenerators,
    ViaConvertRightFactor,
    ViaConvertLeftFactor,
    AlreadyInTarget,
    NoApplicableRoute
  };
  enum class ExpressionToTargetRoute
  {
    ViaSchurTriangularReduction,
    ViaHallLittlewoodTriangularReduction,
    ViaFactorwiseConversion,
    NoApplicableRoute
  };
  enum class BasisElementToPowerSumsRoute
  {
    AlreadyPowerSums,
    ViaCompleteClassicalFormula,
    ViaElementaryClassicalFormula,
    ViaHallLittlewoodGeneratorClassicalFormula,
    ViaSchurCharacters,
    ViaOmegaSchurCharacters,
    ViaMonomialTransition,
    ViaForgottenTransition,
    ViaHallLittlewoodRaisingOperators,
    ViaHallLittlewoodCapitalNormalization,
    ViaSkewSchurJacobiTrudiComplete,
    ViaSkewOmegaSchurJacobiTrudiElementary,
    ViaSkewHallLittlewood,
    NoApplicableRoute
  };

  const Ring *coefficientRing;
  mutable std::map<int, std::string> basisDisplays;
  mutable std::map<int, int> basisOrders;
  mutable std::map<int, bool> multiplicativeBases;
  mutable int powerSumBasisId = -1;
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
  mutable GCMap<int, CoeffMap> powerSumToCompleteMapCache;
  mutable GCMap<int, CoeffMap> powerSumToElementaryMapCache;
  mutable GCMap<int, CoeffMap> powerSumToQGeneratorMapCache;
  mutable GCMap<int, CoeffMap> powerSumToBGeneratorMapCache;
  mutable GCMap<int, GCMap<std::string, ring_elem>> monomialToPowerSumCache;
  mutable GCMap<int, GCMap<std::string, ring_elem>> forgottenToPowerSumCache;
  mutable std::map<std::string, std::vector<SchurConversionRecipeEntry>>
      powerSumsToSchurRecipeCache;
  mutable std::map<std::string, std::vector<LRProductTerm>> littlewoodRichardsonCoefficientProductCache;
  mutable std::map<std::pair<Partition, Partition>, std::vector<LRProductTerm>>
      littlewoodRichardsonTableauProductCache;
  mutable std::map<std::pair<Partition, Partition>, std::vector<LRProductTerm>>
      skewSchurToSchurViaLittlewoodRichardsonCache;
  mutable std::map<std::pair<Partition, Partition>, long> kostkaNumberCache;
  mutable std::map<std::pair<Partition, int>, std::vector<LRProductTerm>>
      schurTimesPowerSumViaBorderStripsCache;
  mutable std::map<std::pair<Partition, Partition>, std::vector<LRProductTerm>>
      monomialProductViaExponentSplittingsCache;
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
  bool isPowerSumBasis(int basisId) const;
  bool isForgottenDisplay(const std::string& display) const;
  void rememberBasis(int basisId,
                       const std::string& display,
                       int order,
                       bool isMultiplicative) const;
  void rememberBasesFrom(const SymmetricEngineRing *R) const;
  std::string displayForBasis(int basisId) const;
  int basisIdForDisplay(const std::string& display) const;
  int basisOrderForId(int basisId) const;
  const CharacterTable& characterTable(int degree) const;
  int characterTableValue(const CharacterTable& table, size_t row, size_t col) const;
  ring_elem rationalCoefficient(long numerator, long denominator) const;
  void clearHallLittlewoodCaches() const;
  ring_elem hallLittlewoodFactor(const Partition& mu) const;
  ring_elem basisElementFromIndex(int basisId,
                                    const std::string& display,
                                    int order,
                                    bool isMultiplicative,
                                    const Partition& index) const;
  ring_elem basisElementFromSkewIndex(int basisId,
                                        const std::string& display,
                                        int order,
                                        bool isMultiplicative,
                                        const Partition& outer,
                                        const Partition& inner) const;
  ring_elem basisPartElement(int basisId,
                               const std::string& display,
                               int order,
                               bool isMultiplicative,
                               int n) const;

  // ============================================================================
  // Topic-Specific Engine Declarations
  // ============================================================================

#include "symmetric-rings/basis-conversion-kernels.hpp"
#include "symmetric-rings/basis-conversion-products.hpp"
#include "symmetric-rings/basis-conversion-dispatch.hpp"
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

 public:
  // ============================================================================
  // Public Engine Ring Methods
  // ============================================================================

  explicit SymmetricEngineRing(const Ring *A);
  static SymmetricEngineRing *create(const Ring *A);
  const Ring *getCoefficientRing() const;
  void rememberBasisMetadata(int basisId,
                               const std::string& display,
                               int order,
                               bool isMultiplicative) const;
  bool setHallLittlewoodParameter(const RingElement *t) const;
  ring_elem fromCoeff(ring_elem coeff) const;
  ring_elem basisElement(int basisId,
                           const std::string& display,
                           int order,
                           bool isMultiplicative,
                           int innerLength,
                           M2_arrayint index) const;
  bool getScalar(const SymmetricRingPoly *f, ring_elem &result) const;
  bool hasPowerSumConversionHook(const SymmetricMonomial& monomial, size_t pos) const;
  std::string displayIndex(const SymmetricMonomial& monomial, size_t pos) const;
  std::string displayBasisElement(const SymmetricMonomial& monomial, size_t pos) const;
  std::string displayMonomial(const SymmetricMonomial& monomial) const;
  std::string elementString(ring_elem f, int maxTerms = -1) const;
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
