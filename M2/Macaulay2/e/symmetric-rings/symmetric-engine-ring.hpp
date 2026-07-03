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
#include <string>
#include <vector>

namespace symmetric_rings {

class SymmetricEngineRing : public Ring
{

 private:
  const Ring *coefficientRing;
  mutable std::map<int, std::string> basisDisplays;
  mutable std::map<int, int> basisOrders;
  mutable std::map<int, bool> multiplicativeBases;
  mutable int powerSumBasisId = -1;
  mutable GCMap<int, ring_elem> hToPowerSumCache;
  mutable GCMap<int, ring_elem> eToPowerSumCache;
  mutable GCMap<int, ring_elem> qToPowerSumCache;
  mutable GCMap<int, ring_elem> bToPowerSumCache;
  mutable std::map<int, CharacterTable> characterTableCache;
  mutable GCMap<int, ring_elem> powerSumToCompleteCache;
  mutable GCMap<int, ring_elem> powerSumToElementaryCache;
  mutable GCMap<int, ring_elem> powerSumToQGeneratorCache;
  mutable GCMap<int, ring_elem> powerSumToBGeneratorCache;
  mutable GCMap<int, GCMap<std::string, ring_elem>> monomialToPowerSumCache;
  mutable std::map<std::string, std::vector<SchurConversionRecipeEntry>>
      powerSumsToSchurRecipeCache;
  mutable std::map<std::string, std::vector<LRProductTerm>> lrProductCache;
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
  // Basis Conversion
  // ============================================================================

  std::string jacobiTrudiCacheKey(int basisId,
                                    const Partition& outer,
                                    const Partition& inner) const;
  int popcountMask(size_t mask) const;
  int selectedGreaterThan(size_t mask, size_t col, size_t n) const;
  ring_elem jacobiTrudi(const Partition& outer,
                          const Partition& inner,
                          int basisId,
                          const std::string& display,
                          int order,
                          bool isMultiplicative) const;
  ring_elem scaled(ring_elem coeff, ring_elem f) const;
  ring_elem coefficientQuotient(ring_elem numerator, ring_elem denominator) const;
  ring_elem hallLittlewoodCFactor(const Partition& lambda) const;
  void addCoeff(CoeffMap& target, const Partition& index, ring_elem coeff) const;
  CoeffMap scaledCoeffMap(ring_elem coeff, const CoeffMap& source) const;
  CoeffMap addCoeffMaps(const CoeffMap& a, const CoeffMap& b) const;
  CoeffMap multiplyCoeffMaps(const CoeffMap& a, const CoeffMap& b) const;
  CoeffMap oneCoeffMap() const;
  bool isPartitionIndex(const Partition& p) const;
  int partitionPart(const Partition& p, size_t i) const;
  bool partitionContains(const Partition& outer, const Partition& inner) const;
  std::string lrProductKey(const Partition& lambda, const Partition& mu) const;
  long lrCoefficient(const Partition& lambda,
                       const Partition& content,
                       const Partition& nu) const;
  void partitionsContainingRec(const Partition& lambda,
                                 int addedWeight,
                                 size_t row,
                                 int previousPart,
                                 Partition& current,
                                 std::vector<Partition>& result) const;
  std::vector<Partition> partitionsContaining(const Partition& lambda,
                                                int addedWeight) const;
  const std::vector<LRProductTerm>& lrProduct(const Partition& a,
                                                const Partition& b) const;
  ring_elem multiplySchurElements(ring_elem f,
                                    ring_elem g,
                                    int schurId,
                                    const std::string& schurDisplay,
                                    int schurOrder) const;
  std::string powerSumsToSchurRecipeKey(int degree,
                                          const std::vector<Partition>& inputPartitions,
                                          bool omegaStyle) const;
  const std::vector<SchurConversionRecipeEntry>&
    powerSumsToSchurRecipe(int degree,
                           const std::vector<Partition>& inputPartitions,
                           bool omegaStyle) const;
  Partition leadingPartition(const CoeffMap& H) const;
  ring_elem hPartToPowerSums(int n) const;
  ring_elem ePartToPowerSums(int n) const;
  ring_elem hallLittlewoodPartToPowerSums(int n, bool omega) const;
  CoeffMap raisingExpansion(const Partition& lambda) const;
  CoeffMap raisingGeneratorMap(const Partition& lambda) const;
  ring_elem hallCapitalToPowerSums(const Partition& lambda, bool omega) const;
  ring_elem hallPToPowerSums(const Partition& lambda, bool omega) const;
  ring_elem schurToPowerSums(const Partition& lambda) const;
  ring_elem powerSumPartToComplete(int n, int hId, int hOrder) const;
  ring_elem powerSumPartToElementary(int n, int eId, int eOrder) const;
  ring_elem powerSumPartToHallGenerator(int n,
                                          int generatorId,
                                          const std::string& display,
                                          int generatorOrder,
                                          bool omega) const;
  CoeffMap powerSumPartToHallGeneratorMap(int n, bool omega) const;
  CoeffMap powerSumIndexToHallGeneratorMap(const Partition& index, bool omega) const;
  CoeffMap powerSumsToHallGeneratorMap(ring_elem f, bool omega) const;
  CoeffMap triangularReduceHallCapital(const CoeffMap& generatorMap,
                                         bool omega) const;
  ring_elem powerSumToSchurLike(const Partition& mu,
                                  int schurId,
                                  int schurOrder,
                                  const std::string& display,
                                  long sign) const;
  ring_elem powerSumToSchur(const Partition& mu, int schurId, int schurOrder) const;
  ring_elem powerSumsToSchurLike(ring_elem f,
                                   int schurId,
                                   int schurOrder,
                                   const std::string& display,
                                   bool omegaStyle) const;
  RingElemVector solveSquareSystem(RingElemMatrix M,
                                     RingElemVector v) const;
  ring_elem omegaPowerSums(ring_elem f) const;
  ring_elem monomialBasisToPowerSums(const Partition& lambda, bool forgotten) const;
  ring_elem powerSumIndexToMonomialTarget(const Partition& lambda,
                                            int targetBasisId,
                                            const std::string& targetDisplay,
                                            int targetDisplayOrder,
                                            bool forgotten) const;
  ring_elem coeffMapToElement(const CoeffMap& H,
                                int targetBasisId,
                                const std::string& targetDisplay,
                                int targetDisplayOrder,
                                bool targetIsMultiplicative) const;
  Partition replaceAdjacentPair(const Partition& alpha,
                                  size_t pos,
                                  int first,
                                  int second) const;
  ring_elem straightenSchurAtom(const Partition& alpha,
                                  const std::string& display) const;
  ring_elem omegaSchurAtomAsSchur(const Partition& alpha) const;
  ring_elem straightenHallCapitalAtom(const Partition& alpha,
                                        const std::string& display) const;
  ring_elem straightenAtom(const SymmetricMonomial& monomial, size_t pos) const;
  ring_elem straightenMonomial(const SymmetricMonomial& monomial) const;
  ring_elem straightenElement(ring_elem f) const;
  int requiredBasisIdForDisplay(const std::string& display) const;
  bool singleBasisIndexFromMonomial(const SymmetricMonomial& monomial,
                                      int basisId,
                                      Partition& index) const;
  CoeffMap coefficientsInBasis(ring_elem f, int basisId) const;
  bool coefficientsInBasisIfPossible(ring_elem f,
                                        int basisId,
                                        CoeffMap& result) const;
  ring_elem coefficientPairing(const CoeffMap& fCoeffs,
                                 const CoeffMap& gCoeffs) const;
  ring_elem powerSumInnerProductFactor(const Partition& lambda) const;
  ring_elem powerSumPairing(const CoeffMap& fCoeffs,
                              const CoeffMap& gCoeffs) const;
  std::map<int, InnerProductTarget> innerProductTargetMap(M2_arrayint innerProductMap) const;
  bool directHallInnerProductFromMetadata(
        ring_elem f,
        ring_elem g,
        const std::map<int, InnerProductTarget>& metadata,
        ring_elem& result) const;
  bool directPowerSumSchurInnerProduct(ring_elem f,
                                         ring_elem g,
                                         ring_elem& result) const;
  ring_elem basisElementForDisplay(const std::string& display,
                                     const Partition& index) const;
  ring_elem hallInnerProductElements(ring_elem f,
                                       ring_elem g,
                                       const std::map<int, InnerProductTarget>& metadata = {}) const;
  ring_elem skewQOrBFunction(const Partition& lambda,
                               const Partition& mu,
                               bool omega) const;
  ring_elem replaceSingleBasis(ring_elem f,
                                 int sourceBasisId,
                                 const std::string& targetDisplay) const;
  ring_elem skewPOrRToPowerSums(const Partition& lambda,
                                  const Partition& mu,
                                  bool omega) const;
  ring_elem skewHallLittlewoodToPowerSums(const Partition& lambda,
                                            const Partition& mu,
                                            const std::string& display) const;
  ring_elem powerSumsToHallCapitalTarget(ring_elem f,
                                           int targetBasisId,
                                           const std::string& targetDisplay,
                                           int targetDisplayOrder) const;
  CoeffMap schurGeneratorMap(const Partition& lambda,
                               bool omegaStyle,
                               int generatorId,
                               const std::string& generatorDisplay) const;
  CoeffMap triangularReduceSchur(const CoeffMap& generatorMap,
                                   bool omegaStyle,
                                   int generatorId,
                                   const std::string& generatorDisplay) const;
  bool directTriangularSchurConversion(ring_elem f,
                                         int targetBasisId,
                                         const std::string& targetDisplay,
                                         int targetDisplayOrder,
                                         ring_elem& result) const;
  bool directTriangularHallConversion(ring_elem f,
                                        int targetBasisId,
                                        const std::string& targetDisplay,
                                        int targetDisplayOrder,
                                        ring_elem& result) const;
  Partition atomIndex(const SymmetricMonomial& monomial, size_t pos) const;
  Partition atomOuterIndex(const SymmetricMonomial& monomial, size_t pos) const;
  Partition atomInnerIndex(const SymmetricMonomial& monomial, size_t pos) const;
  ring_elem atomToPowerSums(const SymmetricMonomial& monomial, size_t pos) const;
  ring_elem monomialToPowerSums(const SymmetricMonomial& monomial) const;
  ring_elem elementToPowerSums(ring_elem f) const;
  ring_elem powerSumElementFromIndex(const Partition& index) const;

  // ============================================================================
  // Plethysm
  // ============================================================================

  ring_elem adamsPowerSums(ring_elem f, int multiplier) const;
  ring_elem plethysmPowerSums(ring_elem fPowerSums, ring_elem gPowerSums) const;
  int singleBasisIdInMonomial(const SymmetricMonomial& monomial) const;
  int singleBasisId(ring_elem f) const;
  bool powerSumIndexFromMonomial(const SymmetricMonomial& monomial,
                                   Partition& index) const;
  ring_elem powerSumMonomialToTarget(const Partition& index,
                                       const std::string& targetDisplay,
                                       int targetBasisId,
                                       int targetDisplayOrder,
                                       bool targetIsMultiplicative) const;
  ring_elem powerSumsToTarget(ring_elem f,
                                int targetBasisId,
                                const std::string& targetDisplay,
                                int targetDisplayOrder,
                                bool targetIsMultiplicative) const;
  bool atomToDirectTarget(const SymmetricMonomial& monomial,
                            size_t pos,
                            int targetBasisId,
                            const std::string& targetDisplay,
                            int targetDisplayOrder,
                            bool targetIsMultiplicative,
                            ring_elem& result) const;
  bool monomialToDirectTarget(const SymmetricMonomial& monomial,
                                int targetBasisId,
                                const std::string& targetDisplay,
                                int targetDisplayOrder,
                                bool targetIsMultiplicative,
                                ring_elem& result) const;
  bool elementToDirectTarget(ring_elem f,
                               int targetBasisId,
                               const std::string& targetDisplay,
                               int targetDisplayOrder,
                               bool targetIsMultiplicative,
                               ring_elem& result) const;
  ring_elem omegaDirectAtom(const SymmetricMonomial& monomial,
                              size_t pos,
                              const std::map<int, OmegaTarget>& omegaTargets,
                              bool useSomega) const;
  ring_elem omegaMonomial(const SymmetricMonomial& monomial,
                            const std::map<int, OmegaTarget>& omegaTargets,
                            bool useSomega) const;
  std::map<int, OmegaTarget> omegaTargetMap(M2_arrayint omegaMap) const;

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
  std::string displayAtom(const SymmetricMonomial& monomial, size_t pos) const;
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
  ring_elem batchSum(engine_RawRingElementArray elements) const;
  ring_elem batchProduct(engine_RawRingElementArray elements) const;
  bool singleSchurPartition(ring_elem f, int schurId, Partition& lambda) const;
  std::string schurCompletePlethysmKey(int n,
                                         int schurId,
                                         const Partition& inner) const;
  ring_elem schurCompletePlethysm(int n,
                                    const Partition& inner,
                                    int powerSumBasisId0,
                                    const std::string& powerSumDisplay,
                                    int powerSumOrder,
                                    bool powerSumIsMultiplicative,
                                    int schurId,
                                    const std::string& schurDisplay,
                                    int schurOrder) const;
  ring_elem schurPlethysmJacobiTrudi(const Partition& outer,
                                       const Partition& inner,
                                       int powerSumBasisId0,
                                       const std::string& powerSumDisplay,
                                       int powerSumOrder,
                                       bool powerSumIsMultiplicative,
                                       int schurId,
                                       const std::string& schurDisplay,
                                       int schurOrder) const;
  bool schurPlethysmToSchur(ring_elem f,
                              ring_elem g,
                              int powerSumBasisId0,
                              const std::string& powerSumDisplay,
                              int powerSumOrder,
                              bool powerSumIsMultiplicative,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetOrder,
                              ring_elem& result) const;
  CoeffMap multiplySchurCoeffMapByRow(const CoeffMap& source, int row) const;
  CoeffMap hToSchurRecTransCoeffs(const CoeffMap& hCoeffs) const;
  ring_elem hToSchurViaRecTrans(ring_elem f,
                                  int hBasisId,
                                  const std::string& hDisplay,
                                  int hOrder,
                                  bool hIsMultiplicative,
                                  int schurId,
                                  const std::string& schurDisplay,
                                  int schurOrder) const;
  ring_elem jacobiTrudiBasis(int basisId,
                               const std::string& display,
                               int order,
                               bool isMultiplicative,
                               const Partition& outer,
                               const Partition& inner) const;
  ring_elem toBasis(ring_elem f,
                      int pBasisId,
                      const std::string& pDisplay,
                      int pOrder,
                      bool pIsMultiplicative,
                      int targetBasisId,
                      const std::string& targetDisplay,
                      int targetOrder,
                      bool targetIsMultiplicative) const;
  ring_elem plethysm(ring_elem f,
                       ring_elem g,
                       int pBasisId,
                       const std::string& pDisplay,
                       int pOrder,
                       bool pIsMultiplicative) const;
  ring_elem plethysmToBasis(ring_elem f,
                              ring_elem g,
                              int pBasisId,
                              const std::string& pDisplay,
                              int pOrder,
                              bool pIsMultiplicative,
                              int targetBasisId,
                              const std::string& targetDisplay,
                              int targetOrder,
                              bool targetIsMultiplicative) const;
  int uniformBasisId(ring_elem f) const;
  ring_elem omegaInvolution(ring_elem f, M2_arrayint omegaMap, bool useSomega) const;
  ring_elem straighten(ring_elem f) const;
  ring_elem hallInnerProduct(ring_elem f,
                               ring_elem g,
                               M2_arrayint innerProductMap) const;
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
