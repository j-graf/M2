// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_INNER_PRODUCT_DISPATCH_HPP_
#define M2_SYMMETRIC_RINGS_INNER_PRODUCT_DISPATCH_HPP_

// Declaration fragment included inside SymmetricEngineRing.

  // ============================================================================
  // Inner-Product Requests And Profiles
  // ============================================================================

  struct InnerProductProfile
  {
    std::optional<int> pureBasis;
    std::optional<int> expandedBasis;
    std::optional<int> homogeneousWeight;
    std::optional<size_t> termCount;
    mutable std::optional<size_t> maximumPartitionLength;
    std::optional<int> singleBasisId;
    Partition singleIndex;
    KnownState singleBasisElement = KnownState::Unknown;
    KnownState noProducts = KnownState::Unknown;
    KnownState normalized = KnownState::Unknown;
    KnownState skewFree = KnownState::Unknown;
    KnownState collected = KnownState::Unknown;
  };

  enum class InnerProductKind
  {
    OrdinaryHall = 0,
    HallLittlewood = 1,
    SchurQ = 2,
    Macdonald = 3
  };

  enum class PowerSumPairingKind
  {
    OrdinaryHall,
    HallLittlewood,
    SchurQ,
    Macdonald
  };

  struct InnerProductContext
  {
    InnerProductKind kind;
    PowerSumPairingKind powerSumPairing;
    std::map<int, InnerProductTarget> pairings;
  };

  struct InnerProductRequest
  {
    ring_elem left;
    ring_elem right;
    InnerProductContext context;
    InnerProductProfile leftProfile;
    InnerProductProfile rightProfile;
  };

  InnerProductProfile inferInnerProductProfile(ring_elem f) const;
  size_t innerProductMaximumPartitionLength(
      ring_elem f,
      const InnerProductProfile& profile) const;
  InnerProductRequest buildInnerProductRequest(
      ring_elem f,
      ring_elem g,
      InnerProductKind kind,
      std::map<int, InnerProductTarget> metadata) const;
  InnerProductKind innerProductKindFromCode(int kindCode) const;
  PowerSumPairingKind powerSumPairingKind(
      InnerProductKind kind) const;
  const char *innerProductKindName(InnerProductKind kind) const;

  // ============================================================================
  // Candidate Cost And Selection
  // ============================================================================

  enum class InnerProductPipeline
  {
    DiagonalBasis,
    SingleBasisElement,
    PowerSumsStructured,
    FallbackPowerSums
  };

  enum class InnerProductRoute
  {
    ViaRegisteredDiagonalPairing,
    ViaDualBasisCoefficient,
    ViaSchurCompleteKostkaNumbers,
    ViaSchurElementaryConjugateKostkaNumbers,
    ViaOmegaSchurCompleteConjugateKostkaNumbers,
    ViaOmegaSchurElementaryKostkaNumbers,
    ViaWeightedSchurCharacters,
    ViaPowerSumDiagonalPairing,
    ViaConvertBothToPowerSums
  };

  enum class InnerProductOrientation { Original, Swapped };

  struct InnerProductCandidate
  {
    InnerProductPipeline pipeline;
    InnerProductRoute route;
    InnerProductOrientation orientation;
    int coefficientBasisId = -1;
    int pairingSourceBasisId = -1;
    int pairingDualBasisId = -1;
    int pairingKind = 0;
    size_t estimatedCost = 0;
    KnownState transitionCached = KnownState::Unknown;
  };

  KnownState innerProductTransitionCacheState(
      const InnerProductCandidate& candidate,
      const InnerProductRequest& request) const;
  size_t estimateInnerProductCandidateCost(
      InnerProductCandidate& candidate,
      const InnerProductRequest& request) const;
  InnerProductCandidate selectInnerProductCandidate(
      std::vector<InnerProductCandidate> candidates,
      const InnerProductRequest& request) const;
  const char *innerProductRouteName(InnerProductRoute route) const;
  void traceInnerProductRouteSelection(
      const InnerProductCandidate& candidate) const;
  ring_elem executeInnerProductRoute(
      const InnerProductCandidate& candidate,
      const InnerProductRequest& request) const;

  // ============================================================================
  // Top-Level Inner-Product Pipeline
  // ============================================================================

  InnerProductPipeline selectInnerProductPipeline(
      const InnerProductRequest& request) const;
  const char *innerProductPipelineName(InnerProductPipeline pipeline) const;
  void traceInnerProductPipelineSelection(
      InnerProductPipeline pipeline,
      const InnerProductRequest& request) const;
  ring_elem executeInnerProductPipeline(
      InnerProductPipeline pipeline,
      const InnerProductRequest& request) const;
  ring_elem hallInnerProductDispatch(const InnerProductRequest& request) const;

  // ============================================================================
  // Diagonal-Basis Pipeline
  // ============================================================================

  InnerProductCandidate selectDiagonalBasisRoute(
      const InnerProductRequest& request) const;
  ring_elem runDiagonalBasisPipeline(
      const InnerProductRequest& request) const;

  // ============================================================================
  // Single-Basis-Element Pipeline
  // ============================================================================

  InnerProductCandidate selectSingleBasisElementRoute(
      const InnerProductRequest& request) const;
  ring_elem runSingleBasisElementPipeline(
      const InnerProductRequest& request) const;

  // ============================================================================
  // Structured Power-Sum Pipeline
  // ============================================================================

  InnerProductCandidate selectPowerSumsStructuredRoute(
      const InnerProductRequest& request) const;
  ring_elem runPowerSumsStructuredPipeline(
      const InnerProductRequest& request) const;

  // ============================================================================
  // Fallback Power-Sum Pipeline
  // ============================================================================

  ring_elem runFallbackPowerSumsPipeline(
      const InnerProductRequest& request) const;

  // ============================================================================
  // Public Inner-Product Entry Points
  // ============================================================================

  std::map<int, InnerProductTarget> innerProductTargetMap(
      M2_arrayint innerProductMap) const;
  ring_elem hallInnerProductElements(
      ring_elem f,
      ring_elem g,
      InnerProductKind kind,
      const std::map<int, InnerProductTarget>& metadata = {}) const;

 public:
  ring_elem hallInnerProduct(ring_elem f,
                             ring_elem g,
                             int kindCode,
                             M2_arrayint innerProductMap) const;

 private:

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
