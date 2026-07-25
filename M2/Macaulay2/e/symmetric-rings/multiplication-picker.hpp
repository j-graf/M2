// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_MULTIPLICATION_PICKER_HPP_
#define M2_SYMMETRIC_RINGS_MULTIPLICATION_PICKER_HPP_

// Declaration fragment included inside SymmetricEngineRing.

// ============================================================================
// Commutative Binary Endpoints
// ============================================================================

  struct UnorderedBasisPair
  {
    BasisKind first = BasisKind::Custom;
    BasisKind second = BasisKind::Custom;

    static UnorderedBasisPair from(
        BasisKind left, BasisKind right);

    bool operator<(const UnorderedBasisPair& other) const
    {
      return std::tie(first, second) <
             std::tie(other.first, other.second);
    }

    bool operator==(const UnorderedBasisPair& other) const
    {
      return first == other.first && second == other.second;
    }

    bool operator!=(const UnorderedBasisPair& other) const
    {
      return !(*this == other);
    }
  };

// ============================================================================
// Inspectable Strict-Pair Conditions
// ============================================================================

  enum class BinaryMultiplicationConditionKind
  {
    Always,
    Otherwise,
    BasisIndexHasAtMostOnePart,
    CombinedWeightAtMost,
    And,
    Or,
    Not
  };

  struct BinaryMultiplicationCondition
  {
    BinaryMultiplicationConditionKind kind =
        BinaryMultiplicationConditionKind::Always;
    BasisKind namedBasis = BasisKind::Custom;
    int integer = 0;
    std::vector<BinaryMultiplicationCondition> operands;
  };

  struct BinaryMultiplicationFacts
  {
    const CanonicalBasisTermView *left = nullptr;
    const CanonicalBasisTermView *right = nullptr;
    UnorderedBasisPair factorBases;
    BasisKind targetBasis = BasisKind::Custom;

    int combinedWeight() const;
    const Partition *uniqueIndexInBasis(BasisKind basis) const;
  };

  static BinaryMultiplicationCondition binaryAlways();
  static BinaryMultiplicationCondition binaryOtherwise();
  static BinaryMultiplicationCondition completeIndexIsOneRow();
  static BinaryMultiplicationCondition elementaryIndexIsOneColumn();
  static BinaryMultiplicationCondition powerSumIndexIsSingleCycle();
  static BinaryMultiplicationCondition combinedWeightAtMost(int weight);
  static BinaryMultiplicationCondition binaryAnd(
      BinaryMultiplicationCondition left,
      BinaryMultiplicationCondition right);
  static BinaryMultiplicationCondition binaryOr(
      BinaryMultiplicationCondition left,
      BinaryMultiplicationCondition right);
  static BinaryMultiplicationCondition binaryNot(
      BinaryMultiplicationCondition condition);
  static bool binaryMultiplicationConditionHolds(
      const BinaryMultiplicationCondition& condition,
      const BinaryMultiplicationFacts& facts);
  static std::string binaryMultiplicationConditionToString(
      const BinaryMultiplicationCondition& condition);
  static void validateBinaryMultiplicationCondition(
      const BinaryMultiplicationCondition& condition,
      const UnorderedBasisPair& endpoint,
      bool allowOtherwise = false);

// ============================================================================
// Kernel Declarations And Performance Picker
// ============================================================================

  struct BinaryMultiplicationKernelDefinition
  {
    std::string identifier;
    BasisKind firstArgumentBasis = BasisKind::Custom;
    BasisKind secondArgumentBasis = BasisKind::Custom;
    BasisKind targetBasis = BasisKind::Custom;
    BinaryMultiplicationCondition applicableWhen;
    BinaryMultiplicationKernel kernel = nullptr;
  };

  struct BinaryMultiplicationPreference
  {
    BinaryMultiplicationCondition condition;
    std::vector<std::string> kernelsInOrder;
    mutable std::vector<
        const BinaryMultiplicationKernelDefinition *> kernelDefinitions;
  };

  struct BinaryMultiplicationPicker
  {
    UnorderedBasisPair factorBases;
    BasisKind targetBasis = BasisKind::Custom;
    std::vector<BinaryMultiplicationPreference> preferences;
    std::vector<std::string> alternativeKernels;
    mutable std::vector<
        const BinaryMultiplicationKernelDefinition *>
        alternativeKernelDefinitions;
  };

  struct SelectedBinaryMultiplicationKernel
  {
    const BinaryMultiplicationKernelDefinition *definition = nullptr;
    RingBasis target;
    bool callableTakesReversedArguments = false;

    bool valid() const { return definition != nullptr; }
  };

  struct BinaryMultiplicationRequest
  {
    std::optional<std::string> forcedKernel;
    bool usePowerSumReference = false;
    bool traceWorkflow = false;

    bool automatic() const
    {
      return !forcedKernel && !usePowerSumReference;
    }
  };

  static const std::vector<BinaryMultiplicationKernelDefinition>&
  binaryMultiplicationKernelDatabase();
  static const std::vector<BinaryMultiplicationPicker>&
  binaryMultiplicationPickerDatabase();
  static const std::map<
      std::string,
      const BinaryMultiplicationKernelDefinition *>&
  binaryMultiplicationKernelsByIdentifier();
  static const std::map<
      std::pair<UnorderedBasisPair, BasisKind>,
      const BinaryMultiplicationPicker *>&
  binaryMultiplicationPickersByEndpoint();
  static void validateBinaryMultiplicationDatabase();

  BinaryMultiplicationFacts binaryMultiplicationFacts(
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right,
      BasisKind targetBasis) const;
  std::optional<SelectedBinaryMultiplicationKernel>
  selectBinaryMultiplicationKernel(
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right,
      int targetBasisId,
      const BinaryMultiplicationRequest& request,
      bool traceSelection) const;
  bool automaticBinaryMultiplicationKernelAvailable(
      const CanonicalBasisTermView& left,
      const CanonicalBasisTermView& right,
      int targetBasisId) const;

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
