-- Named mathematical inputs for the systematic SymmetricRings benchmarks.
-- Keep identifiers stable: result histories and reports refer to them.

benchmarkPartitionRecord = (id, lambda, shapeClass) -> hashTable {
    "ID" => id,
    "Partition" => lambda,
    "ShapeClass" => shapeClass,
    "Weight" => sum lambda,
    "Length" => #lambda,
    "LargestPart" => if #lambda == 0 then 0 else lambda#0,
    "DistinctParts" => #unique lambda
    }

benchmarkPartitions = {
    benchmarkPartitionRecord("w06-row",              {6},             "row"),
    benchmarkPartitionRecord("w06-hook",             {4,1,1},         "hook"),
    benchmarkPartitionRecord("w06-two-balanced",     {3,3},           "two-row-balanced"),
    benchmarkPartitionRecord("w06-staircase",        {3,2,1},         "staircase"),
    benchmarkPartitionRecord("w08-row",              {8},             "row"),
    benchmarkPartitionRecord("w08-hook",             {5,1,1,1},       "hook"),
    benchmarkPartitionRecord("w08-two-balanced",     {4,4},           "two-row-balanced"),
    benchmarkPartitionRecord("w08-three-generic",    {4,3,1},         "three-row-generic"),
    benchmarkPartitionRecord("w08-repeated",         {3,3,2},         "repeated-parts"),
    benchmarkPartitionRecord("w10-row",              {10},            "row"),
    benchmarkPartitionRecord("w10-hook",             {7,1,1,1},       "hook"),
    benchmarkPartitionRecord("w10-two-unbalanced",   {7,3},           "two-row-unbalanced"),
    benchmarkPartitionRecord("w10-three-balanced",   {4,3,3},         "three-row-balanced"),
    benchmarkPartitionRecord("w10-four-generic",     {4,3,2,1},       "four-row-generic"),
    benchmarkPartitionRecord("w12-row",              {12},            "row"),
    benchmarkPartitionRecord("w12-hook",             {8,1,1,1,1},     "hook"),
    benchmarkPartitionRecord("w12-rectangle",        {4,4,4},         "rectangle"),
    benchmarkPartitionRecord("w12-staircase",        {5,4,2,1},       "near-staircase"),
    benchmarkPartitionRecord("w13-four-generic",     {5,4,3,1},       "four-row-generic"),
    benchmarkPartitionRecord("w13-three-generic",    {7,4,2},         "three-row-generic"),
    benchmarkPartitionRecord("w14-three-balanced",   {7,4,3},         "three-row-balanced"),
    benchmarkPartitionRecord("w14-three-unbalanced", {8,4,2},         "three-row-unbalanced"),
    benchmarkPartitionRecord("w14-repeated",         {5,5,2,2},       "repeated-parts"),
    benchmarkPartitionRecord("w16-two-balanced",     {8,8},           "two-row-balanced"),
    benchmarkPartitionRecord("w16-four-generic",     {7,5,3,1},       "four-row-generic"),
    benchmarkPartitionRecord("w18-three-generic",    {9,6,3},         "three-row-generic"),
    benchmarkPartitionRecord("w20-two-unbalanced",   {12,8},          "two-row-unbalanced"),
    benchmarkPartitionRecord("w20-three-generic",    {10,6,4},        "three-row-generic"),
    benchmarkPartitionRecord("w24-four-generic",     {9,7,5,3},       "four-row-generic"),
    benchmarkPartitionRecord("w28-three-generic",    {13,9,6},        "three-row-generic"),
    benchmarkPartitionRecord("w30-row",              {30},            "row"),
    benchmarkPartitionRecord("w30-three-generic",    {14,10,6},       "three-row-generic")
    }

benchmarkPartitionById = hashTable apply(benchmarkPartitions, x -> (x#"ID", x))

benchmarkPartition = id -> (
    if not benchmarkPartitionById#?id then error("unknown benchmark partition: ", id);
    benchmarkPartitionById#id
    )

benchmarkPairRecord = (id, left, right, pairClass, tier) -> hashTable {
    "ID" => id,
    "Left" => left,
    "Right" => right,
    "PairClass" => pairClass,
    "Tier" => tier
    }

-- Schur product pairs emphasize Littlewood-Richardson support profiles.
benchmarkProductPairs = {
    benchmarkPairRecord("prod-small-balanced", "w06-two-balanced", "w06-staircase", "balanced-by-staircase", "Small"),
    benchmarkPairRecord("prod-hook-balanced", "w08-hook", "w08-two-balanced", "hook-by-balanced", "Small"),
    benchmarkPairRecord("prod-three-three", "w08-three-generic", "w08-repeated", "three-row-by-repeated", "Small"),
    benchmarkPairRecord("prod-four-three", "w10-four-generic", "w10-three-balanced", "four-row-by-three-row", "Medium"),
    benchmarkPairRecord("prod-unequal-weight", "w12-rectangle", "w08-hook", "unequal-weight", "Medium"),
    benchmarkPairRecord("prod-large-two", "w16-two-balanced", "w10-two-unbalanced", "large-two-row", "Large"),
    benchmarkPairRecord("prod-large-three", "w13-three-generic", "w14-three-balanced", "large-three-row", "Large"),
    benchmarkPairRecord("prod-large-four", "w16-four-generic", "w14-repeated", "large-four-row", "Stress")
    }

-- Plethysm pairs cover determinant size, Adams dilation, and multirow inners.
benchmarkPlethysmPairs = {
    benchmarkPairRecord("pleth-row-row", "w08-row", "inner-row3", "row-by-row", "Large"),
    benchmarkPairRecord("pleth-two-row2", "w10-two-unbalanced", "inner-row2", "two-row-by-row", "Medium"),
    benchmarkPairRecord("pleth-three-row2", "w13-three-generic", "inner-row2", "three-row-by-row", "Medium"),
    benchmarkPairRecord("pleth-three-balanced-row2", "w14-three-balanced", "inner-row2", "three-balanced-by-row", "Large"),
    benchmarkPairRecord("pleth-four-row2", "w13-four-generic", "inner-row2", "four-row-by-row", "Large"),
    benchmarkPairRecord("pleth-two-multirow", "w10-two-unbalanced", "inner-21", "two-row-by-multirow", "Large"),
    benchmarkPairRecord("pleth-three-multirow", "w08-three-generic", "inner-31", "three-row-by-multirow", "Large"),
    benchmarkPairRecord("pleth-hook-row", "w08-hook", "inner-row2", "hook-by-row", "Medium"),
    benchmarkPairRecord("pleth-repeated-row", "w08-repeated", "inner-row2", "repeated-by-row", "Medium")
    }

-- Dedicated small one-row inner indices make crossover tables readable.
benchmarkPlethysmInnerRows = {
    benchmarkPartitionRecord("inner-row2", {2}, "row"),
    benchmarkPartitionRecord("inner-row3", {3}, "row"),
    benchmarkPartitionRecord("inner-row4", {4}, "row"),
    benchmarkPartitionRecord("inner-row6", {6}, "row"),
    benchmarkPartitionRecord("inner-row8", {8}, "row"),
    benchmarkPartitionRecord("inner-21", {2,1}, "two-row"),
    benchmarkPartitionRecord("inner-31", {3,1}, "two-row")
    }

benchmarkPartitionById = hashTable apply(
    benchmarkPartitions | benchmarkPlethysmInnerRows, x -> (x#"ID", x))
