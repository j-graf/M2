-- Verify the formulas in plethysm stability theorem 2.
--
-- The calculations below are driven entirely by the three input partitions.
-- For a multirow mu they verify part (1), including the asserted upper bound.
-- For a one-row mu they verify both descriptions in part (2), coefficientwise
-- through the user-selected value rowMMax.

restart
needsPackage "SymmetricRings"

A = frac(QQ[t,z])
Sym = symmetricRing A

stability2qBasis = basis(Sym, "HallLittlewoodQGenerator")
stability2QBasis = basis(Sym, "HallLittlewoodQ")

stability2q = index -> stability2qBasis_index
stability2Q = index -> stability2QBasis_index

-- The parentheses in the theorem denote the Hall--Littlewood scalar product.
-- Never rely on the package's context-sensitive "Automatic" choice here.
stability2Pair = (f, g) ->
    hallInnerProduct(f, g, "InnerProduct" => "HallLittlewood")

stability2Plethysm = (f, g) ->
    plethysm(f, g)

stability2Range = (lower, upper) ->
    if lower > upper then {} else toList(lower..upper)

stability2RequireEqual = (label, left, right) -> (
    if left != right then
        error(label
            | "\nleft  = " | toString left
            | "\nright = " | toString right
            | "\ndifference = " | toString(left-right));
    )

-- Print one comparison per line.  Large coefficient polynomials are much
-- easier to inspect this way than in the pretty-printed matrix for a list.
stability2PrintRows = (label, headings, rows) -> (
    print(label | ": " | headings);
    scan(rows, row -> print("  " | toString row));
    )

-- Each row is {exponent, coefficient}.  Printing the terms separately avoids
-- the single rational-function fraction used by frac(QQ[t,z]).
stability2PrintLaurentPolynomial = (label, variableName, rows) -> (
    print(label | ", grouped by powers of " | variableName | ":");
    scan(rows, row ->
        print("  (" | toString(row#1) | ")*"
            | variableName | "^(" | toString(row#0) | ")"));
    )

-- All weak compositions of total weight at most bound with the requested
-- number of coordinates.  In particular, this adapts to ell(lam).
stability2WeakCompositionsOfWeight = (slots, total) -> (
    result := {};
    if slots == 0 then (
        if total == 0 then result = {{}};
        return result;
        );
    scan(stability2Range(0, total), firstPart ->
        scan(stability2WeakCompositionsOfWeight(slots - 1, total - firstPart),
            tail -> result = append(result, prepend(firstPart, tail))));
    result
    )

stability2WeakCompositionsAtMost = (slots, bound) -> (
    result := {};
    scan(stability2Range(0, bound), total ->
        result = result | stability2WeakCompositionsOfWeight(slots, total));
    result
    )

stability2NonzeroLength = gamma -> #select(gamma, part -> part != 0)

-- These are plethystic alphabet substitutions.  Thus p_r[X+1]=p_r+1 and
-- p_r[X-1]=p_r-1, as required by the theorem.
stability2AtXPlusOne = f ->
    stability2Plethysm(f, p_1 + 1)
stability2AtXMinusOne = f ->
    stability2Plethysm(f, p_1 - 1)

-- The coefficient q_r[1-t], returned in the coefficient field rather than as
-- a scalar symmetric function.  Negative subscripts use q_r=0.
stability2qAtOneMinusT = r -> (
    if r < 0 then return 0_A;
    lift(stability2Plethysm(stability2q({r}), (1-t)*1_Sym), A)
    )

-- Ordinary generalized binomial coefficients.  They give an independent
-- expansion of ((1-tu)/(1-u))^(1-t).
stability2GeneralizedBinomial = (a, r) -> (
    if r < 0 then return 0_A;
    if r == 0 then return 1_A;
    value := 1_A;
    scan(0..r-1, i -> value = value * (a-i)/(i+1));
    value
    )

stability2ClosedFactorCoefficient = r -> (
    if r < 0 then return 0_A;
    value := 0_A;
    scan(0..r, a ->
        value = value
            + stability2GeneralizedBinomial(1-t, a) * (-t)^a
            * stability2GeneralizedBinomial(t-1, r-a) * (-1)^(r-a));
    value
    )

stability2VerifyMultirow = (lam, mu, nu, boundaryChecks) -> (
    lamWeight := partitionWeight lam;
    muWeight := partitionWeight mu;
    nuWeight := partitionWeight nu;
    muFirstPart := mu#0;
    lowerM := -lamWeight;
    upperM := floor(nuWeight/(muWeight-muFirstPart)) - lamWeight;
    qMuAtXPlusOne := stability2AtXPlusOne stability2Q(mu);
    qNuAtXMinusOne := stability2AtXMinusOne stability2Q(nu);

    theoremPolynomial := 0_A;
    shiftedRows := {};
    scan(stability2Range(lowerM, upperM), m -> (
        n := (m+lamWeight)*muWeight-nuWeight;
        shiftedCoefficient := stability2Pair(
            stability2Plethysm(
                stability2Q(prepend(m, lam)), qMuAtXPlusOne),
            qNuAtXMinusOne);
        theoremPolynomial = theoremPolynomial
            + shiftedCoefficient*z^n;
        shiftedRows = append(
            shiftedRows, {m, n, shiftedCoefficient});
        ));

    -- Compare the definition with the theorem on its asserted support.
    definitionPolynomial := 0_A;
    definitionRows := {};
    scan(stability2Range(lowerM, upperM), m -> (
        n := (m+lamWeight)*muWeight-nuWeight;
        coefficient := stability2Pair(
            stability2Plethysm(
                stability2Q(prepend(m, lam)), stability2Q(mu)),
            stability2Q(prepend(n, nu)));
        definitionPolynomial = definitionPolynomial + coefficient*z^n;
        definitionRows = append(
            definitionRows, {m, n, coefficient});
        ));

    stability2RequireEqual(
        "part (1): defining Laurent polynomial differs from the X+1/X-1 formula"
            | "; lam=" | toString lam
            | ", mu=" | toString mu | ", nu=" | toString nu,
        definitionPolynomial,
        theoremPolynomial);
    print("part (1) Q_mu[X+1]: " | toString qMuAtXPlusOne);
    print("part (1) Q_nu[X-1]: " | toString qNuAtXMinusOne);
    stability2PrintRows(
        "part (1) coefficients",
        "{m, n, defining coefficient, X+1/X-1 coefficient, difference}",
        apply(#definitionRows, i -> {
            definitionRows#i#0,
            definitionRows#i#1,
            definitionRows#i#2,
            shiftedRows#i#2,
            definitionRows#i#2-shiftedRows#i#2
            }));
    stability2PrintLaurentPolynomial(
        "part (1) g(z;t)",
        "z",
        apply(shiftedRows, row -> {row#1, row#2}));

    -- A finite experiment cannot inspect infinitely many zero coefficients.
    -- Probe a configurable margin beyond each end of the asserted support,
    -- without mixing those probes into the polynomial comparison above.
    lowerSupportRows := {};
    scan(stability2Range(lowerM-boundaryChecks, lowerM-1), m -> (
        n := (m+lamWeight)*muWeight-nuWeight;
        coefficient := stability2Pair(
            stability2Plethysm(
                stability2Q(prepend(m, lam)), stability2Q(mu)),
            stability2Q(prepend(n, nu)));
        stability2RequireEqual(
            "part (1): coefficient below m=-|lam| is nonzero"
                | "; m=" | toString m | ", n=" | toString n,
            coefficient,
            0_A);
        lowerSupportRows = append(
            lowerSupportRows, {m, n, coefficient});
        ));
    upperSupportRows := {};
    scan(stability2Range(upperM+1, upperM+boundaryChecks), m -> (
        n := (m+lamWeight)*muWeight-nuWeight;
        coefficient := stability2Pair(
            stability2Plethysm(
                stability2Q(prepend(m, lam)), stability2Q(mu)),
            stability2Q(prepend(n, nu)));
        stability2RequireEqual(
            "part (1): coefficient above M is nonzero"
                | "; m=" | toString m | ", n=" | toString n,
            coefficient,
            0_A);
        upperSupportRows = append(
            upperSupportRows, {m, n, coefficient});
        ));
    stability2PrintRows(
        "part (1) support probes below -|lam|",
        "{m, n, coefficient}",
        lowerSupportRows);
    stability2PrintRows(
        "part (1) support probes above M",
        "{m, n, coefficient}",
        upperSupportRows);
    print("part (1) passed for lam = " | toString lam
        | ", mu = " | toString mu | ", nu = " | toString nu);
    hashTable {
        "Case" => 1,
        "LowerM" => lowerM,
        "M" => upperM,
        "g" => theoremPolynomial
        }
    )

stability2VerifyOneRow = (lam, mu, nu, rowMMax, boundaryChecks) -> (
    lamWeight := partitionWeight lam;
    nuWeight := partitionWeight nu;
    pValue := mu#0;
    lowerM := -lamWeight;
    if rowMMax < lowerM then
        error "rowMMax must be at least -|lam|";

    -- A_p^+ and A_p are constructed from p=mu_1, not from the sample value.
    stability2RequireEqual(
        "part (2): a one-row Q_mu is not Q_p; mu=" | toString mu,
        stability2Q(mu),
        stability2Q({pValue}));
    aPlus := stability2q({pValue});
    scan(stability2Range(1, pValue-1), j ->
        aPlus = aPlus + (1-t)*stability2q({j}));
    aAlphabet := (1-t)*1_Sym + aPlus;

    -- This is one of the identities used in the theorem's second formula.
    shiftedAAlphabet := stability2AtXPlusOne(stability2Q(mu));
    stability2RequireEqual(
        "part (2): A_p != Q_(p)[X+1]; p=" | toString pValue,
        aAlphabet,
        shiftedAAlphabet);
    aAlphabetInQGenerators := toBasis(aAlphabet, stability2qBasis);
    shiftedAAlphabetInQGenerators := toBasis(
        shiftedAAlphabet, stability2qBasis);
    print("part (2) A_p formula: " | toString aAlphabet);
    stability2PrintRows(
        "part (2) A_p comparison in Hall--Littlewood q generators",
        "{formula, Q_(p)[X+1], difference}",
        {{aAlphabetInQGenerators, shiftedAAlphabetInQGenerators,
                aAlphabetInQGenerators-shiftedAAlphabetInQGenerators}});

    shiftedNu := stability2AtXMinusOne stability2Q(nu);
    gammas := stability2WeakCompositionsAtMost(#lam, lamWeight);
    cValues := new MutableHashTable;

    scan(stability2Range(-nuWeight, lamWeight), k -> (
        cK := 0_A;
        scan(gammas, gamma -> (
            gammaWeight := sum gamma;
            if 0 <= gammaWeight-k and gammaWeight-k <= nuWeight then (
                gammaLength := stability2NonzeroLength gamma;
                lamMinusGamma := apply(#lam,
                    i -> lam#i-gamma#i);
                cK = cK
                    + (-1)^gammaLength
                    * (1-t)^gammaLength
                    * t^(gammaWeight-gammaLength)
                    * stability2Pair(
                        stability2Plethysm(
                            stability2q({gammaWeight-k}), aPlus)
                            * stability2Plethysm(
                                stability2Q(lamMinusGamma), aAlphabet),
                        shiftedNu);
                );
            ));
        cValues#k = cK;
        ));
    stability2PrintRows(
        "part (2) c_k(t)",
        "{k, c_k(t)}",
        apply(stability2Range(-nuWeight, lamWeight),
            k -> {k, cValues#k}));

    -- Independently check the coefficients of the formal factor
    -- ((1-tu)/(1-u))^(1-t) needed by this truncation.
    factorDegree := rowMMax + lamWeight;
    factorRows := {};
    scan(stability2Range(0, factorDegree), r -> (
        qValue := stability2qAtOneMinusT r;
        closedValue := stability2ClosedFactorCoefficient r;
        stability2RequireEqual(
            "part (2): q_r[1-t] differs from the closed factor"
                | "; r=" | toString r,
            qValue,
            closedValue);
        factorRows = append(factorRows, {r, qValue, closedValue});
        ));
    stability2PrintRows(
        "part (2) closed factor coefficients",
        "{r, q_r[1-t], closed coefficient, difference}",
        apply(factorRows, row ->
            row | {row#1-row#2}));

    cLaurentPolynomial := 0_A;
    scan(stability2Range(-nuWeight, lamWeight), k ->
        cLaurentPolynomial =
            cLaurentPolynomial + cValues#k*z^(-k));
    closedFactorTruncation := 0_A;
    scan(stability2Range(0, factorDegree), r ->
        closedFactorTruncation = closedFactorTruncation
            + stability2ClosedFactorCoefficient(r)*z^r);
    stability2PrintLaurentPolynomial(
        "part (2) C(u;t)",
        "u",
        apply(stability2Range(-nuWeight, lamWeight),
            k -> {-k, cValues#k}));
    stability2PrintLaurentPolynomial(
        "part (2) closed factor truncation",
        "u",
        apply(factorRows, row -> {row#0, row#2}));

    -- The theorem starts at m=-|lam|.  Probe a configurable number of lower
    -- indices to verify that no earlier term of the defining double sum occurs.
    lowerSupportRows := {};
    scan(stability2Range(lowerM-boundaryChecks, lowerM-1), m -> (
        n := (m+lamWeight)*pValue-nuWeight;
        coefficient := stability2Pair(
            stability2Plethysm(
                stability2Q(prepend(m, lam)), stability2Q(mu)),
            stability2Q(prepend(n, nu)));
        stability2RequireEqual(
            "part (2): coefficient below m=-|lam| is nonzero"
                | "; m=" | toString m | ", n=" | toString n,
            coefficient,
            0_A);
        lowerSupportRows = append(
            lowerSupportRows, {m, n, coefficient});
        ));
    stability2PrintRows(
        "part (2) support probes below -|lam|",
        "{m, n, coefficient}",
        lowerSupportRows);

    definitionTruncation := 0_A;
    coefficientFormulaTruncation := 0_A;
    closedFormTruncation := 0_A;
    coefficientRows := {};
    scan(stability2Range(lowerM, rowMMax), m -> (
        n := (m+lamWeight)*pValue-nuWeight;
        directCoefficient := stability2Pair(
            stability2Plethysm(
                stability2Q(prepend(m, lam)), stability2Q(mu)),
            stability2Q(prepend(n, nu)));

        formulaCoefficient := 0_A;
        closedCoefficient := 0_A;
        scan(stability2Range(-nuWeight, lamWeight), k -> (
            formulaCoefficient = formulaCoefficient
                + cValues#k * stability2qAtOneMinusT(m+k);
            closedCoefficient = closedCoefficient
                + cValues#k * stability2ClosedFactorCoefficient(m+k);
            ));

        stability2RequireEqual(
            "part (2): defining coefficient differs from the c_k formula"
                | "; m=" | toString m | ", n=" | toString n,
            directCoefficient,
            formulaCoefficient);
        stability2RequireEqual(
            "part (2): c_k formula differs from the closed form"
                | "; m=" | toString m | ", n=" | toString n,
            formulaCoefficient,
            closedCoefficient);
        definitionTruncation = definitionTruncation
            + directCoefficient*z^n;
        coefficientFormulaTruncation = coefficientFormulaTruncation
            + formulaCoefficient*z^n;
        closedFormTruncation = closedFormTruncation
            + closedCoefficient*z^n;
        coefficientRows = append(coefficientRows, {
                m, n, directCoefficient,
                formulaCoefficient, closedCoefficient
                });
        ));

    stability2RequireEqual(
        "part (2): defining truncation differs from the c_k truncation",
        definitionTruncation,
        coefficientFormulaTruncation);
    stability2RequireEqual(
        "part (2): c_k truncation differs from the closed-form truncation",
        coefficientFormulaTruncation,
        closedFormTruncation);
    stability2PrintRows(
        "part (2) coefficients",
        "{m, n, defining coefficient, c_k coefficient, closed coefficient, "
            | "defining-c_k, c_k-closed}",
        apply(coefficientRows, row ->
            row | {row#2-row#3, row#3-row#4}));
    stability2PrintLaurentPolynomial(
        "part (2) g(z;t) truncation",
        "z",
        apply(coefficientRows, row -> {row#1, row#2}));
    print("part (2) passed through m = " | toString rowMMax
        | " for lam = " | toString lam
        | ", mu = " | toString mu | ", nu = " | toString nu);
    hashTable {
        "Case" => 2,
        "Lambda" => lam,
        "Mu" => mu,
        "Nu" => nu,
        "LowerM" => lowerM,
        "MMaxChecked" => rowMMax,
        "c" => new HashTable from cValues,
        "C(u), with u=z" => cLaurentPolynomial,
        "ClosedFactorTruncation, with u=z" => closedFactorTruncation,
        "gTruncation" => definitionTruncation
        }
    )

stability2Verify = (lam, mu, nu, rowMMax, boundaryChecks) -> (
    lam = take(lam, partitionLength lam);
    mu = take(mu, partitionLength mu);
    nu = take(nu, partitionLength nu);
    if partitionLength mu == 0 then
        error "mu must be a nonempty partition";
    if boundaryChecks < 0 then
        error "boundaryChecks must be nonnegative";
    if partitionLength mu > 1 then
        stability2VerifyMultirow(lam, mu, nu, boundaryChecks)
    else
        stability2VerifyOneRow(lam, mu, nu, rowMMax, boundaryChecks)
    )

-- ============================================================================
-- Corollary: binomial formula and integer specializations of P_t(m)
-- ============================================================================

stability2Differences = values -> (
    if #values < 2 then return {};
    apply(0..#values-2, i -> values#(i+1)-values#i)
    )

stability2DifferenceOrder = (values, order) -> (
    result := values;
    scan(stability2Range(1, order), i ->
        result = stability2Differences result);
    result
    )

stability2PrintValues = (label, mValues, values) -> (
    stability2PrintRows(
        label,
        "{m, value}",
        apply(#mValues, i -> {mValues#i, values#i}));
    )

stability2VerifyDegreeBound = (label, values, degree) -> (
    differences := stability2DifferenceOrder(values, degree+1);
    scan(0..#differences-1, i ->
        stability2RequireEqual(
            label | "; finite difference index=" | toString i,
            differences#i,
            0_A));
    )

-- Reconstructs the sampled polynomial in the Newton basis.  The coefficient
-- ring variable z is used here only as a printable stand-in for m.
stability2NewtonPolynomial = (mStart, values, degree) -> (
    result := 0_A;
    differences := values;
    scan(stability2Range(0, degree), order -> (
        result = result
            + differences#0
            * stability2GeneralizedBinomial(z-mStart, order);
        differences = stability2Differences differences;
        ));
    result
    )

stability2PValue = (cValues, lamWeight, nuWeight, m) -> (
    result := 0_A;
    scan(stability2Range(-nuWeight, lamWeight), k ->
        result = result
            + cValues#k * stability2qAtOneMinusT(m+k));
    result
    )

-- The coefficient formula in the corollary, before specializing t.
stability2BinomialPValue = (cValues, lamWeight, nuWeight, m) -> (
    result := 0_A;
    scan(stability2Range(1, m+lamWeight), r -> (
        innerSum := 0_A;
        scan(stability2Range(-nuWeight, lamWeight), k ->
            innerSum = innerSum
                + cValues#k
                * stability2GeneralizedBinomial(m+k-1, r-1));
        result = result
            + stability2GeneralizedBinomial(1-t, r)
            * (1-t)^r * innerSum;
        ));
    result
    )

-- Specialize only after the complete symbolic sum has been collected.  This
-- preserves cancellations at integer t-values where individual rational
-- intermediates may have removable singularities.
stability2SpecializedPValue =
    (cValues, lamWeight, nuWeight, m, tValue) ->
        sub(stability2PValue(cValues, lamWeight, nuWeight, m),
            {t => tValue})

stability2Factorial = n -> (
    if n == 0 then return 1_A;
    product toList(1..n)
    )

stability2SpecializedCMoment =
    (cValues, lamWeight, nuWeight, tValue, base) -> (
        result := 0_A;
        scan(stability2Range(-nuWeight, lamWeight), k ->
            result = result
                + sub(cValues#k, {t => tValue}) * base^k);
        result
        )

-- The explicit polynomial P_(-d)(m), with z used as the variable m.
stability2NegativePolynomial =
    (cValues, lamWeight, nuWeight, d) -> (
        result := 0_A;
        scan(stability2Range(1, d+1), r -> (
            innerSum := 0_A;
            scan(stability2Range(-nuWeight, lamWeight), k ->
                innerSum = innerSum
                    + sub(cValues#k, {t => -d})
                    * stability2GeneralizedBinomial(z+k-1, r-1));
            result = result
                + stability2GeneralizedBinomial(d+1, r)
                * (d+1)^r * innerSum;
            ));
        result
        )

-- The explicit polynomial R_d(m), with z used as the variable m.
stability2PositivePolynomial =
    (cValues, lamWeight, nuWeight, d) -> (
        result := 0_A;
        scan(stability2Range(1, d-1), r -> (
            innerSum := 0_A;
            scan(stability2Range(-nuWeight, lamWeight), k ->
                innerSum = innerSum
                    + sub(cValues#k, {t => d}) * d^k
                    * stability2GeneralizedBinomial(z+k-1, r-1));
            result = result
                + stability2GeneralizedBinomial(d-1, r)
                * ((d-1)/d)^r * innerSum;
            ));
        result
        )

stability2MaximumOrZero = values ->
    if #values == 0 then 0 else max values

stability2VerifyCorollary = (
    lam, mu, nu, theoremResult,
    nonnegativeDs, positiveDs, extraChecks) -> (
    lam = take(lam, partitionLength lam);
    mu = take(mu, partitionLength mu);
    nu = take(nu, partitionLength nu);
    if partitionLength mu != 1 then
        error "the corollary requires mu=(p)";
    if extraChecks < 0 then
        error "extraChecks must be nonnegative";
    if not all(nonnegativeDs, d -> instance(d, ZZ) and d >= 0) then
        error "expected nonnegative integers in nonnegativeDs";
    if not all(positiveDs, d -> instance(d, ZZ) and d >= 2) then
        error "expected integers at least 2 in positiveDs";
    if theoremResult#"Case" != 2
        or theoremResult#"Lambda" != lam
        or theoremResult#"Mu" != mu
        or theoremResult#"Nu" != nu then
        error "theoremResult was computed from different partitions";

    lamWeight := partitionWeight lam;
    nuWeight := partitionWeight nu;
    mStart := nuWeight+1;
    cValues := theoremResult#"c";
    negativePolynomials := new MutableHashTable;
    positivePolynomials := new MutableHashTable;

    -- Verify the unspecialized binomial formula coefficientwise in t.
    generalMValues := stability2Range(mStart, mStart+extraChecks);
    generalRows := apply(generalMValues, m -> (
            pValue := stability2PValue(
                cValues, lamWeight, nuWeight, m);
            binomialValue := stability2BinomialPValue(
                cValues, lamWeight, nuWeight, m);
            stability2RequireEqual(
                "corollary binomial formula failed; m=" | toString m,
                pValue,
                binomialValue);
            {m, pValue, binomialValue, pValue-binomialValue}
            ));
    stability2PrintRows(
        "corollary binomial formula",
        "{m, P_t(m), binomial formula, difference}",
        generalRows);

    -- At t=-d, vanishing (d+1)-st finite differences verify degree at most d.
    scan(nonnegativeDs, d -> (
        mValues := stability2Range(
            mStart, mStart+d+1+extraChecks);
        values := apply(mValues, m ->
            stability2SpecializedPValue(
                cValues, lamWeight, nuWeight, m, -d));
        stability2VerifyDegreeBound(
            "corollary (1) failed for d=" | toString d,
            values,
            d);
        polynomial := stability2NewtonPolynomial(mStart, values, d);
        explicitPolynomial := stability2NegativePolynomial(
            cValues, lamWeight, nuWeight, d);
        stability2RequireEqual(
            "corollary (1) explicit polynomial failed for d="
                | toString d,
            polynomial,
            explicitPolynomial);
        explicitValues := apply(mValues, m ->
            sub(explicitPolynomial, {z => m}));
        scan(0..#values-1, i ->
            stability2RequireEqual(
                "corollary (1) explicit value failed for d="
                    | toString d | ", m=" | toString(mValues#i),
                values#i,
                explicitValues#i));

        factorialD := stability2Factorial d;
        leadingFromValues :=
            (stability2DifferenceOrder(values, d))#0/factorialD;
        leadingFormula := (d+1)^(d+1)/factorialD
            * stability2SpecializedCMoment(
                cValues, lamWeight, nuWeight, -d, 1);
        stability2RequireEqual(
            "corollary (1) leading coefficient failed for d="
                | toString d,
            leadingFromValues,
            leadingFormula);

        negativePolynomials#d = polynomial;
        stability2PrintRows(
            "corollary (1) explicit formula for P_"
                | toString(-d) | "(m)",
            "{m, P_(-d)(m), explicit formula, difference}",
            apply(#mValues, i -> {
                    mValues#i,
                    values#i,
                    explicitValues#i,
                    values#i-explicitValues#i
                    }));
        print("order " | toString(d+1) | " finite differences: "
            | toString stability2DifferenceOrder(values, d+1));
        print("explicit polynomial in m (displayed with m=z): "
            | toString explicitPolynomial);
        stability2PrintRows(
            "coefficient of m^" | toString d,
            "{from sampled polynomial, claimed formula, difference}",
            {{leadingFromValues, leadingFormula,
                    leadingFromValues-leadingFormula}});
        print("corollary (1) passed for d=" | toString d);
        ));

    -- Since m>|nu| and k>=-|nu|, every q_(m+k)[1-t] has positive
    -- subscript.  The complete P_1(m) must therefore vanish.
    oneChecks := stability2MaximumOrZero(nonnegativeDs)
        + stability2MaximumOrZero(positiveDs)
        + extraChecks + 2;
    oneMValues := stability2Range(mStart, mStart+oneChecks);
    oneValues := apply(oneMValues, m ->
        stability2SpecializedPValue(
            cValues, lamWeight, nuWeight, m, 1));
    scan(0..#oneValues-1, i ->
        stability2RequireEqual(
            "corollary (2) failed; m=" | toString(oneMValues#i),
            oneValues#i,
            0_A));
    stability2PrintValues(
        "P_1(m)",
        oneMValues,
        oneValues);
    print("corollary (2) passed through m="
        | toString(mStart+oneChecks));

    -- At t=d, divide by d^m and check that the (d-1)-st finite differences
    -- vanish, which is equivalent to degree at most d-2.
    scan(positiveDs, d -> (
        degree := d-2;
        mValues := stability2Range(
            mStart, mStart+degree+1+extraChecks);
        normalizedValues := apply(mValues, m ->
            stability2SpecializedPValue(
                cValues, lamWeight, nuWeight, m, d)/d^m);
        unnormalizedValues := apply(#mValues, i ->
            d^(mValues#i)*normalizedValues#i);
        stability2VerifyDegreeBound(
            "corollary (3) failed for d=" | toString d,
            normalizedValues,
            degree);
        polynomial := stability2NewtonPolynomial(
            mStart, normalizedValues, degree);
        explicitPolynomial := stability2PositivePolynomial(
            cValues, lamWeight, nuWeight, d);
        stability2RequireEqual(
            "corollary (3) explicit R_d failed for d=" | toString d,
            polynomial,
            explicitPolynomial);
        explicitValues := apply(mValues, m ->
            sub(explicitPolynomial, {z => m}));
        scan(0..#normalizedValues-1, i ->
            stability2RequireEqual(
                "corollary (3) explicit value failed for d="
                    | toString d | ", m=" | toString(mValues#i),
                normalizedValues#i,
                explicitValues#i));

        factorialDegree := stability2Factorial degree;
        leadingFromValues :=
            (stability2DifferenceOrder(normalizedValues, degree))#0
            / factorialDegree;
        -- Only r=d-1 contributes to m^(d-2); the exponent d-1 corrects
        -- the unbound r in the displayed corollary statement.
        leadingFormula := ((d-1)/d)^(d-1)/factorialDegree
            * stability2SpecializedCMoment(
                cValues, lamWeight, nuWeight, d, d);
        stability2RequireEqual(
            "corollary (3) leading coefficient failed for d="
                | toString d,
            leadingFromValues,
            leadingFormula);

        positivePolynomials#d = polynomial;
        stability2PrintValues(
            "P_" | toString d | "(m)",
            mValues,
            unnormalizedValues);
        stability2PrintRows(
            "corollary (3) explicit formula for R_"
                | toString d | "(m)",
            "{m, P_d(m)/d^m, explicit R_d(m), difference}",
            apply(#mValues, i -> {
                    mValues#i,
                    normalizedValues#i,
                    explicitValues#i,
                    normalizedValues#i-explicitValues#i
                    }));
        print("order " | toString(d-1) | " finite differences: "
            | toString stability2DifferenceOrder(
                normalizedValues, d-1));
        print("explicit R_" | toString d
            | "(m) (displayed with m=z): "
            | toString explicitPolynomial);
        stability2PrintRows(
            "coefficient of m^" | toString degree
                | " in R_" | toString d,
            "{from sampled polynomial, claimed formula, difference}",
            {{leadingFromValues, leadingFormula,
                    leadingFromValues-leadingFormula}});
        print("corollary (3) passed for d=" | toString d);
        ));

    hashTable {
        "NegativeSpecializationPolynomials, with m=z" =>
            new HashTable from negativePolynomials,
        "PositiveSpecializationPolynomials R_d, with m=z" =>
            new HashTable from positivePolynomials,
        "SmallestMChecked" => mStart,
        "ExtraChecks" => extraChecks
        }
    )

-- Change these partitions freely.  The verifier selects the applicable part
-- of the theorem from ell(mu); no weights, lengths, bounds, or exponents below
-- are specialized to these examples.
lam = {2}
nu = {3}
muMultirow = {2,1}
muOneRow = {2}

multirowResult = stability2Verify(lam, muMultirow, nu, 4, 2)
oneRowResult = stability2Verify(lam, muOneRow, nu, 4, 2)
corollaryResult = stability2VerifyCorollary(
    lam, muOneRow, nu, oneRowResult,
    {0,1,2}, {2,3,4}, 3)

multirowResult
oneRowResult
corollaryResult
