#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../../../../.." && pwd)
M2_BIN=${M2_BIN:-"$SOURCE_ROOT/BUILD/build/M2"}
TEST_HOME=$(mktemp -d "${TMPDIR:-/tmp}/symmetricrings-plan-test.XXXXXX")
trap 'rm -rf "$TEST_HOME"' EXIT HUP INT TERM

run_m2()
{
    HOME="$TEST_HOME" "$M2_BIN" --no-preload --silent --stop -q -e "$1"
}

assert_trace()
{
    expected=$1
    code=$2
    output=$(
        M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1 \
            run_m2 "$code" 2>&1
    )
    if ! printf '%s\n' "$output" | grep -F -- "$expected" >/dev/null
    then
        printf '%s\n' "missing workflow trace: $expected" >&2
        printf '%s\n' "$output" >&2
        return 1
    fi
}

assert_bench_trace_count()
{
    expected_count=$1
    expected_text=$2
    code=$3
    output=$(
        run_m2 "$code" 2>&1
    )
    actual_count=$(
        printf '%s\n' "$output" |
            grep -F -c -- "$expected_text" || true
    )
    if test "$actual_count" -ne "$expected_count"
    then
        printf '%s\n' "unexpected benchmark trace count: $expected_text" >&2
        printf '%s\n' "$output" >&2
        return 1
    fi
}

assert_forced_multiplication_agrees()
{
    plan=$1
    code=$2
    expected=$(run_m2 "$code")
    actual=$(
        M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN=$plan \
            run_m2 "$code"
    )
    if test "$actual" != "$expected"
    then
        printf '%s\n' "forced multiplication disagrees: $plan" >&2
        return 1
    fi
}

assert_trace \
    'plethysm-to-basis: target=Schur route=S->S:Adams-Jacobi-Trudi' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=S_{4,2,1}@S_2; exit 0'
assert_trace \
    'plethysm-to-basis: target=Schur route=p-materialization->target:default-conversion' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=S_{3,2,1,1}@S_2; exit 0'
assert_trace \
    'basis-coefficient: target=S_2 route=full-basis-conversion' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=basisCoefficient(S_1*S_1,S_2); exit 0'
assert_bench_trace_count \
    1 \
    'conversion-plan: source=PowerSum target=Schur plan=PowerSum->Schur:homogeneous-component-formulas' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasisBench(p_3+2*p_2+3*p_1,S,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'conversion-plan: source=Elementary target=Schur plan=u->v:via-power-sums' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasisBench(e_3+2*e_2+3*e_1,S,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'conversion-plan: source=Schur target=Complete plan=Schur->Complete:Jacobi-Trudi' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasisBench(S_{4,2}+S_{3,2,1},h,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'conversion-plan: source=SchurOmega target=Elementary plan=SchurOmega->Elementary:Jacobi-Trudi' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing(QQ,"NormalizeSomega"=>false); A=toBasisBench(Somega_{4,2}+Somega_{3,2,1},e,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'conversion-plan: source=PowerSum target=HallLittlewoodQ plan=PowerSum->HallLittlewoodQ:single-cycles-via-Green-polynomials' \
    'debug needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; A=toBasisBench(p_3+2*p_2+3*p_1,Q,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'conversion-plan: source=PowerSum target=HallLittlewoodQ plan=PowerSum->HallLittlewoodQ:Green-polynomials-and-duality' \
    'debug needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; A=toBasisBench(p_{2,1},Q,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'conversion-plan: source=PowerSum target=HallLittlewoodQ plan=PowerSum->HallLittlewoodQ:triangular-reduction' \
    'debug needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; A=toBasisBench(p_{2,1}+p_{3,1},Q,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'

# Compare representative direct, composition, term-hybrid, component-hybrid,
# Hall--Littlewood, and generic power-sum plans through the private benchmark
# entry point. Each call verifies exact agreement with automatic selection.
run_m2 '
    debug needsPackage "SymmetricRings";
    comparePlans = (input, target, plans) -> (
        report := toBasisBench(
            input,
            target,
            "Plans" => prepend("Automatic", plans),
            "Repetitions" => 1,
            "Warmups" => 0);
        assert report#"Verified");

    R = symmetricRing QQ;
    schurInput = p_{4,2} + p_{3,2,1};
    comparePlans(
        schurInput,
        S,
        {
            "PowerSum->Schur:Frobenius-character-formula",
            "PowerSum->Schur:Murnaghan-Nakayama",
            "PowerSum->Schur:abacus-rim-hooks",
            "PowerSum->Schur:via-complete-basis"
            });
    comparePlans(
        S_{4,2} + S_{3,2,1},
        h,
        {"Schur->Complete:Jacobi-Trudi"});
    comparePlans(
        h_{3,1},
        p,
        {"Complete->PowerSum:Newton-identities"});
    comparePlans(
        e_{3,1},
        p,
        {"Elementary->PowerSum:Newton-identities"});
    comparePlans(
        S_{3,1},
        p,
        {"Schur->PowerSum:Frobenius-character-formula"});
    comparePlans(
        m_{3,1},
        p,
        {"Monomial->PowerSum:transition-matrix"});
    comparePlans(
        h_{3,1} + h_2,
        S,
        {"u->v:via-power-sums"});

    hybridInput =
        p_{8,1} + p_{7,2} + p_{6,3} + p_{5,4} +
        p_{5,3,1} + p_{4,3,2} + p_{4,2,2,1} +
        p_{3,3,2,1};
    comparePlans(
        hybridInput,
        S,
        {"PowerSum->Schur:short-cycle-hybrid"});
    comparePlans(
        hybridInput,
        Somega,
        {"PowerSum->SchurOmega:short-cycle-hybrid"});

    R = symmetricRing(QQ, "NormalizeSomega" => false);
    comparePlans(
        Somega_{4,2} + Somega_{3,2,1},
        e,
        {"SchurOmega->Elementary:Jacobi-Trudi"});

    E = frac(QQ[t]);
    R = symmetricRing E;
    comparePlans(
        p_3 + 2*p_2 + 3*p_1,
        Q,
        {"PowerSum->HallLittlewoodQ:single-cycles-via-Green-polynomials"});
    exit 0'

assert_bench_trace_count \
    1 \
    'plan=PowerSum->Complete:logarithm-formula' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; F=p_{8,1}+p_{7,2}+p_{6,3}+p_{5,4}+p_{5,3,1}+p_{4,3,2}+p_{4,2,2,1}+p_{3,3,2,1}; A=toBasisBench(F,S,"Plans"=>"PowerSum->Schur:short-cycle-hybrid","Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'plan=PowerSum->Schur:short-cycle-hybrid' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; F=p_{4,4,3,3}+p_{4,4,3,2,1}+p_{4,3,3,2,2}+p_{4,3,2,2,1,1,1}+p_{3,3,3,3,2}+p_{3,3,2,2,2,1,1}+p_{2,2,2,2,2,2,1,1}+p_{2,2,2,2,2,1,1,1,1}+p_14; A=toBasisBench(F,S,"Plans"=>"PowerSum->Schur:homogeneous-component-formulas","Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'

# Keep this list exhaustive over MultiplicationPlanId. Unlike conversion plans,
# multiplication plans depend on operand/target context and cannot be
# enumerated from one endpoint matrix.
assert_forced_multiplication_agrees \
    'Schur-product:Littlewood-Richardson' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,S_2,S); exit 0'
assert_forced_multiplication_agrees \
    'Schur-product:horizontal-Pieri' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,h_2,S); exit 0'
assert_forced_multiplication_agrees \
    'Schur-product:vertical-Pieri' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,e_2,S); exit 0'
assert_forced_multiplication_agrees \
    'Schur-product:Murnaghan-Nakayama' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,p_2,S); exit 0'
assert_forced_multiplication_agrees \
    'Schur-product:mixed-Pieri-and-Murnaghan-Nakayama' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(h_3,e_2,S); exit 0'
assert_forced_multiplication_agrees \
    'monomial-like-product:exponent-splittings' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(m_2,m_1,m); exit 0'
assert_forced_multiplication_agrees \
    'Hall-Littlewood-product:via-generator-basis' \
    'needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; print toExternalString rawTerms multiplyToBasis(Q_2,Q_1,Q); exit 0'
assert_forced_multiplication_agrees \
    'product:in-target-basis' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(h_3,h_2,h); exit 0'
assert_forced_multiplication_agrees \
    'product:via-power-sums' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,S_2,S); exit 0'

run_m2 'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(toBasisBench(p_2,S,"Plans"=>"unknown-plan","Repetitions"=>1,"Warmups"=>0); false) else true); exit 0'
run_m2 'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(toBasisBench(e_2,m,"Plans"=>"Elementary->Monomial:via-power-sums","Repetitions"=>1,"Warmups"=>0); false) else true); exit 0'
M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN='Schur-product:horizontal-Pieri' \
    run_m2 'needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(multiplyToBasis(S_3,S_2,S); false) else true); exit 0'

# Ordinary conversion is invariant under every former conversion diagnostic,
# including forced multiplication encountered while resolving a product.
ordinary_output=$(
    M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN='unknown-plan' \
    M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1 \
    M2_SYMMETRIC_RINGS_CHECK_ALL_CONVERSION_PLANS=1 \
    M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN='product:via-power-sums' \
        run_m2 'needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(toBasis(p_2,S)==S_2-S_{1,1}); assert(toBasis(S_2*S_2,S)==S_4+S_{3,1}+S_{2,2}); exit 0' \
        2>&1
)
if test -n "$ordinary_output"
then
    printf '%s\n' 'ordinary toBasis emitted diagnostic output' >&2
    printf '%s\n' "$ordinary_output" >&2
    exit 1
fi

# Contributor-surface invariants: mathematical kernels, complete plans, and
# endpoint picker policy have exactly one source owner each. Workflow
# infrastructure must not regain a kernel enum, execution switch, plan IDs, or
# a third formula variant for one-plan compositions.
if grep -F 'enum class BasisConversionKernel' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion.hpp" \
    >/dev/null ||
   grep -E \
    'KernelPlanFormula|NamedPlanFormula|ComposedPlansFormula|usePlan[[:space:]]*\(' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion.hpp" \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-plans.cpp" \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion.cpp" \
    >/dev/null ||
   grep -F 'executeBasisConversionKernel' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion.cpp" \
    >/dev/null ||
   grep -F 'canonicalBuiltInExpressionToPowerSumsViaBasisFormulas' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-plans.cpp" \
    >/dev/null ||
   grep -E 'default-policy|diagnosticOnly|powerSumExpressionTo' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-plans.cpp" \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-picker.cpp" \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-kernels.cpp" \
    >/dev/null ||
   grep -E \
    'powerSumIndexToTargetViaTermwiseKernel|powerSumsToTargetViaTermwiseConversion|schurUsesCompleteIntermediate|endpointSpecificPowerSumComposition' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-plans.cpp" \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-kernels.cpp" \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-kernels.hpp" \
    >/dev/null ||
   sed -n \
    '/Plan-Callable Conversion Kernels/,/Shared Conversion Utilities/p' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-kernels.cpp" |
    grep -E 'switch[[:space:]]*\(' >/dev/null ||
   grep -F 'basisExpansionToPowerSumsForNormalization' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-kernels.cpp" \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-kernels.hpp" \
    >/dev/null ||
   grep -E \
    'powerSumsToSchurVia|schurLikeToPowerSumsViaFrobeniusCharacterFormula' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/plethysm.cpp" \
    >/dev/null ||
   grep -E '"[A-Za-z]+->[A-Za-z]+' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion.cpp" \
    >/dev/null
then
    printf '%s\n' \
        'basis-conversion contributor-surface invariant failed' >&2
    exit 1
fi

printf '%s\n' 'Conversion and multiplication plan-selection tests passed'
