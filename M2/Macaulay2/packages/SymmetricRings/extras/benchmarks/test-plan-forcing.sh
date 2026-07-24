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

assert_forced_conversion_agrees()
{
    plan=$1
    code=$2
    expected=$(run_m2 "$code")
    actual=$(
        M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN=$plan \
            run_m2 "$code"
    )
    if test "$actual" != "$expected"
    then
        printf '%s\n' "forced conversion disagrees: $plan" >&2
        return 1
    fi
}

assert_forced_trace_count()
{
    plan=$1
    expected_count=$2
    expected_text=$3
    code=$4
    output=$(
        M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN=$plan \
        M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1 \
            run_m2 "$code" 2>&1
    )
    actual_count=$(
        printf '%s\n' "$output" |
            grep -F -c -- "$expected_text" || true
    )
    if test "$actual_count" -ne "$expected_count"
    then
        printf '%s\n' \
            "unexpected trace count for $plan: $expected_text" >&2
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
    'toBasis: bypass=canonical-metadata target=S' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=S_{3,2,1,1}@S_2; exit 0'
assert_trace \
    'basis-coefficient: target=S_2 route=full-basis-conversion' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=basisCoefficient(S_1*S_1,S_2); exit 0'
assert_trace \
    'conversion-plan: source=PowerSum target=Schur plan=PowerSum->Schur:default-policy' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasis(p_3+2*p_2+3*p_1,S); exit 0'
assert_trace \
    'conversion-plan: source=Elementary target=Schur plan=Elementary->Schur:via-PowerSum-default' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasis(e_3+2*e_2+3*e_1,S); exit 0'
assert_trace \
    'conversion-plan: source=Schur target=Complete plan=Schur->Complete:Jacobi-Trudi' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasis(S_{4,2}+S_{3,2,1},h); exit 0'
assert_trace \
    'conversion-plan: source=SchurOmega target=Elementary plan=SchurOmega->Elementary:Jacobi-Trudi' \
    'needsPackage "SymmetricRings"; R=symmetricRing(QQ,"NormalizeSomega"=>false); A=toBasis(Somega_{4,2}+Somega_{3,2,1},e); exit 0'

assert_forced_conversion_agrees \
    'PowerSum->Schur:grouped-characters' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{4,2}+p_{3,2,1},S); exit 0'
assert_forced_conversion_agrees \
    'Elementary->Schur:via-PowerSum-default' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(e_3+2*e_2+3*e_1,S); exit 0'
assert_forced_conversion_agrees \
    'Schur->Complete:Jacobi-Trudi' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(S_{4,2}+S_{3,2,1},h); exit 0'
assert_forced_conversion_agrees \
    'SchurOmega->Elementary:Jacobi-Trudi' \
    'needsPackage "SymmetricRings"; R=symmetricRing(QQ,"NormalizeSomega"=>false); print toExternalString rawTerms toBasis(Somega_{4,2}+Somega_{3,2,1},e); exit 0'
assert_forced_conversion_agrees \
    'PowerSum->Schur:border-strips' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{4,2}+p_{3,2,1},S); exit 0'
assert_forced_conversion_agrees \
    'PowerSum->Schur:abacus-rim-hooks' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{4,2}+p_{3,2,1},S); exit 0'
assert_forced_conversion_agrees \
    'PowerSum->Schur:via-complete' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{4,2}+p_{3,2,1},S); exit 0'
assert_forced_conversion_agrees \
    'PowerSum->Schur:complete-friendly-hybrid-plan' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{8,1}+p_{7,2}+p_{6,3}+p_{5,4}+p_{5,3,1}+p_{4,3,2}+p_{4,2,2,1}+p_{3,3,2,1},S); exit 0'
assert_forced_conversion_agrees \
    'PowerSum->SchurOmega:complete-friendly-hybrid-plan' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{8,1}+p_{7,2}+p_{6,3}+p_{5,4}+p_{5,3,1}+p_{4,3,2}+p_{4,2,2,1}+p_{3,3,2,1},Somega); exit 0'
assert_forced_trace_count \
    'PowerSum->Schur:complete-friendly-hybrid-plan' \
    1 \
    'plan=PowerSum->Complete:logarithm-formula' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasis(p_{8,1}+p_{7,2}+p_{6,3}+p_{5,4}+p_{5,3,1}+p_{4,3,2}+p_{4,2,2,1}+p_{3,3,2,1},S); exit 0'
assert_forced_trace_count \
    'PowerSum->Schur:default-policy' \
    1 \
    'plan=PowerSum->Schur:complete-friendly-hybrid-plan' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasis(p_{4,4,3,3}+p_{4,4,3,2,1}+p_{4,3,3,2,2}+p_{4,3,2,2,1,1,1}+p_{3,3,3,3,2}+p_{3,3,2,2,2,1,1}+p_{2,2,2,2,2,2,1,1}+p_{2,2,2,2,2,1,1,1,1}+p_14,S); exit 0'
assert_forced_conversion_agrees \
    'Complete->PowerSum:classical-formula' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(h_{3,1},p); exit 0'
assert_forced_conversion_agrees \
    'Elementary->PowerSum:classical-formula' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(e_{3,1},p); exit 0'
assert_forced_conversion_agrees \
    'Schur->PowerSum:characters' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(S_{3,1},p); exit 0'
assert_forced_conversion_agrees \
    'Monomial->PowerSum:transition' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(m_{3,1},p); exit 0'
assert_forced_conversion_agrees \
    'PowerSum->HallLittlewoodQ:single-cycles-Green' \
    'needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; print toExternalString rawTerms toBasis(p_3+2*p_2+3*p_1,Q); exit 0'
assert_forced_conversion_agrees \
    'Elementary->Monomial:via-PowerSum-default' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(e_{3,1}+e_2,m); exit 0'

# Exercise every applicable registered conversion plan, rather than only the
# automatically selected one. The small all-family matrix reaches every broad
# X -> p -> Y composition and every unconditional atomic competitor. The
# additional power-sum inputs cover conditional Hall--Littlewood plans,
# nonhomogeneous component routing, and the term-level Schur hybrid.
M2_SYMMETRIC_RINGS_CHECK_ALL_CONVERSION_PLANS=1 \
    run_m2 '
        needsPackage "SymmetricRings";
        E=frac(QQ[t]);
        R=symmetricRing E;
        smallBasisElements = {
            p_{2,1}, h_{2,1}, e_{2,1}, m_{2,1}, ff_{2,1},
            S_{2,1}, Somega_{2,1}, q_{2,1}, b_{2,1},
            Q_{2,1}, B_{2,1}, P_{2,1}, Pomega_{2,1}};
        allTargets = {p,h,e,m,ff,S,Somega,q,b,Q,B,P,Pomega};
        scan(smallBasisElements, F ->
            scan(allTargets, target -> toBasis(F,target)));
        scan({Q,B,P,Pomega}, target -> toBasis(p_2,target));
        toBasis(p_3+2*p_2+3*p_1,S);
        hybridInput =
            p_{4,4,3,3} + p_{4,4,3,2,1} +
            p_{4,3,3,2,2} + p_{4,3,2,2,1,1,1} +
            p_{3,3,3,3,2} + p_{3,3,2,2,2,1,1} +
            p_{2,2,2,2,2,2,1,1} +
            p_{2,2,2,2,2,1,1,1,1} + p_14;
        toBasis(hybridInput,S);
        toBasis(hybridInput,Somega);
        exit 0'

# Keep this list exhaustive over MultiplicationPlanId. Unlike conversion plans,
# multiplication catalogs depend on operand/target context and cannot be
# enumerated from one endpoint matrix.
assert_forced_multiplication_agrees \
    'Schur-product:littlewood-richardson' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,S_2,S); exit 0'
assert_forced_multiplication_agrees \
    'Schur-product:horizontal-pieri' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,h_2,S); exit 0'
assert_forced_multiplication_agrees \
    'Schur-product:vertical-pieri' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,e_2,S); exit 0'
assert_forced_multiplication_agrees \
    'Schur-product:border-strips' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,p_2,S); exit 0'
assert_forced_multiplication_agrees \
    'Schur-product:compatible-factor-rules' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(h_3,e_2,S); exit 0'
assert_forced_multiplication_agrees \
    'monomial-like-product:exponent-splittings' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(m_2,m_1,m); exit 0'
assert_forced_multiplication_agrees \
    'Hall-Littlewood-product:generators' \
    'needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; print toExternalString rawTerms multiplyToBasis(Q_2,Q_1,Q); exit 0'
assert_forced_multiplication_agrees \
    'product:multiplicative-target' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(h_3,h_2,h); exit 0'
assert_forced_multiplication_agrees \
    'product:power-sums' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(S_3,S_2,S); exit 0'

M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN='unknown-plan' \
    run_m2 'needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(toBasis(p_2,S); false) else true); exit 0'
M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN='Schur-product:horizontal-pieri' \
    run_m2 'needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(multiplyToBasis(S_3,S_2,S); false) else true); exit 0'

printf '%s\n' 'Conversion and multiplication plan-selection tests passed'
