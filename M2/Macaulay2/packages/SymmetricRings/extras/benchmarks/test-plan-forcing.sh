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

assert_forced_conversion_agrees \
    'p->S:grouped-characters' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{4,2}+p_{3,2,1},S); exit 0'
assert_forced_conversion_agrees \
    'p->S:border-strips' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{4,2}+p_{3,2,1},S); exit 0'
assert_forced_conversion_agrees \
    'p->S:abacus-rim-hooks' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms toBasis(p_{4,2}+p_{3,2,1},S); exit 0'
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
    'single-cycle-p-terms->Hall-Littlewood:Green-polynomials' \
    'needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; print toExternalString rawTerms toBasis(p_3+2*p_2+3*p_1,Q); exit 0'

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
    'monomial-like-product:exponent-splittings' \
    'needsPackage "SymmetricRings"; R=symmetricRing QQ; print toExternalString rawTerms multiplyToBasis(m_2,m_1,m); exit 0'
assert_forced_multiplication_agrees \
    'Hall-Littlewood-product:generators' \
    'needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; print toExternalString rawTerms multiplyToBasis(Q_2,Q_1,Q); exit 0'

M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN='unknown-plan' \
    run_m2 'needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(toBasis(p_2,S); false) else true); exit 0'
M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN='Schur-product:horizontal-pieri' \
    run_m2 'needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(multiplyToBasis(S_3,S_2,S); false) else true); exit 0'

printf '%s\n' 'Conversion and multiplication plan-selection tests passed'
