#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../../../../.." && pwd)
M2_BIN=${M2_BIN:-"$SOURCE_ROOT/BUILD/build/M2"}
BENCHMARK_HELPERS="$SCRIPT_DIR/benchmark-helpers.m2"
TEST_HOME=$(mktemp -d "${TMPDIR:-/tmp}/symmetricrings-plan-test.XXXXXX")
trap 'rm -rf "$TEST_HOME"' EXIT HUP INT TERM

run_m2()
{
    HOME="$TEST_HOME" "$M2_BIN" --no-preload --silent --stop -q \
        -e "debug needsPackage \"SymmetricRings\"; load \"$BENCHMARK_HELPERS\"; $1"
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

assert_binary_kernel_agrees()
{
    kernel=$1
    setup=$2
    left=$3
    right=$4
    target=$5
    run_m2 "debug needsPackage \"SymmetricRings\"; $setup; report=multiplyToBasisBench($left,$right,$target,\"Kernels\"=>{\"Automatic\",\"$kernel\",\"PowerSumReference\"},\"Repetitions\"=>1,\"Warmups\"=>0,\"Verify\"=>true); assert(report#\"Verified\"); exit 0"
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

# Force every declared kernel on a mathematically applicable pair and compare
# it with the independent power-sum reference.  This table is deliberately
# exhaustive: adding a kernel requires adding its differential case here.
run_m2 '
    debug needsPackage "SymmetricRings";
    R = symmetricRing QQ;
    kernelCases = {
        {"Schur*Schur->Schur:Littlewood-Richardson-tableaux", S_3, S_2, S},
        {"Schur*Schur->Schur:Littlewood-Richardson-coefficients", S_3, S_2, S},
        {"Schur*Complete->Schur:horizontal-Pieri", S_3, h_2, S},
        {"Schur*Complete->Schur:repeated-horizontal-Pieri", S_3, h_{2,1}, S},
        {"Schur*Elementary->Schur:vertical-Pieri", S_3, e_2, S},
        {"Schur*Elementary->Schur:repeated-vertical-Pieri", S_3, e_{2,1}, S},
        {"Schur*PowerSum->Schur:Murnaghan-Nakayama", S_3, p_2, S},
        {"Schur*PowerSum->Schur:repeated-Murnaghan-Nakayama", S_3, p_{2,1}, S},
        {"Complete*Complete->Schur:repeated-horizontal-Pieri", h_{2,1}, h_2, S},
        {"Complete*Elementary->Schur:repeated-Pieri", h_{2,1}, e_{2,1}, S},
        {"Complete*PowerSum->Schur:Pieri-and-Murnaghan-Nakayama", h_{2,1}, p_{2,1}, S},
        {"Elementary*Elementary->Schur:repeated-vertical-Pieri", e_{2,1}, e_2, S},
        {"Elementary*PowerSum->Schur:Pieri-and-Murnaghan-Nakayama", e_{2,1}, p_{2,1}, S},
        {"PowerSum*PowerSum->Schur:repeated-Murnaghan-Nakayama", p_{2,1}, p_2, S},

        {"Monomial*Monomial->Monomial:exponent-splittings", m_2, m_1, m},
        {"Monomial*Complete->Monomial:exponent-splittings", m_2, h_1, m},
        {"Monomial*Elementary->Monomial:exponent-splittings", m_2, e_1, m},
        {"Monomial*PowerSum->Monomial:exponent-splittings", m_2, p_1, m},
        {"Complete*Complete->Monomial:exponent-splittings", h_2, h_1, m},
        {"Complete*Elementary->Monomial:exponent-splittings", h_2, e_1, m},
        {"Complete*PowerSum->Monomial:exponent-splittings", h_2, p_1, m},
        {"Elementary*Elementary->Monomial:exponent-splittings", e_2, e_1, m},
        {"Elementary*PowerSum->Monomial:exponent-splittings", e_2, p_1, m},
        {"PowerSum*PowerSum->Monomial:exponent-splittings", p_2, p_1, m},

        {"Forgotten*Forgotten->Forgotten:exponent-splittings", ff_2, ff_1, ff},
        {"Forgotten*Complete->Forgotten:exponent-splittings", ff_2, h_1, ff},
        {"Forgotten*Elementary->Forgotten:exponent-splittings", ff_2, e_1, ff},
        {"Forgotten*PowerSum->Forgotten:exponent-splittings", ff_2, p_1, ff},
        {"Complete*Complete->Forgotten:exponent-splittings", h_2, h_1, ff},
        {"Complete*Elementary->Forgotten:exponent-splittings", h_2, e_1, ff},
        {"Complete*PowerSum->Forgotten:exponent-splittings", h_2, p_1, ff},
        {"Elementary*Elementary->Forgotten:exponent-splittings", e_2, e_1, ff},
        {"Elementary*PowerSum->Forgotten:exponent-splittings", e_2, p_1, ff},
        {"PowerSum*PowerSum->Forgotten:exponent-splittings", p_2, p_1, ff}
        };
    scan(kernelCases, case -> (
        report := multiplyToBasisBench(
            case#1,
            case#2,
            case#3,
            "Kernels" => {case#0, "PowerSumReference"},
            "Repetitions" => 1,
            "Warmups" => 0,
            "Verify" => true);
        assert report#"Verified";
        reversedReport := multiplyToBasisBench(
            case#2,
            case#1,
            case#3,
            "Kernels" => {case#0, "PowerSumReference"},
            "Repetitions" => 1,
            "Warmups" => 0,
            "Verify" => true);
        assert reversedReport#"Verified"));

    -- Reversing a mixed pair uses the same stable kernel declaration and
    -- commutative picker endpoint.
    commuted := multiplyToBasisBench(
        h_2, S_3, S,
        "Kernels" => {
            "Automatic",
            "Schur*Complete->Schur:horizontal-Pieri",
            "PowerSumReference"
            },
        "Repetitions" => 1,
        "Warmups" => 0,
        "Verify" => true);
    assert commuted#"Verified";

    E = frac(QQ[t]);
    R = symmetricRing E;
    hallCases = {
        {"HallLittlewoodQ*HallLittlewoodQ->HallLittlewoodQ:generator-triangular-formula", Q_{3,1}, Q_2, Q},
        {"HallLittlewoodB*HallLittlewoodB->HallLittlewoodB:generator-triangular-formula", B_{3,1}, B_2, B},
        {"HallLittlewoodP*HallLittlewoodP->HallLittlewoodP:generator-triangular-formula", P_{3,1}, P_2, P},
        {"HallLittlewoodPOmega*HallLittlewoodPOmega->HallLittlewoodPOmega:generator-triangular-formula", Pomega_{3,1}, Pomega_2, Pomega}
        };
    scan(hallCases, case -> (
        report := multiplyToBasisBench(
            case#1,
            case#2,
            case#3,
            "Kernels" => {case#0, "PowerSumReference"},
            "Repetitions" => 1,
            "Warmups" => 0,
            "Verify" => true);
        assert report#"Verified";
        reversedReport := multiplyToBasisBench(
            case#2,
            case#1,
            case#3,
            "Kernels" => {case#0, "PowerSumReference"},
            "Repetitions" => 1,
            "Warmups" => 0,
            "Verify" => true);
        assert reversedReport#"Verified"));
    exit 0'

# The two strict structural workflows use explicit requests rather than
# process-global route forcing.
run_m2 'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; report=multiplyToBasisBench(h_3,h_2,h,"Kernels"=>{"Automatic","PowerSumReference"},"Repetitions"=>1,"Warmups"=>0); assert(report#"Verified"); exit 0'
run_m2 'debug needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; report=multiplyToBasisBench(Q_2,Q_1,Q,"Kernels"=>{"Automatic","PowerSumReference"},"Repetitions"=>1,"Warmups"=>0); assert(report#"Verified"); exit 0'
run_m2 'debug needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; scan({h,e,p,q,b}, B -> (report := multiplyToBasisBench(S_2,e_1,B,"Kernels"=>{"Automatic","PowerSumReference"},"Repetitions"=>1,"Warmups"=>0); assert(report#"Verified"))); exit 0'

# The outer benchmark boundaries report the actual complete-input or
# complete-term strategy and whether that strategy invokes strict binary
# multiplication.  Trace executions are untimed.
run_m2 'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; BS=basis(R,"Schur"); BP=basis(R,"PowerSum"); A=multiplyExpressionsToBasisBenchExecute(S_2+S_{1,1},S_1,BS,"Automatic",false); K=multiplyExpressionsToBasisBenchExecute(S_2+S_{1,1},S_1,BS,"KernelDistribution",false); F=multiplyExpressionsToBasisBenchExecute(S_2+S_{1,1},S_1,BS,"PowerSumFallback",false); M=multiplyExpressionsToBasisBenchExecute(S_2+S_{1,1},S_1,BP,"MultiplicativeTarget",false); assert(A==K and K==F and F==M); exit 0'
assert_bench_trace_count \
    1 \
    'multiplication-expression: strategy=direct-kernel-distribution strict-binary=yes target=Schur' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; B=basis(R,"Schur"); A=multiplyExpressionsToBasisBenchExecute(S_{3,1}+S_{2,2},h_2+h_{1,1},B,true); exit 0'
assert_bench_trace_count \
    1 \
    'multiplication-expression: strategy=complete-input-power-sums strict-binary=no target=HallLittlewoodQ' \
    'debug needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; BQ=basis(R,"HallLittlewoodQ"); A=multiplyExpressionsToBasisBenchExecute(Q_2+Q_{1,1},B_2+B_{1,1},BQ,true); exit 0'
assert_bench_trace_count \
    1 \
    'multiplication-term: strategy=multiplicative-target-balanced strict-binary=no factors=4 target=PowerSum' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasisBench(S_2*h_2*e_1*S_1,p,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'multiplication-term: strategy=target-closed-fold strict-binary=yes factors=4 family=Schur-compatible target=Schur' \
    'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; A=toBasisBench(h_{2,1}*e_2*p_{2,1}*S_1,S,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'multiplication-term: strategy=target-closed-fold strict-binary=yes factors=3 family=HallLittlewoodQ-same-basis target=HallLittlewoodQ' \
    'debug needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; A=toBasisBench(Q_{2,1}*Q_2*Q_1,Q,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'
assert_bench_trace_count \
    1 \
    'multiplication-term: strategy=complete-term-power-sums strict-binary=no factors=3 target=HallLittlewoodQ' \
    'debug needsPackage "SymmetricRings"; E=frac(QQ[t]); R=symmetricRing E; A=toBasisBench(Q_{2,1}*e_2*B_{1,1},Q,"Repetitions"=>1,"Warmups"=>0,"Track"=>true); exit 0'

run_m2 'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(toBasisBench(p_2,S,"Plans"=>"unknown-plan","Repetitions"=>1,"Warmups"=>0); false) else true); exit 0'
run_m2 'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(toBasisBench(e_2,m,"Plans"=>"Elementary->Monomial:via-power-sums","Repetitions"=>1,"Warmups"=>0); false) else true); exit 0'
run_m2 'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(multiplyToBasisBench(S_3,S_2,S,"Kernels"=>"Schur*Complete->Schur:horizontal-Pieri","Repetitions"=>1,"Warmups"=>0); false) else true); exit 0'
run_m2 'debug needsPackage "SymmetricRings"; R=symmetricRing QQ; assert(try(multiplyToBasisBench(2*S_2,h_1,S,"Repetitions"=>1,"Warmups"=>0); false) else true); assert(try(userSymmetricElement(R,rawSymmetricRingsMultiplyToBasisBench(raw S_2,raw h_1,S#"BasisId","Schur*Complete->Schur:horizontal-Pieri",true,false)); false) else true); exit 0'

# Ordinary conversion is invariant under every former conversion diagnostic,
# including unrelated former ambient diagnostics.
ordinary_output=$(
    M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN='unknown-plan' \
    M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1 \
    M2_SYMMETRIC_RINGS_CHECK_ALL_CONVERSION_PLANS=1 \
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

if grep -R -E \
    'MultiplicationPlan|M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings" \
    --include='*.cpp' --include='*.hpp' >/dev/null ||
   test -e \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-products.cpp" ||
   test -e \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/basis-conversion-products.hpp" ||
   grep -E '"(Schur|Complete|Elementary|PowerSum|Monomial|Forgotten)\*' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/binary-multiplication.cpp" \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/multiplication.cpp" \
    >/dev/null ||
   grep -E \
    'selectBinaryMultiplicationKernel|multiplyCanonicalBasisTerms[[:space:]]*[(]|multiplyTermToBasis[[:space:]]*[(]|multiplyToBasis[[:space:]]*[(]|convertCanonicalExpressionToBasis[[:space:]]*[(]' \
    "$SOURCE_ROOT/Macaulay2/e/symmetric-rings/multiplication-kernels.cpp" \
    >/dev/null
then
    printf '%s\n' \
        'binary-multiplication contributor-surface invariant failed' >&2
    exit 1
fi

printf '%s\n' \
    'Conversion-plan and binary-kernel selection tests passed'
