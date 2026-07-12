#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../../../../.." && pwd)
M2_BIN=${M2_BIN:-$SOURCE_ROOT/BUILD/build/M2}
REPETITIONS=3
FAMILY=
TIER=
CASE_ID=
OUTPUT=
LIST_ONLY=0
VERIFY=1
BASELINES="$SCRIPT_DIR/baselines.tsv"

usage() {
    printf '%s\n' \
      "Usage: $0 [--family NAME] [--tier NAME] [--case ID] [--repetitions N]" \
      "          [--output FILE] [--baselines FILE] [--list] [--no-verify]" \
      "" \
      "Environment: M2_BIN may override the local Macaulay2 executable."
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --family) FAMILY=$2; shift 2 ;;
        --tier) TIER=$2; shift 2 ;;
        --case) CASE_ID=$2; shift 2 ;;
        --repetitions) REPETITIONS=$2; shift 2 ;;
        --output) OUTPUT=$2; shift 2 ;;
        --baselines) BASELINES=$2; shift 2 ;;
        --list) LIST_ONLY=1; shift ;;
        --no-verify) VERIFY=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ ! -x "$M2_BIN" ]; then
    printf 'Macaulay2 executable not found: %s\n' "$M2_BIN" >&2
    exit 2
fi

TMP_HOME=$(mktemp -d "${TMPDIR:-/tmp}/symmetricrings-bench.XXXXXX")
trap 'rm -rf "$TMP_HOME"' EXIT HUP INT TERM

list_cases() {
    env HOME="$TMP_HOME" \
        SYMRINGS_BENCH_MODE=list \
        SYMRINGS_BENCH_FAMILY="$FAMILY" \
        SYMRINGS_BENCH_TIER="$TIER" \
        SYMRINGS_BENCH_CASE="$CASE_ID" \
        "$M2_BIN" --script "$SCRIPT_DIR/runner.m2"
}

if [ "$LIST_ONLY" -eq 1 ]; then
    list_cases
    exit 0
fi

if [ -z "$OUTPUT" ]; then
    stamp=$(date '+%Y%m%d-%H%M%S')
    OUTPUT="$SCRIPT_DIR/results-$stamp.tsv"
fi

system_file="${OUTPUT%.tsv}-system.tsv"
"$SCRIPT_DIR/collect-system-info.sh" "$M2_BIN" "$SOURCE_ROOT" > "$system_file"

if [ ! -e "$OUTPUT" ]; then
    printf 'case_id\tfamily\toperation\ttier\tcoefficient_ring\trepetition\tcpu_seconds\twall_seconds\tresult_terms\tresult_weight\tinput_group\tlambda\tmu\tprobe\tlambda_weight\tlambda_length\tmu_weight\tmu_length\texpected_weight\tpair_class\n' > "$OUTPUT"
fi

case_list=$(list_cases)
if [ -z "$case_list" ]; then
    printf 'No benchmark cases matched the filters.\n' >&2
    exit 2
fi

printf '%s\n' "$case_list" | while IFS= read -r selected_case; do
    [ -n "$selected_case" ] || continue
    repetition=1
    while [ "$repetition" -le "$REPETITIONS" ]; do
        printf '[%s/%s] %s\n' "$repetition" "$REPETITIONS" "$selected_case" >&2
        if ! worker_output=$(env HOME="$TMP_HOME" \
                SYMRINGS_BENCH_MODE=run \
                SYMRINGS_BENCH_CASE="$selected_case" \
                SYMRINGS_BENCH_REPETITION="$repetition" \
                SYMRINGS_BENCH_VERIFY="$VERIFY" \
                "$M2_BIN" --script "$SCRIPT_DIR/runner.m2"); then
            printf 'Benchmark worker failed: %s\n' "$selected_case" >&2
            exit 1
        fi
        benchmark_line=$(printf '%s\n' "$worker_output" |
            awk -F '|' 'BEGIN {OFS="\t"} index($0, "BENCH|") == 1 {
                $1=""; sub(/^\t/, ""); print
            }')
        if [ -z "$benchmark_line" ]; then
            printf 'Benchmark worker produced no result: %s\n' "$selected_case" >&2
            [ -z "$worker_output" ] || printf '%s\n' "$worker_output" >&2
            exit 1
        fi
        printf '%s\n' "$benchmark_line" >> "$OUTPUT"
        repetition=$((repetition + 1))
    done
done

printf 'Results appended to %s\n' "$OUTPUT" >&2

summary_file="${OUTPUT%.tsv}-summary.tsv"
comparison_file="${OUTPUT%.tsv}-comparison.tsv"
report_file="${OUTPUT%.tsv}-report.md"
"$SCRIPT_DIR/summarize-results.awk" "$OUTPUT" > "$summary_file"
"$SCRIPT_DIR/compare-results.awk" "$BASELINES" "$OUTPUT" > "$comparison_file"
"$SCRIPT_DIR/make-report.sh" "$system_file" "$comparison_file" > "$report_file"
printf 'Summary written to %s\n' "$summary_file" >&2
printf 'Baseline comparison written to %s\n' "$comparison_file" >&2
printf 'System information written to %s\n' "$system_file" >&2
printf 'Markdown report written to %s\n' "$report_file" >&2
printf '\nBaseline comparison:\n' >&2
cat "$comparison_file" >&2
