#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../../../../.." && pwd)
M2_BIN=${M2_BIN:-$SOURCE_ROOT/BUILD/build/M2}
DEVICE_PROFILE=${SYMRINGS_BENCH_DEVICE:-}
REPETITIONS=3
FAMILY=
TIER=
CASE_ID=
RUN_NAME=
LIST_ONLY=0
VERIFY=1
RECORDS=
LATEST=
RECORDS_OVERRIDE=
LATEST_OVERRIDE=
VARIED_MODE=
VARIED_LEVEL=
RANDOM_SEED=
NEW_ONLY=0
NEW_REFERENCE_HISTORY=none
ESTIMATE_ONLY=0

usage() {
    printf '%s\n' \
      "Usage: $0 [--device PROFILE] [--family NAME] [--tier NAME] [--case ID]" \
      "          [--repetitions N]" \
      "          [--varied-fixed LEVEL | --varied-random LEVEL] [--seed N]" \
      "          [--new] [--output LABEL] [--records FILE] [--latest FILE]" \
      "          [--list | --estimate] [--no-verify]" \
      "" \
      "Varied levels: light, standard, thorough." \
      "Device profiles are user-chosen slugs such as m4-laptop or workstation." \
      "A profile is required for runs, estimates, and --new selections." \
      "Environment: SYMRINGS_BENCH_DEVICE and M2_BIN may supply defaults."
}

set_varied() {
    requested_mode=$1
    requested_level=$2
    if [ -n "$VARIED_MODE" ]; then
        printf 'Choose only one varied selection mode.\n' >&2
        exit 2
    fi
    case "$requested_level" in
        light|standard|thorough) ;;
        *) printf 'Unknown varied level: %s\n' "$requested_level" >&2; exit 2 ;;
    esac
    VARIED_MODE=$requested_mode
    VARIED_LEVEL=$requested_level
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --family) FAMILY=$2; shift 2 ;;
        --tier) TIER=$2; shift 2 ;;
        --case) CASE_ID=$2; shift 2 ;;
        --repetitions) REPETITIONS=$2; shift 2 ;;
        --varied-fixed) set_varied fixed "$2"; shift 2 ;;
        --varied-random) set_varied random "$2"; shift 2 ;;
        --seed) RANDOM_SEED=$2; shift 2 ;;
        --new) NEW_ONLY=1; shift ;;
        --estimate) ESTIMATE_ONLY=1; shift ;;
        --device) DEVICE_PROFILE=$2; shift 2 ;;
        --output) RUN_NAME=$2; shift 2 ;;
        --records) RECORDS_OVERRIDE=$2; shift 2 ;;
        --latest) LATEST_OVERRIDE=$2; shift 2 ;;
        --list) LIST_ONLY=1; shift ;;
        --no-verify) VERIFY=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ "$LIST_ONLY" -eq 1 ] && [ "$ESTIMATE_ONLY" -eq 1 ]; then
    printf 'Choose only one of --list and --estimate.\n' >&2
    exit 2
fi

if [ "$LIST_ONLY" -eq 0 ] || [ "$NEW_ONLY" -eq 1 ]; then
    if [ -z "$DEVICE_PROFILE" ]; then
        printf '%s\n' \
          'A device profile is required. Use --device PROFILE or set SYMRINGS_BENCH_DEVICE.' >&2
        exit 2
    fi
fi
if [ -n "$DEVICE_PROFILE" ]; then
    case "$DEVICE_PROFILE" in
        [A-Za-z0-9]*)
            case "$DEVICE_PROFILE" in
                *[!A-Za-z0-9._-]*)
                    printf 'Invalid device profile: %s\n' "$DEVICE_PROFILE" >&2
                    exit 2
                    ;;
            esac
            ;;
        *)
            printf 'Invalid device profile: %s\n' "$DEVICE_PROFILE" >&2
            exit 2
            ;;
    esac
    HISTORY_DIR="$SCRIPT_DIR/history/$DEVICE_PROFILE"
    RECORDS=${RECORDS_OVERRIDE:-$HISTORY_DIR/records.tsv}
    LATEST=${LATEST_OVERRIDE:-$HISTORY_DIR/latest.tsv}
    if [ "$RECORDS" = "$LATEST" ]; then
        printf '%s\n' 'Records and latest history must use different files.' >&2
        exit 2
    fi
    if [ -n "$RECORDS_OVERRIDE" ]; then
        RECORDS_HISTORY_LABEL="custom:$(basename -- "$RECORDS")"
    else
        RECORDS_HISTORY_LABEL="history/$DEVICE_PROFILE/records.tsv"
    fi
    if [ -n "$LATEST_OVERRIDE" ]; then
        LATEST_HISTORY_LABEL="custom:$(basename -- "$LATEST")"
    else
        LATEST_HISTORY_LABEL="history/$DEVICE_PROFILE/latest.tsv"
    fi
fi

validate_history_device() {
    history_file=$1
    [ -f "$history_file" ] || return 0
    if ! awk -F '\t' -v expected="$DEVICE_PROFILE" '
        NR == 1 {
            for (i = 1; i <= NF; i++)
                if ($i == "device_profile") deviceColumn = i
            next
        }
        deviceColumn > 0 && $deviceColumn != "" && $deviceColumn != expected {
            exit 1
        }
    ' "$history_file"; then
        printf 'History file belongs to another device profile: %s\n' \
            "$history_file" >&2
        exit 2
    fi
}

if [ -n "$DEVICE_PROFILE" ]; then
    validate_history_device "$RECORDS"
    validate_history_device "$LATEST"
fi
case "$REPETITIONS" in
    *[!0-9]*|'') printf 'Repetitions must be a positive integer: %s\n' "$REPETITIONS" >&2; exit 2 ;;
esac
if [ "$REPETITIONS" -le 0 ]; then
    printf 'Repetitions must be a positive integer: %s\n' "$REPETITIONS" >&2
    exit 2
fi

if [ -n "$RANDOM_SEED" ] && [ "$VARIED_MODE" != random ]; then
    printf '%s\n' '--seed requires --varied-random.' >&2
    exit 2
fi
if [ "$VARIED_MODE" = random ]; then
    if [ -z "$RANDOM_SEED" ]; then RANDOM_SEED=$(date '+%s'); fi
    case "$RANDOM_SEED" in
        *[!0-9]*|'') printf 'Random seed must be a nonnegative integer: %s\n' "$RANDOM_SEED" >&2; exit 2 ;;
    esac
fi

if [ ! -x "$M2_BIN" ]; then
    printf 'Macaulay2 executable not found: %s\n' "$M2_BIN" >&2
    exit 2
fi

TMP_HOME=$(mktemp -d "${TMPDIR:-/tmp}/symmetricrings-bench.XXXXXX")
trap 'rm -rf "$TMP_HOME"' EXIT HUP INT TERM

plan_cases() {
    fixed_level=
    if [ "$VARIED_MODE" = fixed ]; then fixed_level=$VARIED_LEVEL; fi
    env HOME="$TMP_HOME" \
        SYMRINGS_BENCH_MODE=plan \
        SYMRINGS_BENCH_FAMILY="$FAMILY" \
        SYMRINGS_BENCH_TIER="$TIER" \
        SYMRINGS_BENCH_CASE="$CASE_ID" \
        SYMRINGS_BENCH_VARIED_FIXED="$fixed_level" \
        "$M2_BIN" --script "$SCRIPT_DIR/runner.m2"
}

randomize_plan() {
    case "$VARIED_LEVEL" in
        light) picks_per_family=1 ;;
        standard) picks_per_family=2 ;;
        thorough) picks_per_family=4 ;;
    esac
    awk -F '|' -v picks="$picks_per_family" -v seed="$RANDOM_SEED" '
        BEGIN { srand(seed) }
        {
            family = $2
            if (!family_seen[family]++) families[++family_count] = family
            count = ++seen[family]
            if (count <= picks) selected[family, count] = $0
            else {
                slot = int(rand() * count) + 1
                if (slot <= picks) selected[family, slot] = $0
            }
        }
        END {
            for (i = 1; i <= family_count; i++) {
                family = families[i]
                limit = seen[family] < picks ? seen[family] : picks
                for (slot = 1; slot <= limit; slot++)
                    print selected[family, slot]
            }
        }'
}

case_plan=$(plan_cases)
if [ "$NEW_ONLY" -eq 1 ] && [ -n "$case_plan" ]; then
    known_cases_file="$TMP_HOME/known-cases.txt"
    plan_file="$TMP_HOME/case-plan.txt"
    printf '%s\n' __known_case_sentinel__ > "$known_cases_file"
    reference_histories=
    if [ -f "$LATEST" ]; then
        awk -F '\t' 'NR > 1 && $1 != "" {print $1}' \
            "$LATEST" >> "$known_cases_file"
        reference_histories=$LATEST_HISTORY_LABEL
    fi
    if [ -f "$RECORDS" ]; then
        awk -F '\t' 'NR > 1 && $1 != "" {print $1}' \
            "$RECORDS" >> "$known_cases_file"
        if [ -n "$reference_histories" ]; then
            reference_histories="$reference_histories + $RECORDS_HISTORY_LABEL"
        else
            reference_histories=$RECORDS_HISTORY_LABEL
        fi
    fi
    if [ -n "$reference_histories" ]; then
        NEW_REFERENCE_HISTORY=$reference_histories
    fi
    printf '%s\n' "$case_plan" > "$plan_file"
    case_plan=$(awk -F '|' '
        NR == FNR {known[$1] = 1; next}
        !($1 in known) {print}
    ' "$known_cases_file" "$plan_file")
fi
if [ "$VARIED_MODE" = random ] && [ -n "$case_plan" ]; then
    case_plan=$(printf '%s\n' "$case_plan" | randomize_plan)
    printf 'Random varied selection seed: %s\n' "$RANDOM_SEED" >&2
fi

if [ "$LIST_ONLY" -eq 1 ]; then
    if [ -n "$case_plan" ]; then
        printf '%s\n' "$case_plan" | cut -d '|' -f 1
    fi
    exit 0
fi

if [ "$ESTIMATE_ONLY" -eq 1 ]; then
    if [ -z "$case_plan" ]; then
        printf 'No benchmark cases matched the filters.\n'
        exit 0
    fi
    estimate_plan="$TMP_HOME/estimate-plan.txt"
    printf '%s\n' "$case_plan" > "$estimate_plan"
    if [ -f "$LATEST" ]; then
        estimate_latest=$LATEST
        latest_label=$LATEST_HISTORY_LABEL
    else
        estimate_latest=/dev/null
        latest_label=none
    fi
    PROCESS_OVERHEAD_SECONDS=${SYMRINGS_BENCH_PROCESS_OVERHEAD_SECONDS:-0.30}
    CALIBRATION_SECONDS=${SYMRINGS_BENCH_CALIBRATION_SECONDS:-0.19}
    if [ -f "$RECORDS" ]; then
        estimate_records=$RECORDS
    else
        estimate_records=/dev/null
    fi
    "$SCRIPT_DIR/estimate-run.awk" \
        -v recordFile="$estimate_records" \
        -v latestFile="$estimate_latest" \
        -v planFile="$estimate_plan" \
        -v repetitions="$REPETITIONS" \
        -v processOverhead="$PROCESS_OVERHEAD_SECONDS" \
        -v calibrationSeconds="$CALIBRATION_SECONDS" \
        -v latestLabel="$latest_label" \
        -v deviceProfile="$DEVICE_PROFILE" \
        "$estimate_records" "$estimate_latest" "$estimate_plan"
    exit 0
fi

if [ -z "$case_plan" ]; then
    printf 'No benchmark cases matched the filters.\n' >&2
    exit 2
fi

case "$RUN_NAME" in
    .|..|*/*) printf 'Benchmark run name must be one path component: %s\n' "$RUN_NAME" >&2; exit 2 ;;
esac

RUN_TIMESTAMP=$(date '+%Y%m%d-%H%M%S')
if [ -z "$RUN_NAME" ]; then
    RUN_NAME=$RUN_TIMESTAMP
else
    RUN_NAME="$RUN_TIMESTAMP-$RUN_NAME"
fi

RUN_DIR="$SCRIPT_DIR/results/$DEVICE_PROFILE/$RUN_NAME"

if [ -e "$RUN_DIR" ]; then
    printf 'Benchmark result directory already exists: %s\n' "$RUN_DIR" >&2
    exit 2
fi
mkdir -p "$RUN_DIR"

raw_file="$RUN_DIR/raw.tsv"
summary_file="$RUN_DIR/summary.tsv"
comparison_file="$RUN_DIR/comparison.tsv"
system_file="$RUN_DIR/system.tsv"
report_file="$RUN_DIR/report.md"
conditions_file="$RUN_DIR/conditions.tsv"

"$SCRIPT_DIR/collect-system-info.sh" "$M2_BIN" "$SOURCE_ROOT" \
    "$DEVICE_PROFILE" > "$system_file"

printf 'case_id\tfamily\toperation\ttier\tcoefficient_ring\trepetition\tcpu_seconds\twall_seconds\tresult_terms\tresult_weight\tinput_group\tlambda\tmu\tprobe\tlambda_weight\tlambda_length\tmu_weight\tmu_length\texpected_weight\tpair_class\n' > "$raw_file"
printf 'kind\tscope\tphase\tmetric\tvalue\n' > "$conditions_file"
printf 'configuration\trun\tselection\tvaried_mode\t%s\n' \
    "${VARIED_MODE:-none}" >> "$conditions_file"
printf 'configuration\trun\tselection\tvaried_level\t%s\n' \
    "${VARIED_LEVEL:-none}" >> "$conditions_file"
printf 'configuration\trun\tselection\trandom_seed\t%s\n' \
    "${RANDOM_SEED:-none}" >> "$conditions_file"
if [ "$NEW_ONLY" -eq 1 ]; then new_only_value=true; else new_only_value=false; fi
printf 'configuration\trun\tselection\tnew_only\t%s\n' \
    "$new_only_value" >> "$conditions_file"
printf 'configuration\trun\tselection\tnew_reference_history\t%s\n' \
    "$NEW_REFERENCE_HISTORY" >> "$conditions_file"
printf 'configuration\trun\tidentity\trun_name\t%s\n' \
    "$RUN_NAME" >> "$conditions_file"
printf 'configuration\trun\tidentity\tdevice_profile\t%s\n' \
    "$DEVICE_PROFILE" >> "$conditions_file"
printf 'configuration\trun\tidentity\trecords_history\t%s\n' \
    "$RECORDS_HISTORY_LABEL" >> "$conditions_file"
printf 'configuration\trun\tidentity\tlatest_history\t%s\n' \
    "$LATEST_HISTORY_LABEL" >> "$conditions_file"
"$SCRIPT_DIR/collect-run-conditions.sh" before >> "$conditions_file"

run_calibration() {
    calibration_family=$1
    calibration_phase=$2
    printf '[calibration %s] %s\n' "$calibration_phase" "$calibration_family" >&2
    if ! calibration_output=$(env HOME="$TMP_HOME" \
            SYMRINGS_BENCH_MODE=calibration \
            "$M2_BIN" --script "$SCRIPT_DIR/runner.m2"); then
        printf 'Calibration probe failed: %s %s\n' \
            "$calibration_family" "$calibration_phase" >&2
        exit 1
    fi
    calibration_line=$(printf '%s\n' "$calibration_output" |
        awk -F '|' 'index($0, "CALIBRATION|") == 1 {print; exit}')
    if [ -z "$calibration_line" ]; then
        printf 'Calibration probe produced no result: %s %s\n' \
            "$calibration_family" "$calibration_phase" >&2
        exit 1
    fi
    calibration_case=$(printf '%s\n' "$calibration_line" | cut -d '|' -f 2)
    calibration_cpu=$(printf '%s\n' "$calibration_line" | cut -d '|' -f 3)
    calibration_wall=$(printf '%s\n' "$calibration_line" | cut -d '|' -f 4)
    printf 'calibration\t%s\t%s\tcase_id\t%s\n' \
        "$calibration_family" "$calibration_phase" "$calibration_case" >> "$conditions_file"
    printf 'calibration\t%s\t%s\tcpu_seconds\t%s\n' \
        "$calibration_family" "$calibration_phase" "$calibration_cpu" >> "$conditions_file"
    printf 'calibration\t%s\t%s\twall_seconds\t%s\n' \
        "$calibration_family" "$calibration_phase" "$calibration_wall" >> "$conditions_file"
}

family_list=$(printf '%s\n' "$case_plan" |
    awk -F '|' '!seen[$2]++ {print $2}')

for selected_family in $family_list; do
    run_calibration "$selected_family" before
    selected_cases=$(printf '%s\n' "$case_plan" |
        awk -F '|' -v family="$selected_family" '$2 == family {print $1}')
    for selected_case in $selected_cases; do
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
            printf '%s\n' "$benchmark_line" >> "$raw_file"
            repetition=$((repetition + 1))
        done
    done
    run_calibration "$selected_family" after
done

"$SCRIPT_DIR/collect-run-conditions.sh" after >> "$conditions_file"

printf 'Raw results written to %s\n' "$raw_file" >&2

"$SCRIPT_DIR/summarize-results.awk" "$raw_file" > "$summary_file"
if [ -f "$RECORDS" ]; then
    records_for_comparison=$RECORDS
else
    records_for_comparison="$TMP_HOME/empty-records.tsv"
    printf 'case_id\trecorded_at\tmedian_cpu_seconds\tcoefficient_ring\truns\trun_name\tdevice_profile\n' \
        > "$records_for_comparison"
fi
"$SCRIPT_DIR/compare-results.awk" \
    "$records_for_comparison" "$raw_file" > "$comparison_file"
"$SCRIPT_DIR/make-report.sh" "$system_file" "$comparison_file" \
    "$conditions_file" > "$report_file"

records_update="$TMP_HOME/records.tsv"
if [ "$VERIFY" -eq 1 ]; then
    mkdir -p "$(dirname -- "$RECORDS")" "$(dirname -- "$LATEST")"
    recorded_at=$(awk -F '\t' '$1 == "captured_at" {print $2; exit}' "$system_file")
    RECORD_RUN_NAME="$RUN_NAME" RECORD_RECORDED_AT="$recorded_at" \
        RECORD_DEVICE_PROFILE="$DEVICE_PROFILE" \
        "$SCRIPT_DIR/update-records.awk" \
        "$records_for_comparison" "$summary_file" > "$records_update"
    latest_for_update=$LATEST
    if [ ! -f "$latest_for_update" ]; then
        latest_for_update="$TMP_HOME/empty-latest.tsv"
        printf 'case_id\trecorded_at\tmedian_cpu_seconds\tcoefficient_ring\truns\trun_name\tdevice_profile\n' \
            > "$latest_for_update"
    fi
    latest_update="$TMP_HOME/latest.tsv"
    RECORD_RUN_NAME="$RUN_NAME" RECORD_RECORDED_AT="$recorded_at" \
        RECORD_DEVICE_PROFILE="$DEVICE_PROFILE" \
        "$SCRIPT_DIR/update-latest.awk" \
        "$latest_for_update" "$summary_file" > "$latest_update"
    mv "$records_update" "$RECORDS"
    mv "$latest_update" "$LATEST"
fi
printf 'Summary written to %s\n' "$summary_file" >&2
printf 'Performance comparison written to %s\n' "$comparison_file" >&2
printf 'System information written to %s\n' "$system_file" >&2
printf 'Run conditions written to %s\n' "$conditions_file" >&2
printf 'Markdown report written to %s\n' "$report_file" >&2
if [ "$VERIFY" -eq 1 ]; then
    printf 'Fastest qualifying records updated in %s\n' "$RECORDS" >&2
    printf 'Most recent per-case results updated in %s\n' "$LATEST" >&2
else
    printf '%s\n' \
      'Diagnostic --no-verify run: records and latest history were not updated.' >&2
fi
printf '\nPerformance comparison:\n' >&2
cat "$comparison_file" >&2
