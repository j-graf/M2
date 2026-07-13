#!/bin/sh
set -eu

SYSTEM_FILE=$1
COMPARISON_FILE=$2
CONDITIONS_FILE=${3:-}
BASELINES_FILE=${4:-}
REPORT_DIR=$(dirname -- "$0")

printf '# SymmetricRings benchmark report\n\n'

if [ -n "$CONDITIONS_FILE" ] && [ -s "$CONDITIONS_FILE" ]; then
    run_health=$("$REPORT_DIR/summarize-run-conditions.awk" "$CONDITIONS_FILE" |
        awk '/^Overall health:/ {
            sub(/^Overall health: \*\*/, "")
            sub(/\*\*\..*$/, "")
            print
            exit
        }')
    [ -n "$run_health" ] || run_health=unknown
else
    run_health='not assessed'
fi

if [ -n "$BASELINES_FILE" ] && [ -s "$BASELINES_FILE" ]; then
    baseline_system=$(
        "$REPORT_DIR/compare-baseline-system.awk" \
            "$SYSTEM_FILE" "$BASELINES_FILE" "$COMPARISON_FILE")
    baseline_system_status=$(printf '%s\n' "$baseline_system" |
        awk -F '\t' '$1 == "overall" {print $4; exit}')
    [ -n "$baseline_system_status" ] || baseline_system_status=unknown
else
    baseline_system=''
    baseline_system_status=unknown
fi

printf '## Overall summary\n\n'
awk -F '\t' -v run_health="$run_health" \
    -v baseline_system_status="$baseline_system_status" '
function comparisonDetail(case_id, versus_best, versus_recent) {
    return "`" case_id "` (" versus_best "% vs best; " \
        versus_recent "% vs recent)"
}
function appendDetail(current, detail) {
    return current == "" ? detail : current ", " detail
}
NR == 1 { next }
{
    cases++
    if (!seen_family[$2]++) families++
    classification = $NF
    counts[classification]++
    if (classification == "improvement")
        improvements = appendDetail(improvements,
            comparisonDetail($1, $9, $11))
    else if (classification == "regression")
        regressions = appendDetail(regressions,
            comparisonDetail($1, $9, $11))
}
END {
    printf "| Run health | Baseline system | Families | Cases | Improvements | Regressions | Stable | New |\n"
    printf "|---|---|---:|---:|---:|---:|---:|---:|\n"
    printf "| **%s** | **%s** | %d | %d | %d | %d | %d | %d |\n\n", \
        run_health, baseline_system_status, families + 0, cases + 0, \
        counts["improvement"] + 0, counts["regression"] + 0, \
        counts["stable"] + 0, counts["new"] + 0
    if (improvements != "")
        printf "Improvements: %s.\n\n", improvements
    if (regressions != "")
        printf "Regressions: %s.\n\n", regressions
}' "$COMPARISON_FILE"

printf '## Baseline comparison\n\n'
families=$(awk -F '\t' 'NR > 1 { print $2 }' "$COMPARISON_FILE" | LC_ALL=C sort -u)

if [ -z "$families" ]; then
    printf 'No benchmark cases were recorded.\n'
else
  printf '%s\n' "$families" | while IFS= read -r family; do
    [ -n "$family" ] || continue
    awk -F '\t' -v selected_family="$family" '
    function escaped(value) {
        gsub(/\|/, "\\|", value)
        return value
    }
    function comparisonDetail(case_id, versus_best, versus_recent) {
        return "`" case_id "` (" versus_best "% vs best; " \
            versus_recent "% vs recent)"
    }
    function appendDetail(current, detail) {
        return current == "" ? detail : current ", " detail
    }
    NR == 1 {
        field_count = NF
        for (i = 1; i <= NF; i++) header[i] = $i
        next
    }
    $2 == selected_family {
        row_count++
        for (i = 1; i <= NF; i++) rows[row_count, i] = $i
        classification = $NF
        counts[classification]++
        if (classification == "improvement")
            improvements = appendDetail(improvements,
                comparisonDetail($1, $9, $11))
        else if (classification == "regression")
            regressions = appendDetail(regressions,
                comparisonDetail($1, $9, $11))
    }
    END {
        printf "### %s\n\n", selected_family
        printf "| Cases | Improvements | Regressions | Stable | New |\n"
        printf "|---:|---:|---:|---:|---:|\n"
        printf "| %d | %d | %d | %d | %d |\n\n", \
            row_count, counts["improvement"] + 0, \
            counts["regression"] + 0, counts["stable"] + 0, \
            counts["new"] + 0
        if (improvements != "")
            printf "Improvements: %s.\n\n", improvements
        if (regressions != "")
            printf "Regressions: %s.\n\n", regressions

        printf "|"
        for (i = 1; i <= field_count; i++)
            if (i != 2) printf " %s |", escaped(header[i])
        printf "\n|"
        for (i = 1; i <= field_count; i++)
            if (i != 2) printf "---|"
        printf "\n"
        for (row = 1; row <= row_count; row++) {
            printf "|"
            for (i = 1; i <= field_count; i++)
                if (i != 2) printf " %s |", escaped(rows[row, i])
            printf "\n"
        }
        printf "\n"
    }' "$COMPARISON_FILE"
  done
fi

printf '## Baseline system comparison\n\n'
case "$baseline_system_status" in
    same)
        printf 'The major system configuration matches the most recent accepted baseline configuration for the selected cases.\n\n'
        ;;
    different)
        printf 'At least one major system field differs from the most recent accepted baseline configuration for the selected cases. Timing changes may therefore include machine effects.\n\n'
        ;;
    mixed)
        printf 'The selected cases have most recent accepted baselines from more than one major system configuration.\n\n'
        ;;
    *)
        printf 'The major system configuration could not be compared completely, usually because no selected baseline or structured baseline system metadata was available.\n\n'
        ;;
esac
printf '| Field | Current run | Accepted baseline | Comparison |\n'
printf '|---|---|---|---|\n'
if [ -n "$baseline_system" ]; then
    printf '%s\n' "$baseline_system" | awk -F '\t' 'NR > 1 && $1 != "overall" {
        for (i = 1; i <= 4; ++i) gsub(/\|/, "\\|", $i)
        printf "| %s | %s | %s | %s |\n", $1, $2, $3, $4
    }'
else
    printf '| CPU | unknown | unknown | unknown |\n'
    printf '| RAM | unknown | unknown | unknown |\n'
    printf '| Operating system | unknown | unknown | unknown |\n'
fi

printf '\n'

printf '## System\n\n'
printf '| Field | Value |\n|---|---|\n'
awk -F '\t' 'NR > 1 {
    gsub(/\|/, "\\|", $1); gsub(/\|/, "\\|", $2)
    printf "| %s | %s |\n", $1, $2
}' "$SYSTEM_FILE"

if [ -n "$CONDITIONS_FILE" ] && [ -s "$CONDITIONS_FILE" ]; then
    printf '\n'
    "$REPORT_DIR/summarize-run-conditions.awk" "$CONDITIONS_FILE"
else
    printf '\n## Run conditions\n\n'
    printf 'Overall health: **not assessed** (this run predates condition collection).\n\n'
    printf '%s\n' \
      '- **Clean:** conditions were stable and no acceptance concern was detected.' \
      '- **Warning:** moderate drift or resource pressure was detected; review before accepting baselines.' \
      '- **Compromised:** severe drift, thermal pressure, or a performance-limiting mode makes the run unsuitable as a baseline.'
fi
