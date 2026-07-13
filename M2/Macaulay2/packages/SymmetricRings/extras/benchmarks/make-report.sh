#!/bin/sh
set -eu

SYSTEM_FILE=$1
COMPARISON_FILE=$2
CONDITIONS_FILE=${3:-}
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

printf '## Overall summary\n\n'
awk -F '\t' -v run_health="$run_health" '
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
    printf "| Run health | Families | Cases | Improvements | Regressions | Stable | New |\n"
    printf "|---|---:|---:|---:|---:|---:|---:|\n"
    printf "| **%s** | %d | %d | %d | %d | %d | %d |\n\n", \
        run_health, families + 0, cases + 0, \
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
