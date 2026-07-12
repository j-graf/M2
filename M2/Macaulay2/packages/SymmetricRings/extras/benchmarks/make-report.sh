#!/bin/sh
set -eu

SYSTEM_FILE=$1
COMPARISON_FILE=$2

printf '# SymmetricRings benchmark report\n\n'
printf '## System\n\n'
printf '| Field | Value |\n|---|---|\n'
awk -F '\t' 'NR > 1 {
    gsub(/\|/, "\\|", $1); gsub(/\|/, "\\|", $2)
    printf "| %s | %s |\n", $1, $2
}' "$SYSTEM_FILE"

printf '\n## Baseline comparison\n\n'
awk -F '\t' '
NR == 1 {
    printf "|"
    for (i=1; i<=NF; i++) printf " %s |", $i
    printf "\n|"
    for (i=1; i<=NF; i++) printf "---|"
    printf "\n"
    next
}
{
    printf "|"
    for (i=1; i<=NF; i++) {
        gsub(/\|/, "\\|", $i)
        printf " %s |", $i
    }
    printf "\n"
}' "$COMPARISON_FILE"

