#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TEST_DIR=$(mktemp -d "${TMPDIR:-/tmp}/symmetricrings-history-test.XXXXXX")
trap 'rm -rf "$TEST_DIR"' EXIT HUP INT TERM

records="$TEST_DIR/records.tsv"
latest="$TEST_DIR/latest.tsv"
summary="$TEST_DIR/summary.tsv"
updated_records="$TEST_DIR/updated-records.tsv"
updated_latest="$TEST_DIR/updated-latest.tsv"

printf '%s\n' \
  'case_id	recorded_at	median_cpu_seconds	coefficient_ring	runs	run_name	device_profile' \
  'case-a	2026-01-01 00:00:00 UTC	1	QQ	3	old-record	device-one' \
  > "$records"

printf '%s\n' \
  'case_id	recorded_at	median_cpu_seconds	coefficient_ring	runs	run_name	device_profile' \
  'case-a	2026-01-02 00:00:00 UTC	2	QQ	3	old-latest	device-one' \
  'case-b	2026-01-02 00:00:00 UTC	3	QQ	3	old-latest	device-one' \
  > "$latest"

printf '%s\n' \
  'case_id	family	operation	tier	coefficient_ring	runs	min_cpu_seconds	median_cpu_seconds	max_cpu_seconds' \
  'case-a	Family	Operation	Small	QQ	3	1.4	1.5	1.6' \
  'case-c	Family	Operation	Small	QQ	1	0.5	0.5	0.5' \
  > "$summary"

RECORD_RUN_NAME=new-run \
RECORD_RECORDED_AT='2026-01-03 00:00:00 UTC' \
RECORD_DEVICE_PROFILE=device-one \
    "$SCRIPT_DIR/update-records.awk" "$records" "$summary" > "$updated_records"

RECORD_RUN_NAME=new-run \
RECORD_RECORDED_AT='2026-01-03 00:00:00 UTC' \
RECORD_DEVICE_PROFILE=device-one \
    "$SCRIPT_DIR/update-latest.awk" "$latest" "$summary" > "$updated_latest"

awk -F '\t' '
    NR == 1 { next }
    $1 == "case-a" && $3 == 1 && $6 == "old-record" &&
        $7 == "device-one" { foundA = 1 }
    $1 == "case-c" { foundC = 1 }
    END { exit !(foundA && !foundC) }
' "$updated_records"

awk -F '\t' '
    NR == 1 { next }
    $1 == "case-a" && $3 == 1.5 && $5 == 3 && $6 == "new-run" &&
        $7 == "device-one" { foundA = 1 }
    $1 == "case-b" && $3 == 3 && $6 == "old-latest" { foundB = 1 }
    $1 == "case-c" && $3 == 0.5 && $5 == 1 && $6 == "new-run" &&
        $7 == "device-one" { foundC = 1 }
    END { exit !(foundA && foundB && foundC) }
' "$updated_latest"

printf '%s\n' 'history update tests passed'
