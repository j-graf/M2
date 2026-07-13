#!/bin/sh
set -eu

PHASE=${1:?usage: collect-run-conditions.sh before-or-after}
UNKNOWN=unknown

emit() {
    value=$(printf '%s' "$2" | tr '\t\r\n' '   ')
    printf 'snapshot\trun\t%s\t%s\t%s\n' "$PHASE" "$1" "$value"
}

kernel=$(uname -s 2>/dev/null || printf '%s' "$UNKNOWN")
captured_at=$(date '+%Y-%m-%d %H:%M:%S %Z')
snapshot_epoch=$(date '+%s')
power_source=$UNKNOWN
battery_state=$UNKNOWN
low_power_mode=$UNKNOWN
thermal_state=$UNKNOWN
memory_free_percent=$UNKNOWN
pageouts=$UNKNOWN
page_size_bytes=$UNKNOWN
load_averages=$UNKNOWN

if command -v uptime >/dev/null 2>&1; then
    load_averages=$(uptime 2>/dev/null |
        sed -E 's/^.*load averages?:[[:space:]]*//' || printf '%s' "$UNKNOWN")
fi

if [ "$kernel" = Darwin ]; then
    if command -v pmset >/dev/null 2>&1; then
        battery_output=$(pmset -g batt 2>/dev/null || true)
        power_source=$(printf '%s\n' "$battery_output" |
            awk -F "'" 'NR == 1 && NF >= 2 {print $2; exit}')
        battery_state=$(printf '%s\n' "$battery_output" |
            awk -F '\t' 'NR == 2 {print $2; exit}')
        low_power_mode=$(pmset -g custom 2>/dev/null |
            awk '/lowpowermode/ {print $2; exit}')
        case "$low_power_mode" in
            0) low_power_mode=disabled ;;
            1) low_power_mode=enabled ;;
        esac
    fi
    if command -v osascript >/dev/null 2>&1; then
        thermal_code=$(osascript -l JavaScript \
            -e 'ObjC.import("Foundation"); Number($.NSProcessInfo.processInfo.thermalState)' \
            2>/dev/null || true)
        case "$thermal_code" in
            0) thermal_state=nominal ;;
            1) thermal_state=fair ;;
            2) thermal_state=serious ;;
            3) thermal_state=critical ;;
        esac
    fi
    if command -v memory_pressure >/dev/null 2>&1; then
        memory_free_percent=$(memory_pressure 2>/dev/null |
            awk -F ': ' '/System-wide memory free percentage:/ {
                gsub(/%/, "", $2); print $2; exit
            }')
    fi
    if command -v vm_stat >/dev/null 2>&1; then
        vm_output=$(vm_stat 2>/dev/null || true)
        page_size_bytes=$(printf '%s\n' "$vm_output" |
            awk 'NR == 1 && match($0, /page size of [0-9]+ bytes/) {
                value = substr($0, RSTART, RLENGTH)
                gsub(/[^0-9]/, "", value); print value; exit
            }')
        pageouts=$(printf '%s\n' "$vm_output" |
            awk -F ': ' '/^Pageouts:/ {
                gsub(/\./, "", $2); gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
                print $2; exit
            }')
    fi
elif [ "$kernel" = Linux ] && [ -r /proc/meminfo ]; then
    memory_free_percent=$(awk '
        /^MemTotal:/ {total=$2}
        /^MemAvailable:/ {available=$2}
        END {if (total > 0) printf "%.1f", 100 * available / total}
    ' /proc/meminfo)
    pageouts=$(awk '/^pswpout / {print $2; exit}' /proc/vmstat 2>/dev/null || true)
    if command -v getconf >/dev/null 2>&1; then
        page_size_bytes=$(getconf PAGESIZE 2>/dev/null || true)
    fi
fi

for variable in power_source battery_state low_power_mode thermal_state \
                memory_free_percent pageouts page_size_bytes load_averages; do
    eval value=\${$variable}
    [ -n "$value" ] || eval "$variable=\$UNKNOWN"
done

emit captured_at "$captured_at"
emit snapshot_epoch "$snapshot_epoch"
emit power_source "$power_source"
emit battery_state "$battery_state"
emit low_power_mode "$low_power_mode"
emit thermal_state "$thermal_state"
emit memory_free_percent "$memory_free_percent"
emit pageouts "$pageouts"
emit page_size_bytes "$page_size_bytes"
emit load_averages "$load_averages"
