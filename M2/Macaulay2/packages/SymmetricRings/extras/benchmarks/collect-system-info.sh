#!/bin/sh
set -eu

M2_BIN=${1:-M2}
SOURCE_ROOT=${2:-.}
DEVICE_PROFILE=${3:-unknown}

clean_value() {
    printf '%s' "$1" | tr '\t\r\n' '   '
}

emit() {
    printf '%s\t%s\n' "$1" "$(clean_value "$2")"
}

unknown=unknown
kernel_name=$(uname -s 2>/dev/null || printf '%s' "$unknown")
kernel_release=$(uname -r 2>/dev/null || printf '%s' "$unknown")
architecture=$(uname -m 2>/dev/null || printf '%s' "$unknown")

computer_model=$unknown
model_identifier=$unknown
cpu_model=$unknown
physical_cores=$unknown
logical_cores=$unknown
memory=$unknown
memory_bytes=$unknown
os_name=$kernel_name
os_version=$kernel_release
os_build=$unknown

if [ "$kernel_name" = Darwin ]; then
    if command -v system_profiler >/dev/null 2>&1; then
        hardware=$(system_profiler SPHardwareDataType 2>/dev/null || true)
        computer_model=$(printf '%s\n' "$hardware" | awk -F ': ' '/^[[:space:]]*Model Name:/ {print $2; exit}')
        model_identifier=$(printf '%s\n' "$hardware" | awk -F ': ' '/^[[:space:]]*Model Identifier:/ {print $2; exit}')
        cpu_model=$(printf '%s\n' "$hardware" | awk -F ': ' '/^[[:space:]]*Chip:/ {print $2; exit}')
        logical_cores=$(printf '%s\n' "$hardware" | awk -F ': ' '/^[[:space:]]*Total Number of Cores:/ {print $2; exit}')
        physical_cores=$logical_cores
        memory=$(printf '%s\n' "$hardware" | awk -F ': ' '/^[[:space:]]*Memory:/ {print $2; exit}')
    fi
    if command -v sysctl >/dev/null 2>&1; then
        memory_bytes=$(sysctl -n hw.memsize 2>/dev/null || printf '%s' "$unknown")
    fi
    if command -v sw_vers >/dev/null 2>&1; then
        os_name=$(sw_vers -productName 2>/dev/null || printf '%s' macOS)
        os_version=$(sw_vers -productVersion 2>/dev/null || printf '%s' "$unknown")
        os_build=$(sw_vers -buildVersion 2>/dev/null || printf '%s' "$unknown")
    fi
elif [ "$kernel_name" = Linux ]; then
    if [ -r /sys/devices/virtual/dmi/id/product_name ]; then
        computer_model=$(sed -n '1p' /sys/devices/virtual/dmi/id/product_name)
    fi
    if [ -r /sys/devices/virtual/dmi/id/product_version ]; then
        model_identifier=$(sed -n '1p' /sys/devices/virtual/dmi/id/product_version)
    fi
    if [ -r /proc/cpuinfo ]; then
        cpu_model=$(awk -F ': ' '/^(model name|Hardware|Processor)[[:space:]]*:/ {print $2; exit}' /proc/cpuinfo)
        logical_cores=$(awk '/^processor[[:space:]]*:/ {n++} END {print n+0}' /proc/cpuinfo)
        physical_cores=$(awk -F ': ' '
            /^physical id[[:space:]]*:/ {p=$2}
            /^core id[[:space:]]*:/ {c=$2; seen[p ":" c]=1}
            END {n=0; for (x in seen) n++; if (n>0) print n; else print "unknown"}
            ' /proc/cpuinfo)
    fi
    if [ -r /proc/meminfo ]; then
        memory_kib=$(awk '/^MemTotal:/ {print $2; exit}' /proc/meminfo)
        if [ -n "$memory_kib" ]; then
            memory_bytes=$((memory_kib * 1024))
            memory=$(awk -v kib="$memory_kib" 'BEGIN {printf "%.1f GiB", kib/1048576}')
        fi
    fi
    if [ -r /etc/os-release ]; then
        os_name=$(awk -F= '/^NAME=/ {gsub(/^"|"$/, "", $2); print $2; exit}' /etc/os-release)
        os_version=$(awk -F= '/^VERSION_ID=/ {gsub(/^"|"$/, "", $2); print $2; exit}' /etc/os-release)
    fi
fi

for variable in computer_model model_identifier cpu_model physical_cores logical_cores memory memory_bytes os_name os_version os_build; do
    eval current_value=\${$variable}
    if [ -z "$current_value" ]; then eval "$variable=\$unknown"; fi
done

m2_version=$($M2_BIN --version 2>&1 | sed -n '1p' || printf '%s' "$unknown")
git_commit=$(git -C "$SOURCE_ROOT" rev-parse HEAD 2>/dev/null || printf '%s' "$unknown")
git_branch=$(git -C "$SOURCE_ROOT" branch --show-current 2>/dev/null || printf '%s' "$unknown")
if git -C "$SOURCE_ROOT" diff --quiet --ignore-submodules HEAD 2>/dev/null; then
    git_dirty=false
else
    git_dirty=true
fi
build_type=$unknown
if [ -r "$SOURCE_ROOT/BUILD/build/CMakeCache.txt" ]; then
    build_type=$(awk -F= '/^CMAKE_BUILD_TYPE:/ {print $2; exit}' "$SOURCE_ROOT/BUILD/build/CMakeCache.txt")
    [ -n "$build_type" ] || build_type=$unknown
fi

printf 'field\tvalue\n'
emit captured_at "$(date '+%Y-%m-%d %H:%M:%S %Z')"
emit device_profile "$DEVICE_PROFILE"
emit computer_model "$computer_model"
emit model_identifier "$model_identifier"
emit cpu_model "$cpu_model"
emit physical_cores "$physical_cores"
emit logical_cores "$logical_cores"
emit memory "$memory"
emit memory_bytes "$memory_bytes"
emit os_name "$os_name"
emit os_version "$os_version"
emit os_build "$os_build"
emit kernel "$kernel_name $kernel_release"
emit architecture "$architecture"
emit m2_version "$m2_version"
emit m2_executable "$M2_BIN"
emit git_commit "$git_commit"
emit git_branch "$git_branch"
emit git_dirty "$git_dirty"
emit cmake_build_type "$build_type"
