#!/usr/bin/awk -f
# Usage: RECORD_RUN_NAME=... RECORD_RECORDED_AT=... RECORD_DEVICE_PROFILE=... \
#   update-records.awk records.tsv summary.tsv

BEGIN {
    FS = "\t"
    OFS = "\t"
    runName = ENVIRON["RECORD_RUN_NAME"]
    recordedAt = ENVIRON["RECORD_RECORDED_AT"]
    deviceProfile = ENVIRON["RECORD_DEVICE_PROFILE"]
    if (runName == "") runName = "unknown"
    if (recordedAt == "") recordedAt = "unknown"
    if (deviceProfile == "") deviceProfile = "unknown"
}

FILENAME == ARGV[1] {
    if (FNR == 1) next
    if (NF < 6) next
    key = $1 SUBSEP $4
    if (!(key in position)) {
        position[key] = ++recordCount
        orderedKey[recordCount] = key
        caseId[key] = $1
        coefficientRing[key] = $4
        recordMedian[key] = $3 + 0
        recordRuns[key] = $5
        recordDate[key] = $2
        recordRun[key] = $6
        recordDevice[key] = NF >= 7 && $7 != "" ? $7 : deviceProfile
    } else if (($3 + 0) < recordMedian[key]) {
        recordMedian[key] = $3 + 0
        recordRuns[key] = $5
        recordDate[key] = $2
        recordRun[key] = $6
        recordDevice[key] = NF >= 7 && $7 != "" ? $7 : deviceProfile
    }
    next
}

FILENAME == ARGV[2] {
    if (FNR == 1 || NF < 9) next
    runs = $6 + 0
    if (runs < 3) next
    key = $1 SUBSEP $5
    median = $8 + 0
    if (!(key in position)) {
        position[key] = ++recordCount
        orderedKey[recordCount] = key
        caseId[key] = $1
        coefficientRing[key] = $5
        recordMedian[key] = median
        recordRuns[key] = runs
        recordDate[key] = recordedAt
        recordRun[key] = runName
        recordDevice[key] = deviceProfile
    } else if (median < recordMedian[key]) {
        recordMedian[key] = median
        recordRuns[key] = runs
        recordDate[key] = recordedAt
        recordRun[key] = runName
        recordDevice[key] = deviceProfile
    }
    next
}

END {
    print "case_id", "recorded_at", "median_cpu_seconds", \
        "coefficient_ring", "runs", "run_name", "device_profile"
    for (i = 1; i <= recordCount; ++i) {
        key = orderedKey[i]
        print caseId[key], recordDate[key], sprintf("%.9g", recordMedian[key]), \
            coefficientRing[key], recordRuns[key], recordRun[key], recordDevice[key]
    }
}
