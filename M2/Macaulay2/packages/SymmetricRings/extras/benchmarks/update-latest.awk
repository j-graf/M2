#!/usr/bin/awk -f
# Usage: RECORD_RUN_NAME=... RECORD_RECORDED_AT=... RECORD_DEVICE_PROFILE=... \
#   update-latest.awk latest.tsv summary.tsv

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
        position[key] = ++resultCount
        orderedKey[resultCount] = key
    }
    caseId[key] = $1
    coefficientRing[key] = $4
    resultMedian[key] = $3 + 0
    resultRuns[key] = $5
    resultDate[key] = $2
    resultRun[key] = $6
    resultDevice[key] = NF >= 7 && $7 != "" ? $7 : deviceProfile
    next
}

FILENAME == ARGV[2] {
    if (FNR == 1 || NF < 9) next
    key = $1 SUBSEP $5
    if (!(key in position)) {
        position[key] = ++resultCount
        orderedKey[resultCount] = key
    }
    caseId[key] = $1
    coefficientRing[key] = $5
    resultMedian[key] = $8 + 0
    resultRuns[key] = $6 + 0
    resultDate[key] = recordedAt
    resultRun[key] = runName
    resultDevice[key] = deviceProfile
    next
}

END {
    print "case_id", "recorded_at", "median_cpu_seconds", \
        "coefficient_ring", "runs", "run_name", "device_profile"
    for (i = 1; i <= resultCount; ++i) {
        key = orderedKey[i]
        print caseId[key], resultDate[key], sprintf("%.9g", resultMedian[key]), \
            coefficientRing[key], resultRuns[key], resultRun[key], resultDevice[key]
    }
}
