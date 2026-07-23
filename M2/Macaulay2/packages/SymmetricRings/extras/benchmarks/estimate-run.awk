#!/usr/bin/awk -f

function duration(seconds) {
    if (seconds < 1) return sprintf("%.3f s", seconds)
    if (seconds < 60) return sprintf("%.1f s", seconds)
    if (seconds < 3600) return sprintf("%.1f min", seconds / 60)
    return sprintf("%.2f h", seconds / 3600)
}

FILENAME == recordFile {
    if (FNR == 1) next
    split($0, field, "\t")
    key = field[1] SUBSEP field[4]
    recordValue[key] = field[3] + 0
    next
}

FILENAME == latestFile {
    if (FNR == 1) next
    split($0, field, "\t")
    if (length(field) < 6) next
    key = field[1] SUBSEP field[4]
    latestValue[key] = field[3] + 0
    next
}

FILENAME == planFile {
    split($0, field, "|")
    caseCount++
    caseId[caseCount] = field[1]
    caseFamily[caseCount] = field[2]
    caseRing[caseCount] = field[3]
    if (!familySeen[field[2]]++) families[++familyCount] = field[2]
}

END {
    for (i = 1; i <= caseCount; i++) {
        key = caseId[i] SUBSEP caseRing[i]
        family = caseFamily[i]
        familyCases[family]++
        if (key in latestValue) {
            estimate = latestValue[key]
            latestSources++
        } else if (key in recordValue) {
            estimate = recordValue[key]
            recordSources++
        } else {
            unknownSources++
            familyUnknown[family]++
            continue
        }
        familyWork[family] += repetitions * estimate
    }

    printf "Benchmark time estimate (no tests run)\n\n"
    printf "Device profile: %s.\n", deviceProfile
    caseWord = caseCount == 1 ? "case" : "cases"
    familyWord = familyCount == 1 ? "family" : "families"
    printf "Selection: %d %s across %d %s, %d repetitions.\n", \
        caseCount, caseWord, familyCount, familyWord, repetitions
    printf "Latest per-case history: %s.\n\n", latestLabel
    printf "| Family | Cases | Timed work | Process overhead | Calibration | Estimated wall | Unknown |\n"
    printf "|---|---:|---:|---:|---:|---:|---:|\n"
    for (i = 1; i <= familyCount; i++) {
        family = families[i]
        processes = familyCases[family] * repetitions + 2
        overhead = processes * processOverhead
        calibration = 2 * calibrationSeconds
        estimate = familyWork[family] + overhead + calibration
        totalWork += familyWork[family]
        totalOverhead += overhead
        totalCalibration += calibration
        printf "| %s | %d | %s | %s | %s | %s | %d |\n", family, \
            familyCases[family], duration(familyWork[family]), duration(overhead), \
            duration(calibration), duration(estimate), familyUnknown[family] + 0
    }
    total = totalWork + totalOverhead + totalCalibration
    printf "\nEstimated total wall time: **%s**.\n\n", duration(total)
    printf "Timing sources: %d latest per-case medians, %d fastest records, %d unknown.\n\n", \
        latestSources, recordSources, unknownSources
    printf "Model assumptions: %.2f s fresh-process overhead per worker and %.2f s per calibration computation. ", \
        processOverhead, calibrationSeconds
    printf "Latest per-case CPU medians are preferred, with fastest records as fallback. "
    printf "Override the model with SYMRINGS_BENCH_PROCESS_OVERHEAD_SECONDS or SYMRINGS_BENCH_CALIBRATION_SECONDS.\n"
}
