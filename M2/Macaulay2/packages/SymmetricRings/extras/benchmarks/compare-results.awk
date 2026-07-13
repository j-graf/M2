#!/usr/bin/awk -f
# Usage: compare-results.awk baselines.tsv records.tsv raw-results.tsv

BEGIN {
    FS = "\t"; OFS = "\t"
    threshold = ENVIRON["THRESHOLD_PERCENT"] + 0
    if (threshold <= 0) threshold = 10
    inputFile = 0
}
FNR == 1 { ++inputFile; next }
inputFile == 1 && NF >= 5 {
    status = $5
    if (status != "accepted" && status != "valid") next
    key = $1 SUBSEP $4
    value = $3 + 0
    if (!(key in bestBaseline) || value < bestBaseline[key]) bestBaseline[key] = value
    if (!(key in recentDate) || $2 > recentDate[key]) {
        recentDate[key] = $2
        recentBaseline[key] = value
    }
    next
}
inputFile == 2 && NF >= 6 {
    key = $1 SUBSEP $4
    value = $3 + 0
    if (!(key in fastestRecord) || value < fastestRecord[key])
        fastestRecord[key] = value
    next
}
inputFile == 3 && NF >= 11 {
    key = $1 SUBSEP $5
    n = ++counts[key]
    values[key, n] = $7 + 0
    caseId[key] = $1; family[key] = $2; operation[key] = $3
    tier[key] = $4; coefficientRing[key] = $5
}
END {
    print "case_id", "family", "operation", "tier", "coefficient_ring", \
          "runs", "current_median", "best_valid_baseline", \
          "change_vs_best_percent", "most_recent_baseline", \
          "change_vs_recent_percent", "fastest_record", \
          "change_vs_record_percent", "record_status", "classification"
    for (key in counts) {
        n = counts[key]
        for (i = 1; i <= n; ++i) sorted[i] = values[key, i]
        for (i = 2; i <= n; ++i) {
            x = sorted[i]; j = i - 1
            while (j >= 1 && sorted[j] > x) { sorted[j + 1] = sorted[j]; --j }
            sorted[j + 1] = x
        }
        current = n % 2 == 1 ? sorted[(n + 1) / 2] : (sorted[n / 2] + sorted[n / 2 + 1]) / 2
        if (!(key in bestBaseline)) {
            bestText = recentText = bestChangeText = recentChangeText = "NA"
            classification = "new"
        } else {
            bestText = sprintf("%.9g", bestBaseline[key])
            recentText = sprintf("%.9g", recentBaseline[key])
            bestChange = bestBaseline[key] == 0 ? 0 : 100 * (current / bestBaseline[key] - 1)
            recentChange = recentBaseline[key] == 0 ? 0 : 100 * (current / recentBaseline[key] - 1)
            bestChangeText = sprintf("%+.2f", bestChange)
            recentChangeText = sprintf("%+.2f", recentChange)
            classification = bestChange >= threshold ? "regression" : \
                (bestChange <= -threshold ? "improvement" : "stable")
        }
        if (!(key in fastestRecord)) {
            recordText = recordChangeText = "NA"
            recordStatus = n >= 3 ? "first-record" : "ineligible"
        } else {
            recordText = sprintf("%.9g", fastestRecord[key])
            recordChange = fastestRecord[key] == 0 ? 0 : \
                100 * (current / fastestRecord[key] - 1)
            recordChangeText = sprintf("%+.2f", recordChange)
            recordStatus = n < 3 ? "ineligible" : \
                (current < fastestRecord[key] ? "new-record" : \
                    (current == fastestRecord[key] ? "ties-record" : "not-record"))
        }
        print caseId[key], family[key], operation[key], tier[key], coefficientRing[key], \
              n, sprintf("%.9g", current), bestText, bestChangeText, recentText, \
              recentChangeText, recordText, recordChangeText, recordStatus, classification
        for (i = 1; i <= n; ++i) delete sorted[i]
    }
}
