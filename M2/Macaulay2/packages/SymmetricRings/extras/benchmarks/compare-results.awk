#!/usr/bin/awk -f
# Usage: compare-results.awk records.tsv raw-results.tsv

BEGIN {
    FS = "\t"; OFS = "\t"
    threshold = ENVIRON["THRESHOLD_PERCENT"] + 0
    if (threshold <= 0) threshold = 10
    inputFile = 0
}
FNR == 1 { ++inputFile; next }
inputFile == 1 && NF >= 6 {
    key = $1 SUBSEP $4
    value = $3 + 0
    if (!(key in fastestRecord) || value < fastestRecord[key])
        fastestRecord[key] = value
    next
}
inputFile == 2 && NF >= 11 {
    key = $1 SUBSEP $5
    n = ++counts[key]
    values[key, n] = $7 + 0
    caseId[key] = $1; family[key] = $2; operation[key] = $3
    tier[key] = $4; coefficientRing[key] = $5
}
END {
    print "case_id", "family", "operation", "tier", "coefficient_ring", \
          "runs", "current_median", "fastest_record", \
          "change_vs_record_percent", "record_status", "record_classification"
    for (key in counts) {
        n = counts[key]
        for (i = 1; i <= n; ++i) sorted[i] = values[key, i]
        for (i = 2; i <= n; ++i) {
            x = sorted[i]; j = i - 1
            while (j >= 1 && sorted[j] > x) { sorted[j + 1] = sorted[j]; --j }
            sorted[j + 1] = x
        }
        current = n % 2 == 1 ? sorted[(n + 1) / 2] : (sorted[n / 2] + sorted[n / 2 + 1]) / 2
        if (!(key in fastestRecord)) {
            recordText = recordChangeText = "NA"
            recordStatus = n >= 3 ? "first-record" : "ineligible"
            recordClassification = "new"
        } else {
            recordText = sprintf("%.9g", fastestRecord[key])
            recordChange = fastestRecord[key] == 0 ? 0 : \
                100 * (current / fastestRecord[key] - 1)
            recordChangeText = sprintf("%+.2f", recordChange)
            recordClassification = recordChange >= threshold ? "regression" : \
                (recordChange <= -threshold ? "improvement" : "stable")
            recordStatus = n < 3 ? "ineligible" : \
                (current < fastestRecord[key] ? "new-record" : \
                    (current == fastestRecord[key] ? "ties-record" : "not-record"))
        }
        print caseId[key], family[key], operation[key], tier[key], coefficientRing[key], \
              n, sprintf("%.9g", current), recordText, recordChangeText, \
              recordStatus, recordClassification
        for (i = 1; i <= n; ++i) delete sorted[i]
    }
}
