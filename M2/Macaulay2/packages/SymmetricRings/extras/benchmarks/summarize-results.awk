#!/usr/bin/awk -f
BEGIN {
    FS = "\t"
    OFS = "\t"
    print "case_id", "family", "operation", "tier", "coefficient_ring", \
          "runs", "min_cpu_seconds", "median_cpu_seconds", "max_cpu_seconds"
}
NR == 1 { next }
{
    id = $1
    n = ++counts[id]
    values[id, n] = $7 + 0
    family[id] = $2
    operation[id] = $3
    tier[id] = $4
    coefficientRing[id] = $5
}
END {
    for (id in counts) {
        n = counts[id]
        for (i = 1; i <= n; ++i) sorted[i] = values[id, i]
        for (i = 2; i <= n; ++i) {
            x = sorted[i]
            j = i - 1
            while (j >= 1 && sorted[j] > x) {
                sorted[j + 1] = sorted[j]
                --j
            }
            sorted[j + 1] = x
        }
        if (n % 2 == 1) median = sorted[(n + 1) / 2]
        else median = (sorted[n / 2] + sorted[n / 2 + 1]) / 2
        print id, family[id], operation[id], tier[id], coefficientRing[id], \
              n, sorted[1], median, sorted[n]
        for (i = 1; i <= n; ++i) delete sorted[i]
    }
}

