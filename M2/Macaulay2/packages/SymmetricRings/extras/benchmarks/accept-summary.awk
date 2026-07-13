#!/usr/bin/awk -f
# Usage: accept-summary.awk system.tsv summary.tsv
BEGIN {
    FS = "\t"; OFS = "\t"
    if (ARGC < 3) {
        print "usage: accept-summary.awk system.tsv summary.tsv" > "/dev/stderr"
        exit 2
    }
    acceptedDate = ENVIRON["ACCEPTED_DATE"]
    notes = ENVIRON["NOTES"]
    if (acceptedDate == "") {
        print "ACCEPTED_DATE is required" > "/dev/stderr"
        exit 2
    }
}
FILENAME == ARGV[1] {
    if (FNR > 1) systemFields[$1] = $2
    next
}
FILENAME == ARGV[2] && FNR == 1 { next }
FILENAME == ARGV[2] && NF >= 9 {
    print $1, acceptedDate, $8, $5, "accepted", notes, \
        systemValue("cpu_model"), systemValue("memory"), \
        systemValue("os_name"), systemValue("os_version")
}

function systemValue(field) {
    return systemFields[field] == "" ? "unknown" : systemFields[field]
}
