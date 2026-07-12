#!/usr/bin/awk -f
BEGIN {
    FS = "\t"; OFS = "\t"
    acceptedDate = ENVIRON["ACCEPTED_DATE"]
    notes = ENVIRON["NOTES"]
    if (acceptedDate == "") {
        print "ACCEPTED_DATE is required" > "/dev/stderr"
        exit 2
    }
}
NR == 1 { next }
NF >= 9 { print $1, acceptedDate, $8, $5, "accepted", notes }

