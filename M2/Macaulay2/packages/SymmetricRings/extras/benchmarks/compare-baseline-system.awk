#!/usr/bin/awk -f
# Usage: compare-baseline-system.awk system.tsv baselines.tsv comparison.tsv

BEGIN {
    FS = "\t"
    OFS = "\t"
}

FILENAME == ARGV[1] {
    if (FNR > 1) current[$1] = $2
    next
}

FILENAME == ARGV[2] {
    if (FNR == 1) {
        for (i = 1; i <= NF; ++i) baselineColumn[$i] = i
        next
    }
    statusColumn = baselineColumn["status"]
    status = statusColumn ? $(statusColumn) : ""
    if (status != "accepted" && status != "valid") next
    caseColumn = baselineColumn["case_id"]
    ringColumn = baselineColumn["coefficient_ring"]
    dateColumn = baselineColumn["accepted_date"]
    if (!caseColumn || !ringColumn || !dateColumn) next
    key = $(caseColumn) SUBSEP $(ringColumn)
    date = $(dateColumn)
    if (!(key in latestDate) || date > latestDate[key]) {
        latestDate[key] = date
        latestCpu[key] = baselineColumn["cpu_model"] ? \
            $(baselineColumn["cpu_model"]) : "unknown"
        latestMemory[key] = baselineColumn["memory"] ? \
            $(baselineColumn["memory"]) : "unknown"
        latestOsName[key] = baselineColumn["os_name"] ? \
            $(baselineColumn["os_name"]) : "unknown"
        latestOsVersion[key] = baselineColumn["os_version"] ? \
            $(baselineColumn["os_version"]) : "unknown"
    }
    next
}

FILENAME == ARGV[3] {
    if (FNR == 1) {
        for (i = 1; i <= NF; ++i) comparisonColumn[$i] = i
        next
    }
    caseColumn = comparisonColumn["case_id"]
    ringColumn = comparisonColumn["coefficient_ring"]
    if (!caseColumn || !ringColumn) next
    key = $(caseColumn) SUBSEP $(ringColumn)
    if (!(key in latestDate)) next
    selectedBaselineCount++
    cpu = latestCpu[key]
    memory = latestMemory[key]
    osName = latestOsName[key]
    osVersion = latestOsVersion[key]
    if (!known(cpu) || !known(memory) || !known(osName) || !known(osVersion)) {
        missingConfiguration = 1
        next
    }
    configuration = cpu SUBSEP memory SUBSEP osName SUBSEP osVersion
    configurations[configuration] = 1
    representativeCpu = cpu
    representativeMemory = memory
    representativeOsName = osName
    representativeOsVersion = osVersion
    next
}

function known(value) {
    return value != "" && value != "unknown"
}

function fieldResult(currentValue, baselineValue) {
    if (!known(currentValue) || !known(baselineValue)) return "unknown"
    return currentValue == baselineValue ? "same" : "different"
}

END {
    configurationCount = 0
    for (configuration in configurations) configurationCount++

    currentCpu = current["cpu_model"]
    currentMemory = current["memory"]
    currentOsName = current["os_name"]
    currentOsVersion = current["os_version"]
    currentOs = currentOsName " " currentOsVersion

    if (selectedBaselineCount == 0 || missingConfiguration || \
        configurationCount == 0) {
        overall = "unknown"
    } else if (configurationCount > 1) {
        overall = "mixed"
    } else {
        cpuResult = fieldResult(currentCpu, representativeCpu)
        memoryResult = fieldResult(currentMemory, representativeMemory)
        osNameResult = fieldResult(currentOsName, representativeOsName)
        osVersionResult = fieldResult(currentOsVersion, representativeOsVersion)
        if (cpuResult == "unknown" || memoryResult == "unknown" || \
            osNameResult == "unknown" || osVersionResult == "unknown")
            overall = "unknown"
        else if (cpuResult == "same" && memoryResult == "same" && \
                 osNameResult == "same" && osVersionResult == "same")
            overall = "same"
        else
            overall = "different"
    }

    print "field", "current", "baseline", "comparison"
    print "overall", "-", "-", overall

    if (configurationCount > 1) {
        baselineCpu = baselineMemory = baselineOs = \
            "multiple baseline configurations"
        cpuResult = memoryResult = osResult = "mixed"
    } else {
        baselineCpu = known(representativeCpu) ? representativeCpu : "unknown"
        baselineMemory = known(representativeMemory) ? representativeMemory : "unknown"
        baselineOs = known(representativeOsName) && known(representativeOsVersion) ? \
            representativeOsName " " representativeOsVersion : "unknown"
        cpuResult = fieldResult(currentCpu, baselineCpu)
        memoryResult = fieldResult(currentMemory, baselineMemory)
        osResult = fieldResult(currentOs, baselineOs)
    }

    print "CPU", known(currentCpu) ? currentCpu : "unknown", baselineCpu, cpuResult
    print "RAM", known(currentMemory) ? currentMemory : "unknown", \
        baselineMemory, memoryResult
    print "Operating system", \
        known(currentOsName) && known(currentOsVersion) ? currentOs : "unknown", \
        baselineOs, osResult
}

