#!/usr/bin/awk -f

BEGIN {
    FS = "\t"
    status = "clean"
}

function markdown(value) {
    gsub(/\|/, "\\|", value)
    return value
}

function warn(message) {
    if (status == "clean") status = "warning"
    notes[++note_count] = message
}

function compromise(message) {
    status = "compromised"
    notes[++note_count] = message
}

function absolute(value) { return value < 0 ? -value : value }

NR == 1 { next }

$1 == "configuration" {
    configuration[$4] = $5
    next
}

$1 == "snapshot" {
    snapshot[$3, $4] = $5
    next
}

$1 == "calibration" {
    family = $2
    if (!(family in seen_family)) {
        seen_family[family] = 1
        families[++family_count] = family
    }
    probe[family, $3, $4] = $5
}

END {
    for (i = 1; i <= family_count; i++) {
        family = families[i]
        before = probe[family, "before", "cpu_seconds"] + 0
        after = probe[family, "after", "cpu_seconds"] + 0
        if (before <= 0 || after <= 0) {
            drift[family] = "NA"
            warn("A calibration probe was missing for " family ".")
        } else {
            change = 100 * (after / before - 1)
            drift[family] = sprintf("%+.2f", change)
            magnitude = absolute(change)
            if (magnitude >= 20)
                compromise(family " calibration drifted " drift[family] "%.")
            else if (magnitude >= 10)
                warn(family " calibration drifted " drift[family] "%.")
        }
    }

    for (phase_index = 1; phase_index <= 2; phase_index++) {
        phase = phase_index == 1 ? "before" : "after"
        thermal = snapshot[phase, "thermal_state"]
        if (thermal == "serious" || thermal == "critical")
            compromise("Thermal state was " thermal " " phase " the run.")
        else if (thermal == "fair")
            warn("Thermal state was fair " phase " the run.")

        low_power = snapshot[phase, "low_power_mode"]
        if (low_power == "1" || low_power == "true" || low_power == "enabled")
            compromise("Low Power Mode was enabled " phase " the run.")

        free_memory = snapshot[phase, "memory_free_percent"]
        if (free_memory != "" && free_memory != "unknown") {
            free_memory += 0
            if (free_memory < 5)
                compromise("Free-memory percentage was below 5% " phase " the run.")
            else if (free_memory < 15)
                warn("Free-memory percentage was below 15% " phase " the run.")
        }
    }

    before_pageouts = snapshot["before", "pageouts"]
    after_pageouts = snapshot["after", "pageouts"]
    if (before_pageouts != "" && after_pageouts != "" &&
        before_pageouts != "unknown" && after_pageouts != "unknown") {
        pageout_delta = after_pageouts + 0 - (before_pageouts + 0)
        if (pageout_delta < 0) pageout_delta = 0
        page_size = snapshot["after", "page_size_bytes"]
        if (page_size == "" || page_size == "unknown" || page_size + 0 <= 0)
            page_size = snapshot["before", "page_size_bytes"]
        # Historical macOS condition files predate explicit page-size capture.
        if (page_size == "" || page_size == "unknown" || page_size + 0 <= 0)
            page_size = 16384
        pageout_mib = pageout_delta * (page_size + 0) / 1048576

        before_epoch = snapshot["before", "snapshot_epoch"]
        after_epoch = snapshot["after", "snapshot_epoch"]
        elapsed = after_epoch + 0 - (before_epoch + 0)
        if (before_epoch == "" || after_epoch == "" ||
            before_epoch == "unknown" || after_epoch == "unknown" || elapsed <= 0)
            elapsed = 0
        pageout_rate = elapsed > 0 ? pageout_mib / elapsed : 0

        before_free = snapshot["before", "memory_free_percent"]
        after_free = snapshot["after", "memory_free_percent"]
        minimum_free = 100
        if (before_free != "" && before_free != "unknown") minimum_free = before_free + 0
        if (after_free != "" && after_free != "unknown" && after_free + 0 < minimum_free)
            minimum_free = after_free + 0

        pageout_concern = elapsed > 0 \
            ? sprintf("Pageout activity was %.2f MiB at %.2f MiB/s.",
                      pageout_mib, pageout_rate) \
            : sprintf("Pageout activity was %.2f MiB.", pageout_mib)

        if (pageout_mib >= 512 && (pageout_rate >= 5 || minimum_free < 5))
            compromise(pageout_concern)
        else if (pageout_mib >= 64 && (pageout_rate >= 1 || minimum_free < 15))
            warn(pageout_concern)
        if (elapsed > 0)
            pageout_summary = sprintf("%d pages (%.2f MiB, %.2f MiB/s)",
                                      pageout_delta, pageout_mib, pageout_rate)
        else
            pageout_summary = sprintf("%d pages (%.2f MiB)",
                                      pageout_delta, pageout_mib)
    }

    if (family_count == 0)
        warn("No calibration probes were recorded.")

    for (i = 1; i <= 8; i++) {
        key = i == 1 ? "varied_mode" : (i == 2 ? "varied_level" : \
            (i == 3 ? "random_seed" : (i == 4 ? "new_only" : \
            (i == 5 ? "new_reference_history" : \
            (i == 6 ? "device_profile" : \
            (i == 7 ? "records_history" : "latest_history"))))))
        if (configuration[key] == "") configuration[key] = "unknown"
    }

    printf "## Run conditions\n\n"
    printf "Overall health: **%s**.\n\n", status
    printf "- **Clean:** conditions were stable and no comparison concern was detected.\n"
    printf "- **Warning:** moderate drift or resource pressure was detected; review performance conclusions.\n"
    printf "- **Compromised:** severe drift, thermal pressure, or a performance-limiting mode makes the run unsuitable for performance conclusions.\n\n"

    if (note_count == 0)
        printf "No run-condition warnings were detected.\n\n"
    else {
        printf "Observed concerns:\n\n"
        for (i = 1; i <= note_count; i++) printf "- %s\n", notes[i]
        printf "\n"
    }

    if (pageout_summary != "")
        printf "Pageout activity: %s.\n\n", pageout_summary

    printf "### Benchmark selection\n\n"
    printf "| Device | Mode | Level | Random seed | New only | New-case history |\n"
    printf "|---|---|---|---:|---|---|\n"
    printf "| %s | %s | %s | %s | %s | %s |\n\n", \
        configuration["device_profile"], configuration["varied_mode"], \
        configuration["varied_level"], configuration["random_seed"], \
        configuration["new_only"], configuration["new_reference_history"]

    printf "History sources: fastest records `%s`; latest per-case results `%s`.\n\n", \
        configuration["records_history"], configuration["latest_history"]

    printf "### Condition snapshots\n\n"
    printf "| Metric | Before | After |\n|---|---|---|\n"
    split("captured_at snapshot_epoch power_source battery_state low_power_mode thermal_state memory_free_percent pageouts page_size_bytes load_averages", metrics, " ")
    for (i = 1; i <= 10; i++) {
        before_value = snapshot["before", metrics[i]]
        after_value = snapshot["after", metrics[i]]
        if (before_value == "") before_value = "unknown"
        if (after_value == "") after_value = "unknown"
        printf "| %s | %s | %s |\n", metrics[i], \
            markdown(before_value), markdown(after_value)
    }

    printf "\n### Family calibration probes\n\n"
    printf "The fixed `conv-S-three-p` probe runs in a fresh process immediately before and after each family; no monitor runs concurrently with timed cases.\n\n"
    printf "| Family | Before CPU (s) | After CPU (s) | Drift | Before wall (s) | After wall (s) |\n"
    printf "|---|---:|---:|---:|---:|---:|\n"
    for (i = 1; i <= family_count; i++) {
        family = families[i]
        drift_text = drift[family] == "NA" ? "NA" : drift[family] "%"
        printf "| %s | %s | %s | %s | %s | %s |\n", markdown(family), \
            probe[family, "before", "cpu_seconds"], \
            probe[family, "after", "cpu_seconds"], drift_text, \
            probe[family, "before", "wall_seconds"], \
            probe[family, "after", "wall_seconds"]
    }
    printf "\n"
}
