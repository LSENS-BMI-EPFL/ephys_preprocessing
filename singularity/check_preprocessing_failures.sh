#!/bin/bash

# Script to detect Ephys sessions where spike sorting preprocessing failed
# Checks that all expected output folders/files exist and are non-empty
# per probe: session_id/Ephys/catgt_*/*imec*/kilosort4/...

DATA_DIR="/scratch/bisi/data"
OUTPUT_FILE="/home/bisi/code/ephys_preprocessing/config/failed_sessions_axel.txt"
SUFFIX="tcat_corrected"  # Must match detect_ephys_sessions.sh

echo "========================================"
echo "Checking preprocessing outputs in: $DATA_DIR"
echo "========================================"
echo ""

# Temp file used to pass failure messages out of check_probe_outputs
# (avoids subshell variable-scoping issues with process substitution)
_PROBE_FAIL_LOG=$(mktemp)
trap 'rm -f "$_PROBE_FAIL_LOG"' EXIT

# -----------------------------------------------------------------------
# Helpers — print on failure only; write short reason to _PROBE_FAIL_LOG
# -----------------------------------------------------------------------
check_dir_nonempty() {
    local label="$1" path="$2" indent="$3"
    if [ ! -d "$path" ]; then
        echo "${indent}  ✗ MISSING dir  : $label"
        echo "MISSING dir: $label" >> "$_PROBE_FAIL_LOG"
        return 1
    fi
    if [ -z "$(find "$path" -mindepth 1 -maxdepth 1 2>/dev/null | head -1)" ]; then
        echo "${indent}  ✗ EMPTY dir    : $label"
        echo "EMPTY dir: $label" >> "$_PROBE_FAIL_LOG"
        return 1
    fi
    return 0
}

check_file_exists() {
    local label="$1" dir="$2" pattern="$3" indent="$4"
    if [ ! -d "$dir" ]; then
        echo "${indent}  ✗ MISSING dir  : $label (parent dir absent)"
        echo "MISSING dir: $label (parent absent)" >> "$_PROBE_FAIL_LOG"
        return 1
    fi
    local match
    match=$(find "$dir" -maxdepth 1 -name "$pattern" 2>/dev/null | head -1)
    if [ -z "$match" ]; then
        echo "${indent}  ✗ MISSING file : $label"
        echo "MISSING file: $label" >> "$_PROBE_FAIL_LOG"
        return 1
    fi
    return 0
}

check_glob_nonempty() {
    local label="$1" dir="$2" pattern="$3" indent="$4"
    if [ ! -d "$dir" ]; then
        echo "${indent}  ✗ MISSING dir  : $label (parent dir absent)"
        echo "MISSING dir: $label (parent absent)" >> "$_PROBE_FAIL_LOG"
        return 1
    fi
    local match
    match=$(find "$dir" -maxdepth 1 -name "$pattern" 2>/dev/null | head -1)
    if [ -z "$match" ]; then
        echo "${indent}  ✗ MISSING      : $label"
        echo "MISSING: $label" >> "$_PROBE_FAIL_LOG"
        return 1
    fi
    return 0
}

# -----------------------------------------------------------------------
# check_probe_outputs <probe_dir> <indent>
# Runs all 8 checks; appends per-check reasons to _PROBE_FAIL_LOG
# Returns number of failures
# -----------------------------------------------------------------------
check_probe_outputs() {
    local probe_dir="$1" indent="$2"
    local failures=0
    local ks4="$probe_dir/kilosort4"
    local so="$ks4/sorter_output"
    local ibl="$ks4/ibl_format"

    check_dir_nonempty  "kilosort4"                            "$ks4"               "$indent" || (( failures++ ))
    check_dir_nonempty  "kilosort4/sorter_output"              "$so"                "$indent" || (( failures++ ))
    check_dir_nonempty  "kilosort4/ibl_format"                 "$ibl"               "$indent" || (( failures++ ))
    check_file_exists   "sorter_output/cluster_info.tsv"       "$so" "cluster_info.tsv" "$indent" || (( failures++ ))
    check_dir_nonempty  "sorter_output/bombcell"               "$so/bombcell"       "$indent" || (( failures++ ))
    check_dir_nonempty  "sorter_output/bombcell_plots"         "$so/bombcell_plots" "$indent" || (( failures++ ))
    check_dir_nonempty  "sorter_output/cwaves"                 "$so/cwaves"         "$indent" || (( failures++ ))
    check_glob_nonempty "sorter_output/cwaves/waveform_metrics*" "$so/cwaves" "waveform_metrics*" "$indent" || (( failures++ ))

    return $failures
}

# -----------------------------------------------------------------------
# Collect valid sessions (same logic as detect_ephys_sessions.sh)
# -----------------------------------------------------------------------
if [ -z "$SUFFIX" ]; then
    search_pattern="*.ap.bin"
else
    search_pattern="*${SUFFIX}.imec*.ap.bin"
fi

ephys_dirs=$(find "$DATA_DIR" -type d -name "Ephys" 2>/dev/null)

if [ -z "$ephys_dirs" ]; then
    echo "ERROR: No 'Ephys' directories found in $DATA_DIR"
    exit 1
fi

all_sessions=()
while IFS= read -r ephys_dir; do
    ap_files=$(find "$ephys_dir" -name "$search_pattern" -type f 2>/dev/null)
    [ -z "$ap_files" ] && continue
    rel_path="${ephys_dir#$DATA_DIR/}"
    mouse=$(echo "$rel_path" | cut -d'/' -f1)
    [[ "$mouse" == AB* || "$mouse" == MH* ]] || continue
    all_sessions+=("$ephys_dir")
done <<< "$ephys_dirs"

IFS=$'\n' all_sessions=($(sort <<<"${all_sessions[*]}")); unset IFS
total_sessions=${#all_sessions[@]}

if [ "$total_sessions" -eq 0 ]; then
    echo "ERROR: No valid sessions found in $DATA_DIR"
    exit 1
fi

echo "Found $total_sessions valid session(s) to check"
echo ""

# -----------------------------------------------------------------------
# Per-session / per-probe checks
# -----------------------------------------------------------------------
failed_sessions=()        # Ephys paths of failed sessions
session_reasons=()        # Human-readable failure reasons (parallel array)

for ephys_dir in "${all_sessions[@]}"; do
    rel_path="${ephys_dir#$DATA_DIR/}"
    mouse=$(echo "$rel_path" | cut -d'/' -f1)
    session=$(echo "$rel_path" | cut -d'/' -f2)

    echo "========================================"
    echo "Session : $mouse / $session"
    echo "Path    : $ephys_dir"
    echo "========================================"

    catgt_dirs=$(find "$ephys_dir" -maxdepth 1 -type d -name "catgt_*" 2>/dev/null | sort)

    if [ -z "$catgt_dirs" ]; then
        echo "  ✗ MISSING: no catgt_* directory found under Ephys/"
        echo "  ★ Session FAILED"
        echo ""
        failed_sessions+=("$ephys_dir")
        session_reasons+=("no catgt_* directory")
        continue
    fi

    session_total_failures=0
    session_probes_checked=0
    session_probes_failed=0
    session_reason_parts=()   # one entry per failing probe

    while IFS= read -r catgt_dir; do
        catgt_name=$(basename "$catgt_dir")
        echo ""
        echo "  catgt: $catgt_name"

        probe_dirs=$(find "$catgt_dir" -maxdepth 1 -type d -name "*imec*" 2>/dev/null | sort)

        while IFS= read -r probe_dir; do
            [ -z "$probe_dir" ] && continue
            probe_name=$(basename "$probe_dir")
            echo ""
            echo "    Probe: $probe_name"
            (( session_probes_checked++ ))

            # Clear the per-probe failure log before each probe
            > "$_PROBE_FAIL_LOG"

            check_probe_outputs "$probe_dir" "    "
            probe_failures=$?

            if [ "$probe_failures" -gt 0 ]; then
                echo "    → $probe_failures check(s) FAILED"
                (( session_probes_failed++ ))
                (( session_total_failures += probe_failures ))

                # Collect reasons for this probe: "probeName: reason1; reason2"
                probe_reasons=$(paste -sd '; ' "$_PROBE_FAIL_LOG")
                session_reason_parts+=("${probe_name}: ${probe_reasons}")
            else
                echo "    → All checks passed"
            fi
        done <<< "$probe_dirs"

    done <<< "$catgt_dirs"

    echo ""
    echo "  Session summary: $session_probes_checked probe(s) checked, $session_probes_failed failed"

    if [ "$session_total_failures" -gt 0 ]; then
        echo "  ★ Session FAILED ($session_total_failures total check failure(s))"
        failed_sessions+=("$ephys_dir")
        # Join all probe reason strings with " | "
        joined_reasons=$(printf '%s | ' "${session_reason_parts[@]}")
        joined_reasons="${joined_reasons% | }"   # trim trailing separator
        session_reasons+=("$joined_reasons")
    else
        echo "  ★ Session OK"
    fi
    echo ""
done

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
num_failed=${#failed_sessions[@]}
num_ok=$(( total_sessions - num_failed ))

echo "========================================"
echo "SUMMARY"
echo "========================================"
echo "Total sessions checked : $total_sessions"
echo "Passed                 : $num_ok"
echo "Failed                 : $num_failed"
echo ""

if [ "$num_failed" -eq 0 ]; then
    echo "All sessions passed — no preprocessing failures detected."
    exit 0
fi

echo "========================================"
echo "FAILED SESSIONS"
echo "========================================"
for i in "${!failed_sessions[@]}"; do
    printf "%3d. %s\n     %s\n" $((i+1)) "${failed_sessions[$i]}" "${session_reasons[$i]}"
done
echo ""

# -----------------------------------------------------------------------
# Write output file: one session per line, reason appended after a tab
# Format: <ephys_path>\t<reason>
# -----------------------------------------------------------------------
echo "========================================"
echo "Writing failed sessions to: $OUTPUT_FILE"
echo "========================================"

mkdir -p "$(dirname "$OUTPUT_FILE")"
> "$OUTPUT_FILE"
for i in "${!failed_sessions[@]}"; do
    printf "%s\t%s\n" "${failed_sessions[$i]}" "${session_reasons[$i]}" >> "$OUTPUT_FILE"
done

echo "Written $num_failed failed session(s) to:"
echo "  $OUTPUT_FILE"
echo ""
echo "To re-run spike sorting only on failed sessions:"
echo "  # Extract paths only (column 1):"
echo "  cut -f1 $OUTPUT_FILE > /tmp/rerun_sessions.txt"
echo "  sbatch --array=1-${num_failed}%4 slurm_spikesort_array_jed.sh"
