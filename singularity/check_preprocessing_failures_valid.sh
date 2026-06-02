#!/bin/bash

# Script to detect Ephys sessions where spike sorting preprocessing failed.
# Only checks probes that have a valid insertion entry in:
#   /home/bisi/mice_info/joint_probe_insertion_info.xlsx
# Columns expected: mouse_id, date (DD.MM.YYYY as text), probe_id, valid
#
# Session date is extracted as the 2nd underscore-field of the session name
# (YYYYMMDD) and converted to DD.MM.YYYY for comparison.

DATA_DIR="/scratch/bisi/data"
OUTPUT_FILE="/home/bisi/code/ephys_preprocessing/config/failed_sessions_axel.txt"
INSERTION_XLSX="/home/bisi/mice_info/joint_probe_insertion_info.xlsx"
SUFFIX="tcat_corrected"

echo "========================================"
echo "Checking preprocessing outputs in: $DATA_DIR"
echo "========================================"
echo ""

# -----------------------------------------------------------------------
# Step 1 — Parse the insertion table with Python into a temp lookup file
# Format: mouse_id|DD.MM.YYYY|probe_id   (one valid insertion per line)
# -----------------------------------------------------------------------
VALID_INSERTIONS=$(mktemp)
_PROBE_FAIL_LOG=$(mktemp)
trap 'rm -f "$VALID_INSERTIONS" "$_PROBE_FAIL_LOG"' EXIT

python3 - "$INSERTION_XLSX" "$VALID_INSERTIONS" << 'PYEOF'
import sys, subprocess
subprocess.check_call([sys.executable, "-m", "pip", "install", "openpyxl", "--quiet", "--user"])
import pandas as pd

xlsx_path = sys.argv[1]
out_path  = sys.argv[2]

df = pd.read_excel(xlsx_path, dtype=str)

# Normalise column names: lowercase + strip whitespace
df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

# Identify the valid column (accept 'valid', 'is_valid', 'validity', ...)
valid_col = next((c for c in df.columns if "valid" in c), None)
if valid_col is None:
    print(f"ERROR: no 'valid' column found. Columns: {list(df.columns)}", file=sys.stderr)
    sys.exit(1)

# Keep only valid rows (handle 1/0, TRUE/FALSE, yes/no, true/false)
truthy = {"1", "true", "yes", "y"}
df = df[df[valid_col].str.strip().str.lower().isin(truthy)]

with open(out_path, "w") as fh:
    for _, row in df.iterrows():
        mouse    = str(row.get("mouse_name", "")).strip()
        date_str = str(row.get("date", "")).strip()          # DD.MM.YYYY
        probe    = str(int(str(row.get("probe_id", "")).strip()))  # bare int: 0, 1, 2

        if not mouse or not date_str or not probe:
            continue

        # Normalise date: accept both DD.MM.YYYY and DD/MM/YYYY
        date_str = date_str.replace("/", ".")

        fh.write(f"{mouse}|{date_str}|{probe}\n")

print(f"Parsed {len(df)} valid insertion(s) from {xlsx_path}")
PYEOF

python_exit=$?
if [ "$python_exit" -ne 0 ]; then
    echo "ERROR: failed to parse insertion table — aborting."
    exit 1
fi

echo "Valid insertions loaded: $(wc -l < "$VALID_INSERTIONS") entries"
echo ""

# -----------------------------------------------------------------------
# Helper — look up whether a (mouse, session_date_YYYYMMDD, probe_dir_name)
# combination is marked valid in the insertion table.
#
# probe_dir_name example: AB116_g0_imec0  → extract token containing "imec"
# date example: 20240724 → convert to 24.07.2024
# -----------------------------------------------------------------------
is_valid_insertion() {
    local mouse="$1"
    local date_ymd="$2"   # YYYYMMDD
    local probe_dir="$3"  # e.g. AB116_g0_imec0

    # Convert date YYYYMMDD → DD.MM.YYYY
    local dd="${date_ymd:6:2}"
    local mm="${date_ymd:4:2}"
    local yyyy="${date_ymd:0:4}"
    local date_dmy="${dd}.${mm}.${yyyy}"

    # Extract imecN token from probe dir name (last field containing "imec")
    local probe_id
    probe_id=$(echo "$probe_dir" | tr '_' '\n' | grep -i "imec" | tail -1 | sed 's/[^0-9]//g')

    if [ -z "$probe_id" ]; then
        echo "    [warn] could not extract probe number from: $probe_dir"
        return 1
    fi

    # Check lookup file
    grep -qF "${mouse}|${date_dmy}|${probe_id}" "$VALID_INSERTIONS"
}

# -----------------------------------------------------------------------
# Check helpers — print on failure only; write short reason to _PROBE_FAIL_LOG
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

check_probe_outputs() {
    local probe_dir="$1" indent="$2"
    local failures=0
    local ks4="$probe_dir/kilosort4"
    local so="$ks4/sorter_output"
    local ibl="$ks4/ibl_format"

    check_dir_nonempty  "kilosort4"                              "$ks4"               "$indent" || (( failures++ ))
    check_dir_nonempty  "kilosort4/sorter_output"                "$so"                "$indent" || (( failures++ ))
    check_dir_nonempty  "kilosort4/ibl_format"                   "$ibl"               "$indent" || (( failures++ ))
    check_file_exists   "sorter_output/cluster_info.tsv"         "$so" "cluster_info.tsv" "$indent" || (( failures++ ))
    check_dir_nonempty  "sorter_output/bombcell"                 "$so/bombcell"       "$indent" || (( failures++ ))
    check_dir_nonempty  "sorter_output/bombcell_plots"           "$so/bombcell_plots" "$indent" || (( failures++ ))
    check_dir_nonempty  "sorter_output/cwaves"                   "$so/cwaves"         "$indent" || (( failures++ ))
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
failed_sessions=()
session_reasons=()

for ephys_dir in "${all_sessions[@]}"; do
    rel_path="${ephys_dir#$DATA_DIR/}"
    mouse=$(echo "$rel_path" | cut -d'/' -f1)
    session=$(echo "$rel_path" | cut -d'/' -f2)

    # Extract date from session name: e.g. AB116_20240724_102941 → 20240724
    session_date=$(echo "$session" | cut -d'_' -f2)

    echo "========================================"
    echo "Session : $mouse / $session  (date: $session_date)"
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
    session_probes_skipped=0
    session_reason_parts=()

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

            # --- Validity gate ---
            if ! is_valid_insertion "$mouse" "$session_date" "$probe_name"; then
                echo "    → Skipped (not a valid insertion)"
                (( session_probes_skipped++ ))
                continue
            fi

            (( session_probes_checked++ ))

            > "$_PROBE_FAIL_LOG"
            check_probe_outputs "$probe_dir" "    "
            probe_failures=$?

            if [ "$probe_failures" -gt 0 ]; then
                echo "    → $probe_failures check(s) FAILED"
                (( session_probes_failed++ ))
                (( session_total_failures += probe_failures ))
                probe_reasons=$(paste -sd "; " "$_PROBE_FAIL_LOG")
                session_reason_parts+=("${probe_name}: ${probe_reasons}")
            else
                echo "    → All checks passed"
            fi
        done <<< "$probe_dirs"

    done <<< "$catgt_dirs"

    echo ""
    echo "  Session summary: $session_probes_checked probe(s) checked, $session_probes_failed failed, $session_probes_skipped skipped (invalid insertion)"

    if [ "$session_total_failures" -gt 0 ]; then
        echo "  ★ Session FAILED ($session_total_failures total check failure(s))"
        failed_sessions+=("$ephys_dir")
        joined=$(printf "%s | " "${session_reason_parts[@]}"); joined="${joined% | }"
        session_reasons+=("$joined")
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
# Write output file: path TAB reasons
# Extract paths only: cut -f1 failed_sessions_axel.txt
# -----------------------------------------------------------------------
mkdir -p "$(dirname "$OUTPUT_FILE")"
> "$OUTPUT_FILE"
for i in "${!failed_sessions[@]}"; do
    printf "%s\t%s\n" "${failed_sessions[$i]}" "${session_reasons[$i]}" >> "$OUTPUT_FILE"
done

echo "========================================"
echo "Writing failed sessions to: $OUTPUT_FILE"
echo "========================================"
echo "Written $num_failed failed session(s) to:"
echo "  $OUTPUT_FILE"
echo ""
echo "To re-run spike sorting only on failed sessions:"
echo "  cut -f1 $OUTPUT_FILE > /tmp/rerun_sessions.txt"
echo "  sbatch --array=1-${num_failed}%4 slurm_spikesort_array_jed.sh"
