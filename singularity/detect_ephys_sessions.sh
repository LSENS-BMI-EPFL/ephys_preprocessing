#!/bin/bash

# Script to detect all Ephys sessions with ap.bin files
# Searches for mouse_name/session_name/Ephys folders containing .ap.bin files

DATA_DIR="/scratch/bisi/data"
OUTPUT_FILE="/home/bisi/code/ephys_preprocessing/config/inputs_axel_test.txt"
SUFFIX="corrected"  # Set to "" to find any .ap.bin, or "corrected" for files like *corrected.imec*.ap.bin

echo "========================================"
echo "Scanning for Ephys sessions in: $DATA_DIR"
echo "Looking for suffix: ${SUFFIX:-'any'}"
echo "========================================"
echo ""

# Find all sessions with ap.bin files
sessions=()

# Search pattern based on suffix
if [ -z "$SUFFIX" ]; then
    search_pattern="*.ap.bin"
else
    search_pattern="*${SUFFIX}.imec*.ap.bin"
fi

echo "Search pattern: $search_pattern"
echo "Searching..."
echo ""

# Find all Ephys directories
ephys_dirs=$(find "$DATA_DIR" -type d -name "Ephys" 2>/dev/null)

if [ -z "$ephys_dirs" ]; then
    echo "ERROR: No 'Ephys' directories found in $DATA_DIR"
    exit 1
fi

echo "Found $(echo "$ephys_dirs" | wc -l) Ephys directories"
echo ""

# Process each Ephys directory
while IFS= read -r ephys_dir; do
    echo "Checking: $ephys_dir"

    # Check if this Ephys dir contains any .ap.bin files (NO MAXDEPTH LIMIT)
    ap_files=$(find "$ephys_dir" -name "$search_pattern" -type f 2>/dev/null)

    if [ -n "$ap_files" ]; then
        num_files=$(echo "$ap_files" | wc -l)
        echo "  ✓ Found $num_files matching files"
        echo "$ap_files" | head -3 | sed 's/^/    /'
        # Store FULL absolute path
        sessions+=("$ephys_dir")
    else
        echo "  ✗ No matching files"
    fi
done <<< "$ephys_dirs"

echo ""

# Sort sessions
IFS=$'\n' sessions=($(sort <<<"${sessions[*]}"))
unset IFS

# Count unique mice (first directory level after DATA_DIR)
declare -A mice
for session in "${sessions[@]}"; do
    rel_path="${session#$DATA_DIR/}"
    mouse=$(echo "$rel_path" | cut -d'/' -f1)
    mice[$mouse]=1
done

num_mice=${#mice[@]}
num_sessions=${#sessions[@]}

echo "========================================"
echo "SUMMARY"
echo "========================================"
echo "Total mice: $num_mice"
echo "Total sessions: $num_sessions"
echo ""

if [ $num_sessions -eq 0 ]; then
    echo "ERROR: No sessions found with pattern: $search_pattern"
    echo ""
    echo "Debugging: Let's check what .ap.bin files exist..."
    echo "First 20 .ap.bin files in $DATA_DIR:"
    find "$DATA_DIR" -name "*.ap.bin" -type f 2>/dev/null | head -20
    exit 1
fi

echo "========================================"
echo "MICE FOUND"
echo "========================================"
for mouse in "${!mice[@]}"; do
    echo "  - $mouse"
done | sort
echo ""

echo "========================================"
echo "SESSIONS FOUND"
echo "========================================"
for i in "${!sessions[@]}"; do
    session="${sessions[$i]}"
    # Count .ap.bin files in this session
    num_probes=$(find "$session" -name "$search_pattern" -type f 2>/dev/null | wc -l)
    printf "%3d. %s (%d probes)\n" $((i+1)) "$session" "$num_probes"
done
echo ""

# Write to inputs_axel.txt
echo "========================================"
echo "Writing to: $OUTPUT_FILE"
echo "========================================"

# Create directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Write sessions to file (FULL PATHS)
printf "%s\n" "${sessions[@]}" > "$OUTPUT_FILE"

echo "Successfully written $num_sessions sessions to:"
echo "  $OUTPUT_FILE"
echo ""
echo "Content preview:"
head -5 "$OUTPUT_FILE" | sed 's/^/  /'
echo ""
echo "You can now run your SLURM array job with:"
echo "  sbatch --array=1-${num_sessions}%4 slurm_spikesort_array_jed.sh"