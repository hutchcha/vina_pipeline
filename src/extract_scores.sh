#!/usr/bin/env bash

# ============================================================
# Script: extract_model1_scores.sh
# Description:
#   Extracts docking scores from MODEL 1 of each *_out.pdbqt file.
#   Outputs a text file with two columns: ZINCID and DockingScore.
#
# Usage:
#   ./extract_model1_scores.sh -d /path/to/pdbqt_dir -o scores_model1.txt -j 8
#
# Options:
#   -d, --directory     Directory containing *_out.pdbqt files.
#   -o, --output        Output file to store ZINCID and DockingScore.
#   -j, --jobs          Number of parallel jobs (default: number of CPU cores).
#   -h, --help          Display help message.
#
# =============================================================

# Function to display usage
usage() {
    echo "Usage: $0 -d /path/to/pdbqt_dir -o scores_model1.txt -j 8"
    echo ""
    echo "Options:"
    echo "  -d, --directory     Directory containing *_out.pdbqt files."
    echo "  -o, --output        Output file to store ZINCID and DockingScore."
    echo "  -j, --jobs          Number of parallel jobs (default: number of CPU cores)."
    echo "  -h, --help          Display this help message."
    exit 1
}

# Default values
PARALLEL_JOBS=$(nproc)  # Number of CPU cores
OUTPUT_FILE="scores_model1.txt"

# Parse arguments
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -d|--directory)
        DIRECTORY="$2"
        shift # past argument
        shift # past value
        ;;
        -o|--output)
        OUTPUT_FILE="$2"
        shift
        shift
        ;;
        -j|--jobs)
        PARALLEL_JOBS="$2"
        shift
        shift
        ;;
        -h|--help)
        usage
        ;;
        *)    # unknown option
        echo "Unknown option: $1"
        usage
        ;;
    esac
done

# Check if DIRECTORY is provided
if [ -z "$DIRECTORY" ]; then
    echo "Error: Directory not specified."
    usage
fi

# Check if DIRECTORY exists
if [ ! -d "$DIRECTORY" ]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# Check if pv is installed for progress bar
if command -v pv >/dev/null 2>&1; then
    USE_PV=1
else
    USE_PV=0
    echo "Warning: 'pv' not found. Progress bar will not be displayed."
fi

# Inform user of the settings
echo "Directory to search: $DIRECTORY"
echo "Output file: $OUTPUT_FILE"
echo "Parallel jobs: $PARALLEL_JOBS"

# Temporary file to store grep output
TMP_GREP_OUTPUT=$(mktemp)

# Function to clean up temporary files on exit
cleanup() {
    rm -f "$TMP_GREP_OUTPUT"
}
trap cleanup EXIT

# Start extraction
echo "Starting extraction of docking scores from MODEL 1..."

# Step 1: Find all *_out.pdbqt files
# Step 2: Use xargs to run grep in parallel, extracting the first REMARK VINA RESULT per file
# Step 3: Use awk to parse and format the output

# Calculate total number of files for progress bar
TOTAL_FILES=$(find "$DIRECTORY" -type f -name "*_out.pdbqt" | wc -l)

echo "Total files to process: $TOTAL_FILES"

if [ "$USE_PV" -eq 1 ]; then
    # Use pv for progress bar
    find "$DIRECTORY" -type f -name "*_out.pdbqt" -print0 \
    | pv -0 -s "$TOTAL_FILES" -l -c -N "Processing files" \
    | xargs -0 -P "$PARALLEL_JOBS" -I{} sh -c '
        # For each file "{}", grep "REMARK VINA RESULT" and take the first occurrence
        grep -H "REMARK VINA RESULT" "{}" | head -n 1
    ' > "$TMP_GREP_OUTPUT"
else
    # Without pv, just run the commands
    find "$DIRECTORY" -type f -name "*_out.pdbqt" -print0 \
    | xargs -0 -P "$PARALLEL_JOBS" -I{} sh -c '
        grep -H "REMARK VINA RESULT" "{}" | head -n 1
    ' > "$TMP_GREP_OUTPUT"
fi

# Step 4: Parse the grep output with awk
echo "Parsing extracted lines..."

awk -F ':' '
{
    # Example line:
    # /path/to/dir/ZINCui00000000EZ_out.pdbqt:REMARK VINA RESULT:     -10.2      0.000      0.000

    if (NF < 3) {
        next  # skip malformed lines
    }

    filename = $1
    remark = $2
    result = $3

    # Extract docking score
    gsub(/^ +/, "", result)  # Remove leading spaces
    split(result, arr, " +")
    score = arr[1]

    # Extract ZINC ID from filename
    n = split(filename, pathArr, "/")
    fname = pathArr[n]  # e.g., ZINCui00000000EZ_out.pdbqt
    sub("_out.pdbqt$", "", fname)  # Remove suffix

    # Print ZINCID and docking score
    print fname, score
}
' "$TMP_GREP_OUTPUT" > "$OUTPUT_FILE"

echo "Extraction and parsing complete."

# Final output
echo "Scores saved to '$OUTPUT_FILE'."

# Optionally, display a few lines of the output
echo "First 5 lines of '$OUTPUT_FILE':"
head -n 5 "$OUTPUT_FILE"
=
