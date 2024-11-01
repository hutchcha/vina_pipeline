#!/bin/bash

# ======================================
# Serial Batch Virtual Screening Script with Aggregation
# ======================================

# Exit immediately if a command exits with a non-zero status
set -e

# ----------------------------
# Default Configuration Parameters
# ----------------------------

# Path to AutoDock Vina GPU executable
VINA_GPU_EXEC="/mnt/data_SSD/VINA-GPU/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1"

# Log file for recording progress and errors
LOG_FILE="serial_vina_gpu_with_aggregation.log"

# ----------------------------
# Command-Line Arguments Parsing
# ----------------------------

show_help() {
  echo "Usage: $0 -c CONFIG_FILE -l LIGAND_DIR -o OUTPUT_DIR -a AGGREGATE_OUTPUT_DIR [-s START_FROM_DIR] [--single-directory]"
  echo ""
  echo "Required Arguments:"
  echo "  -c CONFIG_FILE             Path to AutoDock Vina configuration file."
  echo "  -l LIGAND_DIR              Directory containing ligand directories (e.g., part_1, part_2, ...) or a single directory."
  echo "  -o OUTPUT_DIR              Directory for output directories."
  echo "  -a AGGREGATE_OUTPUT_DIR    Directory to aggregate all output files."
  echo ""
  echo "Optional Arguments:"
  echo "  -s START_FROM_DIR          Start processing from this directory number (default: process all)."
  echo "  --single-directory         Treat LIGAND_DIR as a single directory rather than multiple subdirectories."
}

# Initialize variables
START_FROM_DIR=1
SINGLE_DIRECTORY_MODE=false

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c) CONFIG_FILE="$2"; shift ;;
        -l) LIGAND_PARENT_DIR="$2"; shift ;;
        -o) OUTPUT_PARENT_DIR="$2"; shift ;;
        -a) AGGREGATE_OUTPUT_DIR="$2"; shift ;;
        -s) START_FROM_DIR="$2"; shift ;;
        --single-directory) SINGLE_DIRECTORY_MODE=true ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown option: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Check for required arguments
if [ -z "$CONFIG_FILE" ] || [ -z "$LIGAND_PARENT_DIR" ] || [ -z "$OUTPUT_PARENT_DIR" ] || [ -z "$AGGREGATE_OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# ----------------------------
# Create Output and Aggregate Directories
# ----------------------------

mkdir -p "$OUTPUT_PARENT_DIR"
mkdir -p "$AGGREGATE_OUTPUT_DIR"

# ----------------------------
# Initialize Log File
# ----------------------------

echo "Serial Batch Virtual Screening with Aggregation Started at $(date)" > "$LOG_FILE"

# ----------------------------
# Determine Ligand Directories to Process
# ----------------------------

if [ "$SINGLE_DIRECTORY_MODE" = true ]; then
    ligand_dirs=("$LIGAND_PARENT_DIR")  # Treat as a single directory
else
    # Find all subdirectories named "part_*"
    ligand_dirs=($(find "$LIGAND_PARENT_DIR" -maxdepth 1 -type d -name "part_*" | sort))
fi

# Check if any directories were found
total_dirs=${#ligand_dirs[@]}
if [ "$total_dirs" -eq 0 ]; then
    echo "No ligand directories found in '$LIGAND_PARENT_DIR'." | tee -a "$LOG_FILE"
    exit 1
fi

echo "Found $total_dirs ligand directories to process." | tee -a "$LOG_FILE"

# ----------------------------
# Loop Through Each Ligand Directory
# ----------------------------

current_dir=1
start_time=$(date +%s)

for ligand_dir in "${ligand_dirs[@]}"; do
    dir_name=$(basename "$ligand_dir")
    output_dir="${OUTPUT_PARENT_DIR}/${OUTPUT_DIR_PREFIX}_${dir_name}"

    if [ "$SINGLE_DIRECTORY_MODE" = false ]; then
        # Extract the numeric part of the directory name for skipping functionality
        dir_num=${dir_name#part_}
        dir_num=$(echo "$dir_num" | sed 's/^0*//')
        if ! [[ "$dir_num" =~ ^[0-9]+$ ]] || (( dir_num < START_FROM_DIR )); then
            echo "Skipping directory $dir_name (number $dir_num)." | tee -a "$LOG_FILE"
            ((current_dir++))
            continue
        fi
    fi

    echo "Processing directory $current_dir of $total_dirs: $dir_name" | tee -a "$LOG_FILE"
    mkdir -p "$output_dir"

    # Run AutoDock Vina GPU
    "$VINA_GPU_EXEC" --config "$CONFIG_FILE" --ligand_directory "$ligand_dir" --output_directory "$output_dir" >> "$LOG_FILE" 2>&1

    echo "----------------------------------------" | tee -a "$LOG_FILE"
    ((current_dir++))
done

# ----------------------------
# Aggregation Step: Move All Output Files to Aggregate Directory
# ----------------------------

echo "Starting aggregation into '$AGGREGATE_OUTPUT_DIR'..." | tee -a "$LOG_FILE"
output_dirs=($(find "$OUTPUT_PARENT_DIR" -maxdepth 1 -type d -name "${OUTPUT_DIR_PREFIX}_part_*" | sort))

for out_dir in "${output_dirs[@]}"; do
    mv "$out_dir"/* "$AGGREGATE_OUTPUT_DIR"/ || echo "Failed to move files from $out_dir." | tee -a "$LOG_FILE"
    rmdir "$out_dir" 2>/dev/null || true
done

end_time=$(date +%s)
total_elapsed=$(( end_time - start_time ))
printf "Total Time Taken: %02d:%02d:%02d (HH:MM:SS)\n" $((total_elapsed/3600)) $(((total_elapsed%3600)/60)) $((total_elapsed%60)) | tee -a "$LOG_FILE"
echo "All results are aggregated in '$AGGREGATE_OUTPUT_DIR/'." | tee -a "$LOG_FILE"
