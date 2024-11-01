#!/bin/bash

# ======================================
# Serial Batch Virtual Screening Script with Aggregation
# ======================================

# Exit immediately if a command exits with a non-zero status
set -e

# ----------------------------
# Configuration Parameters
# ----------------------------

# Path to AutoDock Vina GPU executable
VINA_GPU_EXEC="/mnt/data_SSD/VINA-GPU/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1"

# Path to configuration file
CONFIG_FILE="/mnt/data/virtualscreen/rheb/os1/p2/conf.txt"

# Parent directory containing ligand directories (e.g., part_1, part_2, ...)
LIGAND_PARENT_DIR="/mnt/data/virtualscreen/rheb/os1/molecule_sets/lead_like_sample"

# Parent directory for output directories
OUTPUT_PARENT_DIR="/mnt/data/virtualscreen/rheb/os1/p2/large_set_out"

# Prefix for output directories
OUTPUT_DIR_PREFIX="out_test"

# Directory to aggregate all output files
AGGREGATE_OUTPUT_DIR="/mnt/data/virtualscreen/rheb/os1/p2/aggregate_results"

# Log file for recording progress and errors
LOG_FILE="serial_vina_gpu_with_aggregation.log"

# ----------------------------
# Create Output and Aggregate Directories
# ----------------------------

mkdir -p "$OUTPUT_PARENT_DIR"
mkdir -p "$AGGREGATE_OUTPUT_DIR"
START_FROM_DIR=216
# ----------------------------
# Initialize Log File
# ----------------------------

echo "Serial Batch Virtual Screening with Aggregation Started at $(date)" > "$LOG_FILE"

# ----------------------------
# Find All Ligand Directories
# ----------------------------

# Assuming ligand directories are named "part_1", "part_2", etc.
ligand_dirs=($(find "$LIGAND_PARENT_DIR" -maxdepth 1 -type d -name "part_*" | sort))

# Get the total number of ligand directories
total_dirs=${#ligand_dirs[@]}
if [ "$total_dirs" -eq 0 ]; then
    echo "No ligand directories found in '$LIGAND_PARENT_DIR'." | tee -a "$LOG_FILE"
    exit 1
fi

echo "Found $total_dirs ligand directories to process." | tee -a "$LOG_FILE"

# ----------------------------
# Initialize Counters and Timer
# ----------------------------

current_dir=1
start_time=$(date +%s)

# ----------------------------
# Loop Through Each Ligand Directory
# ----------------------------

for ligand_dir in "${ligand_dirs[@]}"; do
    dir_name=$(basename "$ligand_dir")
    output_dir="${OUTPUT_PARENT_DIR}/${OUTPUT_DIR_PREFIX}_${dir_name}"
  # Extract the numeric part of the directory name
    dir_num=${dir_name#part_}
    dir_num=$(echo "$dir_num" | sed 's/^0*//')  # Remove leading zeros if any

    # Check if dir_num is a valid number
    if ! [[ "$dir_num" =~ ^[0-9]+$ ]]; then
        echo "Directory name $dir_name does not contain a valid number. Skipping." | tee -a "$LOG_FILE"
        ((current_dir++))
        continue
    fi


    # Skip directories before the starting directory
    if (( dir_num < START_FROM_DIR )); then
        echo "Skipping directory $dir_name (number $dir_num) as it's before START_FROM_DIR $START_FROM_DIR." | tee -a "$LOG_FILE"
        ((current_dir++))
        continue
    fi

    echo "Processing directory $current_dir of $total_dirs: $dir_name" | tee -a "$LOG_FILE"
    echo "Output will be saved to: $output_dir" | tee -a "$LOG_FILE"

    # Create the output directory
    mkdir -p "$output_dir"

    # Run AutoDock Vina GPU
    "$VINA_GPU_EXEC" \
        --config "$CONFIG_FILE" \
        --ligand_directory "$ligand_dir" \
        --output_directory "$output_dir" >> "$LOG_FILE" 2>&1

    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed $dir_name." | tee -a "$LOG_FILE"
    else
        echo "ERROR: Failed to process $dir_name. Check log for details." | tee -a "$LOG_FILE"
    fi

    # Increment the directory counter
    ((current_dir++))

    # Optional: Display elapsed time
    current_time=$(date +%s)
    elapsed=$(( current_time - start_time ))
    hours=$((elapsed / 3600))
    minutes=$(( (elapsed % 3600) / 60 ))
    seconds=$((elapsed % 60))
    echo "Elapsed Time: $hours:$minutes:$seconds (HH:MM:SS)" | tee -a "$LOG_FILE"

    echo "----------------------------------------" | tee -a "$LOG_FILE"
done

# ----------------------------
# Aggregation Step: Move All Output Files to Aggregate Directory
# ----------------------------

echo "Starting aggregation of all output files into '$AGGREGATE_OUTPUT_DIR'..." | tee -a "$LOG_FILE"

# Find all output directories
output_dirs=($(find "$OUTPUT_PARENT_DIR" -maxdepth 1 -type d -name "${OUTPUT_DIR_PREFIX}_part_*" | sort))

# Initialize aggregation counter
total_output_dirs=${#output_dirs[@]}
current_output_dir=1

for out_dir in "${output_dirs[@]}"; do
    echo "Aggregating directory $current_output_dir of $total_output_dirs: $(basename "$out_dir")" | tee -a "$LOG_FILE"
    
    # Move all files from the current output directory to the aggregate directory
    # Use 'mv' with wildcard to move all files; adjust patterns if needed
    if mv "$out_dir"/* "$AGGREGATE_OUTPUT_DIR"/; then
        echo "Successfully moved files from $(basename "$out_dir") to aggregate directory." | tee -a "$LOG_FILE"
    else
        echo "ERROR: Failed to move files from $(basename "$out_dir")." | tee -a "$LOG_FILE"
    fi
    
    # Optionally, remove the now-empty output directory
    rmdir "$out_dir" && echo "Removed empty directory: $(basename "$out_dir")" | tee -a "$LOG_FILE" || echo "Directory not empty or failed to remove: $(basename "$out_dir")" | tee -a "$LOG_FILE"
    
    # Increment the aggregation counter
    ((current_output_dir++))
    
    echo "----------------------------------------" | tee -a "$LOG_FILE"
done

# ----------------------------
# Final Summary
# ----------------------------

end_time=$(date +%s)
total_elapsed=$(( end_time - start_time ))
hours=$((total_elapsed / 3600))
minutes=$(( (total_elapsed % 3600) / 60 ))
seconds=$((total_elapsed % 60))

echo "Serial Batch Virtual Screening with Aggregation Completed at $(date)." | tee -a "$LOG_FILE"
printf "Total Time Taken: %02d:%02d:%02d (HH:MM:SS)\n" $hours $minutes $seconds | tee -a "$LOG_FILE"
echo "All results are aggregated in '$AGGREGATE_OUTPUT_DIR/'." | tee -a "$LOG_FILE"
echo "Detailed logs are available in '$LOG_FILE'." | tee -a "$LOG_FILE"
