#!/bin/bash

# ============================================
# Master Script for Virtual Screening Pipeline
# ============================================

# Exit on errors
set -e

# ----------------------------
# Helper Functions
# ----------------------------
function print_usage() {
    echo "Usage: $0 -c CONFIG_FILE -i INPUT_CSV -o OUTPUT_DIR [-s START_FROM_DIR] [--single-dir] [-v]"
    echo ""
    echo "Required Arguments:"
    echo "  -c CONFIG_FILE       Path to the VINA-GPU configuration file."
    echo "  -i INPUT_CSV         Path to the input CSV file containing SMILES and IDs."
    echo "  -o OUTPUT_DIR        Path to the main output directory."
    echo ""
    echo "Optional Arguments:"
    echo "  -s START_FROM_DIR    Start processing ligands from this directory index (default: 1)."
    echo "  --single-dir         Use a single directory instead of splitting ligands into batches."
    echo "  -v                   Enable verbose mode."
    echo ""
    echo "Example:"
    echo "  $0 -c config.txt -i molecules.csv -o ./screen_output -s 2 --single-dir -v"
}

# ----------------------------
# Parse Arguments
# ----------------------------
START_FROM_DIR=1
SINGLE_DIR_MODE=false
VERBOSE=false

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c) CONFIG_FILE="$2"; shift ;;
        -i) INPUT_CSV="$2"; shift ;;
        -o) OUTPUT_DIR="$2"; shift ;;
        -s) START_FROM_DIR="$2"; shift ;;
        --single-dir) SINGLE_DIR_MODE=true ;;
        -v) VERBOSE=true ;;
        -h|--help) print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; print_usage; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [ -z "$CONFIG_FILE" ] || [ -z "$INPUT_CSV" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments."
    print_usage
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Log file for this pipeline
LOG_FILE="$OUTPUT_DIR/pipeline.log"
echo "Virtual Screening Pipeline Log - $(date)" > "$LOG_FILE"

# ----------------------------
# Step 1: Convert SMILES to PDBQT
# ----------------------------
echo "Step 1: Converting SMILES to PDBQT..." | tee -a "$LOG_FILE"
PDBQT_OUTPUT_DIR="$OUTPUT_DIR/ligand_pdbqts"
mkdir -p "$PDBQT_OUTPUT_DIR"

SMILES_TO_PDBQT_CMD="python3 smiles_to_pdbqt.py \
    $INPUT_CSV \
    $PDBQT_OUTPUT_DIR \
    --molecules_per_dir 1000"

# Add optional flags
if $SINGLE_DIR_MODE; then
    SMILES_TO_PDBQT_CMD+=" --single_dir"
fi
if $VERBOSE; then
    SMILES_TO_PDBQT_CMD+=" -v"
fi

# Run the SMILES-to-PDBQT script
echo "Running: $SMILES_TO_PDBQT_CMD" | tee -a "$LOG_FILE"
eval $SMILES_TO_PDBQT_CMD

# ----------------------------
# Step 2: Run VINA-GPU Screening
# ----------------------------
echo "Step 2: Running VINA-GPU screening..." | tee -a "$LOG_FILE"
VINA_OUTPUT_DIR="$OUTPUT_DIR/vina_results"
AGGREGATE_DIR="$OUTPUT_DIR/aggregated_results"
mkdir -p "$VINA_OUTPUT_DIR"
mkdir -p "$AGGREGATE_DIR"

RUN_SCREEN_CMD="./run_large_screen.sh \
    -c $CONFIG_FILE \
    -l $PDBQT_OUTPUT_DIR \
    -o $VINA_OUTPUT_DIR \
    -a $AGGREGATE_DIR \
    -s $START_FROM_DIR"

# Add optional flags
if $SINGLE_DIR_MODE; then
    RUN_SCREEN_CMD+=" --single-directory"
fi

# Run the VINA-GPU screening script
echo "Running: $RUN_SCREEN_CMD" | tee -a "$LOG_FILE"
eval $RUN_SCREEN_CMD

# ----------------------------
# Step 3: Summarize Virtual Screening Results
# ----------------------------
echo "Step 3: Summarizing virtual screening results..." | tee -a "$LOG_FILE"
SUMMARY_CSV="$OUTPUT_DIR/screen_summary.csv"

SUMMARIZE_CMD="python3 summarize_virtual_screen.py \
    -c $INPUT_CSV \
    -d $AGGREGATE_DIR \
    -g -n 15 \
    -i $OUTPUT_DIR/images \
    -f $OUTPUT_DIR/top_hits_figure.png"

# Add optional verbosity
if $VERBOSE; then
    SUMMARIZE_CMD+=" -v"
fi

# Run the summarization script
echo "Running: $SUMMARIZE_CMD" | tee -a "$LOG_FILE"
eval $SUMMARIZE_CMD

# ----------------------------
# Pipeline Complete
# ----------------------------
echo "Pipeline completed successfully!" | tee -a "$LOG_FILE"
echo "Results are in: $OUTPUT_DIR" | tee -a "$LOG_FILE"
