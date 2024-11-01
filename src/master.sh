#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ----------------------------
# Default Configuration Parameters
# ----------------------------

# Required arguments (to be provided by the user)
CONFIG_FILE=""
SMILES_FILE=""
LIGAND_PARENT_DIR=""

# Optional parameters with default values
OUTPUT_DIR="output_directory"
AGGREGATE_DIR="aggregate_results"
SUMMARY_FILE="docking_summary.txt"
IMAGE_DIR="images"
FIGURE_FILE="top_hits_figure.png"
MAX_FILES_PER_DIR=1000
MOLECULES_PER_DIR=1000
NUM_WORKERS=$(nproc)
MAX_MOLECULES=""
TOP_N=15
VERBOSE=false
SINGLE_DIRECTORY=false
START_FROM_DIR=1

# Help function to display usage information
show_help() {
    echo "Usage: bash run_pipeline.sh -c CONFIG_FILE -s SMILES_FILE -l LIGAND_PARENT_DIR [OPTIONS]"
    echo ""
    echo "Required Arguments:"
    echo "  -c CONFIG_FILE             Path to the AutoDock Vina configuration file."
    echo "  -s SMILES_FILE             Path to the SMILES file with ligand information."
    echo "  -l LIGAND_PARENT_DIR       Directory to store converted ligands (PDBQT files)."
    echo ""
    echo "Optional Arguments:"
    echo "  -o OUTPUT_DIR              Directory for docking results (default: 'output_directory')."
    echo "  -a AGGREGATE_DIR           Directory to store aggregated results (default: 'aggregate_results')."
    echo "  -f SUMMARY_FILE            Output summary file (default: 'docking_summary.txt')."
    echo "  -i IMAGE_DIR               Directory for images (default: 'images')."
    echo "  -g FIGURE_FILE             Combined figure file for top hits (default: 'top_hits_figure.png')."
    echo "  -m MAX_FILES_PER_DIR       Max files per subdirectory for SMILES conversion (default: 1000)."
    echo "  -p MOLECULES_PER_DIR       Number of molecules per output subdirectory (default: 1000)."
    echo "  -w NUM_WORKERS             Number of worker processes for SMILES conversion (default: number of CPU cores)."
    echo "  -x MAX_MOLECULES           Maximum number of molecules to process (default: all)."
    echo "  -n TOP_N                   Number of top hits to generate images for (default: 15)."
    echo "  -v                         Enable verbose output."
    echo "  -d                         Treat LIGAND_PARENT_DIR as a single directory (no subdirectories)."
    echo "  -s START_FROM_DIR          Start processing from this directory number (default: 1)."
    echo "  -h                         Show this help message and exit."
}

# Parse command-line options
while getopts "c:s:l:o:a:f:i:g:m:p:w:x:n:s:vdh" opt; do
    case $opt in
        c) CONFIG_FILE="$OPTARG" ;;
        s) SMILES_FILE="$OPTARG" ;;
        l) LIGAND_PARENT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        a) AGGREGATE_DIR="$OPTARG" ;;
        f) SUMMARY_FILE="$OPTARG" ;;
        i) IMAGE_DIR="$OPTARG" ;;
        g) FIGURE_FILE="$OPTARG" ;;
        m) MAX_FILES_PER_DIR="$OPTARG" ;;
        p) MOLECULES_PER_DIR="$OPTARG" ;;
        w) NUM_WORKERS="$OPTARG" ;;
        x) MAX_MOLECULES="$OPTARG" ;;
        n) TOP_N="$OPTARG" ;;
        v) VERBOSE=true ;;
        d) SINGLE_DIRECTORY=true ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

# Check if required arguments are provided
if [ -z "$CONFIG_FILE" ] || [ -z "$SMILES_FILE" ] || [ -z "$LIGAND_PARENT_DIR" ]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# ----------------------------
# Step 1: Convert SMILES to PDBQT format with directory splitting
# ----------------------------

echo "Converting SMILES to PDBQT with molecules per directory limit..."

SMILES_TO_PDBQT_SCRIPT="smiles_to_pdbqt.py"

SMILES_TO_PDBQT_ARGS=("$SMILES_FILE" "$LIGAND_PARENT_DIR" "-n" "$NUM_WORKERS" "--molecules_per_dir" "$MOLECULES_PER_DIR")
if [ -n "$MAX_MOLECULES" ]; then
    SMILES_TO_PDBQT_ARGS+=("-m" "$MAX_MOLECULES")
fi

python "$SMILES_TO_PDBQT_SCRIPT" "${SMILES_TO_PDBQT_ARGS[@]}"

# ----------------------------
# Step 2: Run docking using AutoDock Vina for each ligand batch
# ----------------------------

echo "Running virtual screening..."

RUN_SCREEN_SCRIPT="run_large_screen.sh"

RUN_SCREEN_ARGS=("-c" "$CONFIG_FILE" "-l" "$LIGAND_PARENT_DIR" "-o" "$OUTPUT_DIR" "-a" "$AGGREGATE_DIR" "-s" "$START_FROM_DIR")
if [ "$SINGLE_DIRECTORY" = true ]; then
    RUN_SCREEN_ARGS+=("--single-directory")
fi

bash "$RUN_SCREEN_SCRIPT" "${RUN_SCREEN_ARGS[@]}"

# ----------------------------
# Step 3: Summarize docking results
# ----------------------------

echo "Summarizing results..."

SUMMARIZE_SCRIPT="summarize_virtual_screen.py"

SUMMARIZE_ARGS=("-d" "$AGGREGATE_DIR" "-s" "$SMILES_FILE" "-o" "$SUMMARY_FILE" "-n" "$TOP_N" "-i" "$IMAGE_DIR" "-f" "$FIGURE_FILE")
if [ "$VERBOSE" = true ]; then
    SUMMARIZE_ARGS+=("-v")
fi

# Enable image generation by default
SUMMARIZE_ARGS+=("-g")

python "$SUMMARIZE_SCRIPT" "${SUMMARIZE_ARGS[@]}"

echo "Pipeline complete. Results saved to:"
echo " - Summary File: $SUMMARY_FILE"
echo " - Images Directory: $IMAGE_DIR"
echo " - Combined Figure File: $FIGURE_FILE"
