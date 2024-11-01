# Virtual Screening Pipeline Using AutoDock Vina

## Overview

This repository contains scripts for a complete virtual screening pipeline using AutoDock Vina GPU. The pipeline automates ligand preparation, docking, and result summarization, facilitating high-throughput screening of compounds against target proteins.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Setup](#setup)
  - [Updating Script Paths](#updating-script-paths)
  - [Preparing the Configuration File](#preparing-the-configuration-file)
- [Usage](#usage)
  - [Master Script: `run_pipeline.sh`](#master-script-run_pipelinesh)
  - [Scripts](#scripts)
    - [`smiles_to_pdbqt.py`](#smiles_to_pdbqtpy)
    - [`run_large_screen.sh`](#run_large_screensh)
    - [`summarize_virtual_screen.py`](#summarize_virtual_screenpy)
- [Example Command](#example-command)
- [Notes](#notes)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Prerequisites

Ensure you have the following installed on your system:

- **Python 3.x**
- **AutoDock Vina GPU**: Download and install AutoDock Vina with GPU support.
- **MGLTools**: For ligand preparation (`prepare_ligand4.py` script).
- **RDKit**: For molecular handling and conversions.
- **Pillow (PIL)**: For image processing.
- **Matplotlib**: For plotting images.
- **tqdm**: For progress bars in the scripts.

## Installation

1. **Clone the Repository**

   ```bash
   git clone https://github.com/your_username/virtual_screening_pipeline.git
   cd virtual_screening_pipeline
```
 
1. **Install Python Dependencies** Use the provided `requirements.txt` file to install necessary Python packages:

```bash
pip install -r requirements.txt
```
`requirements.txt`:

```plaintext
rdkit
pillow
matplotlib
tqdm
```

## Setup 

### Updating Script Paths 

Before running the pipeline, update the paths in the scripts to match your system:
 
- **`smiles_to_pdbqt.py`** :Locate and update the paths to `pythonsh` and `prepare_ligand4.py`:

```python
command = [
    "/path/to/mgltools/bin/pythonsh",
    "/path/to/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py",
    '-l', pdb_filename,
    '-o', pdbqt_filename
]
```
 
- **`run_large_screen.sh`** :
Update the path to the AutoDock Vina GPU executable:


```bash
VINA_GPU_EXEC="/path/to/AutoDock-Vina-GPU"
```
 
- **`summarize_virtual_screen.py`** :
Update the font path for image annotations:


```python
try:
    font = ImageFont.truetype("/path/to/fonts/Arial.TTF", size=50)
except IOError:
    print("Arial font not found, using default")
    font = ImageFont.load_default()
```

### Preparing the Configuration File 
Create a configuration file (`conf.txt`) for AutoDock Vina, specifying docking parameters and the search space based on your protein target:

```plaintext
center_x = 10.0
center_y = 25.0
center_z = -5.0

size_x = 20.0
size_y = 20.0
size_z = 20.0
thread = 8000
exhaustiveness = 8
num_modes = 9
```

Modify the values according to your protein's binding site.
**Note** : You must provide your own configuration file and define the search space based on your specific protein pocket and structure file, as automating this step is not straightforward.
## Usage 
Master Script: `master.sh`The `run_pipeline.sh` script orchestrates the entire pipeline.
#### Usage 


```bash
bash master.sh -c CONFIG_FILE -s SMILES_FILE -l LIGAND_PARENT_DIR [OPTIONS]
```

#### Required Arguments 
 
- `-c CONFIG_FILE`: Path to the AutoDock Vina configuration file.
 
- `-s SMILES_FILE`: Path to the SMILES file containing ligand information.
 
- `-l LIGAND_PARENT_DIR`: Directory to store converted ligands (PDBQT files).

#### Optional Arguments 
 
- `-o OUTPUT_DIR`: Directory for docking results (default: `output_directory`).
 
- `-a AGGREGATE_DIR`: Directory to store aggregated results (default: `aggregate_results`).
 
- `-f SUMMARY_FILE`: Output summary file (default: `docking_summary.txt`).
 
- `-i IMAGE_DIR`: Directory for images (default: `images`).
 
- `-g FIGURE_FILE`: Combined figure file for top hits (default: `top_hits_figure.png`).
 
- `-m MAX_FILES_PER_DIR`: Max files per subdirectory for SMILES conversion (default: 1000).
 
- `-p MOLECULES_PER_DIR`: Number of molecules per output subdirectory (default: 1000).
 
- `-w NUM_WORKERS`: Number of worker processes for SMILES conversion (default: number of CPU cores).
 
- `-x MAX_MOLECULES`: Maximum number of molecules to process (default: all).
 
- `-n TOP_N`: Number of top hits to generate images for (default: 15).
 
- `-v`: Enable verbose output.
 
- `-d`: Treat `LIGAND_PARENT_DIR` as a single directory (no subdirectories).
 
- `-s START_FROM_DIR`: Start processing from this directory number (default: 1).
 
- `-h`: Show help message and exit.

### Scripts 
`smiles_to_pdbqt.py`
Converts SMILES strings to PDBQT files, splitting outputs into subdirectories.
**Usage** :

```bash
python smiles_to_pdbqt.py INPUT_FILE OUTPUT_DIR [OPTIONS]
```
**Options** : 
- `-n`, `--num_workers`: Number of worker processes.
 
- `-m`, `--max_molecules`: Maximum number of molecules to process.
 
- `--molecules_per_dir`: Number of molecules per subdirectory.
`run_large_screen.sh`
Runs AutoDock Vina docking on batches of ligands, aggregating results.
**Usage** :

```bash
bash run_large_screen.sh -c CONFIG_FILE -l LIGAND_DIR -o OUTPUT_DIR -a AGGREGATE_OUTPUT_DIR [OPTIONS]
```
**Options** : 
- `-s START_FROM_DIR`: Start processing from this directory number.
 
- `--single-directory`: Treat `LIGAND_DIR` as a single directory.
`summarize_virtual_screen.py`
Summarizes docking results and generates images for top hits.
**Usage** :

```bash
python summarize_virtual_screen.py -d DIRECTORY -s SMILES_FILE [OPTIONS]
```
**Options** : 
- `-o`, `--output`: Path to the output summary file.
 
- `-g`, `--generate_images`: Enable image generation.
 
- `-n`, `--top_n`: Number of top hits to generate images for.
 
- `-i`, `--image_dir`: Directory to save images.
 
- `-f`, `--figure_file`: Filename for the combined figure.
 
- `-v`, `--verbose`: Enable verbose output.

## Example Command 


```bash
bash run_pipeline.sh \
  -c conf.txt \
  -s data/zinc_smiles.txt \
  -l data/ligands \
  -o results/docking_outputs \
  -a results/aggregated_results \
  -f results/docking_summary.csv \
  -i results/images \
  -g top_hits_figure.png \
  -m 500 \
  -p 500 \
  -w 8 \
  -n 10 \
  -v \
  -d
```

This command will:
 
1. **Convert SMILES to PDBQT** : Processes `zinc_smiles.txt`, splitting outputs into subdirectories with 500 molecules each, using 8 worker processes.
 
2. **Run Docking** : Performs virtual screening using AutoDock Vina, treating the ligand directory as a single directory.
 
3. **Summarize Results** : Generates a summary CSV, images for the top 10 hits, and a combined figure.
 
4. **Verbose Output** : Provides detailed output during execution.

## Notes 
 
- **Path Adjustments** : Ensure all paths in the scripts are correctly set to match your environment.
 
- **Provide Your Own Configuration File** : You must supply your own configuration file (`conf.txt`) and define the search space based on your protein pocket and structure file.
 
- **Font Availability** : Update the font path in `summarize_virtual_screen.py` or use a default font if necessary.
 
- **Permissions** : Ensure you have read/write permissions for all specified directories.
 
- **Atomic Symbols Exclusion** : The `smiles_to_pdbqt.py` script excludes molecules containing certain atomic symbols (e.g., Si, P, B, Fe). Modify the `OFF_LIMITS_SYMBOLS` list in the script if needed.

## Troubleshooting 
 
- **Execution Errors** : Verify that all script paths and permissions are correct.
 
- **Module Not Found** : Ensure all Python dependencies are installed.
 
- **Font Errors** : Modify the font path or use the default font in `summarize_virtual_screen.py`.
 
- **Connection of Scripts** : The pipeline is designed to run sequentially through the master script. If you wish to run individual scripts, ensure that the output directories and file names are consistent across scripts.

## License 

This project is licensed under the MIT License.

## Acknowledgments 
 
- **AutoDock Vina**  for molecular docking software.
 
- **RDKit**  for cheminformatics tools.
 
- **MGLTools**  for ligand preparation utilities.
