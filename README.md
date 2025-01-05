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
 
- **AutoDock Vina GPU** : Download and install AutoDock Vina with GPU support.
 
- **MGLTools** : For ligand preparation (`prepare_ligand4.py` script).
 
- **RDKit** : For molecular handling and conversions.
 
- **Pillow (PIL)** : For image processing.
 
- **Matplotlib** : For plotting images.
 
- **tqdm** : For progress bars in the scripts.


---


## Installation 

### 1. Clone the Repository 


```bash
git clone https://github.com/your_username/virtual_screening_pipeline.git
cd virtual_screening_pipeline
```

### 2. Install Python Dependencies 
Use the provided `requirements.txt` file to install necessary Python packages:

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


---


## Setup 

### Updating Script Paths 

Before running the pipeline, update the paths in the scripts to match your system:
`smiles_to_pdbqt.py`** Update the paths to `pythonsh` and `prepare_ligand4.py`:

```python
command = [
    "/path/to/mgltools/bin/pythonsh",
    "/path/to/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py",
    '-l', pdb_filename,
    '-o', pdbqt_filename
]
```
`run_large_screen.sh`** 
Update the path to the AutoDock Vina GPU executable:


```bash
VINA_GPU_EXEC="/path/to/AutoDock-Vina-GPU"
```
`summarize_virtual_screen.py`** 
Update the font path for image annotations:


```python
try:
    font = ImageFont.truetype("/path/to/fonts/Arial.TTF", size=50)
except IOError:
    print("Arial font not found, using default")
    font = ImageFont.load_default()
```


---


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


---


## Usage 
Master Script: `master.sh`The `master.sh`**  script automates the entire pipeline, including ligand preparation, docking, and result summarization. It also ensures that previously processed directories are skipped, allowing for smooth resumption of interrupted runs.
#### Usage 


```bash
bash master.sh -c CONFIG_FILE -s SMILES_FILE -l LIGAND_PARENT_DIR [OPTIONS]
```

#### Required Arguments 
 
- **`-c CONFIG_FILE`** : Path to the AutoDock Vina configuration file.
 
- **`-s SMILES_FILE`** : Path to the SMILES file containing ligand information.
 
- **`-l LIGAND_PARENT_DIR`** : Directory to store converted ligands (PDBQT files).

#### Optional Arguments 
 
- **`-o OUTPUT_DIR`** : Directory for docking results (default: `output_directory`).
 
- **`-a AGGREGATE_DIR`** : Directory to store aggregated results (default: `aggregate_results`).
 
- **`-n TOP_N`** : Number of top hits to generate images for (default: 15).
 
- **`-d`** : Treat `LIGAND_PARENT_DIR` as a single directory (no subdirectories).
 
- **`-s START_FROM_DIR`** : Start processing from this directory number (default: 1).
 
- **`-v`** : Enable verbose output.


---


### Individual Scripts 
`smiles_to_pdbqt.py`
Converts SMILES strings to PDBQT files, splitting outputs into subdirectories.
**Usage** :

```bash
python smiles_to_pdbqt.py INPUT_FILE OUTPUT_DIR [OPTIONS]
```
**Options** : 
- **`-n`** , ****`-n`** , `--num_workers`** : Number of worker processes.
 
- **`-m`** , ****`-m`** , `--max_molecules`** : Maximum number of molecules to process.
 
- **`--molecules_per_dir`** : Number of molecules per subdirectory.


---

`run_large_screen.sh`
Runs AutoDock Vina docking on batches of ligands, skipping directories that already contain output files.
**Usage** :

```bash
bash run_large_screen.sh -c CONFIG_FILE -l LIGAND_DIR -o OUTPUT_DIR -a AGGREGATE_OUTPUT_DIR [OPTIONS]
```
**Options** : 
- **`-s START_FROM_DIR`** : Start processing from this directory number.
 
- **`--single-directory`** : Treat `LIGAND_DIR` as a single directory.


---

`summarize_virtual_screen.py`
Summarizes docking results and generates images for top hits.
**Usage** :

```bash
python summarize_virtual_screen.py -d DIRECTORY -s SMILES_FILE [OPTIONS]
```
**Options** : 
- **`-g`** , ****`-g`** , `--generate_images`** : Enable image generation.
 
- **`-n`** , ****`-n`** , `--top_n`** : Number of top hits to generate images for.
 
- **`-i`** , ****`-i`** , `--image_dir`** : Directory to save images.
 
- **`-f`** , ****`-f`** , `--figure_file`** : Filename for the combined figure.


---


## Example Command 


```bash
bash master.sh \
  -c conf.txt \
  -s data/zinc_smiles.txt \
  -l data/ligands \
  -o results/docking_outputs \
  -a results/aggregated_results \
  -n 10 \
  -v \
  -d
```

This command will:
 
1. **Convert SMILES to PDBQT** : Processes `zinc_smiles.txt` into PDBQT files.
 
2. **Run Docking** : Performs virtual screening with AutoDock Vina, skipping already processed directories.
 
3. **Summarize Results** : Generates a summary CSV and images for the top 10 hits.


---


## Notes 

- Ensure all paths in the scripts are correctly set to match your environment.

- The pipeline is designed to resume interrupted runs by skipping directories that already contain output files.


---


## Troubleshooting 
 
- **Font Errors** : Modify the font path or use the default font in `summarize_virtual_screen.py`.
 
- **Connection of Scripts** : Ensure output directories and file names are consistent across scripts.


---


## License 

This project is licensed under the MIT License.


---


## Acknowledgments 
 
- **AutoDock Vina**  for molecular docking software.
 
- **RDKit**  for cheminformatics tools.
 
- **MGLTools**  for ligand preparation utilities.