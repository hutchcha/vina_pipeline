## Vina Virtual Screening Pipeline

### Purpose

This repository contains a comprehensive and modular pipeline for performing high-throughput virtual screening of chemical compounds using AutoDock Vina. Designed to streamline and automate the preparation, docking, and analysis processes, this pipeline supports researchers in identifying potential ligand binders to target proteins efficiently.

By handling complex tasks such as SMILES-to-PDBQT conversions, subdirectory organization to manage memory constraints, and docking with user-defined configurations, this repository provides an end-to-end solution for large-scale virtual screening studies. The goal is to assist researchers and scientists in narrowing down the chemical space, ultimately identifying compounds with the most promising binding affinities for further investigation.

### Key Features

- **Automated Ligand Preparation**: Convert SMILES to PDBQT format with optional output subdirectory organization, accommodating memory limitations for extensive libraries.
- **Flexible Docking Configuration**: Easily integrate AutoDock Vina for batch processing of ligands, with customizable docking parameters for different receptor sites.
- **Efficient Results Management**: Organize docking results and analyze scores, allowing for straightforward identification of top candidates for further validation.

### Applications

Created for my own research, but could be potentially useful for others in my lab and outside of it.


## SMILES to PDBQT Conversion Script (smiles_to_pdbqt.py)

This script converts a list of SMILES strings from a file into individual PDBQT files using RDKit, with an option to split output files into subdirectories to manage large datasets.

Smiles file must be formatted as:
ZINCID SMILES

this is the typical way to download from the ZINC database.

### Usage

```bash
python src/main.py <input_file> <output_dir> [--max_files_per_dir <number>]

### Example
python src/main.py zinc_smiles.txt output_dir --max_files_per_dir 1000

