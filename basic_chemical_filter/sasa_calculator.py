#!/usr/bin/env python3
"""
sasa_calculator.py

A script to calculate polar, nonpolar, and total SASA for specified binding pockets
from molecular dynamics trajectories.

Usage:
    python sasa_calculator.py \
        --trajectory /path/to/trajectory.xtc \
        --topology /path/to/topology.pdb \
        --pocket pocket1:10,20,30 \
        --pocket pocket2:40,50,60 \
        --output output_sasa.csv

Author: Your Name
Date: YYYY-MM-DD
"""

import argparse
import mdtraj as md
import pandas as pd
import os
import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("sasa_calculator.log"),
        logging.StreamHandler()
    ]
)

# Predefined lists of polar and nonpolar residues based on standard amino acid properties
POLAR_RESIDUES = {
    'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU',
    'HIS', 'LYS', 'SER', 'THR', 'TYR', 'TRP'
}

NONPOLAR_RESIDUES = {
    'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE',
    'PRO', 'GLY', 'TRP', 'TYR', 'CYS'
}

# Note: Some residues like TRP and TYR have both polar and nonpolar characteristics.
# Depending on the context, you might categorize them differently.
# Adjust the lists as needed based on your specific requirements.

def parse_arguments():
    """
    Parse command-line arguments using argparse.

    Returns:
        args: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Calculate polar, nonpolar, and total SASA for specified binding pockets."
    )

    parser.add_argument(
        "--trajectory",
        required=True,
        help="Path to the trajectory file (e.g., .xtc, .dcd)."
    )

    parser.add_argument(
        "--topology",
        required=True,
        help="Path to the topology file (e.g., .pdb, .prmtop)."
    )

    parser.add_argument(
        "--pocket",
        action='append',
        required=True,
        help=(
            "Define a pocket and its residue IDs in the format "
            "'pocket_name:res1,res2,res3'. "
            "Use this argument multiple times for multiple pockets."
        )
    )

    parser.add_argument(
        "--output",
        default="sasa_results.csv",
        help="Output CSV file to save SASA results."
    )

    return parser.parse_args()

def parse_pockets(pocket_args):
    """
    Parse pocket arguments into a dictionary.

    Args:
        pocket_args (list): List of pocket arguments in 'name:res1,res2' format.

    Returns:
        dict: Mapping of pocket names to lists of residue IDs.
    """
    pocket_dict = {}
    for pocket in pocket_args:
        try:
            name, residues = pocket.split(":")
            res_list = [int(res.strip()) for res in residues.split(",") if res.strip().isdigit()]
            if not res_list:
                raise ValueError
            pocket_dict[name.strip()] = res_list
        except ValueError:
            logging.error(f"Error parsing pocket argument: '{pocket}'. Expected format 'name:res1,res2,...' with integer residue IDs.")
            sys.exit(1)
    return pocket_dict

def load_trajectory(trajectory_file, topology_file):
    """
    Load the trajectory using MDTraj.

    Args:
        trajectory_file (str): Path to the trajectory file.
        topology_file (str): Path to the topology file.

    Returns:
        md.Trajectory: Loaded trajectory object.
    """
    logging.info(f"Loading trajectory from {trajectory_file} with topology {topology_file}...")
    try:
        traj = md.load(trajectory_file, top=topology_file)
    except Exception as e:
        logging.error(f"Error loading trajectory or topology files: {e}")
        sys.exit(1)
    logging.info(f"Trajectory loaded: {traj.n_frames} frames, {traj.n_atoms} atoms.")
    return traj

def map_resSeq_to_indices(traj):
    """
    Create a mapping from residue sequence numbers (resSeq) to atom indices.

    Args:
        traj (md.Trajectory): Loaded trajectory.

    Returns:
        dict: Mapping from resSeq to list of atom indices.
    """
    resSeq_to_indices = {}
    for residue in traj.topology.residues:
        resSeq = residue.resSeq
        atom_indices = [atom.index for atom in residue.atoms]
        resSeq_to_indices[resSeq] = atom_indices
    return resSeq_to_indices

def select_atoms(traj, residues, resSeq_to_indices):
    """
    Select atom indices for given residue IDs.

    Args:
        traj (md.Trajectory): Loaded trajectory.
        residues (list): List of residue sequence numbers.
        resSeq_to_indices (dict): Mapping from resSeq to atom indices.

    Returns:
        list: List of atom indices corresponding to the specified residues.
    """
    selected_indices = []
    for res in residues:
        if res in resSeq_to_indices:
            selected_indices.extend(resSeq_to_indices[res])
        else:
            logging.warning(f"Residue ID {res} not found in topology.")
    return selected_indices

def calculate_sasa(traj, polar_residues, nonpolar_residues):
    """
    Calculate polar, nonpolar, and total SASA for selected atoms.

    Args:
        traj (md.Trajectory): Sub-trajectory containing selected atoms.
        polar_residues (set): Set of polar residue names.
        nonpolar_residues (set): Set of nonpolar residue names.

    Returns:
        tuple: (polar_sasa, nonpolar_sasa, total_sasa) as numpy arrays.
    """
    # Calculate SASA using Shrake-Rupley algorithm
    sasa = md.shrake_rupley(traj, mode='atom')  # Shape: (n_frames, n_atoms)

    # Map atom indices to residue names
    atom_to_residue = [atom.residue.name.upper() for atom in traj.topology.atoms]

    # Identify polar and nonpolar atoms
    polar_atoms = [
        idx for idx, res in enumerate(atom_to_residue) if res in polar_residues
    ]
    nonpolar_atoms = [
        idx for idx, res in enumerate(atom_to_residue) if res in nonpolar_residues
    ]

    # Calculate SASA
    polar_sasa = sasa[:, polar_atoms].sum(axis=1) * 100  # Convert to Å²
    nonpolar_sasa = sasa[:, nonpolar_atoms].sum(axis=1) * 100  # Convert to Å²
    total_sasa = sasa.sum(axis=1) * 100  # Convert to Å²

    return polar_sasa, nonpolar_sasa, total_sasa

def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Parse pocket residues
    pocket_residues = parse_pockets(args.pocket)

    # Load trajectory
    traj = load_trajectory(args.trajectory, args.topology)

    # Create resSeq to atom indices mapping
    resSeq_to_indices = map_resSeq_to_indices(traj)

    # Initialize dictionary to store SASA data
    sasa_framewise = {}

    # Process each pocket
    for pocket, residues in pocket_residues.items():
        logging.info(f"\nProcessing pocket '{pocket}' with residues {residues}...")

        # Select atom indices for the specified residues
        selected_atom_indices = select_atoms(traj, residues, resSeq_to_indices)
        logging.info(f"Selected {len(selected_atom_indices)} atoms for pocket '{pocket}'.")

        if len(selected_atom_indices) == 0:
            logging.warning(f"No atoms found for pocket '{pocket}'. Skipping...")
            continue

        # Create a sub-trajectory for the selected atoms
        sub_traj = traj.atom_slice(selected_atom_indices)
        logging.info(f"Sub-trajectory: {sub_traj.n_frames} frames, {sub_traj.n_atoms} atoms.")

        # Calculate SASA
        polar_sasa, nonpolar_sasa, total_sasa = calculate_sasa(
            sub_traj, POLAR_RESIDUES, NONPOLAR_RESIDUES
        )

        # Store SASA data
        sasa_framewise[pocket] = {
            'Polar_SASA': polar_sasa,
            'Nonpolar_SASA': nonpolar_sasa,
            'Total_SASA': total_sasa
        }

    if not sasa_framewise:
        logging.error("No SASA data calculated. Exiting.")
        sys.exit(1)

    # Prepare data for saving
    # Assuming all pockets have the same number of frames
    frames = range(traj.n_frames)
    data = {'Frame': frames}

    for pocket, sasa_values in sasa_framewise.items():
        data[f"{pocket}_Polar_SASA"] = sasa_values['Polar_SASA']
        data[f"{pocket}_Nonpolar_SASA"] = sasa_values['Nonpolar_SASA']
        data[f"{pocket}_Total_SASA"] = sasa_values['Total_SASA']

    df = pd.DataFrame(data)
    logging.info(f"\nSaving SASA results to {args.output}...")
    try:
        df.to_csv(args.output, index=False)
    except Exception as e:
        logging.error(f"Error saving output file: {e}")
        sys.exit(1)
    logging.info("SASA calculation and saving completed successfully.")

if __name__ == "__main__":
    main()
