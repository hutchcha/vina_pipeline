import os
import argparse
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from tqdm import tqdm
import pandas as pd

# List of off-limits atomic symbols
OFF_LIMITS_SYMBOLS = ['Si', 'P', 'B', 'Fe']  # Add any others you need to exclude

def process_smiles(args):
    """
    Convert a single SMILES to PDBQT format, excluding any molecule
    containing off-limits atoms. 
    """
    zinc_id, smiles_string, output_dir = args

    try:
        # Generate RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            print(f"Failed to parse SMILES for {zinc_id}")
            return

        # Check for off-limits atoms
        if any(atom.GetSymbol() in OFF_LIMITS_SYMBOLS for atom in mol.GetAtoms()):
            print(f"Skipping molecule {zinc_id} due to presence of off-limits atoms.")
            return

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)

        # Write to temporary PDB file
        pdb_filename = f"{zinc_id}.pdb"
        pdb_file = os.path.join(output_dir, pdb_filename)
        Chem.MolToPDBFile(mol, pdb_file)

        # Prepare the PDBQT filename
        pdbqt_filename = f"{zinc_id}.pdbqt"
        pdbqt_file = os.path.join(output_dir, pdbqt_filename)

        # Run prepare_ligand4.py script
        command = [
            "/mnt/data_SSD/VINA/autodocktools/bin/pythonsh",
            "/mnt/data_SSD/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py",
            '-l', pdb_filename,
            '-o', pdbqt_filename
        ]
        result = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, cwd=output_dir
        )

        if result.returncode != 0:
            print(f"Failed to convert {pdb_file} to PDBQT format for {zinc_id}")
            print(f"Error message: {result.stderr}")
            os.remove(pdb_file)
            return

        # Remove the temporary PDB file
        os.remove(pdb_file)
        print(f"Saved: {pdbqt_file}")

    except Exception as e:
        print(f"Error processing {zinc_id}: {e}")

def create_output_subdir(base_dir, subdir_index):
    """
    Create a subdirectory named batch_<index> under base_dir.
    """
    subdir_path = os.path.join(base_dir, f"batch_{subdir_index}")
    os.makedirs(subdir_path, exist_ok=True)
    return subdir_path

def main():
    parser = argparse.ArgumentParser(
        description="Convert ZINC SMILES to PDBQT files using multiprocessing"
    )
    parser.add_argument("input_csv", type=str, 
                        help="Path to the input CSV file (must have 'SMILES' and 'ID' columns).")
    parser.add_argument("output_dir", type=str, 
                        help="Path to the output directory.")
    parser.add_argument("-n", "--num_workers", type=int, default=mp.cpu_count(),
                        help="Number of worker processes (default: number of CPU cores)")
    parser.add_argument("-m", "--max_molecules", type=int, default=None,
                        help="Maximum number of molecules to process (default: all)")
    parser.add_argument("--molecules_per_dir", type=int, default=1000,
                        help="Number of molecules per subdirectory (default: 1000).")
    parser.add_argument("--single_dir", action="store_true",
                        help="If set, outputs all PDBQT files to a single directory.")

    args = parser.parse_args()

    # ------------------------------------------------------
    # (1) Read the CSV via pandas, focusing on SMILES and ID
    # ------------------------------------------------------
    df = pd.read_csv(args.input_csv, usecols=["SMILES", "ID"])

    # If you only want to process a subset:
    if args.max_molecules is not None:
        df = df.head(args.max_molecules)

    # Create output base directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------
    # (2) Build list of (ID, SMILES) to pass to workers
    # ------------------------------------------------------
    smiles_data = list(zip(df["ID"], df["SMILES"]))

    # ------------------------------------------------------
    # (3) Handle single or multiple subdirectories
    # ------------------------------------------------------
    if args.single_dir:
        current_dir = args.output_dir
        pool_args = [(zinc_id, smiles_str, current_dir) for zinc_id, smiles_str in smiles_data]

    else:
        current_dir_index = 0
        current_dir = create_output_subdir(args.output_dir, current_dir_index)
        molecules_in_current_dir = 0

        pool_args = []
        for zinc_id, smiles_str in smiles_data:
            if molecules_in_current_dir >= args.molecules_per_dir:
                current_dir_index += 1
                current_dir = create_output_subdir(args.output_dir, current_dir_index)
                molecules_in_current_dir = 0

            pool_args.append((zinc_id, smiles_str, current_dir))
            molecules_in_current_dir += 1

    # ------------------------------------------------------
    # (4) Run multiprocessing pool
    # ------------------------------------------------------
    with mp.Pool(processes=args.num_workers) as pool:
        for _ in pool.imap_unordered(process_smiles, pool_args):
            pass  # The process_smiles function prints status/info

if __name__ == "__main__":
    main()
