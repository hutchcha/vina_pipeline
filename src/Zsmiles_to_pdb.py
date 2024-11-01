import os
import argparse
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from tqdm import tqdm
# Define the list of off-limits atomic symbols
OFF_LIMITS_SYMBOLS = ['Si', 'P', 'B', 'Fe']  # Add any other atoms to exclude

def process_smiles(args):
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

        # Proceed with the rest of the processing
        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        # Optimize geometry
        AllChem.UFFOptimizeMolecule(mol)

        # Write to temporary PDB file in the output directory
        pdb_filename = f"{zinc_id}.pdb"
        pdb_file = os.path.join(output_dir, pdb_filename)
        Chem.MolToPDBFile(mol, pdb_file)

        # Prepare the PDBQT filename
        pdbqt_filename = f"{zinc_id}.pdbqt"
        pdbqt_file = os.path.join(output_dir, pdbqt_filename)

        # Since we're changing the working directory, use filenames without paths
        command = [
            "/mnt/data_SSD/VINA/autodocktools/bin/pythonsh",
            "/mnt/data_SSD/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py",
            '-l', pdb_filename,
            '-o', pdbqt_filename
        ]

        # Set the working directory to output_dir
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=output_dir  # Set working directory
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

def read_smiles_file(file_path, max_molecules=None):
    smiles_list = []
    with open(file_path, 'r') as file:
        for idx, line in tqdm(enumerate(file)):
            if max_molecules is not None and idx >= max_molecules:
                break
            line = line.strip()
            if line:
                parts = line.strip().split()
                if len(parts) >= 3:
                    tranche = parts[0]
                    zinc_id = parts[1]
                    smiles_string = ' '.join(parts[2:])
                    smiles_list.append((zinc_id, smiles_string))
                else:
                    print(f"Invalid line format: {line}")
    return smiles_list

def main():
    parser = argparse.ArgumentParser(
        description="Convert ZINC SMILES to PDBQT files using multiprocessing"
    )
    parser.add_argument("input_file", type=str,
                        help="Path to the input ZINC SMILES file.")
    parser.add_argument("output_dir", type=str,
                        help="Path to the output directory.")
    parser.add_argument("-n", "--num_workers", type=int,
                        default=mp.cpu_count(),
                        help="Number of worker processes "
                             "(default: number of CPU cores)")
    parser.add_argument("-m", "--max_molecules", type=int, default=None,
                        help="Maximum number of molecules to process "
                             "(default: all)")
    args = parser.parse_args()

    # Read SMILES from file with a limit on the number of molecules
    smiles_data = read_smiles_file(args.input_file,
                                   max_molecules=args.max_molecules)

    # Create output directory if t doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Prepare arguments for multiprocessing
    pool_args = [
        (zinc_id, smiles_string, args.output_dir)
        for zinc_id, smiles_string in smiles_data
    ]

    # Start multiprocessing pool
    with mp.Pool(processes=args.num_workers) as pool:
        for _ in pool.imap_unordered(process_smiles, pool_args):
            pass  # The process_smiles function handles printing

if __name__ == "__main__":
    main()

