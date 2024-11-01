import os
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse

def read_smiles_file(file_path):
    smiles_list = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip().split()
            if line:
                smiles_list.append(line)
    return smiles_list

# Function to convert SMILES to PDB files, with subdirectory option
def convert_smiles_to_individual_sdf(smiles_list, output_dir, max_files_per_dir=None):
    os.makedirs(output_dir, exist_ok=True)
    dir_count = 0
    file_count = 0
    current_dir = os.path.join(output_dir, f"batch_{dir_count}")
    os.makedirs(current_dir, exist_ok=True)

    for smiles in smiles_list:
        zinc_id, smiles_string = smiles[0], smiles[1]
        mol = Chem.MolFromSmiles(smiles_string)
        
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            mol.SetProp("_Name", zinc_id)
            
            pdb_file_path = os.path.join(current_dir, f"{zinc_id}.pdb")
            with Chem.rdmolfiles.PDBWriter(pdb_file_path) as writer:
                writer.write(mol)
            
            print(f"Saved: {pdb_file_path}")
            file_count += 1

            # Check if we've reached the max files per directory
            if max_files_per_dir and file_count >= max_files_per_dir:
                dir_count += 1
                current_dir = os.path.join(output_dir, f"batch_{dir_count}")
                os.makedirs(current_dir, exist_ok=True)
                file_count = 0
        else:
            print(f"Failed to convert SMILES for {zinc_id}")

def main():
    parser = argparse.ArgumentParser(description="Convert ZINC SMILES Lists into individual SDF files")
    parser.add_argument("input_file", type=str, help="Path to the input ZINC SMILES file.")
    parser.add_argument("output_dir", type=str, help="Path for output files")
    parser.add_argument("--max_files_per_dir", type=int, default=None, help="Maximum number of files per subdirectory")
    args = parser.parse_args()

    smiles = read_smiles_file(args.input_file)
    convert_smiles_to_individual_sdf(smiles, output_dir=args.output_dir, max_files_per_dir=args.max_files_per_dir)

if __name__ == "__main__":
    main()
