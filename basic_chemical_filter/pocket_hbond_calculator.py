import MDAnalysis as mda
import mdtraj as md
import pandas as pd

# Define donor and acceptor atom names
DONOR_NAMES = {"N", "ND1", "NE", "NE2", "NH1", "NH2", "NZ", "OG", "OG1", "OH", "SD"}
ACCEPTOR_NAMES = {"O", "OD1", "OD2", "OE1", "OE2", "ND1", "NE2", "OG", "OG1", "OH", "SD"}

pocket_residues = {
    'p1': ['182', '183', '184', '203', '224', '225', '226', '227', '231', '232', '233', '334', '337', '338', '341', '342', '345'],
    'p2': ['187', '189', '190', '191', '192', '236', '238', '239', '240', '241', '242', '243', '244', '245', '246', '247', '248', '250', '253', '255', '258', '260', '272', '273', '275', '276', '279', '280', '283', '284', '291'],
    'p3': ['298', '299', '300', '301', '324', '325', '326', '328', '329', '332'],
    'p4': ['184', '255', '256', '257', '284', '287', '288', '289', '290', '291', '319', '320', '321', '343', '344', '345', '346', '347', '348', '349']
}

def identify_hbond_donors_acceptors_by_name(pdb_file, residue_ids):
    """
    Identify hydrogen bond donors and acceptors based on atom names for backbone and side chains
    of specified residues in a PDB file, considering only solvent-accessible atoms.

    Parameters:
    - pdb_file: Path to the PDB file.
    - residue_ids: List of residue IDs to analyze.

    Returns:
    - The number of donors and acceptors.
    """
    # Load the structure into MDAnalysis
    u = mda.Universe(pdb_file)

    # Select atoms using MDAnalysis's `resid`
    selected_atoms = u.select_atoms(f"resid {' '.join(map(str, residue_ids))}")
    
    # Write the selected atoms to a temporary PDB file for MDTraj
    temp_pdb_file = "temp_selected.pdb"
    selected_atoms.write(temp_pdb_file)

    # Load the PDB into MDTraj
    traj = md.load(temp_pdb_file)
    sasa = md.shrake_rupley(traj, mode='atom')  # SASA in nm²

    # Extract solvent-accessible atoms
    solvent_accessible_mask = sasa[0] > 0  # First frame only
    solvent_accessible_atoms = [atom for atom, accessible in zip(traj.topology.atoms, solvent_accessible_mask) if accessible]

    # Classify atoms as donors or acceptors based on names
    donor_count = 0
    acceptor_count = 0

    for atom in solvent_accessible_atoms:
        if atom.name in DONOR_NAMES:
            donor_count += 1
        if atom.name in ACCEPTOR_NAMES:
            acceptor_count += 1

    return donor_count, acceptor_count

# Main script
if __name__ == "__main__":
    pdb_file = "rheb_test.pdb"  # Replace with your PDB file path

    # Collect results for all pockets
    results = []
    for pocket, residues in pocket_residues.items():
        donors, acceptors = identify_hbond_donors_acceptors_by_name(pdb_file, residues)
        results.append({"Pocket": pocket, "Donors": donors, "Acceptors": acceptors})

    # Convert to a DataFrame
    df = pd.DataFrame(results)

    # Save to CSV
    output_csv = "hydrogen_bond_analysis.csv"
    df.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

    # Print the DataFrame
    print(df)
