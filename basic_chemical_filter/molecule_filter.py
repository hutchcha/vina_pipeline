import argparse
import logging
import os
import gzip
import random
import csv
import time
import multiprocessing
from typing import Tuple, Generator
from multiprocessing import Manager

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import QED

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class MoleculeFilter:
    def __init__(self,
                 mw_range: Tuple[float, float] = (150, 500),
                 hba_range: Tuple[int, int] = (0, 10),
                 hbd_range: Tuple[int, int] = (0, 5),
                 tpsa_range: Tuple[float, float] = (60, 560),
                 rotb_range: Tuple[int, int] = (0, 15),
                 drug_score_range: Tuple[float, float] = (0.5, 1),
                 # New: include logP and aromatic ring filters
                 logp_range: Tuple[float, float] = (-5, 5),
                 aromatic_ring_range: Tuple[int, int] = (0, 5)
                ):
        """
        Initialize the MoleculeFilter with specified thresholds.
        """
        self.mw_min, self.mw_max = mw_range
        self.hba_min, self.hba_max = hba_range
        self.hbd_min, self.hbd_max = hbd_range
        self.tpsa_min, self.tpsa_max = tpsa_range
        self.rotb_min, self.rotb_max = rotb_range
        self.drug_score_min, self.drug_score_max = drug_score_range
        
        # New: store logP and aromatic ring cutoff ranges
        self.logp_min, self.logp_max = logp_range
        self.arom_ring_min, self.arom_ring_max = aromatic_ring_range

    def calculate_drug_score(self, mol: Chem.Mol) -> float:
        """
        Calculate the drug score using RDKit's QED as a proxy.
        QED is a measure of drug-likeness.
        """
        try:
            qed_score = QED.qed(mol)
            return qed_score
        except:
            return 0.0  # Return a default score if calculation fails
    
    def compute_descriptors(self, mol: Chem.Mol):
        """
        Compute all necessary descriptors for the molecule.
        """
        descriptors = {
            'MolecularWeight': Descriptors.MolWt(mol),
            'HBA': Descriptors.NumHAcceptors(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'TPSA': rdMolDescriptors.CalcTPSA(mol),
            'RotatableBonds': Descriptors.NumRotatableBonds(mol),
            'DrugScore': self.calculate_drug_score(mol),
            # New: compute LogP and # of aromatic rings
            'LogP': Descriptors.MolLogP(mol),
            'NumAromaticRings': rdMolDescriptors.CalcNumAromaticRings(mol)
        }
        return descriptors
    
    def filter_molecule(self, smiles: str):
        """
        Determine whether a molecule meets the filtering criteria.

        Returns:
            (bool, dict): Pass/Fail and descriptors dictionary.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.warning(f"Invalid SMILES string: {smiles}")
            return False, {}
        
        d = self.compute_descriptors(mol)
        
        # Apply filtering criteria
        if not (self.mw_min <= d['MolecularWeight'] <= self.mw_max):
            return False, d
        if not (self.hba_min <= d['HBA'] <= self.hba_max):
            return False, d
        if not (self.hbd_min <= d['HBD'] <= self.hbd_max):
            return False, d
        if not (self.tpsa_min <= d['TPSA'] <= self.tpsa_max):
            return False, d
        if not (self.rotb_min <= d['RotatableBonds'] <= self.rotb_max):
            return False, d
        if not (self.drug_score_min <= d['DrugScore'] <= self.drug_score_max):
            return False, d
        
        # New: apply logP and aromatic ring filters
        if not (self.logp_min <= d['LogP'] <= self.logp_max):
            return False, d
        if not (self.arom_ring_min <= d['NumAromaticRings'] <= self.arom_ring_max):
            return False, d
        
        return True, d


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Filter molecules based on physicochemical parameters with parallelization."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input directory containing .smi.gz files."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output CSV file."
    )
    parser.add_argument(
        "--mw_min", type=float, default=150,
        help="Minimum molecular weight (default: 150)."
    )
    parser.add_argument(
        "--mw_max", type=float, default=500,
        help="Maximum molecular weight (default: 500)."
    )
    parser.add_argument(
        "--hba_min", type=int, default=0,
        help="Minimum HBA (default: 0)."
    )
    parser.add_argument(
        "--hba_max", type=int, default=10,
        help="Maximum HBA (default: 10)."
    )
    parser.add_argument(
        "--hbd_min", type=int, default=0,
        help="Minimum HBD (default: 0)."
    )
    parser.add_argument(
        "--hbd_max", type=int, default=5,
        help="Maximum HBD (default: 5)."
    )
    parser.add_argument(
        "--tpsa_min", type=float, default=60,
        help="Minimum TPSA (default: 60)."
    )
    parser.add_argument(
        "--tpsa_max", type=float, default=560,
        help="Maximum TPSA (default: 560)."
    )
    parser.add_argument(
        "--rotb_min", type=int, default=0,
        help="Minimum rotatable bonds (default: 0)."
    )
    parser.add_argument(
        "--rotb_max", type=int, default=15,
        help="Maximum rotatable bonds (default: 15)."
    )
    parser.add_argument(
        "--drugscore_min", type=float, default=0.5,
        help="Minimum drug score (QED) (default: 0.5)."
    )
    parser.add_argument(
        "--drugscore_max", type=float, default=1,
        help="Maximum drug score (QED) (default: 1)."
    )
    
    # New: arguments for LogP range
    parser.add_argument(
        "--logp_min", type=float, default=-5,
        help="Minimum LogP (default: -5)."
    )
    parser.add_argument(
        "--logp_max", type=float, default=5,
        help="Maximum LogP (default: 5)."
    )
    
    # New: arguments for aromatic ring range
    parser.add_argument(
        "--arom_ring_min", type=int, default=0,
        help="Minimum number of aromatic rings (default: 0)."
    )
    parser.add_argument(
        "--arom_ring_max", type=int, default=5,
        help="Maximum number of aromatic rings (default: 5)."
    )

    parser.add_argument(
        "--max_count", type=int, default=None,
        help="Optional maximum number of molecules to process/keep overall."
    )
    parser.add_argument(
        "--per_file_limit", type=int, default=10000,
        help="Maximum number of molecules to select from each file (default: 10000)."
    )
    parser.add_argument(
        "--num_processes", type=int, default=4,
        help="Number of processes to use for parallelization (default: 4)."
    )
    parser.add_argument(
        "--progress_interval", type=int, default=10,
        help="Interval in seconds to print progress updates (default: 10)."
    )
    parser.add_argument(
        "--random_seed", type=int, default=None,
        help="Optional random seed for reproducibility. If not set, shuffle is random each run."
    )
    return parser.parse_args()


def find_smi_gz_files(dir_path: str):
    """
    Recursively find all .smi.gz files in the given directory.
    """
    smi_files = []
    for root, dirs, files in os.walk(dir_path):
        for file in files:
            if file.endswith(".smi.gz"):
                smi_files.append(os.path.join(root, file))
    return smi_files


def load_smiles_from_gz(file_path: str) -> Generator[Tuple[str, str], None, None]:
    """
    Load SMILES and associated ID from a gzipped .smi file.
    Format: SMILES <tab> ID
    """
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                smiles = parts[0]
                mol_id = parts[1]
                yield (smiles, mol_id)
            else:
                logging.warning(f"Invalid line format in {file_path}: {line}")


def gzip_file(file_path: str):
    """
    Gzip the specified file to file_path.gz.
    """
    gz_path = file_path + '.gz'
    with open(file_path, 'rb') as f_in, gzip.open(gz_path, 'wb') as f_out:
        f_out.writelines(f_in)
    logging.info(f"Compressed {file_path} to {gz_path}")


def process_files(args_tuple):
    (files, filter_params, output_path, max_count, passed_count, increment_lock, per_file_limit) = args_tuple
    
    logging.info(f"Worker started. Handling {len(files)} files.")
    logging.info(f"Files: {files}")

    # Instantiate the filter object
    filter_obj = MoleculeFilter(
        mw_range=(filter_params['mw_min'], filter_params['mw_max']),
        hba_range=(filter_params['hba_min'], filter_params['hba_max']),
        hbd_range=(filter_params['hbd_min'], filter_params['hbd_max']),
        tpsa_range=(filter_params['tpsa_min'], filter_params['tpsa_max']),
        rotb_range=(filter_params['rotb_min'], filter_params['rotb_max']),
        drug_score_range=(filter_params['drugscore_min'], filter_params['drugscore_max']),
        # Pass new filters to MoleculeFilter
        logp_range=(filter_params['logp_min'], filter_params['logp_max']),
        aromatic_ring_range=(filter_params['arom_ring_min'], filter_params['arom_ring_max'])
    )
    
    # Extended fieldnames for new descriptors
    fieldnames = [
        'SMILES', 'ID', 'MolecularWeight', 'HBA', 'HBD', 'TPSA',
        'RotatableBonds', 'DrugScore', 'LogP', 'NumAromaticRings'
    ]
    
    try:
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for smi_file in files:
                file_count = 0
                file_passed = 0
                logging.info(f"Processing file: {smi_file}")
                try:
                    for smi, zid in load_smiles_from_gz(smi_file):
                        file_count += 1
                        
                        # Check global passed_count for max_count
                        if max_count is not None and passed_count.value >= max_count:
                            logging.info("Max count reached, stopping worker.")
                            return
                        
                        # If we've already reached the per-file limit of selected molecules, stop this file
                        if file_passed >= per_file_limit:
                            logging.info(f"Per-file limit of {per_file_limit} reached for file {smi_file}.")
                            break
                        
                        try:
                            passed, descriptors = filter_obj.filter_molecule(smi)
                        except Exception as e:
                            logging.error(f"Error filtering molecule {smi}: {e}")
                            continue
                        
                        if passed:
                            with increment_lock:
                                if max_count is not None and passed_count.value >= max_count:
                                    logging.info("Max count reached, stopping worker.")
                                    return
                                passed_count.value += 1
                                current_count = passed_count.value
                            
                            molecule_data = {
                                'SMILES': smi,
                                'ID': zid,
                                'MolecularWeight': descriptors['MolecularWeight'],
                                'HBA': descriptors['HBA'],
                                'HBD': descriptors['HBD'],
                                'TPSA': descriptors['TPSA'],
                                'RotatableBonds': descriptors['RotatableBonds'],
                                'DrugScore': descriptors['DrugScore'],
                                'LogP': descriptors['LogP'],
                                'NumAromaticRings': descriptors['NumAromaticRings']
                            }
                            writer.writerow(molecule_data)
                            file_passed += 1
                            
                            if max_count is not None and current_count >= max_count:
                                logging.info("Max count reached during writing, stopping worker.")
                                return
                            
                            # Check per-file limit after writing the molecule
                            if file_passed >= per_file_limit:
                                logging.info(f"Per-file limit of {per_file_limit} reached for file {smi_file}.")
                                break
                except Exception as e:
                    logging.error(f"Error reading file {smi_file}: {e}")
                
                logging.info(f"Finished file {smi_file}: processed {file_count} molecules, {file_passed} passed.")
    except Exception as e:
        logging.error(f"Error in process_files writing to {output_path}: {e}")


def main():
    args = parse_arguments()
    if args.random_seed is not None:
        random.seed(args.random_seed)
        logging.info(f"Random seed set to {args.random_seed}")
    else:
        logging.info("No random seed set, file order will vary each run")
    
    smi_files = find_smi_gz_files(args.input)
    if not smi_files:
        logging.error("No .smi.gz files found in the provided directory.")
        return
    
    random.shuffle(smi_files)
    
    # Adjust number of processes if there are fewer files than processes
    if args.num_processes > len(smi_files):
        args.num_processes = len(smi_files)
    
    # Split files among processes
    chunk_size = len(smi_files) // args.num_processes
    file_chunks = [smi_files[i*chunk_size:(i+1)*chunk_size] for i in range(args.num_processes)]
    remainder = len(smi_files) % args.num_processes
    if remainder > 0:
        file_chunks[-1].extend(smi_files[-remainder:])
    
    manager = Manager()
    passed_count = manager.Value('i', 0)
    increment_lock = manager.Lock() 
    filter_params = {
        'mw_min': args.mw_min, 'mw_max': args.mw_max,
        'hba_min': args.hba_min, 'hba_max': args.hba_max,
        'hbd_min': args.hbd_min, 'hbd_max': args.hbd_max,
        'tpsa_min': args.tpsa_min, 'tpsa_max': args.tpsa_max,
        'rotb_min': args.rotb_min, 'rotb_max': args.rotb_max,
        'drugscore_min': args.drugscore_min, 'drugscore_max': args.drugscore_max,
        # Pass new filters into the process pool
        'logp_min': args.logp_min, 'logp_max': args.logp_max,
        'arom_ring_min': args.arom_ring_min, 'arom_ring_max': args.arom_ring_max
    }
    
    # Temporary output files for each process
    temp_files = [f"{args.output}.part{i}" for i in range(args.num_processes)]
    
    process_args = [
        (
            file_chunks[i],
            filter_params,
            temp_files[i],
            args.max_count,
            passed_count,
            increment_lock,
            args.per_file_limit
        )
        for i in range(args.num_processes)
    ]
    
    pool = multiprocessing.Pool(args.num_processes)
    result = pool.map_async(process_files, process_args)
    
    # Progress monitoring loop
    while not result.ready():
        time.sleep(args.progress_interval)
        current_count = passed_count.value
        if args.max_count:
            logging.info(f"Selected {current_count}/{args.max_count} desired molecules.")
        else:
            logging.info(f"Selected {current_count} molecules so far.")
        
        if args.max_count is not None and current_count >= args.max_count:
            break
    
    pool.close()
    pool.join()
    
    final_count = passed_count.value
    logging.info(f"Total selected molecules: {final_count}")
    
    # Merge partial files
    if final_count > 0:
        fieldnames = [
            'SMILES', 'ID', 'MolecularWeight', 'HBA', 'HBD', 'TPSA',
            'RotatableBonds', 'DrugScore', 'LogP', 'NumAromaticRings'
        ]
        
        with open(args.output, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            for tf in temp_files:
                if os.path.exists(tf):
                    with open(tf, 'r') as infile:
                        reader = csv.DictReader(infile)
                        for row in reader:
                            writer.writerow(row)
                    os.remove(tf)
    else:
        # If no molecules were selected
        for tf in temp_files:
            if os.path.exists(tf):
                os.remove(tf)
        logging.info("No molecules passed the filtering criteria.")
        return
    
    # Gzip the final CSV
    gzip_file(args.output)


if __name__ == "__main__":
    main()
