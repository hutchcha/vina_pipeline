import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from tqdm import tqdm
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
import math

def compute_fingerprints(smiles_file):
    """
    Compute Morgan fingerprints for all molecules in a SMILES file.
    """
    fingerprints = []
    ids = []
    with open(smiles_file, 'r') as f:
        for line in tqdm(f, desc="Generating fingerprints", unit="molecule"):
            parts = line.strip().split()
            if len(parts) != 2:
                continue  # Skip malformed lines
            smiles, zincid = parts[0], parts[1]
            
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                fingerprints.append(fp)
                ids.append(zincid)
    return fingerprints, ids

def pairwise_similarity_chunk(args):
    """
    Compute pairwise similarities for a subset of rows.
    args: (start_i, end_i, fingerprints)
    """
    start_i, end_i, fingerprints = args
    n = len(fingerprints)
    sim_chunk = []
    # Loop over the assigned range of i
    for i in range(start_i, end_i):
        fp_i = fingerprints[i]
        for j in range(i+1, n):
            sim = DataStructs.FingerprintSimilarity(fp_i, fingerprints[j])
            sim_chunk.append(sim)
    return sim_chunk

def compute_pairwise_similarity_parallel(fingerprints, num_processes):
    """
    Compute pairwise similarities in parallel.
    """
    n = len(fingerprints)
    # We will divide the range [0, n) into tasks_count tasks
    # For large n, a tasks_count ~ num_processes*4 or num_processes*10 is often good
    tasks_count = num_processes * 4
    lines_per_task = math.ceil(n / tasks_count)

    tasks = []
    for start in range(0, n, lines_per_task):
        end = min(start + lines_per_task, n)
        tasks.append((start, end, fingerprints))
    
    similarities = []
    # Use a Pool of workers
    with Pool(num_processes) as pool:
        # imap_unordered returns results as soon as they are ready
        # We can use a progress bar by counting the tasks
        for result in tqdm(pool.imap_unordered(pairwise_similarity_chunk, tasks),
                           total=len(tasks),
                           desc="Calculating similarities"):
            similarities.extend(result)

    return similarities

def plot_similarity_histogram(similarities, output_path):
    """
    Plot a histogram of the pairwise similarities.
    """
    plt.figure(figsize=(8, 6))
    plt.hist(similarities, bins=50, alpha=0.75, color='blue')
    plt.title("Pairwise Tanimoto Similarity Distribution")
    plt.xlabel("Tanimoto Similarity")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()
    print(f"Similarity histogram saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Summarize molecular similarity using Morgan fingerprints with parallelization.")
    parser.add_argument("--input", required=True, help="Path to the SMILES text file.")
    parser.add_argument("--output_hist", required=True, help="Path to save the similarity histogram.")
    parser.add_argument("--output_summary", required=True, help="Path to save the similarity summary statistics.")
    parser.add_argument("--num_processes", type=int, default=cpu_count(),
                        help="Number of parallel processes to use (default: all cores).")
    args = parser.parse_args()

    # Step 1: Compute Morgan fingerprints
    print("Computing fingerprints...")
    fingerprints, ids = compute_fingerprints(args.input)
    
    # Step 2: Compute pairwise Tanimoto similarities in parallel
    print("Calculating pairwise similarities...")
    similarities = compute_pairwise_similarity_parallel(fingerprints, args.num_processes)
    
    # Step 3: Summarize similarities
    summary = {
        "Mean Similarity": np.mean(similarities),
        "Median Similarity": np.median(similarities),
        "Min Similarity": np.min(similarities),
        "Max Similarity": np.max(similarities)
    }
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(args.output_summary, index=False)
    print(f"Similarity summary saved to: {args.output_summary}")
    print(summary_df)

    # Step 4: Plot histogram of similarities
    plot_similarity_histogram(similarities, args.output_hist)

if __name__ == "__main__":
    main()
