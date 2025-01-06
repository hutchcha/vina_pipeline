import os
import sys
import argparse
import pandas as pd
import concurrent.futures

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt

# We use tqdm for the progress bar
from tqdm import tqdm


def extract_zinc_id(filename):
    """
    Extracts the ZINC ID from a filename ending with '_out.pdbqt'.
    E.g., 'ZINC12345678_out.pdbqt' -> 'ZINC12345678'
    """
    return filename.split('_out.pdbqt')[0]


def extract_docking_score(pdbqt_file, verbose=False):
    """
    Extracts the docking score from a PDBQT file by scanning
    for the line starting with 'REMARK VINA RESULT'.
    """
    score = None
    with open(pdbqt_file, 'r') as file:
        for line in file:
            if line.startswith("REMARK VINA RESULT"):
                try:
                    score = float(line.split()[3])  # The energy score is the 4th token
                except (IndexError, ValueError):
                    if verbose:
                        print(f"Warning: Could not parse docking score in file {pdbqt_file}")
                break
    return score


def parse_pdbqt_file(pdbqt_file, verbose=False):
    """
    Function suited for parallel execution:
      - Extracts the docking score
      - Extracts the ZINC ID
    Returns a tuple: (zinc_id, docking_score)
    """
    zinc_id = extract_zinc_id(os.path.basename(pdbqt_file))
    score = extract_docking_score(pdbqt_file, verbose=verbose)
    return (zinc_id, score)


def smiles_to_2d_structure(smiles, output_image_name, label_text=None, score=None):
    """
    Converts a SMILES string to a 2D molecular image with an optional
    label (ZINC ID and docking score).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES string: {smiles}")
        return False

    # Sanitize and recalculate aromaticity
    Chem.SanitizeMol(mol)

    # Draw the 2D structure as a PIL Image
    img = Draw.MolToImage(mol, size=(550, 500))
    img = img.convert("RGBA")
    img = img.crop(img.getbbox())

    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("/home/chutchins@uthouston.edu/fonts/Arial.TTF", size=50)
    except IOError:
        print("Arial font not found, using default font.")
        font = ImageFont.load_default()

    if label_text and score is not None:
        text_label = f"    {label_text}\n    {score} kcal/mol"
        text_bbox = draw.textbbox((0, 0), text_label, font=font)
        text_width = text_bbox[2] - text_bbox[0]
        text_height = text_bbox[3] - text_bbox[1]
        position = ((img.width - text_width - 110) // 2, img.height - text_height - 80)
        draw.text(position, text_label, fill="black", font=font)

    img.save(output_image_name)
    print(f"Image saved to {output_image_name}")
    return True


def create_combined_figure(image_paths, figure_file, top_n, verbose=False):
    """
    Creates a combined figure of the individual molecule images in a grid layout.
    """
    import math

    cols = min(5, top_n)
    rows = math.ceil(top_n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(4, 1.875))
    axes = axes.flatten()

    for ax, img_path in zip(axes, image_paths):
        img = Image.open(img_path)
        ax.imshow(img)
        ax.axis('off')

    # Hide any remaining empty subplots
    for ax in axes[len(image_paths):]:
        ax.axis('off')

    plt.tight_layout()
    plt.savefig(figure_file, dpi=1200)
    plt.close()

    if verbose:
        print(f"Combined figure saved to {figure_file}")


def generate_image_parallel(hit_info):
    """
    Helper function to generate an image in parallel.
    Expects a dict or tuple with needed info to build the image.
    """
    zinc_id = hit_info['ID']
    smiles = hit_info['SMILES']
    score = hit_info['DockingScore']
    image_dir = hit_info['image_dir']

    image_filename = f"{zinc_id}.png"
    image_path = os.path.join(image_dir, image_filename)
    success = smiles_to_2d_structure(smiles, image_path, label_text=zinc_id, score=score)
    if success:
        return image_path
    return None


def summarize_docking_results(
    csv_file,
    directory,
    generate_images=False,
    top_n=15,
    image_dir="images",
    figure_file="top_hits_figure.png",
    verbose=False
):
    """
    Reads an existing CSV that contains at least columns:
        - ID
        - SMILES
    Searches (recursively) through `directory` (and subdirectories) for PDBQT files,
    extracts the docking score, and appends it as a new column ('DockingScore') in
    the same CSV. Also supports generating images for the top N hits by docking score.
    """

    # Load the existing CSV into a DataFrame
    if verbose:
        print(f"Reading CSV file: {csv_file}")
    df = pd.read_csv(csv_file)

    # Create the DockingScore column if it doesn't exist
    if 'DockingScore' not in df.columns:
        df['DockingScore'] = None

    # Gather all *_out.pdbqt files from directory + subdirectories
    if verbose:
        print(f"Recursively scanning '{directory}' for PDBQT files...")
    pdbqt_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith("_out.pdbqt"):
                pdbqt_files.append(os.path.join(root, file))

    if verbose:
        print(f"Found {len(pdbqt_files)} files to parse.")

    # ----------------------------------------------------------
    # 1) Parse all PDBQT files in parallel, collecting results
    #    in a dictionary mapping {zinc_id -> docking_score}.
    #    We'll use a tqdm progress bar here.
    # ----------------------------------------------------------
    docking_scores = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        parse_futures = [executor.submit(parse_pdbqt_file, pdbqt_file, verbose)
                         for pdbqt_file in pdbqt_files]

        for future in tqdm(
            concurrent.futures.as_completed(parse_futures),
            total=len(parse_futures),
            desc="Parsing PDBQT files"
        ):
            zinc_id, score = future.result()
            docking_scores[zinc_id] = score

    # ----------------------------------------------------------
    # 2) Update the DataFrame *once* at the end, rather than
    #    row-by-row.
    # ----------------------------------------------------------
    if verbose:
        print("Updating DataFrame with docking scores...")

    # For faster lookups, convert 'ID' to the index if you haven't already
    if 'ID' in df.columns:
        df.set_index('ID', inplace=True)

    # Map the docking_scores dict onto the DataFrame index
    df['DockingScore'] = df.index.map(docking_scores)

    # Reset the index if you prefer to keep 'ID' as a column
    df.reset_index(inplace=True)

    # Save the updated DataFrame back to the CSV
    df.to_csv(csv_file, index=False)
    if verbose:
        print(f"Updated CSV saved to {csv_file}")

    # Optionally generate images for top hits
    if generate_images:
        df_with_scores = df.dropna(subset=['DockingScore'])
        df_sorted = df_with_scores.sort_values(by='DockingScore', ascending=True)
        top_hits = df_sorted.head(top_n)

        if not top_hits.empty:
            os.makedirs(image_dir, exist_ok=True)

            hit_dicts = []
            for _, row in top_hits.iterrows():
                hit_dicts.append({
                    'ID': row['ID'],
                    'SMILES': row['SMILES'],
                    'DockingScore': row['DockingScore'],
                    'image_dir': image_dir
                })

            img_futures = []
            image_paths = []
            with concurrent.futures.ProcessPoolExecutor() as executor:
                for hit in hit_dicts:
                    img_futures.append(executor.submit(generate_image_parallel, hit))

                for future in tqdm(
                    concurrent.futures.as_completed(img_futures),
                    total=len(img_futures),
                    desc="Generating images"
                ):
                    path = future.result()
                    if path:
                        image_paths.append(path)

            if image_paths:
                create_combined_figure(image_paths, figure_file, top_n, verbose=verbose)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Traverse a directory (and subdirectories) for PDBQT files, extract docking scores, "
            "and append them to an existing CSV (which must contain columns for 'ID' and 'SMILES'). "
            "Optionally generate images and a combined figure for the top N hits, with a progress bar."
        )
    )
    parser.add_argument(
        "-c", "--csv_file",
        type=str,
        required=True,
        help="Path to the existing CSV file (with 'ID' and 'SMILES' columns)."
    )
    parser.add_argument(
        "-d", "--directory",
        type=str,
        required=True,
        help="Path to the directory that will be searched recursively for *_out.pdbqt files."
    )
    parser.add_argument(
        "-g", "--generate_images",
        action="store_true",
        help="Enable image generation for the top hits."
    )
    parser.add_argument(
        "-n", "--top_n",
        type=int,
        default=15,
        help="Number of top hits to generate images for (default = 15)."
    )
    parser.add_argument(
        "-i", "--image_dir", 
        type=str,
        default="images",
        help="Directory where individual molecule images will be saved (default = 'images')."
    )
    parser.add_argument(
        "-f", "--figure_file",
        type=str,
        default="top_hits_figure.png",
        help="Filename for the combined figure of top hits (default = 'top_hits_figure.png')."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose mode for detailed output."
    )

    args = parser.parse_args()

    # Check inputs
    if not os.path.exists(args.csv_file):
        print(f"Error: CSV file '{args.csv_file}' not found.")
        sys.exit(1)

    if not os.path.isdir(args.directory):
        print(f"Error: Directory '{args.directory}' not found or not a directory.")
        sys.exit(1)

    summarize_docking_results(
        csv_file=args.csv_file,
        directory=args.directory,
        generate_images=args.generate_images,
        top_n=args.top_n,
        image_dir=args.image_dir,
        figure_file=args.figure_file,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
