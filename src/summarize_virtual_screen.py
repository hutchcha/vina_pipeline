import os
import argparse
import sys
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt

def extract_zinc_id(filename):
    """Extracts the ZINC ID from the filename."""
    # Remove the _out.pdbqt suffix and return the part before it
    zinc_id = filename.split('_out.pdbqt')[0]
    return zinc_id

def load_smiles_file(smiles_filename, verbose=False):
    """Loads the ZINC ID to SMILES mapping from the separate SMILES file."""
    smiles_dict = {}
    with open(smiles_filename, 'r') as file:
        for line_num, line in enumerate(file, 1):
            # Assuming the SMILES file has format: SMILES ZINCID
            parts = line.strip().split()
            if len(parts) == 3:
                tranche, zinc_id, smiles = parts
                smiles_dict[zinc_id] = smiles
                if verbose:
                    print(f"Loaded SMILES for {zinc_id}: {smiles}")
            else:
                print(f"Warning: Line {line_num} in SMILES file is malformed: '{line.strip()}'")
    return smiles_dict

def extract_docking_score(pdbqt_file, verbose=False):
    """Extracts the docking score from the pdbqt file."""
    score = None
    with open(pdbqt_file, 'r') as file:
        for line in file:
            if line.startswith("REMARK VINA RESULT"):
                try:
                    score = float(line.split()[3])  # Extracting the energy score
                    if verbose:
                        print(f"Found docking score {score} in file {pdbqt_file}")
                except (IndexError, ValueError):
                    print(f"Warning: Could not parse docking score in file {pdbqt_file}")
                break
    return score

def smiles_to_2d_structure(smiles, output_image_name, label_text=None, score=None):
    """Converts a SMILES string to a 2D molecular image with an optional label."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES string: {smiles}")
        return False

    # Recalculate aromaticity
    Chem.SanitizeMol(mol)

    # Create the 2D structure image
    img = Draw.MolToImage(mol, size=(550, 500))
    
    # Convert the RDKit image (PIL Image object) to allow drawing text on it
    img = img.convert("RGBA")
    img = img.crop(img.getbbox())
    # Create Pillow drawing object
    draw = ImageDraw.Draw(img)
    
    # Choose the Arial font or default if unavailable
    try:
        # Update the font path as per your system
        font = ImageFont.truetype("/home/chutchins@uthouston.edu/fonts/Arial.TTF", size=50)
    except IOError:
        print("Arial font not found, using default")
        font = ImageFont.load_default()
    
    if label_text:
        # Get the bounding box of the text and calculate its width and height
        text_bbox = draw.textbbox((0, 0), label_text, font=font)
        text_width = text_bbox[2] - text_bbox[0]
        text_height = text_bbox[3] - text_bbox[1]
        
        # Position the text slightly above the bottom of the image
        position = ((img.width - text_width - 110) // 2, img.height - text_height - 80)
        draw.text(position, f"    {label_text} \n     {score} kcal/mol", fill="black", font=font)

    img.save(output_image_name)
    print(f"Image saved to {output_image_name}")
    return True

def summarize_docking_results(directory, smiles_filename, output_summary, generate_images=False, top_n=15, image_dir="images", figure_file="top_hits_figure.png", verbose=False):
    """Summarizes the docking results, including ZINC ID, score, and SMILES formula."""
    if verbose:
        print(f"Loading SMILES from {smiles_filename}...")
    smiles_dict = load_smiles_file(smiles_filename, verbose=verbose)
    
    results = []  # List to store tuples of (ZINC_ID, Docking_Score, SMILES)

    if verbose:
        print(f"Processing docking output files in directory: {directory}")

    for filename in os.listdir(directory):
        if filename.endswith("_out.pdbqt"):
            pdbqt_file = os.path.join(directory, filename)
            zinc_id = extract_zinc_id(filename)
            docking_score = extract_docking_score(pdbqt_file, verbose=verbose)
            
            # Find the SMILES for the ZINC ID
            smiles = smiles_dict.get(zinc_id, "N/A")
            
            # Append the result to the list
            if docking_score is not None:
                results.append((zinc_id, docking_score, smiles))
                if verbose:
                    print(f"Added {zinc_id}: Score={docking_score}, SMILES={smiles}")
            else:
                print(f"Warning: No docking score found for {zinc_id} in file {filename}")

    if verbose:
        print("Sorting results by Docking_Score in ascending order (lowest scores first)...")

    # Sort the results by Docking_Score (ascending order)
    results.sort(key=lambda x: x[1])

    if verbose:
        print(f"Writing sorted results to {output_summary}...")

    # Write the sorted results to the summary file
    with open(output_summary, 'w') as summary_file:
        summary_file.write("ZINC_ID,Docking_Score (kcal/mol),SMILES\n")
        for zinc_id, docking_score, smiles in results:
            summary_file.write(f"{zinc_id},{docking_score},{smiles}\n")

    print(f"Summary saved to {output_summary} with {len(results)} entries.")

    # If image generation is requested, generate images for the top N hits
    if generate_images:
        if verbose:
            print(f"Generating images for the top {top_n} hits...")
        
        # Ensure the image directory exists
        os.makedirs(image_dir, exist_ok=True)
        
        top_hits = results[:top_n]
        image_paths = []

        for idx, (zinc_id, docking_score, smiles) in enumerate(top_hits, 1):
            image_filename = f"{zinc_id}.png"
            image_path = os.path.join(image_dir, image_filename)
            success = smiles_to_2d_structure(smiles, image_path, label_text=zinc_id, score=docking_score)
            if success:
                image_paths.append(image_path)
            else:
                print(f"Failed to generate image for {zinc_id}")

        if image_paths:
            if verbose:
                print(f"Creating combined figure for the top {len(image_paths)} hits...")
            create_combined_figure(image_paths, figure_file, top_n, verbose=verbose)
            print(f"Combined figure saved to {figure_file}")
        else:
            print("No images were generated.")

def create_combined_figure(image_paths, figure_file, top_n, verbose=False):
    """Creates a combined figure of individual molecule images arranged in a grid."""
    import math

    # Determine grid size (e.g., 3x5 for 15 images)
    cols = min(5, top_n)
    rows = math.ceil(top_n / cols)

    # Create a matplotlib figure
    fig, axes = plt.subplots(rows, cols, figsize=(4, 1.875))
    axes = axes.flatten()  # Flatten in case of multiple rows

    for ax, img_path in zip(axes, image_paths):
        # Open the image and display it
        img = Image.open(img_path)
        ax.imshow(img)
        ax.axis('off')  # Hide axes

    # If there are unused subplots, hide them
    for ax in axes[len(image_paths):]:
        ax.axis('off')

    plt.tight_layout()
    plt.savefig(figure_file, dpi=1200)
    plt.close()

def main():
    # Initialize the parser
    parser = argparse.ArgumentParser(
        description="Extract docking scores from PDBQT files, lookup corresponding SMILES, summarize the results, and optionally generate images for top hits."
    )
    
    # Add arguments
    parser.add_argument(
        "-d", "--directory",
        type=str,
        required=True,
        help="Path to the directory containing the docking output PDBQT files."
    )
    
    parser.add_argument(
        "-s", "--smiles_file",
        type=str,
        required=True,
        help="Path to the SMILES file containing ZINC IDs and their corresponding SMILES strings."
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="docking_summary.txt",
        help="Path to the output summary file. Default is 'docking_summary.txt'."
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
        help="Number of top hits to generate images for. Default is 15."
    )
    
    parser.add_argument(
        "-i", "--image_dir",
        type=str,
        default="images",
        help="Directory to save individual molecule images. Default is 'images'."
    )
    
    parser.add_argument(
        "-f", "--figure_file",
        type=str,
        default="top_hits_figure.png",
        help="Filename for the combined figure of top hits. Default is 'top_hits_figure.png'."
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose mode for detailed output."
    )
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.isdir(args.directory):
        print(f"Error: The directory '{args.directory}' does not exist or is not a directory.")
        sys.exit(1)
    
    # Validate SMILES file
    if not os.path.isfile(args.smiles_file):
        print(f"Error: The SMILES file '{args.smiles_file}' does not exist or is not a file.")
        sys.exit(1)
    
    # Run the summarization
    summarize_docking_results(
        directory=args.directory,
        smiles_filename=args.smiles_file,
        output_summary=args.output,
        generate_images=args.generate_images,
        top_n=args.top_n,
        image_dir=args.image_dir,
        figure_file=args.figure_file,
        verbose=args.verbose
    )

if __name__ == "__main__":
    main()
