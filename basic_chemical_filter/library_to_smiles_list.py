import csv
import argparse

def extract_smiles_and_id(input_csv, output_file):
    """
    Extracts SMILES and ID columns from the input CSV file and writes them
    to an output file formatted as 'smiles ZINCID'.
    
    Args:
        input_csv (str): Path to the input CSV file.
        output_file (str): Path to the output file.
    """
    try:
        with open(input_csv, 'r') as infile, open(output_file, 'w') as outfile:
            reader = csv.DictReader(infile)
            
            # Check if required columns exist
            if 'SMILES' not in reader.fieldnames or 'ID' not in reader.fieldnames:
                raise ValueError("Input file does not contain 'SMILES' and 'ID' columns.")
            
            # Write formatted output
            for row in reader:
                outfile.write(f"{row['SMILES']} {row['ID']}\n")
        
        print(f"File successfully written to: {output_file}")
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Extract SMILES and ZINCID from CSV file.")
    parser.add_argument("--input", required=True, help="Path to the input CSV file.")
    parser.add_argument("--output", required=True, help="Path to the output file.")
    args = parser.parse_args()

    extract_smiles_and_id(args.input, args.output)

if __name__ == "__main__":
    main()
