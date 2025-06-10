import os
import argparse
import subprocess
from pathlib import Path


# --- Align Nanopore reads to a reference sequence using minimap2 and samtools ---

if __name__ == "__main__":
    # --- Set up argument parser for input and output folder paths ---
    parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
    parser.add_argument("input", type=str, help="Path to the folder storing the input bam files")
    parser.add_argument("output", type=str, help="Folder path to store the output")
    parser.add_argument("ref_path", type=str, help="Path to the reference sequence fasta")

    args = parser.parse_args()

    # --- Convert input strings to Path objects for safer path handling ---
    input_folder = Path(args.input)
    output_folder = Path(args.output)
    ref_path = Path(args.ref_path)

    # --- Validate that input folder exists ---
    if not input_folder.is_dir():
        raise FileNotFoundError(f"Input folder '{input_folder}' does not exist!")
    if not ref_path.is_file():
        raise FileNotFoundError(f"Reference file '{ref_path}' does not exist!")
    # --- Create the output directory if it doesn't already exist ---
    output_folder.mkdir(parents=True, exist_ok=True)

    print(f"aligning {input_folder} files to {ref_path}...")

    # --- Generate a list of input file basenames without the .fastq.gz suffix ---
    input_files =  [f.name.replace(".fastq.gz","") for f in input_folder.glob("*.fastq.gz")]
    
    # --- Loop through each file and run alignment using minimap2 ---
    for input_file in input_files:
        # --- Define output path for the resulting BAM file ---
        output_file = output_folder / input_file  

        print(f"Aligning {input_file}.fastq.gz -> {output_file}.bam")

        # --- Construct the minimap2 command for Oxford Nanopore data ---
        #     - `--MD`: include the MD tag in the SAM output (used for mismatch information)
        #     - `-ax map-ont`: preset for Oxford Nanopore reads
        #     - Pipe output to `samtools view -bS` to convert SAM to BAM
        #     - Pipe to `samtools sort` to sort the BAM file before saving
        command = f"minimap2 --MD -ax map-ont {ref_path} {input_folder}/{input_file}.fastq.gz | samtools view -bS | samtools sort -o {output_file}.bam"

        # --- Run the full pipeline (alignment -> BAM conversion -> sorting) ---
        subprocess.run(command, shell=True, check=True)
        # --- Index the resulting BAM file so it can be used for downstream tools (e.g. IGV) ---
        subprocess.run(f"samtools index {output_file}.bam", shell=True, check=True)

    print("Processing complete!")
