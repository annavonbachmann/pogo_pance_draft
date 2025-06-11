import pysam
import os
import argparse
from Bio import SeqIO
import pandas as pd
import csv
from scripts.Nanopore_functions import read_cleaning_

# --- Process Nanopore reads, mainly enforcing in-frame alignment of reads (necessary due to high Nanopore sequencing noise) ---

if __name__ == "__main__":
    # --- Set up command-line argument parser to accept input/output paths and reference file ---
    parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
    parser.add_argument("input", type=str, help="Path to the folder storing the input bam files")
    parser.add_argument("output", type=str, help="Folder path to store the output")
    parser.add_argument("ref_path", type=str, help="Path to the reference sequence fasta")

    args = parser.parse_args()

    # --- Assign argument values to local variables ---
    input_folder = args.input
    output_folder = args.output
    ref_path = args.ref_path

    # --- Set how many bases to remove from the beginning of each read ---
    cut_n_bases_from_start = 48 

    # --- Validate that the input folder exists ---
    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"Input folder '{input_folder}' does not exist!")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  # Create output folder if missing

    # --- Check that the reference file exists before proceeding ---
    if not os.path.exists(ref_path):
                    print(f"Reference file {ref_path} does not exist, please check the reference file path")
                    exit()


    # --- Read the reference sequence using Biopython and store it as a string ---
    ref = str(SeqIO.read(ref_path, "fasta").seq)

    # --- Call a custom function to clean and process the reads ---
    #     This function:
    #     - Loads BAM files from input_folder
    #     - Aligns reads to the given reference
    #     - Trims the first N bases
    #     - Returns: cleaned read sequences, a DataFrame of indels, and base quality info
    all_reads, indels, all_qualitities = read_cleaning_(input_folder, ref, cut_n_bases_from_start)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # --- Save the list of cleaned read sequences to a CSV file ---
    with open(f"{output_folder}/cleaned_reads.csv", "w", newline="") as f:
        writer = csv.writer(f)
        for item in all_reads:
            writer.writerow([item]) 

    print(f"Saved cleaned reads to {output_folder}/cleaned_reads.csv")

    # --- Save the reference sequence to a separate CSV file (single row) ---
    with open(f"{output_folder}/ref.csv", "w", newline="") as f:
        writer = csv.writer(f)    
        writer.writerow([ref]) 

    print(f"Saved reference sequence to {output_folder}/ref.csv")

    indels.to_csv(f"{output_folder}/indels.csv") 

    print(f"Saved indels to {output_folder}/indels.csv")
    # --- Save the reference sequence to a separate CSV file (single row) ---
    with open(f"{output_folder}/cleaned_reads_base_qualitities.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(all_qualitities)

    print("Done!")

