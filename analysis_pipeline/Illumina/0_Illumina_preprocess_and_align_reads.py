# --- Import necessary modules ---
import json
import sys
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from scripts.preprocessing_functions import *

# ============================== ARGUMENT PARSING ==============================


# --- Parse command line arguments for input folder and reference saving option ---

parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
parser.add_argument("filepath", type=str, help="Path to the fastq files")

parser.add_argument("--save_ref", action="store_true", help="Enable saving reference (default: False)") # Default `save_ref` to False, enable it only if `--save_ref` is provided


args = parser.parse_args()

base_dir = args.filepath
save_ref = args.save_ref

print(f"Filepath: {base_dir}")
print(f"Save reference: {save_ref}")


# ============================== LOAD CONFIGURATION ==============================

# --- Load configuration file containing processing parameters ---
if not os.path.exists(f"{base_dir}/config.json"):
    raise FileNotFoundError(f"No config file found. Please provide a config file in the directory.")
    
with open(f"{base_dir}/config.json", "r") as file:
    config = json.load(file)

# --- Extract parameters from config ---
catch_left = config["catch_left"]
catch_right = config["catch_right"]
remove_read_qualities = config["remove_read_qualities"]
Barcodes = config["Barcodes"]
Primer_seq = config["Primer_seq"] 
Primer_out_of_frame = config["Primer_out_of_triplets"]
variant = config["variant"]
cut_primer_start = config["cut_primer_start"]
cut_BC_seq = config["cut_BC_seq"]
used_Barcodes = config["used_Barcodes"]
Sections = config["Sections"]
amplicon = config["amplicon"]   
include_only_complete_reads = config["include_only_complete_reads"] if "include_only_complete_reads" in config else False
cutoff_a_read = config["cutoff_a_read"] if "cutoff_a_read" in config else False
cutoff_b_read = config["cutoff_b_read"] if "cutoff_b_read" in config else False

# --- Quality score mapping for filtering low-quality reads ---
quality_score = {
  '!':0, '"':1, '#':2, '$':3, '%':4, '&':5, "'":6, '(':7, ')':8, '*':9,
  '+':10, ',':11, '-':12, '.':13, '/':14, '0':15, '1':16, '2':17, '3':18, '4':19,
  '5':20, '6':21, '7':22, '8':23, '9':24, ':':25, ';':26, '<':27, '=':28, '>':29,
  '?':30, '@':31, 'A':32, 'B':33, 'C':34, 'D':35, 'E':36, 'F':37, 'G':38, 'H':39, 'I':40
}

# ============================== READ AND DEMULTIPLEX ==============================


# --- Read fastq files and filter reads based on low-quality scores ---
# Returns forward (a-) and reverse (b-) reads along with read IDs and (optionally) quality scores
a_seq, b_seq, _, _, a_ids, b_ids = read_sequences(variant = variant, cutoff_a_read = cutoff_a_read, cutoff_b_read = cutoff_b_read, catch_left=catch_left, catch_right=catch_right, return_qualities_ids=True, quality_score=remove_read_qualities, base_dir = base_dir)

# --- Assign reads to barcode groups and align with expected primers (allowing mismatches) ---
all_reads, all_ids = demultiplex_reads(a_seq, b_seq,Barcodes=Barcodes, Primer_seq=Primer_seq, used_Barcodes = used_Barcodes, Sections = Sections, max_mismatch_primerseq = 5, a_ids=a_ids, b_ids=b_ids, Primer_out_of_frame= Primer_out_of_frame, cut_primer_start=True, cut_BC_seq=True, catch_left=catch_left, catch_right=catch_right, include_only_complete_reads=include_only_complete_reads)


# ============================== SAVE READS AND REFERENCES ==============================

# --- Create output directories for preprocessed reads and reference sequences ---
# R2 reads are saved as reverse complemented reads
Path(f"{base_dir}/preprocessed/").mkdir(parents = True, exist_ok=True)
Path(f"{base_dir}/references/").mkdir(parents = True, exist_ok=True)

# --- For each barcode and section, save the reads and optionally the reference sequences ---

for Bc in used_Barcodes: 

    for section in Sections:
        # Generate and store reference sequence if requested
        if save_ref: # Only if save_ref is True, the reference sequences are saved. Otherwise, the reference sequences have to be provided prior to running the script in the references folder, with the correct file names (e.g. {variant}_{Bc}_{section}_Nt_filt_ref.fasta)
            ref = find_reference_seq(ref_gene=amplicon, Primer_seq=Primer_seq, Section=section, Primer_out_of_frame=Primer_out_of_frame) 
            ref_sequences = [SeqRecord(Seq(ref), id = f"{variant}_{section}_ref", description = f"{variant} {section} DNA sequence")]

        for Read_dir in ["R1", "R2"]:
            # Get reads for this barcode/section and reverse-complement R2 reads
            reads = all_reads[f"{Bc}_{section}_{Read_dir}"] if Read_dir == "R1" else [dna_rev_comp(r) for r in all_reads[f"{Bc}_{section}_{Read_dir}"]]

            output_file = f"{base_dir}/preprocessed/{variant}_{Bc}_{section}_Nt_filt_{Read_dir}.fasta"
            sequences = [SeqIO.SeqRecord(Seq(read), id = all_ids[f"{Bc}_{section}_{Read_dir}"][i], description = f"{variant} {Bc} DNA sequence") for i, read in enumerate(reads)]

            # Write sequences to FASTA
            count = len(sequences)
            with open(output_file, "w") as output_handle:
                SeqIO.write(sequences, output_handle, "fasta")
            print("Saved %i records to %s" % (count, output_file))
          
           # Write reference FASTA if save_ref enabled
            if save_ref:
                with open(f"{base_dir}/references/{variant}_{Bc}_{section}_Nt_filt_ref.fasta", "w") as output_handle:
                    SeqIO.write(ref_sequences, output_handle, "fasta")
                print(f"Saved reference sequence for {variant} {section} to {base_dir}/preprocessed/{variant}_{Bc}_{section}_Nt_filt_ref.fasta")



# ============================== RUN BLAST ALIGNMENT ==============================
                
# --- Define key directories for BLAST ---
input_dir =  Path(f"{base_dir}/preprocessed/")      # Contains preprocessed reads
reference_dir = Path(f"{base_dir}/references/")     # Contains reference sequences
output_dir =  Path(f"{base_dir}/blast/alignments/") # Stores output alignment files
blast_db_dir = Path(f"{base_dir}/blast/db")         # Stores BLAST-formatted databases

# --- Ensure necessary directories exist ---
output_dir.mkdir(parents=True, exist_ok=True)
blast_db_dir.mkdir(parents=True, exist_ok=True)

# --- Create BLAST databases for each read file (self-alignment approach) ---
for read_file in input_dir.glob("*.fasta"):

    read_basename = read_file.stem 
    db_path = blast_db_dir / read_basename

    print(f"Creating BLAST database for {read_file}...")
    subprocess.run([
        "makeblastdb",
        "-in", str(read_file),
        "-dbtype", "nucl", 
        "-out", str(db_path)
    ], check=True)

# --- Run BLAST alignment: reference vs. processed reads ---
for read_file in input_dir.glob("*.fasta"):

    read_basename = read_file.stem
    ref_file = reference_dir / f"{read_basename[:-2]}ref.fasta"
    read_basename = read_file.stem
    db_path = blast_db_dir / read_basename
    output_file = output_dir / f"{read_basename}.out"

    print(f"Running BLAST: aligning reads from {read_file} to {ref_file}...")
    subprocess.run([
        "blastn",  
        "-query", str(ref_file),  # Reference as query
        "-db", str(db_path),      # Read database as target
        "-out", str(output_file),
        "-outfmt", "15", # JSON output format
        "-max_target_seqs", "100000",  # Limit to 100k alignments (adjust if needed)
    ], check=True)

print("BLAST pipeline completed!")
