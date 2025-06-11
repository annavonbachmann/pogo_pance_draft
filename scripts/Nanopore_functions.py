# --- Import necessary modules ---
import pandas as pd 
import pysam
import os
from Bio import SeqIO
from scripts.utils import translate_dna2aa
from scripts.functions_ import mask_ref_in_variants_df
import numpy as np
import re


# --- Characterizes variants and indels in Nanopore reads for DMS analysis ---
def characterize_DMS_Nanopore(aligned_reads, ref, data_type = "AA"):
    """
    Characterizes DMS alignments by counting insertions, deletions, and substitutions per position.

    Parameters:
    - aligned_reads: List of aligned read sequences (list of str)
    - ref: Reference DNA sequence (str)
    - data_type: Output format: "AA" for amino acids, "DNA", or "Codons" (str)
    - read_dir: string, "R1" for forward or "R2" for reverse
    - cut_to_same_start: bool, sequences cut at same start position if True, otherwise provide start/end positions as query_form and query_to

    Returns:
    - all_variants: Variant counts per position (dict)
    - enrichment_counts: Variant counts per position with reference masked (pd.DataFrame)
    - enrichment_relative: Relative frequencies of variants with reference masked (pd.DataFrame)
    - indels_freq: Frequency of insertions and deletions per position (pd.DataFrame)
    """

    all_variants = {}
    seq_with_off_target_indels = 0
    included_seq = 0
    indels = pd.DataFrame(columns = range(len(ref)), index = ["insertion"], data = 0)

    # --- Characterizes variants and indels in Nanopore reads for DMS analysis ---
    reference = translate_dna2aa(ref) if data_type == "AA" else ref
    if data_type == "Codons": 
        reference = [reference[i:i+3] for i in range(0, len(reference)//3*3, 3)]

    # --- Initialize variant tracking dictionary ---
    if data_type == "AA":
        for idx in range(len(reference)):
            all_variants[idx] = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 
                                    'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 
                                    'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 
                                    'W':0, 'Y':0, '*':0, 'X':0} ## X: are triplets with missing nucleotides (i.e. in the aligned seq, "-" is present), that would lead to a frameshift
    elif data_type == "DNA": 
        for idx in range(len(reference)):
            all_variants[idx] = {'A':0, 'C':0, 'G':0, 'T':0, "-":0}
    
    else:
        codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT', 'X']
        for idx in range(len(reference)): 
            all_variants[idx] = {codon: 0 for codon in codons}


    # --- Analyze each alignment for indels and substitutions ---
    for alignment in aligned_reads:
        
        if "-" in alignment: 
            seq_with_off_target_indels += 1
            shift = 0 # Count the shift of the position compared to the reference that occurs if there is an insertion in the qseq 

            for idx,nt in enumerate(alignment): 
                pos = idx - shift # Adjust for the shift in the index, due to prior insertions

                if nt == "-":
                    indels.loc["insertion", pos] += 1
                    shift += 1 # Correct for the shift in the index due to the insertion

        if data_type == "Codons":
            alignment = [alignment[i:i+3] for i in range(0, len(alignment)//3*3, 3)]

        elif data_type == "AA":
            alignment = translate_dna2aa(alignment)
        
        for idx, variant in enumerate(alignment): 
            if "-" in variant:
                variant = "X"
            all_variants[idx][variant] += 1
       
    # --- Normalize and mask reference variants ---
    indels_freq = indels/len(aligned_reads)
    print(seq_with_off_target_indels, "sequences have off target indels")
    print(len(aligned_reads), "sequences are included in the enrichment analysis")

    all_variants = pd.DataFrame.from_dict(all_variants)

    enrichment_counts, enrichment_relative = mask_ref_in_variants_df(variant_df=all_variants, ref_seq=reference if data_type=="AA" else ref, data_type=data_type)
    

    return all_variants, enrichment_counts,enrichment_relative, indels_freq


# --- Processes BAM files to clean reads and align them to reference with frame correction ---
def read_cleaning_(input_folder, ref, cut_n_bases_from_start=48):
    """
    Cleans and aligns Nanopore reads to force frame consistency and correct for indels.

    Parameters:
    - input_folder: Path to folder containing .bam files (str)
    - ref: Reference DNA sequence (str)
    - cut_n_bases_from_start: Number of bases to trim from start for frame alignment (int)

    Returns:
    - all_reads: List of processed read sequences (list of str)
    - indels: Insertions and deletions per reference position (pd.DataFrame)
    - all_qualities: Base quality scores per read (list of lists)
    """
    # --- Collect all BAM files in the input folder ---
    bam_files = [f for f in os.listdir(input_folder) if f.endswith('.bam')]
    all_reads = []
    indels = pd.DataFrame(columns=range(len(ref)), index=["I", "D"], data=0)
    all_qualities = []

    # --- Define regex pattern to parse CIGAR strings ---
    cigar_pattern = re.compile(r"(\d+)([MIDNSHP=X])")  

    for file_nr, bamfile_name in enumerate(bam_files):
        bam_path = os.path.join(input_folder, bamfile_name)
        bamfile = pysam.AlignmentFile(bam_path, "rb")

        print("Status:", file_nr + 1, "/", len(bam_files), "done")

        for read in bamfile.fetch():
            if read.is_unmapped or read.query_sequence is None:
                print(f"Skipping read {read.query_name}")
                continue
            # --- Initialize alignment variables ---
            alignment_start = read.reference_start
            seq = read.query_sequence
            qualitities = read.query_qualities
            refined_qualities = []
            refined_seq_list = []
            ref_pos = 0 

            # --- Parse the CIGAR operations ---
            cigar_operations = [(int(length), op) for length, op in cigar_pattern.findall(read.cigarstring)]

            query_pos = 0  # Track position in the read sequence
            ref_pos = alignment_start  # Start position in reference

            for length, operation in cigar_operations:
                if operation == "M":  # Matches (or mismatches)
                    refined_seq_list.extend(seq[query_pos:query_pos + length])
                    refined_qualities.extend(qualitities[query_pos:query_pos + length])
                    query_pos += length
                    ref_pos += length
                elif operation == "I":  # Insertion 
                    indels.loc["I", ref_pos] += 1
                    query_pos += length
                elif operation == "D":  # Deletion
                    refined_qualities.extend([""] * length)
                    refined_seq_list.extend("-" * length)
                    indels.loc["D",ref_pos] += 1
                    ref_pos += length
                elif operation == "N":  # Skipped region in reference
                    ref_pos += length
                elif operation == "S":  # Soft clipping (ignored bases at ends)
                    query_pos += length
                elif operation == "H":  # Hard clipping (ignored bases, not in read)
                    continue
                elif operation == "P":  # Padding (shouldn't appear in Nanopore data)
                    continue

            # --- Join sequence list into string ---
            refined_seq = "".join(refined_seq_list)

            # --- Trim to standardized start position if possible ---
            if alignment_start < cut_n_bases_from_start:
                cut_start = cut_n_bases_from_start - alignment_start
                refined_seq = refined_seq[cut_start:]
                refined_qualities = refined_qualities[cut_start:]

                all_reads.append(refined_seq)
                all_qualities.append(refined_qualities)

        print(f"Processed {bamfile_name}")

    print("Total reads:", len(all_reads))
    return all_reads, indels, all_qualities


# --- Extracts left and right linker regions from aligned Nanopore reads ---
def get_linker_regions(input_folder, ref, cut_site_seq_left, cut_site_seq_right, left_linker_region_len, right_linker_region_len):
    """ 
    Extracts linker regions from aligned Nanopore reads based on cut sites.

    Parameters:
    - input_folder: Path to folder with BAM files (str)
    - ref: Reference DNA sequence (str)
    - cut_site_seq_left: Sequence of left linker at start of insert (str)
    - cut_site_seq_right: Sequence of right linker at end of insert (str)
    - left_linker_region_len: Length of left linker region to extract (int)
    - right_linker_region_len: Length of right linker region to extract (int)

    Returns:
    - all_left_linkers: Left linker regions by read ID (dict)
    - all_right_linkers: Right linker regions by read ID (dict)
    """
    # --- List all BAM files in the input folder ---
    bam_files = [f for f in os.listdir(input_folder) if f.endswith('.bam')]
    # --- Initialize dictionaries to store linker regions ---
    all_left_linkers = {}
    all_right_linkers = {}
    indels = pd.DataFrame(columns=range(len(ref)), index=["I", "D"], data=0)
    
    # --- Track number of reads where linker could not be extracted ---
    left_linker_excluded = 0
    right_linker_excluded = 0
    # --- Compile regex pattern to parse CIGAR strings --- 
    cigar_pattern = re.compile(r"(\d+)([MIDNSHP=X])") 

    id_nr = 0
    # --- Loop through all BAM files ---
    for file_nr, bamfile_name in enumerate(bam_files):
        bam_path = os.path.join(input_folder, bamfile_name)
        bamfile = pysam.AlignmentFile(bam_path, "rb")

        print("Status:", file_nr + 1, "/", len(bam_files), "done")

        for read in bamfile.fetch():
            # --- Skip unmapped reads or reads without sequence ---
            if read.is_unmapped or read.query_sequence is None:
                print(f"Skipping read {read.query_name}")
                continue
            
            alignment_start = read.reference_start
            seq = read.query_sequence
            # --- Initialize aligned read and reference lists ---
            refined_seq_list = []
            refined_ref_list = []
            ref_pos = 0

            # --- Get CIGAR operations and initialize positions ---
            cigar_operations = [(int(length), op) for length, op in cigar_pattern.findall(read.cigarstring)]
            query_pos = 0  # Track position in the read sequence
            ref_pos = alignment_start  # Start position in reference

            # --- Simulate alignment by iterating through CIGAR operations ---
            for length, operation in cigar_operations:
                if operation == "M":  # Matches (or mismatches)
                    refined_seq_list.extend(seq[query_pos:query_pos + length])
                    refined_ref_list.extend(ref[ref_pos:ref_pos+length])
                    query_pos += length
                    ref_pos += length

                elif operation == "I":  # Insertion (extra bases in read)
                    indels.loc["I", ref_pos] += 1# length
                    refined_seq_list.extend(seq[query_pos:query_pos + length])
                    refined_ref_list.extend("-"*length)
                    query_pos += length
                elif operation == "D":  # Deletion (missing bases in read)
                    #refined_qualities.extend([""] * length)
                    refined_seq_list.extend("-" * length)
                    indels.loc["D",ref_pos] += 1#length
                    refined_ref_list.extend(ref[ref_pos:ref_pos+length])
                    ref_pos += length
                elif operation == "N":  # Skipped region in reference
                    ref_pos += length
                elif operation == "S":  # Soft clipping (ignored bases at ends)
                    query_pos += length
                elif operation == "H":  # Hard clipping (ignored bases, not in read)
                    continue
                elif operation == "P":  # Padding (shouldn't appear in Nanopore data)
                    continue
            # --- Build aligned sequences as strings ---
            refined_seq = "".join(refined_seq_list)
            refined_ref = "".join(refined_ref_list)
            # --- Find cut sites in the aligned reference ---
            cut_site_left = refined_ref.find(cut_site_seq_left) 
            cut_site_right = refined_ref.find(cut_site_seq_right) 

            # --- Extract linker regions only if both cut sites are found ---
            if cut_site_left != -1 and cut_site_right != 1:
                # --- Extract left linker sequence preceding the left cut site ---
                all_left_linkers["id"+str(id_nr)] = {"hseq" : refined_seq[cut_site_left-left_linker_region_len:cut_site_left], 
                                                     "qseq" : refined_ref[cut_site_left-left_linker_region_len:cut_site_left]}
                # --- Extract right linker sequence following the right cut site ---
                cut_site_right = cut_site_right + len(cut_site_seq_right)
                all_right_linkers["id"+str(id_nr)] = {"hseq": refined_seq[cut_site_right:cut_site_right+right_linker_region_len],
                                                 "qseq": refined_ref[cut_site_right:cut_site_right+right_linker_region_len]}
                
            else: 
                # --- Count how many reads were excluded due to missing cut sites ---
                left_linker_excluded +=1
                right_linker_excluded +=1

            id_nr += 1 # --- Update read ID counter ---
        print(f"Processed {bamfile_name}")

    print("Total left linkers:", sum([l != "" for l in all_left_linkers]))
    print("Total right linkers:", sum([l != "" for l in all_right_linkers]))
    print(left_linker_excluded, "left linkers are excluded")
    print(right_linker_excluded, "right linkers are excluded")
    return all_left_linkers, all_right_linkers



# --- Generates genotype dictionary from amino acid sequences of the reads ---
def get_genotype_dict_from_AAseqs(all_Aas, ref_AAseq, ref_aa_annot, not_masked_positions = None, combined = False): 
    """
    Extracts genotype information from amino acid (AA) sequences of aligned reads by comparing them to a reference.

    Parameters:
    - all_Aas: List of AA sequences from the reads (list of str)
    - ref_AAseq: Reference AA sequence (str)
    - ref_aa_annot: List of AA position annotations (e.g. 'A23') from the reference (list of str)
    - not_masked_positions: Positions to include in the comparison; others will be ignored (list of int or None)
    - combined: Whether to group mutations by position only (True) or include the substituted AA (False) (bool)

    Returns:
    - genotypes: Genotype labels (e.g. 'A23V_T45M' or 'A23_T45' if combined=True) with counts (dict)
    """
    # --- Initialize dictionary to store genotype counts ---
    genotypes = {}

    # --- Loop through each AA sequence from the reads ---
    for read in all_Aas:
        variant = []
        # --- Loop through each amino acid in the sequence ---
        for idx, Aa in enumerate(read): 
            # --- Skip masked positions if specified ---
            if not_masked_positions and idx not in not_masked_positions: ## if the position is masked, skip
                continue
            else:
                # --- If the amino acid is not a masked 'X' and differs from the reference ---
                if Aa != "X" and Aa != ref_AAseq[idx]: 
                    # --- Create mutation string; either with or without the changed AA, depending on 'combined' flag ---
                    observed_mut = f"{ref_aa_annot[idx]}{Aa}" if not combined else ref_aa_annot[idx]
                    variant.append(observed_mut)
            
        # --- Determine genotype label: 'WT' if no mutations, or joined mutation list ---
        if len(variant) == 0: 
            variant = "WT"
        else: 
            variant = "_".join(variant)

        # --- Count how many times each genotype occurs ---
        genotypes[variant] = genotypes.get(variant, 0) + 1
        
    return genotypes
        