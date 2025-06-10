# --- Import necessary modules ---
import os
from Bio.SeqIO import QualityIO
import numpy as np
from scripts.utils import *
import pandas as pd

# --- Read and quality-filter forward (R1) and reverse (R2) reads from FASTQ files ---
def read_sequences(variant, 
                   catch_left, 
                   catch_right, 
                   base_dir = None, 
                   cutoff_a_read = None, 
                   cutoff_b_read = None, 
                   quality_score = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*','+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5'], 
                   return_qualities_ids = False):
    """
    Reads and quality-filters forward (R1) and reverse (R2) reads from FASTQ files.

    Parameters:
    - variant: Identifier of the FASTQ files (str), expects files named {variant}_R1_001.fastq and {variant}_R2_001.fastq
    - catch_left: Sequence indicating the end of a region to be removed from the start of R1 (str)
    - catch_right: Sequence indicating the end of a region to be removed from the start of R2 (str)
    - base_dir: Path to the directory containing the FASTQ files (str, default: ./data/fastq/)
    - cutoff_a_read: Max length of R1 reads after filtering (int or None)
    - cutoff_b_read: Max length of R2 reads after filtering (int or None)
    - quality_score: List of low-quality symbols (str) used to truncate reads (default: corresponding to >1% error rate)
    - return_qualities_ids: If True, returns read qualities and IDs (bool)

    Returns:
    - a_sequences: Filtered R1 sequences (list of str)
    - b_sequences: Filtered R2 sequences (list of str)
    - (optional) a_qualities, b_qualities, a_ids, b_ids: Corresponding qualities and identifiers (list of str)
    """

    if not base_dir:
        base_dir = os.getcwd() + "/data/fastq/"

    # --- Initialize lists to store sequences, qualities, and IDs ---
    a_sequences, b_sequences = [], []
    a_qualities, b_qualities = [], []
    a_ids, b_ids = [], []

    # --- Open R1 and R2 FASTQ files ---
    with open(f'{base_dir}/{variant}_R1_001.fastq', "rt") as a_file, open(f'{base_dir}/{variant}_R2_001.fastq', "rt") as b_file:

        a_reader = QualityIO.FastqGeneralIterator(a_file)
        b_reader = QualityIO.FastqGeneralIterator(b_file)
        # --- Loop through paired reads ---
        for total_read, (a, b) in enumerate(zip(a_reader, b_reader)):
                
                a_id, a_seq, a_qual = a
                b_id, b_seq, b_qual = b
                # --- Get the first position with low quality ---
                cutoff_a = find_(a_qual, quality_score)
                cutoff_b = find_(b_qual, quality_score)

                # --- Optionally crop after catch sequences ---
                if cutoff_a_read and catch_left in a_seq:
                    if cutoff_a > (a_seq.index(catch_left) + cutoff_a_read):
                        cutoff_a = a_seq.index(catch_left)  + len(catch_left) + cutoff_a_read 
                
                if cutoff_b_read and dna_rev_comp(catch_right) in b_seq: 
                    if cutoff_b > (b_seq.index(dna_rev_comp(catch_right)) + cutoff_b_read):
                        cutoff_b = b_seq.index(dna_rev_comp(catch_right))+ len(catch_right) + cutoff_b_read
                
                # --- Append truncated reads and metadata ---
                a_sequences.append(a_seq[:cutoff_a])
                a_qualities.append(a_qual[:cutoff_a])
                b_sequences.append(b_seq[:cutoff_b])
                b_qualities.append(b_qual[:cutoff_b])
                a_ids.append(a_id)
                b_ids.append(b_id)
                
        print("total reads", total_read+1)

    # --- Return sequences, optionally with qualities and IDs ---
    if return_qualities_ids:
        
        return a_sequences, b_sequences, a_qualities, b_qualities, a_ids, b_ids 
     
    else: 

        return a_sequences, b_sequences

# --- Split pooled R1/R2 reads by barcode and primer section and extract regions of interest ---
def demultiplex_reads(a_seqs:list, 
                      b_seqs:list,
                      Barcodes:dict, 
                      Primer_seq:dict, 
                      Primer_out_of_frame:dict,
                      used_Barcodes:list, 
                      Sections:list, 
                      max_mismatch_primerseq:int = 5, 
                      a_ids:list = None, 
                      b_ids:list = None, 
                      cut_BC_seq = True,
                      cut_primer_start=True,
                      catch_left = "",
                      catch_right = "",
                      include_only_complete_reads = False):
    """
    Splits pooled R1/R2 reads by barcode and primer section and extracts regions of interest.

    Parameters:
    - a_seqs, b_seqs: Forward (R1) and reverse (R2) reads (list of str)
    - Barcodes: Dictionary of barcode sequences (dict), keys like 'BC1_fwd', 'BC1_rev'
    - Primer_seq: Dictionary of primer sequences for each section (dict), keys like 'S1_fwd', 'S1_rev'
    - Primer_out_of_frame: Dictionary specifying frame offset before codon start (dict)
    - used_Barcodes: List of barcode identifiers used in sequencing (list of str)
    - Sections: List of sequenced section identifiers (list of str)
    - max_mismatch_primerseq: Allowed mismatches in primer sequence (default: 5)
    - a_ids, b_ids: Optional list of sequence IDs for R1 and R2 (list of str)
    - cut_BC_seq: Whether to remove barcode and primer start (bool)
    - cut_primer_start: Whether to remove nucleotides before triplet start (bool)
    - catch_left, catch_right: Sequences used to crop start and end (str)
    - include_only_complete_reads: If True, include only reads with both cut sites (bool)

    Returns:
    - read_Dict: Dictionary with demultiplexed reads by barcode and section
    - ids_Dict: (optional) Dictionary with corresponding read IDs
    """

    read_Dict = {}
    ids_Dict = {}

    # --- Loop through each barcode and section to group reads ---
    for Barcode in used_Barcodes: 

        for Section in Sections:

            # --- Construct combined barcode+primer sequences ---
            fwd_BC_Primer_seq = Barcodes[Barcode + "_fwd"] + Primer_seq[Section+"_fwd"]
            rev_BC_Primer_seq = Barcodes[Barcode + "_rev"] + Primer_seq[Section+"_rev"] 

            ### select the reads that contain the forward and reverse BC + primer sequences, thereby allowing for n mismatches in the primer sequences but no errors in BCs
            fwd_idxs = []
            rev_idxs = []
            # --- Identify forward reads matching barcode and allowing mismatches in primers ---
            for a_idx, seq in enumerate(a_seqs):
                a_mismatch_to_primer_seq = sum([sequence!=primer_ref for sequence, primer_ref in zip(seq[len(Barcodes[Barcode + "_fwd"]):len(fwd_BC_Primer_seq)], Primer_seq[Section+"_fwd"])])
                if seq[:len(Barcodes[Barcode + "_fwd"])] == Barcodes[Barcode + "_fwd"] and a_mismatch_to_primer_seq <= max_mismatch_primerseq:
                    fwd_idxs.append(a_idx)
            # --- Identify reverse reads matching barcode and primer ---
            for b_idx, seq in enumerate(b_seqs):
                b_mismatch_to_primer_seq = sum([sequence!=primer_ref for sequence, primer_ref in zip(seq[len(Barcodes[Barcode + "_rev"]):len(rev_BC_Primer_seq)], Primer_seq[Section+"_rev"])])
                if seq[:len(Barcodes[Barcode + "_rev"])] == Barcodes[Barcode + "_rev"] and b_mismatch_to_primer_seq <= max_mismatch_primerseq:
                    rev_idxs.append(b_idx)
            # --- Keep only matching reads that appear in both forward and reverse lists ---
            indexes = set(
                [idx for idx in fwd_idxs if b_seqs[idx][:len(Barcodes[Barcode + "_rev"])] == Barcodes[Barcode + "_rev"]]  +  
                [idx for idx in rev_idxs if a_seqs[idx][:len(Barcodes[Barcode + "_fwd"])] == Barcodes[Barcode + "_fwd"]]) # Only keep reads that match in the fwd and rev BC seqs
                
            print(sum([len(b_seqs[fwd_i]) < len(Barcodes[Barcode + "_rev"]) for fwd_i in fwd_idxs]), "b reads are empty") # Reads that are only in the reverse list
            print(sum([len(a_seqs[rev_i]) < len(Barcodes[Barcode + "_fwd"]) for rev_i in rev_idxs]), "a reads are empty") # Reads that are only in the reverse list
            print(len(indexes), "reads with matching BC and primer seq")
            print(len(set(fwd_idxs+ rev_idxs)) - len(indexes), "reads with index swapping")
            
            # --- Extract reads and IDs ---
            a_seq_Bc_Sec = [a_seqs[i] for i in indexes]
            b_seq_Bc_Sec = [b_seqs[i] for i in indexes]

            print(Barcode, Section, len(a_seq_Bc_Sec), "total fwd reads")

            if a_ids and b_ids:
                a_ids_Bc_Sec = [a_ids[i].split(" ")[0] for i in indexes]
                b_ids_Bc_Sec = [b_ids[i].split(" ")[0]  for i in indexes]

            # --- Optionally cut off barcodes and primer start ---
            if cut_BC_seq: 
                cutoff_a = len(Barcodes[Barcode + "_fwd"]) if not cut_primer_start else len(Barcodes[Barcode + "_fwd"]) + Primer_out_of_frame[Section + "_fwd"]
                cutoff_b = len(Barcodes[Barcode + "_rev"]) if not cut_primer_start else len(Barcodes[Barcode + "_rev"]) + Primer_out_of_frame[Section + "_rev"]

                a_seq_Bc_Sec = [a[cutoff_a:] if len(a)>=cutoff_a else "" for a in a_seq_Bc_Sec]
                b_seq_Bc_Sec = [b[cutoff_b:] if len(b)>=cutoff_b else "" for b in b_seq_Bc_Sec]
            
            # --- Trim reads using catch sequences ---
            if include_only_complete_reads: # Only include reads that contain the full sequence (i.e. catch_left **and** catch_right is present)
                a_seq_Bc_Sec = [read[read.index(catch_left)+len(catch_left):read.index(catch_right)] if catch_left in read and catch_right in read else "" for read in a_seq_Bc_Sec ]

                b_seq_Bc_Sec = [read[read.index(dna_rev_comp(catch_right))+len(catch_right):read.index(dna_rev_comp(catch_left))] if dna_rev_comp(catch_right) in read  and dna_rev_comp(catch_left) in read else "" for read in b_seq_Bc_Sec ]

            else: ## cut sequences at the catch_left and catch_right positions, reads do not have to be complete 
                a_seq_Bc_Sec = [a[a.index(catch_left)+len(catch_left):] if catch_left in a else "" for a in a_seq_Bc_Sec]
                b_seq_Bc_Sec = [b[b.index(dna_rev_comp(catch_right))+len(catch_right):] if dna_rev_comp(catch_right) in b else "" for b in b_seq_Bc_Sec]

            # --- Store reads and optionally IDs ---
            read_Dict[f"{Barcode}_{Section}_R1"] = a_seq_Bc_Sec
            read_Dict[f"{Barcode}_{Section}_R2"] = b_seq_Bc_Sec

            if a_ids and b_ids:
                ids_Dict[f"{Barcode}_{Section}_R1"] = a_ids_Bc_Sec
                ids_Dict[f"{Barcode}_{Section}_R2"] = b_ids_Bc_Sec

            print(f"################# Completed {Barcode} {Section} #################")

        print(f"################# Completed {Barcode} #################")
    if a_ids and b_ids:
        return  read_Dict, ids_Dict
    else:
        return read_Dict 


# --- Extract a section of the reference gene using primer sequences and frame correction ---
def find_reference_seq(ref_gene, 
                       Primer_seq, 
                       Section, 
                       Primer_out_of_frame):
    """
    Extracts a section of the reference gene using primer sequences and frame correction.

    Parameters:
    - ref_gene: Full reference DNA sequence (str)
    - Primer_seq: Dictionary with forward and reverse primers per section (dict)
    - Section: Section name to extract (str)
    - Primer_out_of_frame: Dictionary specifying nt offset before codon start (dict)

    Returns:
    - ref_gene_section: Trimmed reference gene section for the specified primers (str)
    """
    # --- Apply offset to primer sequences to get in-frame triplet start ---
    tripl_st = Primer_out_of_frame[Section+"_fwd"]
    tripl_end = Primer_out_of_frame[Section+"_rev"]
    primer_fwd = Primer_seq[Section + "_fwd"][tripl_st:]
    primer_rev = dna_rev_comp(Primer_seq[Section+"_rev"][tripl_end:])
    # --- Find matching region in reference gene ---
    ref_gene_section = ref_gene[ref_gene.index(primer_fwd):ref_gene.index(primer_rev)+len(primer_rev)]

    return ref_gene_section

