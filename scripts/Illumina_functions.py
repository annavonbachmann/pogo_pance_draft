# --- Import necessary modules ---
import os
import sys
os.chdir("..")
import pandas as pd
import re
from scripts.utils import translate_dna2aa
from scripts.functions_ import mask_ref_in_variants_df
import numpy as np

def divide_alignments(blast_alignments, cut_site_seq, query_seq, read_dir="R1", cut_read_start= 12): 
    """ 
    Divides BLAST alignments into linker and insert (LOV2) regions based on a cut-site sequence.
    
    Parameters:
    - blast_alignments: list of BLAST alignments in JSON format
    - cut_site: DNA sequence used to detect where the LOV2 insert starts (if read_dir is R1) or ends (if read_dir is R2), string
    - read_dir: "R1" or "R2" to indicate read direction, string
    - query_seq: full reference DNA sequence, string
    - cut_read_start: how many bases to trim at read start to isolate the linker region, integer

    Returns:
    - linker_alignments: dictionary of linker sequences by read ID {seq_id : {"qseq": qseq, "hseq": hseq, "midline": midline}, ...}
    - LOV2_alignments: dictionary of insert sequences by read ID {seq_id : {"qseq": qseq, "hseq": hseq, "midline": midline}, ...}
    - coverage_linker : np.array, per-base coverage array across the linker region
    - coverage_LOV2 : np.array, per-base coverage array across the LOV2 region
    """

    linker_alignments = {}
    LOV2_alignments = {}
    LOV2_start_indel_count = 0
    coverages = np.zeros(len(query_seq), dtype = int)

    for alignment in blast_alignments:
        query_from = alignment["hsps"][0]["query_from"]-1  # Convert to 0-based index
        query_to = alignment["hsps"][0]["query_to"]  # Convert to 0-based index
        qseq = alignment["hsps"][0]["qseq"].upper()
        hseq = alignment["hsps"][0]["hseq"].upper()
        seq_id = alignment["description"][0]["title"]
        midline = alignment["hsps"][0]["midline"]

        cut_site = qseq.find(cut_site_seq) # Find LOV2 position, we add len(LOV_endseq) if "R2" later, so that we can first filter out reads that do not contain the sequence of interest (i.e. cut_site = -1) due to insertions at these sites 
        if cut_site != -1: # If -1, there are insertions in start of LOV2, thus sequence not in ref_seq and we do not include these sequences
            
            if read_dir=="R2":
                cut_site += len(cut_site_seq)

                read_out_of_frame = cut_site % 3 # Correct for out-of-frame reads (only for R2, since R1 is always in frame)
                
                linker_alignments[seq_id] = {"qseq": qseq[cut_site:-cut_read_start], "hseq": hseq[cut_site:-cut_read_start], "midline": midline[cut_site:cut_read_start]} # Always in frame since it starts at the cut site
                LOV2_alignments[seq_id] = {"qseq": qseq[read_out_of_frame:cut_site], "hseq": hseq[read_out_of_frame:cut_site], "midline": midline[read_out_of_frame:cut_site]}

            else:    
                
                linker_alignments[seq_id] = {"qseq": qseq[cut_read_start:cut_site], "hseq": hseq[cut_read_start:cut_site], "midline": midline[cut_read_start:cut_site]}
                LOV2_alignments[seq_id] = {"qseq": qseq[cut_site:], "hseq": hseq[cut_site:], "midline": midline[cut_site:]}
        else:
            LOV2_start_indel_count +=1
            continue

        coverages[query_from:query_to] += 1

    print(LOV2_start_indel_count, "Sequences are excluded, since LOV2 start site could not be found in the reference (due to '-' i.e. insertions at the start of LOV2)")

    return linker_alignments, LOV2_alignments, coverages


def restructure_alignments(blast_alignments, query_seq, read_dir = "R1"): 
    """ 
    Filters and formats BLAST alignments to ensure in-frame reads and calculates coverage.

    Returns:
    - alignments: dict with clean qseq/hseq/midline for in-frame reads
    - coverages: array of per-base read depth
    """
    exlude_seqs = 0
    alignments = {}
    coverages = np.zeros(len(query_seq), dtype = int)


    for alignment in blast_alignments:
        query_from = alignment["hsps"][0]["query_from"]-1 # Convert to 0-based index
        query_to = alignment["hsps"][0]["query_to"] # Keep as 1-based index
        qseq = alignment["hsps"][0]["qseq"].upper()
        hseq = alignment["hsps"][0]["hseq"].upper()
        seq_id = alignment["description"][0]["title"]
        midline = alignment["hsps"][0]["midline"]
        
        if read_dir == "R1": 
            if query_from != 0 :
                exlude_seqs += 1
                continue
        else: 
            if query_to != len(query_seq):
                exlude_seqs += 1
                continue
        
        if read_dir == "R1":
            read_out_of_frame = len(hseq) % 3  ## Correct for out of frame reads, due to different lengths of the reads, so that they end on a codon boundary
            alignments[seq_id] = {"qseq": qseq[:-read_out_of_frame], 
                                  "hseq": hseq[:-read_out_of_frame], 
                                  "midline": midline[-read_out_of_frame]}
            
        if read_dir == "R2":
            read_out_of_frame = len(hseq) % 3  # Correct for out of frame reads, due to different lengths of the reads
            alignments[seq_id] = {"qseq": qseq[read_out_of_frame:], 
                                  "hseq": hseq[read_out_of_frame:], 
                                  "midline": midline[read_out_of_frame:]}

        coverages[query_from:query_to] += 1

    print(f"{exlude_seqs} sequences are excluded, since they do not cover the start (R1) or end (R2) of the amplicon sequence.")

    return alignments, coverages


def characterize_DMS_blast_alignment(DMS_alignments, ref, data_type = "AA", read_dir = "R1", exclude_not_covered_regions = True):
    """
    Counts mutations (AA/DNA/Codon), insertions, and deletions per position from DMS BLAST alignment.

    Parameters:
    - DMS_alignments: dictionary with the sequences of insert that is mutated 
    - ref: string, reference DNA sequence 
    - data_type: string, "AA", "DNA" or "Codons
    - read_dir: string, "R1" or "R2"

    Returns:
    - all_variants: dictionary, with the counts of the variants per position
    - indels: pd.DataFrame, with the counts of insertions and deletions per position
    - enrichment_counts: pd.DataFrame, with the counts of the variants per position, with the reference sequence masked
    - enrichment_relative: pd.DataFrame, with the relative counts of the variants per position, with the reference sequence masked
    """

    all_variants = {}
    seq_with_off_target_indels = 0
    included_seq = 0
    indels = pd.DataFrame(columns = range(len(ref)), index = ["insertion", "deletion"], data = 0)


    ref_prot = translate_dna2aa(ref)

    if data_type == "AA":
        for idx in range(len(ref_prot)):
            all_variants[idx] = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 
                                    'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 
                                    'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 
                                    'W':0, 'Y':0, '*':0} ## X: are triplets with missing nucleotides (i.e. in the aligned seq, "-" is present), that would lead to a frameshift
    elif data_type == "DNA": 
        for idx in range(len(ref)):
            all_variants[idx] = {'A':0, 'C':0, 'G':0, 'T':0}

    else:
        codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
        for idx in range(len(ref_prot)): 
            all_variants[idx] = {codon: 0 for codon in codons}

    for alignment in DMS_alignments.values():
        qseq = alignment["qseq"]
        hseq = alignment["hseq"]

        if "-" in hseq or "-" in qseq:
            seq_with_off_target_indels += 1
            shift = 0 # Count the shift of the position compared to the reference that occurs if there is an insertion in the qseq 

            for idx,nt in enumerate(qseq): 
                pos = idx - shift # Adjust for the shift in the index due to prior insertions

                if read_dir == "R2": 
                    pos = len(ref) - (len(hseq)-qseq.count("-")) + pos

                if hseq[idx] == "-":
                    indels.loc["deletion", pos] += 1

                if nt == "-":
                    indels.loc["insertion", pos] += 1
                    shift += 1 # Correct for the shift in the index due to the insertion
                
            continue

        if data_type == "AA":
            hseq = translate_dna2aa(hseq)
            
        elif data_type == "Codons":
            hseq = [hseq[i:i+3] for i in range(0, len(hseq), 3)]

        if read_dir == "R2": 
            hseq = hseq[::-1]
            for idx, variant in enumerate(hseq):
                all_variants[len(all_variants)-idx-1][variant] += 1
            included_seq +=1
        else: 
            for idx, variant in enumerate(hseq): 
                all_variants[idx][variant] += 1
            included_seq +=1
       
    print(seq_with_off_target_indels, "sequences with off target indels are excluded")
    print(included_seq, "sequences are included in the enrichment analysis")

    enrichment_df = pd.DataFrame.from_dict(all_variants)

    if exclude_not_covered_regions: 
        enrichment_df = enrichment_df.loc[:,enrichment_df.sum() > 0]

    enrichment_counts, enrichment_relative = mask_ref_in_variants_df(variant_df=enrichment_df, ref_seq=ref_prot if data_type=="AA" else ref, data_type=data_type, reverse = True if read_dir == "R2" else False)
    

    return all_variants, indels, enrichment_counts, enrichment_relative


def calc_mut_spectrum_from_enrichment(enrichment_df, ref_seq, data_type = "DNA", set_diag_to_NA = True):
    """
    Constructs a mutation matrix: ref (row) -> mutated (column).

    Parameters:
    - enrichment_df: Dataframe with the counts of each AA/Codon/Nt at each position
    - data_type: string, "AA", "DNA" or "Codons"
    - ref_seq: string, reference DNA (if data_type = "DNA" or "Codon") or AA (if data_type = "AA") sequence

    Returns: 
    - mut_spectrum: count matrix, DataFrame
    - mut_spectrum_perc: normalized matrix in %, DataFrame
    """

    if data_type == "DNA":
        variants = ["A", "C", "G", "T"]

    elif data_type == "Codons":
        variants = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
        ref_seq = [ref_seq[i:i+3] for i in range(0,len(ref_seq)//3*3,3)]

    elif data_type == "AA":
        variants = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*']

    mut_spectrum = pd.DataFrame(index = variants, columns = variants, data = 0, dtype = np.float64) # rows = reference, cols = mutated

    for idx, ref_var in enumerate(ref_seq): 
        for mut_nt in enrichment_df.index:
            mut_pos = enrichment_df.iloc[:,idx]
            mut_count = mut_pos[mut_nt]
            mut_spectrum.loc[ref_var, mut_nt] += mut_count if mut_count > 0 else 0
        
    # Optional: mask identity substitutions
    if set_diag_to_NA:
        np.fill_diagonal(mut_spectrum.values, np.nan)
    
    mut_spectrum_perc = mut_spectrum/mut_spectrum.sum().sum()*100
            
    return mut_spectrum, mut_spectrum_perc


