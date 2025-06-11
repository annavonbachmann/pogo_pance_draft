# --- Import necessary libraries ---
import os
from Bio.SeqIO import QualityIO
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import glob
from scripts.utils import dna_rev_comp, translate_dna2aa
import pandas as pd
import seaborn as sns
from scripts.preprocessing_functions import *
import pickle as pkl
import matplotlib.colors as mcolors
import os.path
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

# --- Filters sequencing reads based on mutation count threshold and optionally read length ---

def read_filtering(a_seqs, 
                   b_seqs,
                   ref, 
                   catch_left = "", 
                   catch_right = "", 
                   n_mut_treshold = 10, 
                   filter_for_read_len = None): 
    """
    Filters sequencing reads based on number of mutations and optional read length.

    Parameters:
    - a_seqs (list): Forward reads (R1).
    - b_seqs (list): Reverse reads (R2).
    - ref (str): Reference gene sequence.
    - catch_left (str): Start marker in R1 read.
    - catch_right (str): Start marker in R2 read (reverse complemented).
    - n_mut_treshold (int): Max allowed number of mutations.
    - filter_for_read_len (tuple): Minimum read lengths (R1, R2).

    Returns:
    - tuple: Filtered lists of forward and reverse reads.
    """
    # --- Print number of reads before filtering ---

    print("total forward reads before filtering", sum([a_Seq != "" for a_Seq in a_seqs]))
    print("total reverse reads before filtering", sum([b_Seq != "" for b_Seq in b_seqs]))

    # --- Optional: define read length thresholds ---

    if filter_for_read_len:
        if type(filter_for_read_len) != tuple:
            print("filter_for_read_len should be a tuple (threshold a_reads, threshold b_reads)")
            raise ValueError
        print("filtering for read length", filter_for_read_len)
        a_len_T = filter_for_read_len[0]
        b_len_T = filter_for_read_len[1]

    # --- Initialize output lists ---
    a_sequences = []
    b_sequences = []

    # --- Iterate through paired reads ---

    for a_seq, b_seq in zip(a_seqs, b_seqs):
            # --- Process forward read (R1) --- 
            if catch_left in a_seq:

                if filter_for_read_len and len(a_seq) < a_len_T:
                    a_sequences.append("")
                
                else:
                    index = a_seq.index(catch_left) + len(catch_left)
                    gene_a = a_seq[index:]
                    total_muts_a = sum([ref[idx] != gene_a[idx] for idx in range(len(gene_a))])

                    # --- Filter based on mutation threshold ---
                    if total_muts_a <= n_mut_treshold: 
                        a_sequences.append(a_seq)
                    else: 
                        a_sequences.append("") #append by empty string to keep the index of the reads in the list, i.e. R1 and R2 can be matched later
            else: 
                a_sequences.append("")
            # --- Process reverse read (R2) ---
            if dna_rev_comp(catch_right) in b_seq:
                
                if filter_for_read_len and len(b_seq) < b_len_T:
                    b_sequences.append("")

                else:
                    index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                    gene_b = dna_rev_comp(b_seq[index:])
                    total_muts_b = sum([ref[::-1][idx] != gene_b[::-1][idx] for idx in range(len(gene_b))])
                    # --- Filter based on mutation threshold ---

                    if total_muts_b <= n_mut_treshold:
                        b_sequences.append(b_seq)
                    else:
                        b_sequences.append("")
            else: 
                b_sequences.append("")
    # --- Print number of reads after filtering ---

    print("total forward reads after filtering", sum([a_seq != "" for a_seq in a_sequences]))
    print("total reverse reads after filtering", sum([b_seq != "" for b_seq in b_sequences]))

    return a_sequences, b_sequences

# --- Collects amino acid variants from sequences ---

def gather_AA_variants(a_seq, 
                       b_seq, 
                       ref,
                       catch_left = "", 
                       catch_right = "",  
                       use_rev_read=True, 
                       use_forward_read = True, 
                       mask_ref = False):
    """
    Collects amino acid counts at each position from sequencing reads.

    Parameters:
    - a_seq (list): Forward reads (R1).
    - b_seq (list): Reverse reads (R2).
    - ref (str): Reference protein sequence.
    - catch_left (str): Start sequence in R1 reads.
    - catch_right (str): Start sequence in R2 reads (reverse complemented).
    - use_rev_read (bool): Whether to process R2.
    - use_forward_read (bool): Whether to process R1.
    - mask_ref (bool): Whether to mask wild-type amino acids (not used here).

    Returns:
    - dict: Dictionary of position-wise amino acid counts.
    """
    # --- Initialize mutation dictionary with all AA positions and possible AAs ---
    mutation_dict = {}
    prot_len = len(ref)

    for idx in range(prot_len):
        mutation_dict[idx] = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 
                            'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 
                            'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 
                            'W':0, 'Y':0, '*':0}
    # --- Iterate through each pair of reads ---
    for a_seq, b_seq in zip(a_seq, b_seq):
        # --- Process forward read ---

        if use_forward_read:

            if catch_left in a_seq:
                
                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]
                tr_a = translate_dna2aa(gene_a)

                for idx, pos in enumerate(tr_a):
                    mutation_dict[idx][pos] += 1
        # --- Process reverse read ---

        if use_rev_read: 

            if dna_rev_comp(catch_right) in b_seq:

                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                tr_b = translate_dna2aa(gene_b)
                tr_b = tr_b[::-1]
            
                for idx, pos in enumerate(tr_b):
                    mutation_dict[prot_len-idx-1][pos] += 1

    return mutation_dict

# --- Collects codon variants from sequencing reads ---

def gather_codon_variants(a_seq, 
                          b_seq, 
                          ref,
                          catch_left = "", 
                          catch_right = "",
                          use_rev_read= True,
                          use_forward_read = True):
    """
    Collects codon counts at each position from sequencing reads.

    Parameters:
    - a_seq (list): Forward reads (R1).
    - b_seq (list): Reverse reads (R2).
    - ref (str): Reference gene sequence.
    - catch_left (str): Start sequence in R1 reads.
    - catch_right (str): Start sequence in R2 reads (reverse complemented).
    - use_forward_read (bool): Whether to process R1.
    - use_rev_read (bool): Whether to process R2.

    Returns:
    - dict: Dictionary of position-wise codon counts.
    """

    mutation_dict = {}
    codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
    gene_len = len(ref)

    for idx in range(0, gene_len//3):

        mutation_dict[idx] = {codon: 0 for codon in codons}

    for a_seq, b_seq in zip(a_seq, b_seq):

        if use_forward_read: 

            if catch_left in a_seq:

                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]

                for i in range(0,len(gene_a)//3*3,3): # triplets, exclude last codon if not complete
                    mutation_dict[i//3][gene_a[i:i+3]] += 1

        if use_rev_read:

            if dna_rev_comp(catch_right) in b_seq:

                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index]) # exclude last codon if not complete
                bSeq_codons = [gene_b[i:i+3] for i in range(0,len(gene_b),3)][::-1] # list of codons, start from the end of ref_gene

                if len(bSeq_codons) > 0:
                    for idx, bSeq_codon in enumerate(bSeq_codons):  
                        mutation_dict[len(mutation_dict)-idx-1][bSeq_codon] += 1

    return mutation_dict

# --- Collects nucleotide variants from sequencing reads ---

def gather_nt_variants(a_seq, 
                       b_seq,
                       ref, 
                       catch_left = "", 
                       catch_right = "", 
                       use_rev_read= True,
                       use_forward_read = True):
    """
    Collects nucleotide counts at each position from sequencing reads.

    Parameters:
    - a_seq (list): Forward reads (R1).
    - b_seq (list): Reverse reads (R2).
    - ref (str): Reference gene sequence.
    - catch_left (str): Start sequence in R1 reads.
    - catch_right (str): Start sequence in R2 reads (reverse complemented).
    - use_forward_read (bool): Whether to process R1.
    - use_rev_read (bool): Whether to process R2.

    Returns:
    - dict: Dictionary of position-wise nucleotide counts.
    """

    mutation_dict = {}
    gene_len = len(ref)
    
    for idx in range(gene_len):
        mutation_dict[idx] = {'A':0, 'T':0, 'G':0, 'C':0}

    for a_seq, b_seq in zip(a_seq, b_seq):

        if use_forward_read: 

            if catch_left in a_seq:
                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]

                for idx, pos in enumerate(gene_a):
                    mutation_dict[idx][pos] += 1

        if use_rev_read:

            if dna_rev_comp(catch_right) in b_seq:

                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                gene_b = gene_b[::-1]

                for idx, pos in enumerate(gene_b):
                    mutation_dict[gene_len-idx-1][pos] += 1

    return mutation_dict

# ---  Processes reads for given variants ---
def process_reads(ref_prot, 
                  ref_gene,
                  variants = None, 
                  base_dir = os.getcwd(),
                  catch_left = "", 
                  catch_right = "", 
                  use_rev_read = True, 
                  use_forward_read = True, 
                  arbitrary_cutoff_a = False, 
                  arbitrary_cutoff_b = False, 
                  quality_score = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*','+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5'], 
                  n_mut_treshold = 10,
                  filter_for_read_len = None):
    """
    Process reads for given variants, i.e. generate dictionaries with the counts of each amino acid, codon and nucleotide at each position from fastq-files.

    Paramters:
    - variants: list of variants to process (variant names of the fastq files, which follow this structure {variant}_R1_001.fastq and {variant}_R2_001.fastq), if None, all variants in the fastq folder are processed
    - catch_left, catch_right: start (end) of the sequence in the forward read (R1) (reverse read (R2)), e.g. Barcodes (will not be included in the analysis)
    - use_forward_read, use_rev_read: whether or not to include the forward read (R1) and/or reverse read (R2) in the analysis (default: True)
    - arbitrary_cutoff_a, arbitrary_cutoff_b: where to cut off the forward (rev) sequence ( = maximum length of the reads, otherwise the cutoff is determined by the quality score = 1% (default) error rate)
    - quality_score: list of quality scores, at which the reads should be aborted (default: 1% error rate)
    - n_mut_treshold: number of mutations at which a read is excluded from the analysis (default: 10), if None, no filtering for n_mut
    - filter_for_read_len: tuple (threshold a_reads, treshold b_reads) to filter out reads with a len(read)<treshold (default: None, i.e. no filtering for read length)

    Returns:
    - returns: dict with the counts of each amino acid, codon and nucleotide at each position for each variant
    """

    variants_dict = {}
    path = f'{base_dir}/data/fastq'
    filenames = glob.glob(f'{path}/*')

    if variants is not None: # filter filenames for given variants
        filenames = [path for path in filenames if any(variant in path for variant in variants)]        

    for name in filenames: 

        if '_R1' in name:

            name = name.split('/')[-1].split('_R')[0]
            f1 = name
            a_seq, b_seq = read_sequences(f1, arbitrary_cutoff_a=arbitrary_cutoff_a, arbitrary_cutoff_b=arbitrary_cutoff_b, catch_left=catch_left, catch_right=catch_right, quality_score=quality_score)

            if n_mut_treshold:
                a_seq, b_seq = read_filtering(a_seq, b_seq, ref=ref_gene, catch_left=catch_left, catch_right=catch_right, n_mut_treshold=n_mut_treshold, filter_for_read_len=filter_for_read_len)

            variants_dict[name] = {}
            variants_dict[name] = get_variants(a_seq, b_seq, ref_prot = ref_prot, ref_gene = ref_gene, use_rev_read=use_rev_read, use_forward_read=use_forward_read, catch_left=catch_left, catch_right=catch_right)

            print(f'Done: {name}')

    return variants_dict

# --- Extracts AA, Codon or Nucleotide variants for a given set of sequences ---
def get_variants(a_seq,b_seq, ref_prot, ref_gene ,catch_right , catch_left , use_rev_read=True,use_forward_read=True):
    """
    Get the amino acid, codon and nucleotide variants for a given set of sequences

    Parameters:
    - a_seq: list of forward reads (R1)
    - b_seq: list of reverse reads (R2)
    - ref_prot: reference AA sequence
    - ref_gene: reference DNA sequence
    - catch_left, catch_right: start (end) of the sequence in the forward read (R1) (reverse read (R2)), e.g. Barcodes (will not be included in the analysis)
    - use_forward_read, use_rev_read: whether or not to include the foward read (R1) and/or reverse read (R2) in the analysis (default: True)

    Returns: 
    - dictionary with the counts of each amino acid, codon and nucleotide at each position
    """
    
    variants_dict = {}

    variants_dict["AA"] = gather_AA_variants(a_seq, b_seq, use_rev_read=use_rev_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref=ref_prot)
    variants_dict["DNA"] = gather_nt_variants(a_seq, b_seq, use_rev_read=use_rev_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref=ref_gene)
    variants_dict["Codons"] = gather_codon_variants(a_seq, b_seq, use_rev_read=use_rev_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref=ref_gene)

    return variants_dict

# --- Masks reference Nt/Codon/AA with counts of each Nt/Codon/AA at each position ---
def mask_ref_in_variants_df(variant_df:pd.DataFrame,
                            ref_seq:str, 
                            data_type:str,
                            reverse:bool = False):
    """
    Mask reference Nt/Codon/AA in dataframe with the counts of each Nt/Codon/AA at each position

    Parameters:
    - variants_df: dataframe with the counts of each AA/Codon/Nt at each position
    - ref_seq: reference DNA (if data_type = "DNA" or "Codon") or AA (if data_type = "AA") sequence
    - data_type: "AA", "DNA" or "Codon"
    - reverse: whether the reverse read only is used, i.e. the analysis focuses on the end of the reference sequence 

    Returns: 
    - pd dataframe with the counts, pd.dataframe with relative frequencies
    """
    variant_df = variant_df.copy()
    read_len = variant_df.shape[1]
    total_counts = variant_df.sum() # sum total counts before masking
    
    if data_type in ["DNA", "AA"]:

        if reverse: 
            for idx in range(read_len):
                variant_df.loc[ref_seq[::-1][idx], variant_df.columns[read_len-idx-1]] = np.nan #select column based on idx
            
        else: 
            for idx in range(read_len):
                variant_df.loc[ref_seq[idx], variant_df.columns[idx]] = np.nan

    elif data_type == "Codons":

        codons = [ref_seq[idx:idx+3] for idx in range(0,len(ref_seq),3)]

        if reverse: 
            for idx in range(read_len):
                variant_df.loc[codons[::-1][idx], variant_df.columns[read_len-idx-1]] = np.nan
        
        else:
            for idx in range(read_len):
                variant_df.loc[codons[idx], variant_df.columns[idx]] = np.nan
    
    variant_df_relative = variant_df/total_counts # calculate relative frequencies (columns with total counts = 0 --> NaN)

    return variant_df, variant_df_relative

# --- Computes mutation spectrum (nucleotide substitutions) ---

def mut_spectrum(a_seq, 
                 b_seq, 
                 reference_seq, 
                 use_rev_read = True, 
                 use_forward_read = True, 
                 catch_left = "", 
                 catch_right = "", 
                 set_diag_to_NA = True):
    """
    Computes the nucleotide substitution spectrum from sequencing reads.

    Parameters:
    - a_seq (list): Forward reads (R1).
    - b_seq (list): Reverse reads (R2).
    - reference_seq (str): Reference gene sequence.
    - use_forward_read (bool): Whether to process R1.
    - use_rev_read (bool): Whether to process R2.
    - catch_left (str): Start sequence in R1 reads.
    - catch_right (str): Start sequence in R2 reads (reverse complemented).
    - set_diag_to_NA (bool): Whether to mask self-substitutions.

    Returns:
    - tuple: (Raw substitution count dataframe, relative frequency dataframe).
    """

    ## reference nt : {sequenced (mutated) nt: count}
    Nts = ['A', 'C', 'G', 'T']
    mut_spec = {ref_nt: {nt:0 for nt in Nts} for ref_nt in Nts}
    
    for a_seq, b_seq in zip(a_seq, b_seq):
                
        if use_forward_read: 

            if catch_left in a_seq:
                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]
            
                for idx, nt in enumerate(gene_a): 
                    mut_spec[reference_seq[idx]][nt] += 1
                
        if use_rev_read:        

            if dna_rev_comp(catch_right) in b_seq:

                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                gene_b = gene_b[::-1]

                for idx, nt in enumerate(gene_b):
                    mut_spec[reference_seq[::-1][idx]][nt] += 1
    
    mut_spec_df = pd.DataFrame.from_dict(mut_spec, orient='index', dtype = float)


    if set_diag_to_NA:
        np.fill_diagonal(mut_spec_df.values, np.nan)

    # ## calculate mutagenic spectrum in percentage
    # total_n_muts = sum([sum(value.values()) for value in mut_spec.values()])
    # mut_spec_perc = {ref_base: {mut_base: round(val/total_n_muts*100, 3) for mut_base, val in value.items()} for ref_base, value in mut_spec.items()}
    total_n_muts = mut_spec_df.sum().sum()
    mut_spec_perc = mut_spec_df/total_n_muts*100

    return mut_spec_df, mut_spec_perc


# --- Calculate the mutagenic spectrum codon-wise ---
def mut_spectrum_codons(a_seq,
                        b_seq, 
                        reference_seq, 
                        use_rev_read = False, 
                        use_forward_read = True, 
                        catch_left = "", 
                        catch_right = "",
                        set_diag_to_NA = True):
    """
    Calculates mutagenic spectrum for a given set of sequences

    Parameters:
    - a_seq, b_seq: list of foward (reverse) sequences (R1, R2)
    - reference_seq: reference DNA sequence
    - use_forward_read, use_rev_read: whether or not to include the foward read (R1) and/or reverse read (R2) in the analysis (default: True)
    - catch_left, catch_right: start (end) of the sequence in the forward read (R1) (reverse read (R2)), e.g. Barcodes (will not be included in the analysis)
    - set_diag_to_NA: whether or not to set the diagonal (no change) of the mutagenic spectrum to np.nan (default: True)
    
    Returns:
    - DataFrame with the total counts of the mutagenic spectrum
    - DataFrame with percentages of the mutagenic spectrum (rows = ref codon, columns = read codon)
    """

    codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']

    # Reference codon : {mutated codon: count}
    mut_spec = {ref_codon: {codon:0 for codon in codons} for ref_codon in codons}

    ref_codons = [reference_seq[i:i+3] for i in range(0,len(reference_seq)//3*3,3)]

    for a_seq, b_seq in zip(a_seq, b_seq):
                
        if use_forward_read: 

            if catch_left in a_seq:

                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]
                gene_a_codons = [gene_a[i:i+3] for i in range(0,len(gene_a)//3*3,3)]

                if len(gene_a_codons) > 0:
                    for idx, a_codon in enumerate(gene_a_codons): 
                            mut_spec[ref_codons[idx]][a_codon] += 1
                
        if use_rev_read:   

            if dna_rev_comp(catch_right) in b_seq:
                
                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                gene_b_codons = [gene_b[i:i+3] for i in range(0,len(gene_b)//3*3,3)][::-1] # start from end of ref sequence
                
                if len(gene_b_codons) > 0:
                    for idx, b_codon in enumerate(gene_b_codons):
                        mut_spec[ref_codons[::-1][idx]][b_codon] += 1
    
    
    mut_spec_df = pd.DataFrame.from_dict(mut_spec, orient = "index", dtype = float)

    if set_diag_to_NA:
        np.fill_diagonal(mut_spec_df.values, np.nan)

    # Calculate mutagenic spectrum in percentage
    total_n_muts = mut_spec_df.sum().sum()
    mut_spec_perc = mut_spec_df/total_n_muts*100

    return mut_spec_df, mut_spec_perc

# --- Finds positions with a mutation rate above threshold and position with coverage above coverage threshold ---
def find_mutated_pos(read_dict,
                    Section, 
                    ref_gene, 
                    Primer_seq, 
                    Primer_out_of_triplets, 
                    Bc = None, 
                    seq_include_Primer_start = False,
                    Barcodes = None, 
                    data_type = "AA", 
                    cyclename = "Mutagenesis", 
                    mut_rate_filter_treshold = 0.05, 
                    cov_filter_treshold=50):
    """
    Finds the positions with a mutation rate above the mut_rate_filter_treshold and the positions with a coverage above the cov_filter_treshold

    Parameters:
    - read_dict = dictionary with the reads (following this naming convention: {cyclename}_{Barcode}_{Section}_R1:[read1_a, read2_a], {cyclename}_{Barcode}_{Section}_R2: [read1_b, read2_b],...})
    - Bc = name of the barcode, if None, it is expected that the read does not include the BC seq anymore
    - seq_include_Primer_start = whether or not the sequence includes the primer start, i.e. whether correction for triplet starts already happended during multiplexing (default: False), only used if Bc = None
    - Barcodes = dictionary with the barcode sequences, following the structure {BC1_fwd : seq, BC1_rev : seq, BC2_fwd : seq, ...}, can be None, if the reads do not include the BC seq anymore (e.g. by calling "cut_BC_seq" during demultiplexing)
    - Section = name of the section of interest
    - ref_gene = reference gene sequence
    - Primer_seq = dictionary with primer sequences, following the structure {S1_fwd : seq, S1_rev : seq, S2_fwd : seq, ...}
    - Primer_out_of_triplets = dictionary with the number of nucleotides at the beginning of the primer seq before a triplet starts, following the structure {S1_fwd : int, S1_rev : int, S2_fwd : int, ...}
    - data_type = "AA", "Codons" "DNA"
    - cyclename = name of the cycle
    - filter_treshold = treshold for the mutation rate
    - cov_filter_treshold = treshold for the coverage

    Returns: 
    - list of positions with mutation rate above mut_rate_filter_treshold
    - list of positions with coverage above cov_filter_treshold
    """

    dataType_handler = {"DNA": gather_nt_variants, "Codons": gather_codon_variants, "AA": gather_AA_variants}
    gather_variants = dataType_handler.get(data_type)
    if not gather_variants: 
        print("Data type not found!")
        exit()

    ref_seq_Section = find_reference_seq(ref_gene = ref_gene, Primer_seq = Primer_seq, Section = Section, Primer_out_of_triplets = Primer_out_of_triplets)
    ref = ref_seq_Section if data_type != "AA" else translate_dna2aa(ref_seq_Section)

    tripl_st = Primer_out_of_triplets[Section+"_fwd"]
    tripl_end = Primer_out_of_triplets[Section+"_rev"]

    a_seq = read_dict[f"{cyclename}_{Section}_R1"]
    b_seq = read_dict[f"{cyclename}_{Section}_R2"]
    
    if Bc: 
        catch_left = Barcodes[f"{Bc}_fwd"]+Primer_seq[Section + "_fwd"][:tripl_st]
        catch_right = dna_rev_comp(Barcodes[f"{Bc}_rev"]+Primer_seq[Section+"_rev"][:tripl_end])
    else: 
        catch_left = "" if not seq_include_Primer_start else ""+Primer_seq[Section + "_fwd"][:tripl_st] ## if the primer start was already cut from the sequence, the catch_left is empty, i.e. the whole sequence is used for further analysis (the sequence was already corrected for triplet start)
        catch_right = "" if not seq_include_Primer_start else ""+dna_rev_comp(Primer_seq[Section+"_rev"][:tripl_end])

    seq_variants = gather_variants(a_seq=a_seq, b_seq = b_seq, catch_left=catch_left,catch_right=catch_right, ref=ref)

    seq_variants = pd.DataFrame.from_dict(seq_variants)

    coverages = seq_variants.sum()

    _, seq_variants_freq = mask_ref_in_variants_df(ref_seq = ref, variant_df = seq_variants, data_type = data_type)

    # Combine mutation rates
    seq_variants_freq = seq_variants_freq.sum(axis = 0)
    low_cov_pos = coverages[coverages<cov_filter_treshold].index 
    high_mut_positions = seq_variants_freq[seq_variants_freq > mut_rate_filter_treshold].index

    return list(high_mut_positions), list(low_cov_pos)


# --- Calculate sum of single, double and triple mutants ---
def gather_n_mutations(a_seq, 
                       b_seq, 
                       reference_seq, 
                       catch_left = "", 
                       catch_right = "", 
                       use_rev_read = True, 
                       use_forward_read = True, 
                       use_triplets = False, 
                       return_seqs_pos = False):
    """
    Creates a dictionary with the number of single, double, triple (...) mutants, also a dict with the seqs and the positions of the mutations (if return_seqs_pos = True)
    !! positions are based on the reference sequence (i.e. does not match location in (reverse) b_reads directly)
    !! if use_triplets = True, the positions refer to the codons (AAs), otherwise to the nucleotides

    Parameters:
    - a_seq, b_seq: list of sequences
    - reference_seq: reference DNA sequence
    - catch_left, catch_right: start (end) of the sequence in the forward read (R1) (reverse read (R2)), e.g. Barcodes (will not be included in the analysis)
    - use_forward_read, use_rev_read: whether or not to include the foward read (R1) and/or reverse read (R2) in the analysis (default: True)
    - use_triplets: if True, the analysis is done on codons, otherwise (default) on nucleotides 
    - return_seq_pos: if True, also return the sequences and the positions of the mutations as as second dictionary

    Returns: 
    - dictionary with the number of single, double, triple (...) mutants {n_muts : count}, if return_seqs_pos = True, also a dictionary with the sequences and the positions of the mutations {n_muts : [(aSeq1, bSeq1, mut_pos_aSeq1, mut_pos_bSeq1), ...]}
    """
    mutation_dict = {}  
    if return_seqs_pos: 
        mutation_seq_dict = {}

    ref_codons = [reference_seq[i:i+3] for i in range(0,len(reference_seq)//3*3,3)]  
    ref_len = len(reference_seq) if not use_triplets else len(ref_codons)
    
    for a_seq, b_seq in zip(a_seq, b_seq):
        
        if use_forward_read and catch_left in a_seq:
            index = a_seq.index(catch_left) + len(catch_left)
            gene_a = a_seq[index:]

            if use_triplets: 
                muts_a_seq = [reference_seq[i:i+3] != gene_a[i:i+3] for i in range(0,len(gene_a)//3*3,3)]
            else: 
                muts_a_seq = [reference_seq[i] != gene_a[i] for i in range(len(gene_a)//3*3)]

            n_muts_a_seq = sum(muts_a_seq)  # count number of mutations
            muts_pos_a = [pos for pos, mut_exists in enumerate(muts_a_seq) if mut_exists] # get positions of mutations

        else: 
            n_muts_a_seq = 0
            gene_a = ""
            muts_pos_a = []
            
        if use_rev_read and dna_rev_comp(catch_right) in b_seq:
                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index]) # exclude last codon if not complete
                    
                if use_triplets:
                    gene_b_codons = [gene_b[i:i+3] for i in range(0,len(gene_b),3)]
                    muts_b_seq = [ref_codons[::-1][i] != gene_b_codons[::-1][i] for i in range(len(gene_b_codons))] # start from the end because different len of ref and b_seq
                else:
                    muts_b_seq = [reference_seq[::-1][i] != gene_b[::-1][i] for i in range(len(gene_b))]
                
                n_muts_b_seq = sum(muts_b_seq)  # count number of mutations
                muts_pos_b = [ref_len-pos-1 for pos, mut_exists in enumerate(muts_b_seq) if mut_exists] ## get positions of mutations, based on the reference (not the reverse reads!), -1 because of 0-based index

        else: 
            n_muts_b_seq = 0
            muts_pos_b = []
            gene_b = ""

        n_muts = n_muts_a_seq + n_muts_b_seq

        if (catch_left in a_seq) or (use_rev_read and (catch_right) in b_seq): ## add to dictionary if at least one read contains the catch sequence 
            if n_muts in mutation_dict:
                mutation_dict[n_muts] += 1
                if return_seqs_pos: 
                    mutation_seq_dict[n_muts].append((gene_a, b_seq, muts_pos_a, muts_pos_b))
                
            else:
                mutation_dict[n_muts] = 1
                if return_seqs_pos: 
                    mutation_seq_dict[n_muts] = [(gene_a, gene_b, muts_pos_a, muts_pos_b)]

    if return_seqs_pos: 
        return mutation_dict, mutation_seq_dict
    else: 
        return mutation_dict
        
# --- Gets linker variants ---
def get_linker_variants(reads, 
                        seq_before_linker, 
                        seq_after_linker, 
                        total_seq,  
                        wt_linker, 
                        intended_changes,
                        adaptor_left,
                        rev_reads,
                        include_changes_after_linker = False, 
                        include_deletions = True, 
                        combine_other = True, 
                        filter_treshold = 0.05
                        ):
    """
    Get linker variants

    Parameters:
    - reads: list of sequences (if R2 read, call dna_rev_comp on the sequences prior to calling this function)
    - seq_before_linker: short DNA sequence before the linker (8-10 bp)
    - seq_after_linker: short DNA sequence after the linker
    - wt_linker: linker DNA sequence
    - total_seq_before_linker: sequence from the beginning of the read until the linker sequence
    - intended_changes: indels and mutations (deletions are handeled separately) that are intended to be introduced by the retron library
    - rev_reads: whether or not the reads are from the R2 read (on which dna_rev_comp was called already) (default: False)
    - include_changes_after_linker: whether or not to include reads that do not contain the sequence after linker, but are in principle long enough, in the analysis (default: False) -> are counted within "other" --> probably due to unintended changes (off-target retron editing), mutations or sequencing errors (intended changes are all located prior to or in the linker sequence, i.e. the sequence after the linker should not be affected by the retron editing)
    - include_deletions: whether or not to include deletions in the analysis (default: True) 
    - combine_other: whether or not to combine all "other" (not intended) sequences into one category (default: True)
    - filter_treshold: variants with frequency below filter_treshold (given in %) are filtered out (default: 0.05)

    Returns: 
    - dictionaries with (1) counts, (2) percents of linker variants, (3) percents of linker variants for AAs
    """
    print("total reads",len(reads)) # Number of sequences
    print("reads with target sequence", sum([seq_after_linker in seq for seq in reads])) # Sum of seqences that include the sequence after the linker

    linker_variants = {}

    for seq in reads:   # Indels and mutations
        if seq_before_linker in seq and seq_after_linker in seq: # Only consider reads that contain the linker position
            start_idx = seq.index(seq_before_linker) + len(seq_before_linker)
            stop_idx = seq.index(seq_after_linker)
            linker = seq[start_idx:stop_idx]

            if linker in intended_changes: # Include intended changes
                if linker in linker_variants.keys():
                    linker_variants[linker] += 1
                else: 
                    linker_variants[linker] = 1

            elif linker == wt_linker:
                if "wt" in linker_variants.keys():
                    linker_variants["wt"] += 1  
                else: 
                    linker_variants["wt"] = 1
                
            else:  # Include changes not intended to "other" category
                if combine_other: 
                    if "other" in linker_variants.keys():
                        linker_variants["other"] += 1
                    else: 
                        linker_variants["other"] = 1
                elif "other_" + linker in linker_variants.keys():
                    linker_variants["other_" + linker] += 1
                else: 
                    linker_variants["other_" + linker] = 1

        elif include_deletions: # Deletions
            deleted_seq = seq_before_linker if not rev_reads else seq_after_linker # Sequence that is deleted (before or after depending on the read orientation)
            intact_seq = seq_after_linker if not rev_reads else seq_before_linker # Sequence that is not deleted (before or after depending on the read orientation)

            if intact_seq in seq and deleted_seq not in seq: # Because for deletions, we delete at least three bases from the seq before the linker, i.e. seq_before linker is not anymore completely in seq
                seq_until_deletion = seq[len(adaptor_left):seq.index(intact_seq)] if not rev_reads else seq[seq.index(intact_seq)+len(intact_seq):-len(adaptor_left)] # Sequence with deletion
                seq_before_deletion = seq_until_deletion[:10]  if rev_reads else seq_until_deletion[:-10]
                if seq_before_deletion in total_seq:
                    if not rev_reads:
                        del_len = total_seq.index(seq_before_deletion)+len(seq_before_deletion) - (len(total_seq)+len(wt_linker)) ## length of deletion
                    else: 
                        del_len = -(total_seq.index(seq_before_deletion) + len(wt_linker)) # Length of deletion
                    
                    delname = "del"+str(del_len) if del_len !=0 else "del-0"
                    if del_len <=0: # No insertion but mutations in seq_before_linker, thus seq_before_linker is not in seq but we do not want to consider these reads
                        
                        if delname in linker_variants.keys():
                            linker_variants[delname] += 1
                        else:
                            linker_variants[delname] = 1

                    else: 
                        if del_len == 0:
                            if "other" in linker_variants.keys():
                                linker_variants["other"] += 1
                            else:
                                linker_variants["other"] = 1
                                
                        elif "other_"+delname in linker_variants.keys():
                            linker_variants["other"+delname] += 1
                        else:
                            linker_variants["other"+delname] = 1
                        

        elif include_changes_after_linker: # Other (untargeted changes that effect the sequence after the linker)
            if len(seq) >= len(total_seq)+len(wt_linker)+len(seq_after_linker):
                if "after_linker" not in linker_variants.keys():
                    linker_variants["after_linker"] = 1
                else:
                    linker_variants["after_linker"] += 1

    total_vars = sum(linker_variants.values())
    linker_variants_perc = {seq: count/total_vars*100 for seq,count in linker_variants.items()}
    # Exclude everything below given filter_treshold
    linker_variants_perc = {seq: count for seq,count in linker_variants_perc.items() if count > filter_treshold}
    # Order after value size
    linker_variants_perc = dict(sorted(linker_variants_perc.items(), key = lambda x: x[1], reverse = True))
    # Convert dict keys to AAs
    linker_variants_perc_AA = {(translate_dna2aa(seq) if seq[0] in ["A","C","G","T"] else seq) : count for seq,count in linker_variants_perc.items()}
    linker_variants_perc_AA = dict(sorted(linker_variants_perc_AA.items(), key = lambda x: x[1], reverse = True))

    return linker_variants, linker_variants_perc, linker_variants_perc_AA

# --- Calculates mutagenic spectrum from enrichment dataframes ---

def calc_mut_spectrum_from_enrichment(enrichment_df, ref_seq, data_type = "DNA", set_diag_to_NA = True):
    """
    Calculate mutagenic spectrum from enrichment dataframes

    Parameters:
    - enrichment_df: dataframe with the counts of each AA/Codon/Nt at each position
    - data_type: "AA", "DNA" or "Codons"
    - ref_seq = reference DNA (if data_type = "DNA" or "Codon") or AA (if data_type = "AA") sequence

    Returns:
    - pd dataframe with the counts
    - pd.dataframe with relative frequencies
    """

    if data_type == "DNA":
        variants = ["A", "C", "G", "T"]

    elif data_type == "Codons":
        variants = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
        ref_seq = [ref_seq[i:i+3] for i in range(0,len(ref_seq)//3*3,3)]

    elif data_type == "AA":
        variants = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*']

    mut_spectrum = pd.DataFrame(index = variants, columns = variants, data = 0, dtype = np.float64) ## rows = reference, cols = mutated

    for idx, ref_var in enumerate(ref_seq): 
        for mut_nt in enrichment_df.index:
            mut_pos = enrichment_df.iloc[:,idx]
            mut_count = mut_pos[mut_nt]
            mut_spectrum.loc[ref_var, mut_nt] += mut_count
    
    if set_diag_to_NA:
        np.fill_diagonal(mut_spectrum.values, np.nan)
    
    mut_spectrum_perc = mut_spectrum/mut_spectrum.sum().sum()*100
            
    return mut_spectrum, mut_spectrum_perc
