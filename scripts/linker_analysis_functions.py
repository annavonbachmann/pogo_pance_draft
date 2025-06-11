# --- Import necessary libraries ---
import pandas as pd 
import pysam
import os
from Bio import SeqIO
from scripts.utils import translate_dna2aa
from scripts.functions_ import mask_ref_in_variants_df
import numpy as np
import re

# --- Extract linker variants from BLAST alignments ---
def get_linker_variants(linker_alignments, wt_linker = "SG", read_dir = "R1"):
    """
    Characterizes linker variants from BLAST alignments.

    Parameters:
    - linker_alignments (dict): Contains 'qseq' and 'hseq' alignments.
    - wt_linker (str): Wild-type linker amino acid sequence.
    - read_dir (str): Either 'R1' or 'R2' to distinguish read direction.

    Returns:
    - linker_counts (dict): Counts of each unique linker variant.
    - linker_list (list): All linker sequences parsed.
    """

    frameshifts = 0
    linker_counts = {}
    linker_list = []
    for x in linker_alignments.values():
        linker = ""
        qseq = x["qseq"]
        hseq = x["hseq"]
        
        # --- Adjust for R1 reads if not in frame ---
        if read_dir == "R1" and len(hseq)%3 != 0: 
            hseq = hseq[(len(hseq)%3):]
            qseq = qseq[(len(qseq)%3):]

        # --- Exclude reads that result in frameshifts ---
        is_frameshift_read = (qseq.count("-") - hseq.count("-")) %3 != 0
        if is_frameshift_read:
            frameshifts += 1
            continue

        # --- Wild-type sequence ---
        elif qseq == hseq:
            linker_counts["wt"] = linker_counts.get("wt", 0) + 1
            linker = "wt"

        # --- Deletions (multiples of 3 only) ---
        elif (qseq.count("-") - hseq.count("-")) < 0:  
                correct_by = 3 - (hseq.count("-") % 3) 
                del_count = (hseq.count("-") - qseq.count("-")) //3
                hseq_filt = re.sub("-", "", hseq)
                if read_dir == "R2": 
                    linker = translate_dna2aa(hseq_filt)[:len(wt_linker)-(del_count)] 
                else: 
                    hseq_filt = hseq_filt[correct_by:]
                    linker = translate_dna2aa(hseq_filt)[-len(wt_linker)+(del_count):]

                linker_counts[linker] = linker_counts.get(linker, 0) + 1

        # --- Substitutions only ---
        elif (qseq.count("-") - hseq.count("-")) == 0: 
            hseq_filt = re.sub("-", "", hseq)
            correct_by = 3 - (hseq.count("-") %3) 
            del_count = hseq.count("-")
            if read_dir == "R2":
                if hseq_filt[:len(wt_linker) * 3] == qseq[:len(wt_linker) * 3]:
                    linker_counts["wt"] = linker_counts.get("wt", 0) + 1
                else: 
                    linker = translate_dna2aa(hseq_filt)[:len(wt_linker)]
                    linker_counts[linker] = linker_counts.get(linker, 0) + 1
            else: 
                if hseq_filt[-len(wt_linker) * 3:] == qseq[-len(wt_linker) * 3:]:  
                    linker_counts["wt"] = linker_counts.get("wt", 0) + 1
                else: 
                    hseq_filt = hseq_filt[correct_by:]
                    linker = translate_dna2aa(hseq_filt)[-len(wt_linker):]
                    linker_counts[linker] = linker_counts.get(linker, 0) + 1
        # --- Insertions ---
        elif (qseq.count("-")- hseq.count("-")) > 0:
            insertion_len = (qseq.count("-") - hseq.count("-")) //3
            hseq_filt = re.sub("-", "", hseq)
            correct_by = 3 - (hseq.count("-") % 3)
            if read_dir == "R2": 
                linker = translate_dna2aa(hseq_filt)[:len(wt_linker) + insertion_len]
                linker_counts[linker] = linker_counts.get(linker, 0) + 1
            else: 
                hseq_filt = hseq_filt[correct_by:]
                linker = translate_dna2aa(hseq_filt)[-len(wt_linker) - insertion_len:]
                linker_counts[linker] = linker_counts.get(linker, 0) + 1
        else:
            print("sequence", hseq, "does not meet any criteria")
        
        linker_list.append(linker)
        
    print(frameshifts, "reads excluded due to frameshifts")

    return linker_counts, linker_list



# --- Rename and annotate right linker variants ---
def rename_right_linkers(linkernames, linkerperc_dict): 
    """
    Renames right linker sequences based on known patterns (e.g., terminal deletions).

    Parameters:
    - linkernames (list): Original linker variant names.
    - linkerperc_dict (dict): Original percentages associated with each linker.

    Returns:
    - linkers_perc_filt_corr_names (dict): Corrected names mapped to original values.
    - linker_renaming (dict): Original name → new name mapping.
    """
    linkers_perc_filt_corr_names = {}
    linker_renaming = {}

    for startname in linkernames: 
        name = startname
        end_n_del = 0
        start_n_del = 0
        ind = 0

        # --- Check and strip known suffixes ---
        for suffix in ["P", "P", "H", "L"]:
            if name.endswith(suffix) and end_n_del == 0:
                name = name[:-1]
            else:
                end_n_del += 1

        # --- Check and strip known prefixes ---
        for prefix in ["I", "D", "E", "A", "A", "K"]:
            if name.startswith(prefix):
                name = name[1:]
                ind += 1
            else:
                start_n_del += 1

        # --- Annotate deletions ---
        if end_n_del > 0:
            name = f"{name}(+{end_n_del}del)"
        if start_n_del > 0:
            name = f"(-{start_n_del}del){name}"

        print(startname , "->", name)
        linkers_perc_filt_corr_names[name] = linkerperc_dict[startname]
        linker_renaming[startname] = name

    return linkers_perc_filt_corr_names, linker_renaming
    

# --- Rename and annotate left linker variants ---
def rename_left_linkers(linkernames, linkerperc_dict): 
    """
    Renames left linker sequences based on deletion pattern.

    Parameters:
    - linkernames (list): Original linker variant names.
    - linkerperc_dict (dict): Original percentages associated with each linker.

    Returns:
    - linkers_perc_filt_corr_names (dict): Corrected names mapped to original values.
    - linker_renamings (dict): Original name → new name mapping.
    """
    linkers_perc_filt_corr_names = {} 
    linker_renamings = {}

    for startname in linkernames: 
        name = startname
        start_n_del = 0
        end_n_del = 0

        if name.endswith("L"): 
            name = name[:-1]
        else:
            end_n_del += 1

        for prefix in ["I", "N", "E", "S"]:
            if name.startswith(prefix):
                name = name[1:]
            else:
                start_n_del += 1

        if start_n_del > 0:
            name = f"(-{start_n_del}del){name}"
        if end_n_del > 0:
            name = f"{name}(+{end_n_del}del)"

        print(startname , "->", name)
        linkers_perc_filt_corr_names[name] = linkerperc_dict[startname]
        linker_renamings[startname] = name

    return linkers_perc_filt_corr_names, linker_renamings
    