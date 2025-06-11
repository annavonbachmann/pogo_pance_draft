# --- Import necessary modules ---
import os
import sys
from Bio.SeqIO import QualityIO
import numpy as np
from matplotlib import pyplot as plt
from scripts.utils import dna_rev_comp, translate_dna2aa
import pandas as pd
import os.path
import json
from pathlib import Path
from scripts.plotting import *
from Bio import SeqIO
from scripts.Illumina_functions import *
from matplotlib.colors import LinearSegmentedColormap
from scripts.linker_analysis_functions import *

# ====================== CONFIG & PARAMETERS ======================


# --- Define folder containing `config.json`, `/references`, and `/blast/alignments` ---
data_dir = ""
current = Path(__file__).resolve()

for parent in current.parents:
    if (parent / 'final_output').exists():
        repo_root = parent
        break
else:
    raise FileNotFoundError("Could not find 'final_output' directory in any parent folder.")

# --- Check if directory exists ---
if not data_dir.exists() or not data_dir.is_dir():
    raise FileNotFoundError(f"Provided data directory does not exist: {data_dir}")

with open(f"{data_dir}/config.json", "r") as file:
    config = json.load(file)

read_directions = [ "R1", "R2"] # Whether to analyze forward reads (R1), reverse reads (R2), or both

# Define wild-type linker sequences for reference
wt_left_linker = "INESSGL"
wt_right_linker = "IDEAAKGSLHPP"

roi_R1 = "gccacaa".upper() 
roi_R2 = "tgctgaaaac".upper()
filter_for_reads_with_roi = True  # Only include reads that span the full ROI (region of interest): have to include one of the ROIs, but not necessarily both, depending on the read direction

variant = config["variant"] 
used_Barcodes = config["used_Barcodes"]
Sections = config["Sections"] 

min_coverage = 2000 # Minimum read coverage required to include a position in the downstream analysis
data_type = "AA" # Data type to analyze: "DNA" or "AA"

colors = ["#22577A", "#C7F9CC"] # Light green to light blue
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors , N=256)


print("############# Calculation for", data_type, "#############")

# ====================== MAIN LOOP OVER BARCODES ======================

for Bc in used_Barcodes:
    print("##############", Bc, "##############")
    for Section in Sections: 
        print("##############", Section, "##############")
        blast_stemFilename = f"{variant}_{Bc}_{Section}_Nt_filt_" # Update accordingly, total name should be e.g. f"RL8_BC1_S1_Nt_filt_R1.out", same stem should be used for the reference file and the blast output file

        if not os.path.exists(f"{data_dir}/references/{blast_stemFilename}ref.fasta"):
            print(f"Reference file {data_dir}/references/{blast_stemFilename}ref.fasta does not exist, please check the reference file path")
            exit()
        

        # Read and translate the full amplicon sequence
        amplicon_seq = str(SeqIO.read(f"{data_dir}/references/{blast_stemFilename}ref.fasta", "fasta").seq)
        amplicon_AA = translate_dna2aa(amplicon_seq)

        # Extract the LOV2 insert region based on start and end motifs
        roi_startidx = amplicon_seq.find(roi_R1)
        roi_endidx = amplicon_seq.find(roi_R2)+len(roi_R2)
        insert_DNA = amplicon_seq[roi_startidx:roi_endidx]
        insert_AA = translate_dna2aa(insert_DNA)
        
        # Prepare enrichment containers for each read direction (R1, R2)
        all_enrichments = {read_dir:{} for read_dir in read_directions}

        for read_dir in read_directions:

            FigFolder = repo_root / 'final_output' / variant / read_dir / 'plots' / data_type
            if not os.path.exists(FigFolder):
                os.makedirs(FigFolder)

            OutputFolder = repo_root / 'final_output' / variant / read_dir / 'enrichments' / data_type
            if not os.path.exists(OutputFolder):
                os.makedirs(OutputFolder)


            if not os.path.exists(f"{data_dir}/blast/alignments/{blast_stemFilename}{read_dir}.out"):
                print("Blast output file does not exist, please check the blast output file path")
                exit
        
            print("################",  read_dir,   "################")

            # Load the blast alignment file
            with open(f"{data_dir}/blast/alignments/{blast_stemFilename}{read_dir}.out", "r") as file:
                blast_output = json.load(file)

            all_blast_alignments = blast_output["BlastOutput2"][0]["report"]["results"]["search"]["hits"].copy()

            # Filter reads that fully span the LOV2 insert region
            print(len(all_blast_alignments), "total alignments")

            if filter_for_reads_with_roi:
                print("Filtering for reads with region of interest")
                filter_for_region = roi_startidx if read_dir=="R1" else roi_endidx

                blast_alignments = [alignment for alignment in all_blast_alignments if alignment["hsps"][0]["query_from"] <= filter_for_region-10 and alignment["hsps"][0]["query_to"] >= filter_for_region+10]

                print(len(blast_alignments), "alignments after filtering filtering for reads with region of interest")
            
            # Define sequence to cut at (start or end of LOV2 insert)
            cut_site_seq = roi_R1 if read_dir=="R1" else roi_R2

            print(read_dir, ": dividing reads at site", cut_site_seq)
            linker_alignments, LOV2_alignments, total_coverages = divide_alignments(blast_alignments, cut_site_seq, read_dir=read_dir, query_seq=amplicon_seq)
            total_coverages_LOV2 = total_coverages[roi_startidx:roi_endidx]

            # Characterize the aligned sequences in the LOV2 region
            all_variants, indels,  enrichment_counts, enrichment_relative = characterize_DMS_blast_alignment(LOV2_alignments, ref = insert_DNA , data_type=data_type,read_dir=read_dir, exclude_not_covered_regions=False)
            all_variants = pd.DataFrame.from_dict(all_variants)
            coverages = all_variants.sum() # Includes only the reads included in the enrichment analysis, i.e. reads with indels are **not** included here

            indels_freq = indels/total_coverages_LOV2 # coverages_LOV2 includes all reads, also reads with indels

            # Analyze the enrichment within the insert 
            print("analyze enrichment within the insert for", read_dir)
            if data_type == "AA":
                annot_ref = insert_AA[:enrichment_relative.shape[1]] if read_dir == "R1" else insert_AA[-enrichment_relative.shape[1]:]
            if data_type == "DNA":
                annot_ref = insert_DNA[:enrichment_relative.shape[1]] if read_dir == "R1" else insert_DNA[-enrichment_relative.shape[1]:]

            # Mask low-coverage columns
            all_variants.loc[:,coverages < min_coverage] = 0
            enrichment_counts.loc[:,coverages < min_coverage] = np.nan
            enrichment_relative.loc[:, coverages < min_coverage] = np.nan
            
            # Save results in dictionary
            all_enrichments[read_dir] = {"all_variants":all_variants, 
                                         "indels_freq":indels_freq, 
                                         "indels":indels,
                                         "enrichment_counts":enrichment_counts, "enrichment_relative":enrichment_relative, 
                                         "coverages_LOV2":total_coverages_LOV2
                                         }


            filename = f"{variant}_{Bc}_{Section}_{read_dir}_{data_type}"

            plot_mutation_enrichment(enrichment_relative, ref_seq=annot_ref, samplename=filename, data_type=data_type, FigFolder=FigFolder, vmax=None, cmap = cmap)

            plot_indel_freqs(indels_freq, filename=filename, FigFolder=FigFolder, color1 = colors[0], color2 = colors[1])

            # Save all enrichments
            all_enrichments[read_dir]["enrichment_relative"].to_csv(f"{OutputFolder}/{filename}_enrichment_relative.csv")
            all_enrichments[read_dir]["all_variants"].to_csv(f"{OutputFolder}/{filename}_all_variants.csv")
            all_enrichments[read_dir]["indels_freq"].to_csv(f"{OutputFolder}/{filename}_indels.csv")
            all_enrichments[read_dir]["enrichment_counts"].to_csv(f"{OutputFolder}/{filename}_enrichment_counts.csv")
            
            # Optionally analyze linker sequences
            if data_type == "AA":
                print("analyze linker variants for", read_dir)
                linkers , _ = get_linker_variants(linker_alignments,wt_linker = wt_left_linker if read_dir=="R1" else wt_right_linker,read_dir=read_dir) 

                # Sort linkers by frequency
                linkers_sorted = {k: v for k, v in sorted(linkers.items(), key=lambda item: item[1], reverse=True)}
                total_reads = sum(linkers_sorted.values())
                linkers_sorted_perc = {k: v/total_reads*100 for k, v in linkers_sorted.items()}

                # Exclude wt and linkers with less than 0.05% frequency
                linkers_sorted_perc.pop("wt")
                linkers_perc_filt = {k: v for k, v in linkers_sorted_perc.items() if v > 0.05}

                linkers_perc_filt, linker_renaming = rename_left_linkers(linkers_perc_filt.keys(), linkers_perc_filt) if read_dir=="R1" else rename_right_linkers(linkers_perc_filt.keys(), linkers_perc_filt)

                # Plot linker variants
                fig, ax = plt.subplots(1, 1, figsize=(15, 5))
                plt.bar(linkers_perc_filt.keys(), linkers_perc_filt.values(), color = colors[0])
                plt.xticks(rotation=90)
                plt.ylabel("Percentage of reads")
                plt.text(0.9, 0.93, f"Total mutation rate: {round(sum(linkers_perc_filt.values()),3)}", horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
                plt.title(f"Linker variants for {read_dir}")
                plt.savefig(f"{FigFolder}/{filename}_linker_distribution.pdf", bbox_inches="tight")
                plt.savefig(f"{FigFolder}/{filename}_linker_distribution.png", bbox_inches="tight")
                plt.close()
                
                linkers_perc_filt = pd.DataFrame.from_dict(linkers_perc_filt, orient="index")
                linkers_perc_filt.to_csv(f"{OutputFolder}/{filename}_linker_distribution.csv")

                with open(f"{OutputFolder}/{filename}_linker_renaming.json", "w") as file:
                    json.dump(linker_renaming, file)
        
    
        # Combine enrichment for R1 and R2
        if len(read_directions) == 1:
            continue

        FigFolder = f"{os.getcwd()}/final_output/{variant}/combined/plots/{data_type}"
        if not os.path.exists(FigFolder):
            os.makedirs(FigFolder)


        OutputFolder = f"{os.getcwd()}/final_output/{variant}/combined/enrichments/{data_type}"
        if not os.path.exists(OutputFolder):
            os.makedirs(OutputFolder)

        all_enrichments["combined"] = {}

        # Total variants of R1 and R2
        
        all_enrichments["combined"]["all_variants"] =  all_enrichments["R1"]["all_variants"] + all_enrichments["R2"]["all_variants"]

        # Total enrichments of R1 and R2

        all_enrichments["combined"]["enrichment_counts"], all_enrichments["combined"]["enrichment_relative"] = mask_ref_in_variants_df(all_enrichments["combined"]["all_variants"], ref_seq=insert_DNA if data_type == "DNA" else insert_AA, data_type=data_type, reverse = True if read_dir == "R2" else False)
        
        # Combine indels of R1 and R2
        all_enrichments["combined"]["indels_freq"] = ( all_enrichments["R1"]["indels"] + all_enrichments["R2"]["indels"])/(all_enrichments["R1"]["coverages_LOV2"] + all_enrichments["R2"]["coverages_LOV2"])

        filename = f"{variant}_{Bc}_{Section}_combined_{data_type}"

        print("analyze enrichment within the insert for combined reads")
        plot_mutation_enrichment(all_enrichments["combined"]["enrichment_relative"] , ref_seq=insert_DNA if data_type =="DNA" else insert_AA, samplename=filename, data_type=data_type, FigFolder=FigFolder, vmax=None, cmap = cmap)

        plot_indel_freqs(all_enrichments["combined"]["indels_freq"], filename=filename, FigFolder=FigFolder, color1 = colors[0], color2 = colors[1])

        # Save all enrichments
        all_enrichments["combined"]["enrichment_relative"].to_csv(f"{OutputFolder}/{filename}_enrichment_relative.csv")
        all_enrichments["combined"]["all_variants"].to_csv(f"{OutputFolder}/{filename}_all_variants.csv")
        all_enrichments["combined"]["indels_freq"].to_csv(f"{OutputFolder}/{filename}_indels.csv")
        all_enrichments["combined"]["enrichment_counts"].to_csv(f"{OutputFolder}/{filename}_enrichment_counts.csv")

        print( "DONE!")