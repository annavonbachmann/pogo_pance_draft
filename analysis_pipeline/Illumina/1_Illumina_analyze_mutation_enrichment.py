# ====================== Please run 1_ILLUMINA_PREPROCESS_AND_ALIGN_READS.PY before running this script ======================

# --- Import necessary modules ---
import os
import sys
import json
from Bio.SeqIO import QualityIO
import numpy as np
from scripts.utils import dna_rev_comp, translate_dna2aa
import pandas as pd
from scripts.plotting import *
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from scripts.Illumina_functions import *
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl


# ====================== INPUT PARAMETERS ======================


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

# --- Load experiment configuration from JSON file ---
with open(f"{data_dir}/config.json", "r") as file:
    config = json.load(file)

# --- Define read directions to be analyzed ---
read_directions = [ "R1", "R2"] # ["R1", "R2"] or ["R1"] or ["R2"]

# --- Define data types to analyze: DNA, protein (AA), or Codons ---
datatypes = [ "DNA", "AA", "Codons"]

# --- Define start and end sequences that span the region of interest (ROI), here, start (end) of LOV ---
roi_startseq = "ttagccacaa".upper() 
roi_endseq = "cggccaaa".upper()

# --- Set filtering behavior based on region of interest ---
filter_for_reads_with_roi = True # if True, only reads that include the roi are considered for the analysis
cut_to_roi = True # if True, the reads will be filtered for the region of interest (region before/after the roi will be filtered out), if False, the whole read will be considered for the analysis

# --- Extract additional parameters from config ---
variant = config["variant"] 
used_Barcodes = config["used_Barcodes"]
Sections = config["Sections"] 
full_amplicon = config["amplicon"]#[2:] # May be different to the specific reference sequence, if the amplicon was split into different parts for the analysis, IMPORTANT: if you cut the start of the amplicon to keep the reads in frame, please adjust the full_amplicon accordingly
full_amplicon_AA = translate_dna2aa(full_amplicon)

# --- Minimum read coverage required to include position in downstream analysis ---
min_coverage = 2000 

# --- Define visual appearance of coverage barplots ---
bar_color = "#22577A"  

# --- Define color gradient for mutation enrichment plots ---
cmap = LinearSegmentedColormap.from_list("custom_cmap", [
    "#22577A",  # Deep blue
    "#38A3A5",  # Teal
    "#57CC99",  # Medium green
    "#80ED99",  # Bright green
    "#C7F9CC"   # Light pastel green
] , N=256)

# --- Set matplotlib style parameters ---
fs = 8  # font size
plt.rcParams['svg.fonttype'] = 'none'
mpl.rcParams.update({
    'font.family': 'Avenir Next',
    'font.weight': 'demi', 
    'font.size': fs,
    'text.color': '#231F20',
    'axes.labelcolor': '#231F20',
    'xtick.color': '#231F20',
    'ytick.color': '#231F20',
    'axes.edgecolor': '#231F20',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'pdf.fonttype': 42,
    'text.usetex': False
})

# --- Define codon to amino acid translation table ---
genetic_code = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
codons = list(genetic_code.keys())


# --- Loop through all selected data types (DNA, AA, Codons) ---
for data_type in datatypes:

    print("############# Analysis for", data_type, "#############")
   
   # --- Choose correct reference depending on data type ---
    full_reference = full_amplicon if data_type != "AA" else full_amplicon_AA
    key_of_interest = "combined" if len(read_directions) > 1 else read_directions[0]

    # --- Define output folder for plots ---
    FigFolder = repo_root / 'final_output' / variant / key_of_interest / 'plots' / data_type
    if not os.path.exists(FigFolder):
        os.makedirs(FigFolder)
    
    # --- Define output folder for enrichment data ---
    OutputFolder = repo_root / 'final_output' / variant / key_of_interest / 'enrichments' / data_type
    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder) 

# ====================== Analysis of Mutation Enrichments ======================
    # --- Loop through all used barcodes ---
    for Bc in used_Barcodes:
        print("##############", Bc, "##############")

        # --- Loop through all defined sections ---
        for Section in Sections: 
            print("##############", Section, "##############")
            blast_stemFilename = f"{variant}_{Bc}_{Section}_Nt_filt_" ## update accordingly, total name should be e.g. f"RL8_BC1_S1_Nt_filt_R1.out", same stem should be used for the reference file and the blast output file

            if not os.path.exists(f"{data_dir}/references/{blast_stemFilename}ref.fasta"):
                print(f"Reference file {data_dir}/references/{blast_stemFilename}ref.fasta does not exist, please update the the reference file path or run the blast analysis first")
                exit()

            # --- Load reference sequence and translate if needed ---
            amplicon_seq = str(SeqIO.read(f"{data_dir}/references/{blast_stemFilename}ref.fasta", "fasta").seq)
            amplicon_AA = translate_dna2aa(amplicon_seq)
            reference = amplicon_seq if data_type !="AA" else amplicon_AA

            # --- Determine positions of ROI within the reference ---
            roi_startidx_DNA = amplicon_seq.index(roi_startseq)
            roi_endidx_DNA = amplicon_seq.index(roi_endseq) + len(roi_endseq)

            roi_startidx_AA = roi_startidx_DNA//3
            roi_endidx_AA = roi_endidx_DNA//3


            # --- Create dictionary to store enrichment results for each read direction ---                
            all_enrichments = {read_dir:{} for read_dir in read_directions}

            # --- Process each read direction (R1 or R2) ---
            for read_dir in read_directions:

                # --- Load corresponding BLAST output file ---
                print("################",  read_dir,   "################")

                if not os.path.exists(f"{data_dir}/blast/alignments/{blast_stemFilename}{read_dir}.out"):
                    print("Blast output file does not exist, please check the blast output file path or the input parameters, e.g. Barcodes, Sections, read directions")
                    exit

                with open(f"{data_dir}/blast/alignments/{blast_stemFilename}{read_dir}.out", "r") as file:
                    blast_output = json.load(file)

                blast_alignments = blast_output["BlastOutput2"][0]["report"]["results"]["search"]["hits"].copy()

                # --- Filter blast alignments for regions that include the region of interest (e.g. LOV2 insertion site) (span at least 10 nucleotides before and after the region) ---
                print(len(blast_alignments), "total alignments")
               
                if filter_for_reads_with_roi:
                    print("Filtering for reads with region of interest")
                    filter_for_region = roi_startidx_DNA if read_dir=="R1" else roi_endidx_DNA

                    blast_alignments = [alignment for alignment in blast_alignments if alignment["hsps"][0]["query_from"] <= filter_for_region-10 and alignment["hsps"][0]["query_to"] >= filter_for_region+10]

                    print(len(blast_alignments), "alignments after filtering filtering for reads with region of interest")
                
                
                # --- Convert BLAST hits to aligned sequence data ---
                alignments, all_coverages = restructure_alignments(blast_alignments, query_seq=amplicon_seq, read_dir=read_dir) ## get hseqs, qseqs, midline, all_coverages refers to the total coverage of the read, (**including** reads with indels (filtered in the next analysis), i.e. all reads returned in alignments)

                # --- Characterize mutations and indels from the alignments ---
                all_variants, indels,  enrichment_counts, enrichment_relative = characterize_DMS_blast_alignment(alignments, amplicon_seq, data_type=data_type,read_dir=read_dir, exclude_not_covered_regions=False if len(read_directions) > 1 else True)
                all_variants = pd.DataFrame.from_dict(all_variants)

                # --- Define slicing indices depending on whether to cut to ROI ---
                if cut_to_roi:
                    cut_start_idx = roi_startidx_DNA if data_type =="DNA" else roi_startidx_AA
                    cut_end_idx = roi_endidx_DNA if data_type =="DNA" else roi_endidx_AA
                else: 
                    cut_start_idx = 0
                    cut_end_idx = len(amplicon_seq) if data_type =="DNA" else len(amplicon_AA)
               
                # --- Adjust column indexing to match the full reference ---
                if key_of_interest == "combined": 
                    idxs_in_full_reference = [idx for idx in range(full_reference.index(reference), full_reference.index(reference)+len(reference))] if data_type != "Codons" else [idx for idx in range(full_reference.index(reference)//3, (full_reference.index(reference)+len(reference))//3)] 

                    all_variants.columns = idxs_in_full_reference

                    indels.columns = [idx for idx in range(full_amplicon.index(amplicon_seq), full_amplicon.index(amplicon_seq)+len(amplicon_seq))] ## always DNA level
                    indel_freq = indels/all_coverages
                    enrichment_counts.columns = idxs_in_full_reference
                    enrichment_relative.columns = idxs_in_full_reference

                # --- Store processed data ---
                all_enrichments[read_dir] = {
                    "all_variants": all_variants.iloc[:,cut_start_idx:cut_end_idx],
                    "indels": indels.iloc[:,cut_start_idx:cut_end_idx] if data_type == "DNA" else indels.iloc[:,cut_start_idx*3:cut_end_idx*3], # always DNA level
                    "indel_freqs": indel_freq.iloc[:,cut_start_idx:cut_end_idx] if data_type == "DNA" else indel_freq.iloc[:,cut_start_idx*3:cut_end_idx*3], # always DNA level
                    "all_coverages": all_coverages[cut_start_idx:cut_end_idx] if data_type == "DNA" else all_coverages[cut_start_idx*3:cut_end_idx*3], # always DNA level
                    "enrichment_counts": enrichment_counts.iloc[:,cut_start_idx:cut_end_idx],
                    "enrichment_relative": enrichment_relative.iloc[:,cut_start_idx:cut_end_idx],
                }
               
                # --- Filter out low-coverage positions ---

                read_depth = all_variants.sum() # `all_coverages` also includes reads that are excluded due to indels, here we want to have the coverage after filtering these reads

                all_enrichments[read_dir]["all_variants"].loc[:,read_depth < min_coverage] = 0
                all_enrichments[read_dir]["enrichment_counts"].loc[:,read_depth < min_coverage] = np.nan
                all_enrichments[read_dir]["enrichment_relative"].loc[:, read_depth < min_coverage] = np.nan
            
            # --- Set reference sequence to the region of interest ---
            if cut_to_roi:
                if key_of_interest == "combined":
                    roi_reference = reference[roi_startidx_AA:roi_endidx_AA] if data_type == "AA" else reference[roi_startidx_DNA:roi_endidx_DNA]
                else: 
                    roi_len = all_enrichments[key_of_interest]["enrichment_relative"].shape[1]
                    if data_type == "DNA":
                        roi_reference = reference[roi_startidx_DNA:roi_len+roi_startidx_DNA] if read_directions[0] == "R1" else reference[-roi_len+roi_endidx_DNA:roi_endidx_DNA]
                    elif data_type == "AA":
                        roi_reference = reference[roi_startidx_AA:roi_len+ roi_startidx_AA] if read_directions[0] == "R1" else reference[-roi_len+roi_endidx_AA:roi_endidx_AA]
                    elif data_type == "Codons": 
                        roi_reference = reference[roi_startidx_DNA:roi_len*3+roi_startidx_DNA] if read_directions[0] == "R1" else reference[-roi_len*3+roi_endidx_DNA:roi_endidx_DNA]
            else: 
                read_len = all_enrichments[read_directions[0]]["enrichment_relative"].shape[1] if data_type != "Codons" else all_enrichments[read_directions[0]]["enrichment_relative"].shape[1]*3
                roi_reference = reference[:read_len] if read_directions[0] != "R2" else reference[-read_len:]
           
            # ====================== COMBINE R1 AND R2 RESULTS ======================

            if key_of_interest == "combined": 

                all_enrichments["combined"] = {}

                # --- Combine variant counts from R1 and R2 ---
                all_enrichments["combined"]["all_variants"] =  all_enrichments["R1"]["all_variants"] + all_enrichments["R2"]["all_variants"]

                # --- Mask reference bases and calculate enrichments from combined variants ---
                all_enrichments["combined"]["enrichment_counts"], all_enrichments["combined"]["enrichment_relative"] = mask_ref_in_variants_df(all_enrichments["combined"]["all_variants"], ref_seq=roi_reference, data_type=data_type)
                
                # --- Combine indel counts and frequencies ---
                all_enrichments["combined"]["indels"] = all_enrichments["R1"]["indels"] + all_enrichments["R2"]["indels"]

                all_enrichments["combined"]["indel_freqs"] =  (all_enrichments["R1"]["indels"] + all_enrichments["R2"]["indels"])/all_enrichments["combined"]["all_variants"].sum()

                all_enrichments["combined"]["all_coverages"] = all_enrichments["R1"]["all_coverages"] + all_enrichments["R2"]["all_coverages"]

            # --- Define key to use for the analysis ---
            filename = f"{variant}_{Bc}_{Section}_{key_of_interest}_roi{cut_to_roi}_{data_type}"

            # --- Calculate mutation rates ---
            average_coverage = all_enrichments[key_of_interest]["all_variants"].sum().sum()/(all_enrichments[key_of_interest]["all_variants"].max().max()*all_enrichments[key_of_interest]["all_variants"].shape[1])*100
            print(f'The illumina paired reads cover on average {average_coverage.round(1)} % of the LOV sequence')
            
            # --- Compute total mutation rate per sequence (sum of all enrichment values) ---
            mut_rate = all_enrichments[key_of_interest]["enrichment_relative"].sum().sum()

            print(f'The mutation rate is estimated to be {mut_rate.round(3)} {data_type} mutations per sequence')
            
            # --- Save mutation rate and coverage stats to JSON ---
            mut_rates_dict = {"calculated_on": data_type, 
                            "coverage": average_coverage,
                            "mut_per_sequence": mut_rate,
                            }
            
            with open(f'{OutputFolder}/{filename}_mutation_rates.json', 'w') as file:
                file.write(json.dumps(mut_rates_dict, indent=4))
                
            # ====================== PLOT MUTATION ENRICHMENT ======================

            print("Plotting mutation enrichment...")
            plot_mutation_enrichment(all_enrichments[key_of_interest]["enrichment_relative"] , ref_seq=roi_reference, samplename=filename, data_type=data_type, FigFolder=FigFolder, cmap = cmap)

            # ====================== PLOT COVERAGE ======================
            print("Plotting coverage...")
            coverage_plot(all_enrichments[key_of_interest]["all_variants"].sum(), FigFolder=FigFolder, samplename = filename, color = bar_color)
        
            # ====================== PLOT INDEL FREQUENCIES (IF DNA) ======================
            print("Plotting indel frequencies...")
            if data_type == "DNA":

                plot_indel_freqs(all_enrichments[key_of_interest]["indel_freqs"], FigFolder=FigFolder, filename = filename, roi_start_idx=roi_startidx_DNA, roi_end_idx=roi_endidx_DNA)
                
                # --- Plot raw coverage before removing reads with indels ---
                plt.plot(all_enrichments[key_of_interest]["all_coverages"])
                plt.title("Coverage before filtering out reads with indels")
                plt.xlabel("Position")
                plt.ylabel("Coverage")
                plt.savefig(f"{FigFolder}/{filename}_coverage_before_filtering_indel_reads.pdf")
                plt.close()

            # ====================== MUTATIONAL SPECTRUM ANALYSIS ======================

            # --- Calculate and visualize mutation spectrum ---
            mut_spec, mut_spec_perc = calc_mut_spectrum_from_enrichment(all_enrichments[key_of_interest]["enrichment_relative"], ref_seq=roi_reference, data_type=data_type)
            
            print("Plotting mutagenic spectrum...")
            plot_mutation_spectrum(mut_spec_perc, FigFolder=FigFolder, samplename = filename, data_type=data_type, colormap = cmap)

            # ====================== SAVE ALL OUTPUTS ======================

            # --- Save mutational spectrum ---
            mut_spec_perc.to_csv(f"{OutputFolder}/{filename}_mut_spec.csv")
            # --- Save enrichment data ---
            all_enrichments[key_of_interest]["enrichment_relative"].to_csv(f"{OutputFolder}/{filename}_enrichment_relative.csv")
            all_enrichments[key_of_interest]["all_variants"].to_csv(f"{OutputFolder}/{filename}_all_variants.csv")
            all_enrichments[key_of_interest]["enrichment_counts"].to_csv(f"{OutputFolder}/{filename}_enrichment_counts.csv")
            # --- Save indel-related data (DNA only) ---
            if data_type == "DNA":
                all_enrichments[key_of_interest]["indel_freqs"].to_csv(f"{OutputFolder}/{filename}_indel_freq.csv")
                all_enrichments[key_of_interest]["indels"].to_csv(f"{OutputFolder}/{filename}_indel_counts.csv")
                pd.DataFrame(all_enrichments[key_of_interest]["all_coverages"]).to_csv(f"{OutputFolder}/{filename}_all_coverages.csv")
            