# --- Import necessary modules ---
import os
from Bio.SeqIO import QualityIO
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scripts.utils import dna_rev_comp, translate_dna2aa
import pandas as pd
import seaborn as sns
import pickle as pkl
import matplotlib.colors as mcolors
import os.path
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from scripts.functions_ import *
from scripts.preprocessing_functions import *
import matplotlib

import matplotlib.pyplot as plt
import matplotlib as mpl

# ======================== PLOTTING SETTINGS ========================

# --- Seaborn theme configuration ---
custom_params = {
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 1
}
sns.set_theme(context="paper", style='ticks', palette="Greys_r", rc=custom_params)

# --- General matplotlib settings ---
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
sns.set_context("paper", rc={
    "font.size": fs,
    "axes.titlesize": fs + 1,
    "axes.labelsize": fs,
    "axes.linewidth": 1,
    "xtick.labelsize": fs,
    "ytick.labelsize": fs,
    "legend.fontsize": fs,
    "legend.title_fontsize": fs + 1
})


# ======================== PLOTTING FUNCTIONS FOR ANALYSIS ========================

def coverage_plot(coverage_df, 
                  samplename = "", 
                  FigFolder = None, color = "blue"): 
    """
    Plot the read coverage across positions.

    Parameters:
    - coverage_df (pd.Series): Coverage (read count) for each position.
    - samplename (str): Sample or variant name.
    - FigFolder (str): Optional path to save the plot.
    - color (str): Line color for coverage.
    """

    plt.plot(coverage_df, color = color)
    plt.xlabel("Position")
    plt.ylabel("Read counts")
    plt.title(f'{samplename} read depth')
    plt.xticks(list(range(0,len(coverage_df), 50)))
    # if FigFolder:
    #     if not os.path.exists(FigFolder):
    #         os.makedirs(FigFolder)
    #     plt.savefig(f'{FigFolder}/{samplename}_coverage.pdf')
    #     plt.savefig(f'{FigFolder}/{samplename}_coverage.png')
    plt.show()
    plt.close()



def plot_mutation_spectrum(data, 
                           samplename="" , 
                           FigFolder = None, 
                           colormap = "viridis",
                           data_type = "DNA",
                           figuresize = (4.76,6.69291/2),
                           ticks = [0.0, 0.5, 1.0]):
    """
    Plot mutation spectrum as a heatmap showing base substitutions.

    Parameters:
    - data (pd.DataFrame): Matrix with ref bases as rows and mutated bases as columns.
    - samplename (str): Sample name to include in plot title.
    - FigFolder (str): Optional path to save the figure.
    - colormap (str): Matplotlib colormap to use.
    - data_type (str): 'DNA' or 'AA'; determines if values are annotated.
    - figuresize (tuple): Figure dimensions in inches.
    - ticks (list): Custom ticks for colorbar.
    """
    f, ax = plt.subplots(figsize=figuresize)
    f.subplots_adjust(wspace=0.01)
    sns.heatmap(data, annot=True if data_type=='DNA' else False, linewidths=.5, ax=ax, cbar = False, square = True, linecolor = "black", cmap = colormap, vmin=0)

    im = ax.collections[0]
    cbar_ax = f.add_axes([0.15, -0.05, 0.25, 0.03])
    cbar = f.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_ticks(ticks)
    tick_labels = [f"{round(t) if i == len(ticks) - 1 else round(t, 1):.1f}" for i, t in enumerate(ticks)]
    cbar.ax.set_xticklabels(tick_labels)
    cbar.set_label("Proportion (%)", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    plt.xlabel('Mutated base (%)')
    ax.set_ylabel('Reference base (%)', labelpad=8)

    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize = 8)
    ax.tick_params(axis='y', pad=2) 
    ax.yaxis.set_label_coords(-0.2, 0.5) 
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(.5)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=1)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=1)
    plt.title(f"{samplename}")

    # if FigFolder:
        # if not os.path.exists(FigFolder):
            # os.makedirs(FigFolder)
        # plt.savefig(f"{FigFolder}/{samplename}_mutagenic_spectrum_perc.pdf")
        # plt.savefig(f"{FigFolder}/{samplename}_mutagenic_spectrum_perc.png")

    plt.show()
    plt.close()


def plot_mut_rate_per_pos(data, 
                          ref_seq,
                          data_type = "DNA",
                          samplename = "", 
                          FigFolder = None,
                          vmax = None,
                          figuresize = (20,2)):
    """
    Plot mutation rate per sequence position as a heatmap.

    Parameters:
    - data (pd.DataFrame): Relative frequency of each base/codon/AA at each position.
    - ref_seq (str): DNA or protein reference sequence.
    - data_type (str): 'DNA', 'AA', or 'Codons'.
    - samplename (str): Sample name for plot title.
    - FigFolder (str): Optional path to save the plot.
    - vmax (float): Maximum value for color scaling.
    - figuresize (tuple): Size of the figure in inches.
    """
    data = data.sum(axis = 0)

    ref_seq = ref_seq if data_type != "AA" else translate_dna2aa(ref_seq)

    if data_type == "Codons":
        x_ticklabels = [ref_seq[i:i+3] for i in range(0, len(ref_seq)//3*3, 3)]
    else: 
        x_ticklabels = [ref for ref in ref_seq]

    plt.figure(figsize=figuresize)
    sns.heatmap(pd.DataFrame(data).T, cmap = "viridis", cbar = True, cbar_kws = {"pad": 0.02, "label": "Mutation rate" },linecolor="black", xticklabels=x_ticklabels, yticklabels=False, vmax=vmax)
    plt.xlabel("Position")
    plt.xticks(rotation = 1 if data_type != "Codons" else 90,fontsize=6 if data_type != "DNA" else 3)
    plt.title(samplename)
    # if FigFolder:
        # if not os.path.exists(FigFolder):
            # os.makedirs(FigFolder)
        # plt.savefig(f"{FigFolder}/{samplename}_mutation_rate_per_{data_type}_position.pdf", bbox_inches="tight")
    plt.show()
    plt.close()



def plot_mutation_enrichment(variants_df, 
                             ref_seq, 
                             samplename = "",
                             data_type = "DNA",
                             FigFolder = None, 
                             cmap = "viridis", 
                             cbar_label = "mutation rate", 
                             vmax = None,
                             figuresize = (5.669288000000001,2.625)):
    """
    Plot mutation enrichment per position with heatmap.

    Parameters:
    - variants_df (pd.DataFrame): Relative counts matrix of each base/codon/AA.
    - ref_seq (str): Reference sequence.
    - samplename (str): Sample name for title.
    - data_type (str): 'DNA', 'AA', or 'Codons'.
    - FigFolder (str): Optional save path.
    - cmap (str): Color map.
    - cbar_label (str): Label for colorbar.
    - vmax (float): Max value for colorbar.
    - figuresize (tuple): Size of the figure.
    """

    if data_type in ["DNA", "AA"]:
        seq_pos = list(ref_seq)
    if data_type == "Codons":
        seq_pos = [ref_seq[i:i+3] for i in range(0, len(ref_seq)//3*3, 3)]
        print(seq_pos)

    plt.figure(figsize=figuresize)    
    ax = sns.heatmap(data=variants_df, cmap=cmap, cbar_kws={'label': cbar_label, "pad": 0.02}, yticklabels=True, xticklabels=False, center=0 if cmap == "coolwarm" else None, vmax=0.0024, vmin=-vmax if (cmap=="coolwarm" and vmax) else None)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=1)
    rotation = 90 if data_type == "Codons" else 1
    ax.yaxis.set_tick_params(width=2)
    ax.grid(False)
    plt.title(samplename)
    plt.ylim(top = .0024)
    plt.xlabel("Sequence")

    # if FigFolder:    
    #     if not os.path.exists(FigFolder):
    #         os.makedirs(FigFolder)
    #     plt.savefig(f"{FigFolder}/{samplename}_{data_type}_mutation_enrichment.pdf", bbox_inches="tight")
    #     plt.savefig(f"{FigFolder}/{samplename}_{data_type}_mutation_enrichment.png", bbox_inches="tight")
    
    plt.show()
    plt.close()


def plot_temporal_enrichment(enrichment_df_dict, 
                            ref,
                            combine_mut_rates = True,
                            color_above_vmax_orange = True,
                            vmax = 0.1,
                            show_cbar_for_each = False,
                            plt_titles = None,
                            show_plttitles = False,
                            plot_coverage = True,
                            bias_per_pos = None,
                            show_only_pos = None,
                            ref_annot = None,
                            return_df = True,
                            FigFolder = None,
                            figsize = (20,10),
                            data_type = "AA",
                            enrichment_cmap = "viridis", 
                            set_over_color = "orange", 
                            coverage_cmap = "magma"):
    
    """
    Plot mutation enrichment over time (or rounds), optionally with coverage and bias.

    Parameters:
    - enrichment_df_dict (dict): Dictionary with sample data, each containing 'mut_enrichment' and 'coverage'.
    - ref (str): Reference sequence.
    - combine_mut_rates (bool): Whether to sum over rows before plotting.
    - color_above_vmax_orange (bool): Highlight values above vmax in orange.
    - vmax (float): Max value for colorbar.
    - show_cbar_for_each (bool): Show individual colorbars per subplot.
    - plt_titles (list): Titles for each subplot.
    - show_plttitles (bool): Show titles on subplots.
    - plot_coverage (bool): Plot coverage as bottom panel.
    - bias_per_pos (list): Optional bias values to show.
    - show_only_pos (list): Index list of positions to subset.
    - ref_annot (list): Optional custom X-axis labels.
    - return_df (bool): Return combined dataframe.
    - FigFolder (str): Optional path to save.
    - figsize (tuple): Size of the overall figure.
    - data_type (str): 'DNA', 'AA', or 'Codons'.
    - enrichment_cmap (str): Colormap for enrichment heatmaps.
    - set_over_color (str): Color for out-of-range enrichment.
    - coverage_cmap (str): Colormap for coverage heatmap.
    """

    if return_df: 
        combined_df = pd.DataFrame(index = (enrichment_df_dict.keys()), columns = list(range(len(ref))), data = 0, dtype = np.float64)

    pltsize = len(enrichment_df_dict.keys()) if not plot_coverage else len(enrichment_df_dict.keys()) + 1
    pltsize = pltsize + 1 if bias_per_pos else pltsize

    fig, axes =  plt.subplots(pltsize,1, figsize=figsize)
    fig.subplots_adjust(hspace=0.5 if show_plttitles else 0.1)

    # --- Adjust subplot position if coverage and/or bias is shown ---
    if not show_plttitles and plot_coverage:
        last_ax = axes[-1] if not bias_per_pos else axes[-2] 
        pos = last_ax.get_position() 
        last_ax.set_position([pos.x0, pos.y0 - 0.015, pos.width, pos.height]) 
        if bias_per_pos:
            pos =  axes[-1].get_position() 
            axes[-1].set_position([pos.x0, pos.y0 - 0.015, pos.width, pos.height]) 
    idx = 0 
    # --- Plot mutation enrichment for each sample ---
    my_cmap = plt.get_cmap(enrichment_cmap).copy()
    if color_above_vmax_orange:
        my_cmap.set_over(set_over_color)

    for cycles, sample_data in enrichment_df_dict.items():
        if combine_mut_rates:
            variant_relative_freq = pd.DataFrame(sample_data["mut_enrichment"][1].sum(axis=0)).T
        else: 
            variant_relative_freq = sample_data["mut_enrichment"][1]

        if return_df:
            combined_df.loc[cycles] = variant_relative_freq.values[0]

        coverage = sample_data["coverage"]

        if show_only_pos:
            variant_relative_freq = variant_relative_freq.iloc[:,show_only_pos]
            coverage = coverage.iloc[:,show_only_pos]

        sns.heatmap(variant_relative_freq, annot=False, ax=axes[idx], linecolor="black", cmap=my_cmap,  cbar_kws={'label': "", "pad": 0.02}, vmin=0,vmax = vmax, xticklabels=False if idx != len(enrichment_df_dict.keys())-1 else True, yticklabels=True if combine_mut_rates == False else False, cbar=show_cbar_for_each )

        for _, spine in axes[idx].spines.items():
            spine.set_visible(True)
            spine.set_linewidth(2)
        axes[idx].set_yticklabels(axes[idx].get_yticklabels(), rotation=1, fontsize=5)
        
        if show_plttitles: 
            axes[idx].set_title(plt_titles[idx], fontsize=10)
        
        axes[idx].set_facecolor('gray')
        axes[idx].grid(False)

        labels = ref if not ref_annot else ref_annot
        ref_labels = [aa for aa in ref] if not show_only_pos else [labels[i] for i in show_only_pos]

        if idx == len(enrichment_df_dict.keys())-1:
            axes[idx].set_xticklabels(ref_labels, rotation=1, fontsize=4 if not show_only_pos else 10)
        
        idx += 1

    # --- Add coverage heatmap ---
    if plot_coverage:
        sns.heatmap(pd.DataFrame(coverage), ax=axes[pltsize-1],square=False, cbar_kws={'label': f"coverage pos selection c3", "pad": 0.02}, vmin=0, yticklabels=False, xticklabels=False, vmax=500, cbar=show_cbar_for_each, cmap=coverage_cmap)
        
    # --- Add bias heatmap ---
    if bias_per_pos:
        spec_cmap = sns.light_palette("black", n_colors=30, reverse=False, as_cmap=True)
        vmin_bias = min(bias_per_pos)
        vmax_bias = max(bias_per_pos)

        if show_only_pos:
            bias_per_pos = [bias_per_pos[pos] for pos in show_only_pos] # Filter to position of interest


        sns.heatmap(pd.DataFrame(bias_per_pos).T, ax=axes[pltsize-2], cmap=spec_cmap, square=False, cbar_kws={'label': "chance of codon mutation", "pad": 0.02}, yticklabels=False, xticklabels=False, cbar=show_cbar_for_each, vmin=vmin_bias, vmax=vmax_bias)

    # --- Add global colorbars if needed ---
    if not show_cbar_for_each:
        if not vmax: 
            print("Please set vmax to show one colorbar for all")
            exit()
        cbar_ax = fig.add_axes([0.13, 0.05, 0.15, 0.02])
        cbar = fig.colorbar(axes[0].collections[0], cax=cbar_ax, orientation = "horizontal")
        cbar.set_label("mutation rate", fontsize = 15)
        cbar.ax.tick_params(labelsize=10)
        if plot_coverage:
            cbar_ax = fig.add_axes([0.32, 0.05, 0.15, 0.02])
            cbar = fig.colorbar(axes[pltsize-1].collections[0], cax=cbar_ax, orientation = "horizontal")
            cbar.set_label('read depth', fontsize = 15)
            cbar.ax.tick_params(labelsize=10)

        if bias_per_pos: 
            cbar_ax = fig.add_axes([0.51, 0.05, 0.15, 0.02])
            cbar = fig.colorbar(axes[pltsize-2].collections[0], cax=cbar_ax, orientation = "horizontal")
            cbar.set_label('chance of codon mutation', fontsize = 15)
            cbar.ax.tick_params(labelsize=10)

    # if FigFolder:
    #     if not os.path.exists(FigFolder):
    #         os.makedirs(FigFolder)
        
    # dattype = data_type if not combine_mut_rates else "combined" + data_type
    # name = f"PANCE_mut_enrichment_all_sections_{dattype}.pdf" if not show_only_pos else f"PANCE_mut_enrichment_all_sections_high_mut_pos_{dattype}.pdf"
    # plt.savefig(os.path.join(FigFolder, name), dpi=300)

    if return_df:
        if show_only_pos:
            combined_df = combined_df.iloc[:,show_only_pos]
            if ref_annot:
                combined_df.columns = [ref_annot[i] for i in show_only_pos]
        return combined_df


def plot_indel_freqs(indels, filename, FigFolder = None, roi_start_idx = None, roi_end_idx = None, color1= "#C7F9CC", color2 = "#38A3A5"):
    """
    Plot insertion and deletion frequencies across positions.

    Parameters:
    - indels (pd.DataFrame): DataFrame with rows ['insertion', 'deletion'] and columns = positions.
    - filename (str): Sample name for title.
    - FigFolder (str): Optional save path.
    - roi_start_idx (int): Region-of-interest start index (vertical line).
    - roi_end_idx (int): Region-of-interest end index (vertical line).
    - color1 (str): Color for insertions.
    - color2 (str): Color for deletions.
    """
    fig, axes = plt.subplots(1, figsize=(15,5))
    plt.plot(indels.columns, indels.loc["insertion",:], label = "insertion", color = color1)
    plt.plot( indels.columns, indels.loc["deletion",:], label = "deletion",alpha = 0.5, color = color2)

    if roi_start_idx:
        plt.axvline(x=roi_start_idx, color='grey', linestyle='--', label = "insertion site")
    if roi_end_idx:
        plt.axvline(x=roi_end_idx, color='grey', linestyle='--')
    plt.legend(frameon = False)
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    plt.title(f"Indel frequency {filename}")
    # if FigFolder:
    #     plt.savefig(f"{FigFolder}/{filename}_indel_freq.pdf", bbox_inches='tight')
    #     plt.savefig(f"{FigFolder}/{filename}_indel_freq.png", bbox_inches='tight')
    plt.show()
    plt.close()