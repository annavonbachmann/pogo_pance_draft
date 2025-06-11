#!/bin/bash

# Please make sure to download the input dataset `Nanopore` from the Zenodo repository and place it in a folder named `data` at root repository level.

script_path="$(realpath "$0")"
script_dir="$(dirname "$script_path")"
repo_root="$(realpath "$script_dir/../../")"

echo "Using repo root: $repo_root"

# Path to all barcode folders
barcode_parent="$repo_root/data/Nanopore_P0109/basecalling/pass"

# Loop over experiment names inside 'pass'
for barcode_path in "$barcode_parent"/*
do
    barcode="$(basename "$barcode_path")"

    # Define input/output paths using real barcode name
    raw_fastq_input_folder="$barcode_path"
    output_folder="$repo_root/data/Nanopore_P0109/$barcode/highly_accurate_basecalling/filtered_Q20_maxminlen"
    output_plot_folder="$repo_root/final_output/Nanopore_P0109/$barcode/highly_accurate_basecalling/filtered_Q20_maxminlen/quality_control"
    reference_file="$repo_root/data/Nanopore_P0109/AraC_S170_LOV_R5_ref.fa"

    echo "################## Running filtering for $barcode... ##################"
    python3 "$repo_root/analysis_pipeline/Nanopore/01_Nanopore_read_filtering.py" "$raw_fastq_input_folder" "$output_folder/filtered_fastq"

    echo "################## Running alignment for $barcode... ##################"
    python3 "$repo_root/analysis_pipeline/Nanopore/02_Nanopore_alignment.py" "$output_folder/filtered_fastq" "$output_folder/minimap2_alignment" "$reference_file"

    echo "################## Plotting quality plots for $barcode... ##################"
    python3 "$repo_root/analysis_pipeline/Nanopore/03_Nanopore_quality_control.py" "$output_folder/minimap2_alignment" "$output_plot_folder"

    ## IMPORTANT: when you are analyzing linkers, you can skip the next step, since we are using the minimap2 alignments for further analysis (instead of the reads that are enforced to be in-frame (below), because for linkers, we encode indels in the library and thus expect them to be present in the reads)

    echo "################## Running read processing for $barcode... ##################"
    python3 "$repo_root/analysis_pipeline/Nanopore/04_Nanopore_process_reads.py" "$output_folder/minimap2_alignment" "$output_folder/processed_reads" "$reference_file"

    echo "Pipeline finished for $barcode!"
done