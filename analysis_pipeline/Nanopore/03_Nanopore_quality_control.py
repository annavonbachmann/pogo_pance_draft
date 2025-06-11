import os
import argparse
import subprocess

# --- Use NanoPlot for visualization of read quality and alignment quality ---

if __name__ == "__main__":
    # --- Set up argument parser for input and output folder paths ---
    parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
    parser.add_argument("input", type=str, help="Path to the folder storing the input bam files")
    parser.add_argument("output", type=str, help="Folder path to store the output")

    args = parser.parse_args()

    # --- Assign input arguments to shorter variable names ---
    input_folder = args.input
    output_folder = args.output
    # --- Validate that input folder exists ---
    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"Input folder '{input_folder}' does not exist!")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    print(f"running Nanoplot on {input_folder}...")
    
    # --- Determine file type based on extension of the first file in the input folder ---
    filetype = '.fastq' if '.fastq.gz' in os.listdir(input_folder)[0] else '.bam'
    fileend = '.fastq.gz' if filetype == '.fastq' else '.bam'
    # --- Collect full paths of all relevant files in the input folder ---
    files =  [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(fileend)]

    # --- Run NanoPlot with the selected files ---
    #     - "-t 2": use 2 threads
    #     - "--fastq" or "--bam": specify input file type based on earlier logic
    #     - "-o": specify the output folder
    #     - "-f pdf": generate output in PDF format
    subprocess.run([
        "NanoPlot",  
        "-t", "2",
        f"--{filetype[1:]}", *files, 
        "-o", str(output_folder),
        "-f", "pdf",
    ], check=True)
