import os
import argparse
import subprocess

# --- Use chopper for quality filtering and filtering based on read length ---
if __name__ == "__main__":
    # --- Set up argument parser for input and output folder paths ---
    parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
    parser.add_argument("input", type=str, help="Path to the folder storing the input bam files")
    parser.add_argument("output", type=str, help="Folder path to store the output")

    args = parser.parse_args()

    # --- Assign input and output paths to shorter variable names ---
    input_folder = args.input
    output_folder = args.output

    # --- Validate that input folder exists ---
    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"Input folder '{input_folder}' does not exist!")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print(f"running Nanoplot on {input_folder}...")

    # --- Collect all .fastq.gz files in the input folder ---
    input_files =  [f for f in os.listdir(input_folder) if f.endswith(".fastq.gz")]

    # --- Process each input file using chopper ---
    for input_file in  input_files:
        print(f"Processing {input_file}")
        
        # --- Define the full path for the output file ---
        output_file = f"{output_folder}/{input_file}"

        # --- Build the chopper command:
        #     - `-q 20`: keep only reads with a minimum average quality score of 20 (Phred scale)
        #     - `--minlength 1800`: filter out reads shorter than 1800 bases
        #     - `--maxlength 2200`: filter out reads longer than 2200 bases
        #     - `-i`: specify the input FASTQ file
        #     - Output is piped to `gzip`
                
        command = f"chopper -q 20 --minlength 1800 --maxlength 2200  -i {input_folder}/{input_file} | gzip > {output_file}" 
        print(f"Processing {input_file} -> {output_file}")

        # --- Run the command in the shell and raise error on failure ---
        subprocess.run(command, shell=True, check=True)

    print("Processing complete!")
