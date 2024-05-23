import os
import subprocess

def run_blast(pssm_files, output_dir, database, etfinder_dir):
    # Get the current script directory dynamically

    # Ensure the output directory is also relative to the script directory

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for pssm_file, params in pssm_files.items():
        # Construct the full path to the PSSM file relative to the ETfinder directory
        pssm_path = os.path.join(etfinder_dir, "PSSMdb", pssm_file)

        # Extract the base name of the PSSM file to use as the output file name
        base_name = os.path.splitext(os.path.basename(pssm_file))[0]
        output_file = os.path.join(output_dir, f"{base_name}.txt")

        # Construct the BLAST command with the desired output format
        blast_command = [
            "psiblast",                   # Change to "blastp" or other BLAST tools as needed
            "-in_pssm", pssm_path,
            "-db", database,
            "-out", output_file,
            "-outfmt", "6 sseqid sseq"    # Change the output format as needed
        ]

        if params:
            for param, value in params.items():
                blast_command.extend([param, value])

        try:
            # Execute the BLAST command
            result = subprocess.run(blast_command, check=True, capture_output=True, text=True)
            print(f"BLAST completed for {pssm_file}. Results saved to {output_file}.")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running BLAST for {pssm_file}. Error message: {e.stderr}")

if __name__ == "__main__":
    # List of PSSM files to be processed
    pssm_files = [
        "pfam03837.pssm",
        "pfam09588.pssm",
        "prk09729.pssm",
        "TIGR01913.pssm"
        # Add more PSSM file names here
    ]

    # Directory where the BLAST results will be saved (relative to the script directory)
    output_dir = "temp_results"

    # Database to be used for BLAST (assuming it's in a known location or configured in BLAST environment)
    database = "your_database_name"

    # Run BLAST for each PSSM file
    run_blast(pssm_files, output_dir, database)
