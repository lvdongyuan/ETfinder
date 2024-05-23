from phiblastinlistdb_2_0 import blast_and_extract
from extractdbac_fmt6_2_1 import extract_database_ids
from mergehitfaa_3_0 import merge_faa_files
import os

def main():
       # Define the result directory relative to the current script directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    result_dir = os.path.join(current_dir, 'temp_results')
    
    # Ensure the result directory exists
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    # Step 1: Run PHI-BLAST and extract results
    db_directory = '/Users/lvdongyuan/Desktop/alphadb'   #Your path of Indivadral GCAdbs
    query_file = '/Users/lvdongyuan/Downloads/alphaproteobacteria/alphaproteodb/file/ssb.faa'    #Query your interesting SSB
    blast_results_file = os.path.join(result_dir, 'ssb_lddeipf.txt')
    no_hits_file = os.path.join(result_dir, 'no_hits_lddeipf.txt')
    patternfile = '/Users/lvdongyuan/Downloads/alphaproteobacteria/alphaproteodb/file/pattern_lddeipf.txt' #path to your pattern file
    max_hits = 1
    blast_and_extract(db_directory, query_file, blast_results_file, no_hits_file, patternfile, max_hits)

    # Step 2: Extract database identifiers
    input_data_path = blast_results_file
    output_path = os.path.join(result_dir, 'lddeipfdbac.txt')
    extract_database_ids(input_data_path, output_path)

    # Step 3: Merge .faa files
    file_list_path = output_path
    output_file_path = os.path.join(result_dir, 'lddeipf.faa')
    directory_path = '/Users/lvdongyuan/Downloads/alphaproteobacteria/gzdb'
    merge_faa_files(file_list_path, output_file_path, directory_path)

if __name__ == "__main__":
    main()


    current_dir = os.path.dirname(os.path.abspath(__file__))
    pssm_dir = os.path.join(current_dir, 'pssmdb')
    psiblast -db PRK09729 -in_pssm list -outfmt 6 sseqid sseq -output PSSM list.txt -num iteration x -evalue x 