import os
from extractidfrompsi_5_1 import extract_unique_ids
from pssm_psiblast_5_0 import run_blast
from recTvsredB_5_2 import compare_ids
from sortETid_6_0 import process_ids
from blastdbcmd_fmt6_6_1 import extract_sequences
from cdhit_6_2 import run_cd_hit
from extractid_fasta_6_3 import extract_sequence_ids

def main():
    # 动态获取脚本目录
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 定义ETfinder父文件夹路径
    etfinder_parent_dir = script_dir

    # List of PSSM files to be processed
    pssm_files_with_params = {
    "pfam03837.pssm": {"-evalue": "0.05", "-num_iterations": "5"},
    "pfam09588.pssm": {"-evalue": "0.05", "-num_iterations": "6"},
    "prk09729.pssm": {"-evalue": "0.0005", "-num_iterations": "2"},
    "TIGR01913.pssm": {"-evalue": "0.05", "-num_iterations": "9"}
    # Add more PSSM files and their parameters here
}

    # Directory where the BLAST results will be saved (relative to the ETfinder directory)
    global_output_dir = os.path.join(etfinder_parent_dir, "temp_results")

    # Database to be used for BLAST (assuming it's in a known location or configured in BLAST environment)
    database = "/Volumes/ex/recETresults/lddeipfdb/lddeipf"

    # Run BLAST for each PSSM file
    run_blast(pssm_files_with_params, global_output_dir, database, etfinder_parent_dir)

    # Process each BLAST results file with extract_unique_ids
    
    # 定义文件名列表
    for pssm_file in pssm_files_with_params.keys():
        base_name = os.path.splitext(pssm_file)[0]
        blast_results_path = os.path.join(global_output_dir, f"{base_name}.txt")
        extract_unique_ids(blast_results_path, global_output_dir)
    # 指定输出目录
  
    file1_path1 = os.path.join(global_output_dir, "pfam03837_ids.txt")
    file1_path2 = os.path.join(global_output_dir, "TIGR01913_ids.txt")
    output1_path = os.path.join(global_output_dir, "ssap_ids.txt")
    compare_ids(file1_path1, file1_path2, output1_path)

    file2_path1 = os.path.join(global_output_dir, "pfam09588_ids.txt")
    file2_path2 = os.path.join(global_output_dir, "prk09729_ids.txt")
    output2_path = os.path.join(global_output_dir, "exo_ids.txt")
    compare_ids(file2_path1, file2_path2, output2_path)
    #process id
    process_ids(output1_path, output2_path, global_output_dir)

    id_list_path = os.path.join(global_output_dir, 'paired_recT_ids.txt')
    db_path = '/Volumes/ex/recETresults/lddeipfdb/lddeipf'  # Assuming this is the database path
    fasta_output_path = os.path.join(global_output_dir, 'paired_recTseq.fasta')
    extract_sequences(id_list_path, db_path, fasta_output_path)

    # Run CD-HIT on the extracted sequences
    cd_hit_output_path = os.path.join(global_output_dir, 'paired_recTseq_95.fasta')
    run_cd_hit(fasta_output_path, cd_hit_output_path, c=0.95)

    unique_recT_fasta = os.path.join(global_output_dir, 'paired_recTseq_95.fasta')
    recT_id_output = os.path.join(global_output_dir, 'recT95id.txt')
    extract_sequence_ids(unique_recT_fasta, recT_id_output)



if __name__ == "__main__":
    main()
