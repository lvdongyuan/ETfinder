import os
from recTid_gbac_7_0 import extract_sequence_ids, find_gbff_with_target_ids, extract_16S_rRNA_from_gb_folders
from group16s_7_1 import group_sequences_by_gca_number, process_grouped_sequences





def main():
    # 动态获取脚本目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 定义全局输出目录
    global_output_dir = os.path.join(script_dir, "temp_results")

    if not os.path.exists(global_output_dir):
        os.makedirs(global_output_dir)

    # 定义文件路径
    duprecT_file = os.path.join(global_output_dir, "paired_recTseq_95.fasta")
    duprecTid_file = os.path.join(global_output_dir, "recT95ids.txt")
    accession_file = os.path.join(global_output_dir, "recT95ac_list.txt")
    gbff_directory = "/Volumes/ex/recETresults/lddeipffaa"
    recT16s_output = os.path.join(global_output_dir, "recT95_16S.txt")
    missing_output = os.path.join(global_output_dir, "missing_16S.txt")
    output_file_uni = os.path.join(global_output_dir, "uni16s.txt")
    output_file_man = os.path.join(global_output_dir, "mannualdup16s.txt")

  
    
    # 第一步：从FASTA文件中提取ID
    extract_sequence_ids(duprecT_file, duprecTid_file)

    # 第二步：根据ID文件在GBFF文件中查找目标ID
    find_gbff_with_target_ids(gbff_directory, duprecTid_file, accession_file)

    # 第三步：根据accession和protein_id文件提取16S rRNA序列
    extract_16S_rRNA_from_gb_folders(accession_file, gbff_directory, recT16s_output, missing_output)

    # 第四步：根据16S rRNA序列文件进行分组和处理
    grouped_sequences = group_sequences_by_gca_number(recT16s_output)
    process_grouped_sequences(grouped_sequences, output_file_uni, output_file_man)




if __name__ == "__main__":
    main()