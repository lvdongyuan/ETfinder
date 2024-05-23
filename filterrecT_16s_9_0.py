from Bio import SeqIO
import os
import argparse


def extract_sequence_ids(fasta_file, id_list_file):
    """从FASTA文件中提取ID并写入ID列表文件"""
    with open(fasta_file, 'r') as file, open(id_list_file, 'w') as out:
        for record in SeqIO.parse(file, "fasta"):
            # 提取从大于号后到第一个下划线前的ID部分
            sequence_id = record.id.split('_')[0]
            # 输出ID到文件
            out.write(sequence_id + '\n')

def read_ids_from_file(id_list_file):
    """从文件中读取ID列表并存储在集合中返回"""
    with open(id_list_file, 'r') as file:
        ids = {line.strip() for line in file if line.strip()}
    return ids

def filter_sequences_by_id_list(target_fasta_file, id_list_file, output_file):
    """根据ID列表过滤目标FASTA文件，只保留ID列表中的序列"""
    ids = read_ids_from_file(id_list_file)
    with open(target_fasta_file, 'r') as fasta, open(output_file, 'w') as out:
        for record in SeqIO.parse(fasta, "fasta"):
            # 提取从大于号后到第一个下划线前的ID部分
            sequence_id = record.id.split('_')[0]
            if sequence_id in ids:
                SeqIO.write(record, out, "fasta")

def extract_and_filter_sequences(input_fasta_file, target_fasta_file, output_file):
    """从输入FASTA文件中提取ID并根据这些ID过滤目标FASTA文件"""
    id_list_file = "temp_id_list.txt"  # 临时ID列表文件
    extract_sequence_ids(input_fasta_file, id_list_file)
    filter_sequences_by_id_list(target_fasta_file, id_list_file, output_file)
    # 删除临时ID列表文件
    os.remove(id_list_file)


def main():
    parser = argparse.ArgumentParser(description="Extract IDs from input FASTA and filter target FASTA based on these IDs.")
    parser.add_argument("input_fasta_file", help="Path to the input FASTA file")
    parser.add_argument("target_fasta_file", help="Path to the target FASTA file")
    parser.add_argument("output_fasta_file", help="Path to the output FASTA file")

    args = parser.parse_args()

   
    extract_and_filter_sequences(args.input_fasta_file, args.target_fasta_file, args.output_fasta_file)

    print(f"Filtered sequences have been written to {args.output_fasta_file}")

if __name__ == "__main__":
    main()

    # 
    #input_fasta_file = "/Volumes/ex/recETresults/unique16s_rect.txt"
    #target_fasta_file = "/Volumes/ex/recETresults/demult_paired_recTseq.txt"
    #output_fasta_file = "/Volumes/ex/recETresults/testrect_16s.txt"

