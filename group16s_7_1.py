from Bio import SeqIO
from Bio.Align import PairwiseAligner
from collections import defaultdict

def group_sequences_by_gca_number(fasta_file):
    grouped_sequences = defaultdict(list)
    
    # 读取FASTA文件并分组
    for record in SeqIO.parse(fasta_file, "fasta"):
        # 获取描述信息中的倒数第二个下划线后的部分到结尾
        seq_description_parts = record.description.rsplit('_', 2)
        if len(seq_description_parts) >= 2:
            gca_number = seq_description_parts[-2] + "_" + seq_description_parts[-1]
            print(f"Record description: {record.description}, GCA Number: {gca_number}")  # 调试信息
            grouped_sequences[gca_number].append(record)
        else:
            print(f"Invalid description format: {record.description}")
    
    # 打印分组信息
    for gca_number, sequences in grouped_sequences.items():
        print(f"GCA Number: {gca_number}, Number of sequences: {len(sequences)}")
    
    return grouped_sequences


def calculate_similarity(seq1, seq2):
    aligner = PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    similarity = best_alignment.score / max(len(seq1), len(seq2))
    return similarity

def process_grouped_sequences(grouped_sequences, output_file_one, output_file_two):
    sequences_to_write_one = []
    sequences_to_write_two = []
    
    for gca_number, sequences in grouped_sequences.items():
        if len(sequences) == 1:
            sequences_to_write_one.append(sequences[0])
        elif len(sequences) == 2:
            longest_seq = max(sequences, key=lambda seq: len(seq.seq))
            sequences_to_write_one.append(longest_seq)
        elif len(sequences) > 2:
            sorted_sequences = sorted(sequences, key=lambda seq: len(seq.seq), reverse=True)
            top2_sequences = sorted_sequences[:2]
            similarity = calculate_similarity(str(top2_sequences[0].seq), str(top2_sequences[1].seq))
            if similarity > 0.998:
                sequences_to_write_one.append(top2_sequences[0])
            else:
                sequences_to_write_two.extend(top2_sequences)
    
    # 写入输出文件
    SeqIO.write(sequences_to_write_one, output_file_one, "fasta")
    SeqIO.write(sequences_to_write_two, output_file_two, "fasta")

def main(input_fasta, output_file_one, output_file_two):
    # 分组序列
    grouped_sequences = group_sequences_by_gca_number(input_fasta)
    
    # 处理并写入新文件
    process_grouped_sequences(grouped_sequences, output_file_one, output_file_two)

if __name__ == "__main__":
    input_fasta = "/Volumes/ex/recETresults/recT_16S.txt"
    output_file_one = "/Volumes/ex/recETresults/single_and_longest_16s_rRNA_gca.txt"
    output_file_two = "/Volumes/ex/recETresults/top2_longest_16s_rRNA_gca.txt"

    main(input_fasta, output_file_one, output_file_two)

    print(f"Processed sequences have been written to {output_file_one} and {output_file_two}")
