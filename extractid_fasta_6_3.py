from Bio import SeqIO

def extract_sequence_ids(fasta_file, output_file):
    with open(fasta_file, 'r') as file, open(output_file, 'w') as out:
        for record in SeqIO.parse(file, "fasta"):
            # 提取从大于号后到第一个下划线前的ID部分
            sequence_id = record.id.split('_')[0]
            # 输出ID到文件
            out.write(sequence_id + '\n')

