from Bio import SeqIO
import os

def read_ids_from_file(file_path):
    """从给定的文件中读取ID，并存储在集合中返回"""
    with open(file_path, 'r') as file:
        ids = {line.strip() for line in file if line.strip()}
    return ids

def find_gbff_with_target_ids(root_dir, id_file, output_file):
    target_filename = "genomic.gbff"  # 指定查找的文件名
    ids = read_ids_from_file(id_file)  # 从文件读取ID列表

    with open(output_file, 'w') as f:
        for root, dirs, files in os.walk(root_dir):
            for file in files:
                if file.endswith(target_filename) and not file.startswith('._'):  # 确保处理的是有效的GBFF文件
                    file_path = os.path.join(root, file)
                    print(f"Processing file: {file_path}")
                    try:
                        for record in SeqIO.parse(file_path, "genbank"):
                            for feature in record.features:
                                if feature.type == "CDS" and 'protein_id' in feature.qualifiers:
                                    matched_ids = [pid for pid in feature.qualifiers['protein_id'] if pid in ids]
                                    if matched_ids:
                                        accession = os.path.basename(os.path.dirname(file_path))  # 获取包含GBFF文件的目录名
                                        for pid in matched_ids:
                                            f.write(f"{accession} {pid}\n")  # 写入accession和protein_id
                                        print(f"Found target IDs in {accession} for IDs: {', '.join(matched_ids)}")
                                        break
                    except Exception as e:
                        print(f"Error processing {file_path}: {e}")

    print(f"All matching accessions have been saved to {output_file}")


root_directory = "/Volumes/ex/recETresults/lddeipffaa"  # GBFF文件所在的根目录
id_file_path = "/Volumes/ex/recETresults/demult_paired_recT_ids.txt"  # 包含protein_id的文本文件路径
output_accession_file = "/Volumes/ex/recETresults/accession_list.txt"  # 输出文件路径

find_gbff_with_target_ids(root_directory, id_file_path, output_accession_file)
