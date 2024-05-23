import os

def extract_unique_ids(blast_results_path, output_dir):
    # Extract the base name of the blast results file
    base_name = os.path.splitext(os.path.basename(blast_results_path))[0]
    output_filename = f"{base_name}_ids.txt"
    output_file_path = os.path.join(output_dir, output_filename)

    # 初始化一个空集合来存储唯一的序列ID
    unique_ids = set()

    # 打开BLAST结果文件并读取每一行
    with open(blast_results_path, 'r') as file:
        for line in file:
            # 按'|'分割并提取序列ID（分割后的第二部分）
            parts = line.split('|')
            if len(parts) > 2:
                seq_id = parts[1].strip()  # 移除可能的空白字符
                unique_ids.add(seq_id)  # 添加到集合中去除重复

    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 将唯一的序列ID写入新的文件
    with open(output_file_path, 'w') as file:
        for seq_id in unique_ids:
            file.write(seq_id + '\n')

    print(f"Extracted {len(unique_ids)} unique sequence IDs from {blast_results_path} and saved to {output_file_path}.")
